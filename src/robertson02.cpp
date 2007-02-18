/**
 * @brief Robertson02 algorithm for automatic self-calibration.
 *
 * 
 * This file is a part of PFS CALIBRATION package.
 * ---------------------------------------------------------------------- 
 * Copyright (C) 2004 Grzegorz Krawczyk
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * ---------------------------------------------------------------------- 
 * 
 * @author Grzegorz Krawczyk, <gkrawczyk@users.sourceforge.net>
 *
 * $Id: robertson02.cpp,v 1.7 2006/11/16 15:06:17 gkrawczyk Exp $
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <vector>

#include <math.h>

#include <responses.h>
#include <robertson02.h>

#define PROG_NAME "robertson02"

// maximum iterations after algorithm accepts local minima
#define MAXIT 500

// maximum accepted error
#define MAX_DELTA 1e-5f

#define TEST_NEW_UCHAR
#ifdef  TEST_NEW_UCHAR
#define GET_AS_UCHAR(V) reinterpret_cast<unsigned char&>(V)
#define GET_AS_FLOAT(V) (float)reinterpret_cast<unsigned char&>(V)
#else
#define GET_AS_UCHAR(V) (V)
#define GET_AS_FLOAT(V) (V)
#endif

float normalizeI( float* I, int M );

int robertson02_applyResponse( pfs::Array2D* xj, const ExposureList &imgs,
  const float* I, const float* w, int M )
{
  // number of exposures
  int N = imgs.size();

  // frame size
  int width = xj->getCols();
  int height = xj->getRows();

  // number of saturated pixels
  int saturated_pixels = 0;

  // --- anti saturation: calculate trusted camera output range
  int minM = 0;
  for( int m=0 ; m<M ; m++ )
    if( w[m]>0 )
    {
      minM = m;
      break;
    }
  int maxM = M-1;
  for( int m=M-1 ; m>=0 ; m-- )
    if( w[m]>0 )
    {
      maxM = m;
      break;
    }
  
  // --- anti ghosting: for each image i, find images with
  // the immediately higher and lower exposure times
  int* i_lower = new int[N];
  int* i_upper = new int[N];
  for( int i=0 ; i<N ; i++ )
  {
    i_lower[i]=-1;
    i_upper[i]=-1;
    float ti = imgs[i].ti;
    float ti_upper = imgs[0].ti;
    float ti_lower = imgs[0].ti;

    for( int j=0 ; j<N ; j++ )
      if( i!=j )
      {
        if( imgs[j].ti>ti && imgs[j].ti<ti_upper )
        {
          ti_upper=imgs[j].ti;
          i_upper[i]=j;
        }
        if( imgs[j].ti<ti && imgs[j].ti>ti_lower )
        {
          ti_lower=imgs[j].ti;
          i_lower[i]=j;
        }
      }
    if( i_lower[i]==-1 )
      i_lower[i]=i;
    if( i_upper[i]==-1 )
      i_upper[i]=i;
  }
  
  
  // all pixels
  for( int j=0 ; j<width*height ; j++ )
  {
    // all exposures for each pixel
    float sum = 0.0f;
    float div = 0.0f;

    float maxti = -1e6f;
    float minti = +1e6f;
    
    for( int i=0 ; i<N ; i++ )
    {
      int m = (int) GET_AS_UCHAR((*imgs[i].yi)(j));
      float ti = imgs[i].ti;

      // --- anti saturation: observe minimum exposure time at which
      // saturated value is present, and maximum exp time at which
      // black value is present
      if( m>maxM )
        minti = fminf(minti,ti);
      if( m<minM )
        maxti = fmaxf(maxti,ti);
      
      // --- anti ghosting: monotonous increase in time should result
      // in monotonous increase in intensity; make forward and
      // backward check, ignore value if condition not satisfied
      int m_lower = (int) GET_AS_UCHAR((*imgs[i_lower[i]].yi)(j));
      int m_upper = (int) GET_AS_UCHAR((*imgs[i_upper[i]].yi)(j));
      if( m_lower>m || m_upper<m)
        continue;
      
      sum += w[m] * ti * I[m];
      div += w[m] * ti * ti;
    }

    // --- anti saturation: if a meaningful representation of pixel
    // was not found, replace it with information from observed data
    if( div==0.0f )
      saturated_pixels++;
    if( div==0.0f && maxti>-1e6f )
    {
      sum = I[minM];
      div = maxti;
    }
    if( div==0.0f && minti<+1e6f )
    {
      sum = I[maxM];
      div = minti;
    }
      
    if( div!=0.0f )
      (*xj)(j) = sum/div;
    else
      (*xj)(j) = 0.0f;
  }

  delete[] i_lower;
  delete[] i_upper;
  
  return saturated_pixels;
}


int robertson02_getResponse( pfs::Array2D* xj, const ExposureList &imgs,
  float* I, const float* w, int M )
{
  // number of exposures
  int N = imgs.size();

  // frame size
  int width = imgs[0].yi->getCols();
  int height = imgs[0].yi->getRows();

  // number of saturated pixels
  int saturated_pixels = 0;

  // indexes
  int i,j,m;

  float* Ip = new float[M];	// previous response
  if( Ip==NULL )
  {
    std::cerr << "robertson02: could not allocate memory for camera response" << std::endl;
    exit(1);
  }

  // 0. Initialization
  normalizeI( I, M );
  for( m=0 ; m<M ; m++ )
    Ip[m] = I[m];

  robertson02_applyResponse( xj, imgs, I, w, M );

  // Optimization process
  bool converged=false;
  long* cardEm = new long[M];
  float* sum = new float[M];
  if( sum==NULL || cardEm==NULL )
  {
    std::cerr << "robertson02: could not allocate memory for optimization process" << std::endl;
    exit(1);
  }

  int cur_it = 0;
  float pdelta= 0.0f;
  while( !converged )
  {
    // 1. Minimize with respect to I
    for( m=0 ; m<M ; m++ )
    {
      cardEm[m]=0;
      sum[m]=0.0f;
    }

    for( i=0 ; i<N ; i++ )
    {
      pfs::Array2D* yi = imgs[i].yi;
      float ti = imgs[i].ti;
      for( j=0 ; j<width*height ; j++ )
      {
	m = (int) GET_AS_UCHAR((*yi)(j));
	if( m<M && m>=0 )
	{
	  sum[m] += ti * (*xj)(j);
	  cardEm[m]++;
	}
  	else
  	  std::cerr << "robertson02: m out of range: " << m << std::endl;
      }
    }

    for( m=0 ; m<M ; m++ )
      if( cardEm[m]!=0 )
	I[m] = sum[m] / cardEm[m];
      else
	I[m] = 0.0f;
    
    // 2. Normalize I
    float middle_response = normalizeI( I, M );

    // 3. Apply new response
    saturated_pixels = robertson02_applyResponse( xj, imgs, I, w, M );

    // 4. Check stopping condition
    float delta = 0.0f;
    int hits=0;
    for( m=0 ; m<M ; m++ )
      if( I[m]!=0.0f )
      {
	float diff = I[m]-Ip[m];;
	delta += diff * diff;
	Ip[m] = I[m];
	hits++;
      }
    delta /= hits;

    VERBOSE_STR << " #" << cur_it
                << " delta=" << delta
                << " (coverage: " << 100*hits/M << "%)\n";
    
    if( delta < MAX_DELTA )
      converged=true;
    else if( isnan(delta) || (cur_it>MAXIT && pdelta<delta) )
    {
      VERBOSE_STR << "algorithm failed to converge, too noisy data in range\n";
      break;
    }

    pdelta = delta;
    cur_it++;
  }

  if( converged )
    VERBOSE_STR << " #" << cur_it
                << " delta=" << pdelta << " <- converged\n";

  delete[] Ip;
  delete[] cardEm;
  delete[] sum;
  
  return saturated_pixels;
}


//----------------------------------------------------------
// private part

float normalizeI( float* I, int M )
{
  int Mmin, Mmax;
  // find min max
  for( Mmin=0 ; Mmin<M && I[Mmin]==0 ; Mmin++ );
  for( Mmax=M-1 ; Mmax>0 && I[Mmax]==0 ; Mmax-- );
  
  int Mmid = Mmin+(Mmax-Mmin)/2;
  float mid = I[Mmid];

//   std::cerr << "robertson02: middle response, mid=" << mid
//             << " [" << Mmid << "]"
//             << " " << Mmin << ".." << Mmax << std::endl;
  
  if( mid==0.0f )
  {
    // find first non-zero middle response
    while( Mmid<Mmax && I[Mmid]==0.0f )
      Mmid++;
    mid = I[Mmid];
  }

  if( mid!=0.0f )
    for( int m=0 ; m<M ; m++ )
      I[m] /= mid;
  return mid;
}

