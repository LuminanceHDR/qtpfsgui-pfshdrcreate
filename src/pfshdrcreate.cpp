/**
 * @brief Create an HDR image or calibrate a response curve from a set
 * of differently exposed images supplied in PFS stream
 *
 * 
 * This file is a part of PFS CALIBRATION package.
 * ---------------------------------------------------------------------- 
 * Copyright (C) 2006 Giuseppe Rota
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
 * @author Giuseppe Rota, <grota@users.sourceforge.net>
 * @author Grzegorz Krawczyk, <krawczyk@mpi-sb.mpg.de>
 *
 * $Id: pfshdrcalibrate.cpp,v 1.10 2006/09/13 11:52:56 gkrawczyk Exp $
 */

#include <config.h>

#include <iostream>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include <pfs.h>

#include "responses.h"
#include <robertson02.h>
#include <icip06.h>

using namespace std;

#define PROG_NAME "pfshdrcreate"

#define TEST_NEW_UCHAR
#ifdef  TEST_NEW_UCHAR
#define GET_AS_UCHAR(V) reinterpret_cast<unsigned char&>(V)
#define GET_AS_FLOAT(V) (float)reinterpret_cast<unsigned char&>(V)
#else
#define GET_AS_UCHAR(V) (V)
#define GET_AS_FLOAT(V) (V)
#endif

//This class is an array2d that stores uchars and not floats,
//that reduces the memory footprint by a factor of 4.
namespace pfs {
class Array2DImplUCHAR: public Array2D
{
	unsigned char *data;
	int cols, rows;
	public:
	Array2DImplUCHAR( int cols, int rows ) : cols( cols ), rows( rows )
	{
		data = new unsigned char[cols*rows];
	}
	~Array2DImplUCHAR()
	{
		delete[] data;
	}
	inline int getCols() const { return cols; }
	inline int getRows() const { return rows; }
	
	inline float& operator()( int col, int row ) {
		assert( col >= 0 && col < cols );
		assert( row >= 0 && row < rows );
		return (float &)data[ col+row*cols ];
	}
	inline const float& operator()( int col, int row ) const {
		assert( col >= 0 && col < cols );
		assert( row >= 0 && row < rows );
		return (float &)data[ col+row*cols ];
	}
	inline float& operator()( int index ) {
		assert( index >= 0 && index < rows*cols );
		return (float &)data[index];
	}
	inline const float& operator()( int index ) const {
		assert( index >= 0 && index <= rows*cols );
		return (float &)data[index];
	}
	unsigned char* getRawData() {
		return data;
	}
};
}

inline float max3( float a, float b, float c )
{
  float max = (a>b) ? a : b;
  return (c>max) ? c : max;
}

inline float min3( float a, float b, float c )
{
  // ignore zero values
  if( int(a)==0 ) a=1e8;
  if( int(b)==0 ) b=1e8;
  if( int(c)==0 ) c=1e8;
  
  float min = (a<b) ? a : b;
  return (c<min) ? c : min;
}

//---------------------------------------------------
//--- standard PFS stuff
bool verbose = false;

class QuietException 
{
};


void printHelp()
{
  fprintf( stderr, PROG_NAME ": \n"
    "\t[--weights triangular|gaussian|plateau]\n"
    "\t[--response linear|gamma|log|calibrate] [--response-file <filename.m>] \n"
    "\t[--model robertson|debevec]\n"
    "\t[--save-response <filename.m>] \n"
    "\t[--bpp <val>] \n"
    "\t[--gauss <val>] (meaningful only if you select --weights gaussian)\n"
    "\t[--verbose] [--help]\n"
    "See man page for more information.\n" );
}

void pfshdrcalibrate( int argc, char* argv[] )
{
  pfs::DOMIO pfsio;

  enum TWeight
    { TRIANGULAR, GAUSSIAN, PLATEAU } opt_weight = TRIANGULAR;
  enum TResponse
    { FROM_FILE, LINEAR, GAMMA, LOG10, FROM_ROBERTSON } opt_response = GAMMA;
  enum TModel
    { ROBERTSON, DEBEVEC } opt_model = DEBEVEC;

  float input_multiplier = 1.0f;
  FILE* responseFile = NULL;
  FILE* responseSaveFile = NULL;
  int opt_bpp = 8;
  float opt_gauss = 8.0f;
  int opt_maxresponse = -1;
  int opt_minresponse = -1;

  static struct option cmdLineOptions[] = {
    { "weights", required_argument, NULL, 'w' },  //WEIGHTS PREFERENCES
    { "response", required_argument, NULL, 'r' }, //RESPONSE CURVES PREFERENCES
    { "model", required_argument, NULL, 'm' },    //HDR-GENERATION MODEL PREFERECES
    { "gauss", required_argument, NULL, 'g' },
    { "max-response", required_argument, NULL, 'A' },
    { "min-response", required_argument, NULL, 'S' },
    { "response-file", required_argument, NULL, 'f' },
    { "save-response", required_argument, NULL, 's' },
    { "bpp", required_argument, NULL, 'b' },
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
  };

  int optionIndex = 0;
  while( 1 ) {
    int c = getopt_long (argc, argv, "hvw:r:m:g:f:s:b:", cmdLineOptions, &optionIndex);
    if( c == -1 ) break;
    switch( c ) {
    case 'h':
      printHelp();
      throw QuietException();
    case 'v':
      verbose = true;
      break;
    case 'w': //WEIGHTS PREFERENCES
      if( strcmp(optarg,"triangular")==0 )
        opt_weight = TRIANGULAR;
      else if( strcmp(optarg,"gaussian")==0 )
        opt_weight = GAUSSIAN;
      else if( strcmp(optarg,"plateau")==0 )
        opt_weight = PLATEAU;
      else
        throw pfs::Exception("unknown model");
      break;
    case 'r': //RESPONSE CURVES PREFERENCES
      if( strcmp(optarg,"linear")==0 )
        opt_response = LINEAR;
      else if( strcmp(optarg,"gamma")==0 )
        opt_response = GAMMA;
      else if( strcmp(optarg,"log")==0 )
        opt_response = LOG10;
      else if( strcmp(optarg,"calibrate")==0 )
        opt_response = FROM_ROBERTSON;
      else
        throw pfs::Exception("unknown standard response (check the manpage or use default)");
      break;
    case 'm': //HDR-GENERATION MODEL PREFERECES
      if( strcmp(optarg,"robertson")==0 )
        opt_model = ROBERTSON;
      else if( strcmp(optarg,"debevec")==0 )
        opt_model = DEBEVEC;
      else
        throw pfs::Exception("unknown model");
      break;
    case 'g':
      opt_gauss = atof(optarg);
      if( opt_gauss<=0.0f or opt_gauss>1024.0f)
        throw pfs::Exception("sigma value for Gaussian out of range. accepted range 0:1024");
      break;
    case 'f':
      opt_response = FROM_FILE;
      responseFile = fopen(optarg, "r");
      if( !responseFile )
        throw pfs::Exception("could not open file with response curve");
      break;
    case 's':
      responseSaveFile = fopen(optarg,"w");
      if( !responseSaveFile )
        throw pfs::Exception("could not open file to save response curve");
      break;
    case 'b':
      opt_bpp = atoi(optarg);
      if( opt_bpp<8 || opt_bpp>32)
        throw pfs::Exception("bits per pixel value out of range. accepted range >=8");
      break;
    case 'A':                   // max response
      opt_maxresponse = atoi(optarg);
      if( opt_maxresponse<=opt_minresponse )
        throw pfs::Exception("max response should be higher than min response");
      break;
    case 'S':                   // min response
      opt_minresponse = atoi(optarg);
      if( opt_minresponse<0 )
        throw pfs::Exception("min response should be >0");
      break;
    case '?':
      throw QuietException();
    case ':':
      throw QuietException();
    }
  }

  // in PFS streams, 8bit data are mapped to 0:1 range
  if( opt_bpp == 8 )
    input_multiplier = 255.0f;

//PRINT WEIGHTS PREFERENCES
  switch( opt_weight )
  {
  case TRIANGULAR:
    VERBOSE_STR << "using triangular-shaped weights" << std::endl;
    break;
  case GAUSSIAN:
    VERBOSE_STR << "using gaussian-shaped weights" << std::endl;
    break;
  case PLATEAU:
    VERBOSE_STR << "using plateau-shaped weights" << endl;
    break;
  default:
    throw pfs::Exception("weights not recognized");
  }

//PRINT RESPONSE CURVES PREFERENCES
  switch( opt_response )
  {
  case FROM_FILE:
    VERBOSE_STR << "response curve from file" << endl;
    break;
  case LINEAR:
    VERBOSE_STR << "initial response: linear" << endl;
    break;
  case GAMMA:
    VERBOSE_STR << "initial response: gamma" << endl;
    break;
  case LOG10:
    VERBOSE_STR << "initial response: logarithmic" << endl;
    break;
  case FROM_ROBERTSON:
    VERBOSE_STR << "initial response to be found from image set" << endl;
    break;
  default:
    throw pfs::Exception("undefined standard response");
    break;
  }

//PRINT HDR-GENERATION MODEL PREFERENCES
  switch( opt_model )
  {
  case ROBERTSON:
    VERBOSE_STR << "using robertson model" << std::endl;
    break;
  case DEBEVEC:
    VERBOSE_STR << "using debevec model" << std::endl;
    break;
  default:
    throw pfs::Exception("hdr generation method not set or not supported");
  }

  if( responseSaveFile!=NULL )
    VERBOSE_STR << "saving response curve to a file" << endl;

  // number of input levels
  int M = (int) powf(2.0f,opt_bpp);
  VERBOSE_STR << "number of input levels: " << M << endl;
  VERBOSE_STR << "input multiplier: " << input_multiplier << endl;

  //--- read frames from pfs stream
  int frame_no = 1;
  int width=0, height=0, size=0;
  float minResponse = M;
  float maxResponse = 0.0f;

  // collected exposures
  ExposureList imgsR;
  ExposureList imgsG;
  ExposureList imgsB;

  while( true ) {
    pfs::Frame *framein = pfsio.readFrame( stdin );
    if( framein == NULL ) break; // No more frames

    pfs::Channel *X, *Y, *Z;
    framein->getXYZChannels( X, Y, Z );
    if( X==NULL || Y==NULL || Z==NULL )
      throw pfs::Exception( "missing XYZ channels in the PFS stream (try to preview your files using pfsview)" );

    const char* exp_str = framein->getTags()->getString("BV");
    if( exp_str==NULL )
      throw pfs::Exception( "missing exposure information in the PFS stream (use pfsinhdrgen to input files)" );

    // relate APEX brightness value only as a function of exposure time
    // that is assume aperture=1 and sensitivity=1
    float exp_time = 1.0f / powf(2.0f,atof( exp_str ));

    // absolute calibration: this magic multiplier is a result of my
    // research in internet plus a bit of tweaking with a luminance
    // meter. tested with canon 350d, so might be a bit off for other
    // cameras. to control absolute calibration modify iso values in
    // hdrgen script or use pfsabsolute program.
    exp_time /= 1.0592f * 11.4f / 3.125f;

    // frame size
    width = Y->getCols();
    height = Y->getRows();
    size = width*height;

    // new exposure image
    // collect RGB channels if not in luminance mode or we will apply
    // response to image (not in save response mode)
        pfs::transformColorSpace( pfs::CS_XYZ, X, Y, Z, pfs::CS_RGB, X, Y, Z );
        Exposure eR,eG,eB;
        eR.ti = eG.ti = eB.ti = exp_time;
#ifdef  TEST_NEW_UCHAR
        eR.yi = new pfs::Array2DImplUCHAR(width,height);
        eG.yi = new pfs::Array2DImplUCHAR(width,height);
        eB.yi = new pfs::Array2DImplUCHAR(width,height);
#else
        eR.yi = new pfs::Array2DImpl(width,height);
        eG.yi = new pfs::Array2DImpl(width,height);
        eB.yi = new pfs::Array2DImpl(width,height);
#endif
        if( eR.yi==NULL || eG.yi==NULL || eB.yi==NULL )
          throw pfs::Exception( "could not allocate memory for source exposure" );

        for( int i=0 ; i<size ; i++ )
        {
#ifdef TEST_NEW_UCHAR
          GET_AS_UCHAR((*eR.yi)(i)) = (int)((*X)(i) * input_multiplier);
          GET_AS_UCHAR((*eG.yi)(i)) = (int)((*Y)(i) * input_multiplier);
          GET_AS_UCHAR((*eB.yi)(i)) = (int)((*Z)(i) * input_multiplier);
          float maxval = max3(GET_AS_FLOAT((*eR.yi)(i)),GET_AS_FLOAT((*eG.yi)(i)),GET_AS_FLOAT((*eB.yi)(i)));
          float minval = min3(GET_AS_FLOAT((*eR.yi)(i)),GET_AS_FLOAT((*eG.yi)(i)),GET_AS_FLOAT((*eB.yi)(i)));
#else
          (*eR.yi)(i) = (*X)(i) * input_multiplier;
          (*eG.yi)(i) = (*Y)(i) * input_multiplier;
          (*eB.yi)(i) = (*Z)(i) * input_multiplier;
          float maxval = max3((*eR.yi)(i),(*eG.yi)(i),(*eB.yi)(i));
          float minval = min3((*eR.yi)(i),(*eG.yi)(i),(*eB.yi)(i));
#endif
          maxResponse = (maxResponse>maxval) ? maxResponse : maxval;
          minResponse = (minResponse<minval) ? minResponse : minval;
        }
        // add to exposures list
        imgsR.push_back(eR);
        imgsG.push_back(eG);
        imgsB.push_back(eB);

    VERBOSE_STR << "frame #" << frame_no << ", BV=" << atof(exp_str) << endl;
    frame_no++;

    // we don't need the frames that we got from stdin anymore
    pfsio.freeFrame( framein );

  } // end of while, we finished reading XYZ frames from stdin

  if( frame_no<2 )
    throw pfs::Exception( "at least one image required for calibration (check paths in hdrgen script?)" );

  // some more info on input frames
  VERBOSE_STR << "registered values: min=" << (int) minResponse
              << " max=" << (int) maxResponse << endl;
  if( opt_minresponse==-1 )
    opt_minresponse = (int) minResponse;
  if( opt_maxresponse==-1 )
    opt_maxresponse = (int) maxResponse;

  if( opt_response!=FROM_FILE )
    VERBOSE_STR << "camera response range: min=" << opt_minresponse
                << " max=" << opt_maxresponse << endl;

  if( maxResponse >= M )
    throw pfs::Exception( "input value higher than defined number of input levels (adjust the number of bits per pixel)" );

  // weighting function representing confidence in accuracy of acquisition
  float* w = new float[M];
  if( w==NULL )
    throw pfs::Exception( "could not allocate memory for weighting function" );
  //1) now fill in w based on weights preferences
  switch( opt_weight )
  {
  case TRIANGULAR:
    weights_triangle(w, M);
    break;
  case GAUSSIAN:
    weightsGauss( w, M, opt_minresponse, opt_maxresponse, opt_gauss );
    break;
  case PLATEAU:
    exposure_weights_icip06(w, M);
    break;
  }

  // create channels for output
  pfs::Frame *frameout = pfsio.createFrame( width, height );
  pfs::Channel *Xj, *Yj, *Zj;
  frameout->createXYZChannels( Xj, Yj, Zj );
  if( Xj==NULL || Yj==NULL || Zj==NULL )
    throw pfs::Exception( "could not allocate memory for hdr output" );
  else
    VERBOSE_STR << "memory for hdr output has been allocated" << endl;

  long sp = 0; // saturated pixels
  // camera response functions for each channel
  float* Ir = new float[M];
  float* Ig = new float[M];
  float* Ib = new float[M];
  if( Ir==NULL || Ig==NULL || Ib==NULL )
    throw pfs::Exception( "could not allocate memory for camera responses" );
  //2) response curves, either predefined (log,gamma,lin,from_file) or calibration from the set of images.
  switch( opt_response )
  {
  case FROM_FILE:
    {
    // read camera response from file
    bool loadR_ok = responseLoad(responseFile, Ir, M);
    bool loadG_ok = responseLoad(responseFile, Ig, M);
    bool loadB_ok = responseLoad(responseFile, Ib, M);
    fclose(responseFile);
    if( !loadR_ok || !loadG_ok || !loadB_ok )
      throw pfs::Exception( "could not load response curve from file" );
    }
    break;
  case LINEAR:
    responseLinear( Ir, M );
    responseLinear( Ig, M );
    responseLinear( Ib, M );
    //VERBOSE_STR << "created linear response curves" << endl;
    break;
  case GAMMA:
    responseGamma( Ir, M );
    responseGamma( Ig, M );
    responseGamma( Ib, M );
    break;
  case LOG10:
    responseLog10( Ir, M );
    responseLog10( Ig, M );
    responseLog10( Ib, M );
    break;
  case FROM_ROBERTSON:
    responseLinear( Ir, M );
    responseLinear( Ig, M );
    responseLinear( Ib, M );
    //call robertson02_getResponse method which computes both the Ir,Ig,Ib and the output HDR (i.e. its channels Xj,Yj,Zj).
    VERBOSE_STR << "recovering R channel..." << endl;
    sp = robertson02_getResponse( Xj, imgsR, Ir, w, M);
    VERBOSE_STR << "recovering G channel..." << endl;
    sp += robertson02_getResponse( Yj, imgsG, Ig, w, M);
    VERBOSE_STR << "recovering B channel..." << endl;
    sp += robertson02_getResponse( Zj, imgsB, Ib, w, M);
    sp /= 3;
    break;
  default:
    throw pfs::Exception( "camera response not computed/loaded/initialized" );
    break;
  }

  //3) apply model to generate hdr. If model is robertson and user preference was to compute response curve from image dataset using robertson algorithm, i.e. to calibrate, we are done, the robertson02_getResponse function has also computed the HDR (i.e. its channels).
  switch( opt_model )
  {
  case ROBERTSON:
    if (opt_response==FROM_ROBERTSON)
      break;
    else {
      //apply robertson model
      VERBOSE_STR << "applying response to R channel..." << endl;
      sp =  robertson02_applyResponse( Xj, imgsR, Ir, w, M);
      VERBOSE_STR << "applying reWonse to G channel..." << endl;
      sp += robertson02_applyResponse( Yj, imgsG, Ig, w, M);
      VERBOSE_STR << "applying response to B channel..." << endl;
      sp += robertson02_applyResponse( Zj, imgsB, Ib, w, M);
    }
    break;
  case DEBEVEC:
    //apply devebec model
    devebec_applyResponse(  Xj, imgsR, Ir,
			    Yj, imgsG, Ig,
			    Zj, imgsB, Ib, w, M);
    break;
  default:
    throw pfs::Exception("calibration method not set or not supported");
  }

  if( sp>0 )
  {
    float perc = ceilf(100.0f*sp/size);
    VERBOSE_STR << "saturated pixels found in " << perc << "% of the image!" << endl;
  }

  // save response curve to a given file
  if( responseSaveFile!=NULL )
  {
      responseSave(responseSaveFile, Ir, M, "IR");
      responseSave(responseSaveFile, Ig, M, "IG");
      responseSave(responseSaveFile, Ib, M, "IB");
      weightsSave(responseSaveFile, w, M, "W"); //w is saved as well for compatibility, even if we don't load it at startup
      fclose(responseSaveFile);
  }
  // output PFS stream with calibrated response
  pfs::transformColorSpace( pfs::CS_RGB, Xj, Yj, Zj, pfs::CS_XYZ, Xj, Yj, Zj );
  pfsio.writeFrame( frameout, stdout );

  VERBOSE_STR << "HDR image written to stdout" << endl;

	// clean up memory
	pfsio.freeFrame( frameout );
	delete[] w;
	delete[] Ir;
	delete[] Ig;
	delete[] Ib;

	for( unsigned int i=0 ; i<imgsR.size() ; i++ ) {
		delete imgsR[i].yi;
		delete imgsG[i].yi;
		delete imgsB[i].yi;
	}
}


int main( int argc, char* argv[] )
{
  try {
    pfshdrcalibrate( argc, argv );
  }
  catch( pfs::Exception ex ) {
    fprintf( stderr, PROG_NAME " error: %s\n", ex.getMessage() );
    return EXIT_FAILURE;
  }        
  catch( QuietException  ex ) {
    return EXIT_FAILURE;
  }        
  return EXIT_SUCCESS;
}
