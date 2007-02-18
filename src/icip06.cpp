#include <icip06.h>
#include <iostream>

#define TEST_NEW_UCHAR
#ifdef  TEST_NEW_UCHAR
#define GET_AS_UCHAR(V) reinterpret_cast<unsigned char&>(V)
#define GET_AS_FLOAT(V) (float)reinterpret_cast<unsigned char&>(V)
#else
#define GET_AS_UCHAR(V) (V)
#define GET_AS_FLOAT(V) (V)
#endif

int icip06_applyResponse(    const ExposureList &imgsR,
                             const ExposureList &imgsG,
                             const ExposureList &imgsB,
                             pfs::Array2D* Xj,
                             pfs::Array2D* Yj,
                             pfs::Array2D* Zj,
                             const Array2DList &imgsl,
                             const Array2DList &imgsalpha,
                             const Array2DList &imgsbeta,
                             const float* Ir,
                             const float* Ig,
                             const float* Ib,
                             const float* w, const int iterations,
                             Array2DList &Ptemp,
                             Array2DList &P)
{

int width = Xj->getCols();
int height = Yj->getRows();

// Create P and Ptemp, stacks of Array2DImpl; P will contain weigths computed by reinhard06_anti_ghosting funtion
for (unsigned int i=0;i<imgsR.size();i++) {
	pfs::Array2D* toaddfinal = new pfs::Array2DImpl(width,height);
	pfs::Array2D* toaddtemp  = new pfs::Array2DImpl(width,height);
// 	if (toaddfinal==NULL || toaddtemp==NULL)
// 	throw pfs::Exception( "Cannot allocate memory for anti-ghosting weights" );
// 	else {
	P.push_back(toaddfinal);
	Ptemp.push_back(toaddtemp);
// 	}
}

//perform anti-ghosting
reinhard06_anti_ghosting(imgsl,imgsalpha,imgsbeta,imgsR,imgsG,imgsB,w,iterations,Ptemp,P);

//we don't need these anymore
for(unsigned int i=0 ; i<imgsR.size() ; i++ ) {
	delete Ptemp[i];
	delete imgsl[i];
	delete imgsalpha[i];
	delete imgsbeta[i];
}
//apply devebec model to compute actual hdr image
int sp=devebec_applyResponse(   Xj, imgsR, Ir,
				Yj, imgsG, Ig,
				Zj, imgsB, Ib, P);
return sp;
}




void reinhard06_anti_ghosting(  const Array2DList &imgsl,
                                const Array2DList &imgsalpha,
                                const Array2DList &imgsbeta,
                                const ExposureList &imgsR,
                                const ExposureList &imgsG,
                                const ExposureList &imgsB,
                                const float* w, const int iterations,
                                Array2DList &Ptemp,
                                Array2DList &P )
{

// number of exposures
int N = imgsl.size();

// frame size
int width  = (*imgsR[0].yi).getCols();
int height = (*imgsR[0].yi).getRows();

//First step, initialize WeightsList P and Ptemp to its initial values, obtained via w
//for all pixels
for ( int j=0 ; j<width*height ; j++ )
	// for all exposures for each pixel
	for( int i=0 ; i<N ; i++ ) 
		(*Ptemp[i])(j)=(*P[i])(j)=(float)(w[(int)GET_AS_UCHAR((*imgsR[i].yi)(j))]+w[(int)GET_AS_UCHAR((*imgsG[i].yi)(j))]+w[(int)GET_AS_UCHAR((*imgsB[i].yi)(j))])/3.0f;

// constant used in next step
float k=pow(2.0f*M_PI,-2.5f);

//COMPUTE weightslist P for ITER iterations.
for ( int ITER=0; ITER<iterations ; ITER++ ) { // for all the iterations ITER

	// update ALL the Ptemp stack using previous P*w[rgb]/3
	if (ITER!=0)
	for ( int exposure=0 ; exposure<N ; exposure++ ) { //for all the exposures
		for ( int row=0 ; row<height; row++ ) { //for all the rows row (from 0 to height-1)
		for ( int col=0 ; col<width;  col++ ) { //for all the columns col (from 0 to width-1)
		(*Ptemp[exposure])(col,row) = (*P[exposure])(col,row) * (float)(w[(int)GET_AS_UCHAR((*imgsR[exposure].yi)(col,row))]+w[(int)GET_AS_UCHAR((*imgsG[exposure].yi)(col,row))]+w[(int)GET_AS_UCHAR((*imgsB[exposure].yi)(col,row))])/3.0f;
		} //rows in pixel col,row
		} //cols in pixel col,row
	} //exposures i in iteration IT

	for ( int exposure=0 ; exposure<N ; exposure++ ) { //for all the exposures
		for ( int row=0 ; row<height; row++ ) { //for all the rows row (from 0 to height-1)
		for ( int col=0 ; col<width;  col++ ) { //for all the columns col (from 0 to width-1)

			// cumulative values for (*P[i])(col,row)
			float sum=0;
			float div=0;
			// get the 3 channel values for the current pixel, exposure, location col,row
			float Ch1x=(*imgsl[exposure])(col,row);     //should be a L value
			float Ch2x=(*imgsalpha[exposure])(col,row); //should be a alpha value
			float Ch3x=(*imgsbeta[exposure])(col,row);  //should be a beta value
			float Ch1y,Ch2y,Ch3y;

			// go through neighborhood, deep into the exposures, to update weight
			// use Ptemp to read current values, and write in P -at the end of the for-.
			for ( int I=0 ; I<N ; I++ ) {
				// if top-left is in boundaries
				if (col-1>=0 && row-1>=0) {
					Ch1y=(*imgsl[I])(col-1,row-1);
					Ch2y=(*imgsalpha[I])(col-1,row-1);
					Ch3y=(*imgsbeta[I])(col-1,row-1);
					sum+=(*Ptemp[I])(col-1,row-1)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 2.0f ));
					div+=(*Ptemp[I])(col-1,row-1);
				}
				//if top is in boundaries
				if (row-1>=0) {
					Ch1y=(*imgsl[I])(col,row-1);
					Ch2y=(*imgsalpha[I])(col,row-1);
					Ch3y=(*imgsbeta[I])(col,row-1);
					float a=(*Ptemp[I])(col,row-1)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 1.0f ));
					float b=(*Ptemp[I])(col,row-1);
					sum+=a;
					div+=b;
				}
				// if top-right is in boundaries
				if (col+1<width && row-1>=0) {
					Ch1y=(*imgsl[I])(col+1,row-1);
					Ch2y=(*imgsalpha[I])(col+1,row-1);
					Ch3y=(*imgsbeta[I])(col+1,row-1);
					sum+=(*Ptemp[I])(col+1,row-1)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 2.0f ));
					div+=(*Ptemp[I])(col+1,row-1);
				}
				// if left is in boundaries
				if (col-1>=0) {
					Ch1y=(*imgsl[I])(col-1,row);
					Ch2y=(*imgsalpha[I])(col-1,row);
					Ch3y=(*imgsbeta[I])(col-1,row);
					sum+=(*Ptemp[I])(col-1,row)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 1.0f ));
					div+=(*Ptemp[I])(col-1,row);
				}
				// if right is in boundaries
				if (col+1<width) {
					Ch1y=(*imgsl[I])(col+1,row);
					Ch2y=(*imgsalpha[I])(col+1,row);
					Ch3y=(*imgsbeta[I])(col+1,row);
					sum+=(*Ptemp[I])(col+1,row)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 1.0f ));
					div+=(*Ptemp[I])(col+1,row);
				}
				// if bottom-left is in boundaries
				if (col-1>=0 && row+1<height) {
					Ch1y=(*imgsl[I])(col-1,row+1);
					Ch2y=(*imgsalpha[I])(col-1,row+1);
					Ch3y=(*imgsbeta[I])(col-1,row+1);
					sum+=(*Ptemp[I])(col-1,row+1)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 2.0f ));
					div+=(*Ptemp[I])(col-1,row+1);
				}
				// if bottom is in boundaries
				if (row+1<height) {
					Ch1y=(*imgsl[I])(col,row+1);
					Ch2y=(*imgsalpha[I])(col,row+1);
					Ch3y=(*imgsbeta[I])(col,row+1);
					sum+=(*Ptemp[I])(col,row+1)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 1.0f ));
					div+=(*Ptemp[I])(col,row+1);
				}
				// if bottom-right is in boundaries
				if (col+1<width && row+1<height) {
					Ch1y=(*imgsl[I])(col+1,row+1);
					Ch2y=(*imgsalpha[I])(col+1,row+1);
					Ch3y=(*imgsbeta[I])(col+1,row+1);
					sum+=(*Ptemp[I])(col+1,row+1)*k*expf(-0.5*((Ch1x-Ch1y)*(Ch1x-Ch1y) + (Ch2x-Ch2y)*(Ch2x-Ch2y) + (Ch3x-Ch3y)*(Ch3x-Ch3y) + 2.0f ));
					div+=(*Ptemp[I])(col+1,row+1);
				}
			} // exposure-wise neighborood
			if (div!=0.0f)
				(*P[exposure])(col,row)=(sum/div);
			else
				(*P[exposure])(col,row)=0.0f;
		} //rows in pixel col,row
		} //cols in pixel col,row
	} //exposures in iteration ITER
std::cerr << "finished iteration " << ITER+1 << " of " << iterations << std::endl;
} //iterations

}


//TODO integrate anti-saturation
int devebec_applyResponse( pfs::Array2D* xj, const ExposureList &imgs1, const float* I1,
                           pfs::Array2D* yj, const ExposureList &imgs2, const float* I2,
                           pfs::Array2D* zj, const ExposureList &imgs3, const float* I3,
                           const Array2DList &P )
{
  // number of exposures
  int N = imgs1.size();

  // frame size
  int width = xj->getCols();
  int height = xj->getRows();

  // number of saturated pixels
  int saturated_pixels = 0;

  // for all pixels
  for( int j=0 ; j<width*height ; j++ )
  {
    float sum1 = 0.0f;
    float sum2 = 0.0f;
    float sum3 = 0.0f;

    float div1 = 0.0f;
    float div2 = 0.0f;
    float div3 = 0.0f;

    // for all exposures for each pixel
    for( int i=0 ; i<N ; i++ )
    {
      //pick the 3 channel values
      int m1 = (int)GET_AS_UCHAR((*imgs1[i].yi)(j));
      int m2 = (int)GET_AS_UCHAR((*imgs2[i].yi)(j));
      int m3 = (int)GET_AS_UCHAR((*imgs3[i].yi)(j));

      float ti = imgs1[i].ti;

      sum1 += (*P[i])(j) * I1[m1] / float(ti);
      div1 += (*P[i])(j);
      sum2 += (*P[i])(j) * I2[m2] / float(ti);
      div2 += (*P[i])(j);
      sum3 += (*P[i])(j) * I3[m3] / float(ti);
      div3 += (*P[i])(j);
    } //END for all the exposures

    if( div1==0.0f || div2==0.0f || div3==0.0f ) {
        saturated_pixels++;
    }

    if( div1!=0.0f && div2!=0.0f && div3!=0.0f ) {
      (*xj)(j) = sum1/div1;
      (*yj)(j) = sum2/div2;
      (*zj)(j) = sum3/div3;
    }
    else {
      (*xj)(j) = 0.0f;
      (*yj)(j) = 0.0f;
      (*zj)(j) = 0.0f;
    }
  }
  return saturated_pixels;
}



int devebec_applyResponse( pfs::Array2D* xj, const ExposureList &imgs1, const float* Ir,
                           pfs::Array2D* yj, const ExposureList &imgs2, const float* Ig,
                           pfs::Array2D* zj, const ExposureList &imgs3, const float* Ib,
                           const float* w, int M )
{
  // number of exposures
  int N = imgs1.size();

  // frame size
  int width = xj->getCols();
  int height = xj->getRows();

  // number of saturated pixels
  int saturated_pixels = 0;

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

  // for all pixels
  for( int j=0 ; j<width*height ; j++ )
  {
    float sumR = 0.0f;
    float sumG = 0.0f;
    float sumB = 0.0f;

    float divR = 0.0f;
    float divG = 0.0f;
    float divB = 0.0f;

    float maxti = -1e6f;
    float minti = +1e6f;

    int index_for_whiteR=-1;
    int index_for_whiteG=-1;
    int index_for_whiteB=-1;

    int index_for_blackR=-1;
    int index_for_blackG=-1;
    int index_for_blackB=-1;

    bool saturated_exposure=false;

    // for all exposures
    for( int i=0 ; i<N ; i++ )
    {
      //pick the 3 channel values
      int mR = (int)GET_AS_UCHAR((*imgs1[i].yi)(j));
      int mG = (int)GET_AS_UCHAR((*imgs2[i].yi)(j));
      int mB = (int)GET_AS_UCHAR((*imgs3[i].yi)(j));

      float ti = imgs1[i].ti;

      //if at least one of the color channel's values are in the bright "not-trusted zone" and we have min exposure time
      if ( (mR>maxM || mG>maxM || mB>maxM) && (ti<minti) ) {
        //update the indexes_for_whiteRGB, minti, and set this exposure as saturated (for this pixel)
        //so that we don't increment the sumRGB and divRGB variables (weighted average eq)
        index_for_whiteR=mR;
        index_for_whiteG=mG;
        index_for_whiteB=mB;
        minti=ti;
        saturated_exposure=true;
      }

      //if at least one of the color channel's values are in the dim "not-trusted zone" and we have max exposure time
      if ( (mR<minM || mG<minM || mB<minM) && (ti>maxti) ) {
        //update the indexes_for_blackRGB, maxti, and set this exposure as saturated (for this pixel)
        //so that we don't increment the sumRGB and divRGB variables (weighted average eq)
        index_for_blackR=mR;
        index_for_blackG=mG;
        index_for_blackB=mB;
        maxti=ti;
        saturated_exposure=true;
      }

      float w_average=(w[mR]+w[mG]+w[mB])/3.0f;
      sumR += w_average * Ir[mR] / float(ti);
      divR += w_average;
      sumG += w_average * Ig[mG] / float(ti);
      divG += w_average;
      sumB += w_average * Ib[mB] / float(ti);
      divB += w_average;

    } //END for all the exposures

    if( divR==0.0f || divG==0.0f || divB==0.0f ) {
        saturated_pixels++;
        if (maxti>-1e6f) {
          sumR = Ir[index_for_blackR] / float(maxti);
          sumG = Ig[index_for_blackG] / float(maxti);
          sumB = Ib[index_for_blackB] / float(maxti);
          divR = divG = divB = 1.0f;
        }
        if (minti<+1e6f) {
          sumR = Ir[index_for_whiteR] / float(minti);
          sumG = Ig[index_for_whiteG] / float(minti);
          sumB = Ib[index_for_whiteB] / float(minti);
          divR = divG = divB = 1.0f;
        }
    }

    if( divR!=0.0f && divG!=0.0f && divB!=0.0f ) {
      (*xj)(j) = sumR/divR;
      (*yj)(j) = sumG/divG;
      (*zj)(j) = sumB/divB;
    }
    else {
      //we shouldn't be here anyway...
      (*xj)(j) = 0.0f;
      (*yj)(j) = 0.0f;
      (*zj)(j) = 0.0f;
    }
  }
  return saturated_pixels;
}


