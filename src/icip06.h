/**
 * @brief Reinhard ghost removal and Devebel hdr generation algorithm
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
 * @author Giuseppe Rota, <grota@users.sourceforge.net>
 *
 * $Id: reinhard06.h,v 1.3 2006/09/13 11:52:56 gkrawczyk Exp $
 */

#ifndef _icip06_h_
#define _icip06_h_
#include <array2d.h>
#include <responses.h>
#include <math.h>

/**
 * @brief icip06 anti-ghosting algorithm
 *
 * @param imgsR reference to input vector containing source exposures, channel R
 * @param imgsG reference to input vector containing source exposures, channel G
 * @param imgsB reference to input vector containing source exposures, channel B
 * @param xj [out] HDR image channel 1
 * @param yj [out] HDR image channel 2
 * @param zj [out] HDR image channel 3
 * @param imgsl reference to input vector containing source exposures, channel l
 * @param imgsalpha reference to input vector containing source exposures, channel alpha
 * @param imgsbeta reference to input vector containing source exposures, channel beta
 * @param I1 response curve for channel 1, to be found with robertson02
 * @param I2 response curve for channel 2, to be found with robertson02
 * @param I3 response curve for channel 3, to be found with robertson02
 * @param w initial exposure weights, has to be computed with exposure_weights_icip06()
 * @param iterations number of iterations, from commandline
 * @param Ptemp   width*height*#exposures temporary reference array of weights
 * @param P [out] width*height*#exposures reference array of weights
 */
int icip06_applyResponse(   const ExposureList &imgsR,
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
                             Array2DList &P);

/**
 * @brief Reinhard06 anti-ghosting algorithm
 *
 * @param iteration number of iterations for wich you want the algorithm to run
 * @param imgsl reference to input vector containing source exposures, channel l
 * @param imgsalpha reference to input vector containing source exposures, channel alpha
 * @param imgsbeta reference to input vector containing source exposures, channel beta
 * @param imgsR reference to input vector containing source exposures, channel R
 * @param imgsG reference to input vector containing source exposures, channel G
 * @param imgsB reference to input vector containing source exposures, channel B
 * @param w initial exposure weights, has to be computed with exposure_weights_icip06()
 * @param iterations number of iterations, from commandline
 * @param Ptemp   width*height*#exposures temporary reference array of weights
 * @param P [out] width*height*#exposures reference array of weights
 */
void reinhard06_anti_ghosting(  const Array2DList &imgsl,
                                const Array2DList &imgsalpha,
                                const Array2DList &imgsbeta,
                                const ExposureList &imgsR,
                                const ExposureList &imgsG,
                                const ExposureList &imgsB,
                                const float* w, const int iterations,
                                Array2DList &Ptemp,
                                Array2DList &P);

/**
 * @brief Create HDR image by applying response curve to given images using Devebec model
 *
 * @param xj [out] HDR image channel 1
 * @param yj [out] HDR image channel 2
 * @param zj [out] HDR image channel 3
 * @param imgsR reference to input vector containing source exposures, channel R
 * @param imgsG reference to input vector containing source exposures, channel G
 * @param imgsB reference to input vector containing source exposures, channel B
 * @param I1 response curve for channel 1, to be found with robertson02
 * @param I2 response curve for channel 2, to be found with robertson02
 * @param I3 response curve for channel 3, to be found with robertson02
 * @param P  width*height*#exposures array of weights
 */
int devebec_applyResponse( pfs::Array2D* xj, const ExposureList &imgs1, const float* I1,
                           pfs::Array2D* yj, const ExposureList &imgs2, const float* I2,
                           pfs::Array2D* zj, const ExposureList &imgs3, const float* I3,
                           const Array2DList &P );

int devebec_applyResponse( pfs::Array2D* xj, const ExposureList &imgs1, const float* I1,
                           pfs::Array2D* yj, const ExposureList &imgs2, const float* I2,
                           pfs::Array2D* zj, const ExposureList &imgs3, const float* I3,
                           const float *w, int M );

#endif /* #ifndef _reinhard06_h_ */
