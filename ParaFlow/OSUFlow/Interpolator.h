/*************************************************************************
*						OSU Flow Vector Field							 *
*																		 *
*																		 *
*	Created:	Han-Wei Shen, Liya Li									 *
*				The Ohio State University								 *
*	Date:		06/2005													 *
*																		 *
*	Interpolator														 *
*************************************************************************/

#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "header.h"

enum LerpType {LINEAR_LERP, NEAREST_LERP};

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

// barycentric interpolation
double BaryInterp(double dataValue[4], double coeff[3]);

// trilinear interpolation
double TriLerp(double lll, double hll, double lhl, double hhl, double llh, double hlh, double lhh, double hhh, double coeff[3]);

// bilinear interpolation
double BiLerp(double ll, double hl, double lh, double hh, double coeff[2]);

// linear interpolation
double Lerp(double x, double y, double ratio);

// Gaussian smoothing filter
void operateGaussianLPF(int width, int height, int element, double *pData);

#endif
