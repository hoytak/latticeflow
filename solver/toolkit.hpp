#ifndef _TOOLKIT_H_
#define _TOOLKIT_H_

#include <cmath>

void pinnedLevelSets(double *L, double *X, 
		     size_t nz, size_t ny, size_t nx, 
		     double wx, double wy, double wz);

static const int Kernel2d_4 = 0;
static const int Kernel2d_8 = 1;
static const int Kernel2d_16 = 2;
static const int Kernel2d_24 = 3;

void getRGBImageLevelSets(double *L, double *X, double *predictions, 
			  double *labels, 
			  size_t nx, size_t ny, double e2_balance, double label_value, 
			  int kernel);


#endif /* _TOOLKIT_H_ */
