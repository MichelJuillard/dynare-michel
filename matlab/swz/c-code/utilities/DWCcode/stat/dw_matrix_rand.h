
#ifndef __RANDOM_MATRIX__
#define __RANDOM_MATRIX__

#include "matrix.h"

/* Random Matrices and Vectors */
TVector dw_UniformVector(TVector x);
TMatrix dw_UniformMatrix(TMatrix X);
TVector dw_NormalVector(TVector x);
TMatrix dw_NormalMatrix(TMatrix X);
TVector dw_LogNormalVector(TVector x, PRECISION mean, PRECISION standard_deviation);
TMatrix dw_GammaMatrix(TMatrix X, TMatrix A, TMatrix B);
TMatrix dw_Wishart(TMatrix X, TMatrix S, int nu);
TVector dw_StudentT(TVector x, TMatrix T, int nu);
TMatrix dw_UniformOrthogonal(TMatrix Q);
TVector dw_UniformUnitSphere(TVector x);
TVector dw_UniformUnitBall(TVector x);

#endif
