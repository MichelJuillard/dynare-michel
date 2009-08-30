
#ifndef __MATRIX_ARRAY__
#define __MATRIX_ARRAY__

#include "matrix.h"
#include "dw_array.h"

extern TElementSpecification dw_VectorSpecs;
extern TElementSpecification dw_MatrixSpecs;

#define dw_CreateArray_vector(dim)     (TVector*)dw_CreateArray(&dw_VectorSpecs,dim)
void* dw_InitializeArray_vector(void* x, PRECISION y);
#define dw_CreateArray_matrix(dim)     (TMatrix*)dw_CreateArray(&dw_MatrixSpecs,dim)
void* dw_InitializeArray_matrix(void* X, PRECISION y);

#if (PRECISION_SIZE == 8)
  #define dw_CreateArray_scalar(dim)     (double*)dw_CreateArray(&dw_DoubleSpecs,dim)
  #define dw_CreateMultidimensionalArray_scalar(depth,dim) dw_CreateMultidimensionalArray(&dw_DoubleSpecs,depth,dim) 
  #define dw_CreateMultidimensionalArrayList_scalar  dw_CreateMultidimensionalArrayList_double
  #define dw_CreateRectangularArray_scalar(row,col) (double**)dw_CreateMultidimensionalArrayList_double(2,row,col)
  #define dw_InitializeArray_scalar(a,x)  dw_InitializeArray_double(a,x)
#else
  #define dw_CreateArray_scalar(dim) (float*)dw_CreateArray(&dw_FloatSpecs,dim)
  #define dw_CreateMultidimensionalArray_scalar(depth,dim) dw_CreateMultidimensionalArray(&dw_FloatSpecs,depth,dim)
  #define dw_CreateMultidimensionalArrayList_scalar  dw_CreateMultidimensionalArrayList_float
  #define dw_CreateRectangularArray_scalar(row,col) (float**)dw_CreateMultidimensionalArrayList_float(2,row,col)
  #define dw_InitializeArray_scalar(a,x)  dw_InitializeArray_float(a,x)
#endif

/* Tensor Calculus */
TMatrix MatrixTensor(TMatrix X, TMatrix* Y);
TVector VectorTensor(TVector x, TVector* y);

#endif
