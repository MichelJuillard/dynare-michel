
#include "dw_matrix_array.h"
#include "dw_error.h"
#include "bmatrix.h"

#include <stdlib.h>

#include "modify_for_mex.h"

TElementSpecification dw_VectorSpecs =
  {
    dw_ARRAY_POINTER,
    sizeof(TVector),
    sizeof(int)+sizeof(TElementSpecification*),
    (void (*)(void*))FreeVector,
    DefaultPointerConstructor,
    (void* (*)(void*,void*))EquateVector,
    NULL,
    (int (*)(FILE*,void*,char*))dw_PrintVector,
    (int (*)(FILE*,void*))dw_ReadVector
  };

TElementSpecification dw_MatrixSpecs =
  {
    dw_ARRAY_POINTER,
    sizeof(TMatrix),
    sizeof(int)+sizeof(TElementSpecification*),
    (void (*)(void*))FreeMatrix,
    DefaultPointerConstructor,
    (void* (*)(void*,void*))EquateMatrix,
    NULL,
    (int (*)(FILE*,void*,char*))dw_PrintMatrix,
    (int (*)(FILE*,void*))dw_ReadMatrix
  };

/******************************************************************************/
/******************************* Initializaton ********************************/
/******************************************************************************/
void* dw_InitializeArray_vector(void* x, PRECISION y)
{
  int i;
  if (!x)
    dw_Error(NULL_ERR);
  else
    if (dw_IsArrayA(x))
      for (i=dw_DimA(x)-1; i >= 0; i--) dw_InitializeArray_vector(((void**)x)[i],y);
    else
      for (i=dw_DimA(x)-1; i >= 0; i--) InitializeVector(((void**)x)[i],y);
  return x;
}

void* dw_InitializeArray_matrix(void* X, PRECISION y)
{
  int i;
  if (!X)
    dw_Error(NULL_ERR);
  else
    if (dw_IsArrayA(X))
      for (i=dw_DimA(X)-1; i >= 0; i--) dw_InitializeArray_matrix(((void**)X)[i],y);
    else
      for (i=dw_DimA(X)-1; i >= 0; i--) InitializeMatrix(((void**)X)[i],y);
  return X;
}

/******************************************************************************/
/****************************** Tensor Calculus *******************************/
/******************************************************************************/
/*
   Assumes:
    X - r x s matrix or null pointer
    Y - k dimensional array of matrices

   Returns:
    The the tensor product

        Y[0] x Y[1] x ... x Y[k-1]

    If X is null, then space for the tensor product is allocated.  If X is not
    null, then the dimensions must match.

        r=RowM(Y[0]) x ... x RowM(Y[k-1])
        c=ColM(Y[0]) x ... x ColM(Y[k-1])

    Notes:
     Calls bMatrixTensor().
*/
TMatrix MatrixTensor(TMatrix X, TMatrix* Y)
{
  int i, r=1, c=1;
  PRECISION *Z, *U, *V, *W;
  TMatrix rtrn;
  if (!Y)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  for (i=dw_DimA(Y)-1; i >= 0; i--)
    if (!Y[i])
      {
    dw_Error(NULL_ERR);
    return (TMatrix)NULL;
      }
    else
      {
    r*=RowM(Y[i]);
    c*=ColM(Y[i]);
      }
  if (!X)
    {
      if (!(rtrn=CreateMatrix(r,c)))
    return (TMatrix)NULL;
    }
  else
    if ((r != RowM(X)) || (c != ColM(X)))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }
    else
      rtrn=X;
  if (dw_DimA(Y) > 2)
    {
      if (!(Z=(PRECISION*)swzMalloc(r*c*sizeof(PRECISION))))
    {
      if (!X) FreeMatrix(rtrn);
      return (TMatrix)NULL;
    }
      if (dw_DimA(Y) % 2)
    {
      U=Z;
      V=pElementM(rtrn);
    }
      else
    {
      U=pElementM(rtrn);
      V=Z;
    }
      i=dw_DimA(Y)-2;
      bMatrixTensor(U,pElementM(Y[i]),pElementM(Y[i+1]),RowM(Y[i]),ColM(Y[i]),RowM(Y[i+1]),
            ColM(Y[i+1]),MajorForm(rtrn),MajorForm(Y[i]),MajorForm(Y[i+1]));
      r=RowM(Y[i])*RowM(Y[i+1]);
      c=ColM(Y[i])*ColM(Y[i+1]);
      while (--i >= 0)
    {
      bMatrixTensor(V,pElementM(Y[i]),U,RowM(Y[i]),ColM(Y[i]),r,c,MajorForm(rtrn),MajorForm(Y[i]),MajorForm(rtrn));
      r*=RowM(Y[i]);
      c*=ColM(Y[i]);
      W=U;
      U=V;
      V=W;
    }
      swzFree(Z);
    }
  else
    if (dw_DimA(Y) > 1)
      bMatrixTensor(pElementM(rtrn),pElementM(Y[0]),pElementM(Y[1]),RowM(Y[0]),ColM(Y[0]),RowM(Y[1]),
            ColM(Y[1]),MajorForm(rtrn),MajorForm(Y[0]),MajorForm(Y[1]));
    else
      EquateMatrix(rtrn,Y[0]);
  return rtrn;
}

/*
   Assumes:
    X - d dimensional vector or null pointer
    Y - k dimensional array of vectors

   Returns:
    The the tensor product

        y[0] x y[1] x ... x y[k-1]

    If x is null, then space for the tensor product is allocated.  If x is not
    null, then the dimensions must match.

        d=DimV(Y[0]) x ... x DimV(Y[k-1])

    Notes:
     Calls bVectorTensor().
*/
TVector VectorTensor(TVector x, TVector* y)
{
  int i, d=1;
  PRECISION *z, *u, *v, *w;
  TVector rtrn;
  if (!y)
    {
      dw_Error(NULL_ERR);
      return (TVector)NULL;
    }
  for (i=dw_DimA(y)-1; i >= 0; i--)
    if (!y[i])
      {
    dw_Error(NULL_ERR);
    return (TVector)NULL;
      }
    else
      d*=DimV(y[i]);
  if (!x)
    {
      if (!(rtrn=CreateVector(d)))
    return (TVector)NULL;
    }
  else
    if (d != DimV(x))
      {
    dw_Error(SIZE_ERR);
    return (TVector)NULL;
      }
    else
      rtrn=x;
  if (dw_DimA(y) > 2)
    {
      if (!(z=(PRECISION*)swzMalloc(d*sizeof(PRECISION))))
    {
      if (!x) FreeVector(rtrn);
      return (TVector)NULL;
    }
      if (dw_DimA(y) % 2)
    {
      u=z;
      v=pElementV(rtrn);
    }
      else
    {
      u=pElementV(rtrn);
      v=z;
    }
      i=dw_DimA(y)-2;
      bVectorTensor(u,pElementV(y[i]),pElementV(y[i+1]),DimV(y[i]),DimV(y[i+1]));
      d=DimV(y[i])*DimV(y[i+1]);
      while (--i >= 0)
    {
      bVectorTensor(v,pElementV(y[i]),u,DimV(y[i]),d);
      d*=DimV(y[i]);
      w=u;
      u=v;
      v=w;
    }
      swzFree(z);
    }
  else
    if (dw_DimA(y) > 1)
      bVectorTensor(pElementV(rtrn),pElementV(y[0]),pElementV(y[1]),DimV(y[0]),DimV(y[1]));
    else
      EquateVector(rtrn,y[0]);
  return rtrn;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


