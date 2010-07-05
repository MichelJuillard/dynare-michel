
#include "bmatrix.h"
#include "dw_error.h"

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

/********************/
#include "blas_lapack.h"
/********************/

/********************
#include "mkl.h"
/********************/

#include "modify_for_mex.h"

/******************************************************************************/
/***************************** Uniary Operations ******************************/
/******************************************************************************/
/*
  Assumes:
   x : n-vector
   y : n-vector
   n : positive

  Results:
   x[i] = -y[i] for 0 <= i < n

  Returns:
   0 upon success

  Notes:
   x and y do not have to be distinct
*/
int bNegative(PRECISION *x, PRECISION *y, int n)
{
 while (--n >= 0) x[n]=-y[n];
 return NO_ERR;
}

/*
  Assumes:
   x : n-vector
   y : n-vector
   n : positive

  Results:
   x[i] = fabs(y[i]) for 0 <= i < n

  Returns:
   0 upon success

  Notes:
   x and y do not have to be distinct
*/
int bAbs(PRECISION *x, PRECISION *y, int n)
{
 while (--n >= 0) x[n]=fabs(y[n]);
 return NO_ERR;
}

/*
  Assumes:
   x : array of length m*n
   y : array of length m*n
   m : positive
   n : positive
   t : 0 or 1


  Results:
               x          y
        t   (n x m)    (m x n)    results
       -----------------------------------------
        0  row major  row major    x = y'
        1  col major  col major    x = y'

   or
               x          y
        t   (m x n)    (m x n)    results
       -----------------------------------------
        0  col major  row major    x = y
        1  row major  col major    x = y

   or
               x          y
        t  row major  row major   results
       -----------------------------------------
        0    n x m      m x n      x = y'
        1    m x n      n x m      x = y'
   or
               x          y
        t  col major  col major   results
       -----------------------------------------
        0    m x n      n x m      x = y'
        1    n x m      m x n      x = y'

*/
int bTranspose(PRECISION *x, PRECISION *y, int m, int n, int t)
{
 int i, j, k;
 if (t)
   for (i=k=m*n-1; k >= 0; i--)
     for (j=i; j >= 0; j-=m)
       x[k--]=y[j];
 else
   for (i=k=m*n-1; k >= 0; i--)
     for (j=i; j >= 0; j-=n)
       x[k--]=y[j];
 return NO_ERR;
}

/*
  Assumes:
   x : array of length m*m
   m : positive

  Results:
   x = y'

  Notes:
   The major format (row or column) does not matter.
*/
int bTransposeInPlace(PRECISION *x, int m)
{
 PRECISION tmp;
 int i, j;
 for (j=m*m-2; j > 0; j+=i-1)
   for (i=j-m+1; i >= 0; j--, i-=m)
     {
       tmp=x[i];
       x[i]=x[j];
       x[j]=tmp;
     }
 return NO_ERR;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/***************************** Addition Routines ******************************/
/******************************************************************************/
/*
  Assumes:
   x : n-vector
   y : n-vector
   z : n-vector
   n : positive

  Results:
   x[i] = y[i] + z[i]  for 0 <= i < n

  Returns:
   0 upon success

  Notes:
   x, y and z do not have to be distinct
*/
int bAdd(PRECISION *x, PRECISION *y, PRECISION *z, int n)
{
 while (--n >= 0) x[n]=y[n]+z[n];
 return NO_ERR;
}

/*
  Assumes:
   x : n-vector
   y : n-vector
   z : n-vector
   n : positive

  Results
   x[i] = y[i] - z[i]  for 0 <= i < n

  Returns:
   0 upon success

  Notes:
   x, y and z do not have to be distinct
*/
int bSubtract(PRECISION *x, PRECISION *y, PRECISION *z, int n)
{
 while (--n >= 0) x[n]=y[n]-z[n];
 return NO_ERR;
}

/*
  Assumes:
   x : scalar array of dimension at least m
   y : scalar array of dimension at least m
   a : scalar
   m : positive

  Results:
   x = x + a*y

  Returns:
   NO_ERR

  Notes:
   x and y should be distinct.
*/
int bLinearUpdateScalar(PRECISION *x, PRECISION *y, PRECISION a, int m)
{
  blas_int inc = 1;
  blas_int m2 = m;
#if (PRECISION_SIZE == 4)
  saxpy(&m2,&a,y,&inc,x,&inc);
#else
  daxpy(&m2,&a,y,&inc,x,&inc);
#endif
  return NO_ERR;
}

/*
  Assumes:
   x : m*n
   y : m*n
   z : m*n
   m : positive
   n : positive
   xt: 0 or 1
   yt: 0 or 1
   zt: 0 or 1

  Results:
   x = y + z

  Returns:
   0 upon success

  Notes:
   x, y and z do not have to be distinct
*/
int bMatrixAdd(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int xt, int yt, int zt)
{
  int i, j, k, s;
  if (xt == yt)
    if (yt == zt)
      for (k=m*n-1; k >= 0; k--) x[k]=y[k]+z[k];
    else
      for (s=zt ? m : n, k=i=m*n-1; k >= 0; i--)
    for (j=i; j >= 0; k--, j-=s)
      x[k]=y[k]+z[j];
  else
    if (yt == zt)
      for (s=yt ? m : n, k=i=m*n-1; k >= 0; i--)
    for (j=i; j >= 0; k--, j-=s)
      x[k]=y[j]+z[j];
    else
      for (s=yt ? m : n, k=i=m*n-1; k >= 0; i--)
    for (j=i; j >= 0; k--, j-=s)
      x[k]=y[j]+z[k];
 return NO_ERR;
}

int bMatrixSubtract(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int xt, int yt, int zt)
{
  int i, j, k, s;
  if (xt == yt)
    if (yt == zt)
      for (k=m*n-1; k >= 0; k--) x[k]=y[k]-z[k];
    else
      for (s=zt ? m : n, k=i=m*n-1; k >= 0; i--)
    for (j=i; j >= 0; k--, j-=s)
      x[k]=y[k]-z[j];
  else
    if (yt == zt)
      for (s=yt ? m : n, k=i=m*n-1; k >= 0; i--)
    for (j=i; j >= 0; k--, j-=s)
      x[k]=y[j]-z[j];
    else
      for (s=yt ? m : n, k=i=m*n-1; k >= 0; i--)
    for (j=i; j >= 0; k--, j-=s)
      x[k]=y[j]-z[k];
 return NO_ERR;
}

 /*
  Assumes:
   x : m vector
   y : m vector
   z : m vector
   m : positive
   n : positive

  Results:
   x = a*y + b*z

  Returns:
   0 upon success

  Notes:
   x, y and z do not have to be distinct
*/
int bLinearCombination(PRECISION *x, PRECISION a, PRECISION *y, PRECISION b, PRECISION *z, int m)
{
  while (--m >= 0) x[m]=a*y[m]+b*z[m];
  return NO_ERR;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
  Assumes:
   x : m x n matrix
   y : n-vector
   n : positive

  Results:
   x[i] = s * y[i] for 0 <= i < n

  Returns:
   0 upon success

  Notes:
  x and y do not have to be distinct
*/
int bMultiply(PRECISION *x, PRECISION *y, PRECISION s, int n)
{
  while (--n >= 0) x[n]=s*y[n];
  return NO_ERR;
}


/*
  Assumes:
   x : array of length mn
   y : array of length mp
   z : array of length pn
   m, n and p are positive
   xt, yt, and zt are 0 or 1

  Results:
                          x          y          z
       xt   yt   zt    (m x n)    (m x p)    (p x n)     results
       -----------------------------------------------------------
        0    0    0   row major  row major  row major   x = y * z
        1    0    0   col major  row major  row major   x = y * z
        0    1    0   row major  col major  row major   x = y * z
        1    1    0   col major  col major  row major   x = y * z
        0    0    1   row major  row major  col major   x = y * z
        1    0    1   col major  row major  col major   x = y * z
        0    1    1   row major  col major  col major   x = y * z
        1    1    1   col major  col major  col major   x = y * z

   or
                          x          y          z
       xt   yt   zt   row major  row major  row major    results
       ----------------------------------------------------------
        0    0    0     m x n      m x p      p x n     x = y * z
        1    0    0     n x m      m x p      p x n     x'= y * z
        0    1    0     m x n      p x m      p x n     x = y'* z
        1    1    0     n x m      p x m      p x n     x'= y'* z
        0    0    1     m x n      m x p      n x p     x = y * z'
        1    0    1     n x m      m x p      n x p     x'= y * z'
        0    1    1     m x n      p x m      n x p     x = y'* z'
        1    1    1     n x m      p x m      n x p     x'= y'* z'

   or

                          x          y          z
       xt   yt   zt   col major  col major  col major    results
       -----------------------------------------------------------
        0    0    0     n x m      p x m      n x p     x'= y'* z'
        1    0    0     m x n      p x m      n x p     x = y'* z'
        0    1    0     n x m      m x p      n x p     x'= y * z'
        1    1    0     m x n      m x p      n x p     x = y * z'
        0    0    1     n x m      p x m      p x n     x'= y'* z
        1    0    1     m x n      p x m      p x n     x = y'* z
        0    1    1     n x m      m x p      p x n     x'= y * z
        1    1    1     m x n      m x p      p x n     x = y * z

  Returns:
   0 upon success

  Notes:
   An (n x m) matrix x is in row major format if x[i][j]=x[i*n+j] and is in
   column major format if x[i][j]=x[i+j*m].

*/

int bMatrixMultiply(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int p, int xt, int yt, int zt)
{
  char transy, transz;
  blas_int dy, dz;
  PRECISION beta=0.0, alpha=1.0;

  blas_int m2 = m;
  blas_int n2 = n;
  blas_int p2 = p;

#if PRECISION_SIZE == 4
  if (xt)
    {
      if (yt) {transy='N'; dy=m;} else {transy='T'; dy=p;}
      if (zt) {transz='N'; dz=p;} else {transz='T'; dz=n;}
      sgemm(&transy,&transz,&m2,&n2,&p2,&alpha,y,&dy,z,&dz,&beta,x,&m2);
    }
  else
    {
      if (yt) {transy='T'; dy=m;} else {transy='N'; dy=p;}
      if (zt) {transz='T'; dz=p;} else {transz='N'; dz=n;}
      sgemm(&transz,&transy,&n2,&m2,&p2,&alpha,z,&dz,y,&dy,&beta,x,&n2);
    }
#else
  if (xt)
    {
      if (yt) {transy='N'; dy=m;} else {transy='T'; dy=p;}
      if (zt) {transz='N'; dz=p;} else {transz='T'; dz=n;}
      dgemm(&transy,&transz,&m2,&n2,&p2,&alpha,y,&dy,z,&dz,&beta,x,&m2);
    }
  else
    {
      if (yt) {transy='T'; dy=m;} else {transy='N'; dy=p;}
      if (zt) {transz='T'; dz=p;} else {transz='N'; dz=n;}
      dgemm(&transz,&transy,&n2,&m2,&p2,&alpha,z,&dz,y,&dy,&beta,x,&n2);
    }
#endif
  return NO_ERR;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/***************************** LU Decompositions ******************************/
/******************************************************************************/
/*
   Assumes
     p   : integer array of length at least q=min(m,n)
     x   : array of lenth mn
    m,n  : positive


  Results:

               x          L          U          P
       xt  row major  row major  row major   uses q       results
       --------------------------------------------------------------
        0    m x n      m x q      q x n      m x m    x = P * L * U
        1    n x m      n x q      q x m      m x m    x = L * U * P'

   or
               x          L          U          P
       xt  col major  col major  col major   uses q       results
       --------------------------------------------------------------
        0    n x m      n x q      q x m      m x m    x = L * U * P
        1    m x n      m x q      q x n      m x m    x = P * L * U

   or
               x           L         U          P
       xt   (m x n)     (m x q)   (q x n)    (m x m)     results
       --------------------------------------------------------------
        0  row major  row major  row major    uses q   x = P * L * U
        1  col major  col major  col major    uses q   x = P * L * U

   Results
     Computes the LU decomposition of A with partial pivoting. The LU
     decomposition of a matrix A is

                                 A = P * L * U

     where P is a (m x m) permutation matrix, L is a (m x q) lower triangular
     matrix with ones on the diagonal, U is a (q x n) upper triangular matrix,
     and q=min(m,n).  These matrices are stored as follows.

      U is stored in the upper part of x, including the diagonal.

      L is stored in the lower part of x.  The diagonal of L is not stored.

      The matrix P is defined by

                   P = P(0,p[0])*P(1,p[1])*...*P(q-1,p[q-1])

      where P(r,s) is the (m x m) matrix obtained by permuting the rth and sth
      rows of the (m x m) identity matrix.  It is assumed that i <= p[i] < m.

   Returns
     NO_ERR   - success
     SING_ERR - x was singular to machine precision.  LU decomposition is
                still computed and returned.

   Notes
     Uses partial pivoting.  Does not scale.  An (n x m) matrix x is in row major
     format if x[i][j]=x[i*n+j] and is in column major format if x[i][j]=x[i+j*m].
     Only q elements of p are set.
*/
int bLU(int *p, PRECISION *x, int m, int n, int xt)
{
#if PRECISION_SIZE == 4
  #define getrf sgetrf
#else
  #define getrf dgetrf
#endif

 PRECISION *y;
 int i;

 lapack_int m2 = m;
 lapack_int n2 = n;
 lapack_int info;
 lapack_int *p2;

 int minmn = (m < n) ? m : n;

 if(!(p2 = (lapack_int *)swzCalloc(minmn, sizeof(lapack_int))))
   return MEM_ERR;

 for(i=0; i<minmn; i++)
   p2[i] = p[i];

 if (xt)
   {
     getrf(&m2,&n2,x,&m2,p2,&info);
   }
 else
   {
     if (!( y=(PRECISION*)swzMalloc(m*n*sizeof(PRECISION)))) return MEM_ERR;
     bTranspose(y,x,m,n,0);

     getrf(&m2,&n2,y,&m2,p2,&info);

     bTranspose(x,y,m,n,1);
     free(y);
   }

 for(i=0; i<minmn; i++)
   p[i] = p2[i];
 free(p2);

 for (i=(m < n) ? m-1 : n-1; i >= 0; i--) p[i]--;
 return (info < 0) ? SING_ERR : NO_ERR;

#undef getrf
}

/*
   Assumes
     x : array of length m*m representing a triangular matrix.
     b : array of length m*n

   Results
                           x           b/y
        u   xt   bt     (m x m)      (m x n)      solve
        --------------------------------------------------
        0    0    0   L:row major   row major   L * y = b
        1    0    0   U:row major   row major   U * y = b
        1    1    0   U:col major   row major   U * y = b
        0    1    0   L:col major   row major   L * y = b
        0    0    1   L:row major   col major   L * y = b
        1    0    1   U:row major   col major   U * y = b
        1    1    1   U:col major   col major   U * y = b
        0    1    1   L:col major   col major   L * y = b

      or
                         x
                      (m x m)      b/y                                  u^xt^0
        u   xt   bt  row major  row major    solve                 (u^xt^major_form)
        ----------------------------------------------------------------------------
        0    0    0      L        m x n    L * y = b                       0
        1    0    0      U        m x n    U * y = b                       1
        1    1    0      L        m x n    L'* y = b                       0
        0    1    0      U        m x n    U'* y = b                       1
        0    0    1      L        n x m    L * y'= b'  (y * L'= b)         0
        1    0    1      U        n x m    U * y'= b'  (y * U'= b)         1
        1    1    1      L        n x m    L'* y'= b'  (y * L = b)         0
        0    1    1      U        n x m    U'* y'= b'  (y * U = b)         1

      or
                         x
                      (m x m)      b/y                                  u^xt^1
        u   xt   bt  col major  col major    solve                 (u^xt^major_form)
        ----------------------------------------------------------------------------
        0    0    0      U        m x n    U'* y'= b'  (y * U = b)        1
        1    0    0      L        m x n    L'* y'= b'  (y * L = b)        0
        1    1    0      U        m x n    U * y'= b'  (y * U'= b)        1
        0    1    0      L        m x n    L * y'= b'  (y * L'= b)        0
        0    0    1      U        n x m    U'* y = b                      1
        1    0    1      L        n x m    L'* y = b                      0
        1    1    1      U        n x m    U * y = b                      1
        0    1    1      L        n x m    L * y = b                      0

     The solution y is stored in b.

   Returns
     0 (NO_ERR) - success
     SING_ERR   - x is singular

   Notes
     Because this routines tests using the xor operator, the values of u and xt
     must be either 0 or 1.  The matrix x is assumed to be triangular.  Care
     must be taken that the matrix is either upper or lower triangular in the
     correct format.

     An (n x m) matrix x is in row major format if x[i][j]=x[i*n+j] and is in
     column major format if x[i][j]=x[i+j*m].
*/
int bSolveTriangular(PRECISION *x, PRECISION *b, int m, int n, int u, int xt, int bt)
{
 int i, j, k, bi, bj, xi, xj, mbi;
 PRECISION *pxx, *px, *pb;
 for (j=m+1, i=j*(m-1); i >= 0; i-=j) if (x[i] == 0.0) return SING_ERR;
 if (xt) { xi=1; xj=m; } else { xi=m; xj=1; }
 if (bt) { bi=1; bj=m; } else { bi=n; bj=1; }
 mbi=(m-1)*bi;
 if (u)
   for (x+=(m-1)*xj, j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
    for (pxx=x+(m-1)*xi, i=(m-1)*bi; i >= 0; pxx-=xi, i-=bi)
     {
      px=pxx;
      pb=b+i;
      for (k=mbi-i; k > 0; px-=xj, k-=bi) (*pb)-=(*px)*pb[k];
      *pb/=(*px);
     }
  else
   {
    for (j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
     for (pxx=x, i=0; i <= mbi; pxx+=xi, i+=bi)
      {
       px=pxx;
       pb=b+i;
       for (k=-i; k < 0; px+=xj, k+=bi) (*pb)-=(*px)*pb[k];
       *pb/=(*px);
      }
   }
 return NO_ERR;
}

/*
   Assumes
     x : array of length m*m representing a triangular matrix with unit diagonal
     b : array of length m*n

   Results
                           x           b/y
        u   xt   bt     (m x m)      (m x n)      solve
        --------------------------------------------------
        0    0    0   L:row major   row major   L * y = b
        1    0    0   U:row major   row major   U * y = b
        1    1    0   U:col major   row major   U * y = b
        0    1    0   L:col major   row major   L * y = b
        0    0    1   L:row major   col major   L * y = b
        1    0    1   U:row major   col major   U * y = b
        1    1    1   U:col major   col major   U * y = b
        0    1    1   L:col major   col major   L * y = b

      or
                         x
                      (m x m)      b/y                                  u^xt^0
        u   xt   bt  row major  row major    solve                 (u^xt^major_form)
        ----------------------------------------------------------------------------
        0    0    0      L        m x n    L * y = b                       0
        1    0    0      U        m x n    U * y = b                       1
        1    1    0      L        m x n    L'* y = b                       0
        0    1    0      U        m x n    U'* y = b                       1
        0    0    1      L        n x m    L * y'= b'  (y * L'= b)         0
        1    0    1      U        n x m    U * y'= b'  (y * U'= b)         1
        1    1    1      L        n x m    L'* y'= b'  (y * L = b)         0
        0    1    1      U        n x m    U'* y'= b'  (y * U = b)         1

      or
                         x
                      (m x m)      b/y                                  u^xt^1
        u   xt   bt  col major  col major    solve                 (u^xt^major_form)
        ----------------------------------------------------------------------------
        0    0    0      U        m x n    U'* y'= b'  (y * U = b)        1
        1    0    0      L        m x n    L'* y'= b'  (y * L = b)        0
        1    1    0      U        m x n    U * y'= b'  (y * U'= b)        1
        0    1    0      L        m x n    L * y'= b'  (y * L'= b)        0
        0    0    1      U        n x m    U'* y = b                      1
        1    0    1      L        n x m    L'* y = b                      0
        1    1    1      U        n x m    U * y = b                      1
        0    1    1      L        n x m    L * y = b                      0


     The solution y is stored in b.

   Returns
     0 (NO_ERR) - success

   Notes
     If f is zero for row major format and one for a column major format, then
     passing xt = 1^f (1 xor f) implies an upper triangular matrix is passed,
     and passing xt = 0^f  implies a lower triangular matrix is passed.

     The matrix x is assumed to be triangular with unit diagonal.  Care must be
     taken that the matrix is either upper or lower triangular in the correct
     format.

     An (n x m) matrix x is in row major format if x[i][j]=x[i*n+j] and is in
     column major format if x[i][j]=x[i+j*m].
*/
int bSolveUnitTriangular(PRECISION *x, PRECISION *b, int m, int n, int u, int xt, int bt)
{
 int i, j, k, bi, bj, xi, xj, mbi;
 PRECISION *pxx, *px, *pb;
 if (xt) { xi=1; xj=m; } else { xi=m; xj=1; }
 if (bt) { bi=1; bj=m; } else { bi=n; bj=1; }
 mbi=(m-1)*bi;
 if (u)
   for (x+=(m-1)*xj, j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
    for (pxx=x+(m-1)*xi, i=(m-1)*bi; i >= 0; pxx-=xi, i-=bi)
     {
      px=pxx;
      pb=b+i;
      for (k=mbi-i; k > 0; px-=xj, k-=bi) (*pb)-=(*px)*pb[k];
     }
  else
   {
    for (j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
     for (pxx=x, i=0; i <= mbi; pxx+=xi, i+=bi)
      {
       px=pxx;
       pb=b+i;
       for (k=-i; k < 0; px+=xj, k+=bi) (*pb)-=(*px)*pb[k];
      }
   }
 return NO_ERR;
}

/*
   Assumes
     p : integer array of length q  with 0 <= p[i] < m for all 0 <= i < q.
     y : array of length mn

   Results
                    x/y
        pt   yt  row major    product
        ---------------------------------------------
         0    0    m x n     x = P * y
         1    0    m x n     x = P'* y
         0    1    n x m     x'= P * y'  (x = y * P')
         1    1    n x m     x'= P'* y'  (x = y * P )

     or

                    x/y
        pt   yt  col major    product
        ---------------------------------------------
         0    0    n x m     x'= P * y'  (x = y * P')
         1    0    n x m     x'= P'* y'  (x = y * P )
         0    1    m x n     x = P * y
         1    1    m x n     x = P'* y

     or
                    x/y
        pt   yt   (m x n)     product
        -------------------------------
         0    0  row major   x = P * y
         1    0  row major   x = P'* y
         0    1  col major   x = P * y
         1    1  col major   x = P'* y

     The matrix P is defined by

                    P = P(0,p[0])*P(1,p[1])*...*P(q-1,p[q])

     where P(r,s) is the (m x m) matrix obtained by permuting the rth and sth
     rows of the (m x m) identity matrix.

   Notes:
     An (n x m) matrix x is in row major format if x[i][j]=x[i*n+j] and is in
     column major format if x[i][j]=x[i+j*m].
*/
int bPermutationMultiply(int *p, PRECISION *y, int m, int n, int q, int pt, int yt)
{
 int i, j, k, pk;
 PRECISION tmp;
 if (yt)
   if (pt)
     for (j=0; j < q; j++)
      {
       if (j != p[j])
        for (i=(n-1)*m; i >= 0; i-=m)
         {
          tmp=y[i+j];
          y[i+j]=y[i+p[j]];
          y[i+p[j]]=tmp;
         }
      }
    else
     for (j=q-1; j >= 0; j--)
      {
       if (j != p[j])
        for (i=(n-1)*m; i >= 0; i-=m)
         {
          tmp=y[i+j];
          y[i+j]=y[i+p[j]];
          y[i+p[j]]=tmp;
         }
      }
  else
   if (pt)
     for (i=0; i < q; i++)
      {
       if (i != p[i])
        {
         k=i*n;
         pk=p[i]*n;
         for (j=n-1; j >= 0; j--)
          {
           tmp=y[k+j];
           y[k+j]=y[pk+j];
           y[pk+j]=tmp;
          }
        }
      }
    else
     for (i=q-1; i >= 0; i--)
      {
       if (i != p[i])
        {
         k=i*n;
         pk=p[i]*n;
         for (j=n-1; j >= 0; j--)
          {
           tmp=y[k+j];
           y[k+j]=y[pk+j];
           y[pk+j]=tmp;
          }
        }
      }
 return NO_ERR;
}

/*
   Assumes
     p : integer array of length q with i <= p[i] < m for all 0 <= i < m.
     x : PRECISION array of length m*m

   Results
     If P(i,j) is the identity matrix with the ith and jth columns permuted,
     the the permutation matrix is defined from p is

                    P = P(0,p[0])*P(1,p[1])*...*P(q-1,p[q-1])

     Then
                      x
             xt   row major
            ----------------
             0      x = P
             1      x = P'

     or
                      x
             xt   col major
            ----------------
             0      x = P'
             1      x = P


     if where P(r,s) is the (m x m) matrix obtained by permuting the rth and sth
     rows of the (m x m) identity matrix.

   Notes:
     Permuting the ith and jth columns of an identity matrix is equivalent to
     permuting the ith and jth rows.  An (n x m) matrix x is in row major
     format if x[i][j]=x[i*n+j] and is in column major format if
     x[i][j]=x[i+j*m].
*/
int bPermutation(PRECISION *x, int *p, int m, int q, int xt)
{
  int i, j, k;
  for (k=m*m-1; k >= 0; k--) x[k]=0.0;
  if (xt)
    for (j=m-1; j >= 0; j--)
      {
    if (j < q)
      {
        k=j-1;
        i=p[j];
      }
    else
      {
        k=q-1;
        i=j;
      }
    for ( ; k >= 0; k--) if (i == p[k]) i=k;
    x[i+j*m]=1.0;
      }
  else
    for (j=m-1; j >= 0; j--)
      {
    if (j < q)
      {
        k=j-1;
        i=p[j];
      }
    else
      {
        k=q-1;
        i=j;
      }
    for ( ; k >= 0; k--) if (i == p[k]) i=k;
    x[i*m+j]=1.0;
      }
  return NO_ERR;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/************************ Singular Value Decomposition ************************/
/******************************************************************************/
/*
   Assumes
    U       : array of length m*m (compact=0) or m*q (compact=1) or null
    d       : array of length q=min(m,n)
    V       : array of length n*n (compact=0) or n*q (compact=1) or null
    A       : array of length m*n
    m       : positive
    n       : positive
    ut      : 0 or 1
    vt      : 0 or 1
    at      : 0 or 1
    compact : 0 or 1

   Returns
     NO_ERR     : success
     MEM_ERR    : out of memory

   Results
     Finds matrices U and V with orthonormal columns and a diagonal matrix
     D=diag(d) with non-negative diagonal such that A = U*D*V'.  The matrix D is
     m x n if compact = 0 and is q x q if compact = 1.  The elements of d are in
     descending order.  The flags ut, vt, and at determine the format of U, V,
     and A.  A value of 1 indicates column major format and a value of 0
     indicates row major format.  If either U or V is null, then it is not
     computed.

   Notes
     If A=U, U and A must be of the same size and ut=at.  If A=V, then V and A
     must be of the same size and vt=at.  It cannot be the case that U=V unless
     both are null.
*/
int bSVD_new(PRECISION *U, PRECISION *d, PRECISION *V, PRECISION *A, int m, int n, int ut, int vt, int at, int compact)
{
#if (PRECISION_SIZE == 4)
  #define gesvd sgesvd
#else
  #define gesvd dgesvd
#endif

  char jobu, jobv;
  int qu, qv, err=NO_ERR;
  PRECISION  *A_, *U_, *V_, *work, opt_size;

  lapack_int k, m2, n2, qv2, info;

  if (!(A_=(PRECISION*)swzMalloc(m*n*sizeof(PRECISION)))) return MEM_ERR;
  if (at)
    memcpy(A_,A,m*n*sizeof(PRECISION));
  else
    bTranspose(A_,A,m,n,at);

  if (!U)
    {
      jobu='N';
      qu=0;
      U_=(PRECISION*)NULL;
    }
  else
    {
      if (compact && (m > n))
    {
      jobu='S';
      qu=n;
    }
      else
    {
      jobu='A';
      qu=m;
    }
      if (ut)
    U_=U;
      else
    if (!(U_=(PRECISION*)swzMalloc(m*qu*sizeof(PRECISION))))
      {
        free(A_);
        return MEM_ERR;
      }
    }

  if (!V)
    {
      qv=1;
      jobv='N';
      V_=(PRECISION*)NULL;
    }
  else
    {
      if (compact && (m < n))
    {
      jobv='S';
      qv=m;
    }
      else
    {
      jobv='A';
      qv=n;
    }
      if (!vt)
    V_=V;
      else
    if (!(V_=(PRECISION*)swzMalloc(n*qv*sizeof(PRECISION))))
      {
        free(A_);
        if (U_ && (U_ != U)) free(U_);
        return MEM_ERR;
      }
    }

/*    // compute singular value decomposition   ansi-c*/
  k=-1;
  m2 = m;
  n2 = n;
  qv2 = qv;
  gesvd(&jobu,&jobv,&m2,&n2,A_,&m2,d,U_,&m2,V_,&qv2,&opt_size,&k,&info);
  if (info)
    err=BLAS_LAPACK_ERR;
  else
    if (!(work=(PRECISION*)swzMalloc((k=(int)opt_size)*sizeof(PRECISION))))
      err=MEM_ERR;
    else
      {
    gesvd(&jobu,&jobv,&m2,&n2,A_,&m2,d,U_,&m2,V_,&qv2,work,&k,&info);
    free(work);
    if (info)
      err=BLAS_LAPACK_ERR;
    else
      {
        if (U_ != U) bTranspose(U,U_,m,qu,1);
        if (V_ != V) bTranspose(V,V_,qv,n,1);
        err=NO_ERR;
      }
      }

  free(A_);
  if (U_ && (U_ != U)) free(U_);
  if (V_ && (V_ != V)) free(V_);
  return err;


/*   char jobu, jobv, jobt; */
/*   int k=-1, info, err, m_, n_, qu_, qv_, transpose; */
/*   PRECISION  *A_, *U_, *V_, *work, opt_size; */

/*   A_=(PRECISION*)swzMalloc(m*n*sizeof(PRECISION)); */

/*   jobu=jobv=compact ? 'S' : 'A'; */

/*   if (!U) */
/*     { */
/*       jobu='N'; */
/*       if (!V) */
/*     { */
/*       jobv='N'; */
/*       vt=transpose=1-at; */
/*     } */
/*       else */
/*     transpose=vt; */
/*       ut=1-vt; */
/*     } */
/*   else */
/*     if (!V) */
/*       { */
/*     jobv='N'; */
/*     vt=transpose=1-ut; */
/*       } */
/*     else */
/*       { */
/*     if (ut != vt) */
/*       transpose=vt; */
/*     else */
/*       transpose=1-at; */
/*       } */

/*   if (transpose) */
/*     { */
/*       jobt=jobu; */
/*       jobu=jobv; */
/*       jobv=jobt; */
/*       if (at) */
/*     bTranspose(A_,A,m,n,at); */
/*       else */
/*     memcpy(A_,A,m*n*sizeof(PRECISION)); */
/*       if (compact) */
/*     { */
/*       m_=n; */
/*       n_=m; */
/*       qu_=qv_=(m < n) ? m : n; */
/*     } */
/*       else */
/*     { */
/*       qu_=m_=n; */
/*       qv_=n_=m; */
/*     } */
/*       U_=vt ? V : (PRECISION*)swzMalloc(m_*qu_*sizeof(PRECISION)); */
/*       V_=ut ? (PRECISION*)swzMalloc(qv_*n_*sizeof(PRECISION)) : U;           */
/*     } */
/*   else */
/*     { */
/*       if (at) */
/*     memcpy(A_,A,m*n*sizeof(PRECISION)); */
/*       else */
/*     bTranspose(A_,A,m,n,at); */
/*       if (compact) */
/*     { */
/*       m_=m; */
/*       n_=n; */
/*       qu_=qv_=(m < n) ? m : n; */
/*     } */
/*       else */
/*     { */
/*       qu_=m_=m; */
/*       qv_=n_=n; */
/*     } */
/*       U_=ut ? U : (PRECISION*)swzMalloc(m_*qu_*sizeof(PRECISION)); */
/*       V_=vt ? (PRECISION*)swzMalloc(qv_*n_*sizeof(PRECISION)) : V; */
/*     } */

/*   // compute singular value decomposition */
/*   gesvd(&jobu,&jobv,&m_,&n_,A_,&m_,d,U_,&m_,V_,&qv_,&opt_size,&k,&info); */
/*   if (info || !(work=(PRECISION*)swzMalloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/*     err=info ? BLAS_LAPACK_ERR : MEM_ERR; */
/*   else */
/*     { */
/*       gesvd(&jobu,&jobv,&m_,&n_,A_,&m_,d,U_,&m_,V_,&qv_,work,&k,&info); */
/*       free(work); */
/*       if (info) */
/*     err=BLAS_LAPACK_ERR; */
/*       else */
/*     { */
/*       if (transpose) */
/*         { */
/*           if (U != V_) */
/*         bTranspose(U,V_,qv_,n_,1); */
/*           else */
/*         if (V != U_) */
/*           bTranspose(V,U_,m_,qu_,1); */
/*         } */
/*       else */
/*         { */
/*           if (U != U_) */
/*         bTranspose(U,U_,m_,qu_,1); */
/*           else */
/*         if (V != V_) */
/*           bTranspose(V,V_,qv_,n_,1); */
/*         } */
/*       err=NO_ERR; */
/*     } */
/*     } */

/*   free(A_); */

/*   if (transpose) */
/*     { */
/*       if (U != V_) */
/*     free(V_); */
/*       else */
/*     if (V != U_) */
/*       free(U_); */
/*     } */
/*   else */
/*     { */
/*       if (U != U_) */
/*     free(U_); */
/*       else */
/*     if (V != V_) */
/*       free(V_); */
/*     } */

/*   return err; */

/* #undef gesvd */
}

/*
   Assumes
    U  : array of length m*m
    d  : array of length q=min(m,n)
    V  : array of length n*n
    A  : array of length m*n
    m  : positive
    n  : positive
    ut : 0 or 1
    vt : 0 or 1
    at : 0 or 1

   Returns
     NO_ERR     : success
     MEM_ERR    : out of memory
     ITER_ERR   : maximum number of iterations (MAX_ITER) exceeded - only if
                  numerical recipe routines are used.

   Results
                        A          U          V         d
       ut  vt  at    (m x n)    (m x m)    (n x n)   diagonal     solves
       ---------------------------------------------------------------------
        0   0   0   row major  row major  row major   m x n   A = U * D * V'
        1   0   0   row major  col major  row major   m x n   A = U * D * V'
        1   1   0   row major  col major  col major   m x n   A = U * D * V'
        0   1   0   row major  row major  col major   m x n   A = U * D * V'
        0   0   1   col major  row major  row major   m x n   A = U * D * V'
        1   0   1   col major  col major  row major   m x n   A = U * D * V'
        1   1   1   col major  col major  col major   m x n   A = U * D * V'
        0   1   1   col major  row major  col major   m x n   A = U * D * V'

                        A          U          V         d
       ut  vt  at   row major  row major  row major  diagonal     solves
       ---------------------------------------------------------------------
        0   0   0     m x n      m x m      n x n     m x n   A = U * D * V'
        1   0   0     m x n      m x m      n x n     m x n   A = U'* D * V'
        1   1   0     m x n      m x m      n x n     m x n   A = U'* D * V
        0   1   0     m x n      m x m      n x n     m x n   A = U * D * V
        0   0   1     n x m      m x m      n x n     m x n   A'= U * D * V'
        1   0   1     n x m      m x m      n x n     m x n   A'= U'* D * V'
        1   1   1     n x m      m x m      n x n     m x n   A'= U'* D * V
        0   1   1     n x m      m x m      n x n     m x n   A'= U * D * V

                        A          U          V         d
       ut  vt  at   col major  col major  col major  diagonal     solves
       ---------------------------------------------------------------------
        0   0   0     n x m      m x m      n x n     m x n   A'= U'* D * V
        1   0   0     n x m      m x m      n x n     m x n   A'= U * D * V
        1   1   0     n x m      m x m      n x n     m x n   A'= U * D * V'
        0   1   0     n x m      m x m      n x n     m x n   A'= U'* D * V'
        0   0   1     m x n      m x m      n x n     m x n   A = U'* D * V
        1   0   1     m x n      m x m      n x n     m x n   A = U * D * V
        1   1   1     m x n      m x m      n x n     m x n   A = U * D * V'
        0   1   1     m x n      m x m      n x n     m x n   A = U'* D * V'


     U and V are orthogonal matrices and the elemets of d are non-negative.

   Notes
     The lapack routine is avoids unnecessary transpositions when ut == at and
     vt = 1-at.  When m=n, A can be equal to U or V.  U and V must be distinct.

     The current code uses the
*/
int bSVD(PRECISION *U, PRECISION *d, PRECISION *V, PRECISION *A, int m, int n, int ut, int vt, int at)
{
#if (PRECISION_SIZE == 4)
#define gesvd sgesvd
#else
#define gesvd dgesvd
#endif

  char jobz='A';
  int *iwork;
  PRECISION *X, *work, opt_size;

  lapack_int m2, n2, k, info;

  if (!(X=(PRECISION*)swzMalloc(m*n*sizeof(PRECISION)))) return MEM_ERR;
  memcpy(X,A,m*n*sizeof(PRECISION));
  if (!(iwork=(int*)swzMalloc(8*((m < n) ? m : n)*sizeof(int))))
    {
      free(X);
      return MEM_ERR;
    }
  k=-1;
  m2 = m;
  n2 = n;
  if (at)
    {
      memcpy(X,A,m*n*sizeof(PRECISION));
      k=-1;
      gesvd(&jobz,&jobz,&m2,&n2,X,&m2,d,U,&m2,V,&n2,&opt_size,&k,&info);
      if (info || !(work=(PRECISION*)swzMalloc((k=(int)opt_size)*sizeof(PRECISION))))
    {
      free(iwork);
      free(X);
      return info ? BLAS_LAPACK_ERR : MEM_ERR;
    }
      gesvd(&jobz,&jobz,&m2,&n2,X,&m2,d,U,&m2,V,&n2,work,&k,&info);
      if (info)
    {
      free(work);
      free(iwork);
      free(X);
      return BLAS_LAPACK_ERR;
    }
      if (!ut)
    bTransposeInPlace(U,m);
      if (vt)
    bTransposeInPlace(V,n);
    }
  else
    {
      memcpy(X,A,m*n*sizeof(PRECISION));
      k=-1;
      gesvd(&jobz,&jobz,&n2,&m2,X,&n2,d,V,&n2,U,&m2,&opt_size,&k,&info);
      if (info || !(work=(PRECISION*)swzMalloc((k=(int)opt_size)*sizeof(PRECISION))))
    {
      free(iwork);
      free(X);
      return info ? BLAS_LAPACK_ERR : MEM_ERR;
    }
      gesvd(&jobz,&jobz,&n2,&m2,X,&n2,d,V,&n2,U,&m2,work,&k,&info);
      if (info)
    {
      free(work);
      free(iwork);
      free(X);
      return BLAS_LAPACK_ERR;
    }
      if (!vt)
    bTransposeInPlace(V,n);
      if (ut)
    bTransposeInPlace(U,m);
    }
  free(work);
  free(iwork);
  free(X);
  return NO_ERR;

#undef gesvd

/*    // The following code attempts to use the divide and conquer algorithm   ansi-c*/
/* #if (PRECISION_SIZE == 4)  */
/*   #define gesvd sgesdd */
/* #else */
/*   #define gesdd dgesdd */
/* #endif */
/*   int jobz='A', k, *iwork, info; */
/*   PRECISION *X, *work, opt_size; */
/*   if (!(X=(PRECISION*)swzMalloc(m*n*sizeof(PRECISION)))) return MEM_ERR; */
/*   memcpy(X,A,m*n*sizeof(PRECISION)); */
/*   if (!(iwork=(int*)swzMalloc(8*((m < n) ? m : n)*sizeof(int)))) */
/*     { */
/*       free(X); */
/*       return MEM_ERR; */
/*     } */
/*   k=-1;   */
/*   if (at) */
/*     { */
/*       gesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,&opt_size,&k,iwork,&info); */
/*       if (info || !(work=(PRECISION*)swzMalloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/*     { */
/*       free(iwork); */
/*       free(X); */
/*       return info ? BLAS_LAPACK_ERR : MEM_ERR; */
/*     } */
/*       gesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,work,&k,iwork,&info); */
/*       if (info) */
/*     { */
/*       free(work); */
/*       free(iwork); */
/*       free(X); */
/*       return BLAS_LAPACK_ERR; */
/*     } */
/*       if (!ut)  */
/*     bTransposeInPlace(U,m); */
/*       if (vt)  */
/*     bTransposeInPlace(V,n); */
/*     } */
/*   else */
/*     { */
/*       gesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,&opt_size,&k,iwork,&info); */
/*       if (!(work=(PRECISION*)swzMalloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/*           { */
/*             free(iwork); */
/*             free(X); */
/*             return MEM_ERR; */
/*           } */
/*       gesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,work,&k,iwork,&info); */
/*       if (info) */
/*           { */
/*             free(work); */
/*       free(iwork); */
/*       free(X); */
/*       return BLAS_LAPACK_ERR; */
/*     } */
/*       if (!vt)  */
/*     bTransposeInPlace(V,n); */
/*       if (ut)  */
/*     bTransposeInPlace(U,m); */
/*     } */
/*   free(work); */
/*   free(iwork); */
/*   free(X); */
/*   return NO_ERR; */

/* #undef gesdd    */
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/***************************** QR Decompositions ******************************/
/******************************************************************************/
/*
   Assumes
     Q : array of length m*q or null pointer
     R : array of length q*n
     X : array of length m*n
     n : positive
     m : positive
     q : m or min(m,n)
    qt : 0 or 1
    rt : 0 or 1
    xt : 0 or 1

   Returns
     NO_ERR   : Success
     MEM_ERR  : Out of memory

   Results
     Finds an orthogonal matrix Q and an upper triangular matrix R
     such that

                                X = Q * R

     The matrix Q is computed only if it is not null.

                       X          Q          R
       qt  rt  xt   (m x n)    (m x q)    (q x n)    solves
       ------------------------------------------------------
        0   0   0  row major  row major  row major  X = Q * R
        1   0   0  row major  col major  row major  X = Q * R
        0   1   0  row major  row major  col major  X = Q * R
        1   1   0  row major  col major  col major  X = Q * R
        0   0   1  col major  row major  row major  X = Q * R
        1   0   1  col major  col major  row major  X = Q * R
        0   1   1  col major  row major  col major  X = Q * R
        1   1   1  col major  col major  col major  X = Q * R

     or
                    R/U         Q
       qt   rt   row major  row major   solves
       ------------------------------------------
        0    0     m x n      m x m    R = Q * U
        1    0     m x n      m x m    R = Q'* U
        0    1     n x m      m x m    R'= Q * U'
        1    1     n x m      m x m    R'= Q'* U'

     or
                    R/U         Q
       qt   rt   col major  col major   solves
       ------------------------------------------
        0    0     n x m      m x m    R'= Q'* U'
        1    0     n x m      m x m    R'= Q * U'
        0    1     m x n      m x m    R = Q'* U
        1    1     m x n      m x m    R = Q * U

   Notes
     The matrices X and R do not have to be distinct.  If X == R, then it must
     be the case that m == q and rt == xt.  The QR decomposition is formed using
     Householder matrices without pivoting.
*/
int bQR(PRECISION *Q, PRECISION *R, PRECISION *X, int m, int n, int q, int qt, int rt, int xt)
{
#if (PRECISION_SIZE == 4)
  #define geqrf sgeqrf
  #define orgqr sorgqr
  #define gelqf sgelqf
  #define orglq sorglq
#else
  #define geqrf dgeqrf
  #define orgqr dorgqr
  #define gelqf dgelqf
  #define orglq dorglq
#endif

  int i, j, k, l, p=(m < n) ? m : n;
  PRECISION *tau, *work, *ptr, opt_size;

  lapack_int m2, n2, p2, q2, lwork, info;

  if (!(tau=(PRECISION*)swzMalloc(p*sizeof(PRECISION)))) return MEM_ERR;
  if (xt)
    {
      lwork=-1;
      m2 = m;
      n2 = n;
      geqrf(&m2,&n2,X,&m2,tau,&opt_size,&lwork,&info);

      if (!(work=(PRECISION*)swzMalloc((lwork=(int)opt_size)*sizeof(PRECISION))))
    {
      free(tau);
      return MEM_ERR;
    }

      geqrf(&m2,&n2,X,&m2,tau,work,&lwork,&info);

      free(work);
      if (info)
    {
      free(tau);
      return ARG_ERR;
    }
      if (Q)
    {
      if (qt)
        ptr=Q;
      else
        if (!(ptr=(PRECISION*)swzMalloc(m*q*sizeof(PRECISION))))
          {
        free(tau);
        return MEM_ERR;
          }
      memcpy(ptr,X,m*p*sizeof(PRECISION));
      lwork=-1;

          p2 = p;
          q2 = q;
      orgqr(&m2,&q2,&p2,ptr,&m2,tau,&opt_size,&lwork,&info);

      if (!(work=(PRECISION*)swzMalloc((lwork=(int)opt_size)*sizeof(PRECISION))))
        {
          if (!qt) free(ptr);
          free(tau);
          return MEM_ERR;
        }

      orgqr(&m2,&q2,&p2,ptr,&m2,tau,work,&lwork,&info);

      free(work);
      if (!qt)
        {
          bTranspose(Q,ptr,m,q,1);
          free(ptr);
        }
      free(tau);
      if (info) return ARG_ERR;
    }
      else
    free(tau);
      if (R != X)
    if (rt)
      for (k=q*n, j=n-1; j >= 0; j--)
        {
          for (i=q-1; i > j; i--) R[--k]=0.0;
          for (l=i+j*m; i >= 0; i--) R[--k]=X[l--];
        }
    else
      for (k=q*n, i=q-1; i >= 0; i--)
        {
          for (l=i+n*m, j=n-1; j >= i; j--) R[--k]=X[l-=m];
          for ( ; j >= 0; j--) R[--k]=0.0;
        }
      else
    {
      for (j=p-1; j >= 0; j--)
        for (k=m*(j+1), i=m-1; i > j; i--) X[--k]=0.0;
    }
    }
  else
    {
      lwork=-1;
      m2 = m;
      n2 = n;
      gelqf(&n2,&m2,X,&n2,tau,&opt_size,&lwork,&info);

      if (!(work=(PRECISION*)swzMalloc((lwork=(int)opt_size)*sizeof(PRECISION))))
    {
      free(tau);
      return MEM_ERR;
    }

      gelqf(&n2,&m2,X,&n2,tau,work,&lwork,&info);

      free(work);
      if (info)
    {
      free(tau);
      return ARG_ERR;
    }
      if (Q)
    {
      if (!qt)
        ptr=Q;
      else
        if (!(ptr=(PRECISION*)swzMalloc(m*q*sizeof(PRECISION))))
          {
        free(tau);
        return MEM_ERR;
          }
      if (q == n)
        memcpy(ptr,X,m*n*sizeof(PRECISION));
      else
        if (m < n)
          for (k=q*m, j=m-1; j >= 0; j--)
        for (l=p+j*n, i=p-1; i >= 0; i--)
          ptr[--k]=X[--l];
        else
          for (l=n*m, j=m-1; j >= 0; j--)
        for (k=p+j*q, i=p-1; i >= 0; i--)
          ptr[--k]=X[--l];
      lwork=-1;

          m2 = m;
          p2 = p;
          q2 = q;
      orglq(&q2,&m2,&p2,ptr,&q2,tau,&opt_size,&lwork,&info);

      if (!(work=(PRECISION*)swzMalloc((lwork=(int)opt_size)*sizeof(PRECISION))))
        {
          if (!qt) free(ptr);
          free(tau);
          return MEM_ERR;
        }

      orglq(&q2,&m2,&p2,ptr,&q2,tau,work,&lwork,&info);

      free(work);
      if (qt)
        {
          bTranspose(Q,ptr,q,m,1);
          free(ptr);
        }
      free(tau);
      if (info) return ARG_ERR;
    }
      else
    free(tau);
      if (R != X)
    if (rt)
      for (k=n*q, i=n-1; i >= 0; i--)
        {
          for (j=q-1; j > i; j--) R[--k]=0.0;
          for (l=i+j*n; j >= 0; l-=n, j--) R[--k]=X[l];
        }
    else
          for (k=n*q-1, j=q-1; j >= 0; j--)
        {
          for (i=n-1; i >= j; k--, i--) R[k]=X[k];
          for ( ; i >= 0; k--, i--) R[k]=0.0;
        }
      else
    {
      for (i=p-1; i >= 0; i--)
        for (k=i+m*n, j=m-1; j > i; j--) X[k-=n]=0.0;
    }
    }
  return NO_ERR;

#undef geqrf
#undef orgqr
#undef gelqf
#undef orglq
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/********************** Generalized Schur Decompositions **********************/
/******************************************************************************/
/*
   Assumes
    Q       : array of length n*n or null
    Z       : array of length n*n or null
    S       : array of length n*n
    T       : array of length n*n
    A       : array of length n*n
    B       : array of length n*n
    n       : positive integer
    qt      : 0 or 1
    zt      : 0 or 1
    st      : 0 or 1
    tt      : 0 or 1
    at      : 0 or 1
    bt      : 0 or 1
    alpha_i : array of length n or null
    alpha_i : array of length n or null
    beta    : array of length n or null

   Returns
     NO_ERR          : success
     MEM_ERR         : out of memory
     BLAS_LAPACK_ERR : blas or lapack error

   Results
     Finds orthogonal matrices Q and Z, an block upper triangular matrix S with
     1 x 1 or 2 x 2 blocks along the diagonal, an upper triangular matrix T such
     that

         A = Q*S*Z'   and    B = Q*T*Z

     If either Q or Z is null, then it is not returned.

   Notes
     The flags qt, zt, st, tt, at, and bt control the format of the matrices Q,
     Z, S, T, A, and B.  A value of 1 indicate column major format and a value of
     0 indictes row major format.
*/
int bQZ_real(PRECISION *Q, PRECISION *Z, PRECISION *S, PRECISION *T, PRECISION *A, PRECISION *B, int n, int qt, int zt, int st, int tt, int at, int bt,
         PRECISION *alpha_r, PRECISION *alpha_i, PRECISION *beta)
{
#if (PRECISION_SIZE == 4)
  #define gges sgges
#else
  #define gges dgges
#endif

  char jobvsl, jobvsr, sort='N';
  int rtrn;
  PRECISION *work, size, *palpha_r, *palpha_i, *pbeta;

  lapack_int n2, simd, lwork, info;

  jobvsl=Q ? 'V' : 'N';
  jobvsr=Z ? 'V' : 'N';
  palpha_r=alpha_r ? alpha_r : (PRECISION*)swzMalloc(n*sizeof(PRECISION));
  palpha_i=alpha_i ? alpha_i : (PRECISION*)swzMalloc(n*sizeof(PRECISION));
  pbeta=beta ? beta : (PRECISION*)swzMalloc(n*sizeof(PRECISION));

  if (palpha_r && palpha_i && pbeta)
    {
      if (S != A)
    if (at)
      memcpy(S,A,n*n*sizeof(PRECISION));
    else
      bTranspose(S,A,n,n,0);
      else
    if (!at) bTransposeInPlace(A,n);

      if (T != B)
    if (bt)
      memcpy(T,B,n*n*sizeof(PRECISION));
    else
      bTranspose(T,B,n,n,0);
      else
    if (!bt) bTransposeInPlace(B,n);

      lwork=-1;
      n2 = n;
      gges(&jobvsl,&jobvsr,&sort,(void*)NULL,&n2,S,&n2,T,&n2,&simd,palpha_r,palpha_i,pbeta,Q,&n2,Z,&n2,&size,&lwork,(void*)NULL,&info);
      if (!info)
    if (!(work=swzMalloc((lwork=(int)size)*sizeof(PRECISION))))
      rtrn=MEM_ERR;
    else
      {
        gges(&jobvsl,&jobvsr,&sort,(void*)NULL,&n2,S,&n2,T,&n2,&simd,palpha_r,palpha_i,pbeta,Q,&n2,Z,&n2,work,&lwork,(void*)NULL,&info);
        if (!info)
          {
        if (Q && !qt) bTransposeInPlace(Q,n);
        if (Z && !zt) bTransposeInPlace(Z,n);
        if (!st) bTransposeInPlace(S,n);
        if (!tt) bTransposeInPlace(T,n);
        rtrn=NO_ERR;
          }
        else
          rtrn=BLAS_LAPACK_ERR;
        free(work);
      }
      else
    rtrn=BLAS_LAPACK_ERR;
    }
  else
    rtrn=MEM_ERR;

  if (!alpha_r && palpha_r) free(palpha_r);
  if (!alpha_i && palpha_i) free(palpha_i);
  if (!beta && pbeta) free(pbeta);

  return rtrn;

#undef gges
}

/*
   Assumes
    select  : array of length n
    QQ      : array of length n*n or null
    ZZ      : array of length n*n or null
    SS      : array of length n*n
    TT      : array of length n*n
    Q       : array of length n*n or null
    Z       : array of length n*n or null
    S       : array of length n*n
    T       : array of length n*n
    n       : positive integer
    qqt     : 0 or 1
    zzt     : 0 or 1
    sst     : 0 or 1
    ttt     : 0 or 1
    qt      : 0 or 1
    zt      : 0 or 1
    st      : 0 or 1
    tt      : 0 or 1
    alpha_i : array of length n or null
    alpha_i : array of length n or null
    beta    : array of length n or null

   Returns
     NO_ERR          : success
     MEM_ERR         : out of memory
     BLAS_LAPACK_ERR : blas or lapack error

   Results
     Finds orthogonal matrices QQ and ZZ, an block upper triangular matrix SS
     with 1 x 1 or 2 x 2 blocks along the diagonal, an upper triangular matrix TT
     such that

              Q*S*Z' = QQ*SS*ZZ'   and    Q*T*Z' = QQ*TT*ZZ'

     If either Q or QQ are null, then QQ is not computed and if either Z or ZZ is
     null, then ZZ is not computed.  The matrices S and T are multiplied by
     orthogonal matrices in such a manner that their block triangular structure
     retained and the generalized eigenvalues corresponding to value of select
     equal to one are transformed to the upper part of SS and TT.

   Notes
     The flags qqt, zzt, sst, ttt, qt, zt, st, and tt control the format of the
     matrices QQ, ZZ, SS, TT, Q, Z, S, and T.  A value of 1 indicate column major
     format and a value of 0 indictes row major format.
*/
int bReorderQZ_real(int *select, PRECISION *QQ, PRECISION *ZZ, PRECISION *SS, PRECISION *TT, PRECISION *Q, PRECISION *Z, PRECISION *S, PRECISION *T, int n,
            int qqt, int zzt, int sst, int ttt, int qt, int zt, int st, int tt, PRECISION *alpha_r, PRECISION *alpha_i, PRECISION *beta)
{
#if (PRECISION_SIZE == 4)
  #define tgsen stgsen
#else
  #define tgsen dtgsen
#endif

  int wantq, wantz, m=n, rtrn, i;
  PRECISION size, *palpha_r, *palpha_i, *pbeta, *work;

  lapack_int ijob, wantq2, wantz2, *select2, n2, m2, lwork, iwork, liwork, info;
  ijob = 0;
  liwork=1;

  if(!(select2 = (lapack_int *)swzCalloc(n, sizeof(lapack_int))))
    return MEM_ERR;

  for(i=0; i<n; i++)
    select2[i] = select[i];

  wantq=(QQ && Q) ? 1 : 0;
  wantz=(ZZ && Z) ? 1 : 0;

  palpha_r=alpha_r ? alpha_r : (PRECISION*)swzMalloc(n*sizeof(PRECISION));
  palpha_i=alpha_i ? alpha_i : (PRECISION*)swzMalloc(n*sizeof(PRECISION));
  pbeta=beta ? beta : (PRECISION*)swzMalloc(n*sizeof(PRECISION));

  if (palpha_r && palpha_i && pbeta)
    {
      if (SS != S)
    if (st)
      memcpy(SS,S,n*n*sizeof(PRECISION));
    else
      bTranspose(SS,S,n,n,0);
      else
    if (!st) bTransposeInPlace(S,n);

      if (TT != T)
    if (tt)
      memcpy(TT,T,n*n*sizeof(PRECISION));
    else
      bTranspose(TT,T,n,n,0);
      else
    if (!tt) bTransposeInPlace(T,n);

      if (wantq)
    if (QQ != Q)
      if (qt)
        memcpy(QQ,Q,n*n*sizeof(PRECISION));
      else
        bTranspose(QQ,Q,n,n,0);
    else
      if (!qt) bTransposeInPlace(Q,n);

      if (wantz)
    if (ZZ != Z)
      if (zt)
        memcpy(ZZ,Z,n*n*sizeof(PRECISION));
      else
        bTranspose(ZZ,Z,n,n,0);
    else
      if (!zt) bTransposeInPlace(Z,n);

      lwork=-1;
      wantq2 = wantq;
      wantz2 = wantz;
      n2 = n;
      m2 = m;
      tgsen(&ijob,&wantq2,&wantz2,select2,&n2,SS,&n2,TT,&n2,palpha_r,palpha_i,pbeta,QQ,&n2,ZZ,&n2,&m2,
        (PRECISION*)NULL,(PRECISION*)NULL,(PRECISION*)NULL,&size,&lwork,&iwork,&liwork,&info);
      m = m2;
      if (!info)
    if (!(work=swzMalloc((lwork=(int)size)*sizeof(PRECISION))))
      rtrn=MEM_ERR;
    else
      {
        tgsen(&ijob,&wantq2,&wantz2,select2,&n2,SS,&n2,TT,&n2,palpha_r,palpha_i,pbeta,QQ,&n2,ZZ,&n2,&m2,
          (PRECISION*)NULL,(PRECISION*)NULL,(PRECISION*)NULL,work,&lwork,&iwork,&liwork,&info);
            m = m2;
        if (!info)
          {
        if (wantq && !qqt) bTransposeInPlace(QQ,n);
        if (wantz && !zzt) bTransposeInPlace(ZZ,n);
        if (!sst) bTransposeInPlace(SS,n);
        if (!ttt) bTransposeInPlace(TT,n);
        rtrn=NO_ERR;
          }
        else
          rtrn=BLAS_LAPACK_ERR;
        free(work);
      }
      else
    rtrn=BLAS_LAPACK_ERR;
    }

  if (!alpha_r && palpha_r) free(palpha_r);
  if (!alpha_i && palpha_i) free(palpha_i);
  if (!beta && pbeta) free(pbeta);
  free(select2);

  return rtrn;

#undef tgsen
}

/*
   Assumes
    QQ      : array of length n*n or null
    ZZ      : array of length n*n or null
    SS      : array of length n*n
    TT      : array of length n*n
    Q       : array of length n*n or null
    Z       : array of length n*n or null
    S       : array of length n*n
    T       : array of length n*n
    n       : positive integer
    qqt     : 0 or 1
    zzt     : 0 or 1
    sst     : 0 or 1
    ttt     : 0 or 1
    qt      : 0 or 1
    zt      : 0 or 1
    st      : 0 or 1
    tt      : 0 or 1
    alpha_i : array of length n
    alpha_i : array of length n
    beta    : array of length n

   Returns
     NO_ERR          : success
     MEM_ERR         : out of memory
     BLAS_LAPACK_ERR : blas or lapack error

   Results
     Finds orthogonal matrices QQ and ZZ, an block upper triangular matrix SS
     with 1 x 1 or 2 x 2 blocks along the diagonal, an upper triangular matrix TT
     such that

              Q*S*Z' = QQ*SS*ZZ'   and    Q*T*Z' = QQ*TT*ZZ'

     If either Q or QQ are null, then QQ is not computed and if either Z or ZZ is
     null, then ZZ is not computed.  The matrices S and T are multiplied by
     orthogonal matrices in such a manner that their block triangular structure
     retained and the generalized eigenvalues are sorted in descending order.  So
     upon exit,

           sqrt(alpha_r[i]^2 + alpha_i[i]^2)/beta[i]
                        >= sqrt(alpha_r[i+1]^2 + alpha_i[i+1]^2)/beta[i+1]

   Notes
     The flags qqt, zzt, sst, ttt, qt, zt, st, and tt control the format of the
     matrices QQ, ZZ, SS, TT, Q, Z, S, and T.  A value of 1 indicate column major
     format and a value of 0 indictes row major format.
*/
int bSortQZ_real(int *select, PRECISION *QQ, PRECISION *ZZ, PRECISION *SS, PRECISION *TT, PRECISION *Q, PRECISION *Z, PRECISION *S, PRECISION *T, int n,
         int qqt, int zzt, int sst, int ttt, int qt, int zt, int st, int tt, PRECISION *alpha_r, PRECISION *alpha_i, PRECISION *beta)
{
#if (PRECISION_SIZE == 4)
  #define tgexc stgexc
#else
  #define tgexc dtgexc
#endif

  int wantq, wantz, rtrn, *pairs, i, j, ii, jj;
  PRECISION size, *work, *gev, small, x1, x2;

  lapack_int wantq2, wantz2, n2, i2, j2, ii2, jj2, lwork, info;

  if (n == 1) return NO_ERR;

  wantq=(QQ && Q) ? 1 : 0;
  wantz=(ZZ && Z) ? 1 : 0;

  pairs=(int*)swzMalloc(n*sizeof(int));
  gev=(PRECISION*)swzMalloc(n*sizeof(PRECISION));
  small=SQRT_MACHINE_EPSILON;

  if (pairs && gev)
    {
      if (SS != S)
    if (st)
      memcpy(SS,S,n*n*sizeof(PRECISION));
    else
      bTranspose(SS,S,n,n,0);
      else
    if (!st) bTransposeInPlace(S,n);

      if (TT != T)
    if (tt)
      memcpy(TT,T,n*n*sizeof(PRECISION));
    else
      bTranspose(TT,T,n,n,0);
      else
    if (!tt) bTransposeInPlace(T,n);

      if (wantq)
    if (QQ != Q)
      if (qt)
        memcpy(QQ,Q,n*n*sizeof(PRECISION));
      else
        bTranspose(QQ,Q,n,n,0);
    else
      if (!qt) bTransposeInPlace(Q,n);

      if (wantz)
    if (ZZ != Z)
      if (zt)
        memcpy(ZZ,Z,n*n*sizeof(PRECISION));
      else
        bTranspose(ZZ,Z,n,n,0);
    else
      if (!zt) bTransposeInPlace(Z,n);

      lwork=-1;
      j=2; i=1;
      wantq2 = wantq;
      wantz2 = wantz;
      n2 = n;
      j2 = j;
      i2 = i;
      tgexc(&wantq2,&wantz2,&n2,SS,&n2,TT,&n2,QQ,&n2,ZZ,&n2,&j2,&i2,&size,&lwork,&info);
      i = i2;
      j = j2;
      if (!info)
    if (!(work=swzMalloc((lwork=(int)size)*sizeof(PRECISION))))
      rtrn=MEM_ERR;
    else
      {
/*          // Setup pairs and gev   ansi-c*/
        for (i=n-1; i >= 0; i--)
          {
        gev[i]=sqrt(alpha_r[i]*alpha_r[i] + alpha_i[i]*alpha_i[i]);
        gev[i]=(gev[i]*small > beta[i]) ? 1.0/small : gev[i]/beta[i];
        if ((i > 0) && (SS[n*(i-1)+i] != 0.0))
          {
            i--;
            gev[i]=gev[i+1];
            pairs[i]=1;
            pairs[i+1]=-1;
          }
        else
          pairs[i]=0;
          }

        rtrn=NO_ERR;

/*          // Order generalized eigenvalues   ansi-c*/
        j=pairs[0] ? 2 : 1;
        while (j < n)
          {
        i=j;
        while ((i > 0) && (gev[i-1] < gev[j])) i-=pairs[i-1] ? 2 : 1;

        if (i != j)
          {
            ii=i+1;
            jj=j+1;
                    wantq2 = wantq;
                    wantz2 = wantz;
                    n2 = n;
                    jj2 = jj;
                    ii2 = ii;
            tgexc(&wantq2,&wantz2,&n2,SS,&n2,TT,&n2,QQ,&n2,ZZ,&n2,&jj2,&ii2,work,&lwork,&info);
                    ii = ii2;
                    jj = jj2;
            if (!info)
              if (pairs[j])
            {
              memmove(pairs+i+2,pairs+i,(j-i)*sizeof(int));
              pairs[i]=1; pairs[i+1]=-1;
              x1=gev[j];
              memmove(gev+i+2,gev+i,(j-i)*sizeof(PRECISION));
              gev[i]=gev[i+1]=x1;
              x1=alpha_r[j]; x2=alpha_r[j+1];
              memmove(alpha_r+i+2,alpha_r+i,(j-i)*sizeof(PRECISION));
              alpha_r[i]=x1; alpha_r[i+1]=x2;
              x1=alpha_i[j]; x2=alpha_i[j+1];
              memmove(alpha_i+i+2,alpha_i+i,(j-i)*sizeof(PRECISION));
              alpha_i[i]=x1; alpha_i[i+1]=x2;
              x1=beta[j]; x2=beta[j+1];
              memmove(beta+i+2,beta+i,(j-i)*sizeof(PRECISION));
              beta[i]=x1; beta[i+1]=x2;
              j+=2;
            }
              else
            {
              memmove(pairs+i+1,pairs+i,(j-i)*sizeof(int));
              pairs[i]=0;
              x1=gev[j];
              memmove(gev+i+1,gev+i,(j-i)*sizeof(PRECISION));
              gev[i]=x1;
              x1=alpha_r[j];
              memmove(alpha_r+i+1,alpha_r+i,(j-i)*sizeof(PRECISION));
              alpha_r[i]=x1;
              x1=alpha_i[j];
              memmove(alpha_i+i+1,alpha_i+i,(j-i)*sizeof(PRECISION));
              alpha_i[i]=x1;
              x1=beta[j];
              memmove(beta+i+1,beta+i,(j-i)*sizeof(PRECISION));
              beta[i]=x1;
              j+=1;
            }
            else
              {
            rtrn=BLAS_LAPACK_ERR;
            break;
              }
          }
        else
          j+=pairs[j] ? 2 : 1;

          }

        free(work);

        if (rtrn == NO_ERR)
          {
        if (wantq && !qqt) bTransposeInPlace(QQ,n);
        if (wantz && !zzt) bTransposeInPlace(ZZ,n);
        if (!sst) bTransposeInPlace(SS,n);
        if (!ttt) bTransposeInPlace(TT,n);
          }
      }
    }

  if (pairs) free(pairs);
  if (gev) free(gev);

  return rtrn;

#undef tgexc
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/************************** Cholesky Decompositions ***************************/
/******************************************************************************/
/*
   Assumes
     X : Scalar array of m*m representing an m x m symmetric matrix

   Returns
     0 (NO_ERR)  : success
     POS_DEF_ERR : X not positive definite

   Results

        u    t      X/T        solves
        ----------------------------------
        0    0   row major    X = L'* L
        1    0   row major    X = U'* U
        0    1   col major    X = L'* L
        1    1   col major    X = U'* U

     or
                     X/T
        u    t    row major     solves
        ----------------------------------
        0    0        -        X = L'* L
        1    0        -        X = U'* U
        0    1        -        X = L * L'
        1    1        -        X = U * U'

     or
                     X/T
        u    t    col major     solves
        ----------------------------------
        0    0        -        X = L * L'
        1    0        -        X = U * U'
        0    1        -        X = L'* L
        1    1        -        X = U'* U

        u    t        T         solves
        ----------------------------------   v
        0    0   L:row major   X = T'* T     0 v^t=0
        1    0   U:row major   X = T'* T     1 v^t=1
        0    1   U:col major   X = T'* T     1 v^t=0
        1    1   L:col major   X = T'* T     0 v^t=1

      or
                      X
        u    t    row major     solves
        ----------------------------------
        0    0        L        X = L'* L
        1    0        U        X = U'* U
        0    1        L        X = L * L'
        1    1        U        X = U * U'

     Upon successful exit T is upper triangular with positive diagonal and
     satisfies X = T' * T.  T overwrites X.

   Notes
     Failure usually indicates X is not positive definite.  Only half of X
     is accessed.
*/
int bCholesky(PRECISION *X, int m, int u, int t)
{
 int i, j, k, b;
 PRECISION scale, *pX, *pXi, *pXj;

 if (u^t)
   if (t)
     for (i=m-1, pXi=X+i*m; i >= 0; pXi-=m, i--)
      {
       for (j=m-1, pX=X+i+j*m; j > i; pX-=m, j--) *pX=0.0;

       for (k=i+1; k < m; k++) *pX-=pXi[k]*pXi[k];

       if (*pX <= 0.0) return POSDEF_ERR;
       scale=1.0/(*pX=sqrt(*pX));

       pXj=pXi;
       for (j--; j >= 0; j--)
        {
         pX-=m;
         pXj-=m;
         for (k=i+1; k < m; k++) *pX-=pXi[k]*pXj[k];
         *pX*=scale;
        }
      }
    else
     for (i=0, pXi=X; i < m; pXi++, i++)
      {
       for (j=0, pX=X+i*m; j < i; pX++, j++) *pX=0.0;

       for (k=(i-1)*m; k >= 0; k-=m) *pX-=pXi[k]*pXi[k];

       if (*pX <= 0.0) return POSDEF_ERR;
       scale=1.0/(*pX=sqrt(*pX));

       pXj=pXi;
       for (j++; j < m; j++)
        {
         pXj++;
         pX++;
         for (k=(i-1)*m; k >= 0; k-=m) *pX-=pXi[k]*pXj[k];
         *pX*=scale;
        }
      }
  else
   if (t)
    for (i=0, pXi=X; i < m; pXi+=m, i++)
      {
       for (j=0, pX=X+i; j < i; pX+=m, j++) *pX=0.0;

       for (k=i-1; k >= 0; k--) *pX-=pXi[k]*pXi[k];

       if (*pX <= 0.0) return POSDEF_ERR;
       scale=1.0/(*pX=sqrt(*pX));

       pXj=pXi;
       for (j++; j < m; j++)
        {
         pX+=m;
         pXj+=m;
         for (k=i-1; k >= 0; k--) *pX-=pXi[k]*pXj[k];
         *pX*=scale;
        }
      }
    else
     for (b=m*m, i=m-1, pXi=X+i; i >= 0; pXi--, i--)
      {
       for (j=m-1, pX=X+i*m+j; j > i; pX--, j--) *pX=0.0;

       for (k=(i+1)*m; k < b; k+=m) *pX-=pXi[k]*pXi[k];

       if (*pX <= 0.0) return POSDEF_ERR;
       scale=1.0/(*pX=sqrt(*pX));

       pXj=pXi;
       for (j--; j >= 0; j--)
        {
         pXj--;
         pX--;
         for (k=(i+1)*m; k < b; k+=m)  *pX-=pXi[k]*pXj[k];
         *pX*=scale;
        }
      }
 return NO_ERR;
}

/*
   Assumes
       x     : array of length m*r*n*s
       y     : array of length m*n
       z     : array of length r*s
    m,n,r,s  : positive
    xt,yt,zt : 0 or 1

   Returns
     NO_ERR     : success

   Results
                        x          y          z
       xt  yt  zt   (mr x ns)   (m x n)    (r x s)      computes
       ---------------------------------------------------------------------
        0   0   0   row major  row major  row major   x = y tensor z
        1   0   0   col major  row major  row major   x = y tensor z
        1   1   0   col major  col major  row major   x = y tensor z
        0   1   0   row major  col major  row major   x = y tensor z
        0   0   1   row major  row major  col major   x = y tensor z
        1   0   1   col major  row major  col major   x = y tensor z
        1   1   1   col major  col major  col major   x = y tensor z
        0   1   1   row major  col major  col major   x = y tensor z
*/
int bMatrixTensor(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int r, int s, int xt, int yt, int zt)
{
  int iy, jy, iz, jz, k, l, stride;
  PRECISION t, *pz=z+r*s-1;
  if (xt)
    if (zt)
      {
    stride=m*r;
    for (iy=m-1; iy >= 0; iy--)
      for (jy=n-1; jy >= 0; jy--)
        {
          t=y[yt ? iy+m*jy : n*iy+jy];
          l=(iy+1)*r-1 + ((jy+1)*s-1)*stride;
          z=pz;
          for (jz=s-1; jz >= 0; l-=stride, jz--)
        for (iz=r-1, k=l; iz >= 0; z--, k--, iz--)
          x[k]=t*(*z);
        }
      }
    else
      {
    stride=m*r;
    for (iy=m-1; iy >= 0; iy--)
      for (jy=n-1; jy >= 0; jy--)
        {
          t=y[yt ? iy+m*jy : n*iy+jy];
          l=(iy+1)*r-1 + ((jy+1)*s-1)*stride;
          z=pz;
          for (iz=r-1; iz >= 0; l--, iz--)
        for (jz=s-1, k=l; jz >= 0; z--, k-=stride, jz--)
          x[k]=t*(*z);
        }
      }
  else
    if (zt)
      {
    stride=n*s;
    for (iy=m-1; iy >= 0; iy--)
      for (jy=n-1; jy >= 0; jy--)
        {
          t=y[yt ? iy+m*jy : n*iy+jy];
          l=((iy+1)*r-1)*stride + (jy+1)*s-1;
          z=pz;
          for (jz=s-1; jz >= 0; l--, jz--)
        for (iz=r-1, k=l; iz >= 0; z--, k-=stride, iz--)
          x[k]=t*(*z);

        }
      }
    else
      {
    stride=n*s;
    for (iy=m-1; iy >= 0; iy--)
      for (jy=n-1; jy >= 0; jy--)
        {
          t=y[yt ? iy+m*jy : n*iy+jy];
          l=((iy+1)*r-1)*stride + (jy+1)*s-1;
          z=pz;
          for (iz=r-1; iz >= 0; l-=stride, iz--)
        for (jz=s-1, k=l; jz >= 0; z--, k--, jz--)
          x[k]=t*(*z);

        }
      }
  return NO_ERR;
}

int bVectorTensor(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n)
{
  int j, k;
  PRECISION s;
  for (x+=m*n-1, j=m-1; j >= 0; j--)
    for (s=y[j], k=n-1; k >= 0; x--, k--)
      *x=s*z[k];
  return NO_ERR;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
