  
#include "bmatrix.h"
#include "dw_error.h"

#include <math.h>
#include <stdlib.h>
#include <memory.h>

/********************/
#include "blas_lapack.h"
#define USE_BLAS_LAPACK
/********************/

#include "modify_for_mex.h"

/********************
#include "mkl.h"
#define USE_BLAS_LAPACK
/********************/

/********************
#define USE_INLINE
/********************/

static PRECISION pythag(PRECISION a, PRECISION b);
static int bSVD_NumericalRecipes(PRECISION *U, PRECISION *d, PRECISION *V, int m, int n);
static int bQR_NumericalRecipes(PRECISION *Q, PRECISION *R, int m, int n);

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
#if defined USE_INLINE
  __asm {
        fninit
        mov  edi,x
        mov  esi,y
        mov  eax,n
        dec  eax       
a1:     fld  PRECISION_WORD ptr [esi+PRECISION_SIZE*eax]
        fabs
        fstp PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]                     
        dec  eax
        jge  a1 
       }
 return NO_ERR;
#else
 while (--n >= 0) x[n]=fabs(y[n]);
 return NO_ERR;
#endif
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
#if defined USE_INLINE
 int rtrn;
 __asm {
        fninit
        mov     edi,x
        mov     esi,y
        mov     ebx,z
        mov     eax,n
        jmp     a2   

a1:     fld     PRECISION_WORD ptr [esi+PRECISION_SIZE*eax]
        fadd    PRECISION_WORD ptr [ebx+PRECISION_SIZE*eax]
        fstp    PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]   
        
a2:     dec     eax
        jge     a1
       
        /* check coprocessor for errors */
	fnstsw  ax                                               
	and     ax,0x000D
	je      a3
	mov     rtrn,NO_ERR
        jmp     a4
a3:     mov     rtrn,FLOAT_ERR
a4:
       }      
 return rtrn;
#else
 while (--n >= 0) x[n]=y[n]+z[n];
 return NO_ERR;
#endif
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
#if defined USE_INLINE
 int rtrn;
 __asm {
        fninit
        mov     edi,x
        mov     esi,y
        mov     ebx,z
        mov     eax,n
        jmp     a2   

a1:     fld     PRECISION_WORD ptr [esi+PRECISION_SIZE*eax]
        fsub    PRECISION_WORD ptr [ebx+PRECISION_SIZE*eax]
        fstp    PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]   
        
a2:     dec     eax
        jge     a1
       
        /* check coprocessor for errors */
	fnstsw  ax                                               
	and     ax,0x000D
	je      a3
	mov     rtrn,NO_ERR
        jmp     a4
a3:     mov     rtrn,FLOAT_ERR
a4:
       }      
 return rtrn;
#else
 while (--n >= 0) x[n]=y[n]-z[n];
 return NO_ERR;
#endif
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
#ifdef USE_BLAS_LAPACK
  int inc=1;
  #if (PRECISION_SIZE == 8)
    daxpy(&m,&a,y,&inc,x,&inc);
  #else
    saxpy(&m,&a,y,&inc,x,&inc);
  #endif
#else
  while (--m >= 0) x[m]+=a*y[m];
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
#if defined USE_INLINE
  __asm {
          fninit
          mov  edi,x
          mov  esi,y
          mov  eax,n
          dec  eax
          fld  s
          a1:     fld  PRECISION_WORD ptr [esi+PRECISION_SIZE*eax]
          fmul st,st(1)
          fstp PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]                     
          dec  eax
          jge  a1
        }
#else
  while (--n >= 0) x[n]=s*y[n];
#endif
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
#if defined USE_INLINE
  /*
     ebx = uj
     ecx = vi
     edx = x
     esi = pu
     edi = pv
  */
 int i, j, ui, vj, puj, pvi;
 __asm {
         fninit
         cmp  yt,0
         je   short dest1                                   // if (yt)

         mov  ui,PRECISION_SIZE                             //     ui=PRECISION_SIZE
                                              
         mov  ebx,m         
         shl  ebx,PRECISION_SHIFT                           //     uj=m*PRECISION_SIZE 

         jmp  short dest2                                   //   else 

dest1:   mov  eax,p                                         
         shl  eax,PRECISION_SHIFT
         mov  ui,eax                                        //     ui=p*PRECISION_SIZE 

         mov  ebx,PRECISION_SIZE                            //     uj=PRECISION_SIZE   

dest2:   cmp  zt,0
         je   short dest3                                   // if (zt) 

         mov  ecx,1                                         //     vi=1
                                           
         mov  eax,p  
         shl  eax,PRECISION_SHIFT                                
         mov  vj,eax                                        //     vj=p*PRECISION_SIZE

         jmp  short dest4                                   //   else

dest3:   mov  ecx,n                                         //     vi=n
                                     
         mov  vj,PRECISION_SIZE                             //     vj=PRECISION_SIZE
        
dest4:   mov  eax,p
         dec  eax
         imul eax,ecx
         mov  pvi,eax                                       // pvi=(p-1)*vi

         mov  edx,m
         imul edx,n
         dec  edx
         shl  edx,PRECISION_SHIFT
         add  edx,x                                         // x+=(m*n-1) 
        
         cmp  xt,0
         je   short dest5                                   // if (xt)

         mov  eax,p
         dec  eax
         imul eax,ebx
         add  y,eax                                         //     y+=(p-1)*uj
         
         mov  eax,n
         dec  eax
         imul eax,vj                                        //     j=(n-1)*vj
outer_1: mov  j,eax                                     

         mov  edi,z
         add  edi,eax                                       //     pv=z+j
       
         mov  eax,m
         dec  eax
         imul eax,ui                                        //     i=(m-1)*ui
inner_1: mov  i,eax                                         

         mov  esi,y                   
         add  esi,eax                                       //     pu=y+i

         mov  eax,pvi                                       //     k=pvi
          
         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]   //     st=(*pu)*pv[k];  

         sub  eax,ecx                                       
         jl   short dest2_1                                 //     (k-=vi) >= 0

dest1_1: sub  esi,ebx                                       //     pu-=uj

         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]
         fadd                                               //     st+=(*pu)*pv[k]

         sub  eax,ecx                                       
         jge  short dest1_1                                 //     (k-=vi) >= 0

dest2_1: fstp PRECISION_WORD ptr [edx]                      //     *x=st

         sub  edx,PRECISION_SIZE                            //     x--
         
         mov  eax,i
         sub  eax,ui
         jge  short inner_1                                 //     (i=i-ui) >= 0

         mov  eax,j
         sub  eax,vj
         jge  short outer_1                                 //     (j=j-vj) >= 0
        
         jmp  short dest6                                   //   else

dest5:   mov  eax,p
         dec  eax
         imul eax,ebx
         mov  puj,eax                                       //     puj=(p-1)*uj;

         mov  eax,m
         dec  eax
         imul eax,ui                                        //     i=(m-1)*ui
outer_0: mov  i,eax

         mov  esi,y
         add  esi,eax                                       //     pu=y+i

         mov  eax,n
         dec  eax
         imul eax,vj                                        //     j=(n-1)*vj
inner_0: mov  j,eax

         mov  edi,z
         add  edi,eax                                       //     pv=z+j

         add  esi,puj
         mov  eax,pvi                                       //     pu+=puj
            
         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]   //     st=(*pu)*pv[k] 

         sub  eax,ecx                                       
         jl   short dest2_0                                 //     (k-=vi) >= 0

dest1_0: sub  esi,ebx                                       //     pu-=uj

         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]
         fadd                                               //     st+=(*pu)*pv[k]

         sub  eax,ecx
         jge  short dest1_0                                 //     (k-=vi) >= 0

dest2_0: fstp PRECISION_WORD ptr [edx]                      //     *x=st

         sub  edx,PRECISION_SIZE                            //     x--

         mov  eax,j
         sub  eax,vj
         jge  short inner_0                                 //     (j-=vj) >= 0

         mov  eax,i
         sub  eax,ui
         jge  short outer_0                                 //     (i-=ui) >= 0
dest6:
       }
 return NO_ERR;
#elif defined USE_BLAS_LAPACK
 int transy, transz, dy, dz;
 PRECISION beta=0.0, alpha=1.0;
#if PRECISION_SIZE == 4
 if (xt) 
   {
    if (yt) {transy='N'; dy=m;} else {transy='T'; dy=p;}   
    if (zt) {transz='N'; dz=p;} else {transz='T'; dz=n;}   
    sgemm(&transy,&transz,&m,&n,&p,&alpha,y,&dy,z,&dz,&beta,x,&m); 
   }
  else
   {
    if (yt) {transy='T'; dy=m;} else {transy='N'; dy=p;}   
    if (zt) {transz='T'; dz=p;} else {transz='N'; dz=n;}   
    sgemm(&transz,&transy,&n,&m,&p,&alpha,z,&dz,y,&dy,&beta,x,&n);  
   }
#else
 if (xt) 
   {
    if (yt) {transy='N'; dy=m;} else {transy='T'; dy=p;}   
    if (zt) {transz='N'; dz=p;} else {transz='T'; dz=n;}   
    dgemm(&transy,&transz,&m,&n,&p,&alpha,y,&dy,z,&dz,&beta,x,&m); 
   }
  else
   {
    if (yt) {transy='T'; dy=m;} else {transy='N'; dy=p;}   
    if (zt) {transz='T'; dz=p;} else {transz='N'; dz=n;}   
    dgemm(&transz,&transy,&n,&m,&p,&alpha,z,&dz,y,&dy,&beta,x,&n);  
   }
#endif
 return NO_ERR;
#else
 int i, j, k, ui, uj, vi, vj, puj, pvi;
 PRECISION *pu, *pv;
 if (yt) { ui=1; uj=m; } else { ui=p; uj=1; }
 if (zt) { vi=1; vj=p; } else { vi=n; vj=1; }
 pvi=(p-1)*vi;
 x+=(m*n-1);
 if (xt)
   {
    y+=(p-1)*uj;
    for (j=(n-1)*vj; j >= 0; j-=vj)
     for (pv=z+j, i=(m-1)*ui; i >= 0; i-=ui)
      {
       k=pvi;
       pu=y+i;
       *x=(*pu)*pv[k];
       while ((k-=vi) >= 0) *x+=(*(pu-=uj))*pv[k];
       x--;
      }
   }
  else
   {
    puj=(p-1)*uj;
    for (i=(m-1)*ui; i >= 0; i-=ui)
     for (pu=y+i, j=(n-1)*vj; j >= 0; j-=vj)
      {
       k=pvi;
       pu+=puj;
       pv=z+j;
       *x=(*pu)*pv[k];
       while ((k-=vi) >= 0) *x+=(*(pu-=uj))*pv[k];
       x--;
      }
   }
 return NO_ERR;
#endif
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
#if defined USE_BLAS_LAPACK
 PRECISION *y;
 int i, info;
 if (xt)
   {
#if PRECISION_SIZE == 4
     sgetrf(&m,&n,x,&m,p,&info);
#else
     dgetrf(&m,&n,x,&m,p,&info);
#endif
   }
 else
   {
     if (!( y=(PRECISION*)malloc(m*n*sizeof(PRECISION)))) return MEM_ERR;
     bTranspose(y,x,m,n,0); 

#if PRECISION_SIZE == 4
     sgetrf(&m,&n,y,&m,p,&info);
#else
     dgetrf(&m,&n,y,&m,p,&info);
#endif

     bTranspose(x,y,m,n,1);
     free(y);
   }
 for (i=(m < n) ? m-1 : n-1; i >= 0; i--) p[i]--;
 return (info < 0) ? SING_ERR : NO_ERR;
#else
 int a, b, c, d, i, j, k, q=(m < n) ? m : n, imax, rtrn=NO_ERR;
 PRECISION big, tmp;
 if (xt)
   {
     for (j=0; j < q; j++)
       {
	 if (j == 0)
	   {
	     /* Find pivot */
	     for (big=fabs(x[imax=0]), i=1; i < m; i++)
	       if (fabs(x[i]) > big)
		 {
		   big=fabs(x[i]);
		   imax=i;
		 }
	   }
	 else
	   {
             /* Perform stored row operations */
             for (i=1; i <= j; i++)
	       {
		 tmp=x[i+j*m];
		 for (c=(i-1)+j*m, b=i+(i-1)*m; b >= 0; c--, b-=m) tmp-=x[b]*x[c];
		 x[i+j*m]=tmp;
	       }

	     /* Perform stored row operations and find pivot */
	     for (big=fabs(tmp), imax=j; i < m; i++)
	       {
		 tmp=x[i+j*m];
		 for (c=(j-1)+j*m, b=i+(j-1)*m; b >= 0; c--, b-=m) tmp-=x[b]*x[c];
		 x[i+j*m]=tmp;
		 if (fabs(tmp) > big)
		   {
		     big=fabs(tmp);
		     imax=i;
		   }
	       }
	   }

	 /* Interchange rows if necessary */
	 if (j != imax)
	   {
	     p[j]=imax;
	     for (a=imax+(n-1)*m, b=j+(n-1)*m; b >= 0; a-=m, b-=m)
	       {
		 tmp=x[a];
		 x[a]=x[b];
		 x[b]=tmp;
	       }
	   }
	 else
	   p[j]=j;

	 /* Is pivot zero? */
	 if (x[a=j+j*m] != 0.0)
	   for (tmp=1.0/x[a], b=m-1+j*m; b > a; b--) x[b]*=tmp;
	 else
	   rtrn=SING_ERR;
       }

     /* Perform stored row operations */
     for ( ; j < n; j++)
       for (i=1; i < m; i++)
	 {
	   tmp=x[i+j*m];
	   for (c=(i-1)+j*m, b=i+(i-1)*m; b >= 0; c--, b-=m) tmp-=x[b]*x[c];
	   x[i+j*m]=tmp;
	 }
 
     return rtrn;
   }
 else
   { 
     for (j=0; j < q; j++)
       {
	 if (j == 0)
	   {
	     /* Find pivot */
	     for (big=fabs(x[imax=0]), a=n, i=1; i < m; a+=n, i++)
	       if (fabs(x[a]) > big)
		 {
		   big=fabs(x[a]);
		   imax=i;
		 }
	   }
	 else
	   {
	     /* Perform stored row operations */
	     for (d=n, a=j+n, i=1; i <= j; d+=n, a+=n, i++)
	       {
		 tmp=x[a];
		 for (b=d, c=j, k=i; k > 0; b++, c+=n, k--) tmp-=x[b]*x[c];
		 x[a]=tmp;
	       }

	     /* Perform stored row operations and find pivot */
	     for (big=fabs(tmp), imax=j; i < m; d+=n, a+=n, i++)
	       {
		 tmp=x[a];
		 for (b=d, c=j, k=j; k > 0; b++, c+=n, k--) tmp-=x[b]*x[c];
		 x[a]=tmp;
		 if (fabs(tmp) > big)
		   {
		     big=fabs(tmp);
		     imax=i;
		   }
	       }
	   }

	 /* Interchange rows if necessary */
	 if (j != imax)
	   {
	     p[j]=imax;

	     for (k=n-1, a=imax*n+k, b=j*n+k; k >= 0; a--, b--, k--)
	       {
		 tmp=x[a];
		 x[a]=x[b];
		 x[b]=tmp;
	       }
	   }
	 else
	   p[j]=j;

	 /* Is pivot zero? */
	 if (x[a=j*n+j] != 0.0)
	   for (tmp=1.0/x[a], a+=n, i=m*n; a < i; a+=n) x[a]*=tmp;
	 else
	   rtrn=SING_ERR;
       }

     /* Perform stored row operations */
     for ( ; j < n; j++)
       for (d=n, a=j+n, i=1; i < m; d+=n, a+=n, i++)
	 {
	   tmp=x[a];
	   for (b=d, c=j, k=i; k > 0; b++, c+=n, k--) tmp-=x[b]*x[c];
	   x[a]=tmp;
	 }

     return rtrn;
   }
#endif
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
#if defined USE_INLINE
 int i, j, bi, bj, xi, xj, mbi;
 PRECISION *pxx;
 for (j=m+1, i=m*m-1; i >= 0; i-=j) if (x[i] == 0.0) return SING_ERR;
 if (xt) { xi=1; xj=m; } else { xi=m; xj=1; }
 if (bt) { bi=1; bj=m; } else { bi=n; bj=1; }
 mbi=(m-1)*bi;
 if (u)
   for (x+=(m-1)*xj, j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
    for (pxx=x+(m-1)*xi, i=(m-1)*bi; i >= 0; pxx-=xi, i-=bi)
     __asm {
            fninit
            mov  esi,pxx
            mov  eax,i
            mov  edi,eax
            shl  edi,PRECISION_SHIFT
            add  edi,b

            fld  PRECISION_WORD ptr [edi]                      //     st=*pb

            sub  eax,mbi
            neg  eax
            je   short dest2_1                                 //     (k=mbi-i) > 0

            mov  ebx,xj
            shl  ebx,PRECISION_SHIFT
            mov  ecx,bi

dest1_1:    fld  PRECISION_WORD ptr [esi]
            fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]
            fsub                                               //     st-=(*px)*pb[k]

            sub  esi,ebx                                       //     px-=bi

            sub  eax,ecx
            jg   short dest1_1                                 //     (k-=bi) > 0

dest2_1:    fld  PRECISION_WORD ptr [esi]
            fdiv
            fstp PRECISION_WORD ptr [edi]                      //     *pb=st/(*px)
           }
  else
   for (j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
    for (pxx=x, i=0; i <= mbi; pxx+=xi, i+=bi)
     __asm {
            fninit
            mov  esi,pxx
            mov  eax,i
            mov  edi,eax
            shl  edi,PRECISION_SHIFT
            add  edi,b

            fld  PRECISION_WORD ptr [edi]                      //     st=*pb

            neg  eax
            je   short dest2_0                                 //     (k=-i) < 0

            mov  ebx,xj
            shl  ebx,PRECISION_SHIFT
            mov  ecx,bi

dest1_0:    fld  PRECISION_WORD ptr [esi]
            fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]
            fsub                                               //     st-=(*px)*pb[k]

            add  esi,ebx                                       //     px+=bi

            add  eax,ecx
            jl   short dest1_0                                 //     (k+=bi) < 0

dest2_0:    fld  PRECISION_WORD ptr [esi]
            fdiv
            fstp PRECISION_WORD ptr [edi]                      //     *pb=st/(*px)
           }
  return NO_ERR;
#else
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
#endif
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
#if defined USE_INLINE
 int i, j, bi, bj, xi, xj, mbi;
 PRECISION *pxx;
 if (xt) { xi=1; xj=m; } else { xi=m; xj=1; }
 if (bt) { bi=1; bj=m; } else { bi=n; bj=1; }
 mbi=(m-1)*bi;
 if (u)
   for (x+=(m-1)*xj, j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
    for (pxx=x+(m-1)*xi, i=(m-1)*bi; i >= 0; pxx-=xi, i-=bi) 
     __asm {                                                   //                       eax    ebx    ecx    edx    esi    edi     st(0)     st(1) 
            mov  eax,i                                         // k=i                    k
            mov  edi,eax                                       // pb=i                   k                                  pb                                              
            sub  eax,mbi                                       // k=i-mbi                k                                  pb                                                  
            neg  eax                                           // k=mbi-i                k                                  pb                                               
            je   short a2                                      // jump if k > 0          k                                  pb

            shl  edi,PRECISION_SHIFT                           // pb=i*PRECISION_SIZE    k                                  pb                                                                                                
            add  edi,b                                         // pb=b+i*PRECISION_SIZE  k                                  pb                                                             
            mov  esi,pxx                                       // px=pxx                 k                           px     pb                                                                                                 
            mov  ebx,xj                                        // ebx=xj                 k      xj                   px     pb                                                                                   
            shl  ebx,PRECISION_SHIFT                           // ebx=xj*PRECISION_SIZE  k      xj                   px     pb                                               
            mov  ecx,bi                                        // ecx=bi                 k      xj     bi            px     pb                                                    

            fld  PRECISION_WORD ptr [edi]                      // st=pb[0]               k      xj     bi            px     pb     pb[0]

a1:                                                            // assumes                k      xj     bi            px     pb     pb[0]
            fld  PRECISION_WORD ptr [esi]                      // st=*px                 k      xj     bi            px     pb      *px      pb[0]                                                                          
            fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]   // st=pb[k]*(*px)         k      xj     bi            px     pb   pb[k]*(*px) pb[0]                                                                                                 
            fsub                                               // st=pb[0]-(*px)*pb[k]   k      xj     bi            px     pb     pb[0]                                                        
         
            sub  esi,ebx                                       // px-=xj                 k      xj     bi            px     pb     pb[0]
         
            sub  eax,ecx                                       // k-=bi                  k      xj     bi            px     pb     pb[0]                                                                                    
            jg   short a1                                      // jump if k > 0          k      xj     bi            px     pb     pb[0]

            fstp PRECISION_WORD ptr [edi]                      // pb[0]=st               k      xj     bi            px     pb    
a2:
           } 
  else
   for (j=(n-1)*bj, b+=j; j >= 0; b-=bj, j-=bj)
    for (pxx=x, i=0; i <= mbi; pxx+=xi, i+=bi)  
     __asm {
            mov  eax,i
            mov  edi,eax
            neg  eax
            je   short dest2_0                                 // (k=-i) < 0
                               
            mov  esi,pxx                                        
            shl  edi,PRECISION_SHIFT
            add  edi,b 
            mov  ebx,xj
            shl  ebx,PRECISION_SHIFT
            mov  ecx,bi

            fld  PRECISION_WORD ptr [edi]                      // st=*pb

dest1_0:    fld  PRECISION_WORD ptr [esi]
            fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]
            fsub                                               // st-=(*px)*pb[k]
         
            add  esi,ebx                                       // px+=bi
         
            add  eax,ecx                                       
            jl   short dest1_0                                 // (k+=bi) < 0  

            fstp PRECISION_WORD ptr [edi]                      // *pb=st  
dest2_0:
           }
#else
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
#endif
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
   results
     returns sqrt(a^2 + b^2)
*/
static PRECISION pythag(PRECISION a, PRECISION b)
{
 PRECISION absa=fabs(a), absb=fabs(b), quotient;

 if (absa > absb)
  {
   quotient=absb/absa;
   return absa*sqrt(1.0+quotient*quotient);
  }

 if (absb == 0) return 0.0;

 quotient=absa/absb;
 return absb*sqrt(1.0+quotient*quotient);
}

/*
   Assumes
     U : m x n matrix in row major format
     d : n-vector
     V : n x n matrix in row major format

   Returns
     0 (NO_ERR) : success
     MEM_ERR    : out of memory
     ITER_ERR   : maximum number of iterations (MAX_ITER) exceeded

   Results
     Finds U, V and d such that U = U * diag(d) * V'.  U and V are orthogonal
     matrices and the elemets of d are non-negative.  U is orthogonal in the
     sense that U'*U and U*U' are diagonal matrices with ones and zeros along
     the diagonal.

   Notes
     Based on Numerical Recipes in C routines.
*/
#define MAX_ITER 100
static int bSVD_NumericalRecipes(PRECISION *U, PRECISION *d, PRECISION *V, int m, int n)
{
 int flag, i, its, j, jj, k, l, nm;
 PRECISION  anorm, c, f, g, h, s, scale, x, y, z, *rv1, tmp;

 rv1=(PRECISION*)malloc(n*sizeof(PRECISION));
 if (!rv1) return MEM_ERR;
 g=scale=anorm=0.0;

 for (i=0; i < n; i++)
  {
   l=i+1;
   rv1[i]=scale*g;
   g=s=scale=0.0;

   if (i < m)
    {
     for (k=i; k < m; k++) scale+=fabs(U[k*n+i]);

     if (scale)
      {
       for (k=i; k < m; k++)
        {
         U[k*n+i]/=scale;
         s+=U[k*n+i]*U[k*n+i];
        }

       f=U[i*n+i];
       g=(f >= 0) ? -sqrt(s) : sqrt(s);
       h=f*g-s;
       U[i*n+i]=f-g;

       for (j=l; j < n; j++)
        {
         for (s=0.0, k=i; k < m; k++) s+=U[k*n+i]*U[k*n+j];
         f=s/h;
         for (k=i; k < m; k++) U[k*n+j]+=f*U[k*n+i];
        }

       for (k=i; k < m; k++) U[k*n+i]*=scale;
      }
    }

   d[i]=scale *g;
   g=s=scale=0.0;
   if (i < m && i != n-1)
    {
     for (k=l; k < n; k++) scale+=fabs(U[i*n+k]);

     if (scale)
      {
       for (k=l; k < n; k++)
        {
         U[i*n+k]/=scale;
         s+=U[i*n+k]*U[i*n+k];
        }

       f=U[i*n+l];
       g=(f >= 0) ? -sqrt(s) : sqrt(s);
       h=f*g-s;
       U[i*n+l]=f-g;

       for (k=l; k < n; k++) rv1[k]=U[i*n+k]/h;
       for (j=l; j < m; j++)
        {
         for (s=0.0, k=l; k < n; k++) s+=U[j*n+k]*U[i*n+k];
         for (k=l; k < n; k++) U[j*n+k]+=s*rv1[k];
        }

       for (k=l; k < n; k++) U[i*n+k]*=scale;
      }
    }

   if (anorm < (tmp=fabs(d[i])+fabs(rv1[i]))) anorm=tmp;
  }

 for (i=n-1; i >= 0; i--)
  {
   if (i < n-1)
    {
     if (g)
      {
       for (j=l; j < n; j++) V[j*n+i]=(U[i*n+j]/U[i*n+l])/g;
       for (j=l; j < n; j++)
        {
         for (s=0.0, k=l; k < n; k++) s+=U[i*n+k]*V[k*n+j];
         for (k=l; k < n; k++) V[k*n+j]+=s*V[k*n+i];
        }
      }

     for (j=l; j < n; j++) V[i*n+j]=V[j*n+i]=0.0;
    }

    V[i*n+i]=1.0;
   g=rv1[i];
   l=i;
  }

 for (i=(m < n) ? m-1 : n-1; i >= 0; i--)
  {
   l=i+1;
   g=d[i];

   for (j=l; j < n; j++) U[i*n+j]=0.0;

   if (g)
     {
      g=1.0/g;

      for (j=l; j < n; j++)
       {
        for (s=0.0, k=l; k < m; k++) s+=U[k*n+i]*U[k*n+j];
        f=(s/U[i*n+i])*g;
        for (k=i; k < m; k++) U[k*n+j]+=f*U[k*n+i];
       }

      for (j=i; j < m; j++) U[j*n+i]*=g;
     }
    else
     for (j=i; j < m; j++) U[j*n+i]=0.0;

   ++U[i*n+i];
  }

 for (k=n-1; k >= 0; k--)
  {
   for (its=1; its <= 30; its++)
    {
     flag=1;
     for (l=k; l >= 0; l--)
      {
       nm=l-1;

       if ((PRECISION)(fabs(rv1[l])+anorm) == anorm)
        {
         flag=0;
         break;
        }

       if ((PRECISION)(fabs(d[nm])+anorm) == anorm) break;
      }

     if (flag)
      {
       c=0.0;
       s=1.0;

       for (i=l; i <= k; i++)
        {
         f=s*rv1[i];
         rv1[i]=c*rv1[i];

         if ((PRECISION)(fabs(f)+anorm) == anorm) break;

         g=d[i];
         h=pythag(f,g);
         d[i]=h;
         h=1.0/h;
         c=g*h;
         s=-f*h;

         for (j=0; j < m; j++)
          {
           y=U[j*n+nm];
           z=U[j*n+i];
           U[j*n+nm]=y*c+z*s;
           U[j*n+i]=z*c-y*s;
          }
        }
      }

     z=d[k];

     if (l == k)
      {
       if (z < 0.0)
        {
         d[k]=-z;
         for (j=0; j < n; j++) V[j*n+k]=-V[j*n+k];
        }

       break;
      }

     if (its >= MAX_ITER)
      {
       free(rv1);
       return ITERATION_ERR;
      }

     x=d[l];
     nm=k-1;
     y=d[nm];
     g=rv1[nm];
     h=rv1[k];
     f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
     g=pythag(f,1.0);
     f=((x-z)*(x+z)+h*((y/(f+((f >= 0.0) ? fabs(g) : -fabs(g))))-h))/x;
     c=s=1.0;

     for (j=l; j <= nm; j++)
      {
       i=j+1;
       g=rv1[i];
       y=d[i];
       h=s*g;
       g=c*g;
       z=pythag(f,h);
       rv1[j]=z;
       c=f/z;
       s=h/z;
       f=x*c+g*s;
       g=g*c-x*s;
       h=y*s;
       y*=c;

       for (jj=0; jj < n; jj++)
        {
         x=V[jj*n+j];
         z=V[jj*n+i];
         V[jj*n+j]=x*c+z*s;
         V[jj*n+i]=z*c-x*s;
        }

       z=pythag(f,h);
       d[j]=z;

       if (z)
        {
         z=1.0/z;
         c=f*z;
         s=h*z;
        }

       f=c*g+s*y;
       x=c*y-s*g;

      for (jj=0; jj < m; jj++)
       {
        y=U[jj*n+j];
        z=U[jj*n+i];
        U[jj*n+j]=y*c+z*s;
        U[jj*n+i]=z*c-y*s;
       }
    }

   rv1[l]=0.0;
   rv1[k]=f;
   d[k]=x;
  }
 }

 free(rv1);

 return 0;
}
#undef MAX_ITER

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
     ITER_ERR   : maximum number of iterations (MAX_ITER) exceeded - only if
                  numerical recipe routines are used.

   Results
     Finds matrices U and V with orthonormal columns and a diagonal matrix 
     D=diag(d) with non-negative diagonal such that A = U*D*V'.  The matrix D is 
     m x n if compact = 0 and is q x q if compact = 1.  The flags ut, vt, and at
     determine the format of U, V, and A.  A value of 1 indicates column major 
     format and a value of 0 indicates row major format.  If either U or V is
     null, then it is not computed.

   Notes
     If A=U, U and A must be of the same size and ut=at.  If A=V, then V and A 
     must be of the same size and vt=at.  It cannot be the case that U=V unless
     both are null.
*/
int bSVD_new(PRECISION *U, PRECISION *d, PRECISION *V, PRECISION *A, int m, int n, int ut, int vt, int at, int compact)
{
#if defined USE_BLAS_LAPACK

#if (PRECISION_SIZE == 4)
#define gesdd sgesdd   
#define gesvd sgesvd
#else
#define gesdd dgesdd  
#define gesvd dgesvd
#endif

  int jobu, jobv, jobt, k=-1, info, err, m_, n_, qu_, qv_, transpose;
  PRECISION  *A_, *U_, *V_, *work, opt_size;

  A_=(PRECISION*)malloc(m*n*sizeof(PRECISION));

  jobu=jobv=compact ? 'S' : 'A';

  if (!U)
    {
      jobu='N';
      if (!V)
	{
	  jobv='N';
	  vt=transpose=1-at;
	}
      else
	transpose=vt;
      ut=1-vt;
    }
  else
    if (!V)
      {
	jobv='N';
	vt=transpose=1-ut;
      }
    else
      {
	if (ut != vt)
	  transpose=vt;
	else
	  transpose=1-at;
      }

  if (transpose)
    {
      jobt=jobu;
      jobu=jobv;
      jobv=jobt;
      if (at)
	bTranspose(A_,A,m,n,at);
      else
	memcpy(A_,A,m*n*sizeof(PRECISION));
      if (compact)
	{
	  m_=n;
	  n_=m;
	  qu_=qv_=(m < n) ? m : n;
	}
      else
	{
	  qu_=m_=n;
	  qv_=n_=m;
	}
      U_=vt ? V : (PRECISION*)malloc(m_*qu_*sizeof(PRECISION));
      V_=ut ? (PRECISION*)malloc(qv_*n_*sizeof(PRECISION)) : U;	      
    }
  else
    {
      if (at)
	memcpy(A_,A,m*n*sizeof(PRECISION));
      else
	bTranspose(A_,A,m,n,at);
      if (compact)
	{
	  m_=m;
	  n_=n;
	  qu_=qv_=(m < n) ? m : n;
	}
      else
	{
	  qu_=m_=m;
	  qv_=n_=n;
	}
      U_=ut ? U : (PRECISION*)malloc(m_*qu_*sizeof(PRECISION));
      V_=vt ? (PRECISION*)malloc(qv_*n_*sizeof(PRECISION)) : V;
    }

  // compute singular value decomposition
  gesvd(&jobu,&jobv,&m_,&n_,A_,&m_,d,U_,&m_,V_,&qv_,&opt_size,&k,&info);
  if (info || !(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION))))
    err=info ? BLAS_LAPACK_ERR : MEM_ERR;
  else
    {
      gesvd(&jobu,&jobv,&m_,&n_,A_,&m_,d,U_,&m_,V_,&qv_,work,&k,&info);
      free(work);
      if (info)
	err=BLAS_LAPACK_ERR;
      else
	{
	  if (transpose)
	    {
	      if (U != V_)
		bTranspose(U,V_,qv_,n_,1);
	      else
		if (V != U_)
		  bTranspose(V,U_,m_,qu_,1);
	    }
	  else
	    {
	      if (U != U_)
		bTranspose(U,U_,m_,qu_,1);
	      else
		if (V != V_)
		  bTranspose(V,V_,qv_,n_,1);
	    }
	  err=NO_ERR;
	}
    }

  free(A_);

  if (transpose)
    {
      if (U != V_)
	free(V_);
      else
	if (V != U_)
	  free(U_);
    }
  else
    {
      if (U != U_)
	free(U_);
      else
	if (V != V_)
	  free(V_);
    }

  return err;

#else
  
  PRECISION *NU;
  int i, j, rtrn;

  if (m == n)
    {
      if (at)
	if (U != A)
	  bTranspose(U,A,m,m,at);
	else
	  bTransposeInPlace(U,m);
      else
	if (U != A)
	  memcpy(U,A,m*m*sizeof(PRECISION));
	
      bSVD_NumericalRecipes(U,d,V,m,m);
      rtrn=NO_ERR;
    }
  else
    if (m < n)
      if (!(NU=(PRECISION*)malloc(m*n*sizeof(PRECISION))))
	rtrn=MEM_ERR;
      else
	{
	  if (at)
	    memcpy(NU,A,m*n*sizeof(PRECISION));
	  else
	    bTranspose(NU,A,n,m,1);
	  bSVD_NumericalRecipes(NU,d,U,n,m);
	  bQR_NumericalRecipes(V,NU,n,m);
	  for (i=m-1; i >= 0; i--)
	    if (NU[i*m+i] < 0)
	      for (j=(n-1)*n+i; j >= 0; j-=n) V[j]=-V[j];
	  rtrn=NO_ERR;
/* 	  if (!(nd=(PRECISION*)malloc(n*sizeof(PRECISION)))) */
/* 	    rtrn=MEM_ERR; */
/* 	  else */
/* 	    { */
/* 	      if (at) */
/* 		bTranspose(NU,A,m,n,at); */
/* 	      else */
/* 		memcpy(NU,A,m*n*sizeof(PRECISION)); */
/* 	      bSVD_NumericalRecipes(NU,nd,V,m,n); */
/* 	      for (i=m-1; i >= 0; i--) memcpy(U+i*m,NU+i*n,m*sizeof(PRECISION)); */
/* 	      memcpy(d,nd,m*sizeof(PRECISION)); */
/* 	      rtrn=NO_ERR; */
/* 	      free(nd); */
/* 	    } */
	  free(NU);
	}
    else
      if (!(NU=(PRECISION*)malloc(m*n*sizeof(PRECISION))))
	rtrn=MEM_ERR;
      else
	{
	  if (at)
	    bTranspose(NU,A,m,n,at);
	  else
	    memcpy(NU,A,m*n*sizeof(PRECISION));
	  bSVD_NumericalRecipes(NU,d,V,m,n);
	  bQR_NumericalRecipes(U,NU,m,n);
	  for (i=n-1; i >= 0; i--)
	    if (NU[i*n+i] < 0)
	      for (j=(m-1)*m+i; j >= 0; j-=m) U[j]=-U[j];
	  rtrn=NO_ERR;
	  free(NU);
	}
  if (vt) bTransposeInPlace(V,n);
  if (ut) bTransposeInPlace(U,m);
  return rtrn;

#endif
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
*/
int bSVD(PRECISION *U, PRECISION *d, PRECISION *V, PRECISION *A, int m, int n, int ut, int vt, int at)
{
#if defined USE_BLAS_LAPACK
#if (PRECISION_SIZE == 4)
#define gesdd sgesdd   
#define gesvd sgesvd
#else
#define gesdd dgesdd  
#define gesvd dgesvd
#endif
  int jobz='A', k, *iwork, info;
  PRECISION *X, *work, opt_size;
  if (!(X=(PRECISION*)malloc(m*n*sizeof(PRECISION)))) return MEM_ERR;
  memcpy(X,A,m*n*sizeof(PRECISION));
  if (!(iwork=(int*)malloc(8*((m < n) ? m : n)*sizeof(int))))
    {
      free(X);
      return MEM_ERR;
    }
  k=-1;  
  if (at)
    {
/*       gesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,&opt_size,&k,iwork,&info); */
/*       if (info || !(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/* 	{ */
/* 	  free(iwork); */
/* 	  free(X); */
/* 	  return info ? BLAS_LAPACK_ERR : MEM_ERR; */
/* 	} */
/*       gesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,work,&k,iwork,&info); */
/*       if (info) */
/* 	{ */
/* 	  free(work); */
	  memcpy(X,A,m*n*sizeof(PRECISION));
	  k=-1;
	  gesvd(&jobz,&jobz,&m,&n,X,&m,d,U,&m,V,&n,&opt_size,&k,&info);
	  if (info || !(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION))))
	    {
	      free(iwork);
	      free(X);
	      return info ? BLAS_LAPACK_ERR : MEM_ERR;
	    }
	  gesvd(&jobz,&jobz,&m,&n,X,&m,d,U,&m,V,&n,work,&k,&info);
	  if (info)
	    {
	      free(iwork);
	      free(X);
	      return BLAS_LAPACK_ERR;
	    }
/* 	} */
      if (!ut) 
	bTransposeInPlace(U,m);
      if (vt) 
	bTransposeInPlace(V,n);
    }
  else
    {
/*       gesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,&opt_size,&k,iwork,&info); */
/*       if (!(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/* 	{ */
/* 	  free(iwork); */
/* 	  free(X); */
/* 	  return MEM_ERR; */
/* 	} */
/*       gesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,work,&k,iwork,&info); */
/*       if (info) */
/* 	{ */
/* 	  free(work); */
	  memcpy(X,A,m*n*sizeof(PRECISION));
	  k=-1;
	  gesvd(&jobz,&jobz,&n,&m,X,&n,d,V,&n,U,&m,&opt_size,&k,&info);
	  if (info || !(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION))))
	    {
	      free(iwork);
	      free(X);
	      return info ? BLAS_LAPACK_ERR : MEM_ERR;
	    }
	  gesvd(&jobz,&jobz,&n,&m,X,&n,d,V,&n,U,&m,work,&k,&info);
	  if (info)
	    {
	      free(iwork);
	      free(X);
	      return BLAS_LAPACK_ERR;
	    }
/* 	} */
      if (!vt) 
	bTransposeInPlace(V,n);
      if (ut) 
	bTransposeInPlace(U,m);
    }
  free(work);
  free(iwork);
  free(X);
  return NO_ERR;
#undef gesdd   
#undef gesvd 
#else
  // bSVD
  PRECISION *NU;
  int i, j, rtrn;

  if (m == n)
    {
      if (at)
	if (U != A)
	  bTranspose(U,A,m,m,at);
	else
	  bTransposeInPlace(U,m);
      else
	if (U != A)
	  memcpy(U,A,m*m*sizeof(PRECISION));
	
      bSVD_NumericalRecipes(U,d,V,m,m);
      rtrn=NO_ERR;
    }
  else
    if (m < n)
      if (!(NU=(PRECISION*)malloc(m*n*sizeof(PRECISION))))
	rtrn=MEM_ERR;
      else
	{
	  if (at)
	    memcpy(NU,A,m*n*sizeof(PRECISION));
	  else
	    bTranspose(NU,A,n,m,1);
	  bSVD_NumericalRecipes(NU,d,U,n,m);
	  bQR_NumericalRecipes(V,NU,n,m);
	  for (i=m-1; i >= 0; i--)
	    if (NU[i*m+i] < 0)
	      for (j=(n-1)*n+i; j >= 0; j-=n) V[j]=-V[j];
	  rtrn=NO_ERR;
/* 	  if (!(nd=(PRECISION*)malloc(n*sizeof(PRECISION)))) */
/* 	    rtrn=MEM_ERR; */
/* 	  else */
/* 	    { */
/* 	      if (at) */
/* 		bTranspose(NU,A,m,n,at); */
/* 	      else */
/* 		memcpy(NU,A,m*n*sizeof(PRECISION)); */
/* 	      bSVD_NumericalRecipes(NU,nd,V,m,n); */
/* 	      for (i=m-1; i >= 0; i--) memcpy(U+i*m,NU+i*n,m*sizeof(PRECISION)); */
/* 	      memcpy(d,nd,m*sizeof(PRECISION)); */
/* 	      rtrn=NO_ERR; */
/* 	      free(nd); */
/* 	    } */
	  free(NU);
	}
    else
      if (!(NU=(PRECISION*)malloc(m*n*sizeof(PRECISION))))
	rtrn=MEM_ERR;
      else
	{
	  if (at)
	    bTranspose(NU,A,m,n,at);
	  else
	    memcpy(NU,A,m*n*sizeof(PRECISION));
	  bSVD_NumericalRecipes(NU,d,V,m,n);
	  bQR_NumericalRecipes(U,NU,m,n);
	  for (i=n-1; i >= 0; i--)
	    if (NU[i*n+i] < 0)
	      for (j=(m-1)*m+i; j >= 0; j-=m) U[j]=-U[j];
	  rtrn=NO_ERR;
	  free(NU);
	}
  if (vt) bTransposeInPlace(V,n);
  if (ut) bTransposeInPlace(U,m);
  return rtrn;
#endif
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/***************************** QR Decompositions ******************************/
/******************************************************************************/
/*
   Assumes
     Q : m x m matrix in row major format or null pointer
     R : m x n matrix in row major format

   Returns
     0 (NO_ERR) : Success
     MEM_ERR    : Out of memory

   Results
     Finds an orthogonal matrix Q and an upper triangular matrix U
     such that

                                R = Q * U
     
     The matrix U is returned in R and Q is computed only if it is
     not null.

   Notes
     The QR decomposition is formed using Householder matrices without
     pivoting.  If Q is null, then the matrix Q will have the property 
     that 

                             det(Q) = (-1)^s

     where s=min{m-1,n).
*/
static int bQR_NumericalRecipes(PRECISION *Q, PRECISION *R, int m, int n)
{
 int i, j, k, s=(m <= n) ? m-1 : n;
 PRECISION tmp, scale, sigma, c, sum;
 PRECISION *diag, *norm;

 if (Q)
   {
    if (!(diag=(PRECISION*)malloc(2*s*sizeof(PRECISION)))) return MEM_ERR;
    norm=diag+s;

    for (k=0; k < s; k++)
     {
      for (scale=0.0, i=k; i < m; i++)
       if ((tmp=fabs(R[i*n+k])) > scale) scale=tmp;

      if (scale == 0.0)
       {
        diag[k]=norm[k]=0.0;
        continue;
       }
   
      for (scale=1.0/scale, sigma=0.0, i=k; i < m; i++)
       {
        R[i*n+k]*=scale;
        sigma+=R[i*n+k]*R[i*n+k];
       }

      sigma=sqrt(sigma);

      if (R[k*n+k] < 0.0) sigma=-sigma;

      diag[k]=R[k*n+k]+=sigma;

      norm[k]=c=1.0/(sigma*diag[k]);

      for (j=k+1; j < n; j++)
       {
        for (sum=0.0, i=k; i < m; i++) sum+=R[i*n+k]*R[i*n+j];
        sum*=c;
        for (i=k; i < m; i++) R[i*n+j]-=R[i*n+k]*sum;
       }

      R[k*n+k]=-sigma/scale;
     }

    for (i=m*m-1; i >= 0; i--) Q[i]=0.0;
    for (i=m*m-1; i >= 0; i-=m+1) Q[i]=1.0;

    for (k=s-1; k >= 0; k--)
     {
      c=norm[k];
      for (j=k; j < m; j++)
       {
        for (sum=diag[k]*Q[k*m+j], i=k+1; i < m; i++) sum+=R[i*n+k]*Q[i*m+j];
        for (Q[k*m+j]-=diag[k]*(sum*=c), i=k+1; i < m; i++) Q[i*m+j]-=R[i*n+k]*sum;
       }
     }

    for (k=s-1; k >= 0; k--)
     for (i=k+1; i < m; i++)
      R[i*n+k]=0;

    free(diag);

    return 0;
   }
  else
   {
    for (k=0; k < s; k++)
     {
      /* Find largest element below diagonal in kth column */
      for (scale=0.0, i=m-1; i > k; i--)
       if ((tmp=fabs(R[i*n+k])) > scale) scale=tmp;

      if (scale == 0.0)
        for (j=k; j < n; j++) R[k*n+j]=-R[k*n+j];
       else
        {
         if ((tmp=fabs(R[k*n+k])) > scale) scale=tmp;

         /* Compute normalized sum of squares */
         for (scale=1.0/scale, sigma=0.0, i=k; i < m; i++)
          {
           R[i*n+k]*=scale;
           sigma+=R[i*n+k]*R[i*n+k];
          }

         sigma=sqrt(sigma);

         /* choose sign of sigma to match the diagonal elements sign */
         if (R[k*n+k] < 0.0) sigma=-sigma;

         /* c = sigma*(diagonal + sigma) */
         R[k*n+k]+=sigma;
         c=1.0/(sigma*R[k*n+k]);

         for (j=k+1; j < n; j++)
          {
           for (sum=0.0, i=k; i < m; i++) sum+=R[i*n+k]*R[i*n+j];
           sum*=c;
           for (i=k; i < m; i++) R[i*n+j]-=R[i*n+k]*sum;
          }

         R[k*n+k]=-sigma/scale;
         for (i=k+1; i < m; i++) R[i*n+k]=0;
        }
     }

    return 0;
   }
}


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

     -- Not implemented --
     If Q is not null, then the matrix Q will have the 
     property that 

                             det(Q) = (-1)^s

     where s=min{m-1,n).
*/
int bQR(PRECISION *Q, PRECISION *R, PRECISION *X, int m, int n, int q, int qt, int rt, int xt)
{
#if defined USE_BLAS_LAPACK
  int i, j, k, l, lwork, info, p=(m < n) ? m : n;
  PRECISION *tau, *work, *ptr, opt_size;
  if (!(tau=(PRECISION*)malloc(p*sizeof(PRECISION)))) return MEM_ERR;
  if (xt)
    {
      lwork=-1;
#if (PRECISION_SIZE == 4)
      sgeqrf(&m,&n,X,&m,tau,&opt_size,&lwork,&info);
#else
      dgeqrf(&m,&n,X,&m,tau,&opt_size,&lwork,&info);
#endif
      if (!(work=(PRECISION*)malloc((lwork=(int)opt_size)*sizeof(PRECISION))))
	{
	  free(tau);
	  return MEM_ERR;
	}
#if (PRECISION_SIZE == 4)
      sgeqrf(&m,&n,X,&m,tau,work,&lwork,&info);
#else
      dgeqrf(&m,&n,X,&m,tau,work,&lwork,&info);
#endif
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
	    if (!(ptr=(PRECISION*)malloc(m*q*sizeof(PRECISION))))
	      {
		free(tau);
		return MEM_ERR;
	      }
	  memcpy(ptr,X,m*p*sizeof(PRECISION));
	  lwork=-1;
#if (PRECISION_SIZE == 4)
	  sorgqr(&m,&q,&p,ptr,&m,tau,&opt_size,&lwork,&info);
#else
	  dorgqr(&m,&q,&p,ptr,&m,tau,&opt_size,&lwork,&info);
#endif
	  if (!(work=(PRECISION*)malloc((lwork=(int)opt_size)*sizeof(PRECISION))))
	    {
	      if (!qt) free(ptr);
	      free(tau);
	      return MEM_ERR;
	    }
#if (PRECISION_SIZE == 4)
	  sorgqr(&m,&q,&p,ptr,&m,tau,work,&lwork,&info);
#else
	  dorgqr(&m,&q,&p,ptr,&m,tau,work,&lwork,&info);
#endif
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
#if (PRECISION_SIZE == 4)
      sgelqf(&n,&m,X,&n,tau,&opt_size,&lwork,&info);
#else
      dgelqf(&n,&m,X,&n,tau,&opt_size,&lwork,&info);
#endif
      if (!(work=(PRECISION*)malloc((lwork=(int)opt_size)*sizeof(PRECISION))))
	{
	  free(tau);
	  return MEM_ERR;
	}
#if (PRECISION_SIZE == 4)
      sgelqf(&n,&m,X,&n,tau,work,&lwork,&info);
#else
      dgelqf(&n,&m,X,&n,tau,work,&lwork,&info);
#endif
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
	    if (!(ptr=(PRECISION*)malloc(m*q*sizeof(PRECISION))))
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
#if (PRECISION_SIZE == 4)
	  sorglq(&q,&m,&p,ptr,&q,tau,&opt_size,&lwork,&info);
#else
	  dorglq(&q,&m,&p,ptr,&q,tau,&opt_size,&lwork,&info);
#endif
	  if (!(work=(PRECISION*)malloc((lwork=(int)opt_size)*sizeof(PRECISION))))
	    {
	      if (!qt) free(ptr);
	      free(tau);
	      return MEM_ERR;
	    }
#if (PRECISION_SIZE == 4)
	  sorglq(&q,&m,&p,ptr,&q,tau,work,&lwork,&info);
#else
	  dorglq(&q,&m,&p,ptr,&q,tau,work,&lwork,&info);
#endif
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
#else
  // bQR
  PRECISION *NQ, *NR, *pQ;
  int i, j;

  if (Q && (q != m))
    {
      if (!(NQ=(PRECISION*)malloc(m*m*sizeof(PRECISION))))
	return MEM_ERR;
    }
  else
    NQ=Q;

  if (rt || (q != m))
    {
      if (!(NR=(PRECISION*)malloc(m*n*sizeof(PRECISION))))
	{
	  if (NQ != Q) free(NQ);
	  return MEM_ERR;
	}
      if (xt)
	bTranspose(NR,X,m,n,xt);
      else
	memcpy(NR,X,m*n*sizeof(PRECISION));
    }
  else
    {
      NR=R;
      if (xt)
	bTranspose(NR,X,m,n,xt);
      else
	if (X != NR)
	  memcpy(NR,X,m*n*sizeof(PRECISION));
    }

  bQR_NumericalRecipes(NQ,NR,m,n);
  
  if (Q)
    if (q != m)
      {
	if (qt)
	  for (pQ=Q+m*n-1, j=n-1; j >= 0; j--)
	    for (i=m*(m-1)+j; i >= 0; pQ--, i-=m)
	      *pQ=NQ[i];	  
	else
	  for (i=m-1; i >= 0; i--) memcpy(Q+i*n,NQ+i*m,n*sizeof(PRECISION));
	free(NQ);
      }
    else
      if (qt)
	bTransposeInPlace(Q,m);

  if (rt)
    {
      bTranspose(R,NR,q,n,0);
      free(NR);
    }
  else
    if (q != m)
      {
	memcpy(R,NR,n*n*sizeof(PRECISION));
	free(NR);
      }

  return NO_ERR;
#endif
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
#if defined USE_INLINE
 int i, b, rtrn=NO_ERR;
 __asm{                                                       //    eax        ebx      ecx        edx        esi        edi          st(0)          st(1)          st(2) 
        fninit
        mov    edx,m                                          //     ?          ?        ?          m          ?          ?                        

        cmp    u,0
        //je     short c0                                       //     ?          ?        ?          m          ?          ?
        je     c0

        cmp    t,0
        //je     short b0                                       //     ?          ?        ?          m          ?          ?
        je     b0
 /********************************************************************************************************************************/
                                                              //     ?         ?         ?          m          ?          ?
        mov    eax,edx
        dec    eax                                            //     i         ?         ?          m          ?          ?     

        mov    edi,edx
        imul   edi,eax
        shl    edi,PRECISION_SHIFT
        add    edi,X                                          //     i         ?         ?          m          ?        X+i*m       

a1:                                                           //     i         ?         ?          m          ?        X+i*m 
        mov    i,eax

        mov    ebx,edx
        dec    ebx                                            //     i         j         ?          m          ?        X+i*m

        mov    esi,edx
        imul   esi,ebx
        add    esi,eax
        shl    esi,PRECISION_SHIFT
        add    esi,X                                          //     i         j         ?          m       X+i+j*m     X+i*m
       
a2:                                                           //     i         j         ?          m       X+i+j*m     X+i*m
        cmp    ebx,eax
        jle    short a3

        mov    dword ptr[esi],0

#if (PRECISION_SIZE == 8)
        mov    dword ptr[esi+4],0
#endif

        mov    ecx,m                                          //     i         j         m          m       X+i+j*m     X+i*m
        shl    ecx,PRECISION_SHIFT                            //     i         j  m*PRECISION_SIZE  m       X+i+j*m     X+i*m
        sub    esi,ecx                                        //     i         j  m*PRECISION_SIZE  m       X+i+j*m     X+i*m
        dec    ebx                                            //     i         j         ?          m       X+i+j*m     X+i*m
        jmp    short a2                                       //     i         j         ?          m       X+i+j*m     X+i*m

a3:
        fld    PRECISION_WORD ptr[esi]                        //     i         j         ?          m       X+i+j*m     X+i*m       X[i+j*m]

a4:                                                           //     k         j         ?          m       X+i+j*m     X+i*m       X[i+j*m]
        inc    eax                                            //     k         j         ?          m       X+i+j*m     X+i*m       X[i+j*m]
        cmp    edx,eax
        jle    short a5

        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]     //     k         j         ?          m       X+i+j*m     X+i*m        X[i*m+k]      X[i+j*m]
        fld    st(0)                                          //     k         j         ?          m       X+i+j*m     X+i*m        X[i*m+k]      X[i*m+k]       X[i+j*m]
        fmul                                                  //     k         j         ?          m       X+i+j*m     X+i*m       X[i*m+k]^2     X[i+j*m]
        fsub                                                  //     k         j         ?          m       X+i+j*m     X+i*m        X[i+j*m]

        jmp    short a4                                       //     k         j         ?          m       X+i+j*m     X+i*m        X[i+j*m]

a5:                                                           //     ?         j         ?          m       X+i+j*m     X+i*m        X[i+j*m]
        ftst
        fnstsw  ax                                               
        sahf                                                      
        //jbe     short e1                                      //     ?         j         ?          m       X+i+j*m     X+i*m        X[i+j*m]
        jbe     e1

        fsqrt                                                 //     ?         j         ?          m       X+i+j*m     X+i*m        X[i+j*m] 
        fst     PRECISION_WORD ptr[esi]                       //     ?         j         ?          m       X+i+j*m     X+i*m        X[i+j*m] 
        
        fld1                                                  //     ?         j         ?          m       X+i+j*m     X+i*m            1           
        fdivr                                                 //     ?         j         ?          m       X+i+j*m     X+i*m         scale 

        mov    ecx,edi                                        //     ?         j       X+j*m        m       X+i+j*m     X+i*m         scale

a6:                                                           //     ?         j       X+j*m        m       X+i+j*m     X+i*m         scale
        dec    ebx                                            //     ?         j       X+j*m        m       X+i+j*m     X+i*m         scale
        jl     short a9                                       //     ?         j       X+j*m        m       X+i+j*m     X+i*m         scale

        mov    eax,m                                          //     m         j       X+j*m        m       X+i+j*m     X+i*m         scale
        shl    eax,PRECISION_SHIFT                        // m*PRECISION_SIZE  j       X+j*m        m       X+i+j*m     X+i*m         scale
        sub    ecx,eax                                    // m*PRECISION_SIZE  j       X+j*m        m       X+i+j*m     X+i*m         scale
        sub    esi,eax                                        //     ?         j       X+j*m        m       X+i+j*m     X+i*m         scale

        fld    PRECISION_WORD ptr[esi]                        //     ?         j       X+j*m        m       X+i+j*m     X+i*m        X[i+j*m]        scale

        mov    eax,i                                          //     k         j       X+j*m        m       X+i+j*m     X+i*m        X[i+j*m]        scale

a7:                                                           //     k         j       X+j*m        m       X+i+j*m     X+i*m        X[i+j*m]        scale
        inc    eax                                            
        cmp    edx,eax
        jle    short a8

        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]     //     k         j       X+j*m        m       X+i+j*m     X+i*m        X[i*m+k]       X[i+j*m]        scale
        fmul   PRECISION_WORD ptr[ecx+PRECISION_SIZE*eax]     //     k         j       X+j*m        m       X+i+j*m     X+i*m   X[i*m+k]*X[j*m+k]   X[i+j*m]        scale
        fsub                                                  //     k         j       X+j*m        m       X+i+j*m     X+i*m        X[i+j*m]        scale

        jmp    short a7                                       //     k         j       X+j*m        m       X+i+j*m     X+i*m        X[i+j*m]        scale

a8:
        fmul   st,st(1)                                       //     k         j       X+j*m        m       X+i+j*m     X+i*m        X[i+j*m]        scale
        fstp   PRECISION_WORD ptr[esi]                        //     k         j       X+j*m        m       X+i+j*m     X+i*m         scale
        jmp    short a6                                       //     k         j       X+j*m        m       X+i+j*m     X+i*m         scale

a9:     ffree  st(0)                                          //     ?         ?         ?          m       X+i+j*m     X+i*m

        mov    eax,edx
        shl    eax,PRECISION_SHIFT
        sub    edi,eax                                        //     ?         ?         ?          m       X+i+j*m     X+i*m

        mov    eax,i                                          //     i         ?         ?          m       X+i+j*m     X+i*m
        dec    eax                                            //     i         ?         ?          m       X+i+j*m     X+i*m
        //jge    short a1                                       //     i         ?         ?          m       X+i+j*m     X+i*m
        jge    a1

        //jmp    short e2                                       //     ?         ?         ?          ?          ?          ?
        jmp    e2
/********************************************************************************************************************************/ 
b0:                                                           //     ?          ?        ?          m          ?          ?                        
        mov    eax,0

        mov    edi,X                                          //     i                              m                    X+i

b1:                                                           //     i                              m                    X+i
        mov    i,eax

        mov    ebx,0                                          //     i         j                    m       X+i*m+j      X+i

        mov    esi,eax
        imul   esi,edx
        shl    esi,PRECISION_SHIFT
        add    esi,X                                          //     i         j                    m       X+i*m+j      X+i

b2:       
        cmp    eax,ebx
        jle    short b3

        mov    dword ptr[esi],0

#if (PRECISION_SIZE == 8)
        mov    dword ptr[esi+4],0
#endif

        add    esi,PRECISION_SIZE
        inc    ebx
        jmp    short b2

b3:
        fld    PRECISION_WORD ptr[esi]
        dec    eax                                            
        imul   eax,edx                                        //     k         j                    m       X+i+j*m     X+i

        cmp    eax,0
        jl     short b5

b4:
        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]
        fld    st(0)
        fmul
        fsub

        sub    eax,edx
        jge    short b4

b5:
        ftst
        fnstsw  ax                                               
        sahf
        //jbe     short e1
        jbe     e1

        fsqrt
        fst     PRECISION_WORD ptr[esi]
        
        fld1
        fdivr  

        mov    ecx,edi                                        //               j        X+j         m       X+i*m+j     X+i

b6:                                                           //               j        X+j         m       X+i*m+j     X+i
        inc    ebx
        cmp    edx,ebx
        jle    short b9

        add    ecx,PRECISION_SIZE
        add    esi,PRECISION_SIZE

        fld    PRECISION_WORD ptr[esi]

        mov    eax,i                                          
        dec    eax
        imul   eax,edx                                        //     k         j       X+j*m        m       X+i+j*m     X+i*m

        cmp    eax,0
        jl     short b8

b7:
        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]
        fmul   PRECISION_WORD ptr[ecx+PRECISION_SIZE*eax]
        fsub

        sub    eax,edx                                            
        jge    short b7

b8:
        fmul   st,st(1)
        fstp   PRECISION_WORD ptr[esi]
        jmp    short b6


b9:     ffree   st(0) 

        add     edi,PRECISION_SIZE

        mov     eax,i
        inc     eax
        cmp     edx,eax
        //jg      short b1
        jg      b1

        //jmp     short e2
        jmp     e2
/********************************************************************************************************************************/  
c0:
        mov    eax,t
        cmp    eax,0
        //je     short d0
        je     d0
/********************************************************************************************************************************/ 
        mov    eax,0                                          //     i                              m                   X+i*m

        mov    edi,X                                          //                                    m                   X+i*m

c1:                                                           //     i                              m                   X+i*m
        mov    i,eax

        mov    ebx,0                                          //     i         j                    m       X+i+j*m     X+i*m
       
        mov    esi,eax
        shl    esi,PRECISION_SHIFT
        add    esi,X                                          //     i                              m       X+i+j*m     X+i*m

c2:       
        cmp    eax,ebx
        jle    short c3

        mov    dword ptr[esi],0

#if (PRECISION_SIZE == 8)
        mov    dword ptr[esi+4],0
#endif

        mov    ecx,edx
        shl    ecx,PRECISION_SHIFT
        add    esi,ecx
        inc    ebx
        jmp    short c2

c3:
        fld    PRECISION_WORD ptr[esi]
c4:                                                           //     k         j                    m       X+i+j*m    X+i*m
        dec    eax                                            
        jl     short c5

        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]
        fld    st(0)
        fmul
        fsub

        jmp    short c4

c5:
        ftst
        fnstsw  ax                                               
        sahf                                                      
        //jbe     short e1
        jbe     e1

        fsqrt
        fst     PRECISION_WORD ptr[esi]
        
        fld1
        fdivr  

        mov    ecx,edi                                        //               j        X+j         m       X+i+j*m    X+i*m

c6:                                                           //               j        X+j         m       X+i+j*m    X+i*m
        inc    ebx
        cmp    edx,ebx
        jle    short c9

        mov    eax,edx
        shl    eax,PRECISION_SHIFT
        add    ecx,eax
        add    esi,eax

        fld    PRECISION_WORD ptr[esi]

        mov    eax,i
c7:                                                          //     k         j         X+j        m       X+i+j*m    X+i*m                      
        dec    eax                                           
        jl     short c8

        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]
        fmul   PRECISION_WORD ptr[ecx+PRECISION_SIZE*eax]
        fsub

        jmp    short c7

c8:
        fmul   st,st(1)
        fstp   PRECISION_WORD ptr[esi]
        jmp    short c6

c9:     ffree   st(0) 

        mov    eax,edx
        shl    eax,PRECISION_SHIFT
        add    edi,eax

        mov    eax,i
        inc    eax
        cmp    edx,eax
        //jg     short c1
        jg     c1

        //jmp    short e2
        jmp    e2
/********************************************************************************************************************************/ 
d0:                                                           //     ?          ?        ?          m          ?          ?
        mov    eax,edx
        imul   eax,edx
        mov    b,eax

        mov    eax,edx
        dec    eax                                            //     i         ?         ?          m       X+i*m+j      X+i         

        mov    edi,eax
        shl    edi,PRECISION_SHIFT
        add    edi,X                                          //     ?         ?         ?          m          ?         X+i

        mov    esi,edx
        imul   esi,edx
        shl    esi,PRECISION_SHIFT
        add    esi,X                                          //     ?         j         ?          m       X+i*m+j      X+i   
                
d1:                                                           //     i         ?         ?          m       X+i*m+j      X+i 
        mov    i,eax
        
        mov    ebx,edx                                        //     i         j         ?          m       X+i*m+j      X+i
        dec    ebx

        sub    esi,PRECISION_SIZE
               
d2:                                                           //     i         j         ?          m       X+i*m+j      X+i
        cmp    ebx,eax
        jle    short d3

        mov    dword ptr[esi],0

#if (PRECISION_SIZE == 8)
        mov    dword ptr[esi+4],0
#endif

        sub    esi,PRECISION_SIZE                             //     i         j         ?          m       X+i*m+j      X+i
        dec    ebx
        jmp    short d2                                       //     i         j         ?          m       X+i*m+j      X+i

d3:
        fld    PRECISION_WORD ptr[esi]                        //     i         j         ?          m       X+i*m+j      X+i        X[i+j*m]

        imul   eax,edx
   
d4:                                                           //     k         j         ?          m       X+i*m+j      X+i        X[i+j*m]
        add    eax,edx                                        //     k         j         ?          m       X+i*m+j      X+i        X[i+j*m]
        cmp    b,eax
        jle    short d5

        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]     //     k         j         ?          m       X+i*m+j      X+i         X[i*m+k]      X[i+j*m]
        fld    st(0)                                          //     k         j         ?          m       X+i*m+j      X+i         X[i*m+k]      X[i*m+k]       X[i+j*m]
        fmul                                                  //     k         j         ?          m       X+i*m+j      X+i        X[i*m+k]^2     X[i+j*m]
        fsub                                                  //     k         j         ?          m       X+i*m+j      X+i         X[i+j*m]

        jmp    short d4                                       //     k         j         ?          m       X+i*m+j      X+i         X[i+j*m]
 
d5:                                                           //     ?         j         ?          m       X+i*m+j      X+i         X[i+j*m]
        ftst                                                    
        fnstsw  ax                                               
        sahf                                                      
        jbe     short e1                                      //     ?         j         ?          m       X+i*m+j      X+i         X[i+j*m] 
  
        fsqrt                                                 //     ?         j         ?          m       X+i*m+j      X+i         X[i+j*m] 
        fst     PRECISION_WORD ptr[esi]                       //     ?         j         ?          m       X+i*m+j      X+i         X[i+j*m] 
        
        fld1                                                  //     ?         j         ?          m       X+i*m+j      X+i             1           
        fdivr                                                 //     ?         j         ?          m       X+i*m+j      X+i          scale 

        mov    ecx,edi                                        //     ?         j        X+j         m       X+i*m+j      X+i          scale

d6:                                                           //     ?         j        X+j         m       X+i*m+j      X+i          scale
        dec    ebx                                            //     ?         j        X+j         m       X+i*m+j      X+i          scale
        jl     short d9                                       //     ?         j        X+j         m       X+i*m+j      X+i          scale

        sub    ecx,PRECISION_SIZE                             //     ?         j        X+j         m       X+i*m+j      X+i          scale
        sub    esi,PRECISION_SIZE                             //     ?         j        X+j         m       X+i*m+j      X+i          scale

        fld    PRECISION_WORD ptr[esi]                        //     ?         j        X+j         m       X+i*m+j      X+i         X[i+j*m]        scale

        mov    eax,i                                          //     k         j        X+j         m       X+i*m+j      X+i         X[i+j*m]        scale
        imul   eax,edx

d7:                                                           //     k         j        X+j         m       X+i*m+j      X+i         X[i+j*m]        scale
        add    eax,edx                                            
        cmp    b,eax
        jle    short d8

        fld    PRECISION_WORD ptr[edi+PRECISION_SIZE*eax]     //     k         j        X+j         m       X+i*m+j      X+i         X[i*m+k]       X[i+j*m]        scale
        fmul   PRECISION_WORD ptr[ecx+PRECISION_SIZE*eax]     //     k         j        X+j         m       X+i*m+j      X+i    X[i*m+k]*X[j*m+k]   X[i+j*m]        scale
        fsub                                                  //     k         j        X+j         m       X+i*m+j      X+i         X[i+j*m]        scale

        jmp    short d7                                       //     k         j        X+j         m       X+i*m+j      X+i         X[i+j*m]        scale

d8:
        fmul   st,st(1)                                       //     k         j        X+j         m       X+i*m+j      X+i         X[i+j*m]        scale
        fstp   PRECISION_WORD ptr[esi]                        //     k         j        X+j         m       X+i*m+j      X+i          scale
        jmp    short d6                                       //     k         j        X+j         m       X+i*m+j      X+i          scale

d9:     ffree  st(0)                                          //     ?         ?         ?          m       X+i*m+j      X+i

        sub    edi,PRECISION_SIZE

        mov    eax,i                                          //     i         ?         ?          m       X+i*m+j      X+i
        dec    eax                                            //     i         ?         ?          m       X+i*m+j      X+i
        jge    short d1                                       //     i         ?         ?          m       X+i*m+j      X+i

        jmp    short e2 
/********************************************************************************************************************************/ 
e1:     mov     rtrn,POSDEF_ERR                 
e2:
       }
 return rtrn;
#else
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
#endif
 
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

/*******************************************************************************

  Accumulation:
   u is an m x p matrix
   v is a p x n matrix
   0 <= i < m
   0 <= j < n
   ui=p  uj=1  (u is in row major format)
   ui=1  uj=m  (u is in column major format)
   vi=n  vj=1  (v is in row major format)
   vi=1  vj=p  (v is in column major format)  
   
   Computes the following sum:

     u[i*ui + 0*uj]*v[0*vi + j*vj] + u[i*ui+1*uj]*v[1*vi+j*vj] + ...

                                   ... + u[i*ui + (p-1)*uj]*v[(p-1)*vi + j*vj]
      

  *** C code ***

  int k=(p-1)*vi;
  PRECISION *pu=u + i*ui + (p-1)*uj;
  PRECISION *pv=v + j*vj;

  PRECISION tmp=(*pu)*pv[k];
  while ((k-=vi) >= 0) tmp+=(*(pu-=uj))*pv[k];
  w(i,j)=tmp;

  *** assembly code *** 

  __asm {

         
         // Upon entry:  
         //   eax=(p-1)*kv                           
         //   ebx=ku*PRECISION_SIZE  
         //   ecx=kv  
         //   esi=u(i) + (p-1)*uj*PRECISION_SIZE
         //   edi=v(j)
         // 
         // Upon exit:   
         //   eax=0         
         //   ebx=ku*PRECISION_SIZE  
         //   ecx=kv  
         //   esi=u(i)
         //   edi=v(j)   
             
         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]      // tmp=(*pu)*v(j)[k]

         sub  eax,ecx                                          // k-=kv
         jl   short dest2                                      // k >= 0

dest1  : sub  esi,ebx                                          // pu-=ku

         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]      
         fadd                                                  // tmp+=(*pu)*v(j)[k]

         sub  eax,ecx                                          // k-=kv
         jge  short dest1_0                                    // k >= 0

dest2  : fstp PRECISION_WORD ptr [edx]                         // w(i,j)=tmp
        }

 Notes:
  The register edx is unused.  It is recommended that this register be used to store a pointer
  to the storage position for the accumlated value.  Furthermore, the direction of outer loops 
  should be such that edx-PRECISION_SIZE is the storage position for the next accumlated value.
   
  If vi is known to be one, then ecx is available.

  If uj is known to be one, then ebx is available and u+(i-1)*ui+(p-1)*uj = u+i*ui-1.  If u is
  a matrix in row major format, then uj will be one. 

*******************************************************************************/


/*******************************************************************************
  Back solving
   Given u(i), v(j), ku, kv, and p, computes

    v(j)[0]=(u(i)[(p-1)*ku]*v(j)[(p-1)*kv] + ... + u(i)[ku]*v(j)[kv])/u(i)[0]

   p is assumed to be positive and ku and kv are positive.
  
  

  *** C code ***

  int k=(p-1)*kv;
  PRECISION *pu=u(i)+(p-1)*ku;

  PRECISION tmp=(*pu)*v(j)[k];
  while (k != 0) tmp+=(*(pu-=ku))*v(j)[k-=kv];

  *** assembly code *** 

  __asm {

         
         // Upon entry:  
         //   eax=(p-1)*kv  
         //   ebx=ku*PRECISION_SIZE  
         //   ecx=kv  
         //   esi=u(i) + (p-1)*uj*PRECISION_SIZE
         //   edi=v(j)
         // 
         // Upon exit:   
         //   eax=0         
         //   ebx=ku*PRECISION_SIZE  
         //   ecx=kv  
         //   esi=u(i)
         //   edi=v(j)   
             
         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]

         test eax,eax
         je   short dest2

dest1:   sub  eax,ecx                                       // if kv=1, replace with:  dec  eax 
         sub  esi,ebx                                       // if ku=1, replace with:  sub  esi,PRECISION_SIZE

         fld  PRECISION_WORD ptr [esi]
         fmul PRECISION_WORD ptr [edi+PRECISION_SIZE*eax]
         fadd

         test eax,eax
         jne  dest1

dest2:   fstp PRECISION_WORD ptr [edx]                      // pop accumulated value off stack and store
       
         sub  edx,PRECISION_SIZE                          

        }

 Notes:
  The register edx is unused.  It is recommended that this register be used to store a pointer
  to the storage position for the accumlated value.  Furthermore, the direction of outer loops 
  should be such that edx-PRECISION_SIZE is the storage position for the next accumlated value.
   
  If vi is known to be one, then ecx is available.

  If uj is known to be one, then ebx is available and u+(i-1)*ui+(p-1)*uj = u+i*ui-1.  If u is
  a matrix in row major format, then uj will be one. 

*******************************************************************************/


/*   int jobz='A', k, *iwork, info; */
/*   PRECISION *X, *work, opt_size; */
/*   if (!(X=(PRECISION*)malloc(m*n*sizeof(PRECISION)))) return MEM_ERR; */
/*   memcpy(X,A,m*n*sizeof(PRECISION)); */
/*   if (!(iwork=(int*)malloc(8*((m < n) ? m : n)*sizeof(int)))) */
/*     { */
/*       free(X); */
/*       return MEM_ERR; */
/*     } */
/*   k=-1;   */
/*   if (at) */
/*     { */
/* #if (PRECISION_SIZE == 4) */
/*       sgesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,&opt_size,&k,iwork,&info); */
/* #else */
/*       dgesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,&opt_size,&k,iwork,&info); */
/* #endif */
/*       if (info || !(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/* 	{ */
/* 	  free(iwork); */
/* 	  free(X); */
/* 	  return info ? BLAS_LAPACK_ERR : MEM_ERR; */
/* 	} */
/* #if (PRECISION_SIZE == 4) */
/*       sgesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,work,&k,iwork,&info); */
/* #else */
/*       dgesdd(&jobz,&m,&n,X,&m,d,U,&m,V,&n,work,&k,iwork,&info); */
/* #endif */
/*       if (info) */
/* 	{ */
/* 	  free(work); */
/* 	  memcpy(X,A,m*n*sizeof(PRECISION)); */
/* 	  k=-1; */
/* #if (PRECISION_SIZE == 4) */
/* 	  sgesvd(&jobz,&jobz,&m,&n,X,&m,d,U,&m,V,&n,&opt_size,&k,&info); */
/* #else */
/* 	  dgesvd(&jobz,&jobz,&m,&n,X,&m,d,U,&m,V,&n,&opt_size,&k,&info); */
/* #endif */
/* 	  if (info || !(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/* 	    { */
/* 	      free(iwork); */
/* 	      free(X); */
/* 	      return info ? BLAS_LAPACK_ERR : MEM_ERR; */
/* 	    } */
/* #if (PRECISION_SIZE == 4) */
/* 	  sgesvd(&jobz,&jobz,&m,&n,X,&m,d,U,&m,V,&n,work,&k,&info); */
/* #else */
/* 	  dgesvd(&jobz,&jobz,&m,&n,X,&m,d,U,&m,V,&n,work,&k,&info); */
/* #endif */
/* 	  if (info) */
/* 	    { */
/* 	      free(iwork); */
/* 	      free(X); */
/* 	      return BLAS_LAPACK_ERR; */
/* 	    } */
/* 	} */
/*       if (!ut)  */
/* 	bTransposeInPlace(U,m); */
/*       if (vt)  */
/* 	bTransposeInPlace(V,n); */
/*     } */
/*   else */
/*     { */
/* #if (PRECISION_SIZE == 4) */
/*       sgesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,&opt_size,&k,iwork,&info); */
/* #else */
/*       dgesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,&opt_size,&k,iwork,&info); */
/* #endif */
/*       if (!(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/* 	{ */
/* 	  free(iwork); */
/* 	  free(X); */
/* 	  return MEM_ERR; */
/* 	} */
/* #if (PRECISION_SIZE == 4) */
/*       sgesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,work,&k,iwork,&info); */
/* #else */
/*       dgesdd(&jobz,&n,&m,X,&n,d,V,&n,U,&m,work,&k,iwork,&info); */
/* #endif */
/*       if (info) */
/* 	{ */
/* 	  free(work); */
/* 	  memcpy(X,A,m*n*sizeof(PRECISION)); */
/* 	  k=-1; */
/* #if (PRECISION_SIZE == 4) */
/* 	  sgesvd(&jobz,&jobz,&n,&m,X,&n,d,V,&n,U,&m,&opt_size,&k,&info); */
/* #else */
/* 	  dgesvd(&jobz,&jobz,&n,&m,X,&n,d,V,&n,U,&m,&opt_size,&k,&info); */
/* #endif */
/* 	  if (info || !(work=(PRECISION*)malloc((k=(int)opt_size)*sizeof(PRECISION)))) */
/* 	    { */
/* 	      free(iwork); */
/* 	      free(X); */
/* 	      return info ? BLAS_LAPACK_ERR : MEM_ERR; */
/* 	    } */
/* #if (PRECISION_SIZE == 4) */
/* 	  sgesvd(&jobz,&jobz,&n,&m,X,&n,d,V,&n,U,&m,work,&k,&info); */
/* #else */
/* 	  dgesvd(&jobz,&jobz,&n,&m,X,&n,d,V,&n,U,&m,work,&k,&info); */
/* #endif */
/* 	  if (info) */
/* 	    { */
/* 	      free(iwork); */
/* 	      free(X); */
/* 	      return BLAS_LAPACK_ERR; */
/* 	    } */
/* 	} */
/*       if (!vt)  */
/* 	bTransposeInPlace(V,n); */
/*       if (ut)  */
/* 	bTransposeInPlace(U,m); */
/*     } */
/*   free(work); */
/*   free(iwork); */
/*   free(X); */
/*   return NO_ERR; */
