#include "mathlib.h"
#include "math.h"

#include "modify_for_mex.h"

/*  //=======================================================   ansi-c*/
/*  // LAPACK routines -- all based on Intel MKL (or IMSL C Math library)   ansi-c*/
/*  //=======================================================   ansi-c*/
#if defined (INTELCMATHLIBRARY)
int lurgen(TSdmatrix *lu_dm, TSivector *pivot_dv, TSdmatrix *x_dm) {
/*     // PLU = x_dm from the LU decomposition of the input matrix x_dm where P is a permutation matrix, L is lower triangular with unit   ansi-c*/
/*     //   diagonal elements (lower trapezoidal if nrows>ncols) and U is upper triangular (upper trapezoidal if nrows<ncols).   ansi-c*/
/*     // L: (1) If nrows <= ncols, nrows-by-nrows .   ansi-c*/
/*     //    (2) If nrows > ncols, nrows-by-ncols (lower trapezoidal).   ansi-c*/
/*     // U: (1) If nrows <= ncols, nrows-by-ncols (upper trapezoidal).   ansi-c*/
/*     //    (2) If nrows > ncols, ncols-by-ncols.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  lu_dm: Stack L and U in this nrows-by-ncols matrix where the unit diagonal elements of L are not stored.   ansi-c*/
/*     //  pivot_dv: Optional.  An min(nrows, ncols) vector of index integers such that row i was interchanged with row pivot_dv->v[i].   ansi-c*/
/*     //    When NULL, this output argument is not exported (but computed anyway by the MKL hard-coded function).   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_dm: nrows-by-ncols general real matrix.   ansi-c*/

   int nrows, ncols, mindim, i,
       errflag=2;   /*  errflag=0 implies a successful execution.  But we start with 2 so as to let dgetrf_ export a correct flag.   ansi-c*/
   double *LU;

   lapack_int nrows2, ncols2, errflag2;
   lapack_int *pivot_p=NULL;

/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if ( !lu_dm || !x_dm )  fn_DisplayError(".../mathlib.c/lurgen(): The input arguments lu_dm and x_dm must be cretaed (memory-allocated)");
   else if ( ( (nrows=x_dm->nrows) != lu_dm->nrows) || ( (ncols=x_dm->ncols) != lu_dm->ncols) )
       fn_DisplayError(".../mathlib.c/lurgen(): Make sure the dimensions of the input matricies lu_dm and x_dm are the same");

   if ( !(x_dm->flag & M_GE) )
   {
      if (x_dm->flag & M_SU)   SUtoGE(x_dm);
      else if (x_dm->flag & M_SL)   SLtoGE(x_dm);
      else  fn_DisplayError(".../mathlib.c/lurgen(): Haven't got time to make M_UT, M_LT, and other to a general matrix M_GE");
   }
/*     //else if ( !(x_dm->flag & M_GE) )  fn_DisplayError(".../mathlib.c/lurgen(): The input arguments x_dm must be a general real matrix with the flag M_GE");   ansi-c*/


   mindim = _min(nrows, ncols);
   memcpy((LU=lu_dm->M), x_dm->M, nrows*ncols*sizeof(double));
   lu_dm->flag = M_UT;   /*  To make the lower part of lu_dm available, one must create another matrix and explicitly make it a unit lower matrix.   ansi-c*/

/*     //=== Calling the MKL function.   ansi-c*/
   nrows2 = nrows;
   ncols2 = ncols;
   errflag2 = errflag;
   pivot_p = tzMalloc(mindim, lapack_int);
   if (!pivot_dv)
      dgetrf(&nrows2, &ncols2, LU, &nrows2, pivot_p, &errflag2);
   else {
      if ( pivot_dv->n != mindim) fn_DisplayError("Make sure the dimension of the input vector pivot_dv is the minimum number of row number and column number of the input matrix x_dm");

      dgetrf(&nrows2, &ncols2, LU, &nrows2, pivot_p, &errflag2);
      for(i=0; i<mindim; i++)
        pivot_dv->v[i] = pivot_p[i];
   }
   swzFree(pivot_p);    /*  Frees the memory belonging to this function.   ansi-c*/
   errflag = errflag2;

   return( errflag );   /*  (1) If errflag = 0, success.  (2) If errorflag = -i, the ith parameter has an illegal value.   ansi-c*/
/*                         //(3) If errflag = i, u_{ii}=0.0.  The factorization is completed, but U is exactly singular.  Dividing   ansi-c*/
/*                         //      by 0.0 will occur if you use the factor U for solving a system of linear equations.   ansi-c*/
}
#else
/*  //No default routine yet.   ansi-c*/
#endif


#if defined (INTELCMATHLIBRARY)
int eigrsym(TSdvector *eval_dv, TSdmatrix *eVec_dm, const TSdmatrix *S_dm)
{
/*     // Outputs (dependent on Intel MKL):   ansi-c*/
/*     //   eval_dv:  _n-by-1 eigenvalues in ascending order;   ansi-c*/
/*     //   eVec_dm:  _n-by-_n eigenvalues -- if (eVec_m==NULL), no eigenvectors are computed; otherwise, S_dm = eVec_dm*diag(eval_dv)*inv(eVec_dm).   ansi-c*/
/*     //   errflag:  error flag.   ansi-c*/
/*     //------------   ansi-c*/
/*     // Inputs:   ansi-c*/
/*     //   S_dm:  _n-by_n real symmetric matrix.   ansi-c*/
/*     //   ansi-c*/
/*     // Eigenanalysis of real symmetric square matrix with all eigenvalues and, optionally, eigenvectors.   ansi-c*/
/*     //   Experts' opinion: do NOT use Cuppen's divide-and-conquer algorithm; instead, use QR algorithm, which I guess this algorithm uses.   ansi-c*/

   int n1, _n, errflag=2,   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dsyev_ export a correct flag.   ansi-c*/
       lwork;
   double *tmpd0_m = NULL,
          *work_p = NULL;

   lapack_int n2, lwork2, errflag2;

   if ( !S_dm || !(S_dm->flag & (M_SU | M_SL)) )  fn_DisplayError(".../mathlib.c/eigrsym():  input matrix (1) must be created (memory-alloacted) and (2) must be symmetric (either M_SU or M_SL)");
   if ( !eval_dv )  fn_DisplayError(".../mathlib.c/eigrsym():  input eigenvalue vector must be created (memory-allocated)");
   lwork = (n1=_n=S_dm->nrows)*BLOCKSIZE_FOR_INTEL_MKL;


/*     //=== Memory allocated for this function.   ansi-c*/
   tmpd0_m = tzMalloc(square(_n), double),
   work_p = tzMalloc(lwork, double);


/*     //---------------------------   ansi-c*/
/*     // Obtains eigenvalues and, optionally, eigenvectors.   ansi-c*/
/*     //---------------------------   ansi-c*/
   memcpy(tmpd0_m, S_dm->M, square(_n)*sizeof(double));

   n2 = n1;
   lwork2 = lwork;
   errflag2 = errflag;
   dsyev( (eVec_dm) ? "V" : "N", (S_dm->flag & M_SU) ? "U" : "L", &n2, tmpd0_m, &n2, eval_dv->v, work_p, &lwork2, &errflag2);
   errflag = errflag2;
   if (work_p[0]>lwork) printf("Warning for /mathlib.c/eigrsym(): needs at least %d workspace for good performance "
                                 "but lwork is allocated with only %d space!\n", (int)work_p[0], lwork);
   eval_dv->flag = V_DEF;
   if (eVec_dm) {
      memcpy(eVec_dm->M, tmpd0_m, square(_n)*sizeof(double));
      eVec_dm->flag = M_GE;
   }


/*     //---------------------------   ansi-c*/
/*     // Frees the allocated memory.   ansi-c*/
/*     //---------------------------   ansi-c*/
   tzDestroy(tmpd0_m);
   tzDestroy(work_p);

/*     //if (errflag<0) fn_DisplayError("/Subfolder: Calling eigrsym_decomp -- some element in input matrix has an illegal value");   ansi-c*/
/*     //else if (errflag>0) fn_DisplayError("/Subfolder: Calling eigrsym_decomp -- the factor U is exactly singular, so the solution cannot be computed");   ansi-c*/
   return (errflag);
}
#else
/*  //Not default routine yet.   ansi-c*/
#endif



#if defined (INTELCMATHLIBRARY)
int eigrgen(TSdzvector *vals_dzv, TSdzmatrix *rights_dzm, TSdzmatrix *lefts_dzm, const TSdmatrix *x_dm)
{
/*     //--- Eigenanalysis of real general (non-symmetric) square matrix with all eigenvalues and, optionally, eigenvectors. ---   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs (dependent on Intel MKL):   ansi-c*/
/*     //  vals_dzv->real->v: _n-by-1 real parts of eigenvalues;   ansi-c*/
/*     //  vals_dzv->imag->v: _n-by-1 imaginary parts of eigenvalues -- must be *initialized to zero* in this function;   ansi-c*/
/*     //  rights_dzm->real->M: if (rights_dzm==NULL), no right eigenvectors are computed; otherwise, _n-by-_n corresponding *real* parts of right eigenvectors column by column: A*v(j)=lambda(j)*v(j);   ansi-c*/
/*     //                       if (rights_dzm!=NULL), lefts_dzm->Mi must be *initialized to zero* in this function to get _n-by-_n *imaginary* parts of left eigenvectors corresponding to vals_dzv.   ansi-c*/
/*     //  lefts_dzm->imag->M:  if (lefts_dzm==NULL), no left eigenvectors are computed; otherwise, n-by-n corresponding *real* parts of left eigenvectors column by column: u(j)^H*A=lambda(j)*u(j)^H, where H means conjugate transpose;   ansi-c*/
/*     //                       if (lefts_dzm!=NULL), lefts_dzm->Mi must be *initialized to zero* in this function to get _n-by-_n *imaginary* parts of right eigenvectors corresponding to vals_dzv.   ansi-c*/
/*     //  returned errflag: error flag.  If errflag<0, some element in input matrix has an illegal value.   ansi-c*/
/*     //                                 If errflag>0, the QR algorithm failed to compute all the eigenvalues and no eigenvectors have been computed.   ansi-c*/
/*     //                                 if errflag=0, we have a successful decomposition.   ansi-c*/
/*     //------------   ansi-c*/
/*     // Inputs:   ansi-c*/
/*     //   x_dm:  _n-by_n real general (non-symmetric) matrix.   ansi-c*/

   int errflag=2,   /*  errflag=0 implies successful decomposition.  But we start with 2 so as to let dgeev_ export a correct flag.   ansi-c*/
       _n, lwork, n1, _i, _j;
   double *tmpd0_m=NULL,
          *work_p=NULL,
          *x_m=NULL,
          *evalr_v=NULL,
          *evali_v=NULL,
          *revecr_m=NULL, *reveci_m=NULL,     /*  NULL means that by default we dont' compute eigenvectors.   ansi-c*/
          *levecr_m=NULL, *leveci_m=NULL;

   lapack_int n2, lwork2, errflag2;

/*     //---------------------------   ansi-c*/
/*     // Checking dimensions, etc.   ansi-c*/
/*     //---------------------------   ansi-c*/
   if ( !x_dm || !vals_dzv )
      fn_DisplayError(".../mathlib.c/eigrgen():  Input square matrix x_dm and eigen value vectors vals_dzv must be created (memory allocated) before the call to this function");
   else {
      _n = x_dm->nrows;
      lwork = _n*BLOCKSIZE_FOR_INTEL_MKL;
      n1 = _n;
      tmpd0_m = tzMalloc(square(_n), double),     /*  @@Must be freed in this function.@@   ansi-c*/
      work_p = tzMalloc(lwork, double),           /*  @@Must be freed in this function.@@   ansi-c*/
      InitializeConstantVector_lf(vals_dzv->imag, 0.0);   /*  Imaginary part must be initialized to 0.0 to testing purposes later on.   ansi-c*/
/*        //   ansi-c*/
      x_m = x_dm->M;
      evalr_v = vals_dzv->real->v;
      evali_v = vals_dzv->imag->v;
   }
   if ( _n!=vals_dzv->real->n || _n!=x_dm->ncols ) fn_DisplayError(".../mathlib.c/eigrgen(): (1)input real matrix x_dm must be square; (2) the length of vals_dzv must match the dimension of x_dm");
   if (rights_dzm) {
      if ( _n!=rights_dzm->real->nrows || _n!=rights_dzm->real->ncols ) fn_DisplayError(".../mathlib.c/eigrgen(): rights_dzm must have the same dimension as the input square matrix");
      revecr_m = rights_dzm->real->M;   /*   (rights_dzm) ?  rights_dzm->real->M : NULL,   ansi-c*/
      rights_dzm->real->flag = M_GE;
      InitializeConstantMatrix_lf(rights_dzm->imag, 0.0);
      reveci_m = rights_dzm->imag->M;
   }
   if (lefts_dzm) {
      if ( _n!=lefts_dzm->real->nrows || _n!=lefts_dzm->real->ncols ) fn_DisplayError(".../mathlib.c/eigrgen(): lefts_dzm must have the same dimension as the input square matrix");
      levecr_m = lefts_dzm->real->M;   /*   (lefts_dzm) ?  lefts_dzm->real->M : NULL,   ansi-c*/
      lefts_dzm->real->flag = M_GE;
      InitializeConstantMatrix_lf(lefts_dzm->imag, 0.0);
      leveci_m = lefts_dzm->imag->M;
   }



/*     //---------------------------   ansi-c*/
/*     // Starts with x_m -- the matrix to be decomposed.   ansi-c*/
/*     //---------------------------   ansi-c*/
   memcpy(tmpd0_m, x_m, square(_n)*sizeof(double));

/*     //---------------------------   ansi-c*/
/*     // Obtains eigenvalues and, optionally, eigenvectors.   ansi-c*/
/*     //---------------------------   ansi-c*/
   n2 = n1;
   lwork2 = lwork;
   errflag2 = errflag;
   dgeev( (levecr_m) ? "V" : "N", (revecr_m) ? "V" : "N", &n2, tmpd0_m, &n2, evalr_v, evali_v,
                                                          levecr_m, &n2, revecr_m, &n2, work_p, &lwork2, &errflag2);
   errflag = errflag2;
   vals_dzv->real->flag = V_DEF;

/*     //---------------------------   ansi-c*/
/*     // Frees the allocated memory.   ansi-c*/
/*     //---------------------------   ansi-c*/
   if (work_p[0]>lwork) printf("Warning for /mathlib.c/eigrgen(): needs at least %d workspace for good performance "
                                 "but lwork is allocated with only %d space!\n", (int)work_p[0], lwork);
   tzDestroy(work_p);
   tzDestroy(tmpd0_m);

/*     //---------------------------   ansi-c*/
/*     // Checks error conditions.   ansi-c*/
/*     // Exports final results.   ansi-c*/
/*     //---------------------------   ansi-c*/
   if (errflag) return( errflag );
   else {
      if (revecr_m) {                    /*  Tested against Matlab output.  Works!  10/13/02.   ansi-c*/
         for (_j=0; _j<_n-1; _j++)
            if (evali_v[_j] && (evali_v[_j] == -evali_v[_j+1]))
               for (_i=0; _i<_n; _i++) {
                  reveci_m[_i+(_j+1)*_n] = -(reveci_m[_i+_j*_n]=revecr_m[_i+(_j+1)*_n]);
                  revecr_m[_i+(_j+1)*_n] = revecr_m[_i+_j*_n];
               }
      }
      if (levecr_m) {       /*  !!WARNINGS!!: Hasn't tested against any other established program, but it seems working.  10/13/02.   ansi-c*/
         for (_j=0; _j<_n-1; _j++)
            if (evali_v[_j] && (evali_v[_j] == -evali_v[_j+1]))
               for (_i=0; _i<_n; _i++) {
                  leveci_m[_i+(_j+1)*_n] = -(leveci_m[_i+_j*_n]=levecr_m[_i+(_j+1)*_n]);
                  levecr_m[_i+(_j+1)*_n] = levecr_m[_i+_j*_n];
               }
      }
      return( errflag );
   }
}
#else
/*  //Not default routine yet.   ansi-c*/
#endif

#if defined (INTELCMATHLIBRARY)
int chol(TSdmatrix *D_dm, TSdmatrix *S_dm, const char ul) {
/*     //?????????? Some of options are NOT tested yet.   ansi-c*/
/*     // Choleski decomposition of a symmetric, positive definite matrix S.  Intel MKL Lapack dependent code.   ansi-c*/
/*     // The fastest way for chol() is to let D = S, but D will be replaced by the Choleski factor.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  D: _n-by_n -- if ul=='U' or 'u',  D'*D = S where D is stored only in the upper triangular part;   ansi-c*/
/*     //                if ul=='L' or 'l',  D*D' = S where D is stored only in the lower triangular part.   ansi-c*/
/*     //     If D_dm->M == S_dm->M, D_dm->M (and S_dm->M) will be replaced by the triangular Choleski factor after the call to this function.   ansi-c*/
/*     //  errflag: error flag:  0:  successful;   ansi-c*/
/*     //                        -6: not symmetric (see mklman.pdf for other error return codes on ?potrf().   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  S:  _n-by-_n symmetric, positive definite matrix (whose only triangular part is used by dpotrf_).   ansi-c*/
/*     //  ul: if =='U' or 'u', D (NOT necessarily S unless D == S) is upper triangular; if =='L' or 'l', D (NOT necessarily S unless D == S) is lower triangular.   ansi-c*/

   int errflag=2, loc, nrows, _m, _i, _j;   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dpotrf_ export a correct flag.   ansi-c*/
   double *D, *S;

   lapack_int _m2, errflag2;

   if ( !D_dm || !S_dm )  fn_DisplayError(".../mathlib.c/chol():  L and R input square matricies must be created (memory allocated)");
   else {
      nrows = S_dm->nrows;
      _m = nrows;   /*  Used by Lapack.   ansi-c*/
      D = D_dm->M;
      S = S_dm->M;
   }
   if ( (nrows != D_dm->ncols) || (nrows != S_dm->nrows) || (nrows != S_dm->ncols) )   fn_DisplayError(".../mathlib.c/chol():  Make sure both R and L input matricies are square and have the same dimension");


/*     //=== Fills the triangular part that is used for Choleski decomposition.   ansi-c*/
   if ( D != S) {
      switch (ul) {
         case 'U': case 'u':
            if (S_dm->flag & M_SU) {
               for (_j=0; _j<nrows; _j++)
                  for (_i=0; _i<=_j; _i++) {
                     loc=mos(_i,_j,nrows);
                     D[loc] = S[loc];
                  }
               D_dm->flag = M_UT;
            }
            else if (S_dm->flag & M_SL) {
               for (_j=0; _j<nrows; _j++)
                  for (_i=0; _i<=_j; _i++)
                     D[mos(_i,_j,nrows)] = S[mos(_j,_i,nrows)];
               D_dm->flag = M_UT;
            }
            else
            {
/*                 //fn_DisplayError(".../mathlib.c/chol():  R input square matrix must be symmetric (and positive definite)");   ansi-c*/
               printf("\n ------- .../mathlib.c/chol():  R input square matrix must be symmetric!-------\n");
               return (-6);
            }
            _m2 = _m;
            errflag2 = errflag;
            dpotrf("U", &_m2, D, &_m2, &errflag2);
            errflag = errflag2;
            break;
         case 'L': case 'l':
            if (S_dm->flag & M_SL) {
               for (_j=0; _j<nrows; _j++) {
/*                    //for (_i=0; _i<_j; _i++)  D[_i+_j*nrows] = 0.0;   //Initializes the other part of D to zero so as to make it visible and readable.   ansi-c*/
                  for (_i=_j; _i<nrows; _i++) {
                     loc=mos(_i,_j,nrows);
                     D[loc] = S[loc];
                  }
               }
               D_dm->flag = M_LT;
            }
            else if (S_dm->flag & M_SU) {
/*                 //????????????? NOT teste yet for this option.   ansi-c*/
               for (_j=0; _j<nrows; _j++)
                  for (_i=_j; _i<nrows; _i++)
                     D[mos(_i,_j,nrows)] = S[mos(_j,_i,nrows)];
               D_dm->flag = M_LT;
            }
            else
            {
/*                 //fn_DisplayError(".../mathlib.c/chol():  R input square matrix must be symmetric (and positive definite)");   ansi-c*/
               printf("\n ------- .../mathlib.c/chol():  R input square matrix must be symmetric!-------\n");
               return (-6);
            }
/*              //??????NOT tested yet.   ansi-c*/
            _m2 = _m;
            errflag2 = errflag;
            dpotrf("L", &_m2, D, &_m2, &errflag2);
            errflag = errflag2;
            break;
         default:
            fn_DisplayError(".../mathlib.c/chol():  Input ul must be either 'U' or 'L'");
      }
   }
   else {
      if ( (ul=='U' || ul=='u') && (D_dm->flag & M_SU) ) {
        _m2 = _m;
        errflag2 = errflag;
        dpotrf("U", &_m2, D, &_m2, &errflag2);
        errflag = errflag2;
         D_dm->flag = M_UT;
      }
      else if ( (ul=='L' || ul=='l') && (D_dm->flag & M_SL) ) {
/*           //Tested.  It works!   ansi-c*/
        _m2 = _m;
        errflag2 = errflag;
        dpotrf("L", &_m2, D, &_m2, &errflag2);
        errflag = errflag2;
         D_dm->flag = M_LT;
      }
      else {
         printf("\nFatal Error: The input ul is %c and the flag D_dm->flag is %d", ul, D_dm->flag);
         fn_DisplayError(".../mathlib.c/chol():  When D==S, upper or lower triangular output must be consistent with upper or lower symmetric input; otherwise, use the option with D != S");
      }
   }
/*     //=== Choleski decomposition.   ansi-c*/
/*     // dpotrf_(((ul=='u') || (ul=='U')) ? "U" : "L", &_m, D, &_m, &errflag);   ansi-c*/
/*     //---   ansi-c*/
/*     // if (errflag<0) fn_DisplayError("Some element has an illegal value");   ansi-c*/
/*     // else if (errflag>0) fn_DisplayError("The leadding minor of some order, hence the entire matrix, is not positive definite");   ansi-c*/
   return (errflag);
}
#else
/*  //No default routine yet.   ansi-c*/
#endif


#if defined (INTELCMATHLIBRARY)
int invrtri(TSdmatrix *X_dm, TSdmatrix *A_dm, const char un)
{
/*     //Inverse of a real triangular matrix A.   ansi-c*/
/*     //The fastest way is to let X=A and A (and X) will be replaced by inv(A).   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  X: _n-by_n inverse of A;   ansi-c*/
/*     //  errflag: error flag (=0 means successful).   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A: _n-by-_n real triangular matrix.   ansi-c*/
/*     //  un:  if un=='U' or 'u', A is unit triangular; otherwise (un=='N' or 'n', A is not a unit triangular matrix.   ansi-c*/

   int _n, errflag=2;   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dgetri_ export a correct flag.   ansi-c*/
   double *X, *A;

   lapack_int _n2, errflag2;

   if ( !X_dm || !A_dm )  fn_DisplayError(".../mathlib.c/invrtri(): Both input matrices must be created (memory-allocated)");
   else if ( !(A_dm->flag & (M_UT | M_LT)) )  fn_DisplayError(".../mathlib.c/invrtri(): (1) R input matrix A must be given legal values; (2) A must be a real triangular matrix, i.e., M_UT or M_LT");
   else {
      X = X_dm->M;
      A = A_dm->M;
     _n=A_dm->nrows;
   }
   if ( (_n != A_dm->ncols) || (_n != X_dm->nrows) || (_n != X_dm->ncols) )
      fn_DisplayError(".../mathlib.c/invrtri(): both input and output matrices (1) must be square and (2) must have the same dimension");

   _n2 = _n;
   errflag2 = errflag;
   if (X==A) {
      dtrtri((A_dm->flag & M_UT) ? "U" : "L", (un=='U' || un=='u') ? "U" : "N", &_n2, X, &_n2, &errflag2);
      errflag = errflag2;
      if (errflag)  return (errflag);
   }
   else {
      memcpy(X, A, _n*_n*sizeof(double));
      dtrtri((A_dm->flag & M_UT) ? "U" : "L", (un=='U' || un=='u') ? "U" : "N", &_n2, X, &_n2, &errflag2);
      errflag = errflag2;
      if (errflag)  return (errflag);
      else  X_dm->flag = A_dm->flag;
   }

   return errflag;    /*  (1) If errflag = 0, success.  (2) If errorflag = -i, the ith parameter has an illegal value.   ansi-c*/
/*                       //(3) If errflag = i, the ith diagonal element of A is zero, A is singular, and the inversion   ansi-c*/
/*                       //      could not be completed.   ansi-c*/
}
#else
/*  //No default routine yet.   ansi-c*/
#endif


#if defined (INTELCMATHLIBRARY)
int invspd(TSdmatrix *X_dm, TSdmatrix *A_dm, const char ul)
{
/*     //Inverse of a symmetric positive matrix A.   ansi-c*/
/*     //Fastest way: let X=A.  Then, A (and X) will be replaced by inv(A).   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  X: _n-by_n inverse of A;   ansi-c*/
/*     //  errflag: error flag (=0 means successful).   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A: _n-by-_n symmetric positive matrix.   ansi-c*/
/*     //  ul:  if ul=='U' or 'u', only upper part of A is used; otherwise (un=='L' or 'l', only lower part of A is used.   ansi-c*/

   int _n, errflag=2;   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dgetri_ export a correct flag.   ansi-c*/
   double *X, *A;

   lapack_int _n2, errflag2;

   if ( !X_dm || !A_dm )  fn_DisplayError(".../mathlib.c/invspd(): Both input matrices must be created (memory-allocated)");
   else if ( !(A_dm->flag & (M_SU | M_SL)) )  fn_DisplayError(".../mathlib.c/invspd(): (1) R input matrix A must be given legal values; (2) A must be symmetric, positive-definite, i.e., M_SU or M_SL");
   else {
      X = X_dm->M;
      A = A_dm->M;
     _n=A_dm->nrows;
   }


   if (X==A) {
      if ( (_n != A_dm->ncols) )
         fn_DisplayError(".../mathlib.c/invspd(): input matrix (1) must be square and (2) must have the same dimension");
/*        //=== Choleski decomposition.   ansi-c*/
      _n2 = _n;
      errflag2 = errflag;
      dpotrf(((ul=='U') || (ul=='u')) ? "U" : "L", &_n2, X, &_n2, &errflag2);
      errflag = errflag2;
      if (errflag)  return (errflag);
/*        //=== Takes inverse.   ansi-c*/
      dpotri(((ul=='U') || (ul=='u')) ? "U" : "L", &_n2, X, &_n2, &errflag2);
      errflag = errflag2;
      A_dm->flag = ((ul=='U') || (ul=='u')) ? M_SU : M_SL;
      return (errflag);
/*        //---   ansi-c*/
/*        // if (errflag<0) fn_DisplayError("Some element has an illegal value");   ansi-c*/
/*        // else if (errflag>0) fn_DisplayError("Not symmetric positive definite or matrix inversion cannot be computed");   ansi-c*/
   }
   else {
      if ( (_n != A_dm->ncols) || (_n != X_dm->nrows) || (_n != X_dm->ncols) )
         fn_DisplayError(".../mathlib.c/invspd(): both input and output matrices (1) must be square and (2) must have the same dimension");
      memcpy(X, A, _n*_n*sizeof(double));
/*        //=== Choleski decomposition.   ansi-c*/
      _n2 = _n;
      errflag2 = errflag;
      dpotrf(((ul=='U') || (ul=='u')) ? "U" : "L", &_n2, X, &_n2, &errflag2);
      errflag = errflag2;
      if (errflag)  return (errflag);
/*        //=== Takes inverse.   ansi-c*/
      dpotri(((ul=='U') || (ul=='u')) ? "U" : "L", &_n2, X, &_n2, &errflag2);
      errflag = errflag2;
      X_dm->flag = ((ul=='U') || (ul=='u')) ? M_SU : M_SL;
      return (errflag);
/*        //---   ansi-c*/
/*        // if (errflag<0) fn_DisplayError("Some element has an illegal value");   ansi-c*/
/*        // else if (errflag>0) fn_DisplayError("Not symmetric positive definite or matrix inversion cannot be computed");   ansi-c*/
   }
}
#else
/*  //No default routine yet.   ansi-c*/
#endif



#if defined (INTELCMATHLIBRARY)
int invrgen(TSdmatrix *X_dm, TSdmatrix *A_dm)
{
/*     //Inverse of a general real matrix A.   ansi-c*/
/*     //If X=A, A (and X) will be replaced by inv(A).   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  X: _n-by_n inverse of A;   ansi-c*/
/*     //  errflag: error flag (=0 means successful).   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A: _n-by-_n real general matrix.   ansi-c*/
   int _n, errflag=2,   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dgetri_ export a correct flag.   ansi-c*/
       lwork;   /*  Used when calling LAPACK.   ansi-c*/
   double *X, *A,
          *work;   /*  Used when calling LAPACK.   ansi-c*/

   lapack_int _n2, errflag2, *ipivot, lwork2;

   if ( !X_dm || !A_dm )  fn_DisplayError(".../mathlib.c/invrgen(): Both input matrices must be created (memory-allocated)");
   else if ( !(A_dm->flag & M_GE) )  fn_DisplayError(".../mathlib.c/invrgen(): (1) R input matrix A must be given legal values; (2) A must be a general matrix, i.e., M_GE");
   else {
      X = X_dm->M;
      A = A_dm->M;
      ipivot = tzMalloc((_n=A_dm->nrows), lapack_int);
      work = tzMalloc((lwork=_n*BLOCKSIZE_FOR_INTEL_MKL), double);
   }
   if ( (_n != A_dm->ncols) || (_n != X_dm->nrows) || (_n != X_dm->ncols) )
      fn_DisplayError(".../mathlib.c/invrgen(): both input and output matrices (1) must be square and (2) have the same dimension");


   if (X==A) {
     _n2 = _n;
     errflag2 = errflag;
     dgetrf(&_n2, &_n2, A, &_n2, ipivot, &errflag2);
     errflag = errflag2;
      if (errflag) {
/*  //         A_dm->flag = M_UNDEF;   ansi-c*/
         swzFree(ipivot);
         swzFree(work);
         return errflag;
      }
      lwork2 = lwork;
      dgetri(&_n2, A, &_n2, ipivot, work, &lwork2, &errflag2);
      errflag = errflag2;
      if (work[0]>lwork) printf("Warning for /mathlib.c/invrgen(); when calling MKL dgetri_(), we need at least %d workspace for good performance "
                                    "but lwork is allocated with only %d space!\n", (int)work[0], lwork);
      if (errflag) {
         swzFree(ipivot);
         swzFree(work);
         return (errflag);   /*  A_dm->flag = M_UNDEF;   ansi-c*/
      }
   }
   else {
      memcpy(X, A, _n*_n*sizeof(double));
      _n2 = _n;
      errflag2 = errflag;
      dgetrf(&_n2, &_n2, X, &_n2, ipivot, &errflag2);
      errflag = errflag2;
      if (errflag) {
/*  //         X_dm->flag = M_UNDEF;   ansi-c*/
         swzFree(ipivot);
         swzFree(work);
         return errflag;
      }
      lwork2 = lwork;
      dgetri(&_n2, X, &_n2, ipivot, work, &lwork2, &errflag2);
      errflag = errflag2;
      if (work[0]>lwork) printf("Warning for /mathlib.c/invrgen(); when calling MKL dgetri_(), we need at least %d workspace for good performance "
                                    "but lwork is allocated with only %d space!\n", (int)work[0], lwork);
      if (errflag) {
         swzFree(ipivot);
         swzFree(work);
         return (errflag);   /*  X_dm->flag = M_UNDEF;   ansi-c*/
      }
      else  X_dm->flag = A_dm->flag;
   }
/*     //=== Frees memory allocated in this function.   ansi-c*/
   swzFree(ipivot);
   swzFree(work);

   return errflag;    /*  (1) If errflag = 0, success.  (2) If errorflag = -i, the ith parameter has an illegal value.   ansi-c*/
/*                       //(3) If errflag = i, U_{ii}=0.0.  The LU factorization U is literally singular and the inversion   ansi-c*/
/*                       //      could not be completed.   ansi-c*/
}
#else
/*  //No default routine yet.   ansi-c*/
#endif



#if defined (INTELCMATHLIBRARY)
int BdivA_rrect(TSdmatrix *X_dm, const TSdmatrix *B_dm, const char lr, const TSdmatrix *A_dm)
{
/*     //This routine solves left division (\) problme AX=B (X=A\B).  For XA=B (right division problem: X=B/A), we first transpose   ansi-c*/
/*     //  it to A'*X'=B', then solve out X' = A'\B' as a left-division problem, and finally transpose X' back to X.   ansi-c*/
/*     //It handles cases with _m>=_n.  For _n>_m, we have an infinite solution.  To get one solution, take the _m-by-_m nonsigular   ansi-c*/
/*     //  part of A_dm and get the solution for the _m-by-_r part of X with \ or the _r-by-_m part of X with /.  The (_n-_m)-by-_r part   ansi-c*/
/*     //  or _r-by-(_n-_m) part of X can simply be filled with 0.0.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  X = A\B or B/A where X is _n-by_r if \ (AX=B) or _r-by-_n if / (XA=B).   ansi-c*/
/*     //  Returns info (if info==0, the execution is successful; if info == -i, the ith parameter has an illegal value.)   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A:  _m-by-_n real rectangular (rrect) matrix if \ (AX=B) or   ansi-c*/
/*     //      _n-by-_m real rectangular (rrect) matrix if / (XA=B).   ansi-c*/
/*     //  B:  _m-by-_r real rectangular (rrect) matrix if \ (AX=B) or   ansi-c*/
/*     //      _r-by-_m real rectangular (rrect) matrix if / (XA=B).   ansi-c*/
/*     //  lr:  if lr='\\', left division \ is performed; if lr='/', right division / is performed.   ansi-c*/

   int _m, _n, _r, i,   /*  mn_max, mn_min,   ansi-c*/
       lwork, _i, info = -2,
       *jpvt_p = NULL;
   double *A, *B, *X,
          *qra_p = NULL,    /*  QR decomposion for A_dm.   ansi-c*/
          *qrb_p = NULL,    /*  QR decomposion for B_dm.   ansi-c*/
          *tau_p = NULL,
          *work_p = NULL;

   lapack_int _m2, _n2, _r2, lwork2, info2;
   lapack_int *jpvt_p2 = NULL;

   if (!A_dm || !(A_dm->flag & M_GE) || !B_dm || !(B_dm->flag &M_GE))   fn_DisplayError(".../mathlib.c/BdivA_rrect(): both input matricies A_dm and B_dm must be (a) created (allocated memory) and (b) given legal values for all elements (in other words, the flag M_GE must exist)");
   if (!X_dm)   fn_DisplayError(".../mathlib.c/BdivA_rrect(): output matrix X_dm must be created (allocated memory)");
   if (lr=='/') {
      if ( (_n=A_dm->nrows) != X_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rrect(): 1st dim of A_dm and 2nd dim of X_dm must exactly match for / (right division)");
      if ( (_m=A_dm->ncols) != B_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rrect(): 2nd dim of A_dm and 2nd dim of B_dm must exactly match for / (right division)");
      if ( (_r=B_dm->nrows) != X_dm->nrows )  fn_DisplayError(".../mathlib.c/BdivA_rrect(): 1st dim of B_dm and 1st dim of X_dm must exactly match for / (right division)");
      if ( _m < _n)  fn_DisplayError(".../mathlib.c/BdivA_rrect(): A_dm->nrows must be <= A_dm->ncols for / (right division).  Otherwise, take the nonsigular square part of A_dm and use BdivA_rrect()");
   }
   else if (lr=='\\') {
      if ( (_m=A_dm->nrows) != B_dm->nrows )  fn_DisplayError(".../mathlib.c/BdivA_rrect(): 1st dim of A_dm and 1st dim of B_dm must exactly match for \\ (left division)");
      if ( (_n=A_dm->ncols) != X_dm->nrows )  fn_DisplayError(".../mathlib.c/BdivA_rrect(): 2nd dim of A_dm and 1st dim of X_dm must exactly match for \\ (left division)");
      if ( (_r=B_dm->ncols) != X_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rrect(): 2nd dim of B_dm and 2nd dim of X_dm must exactly match for \\ (left division)");
      if ( _m < _n)  fn_DisplayError(".../mathlib.c/BdivA_rrect(): A_dm->nrows must be >= A_dm->ncols for \\ (left division). Otherwise, take the nonsigular square part of A_dm and use BdivA_rrect()");
   }
   else  fn_DisplayError(".../mathlib.c/BdivA_rrect(): input character lr must be / (right division) or \\\\ (left division)");
   A = A_dm->M;
   B = B_dm->M;
   X = X_dm->M;


/*  //   lwork = ((mn_max = _m>_n ? _m : _n)>_r ? nm_max : _r)*BLOCKSIZE_FOR_INTEL_MKL;   ansi-c*/
/*  //   mn_min = _m<_n ? _m : _n;   ansi-c*/


/*     //=== Memory allocation for this function only.   ansi-c*/
   qra_p = tzMalloc(_m*_n, double);
/*  //   qrb_p = tzMalloc((mn_max = _m>_n?_m:_n)*_r, double); //DDDDebug: seems requiring _m>_n, but this may not be the case.   ansi-c*/
   qrb_p = tzMalloc(_m*_r, double);  /*  Note that _m>=_n.   ansi-c*/
   jpvt_p = tzMalloc(_n, int);
   jpvt_p2 = tzMalloc(_n, lapack_int);
   tau_p = tzMalloc(_n, double);
/*  //   work_p = tzMalloc(lwork, double);   ansi-c*/


/*   //=== Making copies of input matrices A and B */
   if (lr=='/')      /*  Right division.  In this case, qra_p is _m-by-_n (transpose of A_dm), and qrb_p is max(_m,_n)-by-_r (transpose of B_dm).   ansi-c*/
      for (_i=0; _i<_m; _i++) {
         cblas_dcopy(_n, &A[_i*_n], 1, &qra_p[_i], _m);
         cblas_dcopy(_r, &B[_i*_r], 1, &qrb_p[_i], _m);
      }
   else {
      memcpy(qra_p, A, _m*_n*sizeof(double));     /*  qra_p is _m-by-_n.   ansi-c*/
      memcpy(qrb_p, B, _m*_r*sizeof(double));     /*  qrb_p is _m-by-_r.   ansi-c*/
/*  //      for (_i=0; _i<_r; _i++)  cblas_dcopy(_m, B+_i*_m, 1, qrb_p+_i*mn_max, 1);    //qrb_p is max(_m,_n)-by-_r.   ansi-c*/
   }


/*     //=== Computes the QR factorization of a general m by n matrix with column pivoting using Level 3 BLAS.   ansi-c*/
   work_p = tzMalloc(lwork=_n*BLOCKSIZE_FOR_INTEL_MKL, double);
   _m2 = _m;
   _n2 = _n;
   lwork2 = lwork;
   info2 = info;
   for(i=0; i<_n; i++)
     jpvt_p2[i] = jpvt_p[i];
   dgeqp3(&_m2,&_n2,qra_p,&_m2,jpvt_p2,tau_p,work_p,&lwork2,&info2);
   info = info2;
   for(i=0; i<_n; i++)
     jpvt_p[i] = jpvt_p2[i];
   if (work_p[0]>lwork) printf("Warning for /mathlib.c/BdivA_rrect(); when calling MKL dgeqp3_(), we need at least %d workspace for good performance "
                                 "but lwork is allocated with only %d space!\n", (int)work_p[0], lwork);
   tzDestroy(work_p);
   if (info)  return (info);    /*  If info==0, the execution is successful; if info == -i, the ith parameter has an illegal value.   ansi-c*/

/*     //=== Multiplies a real matrix by the orthogonal matrix Q of the QR factorization formed by dgeqp3_.   ansi-c*/
   work_p = tzMalloc(lwork=_r*BLOCKSIZE_FOR_INTEL_MKL, double);
   _m2 = _m;
   _n2 = _n;
   _r2 = _r;
   lwork2 = lwork;
   info2 = info;
   dormqr("L","T",&_m2,&_r2,&_n2,qra_p,&_m2,tau_p,qrb_p,&_m2,work_p,&lwork2,&info2);
   info = info2;
   if (work_p[0]>lwork) printf("Warning for /mathlib.c/BdivA_rrect(); when calling MKL dormqr_(), we need at least %d workspace for good performance "
                                 "but lwork is allocated with only %d space!\n", (int)work_p[0], lwork);
/*     //dormqr_("L","T",&_m,&_r,&mn_min,qra_p,&_m,tau_p,qrb_p,&mn_max,work_p,&lwork,&info);   ansi-c*/
   tzDestroy(work_p)
   if (info)  return (info);    /*  If info==0, the execution is successful; if info == -i, the ith parameter has an illegal value.   ansi-c*/


/*     //=== Solves a matrix equation R*x = C (one matrix operand is triangular).  Note that the dim of X is n-by-r for \ or r-by-n for /   ansi-c*/
   cblas_dtrsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,_n,_r,1.0,qra_p,_m,qrb_p,_m);
   if (lr=='/')      /*  Right division.  In this case, qra_p is _m-by-_n (transpose of A_dm), and qrb_p is max(_m, _n)-by-_r (transpose of B_dm).   ansi-c*/
      for (_i=0; _i<_n; _i++)  cblas_dcopy(_r, &qrb_p[_i], _m, &X[(jpvt_p[_i]-1)*_r], 1);        /*  Copying the transpose of the _n-by-_r leading part of qrb_p.   ansi-c*/
/*  //         cblas_dcopy(_r, &qrb_p[_i], mn_max, &X[(jpvt_p[_i]-1)*_r], 1);       //Copying the transpose of the _n-by-_r leading part of qrb_p.   ansi-c*/

   else      /*qrb_p is max(_m, _n)-by-_r.*/
      for (_i=0; _i<_n; _i++)  cblas_dcopy(_r, &qrb_p[_i], _m, &X[jpvt_p[_i]-1], _n);
/*  //         cblas_dcopy(_r, &qrb_p[_i], mn_max, &X[jpvt_p[_i]-1], _n);   ansi-c*/
   X_dm->flag = M_GE;

/*     //=== Destroyes memory allocated for this function only.   ansi-c*/
   tzDestroy(qra_p);
   tzDestroy(qrb_p);
   tzDestroy(jpvt_p);
   tzDestroy(jpvt_p2);
   tzDestroy(tau_p);
/*  //   tzDestroy(work_p);   ansi-c*/

   return (0);
}
#else
/*  //No default routine yet.  7 Oct 2003   ansi-c*/
#endif




#if defined (INTELCMATHLIBRARY)
int BdivA_rgens(TSdmatrix *X_dm, const TSdmatrix *B_dm, const char lr, const TSdmatrix *A_dm)
{
/*     //Unlike BdivA_rrect(), this routine deals with only the general square matrix A_dm.  It solves left division (\) problme   ansi-c*/
/*     //  AX=B (X=A\B).  For XA=B (right division problem: X=B/A), we first transpose it to A'*X'=B', then solve out X' = A'\B'   ansi-c*/
/*     //  as a left-division problem, and finally transpose X' back to X.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  X = A\B or B/A where X is _m-by_r if \ (AX=B) or _r-by-_m if / (XA=B).   ansi-c*/
/*     //  Returns info (if info==0, the execution is successful; if info == -i, the ith parameter has an illegal value.)   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A:  _m-by-_m real general square (rgens) matrix.   ansi-c*/
/*     //  B:  _m-by-_r real general square (rgens) matrix if \ (AX=B) or   ansi-c*/
/*     //      _r-by-_m real general square (rgens) matrix if / (XA=B).   ansi-c*/
/*     //  lr:  if lr='\\', left division \ is performed; if lr='/', right division / is performed.   ansi-c*/

   int _m, _r, m2,
       _i, info = -2;
   double *A, *B, *X,
          *Atran_p = NULL,    /*  Transpose of A if right division / takes place.   ansi-c*/
          *Btran_p = NULL,     /*  Transpose of B if right division / takes place.   ansi-c*/
          *W = NULL;    /*  Duplicate copy of A when left division \ is used.  This will be replaced by LU decomposition.   ansi-c*/
/*  //          *tau_p = NULL,   ansi-c*/
/*  //          *work_p = NULL;   ansi-c*/


   lapack_int _mlap, _rlap, info2, *ipiv_p = NULL;


   if (!A_dm || !(A_dm->flag & M_GE) || !B_dm || !(B_dm->flag & M_GE))   fn_DisplayError(".../mathlib.c/BdivA_rgens(): both input matricies A_dm and B_dm must be (a) created (allocated memory) and (b) given legal values for all elements (in other words, the flag M_GE must exist)");
   if (!X_dm)   fn_DisplayError(".../mathlib.c/BdivA_rgens(): output matrix X_dm must be created (allocated memory)");
   if (lr=='/') {
      if ( (_m=A_dm->nrows) != X_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): 1st dim of A_dm and 2nd dim of X_dm must exactly match for / (right division)");
      if ( _m != A_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): input matrix A_dm must be square.  For a nonsqaure matrix, use BdivA_rrect()");
      if ( _m != B_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): 2nd dim of A_dm and 2nd dim of B_dm must exactly match for / (right division)");
      if ( (_r=B_dm->nrows) != X_dm->nrows )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): 1st dim of B_dm and 1st dim of X_dm must exactly match for / (right division)");

   }
   else if (lr=='\\') {
      if ( (_m=A_dm->nrows) != B_dm->nrows )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): 1st dim of A_dm and 1st dim of B_dm must exactly match for \\ (left division)");
      if ( _m != A_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): input matrix A_dm must be square.  For a nonsqaure matrix, use BdivA_rrect()");
      if ( _m != X_dm->nrows )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): 2nd dim of A_dm and 1st dim of X_dm must exactly match for \\ (left division)");
      if ( (_r=B_dm->ncols) != X_dm->ncols )  fn_DisplayError(".../mathlib.c/BdivA_rgens(): 2nd dim of B_dm and 2nd dim of X_dm must exactly match for \\ (left division)");
   }
   else  fn_DisplayError(".../mathlib.c/BdivA_rgens(): input character lr must be / (right division) or \\\\ (left division)");
   A = A_dm->M;
   B = B_dm->M;
   X = X_dm->M;



   if (lr=='/') {
/*        //Right divistion /.   ansi-c*/
/*        //=== Memory allocation for this function only.   ansi-c*/
      ipiv_p = tzMalloc(_m, lapack_int);
      Atran_p = tzMalloc(square(_m), double);
      Btran_p = tzMalloc(_m*_r, double);

      for (_i=0; _i<_m; _i++) {
         cblas_dcopy(_m, A+_i*_m, 1, Atran_p+_i, _m);    /*  Copying the transpose of A to Atran.   ansi-c*/
         cblas_dcopy(_r, B+_i*_r, 1, Btran_p+_i, _m);    /*  Copying the transpose of B (_r-by-_m) to Btran (_m-by-_r);   ansi-c*/
      }
      _mlap = _m;
      _rlap = _r;
      info2 = info;
      dgesv(&_mlap, &_rlap, Atran_p, &_mlap, ipiv_p, Btran_p, &_mlap, &info2);
      info = info2;
      for (_i=0; _i<_r; _i++)  cblas_dcopy(_m, Btran_p+_i*_m, 1, X+_i, _r);   /*  Copying the transpose of Btran(_m-by-_r) to X (_r-by-_m);   ansi-c*/
      X_dm->flag = M_GE;

/*        //=== Destroyes memory allocated for this function only.   ansi-c*/
      tzDestroy(ipiv_p);
      tzDestroy(Atran_p);
      tzDestroy(Btran_p);

      return (info);
   }
   else {
/*        //Left division \.   ansi-c*/
/*        //=== Memory allocation for this function only.   ansi-c*/
      ipiv_p = tzMalloc(_m, lapack_int);
      W = tzMalloc(m2=square(_m), double);

      memcpy(X, B, _m*_r*sizeof(double));
      memcpy(W, A, m2*sizeof(double));
      _mlap = _m;
      _rlap = _r;
      info2 = info;
      dgesv(&_mlap, &_rlap, W, &_mlap, ipiv_p, X, &_mlap, &info2);
      info = info2;
      X_dm->flag = M_GE;

/*        //=== Destroyes memory allocated for this function only.   ansi-c*/
      tzDestroy(ipiv_p);
      tzDestroy(W);

      return (info);    /*  If info==0, the execution is successful; if info == -i, the ith parameter has an illegal value;   ansi-c*/
/*                         //  if info==i, U(i,i) in LU decomposition is exactly zero. The factorization has been completed, but   ansi-c*/
/*                         //  the factor U is exactly singular, so the solution could not be computed.   ansi-c*/
   }
}
/*  //---   ansi-c*/
int bdivA_rgens(TSdvector *x_dv, const TSdvector *b_dv, const char lr, const TSdmatrix *A_dm)
{
/*     //Use bdivA_rgens() where x_dv and b_dv are now both vectors and A_dm is a square matrix.   ansi-c*/
/*     //Unlike bdivA_rrect(), this routine deals with only the general square matrix A_dm.  It solves left division (\) problme   ansi-c*/
/*     //  Ax=b (x=A\b).  For xA=b (right division problem: x=B/b), we first transpose it to A'*x'=b', then solve out x' = A'\b'   ansi-c*/
/*     //  as a left-division problem, and finally transpose x' back to x.   ansi-c*/
/*     //   ansi-c*/
/*     //If x_dv->v = b_dv->v.  Then, x_dv->v will be replaced by new values.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs: Solving Ax = b or xA = b.   ansi-c*/
/*     //  x = A\b or b/A where x is _m-by-1 if \ (Ax=b) or 1-by-_m if / (xA=b).   ansi-c*/
/*     //  Returns info (if info==0, the execution is successful; if info == -i, the ith parameter has an illegal value.)   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A:  _m-by-_m real general square (rgens) matrix.   ansi-c*/
/*     //  b:  _m-by-1 vector if \ (Ax=b) or   ansi-c*/
/*     //      1-by-_m vector if / (xA=b).   ansi-c*/
/*     //  lr:  if lr='\\', left division \ is performed; if lr='/', right division / is performed.   ansi-c*/

   int _m, m2,
       _r = 1,
       _i, info = -2;
   double *A, *b, *x,
          *Atran_p = NULL,    /*  Transpose of A if right division / takes place.   ansi-c*/
          *W = NULL;    /*  Duplicate copy of A when left division \ is used.  This will be replaced by LU decomposition.   ansi-c*/

   lapack_int _mlap, _rlap, info2, *ipiv_p = NULL;

   if (!A_dm || !(A_dm->flag & M_GE) || !b_dv || !b_dv->flag)   fn_DisplayError("mathlib.c/bdivA_rgens(): Both A_dm and b_dv must be (a) created (allocated memory) and (b) given legal values for all elements (in other words, the flag M_GE must exist)");
   if (!x_dv)   fn_DisplayError("mathlib.c/bdivA_rgens(): output vector x_dv must be created (allocated memory)");
   if ( b_dv->n != x_dv->n )  fn_DisplayError("mathlib.c/bdivA_rgens(): The dim of b_dv and the dim of x_dv must exactly match");
   if ( (_m=A_dm->nrows) != x_dv->n )  fn_DisplayError("mathlib.c/bdivA_rgens(): Number of rows of A_dm and the dim of x_dv must exactly match");
   if ( _m != A_dm->ncols )  fn_DisplayError("mathlib.c/bdivA_rgens(): Input matrix A_dm must be square");

   A = A_dm->M;
   b = b_dv->v;
   x = x_dv->v;

   if (lr=='/') {
/*        //Right divistion /.   ansi-c*/
/*        //=== Memory allocation for this function only.   ansi-c*/
      ipiv_p = tzMalloc(_m, lapack_int);
      Atran_p = tzMalloc(square(_m), double);

      for (_i=0; _i<_m; _i++)
         cblas_dcopy(_m, A+_i*_m, 1, Atran_p+_i, _m);    /*  Copying the transpose of A to Atran.   ansi-c*/
      if (x_dv != b_dv )  memcpy(x, b, _m*sizeof(double));
      _mlap = _m;
      _rlap = _r;
      info2 = info;
      dgesv(&_mlap, &_rlap, Atran_p, &_mlap, ipiv_p, x, &_mlap, &info2);
      info = info2;
      x_dv->flag = V_DEF;

/*        //=== Destroyes memory allocated for this function only.   ansi-c*/
      tzDestroy(ipiv_p);
      tzDestroy(Atran_p);

      return (info);
   }
   else {
/*        //Left division \.   ansi-c*/
/*        //=== Memory allocation for this function only.   ansi-c*/
      ipiv_p = tzMalloc(_m, lapack_int);
      W = tzMalloc(m2=square(_m), double);

      if (x_dv != b_dv )  memcpy(x, b, _m*sizeof(double));
      memcpy(W, A, m2*sizeof(double));
      _mlap = _m;
      _rlap = _r;
      info2 = info;
      dgesv(&_mlap, &_rlap, W, &_mlap, ipiv_p, x, &_mlap, &info2);
      info = info2;
      x_dv->flag = V_DEF;

/*        //=== Destroyes memory allocated for this function only.   ansi-c*/
      tzDestroy(ipiv_p);
      tzDestroy(W);

      return (info);    /*  If info==0, the execution is successful; if info == -i, the ith parameter has an illegal value;   ansi-c*/
/*                         //  if info==i, U(i,i) in LU decomposition is exactly zero. The factorization has been completed, but   ansi-c*/
/*                         //  the factor U is exactly singular, so the solution could not be computed.   ansi-c*/
   }
}
#else
/*  //No default routine yet.  7 Oct 2003   ansi-c*/
#endif



#if defined (INTELCMATHLIBRARY)
void Aldivb_spd(TSdvector *x_dv, TSdmatrix *A_dm, TSdvector *b_dv, char an) {
/*     //??????? Some options (e.g., whe A_dm is M_SL) are NOT tested yet.   ansi-c*/
/*     //Output x = A\b where x_dv is an _n-by-1 vector.   ansi-c*/
/*     //  Fastest way is to let x_dv->v = b_dv->v.  Then, x_dv->v will be replaced by new values.   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A:  _n-by-_n symmetric, positive definite matrix.   ansi-c*/
/*     //  b:  _n-by-1 vector.   ansi-c*/
/*     //  an: 'N' or 'n': A_dm will NOT be altered and another matrix will be allocated and destroyed in this function.   ansi-c*/
/*     //      Otherwise ('A' or 'a'): A_dm will be altered and A_dm = U if A_dm->flag = M_SU where U is an upper triagular Choleski such that U'*U = A;   ansi-c*/
/*     //                                                    or A_dm = L if A_dm->flag = M_SL (but != M_SU) where L is a lower triangular Choleski such that L*L' = A.   ansi-c*/
/*     //   ansi-c*/
/*     // Note I:   Intel MLK cblas_dtrsv() does not test for singularity or near-singulariy of the system.   ansi-c*/
/*     //   Such tests must be performed before calling this BLAS routine.   ansi-c*/
/*     // Note II:  If x_dv->v = b_dv->v, x_dv->v will be replaced by new values.   ansi-c*/
/*     // Note III: If an != 'N' or 'n', A_dm will be replaced by U if A_dm->M_SU where U'*U = A or by L if A_dm->M_SL (but != M_SU) where L*L'=A.   ansi-c*/

   int errflag=2, nrows, nels;   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dpotrf_ export a correct flag.   ansi-c*/
   double *A, *W=NULL, *x, *b;

   lapack_int nrows2, errflag2;

   if ( !A_dm || !b_dv || !x_dv )  fn_DisplayError(".../mathlib.c/Aldivb_spd():  All input matrices or vectors must be created (memory allocated)");
   nrows = A_dm->nrows;
   nels = square(nrows);
   A = A_dm->M;
   x = x_dv->v;
   b= b_dv->v;
   if ( nrows != A_dm->ncols )  fn_DisplayError(".../mathlib.c/Aldivb_spd():  L input matrix must be square");
   if ( !A_dm->flag || !b_dv->flag )   fn_DisplayError(".../mathlib.c/Aldivb_spd():  L input matrix and vector must be given legal values");
   if ( (an=='N') || (an=='n') )  {
      W = tzMalloc(nels, double);
      memcpy(W, A, nels*sizeof(double));
   }
   else if ( (an=='A') || (an=='a') )  W = A;
   else  fn_DisplayError(".../mathlib.c/Aldivb_spd(): passing charecter an must be A, a, N, or n");


   if (A_dm->flag & M_SU) {
      nrows2 = nrows;
      errflag2 = errflag;
      dpotrf("U", &nrows2, W, &nrows2, &errflag2);   /*  Choleski.  U'*U = W where W will be replaced by upper triangular U.   ansi-c*/
      errflag = errflag2;
      if (errflag)  fn_DisplayError(".../mathlib.c/Aldivb_spd():  Error when calling Choleski dpotrf_().  Check if the L input matrix A_dm is positive definite or has legal values");
      if (x==b)  {
/*           //=== Solving for A*x=b.   ansi-c*/
         cblas_dtrsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, nrows, W, nrows, x, 1);
         cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, nrows, W, nrows, x, 1);
      }
      else {
         memcpy(x, b, nrows*sizeof(double));
         cblas_dtrsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, nrows, W, nrows, x, 1);
         cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, nrows, W, nrows, x, 1);
         x_dv->flag = V_DEF;
      }
      if ( (an!='N') && (an!='n') )  A_dm->flag = M_UT;
   }
   else if (A_dm->flag & M_SL) {    /*  ?????????? Not tested yet.   ansi-c*/
      nrows2 = nrows;
      errflag2 = errflag;
      dpotrf("L", &nrows2, W, &nrows2, &errflag2);   /*  Choleski.  L*L' = W where W will be replaced by lower triangular L.   ansi-c*/
      errflag = errflag2;
      if (errflag)  fn_DisplayError(".../mathlib.c/Aldivb_spd():  Error when calling Choleski dpotrf_().  Check if the L input matrix A_dm is positive definite or has legal values");
      if (x==b)  {
/*           //=== Solving for A*x=b.   ansi-c*/
         cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, nrows, W, nrows, x, 1);
         cblas_dtrsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, nrows, W, nrows, x, 1);
      }
      else {
         memcpy(x, b, nrows*sizeof(double));
         cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, nrows, W, nrows, x, 1);
         cblas_dtrsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, nrows, W, nrows, x, 1);
         x_dv->flag = V_DEF;
      }
      if ( (an!='N') && (an!='n') )  A_dm->flag = M_LT;
   }
   else  fn_DisplayError(".../mathlib.c/Aldivb_spd():  L input matrix A_dm must be symmetric");
/*     //dpotrf_((A_dm->flag & M_SU) ? "U" : "L", &nrows, A, &nrows, &errflag);  //Choleski. If "U", U'*U = A where A will be replaced by upper triangular U.   ansi-c*/
/*     // if (errflag<0) fn_DisplayError("Some element has an illegal value");   ansi-c*/
/*     // else if (errflag>0) fn_DisplayError("The leadding minor of some order, hence the entire matrix, is not positive definite");   ansi-c*/

   if ( (an=='N') || (an=='n') )  swzFree(W);
}
#else
/*  //No default routine yet.   ansi-c*/
#endif



#if defined (INTELCMATHLIBRARY)
double detspd(TSdmatrix *S_dm)
{
/*     //Determinant of symmetric positive definite (SPD) matrix must be positive.   ansi-c*/
/*     //We set the return value to be -1 if this matrix is NOT SPD.   ansi-c*/
   double valuedet;
   int _n, _i;
   int errflag;
/*     //===   ansi-c*/
   TSdmatrix *Work_dm = NULL;

   if (S_dm) {
      _n = S_dm->nrows;
      Work_dm = CreateMatrix_lf(_n, _n);
   }
   else  fn_DisplayError(".../mathlib.c/detspd():  Input matrix must be (1) created, (2) symmetric, (3) positive definite");

   if (S_dm->flag & M_SU)  errflag=chol(Work_dm, S_dm, 'U');
   else if (S_dm->flag & M_SL)  errflag=chol(Work_dm, S_dm, 'L');
   else  fn_DisplayError(".../mathlib.c/detpdf(): Input matrix S_dm must be either M_SU or M_SL");

   if (errflag) {
/*        //printf("\nFatal Error in .../mathlib.c/detspd() when calling chol() with error flag %d\n", errflag);   ansi-c*/
/*        //swzExit(EXIT_FAILURE);   ansi-c*/
      DestroyMatrix_lf(Work_dm);
      return (-1.0);
   }
   else
   {
      for (valuedet=1.0, _i=square(_n)-1; _i>=0;  _i -= _n+1)  valuedet *= Work_dm->M[_i];
/*                       //log(fabs(M[_i]));   ansi-c*/
      if (!isfinite(valuedet))  fn_DisplayError(".../mathlib.c/detspd():  the determinant is overflow.  Use logdetspd() instead");
/*                                            //Done with Work* arrays.   ansi-c*/
      DestroyMatrix_lf(Work_dm);
      return (square(valuedet));    /*  square() because Work_dm is a square root of S_dm.   ansi-c*/
   }
}
#else
/*  //No default routine yet.   ansi-c*/
#endif
/*  //---   ansi-c*/
#if defined (INTELCMATHLIBRARY)
double logdetspd(TSdmatrix *S_dm)
{
/*     //Determinant of symmetric positive definite (SPD) matrix must be positive.   ansi-c*/
/*     //We set the return value to be log(-1.0) (becomeing NaN) if this matrix is NOT SPD.   ansi-c*/
   double logvaluedet;
   int _n, _i;
   int errflag;
/*     //===   ansi-c*/
   TSdmatrix *Work_dm = NULL;

   if (S_dm) {
      _n = S_dm->nrows;
      Work_dm = CreateMatrix_lf(_n, _n);
   }
   else  fn_DisplayError(".../mathlib.c/detspd():  Input matrix must be (1) created, (2) symmetric, (3) positive definite");

   if (S_dm->flag & M_SU)  errflag=chol(Work_dm, S_dm, 'U');
   else if (S_dm->flag & M_SL)  errflag=chol(Work_dm, S_dm, 'L');
   else  fn_DisplayError(".../mathlib.c/logdetspd(): Input matrix S_dm must be either M_SU or M_SL");

   if (errflag) {
/*        //printf("\nFatal Error in .../mathlib.c/logdetspd() when calling chol() with error flag %d\n", errflag);   ansi-c*/
/*        //swzExit(EXIT_FAILURE);   ansi-c*/
      DestroyMatrix_lf(Work_dm);
      printf("----- errflag for chol() when calling logdetspd() in mathlib.c = %d -----", errflag);
      return (log(-1.0));
   }
   else
   {
      for (logvaluedet=0.0, _i=square(_n)-1; _i>=0;  _i -= _n+1)  logvaluedet += log(Work_dm->M[_i]);
/*                                            //Done with Work* arrays.   ansi-c*/
      DestroyMatrix_lf(Work_dm);
      return (2.0*logvaluedet);    /*  2.0* because Work_dm is a square root of S_dm.   ansi-c*/
   }
}
#else
/*  //No default routine yet.   ansi-c*/
#endif



#if defined (INTELCMATHLIBRARY)
double logdeterminant(TSdmatrix *A_dm) {
/*     //Outputs: log|A|.   ansi-c*/
/*     //------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n real general matrix.   ansi-c*/

   double retval;
   int _m, _n,
       errflag = -2;
   TSdmatrix *U_dm;

   if (!A_dm)  fn_DisplayError(".../logdeterminant(): Input matrix must be memory allocated (and make sure it has legal values with the flag M_GE)");
/*                   //NOTE: all properties of A_dm will be double checked again by lurgen() below.   ansi-c*/

/*     //=== Allocates memory used only in this function.   ansi-c*/
   U_dm = CreateMatrix_lf(_m=A_dm->nrows, _n=A_dm->ncols);

   errflag = lurgen(U_dm, (TSivector *)NULL, A_dm);   /*  Gets only the upper part of U_dm.   ansi-c*/
   if (errflag) fn_DisplayError(".../logdeterminant(): Error occurs when calling lurgen()");
   retval = tracelogfabs(U_dm);   /*  tracelogfabs(U) = trace(log(diag(U))).   ansi-c*/

/*     //=== Frees memory allocated only for this function.   ansi-c*/
   DestroyMatrix_lf(U_dm);

   return ( retval );
}
#else
/*  //No default routine yet.   ansi-c*/
#endif


int eigrsym_decomp(double *eval_v, double *evec_m, const double *s_m, const int _n, const char ul)  {   /*  , const char revec_yn, const char levec_yn) {   ansi-c*/
/*     // Outputs (dependent on Intel MKL):   ansi-c*/
/*     //   eval_v:  _n-by-1 eigenvalues in ascending order;   ansi-c*/
/*     //   evec_m:  _n-by-_n eigenvalues -- if (evec_m==NULL), no eigenvectors are computed; otherwise, x_m = evec_m*diag(eval_v)*inv(evec_m).   ansi-c*/
/*     //   errflag: error flag.   ansi-c*/
/*     //------------   ansi-c*/
/*     // Inputs:   ansi-c*/
/*     //   s_m:  _n-by_n real symmetric matrix.   ansi-c*/
/*     //   ul: if =='u' or 'U', s_m is upper triangular; if =='l' or 'L', s_m is lower triangular.   ansi-c*/
/*     //   ansi-c*/
/*     // Eigenanalysis of real symmetric square matrix with all eigenvalues and, optionally, eigenvectors.   ansi-c*/
/*     //   Experts' opinion: do NOT use Cuppen's divide-and-conquer algorithm; instead, use QR algorithm, which I guess this algorithm uses.   ansi-c*/

   int n1=_n, errflag=2,   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dsyev_ export a correct flag.   ansi-c*/
       lwork=_n*BLOCKSIZE_FOR_INTEL_MKL;
   double *tmpd0_m=tzMalloc(square(_n), double),
          *work_p=tzMalloc(lwork, double);

  lapack_int n2, lwork2, errflag2;

/*     //---------------------------   ansi-c*/
/*     // Obtains eigenvalues and, optionally, eigenvectors.   ansi-c*/
/*     //---------------------------   ansi-c*/
   memcpy(tmpd0_m, s_m, square(_n)*sizeof(double));

   n2 = n1;
   lwork2 = lwork;
   errflag2 = errflag;
   dsyev( (evec_m) ? "V" : "N", ((ul=='u') || (ul=='U')) ? "U" : "L", &n2, tmpd0_m, &n2, eval_v, work_p, &lwork2, &errflag2);
   errflag = errflag2;
   if (evec_m) memcpy(evec_m, tmpd0_m, square(_n)*sizeof(double));


/*     //---------------------------   ansi-c*/
/*     // Frees the allocated memory.   ansi-c*/
/*     //---------------------------   ansi-c*/
   if (work_p[0]>lwork) printf("Warning for /mathlib.c/eigrsym_decomp(): needs at least %d workspace for good performance "
                                 "but lwork is allocated with only %d space!\n", (int)work_p[0], lwork);
   tzDestroy(tmpd0_m);
   tzDestroy(work_p);

/*     //if (errflag<0) fn_DisplayError("/Subfolder: Calling eigrsym_decomp -- some element in input matrix has an illegal value");   ansi-c*/
/*     //else if (errflag>0) fn_DisplayError("/Subfolder: Calling eigrsym_decomp -- the factor U is exactly singular, so the solution cannot be computed");   ansi-c*/
   return (errflag);
}

int eigrgen_decomp(double *evalr_v, double *evali_v, double *revecr_m, double *reveci_m,  double *levecr_m, double *leveci_m, const double *x_m, const int _n)  {   /*  , const char revec_yn, const char levec_yn) {   ansi-c*/
/*     // Outputs (dependent on Intel MKL):   ansi-c*/
/*     //   evalr_v:  _n-by-1 real parts of eigenvalues;   ansi-c*/
/*     //   evali_v:  _n-by-1 imaginary parts of eigenvalues;   ansi-c*/
/*     //   revecr_m: if (revecr_m==NULL), no right eigenvectors are computed; otherwise, _n-by-_n corresponding *real* parts of right eigenvectors column by column: A*v(j)=lambda(j)*v(j);   ansi-c*/
/*     //   reveci_m: if (revecr_m!=NULL) -- must be initialized to zero, _n-by-_n *imaginary* parts of right eigenvectors corresponding to revecr_m;   ansi-c*/
/*     //   levecr_m: if (levecr_m==NULL), no left eigenvectors are computed; otherwise, n-by-n corresponding *real* parts of left eigenvectors column by column: u(j)^H*A=lambda(j)*u(j)^H, where H means conjugate transpose;   ansi-c*/
/*     //   leveci_m: if (levecr_m!=NULL) -- must be initialized to zero, _n-by-_n *imaginary* parts of left eigenvectors corresponding to revecr_m;   ansi-c*/
/*     //   errflag: error flag.   ansi-c*/
/*     //------------   ansi-c*/
/*     // Inputs:   ansi-c*/
/*     //   x_m:  _n-by_n real general (non-symmetric) matrix.   ansi-c*/
/*     //   ansi-c*/
/*     // Eigenanalysis of real general (non-symmetric) square matrix with all eigenvalues and, optionally, eigenvectors.   ansi-c*/

   int n1=_n, errflag=2,   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dgeev_ export a correct flag.   ansi-c*/
       lwork=_n*BLOCKSIZE_FOR_INTEL_MKL,
       _i, _j;
   double *tmpd0_m=tzMalloc(square(_n), double),     /*  @@Must be freed in this function.@@   ansi-c*/
          *work_p=tzMalloc(lwork, double);           /*  @@Must be freed in this function.@@   ansi-c*/

   lapack_int n2, lwork2, errflag2;

/*     //---------------------------   ansi-c*/
/*     // Starts with x_m -- the matrix to be decomposed.   ansi-c*/
/*     //---------------------------   ansi-c*/
   memcpy(tmpd0_m, x_m, square(_n)*sizeof(double));

/*     //---------------------------   ansi-c*/
/*     // Obtains eigenvalues and, optionally, eigenvectors.   ansi-c*/
/*     //---------------------------   ansi-c*/
   n2 = n1;
   lwork2 = lwork;
   errflag2 = errflag;
   dgeev( (levecr_m) ? "V" : "N", (revecr_m) ? "V" : "N", &n2, tmpd0_m, &n2, evalr_v, evali_v,
                                                          levecr_m, &n2, revecr_m, &n2, work_p, &lwork2, &errflag2);
   errflag = errflag2;

/*     //---------------------------   ansi-c*/
/*     // Frees the allocated memory.   ansi-c*/
/*     //---------------------------   ansi-c*/
   if (work_p[0]>lwork) printf("Warning for /mathlib.c/eigrgen_decomp(): needs at least %d workspace for good performance "
                                 "but lwork is allocated with only %d space!\n", (int)work_p[0], lwork);
   if (work_p) swzFree(work_p);
   if (tmpd0_m) swzFree(tmpd0_m);

/*     //---------------------------   ansi-c*/
/*     // Checks error conditions.   ansi-c*/
/*     // Exports final results.   ansi-c*/
/*     //---------------------------   ansi-c*/
/*     //if (errflag<0) fn_DisplayError("/Subfolder: Calling eigrgen_decomp -- some element in input matrix has an illegal value");   ansi-c*/
/*     //else if (errflag>0) fn_DisplayError("/Subfolder: Calling eigrgen_decomp -- the QR algorithm failed to compute all the eigenvalues and no eigenvectors have been computed");   ansi-c*/
   if (errflag) return (errflag);
   else {
      if (revecr_m) {                    /*   Tested against Matlab output.  Works!  10/13/02.   ansi-c*/
         for (_j=0; _j<_n-1; _j++)
            if (evali_v[_j] && (evali_v[_j] == -evali_v[_j+1]))
               for (_i=0; _i<_n; _i++) {
                  reveci_m[_i+(_j+1)*_n] = -(reveci_m[_i+_j*_n]=revecr_m[_i+(_j+1)*_n]);
                  revecr_m[_i+(_j+1)*_n] = revecr_m[_i+_j*_n];
               }
      }
      if (levecr_m) {                   /*  Hasn't tested against any other established program, but it seems working.  10/13/02.   ansi-c*/
         for (_j=0; _j<_n-1; _j++)
            if (evali_v[_j] && (evali_v[_j] == -evali_v[_j+1]))
               for (_i=0; _i<_n; _i++) {
                  leveci_m[_i+(_j+1)*_n] = -(leveci_m[_i+_j*_n]=levecr_m[_i+(_j+1)*_n]);
                  levecr_m[_i+(_j+1)*_n] = levecr_m[_i+_j*_n];
               }
      }
      return (errflag);
   }
}



int chol_decomp(double *D, const double *s_m, const int _n, const char ul) {
/*     //Outputs:   ansi-c*/
/*     //  D: _n-by_n -- if ul='u' or 'U',  D'*D = s_m where D is stored only at the upper triangular part;   ansi-c*/
/*     //                if ul='l' or 'L',  D*D' = s_m where D is stored only at the lower triangular part.   ansi-c*/
/*     //  errflag: error flag (=0 means successful).   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  s_m:  _n-by-_n symmetric, positive definite matrix (whose only triangular part is used by dpotrf_).   ansi-c*/
/*     //  ul: if =='u' or 'U', D (as well as s_m) is upper triangular; if =='l' or 'L', D (as well as s_m) is lower triangular.   ansi-c*/
/*     //   ansi-c*/
/*     // Choleski decomposition of s_m.   ansi-c*/
/*     //   For the MATLAB libriary, ul doest not apply and chol_decomp always takes the 'U' form.   ansi-c*/
/*     //     And Matlab 6.5 (R13) has a different number of inputs for mlfChol().   ansi-c*/
/*     //     See ...\extern\include\libmatlbm.h (included by matlab.h) for the definition of mlfChol().   ansi-c*/
/*     //     So R12 version must be used if one chooses to use the MATLAB libriary.   ansi-c*/

   #ifdef INTELCMATHLIBRARY    /*  Intel MKL Lapack dependent code.   ansi-c*/
      int errflag=2, _m=_n, _i, _j, tmpi0;   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dpotrf_ export a correct flag.   ansi-c*/

      lapack_int _m2, errflag2;

/*        //=== Fills the triangular part that is used for Choleski decomposition.   ansi-c*/

      switch (ul) {
         case 'u': case 'U':
            for (_j=0; _j<_n; _j++) {
               for (_i=0; _i<=_j; _i++) {
                  tmpi0 = _i+_j*_n;
                  D[tmpi0] = s_m[tmpi0];
               }
               for (; _i<_n; _i++) D[_i+_j*_n] = 0.0;    /*  Initializes the other part of D to zero so as to make it visible and readable.   ansi-c*/
            }
            break;
         case 'l': case 'L':
            for (_j=0; _j<_n; _j++) {
               for (_i=0; _i<_j; _i++) D[_i+_j*_n] = 0.0;    /*  Initializes the other part of D to zero so as to make it visible and readable.   ansi-c*/
               for (; _i<_n; _i++) {
                  tmpi0 = _i+_j*_n;
                  D[tmpi0] = s_m[tmpi0];
               }
            }
         default:
            return (-1);
      }
/*        //=== Choleski decomposition.   ansi-c*/
      _m2 = _m;
      errflag2 = errflag;
      dpotrf(((ul=='u') || (ul=='U')) ? "U" : "L", &_m2, D, &_m2, &errflag2);
      errflag = errflag2;
/*        //---   ansi-c*/
/*        // if (errflag<0) fn_DisplayError("Some element has an illegal value");   ansi-c*/
/*        // else if (errflag>0) fn_DisplayError("The leadding minor of some order, hence the entire matrix, is not positive definite");   ansi-c*/
      return (errflag);
   #endif
   #ifdef MATLABCMATHLIBRARY   /*  Matlab dependent code.   ansi-c*/
      mxArray *ms_m=mlfDoubleMatrix(_n, _n, s_m, NULL),      /*  @@Must be freed in this function.@@ mx version of s_m.   ansi-c*/
            *mxD=NULL,                                     /*  @@Must be freed in this function.@@ mx version of D.   ansi-c*/
            *mxflag;
      int errflag;

      mxD = mlfChol(&mxflag, ms_m);
      errflag = (int)mxGetScalar(mxflag);
/*        //if (errflag) fn_DisplayError("Function mathlib.c\\chol_decomp():  matrix must be positive definite for choleski decomposition");   ansi-c*/
      memcpy(D, mxGetPr(mxD), square(_n)*sizeof(double));

/*        //=== Frees up allocated mxArray.   ansi-c*/
      mxDestroyArray(mxD);
      mxDestroyArray(ms_m);
      mxDestroyArray(mxflag);

      return errflag;
   #endif
}

int inv_spd(double *D, const double *s_m, const int _n, const char ul) {
/*     // Inverse of symmetric, positive-definite matrix s_m.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  D: _n-by_n inverse of s_m.   ansi-c*/
/*     //  errflag: error flag (=0 means successful).   ansi-c*/
/*     //--------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  s_m:  _n-by-_n symmetric, positive definite matrix (whose only triangular part is used by dpotrf_).   ansi-c*/
/*     //  ul: if =='u' or 'U', D (as well as s_m) is upper triangular; if =='l' or 'L', D (as well as s_m) is lower triangular.   ansi-c*/

   int errflag=2, _m=_n, _i, _j, tmpi0;   /*  errflat=0 implies successful decomposition.  But we start with 2 so as to let dpotrf_ export a correct flag.   ansi-c*/

   lapack_int _m2, errflag2;

/*     //=== Fills the triangular part that is used for Choleski decomposition.   ansi-c*/
   switch (ul) {
      case 'u': case 'U':
         for (_j=0; _j<_n; _j++) {
            for (_i=0; _i<=_j; _i++) {
               tmpi0 = _i+_j*_n;
               D[tmpi0] = s_m[tmpi0];
            }
            for (; _i<_n; _i++) D[_i+_j*_n] = 0.0;    /*  Initializes the other part of D to zero so as to make it visible and readable.   ansi-c*/
         }
         break;
      case 'l': case 'L':
         for (_j=0; _j<_n; _j++) {
            for (_i=0; _i<_j; _i++) D[_i+_j*_n] = 0.0;    /*  Initializes the other part of D to zero so as to make it visible and readable.   ansi-c*/
            for (; _i<_n; _i++) {
               tmpi0 = _i+_j*_n;
               D[tmpi0] = s_m[tmpi0];
            }
         }
         break;
      default:
         return (-1);
   }
/*     //=== Choleski decomposition.   ansi-c*/
   _m2 = _m;
   errflag2 = errflag;
   dpotrf(((ul=='u') || (ul=='U')) ? "U" : "L", &_m2, D, &_m2, &errflag2);
   errflag = errflag2;
   if (errflag) return (errflag);
/*     //=== Takes inverse.   ansi-c*/
   dpotri(((ul=='u') || (ul=='U')) ? "U" : "L", &_m2, D, &_m2, &errflag2);
   errflag = errflag2;
   return (errflag);
/*     //---   ansi-c*/
/*     // if (errflag<0) fn_DisplayError("Some element has an illegal value");   ansi-c*/
/*     // else if (errflag>0) fn_DisplayError("Not symmetric positive definite or matrix inversion cannot be computed");   ansi-c*/
}




/*  //=======================================================   ansi-c*/
/*  // BLAS routines -- all based on Intel MKL (or IMSL C Math library).   ansi-c*/
/*  //=======================================================   ansi-c*/
/*  //void ScalingVectorUpdate(const double _alpha, TSdvector *x_dv) {   ansi-c*/
/*  //   //Output: x = alpha*x;   ansi-c*/
/*  //   //   ansi-c*/
/*  //   call dscal (n, da, DX, incx)   ansi-c*/
/*  //}   ansi-c*/
/*  //Use the scaling vector MKL routine.   ansi-c*/

double VectorDotVector(TSdvector *x1_dv, TSdvector *x2_dv) {
/*     //Output: Return sum(x1[i] * x2[i]) over i=1, ..., n.   ansi-c*/
/*     //  Allows the case x1_dv = x2_dv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x1_dv: _n-by-1 double vector.   ansi-c*/
/*     //  x2_dv: _n-by-1 double vector.   ansi-c*/
   int _n, _i;
   double *x1, *x2,
          sum2 = 0.0;   /*  Cumulative: must be set to 0.0.   ansi-c*/

   if ( !x1_dv || !x2_dv ) fn_DisplayError(".../mathlib.c/VectorDotVector(): All input vectors must be created (memory-allocated)");

   if ( (x1=x1_dv->v) == (x2=x2_dv->v) ) {
      if ( !x1_dv->flag )  fn_DisplayError(".../mathlib.c/VectorDotVector(): Input vectors must be given legal values");
      for (_i=x1_dv->n-1; _i>=0; _i--)  sum2 += square(x1[_i]);
   }
   else {
      if ( (_n=x1_dv->n) != x2_dv->n )  fn_DisplayError(".../mathlib.c/VectorDotVector(): Dimensions of the two input vectors must be same");
      else if ( !x1_dv->flag || !x2_dv->flag )  fn_DisplayError(".../mathlib.c/VectorDotVector(): Both input vectors must be given legal values");
      for (_i=_n-1; _i>=0; _i--)  sum2 += x1[_i]*x2[_i];
   }

   return ( sum2 );
/*     //return cblas_ddot(_n, x1_dv->v, 1, x2_dv->v, 1);   ansi-c*/
}


void ScalarTimesVectorUpdate(TSdvector *x2_dv, const double _alpha, TSdvector *x1_dv) {
/*     //Output:  x2 = alpha * x1 + x2;   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  alpha:  a double scalar;   ansi-c*/
/*     //  x1:  n-by-1 double vector.   ansi-c*/
   int _n;

   if ( !x1_dv || !x2_dv )  fn_DisplayError(".../mathlib.c/ScalarTimesVectorUpdate():  All input vectors must be created (memory-allocated)");
   else _n = x1_dv->n;

   if (_n != x2_dv->n)   fn_DisplayError(".../mathlib.c/ScalarTimesVectorUpdate():  All input vectors must have the same length");
   else if ( !x1_dv->flag )  fn_DisplayError(".../mathlib.c/ScalarTimesVectorUpdate():  R input vector must be given values");
   else {
      if ( x1_dv->v == x2_dv->v )  fn_DisplayError(".../mathlib.c/ScalarTimesVectorUpdate():  Two input vectors cannot be the same.  Instead, use SclarTimesVector() for this option if you are sure this is what you want");
      cblas_daxpy(_n, _alpha, x1_dv->v, 1, x2_dv->v, 1);
      x2_dv->flag = V_DEF;
   }
}

void ScalarTimesVector(TSdvector *x_dv, const double _alpha, TSdvector *a_dv, const double _beta) {
/*     //Output: x_dv = alpha*a_dv + beta*x_dv where x_dv is n-by-1.   ansi-c*/
/*     //  When beta=0.0 and x_dv->v = a_dv->v, x_dv->v will be replaced by new values.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: n-by-1.   ansi-c*/
/*     //  _alpha: a double scalar.   ansi-c*/
/*     //  _beta: a double scalar.   ansi-c*/
   int _i, _n;
   double *x, *a;

   if ( !x_dv || !a_dv )  fn_DisplayError(".../mathlib.c/ScalarTimesVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag )  fn_DisplayError(".../mathlib.c/ScalarTimesVector(): R input vector must be given values");
   else {
      _n = x_dv->n;
      x = x_dv->v;
      a = a_dv->v;
   }

   if ( _n != a_dv->n )  fn_DisplayError(".../mathlib.c/ScalarTimesVector(): Two input vectors must have the same length");
   if (_beta == 0.0) {
      #if defined (INTELCMATHLIBRARY)             /*   define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
         if ( x == a )  cblas_dscal(_n, _alpha, x, 1);
         else {
            memcpy(x_dv->v, a_dv->v, _n*sizeof(double));
            x_dv->flag = V_DEF;
            cblas_dscal(_n, _alpha, x, 1);
         }
      #else   /*  SWITCHTOTZCMATH: use my own C math library (which is faster than MKL sometimes); undef: use others.   ansi-c*/
         for (_i=_n-1; _i>=0; _i--)  x[_i] = _alpha*a[_i];
         x_dv->flag = V_DEF;
      #endif
   }
   else if (_beta == 1.0) {
      if ( x == a )  fn_DisplayError(".../mathlib.c/ScalarTimesVector(): Two input vectors must be different, i.e., pointing to different memory places");
      cblas_daxpy(_n, _alpha, a, 1, x, 1);
      x_dv->flag = V_DEF;
   }
   else {
      if ( x == a )  fn_DisplayError(".../mathlib.c/ScalarTimesVector(): Two input vectors must be different, i.e., pointing to different memory places");
      for (_i=_n-1; _i>=0; _i--)  x[_i] = _alpha*a[_i] + _beta*x[_i];
      x_dv->flag = V_DEF;
   }
}

void VectorPlusMinusVectorUpdate(TSdvector *x_dv, const TSdvector *b_dv, double _alpha)
{
/*     //Output: x_dv =_alpha * b_dv +  x_dv where x_dv is _n-by-1.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
   int _n;


   if ( !x_dv || !b_dv ) fn_DisplayError(".../mathlib.c/VectorPlusMinusVectorUpdate(): All input vectors must be created (memory-allocated)");
   else if ( !b_dv->flag || !x_dv->flag )  fn_DisplayError(".../mathlib.c/VectorPlusMinusVectorUpdate(): All input vectors must be given values");
   else {
      _n = x_dv->n;
   }
   if ( _n != b_dv->n ) fn_DisplayError(".../mathlib.c/VectorPlusMinusVectorUpdate(): Dimensions of all input vectors must be same");

   cblas_daxpy(_n, _alpha, b_dv->v, 1, x_dv->v, 1);
}

void VectorPlusMinusVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv, double _alpha)
{
/*     //???????? Use tz_VectorPlusMinusVector() or VectorPlusVector() or VectorMinusVector().   ansi-c*/
/*     //????????? NOT finished yet.   ansi-c*/
/*     //????????Must add _beta for x_dv = alpha*a_dv + beta*b_dv. If x_dv=b_dv, update.   ansi-c*/
/*     //??????????? NOT fully tested yet.   ansi-c*/
/*     //Output: x_dv = a_dv + _alpha * b_dv where x_dv is _n-by-1.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: _n-by-1 double vector.   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
   int _n;


   if ( !x_dv || !a_dv || !b_dv ) fn_DisplayError(".../mathlib.c/VectorPlusMinusVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag )  fn_DisplayError(".../mathlib.c/VectorPlusMinusVector(): All input vectors must be given values");
   else {
      _n = x_dv->n;
   }
   if ( (_n != a_dv->n) || (_n != b_dv->n) ) fn_DisplayError(".../mathlib.c/VectorPlusMinusVector(): Dimensions of all input vectors must be same");

   memcpy(x_dv->v, a_dv->v, _n*sizeof(double));
   cblas_daxpy(_n, _alpha, b_dv->v, 1, x_dv->v, 1);
}

void VectorTimesSelf(TSdmatrix *C_dm, const TSdvector *a_dv, const double _alpha, const double _beta, const char ul)
{
/*     //Computes C = alpah*a*a' + beta*C where   ansi-c*/
/*     // Output:   ansi-c*/
/*     //   the symmetric matrix C.   ansi-c*/
/*     // Inputs:   ansi-c*/
/*     //   a is m-by-1,   ansi-c*/
/*     //   C is m-by-m symmetric matrix,   ansi-c*/
/*     //   alpha: a double scalar,   ansi-c*/
/*     //   beta: a double scalar.   ansi-c*/
/*     //   ul: if=='U' or 'u', only the upper triangular part of C is to be referenced; otherwise, only the lower triangular part of C is to be referenced;   ansi-c*/
   int _m, _n;
   int _i, _j;
   double *v, *M;

   if ( !C_dm || !a_dv ) fn_DisplayError(".../mathlib.c/VectorTimesSelf(): At least one of the pointer arguments is not created (memory-allocated)");
   else if (!a_dv->flag) fn_DisplayError(".../mathlib.c/VectorTimesSelf():  Input vector must have legal values");
   else {
      _m = C_dm->nrows;
      _n = C_dm->ncols;
   }

   if ( (_m != a_dv->n) || (_m !=_n) ) fn_DisplayError(".../mathlib.c/VectorTimesSelf(): (1) Size of the input matrix and dimensions of the two input vectors do not match.  (2) Output matrix must be square");
   else {
      if ( _beta == 1.0 ) {
         #if defined (INTELCMATHLIBRARY)             /*   define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
/*           //$$$$ cblas_dsyrk is much slower than the following line.  cblas_dsyrk(CblasColMajor, ((ul=='u') || (ul=='U')) ? CblasUpper : CblasLower, CblasNoTrans, _m, 1, _alpha, a_dv->v, _m, _beta, C_dm->M, _m);   ansi-c*/
         cblas_dsyr(CblasColMajor, ((ul=='U') || (ul=='u')) ? CblasUpper : CblasLower, _m, _alpha, a_dv->v, 1, C_dm->M, _m);
         C_dm->flag = ((ul=='U') || (ul=='u')) ? M_SU : M_SL;
         #else  /*  Corresponds to the default: SWITCHTOTZCMATH -- use my own C math library, which is faster than cblas_dsyrk().   ansi-c*/
         v = a_dv->v;
         M = C_dm->M;
         if ( (ul == 'U') || (ul == 'u') )  {
            C_dm->flag = M_SU;
            if ( _alpha==1.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] += v[_i] * v[_j];
                  }
               }
            }
            else {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] += _alpha * v[_i] * v[_j];
                  }
               }
            }
         }
         else {
            C_dm->flag = M_SL;
            if ( _alpha==1.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] += v[_i] * v[_j];
                  }
               }
            }
            else {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] += _alpha * v[_i] * v[_j];
                  }
               }
            }
         }
         #endif
      }
      else {
/*           //Corresponds to the default: SWITCHTOTZCMATH -- use my own C math library (which is faster than MKL sometimes).   ansi-c*/
         v = a_dv->v;
         M = C_dm->M;
         if ( (ul == 'U') || (ul == 'u') )  {
            C_dm->flag = M_SU;
            if ( _alpha==1.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] = v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                  }
               }
            }
            else {
               if ( _beta==0.0 ) {
                  for ( _j=0; _j<_m; _j++ ) {
                     for ( _i=0; _i<=_j; _i++ ) {
                        M[mos(_i, _j, _m)] = _alpha * v[_i] * v[_j];
                     }
                  }
               }
               else {
                  for ( _j=0; _j<_m; _j++ ) {
                     for ( _i=0; _i<=_j; _i++ ) {
                        M[mos(_i, _j, _m)] = _alpha* v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                     }
                  }
               }
            }
         }
         else {
            C_dm->flag = M_SL;
            if ( _alpha==1.0 ) {
               if ( _beta==0.0 ) {
                  for ( _j=0; _j<_m; _j++ ) {
                     for ( _i=_j; _i<_m; _i++ ) {
                        M[mos(_i, _j, _m)] = v[_i] * v[_j];
                     }
                  }
               }
               else {
                  for ( _j=0; _j<_m; _j++ ) {
                     for ( _i=_j; _i<_m; _i++ ) {
                        M[mos(_i, _j, _m)] = v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                     }
                  }
               }
            }
            else {
               if ( _beta==0.0 ) {
                  for ( _j=0; _j<_m; _j++ ) {
                     for ( _i=_j; _i<_m; _i++ ) {
                        M[mos(_i, _j, _m)] = _alpha * v[_i] * v[_j];
                     }
                  }
               }
               else {
                  for ( _j=0; _j<_m; _j++ ) {
                     for ( _i=_j; _i<_m; _i++ ) {
                        M[mos(_i, _j, _m)] = _alpha* v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                     }
                  }
               }
            }
         }
      }
   }
}

#if defined (INTELCMATHLIBRARY)
void VectorTimesVector(TSdmatrix *C_dm, const TSdvector *a_dv, const TSdvector *b_dv, const double _alpha, const double _beta) {
/*     //?????? NOT tested for _beta != 1.0.   ansi-c*/
/*     //Output is the matrix C and all other arguments are inputs.   ansi-c*/
/*     //If beta != 0, always converting C (if symmetric or trianuglar) to a general matrix before the operation.   ansi-c*/
/*     //The fastest way is to let _beta = 1.0.   ansi-c*/
/*     //Computes C = alpah*a*b' + beta*C where   ansi-c*/
/*     //  a is m-by-1,   ansi-c*/
/*     //  b is n-by-1,   ansi-c*/
/*     //  C is m-by-n general matrix,   ansi-c*/
/*     //  alpha: a double scalar,   ansi-c*/
/*     //  beta: a double scalar.   ansi-c*/
   int _m, _n;

   if ( !C_dm || !a_dv || !b_dv ) fn_DisplayError(".../mathlib.c/VectorTimesVector(): At least one of the pointer arguments is not created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag ) fn_DisplayError(".../mathlib.c/VectorTimesVector():  Both input R vectors must be given values");
   else {
      _m = C_dm->nrows;
      _n = C_dm->ncols;
   }
   if (_beta != 0.0) {
      if ( !(C_dm->flag & M_GE) && (_m = _n) ) {
         if (C_dm->flag & M_SU)   SUtoGE(C_dm);
         else if (C_dm->flag & M_SL)   SLtoGE(C_dm);
         else  fn_DisplayError(".../mathlib.c/VectorTimesVector(): (a) make sure C_dm has legal values; (b) for M_UT and M_LT, I have not got time to convert it to a general matrix");
      }
   }



   if ( (_m != a_dv->n) || (_n != b_dv->n) ) fn_DisplayError(".../mathlib.c/VectorTimesVector(): Size of the input matrix and dimensions of the two input vectors do not match");
   else {
      if (_beta==1.0)  {
         cblas_dger(CblasColMajor, _m, _n, _alpha, a_dv->v, 1, b_dv->v, 1, C_dm->M, _m);
         C_dm->flag = M_GE;
      }
      else {
         cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, _m, _n, 1, _alpha, a_dv->v, _m, b_dv->v, 1, _beta, C_dm->M, _m);
/*           //???????????The above is probably too slow.  Try the following two lines instead.  3/10/03.   ansi-c*/
/*           //????? Test to make sure this works for beta=0.   ansi-c*/
/*           //cblas_dscal(_m*_n, _beta, C_dm->M, 1);   ansi-c*/
/*           //cblas_dger(CblasColMajor, _m, _n, _alpha, a_dv->v, 1, b_dv->v, 1, C_dm->M, _m);   ansi-c*/
         C_dm->flag = M_GE;
      }
   }
}
#else
/*  //No default routine yet.   ansi-c*/
#endif


void MatrixPlusMinusMatrixUpdate(TSdmatrix *X_dm, TSdmatrix *A_dm, double _alpha)
{
/*     //$$$$$ If A_dm or X_dm is only upper or lower symmetric, it will be always converted to a general (and symmetric) matrix.  $$$$$$   ansi-c*/
/*     //Output: X =_alpha * A + X where X_dm is an m-by-n general (and possibly symmetric) matrix.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n general or symmetric matrix.   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
   int _m, _n, nels;


   if ( !X_dm || !A_dm )  fn_DisplayError(".../mathlib.c/MatrixPlusMinusMatrixUpdate(): All input matrices must be created (memory-allocated)");
   else if ( !X_dm->flag || !A_dm->flag )   fn_DisplayError(".../mathlib.c/MatrixPlusMinusMatrixUpdate(): Both input matrices must be given values");
   else {
      _m = X_dm->nrows;
      _n = X_dm->ncols;
      nels = _m * _n;
   }

   if ( (_m != A_dm->nrows) || (_n != A_dm->ncols) )  fn_DisplayError(".../mathlib.c/MatrixPlusMinusMatrixUpdate(): Dimensions of all input matrices must be same");

/*     //=== Making both X_dm and A_dm general if not yet.   ansi-c*/
   if ( !(X_dm->flag & M_GE) ) {
      if (X_dm->flag & M_SU)   SUtoGE(X_dm);
      else if (X_dm->flag & M_SL)   SLtoGE(X_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixPlusMinusMatrixUpdate(): Haven't got time to deal with the M_UT and M_LT cases for X_dm");
   }
   if ( !(A_dm->flag & M_GE) ) {
      if (A_dm->flag & M_SU)   SUtoGE(A_dm);
      else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixPlusMinusMatrixUpdate(): Haven't got time to deal with the M_UT and M_LT cases for A_dm");
   }

   cblas_daxpy(nels, _alpha, A_dm->M, 1, X_dm->M, 1);   /*  This operation may be much cheaper than explicitly using SU or SL operations with two for loops and integer multiplications for matrix offsets.   ansi-c*/

   if ( X_dm->flag != A_dm->flag ) {
/*        //printf("WARNING for .../mathlib.c/MatrixPlusMinusMatrixUpdate(): the two input matrices do not have the same matrix type (or flag), so the output matrix is reset to M_GE");   ansi-c*/
      X_dm->flag = M_GE;   /*  Reset to a general matrix only; otherwise, keep the original X_dm->flag.   ansi-c*/
   }
/*     //if ( X_dm->flag != A_dm->flag )  fn_DisplayError(".../mathlib.c/MatrixPlusMinusMatrixUpdate(): both input matrices must have the same matrix type (or flag)");   ansi-c*/
}

void MatrixTimesVector(TSdvector *x_dv, TSdmatrix *A_dm, const TSdvector *b_dv, const double _alpha, const double _beta, const char tn)
{
/*     //?????  This is NOT checked yet: If x_dv = b_dv, x_dv or b_dv will be relaced by alpha*A*x + beta*x.   ansi-c*/
/*     //Output: x_dv = _alpha*A_dm'*b_dv + _beta*x_dv  for tn=='T'; x_dv = _alpha*A_dm*b_dv + _beta*x_dv  for tn=='N'   ansi-c*/
/*     //  where x_dv->v is ncols-by-1 or nrows-by-1 and needs not be initialized outside this function if _beta is set to 0.0.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm->M: nrows-by-ncols;   ansi-c*/
/*     //  b_dv->v: nrows-by-1 or ncols-by-1;   ansi-c*/
/*     //  _alpha: double scalar;   ansi-c*/
/*     //  _beta:  double scalar;   ansi-c*/
/*     //  tn: if =='T' or 't', transpose of A_dm is used; otherwise, A_dm itself (no transpose) is used.   ansi-c*/

   if ( !x_dv || !A_dm || !b_dv) fn_DisplayError(".../mathlib.c/MatrixTimesVector(): At least one of the pointer arguments is not created (memory-allocated)");
   else if ( !A_dm->flag || !b_dv->flag )  fn_DisplayError(".../mathlib.c/MatrixTimesVector(): R input matrix or vector must be given values");
   if ( !(A_dm->flag & M_GE) ) {
      if (A_dm->flag & M_SU)   SUtoGE(A_dm);
      else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixTimesVector(): For M_UT and M_LT, use TrimatrixTimesVector() instead");
   }


   if ((tn=='T' || tn=='t') && (A_dm->nrows==b_dv->n) && (A_dm->ncols==x_dv->n)) {
      cblas_dgemv(CblasColMajor, CblasTrans, A_dm->nrows, A_dm->ncols, _alpha, A_dm->M, A_dm->nrows, b_dv->v, 1, _beta, x_dv->v, 1);
      x_dv->flag = V_DEF;
   }
   else if ( (A_dm->ncols==b_dv->n) && (A_dm->nrows==x_dv->n) ) {
      cblas_dgemv(CblasColMajor, CblasNoTrans, A_dm->nrows, A_dm->ncols, _alpha, A_dm->M, A_dm->nrows, b_dv->v, 1, _beta, x_dv->v, 1);
      x_dv->flag = V_DEF;
   }
/*  //--- The following if clause is wrong because, when tn=='N', A_dm->ncols == b_dv->n, NOT A_dm->nraws==b_dv_>n.   ansi-c*/
/*  //   if ((A_dm->nrows==b_dv->n) && (A_dm->ncols==x_dv->n)) {   ansi-c*/
/*  //      cblas_dgemv(CblasColMajor, (tn=='T' || tn=='t') ? CblasTrans : CblasNoTrans, A_dm->nrows, A_dm->ncols, _alpha, A_dm->M, A_dm->nrows, b_dv->v, 1, _beta, x_dv->v, 1);   ansi-c*/
/*  //      x_dv->flag = V_DEF;   ansi-c*/
/*  //   }   ansi-c*/
   else  fn_DisplayError(".../mathlib.c/MatrixTimesVector(): Size of the input matrix and dimensions of the two input vectors do not match");
}


#if defined (INTELCMATHLIBRARY)
void TrimatrixTimesVector(TSdvector *x_dv, TSdmatrix *A_dm, TSdvector *b_dv, const char tn, const char un)
{
/*     //Output: x_dv = A_dm'*b_dv  for tn=='T'; x_dv = A_dm*b_dv  for tn=='N' where x_dv->v is _n-by-1.   ansi-c*/
/*     //  If x_dv = b_dv (which gives the fastest return, so try to use this option when possible), x_dv will be relaced by A*b or A'*b.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm->M: _n-by-_n triangular matrix.   ansi-c*/
/*     //  b_dv->v: _n-by-1 vector.   ansi-c*/
/*     //  tn: if =='T' or 't', transpose of A_dm is used; otherwise, A_dm itself (no transpose) is used.   ansi-c*/
/*     //  un: if =='U' or 'u', A_dm is unit triangular; otherwise, A_dm is non-unit triangular (i.e., a regular triangular matrix).   ansi-c*/

   int _n;
/*  //   double *x, *b;   ansi-c*/

   if ( !x_dv || !A_dm || !b_dv)  fn_DisplayError(".../mathlib.c/TrimatrixTimesVector(): At least one of the pointer arguments is not created (memory-allocated)");
   else if ( !A_dm->flag || !b_dv->flag )  fn_DisplayError(".../mathlib.c/TrimatrixTimesVector(): R input matrix or vector must be given values");
   else if ( ((_n = A_dm->nrows) != x_dv->n) )  fn_DisplayError(".../mathlib.c/TrimatrixTimesVector(): Size of input matrix and dimension of input vector must match");
   else if ( !(A_dm->flag & (M_UT | M_LT) ) )  fn_DisplayError(".../mathlib.c/TrimatrixTimesVector(): Make sure R matrix is triangular (i.e., M_UT or M_LT)");

/*  //   if ( (x = x_dv->v) == (b = b_dv->v) )  //Commented out on 22 Oct 03.   ansi-c*/
   if ( x_dv == b_dv )
      cblas_dtrmv(CblasColMajor, (A_dm->flag & M_UT) ? CblasUpper : CblasLower, (tn=='T' || tn=='t') ? CblasTrans : CblasNoTrans, (un=='U' || un=='u') ? CblasUnit : CblasNonUnit, A_dm->nrows, A_dm->M, A_dm->nrows, x_dv->v, 1);
   else {
      if ( _n != b_dv->n )  fn_DisplayError(".../mathlib.c/TrimatrixTimesVector(): Two vectors must have the same length");
      memcpy(x_dv->v, b_dv->v, _n*sizeof(double));
      cblas_dtrmv(CblasColMajor, (A_dm->flag & M_UT) ? CblasUpper : CblasLower, (tn=='T' || tn=='t') ? CblasTrans : CblasNoTrans, (un=='U' || un=='u') ? CblasUnit : CblasNonUnit, A_dm->nrows, A_dm->M, A_dm->nrows, x_dv->v, 1);
      x_dv->flag = V_DEF;
   }
}
#else
/*  //No default routine yet.   ansi-c*/
#endif


#if defined (INTELCMATHLIBRARY)
void SymmatrixTimesVector(TSdvector *x_dv, TSdmatrix *A_dm, TSdvector *b_dv, const double _alpha, const double _beta)
{
/*     //?????  This is NOT checked yet: If x_dv = b_dv, x_dv or b_dv will be relaced by alpha*A*x + beta*x.   ansi-c*/
/*     //Output:   ansi-c*/
/*     //  x_dv = alpha*A_dm*b_dv + beta*x_dv  where x_dv->v is _n-by-1.   ansi-c*/
/*     //    When beta=0, there is no need to initialize the value of x_dv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm->M: _n-by-_n triangular matrix.   ansi-c*/
/*     //  b_dv->v: _n-by-1 vector.   ansi-c*/
/*     //  _alpha: double scalar;   ansi-c*/
/*     //  _beta:  double scalar;   ansi-c*/

   int _n;

   if ( !x_dv || !A_dm || !b_dv)  fn_DisplayError(".../mathlib.c/SymmatrixTimesVector(): all input and output arguments must be created (memory-allocated)");
   else if ( !A_dm->flag || !b_dv->flag )  fn_DisplayError(".../mathlib.c/SymmatrixTimesVector(): R input matrix or vector must be given values");
   else if ( ((_n = A_dm->nrows) != x_dv->n) || (_n != b_dv->n) )  fn_DisplayError(".../mathlib.c/SymmatrixTimesVector(): Size of input matrix and dimensions of input and output vectors must all match");
   else if ( !(A_dm->flag & (M_SU | M_SL) ) )  fn_DisplayError(".../mathlib.c/SymmatrixTimesVector(): Make sure R input matrix is symmetric (i.e., M_SU or M_SL)");

   cblas_dsymv(CblasColMajor, (A_dm->flag & M_SU) ? CblasUpper : CblasLower, _n, _alpha, A_dm->M, _n, b_dv->v, 1, _beta, x_dv->v, 1);
   x_dv->flag = V_DEF;
}
#else
/*  //No default routine yet.   ansi-c*/
#endif



void VectorTimesMatrix(TSdvector *x_dv, const TSdvector *b_dv, TSdmatrix *A_dm, const double _alpha, const double _beta, const char tn) {
/*     //Note this function is exactly the opposite of MatrixTimeVector (which is based on the MKL default).   ansi-c*/
/*     //   ansi-c*/
/*     //Output: x_dv->v = _alpha*b_dv*A_dm + _beta*x_dv  for tn=='N'; x_dv = _alpha*b_dv*A_dm' + _beta*x_dv  for tn=='T'   ansi-c*/
/*     //  where x_dv->v is 1-by-ncols or 1-by-nrows and needs not be initialized outside this function if _beta is set to 0.0.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm->M: nrows-by-ncols;   ansi-c*/
/*     //  b_dv->v: 1-by-nrows or 1-by-ncols;   ansi-c*/
/*     //  _alpha: double scalar;   ansi-c*/
/*     //  _beta:  double scalar;   ansi-c*/
/*     //  tn: if =='T' or 't', transpose of A_dm is used; otherwise (=='N' or 'n'), A_dm itself (no transpose) is used.   ansi-c*/

   if ( !x_dv || !A_dm || !b_dv) fn_DisplayError(".../mathlib.c/VectorTimesMatrix(): At least one of the pointer arguments is not created (memory-allocated)");
   else if ( !A_dm->flag || !b_dv->flag )  fn_DisplayError(".../mathlib.c/VectorTimesMatrix(): R input matrix or vector must be given values");

   if ( !(A_dm->flag & M_GE) ) {
      if (A_dm->flag & M_SU)   SUtoGE(A_dm);
      else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
      else  fn_DisplayError(".../mathlib.c/VectorTimesMatrix(): Haven't got time to deal with the M_UT and M_LT cases for A_dm");
   }


   if ( ((tn=='T') || (tn=='t')) && (A_dm->ncols==b_dv->n) && (A_dm->nrows==x_dv->n)) {
      cblas_dgemv(CblasColMajor, CblasNoTrans, A_dm->nrows, A_dm->ncols, _alpha, A_dm->M, A_dm->nrows, b_dv->v, 1, _beta, x_dv->v, 1);
      x_dv->flag = V_DEF;
   }
   else if ( (A_dm->nrows==b_dv->n) && (A_dm->ncols==x_dv->n) ) {
      cblas_dgemv(CblasColMajor, CblasTrans, A_dm->nrows, A_dm->ncols, _alpha, A_dm->M, A_dm->nrows, b_dv->v, 1, _beta, x_dv->v, 1);
      x_dv->flag = V_DEF;
   }
   else {
      fn_DisplayError(".../mathlib.c/VectorTimesMatrix(): Size of the matrix and dimensions of the two vectors do not match");
   }
}

void ScalarTimesMatrix(TSdmatrix *x_dm, const double _alpha, TSdmatrix *a_dm, const double _beta)
{
/*     //$$$$$ If a_dm or x_dm (when _beta!=0) is only upper or lower symmetric, it will be always converted to a general (and symmetric) matrix.  $$$$$$   ansi-c*/
/*     //Output: x_dm = alpha*a_dm + beta*x_dm where x_dm is m-by-n.   ansi-c*/
/*     //  Fastest way is to let beta=0.0 and x_dm->M = a_dm->M.  Then x_dm->M will be replaced by new values.   ansi-c*/
/*     //     However, with beta=0.0, x_dm and a_dm can be different.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dm: m-by-n.   ansi-c*/
/*     //  _alpha: a double scalar.   ansi-c*/
/*     //  _beta: a double scalar.   ansi-c*/
   int _i, _m, _n, nels;
   double *X, *A;

   if ( !x_dm || !a_dm )  fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): All input matrices must be created (memory-allocated)");
   else if ( _beta != 0 && !x_dm->flag )  fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): Input and output matrix must be given legal values because beta != 0");
   else if ( !a_dm->flag )  fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): R input matrix must be given legal values");
   else {
      _m = x_dm->nrows;
      _n = x_dm->ncols;
      nels = _m*_n;
      X = x_dm->M;
      A = a_dm->M;
   }

   if ( !(a_dm->flag & M_GE) ) {
      if (a_dm->flag & M_SU)  SUtoGE(a_dm);
      else if (a_dm->flag & M_SL)  SLtoGE(a_dm);
      else   fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): Haven't got time to deal with the M_UT or M_LT cases for a_dm");
   }
   if ( _beta && !(x_dm->flag & M_GE) ) {
      if (x_dm->flag & M_SU)  SUtoGE(x_dm);
      else if (x_dm->flag & M_SL)  SLtoGE(x_dm);
      else   fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): Haven't got time to deal with the M_UT or M_LT cases for x_dm");
   }

   if ( (_m != a_dm->nrows) || (_n != a_dm->ncols) )  fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): Two input matrices must have the same dimension");
   if (_beta == 0.0) {
      #if defined( SWITCHTOINTELCMATH )             /*   define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
      if ( X == A )  cblas_dscal(nels, _alpha, X, 1);
      else {
         memcpy(X, A, nels*sizeof(double));
         x_dm->flag = a_dm->flag;
         cblas_dscal(nels, _alpha, X, 1);
      }
      #else    /*  My own C math library (which is faster than MKL sometimes); undef: use others.   ansi-c*/
      for (_i=nels-1; _i>=0; _i--)  X[_i] = _alpha*A[_i];
      x_dm->flag = a_dm->flag;
      #endif
   }
   else if (_beta == 1.0) {
      #if defined( SWITCHTOINTELCMATH )             /*   define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
      if ( X == A )  for (_i=nels-1; _i>=0; _i--)  X[_i] += _alpha*A[_i];   /*  fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): to use cblas_daxpy(), the two input matrices must be different, i.e., pointing to different memory places");   ansi-c*/
      else  cblas_daxpy(nels, _alpha, A, 1, X, 1);
      #else
      for (_i=nels-1; _i>=0; _i--)  X[_i] += _alpha*A[_i];
      #endif
      if ( x_dm->flag != a_dm->flag )  x_dm->flag = M_GE;
   }
   else if (_alpha == _beta) {
/*        //if ( X == A )  fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): Two input matrices must be different, i.e., pointing to different memory places");   ansi-c*/
/*        //=== Intel library is at least twice as slow as my own loop for a large matarix like 100-by-100 and can be 20 times slower for a small matrix like 15-by-15.   ansi-c*/
/*        //=== So I don't use Intel library for this.   ansi-c*/
/*        // #ifdef SWITCHTOINTELCMATH             // define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
/*        //    cblas_daxpy(nels, 1.0, A, 1, X, 1);   ansi-c*/
/*        //    cblas_dscal(nels, _alpha, X, 1);   ansi-c*/
/*        // #endif   ansi-c*/
/*        // #ifdef SWITCHTOTZCMATH            // define: use my own C math library (which is faster than MKL sometimes); undef: use others.   ansi-c*/
/*        //    for (_i=nels-1; _i>=0; _i--)  X[_i] = _alpha*(A[_i] + X[_i]);   ansi-c*/
/*        // #endif   ansi-c*/
      if (!(x_dm->flag & M_GE)) fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): x_dm must be M_GE"
                                                "  -- have not got time to convert other matrices to a general matrix");
      for (_i=nels-1; _i>=0; _i--)  X[_i] = _alpha*(A[_i] + X[_i]);
      if ( x_dm->flag != a_dm->flag )  x_dm->flag = M_GE;
   }
   else {
/*        //if ( X == A )  fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): Two input matrices must be different, i.e., pointing to different memory places");   ansi-c*/
      if (!(x_dm->flag & M_GE)) fn_DisplayError(".../mathlib.c/ScalarTimesMatrix(): x_dm must be M_GE"
                                                "  -- have not got time to convert other matrices to a general matrix");
      for (_i=nels-1; _i>=0; _i--)  X[_i] = _alpha*A[_i] + _beta*X[_i];
      if ( x_dm->flag != a_dm->flag )  x_dm->flag = M_GE;
   }
}
/*  //---   ansi-c*/
void ScalarTimesMatrixSquare(TSdmatrix *B_dm, const double _alpha, TSdmatrix *A_dm, const char tn, const double _beta)
{
/*     //Outputs:   ansi-c*/
/*     //  B = alpha*o(A) + beta*B, where o(A) = A' if tn=='T' or 't' or A if tn=='N' or 'n'.   ansi-c*/
/*     //  If A=B, then A is replaced by alpha*o(A) + beta*A.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: n-by-n square matrix.   ansi-c*/
/*     //  B_dm: n-by-n square matrix.   ansi-c*/
/*     //  tn: 'T' (transpose of A) or 'N' (no transpose).   ansi-c*/
/*     //  alpha, beta: double scalars.   ansi-c*/

   int _n = A_dm->nrows;
   TSdmatrix *Atran_dm = NULL;

   if (!A_dm || !B_dm)  fn_DisplayError(".../mathlib.c/ScalarTimesMatrixSquare(): A_dm and B_dm must be created (memory-allocated)");
   if (A_dm->nrows != A_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesMatrixSquare(): A_dm must be square");

   if (A_dm->M == B_dm->M)
   {
      if ((tn == 'T') || (tn == 't'))
      {
         Atran_dm = CreateMatrix_lf(_n,_n);
         TransposeSquare(Atran_dm, A_dm);
         ScalarTimesMatrix(A_dm, _alpha, Atran_dm, _beta);
      }
      else
      {
         ScalarTimesMatrix(A_dm, _alpha, A_dm, _beta);
      }
   }
   else
   {
      if ((B_dm->nrows != B_dm->ncols) || (B_dm->nrows != A_dm->nrows))  fn_DisplayError(".../mathlib.c/ScalarTimesMatrixSquare(): B_dm must be square and B_dm=A_dm");
      if ((tn == 'T') || (tn == 't'))
      {
         Atran_dm = CreateMatrix_lf(_n,_n);
         TransposeSquare(Atran_dm, A_dm);
         ScalarTimesMatrix(B_dm, _alpha, Atran_dm, _beta);
      }
      else
      {
         ScalarTimesMatrix(B_dm, _alpha, A_dm, _beta);
      }
   }

/*     //===   ansi-c*/
   DestroyMatrix_lf(Atran_dm);
}



void MatrixTimesSelf(TSdmatrix *C_dm, const char ul, TSdmatrix *A_dm, const char tn, const double _alpha, const double _beta)
{
/*     //If tn=='N' or 'n', C = alpha*A*A' + beta*C.   ansi-c*/
/*     //If tn=='T' or 't', C = alpha*A'*A + beta*C.   ansi-c*/
/*     //If ul=='U' or 'u', C_dm->flag = M_SU;   ansi-c*/
/*     //If ul=='L' or 'l', C_dm->flag = M_SL;   ansi-c*/
/*     //  C must be different from A.   ansi-c*/
/*     //  C is n-by-n;   ansi-c*/
/*     //  A is n-by-k if tn=='N';   ansi-c*/
/*     //       k-by-n if tn=='T';   ansi-c*/
/*     //  alpha is a double scalar,   ansi-c*/
/*     //  beta is a double scalar.   ansi-c*/
   int _n, _k, lda;

   if ( !C_dm || !A_dm )  fn_DisplayError(".../mathlib.c/MatrixTimesSelf): All input and output matrices must be created (memory-allocated)");
   else if ( !A_dm->flag )  fn_DisplayError(".../mathlib.c/MatrixTimesSelf): Both R input matrices must be given values");
/*     //=== Making this matrix general if not yet.   ansi-c*/
   if ( !(A_dm->flag & M_GE) ) {
      if (A_dm->flag & M_SU)   SUtoGE(A_dm);
      else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixTimesSelf): Haven't got time to deal with the M_UT or M_LT cases for A_dm");
   }


   if ((tn=='T') || (tn=='t')) {
      if (((_n=C_dm->nrows) != A_dm->ncols))  fn_DisplayError(".../mathlib.c/MatrixTimesSelf): Dimensions of A and C do not match where C = alpha*A'*A + beta*C");
      else if (_n != C_dm->ncols)  fn_DisplayError(".../mathlib.c/MatrixTimesSelf): Output matrix C must be a square matrix");
      lda = _k = A_dm->nrows;
   }
   else {
      if ((_n=C_dm->nrows) != (lda=A_dm->nrows))  fn_DisplayError(".../mathlib.c/MatrixTimesSelf): Dimensions of A and C do not match where C = alpha*A*A' + beta*C");
      else if (_n != C_dm->ncols)  fn_DisplayError(".../mathlib.c/MatrixTimesSelf): Output matrix C must be a square matrix");
      _k = A_dm->ncols;
   }

   cblas_dsyrk(CblasColMajor, ((ul=='U') || (ul=='u')) ? CblasUpper : CblasLower, ((tn=='T') || (tn=='t')) ? CblasTrans : CblasNoTrans, _n, _k, _alpha, A_dm->M, lda, _beta, C_dm->M, _n);
   C_dm->flag = ((ul=='U') || (ul=='u')) ? M_SU : M_SL;
}


void MatrixTimesMatrix(TSdmatrix *C_dm, TSdmatrix *A_dm, TSdmatrix *B_dm, const double _alpha, const double _beta, const char tn1, const char tn2) {
/*     //Output is C and all other arguments are inputs.   ansi-c*/
/*     //Computes C = alpah*op(A)*op(B) + beta*C where op() is either transpose or not, depending on 't' or 'n',   ansi-c*/
/*     //  op(A) is m-by-k,   ansi-c*/
/*     //  op(B) is k-by-n,   ansi-c*/
/*     //  C is m-by-n,   ansi-c*/
/*     //  C must be different from A and from B.   ansi-c*/
/*     //  A and B can be the same, however.   ansi-c*/
/*     //  alpha is a double scalar,   ansi-c*/
/*     //  beta is a double scalar.   ansi-c*/
/*     //  tn1: if == 'T' or 't', the transpose of A is used; otherwise (== 'N' or 'n'), A itself (no transpose) is used.   ansi-c*/
/*     //  tn2: if == 'T' or 't', the transpose of B is used; otherwise (== 'N' or 'n'), B itself (no transpose) is used.   ansi-c*/
   int m1, n1, m2, n2, m3, n3;

   if ( !C_dm || !A_dm || !B_dm )  fn_DisplayError(".../mathlib.c/MatrixTimesMatrix(): All input and output matrices must be created (memory-allocated)");
   else if ( !A_dm->flag || !B_dm->flag )  fn_DisplayError(".../mathlib.c/MatrixTimesMatrix(): Both R input matrices must be given values");
   else {
      m1 = A_dm->nrows;
      n1 = A_dm->ncols;
      m2 = B_dm->nrows;
      n2 = B_dm->ncols;
      m3 = C_dm->nrows;
      n3 = C_dm->ncols;
   }
   if ( (_beta != 0.0) && !(C_dm->flag) )  fn_DisplayError(".../mathlib.c/MatrixTimesMatrix(): L input matrix C_dm must be given values when beta !=0.0 (i.e., when C_dm is to be particially updated)");


/*     //=== Making these matrices general if not yet.  For complete symmetric matrix multiplications, do use this function.   ansi-c*/
   if ( !(A_dm->flag & M_GE) ) {
      if (A_dm->flag & M_SU)   SUtoGE(A_dm);
      else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixTimesMatrix(): Haven't got time to deal with the M_UT or M_LT cases for A_dm");
   }
   if ( !(B_dm->flag & M_GE) ) {
      if (B_dm->flag & M_SU)   SUtoGE(B_dm);
      else if (B_dm->flag & M_SL)   SLtoGE(B_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixTimesMatrix(): Haven't got time to deal with the M_UT or M_LT cases for B_dm");
   }



   if ( ((tn1=='N') || (tn1=='n')) && ((tn2=='N') || (tn2=='n')) && (n1 == m2) && (m3 == m1) && (n3 == n2) ) {
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m3, n3, n1, _alpha, A_dm->M, m1, B_dm->M, m2, _beta, C_dm->M, m3);
      C_dm->flag = M_GE;
   }
   else if ( ((tn1=='T') || (tn1=='t')) && ((tn2=='N') || (tn2=='n')) && (m1 == m2) && (m3 == n1) && (n3 == n2) ) {
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m3, n3, m1, _alpha, A_dm->M, m1, B_dm->M, m2, _beta, C_dm->M, m3);
      C_dm->flag = M_GE;
   }
   else if ( ((tn1=='T') || (tn1=='t')) && ((tn2=='T') || (tn2=='t')) && (m1 == n2) && (m3 == n1) && (n3 == m2) ) {
      cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, m3, n3, m1, _alpha, A_dm->M, m1, B_dm->M, m2, _beta, C_dm->M, m3);
      C_dm->flag = M_GE;
   }
   else if ( ((tn1=='N') || (tn1=='n')) && ((tn2=='T') || (tn2=='t')) && (n1 == n2) && (m3 == m1) && (n3 == m2) ) {
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m3, n3, n1, _alpha, A_dm->M, m1, B_dm->M, m2, _beta, C_dm->M, m3);
      C_dm->flag = M_GE;
   }
   else  fn_DisplayError(".../mathlib.c/MatrixTimesMatrix(): (1) Dimensions of both R input matrices must match.  (2) Dimension of L input matrix must compatible with the multiplication of the two R input matrices");
}

void SolveTriSysVector(TSdvector *x_dv, const TSdmatrix *T_dm, TSdvector *b_dv, const char tn, const char un)
{
/*     //Output:  computes x_dv = inv(T_dm)*b_dv by solving a triangular system of equation T_dm * x_dv = b_dv.   ansi-c*/
/*     //  x_dv(_n-by-1) = inv(T_dm)*b_v if tn=='N'; = inv(T_dm')*b_v if tn=='T'.   ansi-c*/
/*     // Fastest way is to let x_dv->v = b_dv->v.  Then, x_dv->v will be replaced by new values.   ansi-c*/
/*     //-------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  T_dm: _n-by-_n upper or lower triangular matrix;   ansi-c*/
/*     //  b_dv: _n-by-1 vector.   ansi-c*/
/*     //  tn: if =='T' or 't', T_dm->M' (transpose), instead of T_m, will be used; otherwise (i.e., =='n' or 'N'), T_dm->M itself (no transpose) will be used.   ansi-c*/
/*     //  un: if =='U' or 'u', T_dm is a unit upper triangular (i.e., the diagonal being 1);   ansi-c*/
/*     //      otherwise (i.e., if =='N' or 'n'), T_dm is a non-unit upper triangular.   ansi-c*/
/*     //   ansi-c*/
/*     // Note I: Intel MLK cblas_dtrsv() does not test for singularity or near-singulariy of the system.   ansi-c*/
/*     //   Such tests must be performed before calling this BLAS routine.   ansi-c*/
/*     // Note II: if x_dv->v = b_dv->v, x_dv->v will be replaced by new values.   ansi-c*/
   int _n, _i;
   double *x, *b;

   if ( !T_dm || !b_dv || !x_dv )  fn_DisplayError(".../mathlib.c/SolveTriSysVector(): All input pointers must be created (memory-allocated)");
   else if ( !( T_dm->flag & (M_UT | M_LT) ) || !b_dv->flag )  fn_DisplayError(".../mathlib.c/SolveTriSysVector(): (1) R input matrix must be triangular.  (2) R input vector must be given legal values");
   else {
      _n = T_dm->nrows;
      x = x_dv->v;
      b = b_dv->v;
   }

   for (_i=square(_n)-1; _i>=0; _i -= _n+1)
      if (fabs(T_dm->M[_i])<=MACHINEZERO)
         fn_DisplayError(".../mathlib.c/SolveTriSysVector():  The input triangular matrix is singular");


   if ( (_n != T_dm->ncols) || (_n != b_dv->n) || (_n != x_dv->n) )   fn_DisplayError(".../mathlib.c/SolveTriSysVector(): (1) R input matrix must be square. (2) Input vectors and the square matrix must have the same dimension");
   else if ( x == b) {
      cblas_dtrsv(CblasColMajor, ( T_dm->flag & M_UT ) ? CblasUpper : CblasLower,
                                 ((tn=='T') || (tn=='t')) ? CblasTrans : CblasNoTrans,
                                 ((un=='U') || (tn=='u')) ? CblasUnit : CblasNonUnit,
                                                                        _n, T_dm->M, _n, x, 1);
   }
   else {
      memcpy(x, b, _n*sizeof(double));
      cblas_dtrsv(CblasColMajor, ( T_dm->flag & M_UT ) ? CblasUpper : CblasLower,
                                 ((tn=='T') || (tn=='t')) ? CblasTrans : CblasNoTrans,
                                 ((un=='U') || (tn=='u')) ? CblasUnit : CblasNonUnit,
                                                                        _n, T_dm->M, _n, x, 1);
      x_dv->flag = V_DEF;
   }
}






void SymmetricMatrixTimesVector(double *x_v, const double a, const double *A_m, const double *a_v, const double b,  const int _n, const char ul) {
/*     //Output: x_v = a*A_m*a_v + b*X_m where x_v (_n-by-1) must be allocated (but needs not be initialized).   ansi-c*/
/*     //Inputs:         X_m??????????????????????   ansi-c*/
/*     //  A_m: _n-by-_n symmetric matrix;   ansi-c*/
/*     //  a_v: _n-by-1;   ansi-c*/
/*     //  a, b: scalars;   ansi-c*/
/*     //  ul: if =='u' or 'U', upper triangular elements in A_m are filled; if =='l' or 'L', lower triangular elements in A_m are filled.   ansi-c*/

   cblas_dsymv(CblasColMajor, ((ul=='u') || (ul=='U')) ? CblasUpper : CblasLower, _n, a, A_m, _n, a_v, 1, b, x_v, 1);
}


void SolveTriangularSystemVector(double *x_v, const double *A_m, const double *b_v, const int _n, const char ul, const char tn, const char un) {
/*     //Outputs:   ansi-c*/
/*     //  x_v(_n-by-1) = inv(A_m)*b_v.  If x_v=b_v, b_v will be overwritten by x_v.   ansi-c*/
/*     //-------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_m: _n-by-_n upper or lower triangular matrix;   ansi-c*/
/*     //  b_v: _n-by-1 vector.   ansi-c*/
/*     //  ul: if =='u' or 'U', A_m is upper triangular; if =='l' or 'L', A_m is lower triangular.   ansi-c*/
/*     //  tn: if =='t' or 'T', A_m' (transpose), instead of A_m, will be used; if =='n', A_m itself (no transpose) will be used.   ansi-c*/
/*     //  un: if =='u' or 'U', A_m is a unit upper triangular (i.e., the diagonal being 1);   ansi-c*/
/*     //      if =='n' or 'N', A_m is a non-unit upper triangular.   ansi-c*/
/*     //   ansi-c*/
/*     // Computes x_v = inv(A_m)*b_v by solving a triangular system of equation A_m * x_v = b_v.   ansi-c*/
/*     // Note I: Intel MLK cblas_dtrsv() does not test for singularity or near-singulariy of the system.   ansi-c*/
/*     //   Such tests must be performed before calling this BLAS routine.   ansi-c*/
/*     // Note II: if x_v=b_v, b_v will be overwritten by x_v.   ansi-c*/

   if (x_v != b_v)  memcpy(x_v, b_v, _n*sizeof(double));
   cblas_dtrsv(CblasColMajor, ((ul=='u') || (ul=='U')) ? CblasUpper : CblasLower,
                              ((tn=='t') || (tn=='T')) ? CblasTrans : CblasNoTrans,
                              ((un=='u') || (tn=='U')) ? CblasNonUnit : CblasNonUnit,
                                                                      _n, A_m, _n, x_v, 1);
}





/*  //=======================================================   ansi-c*/
/*  // MKL Vector Mathematical Library with default using my own routines.   ansi-c*/
/*  //=======================================================   ansi-c*/
void VectorDotDivByVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv) {
/*     //????????? NOT tested yet.  06/13/03.   ansi-c*/
/*     //Output: x_dv = a_dv ./ b_dv (division element by elment) where x_dv, a_dv, and b_dv are all _n-by-1.   ansi-c*/
/*     //  The fastest way is to use MKL VML with x != a and x != b.   ansi-c*/
/*     //  If x_dv = a_dv, x_dv will be replaced by x_dv ./ b_dv.   ansi-c*/
/*     //  If x_dv = b_dv, x_dv will be replaced by a_dv ./ x_dv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: _n-by-1 double vector.   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
   int _i, _n;
   double *x, *a, *b;


   if ( !x_dv || !a_dv || !b_dv)  fn_DisplayError(".../mathlib.c/VectorDotDivideVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag )  fn_DisplayError(".../mathlib.c/VectorDotDivideVector(): R input vectors must be given legal values");
   else if ( ((_n=x_dv->n) != a_dv->n) || (_n != b_dv->n) )  fn_DisplayError(".../mathlib.c/VectorDotDivideVector(): Dimensions of all input vectors must be same");
   else {
      x = x_dv->v;
      a = a_dv->v;
      b = b_dv->v;
   }


   #if !defined (INTELCMATHLIBRARY)
      if ( (x != a) && (x != b) )  {
         vdDiv (_n, a, b, x);
         x_dv->flag = V_DEF;
      }
      else
         for (_i=_n-1; _i>=0; _i--)   x[_i] = a[_i]/b[_i];
   #else     /*  Default to my own routine.   ansi-c*/
      for (_i=_n-1; _i>=0; _i--)   x[_i] = a[_i]/b[_i];
      x_dv->flag = V_DEF;
   #endif
}

void ElementwiseVectorDivideVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv)
{
/*     //The fastest way is to use MKL VML with x != a and x != b.   ansi-c*/
/*     //Output: x_dv = a_dv ./ b_dv (division element by elment) where x_dv, a_dv, and b_dv are all _n-by-1.   ansi-c*/
/*     //  If x_dv = a_dv, x_dv will be replaced by x_dv ./ b_dv.   ansi-c*/
/*     //  If x_dv = b_dv, x_dv will be replaced by a_dv ./ x_dv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: _n-by-1 double vector.   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
   int _i, _n;
   double *x, *a, *b;

   if ( !x_dv || !a_dv || !b_dv)  fn_DisplayError("mathlib.c/ElementwiseVectorDivideVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag )  fn_DisplayError("mathlib.c/ElementwiseVectorDivideVector(): R input vectors must be given legal values");
   else if ( ((_n=x_dv->n) != a_dv->n) || (_n != b_dv->n) )  fn_DisplayError("mathlib.c/ElementwiseVectorDivideVector(): Dimensions of all input vectors must be same");
   else {
      x = x_dv->v;
      a = a_dv->v;
      b = b_dv->v;
   }


   #if !defined (INTELCMATHLIBRARY)
      if ( (x != a) && (x != b) )  {
         vdDiv (_n, a, b, x);
         x_dv->flag = V_DEF;
      }
      else
         for (_i=_n-1; _i>=0; _i--)   x[_i] = a[_i]/b[_i];
   #else     /*  Default to my own routine.   ansi-c*/
      for (_i=_n-1; _i>=0; _i--)   x[_i] = a[_i]/b[_i];
      x_dv->flag = V_DEF;
   #endif
}

void ElementwiseInverseofVector(TSdvector *y_dv, TSdvector *x_dv) {
/*     //The fastest way is to use MKL VML with y_dv != x_dv;   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If y_dv!=x_dv, y_dv = 1 ./ x_dv;   ansi-c*/
/*     //  If y_dv=x_dv, x_dv = 1 ./ x_dv.   ansi-c*/

   int _i;
   #if defined( INTELCMATHLIBRARY )
   int _n;
   double *y;
   #endif
   double *x;

   if ( !x_dv || !x_dv->flag )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofVector(): (1) Input vector must be memory-allocated; (2) Legal values must be given");

   #if !defined( INTELCMATHLIBRARY )
      if ( y_dv == x_dv ) {
         x = x_dv->v;
         for (_i=x_dv->n-1; _i>=0; _i--)  x[_i] = 1.0/x[_i];
      }
      else  {
         if ( !y_dv )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofVector(): Output vector must be memory-allocated (but no need for legal values)");
         if ( x_dv->n != y_dv->n)  fn_DisplayError(".../mathlib.c/ElementwiseInverseofVector(): Lengths of both input and output vectors must be same");
         vdInv (x_dv->n, x_dv->v, y_dv->v);
         y_dv->flag = V_DEF;
      }
   #else
      if ( y_dv == x_dv ) {
         x = x_dv->v;
         for (_i=x_dv->n-1; _i>=0; _i--)  x[_i] = 1.0/x[_i];
      }
      else  {
         if ( !y_dv )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofVector(): Output vector must be memory-allocated (but no need for legal values)");
         if ( (_n=x_dv->n) != y_dv->n)  fn_DisplayError(".../mathlib.c/ElementwiseInverseofVector(): Lengths of both input and output vectors must be same");
         x = x_dv->v;
         y = y_dv->v;
         for (_i=_n-1; _i>=0; _i--)  y[_i] = 1.0/x[_i];
         y_dv->flag = V_DEF;
      }
   #endif
}

void ElementwiseSqrtofVector(TSdvector *y_dv, TSdvector *x_dv)
{
/*     //The fastest way is to use MKL VML with y_dv != x_dv;   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If y_dv!=x_dv, y_dv = sqrt(x_dv);   ansi-c*/
/*     //  If y_dv=x_dv, x_dv = sqrt(x_dv);   ansi-c*/

   int _i;
   #if defined( INTELCMATHLIBRARY )
   int _n;
   double *y;
   #endif
   double *x;

   if ( !x_dv || !x_dv->flag )  fn_DisplayError("mathlib.c/ElementwiseSqrtofVector(): (1) Input vector must be memory-allocated; (2) Legal values must be given");

   #if !defined( INTELCMATHLIBRARY )
      if ( y_dv == x_dv ) {
         x = x_dv->v;
         for (_i=x_dv->n-1; _i>=0; _i--)  x[_i] = sqrt(x[_i]);
      }
      else  {
         if ( !y_dv )  fn_DisplayError("mathlib.c/ElementwiseSqrtofVector(): Output vector must be memory-allocated (but no need for legal values)");
         if ( x_dv->n != y_dv->n)  fn_DisplayError("mathlib.c/ElementwiseSqrtofVector(): Lengths of both input and output vectors must be same");
         vdSqrt(x_dv->n, x_dv->v, y_dv->v);
         y_dv->flag = V_DEF;
      }
   #else
      if ( y_dv == x_dv ) {
         x = x_dv->v;
         for (_i=x_dv->n-1; _i>=0; _i--)  x[_i] = sqrt(x[_i]);
      }
      else  {
         if ( !y_dv )  fn_DisplayError("mathlib.c/ElementwiseSqrtofVector(): Output vector must be memory-allocated (but no need for legal values)");
         if ( (_n=x_dv->n) != y_dv->n)  fn_DisplayError("mathlib.c/ElementwiseSqrtofVector(): Lengths of both input and output vectors must be same");
         x = x_dv->v;
         y = y_dv->v;
         for (_i=_n-1; _i>=0; _i--)  y[_i] = sqrt(x[_i]);
         y_dv->flag = V_DEF;
      }
   #endif
}

void ElementwiseLogofVector(TSdvector *y_dv, TSdvector *x_dv)
{
/*     //The fastest way is to use MKL VML with y_dv != x_dv;   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If y_dv!=x_dv, y_dv = log(x_dv);   ansi-c*/
/*     //  If y_dv=x_dv, x_dv = log(x_dv);   ansi-c*/

   int _i;
   #if defined( INTELCMATHLIBRARY )
   int _n;
   double *y;
   #endif
   double *x;

   if ( !x_dv || !x_dv->flag )  fn_DisplayError("mathlib.c/ElementwiseLogofVector(): (1) Input vector must be memory-allocated; (2) Legal values must be given");

   #if !defined( INTELCMATHLIBRARY )
      if ( y_dv == x_dv ) {
         x = x_dv->v;
         for (_i=x_dv->n-1; _i>=0; _i--)  x[_i] = log(x[_i]);
      }
      else  {
         if ( !y_dv )  fn_DisplayError("mathlib.c/ElementwiseLogofVector(): Output vector must be memory-allocated (but no need for legal values)");
         if ( x_dv->n != y_dv->n)  fn_DisplayError("mathlib.c/ElementwiseLogofVector(): Lengths of both input and output vectors must be same");
         vdLn(x_dv->n, x_dv->v, y_dv->v);
         y_dv->flag = V_DEF;
      }
   #else
      if ( y_dv == x_dv ) {
         x = x_dv->v;
         for (_i=x_dv->n-1; _i>=0; _i--)  x[_i] = log(x[_i]);
      }
      else  {
         if ( !y_dv )  fn_DisplayError("mathlib.c/ElementwiseLogofVector(): Output vector must be memory-allocated (but no need for legal values)");
         if ( (_n=x_dv->n) != y_dv->n)  fn_DisplayError("mathlib.c/ElementwiseLogofVector(): Lengths of both input and output vectors must be same");
         x = x_dv->v;
         y = y_dv->v;
         for (_i=_n-1; _i>=0; _i--)  y[_i] = log(x[_i]);
         y_dv->flag = V_DEF;
      }
   #endif
}


void ElementwiseInverseofMatrix(TSdmatrix *Y_dm, TSdmatrix *X_dm)
{
/*     //The fastest way is to use MKL VML with Y_dm != X_dm;   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If Y_dm!=X_dm, Y_dm = 1 ./ X_dm;   ansi-c*/
/*     //  If Y_dm=X_dm, X_dm = 1 ./ X_dm.   ansi-c*/

   int _i,
       nrows, ncols;
   double *X;
   #if defined( INTELCMATHLIBRARY )
   double *Y;
   #endif



   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofMatrix(): (1) Input matrix must be memory-allocated; (2) Legal values must be given");

   #if !defined( INTELCMATHLIBRARY )
      if ( Y_dm == X_dm ) {
         X = X_dm->M;
         for (_i=X_dm->nrows*X_dm->ncols-1; _i>=0; _i--)  X[_i] = 1.0/X[_i];
      }
      else  {
         if ( !Y_dm )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofMatrix(): Output matrix must be memory-allocated (but no need for legal values)");
         if ( ((nrows=X_dm->nrows) != Y_dm->nrows) || ((ncols=X_dm->ncols) != Y_dm->ncols) )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofMatrix(): Dimensions of both input and output matrices must be same");
         vdInv (nrows*ncols, X_dm->M, Y_dm->M);
         Y_dm->flag = M_GE;
      }
   #else     /*  Default to my own routine.   ansi-c*/
      if ( Y_dm == X_dm ) {
         X = X_dm->M;
         for (_i=X_dm->nrows*X_dm->ncols-1; _i>=0; _i--)  X[_i] = 1.0/X[_i];
      }
      else  {
         if ( !Y_dm )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofMatrix(): Output matrix must be memory-allocated (but no need for legal values)");
         if ( ((nrows=X_dm->nrows) != Y_dm->nrows) || ((ncols=X_dm->ncols) != Y_dm->ncols) )  fn_DisplayError(".../mathlib.c/ElementwiseInverseofMatrix(): Dimensions of both input and output matrices must be same");
         X = X_dm->M;
         Y = Y_dm->M;
         for (_i=nrows*ncols-1; _i>=0; _i--)  Y[_i] = 1.0/X[_i];
         Y_dm->flag = M_GE;
      }
   #endif
}



/*  //=======================================================   ansi-c*/
/*  // Matrix routines (my own).   ansi-c*/
/*  //=======================================================   ansi-c*/
void tz_VectorPlusMinusVector(TSdvector *x_dv, const TSdvector *a_dv, const double _alpha, const TSdvector *b_dv, const double _beta)
{
/*     //Output: x_dv = alpha*a_dv + beta*b_dv where x_dv is _n-by-1.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: _n-by-1 double vector.   ansi-c*/
/*     //  _alpha: double constant.   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
/*     // _beta: double constant.   ansi-c*/
   int _i, _n;
   double *x, *a, *b;


   if ( !x_dv || !a_dv || !b_dv) fn_DisplayError("mathlib.c/tz_VectorPlusMinusVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag ) fn_DisplayError("mathlib.c/tz_VectorPlusMinusVector(): R input vectors must be given values");
   else {
      _n = x_dv->n;
      x = x_dv->v;
      a = a_dv->v;
      b = b_dv->v;
   }
   if ( (_n != a_dv->n) || (_n != b_dv->n) ) fn_DisplayError("mathlib.c/tz_VectorPlusMinusVector(): Dimensions of all input vectors must be same");
   else {
      for (_i=_n-1; _i>=0; _i--)  x[_i] = _alpha*a[_i] + _beta*b[_i];
      x_dv->flag = V_DEF;
   }
}
void VectorPlusVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv) {
/*     //Output: x_dv = a_dv + b_dv where x_dv is _n-by-1.   ansi-c*/
/*     //          If x_dv = a_dv, a_dv will be replaced by x_dv.   ansi-c*/
/*     //          If x_dv = b_dv, b_dv will be replaced by x_dv,   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: _n-by-1 double vector.   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
   int _i, _n;
   double *x, *a, *b;


   if ( !x_dv || !a_dv || !b_dv) fn_DisplayError(".../mathlib.c/VectorPlusVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag ) fn_DisplayError(".../mathlib.c/VectorPlusVector(): R input vectors must be given values");
   else {
      _n = x_dv->n;
      x = x_dv->v;
      a = a_dv->v;
      b = b_dv->v;
   }
   if ( (_n != a_dv->n) || (_n != b_dv->n) ) fn_DisplayError(".../mathlib.c/VectorPlusVector(): Dimensions of all input vectors must be same");
   else {
      for (_i=_n-1; _i>=0; _i--)  x[_i] = a[_i] + b[_i];
      x_dv->flag = V_DEF;
   }
}
void VectorMinusVector(TSdvector *x_dv, const TSdvector *a_dv, const TSdvector *b_dv)
{
/*     //Output: x_dv = a_dv - b_dv where x_dv is _n-by-1.   ansi-c*/
/*     //  If x_dv = a_dv, x_dv will be replaced by x_dv - b_dv.   ansi-c*/
/*     //  If x_dv = b_dv, x_dv will be replaced by a_dv - x_dv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: _n-by-1 double vector.   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
   int _i, _n;
   double *x, *a, *b;


   if ( !x_dv || !a_dv || !b_dv) fn_DisplayError(".../mathlib.c/VectorMinusVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag ) fn_DisplayError(".../mathlib.c/VectorMinusVector(): R input vectors must be given values");
   else {
      x = x_dv->v;
      a = a_dv->v;
      b = b_dv->v;
   }
   if ( ((_n = x_dv->n) != a_dv->n) || (_n != b_dv->n) ) fn_DisplayError(".../mathlib.c/VectorMinusVector(): Dimensions of all input vectors must be same");
   else {
      for (_i=_n-1; _i>=0; _i--)  x[_i] = a[_i] - b[_i];
      x_dv->flag = V_DEF;
   }
}


void VectorPlusVectorUpdate(TSdvector *x_dv, const TSdvector *b_dv) {
/*     //Output: x_dv = b_dv + x_dv where x_dv is _n-by-1.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
   int _n, _i;
   double *x, *b;


   if ( !x_dv || !b_dv )  fn_DisplayError(".../mathlib.c/VectorPlusVectorUpdate(): All input vectors must be created (memory-allocated)");
   if ( !b_dv->flag || !x_dv->flag )  fn_DisplayError(".../mathlib.c/VectorPlusVectorUpdate(): All input vectors must be given values");
   if ( (_n=x_dv->n) != b_dv->n )  fn_DisplayError(".../mathlib.c/VectorPlusVectorUpdate(): Dimensions of all input vectors must be same");

   x = x_dv->v;
   b = b_dv->v;
   for (_i=_n-1; _i>=0; _i--)  x[_i] += b[_i];
}


void VectorDotTimesVector(TSdvector *x_dv, const TSdvector *a_dv, TSdvector *b_dv, const double _alpha, const double _beta) {
/*     //Output:   ansi-c*/
/*     //  x_dv is _n-by-1.   ansi-c*/
/*     //  x_dv = _alpha * a_dv .* b_dv + _beta * x_dv if x_dv != b_dv.   ansi-c*/
/*     //  x_dv = _alpha * a_dv .* x_dv + _beta * x_dv if x_dv = b_dv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  a_dv: _n-by-1 double vector.   ansi-c*/
/*     //  b_dv: _n-by-1 double vector.   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
/*     //  _beta: a double scalar.   ansi-c*/
   int _i, _n;
   double *x, *a, *b;


   if ( !x_dv || !a_dv || !b_dv) fn_DisplayError(".../mathlib.c/VectorDotTimesVector(): All input vectors must be created (memory-allocated)");
   else if ( !a_dv->flag || !b_dv->flag )   fn_DisplayError(".../mathlib.c/VectorDotTimesVector(): Both R input vectors must be given values");
   else {
      _n = x_dv->n;
      x = x_dv->v;
      a = a_dv->v;
      b = b_dv->v;
   }
   if ( (_n != a_dv->n) || (_n != b_dv->n) ) fn_DisplayError(".../mathlib.c/VectorDotTimesVector(): Dimensions of all input vectors must be same");


   if ( _alpha==1.0 ) {
      if ( _beta==0.0 ) {
         for (_i=_n-1; _i>=0; _i--)  x[_i] = a[_i] * b[_i];
         if (!x_dv->flag)  x_dv->flag = V_DEF;
      }
      else if ( _beta==1.0 ) {
         for (_i=_n-1; _i>=0; _i--)  x[_i] += a[_i] * b[_i];
         if (!x_dv->flag)  x_dv->flag = V_DEF;
      }
      else {
         for (_i=_n-1; _i>=0; _i--)   x[_i] = a[_i] * b[_i] + _beta * x[_i];
         if (!x_dv->flag)  x_dv->flag = V_DEF;
      }
   }
   else {
      if ( _beta==0.0 ) {
         for (_i=_n-1; _i>=0; _i--)   x[_i] = _alpha * a[_i] * b[_i];
         if (!x_dv->flag)  x_dv->flag = V_DEF;
      }
      else if ( _beta==1.0 ) {
         for (_i=_n-1; _i>=0; _i--)  x[_i] += _alpha * a[_i] * b[_i];
         if (!x_dv->flag)  x_dv->flag = V_DEF;
      }
      else {
         for (_i=_n-1; _i>=0; _i--)  x[_i] = _alpha * a[_i] * b[_i] + _beta * x[_i];
         if (!x_dv->flag)  x_dv->flag = V_DEF;
      }
   }
}

void SwapColsofMatrix(TSdmatrix *X_dm, int j1, int j2)
{
/*     //??????? NOT tested yet.   ansi-c*/
/*     //Ouputs:   ansi-c*/
/*     //  The j1_th column of X_dm is swapped with the j2_th column of X_dm.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  X_dm:  Memory allocated and legal values given already.   ansi-c*/
/*     //  j1: The j1_th column of X_dm.   ansi-c*/
/*     //  j2: The j2_th column of X_dm.   ansi-c*/

   int nrows;
   double *M1, *M2;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): (1) Input matrix X must be created (memory-allocated); (2) Legal values must be given");
   if (j1 >= X_dm->ncols)  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): The j1_th column specified for the input matrix X exceeds its column dimension");
   if (j2 >= X_dm->ncols)  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): The j2_th column specified for the input matrix X exceeds its column dimension");
   if (j1 == j2)  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): The two columns for swapping must be different");

   #if defined( INTELCMATHLIBRARY )
      M1 = X_dm->M + j1*(nrows=X_dm->nrows);   /*  Points to the beginning of the j1_th column.   ansi-c*/
      M2 = X_dm->M + j2*nrows;   /*  Points to the beginning of the j2_th column.   ansi-c*/
      cblas_dswap(nrows, M1, 1, M2, 1);
   #else
      fn_DisplayError(".../mathlib.c/SwapColsofMatrix():  Haven't got time to write my own for-loop routine to swap two columns");
   #endif
}
void SwapColsofMatrices(TSdmatrix *X1_dm, int j1, TSdmatrix *X2_dm, int j2)
{
/*     //Ouputs:   ansi-c*/
/*     //  The j1_th column of X1_dm is swapped with the j2_th column of X2_dm.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  X1_dm:  Memory allocated and legal values given already.   ansi-c*/
/*     //  X2_dm:  Memory allocated and legal values given already.   ansi-c*/
/*     //  j1: The j1_th column of X1_dm.   ansi-c*/
/*     //  j2: The j2_th column of X2_dm.   ansi-c*/

   int nrows;
   double *M1, *M2;

   if ( !X1_dm || !X1_dm->flag )  fn_DisplayError(".../mathlib.c/SwapColsofMatrices(): (1) Input matrix X1 must be created (memory-allocated); (2) Legal values must be given");
   if ( !X2_dm || !X2_dm->flag )  fn_DisplayError(".../mathlib.c/SwapColsofMatrices(): (1) Input matrix X2 must be created (memory-allocated); (2) Legal values must be given");
   if ( X1_dm == X2_dm )  fn_DisplayError(".../mathlib.c/SwapColsofMatrices(): The two input matrices must be different");
   if (j1 >= X1_dm->ncols)  fn_DisplayError(".../mathlib.c/SwapColsofMatrices(): The jth column specified for the input matrix X1 exceeds its column dimension");
   if (j2 >= X2_dm->ncols)  fn_DisplayError(".../mathlib.c/SwapColsofMatrices(): The jth column specified for the input matrix X1 exceeds its column dimension");
   if ( (nrows=X1_dm->nrows) != X2_dm->nrows )  fn_DisplayError(".../mathlib.c/SwapColsofMatrices(): The number of rows for both input matrices must be the same");


   #if defined (INTELCMATHLIBRARY)
      M1 = X1_dm->M + j1*nrows;   /*  Points to the beginning of the j1_th column.   ansi-c*/
      M2 = X2_dm->M + j2*nrows;   /*  Points to the beginning of the j2_th column.   ansi-c*/
      cblas_dswap(nrows, M1, 1, M2, 1);
   #else
      fn_DisplayError(".../mathlib.c/SwapColsofMatrices():  Haven't got time to write my own for-loop routine to swap two columns");
   #endif
}
void SwapPositionsofMatrix(TSdmatrix *X_dm, int j1, int j2) {
/*     //Ouputs:   ansi-c*/
/*     //  Column operation: first, the j1_th column of X_dm is swapped with the j2_th column of X_dm.   ansi-c*/
/*     //  Row operation: second, the j1_th row of X_dm is swapped with the j2_th row of X_dm.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  X_dm:  Memory allocated and legal values given already.   ansi-c*/
/*     //  j1: The j1_th column and row of X_dm.   ansi-c*/
/*     //  j2: The j2_th column and row of X_dm.   ansi-c*/

   int nrows, ncols;
   double *M1, *M2;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): (1) Input matrix X must be created (memory-allocated); (2) Legal values must be given");
   if (j1 >= (ncols=X_dm->ncols) || j1 >= (nrows=X_dm->nrows) )  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): The j1_th column or row specified for the input matrix X exceeds its column or row dimension");
   if (j2 >= X_dm->ncols || j2 >= X_dm->ncols )  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): The j2_th column or row specified for the input matrix X exceeds its column or row dimension");
   if (j1 == j2)  fn_DisplayError(".../mathlib.c/SwapColsofMatrix(): The two columns for swapping must be different");

   #if defined( INTELCMATHLIBRARY )
      M1 = X_dm->M + j1*nrows;   /*  Points to the beginning of the j1_th column.   ansi-c*/
      M2 = X_dm->M + j2*nrows;   /*  Points to the beginning of the j2_th column.   ansi-c*/
      cblas_dswap(nrows, M1, 1, M2, 1);      /*  Swaps columns.   ansi-c*/
/*        //   ansi-c*/
      M1 = X_dm->M + j1;   /*  Points to the beginning of the j1_th row.   ansi-c*/
      M2 = X_dm->M + j2;   /*  Points to the beginning of the j2_th row.   ansi-c*/
      cblas_dswap(ncols, M1, nrows, M2, nrows);   /*  Swaps corresponding rows.   ansi-c*/
   #else
      fn_DisplayError(".../mathlib.c/SwapPositionsofMatrix():  Haven't got time to write my own for-loop routine to swap two columns and then two corresponding rows");
   #endif
}

void SwapMatricesofCell(TSdcell *A_dc, int c1, int c2)
{
/*     //Ouputs:   ansi-c*/
/*     //  A_dc->C[c1] and A_dc->C[c2] are swapped.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dc:  Memory allocated and legal values given already.   ansi-c*/
/*     //  c1, c2:  Positions of cells.   ansi-c*/

   int _n;
   TSdmatrix *tpnter2dm = NULL;

   if ( !A_dc )  fn_DisplayError(".../mathlib.c/SwapMatricesofCell(): Input cell A_dc must be created (memory-allocated)");
   _n = A_dc->ncells;
   if ( (c1>=_n) || (c1<0) || (c2>=_n) || (c2<0) )  fn_DisplayError(".../mathlib.c/SwapMatricesofCell(): c1 and c2 must be between 0 and A_dc->ncells inclusive");

   tpnter2dm = A_dc->C[c1];
   A_dc->C[c1] = A_dc->C[c2];
   A_dc->C[c2] = tpnter2dm;
}

void SwapVectorsofCellvec(TSdcellvec *x_dcv, int c1, int c2)
{
/*     //Ouputs:   ansi-c*/
/*     //  x_dcv->C[c1] and x_dcv->C[2] are swapped.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_dcv:  Memory allocated and legal values given already.   ansi-c*/
/*     //  c1, c2:  Positions of cells.   ansi-c*/

   int _n;
   TSdvector *tpnter2dv = NULL;

   if ( !x_dcv )  fn_DisplayError(".../mathlib.c/SwapVectorsofCellvec(): Input cell vector x_dcv must be created (memory-allocated)");
   _n = x_dcv->ncells;
   if ( (c1>=_n) || (c1<0) || (c2>=_n) || (c2<0) )  fn_DisplayError(".../mathlib.c/SwapVectorsofCellvec(): c1 and c2 must be between 0 and x_dcv->ncells inclusive");

   tpnter2dv = x_dcv->C[c1];
   x_dcv->C[c1] = x_dcv->C[c2];
   x_dcv->C[c2] = tpnter2dv;
}
/*  //--   ansi-c*/
void SwapVectorsofCellvec_int(TSicellvec *x_icv, int c1, int c2)
{
/*     //Ouputs:   ansi-c*/
/*     //  x_icv->C[c1] and x_icv->C[2] are swapped.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_icv:  Memory allocated and legal values given already.   ansi-c*/
/*     //  c1, c2:  Positions of cells.   ansi-c*/

   int _n;
   TSivector *tpnter2iv = NULL;

   if ( !x_icv )  fn_DisplayError(".../mathlib.c/SwapVectorsofCellvec_int(): Input cell vector x_icv must be created (memory-allocated)");
   _n = x_icv->ncells;
   if ( (c1>=_n) || (c1<0) || (c2>=_n) || (c2<0) )  fn_DisplayError(".../mathlib.c/SwapVectorsofCellvec_int(): c1 and c2 must be between 0 and x_icv->ncells inclusive");

   tpnter2iv = x_icv->C[c1];
   x_icv->C[c1] = x_icv->C[c2];
   x_icv->C[c2] = tpnter2iv;
}


/*  //=== The following is NOT efficient.   ansi-c*/
/*  //#if defined( INTELCMATHLIBRARY )   ansi-c*/
/*  //void SwapMatricesofCell(TSdcell *X_dc, int c1, int c2)   ansi-c*/
/*  //{   ansi-c*/
/*  //   //??????? NOT tested yet.   ansi-c*/
/*  //   //Ouputs:   ansi-c*/
/*  //   //  The c1_th matrix of X_dc is swapped with the c2_th matrix of X_dc.   ansi-c*/
/*  //   //Inputs:   ansi-c*/
/*  //   //  X_dc:  Memory allocated and legal values given already.   ansi-c*/
/*  //   //  c1: The c1_th matrix of X_dc.   ansi-c*/
/*  //   //  c2: The c2_th matrix of X_dc.   ansi-c*/

/*  //   int dim;   ansi-c*/


/*  //   if ( !X_dc )  fn_DisplayError(".../mathlib.c/SwapMatricesofCell(): input cell X_dc must be created (memory-allocated)");   ansi-c*/
/*  //   if ( c1 >= X_dc->ncells || c2 >= X_dc->ncells )  fn_DisplayError(".../mathlib.c/SwapMatricesofCell(): the c1_th or c2_th cell exceeds the cell dimension");   ansi-c*/
/*  //   if ( c1 == c2 )  fn_DisplayError(".../mathlib.c/MatricesofCell(): the two matrices for swapping must be different");   ansi-c*/
/*  //   if ( !X_dc->C[c1]->flag || !X_dc->C[c2]->flag )  fn_DisplayError(".../mathlib.c/MatricesofCell(): both matrices for swapping must have legal values");   ansi-c*/
/*  //   if ( (dim=X_dc->C[c1]->nrows*X_dc->C[c1]->ncols) != (X_dc->C[c2]->nrows*X_dc->C[c2]->ncols) )   ansi-c*/
/*  //      fn_DisplayError(".../mathlib.c/MatricesofCell(): the two matrices for swapping must have the same dimension");   ansi-c*/


/*  //   cblas_dswap(dim, X_dc->C[c1]->M, 1, X_dc->C[c2]->M, 1);   ansi-c*/
/*  //}   ansi-c*/
/*  //#else   ansi-c*/
/*  ////.../mathlib.c/SwapColsofMatrix():  Haven't got time to write my own for-loop routine to swap two columns.  19 Oct. 03   ansi-c*/
/*  //#endif   ansi-c*/


/*  //=== Do NOT know what the following is.  20 Oct. 03.   ansi-c*/
/*  //void SwapPositionsofCell(TSdcell *X_dc, const int c1, const int c2)   ansi-c*/
/*  //{   ansi-c*/
/*  //   //???? Not tested yet.   ansi-c*/
/*  //   int dim;   ansi-c*/


/*  //   if ( !X_dc )  fn_DisplayError(".../mathlib.c/SwapMatricesofCell(): input cell A_dc must be created (memory-allocated)");   ansi-c*/
/*  //   if ( c1 >= X_dc->ncells || c2 >= X_dc->ncells )  fn_DisplayError(".../mathlib.c/SwapMatricesofCell(): the c1_th or c2_th cell exceeds the cell dimension");   ansi-c*/
/*  //   if ( c1 == c2 )  fn_DisplayError(".../mathlib.c/MatricesofCell(): the two matrices for swapping must be different");   ansi-c*/
/*  //   if ( !X_dc->C[c1]->flag || !X_dc->C[c2]->flag )  fn_DisplayError(".../mathlib.c/MatricesofCell(): both matrices for swapping must have legal values");   ansi-c*/
/*  //   if ( (dim=X_dc->C[c1]->nrows*X_dc->C[c1]->ncols) != (X_dc->C[c2]->nrows*X_dc->C[c2]->ncols) )   ansi-c*/
/*  //      fn_DisplayError(".../mathlib.c/MatricesofCell(): the two matrices for swapping must have the same dimension");   ansi-c*/


/*  //   cblas_dswap(dim, X_dc->C[c1]->M, 1, X_dc->C[c2]->M, 1);   ansi-c*/
/*  //}   ansi-c*/


void PermuteColsofMatrix(TSdmatrix *A_dm, const TSivector *indx_iv)
{
/*     //Ouputs:   ansi-c*/
/*     //  A_dm (m-by-n) is replaced by permuted columns only, according to indx_iv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm (m-by-n):  Memory allocated and legal values given already.   ansi-c*/
/*     //  indx_iv (n-by-1): index for columns and rows of A_dm to exchanged simultaneously.  Example: indx_-v->v = {2 0 1} (base 0) for the 3-by-3 matrix   ansi-c*/
/*     //    means that original column 2 is column 0, original column 0 is column 1, etc.   ansi-c*/

   double *B_pd = NULL,
         *A;
   int _j, _n, _m, mn,
       *indx;

   if ( !A_dm || !A_dm->flag )  fn_DisplayError(".../mathlib.c/PermuteColsofMatrix(): input matrix A_dm must (1) be created (memory allocated) and (2) have legal values");
   if ( !indx_iv->flag || (_n=A_dm->ncols) != indx_iv->n )  fn_DisplayError(".../mathlib.c/PermuteColsofMatrix(): (1) sorted index vector, indx_iv, must have legal values; (2) its length must match the number of columns of the input matrix A_dm");


/*     //=== Memory allocated for this function.   ansi-c*/
   B_pd = tzMalloc(mn=_n*(_m=A_dm->nrows), double);

   indx = indx_iv->v;
   memcpy(B_pd, A = A_dm->M, mn*sizeof(double));
   for (_j=_n-1; _j>=0; _j--)
      memcpy(A+_j*_m, B_pd+indx[_j]*_m, _m*sizeof(double));

/*     //=== Destroys memory allocated for this function.   ansi-c*/
   tzDestroy(B_pd);
}

void PermuteRowsofMatrix(TSdmatrix *A_dm, const TSivector *indx_iv)
{
/*     //Ouputs:   ansi-c*/
/*     //  A_dm (n-by-m) is replaced by permuted rows only, according to indx_iv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm (n-by-m):  Memory allocated and legal values given already.   ansi-c*/
/*     //  indx_iv (n-by-1): index for columns and rows of A_dm to exchanged simultaneously.  Example: indx_-v->v = {2 0 1} (base 0) for the 3-by-3 matrix   ansi-c*/
/*     //    means that original column 2 is column 0, original column 0 is column 1, etc.   ansi-c*/

   double *B_pd = NULL,
          *A;
   int _i, _n, _m, mn,
       *indx;
   #if !defined( INTELCMATHLIBRARY )
   int _j;
   #endif

   if ( !A_dm || !A_dm->flag )  fn_DisplayError(".../mathlib.c/PermuteRowsMatrix(): input matrix A_dm must (1) be created (memory allocated) and (2) have legal values");
   if ( !indx_iv->flag || (_n=A_dm->nrows) != indx_iv->n )  fn_DisplayError(".../mathlib.c/PermuteRowsMatrix(): (1) indx_iv must have legal values; (2) number of rows in A_dm must match the length of indx_iv");


/*     //=== Memory allocated for this function.   ansi-c*/
   B_pd = tzMalloc(mn=_n*(_m=A_dm->ncols), double);

   indx = indx_iv->v;
   memcpy(B_pd, A = A_dm->M, mn*sizeof(double));
   #if defined( INTELCMATHLIBRARY )
   for (_i=_n-1; _i>=0; _i--)
      cblas_dcopy(_m, B_pd+indx[_i], _n, A+_i, _n);
   #else   /*  Default to my own routine.   ansi-c*/
   _m = A_dm->ncols;
   for (_j=_m-1; _j>=0; _j--)
     for (_i=_n-1; _i>=0; _i--)
        A[mos(_i, _j, _n)] = B_pd[mos(indx[_i], _j, _n)];
   #endif

/*     //=== Destroys memory allocated for this function.   ansi-c*/
   tzDestroy(B_pd);
}

void PermuteMatrix(TSdmatrix *A_dm, const TSivector *indx_iv)
{
/*     //Ouputs:   ansi-c*/
/*     //  A_dm (n-by-n) is replaced by permuted columns and rows simultaneously.  The permutation is dicated by indx_iv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm (n-by-n):  Memory allocated and legal values given already.   ansi-c*/
/*     //  indx_iv (n-by-1): index for columns and rows of A_dm to exchanged simultaneously.  Example: indx_-v->v = {2 0 1} (base 0) for the 3-by-3 matrix   ansi-c*/
/*     //    means that original column 2 and row 2 are column 0 and row 0, original column 0 and row 0 are column 1 and row 1, etc.   ansi-c*/

   double *B_pd = NULL,
         *A;
   int _i, _j, _n, n2,
       *indx;

   if ( !A_dm || !A_dm->flag )  fn_DisplayError(".../mathlib.c/PermuteMatrix(): input matrix A_dm must (1) be created (memory allocated) and (2) have legal values");
   if ( !indx_iv->flag || (_n=A_dm->nrows) != A_dm->ncols || _n != indx_iv->n )  fn_DisplayError(".../mathlib.c/PermuteMatrix(): (1) indx_iv must have legal values; (2) input matrix A_dm must be square; (3) it dimension must coincide with the length of indx_iv");


/*     //=== Memory allocated for this function.   ansi-c*/
   B_pd = tzMalloc(n2=_n*_n, double);

   indx = indx_iv->v;
   memcpy(B_pd, A = A_dm->M, n2*sizeof(double));
   for (_j=_n-1; _j>=0; _j--)
     for (_i=_n-1; _i>=0; _i--)
        A[mos(_i, _j, _n)] = B_pd[mos(indx[_i], indx[_j], _n)];


/*     //=== Destroys memory allocated for this function.   ansi-c*/
   tzDestroy(B_pd);
}
/*  //=== The following works but may be less efficient and definitely hard to understand.   ansi-c*/
/*  // void PermuteMatrix(TSdmatrix *A_dm, const TSivector *indx_iv)   ansi-c*/
/*  // {   ansi-c*/
/*  //    //Ouputs:   ansi-c*/
/*  //    //  A_dm is replaced by permuted columns and rows simultaneously.  The permutation is dicated by indx_iv.   ansi-c*/
/*  //    //Inputs:   ansi-c*/
/*  //    //  A_dm:  Memory allocated and legal values given already.   ansi-c*/
/*  //    //  indx_iv: index for columns and rows of A_dm to exchanged simultaneously.  Example: indx_-v->v = {2 0 1} (base 0) for the 3-by-3 matrix   ansi-c*/
/*  //    //    means that original column 2 and row 2 are column 0 and row 0, original column 0 and row 0 are column 1 and row 1, etc.   ansi-c*/
/*  //   ansi-c*/
/*  //    double *B_pd = NULL,   ansi-c*/
/*  //          *A;   ansi-c*/
/*  //    int _i, _n, n2,   ansi-c*/
/*  //       *indx;   ansi-c*/
/*  //   ansi-c*/
/*  //    if ( !A_dm || !A_dm->flag )  fn_DisplayError(".../mathlib.c/PermuteMatrix(): input matrix A_dm must (1) be created (memory allocated) and (2) have legal values");   ansi-c*/
/*  //    if ( (_n=A_dm->nrows) != A_dm->ncols || _n != indx_iv->n )  fn_DisplayError(".../mathlib.c/PermuteMatrix(): (1) input matrix A_dm must be square; (2) it dimension must coincide with the length of indx_iv");   ansi-c*/
/*  //   ansi-c*/
/*  //   ansi-c*/
/*  //    //=== Memory allocated for this function.   ansi-c*/
/*  //    B_pd = tzMalloc(n2=_n*_n, double);   ansi-c*/
/*  //   ansi-c*/
/*  //    indx = indx_iv->v;   ansi-c*/
/*  //    memcpy(B_pd, A = A_dm->M, n2*sizeof(double));   ansi-c*/
/*  //    for (_i=0; _i<n2; _i++)  A[_i] = B_pd[indx[_i%_n]+indx[_i/_n]*_n];   ansi-c*/
/*  //   ansi-c*/
/*  //    //=== Destroys memory allocated for this function.   ansi-c*/
/*  //    tzDestroy(B_pd);   ansi-c*/
/*  // }   ansi-c*/


void PermuteMatricesofCell(TSdcell *A_dc, const TSivector *indx_iv)
{
/*     //Ouputs:   ansi-c*/
/*     //  A_dc is replaced by permuted matrices.  The permutation is dicated by indx_iv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dc:  Memory allocated and legal values given already.   ansi-c*/
/*     //  indx_iv: index for matrices of A_dc to exchanged simultaneously.  Example: indx_-v->v = {2 0 1} (base 0) for the 3-by-1 cell   ansi-c*/
/*     //    means that original matrix 2 is matrix 0, original matrix 0 is matrix 1, etc.   ansi-c*/

   int _i, _n,
       *indx_p;
   TSdmatrix **tA_p2dm = NULL;

   if ( !A_dc || !indx_iv || !indx_iv->flag )  fn_DisplayError(".../mathlib.c/PermuteMatricesofCell(): (1) input cell A_dc must be created (memory-allocated); (2) index vector indx_iv must be created and have legal values");
   if ( (_n=A_dc->ncells) != indx_iv->n )  fn_DisplayError(".../mathlib.c/PermuteMatricesofCell(): number of cells must match the length of indx_iv");
   indx_p = indx_iv->v;


/*     //=== Memory allocated for this function.   ansi-c*/
   tA_p2dm = tzMalloc(_n, TSdmatrix *);

   for (_i=_n-1; _i>=0; _i--)  tA_p2dm[_i] = A_dc->C[indx_p[_i]];
/*     //=== This one is less efficient than the following: for (_i=_n-1; _i>=0; _i--)  A_dc->C[_i] = tA_p2dm[_i];   ansi-c*/
   memcpy(A_dc->C, tA_p2dm, _n*sizeof(TSdmatrix *));

/*     //=== Destroys memory allocated for this function.   ansi-c*/
   tzDestroy(tA_p2dm);
}


void ScalarTimesColofMatrix(TSdvector *y_dv, double _alpha, TSdmatrix *X_dm, int _j)
{
/*     //????????? Default option, in the #else, has NOT been tested yet!   ansi-c*/
/*     //Ouputs:   ansi-c*/
/*     //  If y_dv!=NULL, y_dv is the jth column of X_dm is multiplied by _alpha.   ansi-c*/
/*     //  If !y_dv, the jth column of X_dm is replaced by the new value, which will be multiplied by _alpha.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha:  Scalar.   ansi-c*/
/*     //  X_dm:  Memory allocated and legal values given already.   ansi-c*/
/*     //  _j: The jth column of X_dm.   ansi-c*/

   #if !defined( INTELCMATHLIBRARY )
      int _i;
   #endif
   int nrows;
   double *M, *v;

   if ( !X_dm || !X_dm->flag ) fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix(): (1) Input matrix must be created (memory-allocated); (2) Legal values must be given");
   if (_j >= X_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix(): The jth column specified for the input matrix exceeds the column dimension");

   #if defined( INTELCMATHLIBRARY )
      M = X_dm->M + _j*(nrows=X_dm->nrows);   /*  Points to the beginning of the jth column.   ansi-c*/
      if (!y_dv)  cblas_dscal(nrows, _alpha, M, 1);
      else {
         memcpy(v=y_dv->v, M, nrows*sizeof(double));
         cblas_dscal(nrows, _alpha, v, 1);
         y_dv->flag = V_DEF;
      }
   #else
/*        //      Need to be tested for the following.   ansi-c*/
/*        //   ansi-c*/
/*        // M = X_dm->M + (_j+1)*(nrows=X_dm->nrows) - 1;  //Points to the end of the jth column.   ansi-c*/
/*        // if (!y_dv)   ansi-c*/
/*        //    for (_i=nrows-1; _i>=0; _i--, M--)  *M = _alpha * (*M);   ansi-c*/
/*        // else {   ansi-c*/
/*        //    v = y_dv->v;   ansi-c*/
/*        //    for (_i=nrows-1; _i>=0; _i--, M--)  v[_i] = _alpha * (*M);   ansi-c*/
/*        //    y_dv->flag = V_DEF;   ansi-c*/
/*        // }   ansi-c*/
   #endif
}
void ScalarTimesColofMatrix2ColofMatrix(TSdmatrix *Y_dm, int jy, double _alpha, TSdmatrix *X_dm, int jx)
{
/*     //Ouputs:   ansi-c*/
/*     //  If Y_dm!=NULL, the jy_th column of Y_dm is the jx_th column of X_dm multiplied by _alpha.   ansi-c*/
/*     //  If !Y_dm, the jx_th column of X_dm is replaced by the new value, which will be multiplied by _alpha.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha:  Scalar.   ansi-c*/
/*     //  X_dm:  Memory allocated and legal values given already.   ansi-c*/
/*     //  jy: The jy_th column of Y_dm.   ansi-c*/
/*     //  yx: The jx_th column of X_dm.   ansi-c*/

   #if !defined( INTELCMATHLIBRARY )
      int _i;
   #endif
   int nrows_x;
   double *Mx, *My;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): (1) Input matrix must be created (memory-allocated); (2) Legal values must be given");
   if (jx >= X_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The jth column specified for the input matrix exceeds the column dimension");

   #if defined( INTELCMATHLIBRARY )
      Mx = X_dm->M + jx*(nrows_x=X_dm->nrows);   /*  Points to the beginning of the jth column.   ansi-c*/
      if (!Y_dm)  cblas_dscal(nrows_x, _alpha, Mx, 1);
      else {
         if (jy >= Y_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The jth column specified for the output matrix exceeds the column dimension");
         if ( nrows_x != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The number of rows for both input and output matrices must be the same");

         My = Y_dm->M + jy*nrows_x;   /*  Points to the beginning of the jth column.   ansi-c*/
         memcpy(My, Mx, nrows_x*sizeof(double));
         cblas_dscal(nrows_x, _alpha, My, 1);
      }
   #else
      Mx = X_dm->M + (jx+1)*(nrows_x=X_dm->nrows) - 1;   /*  Points to the end of the jth column.   ansi-c*/
      if (!Y_dm)
         for (_i=nrows_x-1; _i>=0; _i--, Mx--)  *Mx = _alpha * (*Mx);
      else {
         if (jy >= Y_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The jth column specified for the output matrix exceeds the column dimension");
         if ( nrows_x != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The number of rows for both input and output matrices must be the same");

         My = Y_dm->M + (jy+1)*nrows_x - 1;   /*  Points to the end of the jth column.   ansi-c*/
         for (_i=nrows_x-1; _i>=0; _i--, Mx--, My--)  *My = _alpha * (*Mx);
      }
   #endif
}


void ScalarTimesColofMatrixPlusVector2ColofMatrix(TSdmatrix *Y_dm, int jy, double _alpha, TSdmatrix *X_dm, int jx, double _beta, TSdvector *x_dv) {
/*     //Ouputs:   ansi-c*/
/*     //  If Y_dm!=NULL, Y(:,jy) = alpha*X(:,jx) + beta*x.   ansi-c*/
/*     //  If !Y_dm, X(:,jx) = alpha*X(:,jx) + beta*x.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha:  Scalar.   ansi-c*/
/*     //  _beta: Scalar.   ansi-c*/
/*     //  X_dm:  Memory allocated and legal values given already.   ansi-c*/
/*     //  jy: The jy_th column of Y_dm.   ansi-c*/
/*     //  yx: The jx_th column of X_dm.   ansi-c*/

   int _i, nrows_x;
   double *Mx, *My, *v;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrixPlusVector2ColofMatrix(): (1) Input matrix must be created (memory-allocated); (2) Legal values must be given");
   if (jx >= X_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrixPlusVector2ColofMatrix(): The jth column specified for the input matrix exceeds the column dimension");
   if ( !x_dv || !x_dv->flag )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrixPlusVector2ColofMatrix(): (1) input vectr must be created (memory-allocated); (2) legal values must be given");

   if (_beta == 0.0) {
      #if defined( INTELCMATHLIBRARY )
         Mx = X_dm->M + jx*(nrows_x=X_dm->nrows);   /*  Points to the beginning of the jth column.   ansi-c*/
         if (!Y_dm)  cblas_dscal(nrows_x, _alpha, Mx, 1);
         else {
            if (jy >= Y_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The jth column specified for the output matrix exceeds the column dimension");
            if ( nrows_x != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The number of rows for both input and output matrices must be the same");

            My = Y_dm->M + jy*nrows_x;   /*  Points to the beginning of the jth column.   ansi-c*/
            memcpy(My, Mx, nrows_x*sizeof(double));
            cblas_dscal(nrows_x, _alpha, My, 1);
         }
      #else
         Mx = X_dm->M + (jx+1)*(nrows_x=X_dm->nrows) - 1;   /*  Points to the end of the jth column.   ansi-c*/
         if (!Y_dm)
            for (_i=nrows_x-1; _i>=0; _i--, Mx--)  *Mx = _alpha * (*Mx);
         else {
            if (jy >= Y_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The jth column specified for the output matrix exceeds the column dimension");
            if ( nrows_x != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrix2ColofMatrix(): The number of rows for both input and output matrices must be the same");

            My = Y_dm->M + (jy+1)*nrows_x - 1;   /*  Points to the end of the jth column.   ansi-c*/
            for (_i=nrows_x-1; _i>=0; _i--, Mx--, My--)  *My = _alpha * (*Mx);
         }
      #endif
   }
   else {
      Mx = X_dm->M + (jx+1)*(nrows_x=X_dm->nrows) - 1;   /*  Points to the end of the jth column.   ansi-c*/
      if ( nrows_x != x_dv->n )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrixPlusVector2ColofMatrix(): The length of the input vector must match the number of rows in the input matrix");

      if (!Y_dm) {
         if ( _alpha == 1.0 && _beta == 1.0 )
            for (_i=nrows_x-1, v=x_dv->v+_i; _i>=0; _i--, Mx--, v--)  *Mx = *Mx + *v;
         else if ( _alpha == 1.0 )
            for (_i=nrows_x-1, v=x_dv->v+_i; _i>=0; _i--, Mx--, v--)  *Mx = _alpha * (*Mx) + *v;
         else
            for (_i=nrows_x-1, v=x_dv->v+_i; _i>=0; _i--, Mx--, v--)  *Mx = _alpha * (*Mx) + _beta * (*v);
      }
      else {
         if (jy >= Y_dm->ncols)  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrixPlusVector2ColofMatrix(): The jth column specified for the output matrix exceeds the column dimension");
         if ( nrows_x != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ScalarTimesColofMatrixPlusVector2ColofMatrix(): The number of rows for both input and output matrices must be the same");

         My = Y_dm->M + (jy+1)*nrows_x - 1;   /*  Points to the end of the jth column.   ansi-c*/
         if ( _alpha == 1.0 && _beta == 1.0 )
            for (_i=nrows_x-1, v=x_dv->v+_i; _i>=0; _i--, Mx--, My--, v--)  *My = *Mx + *v;
         else if ( _alpha == 1.0 )
            for (_i=nrows_x-1, v=x_dv->v+_i; _i>=0; _i--, Mx--, My--, v--)  *My = _alpha * (*Mx) + *v;
         else
            for (_i=nrows_x-1, v=x_dv->v+_i; _i>=0; _i--, Mx--, My--, v--)  *My = _alpha * (*Mx) + _beta * (*v);
      }
   }
}


void MatrixDotDivideVector_row(TSdmatrix *Y_dm, TSdmatrix *X_dm, TSdvector *x_dv, double _alpha, double _beta)
{
/*     //Outputs:   ansi-c*/
/*     //  If (Y_dm != X_dm), Y_dm(ix, :) = _alpha * X_dm(ix, :) ./ x_dv + _beta * X_dm(ix, :), for all ix.   ansi-c*/
/*     //  If (Y_dm = X_dm), X_dm(ix, :) = _alpha * X_dm(ix, :) ./ x_dv + _beta * X_dm(ix, :), for all ix.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
/*     //  _beta: double scalar.   ansi-c*/

   int _i, cnt, ncols_x, nrows_x, nrows_xm1, ix;
   double *X, *x, *Y, xinv;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/MatrixDotDivideVector(): (1) Input matrix must be created (memory-allocated); (2) Legal values must be given");
   if ( !x_dv || !x_dv->flag )  fn_DisplayError(".../mathlib.c/MatrixDotDivideVector(): (1) Input vector must be created (memory-allocated); (2) Legal values must be given");
   if ( (ncols_x=X_dm->ncols) != x_dv->n )  fn_DisplayError(".../mathlib.c/MatrixDotDivideVector(): Number of columns in the input matrix must match the length of the input vector");
   X = X_dm->M;
   x = x_dv->v;
   nrows_xm1 = (nrows_x=X_dm->nrows)-1;


   if ( _beta==0.0 ) {
      if ( Y_dm == X_dm ) {
         for (ix=nrows_x*ncols_x-1, cnt=ncols_x-1; ix>=nrows_xm1; ix -= nrows_x, cnt--) {
/*              //Last row of X_dm   ansi-c*/
            xinv = _alpha/x[cnt];
            for (_i=ix-nrows_x+1; _i<=ix;  _i++)   X[_i] *= xinv;   /*  Must _i<=ix, not _i<ix.  For each column at time.   ansi-c*/
         }
      }
      else {
         Y = Y_dm->M;
         if ( ncols_x != Y_dm->ncols || nrows_x != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/MatrixDotDivideVector(): Dimension of output matrix Y_dm must be the same as that of input matrix X_dm");
         for (ix=nrows_x*ncols_x-1, cnt=ncols_x-1; ix>=nrows_xm1; ix -= nrows_x, cnt--) {
/*              //Last row of X_dm   ansi-c*/
            xinv = _alpha/x[cnt];
            for (_i=ix-nrows_x+1, cnt=0; _i<=ix;  _i++, cnt++)  Y[_i] = X[_i] * xinv;     /*  Must _i<=ix, not _i<ix.  For each column at time.   ansi-c*/
         }
         Y_dm->flag = M_GE;
      }
   }
   else {
      if ( Y_dm == X_dm ) {
         for (ix=nrows_x*ncols_x-1, cnt=ncols_x-1; ix>=nrows_xm1; ix -= nrows_x, cnt--) {
/*              //Last row of X_dm   ansi-c*/
            xinv = _alpha/x[cnt] + _beta;
            for (_i=ix-nrows_x+1; _i<=ix;  _i++)   X[_i] *= xinv;   /*  Must _i<=ix, not _i<ix.  For each column at time.   ansi-c*/
         }
      }
      else {
         Y = Y_dm->M;
         if ( ncols_x != Y_dm->ncols || nrows_x != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/MatrixDotDivideVector(): Dimension of output matrix Y_dm must be the same as that of input matrix X_dm");
         for (ix=nrows_x*ncols_x-1, cnt=ncols_x-1; ix>=nrows_xm1; ix -= nrows_x, cnt--) {
/*              //Last row of X_dm   ansi-c*/
            xinv = _alpha/x[cnt] + _beta;
            for (_i=ix-nrows_x+1, cnt=0; _i<=ix;  _i++, cnt++)  Y[_i] = X[_i] * xinv;     /*  Must _i<=ix, not _i<ix.  For each column at time.   ansi-c*/
         }
         Y_dm->flag = M_GE;
      }
   }
}


void RowofMatrixDotDivideVector(TSdvector *y_dv, TSdmatrix *X_dm, int ix, TSdvector *x_dv, double _alpha, double _beta)
{
/*     //??????? NOT tested yet, 01/02/04.   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If (y_dv), y_dv = _alpha * X_dm(ix, :) ./ x_dv + _beta * X_dm(ix, :).   ansi-c*/
/*     //  If (!y_dv), X_dm(ix, :) = _alpha * X_dm(ix, :) ./ x_dv + _beta * X_dm(ix, :).   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
/*     //  _beta: double scalar.   ansi-c*/

   int _i, cnt, ncols_x, nrows_x;
   double *X, *x, *y;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): (1) Input matrix must be created (memory-allocated); (2) Legal values must be given");
   if ( !x_dv || !x_dv->flag )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): (1) Input vector must be created (memory-allocated); (2) Legal values must be given");
   if ( (ncols_x=X_dm->ncols) != x_dv->n )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): Number of columns in the input matrix must match the length of the input vector");
   if ( ix >= (nrows_x=X_dm->nrows) )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): The specified ix_th row of the input matrix exceeds its row dimension");
   X = X_dm->M;
   x = x_dv->v;


   if ( _alpha==1.0 ) {
      if ( _beta==0.0 ) {
         if ( !y_dv )
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  X[_i] /= x[cnt];
         else {
            y = y_dv->v;
            if ( ncols_x != y_dv->n )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): Number of columns in the input matrix must match the length of the output vector");
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  y[cnt] = X[_i] / x[cnt];
            y_dv->flag = V_DEF;
         }
      }
/*  //      else if ( _beta==1.0 ) {   ansi-c*/
/*  //         if ( !y_dv )   ansi-c*/
/*  //            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  X[_i] *= ( 1.0 / x[cnt] + 1.0);   ansi-c*/
/*  //         else {   ansi-c*/
/*  //            y = y_dv->v;   ansi-c*/
/*  //            if ( ncols_x != y_dv->n )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): Number of columns in the input matrix must match the length of the output vector");   ansi-c*/
/*  //            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  y[cnt] = X[_i] * ( 1.0 / x[cnt] + 1.0 );   ansi-c*/
/*  //            y_dv->flag = V_DEF;   ansi-c*/
/*  //         }   ansi-c*/
/*  //      }   ansi-c*/
      else {
         if ( !y_dv )
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  X[_i] *= 1.0 / x[cnt] + _beta;
         else {
            y = y_dv->v;
            if ( ncols_x != y_dv->n )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): Number of columns in the input matrix must match the length of the output vector");
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  y[cnt] = X[_i] * ( 1.0 / x[cnt] + _beta );
            y_dv->flag = V_DEF;
         }
      }
   }
   else {
      if ( _beta==0.0 ) {
         if ( !y_dv )
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  X[_i] *= _alpha / x[cnt];
         else {
            y = y_dv->v;
            if ( ncols_x != y_dv->n )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): Number of columns in the input matrix must match the length of the output vector");
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  y[cnt] = (_alpha * X[_i]) / x[cnt];
            y_dv->flag = V_DEF;
         }
      }
/*  //      else if ( _beta==1.0 ) {   ansi-c*/
/*  //         if ( !y_dv )   ansi-c*/
/*  //            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  X[_i] *= _alpha / x[cnt] + 1.0;   ansi-c*/
/*  //         else {   ansi-c*/
/*  //            y = y_dv->v;   ansi-c*/
/*  //            if ( ncols_x != y_dv->n )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): Number of columns in the input matrix must match the length of the output vector");   ansi-c*/
/*  //            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  y[cnt] = X[_i] * (_alpha / x[cnt] + 1.0);   ansi-c*/
/*  //            y_dv->flag = V_DEF;   ansi-c*/
/*  //         }   ansi-c*/
/*  //      }   ansi-c*/
      else {
         if ( !y_dv )
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  X[_i] *= _alpha / x[cnt] + _beta;
         else {
            y = y_dv->v;
            if ( ncols_x != y_dv->n )  fn_DisplayError(".../mathlib.c/RowofMatrixDotDivideVector(): Number of columns in the input matrix must match the length of the output vector");
            for (_i=ix + (ncols_x-1)*nrows_x, cnt=ncols_x-1; _i>=ix; _i -= nrows_x, cnt--)  y[cnt] = X[_i] * (_alpha / x[cnt] + _beta);
            y_dv->flag = V_DEF;
         }
      }
   }
}

void ColofMatrixDotTimesVector(TSdvector *y_dv, TSdmatrix *X_dm, int jx, TSdvector *x_dv, double _alpha, double _beta) {
/*     //Outputs:   ansi-c*/
/*     //  If (y_dv), y_dv = _alpha * X_dm(:,jx) .* x_dv + _beta * X_dm(:,jx).   ansi-c*/
/*     //  If (!y_dv), X_dm(:,jx) = _alpha * X_dm(:,jx) .* x_dv + _beta * X_dm(:,jx).   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
/*     //  _beta: double scalar.   ansi-c*/

   int _i, nrows_x;
   double *X, *x, *y;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): (1) Input matrix must be created (memory-allocated); (2) Legal values must be given");
   if ( !x_dv || !x_dv->flag )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): (1) Input vector must be created (memory-allocated); (2) Legal values must be given");
   if ( (nrows_x=X_dm->nrows) != x_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the input vector");
   if ( jx >= X_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jx_th column of the input matrix exceeds its column dimension");


   X = X_dm->M + (jx+1)*nrows_x - 1;   /*  Points to the end of the jx_th column.   ansi-c*/
   x = x_dv->v + nrows_x - 1;   /*  Points to the end of the vector.   ansi-c*/
   if ( _alpha==1.0 ) {
      if ( _beta==0.0 ) {
         if ( !y_dv )
            for (_i=nrows_x-1; _i>=0; _i--, X--, x--)  *X = (*X) * (*x);
         else {
            if ( nrows_x != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows_x-1, y=y_dv->v+_i; _i>=0; _i--, X--, x--, y--)  *y = (*X) * (*x);
            y_dv->flag = V_DEF;
         }
      }
      else if ( _beta==1.0 ) {
         if ( !y_dv )
            for (_i=nrows_x-1; _i>=0; _i--, X--, x--)  *X = (*X) * (*x) + (*X);
         else {
            if ( nrows_x != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows_x-1, y=y_dv->v+_i; _i>=0; _i--, X--, x--, y--)  *y = (*X) * (*x) + (*X);
            y_dv->flag = V_DEF;
         }
      }
      else {
         if ( !y_dv )
            for (_i=nrows_x-1; _i>=0; _i--, X--, x--)  *X = (*X) * (*x) + _beta * (*X);
         else {
            if ( nrows_x != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows_x-1, y=y_dv->v+_i; _i>=0; _i--, X--, x--, y--)  *y = (*X) * (*x) + _beta * (*X);
            y_dv->flag = V_DEF;
         }
      }
   }
   else {
      if ( _beta==0.0 ) {
         if ( !y_dv )
            for (_i=nrows_x-1; _i>=0; _i--, X--, x--)  *X = _alpha * (*X) * (*x);
         else {
            if ( nrows_x != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows_x-1, y=y_dv->v+_i; _i>=0; _i--, X--, x--, y--)  *y = _alpha * (*X) * (*x);
            y_dv->flag = V_DEF;
         }
      }
      else if ( _beta==1.0 ) {
         if ( !y_dv )
            for (_i=nrows_x-1; _i>=0; _i--, X--, x--)  *X = _alpha * (*X) * (*x) + (*X);
         else {
            if ( nrows_x != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows_x-1, y=y_dv->v+_i; _i>=0; _i--, X--, x--, y--)  *y = _alpha * (*X) * (*x) + (*X);
            y_dv->flag = V_DEF;
         }
      }
      else {
         if ( !y_dv )
            for (_i=nrows_x-1; _i>=0; _i--, X--, x--)  *X = _alpha * (*X) * (*x) + _beta * (*X);
         else {
            if ( nrows_x != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows_x-1, y=y_dv->v+_i; _i>=0; _i--, X--, x--, y--)  *y = _alpha * (*X) * (*x) + _beta * (*X);
            y_dv->flag = V_DEF;
         }
      }
   }
}
void ColofMatrixDotTimesColofMatrix(TSdvector *y_dv, TSdmatrix *X1_dm, int jx1, TSdmatrix *X2_dm, int jx2, double _alpha, double _beta) {
/*     //????????? NOT tested yet.   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If y_dv!=NULL, y_dv = _alpha * X1_dm(:,jx1) .* X2_dm(:,jx2) + _beta * X1_dm(:,jx1).   ansi-c*/
/*     //  If !y_dv, X1_dm(:,jx1) = _alpha * X1_dm(:,jx1) .* X2_dm(:,jx2) + _beta * X1_dm(:,jx2).   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
/*     //  _beta: double scalar.   ansi-c*/

   int _i, nrows1;
   double *X1, *X2, *y;

   if ( !X1_dm || !X1_dm->flag )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesColofMatrix(): (1) Input matrix X1 must be created (memory-allocated); (2) Legal values must be given");
   if ( !X2_dm || !X2_dm->flag )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesColofMatrix(): (1) Input matrix X2 must be created (memory-allocated); (2) Legal values must be given");
   if ( (nrows1=X1_dm->nrows) != X2_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesColofMatrix(): Numbers of rows in both input matrices must be the same");
   if ( jx1 >= X1_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jx1_th column of input matrix X1 exceeds its column dimension");
   if ( jx2 >= X2_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jx2_th column of input matrix X2 exceeds its column dimension");

   X1 = X1_dm->M + (jx1+1)*nrows1 - 1;   /*  Points to the end of the jx1_th column.   ansi-c*/
   X2 = X2_dm->M + (jx2+1)*nrows1 - 1;   /*  Points to the end of the jx2_th column.   ansi-c*/
   if ( _alpha==1.0 ) {
      if ( _beta==0.0 ) {
         if ( !y_dv )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = (*X1) * (*X2);
         else {
            if ( nrows1 != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows1-1, y=y_dv->v+_i; _i>=0; _i--, X1--, X2--, y--)  *y = (*X1) * (*X2);
            y_dv->flag = V_DEF;
         }
      }
      else if ( _beta==1.0 ) {
         if ( !y_dv )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = (*X1) * (*X2) + (*X1);
         else {
            if ( nrows1 != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows1-1, y=y_dv->v+_i; _i>=0; _i--, X1--, X2--, y--)  *y = (*X1) * (*X2) + (*X1);
            y_dv->flag = V_DEF;
         }
      }
      else {
         if ( !y_dv )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = (*X1) * (*X2) + _beta * (*X1);
         else {
            if ( nrows1 != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows1-1, y=y_dv->v+_i; _i>=0; _i--, X1--, X2--, y--)  *y = (*X1) * (*X2) + _beta * (*X1);
            y_dv->flag = V_DEF;
         }
      }
   }
   else {
      if ( _beta==0.0 ) {
         if ( !y_dv )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = _alpha * (*X1) * (*X2);
         else {
            if ( nrows1 != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows1-1, y=y_dv->v+_i; _i>=0; _i--, X1--, X2--, y--)  *y = _alpha * (*X1) * (*X2);
            y_dv->flag = V_DEF;
         }
      }
      else if ( _beta==1.0 ) {
         if ( !y_dv )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = _alpha * (*X1) * (*X2) + (*X1);
         else {
            if ( nrows1 != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows1-1, y=y_dv->v+_i; _i>=0; _i--, X1--, X2--, y--)  *y = _alpha * (*X1) * (*X2) + (*X1);
            y_dv->flag = V_DEF;
         }
      }
      else {
         if ( !y_dv )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = _alpha * (*X1) * (*X2) + _beta * (*X1);
         else {
            if ( nrows1 != y_dv->n )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in the input matrix must match the length of the output vector");
            for (_i=nrows1-1, y=y_dv->v+_i; _i>=0; _i--, X1--, X2--, y--)  *y = _alpha * (*X1) * (*X2) + _beta * (*X1);
            y_dv->flag = V_DEF;
         }
      }
   }
}
void ColofMatrixDotTimesColofMatrix2ColofMatrix(TSdmatrix *Y_dm, int jy, TSdmatrix *X1_dm, int jx1, TSdmatrix *X2_dm, int jx2, double _alpha, double _beta) {
/*     //Outputs:   ansi-c*/
/*     //  If Y_dm!=NULL, Y_dm(:,jy) = _alpha * X1_dm(:,jx1) .* X2_dm(:,jx2) + _beta * X1_dm(:,jx1).   ansi-c*/
/*     //  If !Y_dm, X1_dm(:,jx1) = _alpha * X1_dm(:,jx1) .* X2_dm(:,jx2) + _beta * X1_dm(:,jx2).   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
/*     //  _beta: double scalar.   ansi-c*/

   int _i, nrows1;
   double *X1, *X2, *Y;

   if ( !X1_dm || !X1_dm->flag )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesColofMatrix(): (1) Input matrix X1 must be created (memory-allocated); (2) Legal values must be given");
   if ( !X2_dm || !X2_dm->flag )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesColofMatrix(): (1) Input matrix X2 must be created (memory-allocated); (2) Legal values must be given");
   if ( (nrows1=X1_dm->nrows) != X2_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesColofMatrix(): Numbers of rows in both input matrices must be the same");
   if ( jx1 >= X1_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jx1_th column of input matrix X1 exceeds its column dimension");
   if ( jx2 >= X2_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jx2_th column of input matrix X2 exceeds its column dimension");

   X1 = X1_dm->M + (jx1+1)*nrows1 - 1;   /*  Points to the end of the jx1_th column.   ansi-c*/
   X2 = X2_dm->M + (jx2+1)*nrows1 - 1;   /*  Points to the end of the jx2_th column.   ansi-c*/
   if ( _alpha==1.0 ) {
      if ( _beta==0.0 ) {
         if ( !Y_dm )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = (*X1) * (*X2);
         else {
            if ( nrows1 != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in input matrices must match that of the output matrix");
            if ( jy >= Y_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jy_th column of output matrix Y exceeds its column dimension");
            Y = Y_dm->M + (jy+1)*nrows1 - 1;   /*  Points to the end of the jy_th column.   ansi-c*/

            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--, Y--)  *Y = (*X1) * (*X2);
         }
      }
      else if ( _beta==1.0 ) {
         if ( !Y_dm )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = (*X1) * (*X2) + (*X1);
         else {
            if ( nrows1 != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in input matrices must match that of the output matrix");
            if ( jy >= Y_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jy_th column of output matrix Y exceeds its column dimension");
            Y = Y_dm->M + (jy+1)*nrows1 - 1;   /*  Points to the end of the jy_th column.   ansi-c*/

            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--, Y--)  *Y = (*X1) * (*X2) + (*X1);
         }
      }
      else {
         if ( !Y_dm )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = (*X1) * (*X2) + _beta * (*X1);
         else {
            if ( nrows1 != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in input matrices must match that of the output matrix");
            if ( jy >= Y_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jy_th column of output matrix Y exceeds its column dimension");
            Y = Y_dm->M + (jy+1)*nrows1 - 1;   /*  Points to the end of the jy_th column.   ansi-c*/

            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--, Y--)  *Y = (*X1) * (*X2) + _beta * (*X1);
         }
      }
   }
   else {
      if ( _beta==0.0 ) {
         if ( !Y_dm )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = _alpha * (*X1) * (*X2);
         else {
            if ( nrows1 != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in input matrices must match that of the output matrix");
            if ( jy >= Y_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jy_th column of output matrix Y exceeds its column dimension");
            Y = Y_dm->M + (jy+1)*nrows1 - 1;   /*  Points to the end of the jy_th column.   ansi-c*/

            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--, Y--)  *Y = _alpha * (*X1) * (*X2);
         }
      }
      else if ( _beta==1.0 ) {
         if ( !Y_dm )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = _alpha * (*X1) * (*X2) + (*X1);
         else {
            if ( nrows1 != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in input matrices must match that of the output matrix");
            if ( jy >= Y_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jy_th column of output matrix Y exceeds its column dimension");
            Y = Y_dm->M + (jy+1)*nrows1 - 1;   /*  Points to the end of the jy_th column.   ansi-c*/

            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--, Y--)  *Y = _alpha * (*X1) * (*X2) + (*X1);
         }
      }
      else {
         if ( !Y_dm )
            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--)  *X1 = _alpha * (*X1) * (*X2) + _beta * (*X1);
         else {
            if ( nrows1 != Y_dm->nrows )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): Number of rows in input matrices must match that of the output matrix");
            if ( jy >= Y_dm->ncols )  fn_DisplayError(".../mathlib.c/ColofMatrixDotTimesVector(): The specified jy_th column of output matrix Y exceeds its column dimension");
            Y = Y_dm->M + (jy+1)*nrows1 - 1;   /*  Points to the end of the jy_th column.   ansi-c*/

            for (_i=nrows1-1; _i>=0; _i--, X1--, X2--, Y--)  *Y = _alpha * (*X1) * (*X2) + _beta * (*X1);
         }
      }
   }
}


void MatrixPlusMatrixUpdate(TSdmatrix *X_dm, TSdmatrix *A_dm) {
/*     //Output: X = A + X where X_dm is an m-by-n general (and possibly symmetric) matrix.   ansi-c*/
/*     //  If X = A, then X will be replaced by 2*A;   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n general or symmetric matrix.   ansi-c*/
   int _i, _m, _n, nels;
   double *X, *A;


   if ( !X_dm || !A_dm )  fn_DisplayError(".../mathlib.c/MatrixPlusMatrixUpdate(): All input matrices must be created (memory-allocated)");
   else if ( !X_dm->flag || !A_dm->flag )   fn_DisplayError(".../mathlib.c/MatrixPlusMatrixUpdate(): Both input matrices must be given values");
   else {
      _m = X_dm->nrows;
      _n = X_dm->ncols;
      nels = _m * _n;
      X = X_dm->M;
      A = A_dm->M;
   }

   if ( (_m != A_dm->nrows) || (_n != A_dm->ncols) )  fn_DisplayError(".../mathlib.c/MatrixPlusMatrixUpdate(): Dimensions of all input matrices must be same");

/*     //=== Making both X_dm and A_dm general if not yet.   ansi-c*/
   if ( !(X_dm->flag & M_GE) ) {
      if (X_dm->flag & M_SU)   SUtoGE(X_dm);
      else if (X_dm->flag & M_SL)   SLtoGE(X_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixPlusMatrixUpdate(): Haven't got time to deal with the M_UT and M_LT cases for X_dm");
   }
   if ( !(A_dm->flag & M_GE) ) {
      if (A_dm->flag & M_SU)   SUtoGE(A_dm);
      else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
      else  fn_DisplayError(".../mathlib.c/MatrixPlusMatrixUpdate(): Haven't got time to deal with the M_UT and M_LT cases for A_dm");
   }
   for (_i=nels-1; _i>=0; _i--)   X[_i] += A[_i];     /*  This operation may be much cheaper than explicitly using SU or SL operations with two for loops and integer multiplications for matrix offsets.   ansi-c*/

   if ( X_dm->flag != A_dm->flag )  X_dm->flag = M_GE;   /*  Reset to a general matrix only; otherwise, keep the original X_dm->flag.   ansi-c*/
}

void MatrixPlusMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm) {
/*     //Output: X = A + B where X_dm is an m-by-n general matrix.   ansi-c*/
/*     //  If X=A, A will be replaced by X; if X=B, B will be replaced by X.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n general matrix.   ansi-c*/
/*     //  B_dm: m-by-n general matrix.   ansi-c*/
   int _i, _m, _n, nels;
   double *X, *A, *B;


   if ( !X_dm || !A_dm || !B_dm )  fn_DisplayError(".../mathlib.c/MatrixPlusMatrix(): All input matrices must be created (memory-allocated)");
   else if ( !A_dm->flag || !B_dm->flag )   fn_DisplayError(".../mathlib.c/MatrixPlusMatrix(): Two R input matrices must be given values");
   else {
      _m = X_dm->nrows;
      _n = X_dm->ncols;
      nels = _m * _n;
      X = X_dm->M;
      A = A_dm->M;
      B = B_dm->M;
   }


   if ( (_m != A_dm->nrows) || (_m != B_dm->nrows) || (_n != A_dm->ncols) || (_n != B_dm->ncols) )
      fn_DisplayError(".../mathlib.c/MatrixPlusMatrix(): Dimensions of all input matrices must be same");
   else {
      if ( !(A_dm->flag & M_GE) ) {
         if (A_dm->flag & M_SU)   SUtoGE(A_dm);
         else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
         else  fn_DisplayError(".../mathlib.c/MatrixPlusMatrix(): Haven't got time to deal with the M_UT and M_LT cases for A_dm");
      }
      if ( !(B_dm->flag & M_GE) ) {
         if (B_dm->flag & M_SU) SUtoGE(B_dm);
         else if (B_dm->flag & M_SL) SLtoGE(B_dm);
         else  fn_DisplayError(".../mathlib.c/MatrixPlusMatrix(): Haven't got time to deal with the M_UT and M_LT cases for B_dm");
      }

      for (_i=nels-1; _i>=0; _i--)   X[_i] = A[_i] + B[_i];
      if (A_dm->flag == B_dm->flag)  X_dm->flag = A_dm->flag;
      else  X_dm->flag = M_GE;
   }
}

void MatrixMinusMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm) {
/*     //Output: X = A - B where X_dm is an m-by-n general matrix.   ansi-c*/
/*     //  If X=A, A will be replaced by X; if X=B, B will be replaced by X.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n general matrix.   ansi-c*/
/*     //  B_dm: m-by-n general matrix.   ansi-c*/
   int _i, _m, _n, nels;
   double *X, *A, *B;


   if ( !X_dm || !A_dm || !B_dm )  fn_DisplayError(".../mathlib.c/MatrixMinusMatrix(): All input matrices must be created (memory-allocated)");
   else if ( !A_dm->flag || !B_dm->flag )   fn_DisplayError(".../mathlib.c/MatrixMinusMatrix(): Two R input matrices must be given values");
   else {
      _m = X_dm->nrows;
      _n = X_dm->ncols;
      nels = _m * _n;
      X = X_dm->M;
      A = A_dm->M;
      B = B_dm->M;
   }


   if ( (_m != A_dm->nrows) || (_m != B_dm->nrows) || (_n != A_dm->ncols) || (_n != B_dm->ncols) )
      fn_DisplayError(".../mathlib.c/MatrixMinusMatrix(): Dimensions of all input matrices must be same");
   else {
      if ( !(A_dm->flag & M_GE) ) {
         if (A_dm->flag & M_SU)   SUtoGE(A_dm);
         else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
         else  fn_DisplayError(".../mathlib.c/MatrixMinusMatrix(): Haven't got time to deal with the M_UT and M_LT cases for A_dm");
      }
      if ( !(B_dm->flag & M_GE) ) {
         if (B_dm->flag & M_SU) SUtoGE(B_dm);
         else if (B_dm->flag & M_SL) SLtoGE(B_dm);
         else  fn_DisplayError(".../mathlib.c/MatrixMinusMatrix(): Haven't got time to deal with the M_UT and M_LT cases for B_dm");
      }

      for (_i=nels-1; _i>=0; _i--)   X[_i] = A[_i] - B[_i];
      if (A_dm->flag == B_dm->flag)  X_dm->flag = A_dm->flag;
      else  X_dm->flag = M_GE;
   }
}

void Matrix2PlusMinusMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm, TSdmatrix *C_dm, const double _alpha, const double _beta, const double _gamma) {
/*     //????? Not yet exhaust all possibilities of alpha, beta, and gamma to get most efficiency.  Add more as required.  10 February 2003.   ansi-c*/
/*     //Output: X = alpha*A + beta*B + gamma*C where X_dm is an m-by-n general matrix.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n general matrix.   ansi-c*/
/*     //  B_dm: m-by-n general matrix.   ansi-c*/
/*     //  C_dm: m-by-n general matrix.   ansi-c*/
/*     //  _alpha: a double scalar for A_dm.   ansi-c*/
/*     //  _beta: a double scalar for B_dm.   ansi-c*/
/*     //  _gamma: a double scalar for C_dm.   ansi-c*/
   int _i, _m, _n, nels;
   double *X, *A, *B, *C;


   if ( !X_dm || !A_dm || !B_dm || !C_dm )  fn_DisplayError(".../mathlib.c/Matrix2PlusMinusMatrix(): All input matrices must be created (memory-allocated)");
   else if ( !A_dm->flag || !B_dm->flag || !C_dm->flag ) fn_DisplayError(".../mathlib.c/Matrix2PlusMinusMatrix(): Some of R input matrices are not given values");
   else {
      _m = X_dm->nrows;
      _n = X_dm->ncols;
      nels = _m * _n;
      X = X_dm->M;
      A = A_dm->M;
      B = B_dm->M;
      C = C_dm->M;
   }


   if ( (_m != A_dm->nrows) || (_m != B_dm->nrows) || (_m != C_dm->nrows) || (_n != A_dm->ncols) || (_n != B_dm->ncols) || (_n != C_dm->ncols) )
      fn_DisplayError(".../mathlib.c/Matrix2PlusMinusMatrix(): Dimensions of all L and R input matrices must be same");
   else {
      if ( !(A_dm->flag & M_GE) ) {
         if (A_dm->flag & M_SU)   SUtoGE(A_dm);
         else if (A_dm->flag & M_SL)   SLtoGE(A_dm);
         else  fn_DisplayError(".../mathlib.c/Matrix2PlusMinusMatrix(): Haven't got time to deal with the M_UT and M_LT cases for A_dm");
      }
      if ( !(B_dm->flag & M_GE) ) {
         if (B_dm->flag & M_SU) SUtoGE(B_dm);
         else if (B_dm->flag & M_SL) SLtoGE(B_dm);
         else  fn_DisplayError(".../mathlib.c/Matrix2PlusMinusMatrix(): Haven't got time to deal with the M_UT and M_LT cases for B_dm");
      }
      if ( !(C_dm->flag & M_GE) ) {
         if (C_dm->flag & M_SU) SUtoGE(C_dm);
         else if (C_dm->flag & M_SL) SLtoGE(C_dm);
         else  fn_DisplayError(".../mathlib.c/Matrix2PlusMinusMatrix(): Haven't got time to deal with the M_UT and M_LT cases for C_dm");
      }


      if ( (_alpha==1.0) && (_beta==1.0) && (_gamma==1.0)  ) {
         for (_i=nels-1; _i>=0; _i--)  X[_i] = A[_i] + B[_i] + C[_i];
         if ( (A_dm->flag == B_dm->flag) && (A_dm->flag == C_dm->flag) )  X_dm->flag = A_dm->flag;
         else  X_dm->flag = M_GE;
      }
      else if ( (_alpha==1.0) && (_beta==1.0) && (_gamma==-1.0) ) {
         for (_i=nels-1; _i>=0; _i--)  X[_i] = A[_i] + B[_i] - C[_i];
         if ( (A_dm->flag == B_dm->flag) && (A_dm->flag == C_dm->flag) )  X_dm->flag = A_dm->flag;
         else  X_dm->flag = M_GE;
      }
      else if ( (_alpha==1.0) && (_gamma==1.0) ) {
         for (_i=nels-1; _i>=0; _i--)  X[_i] = A[_i] + _beta*B[_i] + C[_i];
         if ( (A_dm->flag == B_dm->flag) && (A_dm->flag == C_dm->flag) )  X_dm->flag = A_dm->flag;
         else  X_dm->flag = M_GE;
      }
      else {
/*           //Default for all cases (thus, may be most inefficient at this point.).   ansi-c*/
         for (_i=nels-1; _i>=0; _i--)  X[_i] = _alpha*A[_i] + _beta*B[_i] + _gamma*C[_i];
         if ( (A_dm->flag == B_dm->flag) && (A_dm->flag == C_dm->flag) )  X_dm->flag = A_dm->flag;
         else  X_dm->flag = M_GE;
      }
   }
}

void MatrixPlusConstantDiagUpdate(TSdmatrix *X_dm, const double _alpha) {
/*     //Output: X = X + diag([_alpha, ..., _alpha]) where X is an n-by-n square real matrix.   ansi-c*/
   int _i, nrows;
   double *M;

   if (!X_dm)  fn_DisplayError(".../mathlib.c/MatrixPlusConstantDiagUpdate():  Input matrix must be created (memory-allocated)");
   else if (!X_dm->flag)   fn_DisplayError(".../mathlib.c/MatrixPlusConstantDiagUpdate():  R input matrix must be given values");

   M = X_dm->M;
   nrows = X_dm->nrows;
   for (_i=square(nrows)-1; _i>=0; _i -= nrows+1)   M[_i] += _alpha;
}


void MatrixDotTimesMatrix(TSdmatrix *X_dm, TSdmatrix *A_dm, TSdmatrix *B_dm, const double _alpha, const double _beta)
{
/*     //$$$$$ If A_dm or B_dm or X_dm (when _beta!=0) is only upper or lower symmetric, it will be always converted to a general (and symmetric) matrix.  $$$$$$   ansi-c*/
/*     //Output:   ansi-c*/
/*     //  X_dm is m-by-n.   ansi-c*/
/*     //  X_dm = _alpha * A_dm .* B_dm + _beta * X_dm if X_dm != B_dm.   ansi-c*/
/*     //  X_dm = _alpha * A_dm .* X_dm + _beta * X_dm if X_dm = B_dm.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n double vector.   ansi-c*/
/*     //  B_dm: m-by-n double vector.   ansi-c*/
/*     //  _alpha: double scalar.   ansi-c*/
/*     //  _beta: a double scalar.   ansi-c*/
   int _i, nrows, ncols;
   double *X, *A, *B;


   if ( !X_dm || !A_dm || !B_dm) fn_DisplayError(".../mathlib.c/MatrixDotTimesMatrix(): All input matrices must be created (memory-allocated)");
   else if ( !A_dm->flag || !B_dm->flag )  fn_DisplayError(".../mathlib.c/MatrixDotTimesMatrix(): Both R input matrices must be given legal values");
   else if ( _beta && !X_dm->flag )  fn_DisplayError(".../mathlib.c/MatrixDotTimesMatrix(): L output matrix, X_dm, must be given legal values");
   else {
      nrows = X_dm->nrows;
      ncols = X_dm->ncols;
      X = X_dm->M;
      A = A_dm->M;
      B = B_dm->M;
   }
   if ( (nrows != A_dm->nrows) || (nrows != B_dm->nrows) || (ncols != A_dm->ncols) || (ncols != B_dm->ncols) )
      fn_DisplayError(".../mathlib.c/MatrixDotTimesMatrix(): Dimensions of all input matrices must be same");
   if ( !(A_dm->flag & M_GE) ) {
      if (A_dm->flag & M_SU)  SUtoGE(A_dm);
      else if (A_dm->flag & M_SL)  SLtoGE(A_dm);
      else   fn_DisplayError(".../mathlib.c/MatrixDotTimesMatrix(): Haven't got time to deal with the M_UT or M_LT cases for A_dm");
   }
   if ( !(B_dm->flag & M_GE) ) {
      if (B_dm->flag & M_SU)  SUtoGE(B_dm);
      else if (B_dm->flag & M_SL)  SLtoGE(B_dm);
      else   fn_DisplayError(".../mathlib.c/MatrixDotTimesMatrix(): Haven't got time to deal with the M_UT or M_LT cases for B_dm");
   }
   if ( _beta && !(X_dm->flag & M_GE) ) {
      if (X_dm->flag & M_SU)  SUtoGE(X_dm);
      else if (X_dm->flag & M_SL)  SLtoGE(X_dm);
      else   fn_DisplayError(".../mathlib.c/MatrixDotTimesMatrix(): Haven't got time to deal with the M_UT or M_LT cases for X_dm");
   }



   if ( _alpha==1.0 ) {
      if ( _beta==0.0 ) {
         for (_i=nrows*ncols-1; _i>=0; _i--)  X[_i] = A[_i] * B[_i];
         if (B_dm->flag != A_dm->flag)  X_dm->flag = M_GE;
         else  X_dm->flag = A_dm->flag;
      }
      else if ( _beta==1.0 ) {
         for (_i=nrows*ncols-1; _i>=0; _i--)  X[_i] += A[_i] * B[_i];
         if (X_dm->flag != A_dm->flag || X_dm->flag != B_dm->flag)  X_dm->flag = M_GE;
      }
      else {
         for (_i=nrows*ncols-1; _i>=0; _i--)   X[_i] = A[_i] * B[_i] + _beta * X[_i];
         if (X_dm->flag != A_dm->flag || X_dm->flag != B_dm->flag)  X_dm->flag = M_GE;
      }
   }
   else {
      if ( _beta==0.0 ) {
         for (_i=nrows*ncols-1; _i>=0; _i--)   X[_i] = _alpha * A[_i] * B[_i];
         if (B_dm->flag != A_dm->flag)  X_dm->flag = M_GE;
         else  X_dm->flag = A_dm->flag;
      }
      else if ( _beta==1.0 ) {
         for (_i=nrows*ncols-1; _i>=0; _i--)  X[_i] += _alpha * A[_i] * B[_i];
         if (X_dm->flag != A_dm->flag || X_dm->flag != B_dm->flag)  X_dm->flag = M_GE;
      }
      else {
         for (_i=nrows*ncols-1; _i>=0; _i--)  X[_i] = _alpha * A[_i] * B[_i] + _beta * X[_i];
         if (X_dm->flag != A_dm->flag || X_dm->flag != B_dm->flag)  X_dm->flag = M_GE;
      }
   }
}



void CopyVector0(TSdvector *x1_dv, const TSdvector *x2_dv) {
/*     //Ouputs:   ansi-c*/
/*     //  x1_dv, whose elements are copied from from x2_dv.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  Copying elements from x2_dv->v to x1_dv.   ansi-c*/
   int _n;

   if ( !x1_dv || !x2_dv || (x1_dv->n != x2_dv->n) )
      fn_DisplayError(".../mathlib.c/CopyVector0(): (1) all input pointers must be created (memory-allocated) and (2) dimensions of x1_dv and x2_dv must be compatible");
   else if ( !x2_dv->flag ) {
/*        // printf("x2_dv->flag is %d, and length is %d, and the vector is: \n", x2_dv->flag, x2_dv->n);   ansi-c*/
/*        // PrintVector(x2_dv, " %g ");   ansi-c*/
      fn_DisplayError(".../mathlib.c/CopyVector0(): R input vector must be given values");
   }
   else _n = x2_dv->n;

   if ( _n != x1_dv->n )  fn_DisplayError(".../mathlib.c/CopyVector0(): Both L and R input vectors must have the same length");
   memcpy(x1_dv->v, x2_dv->v, _n*sizeof(double));
   x1_dv->flag = V_DEF;
}


void CopyMatrix0(TSdmatrix *x1_dm, TSdmatrix *x2_dm) {
/*     //Deals with double matrices.   ansi-c*/
/*     //Copies the entire matrix x2_dm to x1_dm.   ansi-c*/
   int nrows1, ncols1, nrows2, ncols2;

   if ( !x1_dm || !x2_dm )  fn_DisplayError(".../mathlib.c/CopyMatrix0(): All input matrices must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopyMatrix0(): R input matrix must be given values");
   else {
      nrows1=x1_dm->nrows;
      ncols1=x1_dm->ncols;
      nrows2=x2_dm->nrows;
      ncols2=x2_dm->ncols;
   }


   if (nrows2 == nrows1 && ncols2 == ncols1) {
/*        //$$$$$$$$$$ 5/15/2003. At some point, the following if command should be got rid of completely.  For now, we keep this for maintaining the backward compatibility.   ansi-c*/
      if ( !(x2_dm->flag & M_GE) ) {
         if (x2_dm->flag & M_SU)   SUtoGE(x2_dm);
         else if (x2_dm->flag & M_SL)   SLtoGE(x2_dm);
/*  //         else  fn_DisplayError(".../mathlib.c/CopyMatrix0(): Haven't got time to deal with the M_UT and M_LT cases for x2_dm");   ansi-c*/
      }
/*        //Both matrices have the same size.   ansi-c*/
      memcpy(x1_dm->M, x2_dm->M, nrows1 * ncols1 * sizeof(double));
      x1_dm->flag = x2_dm->flag;

/*  //      #ifdef SWITCHTOTZCMATH            // define: use my own C math library ; undef: use others.   ansi-c*/
/*  //         memcpy(x1_dm->M, x2_dm->M, nrows1 * ncols1 * sizeof(double));   ansi-c*/
/*  //      #endif   ansi-c*/
/*  //      #ifdef SWITCHTOINTELCMATH             // define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
/*  //         cblas_dcopy(nrows1*ncols1, x2_dm->M, 1, x1_dm->M, 1);   ansi-c*/
/*  //      #endif   ansi-c*/
/*  //      x1_dm->flag = x2_dm->flag;   ansi-c*/
   }
   else fn_DisplayError(".../mathlib.c/CopyMatrix0(): Copying matrix (x2_m) and copied matrix (x1_dm) must have the same size");

/*  //?????????? The following is good, but should be used in CopySubmatrix0(), which might have already taken this into account.   ansi-c*/
/*  //   else if (nrows2 <= nrows1 && ncols2 <= ncols1) {   ansi-c*/
/*  //      if ( !(x2_dm->flag & M_GE) ) fn_DisplayError(".../mathlib.c/CopyMatrix0(): Haven't got time to deal with the M_UT and M_LT cases for x2_dm");   ansi-c*/
/*  //      //Size of x2_dm is smaller than that of x1_dm.   ansi-c*/
/*  //      for (_i=0; _i<ncols2; _i++) {   ansi-c*/
/*  //         loc1 = _i*nrows1;      //Points to the top of the column in x1_dm.   ansi-c*/
/*  //         loc2 = _i*nrows2;      //Points to the top of the column in x2_dm.   ansi-c*/
/*  //         memcpy((x1_dm->M+loc1), (x2_dm->M+loc2), x2_dm->nrows*sizeof(double));   ansi-c*/
/*  //      }   ansi-c*/
/*  //      x1_dm->flag = M_GE;   ansi-c*/
/*  //   }   ansi-c*/
/*  //   else fn_DisplayError(".../mathlib.c/CopyMatrix0(): number of rows (columns) of the copying matrix (x2) must be no greater than that of the copied matrix (x1)");   ansi-c*/
}

void CopyCellvec0(TSdcellvec *x1_dcv, TSdcellvec *x2_dcv)
{
/*     //Deals with double vectors.   ansi-c*/
/*     //Copies the entire cellvector x2_dcv to x1_dcv.   ansi-c*/
   int _i, ncells;
   if (!x1_dcv || !x2_dcv)  fn_DisplayError(".../mathlib.c/CopyCellvec0(): Both input cellvectors must be created (memory-allocated)");
   else if ( (ncells=x2_dcv->ncells) != x1_dcv->ncells )  fn_DisplayError(".../mathlib.c/CopyCellvec0(): Both input cellvectors must have exactly the same size");
   for (_i=ncells-1; _i>=0; _i--)  CopyVector0(x1_dcv->C[_i], x2_dcv->C[_i]);
}

void CopyCell0(TSdcell *x1_dc, TSdcell *x2_dc)
{
/*     //Deals with double matrices.   ansi-c*/
/*     //Copies the entire cell x2_dc to x1_dc.   ansi-c*/
   int _i, ncells;
   if (!x1_dc || !x2_dc)  fn_DisplayError(".../mathlib.c/CopyCell0(): Both input cells must be created (memory-allocated)");
   else if ( (ncells=x2_dc->ncells) != x1_dc->ncells )  fn_DisplayError(".../mathlib.c/CopyCell0(): Both input cells must have exactly the same size");
   for (_i=ncells-1; _i>=0; _i--)  CopyMatrix0(x1_dc->C[_i], x2_dc->C[_i]);
}


void CopySubmatrix0(TSdmatrix *x1_dm, TSdmatrix *x2_dm, const int br, const int bc, const int nrs, const int ncs)
{
/*     //Copies the nrs-by-ncs submatrix of x2_dm to the most left corner of x1_dm (i.e., at 0).   ansi-c*/
/*     //Note: br means the beginning of the row (*must* be 0 based) for this submatrix of x2_dm, inclusive;   ansi-c*/
/*     //      bc means the beginning of the column (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
   int _j, loc1, loc2,
       nrows1, ncols1, nrows2, ncols2;

   if ( !x1_dm || !x2_dm )  fn_DisplayError(".../mathlib.c/CopySubmatrix0(): All input matrices must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubmatrix0(): R input matrix must be given values");
   else {
      nrows1=x1_dm->nrows;
      ncols1=x1_dm->ncols;
      nrows2=x2_dm->nrows;
      ncols2=x2_dm->ncols;
   }

   if ( !(x2_dm->flag & M_GE) ) {
      if (x2_dm->flag & M_SU)   SUtoGE(x2_dm);
      else if (x2_dm->flag & M_SL)   SLtoGE(x2_dm);
      else  fn_DisplayError(".../mathlib.c/CopySubmatrix0(): Haven't got time to deal with the M_UT and M_LT cases for x2_dm");
   }

/*     //=== Performs the operation.   ansi-c*/
   if ( (bc+ncs)<=ncols2 && (br+nrs)<=nrows2 && ncs<=ncols1 && nrs<=nrows1 ) {
      for (_j=ncs-1; _j>=0; _j--) {
         loc1 = _j*nrows1;       /*  Points to the top of the column in x1_dm.   ansi-c*/
         loc2 = mos(br, bc+_j, nrows2);   /*  Points to the top of the column in the submatrix of x2_dm.   ansi-c*/
         #ifdef SWITCHTOTZCMATH                    /*   define: use my own C math library ; undef: use others.   ansi-c*/
            memcpy(x1_dm->M+loc1, x2_dm->M+loc2, nrs*sizeof(double));
         #endif
         #ifdef SWITCHTOINTELCMATH              /*   define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
            cblas_dcopy(nrs, x2_dm->M+loc2, 1, x1_dm->M+loc1, 1);
         #endif
      }
      x1_dm->flag = M_GE;
   }
   else  fn_DisplayError(".../mathlib.c/CopySubmatrix0(): the submatrix of x2_dm must be within the ranges of both x1_dm and x2_dm");
}

void CopySubmatrix(TSdmatrix *x1_dm, const int br1, const int bc1, TSdmatrix *x2_dm, const int br2, const int bc2, const int nrs, const int ncs) {
/*     //Copies the nrs-by-ncs submatrix of x2_dm to x1_dm at the specified location (br1, bc1).   ansi-c*/
/*     //Note: br1 means the beginning of the row (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      bc1 means the beginning of the column (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      br2 means the beginning of the row (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
/*     //      bc2 means the beginning of the column (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
   int _j, loc1, loc2,
       nrows1, ncols1, nrows2, ncols2;

   if ( !x1_dm || !x2_dm )  fn_DisplayError(".../mathlib.c/CopySubmatrix(): All input matrices must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubmatrix(): R input matrix must be given values");
   else {
      nrows1=x1_dm->nrows;
      ncols1=x1_dm->ncols;
      nrows2=x2_dm->nrows;
      ncols2=x2_dm->ncols;
   }

   if ( (bc2+ncs)<=ncols2 && (br2+nrs)<=nrows2 && (bc1+ncs)<=ncols1 && (br1+nrs)<=nrows1 ) {
      for (_j=ncs-1; _j>=0; _j--) {
         loc1 = mos(br1, bc1+_j, nrows1);       /*  Points to the top of the column in the submatrix of x1_dm.   ansi-c*/
         loc2 = mos(br2, bc2+_j, nrows2);   /*  Points to the top of the column in the submatrix of x2_dm.   ansi-c*/
         memcpy((x1_dm->M+loc1), (x2_dm->M+loc2), nrs*sizeof(double));
      }
      x1_dm->flag = M_GE;
   }
   else fn_DisplayError(".../mathlib.c/CopySubmatrix(): the submatrix of x2_dm must be within the ranges of both x1_dm and x2_dm");
}

#if defined( INTELCMATHLIBRARY )
void CopySubrowmatrix(TSdmatrix *x1_dm, const int br1, const int bc1, TSdmatrix *x2_dm, const int br2, const int bc2, const int nrs, const int ncs)
{
/*     //??????? NOT tested yet.   ansi-c*/
/*     //Copies the nrs-by-ncs submatrix of x2_dm to x1_dm at the specified location (br1, bc1).   ansi-c*/
/*     //Note: br1 means the beginning of the row (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      bc1 means the beginning of the column (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      br2 means the beginning of the row (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
/*     //      bc2 means the beginning of the column (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
   int _i, loc1, loc2,
       nrows1, ncols1, nrows2, ncols2;

   if ( !x1_dm || !x2_dm )  fn_DisplayError(".../mathlib.c/CopySubrowmatrix(): All input matrices must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubrowmatrix(): R input matrix must be given values");
   else {
      nrows1=x1_dm->nrows;
      ncols1=x1_dm->ncols;
      nrows2=x2_dm->nrows;
      ncols2=x2_dm->ncols;
   }

   if ( (bc2+ncs)<=ncols2 && (br2+nrs)<=nrows2 && (bc1+ncs)<=ncols1 && (br1+nrs)<=nrows1 ) {
      for (_i=nrs-1; _i>=0; _i--) {
         loc1 = mos(br1+_i, bc1, nrows1);       /*  Points to the beginning of the row in the submatrix of x1_dm.   ansi-c*/
         loc2 = mos(br2+_i, bc2, nrows2);   /*  Points to the beginning of the row in the submatrix of x2_dm.   ansi-c*/
         cblas_dcopy(ncs, x2_dm->M+loc2, nrows2, x1_dm->M+loc1, nrows1);
      }
   }
   else fn_DisplayError(".../mathlib.c/CopySubrowmatrix(): the submatrix of x2_dm must be within the range of itself as well as x1_dm");
}
#else
/*  //   Havent got time to code up the default using Linux BLAS.   ansi-c*/
#endif


#if defined( INTELCMATHLIBRARY )
void CopySubmatrix2rowmatrix(TSdmatrix *x1_dm, const int br1, const int bc1, TSdmatrix *x2_dm, const int br2, const int bc2, const int nrs, const int ncs)
{
/*     //Column by column operation on x2_dm: coping the transpose of the nrs-by-ncs submatrix of x2_dm to x1_dm at the specified location (br1, bc1).   ansi-c*/
/*     //Note: br1 means the beginning of the row (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      bc1 means the beginning of the column (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      br2 means the beginning of the row (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
/*     //      bc2 means the beginning of the column (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
   int _j, loc1, loc2,
       nrows1, ncols1, nrows2, ncols2;

   if ( !x1_dm || !x2_dm )  fn_DisplayError(".../mathlib.c/CopySubmatrix2rowmatrix(): All input matrices must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubmatrix2rowmatrix(): R input matrix must be given values");
   else {
      nrows1=x1_dm->nrows;
      ncols1=x1_dm->ncols;
      nrows2=x2_dm->nrows;
      ncols2=x2_dm->ncols;
   }

   if ( (bc2+ncs)<=ncols2 && (br2+nrs)<=nrows2 && (bc1+nrs)<=ncols1 && (br1+ncs)<=nrows1 ) {
      for (_j=ncs-1; _j>=0; _j--) {
         loc1 = mos(br1+_j, bc1, nrows1);       /*  Points to the beginning of the row in the submatrix of x1_dm.   ansi-c*/
         loc2 = mos(br2, bc2+_j, nrows2);   /*  Points to the top of the column in the submatrix of x2_dm.   ansi-c*/
         cblas_dcopy(nrs, x2_dm->M+loc2, 1, x1_dm->M+loc1, nrows1);
      }
   }
   else fn_DisplayError(".../mathlib.c/CopySubmatrix2rowmatrix(): the submatrix of x2_dm must be within the range of x2_dm and its transpose must be within the range of x1_dm");
}
#else
/*  //   Havent got time to code up the default using Linux BLAS.   ansi-c*/
#endif


#if defined( INTELCMATHLIBRARY )
void CopySubrowmatrix2matrix(TSdmatrix *x1_dm, const int br1, const int bc1, TSdmatrix *x2_dm, const int br2, const int bc2, const int nrs, const int ncs)
{
/*     //??????? NOT tested yet.   ansi-c*/
/*     //Row by row operation on x2_dm: coping the transpose of the nrs-by-ncs submatrix of x2_dm to x1_dm at the specified location (br1, bc1).   ansi-c*/
/*     //Note: br1 means the beginning of the row (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      bc1 means the beginning of the column (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //      br2 means the beginning of the row (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
/*     //      bc2 means the beginning of the column (*must* be 0 based) for this submatrix of x2_dm, inclusive.   ansi-c*/
   int _i, loc1, loc2,
       nrows1, ncols1, nrows2, ncols2;

   if ( !x1_dm || !x2_dm )  fn_DisplayError(".../mathlib.c/CopySubrowmatrix2matrix(): All input matrices must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubrowmatrix2matrix(): R input matrix must be given values");
   else {
      nrows1=x1_dm->nrows;
      ncols1=x1_dm->ncols;
      nrows2=x2_dm->nrows;
      ncols2=x2_dm->ncols;
   }

   if ( (bc2+ncs)<=ncols2 && (br2+nrs)<=nrows2 && (bc1+nrs)<=ncols1 && (br1+ncs)<=nrows1 ) {
      for (_i=nrs-1; _i>=0; _i--) {
         loc1 = mos(br1, bc1+_i, nrows1);       /*  Points to the top of the column in the submatrix of x1_dm.   ansi-c*/
         loc2 = mos(br2+_i, bc2, nrows2);   /*  Points to the beginning of the row in the submatrix of x2_dm.   ansi-c*/
         cblas_dcopy(ncs, x2_dm->M+loc2, nrows2, x1_dm->M+loc1, 1);
      }
   }
   else fn_DisplayError(".../mathlib.c/CopySubrowmatrix2matrix(): the submatrix of x2_dm must be within the range of itself as well as x1_dm");
}
#else
/*  //   Havent got time to code up the default using Linux BLAS.   ansi-c*/
#endif


void CopySubvector(TSdvector *x1_dv, const int ptrloc1, const TSdvector *x2_dv, const int ptrloc2, const int nels) {
/*     //Ouputs:   ansi-c*/
/*     //  x1_dv, whose elements from x1_dv->v+ptrloc1 to x1_dv->v+ptrloc1+nels-1 are copied from from x2_dv->v.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  Copying elements from x2_dv->v+ptrloc2 to x2_dv->v+ptrloc2+nels-1 (to x1_dv->v).   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
/*     //  ptrloc1: pointer location for x1_dv->v where the copy begins, inclusive.   ansi-c*/
/*     //  ptrloc2: pointer location for x2_dv->v where the copy begins, inclusive.   ansi-c*/

   if ( !x1_dv || !x2_dv )  fn_DisplayError(".../mathlib.c/CopySubvector(): All input vectors must be created (memory-allocated)");
   else if (!x2_dv->flag)  fn_DisplayError(".../mathlib.c/CopySubvector(): R input vector must be given values");

   if ( (ptrloc2+nels)<=x2_dv->n && (ptrloc1+nels)<=x1_dv->n) {
      memcpy(x1_dv->v+ptrloc1, x2_dv->v+ptrloc2, nels*sizeof(double));
      x1_dv->flag = V_DEF;
   }
   else fn_DisplayError(".../mathlib.c/CopySubvector(): Copying (copied) elements are outside the dimension of the copying vector x2_dv (the copied vector x1_dv)");
}
/*  //---   ansi-c*/
void CopySubvector_int(TSivector *x1_iv, const int ptrloc1, const TSivector *x2_iv, const int ptrloc2, const int nels)
{
/*     //Ouputs:   ansi-c*/
/*     //  x1_iv, whose elements from x1_iv->v+ptrloc1 to x1_iv->v+ptrloc1+nels-1 are copied from from x2_iv->v.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  Copying elements from x2_iv->v+ptrloc2 to x2_iv->v+ptrloc2+nels-1 (to x1_iv->v).   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
/*     //  ptrloc1: pointer location for x1_iv->v where the copy begins, inclusive.   ansi-c*/
/*     //  ptrloc2: pointer location for x2_iv->v where the copy begins, inclusive.   ansi-c*/

   if ( !x1_iv || !x2_iv )  fn_DisplayError(".../mathlib.c/CopySubvector_int(): All input vectors must be created (memory-allocated)");
   else if (!x2_iv->flag)  fn_DisplayError(".../mathlib.c/CopySubvector_int(): R input vector must be given values");

   if ( (ptrloc2+nels)<=x2_iv->n && (ptrloc1+nels)<=x1_iv->n) {
      memcpy(x1_iv->v+ptrloc1, x2_iv->v+ptrloc2, nels*sizeof(int));
      x1_iv->flag = V_DEF;
   }
   else fn_DisplayError(".../mathlib.c/CopySubvector_int(): Copying (copied) elements are outside the dimension of the copying vector x2_iv (the copied vector x1_iv)");
}

void CopySubmatrix2vector(TSdvector *x1_dv, const int ptrloc1, TSdmatrix *x2_dm, const int br, const int bc, const int nels)
{
/*     //Ouputs:   ansi-c*/
/*     //  x1_dv whose elements from x1_dv->v+ptrloc1 to x1_dv->v+ptrloc1+nels-1 are copied from x2_dm->M[br,bc] onward (inclusive)   ansi-c*/
/*     //         where copied elements can run column by column all the way to the end of x2_dm->M[end,end].  Thus, this function   ansi-c*/
/*     //         can be used as the Matlab reshape command.   ansi-c*/
/*     //Inputs:  Copying elements from x2_dm->M+ptrloc2 to x2_dm->M+ptrloc2+nels-1 (to x1_dv->v).   ansi-c*/
/*     //  ptrloc1: pointer location for x1_dv->v where the copy begins, inclusive.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
   int ptrloc2;    /*  pointer location for x2_dm->M where the copy begins.   ansi-c*/

   if ( !x1_dv || !x2_dm ) fn_DisplayError(".../mathlib.c/CopySubmatrix2vector():  All input pointers must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubmatrix2vector(): R input matrix must be given values");
   else   ptrloc2 = mos(br,bc,x2_dm->nrows);   /*    ptrloc2: pointer location for x2_dm->M where the copy begins.   ansi-c*/

   if ( !(x2_dm->flag & M_GE) ) {
      if (x2_dm->flag & M_SU)   SUtoGE(x2_dm);
      else if (x2_dm->flag & M_SL)   SLtoGE(x2_dm);
      else  fn_DisplayError(".../mathlib.c/CopySubmatrix2vector(): Haven't got time to deal with the M_UT and M_LT cases for x2_dm");
   }



   if ( (ptrloc2+nels)<=(x2_dm->nrows*x2_dm->ncols) && (ptrloc1+nels)<=x1_dv->n) {
      memcpy(x1_dv->v+ptrloc1, x2_dm->M+ptrloc2, nels*sizeof(double));
      x1_dv->flag = V_DEF;
   }
   else fn_DisplayError(".../mathlib.c/CopySubmatrix2vector(): Copying (copied) elements are outside the dimension of the copying matrix x2_dm (the copied vector x1_dv)");
}

void CopySubmatrix2vector_sub(TSdvector *x1_dv, const int ptrloc1, TSdmatrix *x2_dm, const int br, const int bc, const int nrs, const int ncs) {
/*     //Ouputs: Unlike CopySubmatrix2vector, _sub means a submatrix, NOT just one column, of the copying matrix x2_dm.   ansi-c*/
/*     //          The copying submatrix must start at (br, bc) ane end at (br+nrs-1, bc+ncs-1).   ansi-c*/
/*     //          Copying is done column by column.   ansi-c*/
/*     //  x1_dv, whose elements from x1_dv->v+ptrloc1 to x1_dv->v+ptrloc1+nrs*ncs-1 are copied from the submatrix of of x2_dm->m.   ansi-c*/
/*     //Inputs:  The copying submatrix of x2_dm->M.   ansi-c*/
/*     //  ptrloc1: inclusive pointer location for x1_dv->v where the copy begins.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  nrs:  number of rows to be copied.   ansi-c*/
/*     //  ncs:  number of colums to be copied.   ansi-c*/
   int nrows, ncols, _j, loc1;
   double *v, *M;

   if ( !x1_dv || !x2_dm ) fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_sub():  All input pointers must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_sub(): R input matrix must be given values");
   else {
      v = x1_dv->v;
      M = x2_dm->M;
      nrows = x2_dm->nrows;
      ncols = x2_dm->ncols;
   }

   if ( !(x2_dm->flag & M_GE) ) {
      if (x2_dm->flag & M_SU)   SUtoGE(x2_dm);
      else if (x2_dm->flag & M_SL)   SLtoGE(x2_dm);
      else  fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_sub(): Haven't got time to deal with the M_UT and M_LT cases for x2_dm");
   }

   if ( (bc+ncs)<=ncols && (br+nrs)<=nrows && (ptrloc1+ncs*nrs)<=x1_dv->n ) {
      loc1 = ptrloc1;
      for (_j=0; _j<ncs; _j++) {
         memcpy(v+loc1, M+mos(br, bc+_j, nrows), nrs*sizeof(double));  /*  mos(br, bc+_j, nrows): Points to the top of the column in the submatrix of x2_dm.   ansi-c*/
         loc1 += nrs;   /*  Must be after memcpy().   ansi-c*/
      }
      x1_dv->flag = V_DEF;
   }
   else fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_sub(): the submatrix of x2_dm must be within the ranges of both x1_dv and x2_dm");
}

void CopySubmatrix2vector_int(TSivector *x1_iv, const int ptrloc1, TSimatrix *x2_im, const int br, const int bc, const int nels) {
/*     //Ouputs:   ansi-c*/
/*     //  x1_iv, whose elements from x1_iv->v+ptrloc1 to x1_iv->v+ptrloc1+nels-1 are copied from a column of x2_im->m.   ansi-c*/
/*     //Inputs:  Copying elements from x2_im->M+ptrloc2 to x2_im->M+ptrloc2+nels-1 (to x1_iv->v).   ansi-c*/
/*     //  ptrloc1: pointer location for x1_iv->v where the copy begins, inclusive.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x2_im, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x2_im, inclusive.   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
   int ptrloc2;    /*  pointer location for x2_im->M where the copy begins.   ansi-c*/

   if ( !x1_iv || !x2_im ) fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_int():  All input pointers must be created (memory-allocated)");
   else   ptrloc2 = mos(br,bc,x2_im->nrows);   /*    ptrloc2: pointer location for x2_im->M where the copy begins.   ansi-c*/



   if ( (ptrloc2+nels)<=(x2_im->nrows*x2_im->ncols) && (ptrloc1+nels)<=x1_iv->n) {
      memcpy(x1_iv->v+ptrloc1, x2_im->M+ptrloc2, nels*sizeof(int));
      x1_iv->flag = V_DEF;
   }
   else fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_int(): Copying (copied) elements are outside the dimension of the copying matrix x2_im (the copied vector x1_iv)");
}

void CopySubmatrix2vector_row(TSdvector *x1_dv, const int ptrloc1, TSdmatrix *x2_dm, const int br, const int bc, const int nels) {
/*     //This is much less efficient because we copy a row of x2_dm where x2_dm is a column-major matrix.  But sometimes,   ansi-c*/
/*     //  transposing x2_dm and then using CopySubmatrix2vector() proves more costly if the transpose has to be done in each iteration.   ansi-c*/
/*     //If SWITCHINTELCMATH is activated, it may achieve efficiency.   ansi-c*/
/*     //   ansi-c*/
/*     //Ouputs: copying a row of x2_dm to a vector.   ansi-c*/
/*     //  x1_dv, whose elements from x1_dv->v+ptrloc1 to x1_dv->v+ptrloc1+nels-1 are copied from a row of x2_dm->m.   ansi-c*/
/*     //Inputs:  Copying elements from x2_dm->M(br, bc) to x2_dm->M(br, bc+nels-1) (to x1_dv->v).   ansi-c*/
/*     //  ptrloc1: inclusive pointer location for x1_dv->v where the copy begins.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
   int nrows;
   #if !defined(SWITCHTOINTELCMATH)                   /*   define: use my own C math library ; undef: use others.   ansi-c*/
      int _i;
   #endif
   double *v, *M;

   if ( !x1_dv || !x2_dm ) fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_row():  All input pointers must be created (memory-allocated)");
   else if ( !x2_dm->flag )  fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_row(): R input matrix must be given values");
   else {
      v = x1_dv->v;
      M = x2_dm->M;
      nrows = x2_dm->nrows;
   }

   if ( !(x2_dm->flag & M_GE) ) {
      if (x2_dm->flag & M_SU)   SUtoGE(x2_dm);
      else if (x2_dm->flag & M_SL)   SLtoGE(x2_dm);
      else  fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_row(): Haven't got time to deal with the M_UT and M_LT cases for x2_dm");
   }


   if ( (bc+nels)<=x2_dm->ncols && (ptrloc1+nels)<=x1_dv->n) {
      #if defined (INTELCMATHLIBRARY)              /*   define: use Intek MKL LAPACK library; undef: use others.   ansi-c*/
         cblas_dcopy(nels, M+mos(br, bc, nrows), nrows, v+ptrloc1, 1);   /*  mos(): inclusive pointer location for x2_dm where the copy begins.   ansi-c*/
      #else      /*  Default to  SWITCHTOTZCMATH                      ansi-c*/
/*           //for (_i=0; _i<nels; _i++)   v[ptrloc1+_i] = M[mos(br, bc+_i, nrows)];   ansi-c*/
         for (_i=nels-1; _i>=0; _i--)   v[ptrloc1+_i] = M[mos(br, bc+_i, nrows)];    /*  Changed above to this.  9/2/03.   ansi-c*/
      #endif
      x1_dv->flag = V_DEF;
   }
   else fn_DisplayError(".../mathlib.c/CopySubmatrix2vector_row(): Copying (copied) elements are outside the dimension of the copying matrix x2_dm (the copied vector x1_dv)");
}

void CopySubvector2matrix(TSdmatrix *x1_dm, const int br, const int bc, const TSdvector *x2_dv, const int ptrloc2, const int nels) {
/*     //Ouputs:  only the ``bc''th column of the matrix is copied.  If this is too restrictive, see CopySubvector2matrix_unr().   ansi-c*/
/*     //  Copied elements (x1_dm->M+ptrloc1 to x1_dv->M+ptrloc1+nels-1) from x2_dv->v.   ansi-c*/
/*     //  Always sets x1_dm->flag = M_GE after the call to this function.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  Copying elements from x2_dv->v+ptrloc2 to x2_dv->v+ptrloc2+nels-1 (to x1_dm->M).   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  ptrloc2: pointer location for x2_dv->v where the copy begins, inclusive.   ansi-c*/
   int ptrloc1, nrows, ncols;

   if ( !x1_dm || !x2_dv )  fn_DisplayError(".../mathlib.c/CopySubvector2matrix():  All input pointers must be created (memory-allocated)");
   else if (!x2_dv->flag)  fn_DisplayError(".../mathlib.c/CopySubvector2matrix():  R input vector must be given values");
   else {
      nrows = x1_dm->nrows;
      ncols = x1_dm->ncols;
      ptrloc1 = mos(br, bc, x1_dm->nrows);     /*  ptrloc1: pointer location for x1_dm->M where the copy begins.   ansi-c*/
   }


   if ( (ptrloc2+nels)<=x2_dv->n && (br+nels)<=nrows ) {
      memcpy(x1_dm->M+ptrloc1, x2_dv->v+ptrloc2, nels*sizeof(double));
      x1_dm->flag = M_GE;   /*  Always reset to a general matrix because it will almost surely break, say, the original symmetry of the matrix x1_dm if x1_dm->flag exists in the first place.   ansi-c*/
/*        //-------if (!x1_dm->flag)  x1_dm->flag = M_GE;   //Set to a general matrix if this matrix is not set yet.   ansi-c*/
   }
   else fn_DisplayError(".../mathlib.c/CopySubvector2matrix(): Copying (copied) elements are outside the (row) dimension of the copying vector x2_dv (the copied matrix x1_dm)");
}


#if defined( INTELCMATHLIBRARY )
void CopySubvector2rowmatrix(TSdmatrix *x1_dm, const int br, const int bc, const TSdvector *x2_dv, const int ptrloc2, const int nels)
{
/*     //Ouputs:   ansi-c*/
/*     //  Only the ``br''th row of the matrix x1_dm (starting from the ``bc''th column) is copied from x2_dv->v.   ansi-c*/
/*     //  Always sets x1_dm->flag = M_GE after the call to this function.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  Copying elements from x2_dv->v+ptrloc2 to x2_dv->v+ptrloc2+nels-1 (to x1_dm).   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  ptrloc2: pointer location for x2_dv->v where the copy begins, inclusive.   ansi-c*/
   int ptrloc1, nrows, ncols;

   if ( !x1_dm || !x2_dv )  fn_DisplayError(".../mathlib.c/CopySubvector2rowmatrix():  All input pointers must be created (memory-allocated)");
   else if (!x2_dv->flag)  fn_DisplayError(".../mathlib.c/CopySubvector2rowmatrix():  R input vector must be given values");
   else {
      nrows = x1_dm->nrows;
      ncols = x1_dm->ncols;
      ptrloc1 = mos(br, bc, x1_dm->nrows);     /*  ptrloc1: pointer location for x1_dm->M where the copy begins.   ansi-c*/
   }


   if ( (ptrloc2+nels)<=x2_dv->n && (bc+nels)<=ncols ) {
      cblas_dcopy(nels, x2_dv->v+ptrloc2, 1, x1_dm->M+ptrloc1, nrows);
      x1_dm->flag = M_GE;   /*  Always reset to a general matrix because it will almost surely break, say, the original symmetry of the matrix x1_dm if x1_dm->flag exists in the first place.   ansi-c*/
/*        //-------if (!x1_dm->flag)  x1_dm->flag = M_GE;   //Set to a general matrix if this matrix is not set yet.   ansi-c*/
   }
   else fn_DisplayError(".../mathlib.c/CopySubvector2rowmatrix(): Copying (copied) elements are outside the (row) dimension of the copying vector x2_dv (the copied matrix x1_dm)");
}
#else
/*  //   Havent got time to code up the default using Linux BLAS.   ansi-c*/
#endif


void CopySubvector2matrix_sub(TSdmatrix *x1_dm, const int br, const int bc, const int nrs, const int ncs, TSdvector *x2_dv, const int ptrloc2) {
/*     //Ouputs: Unlike CopySubvector2matrix, _sub means a submatrix, NOT just one column, of the copied matrix x1_dm.   ansi-c*/
/*     //          The copyed submatrix must start at (br, bc) ane end at (br+nrs-1, bc+ncs-1).   ansi-c*/
/*     //  x1_dm: The copyed submatrix of x2_dm->M.   ansi-c*/

/*     //Inputs:   ansi-c*/
/*     //  x2_dv: whose elements from x2_dv->v+ptrloc2 to x2_dv->v+ptrloc2+nrs*ncs-1 are copied to the submatrix of of x1_dm->m.   ansi-c*/
/*     //  ptrloc2: inclusive pointer location for x2_dv->v where the copy begins.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x1_dm, inclusive.   ansi-c*/
/*     //  nrs:  number of rows to be copied.   ansi-c*/
/*     //  ncs:  number of colums to be copied.   ansi-c*/
   int nrows, ncols, _j, loc2;
   double *v, *M;

   if ( !x1_dm || !x2_dv ) fn_DisplayError(".../mathlib.c/CopySubvector2matrix_sub():  All input pointers must be created (memory-allocated)");
   else if ( !x2_dv->flag )  fn_DisplayError(".../mathlib.c/CopySubvector2matrix_sub(): R input vector must be given values");
   else {
      v = x2_dv->v;
      M = x1_dm->M;
      nrows = x1_dm->nrows;
      ncols = x1_dm->ncols;
   }

   if ( (bc+ncs)<=ncols && (br+nrs)<=nrows && (ncs*nrs)<=(x2_dv->n-ptrloc2) ) {
      loc2 = ptrloc2;
      for (_j=0; _j<ncs; _j++) {
         memcpy(M+mos(br, bc+_j, nrows), v+loc2, nrs*sizeof(double));  /*  mos(br, bc+_j, nrows): Points to the top of the column in the submatrix of x2_dm.   ansi-c*/
         loc2 += nrs;   /*  Must be after memcpy().   ansi-c*/
      }
      x1_dm->flag = M_GE;
   }
   else fn_DisplayError(".../mathlib.c/CopySubvector2matrix_sub(): the submatrix of x1_dm must be within the ranges of both x1_dm and x2_dv");
}

void CopySubvector2matrix_unr(TSdmatrix *x1_dm, const int br, const int bc, const TSdvector *x2_dv, const int ptrloc2, const int nels) {
/*     //Ouputs: _unr means that there is no such a restriction that only the ``bc''th column to be copied in the matrix.  Copied   ansi-c*/
/*     //          elements can affect several columns of the matrix other than the ``bc'' column, but one must be aware   ansi-c*/
/*     //          that the bc+1 column may start before the specified br.  Thus, this function can be used as the Matlab reshape function.   ansi-c*/
/*     //  Copied elements (x1_dm->M+ptrloc1 to x1_dv->M+ptrloc1+nels-1) from x2_dv->v.   ansi-c*/
/*     //  Always sets x1_dm->flag = M_GE after the call to this function.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  Copying elements from x2_dv->v+ptrloc2 to x2_dv->v+ptrloc2+nels-1 (to x1_dm->M).   ansi-c*/
/*     //  nels: number of elements to be copied.   ansi-c*/
/*     //  br: beginning of the row (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  bc: beginning of the column (*must* be 0 based) for x2_dm, inclusive.   ansi-c*/
/*     //  ptrloc2: pointer location for x2_dv->v where the copy begins, inclusive.   ansi-c*/
   int ptrloc1, nrows, ncols;

   if ( !x1_dm || !x2_dv )  fn_DisplayError(".../mathlib.c/CopySubvector2matrix_unr():  All input pointers must be created (memory-allocated)");
   else if (!x2_dv->flag)  fn_DisplayError(".../mathlib.c/CopySubvector2matrix_unr():  R input vector must be given values");
   else {
      nrows = x1_dm->nrows;
      ncols = x1_dm->ncols;
      ptrloc1 = mos(br, bc, x1_dm->nrows);     /*  ptrloc1: pointer location for x1_dm->M where the copy begins.   ansi-c*/
   }


   if ( (ptrloc2+nels)<=x2_dv->n && (ptrloc1+nels)<=(nrows*ncols) ) {
      memcpy(x1_dm->M+ptrloc1, x2_dv->v+ptrloc2, nels*sizeof(double));
      x1_dm->flag = M_GE;   /*  Always reset to a general matrix because it will almost surely break, say, the original symmetry of the matrix x1_dm if x1_dm->flag exists in the first place.   ansi-c*/
/*        //-------if (!x1_dm->flag)  x1_dm->flag = M_GE;   //Set to a general matrix if this matrix is not set yet.   ansi-c*/
   }
   else fn_DisplayError(".../mathlib.c/CopySubvector2matrix_unr(): Copying (copied) elements are outside the dimension of the copying vector x2_dv (the copied matrix x1_dm)");
}

void TransposeSquare(TSdmatrix *B_dm, TSdmatrix *A_dm) {
/*     //???????? Some options are NOT test yet.  5/27/03.  ???????????   ansi-c*/
/*     // Transposes the n-by-n matrix A_dm to the n-by-n matrix B_dm.   ansi-c*/
/*     //   If A_dm = B_dm, B_dm will be replaced by transposed values.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  B_dm: n-by-n matrix (memory already allocated outside this function).   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: n-by-n matrix to be transposed.   ansi-c*/
   int _i, _j, _n, Aflag;
   double *A, *B, tmpd;

/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if ( !B_dm || !A_dm )  fn_DisplayError(".../mathlib.c/TransposeSquare():  Input matrices must be created (memory-allocated)");
   if ( ((_n=A_dm->nrows) != B_dm->ncols) || (_n != B_dm->nrows) || (_n != A_dm->ncols) )
      fn_DisplayError(".../mathlib.c/TransposeSquare(): Both input matrices must be square");
/*     //This is too tough by killing the program. if ( ((Aflag=A_dm->flag) & M_SU) && (Aflag & M_SL) )  fn_DisplayError(".../mathlib.c/TransposeSquare():  Matrix is already both SU and SL, so there is no need to transpose.  Check a possible bug in your program");   ansi-c*/
/*           //The above checking is very important even though it takes about 4 clock cycles, because   ansi-c*/
/*           //    (1) one may make a mistake to mis-use this matrix over and over again;   ansi-c*/
/*           //    (2) if there is no mistake, then no need to transpose this matrix.   ansi-c*/
   if ( ((Aflag=A_dm->flag) & M_SU) && (Aflag & M_SL) )
   {
      #if defined (USE_DEBUG_FILE)
      fprintf(FPTR_DEBUG, "\nWARNING: .../mathlib.c/TransposeSquare():  the matrix is already both SU and SL, so there is no need to transpose.\n");
      fflush(FPTR_DEBUG);
      #else
      printf("\nWARNING: .../mathlib.c/TransposeSquare():  the matrix is already both SU and SL, so there is no need to transpose.\n");
      fflush(stdout);
      #endif

      if ( (A=A_dm->M) != (B=B_dm->M) )
      {
         CopyMatrix0(B_dm, A_dm);
         return;
      }
      else  return;
   }



/*     //=== Transposing the square matrix A_dm.   ansi-c*/
   if ( (A=A_dm->M) != (B=B_dm->M) ) {
      if (Aflag & M_GE) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++) {
/*                 //Off-diagonal elements.   ansi-c*/
               B[_i*_n+_j] = A[_j*_n+_i];
               B[_j*_n+_i] = A[_i*_n+_j];
            }
         for (_i=square(_n)-1; _i>=0; _i -= _n+1)  B[_i] = A[_i];    /*  Diagonal elements.   ansi-c*/
         switch (Aflag) {
            case (M_GE | M_UT):
               B_dm->flag = M_GE | M_LT;
               break;
            case (M_GE | M_LT):
               B_dm->flag = M_GE | M_UT;
               break;
            default:
               B_dm->flag = M_GE;
         }
      }
      else if (Aflag & M_SU) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               B[mos(_i, _j, _n)] = A[mos(_j, _i, _n)];   /*  Off-diagonal elements.   ansi-c*/
         for (_i=square(_n)-1; _i>=0; _i -= _n+1)
            B[_i] = A[_i];    /*  Diagonal elements.   ansi-c*/
         B_dm->flag = M_SL;
      }
      else if (Aflag & M_SL) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               B[mos(_j, _i, _n)] = A[mos(_i, _j, _n)];   /*  Off-diagonal elements.   ansi-c*/
         for (_i=square(_n)-1; _i>=0; _i -= _n+1)
            B[_i] = A[_i];    /*  Diagonal elements.   ansi-c*/
         B_dm->flag = M_SU;
      }
      else if (Aflag & M_UT) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               B[mos(_i, _j, _n)] = A[mos(_j, _i, _n)];   /*  Off-diagonal elements.   ansi-c*/
         for (_i=square(_n)-1; _i>=0; _i -= _n+1)
            B[_i] = A[_i];    /*  Diagonal elements.   ansi-c*/
         B_dm->flag = M_LT;
      }
      else if (Aflag & M_LT) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               B[mos(_j, _i, _n)] = A[mos(_i, _j, _n)];   /*  Off-diagonal elements.   ansi-c*/
         for (_i=square(_n)-1; _i>=0; _i -= _n+1)
            B[_i] = A[_i];    /*  Diagonal elements.   ansi-c*/
         B_dm->flag = M_UT;
      }
   }
   else {
/*        // ????? NOT tested yet.  5/27/03.  ????????????   ansi-c*/
      if ( (Aflag & M_GE) && (Aflag & M_UT) ) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++) {
/*                 //Off-diagonal elements.   ansi-c*/
               A[mos(_i, _j, _n)] = A[mos(_j, _i, _n)];
               A[mos(_j, _i, _n)] = 0.0;
            }
         A_dm->flag = M_GE | M_LT;
      }
      else if ( (Aflag & M_GE) && (Aflag & M_LT) ) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++) {
/*                 //Off-diagonal elements.   ansi-c*/
               A[mos(_j, _i, _n)] = A[mos(_i, _j, _n)];
               A[mos(_i, _j, _n)] = 0.0;
            }
         A_dm->flag = M_GE | M_UT;
      }
      else if (Aflag & M_GE) {
/*           //Tested.  Fine.  10 Oct. 03.   ansi-c*/
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               swap(A[mos(_i,_j,_n)], A[mos(_j,_i,_n)], tmpd);  /*  Off-diagonal elements.   ansi-c*/
         A_dm->flag = M_GE;
      }
      else if (Aflag & M_SU) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               A[mos(_i, _j, _n)] = A[mos(_j, _i, _n)];   /*  Off-diagonal elements.   ansi-c*/
         A_dm->flag = M_GE | M_SU | M_SL;
      }
      else if (Aflag & M_SL) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               A[mos(_j, _i, _n)] = A[mos(_i, _j, _n)];   /*  Off-diagonal elements.   ansi-c*/
         A_dm->flag = M_GE | M_SU | M_SL;
      }
      else if (Aflag & M_UT) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               A[mos(_i, _j, _n)] = A[mos(_j, _i, _n)];   /*  Off-diagonal elements.   ansi-c*/
         A_dm->flag = M_LT;
      }
      else if (Aflag & M_LT) {
         for (_j=0; _j<_n; _j++)
            for (_i=_j+1; _i<_n; _i++)
               A[mos(_j, _i, _n)] = A[mos(_i, _j, _n)];   /*  Off-diagonal elements.   ansi-c*/
         A_dm->flag = M_UT;
      }
      else  fn_DisplayError(".../mathlib.c/TransposeSquare():  Input matrix is not M_GE, M_SU, M_SL, M_UT, or M_LT.  Check the matrix to see why transpose is still required");
   }


   /**   OLD code, which has some errors, I belive.  5/28/03.
   //=== Transposing the square matrix A_dm.
   if ( (A=A_dm->M) != (B=B_dm->M) ) {
      for (_j=0; _j<_n; _j++)
         for (_i=_j+1; _i<_n; _i++) {
            //Off-diagonal elements.
            B[_i*_n+_j] = A[_j*_n+_i];
            B[_j*_n+_i] = A[_i*_n+_j];
         }
      for (_i=square(_n)-1; _i>=0; _i -= _n+1)
         B[_i] = A[_i];   //Diagonal elements.

      switch (A_dm->flag) {
         case (M_GE | M_UT):
            B_dm->flag = M_GE | M_LT;
            break;
         case (M_GE | M_LT):
            B_dm->flag = M_GE | M_UT;
            break;
         case M_GE:
            B_dm->flag = M_GE;
            break;
         default:
            fn_DisplayError(".../mathlib.c/TransposeSquare():  (1) If the R input matrix is symmetric, no transpose is needed. (2) If triangular, it may be unnecessary so that no code is written for a triangular case");
      }
   }
   else {
      // ????? NOT tested yet.  2/27/03.  ????????????
      for (_j=0; _j<_n; _j++)
         for (_i=_j+1; _i<_n; _i++) {
            //Off-diagonal elements.
            tmpd =  A[_j*_n+_i];
            A[_j*_n+_i] = A[_i*_n+_j];
            A[_i*_n+_j] = tmpd;
         }

      switch (A_dm->flag) {
         case (M_GE | M_UT):
            B_dm->flag = M_GE | M_LT;
            break;
         case (M_GE | M_LT):
            B_dm->flag = M_GE | M_UT;
            break;
         case M_GE:
            B_dm->flag = M_GE;
            break;
         default:
            fn_DisplayError(".../mathlib.c/TransposeSquare():  (1) If the R input matrix is symmetric, no transpose is needed. (2) If triangular, it may be unnecessary so that no code is written for a triangular case");
      }
   }
   /**/
}

void TransposeRegular(TSdmatrix *B_dm, const TSdmatrix *A_dm) {
/*     // Transposes the m-by-n matrix A_dm to the n-by-m matrix B_dm.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  B_dm: n-by-m general matrix (memory already allocated outside this function).   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n general matrix to be transposed.   ansi-c*/
   int _i, _j, _m, _n;
   double *A, *B;


/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if ( !B_dm || !A_dm ) fn_DisplayError(".../mathlib.c/TransposeRegular():  Input matrices must be created (memory-allocated)");
   else if ( !(A_dm->flag & M_GE) ) fn_DisplayError(".../mathlib.c/TransposeRegular():  (1) Input matrix A_dm must be given values. (2) A_dm must be a general (M_GE) matrix");
   else {
      _m = A_dm->nrows;
      _n = A_dm->ncols;
      A = A_dm->M;
      B = B_dm->M;
   }
   if ( (_m != B_dm->ncols) || (_n != B_dm->nrows) ) fn_DisplayError(".../mathlib.c/TransposeRegular(): Dimension of B_dm must be compatible with those of tranposed A_dm");


/*     //=== Transposing the regular matrix A_dm   ansi-c*/
   if (_m<_n) {
      for (_j=0; _j<_m; _j++)
         for (_i=_j+1; _i<_m; _i++) {
/*              //Off-diagonal elements in the square part.   ansi-c*/
            B[mos(_j, _i, _n)] = A[mos(_i, _j, _m)];
            B[mos(_i, _j, _n)] = A[mos(_j, _i, _m)];
         }

      for (_i=_m-1; _i>=0; _i--)
         B[mos(_i,_i,_n)] = A[mos(_i, _i, _m)];    /*  Diagonal elements.   ansi-c*/

      for (_j=_m; _j<_n; _j++)
         for (_i=0; _i<_m; _i++)
            B[_i*_n+_j] = A[_j*_m+_i];    /*  Off-diagonal elements in the trapozoidal part.   ansi-c*/

      B_dm->flag = M_GE;
   }
   else {
      for (_j=0; _j<_n; _j++)
         for (_i=_j+1; _i<_n; _i++) {
/*              //Off-diagonal elements in the square part.   ansi-c*/
            B[_i*_n+_j] = A[_j*_m+_i];
            B[_j*_n+_i] = A[_i*_m+_j];
         }

      for (_i=0; _i<_n; _i++)
         B[mos(_i,_i,_n)] = A[mos(_i,_i,_m)];    /*  Diagonal elements.   ansi-c*/

      for (_i=_n; _i<_m; _i++)
         for (_j=0; _j<_n; _j++)
            B[_i*_n+_j] = A[_j*_m+_i];   /*  Off-diagonal elements in the trapozoidal part.   ansi-c*/

      B_dm->flag = M_GE;
   }
}
/*  //---   ansi-c*/
TSdmatrix *tz_TransposeRegular(TSdmatrix *B_dm, const TSdmatrix *A_dm)
{
/*     // Transposes the m-by-n matrix A_dm to the n-by-m matrix B_dm.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If B_dm==NULL, B_dm (n-by-m general matrix) is created (memory allocated) and returned (thus, the memory must be destroyed outside this function).   ansi-c*/
/*     //  If B_dm!=NULL, B_dm (n-by-m general matrix)'s memory has already been allocated outside this function and the same B_dm will be returned.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm: m-by-n general matrix to be transposed.   ansi-c*/
   int _i, _j, _m, _n;
   double *A, *B;


/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if (!A_dm)  fn_DisplayError(".../mathlib.c/TransposeRegular():  Input matrix A_dm must be created (memory-allocated)");
   if ( !(A_dm->flag & M_GE) )  fn_DisplayError(".../mathlib.c/TransposeRegular():  (1) Input matrix A_dm must be given values. (2) A_dm must be a general (M_GE) matrix");
   _m = A_dm->nrows;
   _n = A_dm->ncols;
   A = A_dm->M;
/*     //   ansi-c*/
   if (!B_dm)  B_dm = CreateMatrix_lf(_n, _m);
   else
      if ( (_m != B_dm->ncols) || (_n != B_dm->nrows) ) fn_DisplayError(".../mathlib.c/TransposeRegular(): Dimension of B_dm must be compatible with those of tranposed A_dm");
   B = B_dm->M;

/*     //=== Transposing the regular matrix A_dm   ansi-c*/
   if (_m<_n) {
      for (_j=0; _j<_m; _j++)
         for (_i=_j+1; _i<_m; _i++) {
/*              //Off-diagonal elements in the square part.   ansi-c*/
            B[mos(_j, _i, _n)] = A[mos(_i, _j, _m)];
            B[mos(_i, _j, _n)] = A[mos(_j, _i, _m)];
         }

      for (_i=_m-1; _i>=0; _i--)
         B[mos(_i,_i,_n)] = A[mos(_i, _i, _m)];    /*  Diagonal elements.   ansi-c*/

      for (_j=_m; _j<_n; _j++)
         for (_i=0; _i<_m; _i++)
            B[_i*_n+_j] = A[_j*_m+_i];    /*  Off-diagonal elements in the trapozoidal part.   ansi-c*/

      B_dm->flag = M_GE;
   }
   else {
      for (_j=0; _j<_n; _j++)
         for (_i=_j+1; _i<_n; _i++) {
/*              //Off-diagonal elements in the square part.   ansi-c*/
            B[_i*_n+_j] = A[_j*_m+_i];
            B[_j*_n+_i] = A[_i*_m+_j];
         }

      for (_i=0; _i<_n; _i++)
         B[mos(_i,_i,_n)] = A[mos(_i,_i,_m)];    /*  Diagonal elements.   ansi-c*/

      for (_i=_n; _i<_m; _i++)
         for (_j=0; _j<_n; _j++)
            B[_i*_n+_j] = A[_j*_m+_i];   /*  Off-diagonal elements in the trapozoidal part.   ansi-c*/

      B_dm->flag = M_GE;
   }

   return (B_dm);
}


void SUtoGE(TSdmatrix *x_dm) {
/*     //Output: x_dm (nrows<=ncols) becomes a general matrix in addition to being upper symmetric.   ansi-c*/
/*     //Input: x_dm (nrows<=ncols) is upper symmetric.   ansi-c*/
/*     // We do not check if x_dm is upper symmetric because we assume the calling function checks this before this function is called.   ansi-c*/
   int nrows, ncols, _i, _j;
   if (!x_dm)  fn_DisplayError(".../mathlib.c/SUtoGE(): Input upper symmetric matrix must be created (memory-allocated)");
   else if ( !(x_dm->flag & M_SU) )  fn_DisplayError(".../mathlib.c/SUtoGE(): Input matrix must be (1) upper symmetric and (2) given legal values");
   else if ( (nrows=x_dm->nrows) != (ncols=x_dm->ncols) )   fn_DisplayError(".../mathlib.c/SUtoGE(): Upper symmetric matrix must be square.");
   else {
      for (_j=0; _j<nrows; _j++)
         for (_i=_j+1; _i<nrows; _i++)
            x_dm->M[mos(_i, _j, nrows)] = x_dm->M[mos(_j, _i, nrows)];    /*  Off-diagonal elements in the square part.   ansi-c*/
      x_dm->flag = M_GE | M_SL | M_SU;
   }
}

void SLtoGE(TSdmatrix *x_dm) {
/*     //Output: x_dm becomes a general square matrix in addition to being lower symmetric.   ansi-c*/
/*     //Input: x_dm is lower symmetric.   ansi-c*/
/*     // We do not check if x_dm is upper symmetric because we assume the calling function checks this before this function is called.   ansi-c*/
   int nrows, ncols, _i, _j;
   if (!x_dm)  fn_DisplayError(".../mathlib.c/SLtoGE(): Input lower symmetric matrix must be created (memory-allocated)");
   else if ( !(x_dm->flag & M_SL) )  fn_DisplayError(".../mathlib.c/SLtoGE(): Input matrix must be (1) lower symmetric and (2) given legal values");
   else if ( (nrows=x_dm->nrows) != (ncols=x_dm->ncols) )   fn_DisplayError(".../mathlib.c/SLtoGE(): The lower symmetric matrix must be sqaure");
   else {
      for (_j=0; _j<nrows; _j++)
         for (_i=_j+1; _i<nrows; _i++)
            x_dm->M[mos(_j, _i, nrows)] = x_dm->M[mos(_i, _j, nrows)];   /*  Off-diagonal elements in the square part.   ansi-c*/
      x_dm->flag = M_GE | M_SU | M_SL;
   }
}

double SumVector(TSdvector *x_dv) {
   int _i;
   double sum = 0.0,   /*  Cumulative.   ansi-c*/
          *v;
/*     //double *ptrcnt, *endptr;   ansi-c*/

/*     //=== This option may not speed up.   ansi-c*/
/*     // if ( !x_dv ) fn_DisplayError(".../mathlib.c/SumVector(): Input vector must be created (memory-allocated)");   ansi-c*/
/*     // else endptr = (ptrcnt = x_dv->v) + x_dv->n;   ansi-c*/
/*     // for ( ; ptrcnt<endptr; ptrcnt++ ) sum += *ptrcnt;   ansi-c*/

   if ( !x_dv || !x_dv->flag ) fn_DisplayError(".../mathlib.c/SumVector(): input vector must be (a) created (memory-allocated) and (b) assigned legal values");
   for (_i=x_dv->n-1, v=x_dv->v; _i>=0; _i--)  sum += v[_i];

   return( sum );
}

/*  //double SumVector(TSdvector *x_dv) {   ansi-c*/
/*  //   int _i, _n;   ansi-c*/
/*  //   double sum;   ansi-c*/
/*  //   double *v;   ansi-c*/

/*  //   if ( !x_dv ) fn_DisplayError(".../mathlib.c/SumVector(): Input vector must be created (memory-allocated)");   ansi-c*/
/*  //   else {   ansi-c*/
/*  //      v = x_dv->v;   ansi-c*/
/*  //      _n = x_dv->n;   ansi-c*/
/*  //   }   ansi-c*/

/*  //   sum = v[0];   ansi-c*/
/*  //   for ( _i=1; _i<_n; _i++ ) sum += v[_i];   ansi-c*/
/*  //   return( sum );   ansi-c*/
/*  //}   ansi-c*/


double MinVector(TSdvector *x_dv)
{
/*     //Input: no change for x_dv in this function.   ansi-c*/
   int _i, n;
   double minvalue;
   double *v;

   if (!x_dv || !x_dv->flag) fn_DisplayError(".../cstz.c/MinVector_lf():  Input vector x_dv must be (1) allocated memory and (2) assigned legal values");
   n = x_dv->n;
   v = x_dv->v;

   minvalue = v[0];
   for (_i=n-1; _i>0; _i--)
      if (v[_i]<minvalue)  minvalue = v[_i];

   return( minvalue );
}

double MaxVector(TSdvector *x_dv)
{
/*     //Input: no change for x_dv in this function.   ansi-c*/
   int _i, n;
   double maxvalue;
   double *v;

   if (!x_dv || !x_dv->flag) fn_DisplayError(".../cstz.c/MaxVector_lf():  Input vector x_dv must be (1) allocated memory and (2) assigned legal values");
   n = x_dv->n;
   v = x_dv->v;

   maxvalue = v[0];
   for (_i=n-1; _i>0; _i--)
      if (v[_i]>maxvalue)  maxvalue = v[_i];

   return( maxvalue );
}

int MaxVector_int(TSivector *x_iv)
{
/*     //Input: no change for x_iv in this function.   ansi-c*/
   int _i;
   int maxvalue;
   int *v;

   if (!x_iv || !x_iv->flag) fn_DisplayError(".../cstz.c/MaxVector_int():  Input vector x_iv must be (1) allocated memory and (2) assigned legal values");
   v = x_iv->v;

   maxvalue = v[0];
   for (_i=x_iv->n-1; _i>0; _i--)
      if (v[_i]>maxvalue)  maxvalue = v[_i];

   return( maxvalue );
}


void SumMatrix(TSdvector *x_dv, const TSdmatrix *X_dm, const char rc)
{
/*     //Outputs:   ansi-c*/
/*     //  x_dv: if rc=='R' or 'r', x_dv = sum of X_dm across rows; else or if rc=='C' or 'c', x_dv = sum of X_dm across columns.   ansi-c*/
   int _i, _n, _k, nrows, ncols, last;
   double sum;
   double *v, *M;


   if ( !x_dv || !X_dm || !X_dm->flag )  fn_DisplayError(".../mathlib.c/SumMatrix(): (a) output vector must be created (memory-allocated) and (b) input matrix must be created and assigned legal values");


   if (rc=='R' || rc=='r') {
      if ((_n=x_dv->n) != X_dm->ncols)  fn_DisplayError(".../mathlib.c/SumMatrix(): length of the output vector must match the number of columns of the input matrix when summing it across rows");
      v = x_dv->v;
      M = X_dm->M;
      for (_i=_n-1; _i>=0; _i--) {
         sum = 0.0;
         last = (_i+1)*(nrows=X_dm->nrows);
         for (_k=_i*nrows; _k<last; _k++)  sum +=M[_k];
         v[_i] = sum;
      }
   }
   else {
      if ((_n=x_dv->n) != X_dm->nrows)  fn_DisplayError(".../mathlib.c/SumMatrix(): length of the output vector must match the number of rows of the input matrix when summing it across columns");
      v = x_dv->v;
      M = X_dm->M;
      for (_i=_n-1; _i>=0; _i--) {
         sum = 0.0;
         last = _i + ((ncols=X_dm->ncols)-1)*_n;
         for (_k=_i; _k<=last; _k += _n)  sum +=M[_k];    /*  Must _k<=, NOT, _k<.   ansi-c*/
         v[_i] = sum;
      }
   }

   x_dv->flag = V_DEF;
}


void diagdv(TSdvector *x_dv, TSdmatrix *x_dm)
{
/*     //Extract the diagonal elements of x_dm to a vector.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  x_dv: nrows-by-1 vector.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_dm: nrows-by-ncols matrix.   ansi-c*/
   int _i, _n;
   double *v, *M;


/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if ( !x_dv || !x_dm ) fn_DisplayError(".../mathlib.c/diagdv(): Both the input vector and output matrix must be created (memory-allocated)");
   else {
      if (x_dm->nrows < x_dm->ncols) _n = x_dm->nrows;
      else _n = x_dm->ncols;
      v = x_dv->v;
      M = x_dm->M;
   }
   if ( _n != x_dv->n ) fn_DisplayError(".../mathlib.c/diagdv(): Dimensions of input vector and matrix must match");


   for (_i=0; _i<_n; _i++)  v[_i] = M[mos(_i,_i,_n)];
   x_dv->flag = V_DEF;
}
/*  //---   ansi-c*/
TSdmatrix *tz_DiagMatrix(TSdmatrix *X_dm, TSdvector *x_dv)
{
/*     //Converts a vector to a diagonal matrix with the diagonal elements being the input vector.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  X_dm: _n-by-_n diagonal matrix.   ansi-c*/
/*     //  If X_dm = NULL, then Xout_dm is allocated memory and exported (therefore, its memory will be freed outside this function0.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_dv: _n-by-1 vector.   ansi-c*/
   int _i, _n;
   double *v, *M;
   TSdmatrix *Xout_dm;


/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if ( !x_dv || !x_dv->flag)  fn_DisplayError(".../mathlib.c/tz_DiagMatrix(): the input vector must be (1) created (memory-allocated) and (2) given legal values");
   _n = x_dv->n;

   if (X_dm)
   {
      if ((_n != X_dm->nrows) || (_n != X_dm->ncols))  fn_DisplayError(".../mathlib.c/tz_DiagMatrix(): (1) the input matrix must be square; (2) dimensions of input vector and matrix must match");
      if (isdiagonalmatrix(X_dm))  Xout_dm = X_dm;
      else  fn_DisplayError(".../mathlib.c/tz_DiagMatrix(): the input matrix must be diagonal (M_UT | M_LT)");
   }
   else  Xout_dm = CreateConstantMatrix_lf(_n, _n, 0.0);

   M = Xout_dm->M;
   v = x_dv->v;

   for (_i=0; _i<_n; _i++)  M[mos(_i,_i,_n)] = v[_i];
   if ( !X_dm )  Xout_dm->flag = (M_GE | M_UT | M_LT);   /*  Diagonal (i.e., both lower and upper triangular).   ansi-c*/

   return (Xout_dm);
}


double tracefabs(TSdmatrix *x_dm)
{
/*     //Sum of absolute values of the diagonal elements of x_dm.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  y: double value.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_dm: nrows-by-ncols matrix.   ansi-c*/
   int _i, _n;
   double traceval=0.0,   /*  Cumulative.   ansi-c*/
          *M;


/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if (!x_dm) fn_DisplayError(".../mathlib.c/tracefabs(): The input matrix must be created (memory-allocated)");
   else {
      if (x_dm->nrows < x_dm->ncols) _n = x_dm->nrows;
      else _n = x_dm->ncols;
       M = x_dm->M;
   }


   for (_i=0; _i<_n; _i++)  traceval += fabs(M[mos(_i,_i,_n)]);
   return( traceval );
}
double tracelogfabs(TSdmatrix *x_dm)
{
/*     //Sum of logs of absolute values of the diagonal elements of the square x_dm or sum(log(diag(abs(x_dm)))).   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  y: double value.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_dm: nrows-by-ncols matrix.   ansi-c*/
   int _i, _n;
   double traceval=0.0,   /*  Cumulative.   ansi-c*/
          *M;

/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if (!x_dm || !x_dm->flag)  fn_DisplayError(".../mathlib.c/tracelogfabs(): The input matrix must be (1) created (memory-allocated) and (2) gvein legal values");
   else  M = x_dm->M;
   if ((_n = x_dm->nrows) != x_dm->ncols)   fn_DisplayError(".../mathlib.c/tracelogfabs(): The input matrix must be square");
   for (_i=square(_n)-1; _i>=0; _i -= _n+1)   traceval += log(fabs(M[_i]));

/*  //   if (!x_dm)  fn_DisplayError(".../mathlib.c/tracelogfabs(): The input matrix must be created (memory-allocated)");   ansi-c*/
/*  //   else if ( !x_dm->flag )  fn_DisplayError(".../mathlib.c/tracelogfabs(): The input matrix must be given legal values");   ansi-c*/
/*  //   else {   ansi-c*/
/*  //      if (x_dm->nrows < x_dm->ncols) _n = x_dm->nrows;   ansi-c*/
/*  //      else _n = x_dm->ncols;   ansi-c*/
/*  //      M = x_dm->M;   ansi-c*/
/*  //   }   ansi-c*/
/*  //   for (_i=0; _i<_n; _i++)  traceval += log(fabs(M[mos(_i,_i,_n)]));   ansi-c*/

   return( traceval );
}
double tracelog(TSdmatrix *x_dm)
{
/*     //Sum of logs of the diagonal elements of the square x_dm or sum(log(diag(x_dm))).   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  y: double value.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  x_dm: nrows-by-ncols matrix.   ansi-c*/
   int _i, _n;
   double traceval=0.0,   /*  Cumulative.   ansi-c*/
          *M;

/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if (!x_dm || !x_dm->flag)  fn_DisplayError(".../mathlib.c/tracelogfabs(): The input matrix must be (1) created (memory-allocated) and (2) gvein legal values");
   else  M = x_dm->M;
   if ((_n = x_dm->nrows) != x_dm->ncols)   fn_DisplayError(".../mathlib.c/tracelogfabs(): The input matrix must be square");
   for (_i=square(_n)-1; _i>=0; _i -= _n+1)   traceval += log(M[_i]);

   return( traceval );
}
/*  //---   ansi-c*/
double sumoflogvector(TSdvector *x_dv)
{
/*     //Output: sum(log(x_dv)).   ansi-c*/
   int _i;
   double sumlog = 0.0;
   double *v;

   if ( !x_dv || !x_dv->flag )  fn_DisplayError("mathlib.c/sumoflogvector(): Input vector x_dv must be (1) created and (2) given legal values");
   v = x_dv->v;
   for (_i=x_dv->n-1; _i>=0; _i--)  sumlog += log(v[_i]);

   return( sumlog );
}


TSdmatrix *tz_kron(TSdmatrix *C_dm, TSdmatrix *A_dm, TSdmatrix *B_dm)
{
/*     //C = kron(A, B), compatible with Matlab notation.   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  A_dm and B_dm: two real general matrices.   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  If C_dm == NULL, C_dm is created (memory allocated) and returned (thus, the memory must be destroyed outside this function).   ansi-c*/
/*     //  If C_dm != NULL, C_dm's memory has already been allocated outside this function and the same C_dm will be returned.   ansi-c*/
   int _i, _j, ma, na, mb, nb, ki, kj;
   TSdmatrix *Wmb_nb_dm = NULL;

/*     //=== Checking dimensions and memory allocation.   ansi-c*/
   if (!A_dm || !B_dm)  fn_DisplayError("mathlib.c/tz_kron():  Input matrices A_dm and B_dm must be created (memory-allocated)");
   if ( !(A_dm->flag & M_GE) || !(B_dm->flag & M_GE))  fn_DisplayError("mathlib.c/tz_kron():  "
                 "  (1) Input matrices A_dm and B_dm must be given values."
                 "  (2) A_dm and B_dm must be general (M_GE) matrices");
   ma = A_dm->nrows;
   na = A_dm->ncols;
   mb = B_dm->nrows;
   nb = B_dm->ncols;
   Wmb_nb_dm = CreateMatrix_lf(mb,nb);
/*     //   ansi-c*/
   if (!C_dm)  C_dm = CreateZeroMatrix_lf(ma*mb, na*nb);
   else
      if ( (C_dm->nrows != (ma*mb)) || (C_dm->ncols != (na*nb)) )
         fn_DisplayError("mathlib.c/tz_kron(): Dimension of C_dm must be compatible with those of A_dm and B_dm");

   for (_i=ma-1; _i>=0; _i--)
      for (_j=na-1; _j>=0; _j--)
      {
         ki = _i*mb;
         kj = _j*nb;
         ScalarTimesMatrix(Wmb_nb_dm, A_dm->M[mos(_i,_j,ma)], B_dm, 0.0);
         CopySubmatrix(C_dm, ki, kj, Wmb_nb_dm, 0, 0, mb, nb);
      }

/*     //===   ansi-c*/
   DestroyMatrix_lf(Wmb_nb_dm);


   return (C_dm);
}







/*  //=======================================================   ansi-c*/
/*  // Self-written routines.   ansi-c*/
/*  //=======================================================   ansi-c*/
void ergodicp(TSdvector *p_dv, TSdmatrix *P_dm) {
/*     // Computes the ergodic probabilities.  See Hamilton p.681.   ansi-c*/
/*     //   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //  p_dv: n-by-1 vector filled by ergodic probabilities p.   ansi-c*/
/*     //------------   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  P_dm: n-by-n Markovian transition matrix.  Elements in each column sum up to 1.0.   ansi-c*/

   int eigmaxindx,                       /*   Index of the column corresponding to the max eigenvalue.   ansi-c*/
       _i, _j, _n, errflag;
   double gpisum=0.0,
          eigmax, tmpd0;
   double *p_v=NULL,
          *absval_v=NULL,
          *evalr_v=NULL,
          *evali_v=NULL,
          *revecr_m=NULL,
          *reveci_m=NULL,
          *levecr_m=NULL,
          *leveci_m=NULL;
   TSdvector *absvals_dv=NULL;
   TSdzvector *vals_dzv=NULL;
   TSdzmatrix *rights_dzm=NULL, *lefts_dzm=NULL;


   if ( !p_dv || !P_dm || (p_dv->n != P_dm->nrows) || (P_dm->nrows != P_dm->ncols) ) fn_DisplayError(".../mathlib.c/ergodicp(): One of the two pointer arguments is not created (memory-allocated) or sizes of two pointer arguments do not match or input matrix P_dm is not square");
   else if ( !P_dm->flag || !(P_dm->flag & M_GE) )  fn_DisplayError(".../mathlib.c/ergodicp(): (1) R square input matrix (P_dm) must be given values. (2) P_dm must be a general (M_GE) matrix");
   else {
      _n = p_dv->n;
      absvals_dv = CreateVector_lf(_n);
      vals_dzv = CreateVector_dz(_n);
      rights_dzm = CreateMatrix_dz(_n, _n);
      InitializeConstantMatrix_lf(rights_dzm->imag, 0.0);   /*  Imaginary part must be initialized to zero.   ansi-c*/
   }


/*     //=== Obtains eigen values and vectors.   ansi-c*/
/*     //errflag = eigrgen_decomp(evalr_v, evali_v, revecr_m, reveci_m, levecr_m, leveci_m, cp_m, _n);   ansi-c*/
   errflag = eigrgen(vals_dzv, rights_dzm, lefts_dzm, P_dm);
   if (errflag<0) fn_DisplayError("/mathlib.c/eigrgen(): some element in input matrix P_dm has an illegal value");
   else if (errflag>0) fn_DisplayError("/mathlib.c/eigrgen(): the QR algorithm failed to compute all the eigenvalues and no eigenvectors have been computed");

/*     //=== Utilizes old notations because I have no time to polish this function.   ansi-c*/
   p_v = p_dv->v;
   absval_v = absvals_dv->v;
   evalr_v = vals_dzv->real->v;
   evali_v = vals_dzv->imag->v;
   revecr_m = rights_dzm->real->M;
   reveci_m = rights_dzm->imag->M;   /*  Imaginary part must be initialized to zero.   ansi-c*/


   for (_j=0; _j<_n; _j++) {
      if (!(evali_v[_j])) {    /*  No imaginary part (in other words, real solutions).   ansi-c*/
         eigmax = evalr_v[_j];
         eigmaxindx = _j;
         break;
      }
      else {
         eigmax = sqrt(square(evalr_v[_j])+square(evali_v[_j]));
         eigmaxindx = _j;
         break;
      }
   }
/*     //+   ansi-c*/
   for (_j++; _j<_n; _j++) {
      if (!(evali_v[_j]) && (evalr_v[_j] > eigmax)) {
            eigmax = evalr_v[_j];
            eigmaxindx = _j;
      }
      else if (evali_v[_j]) {
         tmpd0 = sqrt(square(evalr_v[_j])+square(evali_v[_j]));
         if (tmpd0 > eigmax) {
            eigmax = tmpd0;
            eigmaxindx = _j;
         }
      }
   }

   if (!(evali_v[eigmaxindx])) {
      for (_i=0;_i<_n;_i++) {
         absval_v[_i] = fabs(revecr_m[_i+_n*eigmaxindx]);
         gpisum += absval_v[_i];   /*   Sum over the eigmaxindx_th column.   ansi-c*/
      }
      tmpd0 = 1.0/gpisum;
      for (_i=0;_i<_n;_i++) p_v[_i] = absval_v[_i]*tmpd0;       /*   Normalized eigmaxindx_th column as ergodic probabilities.   ansi-c*/
   }
   else {
      for (_i=0;_i<_n;_i++) {
         absval_v[_i] = sqrt(square(revecr_m[_i+_n*eigmaxindx])+square(reveci_m[_i+_n*eigmaxindx]));
         gpisum += absval_v[_i];   /*   Sum over the eigmaxindx_th column.   ansi-c*/
      }
      tmpd0 = 1.0/gpisum;
      for (_i=0;_i<_n;_i++) p_v[_i] = absval_v[_i]*tmpd0;       /*   Normalized eigmaxindx_th column as ergodic probabilities.   ansi-c*/
   }


   p_dv->flag = V_DEF;
/*     //=== Frees up allocated memory belonging to this function.   ansi-c*/
   DestroyVector_lf(absvals_dv);
   DestroyVector_dz(vals_dzv);
   DestroyMatrix_dz(rights_dzm);
}



double *alloc_ergodp2(const double *cp_m, const int _n) {
/*     // Output:   ansi-c*/
/*     //   p_v: n-by-1 vector of ergodic probabilities p.   ansi-c*/
/*     //------------   ansi-c*/
/*     // cp_m: n-by-n Markovian transition matrix.   ansi-c*/
/*     // _n:  the order of cp_m.   ansi-c*/
/*     //   ansi-c*/
/*     // Compute the ergodic probabilities.  See Hamilton p.681.   ansi-c*/

   int eigmaxindx,                       /*   Index of the column corresponding to the max eigenvalue.   ansi-c*/
       _i, _j, errflag;
   double gpisum=0.0,
          eigmax, tmpd0;
   double *p_v=NULL,             /*  @@Will be freed outside this function.@@  D: n-by-1.  Erogodic probabilties.   ansi-c*/
          *absval_v=NULL,
          *evalr_v=NULL,
          *evali_v=NULL,
          *revecr_m=NULL,
          *reveci_m=NULL,
          *levecr_m=NULL,
          *leveci_m=NULL;

/*     //=== Allocates memory.   ansi-c*/
   p_v = tzMalloc(_n, double);
   absval_v = tzMalloc(_n, double);
   evalr_v = tzMalloc(_n, double);
   evali_v = tzCalloc(_n, double);    /*  Imaginary part must be initialized to zero.   ansi-c*/
   revecr_m = tzMalloc(square(_n), double);
   reveci_m = tzCalloc(square(_n), double);   /*  Imaginary part must be initialized to zero.   ansi-c*/

/*     //=== Obtains eigen values and vectors.   ansi-c*/
   errflag = eigrgen_decomp(evalr_v, evali_v, revecr_m, reveci_m, levecr_m, leveci_m, cp_m, _n);
   if (errflag<0) fn_DisplayError("/mathlib.c/eigrgen_decomp(): some element in input matrix has an illegal value");
   else if (errflag>0) fn_DisplayError("/mathlib.c/eigrgen_decomp(): the QR algorithm failed to compute all the eigenvalues and no eigenvectors have been computed");

   for (_j=0; _j<_n; _j++) {
      if (!(evali_v[_j])) {    /*  No imaginary part (in other words, real solutions).   ansi-c*/
         eigmax = evalr_v[_j];
         eigmaxindx = _j;
         break;
      }
      else {
         eigmax = sqrt(square(evalr_v[_j])+square(evali_v[_j]));
         eigmaxindx = _j;
         break;
      }
   }
/*     //+   ansi-c*/
   for (_j++; _j<_n; _j++) {
      if (!(evali_v[_j]) && (evalr_v[_j] > eigmax)) {
            eigmax = evalr_v[_j];
            eigmaxindx = _j;
            break;
      }
      else if (evali_v[_j]) {
         tmpd0 = sqrt(square(evalr_v[_j])+square(evali_v[_j]));
         if (tmpd0 > eigmax) {
            eigmax = tmpd0;
            eigmaxindx = _j;
            break;
         }
      }
   }

   if (!(evali_v[eigmaxindx])) {
      for (_i=0;_i<_n;_i++) {
         absval_v[_i] = fabs(revecr_m[_i+_n*eigmaxindx]);
         gpisum += absval_v[_i];   /*   Sum over the eigmaxindx_th column.   ansi-c*/
      }
      tmpd0 = 1.0/gpisum;
      for (_i=0;_i<_n;_i++) p_v[_i] = absval_v[_i]*tmpd0;       /*   Normalized eigmaxindx_th column as ergodic probabilities.   ansi-c*/
   }
   else {
      for (_i=0;_i<_n;_i++) {
         absval_v[_i] = sqrt(square(revecr_m[_i+_n*eigmaxindx])+square(reveci_m[_i+_n*eigmaxindx]));
         gpisum += absval_v[_i];   /*   Sum over the eigmaxindx_th column.   ansi-c*/
      }
      tmpd0 = 1.0/gpisum;
      for (_i=0;_i<_n;_i++) p_v[_i] = absval_v[_i]*tmpd0;       /*   Normalized eigmaxindx_th column as ergodic probabilities.   ansi-c*/
   }


/*     //=== Frees up allocated memory.   ansi-c*/
   if (absval_v) swzFree(absval_v);
   if (evalr_v) swzFree(evalr_v);
   if (evali_v) swzFree(evali_v);
   if (revecr_m) swzFree(revecr_m);
   if (reveci_m) swzFree(reveci_m);
   if (levecr_m) swzFree(levecr_m);
   if (leveci_m) swzFree(leveci_m);

   return (p_v);
}




/**
void eig_rgen_all(double *eval_v, double *evec_m, const double *x_m, const int _n) {
   // Outputs (dependent on MATLAB C math library):  NEED to be fixed about eval_v or mxval_d, which may be *complex*.  10/13/02
   //   eval_v:  n-by-1 eigenvalues;
   //   evec_m: n-by-n corresponding eigenvectors column by column.
   //------------
   // Inputs:
   //   x_m:  _n-by_n real general (non-symmetric) matrix.
   //
   // Eigenanalysis of real general (non-symmetric) square matrix with all eigenvalues and eigenvectors.

   #ifdef MATLABCMATHLIBRARY  //Matlab dependent code.
      mxArray *mxval_d=NULL, *mxvec_m=NULL,   // @@Must be freed in this function.@@  m: n-by-n eigvector matrix; d: n-by-n eigvalue diagonal.
            *mx_m=NULL;                     // @@Must be freed in this function.@@
      double *mxval_d_p;       // _p: pointer to the corresponding mxArray whose name occurs before _p.
      int ki;


      mx_m = mlfDoubleMatrix(_n, _n, x_m, NULL);
      mxvec_m = mlfEig(&mxval_d, mx_m, NULL, NULL);

      memcpy(evec_m, mxGetPr(mxvec_m), square(_n)*sizeof(double));
      //+
      mxval_d_p = mxGetPr(mxval_d);
      for (ki=0; ki<_n; ki++) eval_v[ki] = mxval_d_p[_n*ki+ki];   // Note that n*ki+ki refers to a diagonal location in the n-by-n matrix.

      //=== Frees up allocated mxArray.
      mxDestroyArray(mxval_d);
      mxDestroyArray(mxvec_m);
      mxDestroyArray(mx_m);
  #endif
}


double *fn_ergodp2(const double *cp_m, const int _n) {
   // Output:
   //   p_v: n-by-1 vector of ergodic probabilities p.
   //------------
   // cp_m: n-by-n Markovian transition matrix.
   // _n:  the order of cp_m.
   //
   // Compute the ergodic probabilities.  See Hamilton p.681.

   int eigmaxindx,                      // Index of the column corresponding to the max eigenvalue.
       ki;
   double gpisum=0.0,
          eigmax, tmpd0,
          *p_v,        // @@Will be freed outside this function.@@  D: n-by-1.  Erogodic probabilties.
          *eval_v, *evec_m;             // @@Must be freed in this function.@@  D: n-by-1 and n-by-n.   Eigenvalues and eigenvectors.

   //=== Allocates memory.
   p_v = tzMalloc(_n, double);
   eval_v = tzMalloc(_n, double);
   evec_m = tzMalloc(square(_n), double);


   eig_rgen_all(eval_v, evec_m, cp_m, _n);
   eigmax = *eval_v;
   eigmaxindx = 0;
   if (_n>1) {
      for (ki=1;ki<_n;ki++) {
         if (eval_v[ki] > eigmax) {
            eigmax=eval_v[ki];
            eigmaxindx=ki;
         }                           // Note that n*ki+ki refers to a diagonal location in the n-by-n matrix.
      }
   }
   for (ki=0;ki<_n;ki++) gpisum += evec_m[_n*eigmaxindx+ki];  // Sum over the eigmaxindx_th column.
   tmpd0 = 1.0/gpisum;
   for (ki=0;ki<_n;ki++) p_v[ki] = evec_m[_n*eigmaxindx+ki]*tmpd0;      // Normalized eigmaxindx_th column as ergodic probabilities.

   //=== Frees up allocated memory.
   swzFree(eval_v);
   swzFree(evec_m);


   return p_v;
}
/**/


/**      //Iskander's code.
int eiggen(double *a, int n, double *dr, double *di, double *vr, double *vi) {
   unsigned char msg[101];
   int lwork = -1, info = 0;
   double *work, work1;
   char *jobvr = vr?"V":"N";
   register int i, j;

   // Query dsyev_ on the value of lwork
   dgeev_("N",jobvr,&n,a,&n,dr,di,NULL,&n,vr,&n,&work1,&lwork,&info);

   if (info < 0) {
      sprintf(msg, "Input %d to dgeev_ had an illegal value",-info);
      mexWarnMsgTxt(msg);
      return(info);
   }

   lwork = (int)(work1);
   work = mxCalloc(lwork,sizeof(double));
   dgeev_("N",jobvr,&n,a,&n,dr,di,NULL,&n,vr,&n,work,&lwork,&info);
   mxFree(work);

   if (info < 0) {
      sprintf(msg, "Input %d to dgeev_ had an illegal value",-info);
      mexWarnMsgTxt(msg);
      return(info);
   }

   for (i=0; i<n-1; i++)
      if (di[i] && (di[i]==-di[i+1]))
         for (j=0; j<n; j++) {
            vi[(i+1)*n+j] = -(vi[i*n+j]=vr[(i+1)*n+j]);
            vr[(i+1)*n+j] = vr[i*n+j];
         }

   if (info > 0) {
      sprintf(msg,"The QR algorithm failed to compute all the eigenvalues,\n"
         "and no eigenvectors have been computed, but elements D(%d:N) contain\n"
         "those eigenvalues which have converged.",info+1);
      mexWarnMsgTxt(msg);
      return(info);
   }

   return(info);
}
/**/


/**  SAVE (Dan's code).  Transpose the square matrix.
for (_i=0: _i<_n; _i++)
   for (_j=i+1; _j<_n; _j++) {
      tmp=A[_i+_j*_n];
      A[i+j*n]=A[j+i*n];
      A[j+i*n] = tmp;
   }
/**/

/**  SAVE (Dan's code).  Transpose the regular matrix.
memcpy(tmp, A, ....);
for (_i=0: _i<_n; _i++)
   for (_j=i+1; _j<m; _j++) {
      A[i+j*..] = tmp[j+i*..];
   }
/**/




/**
void VectorTimesSelf(TSdmatrix *C_dm, const TSdvector *a_dv, const double _alpha, const double _beta, const char ul) {
   //No Lapack -- my own function.
   //Output is C and all other arguments are inputs.
   //Computes C = alpah*a*a' + beta*C where
   //  a is m-by-1,
   //  C is m-by-m symmetric matrix,
   //  alpha: a double scalar,
   //  beta: a double scalar,
   //  ul: if == 'u' or 'U', elements in C are stored only in the upper part; otherwise, C is stored only in the lower part.
   int _i, _j, _m, _n;
   double *v, *M;

   if ( !C_dm || !a_dv ) fn_DisplayError(".../mathlib.c/VectorTimesSelf(): At least one of the pointer arguments is not created (memory-allocated)");
   else {
      v = a_dv->v;
      M = C_dm->M;
      _m = C_dm->nrows;
      _n = C_dm->ncols;
   }

   if ( (_m != a_dv->n) || (_m != _n) ) fn_DisplayError(".../mathlib.c/VectorTimesSelf(): (1) Input matrix must square and (2) its size must match the dimension of the input vector");
   else {
      if ( (ul == 'u') || (ul == 'U') )  {
         if ( _alpha==1.0 ) {
            if ( _beta==1.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] += v[_i] * v[_j];
                  }
               }
            }
            else if ( _beta==0.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] = v[_i] * v[_j];
                  }
               }
            }
            else {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] = v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                  }
               }
            }
         }
         else {
            if ( _beta==1.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] += _alpha * v[_i] * v[_j];
                  }
               }
            }
            else if ( _beta==0.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] = _alpha * v[_i] * v[_j];
                  }
               }
            }
            else {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=0; _i<=_j; _i++ ) {
                     M[mos(_i, _j, _m)] = _alpha* v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                  }
               }
            }
         }
      }
      else {
         if ( _alpha==1.0 ) {
            if ( _beta==1.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] += v[_i] * v[_j];
                  }
               }
            }
            else if ( _beta==0.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] = v[_i] * v[_j];
                  }
               }
            }
            else {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] = v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                  }
               }
            }
         }
         else {
            if ( _beta==1.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] += _alpha * v[_i] * v[_j];
                  }
               }
            }
            else if ( _beta==0.0 ) {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] = _alpha * v[_i] * v[_j];
                  }
               }
            }
            else {
               for ( _j=0; _j<_m; _j++ ) {
                  for ( _i=_j; _i<_m; _i++ ) {
                     M[mos(_i, _j, _m)] = _alpha* v[_i] * v[_j] + _beta*M[mos(_i, _j, _m)];
                  }
               }
            }
         }
      }
   }
}
/**/





/**
void swap_lf(double *a, double *b) {
   double tmpd0;

   tmpd0 = *a;
   *a = *b;
   *b = tmpd0;
}
/**/

/**  Creates zeros of the matrix of a given size.
TSdmatrix *CreateZeros_lf(int nrows, int ncols) {
   int _i;
   TSdmatrix *x_im=CreateMatrix_lf(nrows, ncols);
   for (_i=nrows*ncols-1; _i>=0; _i--)
      x_im->M[_i] = 0.0;
   return(x_im);
}
TSdmatrix *CreateIdentity_lf(int nrows, int ncols) {
   int _i;
   TSdmatrix *x_im=CreateMatrix_lf(nrows, ncols);
   for (_i=nrows*ncols-1; _i>=0; _i--)
      x_im->M[_i] = 0.0;
   if (nrows<=ncols)
      for (_i=square(nrows)-1; _i>=0; _i -= nrows+1)
         x_im->M[_i] = 1.0;
   else
      for (_i=(ncols-1)*(nrows+1); _i>=0; _i -= nrows+1)
         x_im->M[_i] = 1.0;
   return(x_im);
}
/**/



/**
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void permute_matrix(double *a, int n, int *indx) {
   double *b;
   int nn=n*n;
   register int i;
   b = swzCalloc(nn,sizeof(double));
   memcpy(b, a, nn*sizeof(double));
   for (i=0; i<nn; i++, a++)
      *a = b[indx[i%n]+indx[i/n]*n];
}

int main() {
   double a[9]={1,2,3,4,5,6,7,8,9};
   int indx[3]={1,2,0};
   permute_matrix(a,3,indx);
   return 0;
}
/**/
