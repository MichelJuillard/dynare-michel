/*******************************************************************
 * [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,psi,pi,div)
 *
 * System given as
 *         g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
 * with z an exogenous variable process and eta being endogenously determined
 * one-step-ahead expectational errors.  Returned system is
 *        y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
 * If z(t) is i.i.d., the last term drops out.
 * If div or stake is omitted from argument list, a div>1 or stake>1 is calculated.
 * eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
 * existence only with not-serially correlated z(t); eu=[-2,-2] for coincident zeros.
 *
 * g0, g1:  n-by-n matrices.
 * c:  n-by-1 constant terms.
 * z(t):  m-by-1 vector of exogenous residuals where m < n.
 * psi:  n-by-m matrix.
 * eta(t):  h-by-1 vector of expectational (endogenous) errors.
 * pi:  n-by-h matrix.
 * div: a real number dividing stable and unstable roots..  If < 1.0, a div>1.0 is calculated mechanically.
 *-------
 * G1 or Theta_dm:  n-by-n matrices.
 * C:  n-by-1 vector of constant terms.
 * impact:  n-by-m matrix.
 * gev:  n-by-2 z vector of stacked generalized eigenvalues where gev(;,2) ./ gev(:,1) = eig(g0, g1).
 * ywt:  n-by-nunstab z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.
 * fmat: nunstab-by-nunstab z matrix where nunstab is the number of non-stable roots.
 * fwt:  nunstab-by-m z matrix.
 *
 * 1996 MATLAB algorithm by Christopher Sims
 * 2002 Mex implementation by Iskander Karibzhanov
 * 2004 Modified to C function by Tao Zha (April), correcting a few bugs of Iskander.
 * 03/01/06  Another modification by Tao Zha, to be consistent with the CAS 3/10/04 correction.
 *
 * Note: Iskander is transforming g0 and g1 to complex matrices and uses zgges() as a qz decomposition.
 *       This is really wasting efficiency.  One should keep g0 and g1 as real matrices and use
 *       dgges() as a qz decomposition.  I don't have time to overhaul this at this point.  04/20/04, T. Zha.
 * Note: 02/22/06.  I take the above note back.  According to DW, it is easy to *order* the
 *       the generalized eigenvalues by using the complex g0 and g1.  In principle, one could
 *       order the roots using the real qz decomposition on real matrices g0 and g1.  But so far
 *       Dan has found it a pain to do it.  Perhaps we should read the MKL Lapack manual more
 *       carefully at a later point.
********************************************************************/
#include "gensys.h"

#include "modify_for_mex.h"

/*  //----- NOTE: We can't replace MKL_Complex16 with a different name because the Intel Lapack uses MKL_Complex16.   ansi-c*/
/*  //-----       The only way to do this is to overhaul the code and put a wrapper function on each Intel Lapack function.   ansi-c*/
static int selctg(MKL_Complex16 *alpha, MKL_Complex16 *beta);
static int qz(MKL_Complex16 *a, MKL_Complex16 *b, MKL_Complex16 *q, MKL_Complex16 *z, int n);
static MKL_Complex16* CreateComplexMatrix5RealMatrix(TSdmatrix *X_dm);
static MKL_Complex16* CreateComplexMatrix5RealVector(TSdvector *x_dv);
static void ComplexMatrix2RealMatrix(TSdmatrix *Y_dm, MKL_Complex16 *Z);
static void ComplexMatrix2RealVector(TSdvector *y_dv, MKL_Complex16 *Z);
static TSdzmatrix *SubComplexMatrix2Zmatrix(TSdzmatrix *X_dzm, MKL_Complex16 *Z, const int nrowsforZ, const int _m, const int _n);
static void copy_eigenvalues(TSdzmatrix *Gev_dzm, MKL_Complex16 *a, MKL_Complex16 *b);
static int compute_svd(MKL_Complex16 *a, MKL_Complex16 **u, double **d, MKL_Complex16 **v, int m, int n);
static int compute_norm(MKL_Complex16 *a, double **d, int m, int n);
/*  //--- 03/01/06 TZ.  Commented out to be consistent with the CAS 3/10/04 correction.   ansi-c*/
/*  // static int compute_normx(MKL_Complex16 *a, MKL_Complex16 *b, MKL_Complex16 *zwt, MKL_Complex16 *ueta, double **normx, int nunstab, int psin, int n, int bigev);   ansi-c*/
static void cblas_zdupe(int m, int n, MKL_Complex16 *a, int lda, MKL_Complex16 *b, int ldb);
static void cblas_zdscali(int n, double *a, int lda, MKL_Complex16 *b, int ldb);
static void cblas_zdscale(int n, double *a, int lda, MKL_Complex16 *b, int ldb);
static void cblas_zdpsb(int m, int n, MKL_Complex16 *a, int lda, MKL_Complex16 *b, int ldb, MKL_Complex16 *c, int ldc);
/*  //   ansi-c*/
static void InitializeConstantMLK_Complex16(MKL_Complex16 *x_clx,  const int _n, const double c);
static void InitializeConstantDouble(double *x_p,  const int _n, const double c);
static void ConverteZeroSquareMatrix2RealDiagonalMLK_Complex16(MKL_Complex16 *x_pc,  const int _n, const double c);


TSgensys *CreateTSgensys(TFlinratexp *func, const int _n, const int _m, const int _k, const double div)
{
/*     //_n is the number of stacked variables (endogenous, Lagurangian multiplier, expected multiplier, etc.).   ansi-c*/
/*     //_m is the number of exogenous shocks.   ansi-c*/
/*     //_k is the number of expectational errors.   ansi-c*/
/*     //div is the dividing number to determine what constitutes an unstable root.  If div<1.0, a div>1.0 is calculated mechanically.   ansi-c*/
   TSgensys *gensys_ps = tzMalloc(1, TSgensys);

/*     //=== Output arguments.   ansi-c*/
   gensys_ps->Theta_dm = CreateMatrix_lf(_n, _n);   /*  n-by-n.   ansi-c*/
   gensys_ps->c_dv = CreateVector_lf(_n);    /*  n-by-1.   ansi-c*/
   gensys_ps->Impact_dm = CreateMatrix_lf(_n, _m);   /*  n-by-m.   ansi-c*/
   gensys_ps->Fmat_dzm = (TSdzmatrix *)NULL;    /*  nunstab-by-nunstab z matrix.  Initialized to NULL and will be dynamically allocated whenever gensys() is called.   ansi-c*/
   gensys_ps->Fwt_dzm = (TSdzmatrix *)NULL;     /*  nunstab-by-m z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.   ansi-c*/
   gensys_ps->Ywt_dzm = (TSdzmatrix *)NULL;     /*  n-by-nunstab z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.   ansi-c*/
   gensys_ps->Gev_dzm = CreateMatrix_dz(_n, 2);   /*  n-by-2 z matrix of possible complex numbers.   ansi-c*/
   gensys_ps->eu_iv = CreateConstantVector_int(2, 0);    /*  2-by-1.   ansi-c*/

/*     //=== Function itself.   ansi-c*/
   gensys_ps->gensys = func;

/*     //=== Input arguments.   ansi-c*/
   gensys_ps->G0_dm = CreateConstantMatrix_lf(_n, _n, 0.0);   /*  n-by-n.   ansi-c*/
   gensys_ps->G0_dm->flag = M_GE;
   gensys_ps->G1_dm = CreateConstantMatrix_lf(_n, _n, 0.0);   /*  n-by-n.   ansi-c*/
   gensys_ps->G1_dm->flag = M_GE;
   gensys_ps->c0_dv = CreateConstantVector_lf(_n, 0.0);   /*  n-by-1.   ansi-c*/
   gensys_ps->Psi_dm = CreateConstantMatrix_lf(_n, _m, 0.0);  /*  n-by-m.   ansi-c*/
   gensys_ps->Psi_dm->flag = M_GE;
   gensys_ps->Pi_dm = CreateConstantMatrix_lf(_n, _k, 0.0);   /*  n-by-k where k is the number of expectational errors.   ansi-c*/
   gensys_ps->Pi_dm->flag = M_GE;
   gensys_ps->div = div;

   return (gensys_ps);
}
/*  //-------   ansi-c*/
TSgensys *DestroyTSgensys(TSgensys *gensys_ps)
{
   if (gensys_ps) {
/*        //=== Output arguments.   ansi-c*/
      DestroyMatrix_lf(gensys_ps->Theta_dm);   /*  n-by-n.   ansi-c*/
      DestroyVector_lf(gensys_ps->c_dv);    /*  n-by-1.   ansi-c*/
      DestroyMatrix_lf(gensys_ps->Impact_dm);   /*  n-by-m.   ansi-c*/
      DestroyMatrix_dz(gensys_ps->Fmat_dzm);    /*  nunstab-by-nunstab z matrix.  Initialized to NULL and will be dynamically allocated whenever gensys() is called.   ansi-c*/
      DestroyMatrix_dz(gensys_ps->Fwt_dzm);    /*  nunstab-by-m z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.   ansi-c*/
      DestroyMatrix_dz(gensys_ps->Ywt_dzm);     /*  n-by-nunstab z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.   ansi-c*/
      DestroyMatrix_dz(gensys_ps->Gev_dzm);   /*  n-by-2 z matrix of possible complex numbers.   ansi-c*/
      DestroyVector_int(gensys_ps->eu_iv);    /*  2-by-1.   ansi-c*/

/*        //=== Input arguments.   ansi-c*/
      DestroyMatrix_lf(gensys_ps->G0_dm);   /*  n-by-n.   ansi-c*/
      DestroyMatrix_lf(gensys_ps->G1_dm);   /*  n-by-n.   ansi-c*/
      DestroyVector_lf(gensys_ps->c0_dv);   /*  n-by-1.   ansi-c*/
      DestroyMatrix_lf(gensys_ps->Psi_dm);  /*  n-by-m.   ansi-c*/
      DestroyMatrix_lf(gensys_ps->Pi_dm);   /*  n-by-k where k is the number of expectational errors.   ansi-c*/

      free(gensys_ps);

      return ((TSgensys *)NULL);
   }
   else  return (gensys_ps);
}


/*  //--------------------------- For the function gensys_sims() ------------------------------------   ansi-c*/
static int fixdiv = 1, zxz = 0;
static double stake = 1.01;
static int nunstab = 0;
static MKL_Complex16 one, minusone, zero;

/* [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensysmkl(g0,g1,c,psi,pi,stake) */
/*  //void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {   ansi-c*/
int gensys_sims(TSgensys *gensys_ps, void *dummy_ps)
{
/*     //Returns 1 if successful and 0 for fatal errors (such as qz or svd fails or all roots are explosive). Added DW and TZ, 03/08/06.   ansi-c*/
   int tmpi;
   int n, psin, pin, nsquare, md, md1, i, bigev, bigev1;   /*  mds, bigevs,    ansi-c*/
   int *eu;
   int exist = 0, existx = 0, unique = 0;
/*     //=== Memory will be allocated to the following.   ansi-c*/
   double *deta = NULL, *deta1 = NULL, *norm = NULL; //*normx = NULL, *dz = NULL,  //03/01/06 TZ.  Commented out to be consistent with the CAS 3/10/04 correction.
   MKL_Complex16 *a = NULL, *b = NULL, *q = NULL, *z = NULL, *pi = NULL, *psi = NULL;
   MKL_Complex16 *tmat = NULL, *g0 = NULL, *g1 = NULL, *dummy = NULL, *tmatq = NULL, *c = NULL, *impact = NULL, *ab = NULL;
   MKL_Complex16 *fmat = NULL, *fwt = NULL, *ywt = NULL;
   MKL_Complex16 *etawt = NULL, *ueta = NULL, *veta = NULL, *etawt1 = NULL, *ueta1 = NULL, *veta1 = NULL;
        // *uz = NULL, *vz = NULL, *zwt = NULL, //03/01/06 TZ.  Commented out to be consistent with the CAS 3/10/04 correction.
   //--- Dimensions.
   n = gensys_ps->G0_dm->nrows;
   psin = gensys_ps->Psi_dm->ncols;
   pin = gensys_ps->Pi_dm->ncols;
   //--- Pointer.
   eu = gensys_ps->eu_iv->v;

   eu[0]=eu[1]=0;  //Must be initialized because gensys_ps->eu_iv->v may have values in repeated loops.

   //=== [a b q z]=qz(g0,g1);
   a = CreateComplexMatrix5RealMatrix(gensys_ps->G0_dm);
   b = CreateComplexMatrix5RealMatrix(gensys_ps->G1_dm);
   q = tzMalloc(nsquare=square(n), MKL_Complex16);
   z = tzMalloc(nsquare, MKL_Complex16);
   InitializeConstantMLK_Complex16(q,  nsquare, 0.0);
   InitializeConstantMLK_Complex16(z,  nsquare, 0.0);

   fixdiv = (gensys_ps->div < 1.0);
   stake = fixdiv ? 1.01 : gensys_ps->div;
   nunstab = 0;
   zxz = 0;

   if (qz(a, b, q, z, n)) {
      printf("WARNING:  QZ factorization failed.\n");
      tzDestroy(a);
      tzDestroy(b);
      tzDestroy(q);
      tzDestroy(z);
      eu[0] = 0;
      return 0;
   }

   nunstab /= 2;

   if (zxz) {
      printf("WARNING: Coincident zeros.  Indeterminacy and/or nonexistence.\n");
      eu[0] = eu[1] = -2;
      tzDestroy(a);
      tzDestroy(b);
      tzDestroy(q);
      tzDestroy(z);
      return 1;
   }
   copy_eigenvalues(gensys_ps->Gev_dzm, a,  b);

   one.real = 1.0;
   one.imag = 0.0;

   minusone.real = -1.0;
   minusone.imag = 0.0;

   zero.real = 0.0;
   zero.imag = 0.0;

   pi = CreateComplexMatrix5RealMatrix(gensys_ps->Pi_dm);
   //=============================================
   // Modified by DW and TZ to deal with the case where nunstab=0 (no explosive roots).  03/08/06.
   //=============================================
   if (nunstab)  //This branch belongs to original CAS code.
   {
      etawt = tzMalloc(tmpi=nunstab*pin, MKL_Complex16);
      InitializeConstantMLK_Complex16(etawt, tmpi, 0.0);    //Must be initialized to 0.0 in order to have legal values of this pointer.
      cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, nunstab, pin, n,
                  &one, q+n*(n-nunstab), n, pi, n, &zero, etawt, nunstab);
      if (compute_svd(etawt, &ueta, &deta, &veta, nunstab, pin)) {
         //Memory is now allocated to ueta, deta, and veta.
         printf("WARNING: SVD failed.\n");
         tzDestroy(pi);
         tzDestroy(ueta);
         tzDestroy(deta);
         tzDestroy(veta);
         tzDestroy(etawt);
         tzDestroy(a);
         tzDestroy(b);
         tzDestroy(q);
         tzDestroy(z);
         eu[0] = 0;
         return 0;
      }
      tzDestroy(etawt);
      md = nunstab<pin?nunstab:pin;
      bigev = md;
      for (i=0; i<md; i++)
         if (deta[i]<=REALSMALL) {
            bigev=i;
            break;
         }
      //------ 03/01/06 TZ: corrected code by CAS, 3/10/04.
      if ((eu[0]=(bigev >= nunstab))==0)  //DW & TZ, 03/08/06
      {
         tzDestroy(pi);
         tzDestroy(ueta);
         tzDestroy(deta);
         tzDestroy(veta);
         tzDestroy(a);
         tzDestroy(b);
         tzDestroy(q);
         tzDestroy(z);
         return 1;
      }
   }
   else   //DW & TZ.  03/08/06.  This is where we deal with the case when nunstab=0.
   {
      eu[0] = 1;  //Existence.
   }


   //---------------------------------
   // ueta = nunstab x bigev
   // deta = bigev x 1
   // veta = bigev x pin, ldveta = md
   // uz = nunstab x bigevs
   // dz = bigevs x 1
   // vz = bigevs x psin, ldvz = mds
   //---------------------------------

   //====== 03/01/06 TZ: the following note is added by CAS 3/10/04.
   //------ Code below allowed "existence" in cases where the initial lagged state was free to take on values
   //------ inconsistent with existence, so long as the state could w.p.1 remain consistent with a stable solution
   //------ if its initial lagged value was consistent with a stable solution.  This is a mistake, though perhaps there
   //------ are situations where we would like to know that this "existence for restricted initial state" situation holds.
   // psi = CreateComplexMatrix5RealMatrix(gensys_ps->Psi_dm);
   // zwt = tzMalloc(tmpi=nunstab*psin, MKL_Complex16);
   // InitializeConstantMLK_Complex16(zwt, tmpi, 0.0);    //Must be initialized to 0.0 in order to have legal values of this pointer.
   // cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, nunstab, psin, n,
   //             &one, q+n*(n-nunstab), n, psi, n, &zero, zwt, nunstab);
   // if (compute_svd(zwt, &uz, &dz, &vz, nunstab, psin)) {
   //    //Memory is now allocated to uz, dz, and vz.
   //    printf("WARNING: SV decomposition failed.\n");
   //    tzDestroy(ueta);
   //    tzDestroy(deta);
   //    tzDestroy(veta);
   //    tzDestroy(uz);
   //    tzDestroy(dz);
   //    tzDestroy(vz);
   //    tzDestroy(zwt);
   //    tzDestroy(a);
   //    tzDestroy(b);
   //    tzDestroy(q);
   //    tzDestroy(z);
   //    return;
   // }
   // tzDestroy(vz);
   // mds = nunstab<psin?nunstab:psin;
   // bigevs = mds;
   // for (i=0; i<mds; i++)
   //    if (dz[i]<=REALSMALL) {
   //       bigevs=i;
   //       break;
   //    }
   // tzDestroy(dz);
   //
   // if (!bigevs) {
   //    exist = 1;
   //    existx = 1;
   // } else {
   //    /* uz-ueta*ueta'*uz */
/*     //    MKL_Complex16 *tmp = tzMalloc(tmpi=nunstab*nunstab, MKL_Complex16);   ansi-c*/
/*     //    InitializeConstantMLK_Complex16(tmp, tmpi, 0.0);    //Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
/*     //    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, nunstab, nunstab,   ansi-c*/
/*     //       bigev, &one, ueta, nunstab, ueta, nunstab, &zero, tmp, nunstab);   ansi-c*/
/*     //    cblas_zhemm(CblasColMajor, CblasLeft, CblasUpper, nunstab,   ansi-c*/
/*     //       bigevs, &minusone, tmp, nunstab, uz, nunstab, &one, uz, nunstab);   ansi-c*/
/*     //    tzDestroy(tmp);   ansi-c*/
/*     //    if (compute_norm(uz, &norm, nunstab, bigevs)) {   ansi-c*/
/*     //       //Memory is now allocated to norm.   ansi-c*/
/*     //       printf("WARNING: SVD failed.\n");   ansi-c*/
/*     //       tzDestroy(norm);   ansi-c*/
/*     //       tzDestroy(ueta);   ansi-c*/
/*     //       tzDestroy(deta);   ansi-c*/
/*     //       tzDestroy(veta);   ansi-c*/
/*     //       tzDestroy(uz);   ansi-c*/
/*     //       tzDestroy(zwt);   ansi-c*/
/*     //       tzDestroy(a);   ansi-c*/
/*     //       tzDestroy(b);   ansi-c*/
/*     //       tzDestroy(q);   ansi-c*/
/*     //       tzDestroy(z);   ansi-c*/
/*     //       return;   ansi-c*/
/*     //    }   ansi-c*/
/*     //    exist = *norm < REALSMALL*n;   ansi-c*/
/*     //    tzDestroy(norm);   ansi-c*/
/*     //    if (compute_normx(a, b, zwt, ueta, &normx, nunstab, psin, n, bigev)) {   ansi-c*/
/*     //       //If 0, memory is now allocated to normx; otherwise, normx is destroyed within the function compute_normx().   ansi-c*/
/*     //       tzDestroy(ueta);   ansi-c*/
/*     //       tzDestroy(deta);   ansi-c*/
/*     //       tzDestroy(veta);   ansi-c*/
/*     //       tzDestroy(uz);   ansi-c*/
/*     //       tzDestroy(zwt);   ansi-c*/
/*     //       tzDestroy(a);   ansi-c*/
/*     //       tzDestroy(b);   ansi-c*/
/*     //       tzDestroy(q);   ansi-c*/
/*     //       tzDestroy(z);   ansi-c*/
/*     //       return;   ansi-c*/
/*     //    }   ansi-c*/
/*     //    existx = *normx < REALSMALL*n;   ansi-c*/
/*     //    tzDestroy(normx);   ansi-c*/
/*     // }   ansi-c*/
/*     //   ansi-c*/
/*     // tzDestroy(uz);   ansi-c*/
/*     // tzDestroy(zwt);   ansi-c*/

/*     //---------------------------------------------------------------------------   ansi-c*/
/*     // Note that existence and uniqueness are not just matters of comparing   ansi-c*/
/*     //   numbers of roots and numbers of endogenous errors.  These counts are   ansi-c*/
/*     //   reported below because usually they point to the source of the problem.   ansi-c*/
/*     //---------------------------------------------------------------------------   ansi-c*/
/*     //=============================================   ansi-c*/
/*     // Modified by DW and TZ to deal with the case   ansi-c*/
/*     //   where nunstab=n (all explosive roots).  03/08/06.   ansi-c*/
/*     //=============================================   ansi-c*/
   if (nunstab == n)
   {
      tzDestroy(pi);
      tzDestroy(ueta);
      tzDestroy(deta);
      tzDestroy(veta);
      tzDestroy(a);
      tzDestroy(b);
      tzDestroy(q);
      tzDestroy(z);

      printf("\n******** Fatal error: All roots are explosive while we have a solution.  But this should NOT happen.***********\n");
      eu[0] = 0;
      return 0;
   }

/*     //=======  Otherwise, returns to CAS's original code.  03/08/06. =======//   ansi-c*/
   etawt1 = tzMalloc(tmpi=(n-nunstab)*pin, MKL_Complex16);
   InitializeConstantMLK_Complex16(etawt1, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
   cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n-nunstab, pin, n,
      &one, q, n, pi, n, &zero, etawt1, n-nunstab);
   tzDestroy(pi);
   if (compute_svd(etawt1, &ueta1, &deta1, &veta1, n-nunstab, pin)) {
/*        //Memory is now allocated to ueta1, deta1, and veta1.   ansi-c*/
      printf("WARNING: SVD failed for compute_svd().\n");
      tzDestroy(ueta1);
      tzDestroy(deta1);
      tzDestroy(veta1);
      tzDestroy(etawt1);
      tzDestroy(ueta);
      tzDestroy(deta);
      tzDestroy(veta);
      tzDestroy(a);
      tzDestroy(b);
      tzDestroy(q);
      tzDestroy(z);
      eu[0] = 0;
      return 0;
   }
   tzDestroy(etawt1);
   md1 = n-nunstab<pin?n-nunstab:pin;
   bigev1 = md1;
   for (i=0; i<md1; i++)
      if (deta1[i]<=REALSMALL) {
         bigev1=i;
         break;
      }

/*     //====== 03/01/06 TZ: the following is commented out by CAS 3/10/04.   ansi-c*/
/*     // if (existx || !nunstab) {   ansi-c*/
/*     //    //=== Solution exists.   ansi-c*/
/*     //    eu[0] = 1;   ansi-c*/
/*     // } else {   ansi-c*/
/*     //    if (exist) {   ansi-c*/
/*     //       printf("WARNING: Solution exists for unforecastable z only\n");   ansi-c*/
/*     //       eu[0] = -1;   ansi-c*/
   //    } /* else
   //       mexPrintf("No solution.  %d unstable roots. %d endog errors.\n",nunstab,bigev1); */
   // /* mexPrintf("Generalized eigenvalues\n");
   //    mexCallMATLAB(0,NULL,1, &plhs[6], "disp"); */
/*     // }   ansi-c*/


/*     //-------------------------------   ansi-c*/
/*     // ueta1 = n-nunstab x bigev1   ansi-c*/
/*     // deta1 = bigev1 x 1   ansi-c*/
/*     // veta1 = bigev1 x pin, ldveta1 = md1   ansi-c*/
/*     //-------------------------------   ansi-c*/
   if (!bigev1)
      unique = 1;
   else {
/*        // veta1-veta1*veta*veta'   ansi-c*/
/*        // veta = bigev x pin, ldveta1 = md   ansi-c*/
/*        // veta1 = bigev1 x pin, ldveta1 = md1   ansi-c*/
      MKL_Complex16 *tmp = tzMalloc(pin*pin, MKL_Complex16);
      MKL_Complex16 *veta1_copy = tzMalloc(pin*bigev1, MKL_Complex16);
      InitializeConstantMLK_Complex16(tmp, pin*pin, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      InitializeConstantMLK_Complex16(veta1_copy, pin*bigev1, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      if (nunstab)
      {
         cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, pin, pin,
            bigev, &one, veta, md, veta, md, &zero, tmp, pin);    /*  tmp=veta'*veta;   ansi-c*/
         cblas_zdupe(bigev1,pin,veta1,md1,veta1_copy,bigev1);
         cblas_zhemm(CblasColMajor, CblasRight, CblasUpper, bigev1, pin,
            &minusone, tmp, pin, veta1_copy, bigev1, &one, veta1_copy, bigev1);
      }
      else //Added by DW & TZ, 03/08/06.
      {
         cblas_zdupe(bigev1,pin,veta1,md1,veta1_copy,bigev1);
      }
      tzDestroy(tmp);
      if (compute_norm(veta1_copy, &norm, bigev1, pin)) {
/*           //Memory is now allocated to norm.   ansi-c*/
         printf("WARNING: SVD failed.\n");
         tzDestroy(norm);
         tzDestroy(ueta1);
         tzDestroy(deta1);
         tzDestroy(veta1);
         tzDestroy(ueta);
         tzDestroy(deta);
         tzDestroy(veta);
         tzDestroy(veta1_copy);
         tzDestroy(a);
         tzDestroy(b);
         tzDestroy(q);
         tzDestroy(z);
         eu[0] = 0;
         return 0;
      }
      tzDestroy(veta1_copy);
      unique = *norm < REALSMALL*n;
      tzDestroy(norm);
   }
   if (unique) {
/*        //=== Unique solution.   ansi-c*/
      eu[1] = 1;
   } else {
      eu[1] = 0;
      #if defined (PRINTWARNINGofSUNSPOT)
      if (nunstab)
         printf("WARNING: Indeterminacy.  %d loose endog errors with eu being [%d, %d].\n",bigev1-bigev, eu[0], eu[1]);
      else
         printf("WARNING: Indeterminacy.  %d loose endog errors with eu being [%d, %d].\n",pin, eu[0], eu[1]);
/*        //printf("WARNING: Indeterminacy.  %d loose endog errors with eu being [%g, %g].\n",bigev1-bigev, gensys_ps->eu_dv->v[0], gensys_ps->eu_dv->v[1]);   ansi-c*/
/*        //printf("WARNING: Indeterminacy.  %d loose endog errors.\n",bigev1-bigev);   ansi-c*/
      #endif
   }

/*     //---------------------------------------------------------//   ansi-c*/
/*     //------------------ Obtaining the outputs. ---------------//   ansi-c*/
/*     //---------------------------------------------------------//   ansi-c*/
   if (nunstab)
   {
/*        //=== All the following lines are used to compute only ONE object tmat, which is used subsequently. ===//   ansi-c*/
      cblas_zdscali(pin,deta,bigev,veta,md);      /* veta' = deta\veta' */
      tzDestroy(deta);
      cblas_zdscale(pin,deta1,bigev1,veta1,md1);  /* veta1' = deta1*veta1' */
      tzDestroy(deta1);
      etawt = tzMalloc(tmpi=nunstab*pin, MKL_Complex16);      /* etawt = ueta*veta' */
      InitializeConstantMLK_Complex16(etawt, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nunstab, pin, bigev,
         &one, ueta, nunstab, veta, md, &zero, etawt, nunstab);
      tzDestroy(ueta);
      tzDestroy(veta);
      etawt1 = tzMalloc(tmpi=(n-nunstab)*pin, MKL_Complex16); /* etawt1 = ueta1*veta1' */
      InitializeConstantMLK_Complex16(etawt1, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n-nunstab, pin, bigev1,
         &one, ueta1, n-nunstab, veta1, md1, &zero, etawt1, n-nunstab);
      tzDestroy(ueta1);
      tzDestroy(veta1);
      tmat = tzMalloc(tmpi=(n-nunstab)*nunstab, MKL_Complex16); /* tmat = etawt1*etawt' */
      InitializeConstantMLK_Complex16(tmat, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, n-nunstab, nunstab, pin,
         &one, etawt1, n-nunstab, etawt, nunstab, &zero, tmat, n-nunstab);
      tzDestroy(etawt1);
      tzDestroy(etawt);

/*        //=== Getting the solution Theta ===//   ansi-c*/
      g0 = tzMalloc(tmpi=n*n, MKL_Complex16);
      InitializeConstantMLK_Complex16(g0, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zdupe(n-nunstab, n, a, n, g0, n);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n-nunstab, nunstab, nunstab,
         &minusone, tmat, n-nunstab, a+(n-nunstab)*(n+1), n, &one, g0+(n-nunstab)*n, n);
      cblas_zcopy(nunstab, &one, 0, g0+(n-nunstab)*(n+1), n+1);

      g1 = tzMalloc(tmpi=n*n, MKL_Complex16);
      InitializeConstantMLK_Complex16(g1, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zdupe(n-nunstab, n, b, n, g1, n);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n-nunstab, nunstab, nunstab,
         &minusone, tmat, n-nunstab, b+(n-nunstab)*(n+1), n, &one, g1+(n-nunstab)*n, n);
      cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
         n, n, &one, g0, n, g1, n);
      dummy = tzMalloc(tmpi=n*n, MKL_Complex16);
      InitializeConstantMLK_Complex16(dummy, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, z, n, g1, n, &zero, dummy, n);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, n, n, n, &one, dummy, n, z, n, &zero, g1, n);
      tzDestroy(dummy);
      ComplexMatrix2RealMatrix(gensys_ps->Theta_dm, g1);   /*  Output.   ansi-c*/
      tzDestroy(g1);

/*        //=== Getting the constant term c ===//   ansi-c*/
      tmatq = tzMalloc(tmpi=n*n, MKL_Complex16);
      InitializeConstantMLK_Complex16(tmatq, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zcopy(n*n, q, 1, tmatq, 1);
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, n, n-nunstab, nunstab,
         &minusone, tmatq+(n-nunstab)*n, n, tmat, n-nunstab, &one, tmatq, n);
      tzDestroy(tmat);

      ab = tzMalloc(tmpi=nunstab*nunstab, MKL_Complex16);
      InitializeConstantMLK_Complex16(ab, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zdpsb(nunstab, nunstab, a+(n-nunstab)*(n+1), n, b+(n-nunstab)*(n+1), n, ab, nunstab);
      cblas_ztrsm(CblasColMajor, CblasRight, CblasUpper, CblasConjTrans, CblasNonUnit,
         n, nunstab, &one, ab, nunstab, tmatq+(n-nunstab)*n, n);
      tzDestroy(ab);

      c = CreateComplexMatrix5RealVector(gensys_ps->c0_dv);
      dummy = tzMalloc(gensys_ps->c0_dv->n, MKL_Complex16);
/*        //$$$$$$$ The following is Iskander's fatal code error.  One cannot use c in the two different places; otherwise, it makes c be zero completely!   ansi-c*/
/*        // cblas_zgemv(CblasColMajor, CblasConjTrans, n, n, &one, tmatq, n, c, 1, &zero, c, 1);   ansi-c*/
/*        // cblas_ztrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, g0, n, c, 1);   ansi-c*/
/*        // cblas_zgemv(CblasColMajor, CblasNoTrans, n, n, &one, z, n, c, 1, &zero, c, 1);   ansi-c*/
      cblas_zgemv(CblasColMajor, CblasConjTrans, n, n, &one, tmatq, n, c, 1, &zero, dummy, 1);
      cblas_ztrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, g0, n, dummy, 1);
      cblas_zgemv(CblasColMajor, CblasNoTrans, n, n, &one, z, n, dummy, 1, &zero, c, 1);
      ComplexMatrix2RealVector(gensys_ps->c_dv, c);    /*  Output.   ansi-c*/
      tzDestroy(c);
      tzDestroy(dummy);

/*        //=== Getting the term Impact ===//   ansi-c*/
      impact = tzMalloc(tmpi=n*psin, MKL_Complex16);
      InitializeConstantMLK_Complex16(impact, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      psi = CreateComplexMatrix5RealMatrix(gensys_ps->Psi_dm);  /*  03/01/06 TZ.  Added to be consistent with the CAS 3/10/04 correction.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n-nunstab, psin, n,
         &one, tmatq, n, psi, n, &zero, impact, n);
      tzDestroy(tmatq);
      cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
         n, psin, &one, g0, n, impact, n);
      dummy = tzMalloc(tmpi=n*psin, MKL_Complex16);
      InitializeConstantMLK_Complex16(dummy, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, psin, n,
         &one, z, n, impact, n, &zero, dummy, n);
      tzDestroy(impact);
      ComplexMatrix2RealMatrix(gensys_ps->Impact_dm, dummy);    /*  Output.   ansi-c*/
      tzDestroy(dummy);

/*        //=== Finishing up the other terms such as Fmat, Fwt, and Ywt. ===//   ansi-c*/
      fmat = a+(n-nunstab)*(n+1);
      cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
         nunstab, nunstab, &one, b+(n-nunstab)*(n+1), n, fmat, n);
      gensys_ps->Fmat_dzm = SubComplexMatrix2Zmatrix(gensys_ps->Fmat_dzm, fmat, n, nunstab, nunstab);
      tzDestroy(a);

      fwt = tzMalloc(tmpi=nunstab*psin, MKL_Complex16);
      InitializeConstantMLK_Complex16(fwt, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_ztrsm(CblasColMajor, CblasRight, CblasUpper, CblasConjTrans, CblasNonUnit,
         n, nunstab, &one, b+(n-nunstab)*(n+1), n, q+(n-nunstab)*n, n);
      tzDestroy(b);
      cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, nunstab, psin, n,
         &minusone, q+(n-nunstab)*n, n, psi, n, &zero, fwt, nunstab);
      tzDestroy(q);
      tzDestroy(psi);
      gensys_ps->Fwt_dzm = SubComplexMatrix2Zmatrix(gensys_ps->Fwt_dzm, fwt, nunstab, nunstab, psin);
      tzDestroy(fwt);

      ywt = tzMalloc(tmpi=n*nunstab, MKL_Complex16);
      InitializeConstantMLK_Complex16(ywt,  tmpi, 0.0);
      cblas_zcopy(nunstab, &one, 0, ywt+n-nunstab, n+1);
      cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
         n, nunstab, &one, g0, n, ywt, n);
      tzDestroy(g0);
      dummy = tzMalloc(tmpi=n*nunstab, MKL_Complex16);
      InitializeConstantMLK_Complex16(dummy, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, nunstab, n,
         &one, z, n, ywt, n, &zero, dummy, n);
      tzDestroy(z);
      tzDestroy(ywt);
      gensys_ps->Ywt_dzm = SubComplexMatrix2Zmatrix(gensys_ps->Ywt_dzm, dummy, n, n, nunstab);
      tzDestroy(dummy);
   }
   else  //This part is added by DW and TZ, 03/08/06.
   {
/*        //======= Getting Theta = real(z*(G0\G1)*z') =======//   ansi-c*/
      cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
         n, n, &one, a, n, b, n);   /*  Note that a is triangular and b = a\b (overwritten).   ansi-c*/
      dummy = tzMalloc(tmpi=n*n, MKL_Complex16);
      InitializeConstantMLK_Complex16(dummy, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
/*        //--- Getting Theta = real(z*b*z');   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, z, n, b, n, &zero, dummy, n);  /*  dummy=z*b;   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, n, n, n, &one, dummy, n, z, n, &zero, b, n);  /*  dummy=dummy*z';   ansi-c*/
      ComplexMatrix2RealMatrix(gensys_ps->Theta_dm, b);   /*  Output.   ansi-c*/
      tzDestroy(dummy);

/*        //======= Getting c = real(z*G0\q*c0) =======//   ansi-c*/
      c = CreateComplexMatrix5RealVector(gensys_ps->c0_dv);
      dummy = tzMalloc(gensys_ps->c0_dv->n, MKL_Complex16);
      cblas_zgemv(CblasColMajor, CblasConjTrans, n, n, &one, q, n, c, 1, &zero, dummy, 1);  /*  dummy = q*c;   ansi-c*/
      cblas_ztrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, a, n, dummy, 1);   /*  dummy=a\dummy where a is triangular.   ansi-c*/
      cblas_zgemv(CblasColMajor, CblasNoTrans, n, n, &one, z, n, dummy, 1, &zero, c, 1);   /*  dummy=z*dummy;   ansi-c*/
      ComplexMatrix2RealVector(gensys_ps->c_dv, c);    /*  Output.   ansi-c*/
      tzDestroy(dummy);

/*        //======= Getting Impact = real(z*G0\q*psi) =======//   ansi-c*/
      impact = tzMalloc(tmpi=n*psin, MKL_Complex16);
      InitializeConstantMLK_Complex16(impact, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      psi = CreateComplexMatrix5RealMatrix(gensys_ps->Psi_dm);  /*  03/01/06 TZ.  Added to be consistent with the CAS 3/10/04 correction.   ansi-c*/
      dummy = tzMalloc(tmpi=n*psin, MKL_Complex16);
      cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, n, psin, n,
         &one, q, n, psi, n, &zero, impact, n);    /*  impact = q*psi;   ansi-c*/
      cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
         n, psin, &one, a, n, impact, n);   /*  impact = a\impact;   ansi-c*/
      InitializeConstantMLK_Complex16(dummy, tmpi, 0.0);     /*  Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
      cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, psin, n,
         &one, z, n, impact, n, &zero, dummy, n);    /*  dummy = z*impact;   ansi-c*/
      ComplexMatrix2RealMatrix(gensys_ps->Impact_dm, dummy);    /*  Output.   ansi-c*/
      tzDestroy(dummy);


/*        //=== Some of destructions may have been done, but it is better to be safe.   ansi-c*/
      tzDestroy(ueta1);
      tzDestroy(deta1);
      tzDestroy(veta1);
      tzDestroy(etawt);
      tzDestroy(etawt1);
/*        //+   ansi-c*/
      tzDestroy(a);
      tzDestroy(b);
      tzDestroy(q);
      tzDestroy(z);
/*        //+   ansi-c*/
      tzDestroy(c);
      tzDestroy(impact);
      tzDestroy(psi);
   }

/*     //=== Save this debugging format -- DDDDDebugging.   ansi-c*/
/*     // if (!nunstab)   ansi-c*/
/*     // {   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "Aind=[\n");   ansi-c*/
/*     //    WriteMatrix(FPTR_DEBUG, gensys_ps->G0_dm, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "Bind=[\n");   ansi-c*/
/*     //    WriteMatrix(FPTR_DEBUG, gensys_ps->G1_dm, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "Consterm=[\n");   ansi-c*/
/*     //    WriteVector(FPTR_DEBUG, gensys_ps->c0_dv, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "]';\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "gUpsiloneind=[\n");   ansi-c*/
/*     //    WriteMatrix(FPTR_DEBUG, gensys_ps->Psi_dm, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "gUpsilonxind=[\n");   ansi-c*/
/*     //    WriteMatrix(FPTR_DEBUG, gensys_ps->Pi_dm, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*     //    fflush(FPTR_DEBUG);   ansi-c*/
/*     //   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "\n********** Output ******************\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "Theta=[\n");   ansi-c*/
/*     //    WriteMatrix(FPTR_DEBUG, gensys_ps->Theta_dm, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "Impact=[\n");   ansi-c*/
/*     //    WriteMatrix(FPTR_DEBUG, gensys_ps->Impact_dm, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "Consterm=[\n");   ansi-c*/
/*     //    WriteVector(FPTR_DEBUG, gensys_ps->c_dv, " %.16e ");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "]';\n");   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "eu=[\n");   ansi-c*/
/*     //    WriteVector_int(FPTR_DEBUG, gensys_ps->eu_iv);   ansi-c*/
/*     //    fprintf(FPTR_DEBUG, "]';\n");   ansi-c*/
/*     //    fflush(FPTR_DEBUG);   ansi-c*/
/*     // }   ansi-c*/

   return 1;
}


/******************************************************************************
 * function selctg orders the eigenvalues so that a selected cluster of       *
 * eigenvalues appears in the leading diagonal blocks of the upper            *
 * quasi-triangular matrix S and the upper triangular matrix T.               *
 ******************************************************************************/

static int selctg(MKL_Complex16 *alpha, MKL_Complex16 *beta)
{
   double absA = sqrt(alpha->real*alpha->real+alpha->imag*alpha->imag),
          absB = fabs(beta->real);
   if (absA) {
      double divhat = absB/absA;
/*        //bug detected by Vasco Curdia and Daria Finocchiaro, 2/25/2004 CAS. A root of   ansi-c*/
/*        //exactly 1.01 and no root between 1 and 1.02, led to div being stuck at 1.01   ansi-c*/
/*        //and the 1.01 root being misclassified as stable.  Changing < to <= below fixes this.   ansi-c*/
      if (fixdiv && 1+REALSMALL<divhat && divhat<=stake)
         stake = (1+divhat)/2;
   }
   if (absA<REALSMALL && absB<REALSMALL)
      zxz = 1;
   if (absB>stake*absA) {
      nunstab++;
      return(0);
   } else
      return(1);
}

/******************************************************************************
 * compute for a pair of N-by-N complex nonsymmetric matrices (A,B),          *
 * the generalized eigenvalues, the generalized complex Schur form (S, T),    *
 * and optionally left and/or right Schur vectors (VSL and VSR)               *
 ******************************************************************************/

static int qz(MKL_Complex16 *a, MKL_Complex16 *b, MKL_Complex16 *q, MKL_Complex16 *z, int n)
{
/*  //   unsigned char msg[101];   ansi-c*/
   int sdim, lwork = -1, info = 0;
   MKL_Complex16 *alpha = tzMalloc(n,MKL_Complex16),
                 *beta = tzMalloc(n,MKL_Complex16),
                 *work, work1;
   double *rwork = tzMalloc(8*n, double);
   int *bwork = tzMalloc(4*n, int);

   /* Query zgges on the value of lwork */
   zgges("V", "V", "S", &selctg, &n, a, &n, b, &n, &sdim, alpha, beta, q,
         &n, z, &n, &work1, &lwork, rwork, bwork, &info);

   if (info < 0) {
      printf("WARNING: Input %d to the Intel MKL function zgges() has an illegal value",-info);
      tzDestroy(bwork);
      tzDestroy(rwork);
      tzDestroy(alpha);
      tzDestroy(beta);
      return(info);
   }

   lwork = (int)(work1.real);
   work = tzMalloc(lwork, MKL_Complex16);
   zgges("V", "V", "S", &selctg, &n, a, &n, b, &n, &sdim, alpha, beta, q,
         &n, z, &n, work, &lwork, rwork, bwork, &info);

   tzDestroy(work);
   tzDestroy(bwork);
   tzDestroy(rwork);
   tzDestroy(alpha);
   tzDestroy(beta);

   if (info < 0) {
      printf("WARNING: Input %d to the Intel MKL function zgges() has an illegal value",-info);
      return(info);
   }

   if (info > 0)
      if (info < n)
         printf("WARNING: The QZ iteration failed.  (A,B) are not in Schur form,\n"
             "but ALPHA(j) and BETA(j) should be correct for j=%d,...,N.",info+1);
      else {
         switch (info-n) {
         case 1:
            printf("WARNING: LAPACK problem: error return from ZGGBAL");
            break;
         case 2:
            printf("WARNING: LAPACK problem: error return from ZGEQRF");
            break;
         case 3:
            printf("WARNING: LAPACK problem: error return from ZUNMQR");
            break;
         case 4:
            printf("WARNING: LAPACK problem: error return from ZUNGQR");
            break;
         case 5:
            printf("WARNING: LAPACK problem: error return from ZGGHRD");
            break;
         case 6:
            printf("WARNING: LAPACK problem: error return from ZHGEQZ (other than failed iteration)");
            break;
         case 7:
            printf("WARNING: LAPACK problem: error return from ZGGBAK (computing VSL)");
            break;
         case 8:
            printf("WARNING: LAPACK problem: error return from ZGGBAK (computing VSR)");
            break;
         case 9:
            printf("WARNING: LAPACK problem: error return from ZLASCL (various places)");
            break;
         default:
            printf("WARNING: LAPACK problem: unknown error.");
            break;
         }
      }
   return(info);
}

/*
 * Convert MATLAB complex matrix to MKL complex storage.
 *
 *          Z = mat2mkl(X,ldz,ndz)
 *
 * converts MATLAB's mxArray X to MKL_Complex16 Z(ldz,ndz).
 * The parameters ldz and ndz determine the storage allocated for Z,
 * while mxGetM(X) and mxGetN(X) determine the amount of data copied.
 */

/*  //MKL_Complex16* mat2mkl(const mxArray *X, int ldz, int ndz) {   ansi-c*/
/*  //   MKL_Complex16 *Z, *zp;   ansi-c*/
/*  //   int m, n, incz, cmplxflag;   ansi-c*/
/*  //   register int i, j;   ansi-c*/
/*  //   double *xr, *xi;   ansi-c*/

/*  //   Z = mxCalloc(ldz*ndz, sizeof(MKL_Complex16));   ansi-c*/
/*  //   xr = mxGetPr(X);   ansi-c*/
/*  //   xi = mxGetPi(X);   ansi-c*/
/*  //   m =  mxGetM(X);   ansi-c*/
/*  //   n =  mxGetN(X);   ansi-c*/
/*  //   zp = Z;   ansi-c*/
/*  //   incz = ldz-m;   ansi-c*/
/*  //   cmplxflag = (xi != NULL);   ansi-c*/
/*  //   for (j = 0; j < n; j++) {   ansi-c*/
/*  //      if (cmplxflag) {   ansi-c*/
/*  //         for (i = 0; i < m; i++) {   ansi-c*/
/*  //            zp->real = *xr++;   ansi-c*/
/*  //            zp->imag = *xi++;   ansi-c*/
/*  //            zp++;   ansi-c*/
/*  //         }   ansi-c*/
/*  //      } else {   ansi-c*/
/*  //         for (i = 0; i < m; i++) {   ansi-c*/
/*  //            zp->real = *xr++;   ansi-c*/
/*  //            zp++;   ansi-c*/
/*  //         }   ansi-c*/
/*  //      }   ansi-c*/
/*  //      zp += incz;   ansi-c*/
/*  //   }   ansi-c*/
/*  //   return(Z);   ansi-c*/
/*  //}   ansi-c*/


/*
 * Convert MKL complex storage to MATLAB real and imaginary parts.
 *
 *          X = mkl2mat(Z,ldz,m,n)
 *
 * copies MKL_Complex16 Z to X, producing a complex mxArray
 * with mxGetM(X) = m and mxGetN(X) = n.
 */

/*  //mxArray* mkl2mat(MKL_Complex16 *Z, int ldz, int m, int n) {   ansi-c*/
/*  //   int i, j, incz;   ansi-c*/
/*  //   double *xr, *xi;   ansi-c*/
/*  //   MKL_Complex16 *zp;   ansi-c*/
/*  //   mxArray *X;   ansi-c*/

/*  //   X = mxCreateDoubleMatrix(m,n,mxCOMPLEX);   ansi-c*/
/*  //   xr = mxGetPr(X);   ansi-c*/
/*  //   xi = mxGetPi(X);   ansi-c*/
/*  //   zp = Z;   ansi-c*/
/*  //   incz = ldz-m;   ansi-c*/
/*  //   for (j = 0; j < n; j++) {   ansi-c*/
/*  //      for (i = 0; i < m; i++) {   ansi-c*/
/*  //         *xr++ = zp->real;   ansi-c*/
/*  //         *xi++ = zp->imag;   ansi-c*/
/*  //         zp++;   ansi-c*/
/*  //      }   ansi-c*/
/*  //      zp += incz;   ansi-c*/
/*  //   }   ansi-c*/
/*  //   return(X);   ansi-c*/
/*  //}   ansi-c*/

/*  //plhs[3] = mkl2mat(fmat, n, nunstab, nunstab)   ansi-c*/

/*
 * Convert MKL complex storage to MATLAB real matrix ignoring imaginary part.
 *
 *          X = mkl2mat(Z,ldz,m,n)
 *
 * copies MKL_Complex16 Z to X, producing a real mxArray
 * with mxGetM(X) = m and mxGetN(X) = n.
 */

/*  //mxArray* mkl2mat_real(MKL_Complex16 *Z, int ldz, int m, int n) {   ansi-c*/
/*  //   int i, j, incz;   ansi-c*/
/*  //   double *xr;   ansi-c*/
/*  //   MKL_Complex16 *zp;   ansi-c*/
/*  //   mxArray *X;   ansi-c*/

/*  //   X = mxCreateDoubleMatrix(m,n,mxREAL);   ansi-c*/
/*  //   xr = mxGetPr(X);   ansi-c*/
/*  //   zp = Z;   ansi-c*/
/*  //   incz = ldz-m;   ansi-c*/
/*  //   for (j = 0; j < n; j++) {   ansi-c*/
/*  //      for (i = 0; i < m; i++) {   ansi-c*/
/*  //         *xr++ = zp->real;   ansi-c*/
/*  //         zp++;   ansi-c*/
/*  //      }   ansi-c*/
/*  //      zp += incz;   ansi-c*/
/*  //   }   ansi-c*/
/*  //   return(X);   ansi-c*/
/*  //}   ansi-c*/

/*  //void copy_eigenvalues(mxArray *gev, MKL_Complex16 *a, MKL_Complex16 *b, int n) {   ansi-c*/
/*  //   double *gevr = mxGetPr(gev),   ansi-c*/
/*  //          *gevi = mxGetPi(gev);   ansi-c*/
/*  //   int i;   ansi-c*/

/*  //   for (i=0; i<n; i++, gevr++, gevi++, a+=n+1) {   ansi-c*/
/*  //      *gevr = a->real;   ansi-c*/
/*  //      *gevi = a->imag;   ansi-c*/
/*  //   }   ansi-c*/

/*  //   for (i=0; i<n; i++, gevr++, gevi++, b+=n+1) {   ansi-c*/
/*  //      *gevr = b->real;   ansi-c*/
/*  //      *gevi = b->imag;   ansi-c*/
/*  //   }   ansi-c*/
/*  //}   ansi-c*/

static void copy_eigenvalues(TSdzmatrix *Gev_dzm, MKL_Complex16 *a, MKL_Complex16 *b)
{
   int n = Gev_dzm->real->nrows;
   double *gevr = Gev_dzm->real->M,
          *gevi = Gev_dzm->imag->M;
   int i;

   for (i=0; i<n; i++, gevr++, gevi++, a+=n+1) {
      *gevr = a->real;
      *gevi = a->imag;
   }

   for (i=0; i<n; i++, gevr++, gevi++, b+=n+1) {
      *gevr = b->real;
      *gevi = b->imag;
   }
}


static int compute_svd(MKL_Complex16 *x, MKL_Complex16 **u, double **d, MKL_Complex16 **v, int m, int n)
{
/*     //$$$Memory allocated to u, d, and v will be destroyed outside this function.$$$   ansi-c*/
   int tmpi;
   int md = m<n?m:n, lwork = -1, info = 0;
   MKL_Complex16 *a, *work, work1;
   double *rwork = tzMalloc(5*md>1?5*md:1, double);

   a = tzMalloc(m*n, MKL_Complex16);
   cblas_zcopy(m*n, x, 1, a, 1);

   *u = tzMalloc(tmpi=m*md,MKL_Complex16);
   InitializeConstantMLK_Complex16(*u, tmpi, 0.0);
   *v = tzMalloc(tmpi=md*n, MKL_Complex16);
   InitializeConstantMLK_Complex16(*v, tmpi, 0.0);
   *d = tzMalloc(md, double);
   InitializeConstantDouble(*d, md, 0.0);

   /* Query zgges on the value of lwork */
   zgesvd("S", "S", &m, &n, a, &m, *d, *u, &m, *v, &md, &work1, &lwork, rwork, &info);

   if (info < 0) {
      printf("WARNING: Input %d to zgesvd had an illegal value",-info);
      tzDestroy(rwork);
      return(info);
   }

   lwork = (int)(work1.real);
   work = tzMalloc(lwork, MKL_Complex16);
   zgesvd("S", "S", &m, &n, a, &m, *d, *u, &m, *v, &md, work, &lwork, rwork, &info);

   tzDestroy(work);
   tzDestroy(rwork);
   tzDestroy(a);

   if (info < 0)
      printf("WARNING: Input %d to zgesvd had an illegal value",-info);

   if (info > 0)
      printf("WARNING: ZBDSQR did not converge.\n%d superdiagonals of an intermediate "
         "bidiagonal form B did not converge to zero.",info);
   return(info);
}

static int compute_norm(MKL_Complex16 *a, double **d, int m, int n)
{
/*     //Memory will be allocated to d, which will be destroyed outside this function.   ansi-c*/
   int md = m<n?m:n, lwork = -1, info = 0;
   MKL_Complex16 *work = NULL, work1;
   double *rwork = tzMalloc(5*md>1?5*md:1, double);

   *d = tzMalloc(md, double);

   /* Query zgges on the value of lwork */
   zgesvd("N", "N", &m, &n, a, &m, *d, NULL, &m, NULL, &md, &work1, &lwork, rwork, &info);

   if (info < 0) {
      printf("WARNING: Input %d to zgesvd had an illegal value",-info);
      tzDestroy(rwork);
      return(info);
   }

   lwork = (int)(work1.real);
   work = tzMalloc(lwork, MKL_Complex16);
   zgesvd("N", "N", &m, &n, a, &m, *d, NULL, &m, NULL, &md, work, &lwork, rwork, &info);

   tzDestroy(work);
   tzDestroy(rwork);

   if (info < 0)
      printf("WARNING: Input %d to zgesvd had an illegal value",-info);

   if (info > 0)
      printf("WARNING: ZBDSQR() in Intel MKL did not converge.\n%d superdiagonals of an intermediate "
         "bidiagonal form B did not converge to zero.",info);

   return(info);
}


/*  //======= 03/01/06 TZ.  Commented out to be consistent with the CAS 3/10/04 correction. =======//   ansi-c*/
/*  //static int compute_normx(MKL_Complex16 *a, MKL_Complex16 *b, MKL_Complex16 *zwt, MKL_Complex16 *ueta, double **normx, int nunstab, int psin, int n, int bigev)   ansi-c*/
/*  //{   ansi-c*/
/*  //   //Memory is allocated to normx, which will be freed outside this function.   ansi-c*/
/*  //   int tmpi;   ansi-c*/
/*  //   int info = 0, i, bigevs;   ansi-c*/
/*  //   //   ansi-c*/
/*  //   MKL_Complex16 *M = NULL, *zwtx = NULL, *ux = NULL, *vx = NULL, *tmp = NULL;   ansi-c*/
/*  //   double *dx = NULL;   ansi-c*/


/*  //   a += (n+1)*(n-nunstab);   ansi-c*/
/*  //   b += (n+1)*(n-nunstab);   ansi-c*/
/*  //   cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,   ansi-c*/
/*  //      nunstab, psin, &one, b, n, zwt, nunstab);   ansi-c*/
/*  //   M = tzMalloc(nunstab*nunstab, MKL_Complex16);   ansi-c*/
/*  //   cblas_zdupe(nunstab, nunstab, a, n, M, nunstab);   ansi-c*/
/*  //   cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,   ansi-c*/
/*  //      nunstab, nunstab, &one, b, n, M, nunstab);   ansi-c*/

/*  //   zwtx = tzMalloc(nunstab*nunstab*psin, MKL_Complex16);   ansi-c*/
/*  //   cblas_zcopy(nunstab*psin, zwt, 1, zwtx, 1);   ansi-c*/
/*  //   for (i=1; i<nunstab; i++) {   ansi-c*/
/*  //      cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, nunstab, psin*i, &one, M, nunstab, zwtx, nunstab);   ansi-c*/
/*  //      cblas_zcopy(nunstab*psin, zwt, 1, zwtx+nunstab*psin*i, 1);   ansi-c*/
/*  //   }   ansi-c*/
/*  //   tzDestroy(M);   ansi-c*/
/*  //   cblas_ztrmm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, nunstab, nunstab*psin, &one, b, n, zwtx, nunstab);   ansi-c*/
/*  //   info = compute_svd(zwtx, &ux, &dx, &vx, nunstab, nunstab*psin);  //Memory is allocated to ux, dx, and vx.   ansi-c*/
/*  //   tzDestroy(vx);   ansi-c*/
/*  //   tzDestroy(zwtx);   ansi-c*/
/*  //   if (info) {   ansi-c*/
/*  //      printf("WARNING: SVD failed.\n");   ansi-c*/
/*  //      tzDestroy(ux);   ansi-c*/
/*  //      tzDestroy(dx);   ansi-c*/
/*  //      return(info);   ansi-c*/
/*  //   }   ansi-c*/
/*  //   bigevs = nunstab;   ansi-c*/
/*  //   for (i=0; i<nunstab; i++)   ansi-c*/
/*  //      if (dx[i]<=REALSMALL) {   ansi-c*/
/*  //         bigevs = i;   ansi-c*/
/*  //         break;   ansi-c*/
/*  //      }   ansi-c*/
/*  //   tzDestroy(dx);   ansi-c*/
//   /* ux-ueta*ueta'*ux */
/*  //   tmp = tzMalloc(tmpi=nunstab*nunstab, MKL_Complex16);   ansi-c*/
/*  //   InitializeConstantMLK_Complex16(tmp, tmpi, 0.0);    //Must be initialized to 0.0 in order to have legal values of this pointer.   ansi-c*/
/*  //   cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, nunstab, nunstab,   ansi-c*/
/*  //      bigev, &one, ueta, nunstab, ueta, nunstab, &zero, tmp, nunstab);   ansi-c*/
/*  //   cblas_zhemm(CblasColMajor, CblasLeft, CblasUpper, nunstab,   ansi-c*/
/*  //      bigevs, &minusone, tmp, nunstab, ux, nunstab, &one, ux, nunstab);   ansi-c*/
/*  //   tzDestroy(tmp);   ansi-c*/
/*  //   info = compute_norm(ux, normx, nunstab, bigevs);  //Memory is allocated to normx.   ansi-c*/
/*  //   if (info) {   ansi-c*/
/*  //      printf("WARNING: SVD failed.\n");   ansi-c*/
/*  //      tzDestroy(normx);   ansi-c*/
/*  //      tzDestroy(ux);   ansi-c*/
/*  //      return(info);   ansi-c*/
/*  //   }   ansi-c*/
/*  //   tzDestroy(ux);   ansi-c*/
/*  //   return(info);   ansi-c*/
/*  //}   ansi-c*/

static void cblas_zdupe(int m, int n, MKL_Complex16 *a, int lda, MKL_Complex16 *b, int ldb)
{
/*     //Copying from a to b.   ansi-c*/
   int i;
   for (i=0; i<m; i++, a++, b++)
      cblas_zcopy(n, a, lda, b, ldb);
}

static void cblas_zdscali(int n, double *a, int lda, MKL_Complex16 *b, int ldb)
{
   int i;
   for (i=0; i<lda; i++, a++, b++)
      cblas_zdscal(n, 1.0/(*a), b, ldb);
}

static void cblas_zdscale(int n, double *a, int lda, MKL_Complex16 *b, int ldb)
{
   int i;
   for (i=0; i<lda; i++, a++, b++)
      cblas_zdscal(n, *a, b, ldb);
}

static void cblas_zdpsb(int m, int n, MKL_Complex16 *a, int lda, MKL_Complex16 *b, int ldb, MKL_Complex16 *c, int ldc)
{
   int i;
   cblas_zdupe(m, n, a, lda, c, ldc);
   for (i=0; i<m; i++, b++, c++)
      cblas_zaxpy(n, &minusone, b, ldb, c, ldc);
}


static MKL_Complex16* CreateComplexMatrix5RealMatrix(TSdmatrix *X_dm)
{
   int mn, k;
   double *M;
/*     //   ansi-c*/
   MKL_Complex16 *Z = NULL;

   if (!X_dm)  fn_DisplayError("CreateComplexMatrix5RealMatrix():  Input matrix X_dm must be allocated memory");
   M = X_dm->M;

   Z = tzMalloc(mn=X_dm->nrows*X_dm->ncols, MKL_Complex16);
   for (k=mn-1; k>=0; k--) {
      Z[k].real = M[k];
      Z[k].imag = 0.0;
   }
   return(Z);
}

static MKL_Complex16* CreateComplexMatrix5RealVector(TSdvector *x_dv)
{
   int n, k;
   double *v;
/*     //   ansi-c*/
   MKL_Complex16 *Z = NULL;

   if (!x_dv)  fn_DisplayError("CreateComplexMatrix5RealVector():  Input vector x_dv must be allocated memory");
   v = x_dv->v;

   Z = tzMalloc(n=x_dv->n, MKL_Complex16);
   for (k=n-1; k>=0; k--) {
      Z[k].real = v[k];
      Z[k].imag = 0.0;
   }
   return(Z);
}


static void ComplexMatrix2RealMatrix(TSdmatrix *Y_dm, MKL_Complex16 *Z)
{
   int _k;
   double *M;

   if (!Y_dm)  fn_DisplayError("ComplexMatrix2RealMatrix():  Output matrix Y_dm must be allocated memory");
   M = Y_dm->M;

   for (_k=Y_dm->nrows*Y_dm->ncols-1; _k>=0; _k--)  M[_k] = Z[_k].real;
   Y_dm->flag = M_GE;
}


static void ComplexMatrix2RealVector(TSdvector *y_dv, MKL_Complex16 *Z)
{
   int _k;
   double *v;

   if (!y_dv)  fn_DisplayError("ComplexMatrix2RealVector():  Output matrix y_dv must be allocated memory");
   v = y_dv->v;

   for (_k=y_dv->n-1; _k>=0; _k--)  v[_k] = Z[_k].real;
   y_dv->flag = V_DEF;
}


static TSdzmatrix *SubComplexMatrix2Zmatrix(TSdzmatrix *X_dzm, MKL_Complex16 *Z, const int nrowsforZ, const int _m, const int _n)
{
/*     //X_dzm is _m-by_n comlex types where nrowsforZ <= _m and Z is nrowsforZ-by-_n.   ansi-c*/
   int _i, _j, incz;
   double *Mreal, *Mimag;
   MKL_Complex16 *zp;
/*     //   ansi-c*/
   TSdzmatrix *Y_dzm = NULL;

   if (!X_dzm || X_dzm->real->nrows != _m || X_dzm->real->ncols != _n) {
      DestroyMatrix_dz(X_dzm);   /*  Destroys Y_dzm if already allocated memory to accommodate a possbible change of its dimension.   ansi-c*/
      Y_dzm = CreateMatrix_dz(_m, _n);
   }
   else  Y_dzm = X_dzm;

   Mreal = Y_dzm->real->M;
   Mimag = Y_dzm->imag->M;
   zp = Z;
   if ((incz=nrowsforZ-_m)<0)  fn_DisplayError("SubComplexMatrix2ZMatrix():  Number of rows for the input complex matrix Z must be greater that of the output Z matrix Y_dzm");

   for (_j=0; _j<_n; _j++) {
      for (_i=0; _i<_m; _i++) {
         *Mreal++ = zp->real;
         *Mimag++ = zp->imag;
         zp++;
      }
      zp += incz;
   }
   return (Y_dzm);
}


static void InitializeConstantMLK_Complex16(MKL_Complex16 *x_clx,  const int _n, const double c)
{
   int _i;

   for (_i=_n-1; _i>=0; _i--)
      x_clx[_i].real = x_clx[_i].imag = c;
}

static void InitializeConstantDouble(double *x_p,  const int _n, const double c)
{
   int _i;

   for (_i=_n-1; _i>=0; _i--)  x_p[_i] = c;
}

static void ConverteZeroSquareMatrix2RealDiagonalMLK_Complex16(MKL_Complex16 *x_pc,  const int _n, const double c)
{
/*     //Written by TZ, 03/08/06.   ansi-c*/
/*     //Output:   ansi-c*/
/*     //  x_pc: _n-by-_n, with the diagonal   ansi-c*/
/*     //Inputs:   ansi-c*/
/*     //  _n: dimension of x_pc so that x_pc is _n-by-_n.   ansi-c*/
/*     //  x_pc: _n-by-_n, all initialized to zeros.   ansi-c*/
   int _i;
   int np1 = _n+1;

   for (_i=_n*_n-1; _i>=0; _i -= np1)
      x_pc[_i].real = x_pc[_i].imag = c;
}

