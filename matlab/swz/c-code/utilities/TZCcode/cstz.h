#ifndef __CSTZ_H__
#define __CSTZ_H__
   #include "tzmatlab.h"
   #include "switch_opt.h"   //DW's Markov-switching routines, only used by gradcd_timet() and ComputeCovarianceFromOuterProduct().


   typedef struct {
           double bound;  //Real bounds to avoid extreme values that may make the P2 algorithm fail.
           double *p;  //5-by-1 probabilities as {0.0, p/2, p, (1+p)/2, 1.0}.
           double *q;  //5-by-1 quantiles.  Only q[2] is used as an estimate of p[2]-quantile or p-quantile.
           int *m;  //5-by-1 markers.
           int cnt;
           int ndeg;  //Number of exceptions such as degenerate numbers like inf.
   } TSdp2m5;
   typedef struct {
           TSdp2m5 **v;
           int n;
   } TSdvectorp2m5;
   typedef struct {
           TSdp2m5 **M;
           int nrows;
           int ncols;
   } TSdmatrixp2m5;
   typedef struct {
           TSdmatrixp2m5 **C;
           int ncells;
   } TSdcellp2m5;
   typedef struct {
           TSdcellp2m5 **F;
           int ndims;
   } TSdfourthp2m5;
   TSdp2m5 *CreateP2m5(const double p, const double bound);
   TSdp2m5 *DestroyP2m5(TSdp2m5 *x_dp2m5);
   TSdvectorp2m5 *CreateVectorP2m5(const int n, const double p, const double bound);
   TSdvectorp2m5 *DestroyVectorP2m5(TSdvectorp2m5 *x_dvp2m5);
   TSdmatrixp2m5 *CreateMatrixP2m5(const int nrows, const int ncols, const double p, const double bound);
   TSdmatrixp2m5 *DestroyMatrixP2m5(TSdmatrixp2m5 *X_dmp2m5);
   TSdcellp2m5 *CreateCellP2m5(const TSivector *rows_iv, const TSivector *cols_iv, const double p, const double bound);
   TSdcellp2m5 *DestroyCellP2m5(TSdcellp2m5 *X_dcp2m5);
   TSdfourthp2m5 *CreateFourthP2m5(const int ndims, const TSivector *rows_iv, const TSivector *cols_iv, const double p, const double bound);
   TSdfourthp2m5 *DestroyFourthP2m5(TSdfourthp2m5 *X_d4p2m5);
   //
   int P2m5Update(TSdp2m5 *x_dp2m5, const double newval);
   void P2m5VectorUpdate(TSdvectorp2m5 *x_dvp2m5, const TSdvector *newval_dv);
   void P2m5MatrixUpdate(TSdmatrixp2m5 *X_dmp2m5, const TSdmatrix *newval_dm);
   void P2m5CellUpdate(TSdcellp2m5 *X_dcp2m5, const TSdcell *newval_dc);
   void P2m5FourthUpdate(TSdfourthp2m5 *X_d4p2m5, const TSdfourth *newval_d4);


   #if defined ( CSMINWEL_OPTIMIZATION )
      void fn_gradcd(double *g, double *x, int n, double grdh,
                     double (*fcn)(double *x, int n, double **args, int *dims),
                     double **args, int *dims);

      void fn_hesscd(double *H, double *x, int n, double grdh,
                     double (*fcn)(double *x, int n, double **args, int *dims),
                     double **args, int *dims);
   #elif defined ( IMSL_OPTIMIZATION )
      void fn_gradcd(double *g, double *x, int n, double grdh,
                  double fcn(int n, double *x));
      void fn_hesscd(double *H, double *x, int n, double grdh,
                  double fcn(int n, double *x));
   #endif

   //=== For the conjugate gradient method I or II
   void gradcd_gen(double *g, double *x, int n, double (*fcn)(double *x, int n), double *grdh, double f0);
   void gradfd_gen(double *g, double *x, int n, double (*fcn)(double *x, int n), double *grdh, double f0);

   //=== For computing inverse Hessian.
   void gradcd_timet(TSdvector *g_dv, TSdvector *x_dv, int t, struct TStateModel_tag *smodel_ps, double (*fcn)(double *x, int t, struct TStateModel_tag *smodel_ps), double grdh, double f0);
   TSdmatrix *ComputeHessianFromOuterProduct(TSdmatrix *Hessian_dm, struct TStateModel_tag *smodel_ps, TSdvector *xhat_dv);
   TSdmatrix *ComputeCovarianceFromOuterProduct(TSdmatrix *Omega_dm, struct TStateModel_tag *smodel_ps, TSdvector *xhat_dv);




   int next_permutation(int *first, int *last);

   //void fn_ergodp(double **aop, int *aod, mxArray *cp);
   void fn_cumsum(double **aos_v, int *aods_v, double *v, int d_v);
   int fn_cumsum_int(int *x_v, const int d_x_v);
   double fn_cumsum_lf(double *x_v, const int d_x_v);
   double fn_mean(const double *a_v, const int _n);


   //=== For sorting according to x_dv.
   void tz_sort(TSdvector *x_dv, char ad);
   void tz_sortindex_lf(TSivector *x_iv, TSdvector *base_dv, char ad);
   void tz_sortindex(TSivector *x_iv, TSvoidvector *base_voidv, char ad);      //??????Not fully tested yet.
   //+
   void tz_sort_matrix(TSdmatrix *X_dm, char ad, char rc);
   TSdvector *tz_prctile_matrix(TSdvector *z_dv, const double prc, TSdmatrix *Z_dm, const char rc);
   TSdvector *tz_mean_matrix(TSdvector *z_dv, TSdmatrix *Z_dm, const char rc);
   //--- The following 3 functions should be hided (static) but are made visible to accomodate the old code that uses these functions.
   void fn_SetBaseArrayForComp(TSdvector *x_dv);
   int fn_compare(const void *i1, const void *i2);
   int fn_compare2(const void *i1, const void *i2);


   //=== Normalization for VARs.
   void fn_wznormalization(TSdvector *wznmlz_dv, TSdmatrix *A0draw_dm, TSdmatrix *A0peak_dm);

   //=== Handling under or over flows with log values.
   typedef struct TSveclogsum_tag {
      //For a recurisve algorithm to compute the log of sum (and therefore log of mean).  See p.81a and p.105 in TVBAR Notes.
      int n;  //Number of sums, which is the dimension for N_iv, Ysum_dv, and ymax_dv.
      TSivector *N_iv;  //(N_1, ..., N_n).
      TSdvector *logsum_dv,     //(logofsum_1, ..., logofsum_n).
                *logmax_dv;  //(logmax_1, ..., logmax_n).
   } TSveclogsum;
   struct TSveclogsum_tag *CreateVeclogsum(int n);
   struct TSveclogsum_tag *DestroyVeclogsum(struct TSveclogsum_tag *);
   //
   void UpdateSumFor1st2ndMoments(TSdvector *x1stsum_dv, TSdmatrix *X2ndsum_dm, const TSdvector *xdraw_dv);
   int tz_update_logofsum(double *Y_N_dp, double *y_Nmax_dp, double ynew, int N);
   int fn_update_logofsum(int N, double ynew, double *Y_N_dp, double *y_Nmax_dp);
   double fn_replace_logofsumsbt(double *yold, double _a, double ynew, double _b);



   //---------------------------- Special functions and densities. ---------------------
   double fn_normalcdf(double x);
   double fn_normalinv(double p);  //Inverse of normal cdf.
   double fn_chi2inv(double p, double df);
        //p = int_{0}^{\infty} chi2pdf(t, df) dt
   double fn_betainv(double p, double _alpha, double _beta);
        //p = int_{0}^{\infty} betapdf(t, _alpha, _beta) dt where betapdf(t,_alpha,_beta) \propt t^{_alpha-1}*(1-t)^(_beta-1}.
   double fn_gammalog(double x);
        //log gamma (x) where gamma(n+1) = n! and gamma(x) = int_0^{\infty} e^{-t} t^{x-1} dt.
   double fn_betalog(double x, double y);
        //log beta(x, y) where beta(x, y) = gamma(x)*gamm(y)/gamma(x+y).
   //+ Density functions
   double tz_lognormalpdf(double _x, double _m, double _s);
   double tz_logbetapdf(double _x, double _a, double _b);
   double tz_loggammapdf(double _x, double _a, double _b);
   double tz_loginversegammapdf(double _x, double _a, double _b);


   //---------------------------- Some high-level VAR functions ---------------------
   void fn_lev2growthanual(TSdmatrix *levgro_dm, const TSdmatrix *levgrominus1_dm, const TSivector *indxlogper_iv);
   void fn_ctfals_givenshocks_sm(TSdmatrix *ctfalstran_dm, TSdvector *xprimeminus1_dv, const int bloc, const int eloc, const TSdmatrix *strshockstran_dm,
        const TSivector *S_Tdraw_iv, const TSdcell *Bsdraw_dc, const TSdcell *A0sdrawinv_dc, const TSivector *noshocks_iv);
   void fn_ctfals_sm(TSdmatrix *ctfalstran_dm, TSdvector *xprimeminus1_dv, const int bloc, const int eloc, const TSdmatrix *strshockstran_dm, const TSivector *Snfores_iv, const TSdcell *Bsdraw_dc, const TSdcell *A0sdrawinv_dc);
   void fn_ctfals_policyonly(TSdmatrix *ctfalstran_dm, TSdvector *xprimeminus1_dv, const int bloc, const int eloc, const TSdmatrix *strshockstran_dm, const TSivector *S_Tdraw_iv, const int statecon, const int selej, const TSdcell *A0sdraw_dc, const TSdcell *Apsdraw_dc);
   void fn_impulse(TSdmatrix *imftran_dm, const TSdmatrix *Bh_dm, const TSdmatrix *swishtran_dm, const int nlags, const int imsteps);
   TSdmatrix *tz_impulse2levels(TSdmatrix *imflev_dm, TSdmatrix *imf_dm, TSivector *vlist2levels_iv);
   //
   void DynamicResponsesAR(TSdvector *resps_dv, const double c0, const TSdvector *a1_dv);
   void DynamicResponsesForStructuralEquation(TSdmatrix *Resps_dm, const int loclv, const int nlags, const TSdvector *a0p_dv);



   //---------------------------- Some regular vector or matrix operations ---------------------
   double MinVector_lf(TSdvector *x_dv);
   TSdvector *ConvertVector2exp(TSdvector *y_dv, TSdvector *x_dv);  //y=exp(x): output; x: input.
   TSdvector *ConvertVector2log(TSdvector *y_dv, TSdvector *x_dv);  //y=log(x): output; x: input.
   double tz_normofvector(TSdvector *x_dv, double p);


   //---------------------------- Old Interface ---------------------
   double gammalog(double x);
        //log gamma (x) where gamma(n+1) = n! and gamma(x) = int_0^{\infty} e^{-t} t^{x-1} dt.




   //----------- Must keep the following forever. -------------
   /**
   typedef struct {
           double *p;  //5-by-1 probabilities as {0.0, p/2, p, (1+p)/2, 1.0}.
           double *q;  //5-by-1 quantiles.  Only q[2] is used as an estimate of p[2]-quantile or p-quantile.
           int *m;  //5-by-1 markers.
           int cnt;
           int ndeg;  //Number of exceptions such as degenerate numbers like inf.
   } TSdp2m5;
   typedef struct {
           TSdp2m5 **v;
           int n;
   } TSdvectorp2m5;
   typedef struct {
           TSdp2m5 **M;
           int nrows;
           int ncols;
   } TSdmatrixp2m5;
   typedef struct {
           TSdmatrixp2m5 **C;
           int ncells;
   } TSdcellp2m5;

   TSdp2m5 *CreateP2m5(const double p);
   TSdp2m5 *DestroyP2m5(TSdp2m5 *x_dp2m5);
   TSdvectorp2m5 *CreateVectorP2m5(const int n, const double p);
   TSdvectorp2m5 *DestroyVectorP2m5(TSdvectorp2m5 *x_dvp2m5);
   TSdmatrixp2m5 *CreateMatrixP2m5(const int nrows, const int ncols, const double p);
   TSdmatrixp2m5 *DestroyMatrixP2m5(TSdmatrixp2m5 *X_dmp2m5);
   TSdcellp2m5 *CreateCellP2m5(const TSivector *rows_iv, const TSivector *cols_iv, const double p);
   TSdcellp2m5 *DestroyCellP2m5(TSdcellp2m5 *X_dcp2m5);
   /**/
#endif
