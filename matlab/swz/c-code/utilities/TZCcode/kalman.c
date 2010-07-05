/*===============================================================================================================
 * Check $$$ for important notes.
 * Check <<>> for updating DW's new switch code or questions for DW.
 *
 *   kalcvf_urw():  the Kalman filter forward prediction specialized for only a univariate random walk (urw) process.
 *
 *   State space model is defined as follows:
 *       z(t+1) = z(t)+eta(t)     (state or transition equation)
 *       y(t) = x(t)'*z(t)+eps(t)     (observation or measurement equation)
 *   where for this function, eta and eps must be uncorrelated; y(t) must be 1-by-1. Note that
 *     x(t): k-by-1;
 *     z(t): k-by-1;
 *     eps(t): 1-by-1 and ~ N(0, sigma^2);
 *     eta(t):  ~ N(0, V) where V is a k-by-k covariance matrix.
 *
 *
 * Written by Tao Zha, May 2004.
 * Revised, May 2008;
=================================================================================================================*/

/**
//=== For debugging purpose.
if (1)
{
   double t_loglht;

   t_loglht = -(0.5*ny)*LOG2PI - 0.5*logdeterminant(Dtdata_dm) - 0.5*VectorDotVector(wny_dv, etdata_dv);
   fprintf(FPTR_DEBUG, " %10.5f\n", t_loglht);

   fprintf(FPTR_DEBUG, "%%st=%d, inpt=%d, and sti=%d\n", st, inpt, sti);

   fprintf(FPTR_DEBUG, "\n wP0_dv:\n");
   WriteVector(FPTR_DEBUG, wP0_dv, " %10.5f ");
   fprintf(FPTR_DEBUG, "\n Vt_dc->C[sti_v=%d]:\n", sti_v);
   WriteMatrix(FPTR_DEBUG, Vt_dc->C[sti_v], " %10.5f ");

   fflush(FPTR_DEBUG);
}
/**/


#include "kalman.h"

#include "modify_for_mex.h"

static int Update_et_Dt_1stapp(int t_1, struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps);
static int Updatekalfilms_1stapp(int inpt, struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, struct TStateModel_tag *smodel_ps);


TSkalcvfurw *CreateTSkalcvfurw(TFlearninguni *func, int T, int k, int tv)   /*  , int storeZ, int storeV)   ansi-c*/
{
   int _i;
/*     //===   ansi-c*/
   TSivector *rows_iv = NULL;
   TSivector *cols_iv = NULL;
/*     //---   ansi-c*/
   TSkalcvfurw *kalcvfurw_ps = tzMalloc(1, TSkalcvfurw);


   kalcvfurw_ps->indx_tvsigmasq = tv;
   kalcvfurw_ps->fss = T;
   kalcvfurw_ps->kx = k;

/*     //===   ansi-c*/
   kalcvfurw_ps->V_dm = CreateMatrix_lf(k, k);
   kalcvfurw_ps->ylhtran_dv = CreateVector_lf(T);
   kalcvfurw_ps->Xrhtran_dm = CreateMatrix_lf(k, T);
   kalcvfurw_ps->z10_dv = CreateVector_lf(k);
   kalcvfurw_ps->P10_dm = CreateMatrix_lf(k, k);

   kalcvfurw_ps->zupdate_dv = CreateVector_lf(k);
   kalcvfurw_ps->Zpredtran_dm = CreateMatrix_lf(k, T);
   kalcvfurw_ps->ylhtranpred_dv = CreateVector_lf(T);
/*     //   ansi-c*/
   rows_iv = CreateVector_int(T);
   cols_iv = CreateVector_int(T);
   for (_i=T-1; _i>=0; _i--)  rows_iv->v[_i] = cols_iv->v[_i] = k;
   kalcvfurw_ps->Ppred_dc = CreateCell_lf(rows_iv, cols_iv);
/*     //   if (!storeZ)  kalcvfurw_ps->Zpredtran_dm = (TSdmatrix *)NULL;   ansi-c*/
/*     //   else  kalcvfurw_ps->Zpredtran_dm = CreateMatrix_lf(k, T);   ansi-c*/
/*     //   if (!storeV)  kalcvfurw_ps->Ppred_dc = (TSdcell *)NULL;   ansi-c*/
/*     //   else {   ansi-c*/
/*     //      rows_iv = CreateVector_int(T);   ansi-c*/
/*     //      cols_iv = CreateVector_int(T);   ansi-c*/
/*     //      for (_i=T; _i>=0; _i--)  rows_iv->v[_i] = cols_iv->v[_i] = k;   ansi-c*/
/*     //      kalcvfurw_ps->Ppred_dc = CreateCell_lf(rows_iv, cols_iv);   ansi-c*/
/*     //   }   ansi-c*/

   DestroyVector_int(rows_iv);
   DestroyVector_int(cols_iv);
   return (kalcvfurw_ps);
}

TSkalcvfurw *DestroyTSkalcvfurw(TSkalcvfurw *kalcvfurw_ps)
{
   if (kalcvfurw_ps) {
      DestroyMatrix_lf(kalcvfurw_ps->V_dm);
      DestroyVector_lf(kalcvfurw_ps->ylhtran_dv);
      DestroyMatrix_lf(kalcvfurw_ps->Xrhtran_dm);
      DestroyVector_lf(kalcvfurw_ps->z10_dv);
      DestroyMatrix_lf(kalcvfurw_ps->P10_dm);

      DestroyVector_lf(kalcvfurw_ps->zupdate_dv);
      DestroyMatrix_lf(kalcvfurw_ps->Zpredtran_dm);
      DestroyCell_lf(kalcvfurw_ps->Ppred_dc);
      DestroyVector_lf(kalcvfurw_ps->ylhtranpred_dv);

      swzFree(kalcvfurw_ps);
      return ((TSkalcvfurw *)NULL);
   }
   else  return (kalcvfurw_ps);
}


void kalcvf_urw(TSkalcvfurw *kalcvfurw_ps, void *dummy_ps)
{
/*     //See the notes of SWZ regarding the government's updating of the parameters in their Phillips-curve equation.   ansi-c*/
/*     //NOTE: make sure that the value of kalcvfurw_ps->sigmasq and other input values are given.   ansi-c*/
   int ti;
   double workd, workdenominv;
/*     //---   ansi-c*/
   int fss, kx;
   double sigmasq_fix = kalcvfurw_ps->sigmasq;
/*  //   double sigmasq;   ansi-c*/
   TSdmatrix *V_dm;
   TSdmatrix *Zpredtran_dm;
   TSdcell *Ppred_dc;
   TSdvector *ylhtran_dv;
   TSdmatrix *Xrhtran_dm;
/*     //===   ansi-c*/
   TSdvector *workkxby1_dv = NULL;   /*  kx-by-1.   ansi-c*/
/*  //   TSdvector *work1kxby1_dv = NULL;  //kx-by-1.   ansi-c*/
   TSdmatrix *workkxbykx_dm = NULL;   /*  kx-by-kx symmetric and positive positive.   ansi-c*/
/*  //   //===   ansi-c*/
/*  //   TSdvector *zbefore_dv = CreateVector_lf(kalcvfurw_ps->kx);   ansi-c*/
/*  //   TSdmatrix *Vbefore_dm = CreateMatrix_lf(kalcvfurw_ps->kx, kalcvfurw_ps->kx);   ansi-c*/
/*  //   TSdvector *zafter_dv = CreateVector_lf(kalcvfurw_ps->kx);   ansi-c*/
/*  //   TSdmatrix *Vafter_dm = CreateMatrix_lf(kalcvfurw_ps->kx, kalcvfurw_ps->kx);   ansi-c*/
   //******* WARNING: Some dangerous pointer movement to gain efficiency *******
//   double *yt_p;
//   double *Vbefore_p;
//   double *Vafter_p;
   TSdvector xt_sdv;
   TSdvector zbefore_sdv;
   //TSdmatrix Vbefore_sdm;
   TSdvector zafter_sdv;
   //TSdmatrix Vafter_sdm;


   if (!kalcvfurw_ps)  fn_DisplayError(".../kalcvf_urw(): the input argument kalcvfurw_ps must be created");
   if (!kalcvfurw_ps->V_dm || !kalcvfurw_ps->ylhtran_dv || !kalcvfurw_ps->Xrhtran_dm || !kalcvfurw_ps->z10_dv || !kalcvfurw_ps->P10_dm)
      fn_DisplayError(".../kalcvf_urw(): input arguments kalcvfurw_ps->V_dm, kalcvfurw_ps->ylhtran_dv, kalcvfurw_ps->Xrhtran_dm, kalcvfurw_ps->z10_dv, kalcvfurw_ps->P10_dm must be given legal values");
   if (!(kalcvfurw_ps->P10_dm->flag & (M_SU | M_SL)))  fn_DisplayError(".../kalcvf_urw(): the input argument kalcvfurw_ps->P10_dm must be symmetric");
   fss = kalcvfurw_ps->fss;
   kx = kalcvfurw_ps->kx;
   V_dm = kalcvfurw_ps->V_dm;
   Zpredtran_dm = kalcvfurw_ps->Zpredtran_dm;
   Ppred_dc = kalcvfurw_ps->Ppred_dc;
   ylhtran_dv = kalcvfurw_ps->ylhtran_dv;
   Xrhtran_dm = kalcvfurw_ps->Xrhtran_dm;
   //---
   xt_sdv.n = kx;
   xt_sdv.flag = V_DEF;
   zbefore_sdv.n = kx;
   zbefore_sdv.flag = V_DEF;
   zafter_sdv.n = kx;
   zafter_sdv.flag = V_DEF;

   //=== Memory allocation.
   workkxby1_dv = CreateVector_lf(kx);
   workkxbykx_dm = CreateMatrix_lf(kx, kx);


   //------- The first period (ti=0). -------
   zbefore_sdv.v = kalcvfurw_ps->z10_dv->v;
   zafter_sdv.v = Zpredtran_dm->M;
   xt_sdv.v = Xrhtran_dm->M;
   //---

   workd = ylhtran_dv->v[0] - (kalcvfurw_ps->ylhtranpred_dv->v[0]=VectorDotVector(&xt_sdv, &zbefore_sdv));   //y_t - x_t'*z_{t-1}.
   SymmatrixTimesVector(workkxby1_dv, kalcvfurw_ps->P10_dm, &xt_sdv, 1.0, 0.0);   //P_{t|t-1} x_t;

   if (!kalcvfurw_ps->indx_tvsigmasq)
      workdenominv = 1.0/(sigmasq_fix + VectorDotVector(&xt_sdv, workkxby1_dv));   //1/[sigma^2 + x_t' P_{t|t-1} x_t]
   else if (kalcvfurw_ps->indx_tvsigmasq == 1)        //See pp.37 and 37a in SWZ Learning NOTES.
      workdenominv = 1.0/(sigmasq_fix*square(kalcvfurw_ps->z10_dv->v[0]) + VectorDotVector(&xt_sdv, workkxby1_dv));   //1/[sigma^2 + x_t' P_{t|t-1} x_t];
   else {
      printf(".../kalman.c/kalcvf_urw(): Have not got time to deal with kalcvfurw_ps->indx_tvsigmasq defined in kalman.h other than 0 or 1");
      exit(EXIT_FAILURE);
   }


   //--- Updating z_{t+1|t}.
   CopyVector0(&zafter_sdv, &zbefore_sdv);
   VectorPlusMinusVectorUpdate(&zafter_sdv, workkxby1_dv, workd*workdenominv);  //z_{t+1|t} = z_{t|t-1} + P_{t|t-1} x_t [y_t - x_t'*z_{t-1}] / [sigma^2 + x_t' P_{t|t-1} x_t];
   //--- Updating P_{t+1|t}.
   CopyMatrix0(workkxbykx_dm, V_dm);
   VectorTimesSelf(workkxbykx_dm, workkxby1_dv, -workdenominv, 1.0, (V_dm->flag & M_SU) ? 'U' : 'L');
                                     // - P_{t|t-1}*x_t * xt'*P_{t|t-1} / [sigma^2 + x_t' P_{t|t-1} x_t] + V;
   MatrixPlusMatrix(Ppred_dc->C[0], kalcvfurw_ps->P10_dm, workkxbykx_dm);
                                     //P_{t|t-1} - P_{t|t-1}*x_t * xt'*P_{t|t-1} / [sigma^2 + x_t' P_{t|t-1} x_t] + V;
   Ppred_dc->C[0]->flag = M_GE | M_SU | M_SL;   //This is necessary because if P10_dm is initialized as diagonal, it will have M_GE | M_SU | M_SL | M_UT | M_LT,
                                          //  which is no longer true for workkxbykx_dm and therefore gives Ppred_dc->C[0] with M_GE only as a result of MatrixPlusMatrix().
                            //Done with all work* arrays.

   //------- The rest of the periods (ti=1:T-1). -------
   for (ti=1; ti<fss; ti++) {
      //NOTE: ti=0 has been taken care of outside of this loop.
      zbefore_sdv.v = Zpredtran_dm->M + (ti-1)*kx;
      zafter_sdv.v = Zpredtran_dm->M + ti*kx;
      xt_sdv.v = Xrhtran_dm->M + ti*kx;
      //---
      workd = ylhtran_dv->v[ti] - (kalcvfurw_ps->ylhtranpred_dv->v[ti]=VectorDotVector(&xt_sdv, &zbefore_sdv));   //y_t - x_t'*z_{t-1}.
      SymmatrixTimesVector(workkxby1_dv, Ppred_dc->C[ti-1], &xt_sdv, 1.0, 0.0);   //P_{t|t-1} x_t;
      if (!kalcvfurw_ps->indx_tvsigmasq)
         workdenominv = 1.0/(sigmasq_fix + VectorDotVector(&xt_sdv, workkxby1_dv));   //1/[sigma^2 + x_t' P_{t|t-1} x_t]
      else if (kalcvfurw_ps->indx_tvsigmasq == 1)    //See pp.37 and 37a in SWZ Learning NOTES.
         workdenominv = 1.0/(sigmasq_fix*square(zbefore_sdv.v[0]) + VectorDotVector(&xt_sdv, workkxby1_dv));   //1/[sigma^2 + x_t' P_{t|t-1} x_t]
      else {
         printf(".../kalman.c/kalcvf_urw(): Have not got time to deal with kalcvfurw_ps->indx_tvsigmasq defined in kalman.h other than 0 or 1");
         exit(EXIT_FAILURE);
      }
      //--- Updating z_{t+1|t}.
      CopyVector0(&zafter_sdv, &zbefore_sdv);
      VectorPlusMinusVectorUpdate(&zafter_sdv, workkxby1_dv, workd*workdenominv);  //z_{t+1|t} = z_{t|t-1} + P_{t|t-1} x_t [y_t - x_t'*z_{t-1}] / [sigma^2 + x_t' P_{t|t-1} x_t];
      //--- Updating P_{t+1|t}.
      CopyMatrix0(workkxbykx_dm, V_dm);
      VectorTimesSelf(workkxbykx_dm, workkxby1_dv, -workdenominv, 1.0, (V_dm->flag & M_SU) ? 'U' : 'L');
                                        // - P_{t|t-1}*x_t * xt'*P_{t|t-1} / [sigma^2 + x_t' P_{t|t-1} x_t] + V;
      MatrixPlusMatrix(Ppred_dc->C[ti], Ppred_dc->C[ti-1], workkxbykx_dm);
                                        //P_{t|t-1} - P_{t|t-1}*x_t * xt'*P_{t|t-1} / [sigma^2 + x_t' P_{t|t-1} x_t] + V;
      Ppred_dc->C[ti]->flag = M_GE | M_SU | M_SL;   //This is necessary because if P10_dm is initialized as diagonal, it will have M_GE | M_SU | M_SL | M_UT | M_LT,
                                             //  which is no longer true for workkxbykx_dm and therefore gives Ppred_dc->C[0] with M_GE only as a result of MatrixPlusMatrix().
                               //Done with all work* arrays.
   }
   CopyVector0(kalcvfurw_ps->zupdate_dv, &zafter_sdv);
   Zpredtran_dm->flag = M_GE;
   kalcvfurw_ps->ylhtranpred_dv->flag = V_DEF;

//   DestroyVector_lf(zbefore_dv);
//   DestroyMatrix_lf(Vbefore_dm);
//   DestroyVector_lf(zafter_dv);
//   DestroyMatrix_lf(Vafter_dm);

   DestroyVector_lf(workkxby1_dv);
//   DestroyVector_lf(work1kxby1_dv);
   DestroyMatrix_lf(workkxbykx_dm);
}



//-----------------------------------------------------------------------------------------------------------------------
//-- General constant (known-time-varying) Kalman filter for DSGE models.
//-----------------------------------------------------------------------------------------------------------------------
struct TSkalfiltv_tag *CreateTSkalfiltv(int ny, int nz, int T)
{
   int _i;
   //===
   TSivector *rows_iv = CreateVector_int(T);
   TSivector *cols_iv = CreateVector_int(T);
   //~~~ Creating the structure and initializing the NULL pointers.
   struct TSkalfiltv_tag *kalfiltv_ps = tzMalloc(1, struct TSkalfiltv_tag);


   //--- Default value.
   kalfiltv_ps->indxIni = 0;   //1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
                               //0: using the unconditional mean for any given regime at time 0.
   //--- Other assignments.
   kalfiltv_ps->ny = ny;
   kalfiltv_ps->nz = nz;
   kalfiltv_ps->T = T;



   //--------- Creates memory and assigns values.  The order matters.
   kalfiltv_ps->yt_dm = CreateMatrix_lf(ny, T);
   kalfiltv_ps->at_dm = CreateMatrix_lf(ny, T);
   //
   for (_i=T-1; _i>=0; _i--)
   {
      rows_iv->v[_i] = ny;
      cols_iv->v[_i] = nz;
   }
   rows_iv->flag = cols_iv->flag = V_DEF;
   kalfiltv_ps->Ht_dc = CreateCell_lf(rows_iv, cols_iv);
   //
   for (_i=T-1; _i>=0; _i--)
   {
      rows_iv->v[_i] = ny;
      cols_iv->v[_i] = ny;
   }
   kalfiltv_ps->Rt_dc = CreateCell_lf(rows_iv, cols_iv);
   //
   for (_i=T-1; _i>=0; _i--)
   {
      rows_iv->v[_i] = nz;
      cols_iv->v[_i] = ny;
   }
   kalfiltv_ps->Gt_dc = CreateCell_lf(rows_iv, cols_iv);
   //
   kalfiltv_ps->bt_dm = CreateMatrix_lf(nz, T);
   //
   for (_i=T-1; _i>=0; _i--)
   {
      rows_iv->v[_i] = nz;
      cols_iv->v[_i] = nz;
   }
   kalfiltv_ps->Ft_dc = CreateCell_lf(rows_iv, cols_iv);
   kalfiltv_ps->Vt_dc = CreateCell_lf(rows_iv, cols_iv);
   //
   kalfiltv_ps->z0_dv = CreateVector_lf(nz);
   kalfiltv_ps->P0_dm = CreateMatrix_lf(nz, nz);


   //---
   kalfiltv_ps->zt_tm1_dm = CreateMatrix_lf(nz, T);
   for (_i=T-1; _i>=0; _i--)
   {
      rows_iv->v[_i] = nz;
      cols_iv->v[_i] = nz;
   }
   kalfiltv_ps->Pt_tm1_dc = CreateCell_lf(rows_iv, cols_iv);


   //===
   DestroyVector_int(rows_iv);
   DestroyVector_int(cols_iv);

   return (kalfiltv_ps);

}
//---
struct TSkalfiltv_tag *DestroyTSkalfiltv(struct TSkalfiltv_tag *kalfiltv_ps)
{
   if (kalfiltv_ps)
   {
      //=== The order matters!
      DestroyMatrix_lf(kalfiltv_ps->yt_dm);
      DestroyMatrix_lf(kalfiltv_ps->at_dm);
      DestroyCell_lf(kalfiltv_ps->Ht_dc);
      DestroyCell_lf(kalfiltv_ps->Rt_dc);
      DestroyCell_lf(kalfiltv_ps->Gt_dc);
      //---
      DestroyMatrix_lf(kalfiltv_ps->bt_dm);
      DestroyCell_lf(kalfiltv_ps->Ft_dc);
      DestroyCell_lf(kalfiltv_ps->Vt_dc);
      //---
      DestroyVector_lf(kalfiltv_ps->z0_dv);
      DestroyMatrix_lf(kalfiltv_ps->P0_dm);
      //---
      DestroyMatrix_lf(kalfiltv_ps->zt_tm1_dm);
      DestroyCell_lf(kalfiltv_ps->Pt_tm1_dc);


      //---
      tzDestroy(kalfiltv_ps);  //Must be freed last!

      return ((struct TSkalfiltv_tag *)NULL);
   }
   else  return (kalfiltv_ps);
};


//-----------------------------------------------------------------------------------------------------------------------
//-- Inputs for filter for Markov-switching DSGE models at any time t.
//-----------------------------------------------------------------------------------------------------------------------
struct TSkalfilmsinputs_1stapp_tag *CreateTSkalfilmsinputs_1stapp(int ny, int nz, int nst, int T)
{
   //~~~ Creating the structure and initializing the NULL pointers.
   struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps = tzMalloc(1, struct TSkalfilmsinputs_1stapp_tag);

   //===
   TSivector *rows_iv = NULL;
   TSivector *cols_iv = NULL;

   //=== Default value.
   kalfilmsinputs_1stapp_ps->indxIni = 0;   //1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
                                            //0: using the unconditional mean for any given regime at time 0.
   kalfilmsinputs_1stapp_ps->indxDiffuse = 1;  //1: using the diffuse condition for z_{1|0} and P_{1|0} (default option), according to Koopman and Durbin, "Filtering and Smoothing of State Vector for Diffuse State-Space Models," J. of Time Series Analysis, Vol 24(1), pp.85-99.
                                               //0: using the unconditional moments.
   kalfilmsinputs_1stapp_ps->DiffuseScale = 100.0;
   kalfilmsinputs_1stapp_ps->ztm1_track = -1;
   kalfilmsinputs_1stapp_ps->dtm1_track = -1;

   //--- Other key assignments.
   kalfilmsinputs_1stapp_ps->ny = ny;
   kalfilmsinputs_1stapp_ps->nz = nz;
   kalfilmsinputs_1stapp_ps->nst = nst;
   kalfilmsinputs_1stapp_ps->T = T;

   //--------- Creates memory and assigns values.  The order matters.
   kalfilmsinputs_1stapp_ps->yt_dm = CreateMatrix_lf(ny, T);
   kalfilmsinputs_1stapp_ps->at_dm = CreateMatrix_lf(ny, nst);
   //
   rows_iv = CreateConstantVector_int(nst, ny);
   cols_iv = CreateConstantVector_int(nst, nz);
   kalfilmsinputs_1stapp_ps->Ht_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nst, ny);
   cols_iv = CreateConstantVector_int(nst, ny);
   kalfilmsinputs_1stapp_ps->Rt_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nst, nz);
   cols_iv = CreateConstantVector_int(nst, ny);
   kalfilmsinputs_1stapp_ps->Gt_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   kalfilmsinputs_1stapp_ps->bt_dm = CreateMatrix_lf(nz, nst);
   //
   rows_iv = CreateConstantVector_int(nst, nz);
   cols_iv = CreateConstantVector_int(nst, nz);
   kalfilmsinputs_1stapp_ps->Ft_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nst, nz);
   cols_iv = CreateConstantVector_int(nst, nz);
   kalfilmsinputs_1stapp_ps->Vt_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   kalfilmsinputs_1stapp_ps->z0_dm = CreateMatrix_lf(nz, nst); //nz-by-nst.
   kalfilmsinputs_1stapp_ps->z0_0_dm = CreateMatrix_lf(nz, nst); //nz-by-nst.
   //
   rows_iv = CreateConstantVector_int(nst, nz);
   cols_iv = CreateConstantVector_int(nst, nz);
   kalfilmsinputs_1stapp_ps->P0_dc = CreateCell_lf(rows_iv, cols_iv);  //nz-by-nz-by-nst.
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);

   //--- For output arguments.
   rows_iv = CreateConstantVector_int(T+1, nz);
   cols_iv = CreateConstantVector_int(T+1, nst);
   kalfilmsinputs_1stapp_ps->zt_tm1_dc = CreateCell_lf(rows_iv, cols_iv); //nz-by-nst-by-(T+1).
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nst, nz);
   cols_iv = CreateConstantVector_int(nst, nz);
   kalfilmsinputs_1stapp_ps->Pt_tm1_d4 = CreateFourth_lf(T+1, rows_iv, cols_iv); //nz-by-nz-by-nst-by-(T+1).
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nst, nz);
   cols_iv = CreateConstantVector_int(nst, ny);
   kalfilmsinputs_1stapp_ps->PHtran_tdata_d4 = CreateFourth_lf(T, rows_iv, cols_iv);  //nz-by-ny-by-nst-T, saved only for updating Kalman filter Updatekalfilms_1stapp().
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(T, ny);
   cols_iv = CreateConstantVector_int(T, nst);
   kalfilmsinputs_1stapp_ps->etdata_dc = CreateCell_lf(rows_iv, cols_iv); //ny-by-nst-by-T, used for updating Kalman filter Updatekalfilms_1stapp().
   rows_iv = CreateConstantVector_int(T, ny);
   cols_iv = CreateConstantVector_int(T, nst);
   //
   rows_iv = CreateConstantVector_int(nst, ny);
   cols_iv = CreateConstantVector_int(nst, ny);
   kalfilmsinputs_1stapp_ps->Dtdata_d4 = CreateFourth_lf(T, rows_iv, cols_iv);  //ny-by-ny-nst-by-T, used for updating Kalman filter Updatekalfilms_1stapp().
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);

   return (kalfilmsinputs_1stapp_ps);
}
//---
struct TSkalfilmsinputs_1stapp_tag *DestroyTSkalfilmsinputs_1stapp(struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps)
{
   if (kalfilmsinputs_1stapp_ps)
   {
      //=== The order matters!
      DestroyMatrix_lf(kalfilmsinputs_1stapp_ps->yt_dm);
      DestroyMatrix_lf(kalfilmsinputs_1stapp_ps->at_dm);
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->Ht_dc);
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->Rt_dc);
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->Gt_dc);
      //---
      DestroyMatrix_lf(kalfilmsinputs_1stapp_ps->bt_dm);
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->Ft_dc);
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->Vt_dc);
      //---
      DestroyMatrix_lf(kalfilmsinputs_1stapp_ps->z0_dm);
      DestroyMatrix_lf(kalfilmsinputs_1stapp_ps->z0_0_dm);
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->P0_dc);
      //---
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->zt_tm1_dc);
      DestroyFourth_lf(kalfilmsinputs_1stapp_ps->Pt_tm1_d4);
      DestroyFourth_lf(kalfilmsinputs_1stapp_ps->PHtran_tdata_d4);
      DestroyCell_lf(kalfilmsinputs_1stapp_ps->etdata_dc);
      DestroyFourth_lf(kalfilmsinputs_1stapp_ps->Dtdata_d4);
      //---
      tzDestroy(kalfilmsinputs_1stapp_ps);  //Must be freed last!

      return ((struct TSkalfilmsinputs_1stapp_tag *)NULL);
   }
   else  return (kalfilmsinputs_1stapp_ps);
};


//-----------------------------------------------------------------------------------------------------------------------
//-- OLD Code: Inputs for filter for Markov-switching DSGE models at any time t.
//-----------------------------------------------------------------------------------------------------------------------
struct TSkalfilmsinputs_tag *CreateTSkalfilmsinputs(int ny, int nz, int nRc, int nRstc, int nRv, int indxIndRegimes, int T)
{
   //~~~ Creating the structure and initializing the NULL pointers.
   struct TSkalfilmsinputs_tag *kalfilmsinputs_ps = tzMalloc(1, struct TSkalfilmsinputs_tag);

   //===
   TSivector *rows_iv = NULL;
   TSivector *cols_iv = NULL;

   //--- Default value.
   kalfilmsinputs_ps->indxIni = 0;   //1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
                                     //0: using the unconditional mean for any given regime at time 0.
   //--- Other assignments.
   kalfilmsinputs_ps->ny = ny;
   kalfilmsinputs_ps->nz = nz;
   kalfilmsinputs_ps->nRc = nRc;
   kalfilmsinputs_ps->nRstc = nRstc;
   kalfilmsinputs_ps->nRv = nRv;
   kalfilmsinputs_ps->indxIndRegimes = indxIndRegimes;
   kalfilmsinputs_ps->T = T;


   //--------- Creates memory and assigns values.  The order matters.
   kalfilmsinputs_ps->yt_dm = CreateMatrix_lf(ny, T);
   kalfilmsinputs_ps->at_dm = CreateMatrix_lf(ny, nRc);
   //
   rows_iv = CreateConstantVector_int(nRc, ny);
   cols_iv = CreateConstantVector_int(nRc, nz);
   kalfilmsinputs_ps->Ht_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nRv, ny);
   cols_iv = CreateConstantVector_int(nRv, ny);
   kalfilmsinputs_ps->Rt_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nRv, nz);
   cols_iv = CreateConstantVector_int(nRv, ny);
   kalfilmsinputs_ps->Gt_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   kalfilmsinputs_ps->bt_dm = CreateMatrix_lf(nz, nRc);
   //
   rows_iv = CreateConstantVector_int(nRc, nz);
   cols_iv = CreateConstantVector_int(nRc, nz);
   kalfilmsinputs_ps->Ft_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   rows_iv = CreateConstantVector_int(nRv, nz);
   cols_iv = CreateConstantVector_int(nRv, nz);
   kalfilmsinputs_ps->Vt_dc = CreateCell_lf(rows_iv, cols_iv);
   rows_iv = DestroyVector_int(rows_iv);
   cols_iv = DestroyVector_int(cols_iv);
   //
   if (indxIndRegimes)
   {
      kalfilmsinputs_ps->z0_dm = CreateMatrix_lf(nz, nRc*nRv); //nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nRv if indxIndRegimes == 0.
      //
      rows_iv = CreateConstantVector_int(nRc*nRv, nz);
      cols_iv = CreateConstantVector_int(nRc*nRv, nz);
      kalfilmsinputs_ps->P0_dc = CreateCell_lf(rows_iv, cols_iv);  //nz-by-nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nz-by-nRv if indxIndRegimes == 0.
      rows_iv = DestroyVector_int(rows_iv);
      cols_iv = DestroyVector_int(cols_iv);
   }
   else
   {
      if (nRstc != nRv)  fn_DisplayError("kalman.c/CreateTSkalfilmsinputs(): nRstc must equal to nRv when indxIndRegimes==0");
      kalfilmsinputs_ps->z0_dm = CreateMatrix_lf(nz, nRv); //nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nRv if indxIndRegimes == 0.
      //
      rows_iv = CreateConstantVector_int(nRv, nz);
      cols_iv = CreateConstantVector_int(nRv, nz);
      kalfilmsinputs_ps->P0_dc = CreateCell_lf(rows_iv, cols_iv);  //nz-by-nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nz-by-nRv if indxIndRegimes == 0.
      rows_iv = DestroyVector_int(rows_iv);
      cols_iv = DestroyVector_int(cols_iv);
   }
   //--- For output arguments.
   if (indxIndRegimes)
   {
      rows_iv = CreateConstantVector_int(T, nz);
      cols_iv = CreateConstantVector_int(T, nRc*nRv);
      kalfilmsinputs_ps->zt_tm1_dc = CreateCell_lf(rows_iv, cols_iv); //nz-by-nRc*nRv-by-T if indxIndRegimes==1, nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
      rows_iv = DestroyVector_int(rows_iv);
      cols_iv = DestroyVector_int(cols_iv);
      //
      rows_iv = CreateConstantVector_int(nRc*nRv, nz);
      cols_iv = CreateConstantVector_int(nRc*nRv, nz);
      kalfilmsinputs_ps->Pt_tm1_d4 = CreateFourth_lf(T, rows_iv, cols_iv); //nz-by-nz-by-nRc*nRv-T if indxIndRegimes==1, nz-by-nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
      rows_iv = DestroyVector_int(rows_iv);
      cols_iv = DestroyVector_int(cols_iv);
   }
   else
   {
      if (nRstc != nRv)  fn_DisplayError("kalman.c/CreateTSkalfilmsinputs(): nRstc must equal to nRv when indxIndRegimes==0");
      rows_iv = CreateConstantVector_int(T, nz);
      cols_iv = CreateConstantVector_int(T, nRv);
      kalfilmsinputs_ps->zt_tm1_dc = CreateCell_lf(rows_iv, cols_iv); //nz-by-nRc*nRv-by-T if indxIndRegimes==1, nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
      rows_iv = DestroyVector_int(rows_iv);
      cols_iv = DestroyVector_int(cols_iv);
      //
      rows_iv = CreateConstantVector_int(nRv, nz);
      cols_iv = CreateConstantVector_int(nRv, nz);
      kalfilmsinputs_ps->Pt_tm1_d4 = CreateFourth_lf(T, rows_iv, cols_iv); //nz-by-nz-by-nRc*nRv-T if indxIndRegimes==1, nz-by-nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
      rows_iv = DestroyVector_int(rows_iv);
      cols_iv = DestroyVector_int(cols_iv);
   }


   //===
   DestroyVector_int(rows_iv);
   DestroyVector_int(cols_iv);

   return (kalfilmsinputs_ps);

}
//---
struct TSkalfilmsinputs_tag *DestroyTSkalfilmsinputs(struct TSkalfilmsinputs_tag *kalfilmsinputs_ps)
{
   if (kalfilmsinputs_ps)
   {
      //=== The order matters!
      DestroyMatrix_lf(kalfilmsinputs_ps->yt_dm);
      DestroyMatrix_lf(kalfilmsinputs_ps->at_dm);
      DestroyCell_lf(kalfilmsinputs_ps->Ht_dc);
      DestroyCell_lf(kalfilmsinputs_ps->Rt_dc);
      DestroyCell_lf(kalfilmsinputs_ps->Gt_dc);
      //---
      DestroyMatrix_lf(kalfilmsinputs_ps->bt_dm);
      DestroyCell_lf(kalfilmsinputs_ps->Ft_dc);
      DestroyCell_lf(kalfilmsinputs_ps->Vt_dc);
      //---
      DestroyMatrix_lf(kalfilmsinputs_ps->z0_dm);
      DestroyCell_lf(kalfilmsinputs_ps->P0_dc);
      //---
      DestroyCell_lf(kalfilmsinputs_ps->zt_tm1_dc);
      DestroyFourth_lf(kalfilmsinputs_ps->Pt_tm1_d4);
      //---
      tzDestroy(kalfilmsinputs_ps);  //Must be freed last!

      return ((struct TSkalfilmsinputs_tag *)NULL);
   }
   else  return (kalfilmsinputs_ps);
};


#define LOG2PI  (1.837877066409345e+000)   //log(2*pi)
//-----------------------------------------------------
//-- Constant-parameters (known-time-varying) Kalman filter
//-----------------------------------------------------
double tz_kalfiltv(struct TSkalfiltv_tag *kalfiltv_ps)
{
   //General constant (known-time-varying) Kalman filter for DSGE models (conditional on all the parameters).
   //  It computes a sequence of one-step predictions and their covariance matrices, and the log likelihood.
   //  The function uses a forward recursion algorithm.  See also the Matlab function fn_kalfil_tv.m
   //
   //   State space model is defined as follows:
   //       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)
   //       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)
   //     where a(t), H(t), b(t), and F(t) depend on s_t that follows a Markov-chain process and are taken as given.
   //
   //   Inputs are as follows:
   //      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].
   //        a is an n_y-by-T matrix of time-varying input vectors in the measurement equation.
   //        H is an n_y-by-n_z-by-T 3-D of time-varying matrices in the measurement equation.
   //        R is an n_y-by-n_y-by-T 3-D of time-varying covariance matrices for the error in the measurement equation.
   //        G is an n_z-by-n_y-by-T 3-D of time-varying E(eta_t * eps_t').
   //        ------
   //        b is an n_z-by-T matrix of time-varying input vectors in the state equation with b(:,1) as an initial condition.
   //        F is an n_z-by-n_z-by-T 3-D of time-varying transition matrices in the state equation with F(:,:,1) as an initial condition.
   //        V is an n_z-by-n_z-by-T 3-D of time-varying covariance matrices for the error in the state equation with V(:,:,1) as an initial condition.
   //        ------
   //        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
   //                 0: using the unconditional mean for any given regime at time 0.
   //        z0 is an n_z-by-1 vector of initial condition when indxIni=1. (Do not enter if indxIni=0.)
   //        P0 is an n_z-by-n_z matrix of initial condition when indxIni=1. (Do not enter if indxIni=0.)
   //
   //   Outputs are as follows:
   //      loglh is a value of the log likelihood function of the state-space model
   //                                under the assumption that errors are multivariate Gaussian.
   //      zt_tm1 is an n_z-by-T matrices of one-step predicted state vectors with z0_0m1 as an initial condition (base-0 first element)
   //                         and with z_{T|T-1} as the last element.  Thus, we can use it as a base-1 vector.
   //      Pt_tm1 is an n_z-by-n_z-by-T 3-D of covariance matrices of zt_tm1 with P0_0m1 as though it were a initial condition
   //                         and with P_{T|T-1} as the last element.  Thus, we can use it as though it were a base-1 cell.
   //
   //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:
   //             z0_0m1 = (I-F(:,:,1))\b(:,1)
   //        vec(P0_0m1) = (I-kron(F(:,:,1),F(:,:,1)))\vec(V(:,:,1))
   //   Note that all eigenvalues of the matrix F(:,:,1) are inside the unit circle when the state-space model is bounded (stationary).
   //
   //   March 2007, written by Tao Zha
   //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.

   int T = kalfiltv_ps->T;
   int Tp1 = T + 1;
   int ny = kalfiltv_ps->ny;
   int nz = kalfiltv_ps->nz;
   int indx_badlh = 0;   //1: bad likelihood with, say, -infinity of the LH value.
   int tdata, ti;
   //--- Work arguments.
   int nz2 = square(nz);
   TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *Wnz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdmatrix *W2nz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdvector *wP0_dv = CreateVector_lf(nz2);
   //+
   TSdvector yt_sdv, at_sdv, zt_tm1_sdv, ztp1_t_sdv, btp1_sdv;  //double loglh_tdata;  //logdetDtdata.
   TSdvector *wny_dv = CreateVector_lf(ny);
   TSdmatrix *Wnzbyny_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *W2nzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *PHtran_tdata_dm = CreateMatrix_lf(nz,ny);
   TSdvector *etdata_dv = CreateVector_lf(ny);
   TSdmatrix *Dtdata_dm = CreateMatrix_lf(ny,ny);
   TSdmatrix *Kt_tdata0_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *Kt_tdata_dm = CreateMatrix_lf(nz,ny);
   //--- For eigenvalue decompositions
   int ki;
   int errflag;
   double eigmax, logdet_Dtdata;
   TSdzvector *evals_dzv = NULL;
   TSdvector *evals_abs_dv = NULL;  //Absolute eigenvalues.
   //--- Input arguments.
   TSdmatrix *yt_dm = kalfiltv_ps->yt_dm;   //ny-by-T.
   TSdmatrix *at_dm = kalfiltv_ps->at_dm;   //ny-by-T.
   TSdcell *Ht_dc = kalfiltv_ps->Ht_dc;   //ny-by-nz-by-T.
   TSdcell *Rt_dc = kalfiltv_ps->Rt_dc;   //ny-by-ny-by-T.  Covariance matrix for the measurement equation.
   TSdcell *Gt_dc = kalfiltv_ps->Gt_dc;   //nz-by-ny-by-T.  Cross-covariance.
   //
   TSdmatrix *bt_dm = kalfiltv_ps->bt_dm;   //nz-by-T.
   TSdcell *Ft_dc = kalfiltv_ps->Ft_dc;   //nz-by-nz-by-T.
   TSdcell *Vt_dc = kalfiltv_ps->Vt_dc;   //nz-by-nz-by-T.  Covariance matrix for the state equation.
   //
   TSdvector *z0_dv = kalfiltv_ps->z0_dv;  //nz-by-1;
   TSdmatrix *P0_dm = kalfiltv_ps->P0_dm;   //nz-by-nz.
   //--- Output arguments.
   double loglh;   //log likelihood.
   TSdmatrix *zt_tm1_dm = kalfiltv_ps->zt_tm1_dm;  //nz-by-T.
   TSdcell *Pt_tm1_dc = kalfiltv_ps->Pt_tm1_dc;   //nz-by-nz-T.



   //=== Initializing.
   if (!kalfiltv_ps->indxIni)
   {
      InitializeDiagonalMatrix_lf(Wnzbynz_dm, 1.0);  //To be used for I(nz) -
      InitializeDiagonalMatrix_lf(Wnz2bynz2_dm, 1.0);  //To be used for I(nz2) -

      //=== Eigenanalysis to determine the roots to ensure boundedness.
      evals_dzv = CreateVector_dz(nz);
      evals_abs_dv = CreateVector_lf(nz);
      errflag = eigrgen(evals_dzv, (TSdzmatrix *)NULL, (TSdzmatrix *)NULL, Ft_dc->C[0]);
      if (errflag)  fn_DisplayError("tz_kalfiltv() in kalman.c: eigen decomposition failed");
      for (ki=nz-1; ki>=0; ki--)  evals_abs_dv->v[ki] = sqrt(square(evals_dzv->real->v[ki]) + square(evals_dzv->imag->v[ki]));
      evals_abs_dv->flag = V_DEF;
      eigmax = MaxVector(evals_abs_dv);
      if (eigmax < (1.0+1.0e-14))
      {
         //--- Getting z0_dv: zt_tm1(:,1) = (eye(n_z)-F(:,:,1))\b(:,1);
         MatrixMinusMatrix(Wnzbynz_dm, Wnzbynz_dm, Ft_dc->C[0]);
         CopySubmatrix2vector(z0_dv, 0, bt_dm, 0, 0, bt_dm->nrows);
         bdivA_rgens(z0_dv, z0_dv, '\\', Wnzbynz_dm);
                   //Done with Wnzbynz_dm.
         //--- Getting P0_dm: Pt_tm1(:,:,1) = reshape((eye(n_z^2)-kron(F(:,:,1),F(:,:,1)))\V1(:),n_z,n_z);
         tz_kron(W2nz2bynz2_dm, Ft_dc->C[0], Ft_dc->C[0]);
         MatrixMinusMatrix(Wnz2bynz2_dm, Wnz2bynz2_dm, W2nz2bynz2_dm);
         CopySubmatrix2vector(wP0_dv, 0, Vt_dc->C[0], 0, 0, nz2);
         bdivA_rgens(wP0_dv, wP0_dv, '\\', Wnz2bynz2_dm);
         CopySubvector2matrix_unr(P0_dm, 0, 0, wP0_dv, 0, nz2);
                    //Done with all w*_dv and W*_dm.
      }
      else
      {
         printf("Fatal error: tz_kalfiltv() in kalman.c: the system is non-stationary solutions\n"
                         "    and the initial conditions must be supplied by, say, input arguments");
         fflush(stdout);
         exit( EXIT_FAILURE );
      }
   }
   CopySubvector2matrix(zt_tm1_dm, 0, 0, z0_dv, 0, z0_dv->n);
   CopyMatrix0(Pt_tm1_dc->C[0], P0_dm);

   //====== See p.002 in LiuWZ. ======
   at_sdv.n = yt_sdv.n = yt_dm->nrows;
   at_sdv.flag = yt_sdv.flag = V_DEF;
   zt_tm1_sdv.n = ztp1_t_sdv.n = zt_tm1_dm->nrows;
   zt_tm1_sdv.flag = ztp1_t_sdv.flag = V_DEF;
   btp1_sdv.n = bt_dm->nrows;
   btp1_sdv.flag = V_DEF;
   loglh = 0.0;
   for (tdata=0; tdata<T; tdata++ )
   {
      //Base-0 timing.
      ti = tdata + 1;  //Next period.

      //--- Setup.
      MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_dc->C[tdata], Ht_dc->C[tdata], 1.0, 0.0, 'N', 'T');

      //--- Data.
      //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata;
      yt_sdv.v = yt_dm->M + tdata*yt_dm->nrows;
      at_sdv.v = at_dm->M + tdata*at_dm->nrows;
      zt_tm1_sdv.v = zt_tm1_dm->M + tdata*zt_tm1_dm->nrows;
      VectorMinusVector(etdata_dv, &yt_sdv, &at_sdv);
      MatrixTimesVector(etdata_dv, Ht_dc->C[tdata], &zt_tm1_sdv, -1.0, 1.0, 'N');
      //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tdata);
      CopyMatrix0(Dtdata_dm, Rt_dc->C[tdata]);
      MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[tdata], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
      ScalarTimesMatrixSquare(Dtdata_dm, 0.5, Dtdata_dm, 'T', 0.5);  //Making it symmetric against some rounding errors.
                         //This making-symmetric is very IMPORTANT; otherwise, we will get the matrix being singular message
                         //    and eigenvalues being negative for the SPD matrix, etc.  Then the likelihood becomes either
                         //    a bad number or a complex number.
      Dtdata_dm->flag = Dtdata_dm->flag | M_SU | M_SL;

      //--- Forming the log likelihood.
      if (!isfinite(logdet_Dtdata=logdetspd(Dtdata_dm)))  return (kalfiltv_ps->loglh = -NEARINFINITY);
      bdivA_rgens(wny_dv, etdata_dv, '/', Dtdata_dm);
      loglh += -(0.5*ny)*LOG2PI - 0.5*logdet_Dtdata - 0.5*VectorDotVector(wny_dv, etdata_dv);
      //loglh += -(0.5*ny)*LOG2PI - 0.5*logdeterminant(Dtdata_dm) - 0.5*VectorDotVector(wny_dv, etdata_dv);
                            //Done with all w*_dv.


      //--- Updating zt_tm1_dm and Pt_tm1_dc by ztp1_t_sdv and Pt_tm1_dc->C[ti].
      if (ti<T)
      {
         //Updating only up to tdata=T-2.  The values at ti=T or tdata=T-1 will not be used in the likelihood function.

         //- Kt_tdata = (Ft*PHtran_tdata+G(:,:,tdata))/Dtdata;
         CopyMatrix0(Kt_tdata0_dm, Gt_dc->C[tdata]);
         MatrixTimesMatrix(Kt_tdata0_dm, Ft_dc->C[ti], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
         BdivA_rrect(Kt_tdata_dm, Kt_tdata0_dm, '/', Dtdata_dm);
         //+ zt_tm1(:,t) = b(:,t) + Ft*zt_tm1(:,tdata) + Kt_tdata*etdata;
         ztp1_t_sdv.v = zt_tm1_dm->M + ti*zt_tm1_dm->nrows;
         MatrixTimesVector(&ztp1_t_sdv, Ft_dc->C[ti], &zt_tm1_sdv, 1.0, 0.0, 'N');
         MatrixTimesVector(&ztp1_t_sdv, Kt_tdata_dm, etdata_dv, 1.0, 1.0, 'N');
         btp1_sdv.v = bt_dm->M + ti*btp1_sdv.n;
         VectorPlusMinusVectorUpdate(&ztp1_t_sdv, &btp1_sdv, 1.0);
         //+ Pt_tm1(:,:,t) = Ft*Ptdata*Fttran - Kt_tdata*Dtdata*Kt_tdatatran + V(:,:,t);
         CopyMatrix0(Pt_tm1_dc->C[ti], Vt_dc->C[ti]);
         MatrixTimesMatrix(Wnzbyny_dm, Kt_tdata_dm, Dtdata_dm, 1.0, 0.0, 'N', 'N');
         MatrixTimesMatrix(Wnzbynz_dm, Wnzbyny_dm, Kt_tdata_dm, 1.0, 0.0, 'N', 'T');
         MatrixPlusMinusMatrixUpdate(Pt_tm1_dc->C[ti], Wnzbynz_dm, -1.0);
                               //Done with all W*_dm.
         MatrixTimesMatrix(Wnzbynz_dm, Ft_dc->C[ti], Pt_tm1_dc->C[tdata], 1.0, 0.0, 'N', 'N');
         MatrixTimesMatrix(W2nzbynz_dm, Wnzbynz_dm, Ft_dc->C[ti], 1.0, 0.0, 'N', 'T');
         MatrixPlusMatrixUpdate(Pt_tm1_dc->C[ti], W2nzbynz_dm);
                               //Done with all W*_dm.
      }
   }
   zt_tm1_dm->flag = M_GE;

   //===
   DestroyVector_dz(evals_dzv);
   DestroyVector_lf(evals_abs_dv);
   DestroyMatrix_lf(Wnzbynz_dm);
   DestroyMatrix_lf(Wnz2bynz2_dm);
   DestroyMatrix_lf(W2nz2bynz2_dm);
   DestroyVector_lf(wP0_dv);
   //
   DestroyVector_lf(wny_dv);
   DestroyMatrix_lf(Wnzbyny_dm);
   DestroyMatrix_lf(W2nzbynz_dm);
   DestroyMatrix_lf(PHtran_tdata_dm);
   DestroyVector_lf(etdata_dv);
   DestroyMatrix_lf(Dtdata_dm);
   DestroyMatrix_lf(Kt_tdata0_dm);
   DestroyMatrix_lf(Kt_tdata_dm);

   return (kalfiltv_ps->loglh = loglh);
}
/**
double tz_kalfiltv(struct TSkalfiltv_tag *kalfiltv_ps)
{
   //This function is used to test tz_logTimetCondLH_kalfiltv().
   int T = kalfiltv_ps->T;
   int tdata;
   double loglh;

   loglh = 0.0;
   for (tdata=0; tdata<T; tdata++)  loglh += tz_logTimetCondLH_kalfiltv(0, tdata+1, kalfiltv_ps);

   return (loglh);
}
/**/
/*  //-----------------------------------------------------   ansi-c*/
/*  //-- Updating Kalman filter at time t for constant-parameters (or known-time-varying) Kalman filter.   ansi-c*/
/*  //-----------------------------------------------------   ansi-c*/
double tz_logTimetCondLH_kalfiltv(int st, int inpt, struct TSkalfiltv_tag *kalfiltv_ps)
{
/*     //st: base-0 grand regime at time t, which is just a dummy for this constant-parameter function in order to use   ansi-c*/
/*     //       Waggoner's automatic functions.   ansi-c*/
/*     //inpt: base-1 in the sense that inpt>=1 to deal with the time series situation where S_T is (T+1)-by-1 and Y_T is T+nlags_max-by-1.   ansi-c*/
/*     //      The 1st element for S_T is S_T[1] while S_T[0] is s_0 (initial condition).   ansi-c*/
/*     //      The 1st element for Y_T, however, is Y_T[nlags_max+1-1].   ansi-c*/
/*     //See (42.3) on p.42 in the SWZII NOTES.   ansi-c*/
/*     //   ansi-c*/
/*     //log LH at time t for constant (known-time-varying) Kalman-filter DSGE models (conditional on all the parameters).   ansi-c*/
/*     //  It computes a sequence of one-step predictions and their covariance matrices, and the log likelihood at time t.   ansi-c*/
/*     //  The function uses a forward recursion algorithm.  See also the Matlab function fn_kalfil_tv.m   ansi-c*/
/*     //   ansi-c*/
/*     //   State space model is defined as follows:   ansi-c*/
/*     //       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)   ansi-c*/
/*     //       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)   ansi-c*/
/*     //     where a(t), H(t), b(t), and F(t) depend on s_t that follows a Markov-chain process and are taken as given.   ansi-c*/
/*     //   ansi-c*/
/*     //   Inputs are as follows:   ansi-c*/
/*     //      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].   ansi-c*/
/*     //        a is an n_y-by-T matrix of time-varying input vectors in the measurement equation.   ansi-c*/
/*     //        H is an n_y-by-n_z-by-T 3-D of time-varying matrices in the measurement equation.   ansi-c*/
/*     //        R is an n_y-by-n_y-by-T 3-D of time-varying covariance matrices for the error in the measurement equation.   ansi-c*/
/*     //        G is an n_z-by-n_y-by-T 3-D of time-varying E(eta_t * eps_t').   ansi-c*/
/*     //        ------   ansi-c*/
/*     //        b is an n_z-by-T matrix of time-varying input vectors in the state equation with b(:,1) as an initial condition.   ansi-c*/
/*     //        F is an n_z-by-n_z-by-T 3-D of time-varying transition matrices in the state equation with F(:,:,1) as an initial condition.   ansi-c*/
/*     //        V is an n_z-by-n_z-by-T 3-D of time-varying covariance matrices for the error in the state equation with V(:,:,1) as an initial condition.   ansi-c*/
/*     //        ------   ansi-c*/
/*     //        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;   ansi-c*/
/*     //                 0: using the unconditional mean for any given regime at time 0.   ansi-c*/
/*     //        z0 is an n_z-by-1 vector of initial condition when indxIni=1. (Value to be assigned if indxIni=0.)   ansi-c*/
/*     //        P0 is an n_z-by-n_z matrix of initial condition when indxIni=1. (Value to be assigned if indxIni=0.)   ansi-c*/
/*     //   ansi-c*/
/*     //   Outputs are as follows:   ansi-c*/
/*     //      loglh is a value of the log likelihood function of the state-space model   ansi-c*/
/*     //                                under the assumption that errors are multivariate Gaussian.   ansi-c*/
/*     //      zt_tm1 is an n_z-by-T matrices of one-step predicted state vectors with z0_0m1 as an initial condition (base-0 first element)   ansi-c*/
/*     //                         and with z_{T|T-1} as the last element.  Thus, we can use it as a base-1 vector.   ansi-c*/
/*     //      Pt_tm1 is an n_z-by-n_z-by-T 3-D of covariance matrices of zt_tm1 with P0_0m1 as though it were a initial condition   ansi-c*/
/*     //                         and with P_{T|T-1} as the last element.  Thus, we can use it as though it were a base-1 cell.   ansi-c*/
/*     //   ansi-c*/
/*     //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:   ansi-c*/
/*     //             z0_0m1 = (I-F(:,:,1))\b(:,1)   ansi-c*/
/*     //        vec(P0_0m1) = (I-kron(F(:,:,1),F(:,:,1)))\vec(V(:,:,1))   ansi-c*/
/*     //   Note that all eigenvalues of the matrix F(:,:,1) are inside the unit circle when the state-space model is bounded (stationary).   ansi-c*/
/*     //   ansi-c*/
/*     //   April 2008, written by Tao Zha   ansi-c*/
/*     //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.   ansi-c*/

/*     //--- Output arguments.   ansi-c*/
   double loglh_timet;   /*  log likelihood at time t.   ansi-c*/
   TSdmatrix *zt_tm1_dm = kalfiltv_ps->zt_tm1_dm;   /*  nz-by-T.   ansi-c*/
   TSdcell *Pt_tm1_dc = kalfiltv_ps->Pt_tm1_dc;    /*  nz-by-nz-T.   ansi-c*/
/*     //--- Input arguments.   ansi-c*/
   int tdata, tp1;
   TSdvector *z0_dv = kalfiltv_ps->z0_dv;   /*  nz-by-1;   ansi-c*/
   TSdmatrix *P0_dm = kalfiltv_ps->P0_dm;    /*  nz-by-nz.   ansi-c*/
   int T = kalfiltv_ps->T;
   int ny = kalfiltv_ps->ny;
   int nz = kalfiltv_ps->nz;
/*     //--- Work arguments.   ansi-c*/
   int nz2 = square(nz);
   TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *Wnz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdmatrix *W2nz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdvector *wP0_dv = CreateVector_lf(nz2);
/*     //+   ansi-c*/
   TSdvector yt_sdv, at_sdv, zt_tm1_sdv, ztp1_t_sdv, btp1_sdv;
   TSdvector *wny_dv = CreateVector_lf(ny);
   TSdmatrix *Wnzbyny_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *W2nzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *PHtran_tdata_dm = CreateMatrix_lf(nz,ny);
   TSdvector *etdata_dv = CreateVector_lf(ny);
   TSdmatrix *Dtdata_dm = CreateMatrix_lf(ny,ny);
   TSdmatrix *Kt_tdata0_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *Kt_tdata_dm = CreateMatrix_lf(nz,ny);
/*     //--- For eigenvalue decompositions   ansi-c*/
   int ki;
   int errflag;
   double eigmax, logdet_Dtdata;
   TSdzvector *evals_dzv = NULL;
   TSdvector *evals_abs_dv = NULL;   /*  Absolute eigenvalues.   ansi-c*/
/*     //--- Input arguments.   ansi-c*/
   TSdmatrix *yt_dm = kalfiltv_ps->yt_dm;    /*  ny-by-T.   ansi-c*/
   TSdmatrix *at_dm = kalfiltv_ps->at_dm;    /*  ny-by-T.   ansi-c*/
   TSdcell *Ht_dc = kalfiltv_ps->Ht_dc;    /*  ny-by-nz-by-T.   ansi-c*/
   TSdcell *Rt_dc = kalfiltv_ps->Rt_dc;    /*  ny-by-ny-by-T.  Covariance matrix for the measurement equation.   ansi-c*/
   TSdcell *Gt_dc = kalfiltv_ps->Gt_dc;    /*  nz-by-ny-by-T.  Cross-covariance.   ansi-c*/
/*     //   ansi-c*/
   TSdmatrix *bt_dm = kalfiltv_ps->bt_dm;    /*  nz-by-T.   ansi-c*/
   TSdcell *Ft_dc = kalfiltv_ps->Ft_dc;    /*  nz-by-nz-by-T.   ansi-c*/
   TSdcell *Vt_dc = kalfiltv_ps->Vt_dc;    /*  nz-by-nz-by-T.  Covariance matrix for the state equation.   ansi-c*/
/*     //   ansi-c*/


   tdata = (tp1=inpt) - 1;  /*  Base-0 time.   ansi-c*/

/*     //======= Initial condition. =======   ansi-c*/
   if (tdata==0)
   {
/*        //=== Initializing.   ansi-c*/
      if (!kalfiltv_ps->indxIni)
      {
         InitializeDiagonalMatrix_lf(Wnzbynz_dm, 1.0);   /*  To be used for I(nz) -   ansi-c*/
         InitializeDiagonalMatrix_lf(Wnz2bynz2_dm, 1.0);   /*  To be used for I(nz2) -   ansi-c*/

/*           //=== Eigenanalysis to determine the roots to ensure boundedness.   ansi-c*/
         evals_dzv = CreateVector_dz(nz);
         evals_abs_dv = CreateVector_lf(nz);
         errflag = eigrgen(evals_dzv, (TSdzmatrix *)NULL, (TSdzmatrix *)NULL, Ft_dc->C[0]);
         if (errflag)  fn_DisplayError("tz_logTimetCondLH_kalfiltv() in kalman.c: eigen decomposition failed");
         for (ki=nz-1; ki>=0; ki--)  evals_abs_dv->v[ki] = sqrt(square(evals_dzv->real->v[ki]) + square(evals_dzv->imag->v[ki]));
         evals_abs_dv->flag = V_DEF;
         eigmax = MaxVector(evals_abs_dv);
         if (eigmax < (1.0+1.0e-14))
         {
/*              //--- Getting z0_dv: zt_tm1(:,1) = (eye(n_z)-F(:,:,1))\b(:,1);   ansi-c*/
            MatrixMinusMatrix(Wnzbynz_dm, Wnzbynz_dm, Ft_dc->C[0]);
            CopySubmatrix2vector(z0_dv, 0, bt_dm, 0, 0, bt_dm->nrows);
            bdivA_rgens(z0_dv, z0_dv, '\\', Wnzbynz_dm);
/*                        //Done with Wnzbynz_dm.   ansi-c*/
/*              //--- Getting P0_dm: Pt_tm1(:,:,1) = reshape((eye(n_z^2)-kron(F(:,:,1),F(:,:,1)))\V1(:),n_z,n_z);   ansi-c*/
            tz_kron(W2nz2bynz2_dm, Ft_dc->C[0], Ft_dc->C[0]);
            MatrixMinusMatrix(Wnz2bynz2_dm, Wnz2bynz2_dm, W2nz2bynz2_dm);
            CopySubmatrix2vector(wP0_dv, 0, Vt_dc->C[0], 0, 0, nz2);
            bdivA_rgens(wP0_dv, wP0_dv, '\\', Wnz2bynz2_dm);
            CopySubvector2matrix_unr(P0_dm, 0, 0, wP0_dv, 0, nz2);
/*                         //Done with all w*_dv and W*_dm.   ansi-c*/
         }
         else
         {
            fprintf(FPTR_DEBUG, "Fatal error: tz_logTimetCondLH_kalfiltv() in kalman.c: the system is non-stationary solutions\n"
                                "   and thus the initial conditions must be supplied by, say, input arguments");
            fflush(FPTR_DEBUG);
            exit( EXIT_FAILURE );
        }
      }
      CopySubvector2matrix(zt_tm1_dm, 0, 0, z0_dv, 0, z0_dv->n);
      CopyMatrix0(Pt_tm1_dc->C[tdata], P0_dm);
   }


/*     //======= Liklihood at time t (see p.002 in LiuWZ). =======   ansi-c*/
   at_sdv.n = yt_sdv.n = yt_dm->nrows;
   at_sdv.flag = yt_sdv.flag = V_DEF;
   zt_tm1_sdv.n = ztp1_t_sdv.n = zt_tm1_dm->nrows;
   zt_tm1_sdv.flag = ztp1_t_sdv.flag = V_DEF;
   btp1_sdv.n = bt_dm->nrows;
   btp1_sdv.flag = V_DEF;

/*     //--- Setup.   ansi-c*/
   MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_dc->C[tdata], Ht_dc->C[tdata], 1.0, 0.0, 'N', 'T');

/*     //--- Data.   ansi-c*/
/*     //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata;   ansi-c*/
   yt_sdv.v = yt_dm->M + tdata*yt_dm->nrows;
   at_sdv.v = at_dm->M + tdata*at_dm->nrows;
   zt_tm1_sdv.v = zt_tm1_dm->M + tdata*zt_tm1_dm->nrows;
   VectorMinusVector(etdata_dv, &yt_sdv, &at_sdv);
   MatrixTimesVector(etdata_dv, Ht_dc->C[tdata], &zt_tm1_sdv, -1.0, 1.0, 'N');
/*     //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tdata);   ansi-c*/
   CopyMatrix0(Dtdata_dm, Rt_dc->C[tdata]);
   MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[tdata], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
   ScalarTimesMatrixSquare(Dtdata_dm, 0.5, Dtdata_dm, 'T', 0.5);   /*  Making it symmetric against some rounding errors.   ansi-c*/
/*                        //This making-symmetric is very IMPORTANT; otherwise, we will get the matrix being singular message   ansi-c*/
/*                        //    and eigenvalues being negative for the SPD matrix, etc.  Then the likelihood becomes either   ansi-c*/
/*                        //    a bad number or a complex number.   ansi-c*/
   Dtdata_dm->flag = Dtdata_dm->flag | M_SU | M_SL;

/*     //--- Forming the log likelihood.   ansi-c*/
   if (!isfinite(logdet_Dtdata=logdetspd(Dtdata_dm)))  return (loglh_timet = -NEARINFINITY);
   bdivA_rgens(wny_dv, etdata_dv, '/', Dtdata_dm);
   loglh_timet = -(0.5*ny)*LOG2PI - 0.5*logdet_Dtdata - 0.5*VectorDotVector(wny_dv, etdata_dv);
/*                           //Done with all w*_dv.   ansi-c*/


/*     //======= Updating for the next period. =======   ansi-c*/
/*     //--- Updating zt_tm1_dm and Pt_tm1_dc by ztp1_t_sdv and Pt_tm1_dc->C[ti].   ansi-c*/
   if (tp1<T)
   {
/*        //Updating only up to tdata=T-2, because the values at tp1=T or tdata=T-1 will NOT be used in the likelihood function.   ansi-c*/

/*        //- Kt_tdata = (Ft*PHtran_tdata+G(:,:,tdata))/Dtdata;   ansi-c*/
      CopyMatrix0(Kt_tdata0_dm, Gt_dc->C[tdata]);
      MatrixTimesMatrix(Kt_tdata0_dm, Ft_dc->C[tp1], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
      BdivA_rrect(Kt_tdata_dm, Kt_tdata0_dm, '/', Dtdata_dm);
/*        //+ zt_tm1(:,t) = b(:,t) + Ft*zt_tm1(:,tdata) + Kt_tdata*etdata;   ansi-c*/
      ztp1_t_sdv.v = zt_tm1_dm->M + tp1*zt_tm1_dm->nrows;
      MatrixTimesVector(&ztp1_t_sdv, Ft_dc->C[tp1], &zt_tm1_sdv, 1.0, 0.0, 'N');
      MatrixTimesVector(&ztp1_t_sdv, Kt_tdata_dm, etdata_dv, 1.0, 1.0, 'N');
      btp1_sdv.v = bt_dm->M + tp1*btp1_sdv.n;
      VectorPlusMinusVectorUpdate(&ztp1_t_sdv, &btp1_sdv, 1.0);
/*        //+ Pt_tm1(:,:,t) = Ft*Ptdata*Fttran - Kt_tdata*Dtdata*Kt_tdatatran + V(:,:,t);   ansi-c*/
      CopyMatrix0(Pt_tm1_dc->C[tp1], Vt_dc->C[tp1]);
      MatrixTimesMatrix(Wnzbyny_dm, Kt_tdata_dm, Dtdata_dm, 1.0, 0.0, 'N', 'N');
      MatrixTimesMatrix(Wnzbynz_dm, Wnzbyny_dm, Kt_tdata_dm, 1.0, 0.0, 'N', 'T');
      MatrixPlusMinusMatrixUpdate(Pt_tm1_dc->C[tp1], Wnzbynz_dm, -1.0);
/*                              //Done with all W*_dm.   ansi-c*/
      MatrixTimesMatrix(Wnzbynz_dm, Ft_dc->C[tp1], Pt_tm1_dc->C[tdata], 1.0, 0.0, 'N', 'N');
      MatrixTimesMatrix(W2nzbynz_dm, Wnzbynz_dm, Ft_dc->C[tp1], 1.0, 0.0, 'N', 'T');
      MatrixPlusMatrixUpdate(Pt_tm1_dc->C[tp1], W2nzbynz_dm);
/*                              //Done with all W*_dm.   ansi-c*/
   }
   zt_tm1_dm->flag = M_GE;

/*     //===   ansi-c*/
   DestroyVector_dz(evals_dzv);
   DestroyVector_lf(evals_abs_dv);
   DestroyMatrix_lf(Wnzbynz_dm);
   DestroyMatrix_lf(Wnz2bynz2_dm);
   DestroyMatrix_lf(W2nz2bynz2_dm);
   DestroyVector_lf(wP0_dv);
/*     //   ansi-c*/
   DestroyVector_lf(wny_dv);
   DestroyMatrix_lf(Wnzbyny_dm);
   DestroyMatrix_lf(W2nzbynz_dm);
   DestroyMatrix_lf(PHtran_tdata_dm);
   DestroyVector_lf(etdata_dv);
   DestroyMatrix_lf(Dtdata_dm);
   DestroyMatrix_lf(Kt_tdata0_dm);
   DestroyMatrix_lf(Kt_tdata_dm);

   return (loglh_timet);
}




/*  //-----------------------------------------------------   ansi-c*/
/*  //- WARNING: bedore using this function, make sure to call the following functions   ansi-c*/
/*  //      Only once in creating lwzmodel_ps: Refresh_kalfilms_*(lwzmodel_ps);   ansi-c*/
/*  //      Everytime when parameters are changed: RefreshEverything(); RefreRunningGensys_allcases(lwzmodel_ps) in particular.   ansi-c*/
/*  //-----------------------------------------------------   ansi-c*/
double logTimetCondLH_kalfilms_1stapp(int st, int inpt, struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, struct TStateModel_tag *smodel_ps)
{
/*     //st: base-0 grand regime -- deals with the cross-section values at time t.   ansi-c*/
/*     //inpt: base-1 in the sense that inpt>=1 to deal with the time series situation where S_T is (T+1)-by-1 and Y_T is T+nlags_max-by-1.   ansi-c*/
/*     //      The 1st element for S_T is S_T[1] while S_T[0] is s_0 (initial condition).   ansi-c*/
/*     //      The 1st element for Y_T, however, is Y_T[nlags_max+1-1].   ansi-c*/
/*     //See (42.3) on p.42 in the SWZII NOTES.   ansi-c*/

/*     //-- Output arguments   ansi-c*/
   double loglh_timet;
/*     //--- Input arguments   ansi-c*/
   TSdcell *etdata_dc = kalfilmsinputs_1stapp_ps->etdata_dc;  /*  ny-by-nst-by-T, save for computing the likelihood.   ansi-c*/
   TSdfourth *Dtdata_d4 = kalfilmsinputs_1stapp_ps->Dtdata_d4;  /*  ny-by-ny-nst-by-T, save for computing the likelihood and updating Kalman filter Updatekalfilms_1stapp().   ansi-c*/
/*     //--- Local variables   ansi-c*/
   int tbase0;
   double logdet_Dtdata;
/*     //--- Accessible variables   ansi-c*/
   int ny = kalfilmsinputs_1stapp_ps->ny;
   TSdvector etdata_sdv;
/*     //=== Work arguments.   ansi-c*/
   TSdvector *wny_dv = CreateVector_lf(ny);



/*     //--- Critical checking.   ansi-c*/
   if (inpt > kalfilmsinputs_1stapp_ps->T)
      fn_DisplayError(".../kalman.c/logTimetCondLH_kalfilms_1stapp(): The time exceeds the\n"
                      "     data sample size allocated the structure TSkalfilmsinputs_1stapp_tag");

/*     //--- The following is for safe guard.  InitializeKalman_z10_P10() should be called in, say, RefreshEverything().   ansi-c*/
   if (kalfilmsinputs_1stapp_ps->ztm1_track < 0)
      if (!InitializeKalman_z10_P10(kalfilmsinputs_1stapp_ps, (TSdmatrix *)NULL, (TSdcell *)NULL))
         fn_DisplayError(".../kalman.c/logTimetCondLH_kalfilms_1stapp(): the system is non-stationary when calling"
                         "     InitializeKalman_z10_P10().  Please call this function in RefreshEverthing() and"
                         "     set the likehood to be -infty for early exit");

   tbase0=inpt-1;

/*     //-------------------  The order matters. Updatekalfilms_1stapp() must be called before Update_et_Dt_1stapp(). -----------------   ansi-c*/
/*     //--- $$$ Critical updating where we MUSt have inpt-1.  If inpt, Updatekalfilms_1stapp() will call this function again   ansi-c*/
/*     //--- $$$   because DW function ProbabilityStateConditionalCurrent() need to access this function at time inpt,   ansi-c*/
/*     //--- $$$   which has not computed before Updatekalfilms_1stapp().  Thus, we'll have an infinite loop.   ansi-c*/
   Updatekalfilms_1stapp(tbase0, kalfilmsinputs_1stapp_ps, smodel_ps);
/*  //   //--- $$$ Critical updating.   ansi-c*/
/*  //   Update_et_Dt_1stapp(tbase0, kalfilmsinputs_1stapp_ps);   ansi-c*/
/*  //             //This function will give Dtdata_d4->F[tbase0], etdata_dc->C[tbase0], and PHtran_tdata_d4->F[tbase0].   ansi-c*/



/*     //======================================================   ansi-c*/
/*     //= Getting the logLH at time tbase0 or time inpt.   ansi-c*/
/*     //======================================================   ansi-c*/
/*     //--- Forming the log conditional likelihood at t.   ansi-c*/
   etdata_sdv.n = ny;
   etdata_sdv.v = etdata_dc->C[tbase0]->M + ny*st;
   etdata_sdv.flag = V_DEF;
   if (!isfinite(logdet_Dtdata=logdetspd(Dtdata_d4->F[tbase0]->C[st])))  return (loglh_timet = -NEARINFINITY);
   bdivA_rgens(wny_dv, &etdata_sdv, '/', Dtdata_d4->F[tbase0]->C[st]);
   loglh_timet = -(0.5*ny)*LOG2PI - 0.5*logdet_Dtdata - 0.5*VectorDotVector(wny_dv, &etdata_sdv);
/*                           //Done with all w*_dv.   ansi-c*/

/*     //===   ansi-c*/
   DestroyVector_lf(wny_dv);

   return (loglh_timet);
}
/*  //======================================================   ansi-c*/
/*  //= Computing z_{1|0} and P_{1|0} for each new parameter values.   ansi-c*/
/*  //======================================================   ansi-c*/
int InitializeKalman_z10_P10(struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, TSdmatrix *z10_dm, TSdcell *P10_dc)
{
/*     //See p.001 and p.004 in LWZ Model II.   ansi-c*/
/*     //Outputs:   ansi-c*/
/*     //   return 1: success in initializing; 0: initializing fails, so the likelihood must be set to -infty outside this function.   ansi-c*/
/*     //   ztm1_track to track the time up to which Kalman filter have been updated.   ansi-c*/
/*     //   z0_dm, zt_tm1_dc->C[0]   ansi-c*/
/*     //   P0_dc, Pt_tm1_d4->F[0]   ansi-c*/

/*     //--- Output arguments   ansi-c*/
   TSdmatrix *z0_0_dm = kalfilmsinputs_1stapp_ps->z0_dm;         /*  nz-by-nst.   ansi-c*/
   TSdmatrix *z0_dm = kalfilmsinputs_1stapp_ps->z0_dm;         /*  nz-by-nst.   ansi-c*/
   TSdcell *P0_dc = kalfilmsinputs_1stapp_ps->P0_dc;           /*  nz-by-nz-by-nst.   ansi-c*/
/*     //+ Used to get zt_tm1_dc->C[0] and Pt_tm1_d4->F[0] only.   ansi-c*/
   TSdcell *zt_tm1_dc = kalfilmsinputs_1stapp_ps->zt_tm1_dc;  /*  nz-by-nst-by-(T+1).   ansi-c*/
   TSdfourth *Pt_tm1_d4 = kalfilmsinputs_1stapp_ps->Pt_tm1_d4;      /*  nz-by-nz-by-nst-by-(T+1).   ansi-c*/
/*     //--- Input arguments   ansi-c*/
   TSdmatrix *yt_dm = kalfilmsinputs_1stapp_ps->yt_dm;          /*  ny-by-T.   ansi-c*/
   TSdmatrix *at_dm = kalfilmsinputs_1stapp_ps->at_dm;          /*  ny-by-nst.   ansi-c*/
   TSdcell *Ht_dc = kalfilmsinputs_1stapp_ps->Ht_dc;            /*  ny-by-nz-by-nst.   ansi-c*/
   TSdcell *Rt_dc = kalfilmsinputs_1stapp_ps->Rt_dc;            /*  ny-by-ny-by-nst.  Covariance matrix for the measurement equation.   ansi-c*/
/*     //+   ansi-c*/
   TSdmatrix *bt_dm = kalfilmsinputs_1stapp_ps->bt_dm;          /*  nz-by-nst.   ansi-c*/
   TSdcell *Ft_dc = kalfilmsinputs_1stapp_ps->Ft_dc;            /*  nz-by-nz-by-nst.   ansi-c*/
   TSdcell *Vt_dc = kalfilmsinputs_1stapp_ps->Vt_dc;            /*  nz-by-nz-by-nst.  Covariance matrix for the state equation.   ansi-c*/
/*     //--- Local variables   ansi-c*/
   int sti;
/*     //--- Accessible variables   ansi-c*/
   int ny = kalfilmsinputs_1stapp_ps->ny;
   int nz = kalfilmsinputs_1stapp_ps->nz;
   int nst = kalfilmsinputs_1stapp_ps->nst;
   TSdvector z0_sdv, z0_0_sdv, bt_sdv;
   TSdvector yt_sdv, at_sdv;
/*     //--- For the initial conditions: eigenvalue decompositions   ansi-c*/
   int ki;
   int errflag;
   double eigmax;
/*     //===   ansi-c*/
   int nz2 = square(nz);
   TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *Wnz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdmatrix *W2nz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdvector *wP0_dv = CreateVector_lf(nz2);
/*     //   ansi-c*/
   TSdzvector *evals_dzv = evals_dzv = CreateVector_dz(nz);
   TSdvector *evals_abs_dv = CreateVector_lf(nz);  /*  Absolute eigenvalues.   ansi-c*/


   if (kalfilmsinputs_1stapp_ps->ztm1_track < 0)
   {
      z0_sdv.n = z0_0_sdv.n = bt_sdv.n = nz;
      z0_sdv.flag = z0_0_sdv.flag = bt_sdv.flag = V_DEF;
      at_sdv.n = yt_sdv.n = ny;
      at_sdv.flag = yt_sdv.flag = V_DEF;


/*        //======= Initial condition. =======   ansi-c*/
      if (!kalfilmsinputs_1stapp_ps->indxIni)
      {
         z0_0_dm->flag = z0_dm->flag = M_GE;
         for (sti=nst-1; sti>=0;  sti--)
         {
            if (kalfilmsinputs_1stapp_ps->DiffuseScale)  /*  Diffuse initial conditions are used.   ansi-c*/
            {
/*                 //--- Diffuse condition for z0_dv.   ansi-c*/
               z0_sdv.v = z0_dm->M + z0_sdv.n*sti;
               z0_0_sdv.v = z0_0_dm->M + z0_0_sdv.n*sti;
               bt_sdv.v = bt_dm->M + bt_sdv.n*sti;
               InitializeConstantVector_lf(&z0_0_sdv, 0.0);
               MatrixTimesVector(&z0_sdv, Ft_dc->C[sti], &z0_0_sdv, 1.0, 0.0, 'N');
               VectorPlusVector(&z0_sdv, &z0_sdv, &bt_sdv);
/*                 //--- Diffuse condition for P0_dm.   ansi-c*/
               InitializeDiagonalMatrix_lf(Wnzbynz_dm, kalfilmsinputs_1stapp_ps->DiffuseScale);   /*  To be used for DiffuseScale*I(nz)   ansi-c*/
               CopyMatrix0(P0_dc->C[sti], Wnzbynz_dm);
/*                             //Done with W*_dm.   ansi-c*/
            }
            else //Unconditional moments for initial conditions are used.
            {
               InitializeDiagonalMatrix_lf(Wnzbynz_dm, 1.0);   /*  To be used for I(nz) -   ansi-c*/
               InitializeDiagonalMatrix_lf(Wnz2bynz2_dm, 1.0);   /*  To be used for I(nz2) -   ansi-c*/

/*                 //=== Eigenanalysis to determine the roots to ensure boundedness.   ansi-c*/
               errflag = eigrgen(evals_dzv, (TSdzmatrix *)NULL, (TSdzmatrix *)NULL, Ft_dc->C[sti]);
               if (errflag)  fn_DisplayError("kalman.c/InitializeKalman_z10_P10(): eigen decomposition failed");
               for (ki=nz-1; ki>=0; ki--)  evals_abs_dv->v[ki] = sqrt(square(evals_dzv->real->v[ki]) + square(evals_dzv->imag->v[ki]));
               evals_abs_dv->flag = V_DEF;
               eigmax = MaxVector(evals_abs_dv);
               if (eigmax < (1.0-SQRTEPSILON))  /*  (1.0+EPSILON))   ansi-c*/
               {
/*                    //--- Getting z0_dv: zt_tm1(:,1) = (eye(n_z)-F(:,:,sti))\b(:,sti);   ansi-c*/
                  z0_0_sdv.v = z0_0_dm->M + z0_0_sdv.n*sti;
                  z0_sdv.v = z0_dm->M + z0_sdv.n*sti;
                  MatrixMinusMatrix(Wnzbynz_dm, Wnzbynz_dm, Ft_dc->C[sti]);
                  CopySubmatrix2vector(&z0_0_sdv, 0, bt_dm, 0, sti, bt_dm->nrows);
                  bdivA_rgens(&z0_0_sdv, &z0_0_sdv, '\\', Wnzbynz_dm);
/*                    //- Under the assumption s_0 = s_1 (this is a short-cut).   ansi-c*/
                  MatrixTimesVector(&z0_sdv, Ft_dc->C[sti], &z0_0_sdv, 1.0, 0.0, 'N');
                  VectorPlusVector(&z0_sdv, &z0_sdv, &bt_sdv);
/*                              //Done with Wnzbynz_dm.   ansi-c*/
/*                    //--- Getting P0_dm: Pt_tm1(:,:,1) = reshape((eye(n_z^2)-kron(F(:,:,sti),F(:,:,sti)))\V1(:),n_z,n_z);   ansi-c*/
                  tz_kron(W2nz2bynz2_dm, Ft_dc->C[sti], Ft_dc->C[sti]);
                  MatrixMinusMatrix(Wnz2bynz2_dm, Wnz2bynz2_dm, W2nz2bynz2_dm);
                  CopySubmatrix2vector(wP0_dv, 0, Vt_dc->C[sti], 0, 0, nz2);
                  bdivA_rgens(wP0_dv, wP0_dv, '\\', Wnz2bynz2_dm);
                  CopySubvector2matrix_unr(P0_dc->C[sti], 0, 0, wP0_dv, 0, nz2);
/*                               //Done with all w*_dv and W*_dm.   ansi-c*/
               }
               else
               {
                  if (0)  /*  0: no printing.   ansi-c*/
                  {
                     #if defined (USE_DEBUG_FILE)
                     fprintf(FPTR_DEBUG, "\n-------WARNING: ----------\n");
                     fprintf(FPTR_DEBUG, "\nIn grand regime sti=%d\n", sti);
                     fprintf(FPTR_DEBUG, ".../kalman.c/InitializeKalman_z10_P10(): the system is non-stationary solutions\n"
                                         "    and see p.003 in LWZ Model II");
                     #else
                     printf("\n-----------------\n");
                     printf("\nIn grand regime sti=%d\n", sti);
                     printf(".../kalman.c/InitializeKalman_z10_P10(): the system is non-stationary solutions\n"
                                     "    and see p.003 in LWZ Model II");
                     #endif
                  }
/*                    //=== See p.000.3 in LWZ Model II.   ansi-c*/
/*                    //=== Do NOT use the following option.  It turns out that this will often generate explosive conditional liklihood   ansi-c*/
/*                    //===   at the end of the sample, because Pt_tm1 shrinks to zero overtime due to the sigularity of   ansi-c*/
/*                    //===   the initila condition P_{1|0}.   ansi-c*/
/*                    //--- Letting z0_dv = 0.0   ansi-c*/
/*                    // z0_sdv.v = z0_dm->M + z0_sdv.n*sti;   ansi-c*/
/*                    // InitializeConstantVector_lf(&z0_sdv, 0.0);   ansi-c*/
/*                    // //--- Letting P0_dm = V   ansi-c*/
/*                    // CopyMatrix0(P0_dc->C[sti], Vt_dc->C[sti]);   ansi-c*/

/*                    //===   ansi-c*/
                  DestroyVector_dz(evals_dzv);
                  DestroyVector_lf(evals_abs_dv);
                  DestroyMatrix_lf(Wnzbynz_dm);
                  DestroyMatrix_lf(Wnz2bynz2_dm);
                  DestroyMatrix_lf(W2nz2bynz2_dm);
                  DestroyVector_lf(wP0_dv);

                  return (0);   /*  Early exit with kalfilmsinputs_1stapp_ps->ztm1_track continues to be -1.   ansi-c*/
               }
            }
         }
      }
      else
      {
         if (!z10_dm)  fn_DisplayError(".../kalman.c/InitializeKalman_z10_P10(): The initial condition z_{1|0}\n"
                                       "     must be supplied as valid input arguments for when indxIni == 1");
         else
            CopyMatrix0(z0_dm, z10_dm);

         if (!P10_dc)  fn_DisplayError(".../kalman.c/InitializeKalman_z10_P10(): The initial condition P_{1|0}\n"
                                       "     must be supplied as valid input arguments for when indxIni == 1");
         else
            CopyCell0(P0_dc, P10_dc);
      }
      CopyMatrix0(zt_tm1_dc->C[0], z0_dm);  /*  At time t-1 = 1.   ansi-c*/
      CopyCell0(Pt_tm1_d4->F[0], P0_dc);  /*  At time t-1 = 1.   ansi-c*/


      kalfilmsinputs_1stapp_ps->ztm1_track = 0;   /*  Must reset to 0, meaning initial setting is done and ready for computing LH at t = 1.   ansi-c*/

      Update_et_Dt_1stapp(0, kalfilmsinputs_1stapp_ps);

/*        //===   ansi-c*/
      DestroyVector_dz(evals_dzv);
      DestroyVector_lf(evals_abs_dv);
      DestroyMatrix_lf(Wnzbynz_dm);
      DestroyMatrix_lf(Wnz2bynz2_dm);
      DestroyMatrix_lf(W2nz2bynz2_dm);
      DestroyVector_lf(wP0_dv);

      return (1);
   }
   else
   {
      fn_DisplayError(".../kalman.c/InitializeKalman_z10_P10(): calling this function makes sense only if"
                         "     kalfilmsinputs_1stapp_ps->ztm1_track is -1.  Please check this value.");

/*        //===   ansi-c*/
      DestroyVector_dz(evals_dzv);
      DestroyVector_lf(evals_abs_dv);
      DestroyMatrix_lf(Wnzbynz_dm);
      DestroyMatrix_lf(Wnz2bynz2_dm);
      DestroyMatrix_lf(W2nz2bynz2_dm);
      DestroyVector_lf(wP0_dv);

      return (0);
   }
}
/*  //======================================================   ansi-c*/
/*  //= Integrating out the lagged regimes in order to   ansi-c*/
/*  //=   updating zt_tm1 and Pt_tm1 for next perid tp1 through Kim-Nelson filter.   ansi-c*/
/*  //= tdata representing base-0 t timing, while inpt represents base-1 t timing.   ansi-c*/
/*  //   ansi-c*/
/*  //= Purpose: for each inpt, we integrate out grand regimes st   ansi-c*/
/*  //=   only ONCE to prevent the dimension of updated zt_tm1 and Pt_tm1 through Kim-Nelson filter.   ansi-c*/
/*  //======================================================   ansi-c*/
static int Updatekalfilms_1stapp(int t_1, struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, struct TStateModel_tag *smodel_ps)
{
/*     //Output:   ansi-c*/
/*     //  tm1update   ansi-c*/
/*     //  z_{t_1+1|t_1}   ansi-c*/
/*     //  P_{t_1+1|t_1}   ansi-c*/
/*     //Input:   ansi-c*/
/*     //  t-1: base-1 t timing.  Thus t-1=inpt-1.   ansi-c*/

/*     //--- Local variables   ansi-c*/
   int stp1i, sti, t_2, t_2p1;
   double prob_previous_regimes;
/*     //-- Output arguments   ansi-c*/
   TSdcell *zt_tm1_dc = kalfilmsinputs_1stapp_ps->zt_tm1_dc;  /*  nz-by-nst-by-(T+1).   ansi-c*/
   TSdfourth *Pt_tm1_d4 = kalfilmsinputs_1stapp_ps->Pt_tm1_d4;      /*  nz-by-nz-by-nst-by-(T+1).   ansi-c*/
/*     //--- Input arguments   ansi-c*/
   TSdcell *Gt_dc = kalfilmsinputs_1stapp_ps->Gt_dc;            /*  nz-by-ny-by-nst.  Cross-covariance.   ansi-c*/
/*     //+   ansi-c*/
   TSdmatrix *bt_dm = kalfilmsinputs_1stapp_ps->bt_dm;          /*  nz-by-nst.   ansi-c*/
   TSdcell *Ft_dc = kalfilmsinputs_1stapp_ps->Ft_dc;            /*  nz-by-nz-by-nst.   ansi-c*/
   TSdcell *Vt_dc = kalfilmsinputs_1stapp_ps->Vt_dc;            /*  nz-by-nz-by-nst.  Covariance matrix for the state equation.   ansi-c*/
/*     //+   ansi-c*/
   TSdfourth *PHtran_tdata_d4 = kalfilmsinputs_1stapp_ps->PHtran_tdata_d4;   /*  nz-by-ny-by-nst-T, saved only for updating Kalman filter Updatekalfilms_1stapp().   ansi-c*/
   TSdcell *etdata_dc = kalfilmsinputs_1stapp_ps->etdata_dc;  /*  ny-by-nst-by-T, save for computing the likelihood.   ansi-c*/
   TSdfourth *Dtdata_d4 = kalfilmsinputs_1stapp_ps->Dtdata_d4;  /*  ny-by-ny-nst-by-T, save for computing the likelihood and updating Kalman filter Updatekalfilms_1stapp().   ansi-c*/
/*     //--- Accessible variables   ansi-c*/
   int ny = kalfilmsinputs_1stapp_ps->ny;
   int nz = kalfilmsinputs_1stapp_ps->nz;
   int nst = kalfilmsinputs_1stapp_ps->nst;
   int T = kalfilmsinputs_1stapp_ps->T;
   TSdvector z0_sdv;
   TSdvector btp1_sdv;
   TSdvector etdata_sdv;
/*     //=== Work arguments.   ansi-c*/
   TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz,nz);
/*     //+   ansi-c*/
   TSdmatrix *Wnzbyny_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *W2nzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *Kt_tdata0_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *Kt_tdata_dm = CreateMatrix_lf(nz,ny);
/*     //=== For updating zt_tm1_dm and Pt_tm1.   ansi-c*/
   TSdvector *ztp1_t_dv = CreateVector_lf(nz);
   TSdmatrix *Ptp1_t_dm = CreateMatrix_lf(nz, nz);
   TSdvector *ztp1_dv = CreateVector_lf(nz);
   TSdmatrix *Ptp1_dm = CreateMatrix_lf(nz, nz);


/*     //--- Critical checking.   ansi-c*/
   if (kalfilmsinputs_1stapp_ps->ztm1_track < 0)
      fn_DisplayError(".../kalman.c/Updatekalfilms_1stapp(): Make sure InitializeKalman_z10_P10() is called in the function RefreshEverthing()");


   z0_sdv.n = nz;
   z0_sdv.flag = V_DEF;
   btp1_sdv.n = nz;
   btp1_sdv.flag = V_DEF;
/*     //+   ansi-c*/
   etdata_sdv.n = ny;
   etdata_sdv.flag = V_DEF;

   for (t_2=kalfilmsinputs_1stapp_ps->ztm1_track; t_2<t_1; t_2++)
   {
/*        //If t_1 <= ztm1_track, no updating.   ansi-c*/
/*        //If t_1 > ztm1_track, updating z_{t|t-1} and P_{t|t-1} up to t-1 = t_1.   ansi-c*/

      zt_tm1_dc->C[t_2p1=t_2+1]->flag = M_GE;
      for (stp1i=nst-1; stp1i>=0;  stp1i--)
      {
         InitializeConstantVector_lf(ztp1_dv, 0.0);   /*  To be summed over sti.   ansi-c*/
         InitializeConstantMatrix_lf(Ptp1_dm, 0.0);   /*  To be summed over sti.   ansi-c*/

         for (sti=nst-1; sti>=0;  sti--)
         {
/*              //=== Updating for next period by integrating out sti..   ansi-c*/
/*              //--- Ktp1_t = (F_tp1*PHtran_t+G(:,:,t))/Dt;   ansi-c*/
/*              //--- Kt_tdata = (Ft*PHtran_tdata+G(:,:,tdata))/Dtdata where t=tp1 and tdata=t.   ansi-c*/
            CopyMatrix0(Kt_tdata0_dm, Gt_dc->C[sti]);
            MatrixTimesMatrix(Kt_tdata0_dm, Ft_dc->C[stp1i], PHtran_tdata_d4->F[t_2]->C[sti], 1.0, 1.0, 'N', 'N');
            BdivA_rrect(Kt_tdata_dm, Kt_tdata0_dm, '/', Dtdata_d4->F[t_2]->C[sti]);
/*              //+ zt_tm1(:,t) = b(:,t) + Ft*zt_tm1(:,tdata) + Kt_tdata*etdata where t=tp1 and tm1=t.   ansi-c*/
            etdata_sdv.v = etdata_dc->C[t_2]->M + ny*sti;
            z0_sdv.v = zt_tm1_dc->C[t_2]->M + nz*sti;   /*  sti: regime at time t_2.   ansi-c*/
            MatrixTimesVector(ztp1_t_dv, Ft_dc->C[stp1i], &z0_sdv, 1.0, 0.0, 'N');
            MatrixTimesVector(ztp1_t_dv, Kt_tdata_dm, &etdata_sdv, 1.0, 1.0, 'N');
            btp1_sdv.v = bt_dm->M + stp1i*btp1_sdv.n;
            VectorPlusMinusVectorUpdate(ztp1_t_dv, &btp1_sdv, 1.0);
/*              //+ Pt_tm1(:,:,t) = Ft*Ptdata*Fttran - Kt_tdata*Dtdata*Kt_tdatatran + V(:,:,t);   ansi-c*/
            CopyMatrix0(Ptp1_t_dm, Vt_dc->C[stp1i]);
            MatrixTimesMatrix(Wnzbyny_dm, Kt_tdata_dm, Dtdata_d4->F[t_2]->C[sti], 1.0, 0.0, 'N', 'N');
            MatrixTimesMatrix(Wnzbynz_dm, Wnzbyny_dm, Kt_tdata_dm, 1.0, 0.0, 'N', 'T');
            MatrixPlusMinusMatrixUpdate(Ptp1_t_dm, Wnzbynz_dm, -1.0);
/*                                    //Done with all W*_dm.   ansi-c*/
            MatrixTimesMatrix(Wnzbynz_dm, Ft_dc->C[stp1i], Pt_tm1_d4->F[t_2]->C[sti], 1.0, 0.0, 'N', 'N');
            MatrixTimesMatrix(W2nzbynz_dm, Wnzbynz_dm, Ft_dc->C[stp1i], 1.0, 0.0, 'N', 'T');
            MatrixPlusMatrixUpdate(Ptp1_t_dm, W2nzbynz_dm);
/*                                    //Done with all W*_dm.   ansi-c*/


/*              //--- Integrating out the state at t_2 using   ansi-c*/
/*              //---    P(s_t_2|Y_{t_2}, theta) = ProbabilityStateConditionalCurrent(sti, t_2, smodel_ps);   ansi-c*/
/*              //--- One can also access to P(s_t_2|Y_{t_2}, theta) by using ElementV(smodel_ps->V[t_2],s_{t_2}i),   ansi-c*/
/*              //---    but this access will not call my function logTimetCondLH(), thus no updating for   ansi-c*/
/*              //---    P(s_t_2|Y_{t_2}, and thus leading to incorrect results.   ansi-c*/
            prob_previous_regimes = ProbabilityStateConditionalCurrent(sti, t_2, smodel_ps);
            ScalarTimesVectorUpdate(ztp1_dv, prob_previous_regimes, ztp1_t_dv);
            ScalarTimesMatrix(Ptp1_dm, prob_previous_regimes, Ptp1_t_dm, 1.0);
            Ptp1_dm->flag = M_GE | M_SU | M_SL;
/*                                        //Done with ztp1_t_dv and Ptp1_t_dm.   ansi-c*/
         }
/*           //--- Filling zt_tm1 and Pt_tm1 for next period.   ansi-c*/
         z0_sdv.v = zt_tm1_dc->C[t_2p1]->M + z0_sdv.n*stp1i;   /*  stp1i: regime at time tp1.   ansi-c*/
         CopyVector0(&z0_sdv, ztp1_dv);
         CopyMatrix0(Pt_tm1_d4->F[t_2p1]->C[stp1i], Ptp1_dm);   /*  stp1i: regime at time tp1.   ansi-c*/
/*                                             //Done with ztp1_dv, z0_sdv, Ptp1_dm.   ansi-c*/
      }
/*        //--- $$$ The following is important because it tells ProbabilityStateConditionalCurrent(), which calls   ansi-c*/
/*        //--- $$$   logTimetCondLH_kalfilms_1stapp(), which calls recursively this function again, that there is no   ansi-c*/
/*        //--- $$$   need to update Kalman filter for the period before kalfilmsinputs_1stapp_ps->ztm1_track.   ansi-c*/
      kalfilmsinputs_1stapp_ps->ztm1_track = t_2p1;  /*  Means that z_{t_2p1+1|t_2p1} and P_{t_2p1+1|t_2p1} are done.   ansi-c*/

/*        //--- $$$ This function must be called after all the above computations are done.   ansi-c*/
      Update_et_Dt_1stapp(t_2p1, kalfilmsinputs_1stapp_ps);
   }


/*     //===   ansi-c*/
   DestroyMatrix_lf(Wnzbynz_dm);
/*     //   ansi-c*/
   DestroyMatrix_lf(Wnzbyny_dm);
   DestroyMatrix_lf(W2nzbynz_dm);
   DestroyMatrix_lf(Kt_tdata0_dm);
   DestroyMatrix_lf(Kt_tdata_dm);
/*     //   ansi-c*/
   DestroyVector_lf(ztp1_t_dv);
   DestroyMatrix_lf(Ptp1_t_dm);
   DestroyVector_lf(ztp1_dv);
   DestroyMatrix_lf(Ptp1_dm);

   return (kalfilmsinputs_1stapp_ps->ztm1_track);
}
/*  //======================================================   ansi-c*/
/*  //= Computes etdata and Dtdata for all grand regimes st at tbase0=inpt-1 or dtm1_track   ansi-c*/
/*  //=   to prevent recomputing this object for different st at given tbase0.   ansi-c*/
/*  //======================================================   ansi-c*/
static int Update_et_Dt_1stapp(int t_1, struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps)
{
/*     //Output:   ansi-c*/
/*     //  dtm1_track is updated in this function.   ansi-c*/
/*     //  PHtran_tdata_d4->F[t-1]   ansi-c*/
/*     //  etdata_dc->C[t-1]   ansi-c*/
/*     //  Dtdata_d4->F[t-1]   ansi-c*/
/*     //Input:   ansi-c*/
/*     //  t_1=inpt-1: base-0 timing for et and Dt before the likelihood at time inpt is computed.   ansi-c*/

/*     //--- Local variables   ansi-c*/
   int sti, tbase0;
/*     //-- Output arguments   ansi-c*/
   TSdfourth *PHtran_tdata_d4 = kalfilmsinputs_1stapp_ps->PHtran_tdata_d4;   /*  nz-by-ny-by-nst-T, saved only for updating Kalman filter Updatekalfilms_1stapp().   ansi-c*/
   TSdcell *etdata_dc = kalfilmsinputs_1stapp_ps->etdata_dc;  /*  ny-by-nst-by-T, save for computing the likelihood.   ansi-c*/
   TSdcell *yt_tm1_dc = kalfilmsinputs_1stapp_ps->yt_tm1_dc;  /*  ny-by-nst-by-T, one-step forecast y_{t|t-1} for t=0 to T-1 (base-0).   ansi-c*/
   TSdfourth *Dtdata_d4 = kalfilmsinputs_1stapp_ps->Dtdata_d4;  /*  ny-by-ny-nst-by-T, save for computing the likelihood and updating Kalman filter Updatekalfilms_1stapp().   ansi-c*/
/*     //--- input arguments   ansi-c*/
   TSdcell *zt_tm1_dc = kalfilmsinputs_1stapp_ps->zt_tm1_dc;  /*  nz-by-nst-by-T.   ansi-c*/
   TSdfourth *Pt_tm1_d4 = kalfilmsinputs_1stapp_ps->Pt_tm1_d4;      /*  nz-by-nz-by-nst-by-T.   ansi-c*/
/*     //+   ansi-c*/
   TSdmatrix *yt_dm = kalfilmsinputs_1stapp_ps->yt_dm;          /*  ny-by-T.   ansi-c*/
   TSdmatrix *at_dm = kalfilmsinputs_1stapp_ps->at_dm;          /*  ny-by-nst.   ansi-c*/
   TSdcell *Ht_dc = kalfilmsinputs_1stapp_ps->Ht_dc;            /*  ny-by-nz-by-nst.   ansi-c*/
   TSdcell *Rt_dc = kalfilmsinputs_1stapp_ps->Rt_dc;            /*  ny-by-ny-by-nst.  Covariance matrix for the measurement equation.   ansi-c*/
/*     //--- Accessible variables   ansi-c*/
   int ny = kalfilmsinputs_1stapp_ps->ny;
   int nz = kalfilmsinputs_1stapp_ps->nz;
   int nst = kalfilmsinputs_1stapp_ps->nst;
   TSdvector z0_sdv;
   TSdvector yt_sdv, at_sdv;
   TSdvector etdata_sdv, yt_tm1_sdv;
/*     //=== Work arguments.   ansi-c*/
   TSdmatrix *PHtran_tdata_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *Dtdata_dm = CreateMatrix_lf(ny,ny);


   z0_sdv.n = nz;
   z0_sdv.flag = V_DEF;
   at_sdv.n = yt_sdv.n = ny;
   at_sdv.flag = yt_sdv.flag = V_DEF;
   etdata_sdv.n = yt_tm1_sdv.n = ny;
   etdata_sdv.flag = yt_tm1_sdv.flag = V_DEF;

   for (tbase0=(kalfilmsinputs_1stapp_ps->dtm1_track+1); tbase0<=t_1; tbase0++)
   {
/*        //Note tbase0<=t_1, NOT tbase0<t_1.   ansi-c*/
/*        //If t_1 < (dtm1_track+1), no updating.   ansi-c*/
/*        //If t_1 >= (dtm1_track+1), updating etdata_dc->C[t-1] and Dtdata_d4->F[t-1] up to t-1=t_1.   ansi-c*/

      for (sti=nst-1; sti>=0;  sti--)
      {
/*           //--- Setup.   ansi-c*/
         MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_d4->F[tbase0]->C[sti], Ht_dc->C[sti], 1.0, 0.0, 'N', 'T');
         CopyMatrix0(kalfilmsinputs_1stapp_ps->PHtran_tdata_d4->F[tbase0]->C[sti], PHtran_tdata_dm);


/*           //--- Data.   ansi-c*/
/*           //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata where tdata = tbase0 = inpt-1.   ansi-c*/
         yt_sdv.v = yt_dm->M + tbase0*yt_dm->nrows;
         at_sdv.v = at_dm->M + sti*at_dm->nrows;   /*  grand regime at time tbase0.   ansi-c*/
         z0_sdv.v = zt_tm1_dc->C[tbase0]->M + z0_sdv.n*sti;   /*  sti: regime at time tbase0.   ansi-c*/
         etdata_sdv.v = etdata_dc->C[tbase0]->M + etdata_sdv.n*sti;
         yt_tm1_sdv.v = etdata_dc->C[tbase0]->M + yt_tm1_sdv.n*sti;
         CopyVector0(&yt_tm1_sdv, &at_sdv);
         MatrixTimesVector(&yt_tm1_sdv, Ht_dc->C[sti], &z0_sdv, 1.0, 1.0, 'N');  /*  a + H*z_{t|t-1}.   ansi-c*/
         VectorMinusVector(&etdata_sdv, &yt_sdv, &yt_tm1_sdv);  /*  y_t - a - H*z_{t|t-1}.   ansi-c*/
/*           //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tbase0);   ansi-c*/
         CopyMatrix0(Dtdata_dm, Rt_dc->C[sti]);
         MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[sti], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
/*                                       //Done with z0_sdv.v.   ansi-c*/
         ScalarTimesMatrixSquare(Dtdata_dm, 0.5, Dtdata_dm, 'T', 0.5);   /*  Making it symmetric against some rounding errors.   ansi-c*/
/*                              //This making-symmetric is very IMPORTANT; otherwise, we will get the matrix being singular message   ansi-c*/
/*                              //    and eigenvalues being negative for the SPD matrix, etc.  Then the likelihood becomes either   ansi-c*/
/*                              //    a bad number or a complex number.   ansi-c*/
         Dtdata_dm->flag = Dtdata_dm->flag | M_SU | M_SL;
         CopyMatrix0(Dtdata_d4->F[tbase0]->C[sti], Dtdata_dm);  /*  Saved to be used for logTimetCondLH_kalfilms_1stapp().   ansi-c*/
      }

/*        //--- $$$ This tracker functions the same way as kalfilmsinputs_1stapp_ps->ztm1_track.   ansi-c*/
      kalfilmsinputs_1stapp_ps->dtm1_track = tbase0;
   }

/*     //===   ansi-c*/
   DestroyMatrix_lf(PHtran_tdata_dm);
   DestroyMatrix_lf(Dtdata_dm);

   return (kalfilmsinputs_1stapp_ps->dtm1_track);
}






/*  //-----------------------------------------------------   ansi-c*/
/*  //------------ OLD Code --------------------------   ansi-c*/
/*  //- Updating or refreshing all Kalman filter at time t for Markov-switching DSGE model.   ansi-c*/
/*  //- WARNING: make sure to call the following functions   ansi-c*/
/*  //      RunningGensys_const7varionly(lwzmodel_ps);   ansi-c*/
/*  //      Refresh_kalfilms_*(lwzmodel_ps);   //Creates or refreshes kalfilmsinputs_ps at new parameter values.   ansi-c*/
/*  //- before using tz_Refresh_z_T7P_T_in_kalfilms_1st_approx().   ansi-c*/
/*  //   ansi-c*/
/*  //- IMPORTANT NOTE: in the Markov-switching input file datainp_markov*.prn, it MUST be that   ansi-c*/
/*  //-                                 the coefficient regime is the 1st state variable, and   ansi-c*/
/*  //-                                 the volatility regime is the 2nd state variable.   ansi-c*/
/*  //-----------------------------------------------------   ansi-c*/
#if defined (NEWVERSIONofDW_SWITCH)
double tz_logTimetCondLH_kalfilms_1st_approx(int st, int inpt, struct TSkalfilmsinputs_tag *kalfilmsinputs_ps, struct TStateModel_tag *smodel_ps)
{
/*     //st, st_c, and st_v: base-0: deals with the cross-section values at time t where   ansi-c*/
/*     //      st is a grand regime, st_c is an encoded coefficient regime, and st_c is an encoded volatility regime.   ansi-c*/
/*     //inpt: base-1 in the sense that inpt>=1 to deal with the time series situation where S_T is (T+1)-by-1 and Y_T is T+nlags_max-by-1.   ansi-c*/
/*     //      The 1st element for S_T is S_T[1] while S_T[0] is s_0 (initial condition).   ansi-c*/
/*     //      The 1st element for Y_T, however, is Y_T[nlags_max+1-1].   ansi-c*/
/*     //See (42.3) on p.42 in the SWZII NOTES.   ansi-c*/


/*     //--- Local variables   ansi-c*/
   int comst_c;   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
   int st_c, stm1_c, st_v;
   int comsti_c;   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
   int sti, sti_c, stm1i_c, sti_v;
   int comstp1i_c;   /*  composite (s_{t+1}c, s_tc)   ansi-c*/
   int stp1i, stp1i_c, stp1i_v;
   int tbase0, tp1;
   double logdet_Dtdata, loglh_timet;
   static int record_tbase1_or_inpt_or_tp1 = 0;
   static int passonce;
   double prob_previous_regimes;
/*     //=== Accessible variables   ansi-c*/
   int ny = kalfilmsinputs_ps->ny;
   int nz = kalfilmsinputs_ps->nz;
   int nRc = kalfilmsinputs_ps->nRc;
   int nRstc = kalfilmsinputs_ps->nRstc;
   int nRv = kalfilmsinputs_ps->nRv;
   int T = kalfilmsinputs_ps->T;
   int indxIndRegimes = kalfilmsinputs_ps->indxIndRegimes;
   int **Index = smodel_ps->sv->index;   /*  Regime-switching states.   ansi-c*/
/*            //smodel_ps->sv->index is for our new code.   ansi-c*/
/*            //  For old code (before 9 April 08 and before dsge_switch is created), use smodel_ps->sv->Index;   ansi-c*/
   TSdvector z0_sdv;
/*     //+ input arguments.   ansi-c*/
   TSdmatrix *yt_dm = kalfilmsinputs_ps->yt_dm;          /*  ny-by-T.   ansi-c*/
   TSdmatrix *at_dm = kalfilmsinputs_ps->at_dm;          /*  ny-by-nRc.   ansi-c*/
   TSdcell *Ht_dc = kalfilmsinputs_ps->Ht_dc;            /*  ny-by-nz-by-nRc.   ansi-c*/
   TSdcell *Rt_dc = kalfilmsinputs_ps->Rt_dc;            /*  ny-by-ny-by-nRv.  Covariance matrix for the measurement equation.   ansi-c*/
   TSdcell *Gt_dc = kalfilmsinputs_ps->Gt_dc;            /*  nz-by-ny-by-nRv.  Cross-covariance.   ansi-c*/
/*     //   ansi-c*/
   TSdmatrix *bt_dm = kalfilmsinputs_ps->bt_dm;          /*  nz-by-nRc.   ansi-c*/
   TSdcell *Ft_dc = kalfilmsinputs_ps->Ft_dc;            /*  nz-by-nz-by-nRc.   ansi-c*/
   TSdcell *Vt_dc = kalfilmsinputs_ps->Vt_dc;            /*  nz-by-nz-by-nRv.  Covariance matrix for the state equation.   ansi-c*/
/*     //   ansi-c*/
   TSdmatrix *z0_dm = kalfilmsinputs_ps->z0_dm;         /*  nz-by-nRc*nRv or nz-by-nRv, depending on indxIndRegimes.   ansi-c*/
   TSdcell *P0_dc = kalfilmsinputs_ps->P0_dc;           /*  nz-by-nz-by-nRc*nRv or nz-by-nRv, depending on indxIndRegimes.   ansi-c*/
/*     //+ Output arguments.   ansi-c*/
   TSdcell *zt_tm1_dc = kalfilmsinputs_ps->zt_tm1_dc;  /*  nz-by-nRc*nRv-by-T if indxIndRegimes==1, nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.   ansi-c*/
   TSdfourth *Pt_tm1_d4 = kalfilmsinputs_ps->Pt_tm1_d4;      /*  nz-by-nz-by-nRc*nRv-T if indxIndRegimes==1, nz-by-nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.   ansi-c*/
/*     //=== Work arguments.   ansi-c*/
   int nz2 = square(nz);
   TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *Wnz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdmatrix *W2nz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdvector *wP0_dv = CreateVector_lf(nz2);
/*     //+   ansi-c*/
   TSdvector yt_sdv, at_sdv, btp1_sdv;   /*  zt_tm1_sdv, ztp1_t_sdv,   ansi-c*/
   TSdvector *wny_dv = CreateVector_lf(ny);
   TSdmatrix *Wnzbyny_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *W2nzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *PHtran_tdata_dm = CreateMatrix_lf(nz,ny);
   TSdvector *etdata_dv = CreateVector_lf(ny);
   TSdmatrix *Dtdata_dm = CreateMatrix_lf(ny,ny);
   TSdmatrix *Kt_tdata0_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *Kt_tdata_dm = CreateMatrix_lf(nz,ny);
/*     //--- For eigenvalue decompositions   ansi-c*/
   int ki;
   int errflag;
   double eigmax;
   TSdzvector *evals_dzv = evals_dzv = CreateVector_dz(nz);
   TSdvector *evals_abs_dv = CreateVector_lf(nz);  /*  Absolute eigenvalues.   ansi-c*/
/*     //--- For updating zt_tm1_dm and Pt_tm1.   ansi-c*/
   TSdvector *ztp1_t_dv = CreateVector_lf(z0_dm->nrows);
   TSdmatrix *Ptp1_t_dm = CreateMatrix_lf(nz, nz);
   TSdvector *ztp1_dv = CreateVector_lf(z0_dm->nrows);
   TSdmatrix *Ptp1_dm = CreateMatrix_lf(nz, nz);



   if (smodel_ps->sv->nstates != z0_dm->ncols)  fn_DisplayError("kalman.c/tz_logTimetLH_kalfilms_1st_approx():\n"
                     "  Make sure that the column dimension of z0_dm is the same as smodel_ps->sv->nstates");
   if (indxIndRegimes && (nRc>1) && (nRv>1))
      if (smodel_ps->sv->n_state_variables != 2)  fn_DisplayError("kalman.c/tz_Refresh_z_T7P_T_in_kalfilms_1st_approx():\n"
           " Number of state variables must be coincide with indxIndRegimes");

   tbase0 = (tp1=inpt) - 1;

   z0_sdv.n = z0_dm->nrows;
   z0_sdv.flag = V_DEF;
/*     //   ansi-c*/
   at_sdv.n = yt_sdv.n = yt_dm->nrows;
   at_sdv.flag = yt_sdv.flag = V_DEF;
   btp1_sdv.n = bt_dm->nrows;
   btp1_sdv.flag = V_DEF;


/*     //======= Initial condition. =======   ansi-c*/
   if (tbase0==0)
   {
      for (sti=smodel_ps->sv->nstates-1; sti>=0;  sti--)
      {
         if (indxIndRegimes)
         {
            if (nRc==1)        /*  Volatility.   ansi-c*/
            {
               comsti_c = sti_c = 0;
               sti_v = sti;
            }
            else if ((nRv>1) && (nRc>nRstc))   /*  Trend inflation, both sc_t and sc_{t-1} enters coefficient regime.   ansi-c*/
            {
               comsti_c = Index[sti][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
               sti_v = Index[sti][1];   /*  volatility state s_tv   ansi-c*/
               sti_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][0];   /*  coefficient regime at t.   ansi-c*/
               stm1i_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
            }
            else if ((nRv==1) && (nRc>nRstc))
            {
               comsti_c = Index[sti][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
               sti_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][0];   /*  coefficient regime at t.   ansi-c*/
               stm1i_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
               sti_v = 0;
            }
            else if ((nRv==1) && (nRc==nRstc))
            {
               comsti_c  = sti_c = sti;
               sti_v = 0;
            }
            else if ((nRv>1) && (nRc==nRstc))   /*  only sc_t enters coefficient regime.   ansi-c*/
            {
               comsti_c = sti_c = Index[sti][0];
               sti_v = Index[sti][1];
            }
         }
         else  //Syncronized regimes.
         {
            if (nRc>nRstc)
            {
               comsti_c = Index[sti][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
               sti_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][0];   /*  coefficient regime at t.   ansi-c*/
               stm1i_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
               sti_v = sti_c;
            }
            else
               comsti_c = sti_c = sti_v = sti;
         }


         if (!kalfilmsinputs_ps->indxIni)
         {
            InitializeDiagonalMatrix_lf(Wnzbynz_dm, 1.0);   /*  To be used for I(nz) -   ansi-c*/
            InitializeDiagonalMatrix_lf(Wnz2bynz2_dm, 1.0);   /*  To be used for I(nz2) -   ansi-c*/

/*              //=== Eigenanalysis to determine the roots to ensure boundedness.   ansi-c*/
            errflag = eigrgen(evals_dzv, (TSdzmatrix *)NULL, (TSdzmatrix *)NULL, Ft_dc->C[comsti_c]);
            if (errflag)  fn_DisplayError("kalman.c/tz_Refresh_z_T7P_T_in_kalfilms_1st_approx(): eigen decomposition failed");
            for (ki=nz-1; ki>=0; ki--)  evals_abs_dv->v[ki] = sqrt(square(evals_dzv->real->v[ki]) + square(evals_dzv->imag->v[ki]));
            evals_abs_dv->flag = V_DEF;
            eigmax = MaxVector(evals_abs_dv);
            if (eigmax < (1.0+1.0e-14))
            {
/*                 //--- Getting z0_dv: zt_tm1(:,1) = (eye(n_z)-F(:,:,1))\b(:,1);   ansi-c*/
               MatrixMinusMatrix(Wnzbynz_dm, Wnzbynz_dm, Ft_dc->C[comsti_c]);
               z0_sdv.v = z0_dm->M + z0_sdv.n*sti;
               CopySubmatrix2vector(&z0_sdv, 0, bt_dm, 0, comsti_c, bt_dm->nrows);
               bdivA_rgens(&z0_sdv, &z0_sdv, '\\', Wnzbynz_dm);
/*                           //Done with Wnzbynz_dm.   ansi-c*/
/*                 //--- Getting P0_dm: Pt_tm1(:,:,1) = reshape((eye(n_z^2)-kron(F(:,:,1),F(:,:,1)))\V1(:),n_z,n_z);   ansi-c*/
               tz_kron(W2nz2bynz2_dm, Ft_dc->C[comsti_c], Ft_dc->C[comsti_c]);
               MatrixMinusMatrix(Wnz2bynz2_dm, Wnz2bynz2_dm, W2nz2bynz2_dm);
               CopySubmatrix2vector(wP0_dv, 0, Vt_dc->C[sti_v], 0, 0, nz2);
/*  //=== ???????? For debugging purpose.   ansi-c*/
/*  //if ((inpt<2) && (st==0))   ansi-c*/
/*  //{   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "%%st=%d, inpt=%d, and sti=%d\n", st, inpt, sti);   ansi-c*/

/*  //   fprintf(FPTR_DEBUG, "wP0_dv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, wP0_dv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "Vt_dc->C[sti_v=%d]:\n", sti_v);   ansi-c*/
/*  //   WriteMatrix(FPTR_DEBUG, Vt_dc->C[sti_v], " %10.5f ");   ansi-c*/

/*  //   fflush(FPTR_DEBUG);   ansi-c*/

/*  //}   ansi-c*/
               bdivA_rgens(wP0_dv, wP0_dv, '\\', Wnz2bynz2_dm);
               CopySubvector2matrix_unr(P0_dc->C[sti], 0, 0, wP0_dv, 0, nz2);
/*                            //Done with all w*_dv and W*_dm.   ansi-c*/
            }
            else
            {
               printf("\n-----------------\n");
               printf("\nIn regime comsti_c=%d and sti_v=%d and at time=%d\n", comsti_c, sti_v, 0);
               fn_DisplayError("kalman.c/tz_Refresh_z_T7P_T_in_kalfilms_1st_approx(): the system is non-stationary solutions\n"
                             "    and the initial conditions must be supplied by, say, input arguments");
               fflush(stdout);
            }
         }
      }
      z0_dm->flag = M_GE;
      CopyMatrix0(zt_tm1_dc->C[0], z0_dm);   /*  At time t=0.   ansi-c*/
      CopyCell0(Pt_tm1_d4->F[0], P0_dc);                               /*  At time t=0.   ansi-c*/
   }


/*     //======================================================   ansi-c*/
/*     //= Getting the logLH at time tbase0 or time inpt.   ansi-c*/
/*     //======================================================   ansi-c*/
   if (indxIndRegimes )
   {
      if (nRc==1)        /*  Volatility.   ansi-c*/
      {
         comst_c = st_c = 0;
         st_v = st;
      }
      else if ((nRv>1) && (nRc>nRstc))   /*  Trend inflation, both sc_t and sc_{t-1} enters coefficient regime.   ansi-c*/
      {
         if (smodel_ps->sv->n_state_variables != 2)  fn_DisplayError("kalman.c/kalfilms_timet_1st_approx():\n"
                         "  Number of state variables must be coincide with indxIndRegimes");

         comst_c = Index[st][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
         st_v = Index[st][1];   /*  volatility state s_tv   ansi-c*/
         st_c = smodel_ps->sv->state_variable[0]->lag_index[comst_c][0];   /*  coefficient regime at t.   ansi-c*/
         stm1_c = smodel_ps->sv->state_variable[0]->lag_index[comst_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
      }
      else if ((nRv==1) && (nRc>nRstc))
      {
         comst_c = Index[st][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
         st_c = smodel_ps->sv->state_variable[0]->lag_index[comst_c][0];   /*  coefficient regime at t.   ansi-c*/
         stm1_c = smodel_ps->sv->state_variable[0]->lag_index[comst_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
         st_v = 0;
      }
      else if ((nRv==1) && (nRc==nRstc))
      {
         comst_c  = st_c = st;
         st_v = 0;
      }
      else if ((nRv>1) && (nRc==nRstc))   /*  only sc_t enters coefficient regime.   ansi-c*/
      {
         if (smodel_ps->sv->n_state_variables != 2)  fn_DisplayError("kalman.c/kalfilms_timet_1st_approx():\n"
                         "  Number of state variables must be coincide with indxIndRegimes");

         comst_c = st_c = Index[st][0];
         st_v = Index[st][1];
      }
   }
   else   //Syncronized regimes
   {
       if (nRc>nRstc)
       {
          comst_c = Index[st][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
          st_c = smodel_ps->sv->state_variable[0]->lag_index[comst_c][0];   /*  coefficient regime at t.   ansi-c*/
          stm1_c = smodel_ps->sv->state_variable[0]->lag_index[comst_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
          st_v = st_c;
       }
       else
          comst_c = st_c = st_v = st;
   }


   z0_sdv.n = zt_tm1_dc->C[0]->nrows;
   z0_sdv.flag = V_DEF;
/*     //   ansi-c*/
   at_sdv.n = yt_sdv.n = yt_dm->nrows;
   at_sdv.flag = yt_sdv.flag = V_DEF;

/*     //====== Computing the conditional LH at time t. ======   ansi-c*/
/*     //--- Setup.   ansi-c*/
   MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_d4->F[tbase0]->C[st], Ht_dc->C[comst_c], 1.0, 0.0, 'N', 'T');

/*     //--- Data.   ansi-c*/
/*     //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata where tdata = tbase0 = inpt-1.   ansi-c*/
   yt_sdv.v = yt_dm->M + tbase0*yt_dm->nrows;
   at_sdv.v = at_dm->M + comst_c*at_dm->nrows;     /*  comst_c: coefficient regime at time tbase0.   ansi-c*/
   z0_sdv.v = zt_tm1_dc->C[tbase0]->M + z0_sdv.n*st;   /*  st: regime at time tbase0 for zt_tm1.   ansi-c*/
   VectorMinusVector(etdata_dv, &yt_sdv, &at_sdv);
   MatrixTimesVector(etdata_dv, Ht_dc->C[comst_c], &z0_sdv, -1.0, 1.0, 'N');
/*     //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tbase0);   ansi-c*/
   CopyMatrix0(Dtdata_dm, Rt_dc->C[st_v]);
   MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[comst_c], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
   ScalarTimesMatrixSquare(Dtdata_dm, 0.5, Dtdata_dm, 'T', 0.5);   /*  Making it symmetric against some rounding errors.   ansi-c*/
/*                        //This making-symmetric is very IMPORTANT; otherwise, we will get the matrix being singular message   ansi-c*/
/*                        //    and eigenvalues being negative for the SPD matrix, etc.  Then the likelihood becomes either   ansi-c*/
/*                        //    a bad number or a complex number.   ansi-c*/
   Dtdata_dm->flag = Dtdata_dm->flag | M_SU | M_SL;


/*     //--- Forming the log conditional likelihood at t.   ansi-c*/
   if (!isfinite(logdet_Dtdata=logdetspd(Dtdata_dm)))  return (loglh_timet = -NEARINFINITY);
   bdivA_rgens(wny_dv, etdata_dv, '/', Dtdata_dm);
/*  //if ((inpt>82) && (inpt<86) )   ansi-c*/
/*  //{   ansi-c*/
/*  //   //Must be declared at the top of this "if" block.   ansi-c*/
/*  //   int kip1;   ansi-c*/
/*  //   double tmp_Dtdata;   ansi-c*/
/*  //   double tmp_expterm;   ansi-c*/

/*  //   fprintf(FPTR_DEBUG, "%%------------------------\n");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "%%st=%d and inpt=%d\n", st, inpt);   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "loglh_timet = %10.5f;\n", loglh_timet);   ansi-c*/


/*  //   fprintf(FPTR_DEBUG, "wny_dv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, wny_dv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "etdata_dv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, etdata_dv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "Dtdata_dm:\n");   ansi-c*/
/*  //   WriteMatrix(FPTR_DEBUG, Dtdata_dm, " %.16e ");   ansi-c*/

/*  //   fflush(FPTR_DEBUG);   ansi-c*/
/*  //}   ansi-c*/
   loglh_timet = -(0.5*ny)*LOG2PI - 0.5*logdet_Dtdata - 0.5*VectorDotVector(wny_dv, etdata_dv);
/*                           //Done with all w*_dv.   ansi-c*/




/*  //=== ???????? For debugging purpose.   ansi-c*/
if (inpt==1)
{
   double wk1, wk2;

   wk1 = logdet_Dtdata;
   wk2 = VectorDotVector(wny_dv, etdata_dv);
   fprintf(FPTR_DEBUG, "logdet_Dtdata = %10.5f\n", wk1);
   fprintf(FPTR_DEBUG, "VectorDotVector(wny_dv, etdata_dv) = %10.5f\n", wk2);
   fprintf(FPTR_DEBUG, "----- etdata_dv: \n");
   WriteVector(FPTR_DEBUG, etdata_dv, " %10.5f ");
   fprintf(FPTR_DEBUG, "----- yt_dv: \n");
   WriteVector(FPTR_DEBUG, &yt_sdv, " %10.5f ");
   fprintf(FPTR_DEBUG, "----- at_dv: \n");
   WriteVector(FPTR_DEBUG, &at_sdv, " %10.5f ");
   fprintf(FPTR_DEBUG, "----- z0_dv: \n");
   WriteVector(FPTR_DEBUG, &z0_sdv, " %10.5f ");
   fprintf(FPTR_DEBUG, "----- Ht_dc->C[comst_c=%d]:\n", comst_c);
   WriteMatrix(FPTR_DEBUG, Ht_dc->C[comst_c], " %10.5f ");

   fprintf(FPTR_DEBUG, "\n\n");

}
/*  //   ansi-c*/
fprintf(FPTR_DEBUG, " %10.5f\n", loglh_timet);
fflush(FPTR_DEBUG);


/*  //=== ???????? For debugging purpose.   ansi-c*/
/*  //fprintf(FPTR_DEBUG, "------------------------\n");   ansi-c*/
/*  //fprintf(FPTR_DEBUG, "st=%d and inpt=%d\n", st, inpt);   ansi-c*/
/*  //fprintf(FPTR_DEBUG, "loglh_timet = %10.5f\n", loglh_timet);   ansi-c*/
/*  //fprintf(FPTR_DEBUG, "&yt_sdv:\n");   ansi-c*/
/*  //WriteVector(FPTR_DEBUG, &yt_sdv, " %10.5f ");   ansi-c*/
/*  ////WriteVector(FPTR_DEBUG, etdata_dv, " %10.5f ");   ansi-c*/
/*  ////fprintf(FPTR_DEBUG, "\n");   ansi-c*/
/*  ////WriteMatrix(FPTR_DEBUG, Dtdata_dm, " %10.5f ");   ansi-c*/
/*  //fflush(FPTR_DEBUG);   ansi-c*/


/*  //=== ???????? For debugging purpose.   ansi-c*/
/*  //if ((inpt>82) && (inpt<86) )   ansi-c*/
/*  //if (inpt<2)   ansi-c*/
/*  //{   ansi-c*/
/*  //   //Must be declared at the top of this "if" block.   ansi-c*/
/*  //   int kip1;   ansi-c*/
/*  //   double tmp_Dtdata;   ansi-c*/
/*  //   double tmp_expterm;   ansi-c*/

/*  //   fprintf(FPTR_DEBUG, "%%------------------------\n");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "%%st=%d and inpt=%d\n", st, inpt);   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "loglh_timet = %10.5f;\n", loglh_timet);   ansi-c*/


/*  //   tmp_Dtdata = logdeterminant(Dtdata_dm);   ansi-c*/
/*  //   tmp_expterm = VectorDotVector(wny_dv, etdata_dv);   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "logdeterminant(Dtdata_dm) = %10.5f;\n", tmp_Dtdata);   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "VectorDotVector(wny_dv, etdata_dv) = %10.5f;\n", tmp_expterm);   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "wny_dv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, wny_dv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "etdata_dv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, etdata_dv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "&yt_sdv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, &yt_sdv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "&at_sdv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, &at_sdv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "&z0_sdv:\n");   ansi-c*/
/*  //   WriteVector(FPTR_DEBUG, &z0_sdv, " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "Ht_dc->C[comst_c=%d]:\n",comst_c);   ansi-c*/
/*  //   WriteMatrix(FPTR_DEBUG, Ht_dc->C[comst_c], " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "Rt_dc->C[st_v=%d]:\n", st_v);   ansi-c*/
/*  //   WriteMatrix(FPTR_DEBUG, Rt_dc->C[st_v], " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "Pt_tm1_d4->F[tbase0]->C[st = %d]:\n",st);   ansi-c*/
/*  //   WriteMatrix(FPTR_DEBUG, Pt_tm1_d4->F[tbase0]->C[st], " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "Dtdata_dm:\n");   ansi-c*/
/*  //   WriteMatrix(FPTR_DEBUG, Dtdata_dm, " %10.5f ");   ansi-c*/




/*  ////   WriteMatrix(FPTR_DEBUG, Dtdata_dm, " %10.5f ");   ansi-c*/
/*  ////   fprintf(FPTR_DEBUG, "zt_tm1_dc->C[tbase0]:\n");   ansi-c*/
/*  ////   WriteMatrix(FPTR_DEBUG, zt_tm1_dc->C[tbase0], " %10.5f ");   ansi-c*/
/*  ////   //WriteVector(FPTR_DEBUG, &z0_sdv, " %10.5f ");   ansi-c*/
/*  ////   //fprintf(FPTR_DEBUG, "\n");   ansi-c*/
/*  ////   fprintf(FPTR_DEBUG, "bt_dm = [\n");   ansi-c*/
/*  ////   WriteMatrix(FPTR_DEBUG, bt_dm, " %10.5f ");   ansi-c*/
/*  ////   fprintf(FPTR_DEBUG, "];\n");   ansi-c*/

/*  ////   fprintf(FPTR_DEBUG, "et:\n");   ansi-c*/
/*  ////   WriteVector(FPTR_DEBUG, etdata_dv, " %10.5f ");   ansi-c*/
/*  ////   fprintf(FPTR_DEBUG, "yt_dv=[\n");   ansi-c*/
/*  ////   WriteVector(FPTR_DEBUG, &yt_sdv, " %10.5f ");   ansi-c*/
/*  ////   fprintf(FPTR_DEBUG, "]';\n");   ansi-c*/

/*  ////   fprintf(FPTR_DEBUG, "at_dv=[\n");   ansi-c*/
/*  ////   WriteVector(FPTR_DEBUG, &at_sdv, " %10.5f ");   ansi-c*/
/*  ////   fprintf(FPTR_DEBUG, "]';\n");   ansi-c*/


/*  ////   for (ki=0; ki<Ht_dc->ncells; ki++)   ansi-c*/
/*  ////   {   ansi-c*/
/*  ////      kip1 = ki+1;   ansi-c*/
/*  ////      fprintf(FPTR_DEBUG, "Ht_dc(:,:,%d)=[\n", kip1);   ansi-c*/
/*  ////      WriteMatrix(FPTR_DEBUG, Ht_dc->C[ki], " %10.5f ");   ansi-c*/
/*  ////      fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*  ////   }   ansi-c*/
/*  ////   for (ki=0; ki<Ft_dc->ncells; ki++)   ansi-c*/
/*  ////   {   ansi-c*/
/*  ////      kip1 = ki+1;   ansi-c*/
/*  ////      fprintf(FPTR_DEBUG, "Ft_dc(:,:,%d)=[\n", kip1);   ansi-c*/
/*  ////      WriteMatrix(FPTR_DEBUG, Ft_dc->C[ki], " %10.5f ");   ansi-c*/
/*  ////      fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*  ////   }   ansi-c*/
/*  ////   for (ki=0; ki<Vt_dc->ncells; ki++)   ansi-c*/
/*  ////   {   ansi-c*/
/*  ////      kip1 = ki+1;   ansi-c*/
/*  ////      fprintf(FPTR_DEBUG, "Vt_dc(:,:,%d)=[\n", kip1);   ansi-c*/
/*  ////      WriteMatrix(FPTR_DEBUG, Vt_dc->C[ki], " %10.5f ");   ansi-c*/
/*  ////      fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*  ////   }   ansi-c*/
/*  //   fflush(FPTR_DEBUG);   ansi-c*/
/*  //}   ansi-c*/


/*     //======================================================   ansi-c*/
/*     //= Updating zt_tm1 and Pt_tm1 for next perid tp1.   ansi-c*/
/*     //= tdata = tbase0 is base-0 timing.   ansi-c*/
/*     //======================================================   ansi-c*/
   if (inpt > record_tbase1_or_inpt_or_tp1)   /*  This condition always satisfies at the 1st period (which is inpt=1).   ansi-c*/
   {
      passonce = 0;
      record_tbase1_or_inpt_or_tp1 = inpt;
   }
   if (!passonce)
   {
      for (stp1i=smodel_ps->sv->nstates-1; stp1i>=0;  stp1i--)
      {
         if (indxIndRegimes)
         {
            if (nRc==1)        /*  Volatility.   ansi-c*/
            {
               comstp1i_c = stp1i_c = 0;
               stp1i_v = stp1i;
            }
            else if ((nRv>1) && (nRc>nRstc))   /*  Trend inflation, both sc_t and sc_{t-1} enters coefficient regime.   ansi-c*/
            {
               comstp1i_c = Index[stp1i][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
               stp1i_v = Index[stp1i][1];   /*  volatility state s_tv   ansi-c*/
               stp1i_c = smodel_ps->sv->state_variable[0]->lag_index[comstp1i_c][0];   /*  coefficient regime at t.   ansi-c*/
/*                 //sti_c = smodel_ps->sv->state_variable[0]->lag_index[comstp1i_c][1];  //coefficient regime at t-1: tm1: t-1;   ansi-c*/
            }
            else if ((nRv==1) && (nRc>nRstc))
            {
               comstp1i_c = Index[stp1i][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
               stp1i_c = smodel_ps->sv->state_variable[0]->lag_index[comstp1i_c][0];   /*  coefficient regime at t.   ansi-c*/
/*                 //sti_c = smodel_ps->sv->state_variable[0]->lag_index[comstp1i_c][1];  //coefficient regime at t-1: tm1: t-1;   ansi-c*/
               stp1i_v = 0;
            }
            else if ((nRv==1) && (nRc==nRstc))
            {
               comstp1i_c  = stp1i_c = stp1i;
               stp1i_v = 0;
            }
            else if ((nRv>1) && (nRc==nRstc))   /*  only sc_t enters coefficient regime.   ansi-c*/
            {
               comstp1i_c = stp1i_c = Index[stp1i][0];
               stp1i_v = Index[stp1i][1];
            }
         }
         else  //Syncronized regimes.
         {
            if (nRc>nRstc)
            {
               comstp1i_c = Index[stp1i][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
               stp1i_c = smodel_ps->sv->state_variable[0]->lag_index[comstp1i_c][0];   /*  coefficient regime at t.   ansi-c*/
/*                 //sti_c = smodel_ps->sv->state_variable[0]->lag_index[comstp1i_c][1];  //coefficient regime at t-1: tm1: t-1;   ansi-c*/
               stp1i_v = stp1i_c;
            }
            else
               comstp1i_c = stp1i_c = stp1i_v = stp1i;
         }


         InitializeConstantVector_lf(ztp1_dv, 0.0);   /*  To be summed over sti.   ansi-c*/
         InitializeConstantMatrix_lf(Ptp1_dm, 0.0);   /*  To be summed over sti.   ansi-c*/

         for (sti=smodel_ps->sv->nstates-1; sti>=0;  sti--)
         {
            if (indxIndRegimes)
            {
               if (nRc==1)        /*  Volatility.   ansi-c*/
               {
                  comsti_c = sti_c = 0;
                  sti_v = sti;
               }
               else if ((nRv>1) && (nRc>nRstc))   /*  Trend inflation, both sc_t and sc_{t-1} enters coefficient regime.   ansi-c*/
               {
                  comsti_c = Index[sti][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
                  sti_v = Index[sti][1];   /*  volatility state s_tv   ansi-c*/
                  sti_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][0];   /*  coefficient regime at t.   ansi-c*/
                  stm1i_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
               }
               else if ((nRv==1) && (nRc>nRstc))
               {
                  comsti_c = Index[sti][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
                  sti_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][0];   /*  coefficient regime at t.   ansi-c*/
                  stm1i_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
                  sti_v = 0;
               }
               else if ((nRv==1) && (nRc==nRstc))
               {
                  comsti_c  = sti_c = sti;
                  sti_v = 0;
               }
               else if ((nRv>1) && (nRc==nRstc))   /*  only sc_t enters coefficient regime.   ansi-c*/
               {
                  comsti_c = sti_c = Index[sti][0];
                  sti_v = Index[sti][1];
               }
            }
            else  //Syncronized regimes.
            {
               if (nRc>nRstc)
               {
                  comsti_c = Index[sti][0];   /*  composite (s_tc, s_{t-1}c)   ansi-c*/
                  sti_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][0];   /*  coefficient regime at t.   ansi-c*/
                  stm1i_c = smodel_ps->sv->state_variable[0]->lag_index[comsti_c][1];   /*  coefficient regime at t-1: tm1: t-1;   ansi-c*/
                  sti_v = sti_c;
               }
               else
                  comsti_c = sti_c = sti_v = sti;
            }


/*              //--- Setup.   ansi-c*/
            MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_d4->F[tbase0]->C[sti], Ht_dc->C[comsti_c], 1.0, 0.0, 'N', 'T');

/*              //--- Data.   ansi-c*/
/*              //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata where tdata = tbase0 = inpt-1.   ansi-c*/
            yt_sdv.v = yt_dm->M + tbase0*yt_dm->nrows;
            at_sdv.v = at_dm->M + comsti_c*at_dm->nrows;   /*  comsti_c: coefficient regime at time tbase0.   ansi-c*/
            z0_sdv.v = zt_tm1_dc->C[tbase0]->M + z0_sdv.n*sti;   /*  sti: regime at time tbase0.   ansi-c*/
            VectorMinusVector(etdata_dv, &yt_sdv, &at_sdv);
            MatrixTimesVector(etdata_dv, Ht_dc->C[comsti_c], &z0_sdv, -1.0, 1.0, 'N');
/*              //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tbase0);   ansi-c*/
            CopyMatrix0(Dtdata_dm, Rt_dc->C[sti_v]);
            MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[comsti_c], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
/*                                          //Done with z0_sdv.v.   ansi-c*/
            ScalarTimesMatrixSquare(Dtdata_dm, 0.5, Dtdata_dm, 'T', 0.5);   /*  Making it symmetric against some rounding errors.   ansi-c*/
/*                                 //This making-symmetric is very IMPORTANT; otherwise, we will get the matrix being singular message   ansi-c*/
/*                                 //    and eigenvalues being negative for the SPD matrix, etc.  Then the likelihood becomes either   ansi-c*/
/*                                 //    a bad number or a complex number.   ansi-c*/
            Dtdata_dm->flag = Dtdata_dm->flag | M_SU | M_SL;


/*              //=== Updating for next period by integrating out sti..   ansi-c*/
            if (tp1<T)
            {
/*                 //Updating only up to tbase0=T-2.  The values at tp1=T or tbase0=T-1 will not be used in the likelihood function.   ansi-c*/

/*                 //--- Ktp1_t = (F_tp1*PHtran_t+G(:,:,t))/Dt;   ansi-c*/
/*                 //--- Kt_tdata = (Ft*PHtran_tdata+G(:,:,tdata))/Dtdata where t=tp1 and tdata=t.   ansi-c*/
               CopyMatrix0(Kt_tdata0_dm, Gt_dc->C[sti_v]);
               MatrixTimesMatrix(Kt_tdata0_dm, Ft_dc->C[stp1i_c], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
               BdivA_rrect(Kt_tdata_dm, Kt_tdata0_dm, '/', Dtdata_dm);
/*                 //+ zt_tm1(:,t) = b(:,t) + Ft*zt_tm1(:,tdata) + Kt_tdata*etdata where t=tp1 and tm1=t.   ansi-c*/
               MatrixTimesVector(ztp1_t_dv, Ft_dc->C[stp1i_c], &z0_sdv, 1.0, 0.0, 'N');
               MatrixTimesVector(ztp1_t_dv, Kt_tdata_dm, etdata_dv, 1.0, 1.0, 'N');
               btp1_sdv.v = bt_dm->M + stp1i_c*btp1_sdv.n;
               VectorPlusMinusVectorUpdate(ztp1_t_dv, &btp1_sdv, 1.0);
/*                 //+ Pt_tm1(:,:,t) = Ft*Ptdata*Fttran - Kt_tdata*Dtdata*Kt_tdatatran + V(:,:,t);   ansi-c*/
               CopyMatrix0(Ptp1_t_dm, Vt_dc->C[stp1i_v]);
               MatrixTimesMatrix(Wnzbyny_dm, Kt_tdata_dm, Dtdata_dm, 1.0, 0.0, 'N', 'N');
               MatrixTimesMatrix(Wnzbynz_dm, Wnzbyny_dm, Kt_tdata_dm, 1.0, 0.0, 'N', 'T');
               MatrixPlusMinusMatrixUpdate(Ptp1_t_dm, Wnzbynz_dm, -1.0);
/*                                       //Done with all W*_dm.   ansi-c*/
               MatrixTimesMatrix(Wnzbynz_dm, Ft_dc->C[stp1i_c], Pt_tm1_d4->F[tbase0]->C[sti], 1.0, 0.0, 'N', 'N');
               MatrixTimesMatrix(W2nzbynz_dm, Wnzbynz_dm, Ft_dc->C[stp1i_c], 1.0, 0.0, 'N', 'T');
               MatrixPlusMatrixUpdate(Ptp1_t_dm, W2nzbynz_dm);
/*                                       //Done with all W*_dm.   ansi-c*/

/*                 //--- Integrating out the state at tbase0 using P(s_t|Y_{t-1}, theta) = ElementV(smodel_ps->Z[inpt],s_{inpt}_i).   ansi-c*/
/*                 //---   Note tbase0 = inpt-1 because the data in DW code (ElementV) is base-1.   ansi-c*/
/*                 //---   Note at this point, we cannot access to P(s_t|Y_t, theta) = ElementV(smodel_ps->V[inpt],s_{inpt}_i)   ansi-c*/
/*                 //---      through DW's code.  But we can modify my own code to do this later.   ansi-c*/
               prob_previous_regimes = ElementV(smodel_ps->Z[inpt],sti);
               ScalarTimesVectorUpdate(ztp1_dv, prob_previous_regimes, ztp1_t_dv);
               ScalarTimesMatrix(Ptp1_dm, prob_previous_regimes, Ptp1_t_dm, 1.0);
               Ptp1_dm->flag = M_GE | M_SU | M_SL;
/*                                           //Done with ztp1_t_dv and Ptp1_t_dm.   ansi-c*/
            }
         }
/*           //--- Filling zt_tm1 and Pt_tm1 for next period   ansi-c*/
         if (tp1<T)
         {
            z0_sdv.v = zt_tm1_dc->C[tp1]->M + z0_sdv.n*stp1i;   /*  stp1i: regime at time tp1.   ansi-c*/
            CopyVector0(&z0_sdv, ztp1_dv);
            CopyMatrix0(Pt_tm1_d4->F[tp1]->C[stp1i], Ptp1_dm);   /*  stp1i: regime at time tp1.   ansi-c*/
/*                                             //Done with ztp1_dv, z0_sdv, Ptp1_dm.   ansi-c*/
         }
      }
      if (tp1<T)
         zt_tm1_dc->C[tp1]->flag = M_GE;
   }


/*  //=== ???????? For debugging purpose.   ansi-c*/
/*  //if ((inpt>60) && (inpt<65) )  //if (inpt<5)   ansi-c*/
/*  //{   ansi-c*/
/*  //   int kip1;  //Must be declared at the top of this "if" block.   ansi-c*/

/*  //   fprintf(FPTR_DEBUG, "zt_tm1t=[\n");   ansi-c*/
/*  //   WriteMatrix(FPTR_DEBUG, zt_tm1_dc->C[tbase0], " %10.5f ");   ansi-c*/
/*  //   fprintf(FPTR_DEBUG, "];\n");   ansi-c*/

/*  //   for (ki=0; ki<Pt_tm1_d4->F[tbase0]->ncells; ki++)   ansi-c*/
/*  //   {   ansi-c*/
/*  //      kip1 = ki+1;   ansi-c*/
/*  //      fprintf(FPTR_DEBUG, "Pt_tm1_d4t(:,:,%d)=[\n", kip1);   ansi-c*/
/*  //      WriteMatrix(FPTR_DEBUG, Pt_tm1_d4->F[tbase0]->C[ki], " %10.5f ");   ansi-c*/
/*  //      fprintf(FPTR_DEBUG, "];\n");   ansi-c*/
/*  //   }   ansi-c*/

/*  //   fflush(FPTR_DEBUG);   ansi-c*/
/*  //}   ansi-c*/


/*  //=== ???????? For debugging purpose.   ansi-c*/
fprintf(FPTR_DEBUG, " loglh_timet = %10.5f\n", loglh_timet);
fflush(FPTR_DEBUG);


/*     //===   ansi-c*/
   DestroyVector_dz(evals_dzv);
   DestroyVector_lf(evals_abs_dv);
   DestroyMatrix_lf(Wnzbynz_dm);
   DestroyMatrix_lf(Wnz2bynz2_dm);
   DestroyMatrix_lf(W2nz2bynz2_dm);
   DestroyVector_lf(wP0_dv);
/*     //   ansi-c*/
   DestroyVector_lf(wny_dv);
   DestroyMatrix_lf(Wnzbyny_dm);
   DestroyMatrix_lf(W2nzbynz_dm);
   DestroyMatrix_lf(PHtran_tdata_dm);
   DestroyVector_lf(etdata_dv);
   DestroyMatrix_lf(Dtdata_dm);
   DestroyMatrix_lf(Kt_tdata0_dm);
   DestroyMatrix_lf(Kt_tdata_dm);
/*     //   ansi-c*/
   DestroyVector_lf(ztp1_t_dv);
   DestroyMatrix_lf(Ptp1_t_dm);
   DestroyVector_lf(ztp1_dv);
   DestroyMatrix_lf(Ptp1_dm);

   return (loglh_timet);
}
#undef LOG2PI
#endif


/**
//----------------------------------------------------------------
//--  Tested OK, but has not use because tz_Refresh_z_T7P_T_in_kalfilms_1st_approx()
//--   cannot access to ElementV(smodel_ps->V[tp1],sti) or ElementV(smodel_ps->V[tbase0],sti)
//--   because no likelihood has been formed at all before this function is called.
//----------------------------------------------------------------
#define LOG2PI  (1.837877066409345e+000)   //log(2*pi)
//-----------------------------------------------------
//- Updating or refreshing all Kalman filter at time t for Markov-switching DSGE model.
//- WARNING: make sure to call the following functions
//      RunningGensys_const7varionly(lwzmodel_ps);
//      Refresh_kalfilms_*(lwzmodel_ps);   //Creates or refreshes kalfilmsinputs_ps at new parameter values.
//- before using tz_Refresh_z_T7P_T_in_kalfilms_1st_approx().
//-----------------------------------------------------
void tz_Refresh_z_T7P_T_in_kalfilms_1st_approx(struct TSkalfilmsinputs_tag *kalfilmsinputs_ps, struct TStateModel_tag *smodel_ps)
{
   double debug1;
   //--- Local variables
   int stp1i, stp1i_c, stp1i_v, sti, sti_c, sti_v, tbase0, tp1;
   //=== Accessible variables
   int ny = kalfilmsinputs_ps->ny;
   int nz = kalfilmsinputs_ps->nz;
   int nRc = kalfilmsinputs_ps->nRc;
   int nRv = kalfilmsinputs_ps->nRv;
   int T = kalfilmsinputs_ps->T;
   int indxIndRegimes = kalfilmsinputs_ps->indxIndRegimes;
   TSdvector z0_sdv;
   //+ input arguments.
   TSdmatrix *yt_dm = kalfilmsinputs_ps->yt_dm;         //ny-by-T.
   TSdmatrix *at_dm = kalfilmsinputs_ps->at_dm;         //ny-by-nRc.
   TSdcell *Ht_dc = kalfilmsinputs_ps->Ht_dc;           //ny-by-nz-by-nRc.
   TSdcell *Rt_dc = kalfilmsinputs_ps->Rt_dc;           //ny-by-ny-by-nRv.  Covariance matrix for the measurement equation.
   TSdcell *Gt_dc = kalfilmsinputs_ps->Gt_dc;           //nz-by-ny-by-nRv.  Cross-covariance.
   //
   TSdmatrix *bt_dm = kalfilmsinputs_ps->bt_dm;         //nz-by-nRc.
   TSdcell *Ft_dc = kalfilmsinputs_ps->Ft_dc;           //nz-by-nz-by-nRc.
   TSdcell *Vt_dc = kalfilmsinputs_ps->Vt_dc;           //nz-by-nz-by-nRv.  Covariance matrix for the state equation.
   //
   TSdmatrix *z0_dm = kalfilmsinputs_ps->z0_dm;        //nz-by-nRc*nRv or nz-by-nRv, depending on indxIndRegimes.
   TSdcell *P0_dc = kalfilmsinputs_ps->P0_dc;          //nz-by-nz-by-nRc*nRv or nz-by-nRv, depending on indxIndRegimes.
   //+ Output arguments.
   TSdcell *zt_tm1_dc = kalfilmsinputs_ps->zt_tm1_dc; //nz-by-nRc*nRv-by-T if indxIndRegimes==1, nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
   TSdfourth *Pt_tm1_d4 = kalfilmsinputs_ps->Pt_tm1_d4;     //nz-by-nz-by-nRc*nRv-T if indxIndRegimes==1, nz-by-nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
   //=== Work arguments.
   int nz2 = square(nz);
   TSdmatrix *Wnzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *Wnz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdmatrix *W2nz2bynz2_dm = CreateMatrix_lf(nz2,nz2);
   TSdvector *wP0_dv = CreateVector_lf(nz2);
   //+
   TSdvector yt_sdv, at_sdv, btp1_sdv;  //zt_tm1_sdv, ztp1_t_sdv,
   TSdvector *wny_dv = CreateVector_lf(ny);
   TSdmatrix *Wnzbyny_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *W2nzbynz_dm = CreateMatrix_lf(nz,nz);
   TSdmatrix *PHtran_tdata_dm = CreateMatrix_lf(nz,ny);
   TSdvector *etdata_dv = CreateVector_lf(ny);
   TSdmatrix *Dtdata_dm = CreateMatrix_lf(ny,ny);
   TSdmatrix *Kt_tdata0_dm = CreateMatrix_lf(nz,ny);
   TSdmatrix *Kt_tdata_dm = CreateMatrix_lf(nz,ny);
   //--- For eigenvalue decompositions
   int ki;
   int errflag;
   double eigmax;
   TSdzvector *evals_dzv = evals_dzv = CreateVector_dz(nz);
   TSdvector *evals_abs_dv = CreateVector_lf(nz); //Absolute eigenvalues.
   //--- For updating zt_tm1_dm and Pt_tm1.
   TSdvector *ztp1_t_dv = CreateVector_lf(z0_dm->nrows);
   TSdmatrix *Ptp1_t_dm = CreateMatrix_lf(nz, nz);
   TSdvector *ztp1_dv = CreateVector_lf(z0_dm->nrows);
   TSdmatrix *Ptp1_dm = CreateMatrix_lf(nz, nz);


   if (smodel_ps->sv->nstates != z0_dm->ncols)  fn_DisplayError("kalman.c/tz_Refresh_z_T7P_T_in_kalfilms_1st_approx():\n"
                     "  Make sure that the column dimension of z0_dm is the same as smodel_ps->sv->nstates");
   if (indxIndRegimes && (nRc>1) && (nRv>1))
      if (smodel_ps->sv->n_state_variables != 2)  fn_DisplayError("kalman.c/tz_Refresh_z_T7P_T_in_kalfilms_1st_approx():\n"
           " Number of state variables must be coincide with indxIndRegimes");


   z0_sdv.n = z0_dm->nrows;
   z0_sdv.flag = V_DEF;
   //
   at_sdv.n = yt_sdv.n = yt_dm->nrows;
   at_sdv.flag = yt_sdv.flag = V_DEF;
   btp1_sdv.n = bt_dm->nrows;
   btp1_sdv.flag = V_DEF;


   //======= Initial condition. =======
   for (sti=smodel_ps->sv->nstates-1; sti>=0;  sti--)
   {
      if (indxIndRegimes && (nRc==1))
      {
         sti_c = 0;
         sti_v = sti;
      }
      else if (indxIndRegimes && (nRv==1))
      {
         sti_c = sti;
         sti_v = 0;
      }
      else if (indxIndRegimes)
      {
         sti_c = smodel_ps->sv->Index[sti][0];
         sti_v = smodel_ps->sv->Index[sti][1];
      }
      else
      {
         sti_c = sti_v = sti;
      }


      if (!kalfilmsinputs_ps->indxIni)
      {
         InitializeDiagonalMatrix_lf(Wnzbynz_dm, 1.0);  //To be used for I(nz) -
         InitializeDiagonalMatrix_lf(Wnz2bynz2_dm, 1.0);  //To be used for I(nz2) -

         //=== Eigenanalysis to determine the roots to ensure boundedness.
         errflag = eigrgen(evals_dzv, (TSdzmatrix *)NULL, (TSdzmatrix *)NULL, Ft_dc->C[sti_c]);
         if (errflag)  fn_DisplayError("kalman.c/tz_Refresh_z_T7P_T_in_kalfilms_1st_approx(): eigen decomposition failed");
         for (ki=nz-1; ki>=0; ki--)  evals_abs_dv->v[ki] = sqrt(square(evals_dzv->real->v[ki]) + square(evals_dzv->imag->v[ki]));
         evals_abs_dv->flag = V_DEF;
         eigmax = MaxVector(evals_abs_dv);
         if (eigmax < (1.0+1.0e-14))
         {
            //--- Getting z0_dv: zt_tm1(:,1) = (eye(n_z)-F(:,:,1))\b(:,1);
            MatrixMinusMatrix(Wnzbynz_dm, Wnzbynz_dm, Ft_dc->C[sti_c]);
            z0_sdv.v = z0_dm->M + z0_sdv.n*sti;
            CopySubmatrix2vector(&z0_sdv, 0, bt_dm, 0, sti_c, bt_dm->nrows);
            bdivA_rgens(&z0_sdv, &z0_sdv, '\\', Wnzbynz_dm);
                      //Done with Wnzbynz_dm.
            //--- Getting P0_dm: Pt_tm1(:,:,1) = reshape((eye(n_z^2)-kron(F(:,:,1),F(:,:,1)))\V1(:),n_z,n_z);
            tz_kron(W2nz2bynz2_dm, Ft_dc->C[sti_c], Ft_dc->C[sti_c]);
            MatrixMinusMatrix(Wnz2bynz2_dm, Wnz2bynz2_dm, W2nz2bynz2_dm);
            CopySubmatrix2vector(wP0_dv, 0, Vt_dc->C[sti_v], 0, 0, nz2);
            bdivA_rgens(wP0_dv, wP0_dv, '\\', Wnz2bynz2_dm);
            CopySubvector2matrix_unr(P0_dc->C[sti], 0, 0, wP0_dv, 0, nz2);
                       //Done with all w*_dv and W*_dm.
         }
         else
         {
            printf("\n-----------------\n");
            printf("\nIn regime sti_c=%d and sti_v=%d and at time=%d\n", sti_c, sti_v, 0);
            fn_DisplayError("kalman.c/tz_Refresh_z_T7P_T_in_kalfilms_1st_approx(): the system is non-stationary solutions\n"
                          "    and the initial conditions must be supplied by, say, input arguments");
            fflush(stdout);
         }
      }
   }
   z0_dm->flag = M_GE;
   CopyMatrix0(zt_tm1_dc->C[0], z0_dm);  //At time t=0.
   CopyCell0(Pt_tm1_d4->F[0], P0_dc);                              //At time t=0.


//   fprintf(FPTR_DEBUG, "\n zt_tm1_dc->C[0]:\n");
//   WriteMatrix(FPTR_DEBUG, zt_tm1_dc->C[0], " %.16e ");
//   fprintf(FPTR_DEBUG, "\n");
//   fprintf(FPTR_DEBUG, "\n Pt_tm1_d4->F[0]->C[0]:\n");
//   WriteMatrix(FPTR_DEBUG, Pt_tm1_d4->F[0]->C[0], " %.16e ");


   //============== Updating zt_tm1 and Pt_tm1. ==================
   for (tbase0=0; tbase0<T; tbase0++ )
   {
      //tdata = tbase0 is base-0 timing.
      tp1 = tbase0 + 1;  //Next period.

      for (stp1i=smodel_ps->sv->nstates-1; stp1i>=0; stp1i--)
      {
         if (indxIndRegimes && (nRc==1))
         {
            stp1i_c = 0;
            stp1i_v = stp1i;
         }
         else if (indxIndRegimes && (nRv==1))
         {
            stp1i_c = stp1i;
            stp1i_v = 0;
         }
         else if (indxIndRegimes)
         {
            stp1i_c = smodel_ps->sv->Index[stp1i][0];
            stp1i_v = smodel_ps->sv->Index[stp1i][1];
         }
         else
         {
            stp1i_c = stp1i_v = stp1i;
         }


         InitializeConstantVector_lf(ztp1_dv, 0.0);  //To be summed over sti.
         InitializeConstantMatrix_lf(Ptp1_dm, 0.0);  //To be summed over sti.
         for (sti=smodel_ps->sv->nstates-1; sti>=0;  sti--)
         {
            if (indxIndRegimes && (nRc==1))
            {
               sti_c = 0;
               sti_v = sti;
            }
            else if (indxIndRegimes && (nRv==1))
            {
               sti_c = sti;
               sti_v = 0;
            }
            else if (indxIndRegimes)
            {
               sti_c = smodel_ps->sv->Index[sti][0];
               sti_v = smodel_ps->sv->Index[sti][1];
            }
            else
            {
               sti_c = sti_v = sti;
            }

            //--- Setup.
            MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_d4->F[tbase0]->C[sti], Ht_dc->C[sti_c], 1.0, 0.0, 'N', 'T');

            //--- Data.
            //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata where tdata = tbase0 = inpt-1.
            yt_sdv.v = yt_dm->M + tbase0*yt_dm->nrows;
            at_sdv.v = at_dm->M + sti_c*at_dm->nrows;  //sti_c: coefficient regime at time tbase0.
            z0_sdv.v = zt_tm1_dc->C[tbase0]->M + z0_sdv.n*sti;  //sti: regime at time tbase0.
            VectorMinusVector(etdata_dv, &yt_sdv, &at_sdv);
            MatrixTimesVector(etdata_dv, Ht_dc->C[sti_c], &z0_sdv, -1.0, 1.0, 'N');
            //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tbase0);
            CopyMatrix0(Dtdata_dm, Rt_dc->C[sti_v]);
            MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[sti_c], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
                                        //Done with z0_sdv.v.


            //=== Updating for next period by integrating out sti..
            if (tp1<T)
            {
               //Updating only up to tbase0=T-2.  The values at tp1=T or tbase0=T-1 will not be used in the likelihood function.

               //--- Ktp1_t = (F_tp1*PHtran_t+G(:,:,t))/Dt;
               //--- Kt_tdata = (Ft*PHtran_tdata+G(:,:,tdata))/Dtdata where t=tp1 and tdata=t.
               CopyMatrix0(Kt_tdata0_dm, Gt_dc->C[sti_v]);
               MatrixTimesMatrix(Kt_tdata0_dm, Ft_dc->C[stp1i_c], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');
               BdivA_rrect(Kt_tdata_dm, Kt_tdata0_dm, '/', Dtdata_dm);
               //+ zt_tm1(:,t) = b(:,t) + Ft*zt_tm1(:,tdata) + Kt_tdata*etdata where t=tp1 and tm1=t.
               MatrixTimesVector(ztp1_t_dv, Ft_dc->C[stp1i_c], &z0_sdv, 1.0, 0.0, 'N');
               MatrixTimesVector(ztp1_t_dv, Kt_tdata_dm, etdata_dv, 1.0, 1.0, 'N');
               btp1_sdv.v = bt_dm->M + stp1i_c*btp1_sdv.n;
               VectorPlusMinusVectorUpdate(ztp1_t_dv, &btp1_sdv, 1.0);
               //+ Pt_tm1(:,:,t) = Ft*Ptdata*Fttran - Kt_tdata*Dtdata*Kt_tdatatran + V(:,:,t);
               CopyMatrix0(Ptp1_t_dm, Vt_dc->C[stp1i]);
               MatrixTimesMatrix(Wnzbyny_dm, Kt_tdata_dm, Dtdata_dm, 1.0, 0.0, 'N', 'N');
               MatrixTimesMatrix(Wnzbynz_dm, Wnzbyny_dm, Kt_tdata_dm, 1.0, 0.0, 'N', 'T');
               MatrixPlusMinusMatrixUpdate(Ptp1_t_dm, Wnzbynz_dm, -1.0);
                                     //Done with all W*_dm.
               MatrixTimesMatrix(Wnzbynz_dm, Ft_dc->C[stp1i_c], Pt_tm1_d4->F[tbase0]->C[sti], 1.0, 0.0, 'N', 'N');
               MatrixTimesMatrix(W2nzbynz_dm, Wnzbynz_dm, Ft_dc->C[stp1i_c], 1.0, 0.0, 'N', 'T');
               MatrixPlusMatrixUpdate(Ptp1_t_dm, W2nzbynz_dm);
                                     //Done with all W*_dm.

               //--- Integrating out the state at tbase0 using P(s_t|Y_t, theta) = ElementV(smodel_ps->V[t+1],s_{t+1}_i).
               //---   Note because the data in DW code (ElementV) is base-1, t+1 is actually tbase0.
               debug1 = ElementV(smodel_ps->V[tp1],sti);  //?????? Debug.
               //ScalarTimesVectorUpdate(ztp1_dv, ElementV(smodel_ps->V[tp1],sti), ztp1_t_dv);
               //ScalarTimesMatrix(Ptp1_dm, ElementV(smodel_ps->V[tp1],sti), Ptp1_t_dm, 1.0);
               ScalarTimesVectorUpdate(ztp1_dv, 0.5, ztp1_t_dv);
               ScalarTimesMatrix(Ptp1_dm, 0.5, Ptp1_t_dm, 1.0);
               Ptp1_dm->flag = M_GE | M_SU | M_SL;
                                         //Done with ztp1_t_dv and Ptp1_t_dm.
            }
         }
         //--- Filling zt_tm1 and Pt_tm1 for next period
         if (tp1<T)
         {
            z0_sdv.v = zt_tm1_dc->C[tp1]->M + z0_sdv.n*stp1i;  //stp1i: regime at time tp1.
            CopyVector0(&z0_sdv, ztp1_dv);
            CopyMatrix0(Pt_tm1_d4->F[tp1]->C[stp1i], Ptp1_dm);  //stp1i: regime at time tp1.
                                           //Done with ztp1_dv, z0_sdv, Ptp1_dm.
         }
      }
      if (tp1<T)
         zt_tm1_dc->C[tp1]->flag = M_GE;

//      fprintf(FPTR_DEBUG, "\n &yt_sdv:\n");
//      WriteMatrix(FPTR_DEBUG, &yt_sdv, " %.16e ");
//      fprintf(FPTR_DEBUG, "\n zt_tm1_dc->C[tp1]:\n");
//      WriteMatrix(FPTR_DEBUG, zt_tm1_dc->C[tp1], " %.16e ");
//      fprintf(FPTR_DEBUG, "\n");
//      fprintf(FPTR_DEBUG, "\n Pt_tm1_d4->F[tp1]->C[0]:\n");
//      WriteMatrix(FPTR_DEBUG, Pt_tm1_d4->F[tp1]->C[0], " %.16e ");
//      fprintf(FPTR_DEBUG, "\n");
//      fflush(FPTR_DEBUG);


   }

   //===
   DestroyVector_dz(evals_dzv);
   DestroyVector_lf(evals_abs_dv);
   DestroyMatrix_lf(Wnzbynz_dm);
   DestroyMatrix_lf(Wnz2bynz2_dm);
   DestroyMatrix_lf(W2nz2bynz2_dm);
   DestroyVector_lf(wP0_dv);
   //
   DestroyVector_lf(wny_dv);
   DestroyMatrix_lf(Wnzbyny_dm);
   DestroyMatrix_lf(W2nzbynz_dm);
   DestroyMatrix_lf(PHtran_tdata_dm);
   DestroyVector_lf(etdata_dv);
   DestroyMatrix_lf(Dtdata_dm);
   DestroyMatrix_lf(Kt_tdata0_dm);
   DestroyMatrix_lf(Kt_tdata_dm);
   //
   DestroyVector_lf(ztp1_t_dv);
   DestroyMatrix_lf(Ptp1_t_dm);
   DestroyVector_lf(ztp1_dv);
   DestroyMatrix_lf(Ptp1_dm);
}
//-----------------------------------------------------
//- Kalman filter at time t for Markov-switching DSGE model.
//- WARNING: make sure to call the following functions
//      (1) RunningGensys_const7varionly(lwzmodel_ps);
//      (2) Refresh_kalfilms_*(lwzmodel_ps);   //Creates or refreshes kalfilmsinputs_ps at new parameter values.
//      (3) tz_Refresh_z_T7P_T_in_kalfilms_1st_approx();
//- before using kalfilms_timet_1st_approx().
//-----------------------------------------------------
double tz_kalfilms_timet_1st_approx(int st, int inpt, struct TSkalfilmsinputs_tag *kalfilmsinputs_ps, struct TStateModel_tag *smodel_ps)
{
   //st, st_c, and st_v: base-0: deals with the cross-section values at time t where
   //      st is a grand regime, st_c is an encoded coefficient regime, and st_c is an encoded volatility regime.
   //inpt: base-1 in the sense that inpt>=1 to deal with the time series situation where S_T is (T+1)-by-1 and Y_T is T+nlags_max-by-1.
   //      The 1st element for S_T is S_T[1] while S_T[0] is s_0.  The same for (T+1)-by-1 gbeta_dv and nlcoefs-by-(T+1) galpha_dm.
   //      The 1st element for Y_T, however, is Y_T[nlags_max+1-1].
   //See (42.3) on p.42 in the SWZII NOTES.


   //--- Local variables
   int st_c, st_v, tbase0;
   double loglh_timet;
   //--- Accessible variables
   int ny = kalfilmsinputs_ps->ny;
   int nz = kalfilmsinputs_ps->nz;
   int nRc = kalfilmsinputs_ps->nRc;
   int nRv = kalfilmsinputs_ps->nRv;
   int indxIndRegimes = kalfilmsinputs_ps->indxIndRegimes;
   TSdvector z0_sdv;
   //+ input arguments.
   TSdmatrix *yt_dm = kalfilmsinputs_ps->yt_dm;         //ny-by-T.
   TSdmatrix *at_dm = kalfilmsinputs_ps->at_dm;         //ny-by-nRc.
   TSdcell *Ht_dc = kalfilmsinputs_ps->Ht_dc;           //ny-by-nz-by-nRc.
   TSdcell *Rt_dc = kalfilmsinputs_ps->Rt_dc;           //ny-by-ny-by-nRv.  Covariance matrix for the measurement equation.
   //+ Output arguments.
   TSdcell *zt_tm1_dc = kalfilmsinputs_ps->zt_tm1_dc; //nz-by-nRc*nRv-by-T if indxIndRegimes==1, nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
   TSdfourth *Pt_tm1_d4 = kalfilmsinputs_ps->Pt_tm1_d4;     //nz-by-nz-by-nRc*nRv-T if indxIndRegimes==1, nz-by-nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
   //=== Work arguments.
   TSdvector yt_sdv, at_sdv;
   TSdvector *wny_dv = CreateVector_lf(ny);
   TSdmatrix *PHtran_tdata_dm = CreateMatrix_lf(nz,ny);
   TSdvector *etdata_dv = CreateVector_lf(ny);
   TSdmatrix *Dtdata_dm = CreateMatrix_lf(ny,ny);


   if (smodel_ps->sv->nstates != zt_tm1_dc->C[0]->ncols)  fn_DisplayError("kalman.c/kalfilms_timet_1st_approx():\n"
                     "  Make sure that the column dimension of zt_tm1_dc->C is the same as smodel_ps->sv->nstates");

   tbase0 = inpt - 1;  //base-0 time t.

   if (indxIndRegimes && (nRc==1))
   {
      st_c = 0;
      st_v = st;
   }
   else if (indxIndRegimes && (nRv==1))
   {
      st_c = st;
      st_v = 0;
   }
   else if (indxIndRegimes)
   {
      if (smodel_ps->sv->n_state_variables != 2)  fn_DisplayError("kalman.c/kalfilms_timet_1st_approx():\n"
                      "  Number of state variables must be coincide with indxIndRegimes");
      st_c = smodel_ps->sv->Index[st][0];
      st_v = smodel_ps->sv->Index[st][1];
   }
   else
   {
      st_c = st_v = st;
   }


   z0_sdv.n = zt_tm1_dc->C[0]->nrows;
   z0_sdv.flag = V_DEF;
   //
   at_sdv.n = yt_sdv.n = yt_dm->nrows;
   at_sdv.flag = yt_sdv.flag = V_DEF;

   //====== Computing the conditional LH at time t. ======
   //--- Setup.
   MatrixTimesMatrix(PHtran_tdata_dm, Pt_tm1_d4->F[tbase0]->C[st], Ht_dc->C[st_c], 1.0, 0.0, 'N', 'T');

   //--- Data.
   //- etdata = Y_T(:,tdata) - a(:,tdata) - Htdata*ztdata where tdata = tbase0 = inpt-1.
   yt_sdv.v = yt_dm->M + tbase0*yt_dm->nrows;
   at_sdv.v = at_dm->M + st_c*at_dm->nrows;    //st_c: coefficient regime at time tbase0.
   z0_sdv.v = zt_tm1_dc->C[tbase0]->M + z0_sdv.n*st;  //st: regime at time tbase0 for zt_tm1.
   VectorMinusVector(etdata_dv, &yt_sdv, &at_sdv);
   MatrixTimesVector(etdata_dv, Ht_dc->C[st_c], &z0_sdv, -1.0, 1.0, 'N');
   //+   Dtdata = Htdata*PHtran_tdata + R(:,:,tbase0);
   CopyMatrix0(Dtdata_dm, Rt_dc->C[st_v]);
   MatrixTimesMatrix(Dtdata_dm, Ht_dc->C[st_c], PHtran_tdata_dm, 1.0, 1.0, 'N', 'N');

   //--- Forming the log conditional likelihood at t.
   bdivA_rgens(wny_dv, etdata_dv, '/', Dtdata_dm);
   loglh_timet = -(0.5*ny)*LOG2PI - 0.5*logdeterminant(Dtdata_dm) - 0.5*VectorDotVector(wny_dv, etdata_dv);
                         //Done with all w*_dv.


   //===
   DestroyVector_lf(wny_dv);
   DestroyMatrix_lf(PHtran_tdata_dm);
   DestroyVector_lf(etdata_dv);
   DestroyMatrix_lf(Dtdata_dm);

   return (loglh_timet);
}
#undef LOG2PI
/**/




