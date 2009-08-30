#ifndef __KALMAN_H__
   #define __KALMAN_H__

   #include "tzmatlab.h"
   #include "mathlib.h"
   #include "switch.h"
   #include "fn_filesetup.h"  //Used to call WriteMatrix(FPTR_DEBUG,....).


   typedef struct TSkalcvfurw_tag {
           //urw: univariate random walk kalman filter.  Desigend specially for the 2006 AER SWZ paper.

           //=== Input arguments.
           int indx_tvsigmasq;  //0: constant siqmasq in Kalman updating (default);
                                //1: Keyensian (project-specific) type of time-varying sigmasq in Kalman updating;  See pp.37 and 37a in SWZ Learning NOTES;
                                //2: project-specific type;
                                //3: another project-specific type.
           double sigmasq;  //Variance for the residual eps(t) of the measurement equation.
           int fss;   //T: effective sample size (excluding lags).
           int kx;    //dimension for x(t).
           TSdmatrix *V_dm;   //kx-by-kx.  Covariance (symmetric and positive definite) matrix for the residual eta(t) of the transition equation.
           TSdvector *ylhtran_dv;   //1-by-T of y(t).  The term lh means lelf hand side and tran means transpose.
           TSdmatrix *Xrhtran_dm;   //kx-by-T of x(t).  The term rh means right hand side and tran means transpose.
           TSdvector *z10_dv;    //kx-by-1.   Initial condition for prediction: z_{1|0}.
           TSdmatrix *P10_dm;    //kx-by-kx symmetric matrix.  Initial condition for the variance of the prediction: P_{1|0}.

           //=== Output arguments.
           TSdvector *zupdate_dv;  //kx-by-1.  z_{T+1|T}.
           TSdmatrix *Zpredtran_dm;  //kx-by-T matrix of one-step predicted values of z(t).  [z_{2|1}, ..., z_{t+1|t}, ..., z_{T+1|T}].
                                 //Set to NULL (no output) if storeZ = 0;
           TSdcell *Ppred_dc;   //T cells and kx-by-kx symmetric and positive definite matrix for each cell.  Mean square errors of predicted state variables.
                                //{P_{2|1}, ..., P{t+1|t}, ..., P{T+1|T}.  Set to NULL (no output if storeV = 0).
           TSdvector *ylhtranpred_dv;   //1-by-T one-step prediction of y(t) or ylhtran_dv.  Added 03/17/05.

           //=== Function itself.
           void (*learning_fnc)(struct TSkalcvfurw_tag *, void *);
   } TSkalcvfurw;   //urw: univariate random walk.
   //
   typedef void TFlearninguni(struct TSkalcvfurw_tag *, void *);  //For linear rational expectations models.


   //=== Better version is TSkalfilmsinputs_1stapp_tag.  Kalman filter for constant or known-time-varying DSGE models.
   typedef struct TSkalfiltv_tag
   {
      //General (known-time-varying) Kalman filter for DSGE models.
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
      //        z0 is an n_z-by-1 vector of initial condition when indxIni=1. (Not used if indxIni=0.)
      //        P0 is an n_z-by-n_z matrix of initial condition when indxIni=1. (Not used if indxIni=0.)
      //
      //   Outputs are as follows:
      //      loglh is a value of the log likelihood function of the state-space model
      //                                under the assumption that errors are multivariate Gaussian.
      //      zt_tm1 is an n_z-by-T matrices of one-step predicted state vectors with z0_0m1 as a initial condition
      //                         and with z_{t+1|t} as the last element.  Thus, we can use it as a base-1 vector.
      //      Pt_tm1 is an n_z-by-n_z-by-T 3-D of covariance matrices of zt_tm1 with P0_0m1 as a initial condition
      //                         and with P_{t+1|t} as the last element.  Thus, we can use it as a base-1 cell.
      //
      //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:
      //             z0_0m1 = (I-F(:,:,1))\b(:,1)
      //        vec(P0_0m1) = (I-kron(F(:,:,1),F(:,:,1)))\vec(V(:,:,1))
      //   Note that all eigenvalues of the matrix F(:,:,1) are inside the unit circle when the state-space model is bounded (stationary).
      //
      //   March 2007, written by Tao Zha
      //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.

      //=== Input arguments.
      int ny;  //number of observables.
      int nz;  //number of state variables.
      int T;   //sample size.
      int indxIni;   //1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
                     //0: using the unconditional mean for any given regime at time 0. (Default value)
      TSdmatrix *yt_dm;   //ny-by-T.
      TSdmatrix *at_dm;   //ny-by-T.
      TSdcell *Ht_dc;   //ny-by-nz-by-T.
      TSdcell *Rt_dc;   //ny-by-ny-by-T.  Covariance matrix for the measurement equation.
      TSdcell *Gt_dc;   //nz-by-ny-by-T.  Cross-covariance.
      //
      TSdmatrix *bt_dm;   //nz-by-T.
      TSdcell *Ft_dc;   //nz-by-nz-by-T.
      TSdcell *Vt_dc;   //nz-by-nz-by-T.  Covariance matrix for the state equation.
      //
      TSdvector *z0_dv;  //nz-by-1;
      TSdmatrix *P0_dm;   //nz-by-nz.

      //=== Output arguments.
      double loglh;  //log likelihood.
      TSdmatrix *zt_tm1_dm;  //nz-by-T.
      TSdcell *Pt_tm1_dc;   //nz-by-nz-T.
   } TSkalfiltv;



   //=== Inputs for filter for Markov-switching DSGE models at any time t.
   typedef struct TSkalfilmsinputs_1stapp_tag
   {
      //Inputs Markov-switching Kalman filter for DSGE models (conditional on all the regimes at time t).
      //  It computes a sequence of one-step predictions and their covariance matrices, and the log likelihood.
      //  The function uses a forward recursion algorithm.  See also the Matlab function fn_kalfil_tv.m
      //
      //   State space model is defined as follows:
      //       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)
      //       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)
      //     where a(t), H(t), b(t), and F(t) depend on the grand regime s_t that follows a Markov-chain process
      //                                                                          and is taken as given.
      //
      //   Inputs at time t are as follows where nst is number of grand regimes (including lagged regime
      //                                           and coefficients and shock variances):
      //      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].
      //        a is an n_y-by-nst matrix of Markov-switching input vectors in the measurement equation.
      //        H is an n_y-by-n_z-by-nst 3-D of Markov-switching matrices in the measurement equation.
      //        R is an n_y-by-n_y-by-nst 3-D of Markov-switching covariance matrices for the error in the measurement equation.
      //        G is an n_z-by-n_y-by-nst 3-D of Markov-switching E(eta_t * eps_t').
      //        ------
      //        b is an n_z-by-nst matrix of Markov-switching input vectors in the state equation with b(:,st) as an initial condition.
      //                (alternatively, with the ergodic weighted b(:,st) as an initial condition).
      //        F is an n_z-by-n_z-by-nst 3-D of Markov-switching transition matrices in the state equation with F(:,:,st)
      //                as an initial condition (alternatively, with the ergodic weighted F(:,:,st) as an initial condition).
      //        V is an n_z-by-n_z-by-nRv 3-D of Markov-switching covariance matrices for the error in the state equation
      //                with V(:,:,st) as an initial condition (alternatively, with the ergodic weighted V(:,:,st) as an initial condition).
      //        ------
      //        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
      //                 0: using the unconditional mean for any given regime at time 0.
      //        z0 is an n_z-by-nst matrix of initial condition (Not used if indxIni=0).
      //        P0 is an n_z-by-n_z-by-nst 3-D of initial condition (Not used if indxIni=0).
      //
      //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:
      //             z0_0m1 = (I-F(:,:,st))\b(:,st)
      //        vec(P0_0m1) = (I-kron(F(:,:,st),F(:,:,st)))\vec(V(:,:,st))
      //   Note that all eigenvalues of the matrix F(:,:,st) are inside the unit circle when the state-space model is bounded (stationary).
      //
      //   November 2007, written by Tao Zha.  Revised, April 2008.
      //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.

      //=== Input arguments.
      int ny;  //number of observables.
      int nz;  //number of state variables.
      int nst; //number of grand composite regimes (current and past regimes, coefficient and volatility regimes).
      int T;   //sample size.
      int indxIni;   //1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0,
                     //0: using the unconditional momnets for any given regime at time 0 (default when indxDiffuse = 0).
      int indxDiffuse;  //1: using the diffuse condition for z_{1|0} and P_{1|0} (default option), according to Koopman and Durbin, "Filtering and Smoothing of State Vector for Diffuse State-Space Models," J. of Time Series Analysis, Vol 24(1), pp.85-99.
                        //0: using the unconditional moments.
      double DiffuseScale; //A large (infinity) number when indxDiffuse = 1.
      int ztm1_track; //t-1 = -1:      no initial conditions z_{1|0} and P_{1|0} has been computed yet, but will be using InitializeKalman_z10_P10(),
                      //t-1 >= 0:T-1:  z_{t|t-1} and P_{t|t-1} are updated up to t-1.
      int dtm1_track; //t-1 = -1:      no etdata_dc->C[0] or Dtdata_d4->F[0] has been computed yet.
                      //t-1 >= 0:T-1:  etdata_dc->C[t-1] and Dtdata_d4->F[t-1] are updated up to t-1.

      TSdmatrix *yt_dm;   //ny-by-T.
      TSdmatrix *at_dm;   //ny-by-nst.
      TSdcell *Ht_dc;     //ny-by-nz-by-nst.
      TSdcell *Rt_dc;     //ny-by-ny-by-nst.  Covariance matrix for the measurement equation.
      TSdcell *Gt_dc;     //nz-by-ny-by-nst.  Cross-covariance.
      //
      TSdmatrix *bt_dm;   //nz-by-nst.
      TSdcell *Ft_dc;     //nz-by-nz-by-nst.
      TSdcell *Vt_dc;     //nz-by-nz-by-nst.  Covariance matrix for the state equation.
      //
      TSdmatrix *z0_0_dm; //nz-by-nst. z_{0|0}.
      TSdmatrix *z0_dm;  //nz-by-nst. z_{1|0}.
      TSdcell *P0_dc;    //nz-by-nz-by-nst. P_{1|0}


      //=== Output arguments only used for 1st order approximation to zt and Pt depending on infinite past regimes.
      TSdcell *zt_tm1_dc;   //nz-by-nst-by-(T+1), where z_{1|0} is an initial condition (1st element with t-1=0 or t=1 for base-1) and
                            //  the terminal condition z_{T+1|T} using Updatekalfilms_1stapp(T, ...) is not computed
                            //  when the likelihood logTimetCondLH_kalfilms_1stapp() is computed.  Thus, z_{T+1|T}
                            //  has not legal value computed in most applications unless in out-of-sample forecasting problems.
      TSdfourth *Pt_tm1_d4; //nz-by-nz-by-nst-by-(T+1), where P_{1|0} is an initial condition (1st element with t-1=0) and
                            //  the terminal condition P_{T+1|T} using Updatekalfilms_1stapp(T, ...) is not computed
                            //  when the likelihood logTimetCondLH_kalfilms_1stapp() is computed.  Thus, P_{T+1|T}
                            //  has not legal value computed in most applications unless in out-of-sample forecasting problems.
      //+ Will be save for updating likelihood and Kalman filter Updatekalfilms_1stapp(), so save time to recompute these objects again.
      TSdfourth *PHtran_tdata_d4;  //nz-by-ny-by-nst-T, P_{t|t-1}*H_t'.  Saved only for updating Kalman filter Updatekalfilms_1stapp().
      TSdcell *etdata_dc; //ny-by-nst-by-T (with base-0 T), forecast errors e_t in the likelihood.
      TSdcell *yt_tm1_dc; //ny-by-nst-by-T, one-step forecast y_{t|t-1} for t=0 to T-1 (base-0). Used to back out structural shocks.
      TSdfourth *Dtdata_d4; //ny-by-ny-nst-by-T, forecast covariance D_t in the likelihood.  Saved for updating Kalman filter Updatekalfilms_1stapp().
   } TSkalfilmsinputs_1stapp;


   //=== OLD Code: Inputs for filter for Markov-switching DSGE models at any time t.
   typedef struct TSkalfilmsinputs_tag
   {
      //Inputs Markov-switching Kalman filter for DSGE models (conditional on all the regimes at time t).
      //  It computes a sequence of one-step predictions and their covariance matrices, and the log likelihood.
      //  The function uses a forward recursion algorithm.  See also the Matlab function fn_kalfil_tv.m
      //
      //   State space model is defined as follows:
      //       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)
      //       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)
      //     where a(t), H(t), b(t), and F(t) depend on s_t that follows a Markov-chain process and are taken as given.
      //
      //   Inputs at time t are as follows where nRc is number of regimes for coefficients
      //                                         nRv is number of regimes for volatility (shock variances):
      //      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].
      //        a is an n_y-by-nRc matrix of Markov-switching input vectors in the measurement equation.
      //        H is an n_y-by-n_z-by-nRc 3-D of Markov-switching matrices in the measurement equation.
      //        R is an n_y-by-n_y-by-nRv 3-D of Markov-switching covariance matrices for the error in the measurement equation.
      //        G is an n_z-by-n_y-by-nRv 3-D of Markov-switching E(eta_t * eps_t').
      //        ------
      //        b is an n_z-by-nRc matrix of Markov-switching input vectors in the state equation with b(:,1) as an initial condition.
      //        F is an n_z-by-n_z-by-nRc 3-D of Markov-switching transition matrices in the state equation with F(:,:,1) as an initial condition.
      //        V is an n_z-by-n_z-by-nRv 3-D of Markov-switching covariance matrices for the error in the state equation with V(:,:,1) as an initial condition.
      //        ------
      //        indxIndRegimes:  1: coefficient regime and volatility regime are independent; 0: these two regimes are synchronized completely.
      //        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
      //                 0: using the unconditional mean for any given regime at time 0.
      //        z0 is an n_z-by-nRc*nRv matrix of initial condition when indxIni=1 and indxIndRegimes=1. (Not used if indxIni=0.)
      //        z0 is an n_z-by-nRv matrix of initial condition when indxIni=1 and indxIndRegimes=0. (Not used if indxIni=0.)
      //        P0 is an n_z-by-n_z-by-nRc*nRv 3-D of initial condition when indxIni=1 and indxIndRegimes=1. (Not used if indxIni=0.)
      //        P0 is an n_z-by-n_z-by-nRv 3-D of initial condition when indxIni=1 and indxIndRegimes=0. (Not used if indxIni=0.)
      //
      //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:
      //             z0_0m1 = (I-F(:,:,1))\b(:,1)
      //        vec(P0_0m1) = (I-kron(F(:,:,1),F(:,:,1)))\vec(V(:,:,1))
      //   Note that all eigenvalues of the matrix F(:,:,1) are inside the unit circle when the state-space model is bounded (stationary).
      //
      //   November 2007, written by Tao Zha
      //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.

      //=== Input arguments.
      int ny;  //number of observables.
      int nz;  //number of state variables.
      int nRc; //number of composite regimes (current and past regimes) for coefficients.
      int nRstc;  //number of coefficient regimes at time t.
      int nRv; //number of regimes for volatility (shock variances).
      int indxIndRegimes; //1: coefficient regime and volatility regime are independent; 0: these two regimes are synchronized completely.
      int T;   //sample size.
      int indxIni;   //1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;
                     //0: using the unconditional mean for any given regime at time 0. (Default value)
      TSdmatrix *yt_dm;   //ny-by-T.
      TSdmatrix *at_dm;   //ny-by-nRc.
      TSdcell *Ht_dc;     //ny-by-nz-by-nRc.
      TSdcell *Rt_dc;     //ny-by-ny-by-nRv.  Covariance matrix for the measurement equation.
      TSdcell *Gt_dc;     //nz-by-ny-by-nRv.  Cross-covariance.
      //
      TSdmatrix *bt_dm;   //nz-by-nRc.
      TSdcell *Ft_dc;     //nz-by-nz-by-nRc.
      TSdcell *Vt_dc;     //nz-by-nz-by-nRv.  Covariance matrix for the state equation.
      //
      TSdmatrix *z0_dm;  //nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nRv if indxIndRegimes == 0.
      TSdcell *P0_dc;    //nz-by-nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nz-by-nRv if indxIndRegimes == 0.


      //=== Output arguments only used for 1st order approximation to zt and Pt depending on infinite past regimes.
      TSdcell *zt_tm1_dc;      //nz-by-nRc*nRv-by-T if indxIndRegimes==1, nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
      TSdfourth *Pt_tm1_d4;    //nz-by-nz-by-nRc*nRv-T if indxIndRegimes==1, nz-by-nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.
   } TSkalfilmsinputs;




   //--- Functions for univariate random walk kalman filter.
   TSkalcvfurw *CreateTSkalcvfurw(TFlearninguni *func, int T, int k, int tv);  //, int storeZ, int storeV);
   TSkalcvfurw *DestroyTSkalcvfurw(TSkalcvfurw *kalcvfurw_ps);
   void kalcvf_urw(TSkalcvfurw *kalcvfurw_ps, void *dummy_ps);

   //--- New Code: Functions for Markov-switching Kalman filter.
   struct TSkalfilmsinputs_1stapp_tag *CreateTSkalfilmsinputs_1stapp(int ny, int nz, int nst, int T);
   struct TSkalfilmsinputs_1stapp_tag *DestroyTSkalfilmsinputs_1stapp(struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps);
   int InitializeKalman_z10_P10(struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, TSdmatrix *z10_dm, TSdcell *P10_dc);
   double logTimetCondLH_kalfilms_1stapp(int st, int inpt, struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, struct TStateModel_tag *smodel_ps);



   //--- OLD Code: Functions for general constant Kalman filter.
   struct TSkalfiltv_tag *CreateTSkalfiltv(int ny, int nz, int T);
   struct TSkalfiltv_tag *DestroyTSkalfiltv(struct TSkalfiltv_tag *kalfiltv_ps);
   //Used to test tz_logTimetCondLH_kalfiltv(). (Done April 08).  double tz_kalfiltv(struct TSkalfiltv_tag *kalfiltv_ps);
   double tz_logTimetCondLH_kalfiltv(int st, int inpt, struct TSkalfiltv_tag *kalfiltv_ps);

   //--- OLD Code: Functions for Markov-switching Kalman filter.
   struct TSkalfilmsinputs_tag *CreateTSkalfilmsinputs(int ny, int nz, int nRc, int nRstc, int nRv, int indxIndRegimes, int T);
   struct TSkalfilmsinputs_tag *DestroyTSkalfilmsinputs(struct TSkalfilmsinputs_tag *kalfilmsinputs_ps);
   double tz_logTimetCondLH_kalfilms_1st_approx(int st, int inpt, struct TSkalfilmsinputs_tag *kalfilmsinputs_ps, struct TStateModel_tag *smodel_ps);
            //IMPORTANT NOTE: in the Markov-switching input file datainp_markov*.prn, it MUST be that
            //                                the coefficient regime is the 1st state variable, and
            //                                the volatility regime is the 2nd state variable.
#endif

