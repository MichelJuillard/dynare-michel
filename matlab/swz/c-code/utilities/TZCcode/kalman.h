#ifndef __KALMAN_H__
   #define __KALMAN_H__

   #include "tzmatlab.h"
   #include "mathlib.h"
   #include "switch.h"
   #include "fn_filesetup.h"   /*  Used to call WriteMatrix(FPTR_DEBUG,....).   ansi-c*/


   typedef struct TSkalcvfurw_tag {
/*             //urw: univariate random walk kalman filter.  Desigend specially for the 2006 AER SWZ paper.   ansi-c*/

/*             //=== Input arguments.   ansi-c*/
           int indx_tvsigmasq;   /*  0: constant siqmasq in Kalman updating (default);   ansi-c*/
/*                                  //1: Keyensian (project-specific) type of time-varying sigmasq in Kalman updating;  See pp.37 and 37a in SWZ Learning NOTES;   ansi-c*/
/*                                  //2: project-specific type;   ansi-c*/
/*                                  //3: another project-specific type.   ansi-c*/
           double sigmasq;   /*  Variance for the residual eps(t) of the measurement equation.   ansi-c*/
           int fss;    /*  T: effective sample size (excluding lags).   ansi-c*/
           int kx;     /*  dimension for x(t).   ansi-c*/
           TSdmatrix *V_dm;    /*  kx-by-kx.  Covariance (symmetric and positive definite) matrix for the residual eta(t) of the transition equation.   ansi-c*/
           TSdvector *ylhtran_dv;    /*  1-by-T of y(t).  The term lh means lelf hand side and tran means transpose.   ansi-c*/
           TSdmatrix *Xrhtran_dm;    /*  kx-by-T of x(t).  The term rh means right hand side and tran means transpose.   ansi-c*/
           TSdvector *z10_dv;     /*  kx-by-1.   Initial condition for prediction: z_{1|0}.   ansi-c*/
           TSdmatrix *P10_dm;     /*  kx-by-kx symmetric matrix.  Initial condition for the variance of the prediction: P_{1|0}.   ansi-c*/

/*             //=== Output arguments.   ansi-c*/
           TSdvector *zupdate_dv;   /*  kx-by-1.  z_{T+1|T}.   ansi-c*/
           TSdmatrix *Zpredtran_dm;   /*  kx-by-T matrix of one-step predicted values of z(t).  [z_{2|1}, ..., z_{t+1|t}, ..., z_{T+1|T}].   ansi-c*/
/*                                   //Set to NULL (no output) if storeZ = 0;   ansi-c*/
           TSdcell *Ppred_dc;    /*  T cells and kx-by-kx symmetric and positive definite matrix for each cell.  Mean square errors of predicted state variables.   ansi-c*/
/*                                  //{P_{2|1}, ..., P{t+1|t}, ..., P{T+1|T}.  Set to NULL (no output if storeV = 0).   ansi-c*/
           TSdvector *ylhtranpred_dv;    /*  1-by-T one-step prediction of y(t) or ylhtran_dv.  Added 03/17/05.   ansi-c*/

/*             //=== Function itself.   ansi-c*/
           void (*learning_fnc)(struct TSkalcvfurw_tag *, void *);
   } TSkalcvfurw;    /*  urw: univariate random walk.   ansi-c*/
/*     //   ansi-c*/
   typedef void TFlearninguni(struct TSkalcvfurw_tag *, void *);   /*  For linear rational expectations models.   ansi-c*/


/*     //=== Better version is TSkalfilmsinputs_1stapp_tag.  Kalman filter for constant or known-time-varying DSGE models.   ansi-c*/
   typedef struct TSkalfiltv_tag
   {
/*        //General (known-time-varying) Kalman filter for DSGE models.   ansi-c*/
/*        //  It computes a sequence of one-step predictions and their covariance matrices, and the log likelihood.   ansi-c*/
/*        //  The function uses a forward recursion algorithm.  See also the Matlab function fn_kalfil_tv.m   ansi-c*/
/*        //   ansi-c*/
/*        //   State space model is defined as follows:   ansi-c*/
/*        //       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)   ansi-c*/
/*        //       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)   ansi-c*/
/*        //     where a(t), H(t), b(t), and F(t) depend on s_t that follows a Markov-chain process and are taken as given.   ansi-c*/
/*        //   ansi-c*/
/*        //   Inputs are as follows:   ansi-c*/
/*        //      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].   ansi-c*/
/*        //        a is an n_y-by-T matrix of time-varying input vectors in the measurement equation.   ansi-c*/
/*        //        H is an n_y-by-n_z-by-T 3-D of time-varying matrices in the measurement equation.   ansi-c*/
/*        //        R is an n_y-by-n_y-by-T 3-D of time-varying covariance matrices for the error in the measurement equation.   ansi-c*/
/*        //        G is an n_z-by-n_y-by-T 3-D of time-varying E(eta_t * eps_t').   ansi-c*/
/*        //        ------   ansi-c*/
/*        //        b is an n_z-by-T matrix of time-varying input vectors in the state equation with b(:,1) as an initial condition.   ansi-c*/
/*        //        F is an n_z-by-n_z-by-T 3-D of time-varying transition matrices in the state equation with F(:,:,1) as an initial condition.   ansi-c*/
/*        //        V is an n_z-by-n_z-by-T 3-D of time-varying covariance matrices for the error in the state equation with V(:,:,1) as an initial condition.   ansi-c*/
/*        //        ------   ansi-c*/
/*        //        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;   ansi-c*/
/*        //                 0: using the unconditional mean for any given regime at time 0.   ansi-c*/
/*        //        z0 is an n_z-by-1 vector of initial condition when indxIni=1. (Not used if indxIni=0.)   ansi-c*/
/*        //        P0 is an n_z-by-n_z matrix of initial condition when indxIni=1. (Not used if indxIni=0.)   ansi-c*/
/*        //   ansi-c*/
/*        //   Outputs are as follows:   ansi-c*/
/*        //      loglh is a value of the log likelihood function of the state-space model   ansi-c*/
/*        //                                under the assumption that errors are multivariate Gaussian.   ansi-c*/
/*        //      zt_tm1 is an n_z-by-T matrices of one-step predicted state vectors with z0_0m1 as a initial condition   ansi-c*/
/*        //                         and with z_{t+1|t} as the last element.  Thus, we can use it as a base-1 vector.   ansi-c*/
/*        //      Pt_tm1 is an n_z-by-n_z-by-T 3-D of covariance matrices of zt_tm1 with P0_0m1 as a initial condition   ansi-c*/
/*        //                         and with P_{t+1|t} as the last element.  Thus, we can use it as a base-1 cell.   ansi-c*/
/*        //   ansi-c*/
/*        //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:   ansi-c*/
/*        //             z0_0m1 = (I-F(:,:,1))\b(:,1)   ansi-c*/
/*        //        vec(P0_0m1) = (I-kron(F(:,:,1),F(:,:,1)))\vec(V(:,:,1))   ansi-c*/
/*        //   Note that all eigenvalues of the matrix F(:,:,1) are inside the unit circle when the state-space model is bounded (stationary).   ansi-c*/
/*        //   ansi-c*/
/*        //   March 2007, written by Tao Zha   ansi-c*/
/*        //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.   ansi-c*/

/*        //=== Input arguments.   ansi-c*/
      int ny;   /*  number of observables.   ansi-c*/
      int nz;   /*  number of state variables.   ansi-c*/
      int T;    /*  sample size.   ansi-c*/
      int indxIni;    /*  1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;   ansi-c*/
/*                       //0: using the unconditional mean for any given regime at time 0. (Default value)   ansi-c*/
      TSdmatrix *yt_dm;    /*  ny-by-T.   ansi-c*/
      TSdmatrix *at_dm;    /*  ny-by-T.   ansi-c*/
      TSdcell *Ht_dc;    /*  ny-by-nz-by-T.   ansi-c*/
      TSdcell *Rt_dc;    /*  ny-by-ny-by-T.  Covariance matrix for the measurement equation.   ansi-c*/
      TSdcell *Gt_dc;    /*  nz-by-ny-by-T.  Cross-covariance.   ansi-c*/
/*        //   ansi-c*/
      TSdmatrix *bt_dm;    /*  nz-by-T.   ansi-c*/
      TSdcell *Ft_dc;    /*  nz-by-nz-by-T.   ansi-c*/
      TSdcell *Vt_dc;    /*  nz-by-nz-by-T.  Covariance matrix for the state equation.   ansi-c*/
/*        //   ansi-c*/
      TSdvector *z0_dv;   /*  nz-by-1;   ansi-c*/
      TSdmatrix *P0_dm;    /*  nz-by-nz.   ansi-c*/

/*        //=== Output arguments.   ansi-c*/
      double loglh;   /*  log likelihood.   ansi-c*/
      TSdmatrix *zt_tm1_dm;   /*  nz-by-T.   ansi-c*/
      TSdcell *Pt_tm1_dc;    /*  nz-by-nz-T.   ansi-c*/
   } TSkalfiltv;



/*     //=== Inputs for filter for Markov-switching DSGE models at any time t.   ansi-c*/
   typedef struct TSkalfilmsinputs_1stapp_tag
   {
/*        //Inputs Markov-switching Kalman filter for DSGE models (conditional on all the regimes at time t).   ansi-c*/
/*        //  It computes a sequence of one-step predictions and their covariance matrices, and the log likelihood.   ansi-c*/
/*        //  The function uses a forward recursion algorithm.  See also the Matlab function fn_kalfil_tv.m   ansi-c*/
/*        //   ansi-c*/
/*        //   State space model is defined as follows:   ansi-c*/
/*        //       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)   ansi-c*/
/*        //       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)   ansi-c*/
/*        //     where a(t), H(t), b(t), and F(t) depend on the grand regime s_t that follows a Markov-chain process   ansi-c*/
/*        //                                                                          and is taken as given.   ansi-c*/
/*        //   ansi-c*/
/*        //   Inputs at time t are as follows where nst is number of grand regimes (including lagged regime   ansi-c*/
/*        //                                           and coefficients and shock variances):   ansi-c*/
/*        //      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].   ansi-c*/
/*        //        a is an n_y-by-nst matrix of Markov-switching input vectors in the measurement equation.   ansi-c*/
/*        //        H is an n_y-by-n_z-by-nst 3-D of Markov-switching matrices in the measurement equation.   ansi-c*/
/*        //        R is an n_y-by-n_y-by-nst 3-D of Markov-switching covariance matrices for the error in the measurement equation.   ansi-c*/
/*        //        G is an n_z-by-n_y-by-nst 3-D of Markov-switching E(eta_t * eps_t').   ansi-c*/
/*        //        ------   ansi-c*/
/*        //        b is an n_z-by-nst matrix of Markov-switching input vectors in the state equation with b(:,st) as an initial condition.   ansi-c*/
/*        //                (alternatively, with the ergodic weighted b(:,st) as an initial condition).   ansi-c*/
/*        //        F is an n_z-by-n_z-by-nst 3-D of Markov-switching transition matrices in the state equation with F(:,:,st)   ansi-c*/
/*        //                as an initial condition (alternatively, with the ergodic weighted F(:,:,st) as an initial condition).   ansi-c*/
/*        //        V is an n_z-by-n_z-by-nRv 3-D of Markov-switching covariance matrices for the error in the state equation   ansi-c*/
/*        //                with V(:,:,st) as an initial condition (alternatively, with the ergodic weighted V(:,:,st) as an initial condition).   ansi-c*/
/*        //        ------   ansi-c*/
/*        //        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;   ansi-c*/
/*        //                 0: using the unconditional mean for any given regime at time 0.   ansi-c*/
/*        //        z0 is an n_z-by-nst matrix of initial condition (Not used if indxIni=0).   ansi-c*/
/*        //        P0 is an n_z-by-n_z-by-nst 3-D of initial condition (Not used if indxIni=0).   ansi-c*/
/*        //   ansi-c*/
/*        //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:   ansi-c*/
/*        //             z0_0m1 = (I-F(:,:,st))\b(:,st)   ansi-c*/
/*        //        vec(P0_0m1) = (I-kron(F(:,:,st),F(:,:,st)))\vec(V(:,:,st))   ansi-c*/
/*        //   Note that all eigenvalues of the matrix F(:,:,st) are inside the unit circle when the state-space model is bounded (stationary).   ansi-c*/
/*        //   ansi-c*/
/*        //   November 2007, written by Tao Zha.  Revised, April 2008.   ansi-c*/
/*        //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.   ansi-c*/

/*        //=== Input arguments.   ansi-c*/
      int ny;   /*  number of observables.   ansi-c*/
      int nz;   /*  number of state variables.   ansi-c*/
      int nst;  /*  number of grand composite regimes (current and past regimes, coefficient and volatility regimes).   ansi-c*/
      int T;    /*  sample size.   ansi-c*/
      int indxIni;    /*  1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0,   ansi-c*/
/*                       //0: using the unconditional momnets for any given regime at time 0 (default when indxDiffuse = 0).   ansi-c*/
      int indxDiffuse;   /*  1: using the diffuse condition for z_{1|0} and P_{1|0} (default option), according to Koopman and Durbin, "Filtering and Smoothing of State Vector for Diffuse State-Space Models," J. of Time Series Analysis, Vol 24(1), pp.85-99.   ansi-c*/
/*                          //0: using the unconditional moments.   ansi-c*/
      double DiffuseScale;  /*  A large (infinity) number when indxDiffuse = 1.   ansi-c*/
      int ztm1_track;  /*  t-1 = -1:      no initial conditions z_{1|0} and P_{1|0} has been computed yet, but will be using InitializeKalman_z10_P10(),   ansi-c*/
/*                        //t-1 >= 0:T-1:  z_{t|t-1} and P_{t|t-1} are updated up to t-1.   ansi-c*/
      int dtm1_track;  /*  t-1 = -1:      no etdata_dc->C[0] or Dtdata_d4->F[0] has been computed yet.   ansi-c*/
/*                        //t-1 >= 0:T-1:  etdata_dc->C[t-1] and Dtdata_d4->F[t-1] are updated up to t-1.   ansi-c*/

      TSdmatrix *yt_dm;    /*  ny-by-T.   ansi-c*/
      TSdmatrix *at_dm;    /*  ny-by-nst.   ansi-c*/
      TSdcell *Ht_dc;      /*  ny-by-nz-by-nst.   ansi-c*/
      TSdcell *Rt_dc;      /*  ny-by-ny-by-nst.  Covariance matrix for the measurement equation.   ansi-c*/
      TSdcell *Gt_dc;      /*  nz-by-ny-by-nst.  Cross-covariance.   ansi-c*/
/*        //   ansi-c*/
      TSdmatrix *bt_dm;    /*  nz-by-nst.   ansi-c*/
      TSdcell *Ft_dc;      /*  nz-by-nz-by-nst.   ansi-c*/
      TSdcell *Vt_dc;      /*  nz-by-nz-by-nst.  Covariance matrix for the state equation.   ansi-c*/
/*        //   ansi-c*/
      TSdmatrix *z0_0_dm;  /*  nz-by-nst. z_{0|0}.   ansi-c*/
      TSdmatrix *z0_dm;   /*  nz-by-nst. z_{1|0}.   ansi-c*/
      TSdcell *P0_dc;     /*  nz-by-nz-by-nst. P_{1|0}   ansi-c*/


/*        //=== Output arguments only used for 1st order approximation to zt and Pt depending on infinite past regimes.   ansi-c*/
      TSdcell *zt_tm1_dc;    /*  nz-by-nst-by-(T+1), where z_{1|0} is an initial condition (1st element with t-1=0 or t=1 for base-1) and   ansi-c*/
/*                              //  the terminal condition z_{T+1|T} using Updatekalfilms_1stapp(T, ...) is not computed   ansi-c*/
/*                              //  when the likelihood logTimetCondLH_kalfilms_1stapp() is computed.  Thus, z_{T+1|T}   ansi-c*/
/*                              //  has not legal value computed in most applications unless in out-of-sample forecasting problems.   ansi-c*/
      TSdfourth *Pt_tm1_d4;  /*  nz-by-nz-by-nst-by-(T+1), where P_{1|0} is an initial condition (1st element with t-1=0) and   ansi-c*/
/*                              //  the terminal condition P_{T+1|T} using Updatekalfilms_1stapp(T, ...) is not computed   ansi-c*/
/*                              //  when the likelihood logTimetCondLH_kalfilms_1stapp() is computed.  Thus, P_{T+1|T}   ansi-c*/
/*                              //  has not legal value computed in most applications unless in out-of-sample forecasting problems.   ansi-c*/
/*        //+ Will be save for updating likelihood and Kalman filter Updatekalfilms_1stapp(), so save time to recompute these objects again.   ansi-c*/
      TSdfourth *PHtran_tdata_d4;   /*  nz-by-ny-by-nst-T, P_{t|t-1}*H_t'.  Saved only for updating Kalman filter Updatekalfilms_1stapp().   ansi-c*/
      TSdcell *etdata_dc;  /*  ny-by-nst-by-T (with base-0 T), forecast errors e_t in the likelihood.   ansi-c*/
      TSdcell *yt_tm1_dc;  /*  ny-by-nst-by-T, one-step forecast y_{t|t-1} for t=0 to T-1 (base-0). Used to back out structural shocks.   ansi-c*/
      TSdfourth *Dtdata_d4;  /*  ny-by-ny-nst-by-T, forecast covariance D_t in the likelihood.  Saved for updating Kalman filter Updatekalfilms_1stapp().   ansi-c*/
   } TSkalfilmsinputs_1stapp;


/*     //=== OLD Code: Inputs for filter for Markov-switching DSGE models at any time t.   ansi-c*/
   typedef struct TSkalfilmsinputs_tag
   {
/*        //Inputs Markov-switching Kalman filter for DSGE models (conditional on all the regimes at time t).   ansi-c*/
/*        //  It computes a sequence of one-step predictions and their covariance matrices, and the log likelihood.   ansi-c*/
/*        //  The function uses a forward recursion algorithm.  See also the Matlab function fn_kalfil_tv.m   ansi-c*/
/*        //   ansi-c*/
/*        //   State space model is defined as follows:   ansi-c*/
/*        //       y(t) = a(t) + H(t)*z(t) + eps(t)     (observation or measurement equation)   ansi-c*/
/*        //       z(t) = b(t) + F(t)*z(t) + eta(t)     (state or transition equation)   ansi-c*/
/*        //     where a(t), H(t), b(t), and F(t) depend on s_t that follows a Markov-chain process and are taken as given.   ansi-c*/
/*        //   ansi-c*/
/*        //   Inputs at time t are as follows where nRc is number of regimes for coefficients   ansi-c*/
/*        //                                         nRv is number of regimes for volatility (shock variances):   ansi-c*/
/*        //      Y_T is a n_y-by-T matrix containing data [y(1), ... , y(T)].   ansi-c*/
/*        //        a is an n_y-by-nRc matrix of Markov-switching input vectors in the measurement equation.   ansi-c*/
/*        //        H is an n_y-by-n_z-by-nRc 3-D of Markov-switching matrices in the measurement equation.   ansi-c*/
/*        //        R is an n_y-by-n_y-by-nRv 3-D of Markov-switching covariance matrices for the error in the measurement equation.   ansi-c*/
/*        //        G is an n_z-by-n_y-by-nRv 3-D of Markov-switching E(eta_t * eps_t').   ansi-c*/
/*        //        ------   ansi-c*/
/*        //        b is an n_z-by-nRc matrix of Markov-switching input vectors in the state equation with b(:,1) as an initial condition.   ansi-c*/
/*        //        F is an n_z-by-n_z-by-nRc 3-D of Markov-switching transition matrices in the state equation with F(:,:,1) as an initial condition.   ansi-c*/
/*        //        V is an n_z-by-n_z-by-nRv 3-D of Markov-switching covariance matrices for the error in the state equation with V(:,:,1) as an initial condition.   ansi-c*/
/*        //        ------   ansi-c*/
/*        //        indxIndRegimes:  1: coefficient regime and volatility regime are independent; 0: these two regimes are synchronized completely.   ansi-c*/
/*        //        indxIni: 1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;   ansi-c*/
/*        //                 0: using the unconditional mean for any given regime at time 0.   ansi-c*/
/*        //        z0 is an n_z-by-nRc*nRv matrix of initial condition when indxIni=1 and indxIndRegimes=1. (Not used if indxIni=0.)   ansi-c*/
/*        //        z0 is an n_z-by-nRv matrix of initial condition when indxIni=1 and indxIndRegimes=0. (Not used if indxIni=0.)   ansi-c*/
/*        //        P0 is an n_z-by-n_z-by-nRc*nRv 3-D of initial condition when indxIni=1 and indxIndRegimes=1. (Not used if indxIni=0.)   ansi-c*/
/*        //        P0 is an n_z-by-n_z-by-nRv 3-D of initial condition when indxIni=1 and indxIndRegimes=0. (Not used if indxIni=0.)   ansi-c*/
/*        //   ansi-c*/
/*        //   The initial state vector and its covariance matrix are computed under the bounded (stationary) condition:   ansi-c*/
/*        //             z0_0m1 = (I-F(:,:,1))\b(:,1)   ansi-c*/
/*        //        vec(P0_0m1) = (I-kron(F(:,:,1),F(:,:,1)))\vec(V(:,:,1))   ansi-c*/
/*        //   Note that all eigenvalues of the matrix F(:,:,1) are inside the unit circle when the state-space model is bounded (stationary).   ansi-c*/
/*        //   ansi-c*/
/*        //   November 2007, written by Tao Zha   ansi-c*/
/*        //   See Hamilton's book ([13.2.13] -- [13.2.22]), Harvey (pp.100-106), and LiuWZ Model I NOTES pp.001-003.   ansi-c*/

/*        //=== Input arguments.   ansi-c*/
      int ny;   /*  number of observables.   ansi-c*/
      int nz;   /*  number of state variables.   ansi-c*/
      int nRc;  /*  number of composite regimes (current and past regimes) for coefficients.   ansi-c*/
      int nRstc;   /*  number of coefficient regimes at time t.   ansi-c*/
      int nRv;  /*  number of regimes for volatility (shock variances).   ansi-c*/
      int indxIndRegimes;  /*  1: coefficient regime and volatility regime are independent; 0: these two regimes are synchronized completely.   ansi-c*/
      int T;    /*  sample size.   ansi-c*/
      int indxIni;    /*  1: using the initial condition with zt_tm1(:,1)=z0 and Pt_tm1(:,:,1)=P0;   ansi-c*/
/*                       //0: using the unconditional mean for any given regime at time 0. (Default value)   ansi-c*/
      TSdmatrix *yt_dm;    /*  ny-by-T.   ansi-c*/
      TSdmatrix *at_dm;    /*  ny-by-nRc.   ansi-c*/
      TSdcell *Ht_dc;      /*  ny-by-nz-by-nRc.   ansi-c*/
      TSdcell *Rt_dc;      /*  ny-by-ny-by-nRv.  Covariance matrix for the measurement equation.   ansi-c*/
      TSdcell *Gt_dc;      /*  nz-by-ny-by-nRv.  Cross-covariance.   ansi-c*/
/*        //   ansi-c*/
      TSdmatrix *bt_dm;    /*  nz-by-nRc.   ansi-c*/
      TSdcell *Ft_dc;      /*  nz-by-nz-by-nRc.   ansi-c*/
      TSdcell *Vt_dc;      /*  nz-by-nz-by-nRv.  Covariance matrix for the state equation.   ansi-c*/
/*        //   ansi-c*/
      TSdmatrix *z0_dm;   /*  nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nRv if indxIndRegimes == 0.   ansi-c*/
      TSdcell *P0_dc;     /*  nz-by-nz-by-nRc*nRv if indxIndRegimes == 1 or nz-by-nz-by-nRv if indxIndRegimes == 0.   ansi-c*/


/*        //=== Output arguments only used for 1st order approximation to zt and Pt depending on infinite past regimes.   ansi-c*/
      TSdcell *zt_tm1_dc;       /*  nz-by-nRc*nRv-by-T if indxIndRegimes==1, nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.   ansi-c*/
      TSdfourth *Pt_tm1_d4;     /*  nz-by-nz-by-nRc*nRv-T if indxIndRegimes==1, nz-by-nz-by-nRv-by-T if indxIndRegimes==0 where nRc=nRv.   ansi-c*/
   } TSkalfilmsinputs;




/*     //--- Functions for univariate random walk kalman filter.   ansi-c*/
   TSkalcvfurw *CreateTSkalcvfurw(TFlearninguni *func, int T, int k, int tv);   /*  , int storeZ, int storeV);   ansi-c*/
   TSkalcvfurw *DestroyTSkalcvfurw(TSkalcvfurw *kalcvfurw_ps);
   void kalcvf_urw(TSkalcvfurw *kalcvfurw_ps, void *dummy_ps);

/*     //--- New Code: Functions for Markov-switching Kalman filter.   ansi-c*/
   struct TSkalfilmsinputs_1stapp_tag *CreateTSkalfilmsinputs_1stapp(int ny, int nz, int nst, int T);
   struct TSkalfilmsinputs_1stapp_tag *DestroyTSkalfilmsinputs_1stapp(struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps);
   int InitializeKalman_z10_P10(struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, TSdmatrix *z10_dm, TSdcell *P10_dc);
   double logTimetCondLH_kalfilms_1stapp(int st, int inpt, struct TSkalfilmsinputs_1stapp_tag *kalfilmsinputs_1stapp_ps, struct TStateModel_tag *smodel_ps);



/*     //--- OLD Code: Functions for general constant Kalman filter.   ansi-c*/
   struct TSkalfiltv_tag *CreateTSkalfiltv(int ny, int nz, int T);
   struct TSkalfiltv_tag *DestroyTSkalfiltv(struct TSkalfiltv_tag *kalfiltv_ps);
/*     //Used to test tz_logTimetCondLH_kalfiltv(). (Done April 08).  double tz_kalfiltv(struct TSkalfiltv_tag *kalfiltv_ps);   ansi-c*/
   double tz_logTimetCondLH_kalfiltv(int st, int inpt, struct TSkalfiltv_tag *kalfiltv_ps);

/*     //--- OLD Code: Functions for Markov-switching Kalman filter.   ansi-c*/
   struct TSkalfilmsinputs_tag *CreateTSkalfilmsinputs(int ny, int nz, int nRc, int nRstc, int nRv, int indxIndRegimes, int T);
   struct TSkalfilmsinputs_tag *DestroyTSkalfilmsinputs(struct TSkalfilmsinputs_tag *kalfilmsinputs_ps);
   double tz_logTimetCondLH_kalfilms_1st_approx(int st, int inpt, struct TSkalfilmsinputs_tag *kalfilmsinputs_ps, struct TStateModel_tag *smodel_ps);
/*              //IMPORTANT NOTE: in the Markov-switching input file datainp_markov*.prn, it MUST be that   ansi-c*/
/*              //                                the coefficient regime is the 1st state variable, and   ansi-c*/
/*              //                                the volatility regime is the 2nd state variable.   ansi-c*/
#endif

