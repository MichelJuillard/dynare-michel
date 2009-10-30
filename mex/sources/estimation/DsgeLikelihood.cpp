/*
* Copyright (C) 2008-2009 Dynare Team
*
* This file is part of Dynare.
*
* Dynare is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Dynare is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "mexutils.h"
#include "DsgeLikelihood.h"
#include "utils.h"
#include "disclyap_fast.h"
#include "ts_exception.h"
#include "dynamic_dll.h"

DsgeLikelihood::DsgeLikelihood( Vector& inA_init, GeneralMatrix& inQ,  GeneralMatrix& inR,
                               GeneralMatrix& inT,  GeneralMatrix& inZ,  GeneralMatrix& inPstar,  GeneralMatrix& inPinf,
                               GeneralMatrix& inH,  const GeneralMatrix&inData,  GeneralMatrix&inY,  
                               const int INnumPeriods, //  const int INnumVarobs, //  const int INnumTimeObs,
                               const int INorder, const int INendo_nbr, const int INexo_nbr, const int INnstatic,  
                               const int INnpred, const int INnfwrd,   const int INnum_of_observations, const bool INno_more_missing_observations,
                               const vector<int>& INorder_var, const vector<int>& INmfys, const vector<int>& INmf1,
                               Vector& INxparam1, const int INnum_dp, Vector& INdeepParams,
                               const Vector& INub, const Vector& INlb, const vector<int>&INpshape,
                               const Vector&INp6, const Vector&INp7, const Vector&INp3, const Vector&INp4,
                               Vector& INSteadyState,   Vector& INconstant,  GeneralParams& INdynareParams, //GeneralParams& parameterDescription, 
                               GeneralParams& INdr, GeneralMatrix& INkstate, GeneralMatrix& INghx,  GeneralMatrix& INghu
                               ,const int jcols, const char *dfExt)//, KordpDynare& kOrdModel, Approximation& INapprox ) 
                               :a_init(inA_init), Q(inQ), R(inR), T(inT), Z(inZ), Pstar(inPstar), Pinf(inPinf), H(inH), data(inData), Y(inY), 
                               numPeriods(INnumPeriods),   numVarobs(inData.numRows()),  numTimeObs(inData.numCols()),
                               order(INorder),  endo_nbr(INendo_nbr ),  exo_nbr(INexo_nbr),  nstatic(INnstatic),  npred(INnpred),
                               nfwrd(INnfwrd), number_of_observations(INnum_of_observations), no_more_missing_observations(INno_more_missing_observations),
                               order_var(INorder_var),   mfys(INmfys), mf(INmf1), xparam1(INxparam1),
                               num_dp(INnum_dp), deepParams(INdeepParams), //num_dp((int)dynareParams.getDoubleField(string("np),// no of deep params
                               param_ub(INub), param_lb(INlb),pshape(INpshape), //pshape(dynareParams.getIntVectorField(string("pshape),
                               p6(INp6), p7(INp7), p3(INp3), p4(INp4), SteadyState(INSteadyState),
                               constant(INconstant), dynareParams(INdynareParams), 
                               dr(INdr), kstate(INkstate), ghx(INghx),ghu(INghu)   
                               //, model(kOrdModel), approx(INapprox )
  {
 
  /*****
  bayestOptions_("caller","bayestopt_");
  options_("caller","options_");
  M_Options_("caller","M_");
  dr_("caller","dr");
  oo_("caller","oo_");
  *********/
  // setting some frequently used common variables that do not need updating
  //std::vector<double>* vll=new std::vector<double> (nper);
  vll=new std::vector<double> (numTimeObs);// vector of likelihoods
  kalman_algo=(int)dynareParams.getDoubleField(string("kalman_algo"));
  presampleStart=1+(int)dynareParams.getDoubleField(string("presample"));
  mode_compute=(int)dr.getDoubleField(string("mode_compute"));

  // Pepare data for Constructing k-order-perturbation classes
  //const char *
  string fname=dynareParams.getStringField(string("fname"));
  fName = (char *)fname.c_str(); 
  double qz_criterium = dynareParams.getDoubleField(string("qz_criterium"));//qz_criterium = 1+1e-6;
  int nMax_lag =(int)dynareParams.getDoubleField(string("maximum_lag"));
  const int nBoth=(int)dr.getDoubleField(string("nboth"));
  const int nPred = npred-nBoth; // correct nPred for nBoth.
  //vector<int> *var_order_vp = &order_var;
  vCov =  new TwoDMatrix(Q);
  // the lag, current and lead blocks of the jacobian respectively
  llincidence = new TwoDMatrix (dynareParams.getMatrixField(string("lead_lag_incidence")));
  charArraySt * casOrdEndoNames=dynareParams.getCharArrayField(string("var_order_endo_names"));
  const char **endoNamesMX=(const char ** )casOrdEndoNames->charArrayPtr;

#ifdef DEBUG
  for (int i = 0; i < endo_nbr; i++)
    {
      mexPrintf("k_ord_perturbation: EndoNameList[%d][0]= %s.\n", i, endoNamesMX[i]);
    }
#endif

  charArraySt * casExoNames=dynareParams.getCharArrayField(string("exo_names"));
  const char **exoNamesMX=(const char ** )casExoNames->charArrayPtr;

  Vector &NNZD =dynareParams.getDoubleVectorField(string("NNZDerivatives"));

  const int nSteps = 0; // Dynare++ solving steps, for time being default to 0 = deterministic steady state
  const double sstol = 1.e-13; //NL solver tolerance from


  // Construct k-order-perturbation classes

    THREAD_GROUP::max_parallel_threads = 2; //params.num_threads;

    try
      {
        // make journal name and journal
        std::string jName(fName); //params.basename);
        jName += ".jnl";
        journal= new Journal(jName.c_str());
#ifdef DEBUG
        mexPrintf("k_order_perturbation: Calling dynamicDLL constructor.\n");
#endif
        dynamicDLLp=new DynamicModelDLL (fName, endo_nbr, jcols, nMax_lag, exo_nbr, dfExt);

        // intiate tensor library
#ifdef DEBUG
        mexPrintf("k_order_perturbation: Call tls init\n");
#endif
        tls.init(order, nstatic+2*nPred+3*nBoth+2*nfwrd+exo_nbr);

#ifdef DEBUG
        mexPrintf("estimation: Calling dynare model constructor .\n");
#endif
        // make KordpDynare object
        model=new KordpDynare (endoNamesMX,  endo_nbr, exoNamesMX,  exo_nbr, num_dp,//nPar, // paramNames,
                           &SteadyState, vCov, &deepParams /*modParams*/, nstatic, nPred, nfwrd, nBoth,
                           jcols, &NNZD, nSteps, order, *journal, *dynamicDLLp, 
                           sstol, &order_var /*var_order_vp*/, llincidence, qz_criterium);

        // construct main K-order approximation class
#ifdef DEBUG
        mexPrintf("estimation: Call Approximation constructor with qz_criterium=%f \n", qz_criterium);
#endif
        approx= new Approximation(*model, *journal,  nSteps, false, qz_criterium);
        // run stochastic steady
#ifdef DEBUG
        mexPrintf("estimation:k_order_perturbation and Approximation created.\n");
#endif
      }
    catch (const KordException &e)
      {
        printf("Caugth Kord exception: ");
        e.print();
        mexPrintf("Caugth Kord exception: %s", e.get_message());
        std::string errfile(fName); //(params.basename);
        errfile += "_error.log";
        FILE *errfd = NULL;
        if (NULL == (errfd = fopen(errfile.c_str(), "wb")))
          {
            fprintf(stderr, "Couldn't open %s for writing.\n", errfile.c_str());
            return; // e.code();
          }
        fprintf(errfd, "Caugth Kord exception: %s", e.get_message());
        fclose(errfd);
        return; // e.code();
      }
    catch (const TLException &e)
      {
        printf("Caugth TL exception: ");
        e.print();
        return; // 255;
      }
    catch (SylvException &e)
      {
        printf("Caught Sylv exception: ");
        e.printMessage();
        return; // 255;
      }
    catch (const DynareException &e)
      {
        printf("Caught KordpDynare exception: %s\n", e.message());
        mexPrintf("Caugth Dynare exception: %s", e.message());
        std::string errfile(fName); //(params.basename);
        errfile += "_error.log";
        FILE *errfd = NULL;
        if (NULL == (errfd = fopen(errfile.c_str(), "wb")))
          {
            fprintf(stderr, "Couldn't open %s for writing.\n", errfile.c_str());
            return; // e.code();
          }
        fprintf(errfd, "Caugth KordDynare  exception: %s", e.message());
        fclose(errfd);
        return; // 255;
      }
    catch (const ogu::Exception &e)
      {
        printf("Caught ogu::Exception: ");
        e.print();
        mexPrintf("Caugth general exception: %s", e.message());
        std::string errfile(fName); //(params.basename);
        errfile += "_error.log";
        FILE *errfd = NULL;
        if (NULL == (errfd = fopen(errfile.c_str(), "wb")))
          {
            fprintf(stderr, "Couldn't open %s for writing.\n", errfile.c_str());
            return; // e.code();
          }
        e.print(errfd);
        fclose(errfd);
        return; // 255;
      }  //catch
  };


DsgeLikelihood::~DsgeLikelihood()
  {
  delete approx;
  delete model;
  delete dynamicDLLp;
  delete journal;
  delete llincidence;
  delete vCov;
  delete vll;
  delete &H;
  delete &Q;
  delete &SteadyState;
  delete &kstate;
  delete &param_ub;
  delete &param_lb;
  delete &pshape;
  delete &p6;
  delete &p7;
  delete &p3;
  delete &p4;
  delete &xparam1;
  delete &deepParams;
  delete &ghx;
  delete &ghu;
  }

double 
DsgeLikelihood::CalcLikelihood(Vector& xparams)
// runs all routines needed to calculate likelihood
  {
  likelihood=0.0;
  xparam1=xparams;
  /*******************************
  * loop for each sub-sample period
  ********************************/
  for (int sslc=0;sslc<numPeriods;++sslc)
    {
    
    /*****************************************************************************--
    % 1. Get the structural parameters & define penalties
    ******************************************************************************-*/
    cost_flag  	= 1;
    int i;
    if (mode_compute != 1)
      {
      //    vector<int>kk(0);
      double qdelta=0;
      for(i=0;i<xparam1.length();++i)
        {
        if(xparam1[i]<param_lb[i])   //       kk.push_back[i+1];
          qdelta+=(xparam1[i]-param_lb[i])*(xparam1[i]-param_lb[i]);
        }
      if ( qdelta>0) // i.e. kk.size()>0)
        {
        //  fval = bayestopt_.penalty+sum((bayestopt_.lb(k)-xparam1(k)).^2);
        likelihood = penalty+qdelta;
        cost_flag = 0;
        info = 41;
        return likelihood;
        }
      qdelta=0;
      //    kk.clear();
      for(i=0;i<xparam1.length();++i)
        {
        if(xparam1[i]>param_ub[i])  //        kk.push_back[i+1];
          qdelta+=(xparam1[i]-param_ub[i])*(xparam1[i]-param_ub[i]);
        }
      if ( qdelta>0) // i.e. kk.size()>0)
        {
        //fval = bayestopt_.penalty+sum((xparam1(k)-bayestopt_.ub(k)).^2);
        likelihood = penalty+qdelta;
        cost_flag = 0;
        info = 42;
        return likelihood;
        }
      } // mode compute
    
    if(info=updateQHparams()) // updates Q and H matrices and deep parameters
      return likelihood;
    
      /*****************************************************************************--
      % 2. call model setup & reduction program and pre-filter data
    ******************************************************************************-*/
    GeneralMatrix& aux = dynareParams.getMatrixField(string("bayestopt_.restrict_aux"));
    vector<int>&iv= dynareParams.getIntVectorField(string("restrict_var_list"));
    vector<int>&ic= dynareParams.getIntVectorField(string("bayestopt_.restrict_columns"));
    
    /*************************************************************
    // dynare_resolve(()  // ... comes here doing:
    //	resol: 
    //    check if ys is steady state and calculate one i not
    //	  dr
    //  kalman_transition_matrix(out: A,B, in dr)
    // IN: bayestopt_.restrict_var_list, bayestopt_.restrict_columns, bayestopt_.restrict_aux, )
    ***************************************************************/
    if (info = dynareResolveDR (iv,ic, aux)) //OUT: [T,R,SteadyState], 
      return likelihood=penalty;
    
      /*****************************************************************************--
      % 2.b   pre-filter and detrend data
    ******************************************************************************-*/
    
    //if options_.noconstant
    if (dynareParams.getCharField(string("noconstant")))
      constant.zeros();
    else
      {
      //if options_.loglinear
      if (dynareParams.getCharField(string("loglinear")))
        {
        for (i =0;i<numVarobs;++i)
          constant[i] = log(SteadyState[mfys[i]]);
        }
      else
        {
        for (i =0;i<numVarobs;++i)
          constant[i] = SteadyState[mfys[i]];
        }
      }
    Vector trend_coeff(numVarobs);
    //trend = repmat(constant,1,gend);
    GeneralMatrix constMx(constant.base(),numVarobs,1);
    GeneralMatrix&trend = constMx.repmat(1,numTimeObs);
    //if bayestopt_.with_trend
    if (dynareParams.getCharField(string("with_trend")))
      {
      trend_coeff.zeros();
      //      GeneralMatrix& mt = dynareParams.getMatrixField(string("trend_coeffs"));
      //      Vector vt(mt.base, MAX(mt.numCols(), mt.numRows()));
      Vector& vtc = dynareParams.getDoubleVectorField(string("trend_coeffs"));
      for (i=0;i<vtc.length();++i)
        {
        if (vtc[i]!=0.0)
          trend_coeff[i] = vtc[i];
        }
      //trend = repmat(constant,1,gend)+trend_coeff*[1:gend];
      GeneralMatrix trend_coefMx(numVarobs, numTimeObs);
      for (i=1;i<=numTimeObs;++i)
        for (int j=0;j<numVarobs;++j)
          trend_coefMx.get(j,i)=trend_coeff[j]*i;
        
        trend.add(1,trend_coefMx);
      }
    presampleStart =(int) dynareParams.getDoubleField(string("presample"))+1;
    int no_missing_data_flag = (number_of_observations==numTimeObs*numVarobs);
    //Y =data-trend;
    Y=data;
    Y.add(-1,trend);
    
    /*****************************************************************************
    % 3. Initial condition of the Kalman filter
    *******************************************************************************/
    if( InitiateKalmanMatrices())
      return likelihood=penalty;
    
    /*****************************************************************************
    // 4. Likelihood evaluation
    // choose and run KF to get likelihood fval
    *****************************************************************************/
    likelihood+=KalmanFilter(0.000001, false);// calls Kalman
    
    
    /****************************************************************************
    // Adds prior if necessary
    ****************************************************************************/
    //likelihood-= priordens(xparam1,pshape,p6,p7,p3,p4);//fval    = (likelihood-lnprior);
    //options_.kalman_algo = kalman_algo;
    } // end sub-sample loop
  }
  
/**************************************************
* lower level, private Member functions definitions
***************************************************/


/*****************************************************************************--
% 1. Get the structural parameters & define penalties
******************************************************************************-*/
int
DsgeLikelihood::updateQHparams()// updates Q and H matrices and deep parameters
  {
  int i=0,  offset=0, nv=0, k, k1, k2, info=0;
  Q =  dynareParams.getMatrixField(string("Sigma_e"));
  H =  dynareParams.getMatrixField(string("H"));
  nv=(int)dynareParams.getDoubleField(string("nvx"));
  if(nv)
    {
    GeneralMatrix&estvx=dynareParams.getMatrixField(string("var_exo"));
    for (i=0;i<nv;++i)
      {
      k =(int)estvx.get(i,0)-1;
      Q.get(k,k) = xparam1[i]*xparam1[i];
      }
    offset = nv;
    }
  
  nv=(int)dynareParams.getDoubleField(string("nvn"));
  if(nv)
    {
    GeneralMatrix&estvn=dynareParams.getMatrixField(string("var_endo"));
    for (i=0;i<nv;++i)
      {
      k =(int)estvn.get(i,0)-1;
      H.get(k,k) = xparam1[i+offset]*xparam1[i+offset];
      }
    offset += nv;
    }
  
  //if estim_params_.ncx
  //for i=1:estim_params_.ncx
  nv=(int)dynareParams.getDoubleField(string("ncx"));
  if(nv)
    {
    GeneralMatrix&corrx=dynareParams.getMatrixField(string("corrx"));
    for (i=0;i<nv;++i)
      {
      k1 =(int)corrx.get(i,0)-1;
      k2 =(int)corrx.get(i,1)-1;
      Q.get(k1,k2) = xparam1[i+offset]*sqrt(Q.get(k1,k1)*Q.get(k2,k2));
      Q.get(k2,k1) = Q.get(k1,k2);
      }
    //   [CholQ,testQ] = chol(Q);
    int testQ=0;
    try
      {
      NormCholesky chol(Q);
      }
    catch(const TSException &e)
      {
      //      if (string(e.getMessage())==sting("The matrix is not positive definite in NormCholesky constructor"))
      if (0==strncmp(e.getMessage(),"The matrix is not positive definite in NormCholesky constructor",35))
        testQ=1;
      else
        {
        printf("Caugth unhandled TS exception with Q matrix: ");
        likelihood=penalty;
        TS_RAISE(e.getMessage());
        }
      }
    if (testQ) 
      {
      // The variance-covariance matrix of the structural innovations is not definite positive.
      // We have to compute the eigenvalues of this matrix in order to build the penalty.
      double delta=0;
      VDVFact eigQ(Q);  // get eigenvalues
      //k = find(a < 0);
      if(eigQ.hasConverged())
        {
        const Vector& evQ=eigQ.getD();
        for (i=0;i<evQ.length();++i)
          if (evQ[i]<0)
            delta-=evQ[i];
        }
      
      likelihood = penalty+delta;// +sum(-a(k));
      cost_flag = 0;
      info = 43;
      return info;
      }
    //offset = offset+estim_params_.ncx;
    offset += nv;
    }//end
  
  //if estim_params_.ncn 
  //for i=1:estim_params_.ncn
  nv=(int)dynareParams.getDoubleField(string("ncn"));
  if(nv)
    {
    GeneralMatrix&corrn=dynareParams.getMatrixField(string("corrn"));
    vector<int>&lgyidx2varobs= dynareParams.getIntVectorField(string("lgyidx2varobs"));
    for (i=0;i<nv;++i)
      {
      //      k1 = options_.lgyidx2varobs(estim_params_.corrn(i,1));
      //      k2 = options_.lgyidx2varobs(estim_params_.corrn(i,2));
      k1 = lgyidx2varobs[(int)corrn.get(i,0)-1];
      k2 = lgyidx2varobs[(int)corrn.get(i,1)-1];
      //      H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
      //      H(k2,k1) = H(k1,k2);
      H.get(k1,k2) = xparam1[i+offset]*sqrt(H.get(k1,k1)*H.get(k2,k2));
      H.get(k2,k1) = H.get(k1,k2);
      }
    
    //[CholH,testH] = chol(H);
    int testH=0;
    try
      {
      NormCholesky chol(H);
      }
    catch(const TSException &e)
      {
      //      if (string(e.getMessage)==sting("The matrix is not positive definite in NormCholesky constructor"))
      if (0==strncmp(e.getMessage(),"The matrix is not positive definite in NormCholesky constructor",35))
        testH=1;
      else
        {
        printf("Caugth unhandled TS exception with H matrix: ");
        likelihood=penalty;
        TS_RAISE((const char*)e.getMessage());
        }
      }
    if (testH)
      {
      //a = diag(eig(H));
      double delta=0;
      VDVFact eigH(H);  // get eigenvalues
      //k = find(a < 0);
      if(eigH.hasConverged())
        {
        const Vector& evH=eigH.getD();
        for (i=0;i<evH.length();++i)
          if (evH[i]<0)
            delta-=evH[i];
        }
      likelihood = penalty+delta; // +sum(-a(k));
      cost_flag = 0;
      info = 44;
      return info;
      }; //   end if
    //offset = offset+estim_params_.ncn;
    offset += nv;
    }
  
  //if estim_params_.np > 0  // i.e. num of deep parameters >0
  //   M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
  if(num_dp > 0)
    {
    if(xparam1.length()>=offset+num_dp) 
      memcpy(deepParams.base(), xparam1.base()+offset*sizeof(double),num_dp*sizeof(double));
    else
      TS_RAISE("Inssuficient length of the xparam1 parameters vector");
    }
  
    /**********
    M_.Sigma_e = Q;
    M_.H = H;
  ******************/
  return info;
  };

/*****************************************************************************--
% 3. Initial condition of the Kalman filter
******************************************************************************-*/
int
DsgeLikelihood::InitiateKalmanMatrices()
  {
  int np = T.numRows();// size(T,1);
  double lyapunov_tol=dynareParams.getDoubleField(string("lyapunov_tol"));
  int lik_init=(int)dynareParams.getDoubleField(string("lik_init"));
  //if options_.lik_init == 1     % Kalman filter
  GeneralMatrix RQRt(R,Q); // R*Q
  RQRt.multRightTrans(R); // R*Q*Rt
  GeneralMatrix Ptmp(np,np);
  //Pstar = lyapunov_symm (T,R*Q*R',options_.qz_criterium,options_.lyapunov_complex_threshold)
  disclyap_fast(T,RQRt,Ptmp, lyapunov_tol, 1); // 1 to check chol 
  Pstar=Ptmp;
  //Pinf=[]
  //Pinf  = *(new GeneralMatrix(np,np));
  Pinf.zeros();
  
  //Z= zeros(size(data,1), size(T,2))
  GeneralMatrix Ztmp=*(new GeneralMatrix(numVarobs,np));
  Ztmp.zeros();
  
  //a=zeros(size(T,1),1);
  Vector atmp(np);
  atmp.zeros();
  a_init=atmp;
  
  //if (lik_init == 2)// Old Diffuse Kalman filter
  //  Pstar = options_.Harvey_scale_factor*eye(np);
  //Pinf = [];
  //else if (lik_init == 3) // Diffuse Kalman filter
  // else ...
  
  for (int i = 0;i<numVarobs;++i)
    Ztmp.get(i,mf[i])=1.0;
  Z=Ztmp;
  delete &Ztmp;
  }  


  /*****************************************************************************
  // 4. Likelihood evaluation
  // choose and run KF to get likelihood fval
*****************************************************************************/
double 
DsgeLikelihood::KalmanFilter(double riccatiTol=0.000001,bool uni = false)
  {
  bool diffuse=false;
  
  try 
    {
    // make input matrices
    int start = presampleStart;
    
    int nper=Y.numCols();
#ifdef DEBUG		
    mexPrintf("kalman_filter: periods=%d start=%d, a.length=%d, uni=%d diffuse=%d\n", nper, start,a_init.length(), uni, diffuse);
#endif		
    
    // make storage for output
    int per;
    int d;
    // create state init
    StateInit* init = NULL;
    
    if (diffuse||uni) 
      {
      if (diffuse) 
        {
        init = new StateInit(Pstar, Pinf, a_init);
        } 
      else 
        {
        init = new StateInit(Pstar, a_init);
        }
      // fork, create objects and do filtering
      KalmanTask kt(Y, Z, H, T, R, Q, *init);
      if (uni) 
        {
        KalmanUniTask kut(kt);
        likelihood = kut.filter(per, d, (start-1), vll);
        per = per / Y.numRows();
        d = d / Y.numRows();
        } 
      else 
        {
#ifdef TIMING_LOOP
        for (int tt=0;tt<1000;++tt)
          {
#endif
          likelihood = kt.filter(per, d, (start-1), vll);
#ifdef TIMING_LOOP
          }
        mexPrintf("kalman_filter: finished 1000 loops");
#endif
        }
      }
    else // basic Kalman
      {
      init = new StateInit(Pstar, a_init);
      BasicKalmanTask bkt(Y, Z, H, T, R, Q, *init, riccatiTol);
#ifdef TIMING_LOOP
      for (int tt=0;tt<1000;++tt)
        {
#endif
        likelihood = bkt.filter( per, d, (start-1), vll);
#ifdef DEBUG		
        mexPrintf("Basickalman_filter: likelihood=%f \n", likelihood);
#endif		
#ifdef TIMING_LOOP
        }
      mexPrintf("Basickalman_filter: finished 1000 loops");
#endif
      
      }
    // destroy init
    delete init;
    } 
  catch (const TSException& e) 
    {
    mexErrMsgTxt(e.getMessage());
    } 
  catch (SylvException& e) 
    {
    char mes[300];
    e.printMessage(mes, 299);
    mexErrMsgTxt(mes);
    }
  return likelihood;
  }

