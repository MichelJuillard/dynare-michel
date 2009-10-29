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
#include "Estimation.h"
#include "k_ord_dynare.h"
#include "kalman.h"
#include "math.h"
#include "disclyap_fast.h"

#include <cstring>

#include <cctype>


class DsgeLikelihood 
  {
  double likelihood; // sum of vector of KF step log likelihoods
  vector<double>* vll; // vector of KF step log likelihoods
                       /*******
                       MexStruct &options_;
                       MexStruct &M_;
                       MexStruct &bayestOptions_;
                       MexStruct &dr_;
                       MexStruct &oo_;
  ***************/
  Vector& a_init;//initial value of the state, usually set to 0.
  GeneralMatrix& Q;// Kalman Matrices
  GeneralMatrix& R;
  GeneralMatrix& T;
  GeneralMatrix& Z;
  GeneralMatrix& Pstar;
  GeneralMatrix& Pinf;
  GeneralMatrix& H;
  const GeneralMatrix& data;
  GeneralMatrix& Y;
  //GeneralMatrix& currentDataSubSample;
  //int periodStart, periodEnd;  // start and end of current sub sample
  const int numPeriods;//=1; number of structural change periods
  const int numVarobs; // number of observed variables in the observation vector at time t.
  const int numTimeObs; // number of obsevations (vectors) in the time series
  const int order;
  const int endo_nbr;
  const int exo_nbr;
  const int nstatic;
  const int npred;
  const int nfwrd;
  char *fName;
  int presampleStart;
  int kalman_algo; // type of  kalman algorithm: multi- or uni-variate
  int mode_compute;
  int info;
  double penalty;
  int cost_flag;
  const int number_of_observations;
  const bool no_more_missing_observations;
  const vector<int>& order_var;
  const vector<int>& mfys;
  const vector<int>& mf; // positions of observed variables in restricted state vector for likelihood computation.
  Vector& xparam1;  // all estimated parameters incl sderrs
  const int num_dp; // number of deep parameters 
  Vector& deepParams; // estimated deep parameters subset of xparam1 only
  const Vector& param_ub;  // upper and lower bounds
  const Vector& param_lb;
  const vector<int>&pshape;
  const Vector& p6;
  const Vector& p7;
  const Vector& p3;
  const Vector& p4;
  Vector& SteadyState;
  Vector& constant;
  GeneralParams& dynareParams;
  //GeneralParams& parameterDescription;
  GeneralParams& dr;
  
  GeneralMatrix& kstate;
  GeneralMatrix& ghx;
  GeneralMatrix& ghu;

  //DynamicModelDLL* dynamicDLLp;
  KordpDynare* model;// to be initialised by high level calling function
  Approximation* approx;
  //friend class BasicKalmanTask;
  //BasicKalmanTask bkt;
  //friend class KalmanUniTask;
  //KalmanUniTask ukt;// univariate
  // member functions
  MexStructParam& SetDRModel(MexStructParam&params);
  void disclyap_fast(const GeneralMatrix &G, const GeneralMatrix & V, GeneralMatrix &X, double tol = 1e-16, int flag_ch=0);
  GeneralMatrix& SolveDRkOrderPert();//calls k-order pert or whatever;
  int dynareResolveDR(vector<int>&iv,vector<int>&ic,GeneralMatrix& aux);  // returns int info, ys, and TT and RR Decision Rule
  int SolveDRModel(const int endo_nbr, const int exo_nbr, const int nstatic, const int npred, int nfwrd);//int dr1();  // returns int info and updated dr 
  int updateQHparams();// updates Q and H matrices and deep parameters
  int InitiateKalmanMatrices();
  void DataPreparation(MexStructParam&params, const GeneralMatrix &data);
  double KalmanFilter(double riccatiTol,bool uni);// calls Kalman
  
  public:
    DsgeLikelihood( Vector& inA_init, GeneralMatrix& inQ,  GeneralMatrix& R,
      GeneralMatrix& inT,  GeneralMatrix& inZ,  GeneralMatrix& inPstar,  GeneralMatrix& inPinf,
      GeneralMatrix& inH,  const GeneralMatrix&inData,  GeneralMatrix&inY,  
      const int INnumPeriods, //  const int INnumVarobs, //  const int INnumTimeObs,
      const int INorder, const int INendo_nbr, const int INexo_nbr, const int INnstatic,  
      const int INnpred, const int INnfwrd,   const int INnum_of_observations, const bool INno_more_missing_observations,
      const vector<int>& INorder_var, const vector<int>& INmfys, const vector<int>& INmf,
      Vector& INxparam1, const int INnum_dp, Vector& INdeepParams,
      const Vector& INub, const Vector& INlb, const vector<int>&INpshape,
      const Vector&INp6, const Vector&INp7, const Vector&INp3, const Vector&INp4,
      Vector& INSteadyState,   Vector& INconstant,  GeneralParams& INdynareParams,
      //GeneralParams& parameterDescription, 
      GeneralParams& INdr, GeneralMatrix& INkstate, GeneralMatrix& INghx,  GeneralMatrix& INghu 
      ,char *dfExt); //, KordpDynare& inModel, Approximation& INapprox );
    DsgeLikelihood( const Vector&params,const GeneralMatrix&data, const vector<int>& data_index, const int gend,
      const int number_of_observations, const bool no_more_missing_observations);//, KordpDynare& model ); // constructor, and
    DsgeLikelihood( GeneralParams& options_,GeneralParams& M_,GeneralParams& bayestopt_, GeneralMatrix& inData,
      KordpDynare& model); // constructor
      ~DsgeLikelihood();// destructor
    double CalcLikelihood(Vector& xparams);// runs all routines needed to calculate likelihood
    double getLik(){return likelihood;}
    int getInfo(){return info;}
    vector<double>& getLikVector() {return *vll;}  // vector of log likelihoods for each Kalman step
    //GeneralMatrix&lyapunov_symm(const GeneralMatrix &G, const GeneralMatrix & V);
  };
  
  
