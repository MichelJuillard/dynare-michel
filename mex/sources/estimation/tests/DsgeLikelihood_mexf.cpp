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


// set_state_space and kstate preamble 
// should be performed before calling DageLikelihood, not repeatedly withing dr1.

/******************************************
* mexFunction: Matlab Inerface point and the main application driver 
* for DsgeLikelihood
*****************************************************************
% function [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data,data_index,
%                        number_of_observations,no_more_missing_observations)
% Evaluates the posterior kernel of a dsge model. 
% 
% INPUTS 
%   xparam1                        [double]   vector of model parameters.
%   gend                           [integer]  scalar specifying the number of observations.
%   data                           [double]   matrix of data
%   data_index                     [cell]     cell of column vectors
%   number_of_observations         [integer]
%   no_more_missing_observations   [integer]
%
% OUTPUTS 
%   fval        :     value of the posterior kernel at xparam1.
%   cost_flag   :     zero if the function returns a penalty, one otherwise.
%   ys          :     steady state of original endogenous variables
%   trend_coeff :
%   info        :     vector of informations about the penalty:
%                     41: one (many) parameter(s) do(es) not satisfied the lower bound
%                     42: one (many) parameter(s) do(es) not satisfied the upper bound
%   vll         :     vector of time-step log-likelihoods at xparam1.
%
*****************************************************************/
#include "DsgeLikelihood.h"

#include "mexutils.h"

extern const char *DynareParamStructsNm []={"M_", "oo_", "options_", "bayestopt_", "estim_params_", "dr"};
extern const char* mexBase[]={"base", "caller", "global"};

extern "C" {

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  if (nrhs < 6)
    mexErrMsgTxt("Must have at least 6 input parameters.\n");
  if (nlhs == 0)
    mexErrMsgTxt("Must have at least 1 output parameter.\n");

  GeneralMatrix params(mxGetPr(prhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
  if(1!= MIN(params.numCols(),params.numRows()))
    throw SYLV_MES_EXCEPTION("Vextor is 2D Matrix!");
  Vector xparam1(params.getData());//(params.base(), MAX(params.numCols(),params.numRows()));
  const int nper = (const int)mxGetScalar(prhs[1]); //gend
  GeneralMatrix data(mxGetPr(prhs[2]), mxGetM(prhs[2]), mxGetN(prhs[2]));
  const int num_of_observations = (const int)mxGetScalar(prhs[4]); 
  const bool no_more_missing_observations= (const bool)mxGetScalar(prhs[5]);



    const char *dfExt = NULL; //Dyanamic file extension, e.g.".dll" or .mexw32;
    if (prhs[5] != NULL)
      {
        const mxArray *mexExt = prhs[6];
        dfExt = mxArrayToString(mexExt);
      }

#ifdef DEBUG
    mexPrintf("estimation: mexExt=%s.\n", dfExt);
#endif
/***********
***************/
  int numPeriods=1;

  //MexStruct dynareParams();
  MexStruct& dynareParams=*(new MexStruct(numParStructs));
#ifdef DEBUG
        mexPrintf("getting dr\n");
#endif
  MexStructParam& dr=dynareParams.getMexStructField("dr");
  vector<int>&mfys=dynareParams.getIntVectorField("mfys");
  vector<int>&mf=dynareParams.getIntVectorField("mf1");
#ifdef DEBUG
        mexPrintf("getting SS\n");
#endif
  Vector& SteadyState=dr.getDoubleVectorField(string("ys"));
#ifdef DEBUG
  int gg;
  for ( gg=0;gg<SteadyState.length();++gg)
        mexPrintf("SteadyState %d = %f\n", gg, SteadyState[gg]);
#endif

  int numVarobs=data.numRows();

  Vector constant(numVarobs);//=*(new Vector(numVarobs));//   = *(new Vector(nobs));
  GeneralMatrix&kstate = dr.getMatrixField(string("kstate"));
  vector<int>&order_var = dr.getIntVectorField(string("order_var"));
#ifdef DEBUG
  for ( gg=0;gg<order_var.size();++gg)
        mexPrintf("order_var %d = %d\n", gg, order_var[gg]);
#endif
  int order=(int)dynareParams.getDoubleField(string("order"));
  int endo_nbr = (int)dynareParams.getDoubleField(string("endo_nbr"));
  int exo_nbr = (int)dynareParams.getDoubleField(string("exo_nbr"));
  int nstatic = (int)dr.getDoubleField(string("nstatic"));
  int npred = (int)dr.getDoubleField(string("npred"));
  int nfwrd = (int)dr.getDoubleField(string("nfwrd"));
  Vector& ub=dynareParams.getDoubleVectorField(string("ub"));
  Vector& lb=dynareParams.getDoubleVectorField(string("lb"));
#ifdef DEBUG
  for ( gg=0;gg<lb.length();++gg)
        mexPrintf("lb %d = %f\n", gg, lb[gg]);
#endif
  int num_dp=(int)dynareParams.getDoubleField(string("np"));// no of deep params
  Vector& deepParams=*(new Vector(num_dp));
  vector<int>&pshape=dynareParams.getIntVectorField(string("pshape"));
  Vector& p6= dynareParams.getDoubleVectorField(string("p6"));
  Vector& p7= dynareParams.getDoubleVectorField(string("p7"));
  Vector& p3= dynareParams.getDoubleVectorField(string("p3"));
  Vector& p4= dynareParams.getDoubleVectorField(string("p4"));

  //const int jcols = nExog+nEndo+nsPred+nsForw; // Num of Jacobian columns
  int nsPred=(int)dr.getDoubleField(string("nspred"));
  int nsForw=(int)dr.getDoubleField(string("nsfwrd"));
  const int jcols = exo_nbr+endo_nbr+nsPred+nsForw;
#ifdef DEBUG
        mexPrintf("jcols = %d, exo_nbr=%d\n", jcols, exo_nbr);
#endif

  GeneralMatrix& aux = dynareParams.getMatrixField(string("restrict_aux"));
  vector<int>&iv= dynareParams.getIntVectorField(string("restrict_var_list"));
  vector<int>&ic= dynareParams.getIntVectorField(string("restrict_columns"));

  int nr=iv.size()+aux.numRows(); // Size of T matrix
  Vector& a_init=*(new Vector(nr));
  a_init.zeros();

  GeneralMatrix& Q = dynareParams.getMatrixField(string("Sigma_e"));
#ifdef DEBUG
  Q.print();
#endif
  GeneralMatrix& Hrtmp = dynareParams.getMatrixField(string("H"));
  GeneralMatrix * Hp;
  if (Hrtmp.numCols()==0 || Hrtmp.numRows()==0)
    {
    delete &Hrtmp;
    Hp = new GeneralMatrix(numVarobs,numVarobs);
    Hp->zeros();
#ifdef DEBUG
    mexPrintf("finished local initialising of H \n");
#endif
    }
  else
    Hp=&Hrtmp;

  GeneralMatrix& H=*Hp;
  GeneralMatrix Y(data.numRows(),data.numCols());
  GeneralMatrix T(nr,nr);
  GeneralMatrix Z(numVarobs,nr);
  Z.zeros();
  for (int i = 0;i<numVarobs;++i)
    Z.get(i,mf[i]-1)=1;
  GeneralMatrix Pstar(nr,nr);
  GeneralMatrix R(nr,exo_nbr);
  GeneralMatrix ghx(endo_nbr,nsPred);
  GeneralMatrix ghu(endo_nbr,exo_nbr);

//Pinf=[]
  GeneralMatrix Pinf (nr,nr);
  Pinf.zeros();
double loglikelihood;
  try
    {

#ifdef DEBUG
        mexPrintf("Try construction of DsgeLikelihood\n");
#endif
    DsgeLikelihood dl( a_init, Q, R,T, Z, Pstar, Pinf, H,data,Y,  
          numPeriods, //  const int INnumVarobs, //  const int INnumTimeObs,
          order, endo_nbr, exo_nbr, nstatic, npred,nfwrd, num_of_observations, 
          no_more_missing_observations, order_var, mfys, mf, xparam1,
          num_dp, deepParams, ub, lb, pshape, p6, p7, p3, p4, SteadyState, constant, 
          dynareParams, dr, kstate, ghx, ghu, aux, iv, ic, jcols, dfExt);

#ifdef DEBUG
        mexPrintf("Try CalcLikelihood\n");
#endif
#ifdef LL_TIMING_LOOP
        mexPrintf("DsgeLikelihood::CalcLikelihood: starting 1000 loops\n");
        for (int tt=0;tt<1000;++tt)
          {
#endif

    loglikelihood=dl.CalcLikelihood(xparam1);
#ifdef LL_TIMING_LOOP
          }
        mexPrintf("DsgeLikelihood::CalcLikelihood: finished 1000 loops\n");
#endif
/*****************************************************************
% OUTPUTS 
%   fval        :     value of the posterior kernel at xparam1.
%   cost_flag   :     zero if the function returns a penalty, one otherwise.
%   ys          :     steady state of original endogenous variables
%   trend_coeff :
%   info        :     vector of informations about the penalty:
%   vll         :     vector of time-step log-likelihoods at xparam1.
*****************************************************************/
#ifdef DEBUG
        mexPrintf("Try Outputs with nper=%d, loglikelihood = %f\n",nper,loglikelihood);
#endif
    if (nlhs >= 1)
      plhs[0] = mxCreateDoubleScalar(loglikelihood);
    if (nlhs >= 2)
      plhs[1] = mxCreateDoubleScalar((double)dl.getCostFlag());
    if (nlhs >= 3)
      {
        plhs[2] = mxCreateDoubleMatrix(endo_nbr, 1, mxREAL);
        Vector vss(mxGetPr(plhs[2]),endo_nbr);

#ifdef DEBUG
        mexPrintf("SteadyState  size %d \n", dl.getSteadyState().length());
       dl.getSteadyState().print() ;
        mexPrintf("Try getSteadyState into vss size %d \n", vss.length());
#endif
        vss= dl.getSteadyState();
      }
/*********************
    if (nlhs >= 4)
      plhs[3] = mxCreateDoubleScalar((double)dl.getCostFlag());
    if (nlhs >= 5)
      plhs[4] = mxCreateDoubleMatrix(numVarobs,1, mxREAL);//dummy trend_coeff
    if (nlhs >= 6)
      {
      // output full log-likelihood array
      // Set the output pointer to the  array of log likelihood. 
      std::vector<double>& vll=dl.getLikVector();
      plhs[5] = mxCreateDoubleMatrix(nper,1, mxREAL);
      double * mxll= mxGetPr(plhs[5]);
      // assign likelihood array
      for (int j=0;j<nper;++j)
        {
        mxll[j]=vll[j];
#ifdef DEBUG
        mexPrintf("mxll[%d]=%f  vll[%d]=%f\n",j, mxll[j], i, vll[j]);
#endif
        }
      }
*********************/
    }
  catch (const KordException &e)
    {
      printf("Caugth Kord exception: ");
      e.print();
      mexPrintf("Caugth Kord exception: %s", e.get_message());
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
      return; // 255;
    }
  catch (const ogu::Exception &e)
    {
      printf("Caught ogu::Exception: ");
      e.print();
      mexPrintf("Caugth general exception: %s", e.message());
      return; // 255;
    }  //catch
  }; // end of mexFunction()
}; // end of extern C
