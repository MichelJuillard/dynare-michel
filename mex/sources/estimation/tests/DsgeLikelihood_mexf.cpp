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
*****************************************************************%function [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)
% function [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)
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

#ifdef MATLAB
#include "mexutils.h"
#endif

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
    mexPrintf("k_order_perturbation: mexExt=%s.\n", dfExt);
#endif
/***********
***************/
  int numPeriods=1;

  //MexStruct dynareParams();
  GeneralParams& dynareParams=*(new MexStruct());
  //MexStructParam& dr=dynareParams.getStructField("dr");
  GeneralParams& dr=dynareParams.getStructField("dr");
  
  vector<int>&mfys=dynareParams.getIntVectorField("mfys");
  vector<int>&mf=dynareParams.getIntVectorField("mf1");
  int numVarobs=data.numRows();
  Vector& SteadyState=dr.getDoubleVectorField(string("ys"));
  Vector constant(numVarobs);//=*(new Vector(numVarobs));//   = *(new Vector(nobs));
  GeneralMatrix&kstate = dr.getMatrixField(string("kstate"));
  vector<int>&order_var = dr.getIntVectorField(string("order_var"));
  int order=(int)dynareParams.getDoubleField(string("order"));
  int endo_nbr = (int)dynareParams.getDoubleField(string("endo_nbr"));
  int exo_nbr = (int)dynareParams.getDoubleField(string("exo_nbr"));
  int nstatic = (int)dr.getDoubleField(string("nstatic"));
  int npred = (int)dr.getDoubleField(string("npred"));
  int nfwrd = (int)dr.getDoubleField(string("nfwrd"));
  Vector& ub=dynareParams.getDoubleVectorField(string("ub"));
  Vector& lb=dynareParams.getDoubleVectorField(string("lb"));
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

  GeneralMatrix& aux= dr.getMatrixField(string("transition_auxiliary_variables"));
  int nr=endo_nbr+aux.numRows();
  Vector& a_init=*(new Vector(numVarobs));

  GeneralMatrix& Q = dynareParams.getMatrixField(string("Sigma_e"));
  GeneralMatrix& H = dynareParams.getMatrixField(string("H"));
  GeneralMatrix Y(data.numRows(),data.numCols());
  GeneralMatrix T(nr,nr);
  GeneralMatrix Z(numVarobs,nr);
  GeneralMatrix Pstar(nr,nr);
  GeneralMatrix R(nr,exo_nbr);
  GeneralMatrix ghx(endo_nbr,jcols-exo_nbr);
  GeneralMatrix ghu(endo_nbr,exo_nbr);

//Pinf=[]
  GeneralMatrix Pinf (nr,nr);
  Pinf.zeros();

  try
    {

    DsgeLikelihood dl( a_init, Q, R,T, Z, Pstar, Pinf, H,data,Y,  
          numPeriods, //  const int INnumVarobs, //  const int INnumTimeObs,
          order, endo_nbr, exo_nbr, nstatic, npred,nfwrd, num_of_observations, 
          no_more_missing_observations, order_var, mfys, mf, xparam1,
          num_dp, deepParams, ub, lb, pshape, p6, p7, p3, p4, SteadyState, constant, 
          dynareParams, dr, kstate, ghx, ghu, jcols, dfExt);

    double loglikelihood=dl.CalcLikelihood(xparam1);

/*****************************************************************
% OUTPUTS 
%   fval        :     value of the posterior kernel at xparam1.
%   cost_flag   :     zero if the function returns a penalty, one otherwise.
%   ys          :     steady state of original endogenous variables
%   trend_coeff :
%   info        :     vector of informations about the penalty:
%   vll         :     vector of time-step log-likelihoods at xparam1.
*****************************************************************/
    if (nlhs >= 1)
      plhs[0] = mxCreateDoubleScalar(loglikelihood);
    if (nlhs >= 2)
      plhs[1] = mxCreateDoubleScalar((double)dl.getCostFlag());
    if (nlhs >= 3)
      {
        plhs[2] = mxCreateDoubleMatrix(endo_nbr, 1, mxREAL);
        Vector vss(mxGetPr(plhs[2]),endo_nbr);
        vss= dl.getSteadyState();
      }
    if (nlhs >= 4)
         plhs[3] = mxCreateDoubleScalar((double)dl.getCostFlag());
    if (nlhs >= 5)
        plhs[4] = mxCreateDoubleMatrix(numVarobs,1, mxREAL);//dummy trend_coeff
    if (nlhs >= 6)
      {
      // output full log-likelihood array
      /* Set the output pointer to the  array of log likelihood. */
      std::vector<double>& vll=dl.getLikVector();
      plhs[5] = mxCreateDoubleMatrix(nper,1, mxREAL);
      double * mxll= mxGetPr(plhs[5]);
      // assign likelihood array
      for (int j=0;j<nper;++j)
        mxll[j]=vll[j];
      }
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
