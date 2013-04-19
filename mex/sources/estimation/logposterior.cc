/*
 * Copyright (C) 2010-2013 Dynare Team
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

#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>

#include "Vector.hh"
#include "Matrix.hh"
#include "LogPosteriorDensity.hh"

#include <dynmex.h>

class LogposteriorMexErrMsgTxtException
{
public:
  std::string errMsg;
  LogposteriorMexErrMsgTxtException(const std::string &msg) : errMsg(msg)
  {
  }
  inline const char *getErrMsg() { return errMsg.c_str(); }
};

void
fillEstParamsInfo(const mxArray *bayestopt_, const mxArray *estim_params_info, EstimatedParameter::pType type,
                  std::vector<EstimatedParameter> &estParamsInfo)
{
  const mxArray *bayestopt_ubp = mxGetField(bayestopt_, 0, "ub"); // upper bound
  const mxArray *bayestopt_lbp = mxGetField(bayestopt_, 0, "lb"); // lower bound
  const mxArray *bayestopt_p1p = mxGetField(bayestopt_, 0, "p1"); // prior mean
  const mxArray *bayestopt_p2p = mxGetField(bayestopt_, 0, "p2"); // prior standard deviation
  const mxArray *bayestopt_p3p = mxGetField(bayestopt_, 0, "p3"); // lower bound
  const mxArray *bayestopt_p4p = mxGetField(bayestopt_, 0, "p4"); // upper bound
  const mxArray *bayestopt_p6p = mxGetField(bayestopt_, 0, "p6"); // first hyper-parameter (\alpha for the BETA and GAMMA distributions, s for the INVERSE GAMMAs, expectation for the GAUSSIAN distribution, lower bound for the UNIFORM distribution).
  const mxArray *bayestopt_p7p = mxGetField(bayestopt_, 0, "p7"); // second hyper-parameter (\beta for the BETA and GAMMA distributions, \nu for the INVERSE GAMMAs, standard deviation for the GAUSSIAN distribution, upper bound for the UNIFORM distribution).
  const mxArray *bayestopt_jscalep = mxGetField(bayestopt_, 0, "jscale"); // MCMC jump scale

  const size_t bayestopt_size = mxGetM(bayestopt_);
  const VectorConstView bayestopt_ub(mxGetPr(bayestopt_ubp), bayestopt_size, 1);
  const VectorConstView bayestopt_lb(mxGetPr(bayestopt_lbp), bayestopt_size, 1);
  const VectorConstView bayestopt_p1(mxGetPr(bayestopt_p1p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p1");
  const VectorConstView bayestopt_p2(mxGetPr(bayestopt_p2p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p2");
  const VectorConstView bayestopt_p3(mxGetPr(bayestopt_p3p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p3");
  const VectorConstView bayestopt_p4(mxGetPr(bayestopt_p4p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p4");
  const VectorConstView bayestopt_p6(mxGetPr(bayestopt_p6p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p6");
  const VectorConstView bayestopt_p7(mxGetPr(bayestopt_p7p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p7");
  const VectorConstView bayestopt_jscale(mxGetPr(bayestopt_jscalep), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "jscale");

  // loop processsing
  size_t m = mxGetM(estim_params_info), n = mxGetN(estim_params_info);
  MatrixConstView epi(mxGetPr(estim_params_info), m, n, m);
  size_t bayestopt_count = estParamsInfo.size();

  for (size_t i = 0; i < m; i++)
    {
      size_t col = 0;
      size_t id1 = (size_t) epi(i, col++) - 1;
      size_t id2 = 0;
      if (type == EstimatedParameter::shock_Corr
          || type == EstimatedParameter::measureErr_Corr)
        id2 = (size_t) epi(i, col++) - 1;
      col++; // Skip init_val #2 or #3
      double par_low_bound =  bayestopt_lb(bayestopt_count); col++; //#3 epi(i, col++);
      double par_up_bound =  bayestopt_ub(bayestopt_count); col++; //#4 epi(i, col++);
      Prior::pShape shape = (Prior::pShape) epi(i, col++);
      double mean = epi(i, col++);
      double std = epi(i, col++);
      double low_bound =  bayestopt_p3(bayestopt_count);
      double up_bound =  bayestopt_p4(bayestopt_count);
      double fhp =  bayestopt_p6(bayestopt_count); // double p3 = epi(i, col++);
      double shp =  bayestopt_p7(bayestopt_count); // double p4 = epi(i, col++);

      Prior *p = Prior::constructPrior(shape, mean, std, low_bound, up_bound, fhp, shp); //1.0,INFINITY);//p3, p4);

      // Only one subsample
      std::vector<size_t> subSampleIDs;
      subSampleIDs.push_back(0);
      estParamsInfo.push_back(EstimatedParameter(type, id1, id2, subSampleIDs,
                                                 par_low_bound, par_up_bound, p));
      bayestopt_count++;
    }
}

template <class VEC1, class VEC2>
double
logposterior(VEC1 &estParams, const MatrixConstView &data,
             const mxArray *options_, const mxArray *M_, const mxArray *estim_params_,
	     const mxArray *bayestopt_, const mxArray *oo_, VEC2 &steadyState, double *trend_coeff,
	     VectorView &deepParams, Matrix &H, MatrixView &Q)
{
  double loglinear = *mxGetPr(mxGetField(options_, 0, "loglinear"));
  if (loglinear == 1)
    throw LogposteriorMexErrMsgTxtException("Option loglinear is not supported");

  if (*mxGetPr(mxGetField(options_, 0, "endogenous_prior")) == 1)
    throw LogposteriorMexErrMsgTxtException("Option endogenous_prior is not supported");

  double with_trend = *mxGetPr(mxGetField(bayestopt_, 0, "with_trend"));
  if (with_trend == 1)
    throw LogposteriorMexErrMsgTxtException("Observation trends are not supported");

  // Construct arguments of constructor of LogLikelihoodMain
  char *fName = mxArrayToString(mxGetField(M_, 0, "fname"));
  std::string basename(fName);
  mxFree(fName);

  size_t n_endo = (size_t) *mxGetPr(mxGetField(M_, 0, "endo_nbr"));
  size_t n_exo = (size_t) *mxGetPr(mxGetField(M_, 0, "exo_nbr"));
  size_t n_param = (size_t) *mxGetPr(mxGetField(M_, 0, "param_nbr"));
  size_t n_estParams = estParams.getSize();

  std::vector<size_t> zeta_fwrd, zeta_back, zeta_mixed, zeta_static;
  const mxArray *lli_mx = mxGetField(M_, 0, "lead_lag_incidence");
  MatrixConstView lli(mxGetPr(lli_mx), mxGetM(lli_mx), mxGetN(lli_mx), mxGetM(lli_mx));
  if (lli.getRows() != 3)
    throw LogposteriorMexErrMsgTxtException("Purely backward or purely forward models are not supported");
  if (lli.getCols() != n_endo)
    throw LogposteriorMexErrMsgTxtException("Incorrect lead/lag incidence matrix");

  for (size_t i = 0; i < n_endo; i++)
    {
      if (lli(0, i) == 0 && lli(2, i) == 0)
        zeta_static.push_back(i);
      else if (lli(0, i) != 0 && lli(2, i) == 0)
        zeta_back.push_back(i);
      else if (lli(0, i) == 0 && lli(2, i) != 0)
        zeta_fwrd.push_back(i);
      else
        zeta_mixed.push_back(i);
    }

  double qz_criterium = *mxGetPr(mxGetField(options_, 0, "qz_criterium"));
  double lyapunov_tol = *mxGetPr(mxGetField(options_, 0, "lyapunov_complex_threshold"));
  double riccati_tol = *mxGetPr(mxGetField(options_, 0, "riccati_tol"));
  size_t presample = (size_t) *mxGetPr(mxGetField(options_, 0, "presample"));

  std::vector<size_t> varobs;
  const mxArray *varobs_mx = mxGetField(options_, 0, "varobs_id");
  if (mxGetM(varobs_mx) != 1)
    throw LogposteriorMexErrMsgTxtException("options_.varobs_id must be a row vector");

  size_t n_varobs = mxGetN(varobs_mx);
  // substract 1.0 from obsverved variables index
  std::transform(mxGetPr(varobs_mx), mxGetPr(varobs_mx) + n_varobs, back_inserter(varobs),
                 std::bind2nd(std::minus<size_t>(), 1));

  if (data.getRows() != n_varobs)
    throw LogposteriorMexErrMsgTxtException("Data does not have as many rows as there are observed variables");

  std::vector<EstimationSubsample> estSubsamples;
  estSubsamples.push_back(EstimationSubsample(0, data.getCols() - 1));

  std::vector<EstimatedParameter> estParamsInfo;
  fillEstParamsInfo(bayestopt_, mxGetField(estim_params_, 0, "var_exo"), EstimatedParameter::shock_SD,
                    estParamsInfo);
  fillEstParamsInfo(bayestopt_, mxGetField(estim_params_, 0, "var_endo"), EstimatedParameter::measureErr_SD,
                    estParamsInfo);
  fillEstParamsInfo(bayestopt_, mxGetField(estim_params_, 0, "corrx"), EstimatedParameter::shock_Corr,
                    estParamsInfo);
  fillEstParamsInfo(bayestopt_, mxGetField(estim_params_, 0, "corrn"), EstimatedParameter::measureErr_Corr,
                    estParamsInfo);
  fillEstParamsInfo(bayestopt_, mxGetField(estim_params_, 0, "param_vals"), EstimatedParameter::deepPar,
                    estParamsInfo);

  EstimatedParametersDescription epd(estSubsamples, estParamsInfo);

  bool noconstant = (bool) *mxGetPr(mxGetField(options_, 0, "noconstant"));

  // Allocate LogPosteriorDensity object
  LogPosteriorDensity lpd(basename, epd, n_endo, n_exo, zeta_fwrd, zeta_back, zeta_mixed, zeta_static,
                          qz_criterium, varobs, riccati_tol, lyapunov_tol, noconstant);

  // Construct arguments of compute() method

  // Compute the posterior
  double logPD = lpd.compute(steadyState, estParams, deepParams, data, Q, H, presample);

  // Cleanups
  for (std::vector<EstimatedParameter>::iterator it = estParamsInfo.begin();
       it != estParamsInfo.end(); it++)
    delete it->prior;

  return logPD;
}

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  if (nrhs != 7 )
    DYN_MEX_FUNC_ERR_MSG_TXT("logposterior: exactly 7 input arguments are required.");

  if (nlhs > 9 )
    DYN_MEX_FUNC_ERR_MSG_TXT("logposterior returns 8 output arguments at the most.");

  // Check and retrieve the RHS arguments

  if (!mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 1)
    DYN_MEX_FUNC_ERR_MSG_TXT("logposterior: First argument must be a column vector of double-precision numbers");

  VectorConstView estParams(mxGetPr(prhs[0]), mxGetM(prhs[0]), 1);

  for (int i = 1; i < 7; ++i)
    if (!mxIsStruct(prhs[i]))
      {
	std::stringstream msg;
	msg << "logposterior: argument " << i+1 << " must be a Matlab structure";
	DYN_MEX_FUNC_ERR_MSG_TXT(msg.str().c_str());
      }

  const mxArray *dataset = prhs[1];
  const mxArray *options_ = prhs[2];
  const mxArray *M_ = prhs[3];
  const mxArray *estim_params_ = prhs[4];
  const mxArray *bayestopt_ = prhs[5];
  const mxArray *oo_ = prhs[6];

  const mxArray *dataset_data = mxGetField(dataset,0,"data");
  MatrixConstView data(mxGetPr(dataset_data), mxGetM(dataset_data), mxGetN(dataset_data), mxGetM(dataset_data));

  // Creaete LHS arguments

  size_t endo_nbr = (size_t) *mxGetPr(mxGetField(M_, 0, "endo_nbr"));
  size_t exo_nbr = (size_t) *mxGetPr(mxGetField(M_, 0, "exo_nbr"));
  size_t param_nbr = (size_t) *mxGetPr(mxGetField(M_, 0, "param_nbr"));
  size_t varobs_nbr = mxGetM(mxGetField(options_, 0, "varobs"));
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(endo_nbr, 1, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(varobs_nbr, 1, mxREAL);
  plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
  plhs[5] = mxCreateDoubleMatrix(param_nbr, 1, mxREAL);
  plhs[6] = mxCreateDoubleMatrix(varobs_nbr, varobs_nbr, mxREAL);
  plhs[7] = mxCreateDoubleMatrix(exo_nbr, exo_nbr, mxREAL);
  double *lik = mxGetPr(plhs[0]);
  double *exit_flag = mxGetPr(plhs[1]);

  VectorView steadyState(mxGetPr(mxGetField(oo_,0,"steady_state")),endo_nbr, 1);
  VectorView deepParams(mxGetPr(mxGetField(M_, 0, "params")),param_nbr,1);

  MatrixView Q(mxGetPr(mxGetField(M_, 0, "Sigma_e")), exo_nbr, exo_nbr, exo_nbr);

  Matrix H(varobs_nbr,varobs_nbr);
  const mxArray *H_mx = mxGetField(M_, 0, "H");
  if (mxGetM(H_mx) == 1 && mxGetN(H_mx) == 1 && *mxGetPr(H_mx) == 0)
    H.setAll(0.0);
  else
    H = MatrixConstView(mxGetPr(H_mx), varobs_nbr, varobs_nbr, varobs_nbr);

  double *trend_coeff  = mxGetPr(plhs[3]);
  double *info_mx  = mxGetPr(plhs[4]);

  // Compute and return the value
  try
    {
      *lik = logposterior(estParams, data, options_, M_, estim_params_, bayestopt_, oo_,
				steadyState, trend_coeff, deepParams, H, Q);
      *info_mx = 0;
      *exit_flag = 0;
    }
  catch (LogposteriorMexErrMsgTxtException e)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(e.getErrMsg());
    }
  catch (SteadyStateSolver::SteadyStateException e)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(e.message.c_str());
    }
}
