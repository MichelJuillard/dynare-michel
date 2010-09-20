/*
 * Copyright (C) 2010 Dynare Team
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

#include "Vector.hh"
#include "Matrix.hh"
#include "LogPosteriorDensity.hh"

#include "dynmex.h"
#include "mex.h"

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
fillEstParamsInfo(const mxArray *estim_params_info, EstimatedParameter::pType type,
                  std::vector<EstimatedParameter> &estParamsInfo)
{
  // execute once only
  static const mxArray *bayestopt_ = mexGetVariablePtr("global", "bayestopt_");
  static const mxArray *bayestopt_ubp = mxGetField(bayestopt_, 0, "ub"); // upper bound
  static const mxArray *bayestopt_lbp = mxGetField(bayestopt_, 0, "lb"); // lower bound
  static const mxArray *bayestopt_p1p = mxGetField(bayestopt_, 0, "p1"); // prior mean
  static const mxArray *bayestopt_p2p = mxGetField(bayestopt_, 0, "p2"); // prior standard deviation
  static const mxArray *bayestopt_p3p = mxGetField(bayestopt_, 0, "p3"); // lower bound
  static const mxArray *bayestopt_p4p = mxGetField(bayestopt_, 0, "p4"); // upper bound
  static const mxArray *bayestopt_p6p = mxGetField(bayestopt_, 0, "p6"); // first hyper-parameter (\alpha for the BETA and GAMMA distributions, s for the INVERSE GAMMAs, expectation for the GAUSSIAN distribution, lower bound for the UNIFORM distribution).
  static const mxArray *bayestopt_p7p = mxGetField(bayestopt_, 0, "p7"); // second hyper-parameter (\beta for the BETA and GAMMA distributions, \nu for the INVERSE GAMMAs, standard deviation for the GAUSSIAN distribution, upper bound for the UNIFORM distribution).
  static const mxArray *bayestopt_jscalep = mxGetField(bayestopt_, 0, "jscale"); // MCMC jump scale

  static const size_t bayestopt_size = mxGetM(bayestopt_);
  static const VectorConstView bayestopt_ub(mxGetPr(bayestopt_ubp), bayestopt_size, 1);
  static const VectorConstView bayestopt_lb(mxGetPr(bayestopt_lbp), bayestopt_size, 1);
  static const VectorConstView bayestopt_p1(mxGetPr(bayestopt_p1p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p1");
  static const VectorConstView bayestopt_p2(mxGetPr(bayestopt_p2p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p2");
  static const VectorConstView bayestopt_p3(mxGetPr(bayestopt_p3p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p3");
  static const VectorConstView bayestopt_p4(mxGetPr(bayestopt_p4p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p4");
  static const VectorConstView bayestopt_p6(mxGetPr(bayestopt_p6p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p6");
  static const VectorConstView bayestopt_p7(mxGetPr(bayestopt_p7p), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "p7");
  static const VectorConstView bayestopt_jscale(mxGetPr(bayestopt_jscalep), bayestopt_size, 1); //=mxGetField(bayestopt_, 0, "jscale");

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

double
logposterior(const VectorConstView &estParams, const MatrixConstView &data,
             const std::string &mexext)
{
  // Retrieve pointers to global variables
  const mxArray *M_ = mexGetVariablePtr("global", "M_");
  const mxArray *oo_ = mexGetVariablePtr("global", "oo_");
  const mxArray *options_ = mexGetVariablePtr("global", "options_");
  const mxArray *estim_params_ = mexGetVariablePtr("global", "estim_params_");

  // Construct arguments of constructor of LogLikelihoodMain
  char *fName = mxArrayToString(mxGetField(M_, 0, "fname"));
  std::string dynamicDllFile(fName);
  mxFree(fName);
  dynamicDllFile += "_dynamic." + mexext;

  size_t n_endo = (size_t) *mxGetPr(mxGetField(M_, 0, "endo_nbr"));
  size_t n_exo = (size_t) *mxGetPr(mxGetField(M_, 0, "exo_nbr"));
  size_t n_param = (size_t) *mxGetPr(mxGetField(M_, 0, "param_nbr"));
  size_t n_estParams = estParams.getSize();

  std::vector<size_t> zeta_fwrd, zeta_back, zeta_mixed, zeta_static;
  const mxArray *lli_mx = mxGetField(M_, 0, "lead_lag_incidence");
  MatrixConstView lli(mxGetPr(lli_mx), mxGetM(lli_mx), mxGetN(lli_mx), mxGetM(lli_mx));
  if (lli.getRows() != 3 || lli.getCols() != n_endo)
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
  std::transform(mxGetPr(varobs_mx), mxGetPr(varobs_mx) + n_varobs, back_inserter(varobs),
                 std::bind2nd(std::minus<size_t>(), 1));

  if (data.getRows() != n_varobs)
    throw LogposteriorMexErrMsgTxtException("Data has not as many rows as there are observed variables");

  std::vector<EstimationSubsample> estSubsamples;
  estSubsamples.push_back(EstimationSubsample(0, data.getCols() - 1));

  std::vector<EstimatedParameter> estParamsInfo;
  fillEstParamsInfo(mxGetField(estim_params_, 0, "var_exo"), EstimatedParameter::shock_SD,
                    estParamsInfo);
  fillEstParamsInfo(mxGetField(estim_params_, 0, "var_endo"), EstimatedParameter::measureErr_SD,
                    estParamsInfo);
  fillEstParamsInfo(mxGetField(estim_params_, 0, "corrx"), EstimatedParameter::shock_Corr,
                    estParamsInfo);
  fillEstParamsInfo(mxGetField(estim_params_, 0, "corrn"), EstimatedParameter::measureErr_Corr,
                    estParamsInfo);
  fillEstParamsInfo(mxGetField(estim_params_, 0, "param_vals"), EstimatedParameter::deepPar,
                    estParamsInfo);

  EstimatedParametersDescription epd(estSubsamples, estParamsInfo);

  // Allocate LogPosteriorDensity object
  int info;
  LogPosteriorDensity lpd(dynamicDllFile, epd, n_endo, n_exo, zeta_fwrd, zeta_back, zeta_mixed, zeta_static,
      qz_criterium, varobs, riccati_tol, lyapunov_tol, info);

  // Construct arguments of compute() method
  Matrix steadyState(n_endo, 1);
  mat::get_col(steadyState, 0) = VectorConstView(mxGetPr(mxGetField(oo_, 0, "steady_state")), n_endo, 1);

  Vector estParams2(n_estParams);
  estParams2 = estParams;
  Vector deepParams(n_param);
  deepParams = VectorConstView(mxGetPr(mxGetField(M_, 0, "params")), n_param, 1);
  Matrix Q(n_exo);
  Q = MatrixConstView(mxGetPr(mxGetField(M_, 0, "Sigma_e")), n_exo, n_exo, n_exo);
  
  Matrix H(n_varobs);
  const mxArray *H_mx = mxGetField(M_, 0, "H");
  if (mxGetM(H_mx) == 1 && mxGetN(H_mx) == 1 && *mxGetPr(H_mx) == 0)
    H.setAll(0.0);
  else
    H = MatrixConstView(mxGetPr(mxGetField(M_, 0, "H")), n_varobs, n_varobs, n_varobs);

  // Compute the posterior
  double logPD =lpd.compute(steadyState, estParams2, deepParams, data, Q, H, presample, info);

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
  if (nrhs != 3 || nlhs !=2)
    DYN_MEX_FUNC_ERR_MSG_TXT("logposterior: exactly three input arguments and two output arguments are required.");

  // Check and retrieve the arguments

  if (!mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 1)
    DYN_MEX_FUNC_ERR_MSG_TXT("logposterior: First argument must be a column vector of double-precision numbers");

  VectorConstView estParams(mxGetPr(prhs[0]), mxGetM(prhs[0]), 1);

  if (!mxIsDouble(prhs[1]))
    DYN_MEX_FUNC_ERR_MSG_TXT("logposterior: Second argument must be a matrix of double-precision numbers");

  MatrixConstView data(mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]), mxGetM(prhs[1]));

  if (!mxIsChar(prhs[2]))
    DYN_MEX_FUNC_ERR_MSG_TXT("logposterior: Third argument must be a character string");

  char *mexext_mx = mxArrayToString(prhs[2]);
  std::string
  mexext(mexext_mx);
  mxFree(mexext_mx);

  // Compute and return the value
  try
    {
      double lik = logposterior(estParams, data, mexext);
      plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
      *mxGetPr(plhs[1]) = lik;
    }
  catch (LogposteriorMexErrMsgTxtException e)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(e.getErrMsg());
    }
  plhs[0] = mxCreateDoubleScalar(0);
}
