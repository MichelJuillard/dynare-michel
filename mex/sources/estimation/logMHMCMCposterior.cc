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
#include "RandomWalkMetropolisHastings.hh"

#include "mex.h"
#include "mat.h"

#if defined(_WIN32) || defined(__CYGWIN32__) || defined(WINDOWS)
# define DIRECTORY_SEPARATOR "\\"
#else
# define DIRECTORY_SEPARATOR "/"
#endif

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

size_t
sampleMHMC(LogPosteriorDensity &lpd, RandomWalkMetropolisHastings &rwmh,
           Matrix &steadyState, Vector &estParams, Vector &deepParams, const MatrixConstView &data, Matrix &Q, Matrix &H,
           const size_t presampleStart, int &info, const VectorConstView &nruns, const size_t fblock, const size_t nBlocks,
           const Matrix &Jscale, Prior &drawDistribution, EstimatedParametersDescription &epd, 
           const std::string &resultsFileStem, const size_t console_mode)
{
  enum{iMin,iMax};
  std::vector<size_t> OpenOldFile(nBlocks, 1);
  size_t jloop = 0, irun, j; // counters
  double dsum, dmax, dmin, sux = 0, jsux=0;
  std::string mhFName;
  std::stringstream ssFName;
  MATFile *drawmat; // MCMC draws output file pointer
  FILE *fidlog;  // log file
  int matfStatus;
  size_t npar = estParams.getSize();
  Matrix MinMax(npar,2);
  const mxArray *myinputs = mexGetVariablePtr("caller", "myinputs");
  const mxArray *InitSizeArrayPtr = mxGetField(myinputs, 0, "InitSizeArray");
  const VectorConstView InitSizeArrayVw(mxGetPr(InitSizeArrayPtr), nBlocks, 1);
  Vector InitSizeArray(InitSizeArrayVw.getSize());
  InitSizeArray = InitSizeArrayVw;
  const mxArray *flinePtr = mxGetField(myinputs, 0, "fline");
  const VectorConstView fline(mxGetPr(flinePtr), nBlocks, 1);

  const mxArray *NewFileArrayPtr = mxGetField(myinputs, 0, "NewFile");
  VectorView NewFileVw(mxGetPr(NewFileArrayPtr), nBlocks, 1);
  Vector NewFile(NewFileVw.getSize());
  NewFile = NewFileVw;

  const mxArray *MAX_nrunsPtr = mxGetField(myinputs, 0, "MAX_nruns");
  const size_t MAX_nruns = (size_t) mxGetScalar(MAX_nrunsPtr);

  const mxArray *record = mxGetField(myinputs, 0, "record");
  const mxArray *AcceptationRatesPtr = mxGetField(record, 0, "AcceptationRates");
  VectorView AcceptationRates(mxGetPr(AcceptationRatesPtr), nBlocks, 1);

  mxArray *mxLastParametersPtr= mxGetField(record, 0, "LastParameters");
  MatrixView LastParameters(mxGetPr(mxLastParametersPtr), nBlocks, npar, nBlocks);

  mxArray *mxLastLogLikPtr = mxGetField(record, 0, "LastLogLiK");
  VectorView LastLogLiK(mxGetPr(mxLastLogLikPtr), nBlocks, 1);

  mxArray *mxMhLogPostDensPtr = 0;
  mxArray *mxMhParamDrawsPtr = 0;
  size_t currInitSizeArray = 0;

  mexPrintf("\n Starting MH Block Loop\n\n");

  for (size_t b = fblock; b <= nBlocks; ++b)
    {
      jloop = jloop+1;

      sux = 0.0;
      jsux = 0;
      irun = (size_t) fline(b-1);
      j = 0;//1;
      while (j < nruns(b-1))
        {
        if (currInitSizeArray != (size_t) InitSizeArray(b-1))
          {
            // new or different size result arrays/matrices
            currInitSizeArray = (size_t) InitSizeArray(b-1);
            if (mxMhLogPostDensPtr) 
              mxDestroyArray(mxMhLogPostDensPtr);  // log post density array
            mxMhLogPostDensPtr = mxCreateDoubleMatrix(currInitSizeArray, 1, mxREAL);
            if (mxMhParamDrawsPtr) 
              mxDestroyArray(mxMhParamDrawsPtr);  // accepted MCMC MH draws
            mxMhParamDrawsPtr =  mxCreateDoubleMatrix( currInitSizeArray, npar,  mxREAL);
          }
        VectorView mhLogPostDens(mxGetPr(mxMhLogPostDensPtr), currInitSizeArray, (size_t) 1);
        MatrixView mhParamDraws(mxGetPr(mxMhParamDrawsPtr), currInitSizeArray, npar, currInitSizeArray);
        jsux= rwmh.compute(mhLogPostDens, mhParamDraws, steadyState, estParams, deepParams, data, Q, H,
                   presampleStart, info, currInitSizeArray, Jscale, lpd, drawDistribution, epd);
        sux+=jsux*currInitSizeArray;
        j += currInitSizeArray; //j=j+1;
        irun += currInitSizeArray;
        
        if(console_mode)
             mexPrintf("   MH: Computing Metropolis-Hastings (chain %d/%d): %3.f \b%% done, acceptance rate: %3.f \b%%\r", b, nBlocks, 100 * j/nruns(b-1), 100 * sux / j);
        // % Now I save the simulations
        // save draw  2 mat file ([MhDirectoryName '/' ModelName '_mh' int2str(NewFile(b)) '_blck' int2str(b) '.mat'],'x2','logpo2');
        ssFName.clear();
        ssFName.str("");
        ssFName << resultsFileStem << DIRECTORY_SEPARATOR << "metropolis" << DIRECTORY_SEPARATOR << resultsFileStem << "_mh" << (size_t) NewFile(b-1) << "_blck" << b << ".mat" ;
        mhFName = ssFName.str();
        drawmat = matOpen(mhFName.c_str(), "w");
        if (drawmat==0)
          {
            mexPrintf("Error in MH: Can not open draws Mat file for writing:  %s \n", mhFName.c_str());
            exit(1);
          }
        matfStatus = matPutVariable(drawmat, "x2", mxMhParamDrawsPtr);
        if (matfStatus)
          {
            mexPrintf("Error in MH: Can not use draws Mat file for writing:  %s \n", mhFName.c_str());
            exit(1);
          }
        matfStatus = matPutVariable(drawmat, "logpo2", mxMhLogPostDensPtr);
        if (matfStatus)
          {
            mexPrintf("Error in MH: Can not usee draws Mat file for writing:  %s \n", mhFName.c_str());
            exit(1);
          }
        matClose(drawmat);

        // save log to fidlog = fopen([MhDirectoryName '/metropolis.log'],'a');
        ssFName.str("");
        ssFName << resultsFileStem << DIRECTORY_SEPARATOR << "metropolis" << DIRECTORY_SEPARATOR << "metropolis.log" ;
        mhFName = ssFName.str();
        fidlog = fopen(mhFName.c_str(), "a");
        fprintf(fidlog,"\n");
        fprintf(fidlog,"%% Mh%dBlck%d ( %s %s )\n", (int) NewFile(b-1),b , __DATE__ , __TIME__  );
        fprintf(fidlog," \n");
        fprintf(fidlog,"  Number of simulations.: %d \n", currInitSizeArray);// (length(logpo2)) ');
        fprintf(fidlog,"  Acceptation rate......: %f \n", jsux/currInitSizeArray);
        fprintf(fidlog,"  Posterior mean........:\n");
        for (size_t i=0; i<npar; ++i)
          {
            VectorView mhpdColVw=mat::get_col(mhParamDraws,i); 
            fprintf(fidlog,"    params: %d : %f \n", i, vec::meanSumMinMax(dsum, dmin, dmax, mhpdColVw));
            MinMax(i,iMin)=dmin;
            MinMax(i,iMax)=dmax;
          } // end
        fprintf(fidlog,"    log2po: %f \n", vec::meanSumMinMax(dsum, dmin, dmax, mhLogPostDens));
        fprintf(fidlog,"  Minimum value.........:\n");;
        for (size_t i=0; i<npar; ++i)
          fprintf(fidlog,"    params: %d : %f \n", i, MinMax(i,iMin));
        fprintf(fidlog,"    log2po: %f \n", dmin);
        fprintf(fidlog,"  Maximum value.........:\n");
        for (size_t i=0; i<npar; ++i)
          fprintf(fidlog,"    params: %d : %f \n", i, MinMax(i,iMax));
        fprintf(fidlog,"    log2po: %f \n", dmax);
        fprintf(fidlog," \n");
        fclose(fidlog);

        jsux = 0;
  //            if (j == nruns(b-1)) // % I record the last draw...
  //                record.LastParameters(b,:) = x2(end,:);
  //                record.LastLogLiK(b) = logpo2(end);
  //            } // end
        if (j == nruns(b-1)) // % I record the last draw...
          {
            VectorView LastParametersRow= mat::get_row(LastParameters,b-1); 
            LastParametersRow= mat::get_row(mhParamDraws, currInitSizeArray-1);//x2(end,:);
            LastLogLiK(b-1) = mhLogPostDens(currInitSizeArray-1); //logpo2(end);
          } // end
        // size of next file in chain b
        InitSizeArray(b-1) = std::min((size_t) nruns(b-1)-j, MAX_nruns);
        // initialization of next file if necessary
        if (InitSizeArray(b-1))
          {
            NewFile(b-1) = NewFile(b-1) + 1;
            irun = 0;
          } // end

        } // end  while % End of the simulations for one mh-block.
      //set record.
      AcceptationRates(b-1) = sux/j;
    } // end % End of the loop over the mh-blocks.
    if( mexPutVariable("caller", "record.AcceptationRates", AcceptationRatesPtr))
      mexPrintf("MH Warning: due to error record.AcceptationRates is NOT set !! \n") ;

    if( mexPutVariable("caller", "record.LastParameters", mxLastParametersPtr))
      mexPrintf("MH Warning: due to error record.MhParamDraw is NOT set !! \n") ;

    if( mexPutVariable("caller", "record.LastLogLiK", mxLastLogLikPtr))
      mexPrintf("MH Warning: due to error record.LastLogLiK is NOT set !! \n") ;

    NewFileVw=NewFile;
    if( mexPutVariable("caller", "NewFile", NewFileArrayPtr))
      mexPrintf("MH Warning: due to error NewFile is NOT set !! \n") ;

    if (mxMhLogPostDensPtr) 
      mxDestroyArray(mxMhLogPostDensPtr);  // delete log post density array
    if (mxMhParamDrawsPtr) 
      mxDestroyArray(mxMhParamDrawsPtr);  // delete accepted MCMC MH draws

    // return last line run in the last MH block sub-array 
    return currInitSizeArray;

}

size_t
logMCMCposterior(const VectorConstView &estParams, const MatrixConstView &data, const std::string &mexext,
                 const size_t fblock, const size_t nBlocks, const VectorConstView &nMHruns, const MatrixConstView &D)
{
  // Retrieve pointers to global variables
  const mxArray *M_ = mexGetVariablePtr("global", "M_");
  const mxArray *oo_ = mexGetVariablePtr("global", "oo_");
  const mxArray *options_ = mexGetVariablePtr("global", "options_");
  const mxArray *estim_params_ = mexGetVariablePtr("global", "estim_params_");

  // Construct arguments of constructor of LogLikelihoodMain
  char *fName = mxArrayToString(mxGetField(M_, 0, "fname"));
  std::string resultsFileStem(fName);
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
    mexErrMsgTxt("Incorrect lead/lag incidence matrix");
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
  size_t console_mode = (size_t) *mxGetPr(mxGetField(options_, 0, "console_mode"));

  std::vector<size_t> varobs;
  const mxArray *varobs_mx = mxGetField(options_, 0, "varobs_id");
  if (mxGetM(varobs_mx) != 1)
    mexErrMsgTxt("options_.varobs_id must be a row vector");
  size_t n_varobs = mxGetN(varobs_mx);
  std::transform(mxGetPr(varobs_mx), mxGetPr(varobs_mx) + n_varobs, back_inserter(varobs),
                 std::bind2nd(std::minus<size_t>(), 1));

  if (data.getRows() != n_varobs)
    mexErrMsgTxt("Data has not as many rows as there are observed variables");

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

  // Construct MHMCMC Sampler
  RandomWalkMetropolisHastings rwmh(estParams2.getSize());
  // Construct GaussianPrior drawDistribution m=0, sd=1
  GaussianPrior drawGaussDist01(0.0, 1.0, -INFINITY, INFINITY, 0.0, 1.0);
  // get Jscale = diag(bayestopt_.jscale);
  const mxArray *bayestopt_ = mexGetVariablePtr("global", "bayestopt_");
  Matrix Jscale(n_estParams);
  Matrix Dscale(n_estParams);
  //Vector vJscale(n_estParams);
  Jscale.setAll(0.0);
  VectorConstView vJscale(mxGetPr(mxGetField(bayestopt_, 0, "jscale")), n_estParams, 1);
  for (size_t i = 0; i < n_estParams; i++)
    Jscale(i, i) = vJscale(i);
  blas::gemm("N", "N", 1.0, D, Jscale, 0.0, Dscale);

//  Matrix mh_bounds(n_estParams,2);

  // Compute the MHMCMC loop draws
  // and get get last line run in the last MH block sub-array 
  size_t lastMHblockArrayLine = sampleMHMC(lpd, rwmh, steadyState, estParams2, deepParams, data, Q, H, presample, info,
             nMHruns, fblock, nBlocks, Dscale, drawGaussDist01, epd, resultsFileStem, console_mode); 

  // Cleanups
  for (std::vector<EstimatedParameter>::iterator it = estParamsInfo.begin();
       it != estParamsInfo.end(); it++)
    delete it->prior;

  return lastMHblockArrayLine;
  
}

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  if (nrhs != 7)
    mexErrMsgTxt("logposterior: exactly seven arguments are required.");
  if (nlhs != 1)
    mexErrMsgTxt("logposterior: exactly one return argument is required.");

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
// Check and retrieve the arguments

  if (!mxIsDouble(prhs[0]) || mxGetN(prhs[0]) != 1)
    mexErrMsgTxt("logposterior: First argument must be a column vector of double-precision numbers");

  VectorConstView estParams(mxGetPr(prhs[0]), mxGetM(prhs[0]), 1);

  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("logposterior: Second argument must be a matrix of double-precision numbers");

  MatrixConstView data(mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]), mxGetM(prhs[1]));

  if (!mxIsChar(prhs[2]))
    mexErrMsgTxt("logposterior: Third argument must be a character string");

  char *mexext_mx = mxArrayToString(prhs[2]);
  std::string mexext(mexext_mx);
  mxFree(mexext_mx);

  size_t fblock = (size_t) mxGetScalar(prhs[3]);
  size_t nBlocks = (size_t) mxGetScalar(prhs[4]);
  VectorConstView nMHruns(mxGetPr(prhs[5]), mxGetM(prhs[5]), 1);
  assert(nMHruns.getSize() == nBlocks);

  MatrixConstView D(mxGetPr(prhs[6]), mxGetM(prhs[6]), mxGetN(prhs[6]), mxGetM(prhs[6]));

  // Compute MCMC MH Draws and get last line run in the last MH block sub-array 
  size_t lastMHblockArrayLine = logMCMCposterior(estParams, data, mexext, fblock, nBlocks, nMHruns, D);

  *mxGetPr(plhs[0]) = (double) lastMHblockArrayLine;
}
