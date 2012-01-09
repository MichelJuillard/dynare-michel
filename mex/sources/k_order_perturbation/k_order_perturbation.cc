/*
 * Copyright (C) 2008-2012 Dynare Team
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

/*
  Defines the entry point for the k-order perturbation application DLL.

  Inputs:
  1) dr
  2) M_
  3) options

  Outputs:
  - if order == 1: only g_1
  - if order == 2: g_0, g_1, g_2
  - if order == 3: g_0, g_1, g_2, g_3
*/

#include "dynamic_m.hh"
#include "dynamic_dll.hh"

#include <cmath>
#include <cstring>
#include <cctype>
#include <cassert>

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)  // exclude mexFunction for other applications

# include "dynmex.h"

//////////////////////////////////////////////////////
// Convert MATLAB Dynare endo and exo names array to a vector<string> array of string pointers
// Poblem is that Matlab mx function returns a long string concatenated by columns rather than rows
// hence a rather low level approach is needed
///////////////////////////////////////////////////////
void
DynareMxArrayToString(const mxArray *mxFldp, const int len, const int width, vector<string> &out)
{
  char *cNamesCharStr = mxArrayToString(mxFldp);

  out.resize(len);

  for (int i = 0; i < width; i++)
    for (int j = 0; j < len; j++)
      // Allow alphanumeric and underscores "_" only:
      if (isalnum(cNamesCharStr[j+i*len]) || (cNamesCharStr[j+i*len] == '_'))
        out[j] += cNamesCharStr[j+i*len];
}

extern "C" {

  void
  mexFunction(int nlhs, mxArray *plhs[],
              int nrhs, const mxArray *prhs[])
  {
    if (nrhs != 3 || nlhs < 2)
      DYN_MEX_FUNC_ERR_MSG_TXT("Must have exactly 3 input parameters and take at least 2 output parameters.");

    const mxArray *dr = prhs[0];
    const mxArray *M_ = prhs[1];
    const mxArray *options_ = prhs[2];
    int use_dll = (int) mxGetScalar(mxGetField(options_, 0, "use_dll"));

    mxArray *mFname = mxGetField(M_, 0, "fname");
    if (!mxIsChar(mFname))
      DYN_MEX_FUNC_ERR_MSG_TXT("Input must be of type char.");

    string fName = mxArrayToString(mFname);

    int kOrder;
    mxArray *mxFldp = mxGetField(options_, 0, "order");
    if (mxIsNumeric(mxFldp))
      kOrder = (int) mxGetScalar(mxFldp);
    else
      kOrder = 1;

    if (kOrder == 1 && nlhs != 2)
      DYN_MEX_FUNC_ERR_MSG_TXT("k_order_perturbation at order 1 requires exactly 2 arguments in output");
    else if (kOrder > 1 && nlhs != kOrder+2)
      DYN_MEX_FUNC_ERR_MSG_TXT("k_order_perturbation at order > 1 requires exactly order+2 arguments in output");

    double qz_criterium = 1+1e-6;
    mxFldp = mxGetField(options_, 0, "qz_criterium");
    if (mxIsNumeric(mxFldp))
      qz_criterium = (double) mxGetScalar(mxFldp);

    mxFldp = mxGetField(M_, 0, "params");
    double *dparams = mxGetPr(mxFldp);
    int npar = (int) mxGetM(mxFldp);
    Vector modParams(dparams, npar);
    if (!modParams.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("The parameters vector contains NaN or Inf");

    mxFldp = mxGetField(M_, 0, "Sigma_e");
    dparams = mxGetPr(mxFldp);
    npar = (int) mxGetN(mxFldp);
    TwoDMatrix vCov(npar, npar, dparams);
    if (!vCov.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("The covariance matrix of shocks contains NaN or Inf");

    mxFldp = mxGetField(dr, 0, "ys");  // and not in order of dr.order_var
    dparams = mxGetPr(mxFldp);
    const int nSteady = (int) mxGetM(mxFldp);
    Vector ySteady(dparams, nSteady);
    if (!ySteady.isFinite())
      DYN_MEX_FUNC_ERR_MSG_TXT("The steady state vector contains NaN or Inf");

    mxFldp = mxGetField(dr, 0, "nstatic");
    const int nStat = (int) mxGetScalar(mxFldp);
    mxFldp = mxGetField(dr, 0, "npred");
    int nPred = (int) mxGetScalar(mxFldp);
    mxFldp = mxGetField(dr, 0, "nspred");
    const int nsPred = (int) mxGetScalar(mxFldp);
    mxFldp = mxGetField(dr, 0, "nboth");
    const int nBoth = (int) mxGetScalar(mxFldp);
    mxFldp = mxGetField(dr, 0, "nfwrd");
    const int nForw = (int) mxGetScalar(mxFldp);
    mxFldp = mxGetField(dr, 0, "nsfwrd");
    const int nsForw = (int) mxGetScalar(mxFldp);

    mxFldp = mxGetField(M_, 0, "exo_nbr");
    const int nExog = (int) mxGetScalar(mxFldp);
    mxFldp = mxGetField(M_, 0, "endo_nbr");
    const int nEndo = (int) mxGetScalar(mxFldp);
    mxFldp = mxGetField(M_, 0, "param_nbr");
    const int nPar = (int) mxGetScalar(mxFldp);

    nPred -= nBoth; // correct nPred for nBoth.

    mxFldp = mxGetField(dr, 0, "order_var");
    dparams = mxGetPr(mxFldp);
    npar = (int) mxGetM(mxFldp);
    if (npar != nEndo)
      DYN_MEX_FUNC_ERR_MSG_TXT("Incorrect number of input var_order vars.");

    vector<int> var_order_vp(nEndo);
    for (int v = 0; v < nEndo; v++)
      var_order_vp[v] = (int) (*(dparams++));

    // the lag, current and lead blocks of the jacobian respectively
    mxFldp = mxGetField(M_, 0, "lead_lag_incidence");
    dparams = mxGetPr(mxFldp);
    npar = (int) mxGetN(mxFldp);
    int nrows = (int) mxGetM(mxFldp);

    TwoDMatrix llincidence(nrows, npar, dparams);
    if (npar != nEndo)
      {
        ostringstream strstrm;
        strstrm << "dynare:k_order_perturbation " << "Incorrect length of lead lag incidences: ncol=" << npar << " != nEndo=" << nEndo;
        DYN_MEX_FUNC_ERR_MSG_TXT(strstrm.str().c_str());
      }
    //get NNZH =NNZD(2) = the total number of non-zero Hessian elements
    mxFldp = mxGetField(M_, 0, "NNZDerivatives");
    dparams = mxGetPr(mxFldp);
    Vector NNZD(dparams, (int) mxGetM(mxFldp));
    if (NNZD[kOrder-1] == -1)
      DYN_MEX_FUNC_ERR_MSG_TXT("The derivatives were not computed for the required order. Make sure that you used the right order option inside the stoch_simul command");

    const int jcols = nExog+nEndo+nsPred+nsForw; // Num of Jacobian columns

    mxFldp = mxGetField(M_, 0, "var_order_endo_names");
    const int nendo = (int) mxGetM(mxFldp);
    const int widthEndo = (int) mxGetN(mxFldp);
    vector<string> endoNames;
    DynareMxArrayToString(mxFldp, nendo, widthEndo, endoNames);

    mxFldp = mxGetField(M_, 0, "exo_names");
    const int nexo = (int) mxGetM(mxFldp);
    const int widthExog = (int) mxGetN(mxFldp);
    vector<string> exoNames;
    DynareMxArrayToString(mxFldp, nexo, widthExog, exoNames);

    if ((nEndo != nendo) || (nExog != nexo))
      DYN_MEX_FUNC_ERR_MSG_TXT("Incorrect number of input parameters.");

    /* Fetch time index */

    const int nSteps = 0; // Dynare++ solving steps, for time being default to 0 = deterministic steady state
    const double sstol = 1.e-13; //NL solver tolerance from

    THREAD_GROUP::max_parallel_threads = 2; //params.num_threads;

    try
      {
        // make journal name and journal
        std::string jName(fName); //params.basename);
        jName += ".jnl";
        Journal journal(jName.c_str());

        DynamicModelAC *dynamicModelFile;
        if (use_dll == 1)
          dynamicModelFile = new DynamicModelDLL(fName);
        else
          dynamicModelFile = new DynamicModelMFile(fName);

        // intiate tensor library
        tls.init(kOrder, nStat+2*nPred+3*nBoth+2*nForw+nExog);

        // make KordpDynare object
        KordpDynare dynare(endoNames, nEndo, exoNames, nExog, nPar,
                           ySteady, vCov, modParams, nStat, nPred, nForw, nBoth,
                           jcols, NNZD, nSteps, kOrder, journal, dynamicModelFile,
                           sstol, var_order_vp, llincidence, qz_criterium);

        // construct main K-order approximation class

        Approximation app(dynare, journal,  nSteps, false, qz_criterium);
        // run stochastic steady
        app.walkStochSteady();

        /* Write derivative outputs into memory map */
        map<string, ConstTwoDMatrix> mm;
        app.getFoldDecisionRule().writeMMap(mm, string());

        // get latest ysteady
        ySteady = dynare.getSteady();

        if (kOrder == 1)
          {
            /* Set the output pointer to the output matrix ysteady. */
            map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
            ++cit;
            plhs[1] = mxCreateDoubleMatrix((*cit).second.numRows(), (*cit).second.numCols(), mxREAL);

            // Copy Dynare++ matrix into MATLAB matrix
            const ConstVector &vec = (*cit).second.getData();
            assert(vec.skip() == 1);
            memcpy(mxGetPr(plhs[1]), vec.base(), vec.length() * sizeof(double));
          }
        if (kOrder >= 2)
          {
            int ii = 1;
            for (map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
                 ((cit != mm.end()) && (ii < nlhs)); ++cit)
              {
                {
                  plhs[ii] = mxCreateDoubleMatrix((*cit).second.numRows(), (*cit).second.numCols(), mxREAL);

                  // Copy Dynare++ matrix into MATLAB matrix
                  const ConstVector &vec = (*cit).second.getData();
                  assert(vec.skip() == 1);
                  memcpy(mxGetPr(plhs[ii]), vec.base(), vec.length() * sizeof(double));

                  ++ii;
                }
              }
          }
      }
    catch (const KordException &e)
      {
        e.print();
        ostringstream strstrm;
        strstrm << "dynare:k_order_perturbation: Caught Kord exception: " << e.get_message();
        DYN_MEX_FUNC_ERR_MSG_TXT(strstrm.str().c_str());
      }
    catch (const TLException &e)
      {
        e.print();
        DYN_MEX_FUNC_ERR_MSG_TXT("dynare:k_order_perturbation: Caught TL exception");
      }
    catch (SylvException &e)
      {
        e.printMessage();
        DYN_MEX_FUNC_ERR_MSG_TXT("dynare:k_order_perturbation: Caught Sylv exception");
      }
    catch (const DynareException &e)
      {
        ostringstream strstrm;
        strstrm << "dynare:k_order_perturbation: Caught KordDynare exception: " << e.message();
        DYN_MEX_FUNC_ERR_MSG_TXT(strstrm.str().c_str());
      }
    catch (const ogu::Exception &e)
      {
        ostringstream strstrm;
        strstrm << "dynare:k_order_perturbation: Caught general exception: " << e.message();
        DYN_MEX_FUNC_ERR_MSG_TXT(strstrm.str().c_str());
      }
    plhs[0] = mxCreateDoubleScalar(0);
  } // end of mexFunction()
} // end of extern C

#endif // ifdef MATLAB_MEX_FILE  to exclude mexFunction for other applications
