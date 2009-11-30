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

/******************************************************
   // k_order_perturbation.cpp : Defines the entry point for the k-order perturbation application DLL.
   //
   // called from Dynare dr1_k_order.m, (itself called form resol.m instead of regular dr1.m)
   //            if options_.order < 2 % 1st order only
   //                [ysteady, ghx_u]=k_order_perturbation(dr,task,M_,options_, oo_ , ['.' mexext]);
   //            else % 2nd order
   //                [ysteady, ghx_u, g_2]=k_order_perturbation(dr,task,M_,options_, oo_ , ['.' mexext]);
   // inputs:
   //			dr,		- Dynare structure
   //			task,  - check or not, not used
   //			M_		- Dynare structure
   //			options_ - Dynare structure
   //			oo_		- Dynare structure
   //			['.' mexext] Matlab dll extension
   // returns:
   //			 ysteady steady state
   //			ghx_u - first order rules packed in one matrix
   //			g_2 - 2nd order rules packed in one matrix
 **********************************************************/

#include "k_ord_dynare.h"
#include "dynamic_dll.h"

#include <cmath>
#include <cstring>
#include <cctype>

#ifdef _MSC_VER

BOOL APIENTRY
DllMain(HANDLE hModule,
        DWORD  ul_reason_for_call,
        LPVOID lpReserved
        )
{
  switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
      break;
    }
  return TRUE;
}

// Some MS Windows preambles
// This is an example of an exported variable
K_ORDER_PERTURBATION_API int nK_order_perturbation = 0;

// This is an example of an exported function.
K_ORDER_PERTURBATION_API int
fnK_order_perturbation(void)
{
  return 42;
}

#endif // _MSC_VER

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)  // exclude mexFunction for other applications

const char **
DynareMxArrayToString(const char *cNamesCharStr, const int len, const int width)
{
  char **cNamesMX;
  cNamesMX = (char **) calloc(len, sizeof(char *));
  for(int i = 0; i < len; i++)
    cNamesMX[i] = (char *) calloc(width+1, sizeof(char));

  for (int i = 0; i < width; i++)
    {
      for (int j = 0; j < len; j++)
        {
          // Allow alphanumeric and underscores "_" only:
          if (isalnum(cNamesCharStr[j+i*len]) || ('_' == cNamesCharStr[j+i*len]))
            {
              cNamesMX[j][i] = cNamesCharStr[j+i*len];
            }
          else cNamesMX[j][i] = '\0';
        }
    }
  const char **ret = (const char **) calloc(len, sizeof(char *));
  for (int j = 0; j < len; j++)
    {
      cNamesMX[j][width] = '\0';
      char *token = (char *) calloc(strlen(cNamesMX[j])+1, sizeof(char));
      strcpy(token, cNamesMX[j]);
      ret[j] = token;
    }
  mxFree(cNamesMX);
  return ret;
}

//////////////////////////////////////////////////////
// Convert Matlab Dynare endo and exo names array to C type array of string pointers
// Poblem is that Matlab mx function returns a long string concatenated by columns rather than rows
// hence a rather low level approach is needed
///////////////////////////////////////////////////////
const char **
DynareMxArrayToString(const mxArray *mxFldp, const int len, const int width)
{
  char *cNamesCharStr = mxArrayToString(mxFldp);
  const char **ret = DynareMxArrayToString(cNamesCharStr, len, width);
  return ret;
}

extern "C" {

  // mexFunction: Matlab Inerface point and the main application driver
  void
  mexFunction(int nlhs, mxArray *plhs[],
              int nrhs, const mxArray *prhs[])
  {
    if (nrhs < 5)
      mexErrMsgTxt("Must have at least 5 input parameters.");
    if (nlhs == 0)
      mexErrMsgTxt("Must have at least 1 output parameter.");

    const mxArray *dr = prhs[0];
    const int check_flag = (int) mxGetScalar(prhs[1]);
    const mxArray *M_ = prhs[2];
    const mxArray *options_ = prhs[3];
    const mxArray *oo_ = prhs[4];

    mxArray *mFname = mxGetField(M_, 0, "fname");
    if (!mxIsChar(mFname))
        mexErrMsgTxt("Input must be of type char.");
    string fName = mxArrayToString(mFname);
    const mxArray *mexExt = prhs[5];
    string dfExt = mxArrayToString(mexExt); //Dynamic file extension, e.g.".dll" or .mexw32;

    int kOrder;
    mxArray *mxFldp = mxGetField(options_, 0, "order");
    if (mxIsNumeric(mxFldp))
      kOrder = (int) mxGetScalar(mxFldp);
    else
      kOrder = 1;
    
    if (kOrder == 1 && nlhs != 1)
      mexErrMsgTxt("k_order_perturbation at order 1 requires exactly 1 argument in output");
    else if (kOrder > 1 && nlhs != kOrder+1)
      mexErrMsgTxt("k_order_perturbation at order > 1 requires exactly order + 1 argument in output");
    
    double qz_criterium = 1+1e-6;
    mxFldp = mxGetField(options_, 0, "qz_criterium");
    if (mxIsNumeric(mxFldp))
      qz_criterium = (double) mxGetScalar(mxFldp);

    mxFldp      = mxGetField(M_, 0, "params");
    double *dparams = (double *) mxGetData(mxFldp);
    int npar = (int) mxGetM(mxFldp);
    Vector *modParams =  new Vector(dparams, npar);

    mxFldp      = mxGetField(M_, 0, "Sigma_e");
    dparams = (double *) mxGetData(mxFldp);
    npar = (int) mxGetN(mxFldp);
    TwoDMatrix *vCov =  new TwoDMatrix(npar, npar, dparams);

    mxFldp      = mxGetField(dr, 0, "ys");  // and not in order of dr.order_var
    dparams = (double *) mxGetData(mxFldp);
    const int nSteady = (int) mxGetM(mxFldp);
    Vector *ySteady =  new Vector(dparams, nSteady);

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
    // it_ should be set to M_.maximum_lag
    mxFldp = mxGetField(M_, 0, "maximum_lag");
    const int nMax_lag = (int) mxGetScalar(mxFldp);

    nPred -= nBoth; // correct nPred for nBoth.

    mxFldp      = mxGetField(dr, 0, "order_var");
    dparams = (double *) mxGetData(mxFldp);
    npar = (int) mxGetM(mxFldp);
    if (npar != nEndo)    //(nPar != npar)
        mexErrMsgTxt("Incorrect number of input var_order vars.");
    vector<int> *var_order_vp = (new vector<int>(nEndo));
    for (int v = 0; v < nEndo; v++)
        (*var_order_vp)[v] = (int)(*(dparams++));

    // the lag, current and lead blocks of the jacobian respectively
    mxFldp      = mxGetField(M_, 0, "lead_lag_incidence");
    dparams = (double *) mxGetData(mxFldp);
    npar = (int) mxGetN(mxFldp);
    int nrows = (int) mxGetM(mxFldp);


    TwoDMatrix *llincidence =  new TwoDMatrix(nrows, npar, dparams);
    if (npar != nEndo)
      mexErrMsgIdAndTxt("dynare:k_order_perturbation", "Incorrect length of lead lag incidences: ncol=%d != nEndo=%d.", npar, nEndo);

    //get NNZH =NNZD(2) = the total number of non-zero Hessian elements 
    mxFldp = mxGetField(M_, 0, "NNZDerivatives");
    dparams = (double *) mxGetData(mxFldp);
    Vector *NNZD =  new Vector (dparams, (int) mxGetM(mxFldp));

    const int jcols = nExog+nEndo+nsPred+nsForw; // Num of Jacobian columns

    mxFldp = mxGetField(M_, 0, "var_order_endo_names");
    const int nendo = (int) mxGetM(mxFldp);
    const int widthEndo = (int) mxGetN(mxFldp);
    const char **endoNamesMX = DynareMxArrayToString(mxFldp, nendo, widthEndo);

    mxFldp      = mxGetField(M_, 0, "exo_names");
    const int nexo = (int) mxGetM(mxFldp);
    const int widthExog = (int) mxGetN(mxFldp);
    const char **exoNamesMX = DynareMxArrayToString(mxFldp, nexo, widthExog);

    if ((nEndo != nendo) || (nExog != nexo))    //(nPar != npar)
      mexErrMsgTxt("Incorrect number of input parameters.");

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
        DynamicModelDLL dynamicDLL(fName, nEndo, jcols, nMax_lag, nExog, dfExt);

        // intiate tensor library
        tls.init(kOrder, nStat+2*nPred+3*nBoth+2*nForw+nExog);

        // make KordpDynare object
        KordpDynare dynare(endoNamesMX,  nEndo, exoNamesMX,  nExog, nPar, // paramNames,
                           ySteady, vCov, modParams, nStat, nPred, nForw, nBoth,
                           jcols, NNZD, nSteps, kOrder, journal, dynamicDLL, 
                           sstol, var_order_vp, llincidence, qz_criterium);

        // construct main K-order approximation class

        Approximation app(dynare, journal,  nSteps, false, qz_criterium);
        // run stochastic steady
        app.walkStochSteady();

        ConstTwoDMatrix ss(app.getSS());
        // open mat file
        std::string matfile(fName); //(params.basename);
        matfile += ".mat";
        FILE *matfd = NULL;
        if (NULL == (matfd = fopen(matfile.c_str(), "wb")))
          mexErrMsgIdAndTxt("dynare:k_order_perturbation", "Couldn't open %s for writing.", matfile.c_str());

        std::string ss_matrix_name(fName);
        ss_matrix_name += "_steady_states";
        ss.writeMat4(matfd, ss_matrix_name.c_str());

        // write the folded decision rule to the Mat-4 file
        app.getFoldDecisionRule().writeMat4(matfd, fName.c_str()); //params.prefix);

        fclose(matfd);

        /* Write derivative outputs into memory map */
        map<string, ConstTwoDMatrix> mm;
        app.getFoldDecisionRule().writeMMap(mm, string());

        // get latest ysteady
        double *dYsteady = (dynare.getSteady().base());
        ySteady = (Vector *)(&dynare.getSteady());

        // developement of the output.
        double  *dgy, *dgu, *ysteady;
        int nb_row_x;

        ysteady = NULL;
        if (kOrder == 1)
          {
            /* Set the output pointer to the output matrix ysteady. */
	    map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
	    ++cit;
            plhs[0] = mxCreateDoubleMatrix((*cit).second.numRows(), (*cit).second.numCols(), mxREAL);
	    TwoDMatrix g((*cit).second.numRows(), (*cit).second.numCols(), mxGetPr(plhs[0]));
	    g = (const TwoDMatrix &)(*cit).second;
          }
        if (kOrder >= 2)
          {
            int ii = 0;
            for (map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
                 ((cit != mm.end()) && (ii < nlhs)); ++cit)
              {
                {
                  plhs[ii] = mxCreateDoubleMatrix((*cit).second.numRows(), (*cit).second.numCols(), mxREAL);
                  TwoDMatrix g_ii((*cit).second.numRows(), (*cit).second.numCols(), mxGetPr(plhs[ii]));
                  g_ii = (const TwoDMatrix &)(*cit).second;
                  ++ii;
                }
              }
          }
      }
    catch (const KordException &e)
      {
        e.print();
        mexErrMsgIdAndTxt("dynare:k_order_perturbation", "Caught Kord exception: %s", e.get_message());
      }
    catch (const TLException &e)
      {
        e.print();
        mexErrMsgIdAndTxt("dynare:k_order_perturbation", "Caught TL exception");
      }
    catch (SylvException &e)
      {
        e.printMessage();
        mexErrMsgIdAndTxt("dynare:k_order_perturbation", "Caught Sylv exception");
      }
    catch (const DynareException &e)
      {
        mexErrMsgIdAndTxt("dynare:k_order_perturbation", "Caught KordDynare exception: %s", e.message());
      }
    catch (const ogu::Exception &e)
      {
        mexErrMsgIdAndTxt("dynare:k_order_perturbation", "Caught general exception: %s", e.message());
      }
  } // end of mexFunction()
} // end of extern C

#endif // ifdef MATLAB_MEX_FILE  to exclude mexFunction for other applications
