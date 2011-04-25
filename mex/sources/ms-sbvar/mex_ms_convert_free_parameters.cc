/*
 * Copyright (C) 2011 Dynare Team
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

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mex_ms_sbvar.h"

extern "C"  {
#include "modify_for_mex.h"
#include "switch.h"
#include "switchio.h"
#include "VARio.h"
}

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  using namespace std;

  char *input_buf;
  double *free_parameters;
  int nvars, npre, nstates, nfree;
  double *aplus = NULL, *a0 = NULL, *zeta = NULL, *q = NULL;
  mxArray *Q = NULL;
  TStateModel *model;
  SbvarOption *options = NULL;

  /* input must be a string */
  if (mxIsChar(prhs[0]) != 1)
    DYN_MEX_FUNC_ERR_MSG_TXT("First argument has to be a string to the init_filename.");
  if (!mxIsDouble(prhs[1]))
    DYN_MEX_FUNC_ERR_MSG_TXT("Second argument is a vector of free parameters");

  if (nlhs < 3)
    DYN_MEX_FUNC_ERR_MSG_TXT("You must specify at least three output arguments [A0,Aplus,Zeta]");

  // copy the string data from prhs[0] into a C string input_ buf.    */
  input_buf = mxArrayToString(prhs[0]);
  if (input_buf == NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("Could not convert input to string.");

  // second element should be vector of free parameters */
  free_parameters = mxGetPr(prhs[1]);

  // copy the string data from prhs[0] into a C string input_ buf.    */
  input_buf = mxArrayToString(prhs[0]);
  if (input_buf == NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("Could not convert input to string.");
  model =  initialize_model_and_options(input_buf, &options, nrhs, prhs, &nstates, &nvars, &npre, &nfree);
  if (model == NULL || options == NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("There was a problem initializing the model, can not continue");

  // Check the size of the parameters that were passed in
  size_t npars_passed =  mxGetN(prhs[1]) > mxGetM(prhs[1]) ? mxGetN(prhs[1])  : mxGetM(prhs[1]);
  // If it is 2 longer than the number of free parameters assume that the first two numbers should be ignored (log posterior, log likelihood)
  if (nfree + 2 == (int) npars_passed)
    free_parameters = free_parameters + 2;
  else if (nfree != (int) npars_passed && nfree + 2 != (int) npars_passed)
    {
      printf("\n\nThe model requires a free paramter vector of length %d, you passed one of length %d\n\n",
             (int) nfree, (int) npars_passed);
      DYN_MEX_FUNC_ERR_MSG_TXT("The Free Parameter Array is the wrong size for this model\n");
    }

  // Ao (nstates x nvars x nvars)
  mwSize dims[3] = {nstates, nvars, nvars};
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  a0 = mxGetPr(plhs[0]);

  // Zeta (nstates x nvars x nvars)
  plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  zeta = mxGetPr(plhs[2]);

  // Aplus (nstates x (nlags*nvars +nconstant) x nvars)
  dims[1] = npre;
  plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  aplus = mxGetPr(plhs[1]);

  // Grand Transition Matrix
  Q = mxCreateDoubleMatrix(nstates, nstates, mxREAL);
  q = mxGetPr(Q);

  // get the matrices
  int ret = convert_free_parameters_to_VAR(model, free_parameters, a0, aplus, zeta, q);
  if (ret > 0)
    {
      char *error_msg;
      sprintf(error_msg = (char *) mxMalloc(255*sizeof(char)), "There was a problem converting the free parameters for the given model and inputs\nError Code: %d", ret);
      DYN_MEX_FUNC_ERR_MSG_TXT(error_msg);
    }

  // if they have passed 4 arguments then pass back Q
  if (nlhs > 3)
    plhs[3] = Q;

  mxFree(model);

}

#endif
