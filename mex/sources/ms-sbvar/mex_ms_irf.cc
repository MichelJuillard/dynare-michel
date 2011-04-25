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
#include "dynmex.h"

#include "mex_ms_sbvar.h"

extern "C"  {
#include "modify_for_mex.h"
#include "switch.h"
#include "switchio.h"
#include "VARio.h"
#include "dw_rand.h"
#include "dw_histogram.h"
#include "sbvar_impulse_responses.h"
}

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{

  char *input_buf;
  double *out_buf;
  int i, j, k, s, nfree, nstates, nvars, npre;

  TStateModel *model;
  SbvarOption *options = (SbvarOption *) NULL;
  T_VAR_Parameters *p = (T_VAR_Parameters *) NULL;
  TMatrixHistogram *histogram = (TMatrixHistogram *) NULL;
  TMatrix ir;

  int type = F_FREE, ergodic = 1;

  // Check the left hand and right hand side arguments to make sure they conform
  if (mxIsChar(prhs[0]) != 1)
    DYN_MEX_FUNC_ERR_MSG_TXT("First argument has to be a string to the init_filename.");
  if (nlhs != 1)
    DYN_MEX_FUNC_ERR_MSG_TXT("You must specify one output argument");

  // copy the string data from prhs[0] into a C string input_ buf.    */
  input_buf = mxArrayToString(prhs[0]);
  if (input_buf == NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("Could not convert input to string.");
  model =  initialize_model_and_options(input_buf, &options, nrhs, prhs, &nstates, &nvars, &npre, &nfree);
  if (model == NULL || options == NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("There was a problem initializing the model, can not continue");

  p = (T_VAR_Parameters *) (model->theta);

  // Check to make sure that there is a simulation file present if we are
  // using parameter uncertainty
  if (options->parameter_uncertainty && options->simulation_file == NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("Paramter Uncertainty Was Specified but the simulation file was not found, please specify the simulation file to use with: 'simulation_file',<filename>");
  if (!options->parameter_uncertainty)
    options->simulation_file = (FILE *) NULL;

  // Allocate the output matrix
  if (options->regimes)
    if (options->num_percentiles > 1)
      {
        // regimes x percentile x horizon x (nvar*nvar)
        mwSize dims[4] = {nstates, options->num_percentiles, options->horizon, (nvars*nvars)};
        plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
        out_buf = mxGetPr(plhs[0]);
      }
    else
      {
        // regimes x horizon x (nvar*nvar)
        mwSize dims[3] = {nstates, options->horizon, (nvars*nvars)};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        out_buf = mxGetPr(plhs[0]);
      }
  else
  if (options->num_percentiles > 1)
    {
      // percentile x horizon x (nvar*nvar)
      mwSize dims[3] = {options->num_percentiles, options->horizon, (nvars*nvars)};
      plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      out_buf = mxGetPr(plhs[0]);
    }
  else
    {
      // horizon x (nvar*nvar)
      mwSize dims[2] = {options->horizon, (nvars*nvars)};
      plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      out_buf = mxGetPr(plhs[0]);
    }

  // Use filter probabilities?
  if (options->filtered_probabilities)
    ergodic = 0;

  try
    {
      if (options->regimes)
        for (s = 0; s < nstates; s++)
          {
            printf("Constructing percentiles for impulse responses - regime %d\n", s);
            if (histogram = impulse_response_percentile_regime(options->simulation_file, options->thin, s, options->horizon, model, type))
              {
                for (k = 0; k < options->num_percentiles; k++)
                  if (ir = IRHistogramToPercentile((TMatrix) NULL, options->horizon, options->percentiles[k], histogram, model))
                    {
                      // the case where we are only dealing with one percentile so we can reduce dimensionality of output argument
                      if (options->num_percentiles == 1)
                        for (i = 0; i < options->horizon; i++)
                          for (j = 0; j < (nvars*nvars); j++)
                            out_buf[(s) + ((i) + (j)*options->horizon)*nstates] = ElementM(ir, i, j);
                      else
                        for (i = 0; i < options->horizon; i++)
                          for (j = 0; j < (nvars*nvars); j++)
                            out_buf[(s) + ((k) + ((i) + (j)*options->horizon)*options->num_percentiles)*nstates] = ElementM(ir, i, j);
                      mxFree(ir);
                    }
              }
            mxFree(histogram);
          }

      else if (options->regime >= 0)
        {
          s = options->regime;
          printf("Constructing percentiles for impulse responses - regime %d\n", s);
          if (histogram = impulse_response_percentile_regime(options->simulation_file, options->thin, options->regime, options->horizon, model, type))
            {
              for (k = 0; k < options->num_percentiles; k++)
                if (ir = IRHistogramToPercentile((TMatrix) NULL, options->horizon, options->percentiles[k], histogram, model))
                  {
                    // the case where we are only dealing with one percentile so we can reduce dimensionality of output argument
                    if (options->num_percentiles == 1)
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < (nvars*nvars); j++)
                          out_buf[(i) + (j)*options->horizon] = ElementM(ir, i, j);
                    else
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < (nvars*nvars); j++)
                          out_buf[(k) + ((i) + (j)*options->horizon)*options->num_percentiles] = ElementM(ir, i, j);
                    mxFree(ir);
                  }
            }
          mxFree(histogram);
        }
      else
        {
          printf(ergodic ? "Constructing percentiles for ergodic impulse responses - %d draws of regimes per posterior value\n"
                 : "Constructing percentiles for filtered impulse responses - %d draws of regimes per posterior value\n", options->shocks);
          if (histogram = impulse_response_percentile_ergodic(options->shocks, options->simulation_file, options->thin, ergodic, options->horizon, model, type))
            {
              for (k = 0; k < options->num_percentiles; k++)
                if (ir = IRHistogramToPercentile((TMatrix) NULL, options->horizon, options->percentiles[k], histogram, model))
                  {
                    // the case where we are only dealing with one percentile so we can reduce dimensionality of output argument
                    if (options->num_percentiles == 1)
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < (nvars*nvars); j++)
                          out_buf[(i) + (j)*options->horizon] = ElementM(ir, i, j);
                    else
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < (nvars*nvars); j++)
                          out_buf[(k) + ((i) + (j)*options->horizon)*options->num_percentiles] = ElementM(ir, i, j);
                    mxFree(ir);
                  }
            }
          mxFree(histogram);
        }
    }
  catch (const char *s)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("Exception Thrown in IRF: \n");
    }

  mxFree(model);
  mxFree(p);

}

#endif
