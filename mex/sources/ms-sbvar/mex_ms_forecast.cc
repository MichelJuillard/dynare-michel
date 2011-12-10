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

#include <iostream>

extern "C"  {
#include "modify_for_mex.h"
#include "switch.h"
#include "switchio.h"
#include "VARio.h"
#include "dw_histogram.h"
#include "sbvar_forecast.h"
}

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  double *out_buf;
  int i, j, k, s, nfree, nstates, nvars, npre, T;

  TStateModel *model;
  SbvarOption *options = (SbvarOption *) NULL;
  T_VAR_Parameters *p = (T_VAR_Parameters *) NULL;
  TMatrixHistogram *histogram = (TMatrixHistogram *) NULL;
  TMatrix forecast;

  int type = F_FREE;

  // Check the left hand and right hand side arguments to make sure they conform
  if (nrhs != 1)
    DYN_MEX_FUNC_ERR_MSG_TXT("ms_forecast takes one cell array as an input argument.");
  if (nlhs != 2)
    DYN_MEX_FUNC_ERR_MSG_TXT("You must specify two output arguments.");

  model = initialize_model_and_options(&options, prhs, &nstates, &nvars, &npre, &nfree);
  if (model == NULL || options == NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("There was a problem initializing the model, can not continue");

  p = (T_VAR_Parameters *) (model->theta);
  T = p->nobs;

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
        // regimes x percentile x horizon x nvar
        mwSize dims[4] = {nstates, options->num_percentiles, options->horizon, nvars};
        plhs[1] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
        out_buf = mxGetPr(plhs[1]);
      }
    else
      {
        // regimes x horizon x (nvar*nvar)
        mwSize dims[3] = {nstates, options->horizon, nvars};
        plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        out_buf = mxGetPr(plhs[1]);
      }
  else
  if (options->num_percentiles > 1)
    {
      // percentile x horizon x (nvar*nvar)
      mwSize dims[3] = {options->num_percentiles, options->horizon, nvars};
      plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
      out_buf = mxGetPr(plhs[1]);
    }
  else
    {
      // horizon x (nvar*nvar)
      mwSize dims[2] = {options->horizon, nvars};
      plhs[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      out_buf = mxGetPr(plhs[1]);
    }

  // hard coding dates for now to be off
  int dates = 0;
  try
    {
      if (options->regimes)
        for (s = 0; s < nstates; s++)
          {
            printf("Constructing percentiles for forecast - regime %d\n", s);
            if (options->simulation_file) rewind(options->simulation_file);
            if (histogram = forecast_percentile_regime(options->shocks, options->simulation_file, options->thin, s, T, options->horizon, model, type))
              {
                for (k = 0; k < options->num_percentiles; k++)
                  if (forecast = ForecastHistogramToPercentile((TMatrix) NULL, options->number_observations, options->horizon, T, options->percentiles[k], histogram, model, dates))
                    {
                      // the case where we are only dealing with one percentile so we can reduce dimensionality of output argument
                      if (options->num_percentiles == 1)
                        for (i = 0; i < options->horizon; i++)
                          for (j = 0; j < nvars; j++)
                            out_buf[(s) + ((i) + (j)*options->horizon)*nstates] = ElementM(forecast, i, j);
                      else
                        for (i = 0; i < options->horizon; i++)
                          for (j = 0; j < nvars; j++)
                            out_buf[(s) + ((k) + ((i) + (j)*options->horizon)*options->num_percentiles)*nstates] = ElementM(forecast, i, j);
                      mxFree(forecast);
                    }
                mxFree(histogram);
              }
          }

      else if (options->regime >= 0)
        {
          s = options->regime;
          printf("Constructing percentiles for forecast - regime %d\n", s);
          if (histogram = forecast_percentile_regime(options->shocks, options->simulation_file, options->thin, s, T, options->horizon, model, type))
            {
              for (k = 0; k < options->num_percentiles; k++)
                if (forecast = ForecastHistogramToPercentile((TMatrix) NULL, options->number_observations, options->horizon, T, options->percentiles[k], histogram, model, dates))
                  {
                    // the case where we are only dealing with one percentile so we can reduce dimensionality of output argument
                    if (options->num_percentiles == 1)
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < nvars; j++)
                          out_buf[(i) + (j)*options->horizon] = ElementM(forecast, i, j);
                    else
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < nvars; j++)
                          out_buf[(k) + ((i) + (j)*options->horizon)*options->num_percentiles] = ElementM(forecast, i, j);
                    mxFree(forecast);
                  }
              mxFree(histogram);
            }
        }
      else
        {
          printf("Constructing percentiles for forecasts - %d draws of shocks/regimes per posterior value\n", options->shocks);
          if (histogram = forecast_percentile(options->shocks, options->simulation_file, options->thin, T, options->horizon, model, type))
            {
              for (k = 0; k < options->num_percentiles; k++)
                if (forecast = ForecastHistogramToPercentile((TMatrix) NULL, options->number_observations, options->horizon, T, options->percentiles[k], histogram, model, dates))
                  {
                    // the case where we are only dealing with one percentile so we can reduce dimensionality of output argument
                    if (options->num_percentiles == 1)
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < nvars; j++)
                          out_buf[(i) + (j)*options->horizon] = ElementM(forecast, i, j);
                    else
                      for (i = 0; i < options->horizon; i++)
                        for (j = 0; j < nvars; j++)
                          out_buf[(k) + ((i) + (j)*options->horizon)*options->num_percentiles] = ElementM(forecast, i, j);
                    mxFree(forecast);
                  }
            }
          mxFree(histogram);
        }
    }
  catch (const char *s)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("Exception Thrown in Forecast: \n");
    }

  mxFree(p);
  mxFree(model);
  plhs[0] = mxCreateDoubleScalar(0);
}

#endif
