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
#include <fstream>

extern "C"  {
#include "modify_for_mex.h"
#include "switch.h"
#include "switchio.h"
#include "VARio.h"
#include "dw_rand.h"
}

using namespace std;

int
file_exist(char *filename)
{
  int ret = 0;
  char ch;

  ifstream fin;
  fin.open(filename, ios_base::in);
  if (fin.is_open() && fin.get(ch))
    ret = 1;

  fin.close();
  return ret;
}

char *
CreateFilenameFromTag(const char *fmt, const char *tag, const char *dir)
{
  char *filename;
  if (!tag)
    tag = "";
  if (!dir)
    dir = "";
  sprintf(filename = (char *) mxMalloc(strlen(dir) + strlen(fmt) + strlen(tag) - 3), fmt, dir, tag);
  return filename;
}

/*
 *  initialize_ms_model
 *  - initiliazes a TstateModel given the filename or tagname of initialization file
 */
TStateModel *
initialize_ms_model(char *name)
{
  TStateModel *model = (TStateModel *) NULL;

  char *filename = (char *) NULL, *est_filename = (char *) NULL, *init_filename = (char *) NULL;
  est_filename = CreateFilenameFromTag("%sest_final_%s.out", name, "");
  init_filename = CreateFilenameFromTag("%sinit_%s.dat", name, "");

  if (file_exist(name))
    filename = name;
  else if (file_exist(est_filename))
    filename = est_filename;
  else if (file_exist(init_filename))
    filename = init_filename;
  else
    {
      printf("Can not open initialization or estimation file!\n");
      return (TStateModel *) NULL;
    }

  try
    {
      if (!(model = Read_VAR_Specification((FILE *) NULL, filename)))
        {
          mxFree(filename);
          printf("Can not initialize model with given name!\n");
          return (TStateModel *) NULL;
        }

      if (est_filename != NULL)
        {
          if (!strstr(est_filename, filename))
            {
              printf("Initializing with empty parameters %s\n", filename);
              // If no estimation file is around then just initialize the sizes

              double *theta = (double *) mxMalloc(sizeof(double)*NumberFreeParametersTheta(model));
              double *q = (double *) mxMalloc(sizeof(double)*NumberFreeParametersQ(model));
              ConvertFreeParametersToQ(model, q);
              ConvertFreeParametersToTheta(model, theta);
            }
        }
    }
  catch (const char *s)
    {
      printf("Error Initializing TStateModel: %s\n", s);
      model = (TStateModel *) NULL;
    }

  mxFree(filename);
  return model;
}

int
get_var_dimensions(TStateModel *model, int *nstates, int *nvars, int *npre, int *nfree)
{

  if (model)
    {
      try
        {
          *nfree = NumberFreeParametersTheta(model)+NumberFreeParametersQ(model);

          T_VAR_Parameters *p = (T_VAR_Parameters *) (model->theta);

          *nstates = model->sv->nstates;
          *nvars =  p->nvars;
          *npre = p->npre;
          return 0;
        }
      catch (const char *s)
        {
          printf("Error Getting VAR Dimensions");
          return 1;
        }

    }
  else
    {
      printf("The model passed was null!\n");
      return 1;
    }
}

int
set_parameters_in_VAR(TStateModel *model, double *free_parameters)
{
  int nfree;
  try
    {
      nfree = NumberFreeParametersTheta(model)+NumberFreeParametersQ(model);
      ConvertFreeParametersToQ(model, free_parameters+NumberFreeParametersTheta(model));
      ConvertFreeParametersToTheta(model, free_parameters);
      ComputeTransitionMatrix(0, model);
    }
  catch (const char *s)
    {
      return 1;
    }
  return 0;
}

int
convert_free_parameters_to_VAR(TStateModel *model, double *free_parameters, double *A0, double *Aplus, double *Zeta, double *Q)
{
  int nstates = 0, nvars = 0, npre = 0;
  int i = 0, j = 0, s = 0;

  if (model == NULL)
    {
      printf("The model passed was null in Convert Free Parameters!\n");
      return 1;
    }

  if (set_parameters_in_VAR(model, free_parameters) > 0)
    {
      printf("Could not set Parameters in model when converting parameters to var\n");
      return 2;
    }

  try
    {
      T_VAR_Parameters *p = (T_VAR_Parameters *) (model->theta);
      nstates = model->sv->nstates;
      nvars = p->nvars;
      npre = p->npre;

      TMatrix a = NULL, aplus = NULL;
      a = CreateMatrix(p->nvars, p->nvars);
      aplus = CreateMatrix(p->npre, p->nvars);

      // Is this is a good idea in matlab?
      if (A0 == NULL)
        A0 = new double [nstates*nvars*nvars];
      //A0 = malloc(sizeof(double)*(nstates*nvars*nvars));
      if (Aplus == NULL)
        Aplus = new double [nstates*npre*nvars];
      //Aplus = malloc(sizeof(double)*(nstates*npre*nvars));
      if (Zeta == NULL)
        Zeta = new double [nstates*nvars*nvars];
      //Zeta = malloc(sizeof(double)*(nstates*nvars*nvars));
      if (Q == NULL)
        Q = new double [nstates*nstates];
      //Q = malloc(sizeof(double)*(nstates*nstates));

      for (s = 0; s < nstates; s++)
        {
          MakeA0(a, s, p);
          for (i = 0; i < nvars; i++)
            {
              for (j = 0; j < nvars; j++)
                A0[(s)+((i+j*nvars)*nstates)] = ElementM(a, i, j);
            }
          MakeAplus(aplus, s, p);
          for (i = 0; i < npre; i++)
            {
              for (j = 0; j < nvars; j++)
                Aplus[(s)+((i+j*npre)*nstates)] = ElementM(aplus, i, j);
            }
          MakeZeta(a, s, p);
          for (i = 0; i < nvars; i++)
            {
              for (j = 0; j < nvars; j++)
                Zeta[(s)+((i+j*nvars)*nstates)] = ElementM(a, i, j);
            }
        }

      // Grand Transition Matrix
      for (i = 0; i < nstates; i++)
        for (j = 0; j < nstates; j++)
          Q[i + j*nstates] = ElementM(model->sv->Q, i, j);

      FreeMatrix(a);
      FreeMatrix(aplus);
    }
  catch (const char *s)
    {
      printf("Exception: convert_free_parameters\n");
      return 3;
    }
  return 0;
}

/*
   Additional command line parameters

   'horizon', <integer>
      If this argument exists, then the forecast horizon is given by the passed
      integer.  The default value is 12.

   'filtered'
      Uses filtered probabilities at the end of the sample as initial conditions
      for regime probabilities.  The default behavior is to us the erogdic
      distribution for the initial conditions.  This flag only applies if neither
      -regimes nor -regime is specified.

   'error_bands'
      Output error bands.  (default = off - only median is computed)

   'percentiles' n p_1 p_2 ... p_n
      Percentiles to compute. The first parameter after percentiles must be the
      number of percentiles and the following values are the actual percentiles.
      default = 3  0.16  0.50  0.84   if error_bands flag is set
              = 1  0.50               otherwise

   'parameter_uncertainty'
      Apply parameter uncertainty when computing error bands or median. When set,
      will default shocks = 1

   'shocks_per_parameter', <integer>
      Number of regime paths to draw for each parameter draw.  The default value
      is 1 if parameter_uncertainty is set and 10,000 otherwise.

   'thin'
      Thinning factor.  Only 1/thin of the draws in posterior draws file are
      used. The default value is 1.

   'regimes'
      Produces forecasts as if each regime were permanent. (default = off)

   'regime' <integer>
      Produces forecasts as if regime were permanent.  Regime numbers are zero
      based.  (default = off)

   'simulation_file', <string>
      name of the file containing the model's simulated free values

   'number_observations', <integer>
         If this argument exists, then the number of data points included in the
         output is given by the passed integer int the forecast output.  The default value is 0.

   'free_parameters' <vector>:
        Vector of free paramters to initialize the model with.

   'median':
      Shortcut for setting 'percentiles',[0.5]

   'mean':
      Compute the mean, instead of percentiles

   'seed':
      Set the Seed for random number generation, default=0 (random)
 */

SbvarOption *
initialize_sbvar_options(char *file_tag)
{
  SbvarOption *options = new SbvarOption;
  options->shocks = 10000;
  options->thin = 1;
  options->horizon = 20;
  options->number_observations = 0;
  options->regime = -1;
  options->regimes = false;
  options->parameter_uncertainty = false;
  options->num_percentiles = 3;
  options->percentiles = new double[3];
  options->percentiles[0] = 0.16;
  options->percentiles[1] = 0.5;
  options->percentiles[2] = 0.84;
  options->filtered_probabilities = false;
  options->num_parameters = -1;
  options->free_parameters = (double *) NULL;
  options->mean = false;
  options->seed = 0;

  if (file_tag != NULL)
    {
      options->simulation_filename = (char *) CreateFilenameFromTag("%ssimulation_%s.out", file_tag, "");
      if (file_exist(options->simulation_filename))
        options->simulation_file = fopen(options->simulation_filename, "r");
    }
  else
    {
      options->simulation_filename = (char *) NULL;
      options->simulation_file = (FILE *) NULL;
    }
  return options;
}

int
set_options(SbvarOption *options, const mxArray *prhs[])
{
  if (options == NULL)
    options = initialize_sbvar_options((char *) NULL);

  double *temp_buf;
  bool shocks_passed = false;
  int num_options = mxGetN(prhs[0]);
  for (int i = 1; i < num_options; i++)
    {
      mxArray *this_option = mxGetCell(prhs[0],i);
      char *option_name_c = mxArrayToString(mxGetCell(this_option,0));
      string option_name (option_name_c);
      mxArray *this_option_value = NULL;
      if (mxGetN(this_option) > 1)
        this_option_value = mxGetCell(this_option,1);

      if (option_name == "horizon")
        if (this_option_value && mxIsNumeric(this_option_value))
          {
            temp_buf = (double *) mxGetData(this_option_value);
            options->horizon = (int) temp_buf[0];
          }
        else
          {
            cout << "You must pass an integer after specifying the 'horizon' option" << endl;
            return 1;
          }
      else if (option_name == "filtered")
        options->filtered_probabilities = true;
      else if (option_name == "error_bands")
        {
          free(options->percentiles);
          options->num_percentiles = 3;
          options->percentiles = new double[3];
          options->percentiles[0] = 0.16;
          options->percentiles[1] = 0.5;
          options->percentiles[2] = 0.84;

          // Check if the user specified to turn off error bands
          if (this_option_value && mxIsNumeric(this_option_value))
            {
              temp_buf = (double *) mxGetData(this_option_value);
              if (temp_buf[0] == 0)
                {
                  options->num_percentiles = 1;
                  options->percentiles = new double[1];
                  options->percentiles[0] = 0.50;
                }
            }
        }
      else if (option_name == "median")
        {
          free(options->percentiles);
          options->num_percentiles = 1;
          options->percentiles = new double[1];
          options->percentiles[0] = 0.5;
        }
      else if (option_name == "percentiles")
        if (this_option_value)
          {
            options->num_percentiles = mxGetN(this_option_value)
              > mxGetM(this_option_value) ? mxGetN(this_option_value)
              : mxGetM(this_option_value);
            options->percentiles = mxGetPr(this_option_value);
          }
        else
          {
            cout << "You must pass a vector after the 'percentiles' argument with the "
                 << "percentiles that you want to have computed, ex "
                 << "'percentiles',[.16 .5 .84]" << endl;
            return 1;
          }
      else if (option_name == "parameter_uncertainty")
        {
          options->parameter_uncertainty = true;
          if (shocks_passed == false)
            options->shocks = 1;
        }
      else if (option_name == "shocks_per_parameter")
        if (this_option_value && mxIsNumeric(this_option_value))
          {
            temp_buf = (double *) mxGetData(this_option_value);
            options->shocks = (int) temp_buf[0];
            shocks_passed = true;
          }
        else
          {
            cout << "You must pass an integer after specifying the 'shocks_per_parameter' option" << endl;
            return 1;
          }
      else if (option_name == "thin")
        if (this_option_value && mxIsNumeric(this_option_value))
          {
            temp_buf = (double *) mxGetData(this_option_value);
            options->thin = (int) temp_buf[0];
          }
        else
          {
            cout << "You must pass an integer after specifying the 'thin' option" << endl;
            return 1;
          }
      else if (option_name == "simulation_file")
        {
          char *posterior_filename = mxArrayToString(this_option_value);
          strcpy(options->simulation_filename, posterior_filename);
          if (!(options->simulation_file = fopen(posterior_filename, "rt")))
            {
              cout << "Can not open posterior file " << options->simulation_file
                   << " for reading. " << endl;
              return 1;
            }
          mxFree(posterior_filename);
        }
      else if (option_name == "regimes")
        options->regimes = true;
      else if (option_name == "regime")
        if (this_option_value && mxIsNumeric(this_option_value))
          {
            temp_buf = (double *) mxGetData(this_option_value);
            options->regime = (int) temp_buf[0];
          }
        else
          {
            cout << "You must pass an integer after specifying the 'regime' "
                 << "option, or alternatively you can specify 'regimes'" << endl;
            return 1;
          }
      else if (option_name == "number_observations")
        if (this_option_value && mxIsNumeric(this_option_value))
          {
            temp_buf = (double *) mxGetData(this_option_value);
            options->number_observations = (int) temp_buf[0];
          }
        else
          {
            cout << "You must pass an integer after specifying the "
                 << "'number_observations' option" << endl;
            return 1;
          }
      else if (option_name == "free_parameters")
        if (this_option_value && mxIsNumeric(this_option_value))
          {
            options->num_parameters = mxGetM(this_option_value);
            options->free_parameters = mxGetPr(this_option_value);
          }
        else
          {
            cout << "You must pass a vector of free parameters after "
                 << "specifying 'free_parameters'" << endl;
            return 1;
          }
      else if (option_name == "mean")
        {
          options->mean = true;
          options->num_percentiles = 0;
          options->percentiles = (double *) NULL;

        }
      else if (option_name == "seed")
        if (this_option_value && mxIsNumeric(this_option_value))
          {
            temp_buf = (double *) mxGetData(this_option_value);
            options->seed = (long) temp_buf[0];
          }
        else
          {
            cout << "You must pass an integer after specifying the 'seed' option" << endl;
            return 1;
          }
      else
          {
            cout << "set_options error: option '" << option_name << "' not matched" << endl;
            return 1;
          }
      mxFree(option_name_c);
    } // End Optional Arguments
  return 0;
}

int
print_sbvar_options(SbvarOption *options)
{
  using namespace std;
  int i;

  if (options == NULL)
    return 1;

  cout << "SBVAR OPTIONS" << endl;
  cout << "Number of Shocks: " << options->shocks << endl;
  cout << "Thinning Factor: " << options->thin << endl;
  cout << "Horizon: " << options->horizon << endl;
  cout << "Number of Observations to Show: " << options->number_observations << endl;
  cout << "Regime: " << options->regime << endl;
  cout << "Regimes: " << options->regimes << endl;
  cout << "Number Free parameters: " << options->num_parameters << endl;
  cout << "Parameter Uncertainty: " << options->parameter_uncertainty << endl;
  cout << "Mean: " << options->mean << endl;
  cout << "Number Percentiles: " << options->num_percentiles << endl;
  cout << "Percentiles: ";
  for (i = 0; i < options->num_percentiles; i++)
    cout << options->percentiles[i] << " ";
  cout << endl;
  cout << "Filtered Probabilities: " << options->filtered_probabilities << endl;
  cout << "Simulation Filename: " << options->simulation_filename << endl;
  cout << "Random Number Seed: " << options->seed << endl;

  return 0;
}

TStateModel *
initialize_model_and_options(SbvarOption **options, const mxArray *prhs[], int *nstates, int *nvars, int *npre, int *nfree)
{
  char *name;
  TStateModel *model = (TStateModel *) NULL;
  SbvarOption *opt;

  // Initialize the StateSpace Model with the initialization file
  name = mxArrayToString(mxGetCell(mxGetCell(prhs[0],0),1));
  model = initialize_ms_model(name);
  if (model == NULL)
    {
      cout << "Could not initialize State Space Switching model with the given tag." << endl;
      return (TStateModel *) NULL;
    }
  if (get_var_dimensions(model, nstates, nvars, npre, nfree) > 0)
    {
      cout << "Problems Determining the size of the Initialized model." << endl;
      return (TStateModel *) NULL;
    }

  // Process the rest of the options
  opt = initialize_sbvar_options(name);
  mxFree(name);
  if (set_options(opt, prhs) > 0)
    {
      cout << "There was a problem with the options passed." << endl;
      return (TStateModel *) NULL;
    }

  // Set seed value
  try
    {
      dw_initialize_generator(opt->seed);
    }
  catch (const char *s)
    {
      cout << "Exception: " << s << endl;
      cout << "Exception thrown initializing Random Seed." << endl;
      return (TStateModel *) NULL;
    }

  // If free parameters have been passed, then set them in the model
  if (opt->num_parameters > 0)
    {
      if (opt->num_parameters == *nfree + 2)
        {
          opt->free_parameters = opt->free_parameters+2;
          opt->num_parameters = opt->num_parameters - 2;
        }
      if (opt->num_parameters != *nfree)
        {
          cout << "The Free Parameter vector passed is the wrong size for the given model" << endl;
          return (TStateModel *) NULL;
        }

      // Set the paramters as the current 'theta'
      if (set_parameters_in_VAR(model, opt->free_parameters) > 0)
        {
          cout << "Problem with the free parameters that were passed being set in MS-SBVAR Model." << endl;
          return (TStateModel *) NULL;
        }
    }

  *options = opt;

  return model;
}

#endif
