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
 
#ifndef _MEXMSSBVAR_H_
#define _MEXMSSBVAR_H_

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C"  {
  #include "modify_for_mex.h"
  #include "switch.h"
  #include "switchio.h"
  #include "VARio.h"
}

typedef struct sbvar_options_t {
  int shocks;
  int thin;
  int horizon;
  int number_observations;
  int regime;
  bool regimes;
  bool parameter_uncertainty;
  int num_percentiles;
  double *percentiles;
  bool filtered_probabilities;
  char *simulation_filename;
  int num_parameters;
  double *free_parameters; 
  FILE *simulation_file;
  bool mean;
  long seed;
} SbvarOption;

#define printf mexPrintf

int file_exist(char *filename);
char* CreateFilenameFromTag(const char *fmt,const  char *tag,const  char *dir);
SbvarOption * initialize_sbvar_options(char *file_tag);
int set_options(SbvarOption *options, const mxArray *prhs[]);
int print_sbvar_options(SbvarOption *options);
TStateModel * initialize_ms_model(char *filename);
int get_var_dimensions(TStateModel *model, int *nstates, int *nvars, int *npre, int *nfree); 
int set_parameters_in_VAR(TStateModel *model, double *free_parameters); 
int convert_free_parameters_to_VAR(TStateModel *model, double *free_parameters, double *A0, double *Aplus, double *Zeta, double *Q );
TStateModel * initialize_model_and_options(SbvarOption **options, const mxArray *prhs[], int *nstates, int *nvars, int *npre, int *nfree);

#endif
#endif
