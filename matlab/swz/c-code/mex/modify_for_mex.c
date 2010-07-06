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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdarg.h>
#include <string.h>

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#include <dynmex.h>
#endif

int constant_seed;

void
swz_fprintf_err(char *str, ...)
{
  va_list ap;
  va_start(ap, str);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexPrintf(str, ap);
#else
  vfprintf(stderr, str, ap);
#endif

  va_end(ap);
}

void
swzExit(int status)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexErrMsgTxt("Error in mexfile.\n");
#else
  exit(status);
#endif
}
