/*
 * Copyright (C) 2010-2011 Dynare Team
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

#ifndef _MEXMOD
#define _MEXMOD

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)

#include <dynmex.h>
#include <dynblas.h>
#include <dynlapack.h>

#define dw_malloc mxMalloc
#define dw_calloc mxCalloc
#define dw_realloc mxRealloc
#define dw_free mxFree
#define dw_exit msExit

void msExit(int status);
extern int constant_seed;

/* Write Matlab Output */
mxArray *globalMatlabStruct;
void mex_write_to_matlab_matfile(double *, int, int, const char *, const char *);
void mex_write_to_matlab_global_struct(double *, int, int, const char *);
mxArray *getMxArray(double *, int, int);

#endif
#endif
