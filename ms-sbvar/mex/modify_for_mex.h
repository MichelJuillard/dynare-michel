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

#ifndef _MEXMOD
#define _MEXMOD
void swz_exit(int status);
void swz_fprintf_err(const char * str, ...);
extern int constant_seed;
#endif



#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)


#include <dynmex.h>
#include <dynblas.h>
#include <dynlapack.h>

#define swzMalloc mxMalloc
#define swzCalloc mxCalloc
#define swzRealloc mxRealloc
#define swzFree mxFree


#else


#define swz_fprintf_stdout printf
#define swzMalloc malloc
#define swzCalloc calloc
#define swzRealloc realloc
#define swzFree free


#endif
