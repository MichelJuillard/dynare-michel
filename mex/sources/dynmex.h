/*
 * Copyright (C) 2009-2011 Dynare Team
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

#ifndef _DYNMEX_H
#define _DYNMEX_H

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE)
# error You must define either MATLAB_MEX_FILE or OCTAVE_MEX_FILE
#endif

#include <mex.h>

/* mwSize, mwIndex and mwSignedIndex appeared in MATLAB 7.3 */
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION < 0x0703
typedef int mwIndex;
typedef int mwSize;
#endif

/*
 * Fix for trac ticket Ticket #137
 */
#if !defined(DYN_MEX_FUNC_ERR_MSG_TXT)
# define DYN_MEX_FUNC_ERR_MSG_TXT(str)          \
  do {                                          \
    mexPrintf("%s\n", str);                     \
    int i;                                      \
    for (i = 0; i < nlhs; i++)                  \
      plhs[i] = mxCreateDoubleScalar(1);        \
    return;                                     \
  } while (0)
#endif

#endif
