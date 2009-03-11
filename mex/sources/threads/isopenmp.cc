/*
 * Copyright (C) 2009 Dynare Team
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

#include "mex.h"

#ifdef _OPENMP
#define OpenMp 1
#else
#define OpenMp 0
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  // Check input and output:
  if ( (nrhs > 0) )
    {
      mexErrMsgTxt("This function has no input arguments!");
    }
  if (nlhs>1)
    {
      mexErrMsgTxt("Too many output arguments.");
    }
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *z ;
  z = mxGetPr(plhs[0]);
  if (OpenMp)
    {
      *z = 1.0;
    }
  else
    {
      *z = 0.0;
    }
}
