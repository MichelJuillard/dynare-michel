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

#include <stdio.h>
#include <stdlib.h>

#include "modify_for_mex.h"

void
mex_write_to_matlab_matfile(double *data, int rows, int cols, const char *varname, const char *filename)
{
  mxArray *toWrite;
  toWrite = getMxArray(data, rows, cols);
  MATFile *matfile = matOpen(filename, "w");
  if (matfile == NULL)
    mexErrMsgTxt("Error encountered in mex when opening a .mat file");

  if (matPutVariable(matfile, varname, toWrite) != 0 ||
      ferror(matGetFp(matfile)) > 0 ||
      feof(matGetFp(matfile)) > 0)
    mexErrMsgTxt("Error encountered in mex when writing a .mat file");

  if (matClose(matfile) != 0)
    mexErrMsgTxt("Error encountered in mex when closing a .mat file");
}

void
mex_write_to_matlab_global_struct(double *data, int rows, int cols, const char *varname)
{
  mxArray *toWrite;
  toWrite = getMxArray(data, rows, cols);

  int fieldnumber = mxAddField(globalMatlabStruct, varname);
  if (fieldnumber < 0)
    mexErrMsgTxt("Error encountered in mex when assigining to global structure");
  mxSetFieldByNumber(globalMatlabStruct, 0, fieldnumber, toWrite);
}

mxArray *
getMxArray(double *data, int rows, int cols)
{
  mxArray *mat;
  double *pind;
  int i;

  if (rows == 1 && cols == 1)
    mat = mxCreateDoubleScalar(data[0]);
  else
    {
      mat = mxCreateDoubleMatrix(rows, cols, mxREAL);
      if (mat != NULL)
        {
          pind = mxGetPr(mat);
          for (i = 0; i < rows*cols; i++)
            pind[i] = data[i];
        }
    }
  if (mat == NULL)
    mexErrMsgTxt("Error encountered in mex: not enough heap space to allocate memory for mxArray");
  return mat;
}
