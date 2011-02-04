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

#include "dynamic_abstract_class.hh"

void
DynamicModelAC::copyDoubleIntoTwoDMatData(double *dm, TwoDMatrix *tdm, int rows, int cols)
{
  int dmIdx = 0;
  for (int j = 0; j < cols; j++)
    for (int i = 0; i < rows; i++)
      tdm->get(i, j) = dm[dmIdx++];
}

double *
DynamicModelAC::unpackSparseMatrix(mxArray *sparseMat)
{
  int totalCols = mxGetN(sparseMat);
  mwIndex *rowIdxVector = mxGetIr(sparseMat);
  mwSize sizeRowIdxVector = mxGetNzmax(sparseMat);
  mwIndex *colIdxVector = mxGetJc(sparseMat);

  double *ptr = mxGetPr(sparseMat);
  double *newMat = (double *) malloc(sizeRowIdxVector*3*sizeof(double));

  int rind = 0;
  int retvalind0 = 0;
  int retvalind1 = sizeRowIdxVector;
  int retvalind2 = sizeRowIdxVector*2;

  for (int i = 0; i < totalCols; i++)
    for (int j = 0; j < (int) (colIdxVector[i+1]-colIdxVector[i]); j++, rind++)
      {
        newMat[retvalind0++] = rowIdxVector[rind] + 1;
        newMat[retvalind1++] = i + 1;
        newMat[retvalind2++] = ptr[rind];
      }

  return newMat;
}
