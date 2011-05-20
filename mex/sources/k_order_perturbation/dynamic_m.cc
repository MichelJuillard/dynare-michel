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

#include "dynamic_m.hh"

DynamicModelMFile::DynamicModelMFile(const string &modName) throw (DynareException) :
  DynamicMFilename(modName + "_dynamic")
{
}

DynamicModelMFile::~DynamicModelMFile()
{
}

void
DynamicModelMFile::eval(const Vector &y, const Vector &x, const Vector &modParams, const Vector &ySteady,
                        Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException)
{
  mxArray *prhs[nrhs_dynamic], *plhs[nlhs_dynamic];

  prhs[0] = mxCreateDoubleMatrix(y.length(), 1, mxREAL);
  prhs[1] = mxCreateDoubleMatrix(1, x.length(), mxREAL);
  prhs[2] = mxCreateDoubleMatrix(modParams.length(), 1, mxREAL);
  prhs[3] = mxCreateDoubleMatrix(ySteady.length(), 1, mxREAL);
  prhs[4] = mxCreateDoubleScalar(1.0);

  memcpy(mxGetData(prhs[0]), (void *) y.base(), y.length()*sizeof(double));
  memcpy(mxGetData(prhs[1]), (void *) x.base(), x.length()*sizeof(double));
  memcpy(mxGetData(prhs[2]), (void *) modParams.base(), modParams.length()*sizeof(double));
  memcpy(mxGetData(prhs[3]), (void *) ySteady.base(), ySteady.length()*sizeof(double));

  int retVal = mexCallMATLAB(nlhs_dynamic, plhs, nrhs_dynamic, prhs, DynamicMFilename.c_str());
  if (retVal != 0)
    throw DynareException(__FILE__, __LINE__, "Trouble calling " + DynamicMFilename);

  residual = Vector(mxGetPr(plhs[0]), residual.skip(), (int) mxGetM(plhs[0]));
  copyDoubleIntoTwoDMatData(mxGetPr(plhs[1]), g1, (int) mxGetM(plhs[1]), (int) mxGetN(plhs[1]));
  if (g2 != NULL)
    copyDoubleIntoTwoDMatData(unpackSparseMatrix(plhs[2]), g2, (int) mxGetNzmax(plhs[2]), 3);
  if (g3 != NULL)
    copyDoubleIntoTwoDMatData(unpackSparseMatrix(plhs[3]), g3, (int) mxGetNzmax(plhs[3]), 3);

  for (int i = 0; i < nrhs_dynamic; i++)
    mxDestroyArray(prhs[i]);
  for (int i = 0; i < nlhs_dynamic; i++)
    mxDestroyArray(plhs[i]);
}
