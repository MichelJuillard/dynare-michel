/*
 * Copyright (C) 2008-2010 Dynare Team
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

#if defined(_WIN32) || defined(__CYGWIN32__)
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

#include <string>
#include <dynmex.h>

#include "dynare_exception.h"

// <model>_Dynamic DLL pointer
typedef void  (*DynamicFn)
(double *y, double *x, int nb_row_x, double *params, double *steady_state,
 int it_, double *residual, double *g1, double *g2, double *g3);

/**
 * creates pointer to Dynamic function inside <model>_dynamic.dll
 * and handles calls to it.
 **/
class DynamicModelDLL
{
private:
  DynamicFn  Dynamic; // pointer to the Dynamic function in DLL

  const int length;  // tot num vars = Num of Jacobian rows
  const int jcols;  // tot num var t-1, t and t+1 instances + exogs = Num of Jacobian columns
  const int nMax_lag; // no of lags
  const int nExog; // no of exogenous
  const Vector &ySteady;
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE dynamicHinstance;  // DLL instance pointer in Windows
#else
  void *dynamicHinstance; // and in Linux or Mac
#endif

public:
  // construct and load Dynamic model DLL
  DynamicModelDLL(const string &fname, const int length, const int jcols,
                  const int nMax_lag, const int nExog, const Vector &ySteady_arg, const string &sExt) throw (DynareException);
  virtual ~DynamicModelDLL();

  // evaluate Dynamic model DLL
  void eval(double *y, double *x, int nb_row_x, double *params,
            int it_, double *residual, double *g1, double *g2, double *g3);
  void eval(const Vector &y, const Vector &x,  const Vector *params,
            Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException);
  void eval(const Vector &y, const TwoDMatrix &x,  const Vector *params,
            int it_, Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException);
  void eval(const Vector &y, const TwoDMatrix &x,  const Vector *params,
            Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException);
};
