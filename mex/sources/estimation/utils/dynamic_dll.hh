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

// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the K_ORDER_PERTURBATION_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// K_ORDER_PERTURBATION_API functions as being imported from a DLL, wheras this DLL sees symbols
// defined with this macro as being exported.

#if defined(_WIN32) || defined(__CYGWIN32__)
# include <windows.h>
# ifdef _MSC_VER
#  define K_ORDER_PERTURBATION_API __declspec(dllexport)
# endif
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

#include <string>
#include "Matrix.hh"

#include "ts_exception.h"

// <model>_Dynamic DLL pointer
typedef void  (*DynamicFn)
(double *y, double *x, int nb_row_x, double *params,
 int it_, double *residual, double *g1, double *g2, double *g3);

/**
* creates pointer to Dynamic function inside <model>_dynamic.dll
* and handles calls to it.
**/
class DynamicModelDLL
{
private:
  DynamicFn  Dynamic; // pointer to the Dynamic function in DLL

  const size_t length;  // tot num vars = Num of Jacobian rows
  const size_t jcols;  // tot num var t-1, t and t+1 instances + exogs = Num of Jacobian columns
  const int nMax_lag; // no of lags
  const size_t nExog; // no of exogenous
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE dynamicHinstance;  // DLL instance pointer in Windows
#else
  void *dynamicHinstance; // and in Linux or Mac
#endif

public:
  // construct and load Dynamic model DLL
  DynamicModelDLL(const std::string &dynamicDllFile, const size_t length, const size_t jcols,
    const int nMax_lag, const size_t nExog) throw (TSException);
  virtual ~DynamicModelDLL();

  // evaluate Dynamic model DLL
  void eval(double *y, double *x, int nb_row_x, double *params,
    int it_, double *residual, double *g1, double *g2, double *g3);
  void eval(const Vector &y, const Vector &x,  const Vector *params,
    Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException);
  void eval(const Vector &y, const Matrix &x,  const Vector *params,
    int it_, Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException);
  void eval(const Vector &y, const Matrix &x,  const Vector *params,
    Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException);
};
