/*
 * Copyright (C) 2010-2012 Dynare Team
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
# ifndef NOMINMAX
#  define NOMINMAX // Do not define "min" and "max" macros
# endif
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

#include <string>
#include "Matrix.hh"

#include "ts_exception.h"

// <model>_Dynamic DLL pointer
typedef void (*DynamicFn)
(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state,
 int it_, double *residual, double *g1, double *g2, double *g3);

/**
 * creates pointer to Dynamic function inside <model>_dynamic.dll
 * and handles calls to it.
 **/
class DynamicModelDLL
{
private:
  DynamicFn Dynamic; // pointer to the Dynamic function in DLL
  const size_t n_exog; // no of exogenous
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE dynamicHinstance;  // DLL instance pointer in Windows
#else
  void *dynamicHinstance; // and in Linux or Mac
#endif

public:
  // construct and load Dynamic model DLL
  DynamicModelDLL(const std::string &dynamicDllFile, size_t n_exog_arg) throw (TSException);
  virtual ~DynamicModelDLL();

  //! evaluate Dynamic model DLL
  void eval(const Vector &y, const Matrix &x, const Vector &params, VectorView &ySteady,
            Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException);
  template<class VEC>
  void eval(const Vector &y, const Matrix &x, const VectorView &modParams, VEC &ySteady,
                      Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException)
  {
    Dynamic(y.getData(), x.getData(), 1, modParams.getData(), ySteady.getData(), 0, residual.getData(),
	    g1 == NULL ? NULL : g1->getData(), g2 == NULL ? NULL : g2->getData(), g3 == NULL ? NULL : g3->getData());
  };
};
