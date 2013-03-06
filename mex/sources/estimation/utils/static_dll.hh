/*
 * Copyright (C) 2010-2013 Dynare Team
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

// Pointer to the Static function in the MEX
typedef void (*StaticFn)
(const double *y, const double *x, int nb_row_x, const double *params, double *residual, double *g1, double *v2);

/**
 * creates pointer to Dynamic function inside <model>_static.dll
 * and handles calls to it.
 **/
class StaticModelDLL
{
private:
  StaticFn Static; // pointer to the Dynamic function in DLL
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE staticHinstance;  // DLL instance pointer in Windows
#else
  void *staticHinstance; // and in Linux or Mac
#endif

public:
  // construct and load Static model DLL
  StaticModelDLL(const std::string &basename) throw (TSException);
  virtual ~StaticModelDLL();

  //! evaluate Static model DLL
  template<class Vec1, class Vec2, class Vec3, class Mat1>
  void eval(const Vec1 &y, const Mat1 &x, const Vec2 &modParams,
            Vec3 &residual, Matrix *g1, Matrix *v2) throw (TSException)
  {
    assert(y.getStride() == 1);
    assert(x.getLd() == x.getRows());
    assert(modParams.getStride() == 1);
    assert(residual.getStride() == 1);
    assert(g1->getLd() == g1->getRows());
    assert(v2->getLd() == v2->getRows());

    Static(y.getData(), x.getData(), 1, modParams.getData(), residual.getData(),
	    g1 == NULL ? NULL : g1->getData(), v2 == NULL ? NULL : v2->getData());
  };
};
