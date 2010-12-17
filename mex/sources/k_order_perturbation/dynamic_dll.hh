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

#ifndef _DYNAMIC_DLL_HH
#define _DYNAMIC_DLL_HH

#if defined(_WIN32) || defined(__CYGWIN32__)
# define NOMINMAX // Do not define "min" and "max" macros
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

#include <string>

#include "dynamic_abstract_class.hh"
#include "dynare_exception.h"

// <model>_Dynamic DLL pointer
typedef void (*DynamicDLLFn)
(const double *y, const double *x, int nb_row_x, const double *params, const double *steady_state,
 int it_, double *residual, double *g1, double *g2, double *g3);

/**
 * creates pointer to Dynamic function inside <model>_dynamic.dll
 * and handles calls to it.
 **/
class DynamicModelDLL : public DynamicModelAC
{
private:
  DynamicDLLFn Dynamic; // pointer to the Dynamic function in DLL
#if defined(_WIN32) || defined(__CYGWIN32__)
  HINSTANCE dynamicHinstance;  // DLL instance pointer in Windows
#else
  void *dynamicHinstance; // and in Linux or Mac
#endif

public:
  // construct and load Dynamic model DLL
  DynamicModelDLL(const string &fname, const string &sExt) throw (DynareException);
  virtual ~DynamicModelDLL();

  void eval(const Vector &y, const Vector &x, const Vector &params, const Vector &ySteady,
            Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException);
};
#endif
