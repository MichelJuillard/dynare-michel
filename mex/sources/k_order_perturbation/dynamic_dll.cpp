/*
 * Copyright (C) 2008-2009 Dynare Team
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

#include "k_ord_dynare.h"
#include "dynamic_dll.h"

#include <sstream>

/***********************************
 * Members of DynamicModelDLL for handling loading and calling
 * <model>_dynamic () function
 **************************************/
DynamicModelDLL::DynamicModelDLL(const string &modName, const int y_length, const int j_cols,
                                 const int n_max_lag, const int n_exog, const string &sExt) throw (DynareException)
  : length(y_length), jcols(j_cols), nMax_lag(n_max_lag), nExog(n_exog)
{
  string fName;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  fName = "./";
#endif
  fName += modName + "_dynamic" + sExt;

  try
    {
#if defined(__CYGWIN32__) || defined(_WIN32)
      HINSTANCE dynamicHinstance;
      dynamicHinstance = ::LoadLibrary(fName.c_str());
      if (dynamicHinstance == NULL)
        throw 1;
      Dynamic = (DynamicFn *) ::GetProcAddress(dynamicHinstance, "Dynamic");

#else // Linux or Mac
      dynamicHinstance = dlopen(fName.c_str(), RTLD_NOW);
      if ((dynamicHinstance == NULL) || dlerror())
        {
          cerr << dlerror() << endl;
          throw 1;
        }
      Dynamic = (DynamicFn) dlsym(dynamicHinstance, "Dynamic");
      if ((Dynamic  == NULL) || dlerror())
        {
          cerr << dlerror() << endl;
          throw 2;
        }
#endif

    }
  catch (int i)
    {
      ostringstream msg;
      msg << "Can't load " << fName << " (error code: " << i;
      throw DynareException(__FILE__, __LINE__, msg.str());
    }
  catch (...)
    {
      throw DynareException(__FILE__, __LINE__, string("Can't load ") + fName);
    }
}

DynamicModelDLL::~DynamicModelDLL()
{
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  FreeLibrary(dynamicHinstance);
#else
  dlclose(dynamicHinstance);
#endif
}

void
DynamicModelDLL::eval(double *y, double *x, int nb_row_x, double *params,
                      int it_, double *residual, double *g1, double *g2, double *g3)
{
  Dynamic(y, x, nb_row_x, params, it_, residual, g1, g2, g3);
}

void
DynamicModelDLL::eval(const Vector &y, const TwoDMatrix &x, const  Vector *modParams,
                      int it_, Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException)
{
  double  *dresidual, *dg1 = NULL, *dg2 = NULL, *dg3 = NULL;

  if ((jcols-nExog) != y.length())
    throw DynareException(__FILE__, __LINE__, "DLL Error: (jcols-nExog)!=ys.length()");

  if (g1 != NULL)
    {
      if (g1->nrows() != length)  // dummy
        {
          delete g1;
          g1 =     new TwoDMatrix(length, jcols); // and get a new one
          g1->zeros();
        }
      dg1 = const_cast<double *>(g1->base());
    }
  if (g2 != NULL)
      dg2 = const_cast<double *>(g2->base());
  dresidual = const_cast<double *>(residual.base());
  if (g3 != NULL)
      dg3 = const_cast<double *>(g3->base());
  dresidual = const_cast<double *>(residual.base());
  double *dy = const_cast<double *>(y.base());
  double *dx = const_cast<double *>(x.base());
  double *dbParams = const_cast<double *>(modParams->base());

      Dynamic(dy, dx, nExog, dbParams, it_, dresidual, dg1, dg2, dg3);
}

void
DynamicModelDLL::eval(const Vector &y, const TwoDMatrix &x, const Vector *modParams,
                      Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException)
{
  eval(y, x, modParams, nMax_lag, residual, g1, g2, g3);
}

void
DynamicModelDLL::eval(const Vector &y, const Vector &x, const Vector *modParams,
                      Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException)
{
  /** ignore given exogens and create new 2D x matrix since
   * when calling <model>_dynamic(z,x,params,it_) x must be equal to
   * zeros(M_.maximum_lag+1,M_.exo_nbr)
   **/
  TwoDMatrix &mx = *(new TwoDMatrix(nMax_lag+1, nExog));
  mx.zeros(); // initialise shocks to 0s

  eval(y, mx, modParams, nMax_lag, residual, g1, g2, g3);
}
