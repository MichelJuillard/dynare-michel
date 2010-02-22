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

//#include "k_ord_dynare.hh"
#include "dynamic_dll.hh"

#include <sstream>

using namespace std;

/***********************************
* Members of DynamicModelDLL for handling loading and calling
* <model>_dynamic () function
**************************************/
DynamicModelDLL::DynamicModelDLL(const std::string &modName, const int y_length, const int j_cols,
                                 const int n_max_lag, const int n_exog, const std::string &sExt) throw (TSException) :
length(y_length), jcols(j_cols), nMax_lag(n_max_lag), nExog(n_exog)
{
  std::string fName;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  fName = "./";
#endif
  if (sExt.size()>0) //construct modelNmae_dynamic file name with the given extension
    fName += modName + "_dynamic" + sExt;
  else // i.e. an already pre-constructed modelNmae_dynamic file name together with the appropriate extension
    fName += modName;
  // end if

  try
  {
#if defined(__CYGWIN32__) || defined(_WIN32)
    dynamicHinstance = LoadLibrary(fName.c_str());
    if (dynamicHinstance == NULL)
      throw 1;
    Dynamic = (DynamicFn) GetProcAddress(dynamicHinstance, "Dynamic");
    if (Dynamic == NULL)
    {
      FreeLibrary(dynamicHinstance); // Free the library
      throw 2;
    }
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
      dlclose(dynamicHinstance); // Free the library
      cerr << dlerror() << endl;
      throw 2;
    }
#endif

  }
  catch (int i)
  {
    std::ostringstream msg;
    msg << "Error when loading " << fName << " (";
    if (i == 1)
      msg << "can't dynamically load the file";
    if (i == 2)
      msg << "can't locate the 'Dynamic' symbol";
    msg << ")";
    throw TSException(__FILE__, __LINE__, msg.str());
  }
  catch (...)
  {
    throw TSException(__FILE__, __LINE__, std::string("Can't find Dynamic function in ") + fName);
  }
}

DynamicModelDLL::~DynamicModelDLL()
{
#if defined(__CYGWIN32__) || defined(_WIN32)
  bool result = FreeLibrary(dynamicHinstance);
  if (result == 0)
    throw TSException(__FILE__, __LINE__, std::string("Can't free the *_dynamic DLL"));
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
DynamicModelDLL::eval(const Vector &y, const Matrix &x, const  Vector *modParams,
                      int it_, Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException)
{
  double  *dresidual, *dg1 = NULL, *dg2 = NULL, *dg3 = NULL;

  if ((jcols-nExog) != y.getSize())
    throw TSException(__FILE__, __LINE__, "DLL Error: (jcols-nExog)!=ys.length()");

  if (g1 != NULL)
  {
    if (g1->getRows() != length)  // dummy
    {
      delete g1;
      g1 =     new Matrix(length, jcols); // and get a new one
      g1->setAll(0.0); //zero
    }
    dg1 = const_cast<double *>(g1->getData());
  }
  if (g2 != NULL)
    dg2 = const_cast<double *>(g2->getData());
  //dresidual = const_cast<double *>(residual.getData());
  if (g3 != NULL)
    dg3 = const_cast<double *>(g3->getData());
  dresidual = const_cast<double *>(residual.getData());
  double *dy = const_cast<double *>(y.getData());
  double *dx = const_cast<double *>(x.getData());
  double *dbParams = const_cast<double *>(modParams->getData());

  Dynamic(dy, dx, nExog, dbParams, it_, dresidual, dg1, dg2, dg3);
}

void
DynamicModelDLL::eval(const Vector &y, const Matrix &x, const Vector *modParams,
                      Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException)
{
  eval(y, x, modParams, nMax_lag, residual, g1, g2, g3);
}

void
DynamicModelDLL::eval(const Vector &y, const Vector &x, const Vector *modParams,
                      Vector &residual, Matrix *g1, Matrix *g2, Matrix *g3) throw (TSException)
{
  /** ignore given exogens and create new 2D x matrix since
  * when calling <model>_dynamic(z,x,params,it_) x must be equal to
  * zeros(M_.maximum_lag+1,M_.exo_nbr)
  **/
  Matrix &mx = *(new Matrix(nMax_lag+1, nExog));
  mx.setAll(0); // initialise shocks to 0s

  eval(y, mx, modParams, nMax_lag, residual, g1, g2, g3);
}
