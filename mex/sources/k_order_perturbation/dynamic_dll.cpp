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
#include "math.h"
#include <cstring>

//////////////////////////////////////////////////////
// Convert Matlab Dynare endo and exo names array to C type array of string pointers
// Poblem is that Matlab mx function returns a long string concatenated by columns rather than rows
// hence a rather low level approach is needed
///////////////////////////////////////////////////////
const char **
DynareMxArrayToString(const mxArray *mxFldp, const int len, const int width)
{
  char *cNamesCharStr = mxArrayToString(mxFldp);
  const char **ret = DynareMxArrayToString(cNamesCharStr, len, width);
  return ret;
}

const char **
DynareMxArrayToString(const char *cNamesCharStr, const int len, const int width)
{
  char cNamesMX[len][width+1]; //
#ifdef DEBUG
  mexPrintf("loop DynareMxArrayToString cNamesCharStr = %s \n", cNamesCharStr);
#endif
  for (int i = 0; i < width; i++)
    {
      for (int j = 0; j < len; j++)
        {
          // Allow alphanumeric and underscores "_" only:
          if (isalnum(cNamesCharStr[j+i*len]) || ('_' == cNamesCharStr[j+i*len]))
            {
              cNamesMX[j][i] = cNamesCharStr[j+i*len];
            }
          else cNamesMX[j][i] = '\0';
        }
    }
  const char **ret = (const char **) mxCalloc(len, sizeof(char *));
  for (int j = 0; j < len; j++)
    {
      cNamesMX[j][width] = '\0';
#ifdef DEBUG
      //				mexPrintf("String [%d]= %s \n", j, cNamesMX[j]);
#endif
      char *token = (char *) mxCalloc(strlen(cNamesMX[j])+1, sizeof(char));
      strcpy(token, cNamesMX[j]);
      ret[j] = token;
#ifdef DEBUG
      mexPrintf("ret [%d]= %s \n", j, ret[j]);
#endif
    }
  return ret;
}

/***********************************
 * Members of DynamicModelDLL for handling loading and calling
 * <model>_dynamic () function
 **************************************/
DynamicModelDLL::DynamicModelDLL(const char *modName, const int y_length, const int j_cols,
                                 const int n_max_lag, const int n_exog, const char *sExt)
  : length(y_length), jcols(j_cols), nMax_lag(n_max_lag), nExog(n_exog)
{
  char fName[MAX_MODEL_NAME];
  strcpy(fName, modName);
  using namespace std;
  strcat(fName, "_dynamic");
#ifdef DEBUG
  mexPrintf("MexPrintf: Call Load run DLL %s .\n", fName);
#endif

  try
    {
      if (sExt == NULL)
        sExt = MEXEXT;
#ifdef _WIN32
      HINSTANCE dynamicHinstance;
      //		dynamicHinstance=::LoadLibraryEx(strcat(fNname,"_.dll"),NULL,DONT_RESOLVE_DLL_REFERENCES);//sExt); //"_.dll");
      dynamicHinstance = ::LoadLibrary(strcat(fName, sExt)); //.dll); //"_.dll");
      if (dynamicHinstance == NULL)
        throw 1;                  //alt: return;
      //		(DynamicFn*)	typedef void * (__stdcall *DynamicFn)();
# ifdef DEBUG
      mexPrintf("MexPrintf: Call GetProcAddress  %s .\n", fName);
# endif
      Dynamic = (DynamicFn *) ::GetProcAddress(dynamicHinstance, "Dynamic");

#else // __linux__
      dynamicHinstance = dlopen(strcat(fName, sExt), RTLD_NOW);
      if ((dynamicHinstance == NULL) || dlerror())
        {
          cerr << dlerror() << endl;
          mexPrintf("MexPrintf:Error loading DLL");
          throw 1;
        }
      Dynamic = (DynamicFn) dlsym(dynamicHinstance, "Dynamic");
      if ((Dynamic  == NULL) || dlerror())
        {
          cerr << dlerror() << endl;
          mexPrintf("MexPrintf:Error finding DLL function");
          throw 2;
        }
#endif

    }
  catch (int i)
    {
      mexPrintf("MexPrintf: error in Load and run DLL %s , %d.\n", fName, i);
      mexErrMsgTxt("Err: An error in Load and run DLL  .\n");
      return;

    }
  catch (...)
    {
      mexPrintf("MexPrintf: Unknown error in Call MATLAB function %s.\n", fName);
      mexErrMsgTxt("Err: Unknown error in Load and run DLL  .\n");
      return;
    }
}

// close DLL: If the referenced object was successfully closed,
// close() returns 0, non 0 otherwise
int
DynamicModelDLL::close()
{
#ifdef _WIN32
  // MS FreeLibrary returns non 0 if OK, 0 if fails.
  bool rb = FreeLibrary(dynamicHinstance);
  if (rb)
    return 0;
  else
    return 1;
#else // linux
  //If OK, dlclose() returns 0, non 0 otherwise
  return dlclose(dynamicHinstance);
#endif
};

void
DynamicModelDLL::eval(const Vector &y, const TwoDMatrix &x, const  Vector *modParams,
                      int it_, Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2)
{

  double  *dresidual, *dg1 = NULL, *dg2 = NULL;
  //int length=y.length(); // not!
  if ((jcols-nExog) != y.length())
    {
      // throw DLL Error
      mexPrintf(" DLL Error: (jcols-nExog)!=ys.length() \n");
      return;
    }
  if (residual.length() < length)  // dummy or insufficient
    {
      Vector *tempv = new Vector(length);
      residual = *tempv;
      delete tempv;
      residual.zeros();
    }
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
    {
      dg2 = const_cast<double *>(g2->base());
    }
  dresidual = const_cast<double *>(residual.base());
  double *dy = const_cast<double *>(y.base());
  double *dx = const_cast<double *>(x.base());
  double *dbParams = const_cast<double *>(modParams->base());
#ifdef DEBUG
  mexPrintf(" try eval Dynamic with ne g1: cols=%d , rows=%d\n",
            g1->ncols(), g1->nrows());
  for (int i = 0; i < modParams->length(); i++)
    {
      mexPrintf("k_ord_perturbation: Params[%d]= %g.\n", i, (*modParams)[i]);
    }
  for (int i = 0; i < jcols-nExog; i++)
    {
      mexPrintf("k_ord_perturbation: Ys[%d]= %g.\n", i, dy[i]);
    }
  mexPrintf("k_order_perturbation: call <model> Dynamic dParams= %g ,  , dy = %g dx = %f .\n",
            dbParams[0], dy[0], dx[0]);

#endif
  try
    {
      Dynamic(dy, dx, nExog, dbParams, it_, dresidual, dg1, dg2);
    }
  catch (...)
    {
      mexPrintf("MexPrintf: error in run Dynamic DLL \n");
    }
};

void
DynamicModelDLL::eval(const Vector &y, const TwoDMatrix &x, const Vector *modParams,
                      Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2)
{

  eval(y, x, modParams, nMax_lag, residual, g1, g2);
};

void
DynamicModelDLL::eval(const Vector &y, const Vector &x, const Vector *modParams,
                      Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2)
{

  /** ignore given exogens and create new 2D x matrix since
   * when calling <model>_dynamic(z,x,params,it_) x must be equal to
   * zeros(M_.maximum_lag+1,M_.exo_nbr)
   **/
  TwoDMatrix &mx = *(new TwoDMatrix(nMax_lag+1, nExog));
  mx.zeros(); // initialise shocks to 0s

  eval(y, mx, modParams, nMax_lag, residual, g1, g2);
};


