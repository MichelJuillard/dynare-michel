/*
 * Copyright (C) 2008-2011 Dynare Team
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

// GP, based on work by O.Kamenik

#include <vector>
#include "first_order.h"
#include "dynamic_abstract_class.hh"

#include <cmath>
#include <sstream>

#include "memory_file.h"

#include <iostream>
#include <fstream>

/**************************************************************************************/
/*       Dynare DynamicModel class                                                                 */
/**************************************************************************************/

KordpDynare::KordpDynare(const vector<string> &endo, int num_endo,
                         const vector<string> &exo, int nexog, int npar,
                         Vector &ysteady, TwoDMatrix &vcov, Vector &inParams, int nstat,
                         int npred, int nforw, int nboth, const int jcols, const Vector &nnzd,
                         const int nsteps, int norder,
                         Journal &jr, DynamicModelAC *dynamicModelFile_arg, double sstol,
                         const vector<int> &var_order, const TwoDMatrix &llincidence, double criterium) throw (TLException) :
  nStat(nstat), nBoth(nboth), nPred(npred), nForw(nforw), nExog(nexog), nPar(npar),
  nYs(npred + nboth), nYss(nboth + nforw), nY(num_endo), nJcols(jcols), NNZD(nnzd), nSteps(nsteps),
  nOrder(norder), journal(jr), ySteady(ysteady), params(inParams), vCov(vcov),
  md(1), dnl(*this, endo), denl(*this, exo), dsnl(*this, dnl, denl), ss_tol(sstol), varOrder(var_order),
  ll_Incidence(llincidence), qz_criterium(criterium), dynamicModelFile(dynamicModelFile_arg)
{
  ReorderDynareJacobianIndices();

  //	Initialise ModelDerivativeContainer(*this, this->md, nOrder);
  for (int iord = 1; iord <= nOrder; iord++)
    md.insert(new FSSparseTensor(iord, nY+nYs+nYss+nExog, nY));
}

KordpDynare::~KordpDynare()
{
  // No need to manually delete tensors in "md", they are deleted by the TensorContainer destructor
}

void
KordpDynare::solveDeterministicSteady()
{
  JournalRecordPair pa(journal);
  pa << "Non-linear solver for deterministic steady state By-passed " << endrec;
}

void
KordpDynare::evaluateSystem(Vector &out, const Vector &yy, const Vector &xx) throw (DynareException)
{
  // This method is only called when checking the residuals at steady state (Approximation::check), so return zero residuals
  out.zeros();
}

void
KordpDynare::evaluateSystem(Vector &out, const Vector &yym, const Vector &yy,
                            const Vector &yyp, const Vector &xx) throw (DynareException)
{
  // This method is only called when checking the residuals at steady state (Approximation::check), so return zero residuals
  out.zeros();
}

/************************************************
 * this is main derivative calculation functin that indirectly calls dynamic.dll
 * which performs actual calculation and reorders
 ***************************************************/
void
KordpDynare::calcDerivativesAtSteady() throw (DynareException)
{
  g1p = new TwoDMatrix(nY, nJcols);
  g1p->zeros();

  if (nOrder > 1)
    {
      // allocate space for sparse Hessian
      g2p = new TwoDMatrix((int) NNZD[1], 3);
      g2p->zeros();
    }

  if (nOrder > 2)
    {
      g3p = new TwoDMatrix((int) NNZD[2], 3);
      g3p->zeros();
    }

  Vector xx(nexog());
  xx.zeros();

  Vector out(nY);
  out.zeros();
  Vector llxSteady(nJcols-nExog);
  LLxSteady(ySteady, llxSteady);

  dynamicModelFile->eval(llxSteady, xx, params, ySteady, out, g1p, g2p, g3p);

  populateDerivativesContainer(*g1p, 1, JacobianIndices);
  delete g1p;

  if (nOrder > 1)
    {
      populateDerivativesContainer(*g2p, 2, JacobianIndices);
      delete g2p;
    }
  if (nOrder > 2)
    {
      populateDerivativesContainer(*g3p, 3, JacobianIndices);
      delete g3p;
    }
}

/*******************************************************************************
 * populateDerivatives to sparse Tensor and fit it in the Derivatives Container
 *******************************************************************************/
void
KordpDynare::populateDerivativesContainer(const TwoDMatrix &g, int ord, const vector<int> &vOrder)
{
  // model derivatives FSSparseTensor instance
  FSSparseTensor *mdTi = (new FSSparseTensor(ord, nJcols, nY));

  IntSequence s(ord, 0);

  if (ord == 1)
    {
      for (int i = 0; i < g.ncols(); i++)
        {
          for (int j = 0; j < g.nrows(); j++)
            {
              double x;
              if (s[0] < nJcols-nExog)
                x = g.get(j, vOrder[s[0]]);
              else
                x = g.get(j, s[0]);
              if (x != 0.0)
                mdTi->insert(s, j, x);
            }
          s[0]++;
        }
    }
  else if (ord == 2)
    {
      int nJcols1 = nJcols-nExog;
      vector<int> revOrder(nJcols1);
      for (int i = 0; i < nJcols1; i++)
        revOrder[vOrder[i]] = i;
      for (int i = 0; i < g.nrows(); i++)
        {
          int j = (int) g.get(i, 0)-1; // hessian indices start with 1
          int i1 = (int) g.get(i, 1) -1;
          int s0 = i1 / nJcols;
          int s1 = i1 % nJcols;
          if (s0 < nJcols1)
            s[0] = revOrder[s0];
          else
            s[0] = s0;
          if (s1 < nJcols1)
            s[1] = revOrder[s1];
          else
            s[1] = s1;
          if (s[1] >= s[0])
            {
              double x = g.get(i, 2);
              mdTi->insert(s, j, x);
            }
        }
    }
  else if (ord == 3)
    {
      int nJcols1 = nJcols-nExog;
      int nJcols2 = nJcols*nJcols;
      vector<int> revOrder(nJcols1);
      for (int i = 0; i < nJcols1; i++)
        revOrder[vOrder[i]] = i;
      for (int i = 0; i < g.nrows(); i++)
        {
          int j = (int) g.get(i, 0)-1;
          int i1 = (int) g.get(i, 1) -1;
          int s0 = i1 / nJcols2;
          int i2 = i1 % nJcols2;
          int s1 = i2 / nJcols;
          int s2 = i2 % nJcols;
          if (s0 < nJcols1)
            s[0] = revOrder[s0];
          else
            s[0] = s0;
          if (s1 < nJcols1)
            s[1] = revOrder[s1];
          else
            s[1] = s1;
          if (s2 < nJcols1)
            s[2] = revOrder[s2];
          else
            s[2] = s2;
          double x = g.get(i, 2);
          if (s.isSorted())
            mdTi->insert(s, j, x);
        }
    }

  // md container
  md.remove(Symmetry(ord));
  md.insert(mdTi);
  // No need to delete mdTi, it will be deleted by TensorContainer destructor
}

/*********************************************************
 * LLxSteady()
 * returns ySteady extended with leads and lags suitable for
 * passing to <model>_dynamic DLL
 *************************************************************/
void
KordpDynare::LLxSteady(const Vector &yS, Vector &llxSteady) throw (DynareException, TLException)
{
  if ((nJcols-nExog) == yS.length())
    throw DynareException(__FILE__, __LINE__, "ySteady already of right size");

  // create temporary square 2D matrix size nEndo x nEndo (sparse)
  // for the lag, current and lead blocks of the jacobian
  if (llxSteady.length() != nJcols-nExog)
    throw DynareException(__FILE__, __LINE__, "llxSteady has wrong size");

  for (int ll_row = 0; ll_row < ll_Incidence.nrows(); ll_row++)
    {
      // populate (non-sparse) vector with ysteady values
      for (int i = 0; i < nY; i++)
        {
          if (ll_Incidence.get(ll_row, i))
            llxSteady[((int) ll_Incidence.get(ll_row, i))-1] = yS[i];
        }
    }
}

/************************************
 * Reorder DynareJacobianIndices of variables in a vector according to
 * given int * varOrder together with lead & lag incidence matrix and
 * any the extra columns for exogenous vars, and then,
 * reorders its blocks given by the varOrder and the Dynare++ expectations:

 * extra	nboth+ npred (t-1) lags
 * varOrder
                static:
    pred
    both
    forward
    * extra both + nforw (t+1) leads, and
    * extra exogen

    * so to match the jacobian organisation expected by the Appoximation class
      both + nforw (t+1) leads
      static
      pred
      both
      forward
      nboth+ npred  (t-1) lags
      exogen
************************************/

void
KordpDynare::ReorderDynareJacobianIndices() throw (TLException)
{
  // create temporary square 2D matrix size nEndo x nEndo (sparse)
  // for the lag, current and lead blocks of the jacobian
  JacobianIndices.resize(nJcols);
  vector <int> tmp(nY);
  int i, j, rjoff = nJcols-nExog-1;

  for (int ll_row = 0; ll_row < ll_Incidence.nrows(); ll_row++)
    {
      // reorder in orde-var order & populate temporary nEndo (sparse) vector with
      // the lag, current and lead blocks of the jacobian respectively
      for (i = 0; i < nY; i++)
        tmp[i] = ((int) ll_Incidence.get(ll_row, varOrder[i]-1));
      // write the reordered blocks back to the jacobian
      // in reverse order
      for (j = nY-1; j >= 0; j--)
        if (tmp[j])
          {
            JacobianIndices[rjoff] = tmp[j] -1;
            rjoff--;
            if (rjoff < 0)
              break;
          }
    }

  //add the indices for the nExog exogenous jacobians
  for (j = nJcols-nExog; j < nJcols; j++)
    JacobianIndices[j] = j;
}

/**************************************************************************************/
/*       DynareNameList class                                                         */
/**************************************************************************************/

DynareNameList::DynareNameList(const KordpDynare &dynare, const vector<string> &names_arg) : names(names_arg)
{
}

DynareStateNameList::DynareStateNameList(const KordpDynare &dynare, const DynareNameList &dnl,
                                         const DynareNameList &denl)
{
  for (int i = 0; i < dynare.nys(); i++)
    names.push_back(string(dnl.getName(i+dynare.nstat())));
  for (int i = 0; i < dynare.nexog(); i++)
    names.push_back(string(denl.getName(i)));
}
