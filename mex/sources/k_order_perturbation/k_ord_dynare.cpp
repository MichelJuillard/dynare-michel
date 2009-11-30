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

// GP, based on work by O.Kamenik

#include <vector>
#include "first_order.h"
#include "k_ord_dynare.h"
#include "dynamic_dll.h"

#include <cmath>
#include <sstream>

#include "memory_file.h"

/**************************************************************************************/
/*       Dynare DynamicModel class                                                                 */
/**************************************************************************************/

KordpDynare::KordpDynare(const char **endo,  int num_endo,
                         const char **exo, int nexog, int npar, //const char** par,
                         Vector *ysteady, TwoDMatrix *vcov, Vector *inParams, int nstat,
                         int npred, int nforw, int nboth, const int jcols, const Vector *nnzd, 
                         const int nsteps, int norder, //const char* modName,
                         Journal &jr, DynamicModelDLL &dynamicDLL, double sstol, 
                         const vector<int> *var_order, const TwoDMatrix *llincidence, double criterium) throw (TLException)
  : nStat(nstat), nBoth(nboth), nPred(npred), nForw(nforw), nExog(nexog), nPar(npar),
  nYs(npred + nboth), nYss(nboth + nforw), nY(num_endo), nJcols(jcols), NNZD(nnzd), nSteps(nsteps), 
  nOrder(norder), journal(jr),  dynamicDLL(dynamicDLL), ySteady(ysteady), vCov(vcov), params(inParams),
  md(1), dnl(NULL), denl(NULL), dsnl(NULL), ss_tol(sstol), varOrder(var_order),
    ll_Incidence(llincidence), qz_criterium(criterium)
{
      dnl = new DynareNameList(*this, endo);
      denl = new DynareExogNameList(*this, exo);
      dsnl = new DynareStateNameList(*this, *dnl, *denl);

      JacobianIndices = ReorderDynareJacobianIndices(varOrder);

      //	Initialise ModelDerivativeContainer(*this, this->md, nOrder);
      for (int iord = 1; iord <= nOrder; iord++)
        {
          FSSparseTensor *t = new FSSparseTensor(iord, nY+nYs+nYss+nExog, nY);
          md.insert(t);
        }
}

KordpDynare::KordpDynare(const KordpDynare &dynare)
  : nStat(dynare.nStat), nBoth(dynare.nBoth),       nPred(dynare.nPred),
  nForw(dynare.nForw), nExog(dynare.nExog),  nPar(dynare.nPar),
  nYs(dynare.nYs), nYss(dynare.nYss), nY(dynare.nY), nJcols(dynare.nJcols),
  NNZD(dynare.NNZD), nSteps(dynare.nSteps), nOrder(dynare.nOrder), 
  journal(dynare.journal), dynamicDLL(dynare.dynamicDLL),
  ySteady(NULL), params(NULL), vCov(NULL), md(dynare.md),
  dnl(NULL), denl(NULL), dsnl(NULL), ss_tol(dynare.ss_tol),
  varOrder(dynare.varOrder), ll_Incidence(dynare.ll_Incidence),
  JacobianIndices(dynare.JacobianIndices), qz_criterium(dynare.qz_criterium)
{
  ySteady = new Vector(*(dynare.ySteady));
  params = new Vector(*(dynare.params));
  vCov = new TwoDMatrix(*(dynare.vCov));
  dnl = new DynareNameList(dynare);
  denl = new DynareExogNameList(dynare);
  dsnl = new DynareStateNameList(*this, *dnl, *denl);
}

KordpDynare::~KordpDynare()
{
  if (ySteady)
    delete ySteady;
  if (params)
    delete params;
  if (vCov)
    delete vCov;
  if (dnl)
    delete dnl;
  if (dsnl)
    delete dsnl;
  if (denl)
    delete denl;
}

/** This clears the container of model derivatives and initializes it
 * inserting empty sparse tensors up to the given order. */
ModelDerivativeContainer::ModelDerivativeContainer(const KordpDynare &model,
                                                   TensorContainer<FSSparseTensor> &mod_ders, int order) : md(mod_ders)
{
  md.clear();
  for (int iord = 1; iord <= order; iord++)
    {
      FSSparseTensor *t = new FSSparseTensor(iord, model.ny()+model.nys()+model.nyss()+model.nexog(), model.ny());
      md.insert(t);
    }
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
  // Nothing to do, this method is only implemented to complete the interface of DynamicModel
}

void
KordpDynare::evaluateSystem(Vector &out, const Vector &yym, const Vector &yy,
                            const Vector &yyp, const Vector &xx) throw (DynareException)
{
  // Nothing to do, this method is only implemented to complete the interface of DynamicModel
}

/************************************************
 * this is main derivative calculation functin that indirectly calls dynamic.dll
 * which performs actual calculation and reorders
 ***************************************************/
void
KordpDynare::calcDerivatives(const Vector &yy, const Vector &xx) throw (DynareException)
{
  TwoDMatrix *g2 = NULL;
  TwoDMatrix *g3 = NULL;
  TwoDMatrix *g1 = new TwoDMatrix(nY, nJcols); // generate g1 for jacobian
  g1->zeros();

  if ((nJcols != g1->ncols()) || (nY != g1->nrows()))
    throw DynareException(__FILE__, __LINE__, "Error in calcDerivatives: Created wrong jacobian");

  if (nOrder > 1)
    {
      // generate g2 space for sparse Hessian 3x NNZH = 3x NNZD[1]
      g2 = new TwoDMatrix((int) (*NNZD)[1],3);
      g2->zeros();
    }
  if (nOrder > 2)
    {
      g3 = new TwoDMatrix((int) (*NNZD)[2],3);
      g3->zeros();
    }
  Vector out(nY);
  out.zeros();
  const Vector *llxYYp; // getting around the constantness
  if ((nJcols - nExog) > yy.length())
      llxYYp =  (LLxSteady(yy));
  else
      llxYYp = &yy;
  const Vector &llxYY = *(llxYYp);

  dynamicDLL.eval(llxYY,  xx, params, out, g1, g2, g3);

  if ((nJcols != g1->ncols()) || (nY != g1->nrows()))
    throw DynareException(__FILE__, __LINE__, "Error in calcDerivatives: dynamicDLL.eval returned wrong jacobian");

  populateDerivativesContainer(g1, 1, JacobianIndices);
  if (nOrder > 1)
      populateDerivativesContainer(g2, 2, JacobianIndices);
  if (nOrder > 2)
      populateDerivativesContainer(g3, 3, JacobianIndices);
}

void
KordpDynare::calcDerivativesAtSteady() throw (DynareException)
{
  Vector xx(nexog());
  xx.zeros();
  calcDerivatives(*ySteady, xx);
}

/*******************************************************************************
* populateDerivatives to sparse Tensor and fit it in the Derivatives Container
*******************************************************************************/
void
KordpDynare::populateDerivativesContainer(TwoDMatrix *g, int ord, const vector<int> *vOrder)
{
  // model derivatives FSSparseTensor instance
  FSSparseTensor *mdTi = (new FSSparseTensor(ord, nJcols, nY));

  IntSequence s(ord, 0);

  if (ord == 1)
    {
      for (int i = 0; i < g->ncols(); i++)
  	  {
	    for (int j = 0; j < g->nrows(); j++)
	      {
	      double x;
	      if (s[0] < nJcols-nExog)
      		x = g->get(j, (*vOrder)[s[0]]);
	      else
		      x = g->get(j, s[0]);
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
      revOrder[(*vOrder)[i]] = i;
    for (int i = 0; i < g->nrows(); i++)
    	{
	    int j = (int)g->get(i,0)-1; // hessian indices start with 1
	    int i1 = (int)g->get(i,1) -1;
	    int s0 = (int)floor(((double) i1)/((double) nJcols));
      int s1 = i1- (nJcols*s0);
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
	      double x = g->get(i,2);
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
	revOrder[(*vOrder)[i]] = i;
      for (int i = 0; i < g->nrows(); i++)
	{
	  int j = (int)g->get(i,0)-1; 
	  int i1 = (int)g->get(i,1) -1;
	  int s0 = (int)floor(i1/nJcols2);
	  int i2 = i1 - nJcols2*s0;
	  int s1 = (int)floor(i2/nJcols);
	  int s2 = i2 - nJcols*s1;
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
	  if ((s[2] >= s[1]) && (s[1] >= s[0]))
	    {
	      double x = g->get(i,2);
	      mdTi->insert(s, j, x);
	    }
    	}
    }

  // md container
  md.remove(Symmetry(ord));
  md.insert(mdTi);
}

void
KordpDynare::writeModelInfo(Journal &jr) const
{
  // write info on variables
    JournalRecordPair rp(journal);
    rp << "Information on variables" << endrec;
    JournalRecord rec1(journal);
    rec1 << "Number of endogenous:            " << ny() << endrec;
    JournalRecord rec2(journal);
    rec2 << "Number of exogenous:             " << nexog() << endrec;
    JournalRecord rec3(journal);
    rec3 << "Number of static:                " << nstat() << endrec;
    JournalRecord rec4(journal);
    rec4 << "Number of predetermined:         " << npred()+nboth() << endrec;
    JournalRecord rec5(journal);
    rec5 << "Number of forward looking:       " << nforw()+nboth() << endrec;
    JournalRecord rec6(journal);
    rec6 << "Number of both:                  " << nboth() << endrec;
}

/*********************************************************
 * LLxSteady()
 * returns ySteady extended with leads and lags suitable for
 * passing to <model>_dynamic DLL
 *************************************************************/
Vector *
KordpDynare::LLxSteady(const Vector &yS) throw (DynareException, TLException)
{
  if ((nJcols-nExog) == yS.length())
    throw DynareException(__FILE__, __LINE__, "ySteady already of right size");

  // create temporary square 2D matrix size nEndo x nEndo (sparse)
  // for the lag, current and lead blocks of the jacobian
  Vector *llxSteady = new Vector(nJcols-nExog);
      for (int ll_row = 0; ll_row < ll_Incidence->nrows(); ll_row++)
        {
          // populate (non-sparse) vector with ysteady values
          for (int i = 0; i < nY; i++)
            {
              if (ll_Incidence->get(ll_row, i))
                (*llxSteady)[((int) ll_Incidence->get(ll_row, i))-1] = yS[i];
            }
        }

  return llxSteady;
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

vector<int> *
KordpDynare::ReorderDynareJacobianIndices(const vector<int> *varOrder) throw (TLException)
{
  // create temporary square 2D matrix size nEndo x nEndo (sparse)
  // for the lag, current and lead blocks of the jacobian
  vector<int> *JacobianIndices = new vector<int>(nJcols);
  vector <int> tmp(nY);
  int i, j, rjoff = nJcols-nExog-1;

      for (int ll_row = 0; ll_row < ll_Incidence->nrows(); ll_row++)
        {
          // reorder in orde-var order & populate temporary nEndo (sparse) vector with
          // the lag, current and lead blocks of the jacobian respectively
          for (i = 0; i < nY; i++)
              tmp[i] = ((int) ll_Incidence->get(ll_row, (*varOrder)[i]-1));
          // write the reordered blocks back to the jacobian
          // in reverse order
          for (j = nY-1; j >= 0; j--)
              if (tmp[j])
                {
                  (*JacobianIndices)[rjoff] = tmp[j] -1;
                  rjoff--;
                  if (rjoff < 0)
                      break;
                }
        }

  //add the indices for the nExog exogenous jacobians
  for (j = nJcols-nExog; j < nJcols; j++)
      (*JacobianIndices)[j] = j;

  return JacobianIndices;
}

/************************************
* Reorder first set of columns of variables in a (jacobian) matrix
* according to order given in  varsOrder together with the extras
* assuming tdx ncols() - nExog is eaqual or less than length of varOrder and
* of any of its elements too.
************************************/

void
KordpDynare::ReorderCols(TwoDMatrix *tdx, const vector<int> *vOrder) throw (DynareException, TLException)
{

  if (tdx->ncols() > vOrder->size())
    throw DynareException(__FILE__, __LINE__, "Size of order var is too small");

  TwoDMatrix tmp(*tdx); // temporary 2D matrix
  TwoDMatrix &tmpR = tmp;
  tdx->zeros(); // empty original matrix
  // reorder the columns

      for (int i = 0; i < tdx->ncols(); i++)
        tdx->copyColumn(tmpR, (*vOrder)[i], i);
}

void
KordpDynare::ReorderCols(TwoDMatrix *tdx, const int *vOrder) throw (TLException)
{

  TwoDMatrix tmp(*tdx); // temporary 2D matrix
  TwoDMatrix &tmpR = tmp;
  tdx->zeros(); // empty original matrix
  // reorder the columns
      for (int i = 0; i < tdx->ncols(); i++)
        tdx->copyColumn(tmpR, vOrder[i], i);
}

/***********************************************************************
* Recursive hierarchical block reordering of the higher order, input model
*	derivatives inc. Hessian
* This is now obsolete but kept in in case it is needed
***********************************************************************/

void
KordpDynare::ReorderBlocks(TwoDMatrix *tdx, const vector<int> *vOrder) throw (DynareException, TLException)
{
  // determine order of the matrix

  double dbOrder = log((double) tdx->ncols())/log((double) nJcols);
  int ibOrder = (int) dbOrder;
  if ((double) ibOrder != dbOrder || ibOrder > nOrder)
    {
      ostringstream msg;
      msg << "Wrong order " << dbOrder;
      throw DynareException(__FILE__, __LINE__, msg.str());
    }

  TwoDMatrix tmp(*tdx); // temporary 2D matrix
  TwoDMatrix &tmpR = tmp;
  tdx->zeros(); // empty original matrix

  if (ibOrder > 1)
    {
      int nBlocks = tmp.ncols()/ nJcols;
      int bSize = tmp.ncols()/nBlocks;
      for (int j = 0; j < nBlocks;  ++j)
        {
          TwoDMatrix subtdx(tmpR, bSize*((*vOrder)[j]), bSize);
          ReorderBlocks(&subtdx, vOrder);
          tdx->place(subtdx, 0, bSize*j);
        }
    }
  else
    {
      if (tdx->ncols() > vOrder->size())
        throw DynareException(__FILE__, __LINE__, "Size of order var is too small");

      // reorder the columns
          for (int i = 0; i < tdx->ncols(); i++)
            tdx->copyColumn(tmpR, (*vOrder)[i], i);
    }
}

void
KordpVectorFunction::eval(const ConstVector &in, Vector &out) throw (DynareException)
{
  check_for_eval(in, out);
  Vector xx(d.nexog());
  xx.zeros();
  d.evaluateSystem(out, in, xx);
}

/**************************************************************************************/
/*       DynareNameList class                                                         */
/**************************************************************************************/
vector<int>
DynareNameList::selectIndices(const vector<const char *> &ns) const throw (DynareException)
{
  vector<int> res;
  for (unsigned int i = 0; i < ns.size(); i++)
    {
      int j = 0;
      while (j < getNum() && strcmp(getName(j), ns[i]) != 0)
        j++;
      if (j == getNum())
        throw DynareException(__FILE__, __LINE__,
                              string("Couldn't find name for ") + ns[i]
                              +" in DynareNameList::selectIndices");
      res.push_back(j);
    }
  return res;
}

DynareNameList::DynareNameList(const  KordpDynare &dynare)
{
  for (int i = 0; i < dynare.ny(); i++)
      names.push_back(dynare.dnl->getName(i));
}
DynareNameList::DynareNameList(const KordpDynare &dynare, const char **namesp)
{
  for (int i = 0; i < dynare.ny(); i++)
    names.push_back(namesp[i]);
}

DynareExogNameList::DynareExogNameList(const KordpDynare &dynare)
{
  for (int i = 0; i < dynare.nexog(); i++)
    names.push_back(dynare.denl->getName(i));
}

DynareExogNameList::DynareExogNameList(const KordpDynare &dynare, const char **namesp)
{
  for (int i = 0; i < dynare.nexog(); i++)
    names.push_back(namesp[i]);
}

DynareStateNameList::DynareStateNameList(const KordpDynare &dynare, const DynareNameList &dnl,
                                         const DynareExogNameList &denl)
{
  for (int i = 0; i < dynare.nys(); i++)
      names.push_back(dnl.getName(i+dynare.nstat()));
  for (int i = 0; i < dynare.nexog(); i++)
      names.push_back(denl.getName(i));
}
