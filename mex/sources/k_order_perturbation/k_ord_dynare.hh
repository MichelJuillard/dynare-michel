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

#ifndef K_ORD_DYNARE3_H
#define K_ORD_DYNARE3_H
#include <vector>
#include "t_container.h"
#include "sparse_tensor.h"
#include "decision_rule.h"
#include "dynamic_model.h"

#include "exception.h"
#include "dynare_exception.h"
#include "fs_tensor.h"
#include "SylvException.h"
#include "tl_exception.h"
#include "kord_exception.h"
#include "nlsolve.h"
#include "approximation.h"

class KordpDynare;

/*////////////////////////////////////////////*/
// instantiations of pure abstract class NameList in dynamic_model.h:
/*////////////////////////////////////////////*/
class DynareNameList : public NameList
{
  vector<string> names;
public:
  DynareNameList(const KordpDynare &dynare, const vector<string> &names_arg);
  int
  getNum() const
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const
  {
    return names[i].c_str();
  }
};

class DynareStateNameList : public NameList
{
  vector<string> names;
public:
  DynareStateNameList(const KordpDynare &dynare, const DynareNameList &dnl,
                      const DynareNameList &denl);
  int
  getNum() const
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const
  {
    return names[i].c_str();
  }
};
/*********************************************/
// The following only implements DynamicModel with help of ogdyn::DynareModel
// instantiation of pure abstract DynamicModel decl. in dynamic_model.h
class DynamicModelAC;
class DynamicModelDLL;
class DynamicModelMFile;

class KordpDynare : public DynamicModel
{
  friend class DynareNameList;
  friend class DynareStateNameList;
  friend class DynamicModelDLL;
  friend class DynamicModelMFile;

  const int nStat;
  const int nBoth;
  const int nPred;
  const int nForw;
  const int nExog;
  const int nPar;
  const int nYs; // ={npred + nboth ; }
  const int nYss; // nyss ={ nboth + nforw ; }
  const int nY;  // = num_endo={ nstat + npred + nboth + nforw ; }
  const int nJcols; // no of jacobian columns= nExog+nEndo+nsPred+nsForw
  const Vector &NNZD;  /* the total number of non-zero derivative elements
                          where hessian is 2nd : NZZD(order=2) */
  const int nSteps;
  const int nOrder;
  Journal &journal;
  Vector &ySteady;
  Vector &params;
  TwoDMatrix &vCov;
  TensorContainer<FSSparseTensor> md; // ModelDerivatives
  DynareNameList dnl, denl;
  DynareStateNameList dsnl;
  const double ss_tol;
  const vector<int> &varOrder;
  const TwoDMatrix &ll_Incidence;
  double qz_criterium;
  vector<int> JacobianIndices;

  TwoDMatrix *g1p;
  TwoDMatrix *g2p;
  TwoDMatrix *g3p;
public:
  KordpDynare(const vector<string> &endo, int num_endo,
              const vector<string> &exo, int num_exo, int num_par,
              Vector &ySteady, TwoDMatrix &vCov, Vector &params, int nstat, int nPred,
              int nforw, int nboth, const int nJcols, const Vector &NNZD,
              const int nSteps, const int ord,
              Journal &jr, DynamicModelAC *dynamicModelFile_arg, double sstol,
              const vector<int> &varOrder, const TwoDMatrix &ll_Incidence,
              double qz_criterium) throw (TLException);

  virtual ~KordpDynare();
  int
  nstat() const
  {
    return nStat;
  }
  int
  nboth() const
  {
    return nBoth;
  }
  int
  npred() const
  {
    return nPred;
  }
  int
  nforw() const
  {
    return nForw;
  }
  int
  nexog() const
  {
    return nExog;
  }
  int
  nys() const
  {
    return nYs;
  }
  int
  nyss() const
  {
    return nYss;
  }
  int
  ny() const
  {
    return nY;
  }
  int
  steps() const
  {
    return nSteps;
  }
  int
  order() const
  {
    return nOrder;
  }
  const NameList &
  getAllEndoNames() const
  {
    return dnl;
  }
  const NameList &
  getStateNames() const
  {
    return dsnl;
  }
  const NameList &
  getExogNames() const
  {
    return denl;
  }
  const TwoDMatrix &
  getVcov() const
  {
    return vCov;
  }
  Vector &
  getParams()
  {
    return params;
  }

  const TensorContainer<FSSparseTensor> &
  getModelDerivatives() const
  {
    return md;
  }
  const Vector &
  getSteady() const
  {
    return ySteady;
  }
  Vector &
  getSteady()
  {
    return ySteady;
  }

  void solveDeterministicSteady();
  void evaluateSystem(Vector &out, const Vector &yy, const Vector &xx) throw (DynareException);
  void evaluateSystem(Vector &out, const Vector &yym, const Vector &yy,
                      const Vector &yyp, const Vector &xx) throw (DynareException);
  void calcDerivativesAtSteady() throw (DynareException);
  DynamicModelAC *dynamicModelFile;
  DynamicModel *
  clone() const
  {
    std::cerr << "KordpDynare::clone() not implemented" << std::endl;
    exit(EXIT_FAILURE);
  }
  void LLxSteady(const Vector &yS, Vector &llxSteady) throw (DynareException, TLException); // Given the steady state in yS, returns in llxSteady the steady state extended with leads and lags

private:
  void ReorderDynareJacobianIndices() throw (TLException);
  void populateDerivativesContainer(const TwoDMatrix &g, int ord, const vector<int> &vOrder);
};

#endif
