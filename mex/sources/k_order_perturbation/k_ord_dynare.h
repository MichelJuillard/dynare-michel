/*
 * Copyright (C) 2005 Ondra Kamenik
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

// based on: $Id: dynare3.h 1764 2008-03-31 14:30:55Z kamenik $
// by 2005, Ondra Kamenik

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
  vector<const char *> names;
public:
  DynareNameList(const KordpDynare &dynare);
  DynareNameList(const KordpDynare &dynare, const char **names);
  int
  getNum() const
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const
  {
    return names[i];
  }
  /** This for each string of the input vector calculates its index
   * in the names. And returns the resulting vector of indices. If
   * the name cannot be found, then an exception is raised. */
  vector<int> selectIndices(const vector<const char *> &ns) const;
};

class DynareExogNameList : public NameList
{
  vector<const char *> names;
public:
  DynareExogNameList(const  KordpDynare &dynare);
  DynareExogNameList(const KordpDynare &dynare, const char **names);
  int
  getNum() const
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const
  {
    return names[i];
  }
};

class DynareStateNameList : public NameList
{
  vector<const char *> names;
public:
  DynareStateNameList(const KordpDynare &dynare, const DynareNameList &dnl,
                      const DynareExogNameList &denl);
  int
  getNum() const
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const
  {
    return names[i];
  }
};
/*********************************************/
// The following only implements DynamicModel with help of ogdyn::DynareModel
// instantiation of pure abstract DynamicModel decl. in dynamic_model.h
class DynamicModelDLL;
class KordpJacobian;
class KordpDynare : public DynamicModel
{
  friend class DynareNameList;
  friend class DynareExogNameList;
  friend class DynareStateNameList;
  friend class KordpDynareJacobian;
  friend class DynamicModelDLL;
  //////////
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
  const Vector *NNZD;  //the total number of non-zero derivative elements
                      //  where hessian is 2nd : NZZD(order=2)
  const int nSteps;
  const int nOrder;
  Journal &journal;
  ///	DynamicModel* model;
  ///const char* modName;
  Vector *ySteady;
  Vector *params;
  TwoDMatrix *vCov;
  TensorContainer<FSSparseTensor> md; // ModelDerivatives
  DynareNameList *dnl;
  DynareExogNameList *denl;
  DynareStateNameList *dsnl;
  const double ss_tol;
  const vector<int> *varOrder;
  const TwoDMatrix *ll_Incidence;
  double qz_criterium;
  vector<int> *JacobianIndices;
public:
  KordpDynare(const char **endo, int num_endo,
              const char **exo, int num_exo, int num_par, //const char** par,
              Vector *ySteady, TwoDMatrix *vCov, Vector *params, int nstat, int nPred,
              int nforw, int nboth, const int nJcols, const Vector *NNZD, 
              const int nSteps, const int ord,    //const char* modName,
              Journal &jr, DynamicModelDLL &dynamicDLL, double sstol,
              const vector<int> *varOrder, const TwoDMatrix *ll_Incidence, 
              double qz_criterium);

  /** Makes a deep copy of the object. */
  KordpDynare(const KordpDynare &dyn);
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
    return *dnl;
  }
  const NameList &
  getStateNames() const
  {
    return *dsnl;
  }
  const NameList &
  getExogNames() const
  {
    return *denl;
  }
  const TwoDMatrix &
  getVcov() const
  {
    return *vCov;
  }
  Vector &
  getParams()
  {
    return *params;
  }

  const TensorContainer<FSSparseTensor> &
  getModelDerivatives() const
  {
    return md;
  }
  const Vector &
  getSteady() const
  {
    return *ySteady;
  }
  Vector &
  getSteady()
  {
    return *ySteady;
  }

  // here is true public interface
  void solveDeterministicSteady(Vector &steady);
  void
  solveDeterministicSteady()
  {
    solveDeterministicSteady(*ySteady);
  }
  void evaluateSystem(Vector &out, const Vector &yy, const Vector &xx);
  void evaluateSystem(Vector &out, const Vector &yym, const Vector &yy,
                      const Vector &yyp, const Vector &xx);
  void calcDerivatives(const Vector &yy, const Vector &xx);
  //void calcDerivatives(const Vector& yy, TwoDMatrix& jj);
  void calcDerivatives(const Vector &yy, ogu::Jacobian &jacob);
  void calcDerivativesAtSteady();
  DynamicModelDLL &dynamicDLL;
  ///	void writeMat4(FILE* fd, const char* prefix) const;
  ///	void writeDump(const std::string& basename) const;
  DynamicModel *
  clone() const
  {
    return new KordpDynare(*this);
  }
  void ReorderCols(TwoDMatrix *tdx, const int *varOrder);
  void ReorderCols(TwoDMatrix *tdx, const vector<int> *varOrder);
  Vector *LLxSteady(const Vector &yS); // returns ySteady extended with leads and lags

private:
  void writeModelInfo(Journal &jr) const;
  int *ReorderDynareJacobianIndices(const int *varOrder);
  vector<int> *ReorderDynareJacobianIndices(const vector<int> *varOrder);
  void ReorderBlocks(TwoDMatrix *tdx, const int *varOrder);
  void ReorderBlocks(TwoDMatrix *tdx, const vector<int> *vOrder);
  void populateDerivativesContainer(TwoDMatrix *g, int ord, const vector<int> *vOrder);
};

/****************************
 * ModelDerivativeContainer manages derivatives container
 ************************************/

class ModelDerivativeContainer //: public ogp::FormulaDerEvalLoader
{
protected:
  //	const ogp::FineAtoms& atoms;
  TensorContainer<FSSparseTensor> &md;
public:
  ModelDerivativeContainer(const KordpDynare &model, TensorContainer<FSSparseTensor> &mod_ders,
                           int order);
  void load(int i, int iord, const int *vars, double res);
};

/****************************
 *  K-Order Perturbation instance of Jacobian:
 ************************************/

class KordpJacobian : public ogu::Jacobian ///, public ogp::FormulaDerEvalLoader
{
protected:
  KordpDynare &dyn;
public:
  KordpJacobian(KordpDynare &dyn);
  virtual ~KordpJacobian()
  {
  }
  // Load <mod>_dynamic.DLL
  //	  void load(const char** modName);
  void eval(const Vector &in);
};

/****************************
 *  K-Order Perturbation instance of VectorFunction:
 ************************************/

class KordpVectorFunction : public ogu::VectorFunction
{
protected:
  KordpDynare &d;
public:
  KordpVectorFunction(KordpDynare &dyn)
    : d(dyn)
  {
  }
  virtual ~KordpVectorFunction()
  {
  }
  int
  inDim() const
  {
    return d.ny();
  }
  int
  outDim() const
  {
    return d.ny();
  }
  void eval(const ConstVector &in, Vector &out);
};

#endif
