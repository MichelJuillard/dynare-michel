/* Dynare is free software: you can redistribute it and/or modify
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
// Copyright 2005, Ondra Kamenik

#ifndef K_ORD_DYNARE3_H
#define K_ORD_DYNARE3_H

#include "Dynare_pp/tl/cc/t_container.h"
#include "Dynare_pp/tl/cc/sparse_tensor.h"
#include "Dynare_pp/kord/decision_rule.h"
#include "Dynare_pp/kord/dynamic_model.h"

#include "nlsolve.h"

#include "k_ord_dynare_model.h"

#include <vector>

class Dynare;


// instantiations of pure abstract class NameList in dynamic_model.h:
/*////////////////////////////////////////////
class DynareNameList : public NameList {
	vector<const char*> names;
public:
	DynareNameList(const Dynare& dynare);
	int getNum() const
		{return (int)names.size();}
	const char* getName(int i) const
		{return names[i];}
	/** This for each string of the input vector calculates its index
	 * in the names. And returns the resulting vector of indices. If
	 * the name cannot be found, then an exception is raised. *
	vector<int> selectIndices(const vector<const char*>& ns) const;
};

class DynareExogNameList : public NameList {
	vector<const char*> names;
public:
	DynareExogNameList(const Dynare& dynare);
	int getNum() const
		{return (int)names.size();}
	const char* getName(int i) const
		{return names[i];}
};

class DynareStateNameList : public NameList {
	vector<const char*> names;
public:
	DynareStateNameList(const Dynare& dynare, const DynareNameList& dnl,
						const DynareExogNameList& denl);
	int getNum() const
		{return (int)names.size();}
	const char* getName(int i) const
		{return names[i];}
};
*********************************************/
// The following only implements DynamicModel with help of ogdyn::DynareModel
// instantiation of pure abstract DynamicModel decl. in dynamic_model.h
class DynareJacobian;
class Dynare : public DynamicModel {
/*	friend class DynareNameList;
	friend class DynareExogNameList;
	friend class DynareStateNameList;
	friend class DynareJacobian;
*//////////
	const int nStat;
	const int nBoth;
	const int nPred;
	const int nForw;
	const int nExog;
	const int nPar;
	const int nYs;
	const int nYss;
	const int nY;
	const int nOrder;
///	Journal& journal;
///	ogdyn::DynareModel* model;
///	DynamicModel* model;
	const char* modName;
	Vector* ysteady;
	Vector* params;
	TensorContainer<FSSparseTensor> md;
///	DynareNameList* dnl;  
///	DynareExogNameList* denl;
///	DynareStateNameList* dsnl;

	TwoDMatrix* Vcov;
	TensorContainer<FSSparseTensor>* ModelDerivatives;

///	ogp::FormulaEvaluator* fe;
///	ogp::FormulaDerEvaluator* fde;
	const double ss_tol;
public:
	/** Parses the given model file and uses the given order to
	 * override order from the model file (if it is != -1). */
///	Dynare(const char* modname, int ord, double sstol, Journal& jr);
	/** Parses the given equations with explicitly given names. */
/*///Dynare(const char** endo, int num_endo,
		   const char** exo, int num_exo,
		   const char** par, int num_par,
		   const char* equations, int len, int ord,
		   double sstol, Journal& jr);
*/	Dynare(const char** endo, int num_endo,
		   const char** exo, int num_exo,
		   const char** par, int num_par,
		   const char** modName, int len, int ord,
		   double sstol, Journal& jr);
	/** Makes a deep copy of the object. */
	Dynare(const Dynare& dyn);
	DynamicModel* clone() const
		{return new Dynare(*this);}
	virtual ~Dynare();
	int nstat() const
		{return nStat;}
	int nboth() const
		{return nBoth;}
	int npred() const
		{return nPred;}
	int nforw() const
		{return nForw;}
	int nexog() const
		{return nExo;}
	int nys() const
		{return nYs;}
	int nyss() const
		{return nYss;}
	int ny() const
		{return nY;}
	int order() const
		{return nOrder;}
/****
	const NameList& getAllEndoNames() const
		{return *dnl;}
	const NameList& getStateNames() const
		{return *dsnl;}
	const NameList& getExogNames() const
		{return *denl;}
********/
	TwoDMatrix& getVcov()
		{return Vcov;}
///	const TwoDMatrix& getVcov() const
///		{return model->getVcov();}
	Vector& getParams()
		{return params;}
///	const Vector& getParams() const
///		{return model->getParams();}
//	void setInitOuter(const Vector& x)
//		{model->setInitOuter(x);}

	const TensorContainer<FSSparseTensor>& getModelDerivatives() const
		{return md;}
	const Vector& getSteady() const
		{return *ysteady;}
	Vector& getSteady()
		{return *ysteady;}
///	const ogdyn::DynareModel& getModel() const
///		{return *model;}

	// here is true public interface
	void solveDeterministicSteady(Vector& steady);
	void solveDeterministicSteady()
		{solveDeterministicSteady(*ysteady);}
	void evaluateSystem(Vector& out, const Vector& yy, const Vector& xx);
	void evaluateSystem(Vector& out, const Vector& yym, const Vector& yy,
						const Vector& yyp, const Vector& xx);
	void calcDerivatives(const Vector& yy, const Vector& xx);
	void calcDerivativesAtSteady();

///	void writeMat4(FILE* fd, const char* prefix) const;
///	void writeDump(const std::string& basename) const;
private:
	void writeModelInfo(Journal& jr) const;
};

/*********************

class DynareJacobian : public ogu::Jacobian ///, public ogp::FormulaDerEvalLoader 
{
protected:
	Dynare& d;
public:
	DynareJacobian(Dynare& dyn);
	virtual ~DynareJacobian() {}
	void load(int i, int iord, const int* vars, double res);
	void eval(const Vector& in);
};
************//////////
/****************************
*  K-Order Perturbation instance of Jacobian:
************************************/

class KordpJacobian : public ogu::Jacobian ///, public ogp::FormulaDerEvalLoader 
{
protected:
	KordpDynare& kdyn;
public:
	KordpJacobian(	KordpDynare& dyn);
	virtual ~KordpJacobian() {}
	// Load <mod>_dynamic.DLL 
	void load(const char** modName);
	void eval(const Vector& in);
};


/****************************
*  K-Order Perturbation instance of VectorFunction:
************************************/

class 	KordpVectorFunction : public ogu::VectorFunction {
protected:
		KordpDynare& d;
public:
	KordpVectorFunction(	KordpDynare& dyn)
		: d(dyn) {}
	virtual ~KordpVectorFunction() {}
	int inDim() const
		{return d.ny();}
	int outDim() const
		{return d.ny();}
	void eval(const ConstVector& in, Vector& out);
};

#endif

// Local Variables:
// mode:C++
// End:
