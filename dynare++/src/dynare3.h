// $Id: dynare3.h 1764 2008-03-31 14:30:55Z kamenik $
// Copyright 2005, Ondra Kamenik

#ifndef DYNARE3_H
#define DYNARE3_H

#include "../tl/cc/t_container.h"
#include "../tl/cc/sparse_tensor.h"
#include "../kord/decision_rule.h"
#include "../kord/dynamic_model.h"

#include "dynare_model.h"
#include "nlsolve.h"

#include <vector>

#include <matio.h>

class Dynare;

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
	 * the name cannot be found, then an exception is raised. */
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

// The following only implements DynamicModel with help of ogdyn::DynareModel

class DynareJacobian;
class Dynare : public DynamicModel {
	friend class DynareNameList;
	friend class DynareExogNameList;
	friend class DynareStateNameList;
	friend class DynareJacobian;
	Journal& journal;
	ogdyn::DynareModel* model;
	Vector* ysteady;
	TensorContainer<FSSparseTensor> md;
	DynareNameList* dnl;
	DynareExogNameList* denl;
	DynareStateNameList* dsnl;
	ogp::FormulaEvaluator* fe;
	ogp::FormulaDerEvaluator* fde;
	const double ss_tol;
public:
	/** Parses the given model file and uses the given order to
	 * override order from the model file (if it is != -1). */
	Dynare(const char* modname, int ord, double sstol, Journal& jr);
	/** Parses the given equations with explicitly given names. */
	Dynare(const char** endo, int num_endo,
		   const char** exo, int num_exo,
		   const char** par, int num_par,
		   const char* equations, int len, int ord,
		   double sstol, Journal& jr);
	/** Makes a deep copy of the object. */
	Dynare(const Dynare& dyn);
	DynamicModel* clone() const
		{return new Dynare(*this);}
	virtual ~Dynare();
	int nstat() const
		{return model->getAtoms().nstat();}
	int nboth() const
		{return model->getAtoms().nboth();}
	int npred() const
		{return model->getAtoms().npred();}
	int nforw() const
		{return model->getAtoms().nforw();}
	int nexog() const
		{return model->getAtoms().nexo();}
	int nys() const
		{return model->getAtoms().nys();}
	int nyss() const
		{return model->getAtoms().nyss();}
	int ny() const
		{return model->getAtoms().ny();}
	int order() const
		{return model->getOrder();}

	const NameList& getAllEndoNames() const
		{return *dnl;}
	const NameList& getStateNames() const
		{return *dsnl;}
	const NameList& getExogNames() const
		{return *denl;}

	TwoDMatrix& getVcov()
		{return model->getVcov();}
	const TwoDMatrix& getVcov() const
		{return model->getVcov();}
	Vector& getParams()
		{return model->getParams();}
	const Vector& getParams() const
		{return model->getParams();}
	void setInitOuter(const Vector& x)
		{model->setInitOuter(x);}

	const TensorContainer<FSSparseTensor>& getModelDerivatives() const
		{return md;}
	const Vector& getSteady() const
		{return *ysteady;}
	Vector& getSteady()
		{return *ysteady;}
	const ogdyn::DynareModel& getModel() const
		{return *model;}

	// here is true public interface
	void solveDeterministicSteady(Vector& steady);
	void solveDeterministicSteady()
		{solveDeterministicSteady(*ysteady);}
	void evaluateSystem(Vector& out, const Vector& yy, const Vector& xx);
	void evaluateSystem(Vector& out, const Vector& yym, const Vector& yy,
						const Vector& yyp, const Vector& xx);
	void calcDerivatives(const Vector& yy, const Vector& xx);
	void calcDerivativesAtSteady();

	void writeMat(mat_t *fd, const char* prefix) const;
	void writeDump(const std::string& basename) const;
private:
	void writeModelInfo(Journal& jr) const;
};

class DynareEvalLoader : public ogp::FormulaEvalLoader, public Vector {
public:
	DynareEvalLoader(const ogp::FineAtoms& a, Vector& out);
	void load(int i, double res)
		{operator[](i) = res;}
};

class DynareDerEvalLoader : public ogp::FormulaDerEvalLoader {
protected:
	const ogp::FineAtoms& atoms;
	TensorContainer<FSSparseTensor>& md;
public:
	DynareDerEvalLoader(const ogp::FineAtoms& a, TensorContainer<FSSparseTensor>& mod_ders,
						int order);
	void load(int i, int iord, const int* vars, double res);
};

class DynareJacobian : public ogu::Jacobian, public ogp::FormulaDerEvalLoader {
protected:
	Dynare& d;
public:
	DynareJacobian(Dynare& dyn);
	virtual ~DynareJacobian() {}
	void load(int i, int iord, const int* vars, double res);
	void eval(const Vector& in);
};

class DynareVectorFunction : public ogu::VectorFunction {
protected:
	Dynare& d;
public:
	DynareVectorFunction(Dynare& dyn)
		: d(dyn) {}
	virtual ~DynareVectorFunction() {}
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
