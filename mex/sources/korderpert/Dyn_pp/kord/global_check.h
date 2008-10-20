/*1:*/

#ifndef GLOBAL_CHECK_H
#define GLOBAL_CHECK_H

#include "vector_function.h"
#include "quadrature.h"

#include "dynamic_model.h"
#include "journal.h"
#include "approximation.h"

/*2:*/

class ResidFunction:public VectorFunction{
protected:
	const Approximation&approx;
	DynamicModel*model;
	Vector*yplus;
	Vector*ystar;
	Vector*u;
	FTensorPolynomial*hss;
public:
	ResidFunction(const Approximation&app);
	ResidFunction(const ResidFunction&rf);
	virtual~ResidFunction();
	virtual VectorFunction*clone()const
	{return new ResidFunction(*this);}
	virtual void eval(const Vector&point,const ParameterSignal&sig,Vector&out);
	void setYU(const Vector&ys,const Vector&xx);
};

/*:2*/
;
/*3:*/

class GResidFunction:public GaussConverterFunction{
public:
	GResidFunction(const Approximation&app)
		:GaussConverterFunction(new ResidFunction(app),app.getModel().getVcov()){}
	GResidFunction(const GResidFunction&rf)
		:GaussConverterFunction(rf){}
	virtual~GResidFunction(){}
	virtual VectorFunction*clone()const
	{return new GResidFunction(*this);}
	void setYU(const Vector&ys,const Vector&xx)
	{((ResidFunction*)func)->setYU(ys,xx);}
};


/*:3*/
;
/*4:*/

class GlobalChecker{
	const Approximation&approx;
	const DynamicModel&model;
	Journal&journal;
	GResidFunction rf;
	VectorFunctionSet vfs;
public:
	GlobalChecker(const Approximation&app,int n,Journal&jr)
		:approx(app),model(approx.getModel()),journal(jr),
		rf(approx),vfs(rf,n){}
	void check(int max_evals,const ConstTwoDMatrix&y,
		const ConstTwoDMatrix&x,TwoDMatrix&out);
	void checkAlongShocksAndSave(FILE*fd,const char*prefix,
		int m,double mult,int max_evals);
	void checkOnEllipseAndSave(FILE*fd,const char*prefix,
		int m,double mult,int max_evals);
	void checkAlongSimulationAndSave(FILE*fd,const char*prefix,
		int m,int max_evals);
	void checkUnconditionalAndSave(FILE*fd,const char*prefix,
		int m,int max_evals);
protected:
	void check(const Quadrature&quad,int level,
		const ConstVector&y,const ConstVector&x,Vector&out);
};


/*:4*/
;
/*5:*/

class ResidFunctionSig:public ResidFunction{
public:
	ResidFunctionSig(const Approximation&app,const Vector&ys,const Vector&xx);
};

/*:5*/
;

#endif

/*:1*/
