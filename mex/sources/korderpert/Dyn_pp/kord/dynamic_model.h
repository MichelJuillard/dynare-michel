/*1:*/
#line 13 "./dynamic_model.hweb"

#ifndef DYNAMIC_MODEL_H
#define DYNAMIC_MODEL_H

#include "t_container.h"
#include "sparse_tensor.h"

#include "Vector.h"

/*2:*/
#line 29 "./dynamic_model.hweb"

class NameList{
public:
virtual~NameList(){}
virtual int getNum()const= 0;
virtual const char*getName(int i)const= 0;
void print()const;
void writeMat4(FILE*fd,const char*vname)const;
void writeMat4Indices(FILE*fd,const char*prefix)const;
};

/*:2*/
#line 22 "./dynamic_model.hweb"
;
/*3:*/
#line 89 "./dynamic_model.hweb"

class DynamicModel{
public:
virtual DynamicModel*clone()const= 0;
virtual~DynamicModel(){}

virtual int nstat()const= 0;
virtual int nboth()const= 0;
virtual int npred()const= 0;
virtual int nforw()const= 0;
virtual int nexog()const= 0;
virtual int order()const= 0;
int numeq()const
{return nstat()+nboth()+npred()+nforw();}

virtual const NameList&getAllEndoNames()const= 0;
virtual const NameList&getStateNames()const= 0;
virtual const NameList&getExogNames()const= 0;
virtual const TwoDMatrix&getVcov()const= 0;
virtual const TensorContainer<FSSparseTensor> &getModelDerivatives()const= 0;
virtual const Vector&getSteady()const= 0;
virtual Vector&getSteady()= 0;

virtual void solveDeterministicSteady()= 0;
virtual void evaluateSystem(Vector&out,const Vector&yy,const Vector&xx)= 0;
virtual void evaluateSystem(Vector&out,const Vector&yym,const Vector&yy,
const Vector&yyp,const Vector&xx)= 0;
virtual void calcDerivativesAtSteady()= 0;
};


/*:3*/
#line 23 "./dynamic_model.hweb"
;

#endif

/*:1*/
