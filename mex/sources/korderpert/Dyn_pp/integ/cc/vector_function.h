/*1:*/

#ifndef VECTOR_FUNCTION_H
#define VECTOR_FUNCTION_H

#include "Vector.h"
#include "GeneralMatrix.h"

#include <vector> 

/*2:*/

class ParameterSignal{
protected:
bool*data;
int num;
public:
ParameterSignal(int n);
ParameterSignal(const ParameterSignal&sig);
~ParameterSignal()
{delete[]data;}
void signalAfter(int l);
const bool&operator[](int i)const
{return data[i];}
bool&operator[](int i)
{return data[i];}
};

/*:2*/
;
/*3:*/

class VectorFunction{
protected:
int in_dim;
int out_dim;
public:
VectorFunction(int idim,int odim)
:in_dim(idim),out_dim(odim){}
VectorFunction(const VectorFunction&func)
:in_dim(func.in_dim),out_dim(func.out_dim){}
virtual~VectorFunction(){}
virtual VectorFunction*clone()const= 0;
virtual void eval(const Vector&point,const ParameterSignal&sig,Vector&out)= 0;
int indim()const
{return in_dim;}
int outdim()const
{return out_dim;}
};

/*:3*/
;
/*4:*/

class VectorFunctionSet{
protected:
std::vector<VectorFunction*> funcs;
bool first_shallow;
public:
VectorFunctionSet(const VectorFunction&f,int n);
VectorFunctionSet(VectorFunction&f,int n);
~VectorFunctionSet();
VectorFunction&getFunc(int i)
{return*(funcs[i]);}
int getNum()const
{return funcs.size();}
};

/*:4*/
;
/*5:*/

class GaussConverterFunction:public VectorFunction{
protected:
VectorFunction*func;
bool delete_flag;
GeneralMatrix A;
double multiplier;
public:
GaussConverterFunction(VectorFunction&f,const GeneralMatrix&vcov);
GaussConverterFunction(VectorFunction*f,const GeneralMatrix&vcov);
GaussConverterFunction(const GaussConverterFunction&f);
virtual~GaussConverterFunction()
{if(delete_flag)delete func;}
virtual VectorFunction*clone()const
{return new GaussConverterFunction(*this);}
virtual void eval(const Vector&point,const ParameterSignal&sig,Vector&out);
private:
double calcMultiplier()const;
void calcCholeskyFactor(const GeneralMatrix&vcov);
};

/*:5*/
;

#endif

/*:1*/
