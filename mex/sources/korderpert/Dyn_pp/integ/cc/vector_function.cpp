/*1:*/


#include "vector_function.h"
#include "cpplapack.h"

#include <math.h> 

#include <string.h> 
#include <algorithm> 

/*2:*/

ParameterSignal::ParameterSignal(int n)
:data(new bool[n]),num(n)
{
for(int i= 0;i<num;i++)
data[i]= true;
}

/*:2*/
;
/*3:*/

ParameterSignal::ParameterSignal(const ParameterSignal&sig)
:data(new bool[sig.num]),num(sig.num)
{
memcpy(data,sig.data,num);
}

/*:3*/
;
/*4:*/

void ParameterSignal::signalAfter(int l)
{
for(int i= 0;i<std::min(l,num);i++)
data[i]= false;
for(int i= l;i<num;i++)
data[i]= true;
}

/*:4*/
;
/*5:*/

VectorFunctionSet::VectorFunctionSet(const VectorFunction&f,int n)
:funcs(n),first_shallow(false)
{
for(int i= 0;i<n;i++)
funcs[i]= f.clone();
}

/*:5*/
;
/*6:*/

VectorFunctionSet::VectorFunctionSet(VectorFunction&f,int n)
:funcs(n),first_shallow(true)
{
if(n> 0)
funcs[0]= &f;
for(int i= 1;i<n;i++)
funcs[i]= f.clone();
}

/*:6*/
;
/*7:*/

VectorFunctionSet::~VectorFunctionSet()
{
unsigned int start= first_shallow?1:0;
for(unsigned int i= start;i<funcs.size();i++)
delete funcs[i];
}

/*:7*/
;
/*8:*/

GaussConverterFunction::GaussConverterFunction(VectorFunction&f,const GeneralMatrix&vcov)
:VectorFunction(f),func(&f),delete_flag(false),A(vcov.numRows(),vcov.numRows()),
multiplier(calcMultiplier())
{

calcCholeskyFactor(vcov);
}

/*:8*/
;
/*9:*/

GaussConverterFunction::GaussConverterFunction(VectorFunction*f,const GeneralMatrix&vcov)
:VectorFunction(*f),func(f),delete_flag(true),A(vcov.numRows(),vcov.numRows()),
multiplier(calcMultiplier())
{

calcCholeskyFactor(vcov);
}


/*:9*/
;
/*10:*/

GaussConverterFunction::GaussConverterFunction(const GaussConverterFunction&f)
:VectorFunction(f),func(f.func->clone()),delete_flag(true),A(f.A),
multiplier(f.multiplier)
{
}

/*:10*/
;
/*11:*/

void GaussConverterFunction::eval(const Vector&point,const ParameterSignal&sig,Vector&out)
{
ParameterSignal s(sig);
int i= 0;
while(i<indim()&&!sig[i])
i++;
s.signalAfter(i);

Vector x(indim());
x.zeros();
A.multaVec(x,point);
x.mult(sqrt(2.0));

func->eval(x,s,out);

out.mult(multiplier);
}

/*:11*/
;
/*12:*/

double GaussConverterFunction::calcMultiplier()const
{
return sqrt(pow(M_PI,-1*indim()));
}

/*:12*/
;
/*13:*/

void GaussConverterFunction::calcCholeskyFactor(const GeneralMatrix&vcov)
{
A= vcov;

int rows= A.numRows();
for(int i= 0;i<rows;i++)
for(int j= i+1;j<rows;j++)
A.get(i,j)= 0.0;

int info;
LAPACK_dpotrf("L",&rows,A.base(),&rows,&info);

}


/*:13*/
;

/*:1*/
