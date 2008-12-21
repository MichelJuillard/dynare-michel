/*1:*/

#ifndef PRODUCT_H
#define PRODUCT_H

#include "int_sequence.h"
#include "vector_function.h"
#include "quadrature.h"

/*2:*/

class ProductQuadrature;

class prodpit{
protected:
const ProductQuadrature*prodq;
int level;
int npoints;
IntSequence*jseq;
bool end_flag;
ParameterSignal*sig;
Vector*p;
double w;
public:
prodpit();
prodpit(const ProductQuadrature&q,int j0,int l);
prodpit(const prodpit&ppit);
~prodpit();
bool operator==(const prodpit&ppit)const;
bool operator!=(const prodpit&ppit)const
{return!operator==(ppit);}
const prodpit&operator= (const prodpit&spit);
prodpit&operator++();
const ParameterSignal&signal()const
{return*sig;}
const Vector&point()const
{return*p;}
double weight()const
{return w;}
void print()const;
protected:
void setPointAndWeight();
};

/*:2*/
;
/*3:*/

class ProductQuadrature:public QuadratureImpl<prodpit> {
friend class prodpit;
const OneDQuadrature&uquad;
public:
ProductQuadrature(int d,const OneDQuadrature&uq);
virtual~ProductQuadrature(){}
int numEvals(int l)const
{
int res= 1;
for(int i= 0;i<dimen();i++)
res*= uquad.numPoints(l);
return res;
}
void designLevelForEvals(int max_eval,int&lev,int&evals)const;
protected:
prodpit begin(int ti,int tn,int level)const;
};

/*:3*/
;

#endif

/*:1*/
