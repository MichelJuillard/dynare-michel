/*1:*/

#ifndef QUASI_MCARLO_H
#define QUASI_MCARLO_H

#include "int_sequence.h"
#include "quadrature.h"

#include "Vector.h"

#include <vector> 

/*2:*/

class PermutationScheme{
public:
PermutationScheme(){}
virtual~PermutationScheme(){}
virtual int permute(int i,int base,int c)const= 0;
};


/*:2*/
;
/*3:*/

class RadicalInverse{
int num;
int base;
int maxn;
int j;
IntSequence coeff;
public:
RadicalInverse(int n,int b,int mxn);
RadicalInverse(const RadicalInverse&ri)
:num(ri.num),base(ri.base),maxn(ri.maxn),j(ri.j),coeff(ri.coeff){}
const RadicalInverse&operator= (const RadicalInverse&radi)
{
num= radi.num;base= radi.base;maxn= radi.maxn;
j= radi.j;coeff= radi.coeff;
return*this;
}
double eval(const PermutationScheme&p)const;
void increase();
void print()const;
};

/*:3*/
;
/*4:*/

class HaltonSequence{
private:
static int primes[];
static int num_primes;
protected:
int num;
int maxn;
vector<RadicalInverse> ri;
const PermutationScheme&per;
Vector pt;
public:
HaltonSequence(int n,int mxn,int dim,const PermutationScheme&p);
HaltonSequence(const HaltonSequence&hs)
:num(hs.num),maxn(hs.maxn),ri(hs.ri),per(hs.per),pt(hs.pt){}
const HaltonSequence&operator= (const HaltonSequence&hs);
void increase();
const Vector&point()const
{return pt;}
const int getNum()const
{return num;}
void print()const;
protected:
void eval();
};

/*:4*/
;
/*5:*/

class QMCSpecification{
protected:
int dim;
int lev;
const PermutationScheme&per_scheme;
public:
QMCSpecification(int d,int l,const PermutationScheme&p)
:dim(d),lev(l),per_scheme(p){}
virtual~QMCSpecification(){}
int dimen()const
{return dim;}
int level()const
{return lev;}
const PermutationScheme&getPerScheme()const
{return per_scheme;}
};


/*:5*/
;
/*6:*/

class qmcpit{
protected:
const QMCSpecification*spec;
HaltonSequence*halton;
ParameterSignal*sig;
public:
qmcpit();
qmcpit(const QMCSpecification&s,int n);
qmcpit(const qmcpit&qpit);
~qmcpit();
bool operator==(const qmcpit&qpit)const;
bool operator!=(const qmcpit&qpit)const
{return!operator==(qpit);}
const qmcpit&operator= (const qmcpit&qpit);
qmcpit&operator++();
const ParameterSignal&signal()const
{return*sig;}
const Vector&point()const
{return halton->point();}
double weight()const;
void print()const
{halton->print();}
};

/*:6*/
;
/*7:*/

class QMCarloCubeQuadrature:public QuadratureImpl<qmcpit> ,public QMCSpecification{
public:
QMCarloCubeQuadrature(int d,int l,const PermutationScheme&p)
:QuadratureImpl<qmcpit> (d),QMCSpecification(d,l,p){}
virtual~QMCarloCubeQuadrature(){}
int numEvals(int l)const
{return l;}
protected:
qmcpit begin(int ti,int tn,int lev)const
{return qmcpit(*this,ti*level()/tn+1);}
};

/*:7*/
;
/*8:*/

class qmcnpit:public qmcpit{
protected:
Vector*pnt;
public:
qmcnpit();
qmcnpit(const QMCSpecification&spec,int n);
qmcnpit(const qmcnpit&qpit);
~qmcnpit();
bool operator==(const qmcnpit&qpit)const
{return qmcpit::operator==(qpit);}
bool operator!=(const qmcnpit&qpit)const
{return!operator==(qpit);}
const qmcnpit&operator= (const qmcnpit&qpit);
qmcnpit&operator++();
const ParameterSignal&signal()const
{return*sig;}
const Vector&point()const
{return*pnt;}
void print()const
{halton->print();pnt->print();}
};

/*:8*/
;
/*9:*/

class QMCarloNormalQuadrature:public QuadratureImpl<qmcnpit> ,public QMCSpecification{
public:
QMCarloNormalQuadrature(int d,int l,const PermutationScheme&p)
:QuadratureImpl<qmcnpit> (d),QMCSpecification(d,l,p){}
virtual~QMCarloNormalQuadrature(){}
int numEvals(int l)const
{return l;}
protected:
qmcnpit begin(int ti,int tn,int lev)const
{return qmcnpit(*this,ti*level()/tn+1);}
};

/*:9*/
;
/*10:*/

class WarnockPerScheme:public PermutationScheme{
public:
int permute(int i,int base,int c)const;
};

/*:10*/
;
/*11:*/

class ReversePerScheme:public PermutationScheme{
public:
int permute(int i,int base,int c)const;
};

/*:11*/
;
/*12:*/

class IdentityPerScheme:public PermutationScheme{
public:
int permute(int i,int base,int c)const
{return c;}
};

/*:12*/
;

#endif

/*:1*/
