/*1:*/

#ifndef SMOLYAK_H
#define SMOLYAK_H

#include "int_sequence.h"
#include "tl_static.h"
#include "vector_function.h"
#include "quadrature.h"

/*2:*/

class SmolyakQuadrature;

class smolpit{
protected:
const SmolyakQuadrature*smolq;
unsigned int isummand;
IntSequence*jseq;
ParameterSignal*sig;
Vector*p;
double w;
public:
smolpit();
smolpit(const SmolyakQuadrature&q,unsigned int isum);
smolpit(const smolpit&spit);
~smolpit();
bool operator==(const smolpit&spit)const;
bool operator!=(const smolpit&spit)const
{return!operator==(spit);}
const smolpit&operator= (const smolpit&spit);
smolpit&operator++();
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

class SmolyakQuadrature:public QuadratureImpl<smolpit> {
friend class smolpit;
int level;
const OneDQuadrature&uquad;
vector<IntSequence> levels;
vector<IntSequence> levpoints;
vector<int> cumevals;
PascalTriangle psc;
public:
SmolyakQuadrature(int d,int l,const OneDQuadrature&uq);
virtual~SmolyakQuadrature(){}
virtual int numEvals(int level)const;
void designLevelForEvals(int max_eval,int&lev,int&evals)const;
protected:
smolpit begin(int ti,int tn,int level)const;
unsigned int numSummands()const
{return levels.size();}
private:
int calcNumEvaluations(int level)const;
};

/*:3*/
;

#endif

/*:1*/
