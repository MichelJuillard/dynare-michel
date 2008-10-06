/*1:*/
#line 37 "./normal_conjugate.hweb"

#ifndef NORMAL_CONJUGATE_H
#define NORMAL_CONJUGATE_H

#include "twod_matrix.h"

/*2:*/
#line 49 "./normal_conjugate.hweb"

class NormalConj{
protected:
Vector mu;
int kappa;
int nu;
TwoDMatrix lambda;
public:
/*3:*/
#line 76 "./normal_conjugate.hweb"

NormalConj(int d);
NormalConj(const ConstTwoDMatrix&ydata);
NormalConj(const NormalConj&nc);


/*:3*/
#line 57 "./normal_conjugate.hweb"
;
virtual~NormalConj(){}
void update(const ConstVector&y);
void update(const ConstTwoDMatrix&ydata);
void update(const NormalConj&nc);
int getDim()const
{return mu.length();}
const Vector&getMean()const
{return mu;}
void getVariance(TwoDMatrix&v)const;
};

/*:2*/
#line 43 "./normal_conjugate.hweb"
;

#endif

/*:1*/
