/*1:*/

#ifndef NORMAL_CONJUGATE_H
#define NORMAL_CONJUGATE_H

#include "twod_matrix.h"

/*2:*/

class NormalConj{
protected:
	Vector mu;
	int kappa;
	int nu;
	TwoDMatrix lambda;
public:
	/*3:*/
	
	NormalConj(int d);
	NormalConj(const ConstTwoDMatrix&ydata);
	NormalConj(const NormalConj&nc);
	
	
	/*:3*/
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
;

#endif

/*:1*/
