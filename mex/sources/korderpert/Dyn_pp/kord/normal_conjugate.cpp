/*1:*/


#include "normal_conjugate.h"
#include "kord_exception.h"

/*2:*/

NormalConj::NormalConj(int d)
:mu(d),kappa(0),nu(-1),lambda(d,d)
{
	mu.zeros();
	lambda.zeros();
}

/*:2*/
;
/*3:*/

NormalConj::NormalConj(const ConstTwoDMatrix&ydata)
:mu(ydata.numRows()),kappa(ydata.numCols()),nu(ydata.numCols()-1),
lambda(ydata.numRows(),ydata.numRows())
{
	mu.zeros();
	for(int i= 0;i<ydata.numCols();i++)
		mu.add(1.0/ydata.numCols(),ConstVector(ydata,i));
	
	lambda.zeros();
	for(int i= 0;i<ydata.numCols();i++){
		Vector diff(ConstVector(ydata,i));
		diff.add(-1,mu);
		lambda.addOuter(diff);
	}
}

/*:3*/
;
/*4:*/

NormalConj::NormalConj(const NormalConj&nc)
:mu(nc.mu),kappa(nc.kappa),nu(nc.nu),lambda(nc.lambda)
{
}

/*:4*/
;
/*5:*/

void NormalConj::update(const ConstVector&y)
{
	KORD_RAISE_IF(y.length()!=mu.length(),
		"Wrong length of a vector in NormalConj::update");
	
	mu.mult(kappa/(1.0+kappa));
	mu.add(1.0/(1.0+kappa),y);
	
	Vector diff(y);
	diff.add(-1,mu);
	lambda.addOuter(diff,kappa/(1.0+kappa));
	
	kappa++;
	nu++;
}

/*:5*/
;
/*6:*/

void NormalConj::update(const ConstTwoDMatrix&ydata)
{
	NormalConj nc(ydata);
	update(nc);
}


/*:6*/
;
/*7:*/

void NormalConj::update(const NormalConj&nc)
{
	double wold= ((double)kappa)/(kappa+nc.kappa);
	double wnew= 1-wold;
	
	mu.mult(wold);
	mu.add(wnew,nc.mu);
	
	Vector diff(nc.mu);
	diff.add(-1,mu);
	lambda.add(1.0,nc.lambda);
	lambda.addOuter(diff);
	
	kappa= kappa+nc.kappa;
	nu= nu+nc.kappa;
}


/*:7*/
;
/*8:*/

void NormalConj::getVariance(TwoDMatrix&v)const
{
	if(nu> getDim()+1){
		v= (const TwoDMatrix&)lambda;
		v.mult(1.0/(nu-getDim()-1));
	}else
		v.nans();
}


/*:8*/
;

/*:1*/
