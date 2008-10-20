/*1:*/

#ifndef NORMAL_MOMENTS_H
#define NORMAL_MOMENTS_H

#include "t_container.h"

/*2:*/

class UNormalMoments:public TensorContainer<URSingleTensor> {
public:
	UNormalMoments(int maxdim,const TwoDMatrix&v);
private:
	void generateMoments(int maxdim,const TwoDMatrix&v);
	static bool selectEquiv(const Equivalence&e);
};

/*:2*/
;
/*3:*/

class FNormalMoments:public TensorContainer<FRSingleTensor> {
public:
	FNormalMoments(const UNormalMoments&moms);
};


/*:3*/
;

#endif

/*:1*/
