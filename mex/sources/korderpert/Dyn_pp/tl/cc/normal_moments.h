/*1:*/
#line 110 "./normal_moments.hweb"

#ifndef NORMAL_MOMENTS_H
#define NORMAL_MOMENTS_H

#include "t_container.h"

/*2:*/
#line 122 "./normal_moments.hweb"

class UNormalMoments:public TensorContainer<URSingleTensor> {
public:
UNormalMoments(int maxdim,const TwoDMatrix&v);
private:
void generateMoments(int maxdim,const TwoDMatrix&v);
static bool selectEquiv(const Equivalence&e);
};

/*:2*/
#line 116 "./normal_moments.hweb"
;
/*3:*/
#line 132 "./normal_moments.hweb"

class FNormalMoments:public TensorContainer<FRSingleTensor> {
public:
FNormalMoments(const UNormalMoments&moms);
};


/*:3*/
#line 117 "./normal_moments.hweb"
;

#endif

/*:1*/
