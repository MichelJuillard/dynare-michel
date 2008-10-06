/*1:*/
#line 69 "./pyramid_prod2.hweb"

#ifndef PYRAMID_PROD2_H
#define PYRAMID_PROD2_H

#include "permutation.h"
#include "tensor.h"
#include "tl_exception.h"
#include "rfs_tensor.h"
#include "stack_container.h"

#include "Vector.h"

/*2:*/
#line 101 "./pyramid_prod2.hweb"

class IrregTensor;
class IrregTensorHeader{
friend class IrregTensor;
int nv;
IntSequence unit_flag;
Vector**const cols;
IntSequence end_seq;
public:
IrregTensorHeader(const StackProduct<FGSTensor> &sp,const IntSequence&c);
~IrregTensorHeader();
int dimen()const
{return unit_flag.size();}
void increment(IntSequence&v)const;
int calcMaxOffset()const;
private:
IrregTensorHeader(const IrregTensorHeader&);
};


/*:2*/
#line 81 "./pyramid_prod2.hweb"
;
/*3:*/
#line 137 "./pyramid_prod2.hweb"

class IrregTensor:public Tensor{
const IrregTensorHeader&header;
public:
IrregTensor(const IrregTensorHeader&h);
void addTo(FRSingleTensor&out)const;
void increment(IntSequence&v)const
{header.increment(v);}
void decrement(IntSequence&v)const
{TL_RAISE("Not implemented error in IrregTensor::decrement");}
int getOffset(const IntSequence&v)const
{TL_RAISE("Not implemented error in IrregTensor::getOffset");return 0;}
};

/*:3*/
#line 82 "./pyramid_prod2.hweb"
;

#endif

/*:1*/
