/*1:*/
#line 13 "./faa_di_bruno.hweb"

#ifndef FAA_DI_BRUNO_H
#define FAA_DI_BRUNO_H

#include "journal.h"
#include "stack_container.h"
#include "t_container.h"
#include "sparse_tensor.h"
#include "gs_tensor.h"

/*2:*/
#line 30 "./faa_di_bruno.hweb"

class FaaDiBruno{
Journal&journal;
public:
FaaDiBruno(Journal&jr)
:journal(jr){}
void calculate(const StackContainer<FGSTensor> &cont,const TensorContainer<FSSparseTensor> &f,
FGSTensor&out);
void calculate(const FoldedStackContainer&cont,const FGSContainer&g,
FGSTensor&out);
void calculate(const StackContainer<UGSTensor> &cont,const TensorContainer<FSSparseTensor> &f,
UGSTensor&out);
void calculate(const UnfoldedStackContainer&cont,const UGSContainer&g,
UGSTensor&out);
protected:
int estimRefinment(const TensorDimens&tdims,int nr,int l,int&avmem_mb,int&tmpmem_mb);
static double magic_mult;
};

/*:2*/
#line 23 "./faa_di_bruno.hweb"
;

#endif

/*:1*/
