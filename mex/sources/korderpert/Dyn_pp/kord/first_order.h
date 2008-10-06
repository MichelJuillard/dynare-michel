/*1:*/
#line 10 "./first_order.hweb"


#ifndef FIRST_ORDER_H
#define FIRST_ORDER_H

#include "korder.h"

/*2:*/
#line 23 "./first_order.hweb"

template<int> class FirstOrderDerivs;
class FirstOrder{
template<int> friend class FirstOrderDerivs;
PartitionY ypart;
int nu;
TwoDMatrix gy;
TwoDMatrix gu;
bool bk_cond;
double b_error;
int sdim;
Vector alphar;
Vector alphai;
Vector beta;
Journal&journal;
public:
FirstOrder(int num_stat,int num_pred,int num_both,int num_forw,
int num_u,const FSSparseTensor&f,Journal&jr)
:ypart(num_stat,num_pred,num_both,num_forw),
nu(num_u),
gy(ypart.ny(),ypart.nys()),
gu(ypart.ny(),nu),
alphar(ypart.ny()+ypart.nboth),
alphai(ypart.ny()+ypart.nboth),
beta(ypart.ny()+ypart.nboth),
journal(jr)
{solve(FFSTensor(f));}
bool isStable()const
{return bk_cond;}
const TwoDMatrix&getGy()const
{return gy;}
const TwoDMatrix&getGu()const
{return gu;}
protected:
void solve(const TwoDMatrix&f);
void journalEigs();
};

/*:2*/
#line 17 "./first_order.hweb"
;
/*3:*/
#line 64 "./first_order.hweb"

template<int t> 
class FirstOrderDerivs:public ctraits<t> ::Tg{
public:
FirstOrderDerivs(const FirstOrder&fo)
:ctraits<t> ::Tg(4)
{
IntSequence nvs(4);
nvs[0]= fo.ypart.nys();nvs[1]= fo.nu;nvs[2]= fo.nu;nvs[3]= 1;
_Ttensor*ten= new _Ttensor(fo.ypart.ny(),TensorDimens(Symmetry(1,0,0,0),nvs));
ten->zeros();ten->add(1.0,fo.gy);
insert(ten);
ten= new _Ttensor(fo.ypart.ny(),TensorDimens(Symmetry(0,1,0,0),nvs));
ten->zeros();ten->add(1.0,fo.gu);
insert(ten);
}
};


/*:3*/
#line 18 "./first_order.hweb"
;

#endif

/*:1*/
