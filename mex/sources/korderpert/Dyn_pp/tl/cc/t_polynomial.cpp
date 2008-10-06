/*1:*/
#line 6 "./t_polynomial.cweb"

#include "t_polynomial.h"
#include "kron_prod.h"

/*2:*/
#line 20 "./t_polynomial.cweb"

const URSingleTensor&PowerProvider::getNext(const URSingleTensor*dummy)
{
if(ut){
URSingleTensor*ut_new= new URSingleTensor(nv,ut->dimen()+1);
KronProd::kronMult(ConstVector(origv),ConstVector(ut->getData()),ut_new->getData());
delete ut;
ut= ut_new;
}else{
ut= new URSingleTensor(nv,1);
ut->getData()= origv;
}
return*ut;
}

/*:2*/
#line 10 "./t_polynomial.cweb"
;
/*3:*/
#line 38 "./t_polynomial.cweb"

const FRSingleTensor&PowerProvider::getNext(const FRSingleTensor*dummy)
{
getNext(ut);
if(ft)
delete ft;
ft= new FRSingleTensor(*ut);
return*ft;
}

/*:3*/
#line 11 "./t_polynomial.cweb"
;
/*4:*/
#line 49 "./t_polynomial.cweb"

PowerProvider::~PowerProvider()
{
if(ut)
delete ut;
if(ft)
delete ft;
}

/*:4*/
#line 12 "./t_polynomial.cweb"
;
/*5:*/
#line 59 "./t_polynomial.cweb"

UTensorPolynomial::UTensorPolynomial(const FTensorPolynomial&fp)
:TensorPolynomial<UFSTensor,UGSTensor,URSingleTensor> (fp.nrows(),fp.nvars())
{
for(FTensorPolynomial::const_iterator it= fp.begin();
it!=fp.end();++it){
insert(new UFSTensor(*((*it).second)));
}
}

/*:5*/
#line 13 "./t_polynomial.cweb"
;
/*6:*/
#line 70 "./t_polynomial.cweb"

FTensorPolynomial::FTensorPolynomial(const UTensorPolynomial&up)
:TensorPolynomial<FFSTensor,FGSTensor,FRSingleTensor> (up.nrows(),up.nvars())
{
for(UTensorPolynomial::const_iterator it= up.begin();
it!=up.end();++it){
insert(new FFSTensor(*((*it).second)));
}
}

/*:6*/
#line 14 "./t_polynomial.cweb"
;


/*:1*/
