/*1:*/

#include "t_polynomial.h"
#include "kron_prod.h"

/*2:*/

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
;
/*3:*/

const FRSingleTensor&PowerProvider::getNext(const FRSingleTensor*dummy)
{
	getNext(ut);
	if(ft)
		delete ft;
	ft= new FRSingleTensor(*ut);
	return*ft;
}

/*:3*/
;
/*4:*/

PowerProvider::~PowerProvider()
{
	if(ut)
		delete ut;
	if(ft)
		delete ft;
}

/*:4*/
;
/*5:*/

UTensorPolynomial::UTensorPolynomial(const FTensorPolynomial&fp)
:TensorPolynomial<UFSTensor,UGSTensor,URSingleTensor> (fp.nrows(),fp.nvars())
{
	for(FTensorPolynomial::const_iterator it= fp.begin();
	it!=fp.end();++it){
		insert(new UFSTensor(*((*it).second)));
	}
}

/*:5*/
;
/*6:*/

FTensorPolynomial::FTensorPolynomial(const UTensorPolynomial&up)
:TensorPolynomial<FFSTensor,FGSTensor,FRSingleTensor> (up.nrows(),up.nvars())
{
	for(UTensorPolynomial::const_iterator it= up.begin();
	it!=up.end();++it){
		insert(new FFSTensor(*((*it).second)));
	}
}

/*:6*/
;


/*:1*/
