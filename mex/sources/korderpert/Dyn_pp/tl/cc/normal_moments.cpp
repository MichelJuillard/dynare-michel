/*1:*/

#include "normal_moments.h"
#include "permutation.h"
#include "kron_prod.h"
#include "tl_static.h"

/*2:*/

UNormalMoments::UNormalMoments(int maxdim,const TwoDMatrix&v)
:TensorContainer<URSingleTensor> (1)
{
	if(maxdim>=2)
		generateMoments(maxdim,v);
}


/*:2*/
;
/*3:*/

void UNormalMoments::generateMoments(int maxdim,const TwoDMatrix&v)
{
	TL_RAISE_IF(v.nrows()!=v.ncols(),
		"Variance-covariance matrix is not square in UNormalMoments constructor");
	
	int nv= v.nrows();
	URSingleTensor*mom2= new URSingleTensor(nv,2);
	mom2->getData()= v.getData();
	insert(mom2);
	URSingleTensor*kronv= new URSingleTensor(nv,2);
	kronv->getData()= v.getData();
	for(int d= 4;d<=maxdim;d+= 2){
		URSingleTensor*newkronv= new URSingleTensor(nv,d);
		KronProd::kronMult(ConstVector(v.getData()),
			ConstVector(kronv->getData()),
			newkronv->getData());
		delete kronv;
		kronv= newkronv;
		URSingleTensor*mom= new URSingleTensor(nv,d);
		/*4:*/
		
		mom->zeros();
		const EquivalenceSet eset= ebundle.get(d);
		for(EquivalenceSet::const_iterator cit= eset.begin();
		cit!=eset.end();cit++){
			if(selectEquiv(*cit)){
				Permutation per(*cit);
				per.inverse();
				for(Tensor::index it= kronv->begin();it!=kronv->end();++it){
					IntSequence ind(kronv->dimen());
					per.apply(it.getCoor(),ind);
					Tensor::index it2(mom,ind);
					mom->get(*it2,0)+= kronv->get(*it,0);
				}
			}
		}
		
		/*:4*/
		;
		insert(mom);
	}
	delete kronv;
}

/*:3*/
;
/*5:*/

bool UNormalMoments::selectEquiv(const Equivalence&e)
{
	if(2*e.numClasses()!=e.getN())
		return false;
	for(Equivalence::const_seqit si= e.begin();
	si!=e.end();++si){
		if((*si).length()!=2)
			return false;
	}
	return true;
}

/*:5*/
;
/*6:*/

FNormalMoments::FNormalMoments(const UNormalMoments&moms)
:TensorContainer<FRSingleTensor> (1)
{
	for(UNormalMoments::const_iterator it= moms.begin();
	it!=moms.end();++it){
		FRSingleTensor*fm= new FRSingleTensor(*((*it).second));
		insert(fm);
	}
}


/*:6*/
;

/*:1*/
