/*1:*/

#include "t_container.h"
#include "kron_prod.h"
#include "ps_tensor.h"
#include "pyramid_prod.h"

const int FGSContainer::num_one_time= 10;
/*2:*/

UGSContainer::UGSContainer(const FGSContainer&c)
:TensorContainer<UGSTensor> (c.num())
{
	for(FGSContainer::const_iterator it= c.begin();
	it!=c.end();++it){
		UGSTensor*unfolded= new UGSTensor(*((*it).second));
		insert(unfolded);
	}
}

/*:2*/
;
/*3:*/

void UGSContainer::multAndAdd(const UGSTensor&t,UGSTensor&out)const
{
	int l= t.dimen();
	int k= out.dimen();
	const EquivalenceSet&eset= ebundle.get(k);
	
	for(EquivalenceSet::const_iterator it= eset.begin();
	it!=eset.end();++it){
		if((*it).numClasses()==l){
			vector<const UGSTensor*> ts= 
				fetchTensors(out.getSym(),*it);
			KronProdAllOptim kp(l);
			for(int i= 0;i<l;i++)
				kp.setMat(i,*(ts[i]));
			kp.optimizeOrder();
			UPSTensor ups(out.getDims(),*it,t,kp);
			ups.addTo(out);
		}
	}
}

/*:3*/
;
/*4:*/

FGSContainer::FGSContainer(const UGSContainer&c)
:TensorContainer<FGSTensor> (c.num())
{
	for(UGSContainer::const_iterator it= c.begin();
	it!=c.end();++it){
		FGSTensor*folded= new FGSTensor(*((*it).second));
		insert(folded);
	}
}


/*:4*/
;
/*5:*/

void FGSContainer::multAndAdd(const FGSTensor&t,FGSTensor&out)const
{
	UGSTensor ut(t);
	multAndAdd(ut,out);
}

/*:5*/
;
/*6:*/

void FGSContainer::multAndAdd(const UGSTensor&t,FGSTensor&out)const
{
	int l= t.dimen();
	int k= out.dimen();
	const EquivalenceSet&eset= ebundle.get(k);
	
	for(EquivalenceSet::const_iterator it= eset.begin();
	it!=eset.end();++it){
		if((*it).numClasses()==l){
			vector<const FGSTensor*> ts= 
				fetchTensors(out.getSym(),*it);
			KronProdAllOptim kp(l);
			for(int i= 0;i<l;i++)
				kp.setMat(i,*(ts[i]));
			kp.optimizeOrder();
			FPSTensor fps(out.getDims(),*it,t,kp);
			fps.addTo(out);
		}
	}
}


/*:6*/
;
/*7:*/

Tensor::index
FGSContainer::getIndices(int num,vector<IntSequence> &out,
						 const Tensor::index&start,
						 const Tensor::index&end)
{
	out.clear();
	int i= 0;
	Tensor::index run= start;
	while(i<num&&run!=end){
		out.push_back(run.getCoor());
		i++;
		++run;
	}
	return run;
}


/*:7*/
;

/*:1*/
