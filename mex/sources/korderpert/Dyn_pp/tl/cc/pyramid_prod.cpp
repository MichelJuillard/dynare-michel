/*1:*/
#line 5 "./pyramid_prod.cweb"


#include "pyramid_prod.h"
#include "permutation.h"
#include "tl_exception.h"

/*2:*/
#line 27 "./pyramid_prod.cweb"

USubTensor::USubTensor(const TensorDimens&bdims,
const TensorDimens&hdims,
const FGSContainer&cont,
const vector<IntSequence> &lst)
:URTensor(lst.size(),hdims.getNVX()[0],hdims.dimen())
{
TL_RAISE_IF(!hdims.getNVX().isConstant(),
"Tensor has not full symmetry in USubTensor()");
const EquivalenceSet&eset= cont.getEqBundle().get(bdims.dimen());
zeros();
for(EquivalenceSet::const_iterator it= eset.begin();
it!=eset.end();++it){
if((*it).numClasses()==hdims.dimen()){
Permutation per(*it);
vector<const FGSTensor*> ts= 
cont.fetchTensors(bdims.getSym(),*it);
for(int i= 0;i<(int)lst.size();i++){
IntSequence perindex(lst[i].size());
per.apply(lst[i],perindex);
addKronColumn(i,ts,perindex);
}
}
}
}

/*:2*/
#line 11 "./pyramid_prod.cweb"
;
/*3:*/
#line 67 "./pyramid_prod.cweb"

void USubTensor::addKronColumn(int i,const vector<const FGSTensor*> &ts,
const IntSequence&pindex)
{
vector<ConstVector> tmpcols;
int lastdim= 0;
for(unsigned int j= 0;j<ts.size();j++){
IntSequence ind(pindex,lastdim,lastdim+ts[j]->dimen());
lastdim+= ts[j]->dimen();
index in(ts[j],ind);
tmpcols.push_back(ConstVector(*(ts[j]),*in));
}

URSingleTensor kronmult(tmpcols);
Vector coli(*this,i);
coli.add(1.0,kronmult.getData());
}


/*:3*/
#line 12 "./pyramid_prod.cweb"
;


/*:1*/
