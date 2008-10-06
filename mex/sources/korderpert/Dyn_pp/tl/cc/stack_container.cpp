/*1:*/
#line 6 "./stack_container.cweb"

#include "stack_container.h"
#include "pyramid_prod2.h"
#include "ps_tensor.h"

double FoldedStackContainer::fill_threshold= 0.00005;
double UnfoldedStackContainer::fill_threshold= 0.00005;
/*2:*/
#line 59 "./stack_container.cweb"

void FoldedStackContainer::multAndAdd(const FSSparseTensor&t,
FGSTensor&out)const
{
TL_RAISE_IF(t.nvar()!=getAllSize(),
"Wrong number of variables of tensor for FoldedStackContainer::multAndAdd");
multAndAddSparse2(t,out);
}

/*:2*/
#line 13 "./stack_container.cweb"
;
/*3:*/
#line 73 "./stack_container.cweb"

void FoldedStackContainer::multAndAdd(int dim,const FGSContainer&c,FGSTensor&out)const
{
TL_RAISE_IF(c.num()!=numStacks(),
"Wrong symmetry length of container for FoldedStackContainer::multAndAdd");

THREAD_GROUP gr;
SymmetrySet ss(dim,c.num());
for(symiterator si(ss);!si.isEnd();++si){
if(c.check(*si)){
THREAD*worker= new WorkerFoldMAADense(*this,*si,c,out);
gr.insert(worker);
}
}
gr.run();
}

/*:3*/
#line 14 "./stack_container.cweb"
;
/*4:*/
#line 93 "./stack_container.cweb"

void WorkerFoldMAADense::operator()()
{
Permutation iden(dense_cont.num());
IntSequence coor(sym,iden.getMap());
const FGSTensor*g= dense_cont.get(sym);
cont.multAndAddStacks(coor,*g,out,&out);
}

/*:4*/
#line 15 "./stack_container.cweb"
;
/*5:*/
#line 103 "./stack_container.cweb"

WorkerFoldMAADense::WorkerFoldMAADense(const FoldedStackContainer&container,
const Symmetry&s,
const FGSContainer&dcontainer,
FGSTensor&outten)
:cont(container),sym(s),dense_cont(dcontainer),out(outten)
{}

/*:5*/
#line 16 "./stack_container.cweb"
;
/*6:*/
#line 113 "./stack_container.cweb"

void FoldedStackContainer::multAndAddSparse1(const FSSparseTensor&t,
FGSTensor&out)const
{
THREAD_GROUP gr;
UFSTensor dummy(0,numStacks(),t.dimen());
for(Tensor::index ui= dummy.begin();ui!=dummy.end();++ui){
THREAD*worker= new WorkerFoldMAASparse1(*this,t,out,ui.getCoor());
gr.insert(worker);
}
gr.run();
}

/*:6*/
#line 17 "./stack_container.cweb"
;
/*7:*/
#line 135 "./stack_container.cweb"

void WorkerFoldMAASparse1::operator()()
{
const EquivalenceSet&eset= ebundle.get(out.dimen());
const PermutationSet&pset= tls.pbundle->get(t.dimen());
Permutation iden(t.dimen());

UPSTensor slice(t,cont.getStackSizes(),coor,
PerTensorDimens(cont.getStackSizes(),coor));
for(int iper= 0;iper<pset.getNum();iper++){
const Permutation&per= pset.get(iper);
IntSequence percoor(coor.size());
per.apply(coor,percoor);
for(EquivalenceSet::const_iterator it= eset.begin();
it!=eset.end();++it){
if((*it).numClasses()==t.dimen()){
StackProduct<FGSTensor> sp(cont,*it,out.getSym());
if(!sp.isZero(percoor)){
KronProdStack<FGSTensor> kp(sp,percoor);
kp.optimizeOrder();
const Permutation&oper= kp.getPer();
if(Permutation(oper,per)==iden){
FPSTensor fps(out.getDims(),*it,slice,kp);
{
SYNCHRO syn(&out,"WorkerUnfoldMAASparse1");
fps.addTo(out);
}
}
}
}
}
}
}

/*:7*/
#line 18 "./stack_container.cweb"
;
/*8:*/
#line 170 "./stack_container.cweb"

WorkerFoldMAASparse1::WorkerFoldMAASparse1(const FoldedStackContainer&container,
const FSSparseTensor&ten,
FGSTensor&outten,const IntSequence&c)
:cont(container),t(ten),out(outten),coor(c),ebundle(*(tls.ebundle)){}


/*:8*/
#line 19 "./stack_container.cweb"
;
/*9:*/
#line 183 "./stack_container.cweb"

void FoldedStackContainer::multAndAddSparse2(const FSSparseTensor&t,
FGSTensor&out)const
{
THREAD_GROUP gr;
FFSTensor dummy_f(0,numStacks(),t.dimen());
for(Tensor::index fi= dummy_f.begin();fi!=dummy_f.end();++fi){
THREAD*worker= new WorkerFoldMAASparse2(*this,t,out,fi.getCoor());
gr.insert(worker);
}
gr.run();
}

/*:9*/
#line 20 "./stack_container.cweb"
;
/*10:*/
#line 207 "./stack_container.cweb"

void WorkerFoldMAASparse2::operator()()
{
GSSparseTensor slice(t,cont.getStackSizes(),coor,
TensorDimens(cont.getStackSizes(),coor));
if(slice.getNumNonZero()){
if(slice.getUnfoldIndexFillFactor()> FoldedStackContainer::fill_threshold){
FGSTensor dense_slice(slice);
int r1= slice.getFirstNonZeroRow();
int r2= slice.getLastNonZeroRow();
FGSTensor dense_slice1(r1,r2-r1+1,dense_slice);
FGSTensor out1(r1,r2-r1+1,out);
cont.multAndAddStacks(coor,dense_slice1,out1,&out);
}else
cont.multAndAddStacks(coor,slice,out,&out);
}
}

/*:10*/
#line 21 "./stack_container.cweb"
;
/*11:*/
#line 226 "./stack_container.cweb"

WorkerFoldMAASparse2::WorkerFoldMAASparse2(const FoldedStackContainer&container,
const FSSparseTensor&ten,
FGSTensor&outten,const IntSequence&c)
:cont(container),t(ten),out(outten),coor(c)
{}


/*:11*/
#line 22 "./stack_container.cweb"
;
/*12:*/
#line 248 "./stack_container.cweb"

void FoldedStackContainer::multAndAddSparse3(const FSSparseTensor&t,
FGSTensor&out)const
{
const EquivalenceSet&eset= ebundle.get(out.dimen());
for(Tensor::index run= out.begin();run!=out.end();++run){
Vector outcol(out,*run);
FRSingleTensor sumcol(t.nvar(),t.dimen());
sumcol.zeros();
for(EquivalenceSet::const_iterator it= eset.begin();
it!=eset.end();++it){
if((*it).numClasses()==t.dimen()){
StackProduct<FGSTensor> sp(*this,*it,out.getSym());
IrregTensorHeader header(sp,run.getCoor());
IrregTensor irten(header);
irten.addTo(sumcol);
}
}
t.multColumnAndAdd(sumcol,outcol);
}
}

/*:12*/
#line 23 "./stack_container.cweb"
;
/*13:*/
#line 278 "./stack_container.cweb"

void FoldedStackContainer::multAndAddSparse4(const FSSparseTensor&t,FGSTensor&out)const
{
THREAD_GROUP gr;
FFSTensor dummy_f(0,numStacks(),t.dimen());
for(Tensor::index fi= dummy_f.begin();fi!=dummy_f.end();++fi){
THREAD*worker= new WorkerFoldMAASparse4(*this,t,out,fi.getCoor());
gr.insert(worker);
}
gr.run();
}

/*:13*/
#line 24 "./stack_container.cweb"
;
/*14:*/
#line 294 "./stack_container.cweb"

void WorkerFoldMAASparse4::operator()()
{
GSSparseTensor slice(t,cont.getStackSizes(),coor,
TensorDimens(cont.getStackSizes(),coor));
if(slice.getNumNonZero())
cont.multAndAddStacks(coor,slice,out,&out);
}

/*:14*/
#line 25 "./stack_container.cweb"
;
/*15:*/
#line 304 "./stack_container.cweb"

WorkerFoldMAASparse4::WorkerFoldMAASparse4(const FoldedStackContainer&container,
const FSSparseTensor&ten,
FGSTensor&outten,const IntSequence&c)
:cont(container),t(ten),out(outten),coor(c)
{}


/*:15*/
#line 26 "./stack_container.cweb"
;
/*16:*/
#line 321 "./stack_container.cweb"

void FoldedStackContainer::multAndAddStacks(const IntSequence&coor,
const FGSTensor&g,
FGSTensor&out,const void*ad)const
{
const EquivalenceSet&eset= ebundle.get(out.dimen());

UGSTensor ug(g);
UFSTensor dummy_u(0,numStacks(),g.dimen());
for(Tensor::index ui= dummy_u.begin();ui!=dummy_u.end();++ui){
IntSequence tmp(ui.getCoor());
tmp.sort();
if(tmp==coor){
Permutation sort_per(ui.getCoor());
sort_per.inverse();
for(EquivalenceSet::const_iterator it= eset.begin();
it!=eset.end();++it){
if((*it).numClasses()==g.dimen()){
StackProduct<FGSTensor> sp(*this,*it,sort_per,out.getSym());
if(!sp.isZero(coor)){
KronProdStack<FGSTensor> kp(sp,coor);
if(ug.getSym().isFull())
kp.optimizeOrder();
FPSTensor fps(out.getDims(),*it,sort_per,ug,kp);
{
SYNCHRO syn(ad,"multAndAddStacks");
fps.addTo(out);
}
}
}
}
}
}
}

/*:16*/
#line 27 "./stack_container.cweb"
;
/*17:*/
#line 362 "./stack_container.cweb"

void FoldedStackContainer::multAndAddStacks(const IntSequence&coor,
const GSSparseTensor&g,
FGSTensor&out,const void*ad)const
{
const EquivalenceSet&eset= ebundle.get(out.dimen());
UFSTensor dummy_u(0,numStacks(),g.dimen());
for(Tensor::index ui= dummy_u.begin();ui!=dummy_u.end();++ui){
IntSequence tmp(ui.getCoor());
tmp.sort();
if(tmp==coor){
Permutation sort_per(ui.getCoor());
sort_per.inverse();
for(EquivalenceSet::const_iterator it= eset.begin();
it!=eset.end();++it){
if((*it).numClasses()==g.dimen()){
StackProduct<FGSTensor> sp(*this,*it,sort_per,out.getSym());
if(!sp.isZero(coor)){
KronProdStack<FGSTensor> kp(sp,coor);
FPSTensor fps(out.getDims(),*it,sort_per,g,kp);
{
SYNCHRO syn(ad,"multAndAddStacks");
fps.addTo(out);
}
}
}
}
}
}
}

/*:17*/
#line 28 "./stack_container.cweb"
;

/*18:*/
#line 397 "./stack_container.cweb"

void UnfoldedStackContainer::multAndAdd(const FSSparseTensor&t,
UGSTensor&out)const
{
TL_RAISE_IF(t.nvar()!=getAllSize(),
"Wrong number of variables of tensor for UnfoldedStackContainer::multAndAdd");
multAndAddSparse2(t,out);
}

/*:18*/
#line 30 "./stack_container.cweb"
;
/*19:*/
#line 419 "./stack_container.cweb"

void UnfoldedStackContainer::multAndAdd(int dim,const UGSContainer&c,
UGSTensor&out)const
{
TL_RAISE_IF(c.num()!=numStacks(),
"Wrong symmetry length of container for UnfoldedStackContainer::multAndAdd");

THREAD_GROUP gr;
SymmetrySet ss(dim,c.num());
for(symiterator si(ss);!si.isEnd();++si){
if(c.check(*si)){
THREAD*worker= new WorkerUnfoldMAADense(*this,*si,c,out);
gr.insert(worker);
}
}
gr.run();
}

/*:19*/
#line 31 "./stack_container.cweb"
;
/*20:*/
#line 438 "./stack_container.cweb"

void WorkerUnfoldMAADense::operator()()
{
Permutation iden(dense_cont.num());
IntSequence coor(sym,iden.getMap());
const UGSTensor*g= dense_cont.get(sym);
cont.multAndAddStacks(coor,*g,out,&out);
}

/*:20*/
#line 32 "./stack_container.cweb"
;
/*21:*/
#line 448 "./stack_container.cweb"

WorkerUnfoldMAADense::WorkerUnfoldMAADense(const UnfoldedStackContainer&container,
const Symmetry&s,
const UGSContainer&dcontainer,
UGSTensor&outten)
:cont(container),sym(s),dense_cont(dcontainer),out(outten){}


/*:21*/
#line 33 "./stack_container.cweb"
;
/*22:*/
#line 473 "./stack_container.cweb"

void UnfoldedStackContainer::multAndAddSparse1(const FSSparseTensor&t,
UGSTensor&out)const
{
THREAD_GROUP gr;
UFSTensor dummy(0,numStacks(),t.dimen());
for(Tensor::index ui= dummy.begin();ui!=dummy.end();++ui){
THREAD*worker= new WorkerUnfoldMAASparse1(*this,t,out,ui.getCoor());
gr.insert(worker);
}
gr.run();
}

/*:22*/
#line 34 "./stack_container.cweb"
;
/*23:*/
#line 512 "./stack_container.cweb"

void WorkerUnfoldMAASparse1::operator()()
{
const EquivalenceSet&eset= ebundle.get(out.dimen());
const PermutationSet&pset= tls.pbundle->get(t.dimen());
Permutation iden(t.dimen());

UPSTensor slice(t,cont.getStackSizes(),coor,
PerTensorDimens(cont.getStackSizes(),coor));
for(int iper= 0;iper<pset.getNum();iper++){
const Permutation&per= pset.get(iper);
IntSequence percoor(coor.size());
per.apply(coor,percoor);
for(EquivalenceSet::const_iterator it= eset.begin();
it!=eset.end();++it){
if((*it).numClasses()==t.dimen()){
StackProduct<UGSTensor> sp(cont,*it,out.getSym());
if(!sp.isZero(percoor)){
KronProdStack<UGSTensor> kp(sp,percoor);
kp.optimizeOrder();
const Permutation&oper= kp.getPer();
if(Permutation(oper,per)==iden){
UPSTensor ups(out.getDims(),*it,slice,kp);
{
SYNCHRO syn(&out,"WorkerUnfoldMAASparse1");
ups.addTo(out);
}
}
}
}
}
}
}

/*:23*/
#line 35 "./stack_container.cweb"
;
/*24:*/
#line 547 "./stack_container.cweb"

WorkerUnfoldMAASparse1::WorkerUnfoldMAASparse1(const UnfoldedStackContainer&container,
const FSSparseTensor&ten,
UGSTensor&outten,const IntSequence&c)
:cont(container),t(ten),out(outten),coor(c),ebundle(*(tls.ebundle)){}


/*:24*/
#line 36 "./stack_container.cweb"
;
/*25:*/
#line 575 "./stack_container.cweb"

void UnfoldedStackContainer::multAndAddSparse2(const FSSparseTensor&t,
UGSTensor&out)const
{
THREAD_GROUP gr;
FFSTensor dummy_f(0,numStacks(),t.dimen());
for(Tensor::index fi= dummy_f.begin();fi!=dummy_f.end();++fi){
THREAD*worker= new WorkerUnfoldMAASparse2(*this,t,out,fi.getCoor());
gr.insert(worker);
}
gr.run();
}

/*:25*/
#line 37 "./stack_container.cweb"
;
/*26:*/
#line 595 "./stack_container.cweb"

void WorkerUnfoldMAASparse2::operator()()
{
GSSparseTensor slice(t,cont.getStackSizes(),coor,
TensorDimens(cont.getStackSizes(),coor));
if(slice.getNumNonZero()){
FGSTensor fslice(slice);
UGSTensor dense_slice(fslice);
int r1= slice.getFirstNonZeroRow();
int r2= slice.getLastNonZeroRow();
UGSTensor dense_slice1(r1,r2-r1+1,dense_slice);
UGSTensor out1(r1,r2-r1+1,out);

cont.multAndAddStacks(coor,dense_slice1,out1,&out);
}
}

/*:26*/
#line 38 "./stack_container.cweb"
;
/*27:*/
#line 613 "./stack_container.cweb"

WorkerUnfoldMAASparse2::WorkerUnfoldMAASparse2(const UnfoldedStackContainer&container,
const FSSparseTensor&ten,
UGSTensor&outten,const IntSequence&c)
:cont(container),t(ten),out(outten),coor(c){}


/*:27*/
#line 39 "./stack_container.cweb"
;
/*28:*/
#line 636 "./stack_container.cweb"

void UnfoldedStackContainer::multAndAddStacks(const IntSequence&fi,
const UGSTensor&g,
UGSTensor&out,const void*ad)const
{
const EquivalenceSet&eset= ebundle.get(out.dimen());

UFSTensor dummy_u(0,numStacks(),g.dimen());
for(Tensor::index ui= dummy_u.begin();ui!=dummy_u.end();++ui){
IntSequence tmp(ui.getCoor());
tmp.sort();
if(tmp==fi){
Permutation sort_per(ui.getCoor());
sort_per.inverse();
for(EquivalenceSet::const_iterator it= eset.begin();
it!=eset.end();++it){
if((*it).numClasses()==g.dimen()){
StackProduct<UGSTensor> sp(*this,*it,sort_per,out.getSym());
if(!sp.isZero(fi)){
KronProdStack<UGSTensor> kp(sp,fi);
if(g.getSym().isFull())
kp.optimizeOrder();
UPSTensor ups(out.getDims(),*it,sort_per,g,kp);
{
SYNCHRO syn(ad,"multAndAddStacks");
ups.addTo(out);
}
}
}
}
}
}
}

/*:28*/
#line 40 "./stack_container.cweb"
;


/*:1*/
