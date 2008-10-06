/*1:*/
#line 6 "./gs_tensor.cweb"

#include "gs_tensor.h"
#include "sparse_tensor.h"
#include "tl_exception.h"
#include "kron_prod.h"

/*2:*/
#line 36 "./gs_tensor.cweb"

TensorDimens::TensorDimens(const IntSequence&ss,const IntSequence&coor)
:nvs(ss),
sym(ss.size(),""),
nvmax(coor.size(),0)
{
TL_RAISE_IF(!coor.isSorted(),
"Coordinates not sorted in TensorDimens slicing constructor");
TL_RAISE_IF(coor[0]<0||coor[coor.size()-1]>=ss.size(),
"A coordinate out of stack range in TensorDimens slicing constructor");

for(int i= 0;i<coor.size();i++){
sym[coor[i]]++;
nvmax[i]= ss[coor[i]];
}
}


/*:2*/
#line 12 "./gs_tensor.cweb"
;
/*3:*/
#line 55 "./gs_tensor.cweb"

int TensorDimens::calcUnfoldMaxOffset()const
{
return nvmax.mult();
}

/*:3*/
#line 13 "./gs_tensor.cweb"
;
/*4:*/
#line 64 "./gs_tensor.cweb"

int TensorDimens::calcFoldMaxOffset()const
{
int res= 1;
for(int i= 0;i<nvs.size();i++){
if(nvs[i]==0&&sym[i]> 0)
return 0;
if(sym[i]> 0)
res*= Tensor::noverk(nvs[i]+sym[i]-1,sym[i]);
}
return res;
}

/*:4*/
#line 14 "./gs_tensor.cweb"
;
/*5:*/
#line 95 "./gs_tensor.cweb"

int TensorDimens::calcFoldOffset(const IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input vector size in TensorDimens::getFoldOffset");

int res= 0;
int pow= 1;
int blstart= v.size();
for(int ibl= getSym().num()-1;ibl>=0;ibl--){
int bldim= getSym()[ibl];
if(bldim> 0){
blstart-= bldim;
int blnvar= getNVX()[blstart];
IntSequence subv(v,blstart,blstart+bldim);
res+= FTensor::getOffset(subv,blnvar)*pow;
pow*= FFSTensor::calcMaxOffset(blnvar,bldim);
}
}
TL_RAISE_IF(blstart!=0,
"Error in tracing symmetry in TensorDimens::getFoldOffset");
return res;
}

/*:5*/
#line 15 "./gs_tensor.cweb"
;
/*6:*/
#line 136 "./gs_tensor.cweb"

void TensorDimens::decrement(IntSequence&v)const
{
TL_RAISE_IF(getNVX().size()!=v.size(),
"Wrong size of input/output sequence in TensorDimens::decrement");

int iblock= getSym().num()-1;
int block_last= v.size();
int block_first= block_last-getSym()[iblock];
/*7:*/
#line 150 "./gs_tensor.cweb"

while(iblock> 0&&v[block_last-1]==0){
for(int i= block_first;i<block_last;i++)
v[i]= getNVX(i);
iblock--;
block_last= block_first;
block_first-= getSym()[iblock];
}

/*:7*/
#line 145 "./gs_tensor.cweb"
;
/*8:*/
#line 160 "./gs_tensor.cweb"

IntSequence vtmp(v,block_first,block_last);
FTensor::decrement(vtmp,getNVX(block_first));



/*:8*/
#line 146 "./gs_tensor.cweb"
;
}

/*:6*/
#line 16 "./gs_tensor.cweb"
;
/*9:*/
#line 169 "./gs_tensor.cweb"

FGSTensor::FGSTensor(const UGSTensor&ut)
:FTensor(along_col,ut.tdims.getNVX(),ut.nrows(),
ut.tdims.calcFoldMaxOffset(),ut.dimen()),
tdims(ut.tdims)
{
for(index ti= begin();ti!=end();++ti){
index ui(&ut,ti.getCoor());
copyColumn(ut,*ui,*ti);
}
}

/*:9*/
#line 17 "./gs_tensor.cweb"
;
/*10:*/
#line 190 "./gs_tensor.cweb"

FGSTensor::FGSTensor(const FSSparseTensor&t,const IntSequence&ss,
const IntSequence&coor,const TensorDimens&td)
:FTensor(along_col,td.getNVX(),t.nrows(),
td.calcFoldMaxOffset(),td.dimen()),
tdims(td)
{
/*11:*/
#line 221 "./gs_tensor.cweb"

IntSequence s_offsets(ss.size(),0);
for(int i= 1;i<ss.size();i++)
s_offsets[i]= s_offsets[i-1]+ss[i-1];

IntSequence lb(coor.size());
IntSequence ub(coor.size());
for(int i= 0;i<coor.size();i++){
lb[i]= s_offsets[coor[i]];
ub[i]= s_offsets[coor[i]]+ss[coor[i]]-1;
}


/*:11*/
#line 197 "./gs_tensor.cweb"
;

zeros();
FSSparseTensor::const_iterator lbi= t.getMap().lower_bound(lb);
FSSparseTensor::const_iterator ubi= t.getMap().upper_bound(ub);
for(FSSparseTensor::const_iterator run= lbi;run!=ubi;++run){
if(lb.lessEq((*run).first)&&(*run).first.lessEq(ub)){
IntSequence c((*run).first);
c.add(-1,lb);
Tensor::index ind(this,c);
TL_RAISE_IF(*ind<0||*ind>=ncols(),
"Internal error in slicing constructor of FGSTensor");
get((*run).second.first,*ind)= (*run).second.second;
}
}
}

/*:10*/
#line 18 "./gs_tensor.cweb"
;
/*12:*/
#line 235 "./gs_tensor.cweb"

FGSTensor::FGSTensor(const FFSTensor&t,const IntSequence&ss,
const IntSequence&coor,const TensorDimens&td)
:FTensor(along_col,td.getNVX(),t.nrows(),
td.calcFoldMaxOffset(),td.dimen()),
tdims(td)
{
if(ncols()==0)
return;

/*11:*/
#line 221 "./gs_tensor.cweb"

IntSequence s_offsets(ss.size(),0);
for(int i= 1;i<ss.size();i++)
s_offsets[i]= s_offsets[i-1]+ss[i-1];

IntSequence lb(coor.size());
IntSequence ub(coor.size());
for(int i= 0;i<coor.size();i++){
lb[i]= s_offsets[coor[i]];
ub[i]= s_offsets[coor[i]]+ss[coor[i]]-1;
}


/*:11*/
#line 245 "./gs_tensor.cweb"
;

zeros();
Tensor::index lbi(&t,lb);
Tensor::index ubi(&t,ub);
++ubi;
for(Tensor::index run= lbi;run!=ubi;++run){
if(lb.lessEq(run.getCoor())&&run.getCoor().lessEq(ub)){
IntSequence c(run.getCoor());
c.add(-1,lb);
Tensor::index ind(this,c);
TL_RAISE_IF(*ind<0||*ind>=ncols(),
"Internal error in slicing constructor of FGSTensor");
copyColumn(t,*run,*ind);
}
}
}

/*:12*/
#line 19 "./gs_tensor.cweb"
;
/*13:*/
#line 264 "./gs_tensor.cweb"

FGSTensor::FGSTensor(const GSSparseTensor&t)
:FTensor(along_col,t.getDims().getNVX(),t.nrows(),
t.getDims().calcFoldMaxOffset(),t.dimen()),tdims(t.getDims())
{
zeros();
for(FSSparseTensor::const_iterator it= t.getMap().begin();
it!=t.getMap().end();++it){
index ind(this,(*it).first);
get((*it).second.first,*ind)= (*it).second.second;
}
}

/*:13*/
#line 20 "./gs_tensor.cweb"
;
/*14:*/
#line 281 "./gs_tensor.cweb"

void FGSTensor::increment(IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input/output vector size in FGSTensor::increment");

UTensor::increment(v,tdims.getNVX());
v.pmonotone(tdims.getSym());
}




/*:14*/
#line 21 "./gs_tensor.cweb"
;
/*15:*/
#line 295 "./gs_tensor.cweb"

UTensor&FGSTensor::unfold()const
{
return*(new UGSTensor(*this));
}


/*:15*/
#line 22 "./gs_tensor.cweb"
;
/*16:*/
#line 321 "./gs_tensor.cweb"

void FGSTensor::contractAndAdd(int i,FGSTensor&out,
const FRSingleTensor&col)const
{
TL_RAISE_IF(i<0||i>=getSym().num(),
"Wrong index for FGSTensor::contractAndAdd");

TL_RAISE_IF(getSym()[i]!=col.dimen()||tdims.getNVS()[i]!=col.nvar(),
"Wrong dimensions for FGSTensor::contractAndAdd");

/*17:*/
#line 349 "./gs_tensor.cweb"

Symmetry sym_left(getSym());
Symmetry sym_right(getSym());
for(int j= 0;j<getSym().num();j++){
if(j<=i)
sym_right[j]= 0;
if(j>=i)
sym_left[j]= 0;
}


/*:17*/
#line 331 "./gs_tensor.cweb"
;
int dleft= TensorDimens(sym_left,tdims.getNVS()).calcFoldMaxOffset();
int dright= TensorDimens(sym_right,tdims.getNVS()).calcFoldMaxOffset();
KronProdAll kp(3);
kp.setUnit(0,dleft);
kp.setMat(1,col);
kp.setUnit(2,dright);
FGSTensor tmp(out.nrows(),out.getDims());
kp.mult(*this,tmp);
out.add(1.0,tmp);
}

/*:16*/
#line 23 "./gs_tensor.cweb"
;
/*18:*/
#line 364 "./gs_tensor.cweb"

UGSTensor::UGSTensor(const FGSTensor&ft)
:UTensor(along_col,ft.tdims.getNVX(),ft.nrows(),
ft.tdims.calcUnfoldMaxOffset(),ft.dimen()),
tdims(ft.tdims)
{
for(index fi= ft.begin();fi!=ft.end();++fi){
index ui(this,fi.getCoor());
copyColumn(ft,*fi,*ui);
}
unfoldData();
}

/*:18*/
#line 24 "./gs_tensor.cweb"
;
/*19:*/
#line 378 "./gs_tensor.cweb"

UGSTensor::UGSTensor(const FSSparseTensor&t,const IntSequence&ss,
const IntSequence&coor,const TensorDimens&td)
:UTensor(along_col,td.getNVX(),t.nrows(),
td.calcUnfoldMaxOffset(),td.dimen()),
tdims(td)
{
if(ncols()==0)
return;

FGSTensor ft(t,ss,coor,td);
for(index fi= ft.begin();fi!=ft.end();++fi){
index ui(this,fi.getCoor());
copyColumn(ft,*fi,*ui);
}
unfoldData();
}

/*:19*/
#line 25 "./gs_tensor.cweb"
;
/*20:*/
#line 397 "./gs_tensor.cweb"

UGSTensor::UGSTensor(const UFSTensor&t,const IntSequence&ss,
const IntSequence&coor,const TensorDimens&td)
:UTensor(along_col,td.getNVX(),t.nrows(),
td.calcUnfoldMaxOffset(),td.dimen()),
tdims(td)
{
FFSTensor folded(t);
FGSTensor ft(folded,ss,coor,td);
for(index fi= ft.begin();fi!=ft.end();++fi){
index ui(this,fi.getCoor());
copyColumn(ft,*fi,*ui);
}
unfoldData();
}


/*:20*/
#line 26 "./gs_tensor.cweb"
;
/*21:*/
#line 415 "./gs_tensor.cweb"

void UGSTensor::increment(IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input/output vector size in UGSTensor::increment");

UTensor::increment(v,tdims.getNVX());
}

void UGSTensor::decrement(IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input/output vector size in UGSTensor::decrement");

UTensor::decrement(v,tdims.getNVX());
}


/*:21*/
#line 27 "./gs_tensor.cweb"
;
/*22:*/
#line 434 "./gs_tensor.cweb"

FTensor&UGSTensor::fold()const
{
return*(new FGSTensor(*this));
}

/*:22*/
#line 28 "./gs_tensor.cweb"
;
/*23:*/
#line 441 "./gs_tensor.cweb"

int UGSTensor::getOffset(const IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input vector size in UGSTensor::getOffset");

return UTensor::getOffset(v,tdims.getNVX());
}

/*:23*/
#line 29 "./gs_tensor.cweb"
;
/*24:*/
#line 453 "./gs_tensor.cweb"

void UGSTensor::unfoldData()
{
for(index in= begin();in!=end();++in)
copyColumn(*(getFirstIndexOf(in)),*in);
}

/*:24*/
#line 30 "./gs_tensor.cweb"
;
/*25:*/
#line 464 "./gs_tensor.cweb"

Tensor::index UGSTensor::getFirstIndexOf(const index&in)const
{
IntSequence v(in.getCoor());
int last= 0;
for(int i= 0;i<tdims.getSym().num();i++){
IntSequence vtmp(v,last,last+tdims.getSym()[i]);
vtmp.sort();
last+= tdims.getSym()[i];
}
return index(this,v);
}

/*:25*/
#line 31 "./gs_tensor.cweb"
;
/*26:*/
#line 480 "./gs_tensor.cweb"

void UGSTensor::contractAndAdd(int i,UGSTensor&out,
const URSingleTensor&col)const
{
TL_RAISE_IF(i<0||i>=getSym().num(),
"Wrong index for UGSTensor::contractAndAdd");
TL_RAISE_IF(getSym()[i]!=col.dimen()||tdims.getNVS()[i]!=col.nvar(),
"Wrong dimensions for UGSTensor::contractAndAdd");

/*17:*/
#line 349 "./gs_tensor.cweb"

Symmetry sym_left(getSym());
Symmetry sym_right(getSym());
for(int j= 0;j<getSym().num();j++){
if(j<=i)
sym_right[j]= 0;
if(j>=i)
sym_left[j]= 0;
}


/*:17*/
#line 489 "./gs_tensor.cweb"
;
int dleft= TensorDimens(sym_left,tdims.getNVS()).calcUnfoldMaxOffset();
int dright= TensorDimens(sym_right,tdims.getNVS()).calcUnfoldMaxOffset();
KronProdAll kp(3);
kp.setUnit(0,dleft);
kp.setMat(1,col);
kp.setUnit(2,dright);
UGSTensor tmp(out.nrows(),out.getDims());
kp.mult(*this,tmp);
out.add(1.0,tmp);
}

/*:26*/
#line 32 "./gs_tensor.cweb"
;

/*:1*/
