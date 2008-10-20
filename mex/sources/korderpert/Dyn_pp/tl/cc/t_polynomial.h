/*1:*/

#include "t_container.h"
#include "fs_tensor.h"
#include "rfs_tensor.h"
#include"tl_static.h"

/*2:*/

class PowerProvider{
Vector origv;
URSingleTensor*ut;
FRSingleTensor*ft;
int nv;
public:
PowerProvider(const ConstVector&v)
:origv(v),ut(NULL),ft(NULL),nv(v.length()){}
~PowerProvider();
const URSingleTensor&getNext(const URSingleTensor*dummy);
const FRSingleTensor&getNext(const FRSingleTensor*dummy);
};

/*:2*/
;
/*3:*/

template<class _Ttype,class _TGStype,class _Stype> 
class TensorPolynomial:public TensorContainer<_Ttype> {
int nr;
int nv;
int maxdim;
typedef TensorContainer<_Ttype> _Tparent;
typedef typename _Tparent::_ptr _ptr;
public:
TensorPolynomial(int rows,int vars)
:TensorContainer<_Ttype> (1),
nr(rows),nv(vars),maxdim(0){}
TensorPolynomial(const TensorPolynomial<_Ttype,_TGStype,_Stype> &tp,int k)
:TensorContainer<_Ttype> (tp),
nr(tp.nr),nv(tp.nv),maxdim(0){derivative(k);}
TensorPolynomial(int first_row,int num,TensorPolynomial<_Ttype,_TGStype,_Stype> &tp)
:TensorContainer<_Ttype> (first_row,num,tp),
nr(num),nv(tp.nv),maxdim(tp.maxdim){}
/*4:*/

TensorPolynomial(const TensorPolynomial<_Ttype,_TGStype,_Stype> &tp,const Vector&xval)
:TensorContainer<_Ttype> (1),
nr(tp.nrows()),nv(tp.nvars()-xval.length()),maxdim(0)
{
TL_RAISE_IF(nvars()<0,
"Length of xval too big in TensorPolynomial contract constructor");
IntSequence ss(2);ss[0]= xval.length();ss[1]= nvars();
IntSequence pp(2);pp[0]= 0;pp[1]= 1;

/*5:*/

PowerProvider pwp(xval);
for(int i= 1;i<=tp.maxdim;i++){
const _Stype&xpow= pwp.getNext((const _Stype*)NULL);
for(int j= 0;j<=tp.maxdim-i;j++){
if(tp.check(Symmetry(i+j))){
/*7:*/

_Ttype*ten;
if(_Tparent::check(Symmetry(j))){
ten= _Tparent::get(Symmetry(j));
}else{
ten= new _Ttype(nrows(),nvars(),j);
ten->zeros();
insert(ten);
}


/*:7*/
;
Symmetry sym(i,j);
IntSequence coor(sym,pp);
_TGStype slice(*(tp.get(Symmetry(i+j))),ss,coor,TensorDimens(sym,ss));
slice.mult(Tensor::noverk(i+j,j));
_TGStype tmp(*ten);
slice.contractAndAdd(0,tmp,xpow);
}
}
}

/*:5*/
;
/*6:*/

for(int j= 0;j<=tp.maxdim;j++){
if(tp.check(Symmetry(j))){
/*7:*/

_Ttype*ten;
if(_Tparent::check(Symmetry(j))){
ten= _Tparent::get(Symmetry(j));
}else{
ten= new _Ttype(nrows(),nvars(),j);
ten->zeros();
insert(ten);
}


/*:7*/
;
Symmetry sym(0,j);
IntSequence coor(sym,pp);
_TGStype slice(*(tp.get(Symmetry(j))),ss,coor,TensorDimens(sym,ss));
ten->add(1.0,slice);
}
}


/*:6*/
;
}

/*:4*/
;
TensorPolynomial(const TensorPolynomial&tp)
:TensorContainer<_Ttype> (tp),nr(tp.nr),nv(tp.nv),maxdim(tp.maxdim){}
int nrows()const
{return nr;}
int nvars()const
{return nv;}
/*8:*/

void evalTrad(Vector&out,const ConstVector&v)const
{
if(_Tparent::check(Symmetry(0)))
out= _Tparent::get(Symmetry(0))->getData();
else
out.zeros();

PowerProvider pp(v);
for(int d= 1;d<=maxdim;d++){
const _Stype&p= pp.getNext((const _Stype*)NULL);
Symmetry cs(d);
if(_Tparent::check(cs)){
const _Ttype*t= _Tparent::get(cs);
t->multaVec(out,p.getData());
}
}
}

/*:8*/
;
/*9:*/

void evalHorner(Vector&out,const ConstVector&v)const
{
if(_Tparent::check(Symmetry(0)))
out= _Tparent::get(Symmetry(0))->getData();
else
out.zeros();

if(maxdim==0)
return;

_Ttype*last;
if(maxdim==1)
last= new _Ttype(*(_Tparent::get(Symmetry(1))));
else
last= new _Ttype(*(_Tparent::get(Symmetry(maxdim))),v);
for(int d= maxdim-1;d>=1;d--){
Symmetry cs(d);
if(_Tparent::check(cs)){
const _Ttype*nt= _Tparent::get(cs);
last->add(1.0,ConstTwoDMatrix(*nt));
}
if(d> 1){
_Ttype*new_last= new _Ttype(*last,v);
delete last;
last= new_last;
}
}
last->multaVec(out,v);
delete last;
}

/*:9*/
;
/*10:*/

void insert(_ptr t)
{
TL_RAISE_IF(t->nrows()!=nr,
"Wrong number of rows in TensorPolynomial::insert");
TL_RAISE_IF(t->nvar()!=nv,
"Wrong number of variables in TensorPolynomial::insert");
TensorContainer<_Ttype> ::insert(t);
if(maxdim<t->dimen())
maxdim= t->dimen();
}

/*:10*/
;
/*11:*/

void derivative(int k)
{
for(int d= 1;d<=maxdim;d++){
if(_Tparent::check(Symmetry(d))){
_Ttype*ten= _Tparent::get(Symmetry(d));
ten->mult((double)max((d-k),0));
}
}
}

/*:11*/
;
/*12:*/

_Ttype*evalPartially(int s,const ConstVector&v)
{
TL_RAISE_IF(v.length()!=nvars(),
"Wrong length of vector for TensorPolynomial::evalPartially");

_Ttype*res= new _Ttype(nrows(),nvars(),s);
res->zeros();

int sfac= 1;
for(int i= 1;i<=s;i++)
sfac*= i;

if(_Tparent::check(Symmetry(s)))
res->add(1.0/sfac,*(_Tparent::get(Symmetry(s))));

int dfac= sfac*(s+1);
for(int d= s+1;d<=maxdim;d++,dfac*= d){
if(_Tparent::check(Symmetry(d))){
const _Ttype&ltmp= *(_Tparent::get(Symmetry(d)));
_Ttype*last= new _Ttype(ltmp);
for(int j= 0;j<d-s;j++){
_Ttype*newlast= new _Ttype(*last,v);
delete last;
last= newlast;
}
res->add(1.0/dfac,*last);
delete last;
}
}

return res;
}

/*:12*/
;
};


/*:3*/
;
/*13:*/

class FTensorPolynomial;
class UTensorPolynomial:public TensorPolynomial<UFSTensor,UGSTensor,URSingleTensor> {
public:
UTensorPolynomial(int rows,int vars)
:TensorPolynomial<UFSTensor,UGSTensor,URSingleTensor> (rows,vars){}
UTensorPolynomial(const UTensorPolynomial&up,int k)
:TensorPolynomial<UFSTensor,UGSTensor,URSingleTensor> (up,k){}
UTensorPolynomial(const FTensorPolynomial&fp);
UTensorPolynomial(const UTensorPolynomial&tp,const Vector&xval)
:TensorPolynomial<UFSTensor,UGSTensor,URSingleTensor> (tp,xval){}
UTensorPolynomial(int first_row,int num,UTensorPolynomial&tp)
:TensorPolynomial<UFSTensor,UGSTensor,URSingleTensor> (first_row,num,tp){}
};

/*:13*/
;
/*14:*/

class FTensorPolynomial:public TensorPolynomial<FFSTensor,FGSTensor,FRSingleTensor> {
public:
FTensorPolynomial(int rows,int vars)
:TensorPolynomial<FFSTensor,FGSTensor,FRSingleTensor> (rows,vars){}
FTensorPolynomial(const FTensorPolynomial&fp,int k)
:TensorPolynomial<FFSTensor,FGSTensor,FRSingleTensor> (fp,k){}
FTensorPolynomial(const UTensorPolynomial&up);
FTensorPolynomial(const FTensorPolynomial&tp,const Vector&xval)
:TensorPolynomial<FFSTensor,FGSTensor,FRSingleTensor> (tp,xval){}
FTensorPolynomial(int first_row,int num,FTensorPolynomial&tp)
:TensorPolynomial<FFSTensor,FGSTensor,FRSingleTensor> (first_row,num,tp){}
};

/*:14*/
;
/*15:*/

template<class _Ttype,class _TGStype,class _Stype> 
class CompactPolynomial:public _Ttype{
public:
/*16:*/

CompactPolynomial(const TensorPolynomial<_Ttype,_TGStype,_Stype> &pol)
:_Ttype(pol.nrows(),pol.nvars()+1,pol.getMaxDim())
{
_Ttype::zeros();

IntSequence dumnvs(2);
dumnvs[0]= 1;
dumnvs[1]= pol.nvars();

int offset= 0;
_Ttype dum(0,2,_Ttype::dimen());
for(Tensor::index i= dum.begin();i!=dum.end();++i){
int d= i.getCoor().sum();
Symmetry symrun(_Ttype::dimen()-d,d);
_TGStype dumgs(0,TensorDimens(symrun,dumnvs));
if(pol.check(Symmetry(d))){
TwoDMatrix subt(*this,offset,dumgs.ncols());
subt.add(1.0,*(pol.get(Symmetry(d))));
}
offset+= dumgs.ncols();
}
}


/*:16*/
;
/*17:*/

void eval(Vector&out,const ConstVector&v)const
{
TL_RAISE_IF(v.length()+1!=_Ttype::nvar(),
"Wrong input vector length in CompactPolynomial::eval");
TL_RAISE_IF(out.length()!=_Ttype::nrows(),
"Wrong output vector length in CompactPolynomial::eval");

Vector x1(v.length()+1);
Vector x1p(x1,1,v.length());
x1p= v;
x1[0]= 1.0;

if(_Ttype::dimen()==0)
out= ConstVector(*this,0);
else{
PowerProvider pp(x1);
const _Stype&xpow= pp.getNext((const _Stype*)NULL);
for(int i= 1;i<_Ttype::dimen();i++)
xpow= pp.getNext((const _Stype*)NULL);
multVec(0.0,out,1.0,xpow);
}
}

/*:17*/
;
};

/*:15*/
;
/*18:*/

class UCompactPolynomial:public CompactPolynomial<UFSTensor,UGSTensor,URSingleTensor> {
public:
UCompactPolynomial(const UTensorPolynomial&upol)
:CompactPolynomial<UFSTensor,UGSTensor,URSingleTensor> (upol){}
};

/*:18*/
;
/*19:*/

class FCompactPolynomial:public CompactPolynomial<FFSTensor,FGSTensor,FRSingleTensor> {
public:
FCompactPolynomial(const FTensorPolynomial&fpol)
:CompactPolynomial<FFSTensor,FGSTensor,FRSingleTensor> (fpol){}
};



/*:19*/
;

/*:1*/
