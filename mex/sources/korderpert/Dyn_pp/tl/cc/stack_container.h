/*1:*/
#line 80 "./stack_container.hweb"

#ifndef STACK_CONTAINER_H
#define STACK_CONTAINER_H

#include "int_sequence.h"
#include "equivalence.h"
#include "tl_static.h"
#include "t_container.h"
#include "kron_prod.h"
#include "permutation.h"
#include "sthread.h"

/*2:*/
#line 135 "./stack_container.hweb"

template<class _Ttype> 
class StackContainerInterface{
public:
typedef TensorContainer<_Ttype> _Ctype;
typedef enum{matrix,unit,zero}itype;
protected:
const EquivalenceBundle&ebundle;
public:
StackContainerInterface()
:ebundle(*(tls.ebundle)){}
virtual~StackContainerInterface(){}
virtual const IntSequence&getStackSizes()const= 0;
virtual IntSequence&getStackSizes()= 0;
virtual const IntSequence&getStackOffsets()const= 0;
virtual IntSequence&getStackOffsets()= 0;
virtual int numConts()const= 0;
virtual const _Ctype*getCont(int i)const= 0;
virtual itype getType(int i,const Symmetry&s)const= 0;
virtual int numStacks()const= 0;
virtual bool isZero(int i,const Symmetry&s)const= 0;
virtual const _Ttype*getMatrix(int i,const Symmetry&s)const= 0;
virtual int getLengthOfMatrixStacks(const Symmetry&s)const= 0;
virtual int getUnitPos(const Symmetry&s)const= 0;
virtual Vector*createPackedColumn(const Symmetry&s,
const IntSequence&coor,
int&iu)const= 0;
int getAllSize()const
{return getStackOffsets()[numStacks()-1]
+getStackSizes()[numStacks()-1];}
};

/*:2*/
#line 92 "./stack_container.hweb"
;
/*3:*/
#line 171 "./stack_container.hweb"

template<class _Ttype> 
class StackContainer:virtual public StackContainerInterface<_Ttype> {
public:
typedef StackContainerInterface<_Ttype> _Stype;
typedef typename StackContainerInterface<_Ttype> ::_Ctype _Ctype;
typedef typename StackContainerInterface<_Ttype> ::itype itype;
protected:
int num_conts;
IntSequence stack_sizes;
IntSequence stack_offsets;
const _Ctype**const conts;
public:
StackContainer(int ns,int nc)
:num_conts(nc),stack_sizes(ns,0),stack_offsets(ns,0),
conts(new const _Ctype*[nc]){}
virtual~StackContainer(){delete[]conts;}
const IntSequence&getStackSizes()const
{return stack_sizes;}
IntSequence&getStackSizes()
{return stack_sizes;}
const IntSequence&getStackOffsets()const
{return stack_offsets;}
IntSequence&getStackOffsets()
{return stack_offsets;}
int numConts()const
{return num_conts;}
const _Ctype*getCont(int i)const
{return conts[i];}
virtual itype getType(int i,const Symmetry&s)const= 0;
int numStacks()const
{return stack_sizes.size();}
/*4:*/
#line 213 "./stack_container.hweb"

bool isZero(int i,const Symmetry&s)const
{
TL_RAISE_IF(i<0||i>=numStacks(),
"Wrong index to stack in StackContainer::isZero.");
return(getType(i,s)==_Stype::zero||
(getType(i,s)==_Stype::matrix&&!conts[i]->check(s)));
}

/*:4*/
#line 203 "./stack_container.hweb"
;
/*5:*/
#line 223 "./stack_container.hweb"

const _Ttype*getMatrix(int i,const Symmetry&s)const
{
TL_RAISE_IF(isZero(i,s)||getType(i,s)==_Stype::unit,
"Matrix is not returned in StackContainer::getMatrix");
return conts[i]->get(s);
}

/*:5*/
#line 204 "./stack_container.hweb"
;
/*6:*/
#line 232 "./stack_container.hweb"

int getLengthOfMatrixStacks(const Symmetry&s)const
{
int res= 0;
int i= 0;
while(i<numStacks()&&getType(i,s)==_Stype::matrix)
res+= stack_sizes[i++];
return res;
}


/*:6*/
#line 205 "./stack_container.hweb"
;
/*7:*/
#line 244 "./stack_container.hweb"

int getUnitPos(const Symmetry&s)const
{
if(s.dimen()!=1)
return-1;
int i= numStacks()-1;
while(i>=0&&getType(i,s)!=_Stype::unit)
i--;
return i;
}


/*:7*/
#line 206 "./stack_container.hweb"
;
/*8:*/
#line 257 "./stack_container.hweb"

Vector*createPackedColumn(const Symmetry&s,
const IntSequence&coor,int&iu)const
{
TL_RAISE_IF(s.dimen()!=coor.size(),
"Incompatible coordinates for symmetry in StackContainer::createPackedColumn");

int len= getLengthOfMatrixStacks(s);
iu= -1;
int i= 0;
if(-1!=(i= getUnitPos(s))){
iu= stack_offsets[i]+coor[0];
len++;
}

Vector*res= new Vector(len);
i= 0;
while(i<numStacks()&&getType(i,s)==_Stype::matrix){
const _Ttype*t= getMatrix(i,s);
Tensor::index ind(t,coor);
Vector subres(*res,stack_offsets[i],stack_sizes[i]);
subres= ConstVector(ConstGeneralMatrix(*t),*ind);
i++;
}
if(iu!=-1)
(*res)[len-1]= 1;

return res;
}

/*:8*/
#line 207 "./stack_container.hweb"
;
protected:
/*9:*/
#line 288 "./stack_container.hweb"

void calculateOffsets()
{
stack_offsets[0]= 0;
for(int i= 1;i<stack_offsets.size();i++)
stack_offsets[i]= stack_offsets[i-1]+stack_sizes[i-1];
}

/*:9*/
#line 209 "./stack_container.hweb"
;
};

/*:3*/
#line 93 "./stack_container.hweb"
;
/*10:*/
#line 297 "./stack_container.hweb"

class WorkerFoldMAADense;
class WorkerFoldMAASparse1;
class WorkerFoldMAASparse2;
class WorkerFoldMAASparse4;
class FoldedStackContainer:virtual public StackContainerInterface<FGSTensor> {
friend class WorkerFoldMAADense;
friend class WorkerFoldMAASparse1;
friend class WorkerFoldMAASparse2;
friend class WorkerFoldMAASparse4;
public:
static double fill_threshold;
void multAndAdd(int dim,const TensorContainer<FSSparseTensor> &c,
FGSTensor&out)const
{if(c.check(Symmetry(dim)))multAndAdd(*(c.get(Symmetry(dim))),out);}
void multAndAdd(const FSSparseTensor&t,FGSTensor&out)const;
void multAndAdd(int dim,const FGSContainer&c,FGSTensor&out)const;
protected:
void multAndAddSparse1(const FSSparseTensor&t,FGSTensor&out)const;
void multAndAddSparse2(const FSSparseTensor&t,FGSTensor&out)const;
void multAndAddSparse3(const FSSparseTensor&t,FGSTensor&out)const;
void multAndAddSparse4(const FSSparseTensor&t,FGSTensor&out)const;
void multAndAddStacks(const IntSequence&fi,const FGSTensor&g,
FGSTensor&out,const void*ad)const;
void multAndAddStacks(const IntSequence&fi,const GSSparseTensor&g,
FGSTensor&out,const void*ad)const;
};


/*:10*/
#line 94 "./stack_container.hweb"
;
/*11:*/
#line 327 "./stack_container.hweb"

class WorkerUnfoldMAADense;
class WorkerUnfoldMAASparse1;
class WorkerUnfoldMAASparse2;
class UnfoldedStackContainer:virtual public StackContainerInterface<UGSTensor> {
friend class WorkerUnfoldMAADense;
friend class WorkerUnfoldMAASparse1;
friend class WorkerUnfoldMAASparse2;
public:
static double fill_threshold;
void multAndAdd(int dim,const TensorContainer<FSSparseTensor> &c,
UGSTensor&out)const
{if(c.check(Symmetry(dim)))multAndAdd(*(c.get(Symmetry(dim))),out);}
void multAndAdd(const FSSparseTensor&t,UGSTensor&out)const;
void multAndAdd(int dim,const UGSContainer&c,UGSTensor&out)const;
protected:
void multAndAddSparse1(const FSSparseTensor&t,UGSTensor&out)const;
void multAndAddSparse2(const FSSparseTensor&t,UGSTensor&out)const;
void multAndAddStacks(const IntSequence&fi,const UGSTensor&g,
UGSTensor&out,const void*ad)const;
};

/*:11*/
#line 95 "./stack_container.hweb"
;
/*12:*/
#line 360 "./stack_container.hweb"

template<class _Ttype> 
class ZContainer:public StackContainer<_Ttype> {
public:
typedef StackContainer<_Ttype> _Tparent;
typedef StackContainerInterface<_Ttype> _Stype;
typedef typename _Tparent::_Ctype _Ctype;
typedef typename _Tparent::itype itype;
ZContainer(const _Ctype*gss,int ngss,const _Ctype*g,int ng,
int ny,int nu)
:_Tparent(4,2)
{
_Tparent::stack_sizes[0]= ngss;_Tparent::stack_sizes[1]= ng;
_Tparent::stack_sizes[2]= ny;_Tparent::stack_sizes[3]= nu;
_Tparent::conts[0]= gss;
_Tparent::conts[1]= g;
_Tparent::calculateOffsets();
}

/*13:*/
#line 385 "./stack_container.hweb"

itype getType(int i,const Symmetry&s)const
{
if(i==0)
return _Stype::matrix;
if(i==1)
if(s[2]> 0)
return _Stype::zero;
else
return _Stype::matrix;
if(i==2)
if(s==Symmetry(1,0,0,0))
return _Stype::unit;
else
return _Stype::zero;
if(i==3)
if(s==Symmetry(0,1,0,0))
return _Stype::unit;
else
return _Stype::zero;

TL_RAISE("Wrong stack index in ZContainer::getType");
return _Stype::zero;
}

/*:13*/
#line 379 "./stack_container.hweb"
;
};

/*:12*/
#line 96 "./stack_container.hweb"
;
/*14:*/
#line 411 "./stack_container.hweb"

class FoldedZContainer:public ZContainer<FGSTensor> ,
public FoldedStackContainer{
public:
typedef TensorContainer<FGSTensor> _Ctype;
FoldedZContainer(const _Ctype*gss,int ngss,const _Ctype*g,int ng,
int ny,int nu)
:ZContainer<FGSTensor> (gss,ngss,g,ng,ny,nu){}
};

/*:14*/
#line 97 "./stack_container.hweb"
;
/*15:*/
#line 422 "./stack_container.hweb"

class UnfoldedZContainer:public ZContainer<UGSTensor> ,
public UnfoldedStackContainer{
public:
typedef TensorContainer<UGSTensor> _Ctype;
UnfoldedZContainer(const _Ctype*gss,int ngss,const _Ctype*g,int ng,
int ny,int nu)
:ZContainer<UGSTensor> (gss,ngss,g,ng,ny,nu){}
};

/*:15*/
#line 98 "./stack_container.hweb"
;
/*16:*/
#line 442 "./stack_container.hweb"

template<class _Ttype> 
class GContainer:public StackContainer<_Ttype> {
public:
typedef StackContainer<_Ttype> _Tparent;
typedef StackContainerInterface<_Ttype> _Stype;
typedef typename StackContainer<_Ttype> ::_Ctype _Ctype;
typedef typename StackContainer<_Ttype> ::itype itype;
GContainer(const _Ctype*gs,int ngs,int nu)
:StackContainer<_Ttype> (4,1)
{
_Tparent::stack_sizes[0]= ngs;_Tparent::stack_sizes[1]= nu;
_Tparent::stack_sizes[2]= nu;_Tparent::stack_sizes[3]= 1;
_Tparent::conts[0]= gs;
_Tparent::calculateOffsets();
}

/*17:*/
#line 467 "./stack_container.hweb"

itype getType(int i,const Symmetry&s)const
{
if(i==0)
if(s[2]> 0||s==Symmetry(0,0,0,1))
return _Stype::zero;
else
return _Stype::matrix;
if(i==1)
if(s==Symmetry(0,0,1,0))
return _Stype::unit;
else
return _Stype::zero;
if(i==2)
return _Stype::zero;
if(i==3)
if(s==Symmetry(0,0,0,1))
return _Stype::unit;
else
return _Stype::zero;

TL_RAISE("Wrong stack index in GContainer::getType");
return _Stype::zero;
}


/*:17*/
#line 459 "./stack_container.hweb"
;
};

/*:16*/
#line 99 "./stack_container.hweb"
;
/*18:*/
#line 494 "./stack_container.hweb"

class FoldedGContainer:public GContainer<FGSTensor> ,
public FoldedStackContainer{
public:
typedef TensorContainer<FGSTensor> _Ctype;
FoldedGContainer(const _Ctype*gs,int ngs,int nu)
:GContainer<FGSTensor> (gs,ngs,nu){}
};

/*:18*/
#line 100 "./stack_container.hweb"
;
/*19:*/
#line 504 "./stack_container.hweb"

class UnfoldedGContainer:public GContainer<UGSTensor> ,
public UnfoldedStackContainer{
public:
typedef TensorContainer<UGSTensor> _Ctype;
UnfoldedGContainer(const _Ctype*gs,int ngs,int nu)
:GContainer<UGSTensor> (gs,ngs,nu){}
};


/*:19*/
#line 101 "./stack_container.hweb"
;
/*20:*/
#line 520 "./stack_container.hweb"

template<class _Ttype> 
class StackProduct{
public:
typedef StackContainerInterface<_Ttype> _Stype;
typedef typename _Stype::_Ctype _Ctype;
typedef typename _Stype::itype itype;
protected:
const _Stype&stack_cont;
InducedSymmetries syms;
Permutation per;
public:
StackProduct(const _Stype&sc,const Equivalence&e,
const Symmetry&os)
:stack_cont(sc),syms(e,os),per(e){}
StackProduct(const _Stype&sc,const Equivalence&e,
const Permutation&p,const Symmetry&os)
:stack_cont(sc),syms(e,p,os),per(e,p){}
int dimen()const
{return syms.size();}
int getAllSize()const
{return stack_cont.getAllSize();}
const Symmetry&getProdSym(int ip)const
{return syms[ip];}
/*21:*/
#line 553 "./stack_container.hweb"

bool isZero(const IntSequence&istacks)const
{
TL_RAISE_IF(istacks.size()!=dimen(),
"Wrong istacks coordinates for StackProduct::isZero");

bool res= false;
int i= 0;
while(i<dimen()&&!(res= stack_cont.isZero(istacks[i],syms[i])))
i++;
return res;
}

/*:21*/
#line 544 "./stack_container.hweb"
;
/*22:*/
#line 567 "./stack_container.hweb"

itype getType(int is,int ip)const
{
TL_RAISE_IF(is<0||is>=stack_cont.numStacks(),
"Wrong index to stack in StackProduct::getType");
TL_RAISE_IF(ip<0||ip>=dimen(),
"Wrong index to stack container in StackProduct::getType");
return stack_cont.getType(is,syms[ip]);
}

/*:22*/
#line 545 "./stack_container.hweb"
;
/*23:*/
#line 578 "./stack_container.hweb"

const _Ttype*getMatrix(int is,int ip)const
{
return stack_cont.getMatrix(is,syms[ip]);
}

/*:23*/
#line 546 "./stack_container.hweb"
;
/*24:*/
#line 585 "./stack_container.hweb"

void createPackedColumns(const IntSequence&coor,
Vector**vs,IntSequence&iu)const
{
TL_RAISE_IF(iu.size()!=dimen(),
"Wrong storage length for unit flags in StackProduct::createPackedColumn");
TL_RAISE_IF(coor.size()!=per.size(),
"Wrong size of index coor in StackProduct::createPackedColumn");
IntSequence perindex(coor.size());
per.apply(coor,perindex);
int off= 0;
for(int i= 0;i<dimen();i++){
IntSequence percoor(perindex,off,syms[i].dimen()+off);
vs[i]= stack_cont.createPackedColumn(syms[i],percoor,iu[i]);
off+= syms[i].dimen();
}
}

/*:24*/
#line 547 "./stack_container.hweb"
;
/*25:*/
#line 604 "./stack_container.hweb"

int getSize(int is)const
{
return stack_cont.getStackSizes()[is];
}


/*:25*/
#line 548 "./stack_container.hweb"
;
/*26:*/
#line 612 "./stack_container.hweb"

int numMatrices(const IntSequence&istacks)const
{
TL_RAISE_IF(istacks.size()!=dimen(),
"Wrong size of stack coordinates in StackContainer::numMatrices");
int ret= 0;
int ip= 0;
while(ip<dimen()&&getType(istacks[ip],ip)==_Stype::matrix){
ret++;
ip++;
}
return ret;
}

/*:26*/
#line 549 "./stack_container.hweb"
;
};

/*:20*/
#line 102 "./stack_container.hweb"
;
/*27:*/
#line 629 "./stack_container.hweb"

template<class _Ttype> 
class KronProdStack:public KronProdAllOptim{
public:
typedef StackProduct<_Ttype> _Ptype;
typedef StackContainerInterface<_Ttype> _Stype;
/*28:*/
#line 646 "./stack_container.hweb"

KronProdStack(const _Ptype&sp,const IntSequence&istack)
:KronProdAllOptim(sp.dimen())
{
TL_RAISE_IF(sp.dimen()!=istack.size(),
"Wrong stack product dimension for KronProdStack constructor");

for(int i= 0;i<sp.dimen();i++){
TL_RAISE_IF(sp.getType(istack[i],i)==_Stype::zero,
"Attempt to construct KronProdStack from zero matrix");
if(sp.getType(istack[i],i)==_Stype::unit)
setUnit(i,sp.getSize(istack[i]));
if(sp.getType(istack[i],i)==_Stype::matrix){
const TwoDMatrix*m= sp.getMatrix(istack[i],i);
TL_RAISE_IF(m->nrows()!=sp.getSize(istack[i]),
"Wrong size of returned matrix in KronProdStack constructor");
setMat(i,*m);
}
}
}


/*:28*/
#line 635 "./stack_container.hweb"
;
};

/*:27*/
#line 103 "./stack_container.hweb"
;
/*29:*/
#line 669 "./stack_container.hweb"

class WorkerFoldMAADense:public THREAD{
const FoldedStackContainer&cont;
Symmetry sym;
const FGSContainer&dense_cont;
FGSTensor&out;
public:
WorkerFoldMAADense(const FoldedStackContainer&container,
const Symmetry&s,
const FGSContainer&dcontainer,
FGSTensor&outten);
void operator()();
};

/*:29*/
#line 104 "./stack_container.hweb"
;
/*30:*/
#line 684 "./stack_container.hweb"

class WorkerFoldMAASparse1:public THREAD{
const FoldedStackContainer&cont;
const FSSparseTensor&t;
FGSTensor&out;
IntSequence coor;
const EquivalenceBundle&ebundle;
public:
WorkerFoldMAASparse1(const FoldedStackContainer&container,
const FSSparseTensor&ten,
FGSTensor&outten,const IntSequence&c);
void operator()();
};

/*:30*/
#line 105 "./stack_container.hweb"
;
/*31:*/
#line 699 "./stack_container.hweb"

class WorkerFoldMAASparse2:public THREAD{
const FoldedStackContainer&cont;
const FSSparseTensor&t;
FGSTensor&out;
IntSequence coor;
public:
WorkerFoldMAASparse2(const FoldedStackContainer&container,
const FSSparseTensor&ten,
FGSTensor&outten,const IntSequence&c);
void operator()();
};

/*:31*/
#line 106 "./stack_container.hweb"
;
/*32:*/
#line 713 "./stack_container.hweb"

class WorkerFoldMAASparse4:public THREAD{
const FoldedStackContainer&cont;
const FSSparseTensor&t;
FGSTensor&out;
IntSequence coor;
public:
WorkerFoldMAASparse4(const FoldedStackContainer&container,
const FSSparseTensor&ten,
FGSTensor&outten,const IntSequence&c);
void operator()();
};

/*:32*/
#line 107 "./stack_container.hweb"
;
/*33:*/
#line 727 "./stack_container.hweb"

class WorkerUnfoldMAADense:public THREAD{
const UnfoldedStackContainer&cont;
Symmetry sym;
const UGSContainer&dense_cont;
UGSTensor&out;
public:
WorkerUnfoldMAADense(const UnfoldedStackContainer&container,
const Symmetry&s,
const UGSContainer&dcontainer,
UGSTensor&outten);
void operator()();
};

/*:33*/
#line 108 "./stack_container.hweb"
;
/*34:*/
#line 742 "./stack_container.hweb"

class WorkerUnfoldMAASparse1:public THREAD{
const UnfoldedStackContainer&cont;
const FSSparseTensor&t;
UGSTensor&out;
IntSequence coor;
const EquivalenceBundle&ebundle;
public:
WorkerUnfoldMAASparse1(const UnfoldedStackContainer&container,
const FSSparseTensor&ten,
UGSTensor&outten,const IntSequence&c);
void operator()();
};

/*:34*/
#line 109 "./stack_container.hweb"
;
/*35:*/
#line 757 "./stack_container.hweb"

class WorkerUnfoldMAASparse2:public THREAD{
const UnfoldedStackContainer&cont;
const FSSparseTensor&t;
UGSTensor&out;
IntSequence coor;
public:
WorkerUnfoldMAASparse2(const UnfoldedStackContainer&container,
const FSSparseTensor&ten,
UGSTensor&outten,const IntSequence&c);
void operator()();
};


/*:35*/
#line 110 "./stack_container.hweb"
;

#endif

/*:1*/
