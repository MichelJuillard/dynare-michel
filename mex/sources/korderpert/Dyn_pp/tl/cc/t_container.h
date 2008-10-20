/*1:*/

#ifndef T_CONTAINER_H
#define T_CONTAINER_H

#include "symmetry.h"
#include "gs_tensor.h"
#include "tl_exception.h"
#include "tl_static.h"
#include "sparse_tensor.h"
#include "equivalence.h"
#include "rfs_tensor.h"
#include "Vector.h"

#include <map> 

/*2:*/

struct ltsym{
bool operator()(const Symmetry&s1,const Symmetry&s2)const
{return s1<s2;}
};

/*:2*/
;
/*3:*/

template<class _Ttype> class TensorContainer{
protected:
typedef const _Ttype*_const_ptr;
typedef _Ttype*_ptr;
typedef map<Symmetry,_ptr,ltsym> _Map;
typedef typename _Map::value_type _mvtype;
public:
typedef typename _Map::iterator iterator;
typedef typename _Map::const_iterator const_iterator;
private:
int n;
_Map m;
protected:
const EquivalenceBundle&ebundle;
public:
TensorContainer(int nn)
:n(nn),ebundle(*(tls.ebundle)){}
/*5:*/

TensorContainer(const TensorContainer<_Ttype> &c)
:n(c.n),m(),ebundle(c.ebundle)
{
for(const_iterator it= c.m.begin();it!=c.m.end();++it){
_Ttype*ten= new _Ttype(*((*it).second));
insert(ten);
}
}

/*:5*/
;
/*6:*/

TensorContainer(int first_row,int num,TensorContainer<_Ttype> &c)
:n(c.n),ebundle(*(tls.ebundle))
{
for(iterator it= c.m.begin();it!=c.m.end();++it){
_Ttype*t= new _Ttype(first_row,num,*((*it).second));
insert(t);
}
}


/*:6*/
;
/*7:*/

_const_ptr get(const Symmetry&s)const
{
TL_RAISE_IF(s.num()!=num(),
"Incompatible symmetry lookup in TensorContainer::get");
const_iterator it= m.find(s);
if(it==m.end()){
TL_RAISE("Symmetry not found in TensorContainer::get");
return NULL;
}else{
return(*it).second;
}
}


_ptr get(const Symmetry&s)
{
TL_RAISE_IF(s.num()!=num(),
"Incompatible symmetry lookup in TensorContainer::get");
iterator it= m.find(s);
if(it==m.end()){
TL_RAISE("Symmetry not found in TensorContainer::get");
return NULL;
}else{
return(*it).second;
}
}

/*:7*/
;
/*8:*/

bool check(const Symmetry&s)const
{
TL_RAISE_IF(s.num()!=num(),
"Incompatible symmetry lookup in TensorContainer::check");
const_iterator it= m.find(s);
return it!=m.end();
}

/*:8*/
;
/*9:*/

void insert(_ptr t)
{
TL_RAISE_IF(t->getSym().num()!=num(),
"Incompatible symmetry insertion in TensorContainer::insert");
TL_RAISE_IF(check(t->getSym()),
"Tensor already in container in TensorContainer::insert");
m.insert(_mvtype(t->getSym(),t));
if(!t->isFinite()){
throw TLException(__FILE__,__LINE__,"NaN or Inf asserted in TensorContainer::insert");
}
}

/*:9*/
;
/*10:*/

void remove(const Symmetry&s)
{
iterator it= m.find(s);
if(it!=m.end()){
_ptr t= (*it).second;
m.erase(it);
delete t;
}
}


/*:10*/
;
/*11:*/

void clear()
{
while(!m.empty()){
delete(*(m.begin())).second;
m.erase(m.begin());
}
}

/*:11*/
;
/*15:*/

vector<_const_ptr> 
fetchTensors(const Symmetry&rsym,const Equivalence&e)const
{
vector<_const_ptr> res(e.numClasses());
int i= 0;
for(Equivalence::const_seqit it= e.begin();
it!=e.end();++it,i++){
Symmetry s(rsym,*it);
res[i]= get(s);
}
return res;
}

/*:15*/
;
/*12:*/

int getMaxDim()const
{
int res= -1;
for(const_iterator run= m.begin();run!=m.end();++run){
int dim= (*run).first.dimen();
if(dim> res)
res= dim;
}
return res;
}


/*:12*/
;
/*13:*/

void print()const
{
printf("Tensor container: nvars=%d, tensors=%d\n",n,m.size());
for(const_iterator it= m.begin();it!=m.end();++it){
printf("Symmetry: ");
(*it).first.print();
((*it).second)->print();
}
}

/*:13*/
;
/*14:*/

void writeMat4(FILE*fd,const char*prefix)const
{
for(const_iterator it= begin();it!=end();++it){
char lname[100];
sprintf(lname,"%s_g",prefix);
const Symmetry&sym= (*it).first;
for(int i= 0;i<sym.num();i++){
char tmp[10];
sprintf(tmp,"_%d",sym[i]);
strcat(lname,tmp);
}
ConstTwoDMatrix m(*((*it).second));
m.writeMat4(fd,lname);
}
}


/*:14*/
;

virtual~TensorContainer()
{clear();}

/*4:*/

int num()const
{return n;}
const EquivalenceBundle&getEqBundle()const
{return ebundle;}

const_iterator begin()const
{return m.begin();}
const_iterator end()const
{return m.end();}
iterator begin()
{return m.begin();}
iterator end()
{return m.end();}

/*:4*/
;
};

/*:3*/
;
/*16:*/

class FGSContainer;
class UGSContainer:public TensorContainer<UGSTensor> {
public:
UGSContainer(int nn)
:TensorContainer<UGSTensor> (nn){}
UGSContainer(const UGSContainer&uc)
:TensorContainer<UGSTensor> (uc){}
UGSContainer(const FGSContainer&c);
void multAndAdd(const UGSTensor&t,UGSTensor&out)const;
};


/*:16*/
;
/*17:*/

class FGSContainer:public TensorContainer<FGSTensor> {
static const int num_one_time;
public:
FGSContainer(int nn)
:TensorContainer<FGSTensor> (nn){}
FGSContainer(const FGSContainer&fc)
:TensorContainer<FGSTensor> (fc){}
FGSContainer(const UGSContainer&c);
void multAndAdd(const FGSTensor&t,FGSTensor&out)const;
void multAndAdd(const UGSTensor&t,FGSTensor&out)const;
private:
static Tensor::index
getIndices(int num,vector<IntSequence> &out,
const Tensor::index&start,
const Tensor::index&end);
};


/*:17*/
;

#endif

/*:1*/
