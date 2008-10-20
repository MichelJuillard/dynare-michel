/*1:*/


#ifndef TENSOR_H
#define TENSOR_H

#include "int_sequence.h"
#include "twod_matrix.h"

/*2:*/

template<class _Tptr> class _index{
typedef _index<_Tptr> _Self;
_Tptr tensor;
int offset;
IntSequence coor;
public:
_index(_Tptr t,int n)
:tensor(t),offset(0),coor(n,0){}
_index(_Tptr t,const IntSequence&cr,int c)
:tensor(t),offset(c),coor(cr){}
_index(_Tptr t,const IntSequence&cr)
:tensor(t),offset(tensor->getOffset(cr)),coor(cr){}
_index(const _index&ind)
:tensor(ind.tensor),offset(ind.offset),coor(ind.coor){}
const _Self&operator= (const _Self&in)
{tensor= in.tensor;offset= in.offset;coor= in.coor;
return*this;}
_Self&operator++()
{tensor->increment(coor);offset++;return*this;}
_Self&operator--()
{tensor->decrement(coor);offset--;return*this;}
int operator*()const
{return offset;}
bool operator==(const _index&n)const
{return offset==n.offset;}
bool operator!=(const _index&n)const
{return offset!=n.offset;}
const IntSequence&getCoor()const
{return coor;}
void print()const
{printf("%4d: ",offset);coor.print();}
};

/*:2*/
;
/*3:*/

class Tensor:public TwoDMatrix{
public:
enum indor{along_row,along_col};
typedef _index<const Tensor*> index;
protected:
const index in_beg;
const index in_end;
int dim;
public:
Tensor(indor io,const IntSequence&last,int r,int c,int d)
:TwoDMatrix(r,c),
in_beg(this,d),
in_end(this,last,(io==along_row)?r:c),
dim(d){}
Tensor(indor io,const IntSequence&first,const IntSequence&last,
int r,int c,int d)
:TwoDMatrix(r,c),
in_beg(this,first,0),
in_end(this,last,(io==along_row)?r:c),
dim(d){}
Tensor(int first_row,int num,Tensor&t)
:TwoDMatrix(first_row,num,t),
in_beg(t.in_beg),
in_end(t.in_end),
dim(t.dim){}
Tensor(const Tensor&t)
:TwoDMatrix(t),
in_beg(this,t.in_beg.getCoor(),*(t.in_beg)),
in_end(this,t.in_end.getCoor(),*(t.in_end)),
dim(t.dim){}
virtual~Tensor(){}
virtual void increment(IntSequence&v)const= 0;
virtual void decrement(IntSequence&v)const= 0;
virtual int getOffset(const IntSequence&v)const= 0;
int dimen()const
{return dim;}

const index&begin()const
{return in_beg;}
const index&end()const
{return in_end;}

static int noverk(int n,int k);
static int power(int a,int b);
static int noverseq(const IntSequence&s)
{
IntSequence seq(s);
return noverseq_ip((IntSequence&)s);
}
private:
static int noverseq_ip(IntSequence&s);
};

/*:3*/
;
/*4:*/

class FTensor;
class UTensor:public Tensor{
public:
UTensor(indor io,const IntSequence&last,int r,int c,int d)
:Tensor(io,last,r,c,d){}
UTensor(const UTensor&ut)
:Tensor(ut){}
UTensor(int first_row,int num,UTensor&t)
:Tensor(first_row,num,t){}
virtual~UTensor(){}
virtual FTensor&fold()const= 0;

static void increment(IntSequence&v,int nv);
static void decrement(IntSequence&v,int nv);
static void increment(IntSequence&v,const IntSequence&nvmx);
static void decrement(IntSequence&v,const IntSequence&nvmx);
static int getOffset(const IntSequence&v,int nv);
static int getOffset(const IntSequence&v,const IntSequence&nvmx);
};

/*:4*/
;
/*5:*/

class FTensor:public Tensor{
public:
FTensor(indor io,const IntSequence&last,int r,int c,int d)
:Tensor(io,last,r,c,d){}
FTensor(const FTensor&ft)
:Tensor(ft){}
FTensor(int first_row,int num,FTensor&t)
:Tensor(first_row,num,t){}
virtual~FTensor(){}
virtual UTensor&unfold()const= 0;

static void decrement(IntSequence&v,int nv);
static int getOffset(const IntSequence&v,int nv)
{IntSequence vtmp(v);return getOffsetRecurse(vtmp,nv);}
private:
static int getOffsetRecurse(IntSequence&v,int nv);
};

/*:5*/
;

#endif

/*:1*/
