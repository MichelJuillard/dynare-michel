/*1:*/
#line 29 "./sparse_tensor.hweb"

#ifndef SPARSE_TENSOR_H
#define SPARSE_TENSOR_H

#include "symmetry.h"
#include "tensor.h"
#include "gs_tensor.h"
#include "Vector.h"

#include <map> 

using namespace std;

/*2:*/
#line 50 "./sparse_tensor.hweb"

struct ltseq{
bool operator()(const IntSequence&s1,const IntSequence&s2)const
{return s1<s2;}
};

/*:2*/
#line 42 "./sparse_tensor.hweb"
;
/*3:*/
#line 60 "./sparse_tensor.hweb"

class SparseTensor{
public:
typedef pair<int,double> Item;
typedef multimap<IntSequence,Item,ltseq> Map;
typedef Map::const_iterator const_iterator;
protected:
typedef Map::iterator iterator;

Map m;
const int dim;
const int nr;
const int nc;
int first_nz_row;
int last_nz_row;
public:
SparseTensor(int d,int nnr,int nnc)
:dim(d),nr(nnr),nc(nnc),first_nz_row(nr),last_nz_row(-1){}
SparseTensor(const SparseTensor&t)
:m(t.m),dim(t.dim),nr(t.nr),nc(t.nc){}
virtual~SparseTensor(){}
void insert(const IntSequence&s,int r,double c);
const Map&getMap()const
{return m;}
int dimen()const
{return dim;}
int nrows()const
{return nr;}
int ncols()const
{return nc;}
double getFillFactor()const
{return((double)m.size())/(nrows()*ncols());}
double getFoldIndexFillFactor()const;
double getUnfoldIndexFillFactor()const;
int getNumNonZero()const
{return m.size();}
int getFirstNonZeroRow()const
{return first_nz_row;}
int getLastNonZeroRow()const
{return last_nz_row;}
virtual const Symmetry&getSym()const= 0;
void print()const;
bool isFinite()const;
}

/*:3*/
#line 43 "./sparse_tensor.hweb"
;
/*4:*/
#line 109 "./sparse_tensor.hweb"

class FSSparseTensor:public SparseTensor{
public:
typedef SparseTensor::const_iterator const_iterator;
private:
const int nv;
const Symmetry sym;
public:
FSSparseTensor(int d,int nvar,int r);
FSSparseTensor(const FSSparseTensor&t);
void insert(const IntSequence&s,int r,double c);
void multColumnAndAdd(const Tensor&t,Vector&v)const;
const Symmetry&getSym()const
{return sym;}
int nvar()const
{return nv;}
void print()const;
};


/*:4*/
#line 44 "./sparse_tensor.hweb"
;
/*5:*/
#line 134 "./sparse_tensor.hweb"

class GSSparseTensor:public SparseTensor{
public:
typedef SparseTensor::const_iterator const_iterator;
private:
const TensorDimens tdims;
public:
GSSparseTensor(const FSSparseTensor&t,const IntSequence&ss,
const IntSequence&coor,const TensorDimens&td);
GSSparseTensor(const GSSparseTensor&t)
:SparseTensor(t),tdims(t.tdims){}
void insert(const IntSequence&s,int r,double c);
const Symmetry&getSym()const
{return tdims.getSym();}
const TensorDimens&getDims()const
{return tdims;}
void print()const;

};

/*:5*/
#line 45 "./sparse_tensor.hweb"
;

#endif

/*:1*/
