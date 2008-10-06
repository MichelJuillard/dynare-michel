/*1:*/
#line 37 "./rfs_tensor.hweb"

#ifndef RFS_TENSOR_H
#define RFS_TENSOR_H

#include "tensor.h"
#include "fs_tensor.h"
#include "symmetry.h"

/*2:*/
#line 53 "./rfs_tensor.hweb"

class FRTensor;
class URTensor:public UTensor{
int nv;
public:
/*3:*/
#line 73 "./rfs_tensor.hweb"

URTensor(int c,int nvar,int d)
:UTensor(along_row,IntSequence(d,nvar),
UFSTensor::calcMaxOffset(nvar,d),c,d),nv(nvar){}
URTensor(const URTensor&ut)
:UTensor(ut),nv(ut.nv){}
URTensor(const FRTensor&ft);

/*:3*/
#line 58 "./rfs_tensor.hweb"
;
virtual~URTensor(){}

void increment(IntSequence&v)const;
void decrement(IntSequence&v)const;
FTensor&fold()const;

int getOffset(const IntSequence&v)const;
int nvar()const
{return nv;}
Symmetry getSym()const
{return Symmetry(dimen());}
};

/*:2*/
#line 45 "./rfs_tensor.hweb"
;
/*4:*/
#line 82 "./rfs_tensor.hweb"

class FRTensor:public FTensor{
int nv;
public:
/*5:*/
#line 102 "./rfs_tensor.hweb"

FRTensor(int c,int nvar,int d)
:FTensor(along_row,IntSequence(d,nvar),
FFSTensor::calcMaxOffset(nvar,d),c,d),nv(nvar){}
FRTensor(const FRTensor&ft)
:FTensor(ft),nv(ft.nv){}
FRTensor(const URTensor&ut);

/*:5*/
#line 86 "./rfs_tensor.hweb"
;
virtual~FRTensor(){}

void increment(IntSequence&v)const;
void decrement(IntSequence&v)const;
UTensor&unfold()const;

int nvar()const
{return nv;}
int getOffset(const IntSequence&v)const
{return FTensor::getOffset(v,nv);}
Symmetry getSym()const
{return Symmetry(dimen());}
};

/*:4*/
#line 46 "./rfs_tensor.hweb"
;
/*6:*/
#line 117 "./rfs_tensor.hweb"

class URSingleTensor:public URTensor{
public:
URSingleTensor(int nvar,int d)
:URTensor(1,nvar,d){}
URSingleTensor(const vector<ConstVector> &cols);
URSingleTensor(const ConstVector&v,int d);
URSingleTensor(const URSingleTensor&ut)
:URTensor(ut){}
virtual~URSingleTensor(){}
FTensor&fold()const;
};

/*:6*/
#line 47 "./rfs_tensor.hweb"
;
/*7:*/
#line 136 "./rfs_tensor.hweb"

class FRSingleTensor:public FRTensor{
public:
FRSingleTensor(int nvar,int d)
:FRTensor(1,nvar,d){}
FRSingleTensor(const URSingleTensor&ut);
FRSingleTensor(const FRSingleTensor&ft)
:FRTensor(ft){}
virtual~FRSingleTensor(){}
};


/*:7*/
#line 48 "./rfs_tensor.hweb"
;

#endif

/*:1*/
