/*1:*/

#ifndef RFS_TENSOR_H
#define RFS_TENSOR_H

#include "tensor.h"
#include "fs_tensor.h"
#include "symmetry.h"

/*2:*/

class FRTensor;
class URTensor:public UTensor{
	int nv;
public:
	/*3:*/
	
	URTensor(int c,int nvar,int d)
		:UTensor(along_row,IntSequence(d,nvar),
		UFSTensor::calcMaxOffset(nvar,d),c,d),nv(nvar){}
	URTensor(const URTensor&ut)
		:UTensor(ut),nv(ut.nv){}
	URTensor(const FRTensor&ft);
	
	/*:3*/
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
;
/*4:*/

class FRTensor:public FTensor{
	int nv;
public:
	/*5:*/
	
	FRTensor(int c,int nvar,int d)
		:FTensor(along_row,IntSequence(d,nvar),
		FFSTensor::calcMaxOffset(nvar,d),c,d),nv(nvar){}
	FRTensor(const FRTensor&ft)
		:FTensor(ft),nv(ft.nv){}
	FRTensor(const URTensor&ut);
	
	/*:5*/
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
;
/*6:*/

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
;
/*7:*/

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
;

#endif

/*:1*/
