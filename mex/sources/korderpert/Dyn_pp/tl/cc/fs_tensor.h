/*1:*/

#ifndef FS_TENSOR_H
#define FS_TENSOR_H

#include "tensor.h"
#include "symmetry.h"

class FGSTensor;
class UGSTensor;
class FRSingleTensor;
class FSSparseTensor;
/*2:*/

class UFSTensor;
class FFSTensor:public FTensor{
	int nv;
public:
	/*3:*/
	
	FFSTensor(int r,int nvar,int d)
		:FTensor(along_col,IntSequence(d,nvar),
		r,calcMaxOffset(nvar,d),d),nv(nvar){}
	FFSTensor(const FFSTensor&t,const ConstVector&x);
	FFSTensor(const FSSparseTensor&t);
	FFSTensor(const FFSTensor&ft)
		:FTensor(ft),nv(ft.nv){}
	FFSTensor(const UFSTensor&ut);
	FFSTensor(int first_row,int num,FFSTensor&t)
		:FTensor(first_row,num,t),nv(t.nv){}
	
	
	/*:3*/
	;
	
	void increment(IntSequence&v)const;
	void decrement(IntSequence&v)const;
	UTensor&unfold()const;
	Symmetry getSym()const
	{return Symmetry(dimen());}
	
	int getOffset(const IntSequence&v)const;
	void addSubTensor(const FGSTensor&t);
	int nvar()const
	{return nv;}
	static int calcMaxOffset(int nvar,int d);
};

/*:2*/
;
/*4:*/

class UFSTensor:public UTensor{
	int nv;
public:
	/*5:*/
	
	UFSTensor(int r,int nvar,int d)
		:UTensor(along_col,IntSequence(d,nvar),
		r,calcMaxOffset(nvar,d),d),nv(nvar){}
	UFSTensor(const UFSTensor&t,const ConstVector&x);
	UFSTensor(const UFSTensor&ut)
		:UTensor(ut),nv(ut.nv){}
	UFSTensor(const FFSTensor&ft);
	UFSTensor(int first_row,int num,UFSTensor&t)
		:UTensor(first_row,num,t),nv(t.nv){}
	
	/*:5*/
	;
	
	void increment(IntSequence&v)const;
	void decrement(IntSequence&v)const;
	FTensor&fold()const;
	Symmetry getSym()const
	{return Symmetry(dimen());}
	
	int getOffset(const IntSequence&v)const;
	void addSubTensor(const UGSTensor&t);
	int nvar()const
	{return nv;}
	static int calcMaxOffset(int nvar,int d)
	{return power(nvar,d);}
private:
	void unfoldData();
};

/*:4*/
;

#endif


/*:1*/
