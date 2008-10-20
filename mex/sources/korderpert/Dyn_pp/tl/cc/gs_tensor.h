/*1:*/

#ifndef GS_TENSOR_H
#define GS_TENSOR_H

#include "tensor.h"
#include "fs_tensor.h"
#include "symmetry.h"
#include "rfs_tensor.h"

class FGSTensor;
class UGSTensor;
class FSSparseTensor;

/*2:*/

class TensorDimens{
protected:
	IntSequence nvs;
	Symmetry sym;
	IntSequence nvmax;
public:
	TensorDimens(const Symmetry&s,const IntSequence&nvars)
		:nvs(nvars),sym(s),nvmax(sym,nvs){}
	TensorDimens(int nvar,int dimen)
		:nvs(1),sym(dimen),nvmax(dimen,nvar)
	{nvs[0]= nvar;}
	TensorDimens(const TensorDimens&td)
		:nvs(td.nvs),sym(td.sym),nvmax(td.nvmax){}
	virtual~TensorDimens(){}
	TensorDimens(const IntSequence&ss,const IntSequence&coor);
	const TensorDimens&operator= (const TensorDimens&td)
	{nvs= td.nvs;sym= td.sym;nvmax= td.nvmax;return*this;}
	bool operator==(const TensorDimens&td)const
	{return nvs==td.nvs&&sym==td.sym;}
	bool operator!=(const TensorDimens&td)const
	{return!operator==(td);}
	
	int dimen()const
	{return sym.dimen();}
	int getNVX(int i)const
	{return nvmax[i];}
	const IntSequence&getNVS()const
	{return nvs;}
	const IntSequence&getNVX()const
	{return nvmax;}
	const Symmetry&getSym()const
	{return sym;}
	
	int calcUnfoldMaxOffset()const;
	int calcFoldMaxOffset()const;
	int calcFoldOffset(const IntSequence&v)const;
	void decrement(IntSequence&v)const;
};

/*:2*/
;
/*3:*/

class GSSparseTensor;
class FGSTensor:public FTensor{
	friend class UGSTensor;
	
	const TensorDimens tdims;
public:
	/*4:*/
	
	FGSTensor(int r,const TensorDimens&td)
		:FTensor(along_col,td.getNVX(),r,
		td.calcFoldMaxOffset(),td.dimen()),tdims(td){}
	FGSTensor(const FGSTensor&ft)
		:FTensor(ft),tdims(ft.tdims){}
	FGSTensor(const UGSTensor&ut);
	FGSTensor(int first_row,int num,FGSTensor&t)
		:FTensor(first_row,num,t),tdims(t.tdims){}
	FGSTensor(const FSSparseTensor&t,const IntSequence&ss,
		const IntSequence&coor,const TensorDimens&td);
	FGSTensor(const FFSTensor&t,const IntSequence&ss,
		const IntSequence&coor,const TensorDimens&td);
	FGSTensor(const GSSparseTensor&sp);
	FGSTensor(FFSTensor&t)
		:FTensor(0,t.nrows(),t),tdims(t.nvar(),t.dimen()){}
	
	
	/*:4*/
	;
	virtual~FGSTensor(){}
	
	void increment(IntSequence&v)const;
	void decrement(IntSequence&v)const
	{tdims.decrement(v);}
	UTensor&unfold()const;
	const TensorDimens&getDims()const
	{return tdims;}
	const Symmetry&getSym()const
	{return getDims().getSym();}
	
	void contractAndAdd(int i,FGSTensor&out,
		const FRSingleTensor&col)const;
	int getOffset(const IntSequence&v)const
	{return tdims.calcFoldOffset(v);}
};

/*:3*/
;
/*5:*/

class UGSTensor:public UTensor{
	friend class FGSTensor;
	
	const TensorDimens tdims;
public:
	/*6:*/
	
	UGSTensor(int r,const TensorDimens&td)
		:UTensor(along_col,td.getNVX(),r,
		td.calcUnfoldMaxOffset(),td.dimen()),tdims(td){}
	UGSTensor(const UGSTensor&ut)
		:UTensor(ut),tdims(ut.tdims){}
	UGSTensor(const FGSTensor&ft);
	UGSTensor(int first_row,int num,UGSTensor&t)
		:UTensor(first_row,num,t),tdims(t.tdims){}
	UGSTensor(const FSSparseTensor&t,const IntSequence&ss,
		const IntSequence&coor,const TensorDimens&td);
	UGSTensor(const UFSTensor&t,const IntSequence&ss,
		const IntSequence&coor,const TensorDimens&td);
	UGSTensor(UFSTensor&t)
		:UTensor(0,t.nrows(),t),tdims(t.nvar(),t.dimen()){}
	
	
	/*:6*/
	;
	virtual~UGSTensor(){}
	
	void increment(IntSequence&v)const;
	void decrement(IntSequence&v)const;
	FTensor&fold()const;
	const TensorDimens&getDims()const
	{return tdims;}
	const Symmetry&getSym()const
	{return getDims().getSym();}
	
	void contractAndAdd(int i,UGSTensor&out,
		const URSingleTensor&col)const;
	int getOffset(const IntSequence&v)const;
private:
	void unfoldData();
public:
	index getFirstIndexOf(const index&in)const;
};


/*:5*/
;

#endif

/*:1*/
