/*1:*/


#ifndef PS_TENSOR_H
#define PS_TENSOR_H

#include "tensor.h"
#include "gs_tensor.h"
#include "equivalence.h"
#include "permutation.h"
#include "kron_prod.h"
#include "sparse_tensor.h"

/*2:*/

class SortIntSequence:public IntSequence{
public:
	SortIntSequence(const IntSequence&s)
		:IntSequence(s){sort();}
};


/*:2*/
;
/*3:*/

class PerTensorDimens:public TensorDimens{
protected:
	Permutation per;
public:
	PerTensorDimens(const Symmetry&s,const IntSequence&nvars,
		const Equivalence&e)
		:TensorDimens(s,nvars),per(e)
	{per.apply(nvmax);}
	PerTensorDimens(const TensorDimens&td,const Equivalence&e)
		:TensorDimens(td),per(e)
	{per.apply(nvmax);}
	PerTensorDimens(const TensorDimens&td,const Permutation&p)
		:TensorDimens(td),per(p)
	{per.apply(nvmax);}
	PerTensorDimens(const IntSequence&ss,const IntSequence&coor)
		:TensorDimens(ss,SortIntSequence(coor)),per(coor)
	{per.apply(nvmax);}
	PerTensorDimens(const PerTensorDimens&td)
		:TensorDimens(td),per(td.per){}
	const PerTensorDimens&operator= (const PerTensorDimens&td)
	{TensorDimens::operator= (td);per= td.per;return*this;}
	bool operator==(const PerTensorDimens&td)
	{return TensorDimens::operator==(td)&&per==td.per;}
	int tailIdentity()const
	{return per.tailIdentity();}
	const Permutation&getPer()const
	{return per;}
};

/*:3*/
;
/*4:*/

class UPSTensor:public UTensor{
	const PerTensorDimens tdims;
public:
	/*5:*/
	
	UPSTensor(const TensorDimens&td,const Equivalence&e,
		const ConstTwoDMatrix&a,const KronProdAll&kp)
		:UTensor(along_col,PerTensorDimens(td,e).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,e)
	{kp.mult(a,*this);}
	UPSTensor(const TensorDimens&td,const Equivalence&e,
		const ConstTwoDMatrix&a,const KronProdAllOptim&kp)
		:UTensor(along_col,PerTensorDimens(td,Permutation(e,kp.getPer())).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,Permutation(e,kp.getPer()))
	{kp.mult(a,*this);}
	UPSTensor(const TensorDimens&td,const Equivalence&e,const Permutation&p,
		const ConstTwoDMatrix&a,const KronProdAll&kp)
		:UTensor(along_col,PerTensorDimens(td,Permutation(e,p)).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,Permutation(e,p))
	{kp.mult(a,*this);}
	UPSTensor(const TensorDimens&td,const Equivalence&e,const Permutation&p,
		const ConstTwoDMatrix&a,const KronProdAllOptim&kp)
		:UTensor(along_col,PerTensorDimens(td,Permutation(e,Permutation(p,kp.getPer()))).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,Permutation(e,Permutation(p,kp.getPer())))
	{kp.mult(a,*this);}
	
	/*:5*/
	;
	UPSTensor(const FSSparseTensor&t,const IntSequence&ss,
		const IntSequence&coor,const PerTensorDimens&ptd);
	UPSTensor(const UPSTensor&ut)
		:UTensor(ut),tdims(ut.tdims){}
	
	void increment(IntSequence&v)const;
	void decrement(IntSequence&v)const;
	FTensor&fold()const;
	
	int getOffset(const IntSequence&v)const;
	void addTo(FGSTensor&out)const;
	void addTo(UGSTensor&out)const;
	
	enum fill_method{first,second};
	static fill_method decideFillMethod(const FSSparseTensor&t);
private:
	int tailIdentitySize()const;
	void fillFromSparseOne(const FSSparseTensor&t,const IntSequence&ss,
		const IntSequence&coor);
	void fillFromSparseTwo(const FSSparseTensor&t,const IntSequence&ss,
		const IntSequence&coor);
};

/*:4*/
;
/*6:*/

class PerTensorDimens2:public PerTensorDimens{
	InducedSymmetries syms;
	IntSequence ds;
public:
	PerTensorDimens2(const TensorDimens&td,const Equivalence&e,
		const Permutation&p)
		:PerTensorDimens(td,Permutation(e,p)),
		syms(e,p,td.getSym()),
		ds(syms.size())
	{setDimensionSizes();}
	PerTensorDimens2(const TensorDimens&td,const Equivalence&e)
		:PerTensorDimens(td,e),
		syms(e,td.getSym()),
		ds(syms.size())
	{setDimensionSizes();}
	int numSyms()const
	{return(int)syms.size();}
	const Symmetry&getSym(int i)const
	{return syms[i];}
	int calcMaxOffset()const
	{return ds.mult();}
	int calcOffset(const IntSequence&coor)const;
	void print()const;
protected:
	void setDimensionSizes();
};

/*:6*/
;
/*7:*/

template<typename _Ttype> class StackProduct;

class FPSTensor:public FTensor{
	const PerTensorDimens2 tdims;
public:
	/*8:*/
	
	FPSTensor(const TensorDimens&td,const Equivalence&e,
		const ConstTwoDMatrix&a,const KronProdAll&kp)
		:FTensor(along_col,PerTensorDimens(td,e).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,e)
	{kp.mult(a,*this);}
	FPSTensor(const TensorDimens&td,const Equivalence&e,
		const ConstTwoDMatrix&a,const KronProdAllOptim&kp)
		:FTensor(along_col,PerTensorDimens(td,Permutation(e,kp.getPer())).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,e,kp.getPer())
	{kp.mult(a,*this);}
	FPSTensor(const TensorDimens&td,const Equivalence&e,const Permutation&p,
		const ConstTwoDMatrix&a,const KronProdAll&kp)
		:FTensor(along_col,PerTensorDimens(td,Permutation(e,p)).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,e,p)
	{kp.mult(a,*this);}
	FPSTensor(const TensorDimens&td,const Equivalence&e,const Permutation&p,
		const ConstTwoDMatrix&a,const KronProdAllOptim&kp)
		:FTensor(along_col,PerTensorDimens(td,Permutation(e,Permutation(p,kp.getPer()))).getNVX(),
		a.nrows(),kp.ncols(),td.dimen()),tdims(td,e,Permutation(p,kp.getPer()))
	{kp.mult(a,*this);}
	
	FPSTensor(const TensorDimens&td,const Equivalence&e,const Permutation&p,
		const GSSparseTensor&t,const KronProdAll&kp);
	
	FPSTensor(const FPSTensor&ft)
		:FTensor(ft),tdims(ft.tdims){}
	
	/*:8*/
	;
	
	void increment(IntSequence&v)const;
	void decrement(IntSequence&v)const;
	UTensor&unfold()const;
	
	int getOffset(const IntSequence&v)const;
	void addTo(FGSTensor&out)const;
};

/*:7*/
;

#endif

/*:1*/
