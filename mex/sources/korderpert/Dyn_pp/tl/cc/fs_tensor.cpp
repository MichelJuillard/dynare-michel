/*1:*/

#include "fs_tensor.h"
#include "gs_tensor.h"
#include "sparse_tensor.h"
#include "rfs_tensor.h"
#include "tl_exception.h"

/*2:*/

FFSTensor::FFSTensor(const FFSTensor&t,const ConstVector&x)
:FTensor(along_col,IntSequence(t.dimen()-1,t.nvar()),
		 t.nrows(),calcMaxOffset(t.nvar(),t.dimen()-1),t.dimen()-1),
		 nv(t.nvar())
{
	TL_RAISE_IF(t.dimen()<1,
		"Wrong dimension for tensor contraction of FFSTensor");
	TL_RAISE_IF(t.nvar()!=x.length(),
		"Wrong number of variables for tensor contraction of FFSTensor");
	
	zeros();
	
	for(Tensor::index to= begin();to!=end();++to){
		for(int i= 0;i<nvar();i++){
			IntSequence from_ind(i,to.getCoor());
			Tensor::index from(&t,from_ind);
			addColumn(x[i],t,*from,*to);
		}
	}
}


/*:2*/
;
/*3:*/

int FFSTensor::calcMaxOffset(int nvar,int d)
{
	if(nvar==0&&d==0)
		return 1;
	if(nvar==0&&d> 0)
		return 0;
	return noverk(nvar+d-1,d);
}

/*:3*/
;
/*4:*/

FFSTensor::FFSTensor(const FSSparseTensor&t)
:FTensor(along_col,IntSequence(t.dimen(),t.nvar()),
		 t.nrows(),calcMaxOffset(t.nvar(),t.dimen()),t.dimen()),
		 nv(t.nvar())
{
	zeros();
	for(FSSparseTensor::const_iterator it= t.getMap().begin();
	it!=t.getMap().end();++it){
		index ind(this,(*it).first);
		get((*it).second.first,*ind)= (*it).second.second;
	}
}


/*:4*/
;
/*5:*/

FFSTensor::FFSTensor(const UFSTensor&ut)
:FTensor(along_col,IntSequence(ut.dimen(),ut.nvar()),
		 ut.nrows(),calcMaxOffset(ut.nvar(),ut.dimen()),ut.dimen()),
		 nv(ut.nvar())
{
	for(index in= begin();in!=end();++in){
		index src(&ut,in.getCoor());
		copyColumn(ut,*src,*in);
	}
}

/*:5*/
;
/*6:*/

UTensor&FFSTensor::unfold()const
{
	return*(new UFSTensor(*this));
}

/*:6*/
;
/*7:*/

void FFSTensor::increment(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in FFSTensor::increment");
	
	UTensor::increment(v,nv);
	v.monotone();
}

/*:7*/
;
/*8:*/

void FFSTensor::decrement(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in FFSTensor::decrement");
	
	FTensor::decrement(v,nv);
}

/*:8*/
;
/*9:*/

int FFSTensor::getOffset(const IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input vector size in FFSTensor::getOffset");
	
	return FTensor::getOffset(v,nv);
}

/*:9*/
;
/*10:*/

void FFSTensor::addSubTensor(const FGSTensor&t)
{
	TL_RAISE_IF(dimen()!=t.getDims().dimen(),
		"Wrong dimensions for FFSTensor::addSubTensor");
	TL_RAISE_IF(nvar()!=t.getDims().getNVS().sum(),
		"Wrong nvs for FFSTensor::addSubTensor");
	
	/*11:*/
	
	IntSequence shift_pre(t.getSym().num(),0);
	for(int i= 1;i<t.getSym().num();i++)
		shift_pre[i]= shift_pre[i-1]+t.getDims().getNVS()[i-1];
	IntSequence shift(t.getSym(),shift_pre);
	
	/*:11*/
	;
	for(Tensor::index ind= t.begin();ind!=t.end();++ind){
		IntSequence c(ind.getCoor());
		c.add(1,shift);
		c.sort();
		Tensor::index tar(this,c);
		addColumn(t,*ind,*tar);
	}
}

/*:10*/
;
/*12:*/

UFSTensor::UFSTensor(const UFSTensor&t,const ConstVector&x)
:UTensor(along_col,IntSequence(t.dimen()-1,t.nvar()),
		 t.nrows(),calcMaxOffset(t.nvar(),t.dimen()-1),t.dimen()-1),
		 nv(t.nvar())
{
	TL_RAISE_IF(t.dimen()<1,
		"Wrong dimension for tensor contraction of UFSTensor");
	TL_RAISE_IF(t.nvar()!=x.length(),
		"Wrong number of variables for tensor contraction of UFSTensor");
	
	zeros();
	
	for(int i= 0;i<ncols();i++){
		ConstTwoDMatrix tpart(t,i*nvar(),nvar());
		Vector outcol(*this,i);
		tpart.multaVec(outcol,x);
	}
}

/*:12*/
;
/*13:*/

UFSTensor::UFSTensor(const FFSTensor&ft)
:UTensor(along_col,IntSequence(ft.dimen(),ft.nvar()),
		 ft.nrows(),calcMaxOffset(ft.nvar(),ft.dimen()),ft.dimen()),
		 nv(ft.nvar())
{
	for(index src= ft.begin();src!=ft.end();++src){
		index in(this,src.getCoor());
		copyColumn(ft,*src,*in);
	}
	unfoldData();
}

/*:13*/
;
/*14:*/

FTensor&UFSTensor::fold()const
{
	return*(new FFSTensor(*this));
}

/*:14*/
;
/*15:*/

void UFSTensor::increment(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in UFSTensor::increment");
	
	UTensor::increment(v,nv);
}

void UFSTensor::decrement(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in UFSTensor::decrement");
	
	UTensor::decrement(v,nv);
}

/*:15*/
;
/*16:*/

int UFSTensor::getOffset(const IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input vector size in UFSTensor::getOffset");
	
	return UTensor::getOffset(v,nv);
}

/*:16*/
;
/*17:*/

void UFSTensor::addSubTensor(const UGSTensor&t)
{
	TL_RAISE_IF(dimen()!=t.getDims().dimen(),
		"Wrong dimensions for UFSTensor::addSubTensor");
	TL_RAISE_IF(nvar()!=t.getDims().getNVS().sum(),
		"Wrong nvs for UFSTensor::addSubTensor");
	
	/*11:*/
	
	IntSequence shift_pre(t.getSym().num(),0);
	for(int i= 1;i<t.getSym().num();i++)
		shift_pre[i]= shift_pre[i-1]+t.getDims().getNVS()[i-1];
	IntSequence shift(t.getSym(),shift_pre);
	
	/*:11*/
	;
	for(Tensor::index tar= begin();tar!=end();++tar){
		IntSequence c(tar.getCoor());
		c.sort();
		c.add(-1,shift);
		if(c.isPositive()&&c.less(t.getDims().getNVX())){
			Tensor::index from(&t,c);
			addColumn(t,*from,*tar);
		}
	}
}


/*:17*/
;
/*18:*/

void UFSTensor::unfoldData()
{
	for(index in= begin();in!=end();++in){
		IntSequence v(in.getCoor());
		v.sort();
		index tmp(this,v);
		copyColumn(*tmp,*in);
	}
}


/*:18*/
;

/*:1*/
