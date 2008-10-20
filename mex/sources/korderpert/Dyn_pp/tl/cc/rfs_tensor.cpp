/*1:*/

#include "rfs_tensor.h"
#include "kron_prod.h"
#include "tl_exception.h"

/*2:*/

FRTensor::FRTensor(const URTensor&ut)
:FTensor(along_row,IntSequence(ut.dimen(),ut.nvar()),
		 FFSTensor::calcMaxOffset(ut.nvar(),ut.dimen()),ut.ncols(),
		 ut.dimen()),
		 nv(ut.nvar())
{
	zeros();
	for(index in= ut.begin();in!=ut.end();++in){
		IntSequence vtmp(in.getCoor());
		vtmp.sort();
		index tar(this,vtmp);
		addRow(ut,*in,*tar);
	}
}

/*:2*/
;
/*3:*/

UTensor&FRTensor::unfold()const
{
	return*(new URTensor(*this));
}

/*:3*/
;
/*4:*/

void FRTensor::increment(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in FRTensor::increment");
	
	UTensor::increment(v,nv);
	v.monotone();
}

/*:4*/
;
/*5:*/

void FRTensor::decrement(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in FRTensor::decrement");
	
	FTensor::decrement(v,nv);
}


/*:5*/
;
/*6:*/

URTensor::URTensor(const FRTensor&ft)
:UTensor(along_row,IntSequence(ft.dimen(),ft.nvar()),
		 UFSTensor::calcMaxOffset(ft.nvar(),ft.dimen()),ft.ncols(),
		 ft.dimen()),
		 nv(ft.nvar())
{
	zeros();
	for(index src= ft.begin();src!=ft.end();++src){
		index in(this,src.getCoor());
		copyRow(ft,*src,*in);
	}
}

/*:6*/
;
/*7:*/

FTensor&URTensor::fold()const
{
	return*(new FRTensor(*this));
}

/*:7*/
;
/*8:*/

void URTensor::increment(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in URTensor::increment");
	
	UTensor::increment(v,nv);
}

void URTensor::decrement(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in URTensor::decrement");
	
	UTensor::decrement(v,nv);
}

/*:8*/
;
/*9:*/

int URTensor::getOffset(const IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input vector size in URTensor::getOffset");
	
	return UTensor::getOffset(v,nv);
}

/*:9*/
;
/*10:*/

URSingleTensor::URSingleTensor(const vector<ConstVector> &cols)
:URTensor(1,cols[0].length(),cols.size())
{
	if(dimen()==1){
		getData()= cols[0];
		return;
	}
	
	Vector*last= new Vector(cols[cols.size()-1]);
	for(int i= cols.size()-2;i> 0;i--){
		Vector*newlast= new Vector(Tensor::power(nvar(),cols.size()-i));
		KronProd::kronMult(cols[i],ConstVector(*last),*newlast);
		delete last;
		last= newlast;
	}
	KronProd::kronMult(cols[0],ConstVector(*last),getData());
	delete last;
}

/*:10*/
;
/*11:*/

URSingleTensor::URSingleTensor(const ConstVector&v,int d)
:URTensor(1,v.length(),d)
{
	if(d==1){
		getData()= v;
		return;
	}
	
	Vector*last= new Vector(v);
	for(int i= d-2;i> 0;i--){
		Vector*newlast= new Vector(last->length()*v.length());
		KronProd::kronMult(v,ConstVector(*last),*newlast);
		delete last;
		last= newlast;
	}
	KronProd::kronMult(v,ConstVector(*last),getData());
	delete last;
}

/*:11*/
;
/*12:*/

FTensor&URSingleTensor::fold()const
{
	return*(new FRSingleTensor(*this));
}



/*:12*/
;
/*13:*/

FRSingleTensor::FRSingleTensor(const URSingleTensor&ut)
:FRTensor(1,ut.nvar(),ut.dimen())
{
	zeros();
	for(index in= ut.begin();in!=ut.end();++in){
		IntSequence vtmp(in.getCoor());
		vtmp.sort();
		index tar(this,vtmp);
		get(*tar,0)+= ut.get(*in,0);
	}
}


/*:13*/
;

/*:1*/
