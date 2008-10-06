/*1:*/
#line 6 "./rfs_tensor.cweb"

#include "rfs_tensor.h"
#include "kron_prod.h"
#include "tl_exception.h"

/*2:*/
#line 29 "./rfs_tensor.cweb"

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
#line 11 "./rfs_tensor.cweb"
;
/*3:*/
#line 46 "./rfs_tensor.cweb"

UTensor&FRTensor::unfold()const
{
return*(new URTensor(*this));
}

/*:3*/
#line 12 "./rfs_tensor.cweb"
;
/*4:*/
#line 54 "./rfs_tensor.cweb"

void FRTensor::increment(IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input/output vector size in FRTensor::increment");

UTensor::increment(v,nv);
v.monotone();
}

/*:4*/
#line 13 "./rfs_tensor.cweb"
;
/*5:*/
#line 66 "./rfs_tensor.cweb"

void FRTensor::decrement(IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input/output vector size in FRTensor::decrement");

FTensor::decrement(v,nv);
}


/*:5*/
#line 14 "./rfs_tensor.cweb"
;
/*6:*/
#line 81 "./rfs_tensor.cweb"

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
#line 15 "./rfs_tensor.cweb"
;
/*7:*/
#line 96 "./rfs_tensor.cweb"

FTensor&URTensor::fold()const
{
return*(new FRTensor(*this));
}

/*:7*/
#line 16 "./rfs_tensor.cweb"
;
/*8:*/
#line 103 "./rfs_tensor.cweb"

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
#line 17 "./rfs_tensor.cweb"
;
/*9:*/
#line 121 "./rfs_tensor.cweb"

int URTensor::getOffset(const IntSequence&v)const
{
TL_RAISE_IF(v.size()!=dimen(),
"Wrong input vector size in URTensor::getOffset");

return UTensor::getOffset(v,nv);
}

/*:9*/
#line 18 "./rfs_tensor.cweb"
;
/*10:*/
#line 133 "./rfs_tensor.cweb"

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
#line 19 "./rfs_tensor.cweb"
;
/*11:*/
#line 156 "./rfs_tensor.cweb"

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
#line 20 "./rfs_tensor.cweb"
;
/*12:*/
#line 179 "./rfs_tensor.cweb"

FTensor&URSingleTensor::fold()const
{
return*(new FRSingleTensor(*this));
}



/*:12*/
#line 21 "./rfs_tensor.cweb"
;
/*13:*/
#line 191 "./rfs_tensor.cweb"

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
#line 22 "./rfs_tensor.cweb"
;

/*:1*/
