/*1:*/

#include "pyramid_prod2.h"
#include "rfs_tensor.h"

/*2:*/

IrregTensorHeader::IrregTensorHeader(const StackProduct<FGSTensor> &sp,
									 const IntSequence&c)
									 :nv(sp.getAllSize()),
									 unit_flag(sp.dimen()),
									 cols(new Vector*[sp.dimen()]),
									 end_seq(sp.dimen())
{
	sp.createPackedColumns(c,cols,unit_flag);
	for(int i= 0;i<sp.dimen();i++){
		end_seq[i]= cols[i]->length();
		if(unit_flag[i]!=-1)
			end_seq[i]= unit_flag[i]+1;
	}
}


/*:2*/
;
/*3:*/

void IrregTensorHeader::increment(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong size of coordinates in IrregTensorHeader::increment");
	
	if(v.size()==0)
		return;
	int i= v.size()-1;
	/*4:*/
	
	v[i]++;
	if(unit_flag[i]!=-1&&v[i]==cols[i]->length()-1)
		v[i]= unit_flag[i];
	
	
	/*:4*/
	;
	while(i> 0&&v[i]==end_seq[i]){
		v[i]= 0;
		i--;
		/*4:*/
		
		v[i]++;
		if(unit_flag[i]!=-1&&v[i]==cols[i]->length()-1)
			v[i]= unit_flag[i];
		
		
		/*:4*/
		;
	}
}

/*:3*/
;
/*5:*/

IrregTensorHeader::~IrregTensorHeader()
{
	for(int i= 0;i<dimen();i++)
		delete cols[i];
	delete[]cols;
}

/*:5*/
;
/*6:*/

int IrregTensorHeader::calcMaxOffset()const
{
	int res= 1;
	for(int i= 0;i<dimen();i++)
		res*= cols[i]->length();
	return res;
}


/*:6*/
;
/*7:*/

IrregTensor::IrregTensor(const IrregTensorHeader&h)
:Tensor(along_row,IntSequence(h.dimen(),0),h.end_seq,
		h.calcMaxOffset(),1,h.dimen()),
		header(h)
{
	if(header.dimen()==1){
		getData()= *(header.cols[0]);
		return;
	}
	
	Vector*last= new Vector(*(header.cols[header.dimen()-1]));
	for(int i= header.dimen()-2;i> 0;i--){
		Vector*newlast= new Vector(last->length()*header.cols[i]->length());
		KronProd::kronMult(ConstVector(*(header.cols[i])),
			ConstVector(*last),*newlast);
		delete last;
		last= newlast;
	}
	KronProd::kronMult(ConstVector(*(header.cols[0])),
		ConstVector(*last),getData());
	delete last;
}

/*:7*/
;
/*8:*/

void IrregTensor::addTo(FRSingleTensor&out)const
{
	for(index it= begin();it!=end();++it){
		IntSequence tmp(it.getCoor());
		tmp.sort();
		Tensor::index ind(&out,tmp);
		out.get(*ind,0)+= get(*it,0);
	}
}


/*:8*/
;

/*:1*/
