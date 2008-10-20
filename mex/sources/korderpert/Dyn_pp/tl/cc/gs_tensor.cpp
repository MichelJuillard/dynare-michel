/*1:*/

#include "gs_tensor.h"
#include "sparse_tensor.h"
#include "tl_exception.h"
#include "kron_prod.h"

/*2:*/

TensorDimens::TensorDimens(const IntSequence&ss,const IntSequence&coor)
:nvs(ss),
sym(ss.size(),""),
nvmax(coor.size(),0)
{
	TL_RAISE_IF(!coor.isSorted(),
		"Coordinates not sorted in TensorDimens slicing constructor");
	TL_RAISE_IF(coor[0]<0||coor[coor.size()-1]>=ss.size(),
		"A coordinate out of stack range in TensorDimens slicing constructor");
	
	for(int i= 0;i<coor.size();i++){
		sym[coor[i]]++;
		nvmax[i]= ss[coor[i]];
	}
}


/*:2*/
;
/*3:*/

int TensorDimens::calcUnfoldMaxOffset()const
{
	return nvmax.mult();
}

/*:3*/
;
/*4:*/

int TensorDimens::calcFoldMaxOffset()const
{
	int res= 1;
	for(int i= 0;i<nvs.size();i++){
		if(nvs[i]==0&&sym[i]> 0)
			return 0;
		if(sym[i]> 0)
			res*= Tensor::noverk(nvs[i]+sym[i]-1,sym[i]);
	}
	return res;
}

/*:4*/
;
/*5:*/

int TensorDimens::calcFoldOffset(const IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input vector size in TensorDimens::getFoldOffset");
	
	int res= 0;
	int pow= 1;
	int blstart= v.size();
	for(int ibl= getSym().num()-1;ibl>=0;ibl--){
		int bldim= getSym()[ibl];
		if(bldim> 0){
			blstart-= bldim;
			int blnvar= getNVX()[blstart];
			IntSequence subv(v,blstart,blstart+bldim);
			res+= FTensor::getOffset(subv,blnvar)*pow;
			pow*= FFSTensor::calcMaxOffset(blnvar,bldim);
		}
	}
	TL_RAISE_IF(blstart!=0,
		"Error in tracing symmetry in TensorDimens::getFoldOffset");
	return res;
}

/*:5*/
;
/*6:*/

void TensorDimens::decrement(IntSequence&v)const
{
	TL_RAISE_IF(getNVX().size()!=v.size(),
		"Wrong size of input/output sequence in TensorDimens::decrement");
	
	int iblock= getSym().num()-1;
	int block_last= v.size();
	int block_first= block_last-getSym()[iblock];
	/*7:*/
	
	while(iblock> 0&&v[block_last-1]==0){
		for(int i= block_first;i<block_last;i++)
			v[i]= getNVX(i);
		iblock--;
		block_last= block_first;
		block_first-= getSym()[iblock];
	}
	
	/*:7*/
	;
	/*8:*/
	
	IntSequence vtmp(v,block_first,block_last);
	FTensor::decrement(vtmp,getNVX(block_first));
	
	
	
	/*:8*/
	;
}

/*:6*/
;
/*9:*/

FGSTensor::FGSTensor(const UGSTensor&ut)
:FTensor(along_col,ut.tdims.getNVX(),ut.nrows(),
		 ut.tdims.calcFoldMaxOffset(),ut.dimen()),
		 tdims(ut.tdims)
{
	for(index ti= begin();ti!=end();++ti){
		index ui(&ut,ti.getCoor());
		copyColumn(ut,*ui,*ti);
	}
}

/*:9*/
;
/*10:*/

FGSTensor::FGSTensor(const FSSparseTensor&t,const IntSequence&ss,
					 const IntSequence&coor,const TensorDimens&td)
					 :FTensor(along_col,td.getNVX(),t.nrows(),
					 td.calcFoldMaxOffset(),td.dimen()),
					 tdims(td)
{
	/*11:*/
	
	IntSequence s_offsets(ss.size(),0);
	for(int i= 1;i<ss.size();i++)
		s_offsets[i]= s_offsets[i-1]+ss[i-1];
	
	IntSequence lb(coor.size());
	IntSequence ub(coor.size());
	for(int i= 0;i<coor.size();i++){
		lb[i]= s_offsets[coor[i]];
		ub[i]= s_offsets[coor[i]]+ss[coor[i]]-1;
	}
	
	
	/*:11*/
	;
	
	zeros();
	FSSparseTensor::const_iterator lbi= t.getMap().lower_bound(lb);
	FSSparseTensor::const_iterator ubi= t.getMap().upper_bound(ub);
	for(FSSparseTensor::const_iterator run= lbi;run!=ubi;++run){
		if(lb.lessEq((*run).first)&&(*run).first.lessEq(ub)){
			IntSequence c((*run).first);
			c.add(-1,lb);
			Tensor::index ind(this,c);
			TL_RAISE_IF(*ind<0||*ind>=ncols(),
				"Internal error in slicing constructor of FGSTensor");
			get((*run).second.first,*ind)= (*run).second.second;
		}
	}
}

/*:10*/
;
/*12:*/

FGSTensor::FGSTensor(const FFSTensor&t,const IntSequence&ss,
					 const IntSequence&coor,const TensorDimens&td)
					 :FTensor(along_col,td.getNVX(),t.nrows(),
					 td.calcFoldMaxOffset(),td.dimen()),
					 tdims(td)
{
	if(ncols()==0)
		return;
	
	/*11:*/
	
	IntSequence s_offsets(ss.size(),0);
	for(int i= 1;i<ss.size();i++)
		s_offsets[i]= s_offsets[i-1]+ss[i-1];
	
	IntSequence lb(coor.size());
	IntSequence ub(coor.size());
	for(int i= 0;i<coor.size();i++){
		lb[i]= s_offsets[coor[i]];
		ub[i]= s_offsets[coor[i]]+ss[coor[i]]-1;
	}
	
	
	/*:11*/
	;
	
	zeros();
	Tensor::index lbi(&t,lb);
	Tensor::index ubi(&t,ub);
	++ubi;
	for(Tensor::index run= lbi;run!=ubi;++run){
		if(lb.lessEq(run.getCoor())&&run.getCoor().lessEq(ub)){
			IntSequence c(run.getCoor());
			c.add(-1,lb);
			Tensor::index ind(this,c);
			TL_RAISE_IF(*ind<0||*ind>=ncols(),
				"Internal error in slicing constructor of FGSTensor");
			copyColumn(t,*run,*ind);
		}
	}
}

/*:12*/
;
/*13:*/

FGSTensor::FGSTensor(const GSSparseTensor&t)
:FTensor(along_col,t.getDims().getNVX(),t.nrows(),
		 t.getDims().calcFoldMaxOffset(),t.dimen()),tdims(t.getDims())
{
	zeros();
	for(FSSparseTensor::const_iterator it= t.getMap().begin();
	it!=t.getMap().end();++it){
		index ind(this,(*it).first);
		get((*it).second.first,*ind)= (*it).second.second;
	}
}

/*:13*/
;
/*14:*/

void FGSTensor::increment(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in FGSTensor::increment");
	
	UTensor::increment(v,tdims.getNVX());
	v.pmonotone(tdims.getSym());
}




/*:14*/
;
/*15:*/

UTensor&FGSTensor::unfold()const
{
	return*(new UGSTensor(*this));
}


/*:15*/
;
/*16:*/

void FGSTensor::contractAndAdd(int i,FGSTensor&out,
							   const FRSingleTensor&col)const
{
	TL_RAISE_IF(i<0||i>=getSym().num(),
		"Wrong index for FGSTensor::contractAndAdd");
	
	TL_RAISE_IF(getSym()[i]!=col.dimen()||tdims.getNVS()[i]!=col.nvar(),
		"Wrong dimensions for FGSTensor::contractAndAdd");
	
	/*17:*/
	
	Symmetry sym_left(getSym());
	Symmetry sym_right(getSym());
	for(int j= 0;j<getSym().num();j++){
		if(j<=i)
			sym_right[j]= 0;
		if(j>=i)
			sym_left[j]= 0;
	}
	
	
	/*:17*/
	;
	int dleft= TensorDimens(sym_left,tdims.getNVS()).calcFoldMaxOffset();
	int dright= TensorDimens(sym_right,tdims.getNVS()).calcFoldMaxOffset();
	KronProdAll kp(3);
	kp.setUnit(0,dleft);
	kp.setMat(1,col);
	kp.setUnit(2,dright);
	FGSTensor tmp(out.nrows(),out.getDims());
	kp.mult(*this,tmp);
	out.add(1.0,tmp);
}

/*:16*/
;
/*18:*/

UGSTensor::UGSTensor(const FGSTensor&ft)
:UTensor(along_col,ft.tdims.getNVX(),ft.nrows(),
		 ft.tdims.calcUnfoldMaxOffset(),ft.dimen()),
		 tdims(ft.tdims)
{
	for(index fi= ft.begin();fi!=ft.end();++fi){
		index ui(this,fi.getCoor());
		copyColumn(ft,*fi,*ui);
	}
	unfoldData();
}

/*:18*/
;
/*19:*/

UGSTensor::UGSTensor(const FSSparseTensor&t,const IntSequence&ss,
					 const IntSequence&coor,const TensorDimens&td)
					 :UTensor(along_col,td.getNVX(),t.nrows(),
					 td.calcUnfoldMaxOffset(),td.dimen()),
					 tdims(td)
{
	if(ncols()==0)
		return;
	
	FGSTensor ft(t,ss,coor,td);
	for(index fi= ft.begin();fi!=ft.end();++fi){
		index ui(this,fi.getCoor());
		copyColumn(ft,*fi,*ui);
	}
	unfoldData();
}

/*:19*/
;
/*20:*/

UGSTensor::UGSTensor(const UFSTensor&t,const IntSequence&ss,
					 const IntSequence&coor,const TensorDimens&td)
					 :UTensor(along_col,td.getNVX(),t.nrows(),
					 td.calcUnfoldMaxOffset(),td.dimen()),
					 tdims(td)
{
	FFSTensor folded(t);
	FGSTensor ft(folded,ss,coor,td);
	for(index fi= ft.begin();fi!=ft.end();++fi){
		index ui(this,fi.getCoor());
		copyColumn(ft,*fi,*ui);
	}
	unfoldData();
}


/*:20*/
;
/*21:*/

void UGSTensor::increment(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in UGSTensor::increment");
	
	UTensor::increment(v,tdims.getNVX());
}

void UGSTensor::decrement(IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input/output vector size in UGSTensor::decrement");
	
	UTensor::decrement(v,tdims.getNVX());
}


/*:21*/
;
/*22:*/

FTensor&UGSTensor::fold()const
{
	return*(new FGSTensor(*this));
}

/*:22*/
;
/*23:*/

int UGSTensor::getOffset(const IntSequence&v)const
{
	TL_RAISE_IF(v.size()!=dimen(),
		"Wrong input vector size in UGSTensor::getOffset");
	
	return UTensor::getOffset(v,tdims.getNVX());
}

/*:23*/
;
/*24:*/

void UGSTensor::unfoldData()
{
	for(index in= begin();in!=end();++in)
		copyColumn(*(getFirstIndexOf(in)),*in);
}

/*:24*/
;
/*25:*/

Tensor::index UGSTensor::getFirstIndexOf(const index&in)const
{
	IntSequence v(in.getCoor());
	int last= 0;
	for(int i= 0;i<tdims.getSym().num();i++){
		IntSequence vtmp(v,last,last+tdims.getSym()[i]);
		vtmp.sort();
		last+= tdims.getSym()[i];
	}
	return index(this,v);
}

/*:25*/
;
/*26:*/

void UGSTensor::contractAndAdd(int i,UGSTensor&out,
							   const URSingleTensor&col)const
{
	TL_RAISE_IF(i<0||i>=getSym().num(),
		"Wrong index for UGSTensor::contractAndAdd");
	TL_RAISE_IF(getSym()[i]!=col.dimen()||tdims.getNVS()[i]!=col.nvar(),
		"Wrong dimensions for UGSTensor::contractAndAdd");
	
	/*17:*/
	
	Symmetry sym_left(getSym());
	Symmetry sym_right(getSym());
	for(int j= 0;j<getSym().num();j++){
		if(j<=i)
			sym_right[j]= 0;
		if(j>=i)
			sym_left[j]= 0;
	}
	
	
	/*:17*/
	;
	int dleft= TensorDimens(sym_left,tdims.getNVS()).calcUnfoldMaxOffset();
	int dright= TensorDimens(sym_right,tdims.getNVS()).calcUnfoldMaxOffset();
	KronProdAll kp(3);
	kp.setUnit(0,dleft);
	kp.setMat(1,col);
	kp.setUnit(2,dright);
	UGSTensor tmp(out.nrows(),out.getDims());
	kp.mult(*this,tmp);
	out.add(1.0,tmp);
}

/*:26*/
;

/*:1*/
