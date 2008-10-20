/*1:*/

#ifndef FINE_CONTAINER_H
#define FINE_CONTAINER_H

#include "stack_container.h"

#include <vector> 

/*2:*/

class SizeRefinement{
	vector<int> rsizes;
	vector<int> ind_map;
	int new_nc;
public:
	SizeRefinement(const IntSequence&s,int nc,int max);
	int getRefSize(int i)const
	{return rsizes[i];}
	int numRefinements()const
	{return rsizes.size();}
	int getOldIndex(int i)const
	{return ind_map[i];}
	int getNC()const
	{return new_nc;}
};


/*:2*/
;
/*3:*/

template<class _Ttype> 
class FineContainer:public SizeRefinement,public StackContainer<_Ttype> {
protected:
	typedef StackContainer<_Ttype> _Stype;
	typedef typename StackContainerInterface<_Ttype> ::_Ctype _Ctype;
	typedef typename StackContainerInterface<_Ttype> ::itype itype;
	_Ctype**const ref_conts;
	const _Stype&stack_cont;
public:
	/*4:*/
	
	FineContainer(const _Stype&sc,int max)
		:SizeRefinement(sc.getStackSizes(),sc.numConts(),max),
		StackContainer<_Ttype> (numRefinements(),getNC()),
		ref_conts(new _Ctype*[getNC()]),
		stack_cont(sc)
	{
		for(int i= 0;i<numRefinements();i++)
			_Stype::stack_sizes[i]= getRefSize(i);
		_Stype::calculateOffsets();
		
		int last_cont= -1;
		int last_row= 0;
		for(int i= 0;i<getNC();i++){
			if(getOldIndex(i)!=last_cont){
				last_cont= getOldIndex(i);
				last_row= 0;
			}
			union{const _Ctype*c;_Ctype*n;}convert;
			convert.c= stack_cont.getCont(last_cont);
			ref_conts[i]= new _Ctype(last_row,_Stype::stack_sizes[i],
				*(convert.n));
			_Stype::conts[i]= ref_conts[i];
			last_row+= _Stype::stack_sizes[i];
		}
	}
	
	/*:4*/
	;
	/*5:*/
	
	virtual~FineContainer()
	{
		for(int i= 0;i<_Stype::numConts();i++)
			delete ref_conts[i];
		delete[]ref_conts;
	}
	
	
	
	/*:5*/
	;
	itype getType(int i,const Symmetry&s)const
	{return stack_cont.getType(getOldIndex(i),s);}
	
};


/*:3*/
;
/*6:*/

class FoldedFineContainer:public FineContainer<FGSTensor> ,public FoldedStackContainer{
public:
	FoldedFineContainer(const StackContainer<FGSTensor> &sc,int max)
		:FineContainer<FGSTensor> (sc,max){}
};

/*:6*/
;
/*7:*/

class UnfoldedFineContainer:public FineContainer<UGSTensor> ,public UnfoldedStackContainer{
public:
	UnfoldedFineContainer(const StackContainer<UGSTensor> &sc,int max)
		:FineContainer<UGSTensor> (sc,max){}
};


/*:7*/
;

#endif

/*:1*/
