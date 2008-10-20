/*1:*/

#include "symmetry.h"
#include "permutation.h"

#include <stdio.h> 

/*2:*/

Symmetry::Symmetry(const IntSequence&s)
:IntSequence(s.getNumDistinct(),0)
{
	int p= 0;
	if(s.size()> 0)
		operator[](p)= 1;
	for(int i= 1;i<s.size();i++){
		if(s[i]!=s[i-1])
			p++;
		operator[](p)++;
	}
}

/*:2*/
;
/*3:*/

int Symmetry::findClass(int i)const
{
	int j= 0;
	int sum= 0;
	do{
		sum+= operator[](j);
		j++;
	}while(j<size()&&sum<=i);
	
	return j-1;
}

/*:3*/
;
/*4:*/

bool Symmetry::isFull()const
{
	int count= 0;
	for(int i= 0;i<num();i++)
		if(operator[](i)!=0)
			count++;
		return count<=1;
}


/*:4*/
;
/*5:*/

symiterator::symiterator(SymmetrySet&ss)
:s(ss),subit(NULL),subs(NULL),end_flag(false)
{
	s.sym()[0]= 0;
	if(s.size()==2){
		s.sym()[1]= s.dimen();
	}else{
		subs= new SymmetrySet(s,s.dimen());
		subit= new symiterator(*subs);
	}
}


/*:5*/
;
/*6:*/

symiterator::~symiterator()
{
	if(subit)
		delete subit;
	if(subs)
		delete subs;
}

/*:6*/
;
/*7:*/

symiterator&symiterator::operator++()
{
	if(!end_flag){
		if(s.size()==2){
			s.sym()[0]++;
			s.sym()[1]--;
		}else{
			++(*subit);
			if(subit->isEnd()){
				delete subit;
				delete subs;
				s.sym()[0]++;
				subs= new SymmetrySet(s,s.dimen()-s.sym()[0]);
				subit= new symiterator(*subs);
			}
		}
		if(s.sym()[0]==s.dimen()+1)
			end_flag= true;
	}
	return*this;
}

/*:7*/
;
/*8:*/

InducedSymmetries::InducedSymmetries(const Equivalence&e,const Symmetry&s)
{
	for(Equivalence::const_seqit i= e.begin();i!=e.end();++i){
		push_back(Symmetry(s,*i));
	}
}

/*:8*/
;
/*9:*/

InducedSymmetries::InducedSymmetries(const Equivalence&e,const Permutation&p,
									 const Symmetry&s)
{
	for(int i= 0;i<e.numClasses();i++){
		Equivalence::const_seqit it= e.find(p.getMap()[i]);
		push_back(Symmetry(s,*it));
	}
}

/*:9*/
;
/*10:*/

void InducedSymmetries::print()const
{
	printf("Induced symmetries: %d\n",size());
	for(unsigned int i= 0;i<size();i++)
		operator[](i).print();
}

/*:10*/
;

/*:1*/
