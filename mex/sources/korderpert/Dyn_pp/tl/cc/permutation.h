/*1:*/

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "int_sequence.h"
#include "equivalence.h"

#include <vector> 

/*2:*/

class Permutation{
protected:
	IntSequence permap;
public:
	Permutation(int len)
		:permap(len){for(int i= 0;i<len;i++)permap[i]= i;}
	Permutation(const Equivalence&e)
		:permap(e.getN()){e.trace(permap);}
	Permutation(const Equivalence&e,const Permutation&per)
		:permap(e.getN()){e.trace(permap,per);}
	Permutation(const IntSequence&s)
		:permap(s.size()){computeSortingMap(s);};
	Permutation(const Permutation&p)
		:permap(p.permap){}
	Permutation(const Permutation&p1,const Permutation&p2)
		:permap(p2.permap){p1.apply(permap);}
	Permutation(const Permutation&p,int i)
		:permap(p.size(),p.permap,i){}
	const Permutation&operator= (const Permutation&p)
	{permap= p.permap;return*this;}
	bool operator==(const Permutation&p)
	{return permap==p.permap;}
	int size()const
	{return permap.size();}
	void print()const
	{permap.print();}
	void apply(const IntSequence&src,IntSequence&tar)const;
	void apply(IntSequence&tar)const;
	void inverse();
	int tailIdentity()const;
	const IntSequence&getMap()const
	{return permap;}
	IntSequence&getMap()
	{return permap;}
protected:
	void computeSortingMap(const IntSequence&s);
};


/*:2*/
;
/*3:*/

class PermutationSet{
	int order;
	int size;
	const Permutation**const pers;
public:
	PermutationSet();
	PermutationSet(const PermutationSet&ps,int n);
	~PermutationSet();
	int getNum()const
	{return size;}
	const Permutation&get(int i)const
	{return*(pers[i]);}
	vector<const Permutation*> getPreserving(const IntSequence&s)const;
};


/*:3*/
;
/*4:*/

class PermutationBundle{
	vector<PermutationSet*> bundle;
public:
	PermutationBundle(int nmax);
	~PermutationBundle();
	const PermutationSet&get(int n)const;
	void generateUpTo(int nmax);
};

/*:4*/
;

#endif

/*:1*/
