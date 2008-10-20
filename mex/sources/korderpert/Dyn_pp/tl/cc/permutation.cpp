/*1:*/


#include "permutation.h"
#include "tl_exception.h"

/*2:*/

void Permutation::apply(const IntSequence&src,IntSequence&tar)const
{
	TL_RAISE_IF(src.size()!=permap.size()||tar.size()!=permap.size(),
		"Wrong sizes of input or output in Permutation::apply");
	for(int i= 0;i<permap.size();i++)
		tar[i]= src[permap[i]];
}


void Permutation::apply(IntSequence&tar)const
{
	IntSequence tmp(tar);
	apply(tmp,tar);
}

/*:2*/
;
/*3:*/

void Permutation::inverse()
{
	IntSequence former(permap);
	for(int i= 0;i<size();i++)
		permap[former[i]]= i;
}


/*:3*/
;
/*4:*/

int Permutation::tailIdentity()const
{
	int i= permap.size();
	while(i> 0&&permap[i-1]==i-1)
		i--;
	return permap.size()-i;
}

/*:4*/
;
/*5:*/

void Permutation::computeSortingMap(const IntSequence&s)
{
	IntSequence srt(s);
	srt.sort();
	IntSequence flags(s.size(),0);
	
	for(int i= 0;i<s.size();i++){
		int j= 0;
		while(j<s.size()&&(flags[j]||srt[j]!=s[i]))
			j++;
		TL_RAISE_IF(j==s.size(),
			"Internal algorithm error in Permutation::computeSortingMap");
		flags[j]= 1;
		permap[i]= j;
	}
}

/*:5*/
;
/*6:*/

PermutationSet::PermutationSet()
:order(1),size(1),pers(new const Permutation*[size])
{
	pers[0]= new Permutation(1);
}

/*:6*/
;
/*7:*/

PermutationSet::PermutationSet(const PermutationSet&sp,int n)
:order(n),size(n*sp.size),
pers(new const Permutation*[size])
{
	for(int i= 0;i<size;i++)
		pers[i]= NULL;
	
	TL_RAISE_IF(n!=sp.order+1,
		"Wrong new order in PermutationSet constructor");
	
	int k= 0;
	for(int i= 0;i<sp.size;i++){
		for(int j= 0;j<order;j++,k++){
			pers[k]= new Permutation(*(sp.pers[i]),j);
		}
	}
}

/*:7*/
;
/*8:*/

PermutationSet::~PermutationSet()
{
	for(int i= 0;i<size;i++)
		if(pers[i])
			delete pers[i];
		delete[]pers;
}

/*:8*/
;
/*9:*/

vector<const Permutation*> PermutationSet::getPreserving(const IntSequence&s)const
{
	TL_RAISE_IF(s.size()!=order,
		"Wrong sequence length in PermutationSet::getPreserving");
	
	vector<const Permutation*> res;
	IntSequence tmp(s.size());
	for(int i= 0;i<size;i++){
		pers[i]->apply(s,tmp);
		if(s==tmp){
			res.push_back(pers[i]);
		}
	}
	
	return res;
}

/*:9*/
;
/*10:*/

PermutationBundle::PermutationBundle(int nmax)
{
	nmax= max(nmax,1);
	generateUpTo(nmax);
}

/*:10*/
;
/*11:*/

PermutationBundle::~PermutationBundle()
{
	for(unsigned int i= 0;i<bundle.size();i++)
		delete bundle[i];
}

/*:11*/
;
/*12:*/

const PermutationSet&PermutationBundle::get(int n)const
{
	if(n> (int)(bundle.size())||n<1){
		TL_RAISE("Permutation set not found in PermutationSet::get");
		return*(bundle[0]);
	}else{
		return*(bundle[n-1]);
	}
}

/*:12*/
;
/*13:*/

void PermutationBundle::generateUpTo(int nmax)
{
	if(bundle.size()==0)
		bundle.push_back(new PermutationSet());
	
	int curmax= bundle.size();
	for(int n= curmax+1;n<=nmax;n++){
		bundle.push_back(new PermutationSet(*(bundle.back()),n));
	}
}

/*:13*/
;


/*:1*/
