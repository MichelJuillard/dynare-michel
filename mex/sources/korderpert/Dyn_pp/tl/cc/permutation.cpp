/*1:*/
#line 5 "./permutation.cweb"


#include "permutation.h"
#include "tl_exception.h"

/*2:*/
#line 25 "./permutation.cweb"

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
#line 10 "./permutation.cweb"
;
/*3:*/
#line 42 "./permutation.cweb"

void Permutation::inverse()
{
IntSequence former(permap);
for(int i= 0;i<size();i++)
permap[former[i]]= i;
}


/*:3*/
#line 11 "./permutation.cweb"
;
/*4:*/
#line 54 "./permutation.cweb"

int Permutation::tailIdentity()const
{
int i= permap.size();
while(i> 0&&permap[i-1]==i-1)
i--;
return permap.size()-i;
}

/*:4*/
#line 12 "./permutation.cweb"
;
/*5:*/
#line 72 "./permutation.cweb"

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
#line 13 "./permutation.cweb"
;
/*6:*/
#line 91 "./permutation.cweb"

PermutationSet::PermutationSet()
:order(1),size(1),pers(new const Permutation*[size])
{
pers[0]= new Permutation(1);
}

/*:6*/
#line 14 "./permutation.cweb"
;
/*7:*/
#line 99 "./permutation.cweb"

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
#line 15 "./permutation.cweb"
;
/*8:*/
#line 119 "./permutation.cweb"

PermutationSet::~PermutationSet()
{
for(int i= 0;i<size;i++)
if(pers[i])
delete pers[i];
delete[]pers;
}

/*:8*/
#line 16 "./permutation.cweb"
;
/*9:*/
#line 129 "./permutation.cweb"

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
#line 17 "./permutation.cweb"
;
/*10:*/
#line 148 "./permutation.cweb"

PermutationBundle::PermutationBundle(int nmax)
{
nmax= max(nmax,1);
generateUpTo(nmax);
}

/*:10*/
#line 18 "./permutation.cweb"
;
/*11:*/
#line 156 "./permutation.cweb"

PermutationBundle::~PermutationBundle()
{
for(unsigned int i= 0;i<bundle.size();i++)
delete bundle[i];
}

/*:11*/
#line 19 "./permutation.cweb"
;
/*12:*/
#line 164 "./permutation.cweb"

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
#line 20 "./permutation.cweb"
;
/*13:*/
#line 176 "./permutation.cweb"

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
#line 21 "./permutation.cweb"
;


/*:1*/
