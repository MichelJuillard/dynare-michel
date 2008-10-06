/*1:*/
#line 6 "./int_sequence.cweb"

#include "int_sequence.h"
#include "symmetry.h"
#include "tl_exception.h"

#include <stdio.h> 
#include <limits.h> 

/*2:*/
#line 42 "./int_sequence.cweb"

IntSequence::IntSequence(const Symmetry&sy,const IntSequence&se)
:data(new int[sy.dimen()]),length(sy.dimen()),destroy(true)
{
int k= 0;
for(int i= 0;i<sy.num();i++)
for(int j= 0;j<sy[i];j++,k++)
operator[](k)= se[i];
}


/*:2*/
#line 14 "./int_sequence.cweb"
;
/*3:*/
#line 63 "./int_sequence.cweb"

IntSequence::IntSequence(const Symmetry&sy,const vector<int> &se)
:data(new int[sy.num()]),length(sy.num()),destroy(true)
{
TL_RAISE_IF(sy.dimen()<=se[se.size()-1],
"Sequence is not reachable by symmetry in IntSequence()");
for(int i= 0;i<length;i++)
operator[](i)= 0;

for(unsigned int i= 0;i<se.size();i++)
operator[](sy.findClass(se[i]))++;
}

/*:3*/
#line 15 "./int_sequence.cweb"
;
/*4:*/
#line 79 "./int_sequence.cweb"

IntSequence::IntSequence(int i,const IntSequence&s)
:data(new int[s.size()+1]),length(s.size()+1),destroy(true)
{
int j= 0;
while(j<s.size()&&s[j]<i)
j++;
for(int jj= 0;jj<j;jj++)
operator[](jj)= s[jj];
operator[](j)= i;
for(int jj= j;jj<s.size();jj++)
operator[](jj+1)= s[jj];
}

/*:4*/
#line 16 "./int_sequence.cweb"
;
/*5:*/
#line 94 "./int_sequence.cweb"

IntSequence::IntSequence(int i,const IntSequence&s,int pos)
:data(new int[s.size()+1]),length(s.size()+1),destroy(true)
{
TL_RAISE_IF(pos<0||pos> s.size(),
"Wrong position for insertion IntSequence constructor");
for(int jj= 0;jj<pos;jj++)
operator[](jj)= s[jj];
operator[](pos)= i;
for(int jj= pos;jj<s.size();jj++)
operator[](jj+1)= s[jj];
}

/*:5*/
#line 17 "./int_sequence.cweb"
;
/*6:*/
#line 108 "./int_sequence.cweb"

const IntSequence&IntSequence::operator= (const IntSequence&s)
{
TL_RAISE_IF(!destroy&&length!=s.length,
"Wrong length for in-place IntSequence::operator=");
if(destroy&&length!=s.length){
delete[]data;
data= new int[s.length];
destroy= true;
length= s.length;
}
memcpy(data,s.data,sizeof(int)*length);
return*this;
}


/*:6*/
#line 18 "./int_sequence.cweb"
;
/*7:*/
#line 125 "./int_sequence.cweb"

bool IntSequence::operator==(const IntSequence&s)const
{
if(size()!=s.size())
return false;

int i= 0;
while(i<size()&&operator[](i)==s[i])
i++;
return i==size();
}

/*:7*/
#line 19 "./int_sequence.cweb"
;
/*8:*/
#line 139 "./int_sequence.cweb"

bool IntSequence::operator<(const IntSequence&s)const
{
int len= min(size(),s.size());

int i= 0;
while(i<len&&operator[](i)==s[i])
i++;
return(i<s.size()&&(i==size()||operator[](i)<s[i]));
}

/*:8*/
#line 20 "./int_sequence.cweb"
;
/*9:*/
#line 151 "./int_sequence.cweb"

bool IntSequence::lessEq(const IntSequence&s)const
{
TL_RAISE_IF(size()!=s.size(),
"Sequence with different lengths in IntSequence::lessEq");

int i= 0;
while(i<size()&&operator[](i)<=s[i])
i++;
return(i==size());
}

/*:9*/
#line 21 "./int_sequence.cweb"
;
/*10:*/
#line 164 "./int_sequence.cweb"

bool IntSequence::less(const IntSequence&s)const
{
TL_RAISE_IF(size()!=s.size(),
"Sequence with different lengths in IntSequence::less");

int i= 0;
while(i<size()&&operator[](i)<s[i])
i++;
return(i==size());
}

/*:10*/
#line 22 "./int_sequence.cweb"
;
/*11:*/
#line 179 "./int_sequence.cweb"

void IntSequence::sort()
{
for(int i= 0;i<length;i++){
int swaps= 0;
for(int j= 0;j<length-1;j++){
if(data[j]> data[j+1]){
int s= data[j+1];
data[j+1]= data[j];
data[j]= s;
swaps++;
}
}
if(swaps==0)
return;
}
}

/*:11*/
#line 23 "./int_sequence.cweb"
;
/*12:*/
#line 200 "./int_sequence.cweb"

void IntSequence::monotone()
{
for(int i= 1;i<length;i++)
if(data[i-1]> data[i])
data[i]= data[i-1];
}

/*:12*/
#line 24 "./int_sequence.cweb"
;
/*13:*/
#line 213 "./int_sequence.cweb"

void IntSequence::pmonotone(const Symmetry&s)
{
int cum= 0;
for(int i= 0;i<s.num();i++){
for(int j= cum+1;j<cum+s[i];j++)
if(data[j-1]> data[j])
data[j]= data[j-1];
cum+= s[i];
}
}

/*:13*/
#line 25 "./int_sequence.cweb"
;
/*14:*/
#line 226 "./int_sequence.cweb"

int IntSequence::sum()const
{
int res= 0;
for(int i= 0;i<length;i++)
res+= operator[](i);
return res;
}

/*:14*/
#line 26 "./int_sequence.cweb"
;
/*15:*/
#line 238 "./int_sequence.cweb"

int IntSequence::mult(int i1,int i2)const
{
int res= 1;
for(int i= i1;i<i2;i++)
res*= operator[](i);
return res;
}

/*:15*/
#line 27 "./int_sequence.cweb"
;
/*16:*/
#line 248 "./int_sequence.cweb"

int IntSequence::getPrefixLength()const
{
int i= 0;
while(i+1<size()&&operator[](i+1)==operator[](0))
i++;
return i+1;
}

/*:16*/
#line 28 "./int_sequence.cweb"
;
/*17:*/
#line 260 "./int_sequence.cweb"

int IntSequence::getNumDistinct()const
{
int res= 0;
if(size()> 0)
res++;
for(int i= 1;i<size();i++)
if(operator[](i)!=operator[](i-1))
res++;
return res;
}

/*:17*/
#line 29 "./int_sequence.cweb"
;
/*18:*/
#line 275 "./int_sequence.cweb"

int IntSequence::getMax()const
{
int res= INT_MIN;
for(int i= 0;i<size();i++)
if(operator[](i)> res)
res= operator[](i);
return res;
}

/*:18*/
#line 30 "./int_sequence.cweb"
;
/*19:*/
#line 286 "./int_sequence.cweb"

void IntSequence::add(int i)
{
for(int j= 0;j<size();j++)
operator[](j)+= i;
}

/*:19*/
#line 31 "./int_sequence.cweb"
;
/*20:*/
#line 294 "./int_sequence.cweb"

void IntSequence::add(int f,const IntSequence&s)
{
TL_RAISE_IF(size()!=s.size(),
"Wrong sequence length in IntSequence::add");
for(int j= 0;j<size();j++)
operator[](j)+= f*s[j];
}

/*:20*/
#line 32 "./int_sequence.cweb"
;
/*21:*/
#line 304 "./int_sequence.cweb"

bool IntSequence::isPositive()const
{
int i= 0;
while(i<size()&&operator[](i)>=0)
i++;
return(i==size());
}

/*:21*/
#line 33 "./int_sequence.cweb"
;
/*22:*/
#line 314 "./int_sequence.cweb"

bool IntSequence::isConstant()const
{
bool res= true;
int i= 1;
while(res&&i<size()){
res= res&&operator[](0)==operator[](i);
i++;
}
return res;
}

/*:22*/
#line 34 "./int_sequence.cweb"
;
/*23:*/
#line 327 "./int_sequence.cweb"

bool IntSequence::isSorted()const
{
bool res= true;
int i= 1;
while(res&&i<size()){
res= res&&operator[](i-1)<=operator[](i);
i++;
}
return res;
}



/*:23*/
#line 35 "./int_sequence.cweb"
;
/*24:*/
#line 342 "./int_sequence.cweb"

void IntSequence::print()const
{
printf("[");
for(int i= 0;i<size();i++)
printf("%2d ",operator[](i));
printf("]\n");
}

/*:24*/
#line 36 "./int_sequence.cweb"
;

/*:1*/
