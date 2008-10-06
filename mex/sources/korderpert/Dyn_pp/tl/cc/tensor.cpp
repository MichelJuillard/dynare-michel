/*1:*/
#line 6 "./tensor.cweb"

#include "tensor.h"
#include "tl_exception.h"
#include "tl_static.h"

/*2:*/
#line 27 "./tensor.cweb"

int Tensor::noverk(int n,int k)
{
return tls.ptriang->noverk(n,k);
}

int Tensor::power(int a,int b)
{
int res= 1;
for(int i= 0;i<b;i++)
res*= a;
return res;
}

/*:2*/
#line 11 "./tensor.cweb"
;
/*3:*/
#line 51 "./tensor.cweb"

int Tensor::noverseq_ip(IntSequence&s)
{
if(s.size()==0||s.size()==1)
return 1;
s[1]+= s[0];
return noverk(s[1],s[0])*noverseq(IntSequence(s,1,s.size()));
}

/*:3*/
#line 12 "./tensor.cweb"
;
/*4:*/
#line 65 "./tensor.cweb"

void UTensor::increment(IntSequence&v,int nv)
{
if(v.size()==0)
return;
int i= v.size()-1;
v[i]++;
while(i> 0&&v[i]==nv){
v[i]= 0;
v[--i]++;
}
}

/*:4*/
#line 13 "./tensor.cweb"
;
/*5:*/
#line 80 "./tensor.cweb"

void UTensor::decrement(IntSequence&v,int nv)
{
if(v.size()==0)
return;
int i= v.size()-1;
v[i]--;
while(i> 0&&v[i]==-1){
v[i]= nv-1;
v[--i]--;
}
}

/*:5*/
#line 14 "./tensor.cweb"
;
/*6:*/
#line 99 "./tensor.cweb"

void UTensor::increment(IntSequence&v,const IntSequence&nvmx)
{
if(v.size()==0)
return;
int i= v.size()-1;
v[i]++;
while(i> 0&&v[i]==nvmx[i]){
v[i]= 0;
v[--i]++;
}
}

/*:6*/
#line 15 "./tensor.cweb"
;
/*7:*/
#line 115 "./tensor.cweb"

void UTensor::decrement(IntSequence&v,const IntSequence&nvmx)
{
if(v.size()==0)
return;
int i= v.size()-1;
v[i]--;
while(i> 0&&v[i]==-1){
v[i]= nvmx[i]-1;
v[--i]--;
}
}

/*:7*/
#line 16 "./tensor.cweb"
;
/*8:*/
#line 131 "./tensor.cweb"

int UTensor::getOffset(const IntSequence&v,int nv)
{
int pow= 1;
int res= 0;
for(int i= v.size()-1;i>=0;i--){
res+= v[i]*pow;
pow*= nv;
}
return res;
}

/*:8*/
#line 17 "./tensor.cweb"
;
/*9:*/
#line 144 "./tensor.cweb"

int UTensor::getOffset(const IntSequence&v,const IntSequence&nvmx)
{
int pow= 1;
int res= 0;
for(int i= v.size()-1;i>=0;i--){
res+= v[i]*pow;
pow*= nvmx[i];
}
return res;
}


/*:9*/
#line 18 "./tensor.cweb"
;
/*10:*/
#line 164 "./tensor.cweb"

void FTensor::decrement(IntSequence&v,int nv)
{
int i= v.size()-1;
while(i> 0&&v[i-1]==v[i])
i--;
v[i]--;
for(int j= i+1;j<v.size();j++)
v[j]= nv-1;
}

/*:10*/
#line 19 "./tensor.cweb"
;
/*11:*/
#line 215 "./tensor.cweb"

int FTensor::getOffsetRecurse(IntSequence&v,int nv)
{
if(v.size()==0)return 0;
int prefix= v.getPrefixLength();
int m= v[0];
int k= v.size();
int s1= noverk(nv+k-1,k)-noverk(nv-m+k-1,k);
IntSequence subv(v,prefix,k);
subv.add(-m);
int s2= getOffsetRecurse(subv,nv-m);
return s1+s2;
}

/*:11*/
#line 20 "./tensor.cweb"
;

/*:1*/
