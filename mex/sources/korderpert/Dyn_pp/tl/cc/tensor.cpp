/*1:*/

#include "tensor.h"
#include "tl_exception.h"
#include "tl_static.h"

/*2:*/

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
;
/*3:*/

int Tensor::noverseq_ip(IntSequence&s)
{
	if(s.size()==0||s.size()==1)
		return 1;
	s[1]+= s[0];
	return noverk(s[1],s[0])*noverseq(IntSequence(s,1,s.size()));
}

/*:3*/
;
/*4:*/

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
;
/*5:*/

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
;
/*6:*/

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
;
/*7:*/

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
;
/*8:*/

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
;
/*9:*/

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
;
/*10:*/

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
;
/*11:*/

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
;

/*:1*/
