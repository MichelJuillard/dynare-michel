/*1:*/

#include "tl_static.h"
#include "tl_exception.h"

TLStatic tls;
/*2:*/

TLStatic::TLStatic()
{
	ebundle= NULL;
	pbundle= NULL;
	ptriang= NULL;
}

TLStatic::~TLStatic()
{
	if(ebundle)
		delete ebundle;
	if(pbundle)
		delete pbundle;
	if(ptriang)
		delete ptriang;
}

void TLStatic::init(int dim,int nvar)
{
	if(ebundle)
		ebundle->generateUpTo(dim);
	else
		ebundle= new EquivalenceBundle(dim);
	
	if(pbundle)
		pbundle->generateUpTo(dim);
	else
		pbundle= new PermutationBundle(dim);
	
	if(ptriang)
		delete ptriang;
	ptriang= new PascalTriangle(nvar,dim);
}

/*:2*/
;
/*3:*/

PascalTriangle::PascalTriangle(int n,int k)
:data(new int[(n+1)*(k+1)]),kmax(k),nmax(n)
{
	for(int i= 0;i<=n;i++)
		data[i]= 1;
	for(int j= 1;j<=k;j++){
		data[j*(nmax+1)]= 1;
		for(int i= 1;i<=n;i++)
			data[j*(nmax+1)+i]= noverk(i+j-1,j)+noverk(i+j-1,j-1);
	}
}

/*:3*/
;
/*4:*/

int PascalTriangle::noverk(int n,int k)const
{
	TL_RAISE_IF(k> n||n<0,
		"Wrong arguments for PascalTriangle::noverk");
	
	if(k<=kmax&&n-k<=nmax)
		return data[k*(nmax+1)+n-k];
	
	if(n-k<=kmax&&k<=nmax)
		return data[(n-k)*(nmax+1)+k];
	
	TL_RAISE("n or k out of range in PascalTriangle::noverk");
	return 0;
}

/*:4*/
;

/*:1*/
