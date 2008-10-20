/*1:*/

#ifndef TL_STATIC_H
#define TL_STATIC_H

#include "equivalence.h"
#include "permutation.h"

/*2:*/

class PascalTriangle{
	int*data;
	int kmax;
	int nmax;
public:
	PascalTriangle(int n,int k);
	~PascalTriangle()
	{delete[]data;}
	int noverk(int n,int k)const;
};


/*:2*/
;
/*3:*/

struct TLStatic{
	EquivalenceBundle*ebundle;
	PermutationBundle*pbundle;
	PascalTriangle*ptriang;
	
	TLStatic();
	~TLStatic();
	void init(int dim,int nvar);
};


/*:3*/
;
extern TLStatic tls;

#endif

/*:1*/
