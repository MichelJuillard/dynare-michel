/*1:*/
#line 19 "./tl_static.hweb"

#ifndef TL_STATIC_H
#define TL_STATIC_H

#include "equivalence.h"
#include "permutation.h"

/*2:*/
#line 36 "./tl_static.hweb"

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
#line 26 "./tl_static.hweb"
;
/*3:*/
#line 50 "./tl_static.hweb"

struct TLStatic{
EquivalenceBundle*ebundle;
PermutationBundle*pbundle;
PascalTriangle*ptriang;

TLStatic();
~TLStatic();
void init(int dim,int nvar);
};


/*:3*/
#line 27 "./tl_static.hweb"
;
extern TLStatic tls;

#endif

/*:1*/
