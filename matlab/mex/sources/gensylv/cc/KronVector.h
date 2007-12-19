/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/KronVector.h,v 1.1.1.1 2004/06/04 13:00:31 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef KRON_VECTOR_H
#define KRON_VECTOR_H

#include "Vector.h"

class ConstKronVector;

class KronVector : public Vector {
protected:
	int m;
	int n;
	int depth;
public:
	KronVector() : Vector((double*)0, 0), m(0), n(0), depth(0)  {}
	KronVector(int mm, int nn, int dp); // new instance
	KronVector(Vector& v, int mm, int nn, int dp); // conversion
	KronVector(KronVector&, int i); // picks i-th subvector
	KronVector(const ConstKronVector& v); // new instance and copy
	const KronVector& operator=(KronVector& v)
		{Vector::operator=(v); m=v.m; n=v.n; depth = v.depth; return *this;}
	const KronVector& operator=(const KronVector& v)
		{Vector::operator=(v); m=v.m; n=v.n; depth = v.depth; return *this;}
	const KronVector& operator=(const ConstKronVector& v);
	const KronVector& operator=(const Vector& v);
	int getM() const {return m;}
	int getN() const {return n;}
	int getDepth() const {return depth;}
};

class ConstKronVector : public ConstVector
{
protected:
	int m;
	int n;
	int depth;
public:
	ConstKronVector(const KronVector& v);
	ConstKronVector(const ConstKronVector& v);
	ConstKronVector(const Vector& v, int mm, int nn, int dp);
	ConstKronVector(const ConstVector& v, int mm, int nn, int dp);
	ConstKronVector(const KronVector& v, int i);
	ConstKronVector(const ConstKronVector& v, int i);
	int getM() const {return m;}
	int getN() const {return n;}
	int getDepth() const {return depth;}
};

int power(int m, int depth);

#endif /* KRON_VECTOR */

// Local Variables:
// mode:C++
// End:
