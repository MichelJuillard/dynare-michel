/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/KronVector.cpp,v 1.1.1.1 2004/06/04 13:00:31 kamenik Exp $ */

/* Tag $Name:  $ */

#include "KronVector.h"
#include "SylvException.h"

int power(int m, int depth)
{
	int p = 1;
	for (int i = 0; i < depth; i++) {
		p *= m;
	}
	return p;
}

KronVector::KronVector(int mm, int nn, int dp)
	: Vector(power(mm,dp)*nn), m(mm), n(nn), depth(dp)
{}

KronVector::KronVector(Vector& v, int mm, int nn, int dp)
	: Vector(v), m(mm), n(nn), depth(dp)
{
	len = power(m,depth)*n;
	if (v.length() != length()) {
		throw SYLV_MES_EXCEPTION("Bad conversion KronVector from Vector.");
	}
}

KronVector::KronVector(KronVector& v, int i)
	: Vector(v, i*power(v.m,v.depth-1)*v.n, power(v.m, v.depth-1)*v.n), m(v.m), n(v.n),
	  depth(v.depth-1)
{
	if (depth < 0) {
		throw SYLV_MES_EXCEPTION("Bad KronVector pick, depth < 0.");
	}
}

KronVector::KronVector(const ConstKronVector& v)
	: Vector(v.length()), m(v.getM()), n(v.getN()), depth(v.getDepth())
{
	Vector::operator=(v);
}

const KronVector& KronVector::operator=(const ConstKronVector& v)
{
	Vector::operator=(v);
	m=v.getM();
	n=v.getN();
	depth = v.getDepth();
	return *this;
}

const KronVector& KronVector::operator=(const Vector& v)
{
	if (length() != v.length()) {
		throw SYLV_MES_EXCEPTION("Wrong lengths for vector operator =.");
	}
	Vector::operator=(v);
	return *this;
}



ConstKronVector::ConstKronVector(const KronVector& v)
	: ConstVector(v), m(v.getM()), n(v.getN()), depth(v.getDepth())
{}

ConstKronVector::ConstKronVector(const ConstKronVector& v)
	: ConstVector(power(v.getM(),v.getDepth())*v.getN()), m(v.getM()), n(v.getN()),
	  depth(v.getDepth())	  
{}

ConstKronVector::ConstKronVector(const Vector& v, int mm, int nn, int dp)
	: ConstVector(v), m(mm), n(nn), depth(dp)
{
	len = power(m,depth)*n;
	if (v.length() != length()) {
		throw SYLV_MES_EXCEPTION("Bad conversion KronVector from Vector.");
	}
}

ConstKronVector::ConstKronVector(const ConstVector& v, int mm, int nn, int dp)
	: ConstVector(v), m(mm), n(nn), depth(dp)
{
	len = power(m,depth)*n;
	if (v.length() != length()) {
		throw SYLV_MES_EXCEPTION("Bad conversion KronVector from Vector.");
	}
}

ConstKronVector::ConstKronVector(const KronVector& v, int i)
	: ConstVector(v, i*power(v.getM(),v.getDepth()-1)*v.getN(),
				  power(v.getM(),v.getDepth()-1)*v.getN()),
	  m(v.getM()), n(v.getN()), depth(v.getDepth()-1)
{
	if (depth < 0) {
		throw SYLV_MES_EXCEPTION("Bad KronVector pick, depth < 0.");
	}
}

ConstKronVector::ConstKronVector(const ConstKronVector& v, int i)
	: ConstVector(v, i*power(v.m,v.depth-1)*v.n, power(v.m,v.depth-1)*v.n),
	  m(v.getM()), n(v.getN()), depth(v.getDepth()-1)
{
	if (depth < 0) {
		throw SYLV_MES_EXCEPTION("Bad KronVector pick, depth < 0.");
	}
}
