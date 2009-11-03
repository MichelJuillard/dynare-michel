/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/QuasiTriangularZero.cpp,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#include "QuasiTriangularZero.h"
#include "SchurDecomp.h"
#include "SylvMatrix.h"
#include "SylvException.h"

#include <cstdio>

QuasiTriangularZero::QuasiTriangularZero(int num_zeros, const double* d,
										 int d_size)
	: QuasiTriangular(SqSylvMatrix(GeneralMatrix(d, num_zeros+d_size, d_size),
								   num_zeros, 0, d_size).getData().base(),
					  d_size),
	  nz(num_zeros),
	  ru(GeneralMatrix(d, num_zeros+d_size, d_size), 0, 0, num_zeros, d_size)
{
}

QuasiTriangularZero::QuasiTriangularZero(double r,
										 const QuasiTriangularZero& t)
	: QuasiTriangular(r, t),
	  nz(t.nz),
	  ru(t.ru)
{
	ru.mult(r);
}

QuasiTriangularZero::QuasiTriangularZero(double r,
										 const QuasiTriangularZero& t,
										 double rr,
										 const QuasiTriangularZero& tt)
	: QuasiTriangular(r, t, rr, tt),
	  nz(t.nz),
	  ru(t.ru)
{
	ru.mult(r);
	ru.add(rr, tt.ru);
}

QuasiTriangularZero::QuasiTriangularZero(int p, const QuasiTriangularZero& t)
	: QuasiTriangular(p, t),
	  nz(t.nz),
	  ru(t.ru)
{
	ru.multRight(t);
}

QuasiTriangularZero::QuasiTriangularZero(const SchurDecompZero& decomp)
	: QuasiTriangular(decomp.getT().getData().base(),
					  decomp.getT().numRows()),
	  nz(decomp.getZeroCols()),
	  ru(decomp.getRU())
{
}

QuasiTriangularZero::QuasiTriangularZero(const QuasiTriangular& t)
	: QuasiTriangular(t),
	  nz(0), ru(0, t.getDiagonal().getSize())
{
}

QuasiTriangularZero::~QuasiTriangularZero()
{
}

void QuasiTriangularZero::solvePre(Vector& x, double& eig_min)
{
	Vector xu(x, 0, nz);
	Vector xl(x, nz, x.length()-nz);
	QuasiTriangular::solvePre(xl, eig_min);
	ru.multsVec(xu, xl);
	if (nz > 0)
		eig_min = (eig_min > 1.0)? 1.0 : eig_min;
}

void QuasiTriangularZero::solvePreTrans(Vector& x, double& eig_min)
{
	Vector xu(x, 0, nz);
	Vector xl(x, nz, x.length()-nz);
	ru.multsVecTrans(xl, xu);
	QuasiTriangular::solvePreTrans(xl, eig_min);
	if (nz > 0)
		eig_min = (eig_min > 1.0)? 1.0 : eig_min;
}

void QuasiTriangularZero::multVec(Vector& x, const ConstVector& b) const
{
	x.zeros();
	multaVec(x, b);
}

void QuasiTriangularZero::multVecTrans(Vector& x, const ConstVector& b) const
{
	x.zeros();
	multaVecTrans(x, b);
}

void QuasiTriangularZero::multaVec(Vector& x, const ConstVector& b) const
{
	ConstVector bl(b, nz, b.length()-nz);
	Vector xu(x, 0, nz);
	Vector xl(x, nz, x.length()-nz);
	xu.zeros();
	ru.multaVec(xu, bl);
	QuasiTriangular::multVec(xl, bl);
}

void QuasiTriangularZero::multaVecTrans(Vector& x, const ConstVector& b) const
{
	ConstVector bu(b, 0, b.length());
	ConstVector bl(b, nz, b.length()-nz);
	Vector xu(x, 0, nz);
	Vector xl(x, nz, x.length()-nz);
	xu.zeros();
	QuasiTriangular::multVecTrans(xl, bl);
	ru.multaVecTrans(xl, bu);
}

void QuasiTriangularZero::multLeftOther(GeneralMatrix& a) const
{
	GeneralMatrix a1(a, 0, 0, nz, a.numCols());
	GeneralMatrix a2(a, nz, 0, a.numRows()-nz, a.numCols());
	a1.mult(ru, a2);
	QuasiTriangular::multLeftOther(a2);
}

void QuasiTriangularZero::print() const
{
	printf("super=\n");
	QuasiTriangular::print();
	printf("nz=%d\n",nz);
	printf("ru=\n");
	ru.print();
}

void QuasiTriangularZero::multKron(KronVector& x) const
{
	throw SYLV_MES_EXCEPTION("Attempt to run QuasiTriangularZero::multKron.");
}

void QuasiTriangularZero::multKronTrans(KronVector& x) const
{
	throw SYLV_MES_EXCEPTION("Attempt to run QuasiTriangularZero::multKronTrans.");
}

