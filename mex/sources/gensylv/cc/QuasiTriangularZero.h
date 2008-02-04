/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/QuasiTriangularZero.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef QUASI_TRIANGULAR_ZERO_H
#define QUASI_TRIANGULAR_ZERO_H

#include "QuasiTriangular.h"
#include "GeneralMatrix.h"

class QuasiTriangularZero : public QuasiTriangular {
	int nz; // number of zero columns
	GeneralMatrix ru; // data in right upper part (nz,d_size)
public:
	QuasiTriangularZero(int num_zeros, const double* d, int d_size);
	QuasiTriangularZero(double r, const QuasiTriangularZero& t);
	QuasiTriangularZero(double r, const QuasiTriangularZero& t,
						double rr, const QuasiTriangularZero& tt);
	QuasiTriangularZero(int p, const QuasiTriangularZero& t);
	QuasiTriangularZero(const QuasiTriangular& t);
	QuasiTriangularZero(const SchurDecompZero& decomp);
	~QuasiTriangularZero();
	void solvePre(Vector& x, double& eig_min);
	void solvePreTrans(Vector& x, double& eig_min);
	void multVec(Vector& x, const ConstVector& b) const;
	void multVecTrans(Vector& x, const ConstVector& b) const;
	void multaVec(Vector& x, const ConstVector& b) const;
	void multaVecTrans(Vector& x, const ConstVector& b) const;
	void multKron(KronVector& x) const;
	void multKronTrans(KronVector& x) const;
	void multLeftOther(GeneralMatrix& a) const;
	/* clone */
	virtual QuasiTriangular* clone() const
		{return new QuasiTriangularZero(*this);}
	virtual QuasiTriangular* clone(int p, const QuasiTriangular& t) const
		{return new QuasiTriangularZero(p, (const QuasiTriangularZero&)t);}
	virtual QuasiTriangular* clone(double r) const
		{return new QuasiTriangularZero(r, *this);}
	virtual QuasiTriangular* clone(double r, double rr, const QuasiTriangular& tt) const
		{return new QuasiTriangularZero(r, *this, rr, (const QuasiTriangularZero&)tt);}
	void print() const;
};

#endif /* QUASI_TRIANGULAR_ZERO_H */

// Local Variables:
// mode:C++
// End:
