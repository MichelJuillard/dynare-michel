/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SchurDecomp.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYM_SCHUR_DECOMP_H
#define SYM_SCHUR_DECOMP_H

#include "SylvMatrix.h"

class SymSchurDecomp {
protected:
	Vector lambda;
	SqSylvMatrix q;
public:
	/** Calculates A = Q*Lambda*Q^T, where A is assummed to be
	 * symmetric and Lambda real diagonal, hence a vector. */
	SymSchurDecomp(const GeneralMatrix& a);
	SymSchurDecomp(const SymSchurDecomp& ssd)
		: lambda(ssd.lambda), q(ssd.q) {}
	virtual ~SymSchurDecomp() {}
	const Vector& getLambda() const
		{return lambda;}
	const SqSylvMatrix& getQ() const
		{return q;}
	/** Return factor F*F^T = A, raises and exception if A is not
	 * positive semidefinite, F must be square. */
	void getFactor(GeneralMatrix& f) const;
	/** Returns true if A is positive semidefinite. */
	bool isPositiveSemidefinite() const;
	/** Correct definitness. This sets all eigenvalues between minus
	 * tolerance and zero to zero. */
	void correctDefinitness(double tol);

};

#endif


// Local Variables:
// mode:C++
// End:
