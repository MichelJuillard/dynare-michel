/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SimilarityDecomp.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SIMILARITY_DECOMP_H
#define SIMILARITY_DECOMP_H

#include "SylvMatrix.h"
#include "BlockDiagonal.h"
#include "SylvParams.h"

class SimilarityDecomp {
	SqSylvMatrix* q;
	BlockDiagonal* b;
	SqSylvMatrix* invq;
	typedef BlockDiagonal::diag_iter diag_iter;
public:
	SimilarityDecomp(const double* d, int d_size, double log10norm = 3.0);
	virtual ~SimilarityDecomp();
	const SqSylvMatrix& getQ() const
		{return *q;}
	const SqSylvMatrix& getInvQ() const
		{return *invq;}
	const BlockDiagonal& getB() const
		{return *b;}
	void check(SylvParams& pars, const GeneralMatrix& m) const;
	void infoToPars(SylvParams& pars) const;
protected:
	void getXDim(diag_iter start, diag_iter end, int& rows, int& cols) const;
	bool solveX(diag_iter start, diag_iter end, GeneralMatrix& X, double norm) const;
	void updateTransform(diag_iter start, diag_iter end, GeneralMatrix& X);
	void bringGuiltyBlock(diag_iter start, diag_iter& end);
	void diagonalize(double norm);
};

#endif /* SIMILARITY_DECOMP_H */


// Local Variables:
// mode:C++
// End:
