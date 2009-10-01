/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SimilarityDecomp.cpp,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SimilarityDecomp.h"
#include "SchurDecomp.h"
#include "SchurDecompEig.h"
#include "SylvException.h"

#include <dynlapack.h>

#include <cmath>

SimilarityDecomp::SimilarityDecomp(const double* d, int d_size, double log10norm)
{
	SchurDecomp sd(SqSylvMatrix(d, d_size));
	q = new SqSylvMatrix(sd.getQ());
	b = new BlockDiagonal(sd.getT());
	invq = new SqSylvMatrix(d_size);
	invq->setUnit();
	invq->multLeftTrans(sd.getQ());
	double norm = pow(10.0, log10norm);
	diagonalize(norm);
}

SimilarityDecomp::~SimilarityDecomp()
{
	delete invq;
	delete b;
	delete q;
}

void SimilarityDecomp::getXDim(diag_iter start, diag_iter end,
							   int &rows, int& cols) const
{
	int si = (*start).getIndex();
	int ei = (*end).getIndex();
	cols = b->numRows() - ei;
	rows = ei - si;
}

/* find solution of X for diagonal block given by start(incl.) and
 * end(excl.). If the solution cannot be found, or it is greater than
 * norm, X is not changed and flase is returned.
 */
bool SimilarityDecomp::solveX(diag_iter start, diag_iter end,
							  GeneralMatrix& X, double norm) const
{
	int si = (*start).getIndex();
	int ei = (*end).getIndex();

	SqSylvMatrix A((const GeneralMatrix&)*b, si, si, X.numRows());
	SqSylvMatrix B((const GeneralMatrix&)*b, ei, ei, X.numCols());
	GeneralMatrix C((const GeneralMatrix&)*b, si, ei, X.numRows(), X.numCols());

	lapack_int isgn = -1;
	lapack_int m = A.numRows();
	lapack_int n = B.numRows();
	double scale;
	lapack_int info;
	dtrsyl("N", "N", &isgn, &m, &n, A.base(), &m, B.base(), &n,
				  C.base(), &m, &scale, &info);
	if (info < -1)
		throw SYLV_MES_EXCEPTION("Wrong parameter to LAPACK dtrsyl.");

	if (info == 1 || scale < 1)
		return false;
	if (C.getData().getMax() > norm)
		return false;
	
	X = C;
	return true;
}

/* multiply Q and invQ with (I -X; 0 I), and (I X; 0 I). This also sets X=-X. */
void SimilarityDecomp::updateTransform(diag_iter start, diag_iter end,
									   GeneralMatrix& X)
{
	int si = (*start).getIndex();
	int ei = (*end).getIndex();
	
	SqSylvMatrix iX(q->numRows());
	iX.setUnit();
	iX.place(X, si, ei);
	invq->GeneralMatrix::multLeft(iX);

	iX.setUnit();
	X.mult(-1.0);
	iX.place(X, si, ei);
	q->multRight(iX);
}

void SimilarityDecomp::bringGuiltyBlock(diag_iter start, diag_iter& end)
{
	double av = b->getAverageDiagSize(start, end);
	diag_iter guilty = b->findClosestDiagBlock(end, b->diag_end(), av);
	SchurDecompEig sd((QuasiTriangular&)*b); // works on b including diagonal structure
	end = sd.bubbleEigen(guilty, end); // iterators are valid
	++end;
	q->multRight(sd.getQ());
	invq->multLeftTrans(sd.getQ());
}

void SimilarityDecomp::diagonalize(double norm)
{
	diag_iter start = b->diag_begin();
	diag_iter end = start;
	++end;

	while (end != b->diag_end()) {
		int xrows;
		int xcols;
		getXDim(start, end, xrows, xcols);
		GeneralMatrix X(xrows, xcols);
		if (solveX(start, end, X, norm)) {
			updateTransform(start, end, X);
			b->setZeroBlockEdge(end);
			start = end;
			++end;
		} else {
			bringGuiltyBlock(start, end); // moves with end
		}
	}
}

void SimilarityDecomp::check(SylvParams& pars, const GeneralMatrix& m) const
{
	// M - Q*B*inv(Q)
	SqSylvMatrix c(getQ(), getB());
	c.multRight(getInvQ());
	c.add(-1.0, m);
	pars.f_err1 = c.getNorm1();
	pars.f_errI = c.getNormInf();

	// I - Q*inv(Q)
	c.setUnit();
	c.mult(-1);
	c.multAndAdd(getQ(), getInvQ());
	pars.viv_err1 = c.getNorm1();
	pars.viv_errI = c.getNormInf();

	// I - inv(Q)*Q
	c.setUnit();
	c.mult(-1);
	c.multAndAdd(getInvQ(), getQ());
	pars.ivv_err1 = c.getNorm1();
	pars.ivv_errI = c.getNormInf();	
}

void SimilarityDecomp::infoToPars(SylvParams& pars) const
{
	pars.f_blocks = getB().getNumBlocks();
	pars.f_largest = getB().getLargestBlock();
	pars.f_zeros = getB().getNumZeros();
	pars.f_offdiag = getB().getNumOffdiagonal();
}

// Local Variables:
// mode:C++
// End:
