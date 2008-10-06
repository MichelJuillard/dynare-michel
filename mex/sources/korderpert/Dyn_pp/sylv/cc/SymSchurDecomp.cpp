/* $Header$ */

/* Tag $Name$ */

#include "SymSchurDecomp.h"
#include "SylvException.h"

#include "cpplapack.h"

#include <algorithm>
#include <cmath>

SymSchurDecomp::SymSchurDecomp(const GeneralMatrix& mata)
	: lambda(mata.numRows()), q(mata.numRows())
{
	// check mata is square
	if (mata.numRows() != mata.numCols())
		throw SYLV_MES_EXCEPTION("Matrix is not square in SymSchurDecomp constructor");

	// prepare for dsyevr
	const char* jobz = "V";
	const char* range = "A";
	const char* uplo = "U";
	int n = mata.numRows();
	GeneralMatrix tmpa(mata);
	double* a = tmpa.base();
	int lda = tmpa.getLD();
	double dum;
	double* vl = &dum;
	double* vu = &dum;
	int idum;
	int* il = &idum;
	int* iu = &idum;
	double abstol = 0.0;
	int m = n;
	double* w = lambda.base();
	double* z = q.base();
	int ldz = q.getLD();
	int* isuppz = new int[2*std::max(1,m)];
	double tmpwork;
	int lwork = -1;
	int tmpiwork;
	int liwork = -1;
	int info;

	// query for lwork and liwork
	LAPACK_dsyevr(jobz, range, uplo, &n, a, &lda, vl, vu, il, iu, &abstol,
				  &m, w, z, &ldz, isuppz, &tmpwork, &lwork, &tmpiwork, &liwork, &info);
	lwork = (int)tmpwork;
	liwork = tmpiwork;
	// allocate work arrays
	double* work = new double[lwork];
	int* iwork = new int[liwork];
	
	// do the calculation
	LAPACK_dsyevr(jobz, range, uplo, &n, a, &lda, vl, vu, il, iu, &abstol,
				  &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

	if (info < 0)
		throw SYLV_MES_EXCEPTION("Internal error in SymSchurDecomp constructor");
	if (info > 0)
		throw SYLV_MES_EXCEPTION("Internal LAPACK error in DSYEVR");

	delete [] work;
	delete [] iwork;
	delete [] isuppz;
}

void SymSchurDecomp::getFactor(GeneralMatrix& f) const
{
	if (f.numRows() != q.numRows())
		throw SYLV_MES_EXCEPTION("Wrong dimension of factor matrix in SymSchurDecomp::getFactor");
	if (f.numRows() != f.numCols())
		throw SYLV_MES_EXCEPTION("Factor matrix is not square in SymSchurDecomp::getFactor");
	if (! isPositiveSemidefinite())
		throw SYLV_MES_EXCEPTION("Symmetric decomposition not positive semidefinite in SymSchurDecomp::getFactor");

	f = q;
	for (int i = 0; i < f.numCols(); i++) {
		Vector fi(f, i);
		fi.mult(std::sqrt(lambda[i]));
	}
}


// LAPACK says that eigenvalues are ordered in ascending order, but we
// do not rely her on it
bool SymSchurDecomp::isPositiveSemidefinite() const
{
	for (int i = 0; i < lambda.length(); i++)
		if (lambda[i] < 0)
			return false;
	return true;
}

void SymSchurDecomp::correctDefinitness(double tol)
{
	for (int i = 0; i < lambda.length(); i++)
		if (lambda[i] < 0 && lambda[i] > - tol)
			lambda[i] = 0.0;
}
