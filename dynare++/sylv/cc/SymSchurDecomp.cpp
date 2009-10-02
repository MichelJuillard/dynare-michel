/* $Header$ */

/* Tag $Name$ */

#include "SymSchurDecomp.h"
#include "SylvException.h"

#include <dynlapack.h>

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
	lapack_int n = mata.numRows();
	GeneralMatrix tmpa(mata);
	double* a = tmpa.base();
	lapack_int lda = tmpa.getLD();
	double dum;
	double* vl = &dum;
	double* vu = &dum;
	lapack_int idum;
	lapack_int* il = &idum;
	lapack_int* iu = &idum;
	double abstol = 0.0;
	lapack_int m = n;
	double* w = lambda.base();
	double* z = q.base();
	lapack_int ldz = q.getLD();
	lapack_int* isuppz = new lapack_int[2*std::max(1,(int) m)];
	double tmpwork;
	lapack_int lwork = -1;
	lapack_int tmpiwork;
	lapack_int liwork = -1;
	lapack_int info;

	// query for lwork and liwork
	dsyevr(jobz, range, uplo, &n, a, &lda, vl, vu, il, iu, &abstol,
				  &m, w, z, &ldz, isuppz, &tmpwork, &lwork, &tmpiwork, &liwork, &info);
	lwork = (int)tmpwork;
	liwork = tmpiwork;
	// allocate work arrays
	double* work = new double[lwork];
	lapack_int* iwork = new lapack_int[liwork];
	
	// do the calculation
	dsyevr(jobz, range, uplo, &n, a, &lda, vl, vu, il, iu, &abstol,
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
