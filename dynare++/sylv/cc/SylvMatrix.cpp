/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvMatrix.cpp,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvException.h"
#include "SylvMatrix.h"

#include <dynblas.h>
#include <dynlapack.h>

#include <cstdio>
#include <cstring>
#include <cmath>

void SylvMatrix::multLeftI(const SqSylvMatrix& m)
{
	int off = rows - m.numRows();
	if (off < 0) {
		throw SYLV_MES_EXCEPTION("Wrong matrix dimensions for multLeftI.");
	}
	GeneralMatrix subtmp(*this, off, 0, m.numRows(), cols);
	subtmp.multLeft(m);
}

void SylvMatrix::multLeftITrans(const SqSylvMatrix& m)
{
	int off = rows - m.numRows();
	if (off < 0) {
		throw SYLV_MES_EXCEPTION("Wrong matrix dimensions for multLeftITrans.");
	}
	GeneralMatrix subtmp(*this, off, 0, m.numRows(), cols);
	subtmp.multLeftTrans(m);
}


void SylvMatrix::multLeft(int zero_cols, const GeneralMatrix& a, const GeneralMatrix& b)
{
	int off = a.numRows() - a.numCols();
	if (off < 0 || a.numRows() != rows || off != zero_cols ||
		rows != b.numRows() || cols != b.numCols()) {
		throw SYLV_MES_EXCEPTION("Wrong matrix dimensions for multLeft.");
	}
	// here we cannot call SylvMatrix::gemm since it would require
	// another copy of (usually big) b (we are not able to do inplace
	// submatrix of const GeneralMatrix)
	if (a.getLD() > 0 && ld > 0) {
		blas_int mm = a.numRows();
		blas_int nn = cols;
		blas_int kk = a.numCols();
		double alpha = 1.0;
		blas_int lda = a.getLD();
		blas_int ldb = ld;
		double beta = 0.0;
		blas_int ldc = ld;
		dgemm("N", "N", &mm, &nn, &kk, &alpha, a.getData().base(), &lda,
				   b.getData().base()+off, &ldb, &beta, data.base(), &ldc);
	}
}


void SylvMatrix::multRightKron(const SqSylvMatrix& m, int order)
{
	if (power(m.numRows(), order) != cols) {
		throw SYLV_MES_EXCEPTION("Wrong number of cols for right kron multiply.");
	}
	KronVector auxrow(m.numRows(), m.numRows(), order-1);
	for (int i = 0; i < rows; i++) {
		Vector rowi(data.base()+i, rows, cols);
		KronVector rowikron(rowi, m.numRows(), m.numRows(), order-1);
		auxrow = rowi; // copy data
		m.multVecKronTrans(rowikron, auxrow);
	}
}

void SylvMatrix::multRightKronTrans(const SqSylvMatrix& m, int order)
{
	if (power(m.numRows(), order) != cols) {
		throw SYLV_MES_EXCEPTION("Wrong number of cols for right kron multiply.");
	}
	
	KronVector auxrow(m.numRows(), m.numRows(), order-1);
	for (int i = 0; i < rows; i++) {
		Vector rowi(data.base()+i, rows, cols);
		KronVector rowikron(rowi, m.numRows(), m.numRows(), order-1);
		auxrow = rowi; // copy data
		m.multVecKron(rowikron, auxrow);
	}
}

void SylvMatrix::eliminateLeft(int row, int col, Vector& x)
{
	double d = get(col, col);
	double e = get(row, col);
	if (std::abs(d) > std::abs(e)) {
		get(row, col) = 0.0;
		double mult = e/d;
		for (int i = col + 1; i < numCols(); i++) {
			get(row, i) = get(row, i) - mult*get(col, i);
		}
		x[row] = x[row] - mult*x[col];
	} else if (std::abs(e) > std::abs(d)) {
		get(row, col) = 0.0;
		get(col, col) = e;
		double mult = d/e;
		for (int i = col + 1; i < numCols(); i++) {
			double tx = get(col, i);
			double ty = get(row, i);
			get(col, i) = ty;
			get(row, i) = tx - mult*ty;
		}
		double tx = x[col];
		double ty = x[row];
		x[col] = ty;
		x[row] = tx - mult*ty;
	}
}

void SylvMatrix::eliminateRight(int row, int col, Vector& x)
{
	double d = get(row, row);
	double e = get(row, col);
	
	if (std::abs(d) > std::abs(e)) {
		get(row, col) = 0.0;
		double mult = e/d;
		for (int i = 0; i < row; i++) {
			get(i, col) = get(i, col) - mult*get(i, row);
		}
		x[col] = x[col] - mult*x[row];
	} else if (std::abs(e) > std::abs(d)) {
		get(row, col) = 0.0;
		get(row, row) = e;
		double mult = d/e;
		for (int i = 0; i < row; i++) {
			double tx = get(i, row);
			double ty = get(i, col);
			get(i, row) = ty;
			get(i, col) = tx - mult*ty;
		}
		double tx = x[row];
		double ty = x[col];
		x[row] = ty;
		x[col] = tx - mult*ty;
	}
}



SqSylvMatrix::SqSylvMatrix(const GeneralMatrix& a, const GeneralMatrix& b)
	: SylvMatrix(a,b)
{
	if (rows != cols)
		throw SYLV_MES_EXCEPTION("Wrong matrix dimensions in multiplication constructor of square matrix.");
}

void SqSylvMatrix::multVecKron(KronVector& x, const KronVector& d) const
{
	x.zeros();
	if (d.getDepth() == 0) {
		multaVec(x, d);
	} else {
		KronVector aux(x.getM(), x.getN(), x.getDepth());
		for (int i = 0; i < x.getM(); i++) {
			KronVector auxi(aux, i);
			ConstKronVector di(d, i);
			multVecKron(auxi, di);
		}
		for (int i = 0; i < rows; i++) {
			KronVector xi(x, i);
			for (int j = 0; j < cols; j++) {
				KronVector auxj(aux, j);
				xi.add(get(i,j),auxj);
			}
		}
	}
}


void SqSylvMatrix::multVecKronTrans(KronVector& x, const KronVector& d) const
{
	x.zeros();
	if (d.getDepth() == 0) {
		multaVecTrans(x, d);
	} else {
		KronVector aux(x.getM(), x.getN(), x.getDepth());
		for (int i = 0; i < x.getM(); i++) {
			KronVector auxi(aux, i);
			ConstKronVector di(d, i);
			multVecKronTrans(auxi, di);
		}
		for (int i = 0; i < rows; i++) {
			KronVector xi(x, i);
			for (int j = 0; j < cols; j++) {
				KronVector auxj(aux, j);
				xi.add(get(j,i), auxj);
			}
		}
	}
}

void SqSylvMatrix::multInvLeft2(GeneralMatrix& a, GeneralMatrix& b,
								double& rcond1, double& rcondinf) const
{
	if (rows != a.numRows() || rows != b.numRows()) {
		throw SYLV_MES_EXCEPTION("Wrong dimensions for multInvLeft2.");
	}
	// PLU factorization
	Vector inv(data);
	lapack_int * const ipiv = new lapack_int[rows];
	lapack_int info;
	lapack_int rows2 = rows;
	dgetrf(&rows2, &rows2, inv.base(), &rows2, ipiv, &info);
	// solve a
	lapack_int acols = a.numCols();
	double* abase = a.base();
	dgetrs("N", &rows2, &acols, inv.base(), &rows2, ipiv,
				  abase, &rows2, &info);
	// solve b
	lapack_int bcols = b.numCols();
	double* bbase = b.base();
	dgetrs("N", &rows2, &bcols, inv.base(), &rows2, ipiv,
				  bbase, &rows2, &info);
	delete [] ipiv;

	// condition numbers
	double* const work = new double[4*rows];
	lapack_int* const iwork = new lapack_int[rows];
	double norm1 = getNorm1();
	dgecon("1", &rows2, inv.base(), &rows2, &norm1, &rcond1, 
				  work, iwork, &info);
	double norminf = getNormInf();
	dgecon("I", &rows2, inv.base(), &rows2, &norminf, &rcondinf, 
				  work, iwork, &info);
	delete [] iwork;
	delete [] work;
}

void SqSylvMatrix::setUnit()
{
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (i==j)
				get(i,j) = 1.0;
			else 
				get(i,j) = 0.0;
		}
	}
}

// Local Variables:
// mode:C++
// End:
