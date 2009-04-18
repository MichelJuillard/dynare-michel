/*
 * Copyright (C) 2003-2005 Ondra Kamenik
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/QuasiTriangular.cpp,v 1.1.1.1 2004/06/04 13:00:31 kamenik Exp $ */

/* Tag $Name:  $ */

#include "QuasiTriangular.h"
#include "SylvException.h"
#include "SchurDecomp.h"

#include "cppblas.h"

#include <stdio.h>
#include <cmath>

using namespace std;

double DiagonalBlock::getDeterminant() const
{
	return (*alpha)*(*alpha) + getSBeta();
}

double DiagonalBlock::getSBeta() const
{
	return -(*beta1)*(*beta2);
}

double DiagonalBlock::getSize() const
{
	if (real)
		return abs(*alpha);
	else
		return sqrt(getDeterminant());
}

// this function makes Diagonal inconsistent, it should only be used
// on temorary matrices, which will not be used any more, e.g. in
// QuasiTriangular::solve (we need fast performance)
void DiagonalBlock::setReal()
{
	*beta1 = 0;
	*beta2 = 0;
	real = true;
}

void DiagonalBlock::checkBlock(const double* d, int d_size)
{
	const double* a1 = d + jbar*d_size+jbar;
	const double* b1 = a1 + d_size;
	const double* b2 = a1 + 1;
	const double* a2 = b1 + 1;
	if (a1 != alpha.a1)
		throw SYLV_MES_EXCEPTION("Bad alpha1.");
	if (!real && b1 != beta1)
		throw SYLV_MES_EXCEPTION("Bad beta1.");
	if (!real && b2 != beta2)
		throw SYLV_MES_EXCEPTION("Bad beta2.");
	if (!real && a2 != alpha.a2)
		throw SYLV_MES_EXCEPTION("Bad alpha2.");
}

Diagonal::Diagonal(double* data, int d_size)
{
	int nc = getNumComplex(data, d_size); // return nc <= d_size/2
	num_all = d_size - nc;
	num_real = d_size - 2*nc;

	int jbar = 0;
	int j = 0;
	while (j < num_all) {		
		int id = jbar*d_size + jbar; // index of diagonal block in data
		int ill = id + 1; // index of element below the diagonal
		int iur = id + d_size; // index of element right to diagonal
		int idd = id + d_size + 1; // index of element next on diagonal
		if ((jbar < d_size-1) && !isZero(data[ill])) {
			// it is not last column and we have nonzero below diagonal
			DiagonalBlock b(jbar, false, &data[id], &data[idd],
							&data[iur], &data[ill]);
			blocks.push_back(b);
			jbar++;
		} else {
			// it is last column or we have zero below diagonal
			DiagonalBlock b(jbar, true, &data[id], &data[id], NULL, NULL);
			blocks.push_back(b);
		}
		jbar++;
		j++;
	}
}


Diagonal::Diagonal(double* data, const Diagonal& d)
{
	num_all = d.num_all;
	num_real = d.num_real;
	int d_size = d.getSize();
	for (const_diag_iter it = d.begin(); it != d.end(); ++it) {
		const DiagonalBlock& dit = *it;
		double* beta1 = NULL;
		double* beta2 = NULL;
		int id = dit.getIndex()*(d_size+1);
		int idd = id;
		if (! dit.isReal()) {
			beta1 = &data[id+d_size];
			beta2 = &data[id+1];
			idd = id + d_size + 1;
		}
		DiagonalBlock b(dit.getIndex(), dit.isReal(),
						&data[id], &data[idd], beta1, beta2);
		blocks.push_back(b);
	}
}


void Diagonal::copy(const Diagonal& d)
{
	num_all = d.num_all;
	num_real = d.num_real;
	blocks = d.blocks;
}

int Diagonal::getNumComplex(const double* data, int d_size)
{
	int num_complex = 0;
	int in = 1;
	for (int i = 0; i < d_size-1; i++, in = in + d_size + 1) {
		if (! isZero(data[in])) {
			num_complex++;
			if (in < d_size - 2 && ! isZero(data[in + d_size +1])) {
				throw SYLV_MES_EXCEPTION("Matrix is not quasi-triangular");
			}
		}
	}
	return num_complex;
}

void Diagonal::changeBase(double* p)
{
	int d_size = getSize();
	for (diag_iter it = begin(); it != end(); ++it) {
		const DiagonalBlock& b = *it;
		int jbar = b.getIndex();
		int base = d_size*jbar + jbar;
		if (b.isReal()) {
			DiagonalBlock bnew(jbar, true, &p[base], &p[base],
							   NULL, NULL);
			*it = bnew;
		} else {
			DiagonalBlock bnew(jbar, false, &p[base], &p[base+d_size+1],
							   &p[base+d_size], &p[base+1]);
			*it = bnew;			
		}
	}
}

void Diagonal::getEigenValues(Vector& eig) const
{
	int d_size = getSize();
	if (eig.length() != 2*d_size) {
		char mes[500];
		sprintf(mes, "Wrong length of vector for eigenvalues len=%d, should be=%d.\n",
				eig.length(), 2*d_size);
		throw SYLV_MES_EXCEPTION(mes);
	}
	for (const_diag_iter it = begin(); it != end(); ++it) {
		const DiagonalBlock& b = *it;
		int ind = b.getIndex();
		eig[2*ind] = *(b.getAlpha());
		if (b.isReal()) {
			eig[2*ind+1] = 0.0;
		} else {
			double beta = sqrt(b.getSBeta());
			eig[2*ind+1] = beta;
			eig[2*ind+2] = eig[2*ind];
			eig[2*ind+3] = -beta;
		}
	}
}

/* swaps logically blocks 'it', and '++it'. remember to move also
 * addresses, alpha, beta1, beta2. This is a dirty (but most
 * effective) way how to do it. */
void Diagonal::swapLogically(diag_iter it)
{
	diag_iter itp = it;
	++itp;

	if ((*it).isReal() && !(*itp).isReal()) {
		// first is real, second is complex
		double* d1 = (*it).alpha.a1;
		double* d2 = (*itp).alpha.a1;
		double* d3 = (*itp).alpha.a2;
		// swap
		DiagonalBlock new_it((*it).jbar, d1, d2);
		*it = new_it;
		DiagonalBlock new_itp((*itp).jbar+1, d3);
		*itp = new_itp;
	} else if (!(*it).isReal() && (*itp).isReal()) {
		// first is complex, second is real
		double* d1 = (*it).alpha.a1;
		double* d2 = (*it).alpha.a2;
		double* d3 = (*itp).alpha.a1;
		// swap
		DiagonalBlock new_it((*it).jbar, d1);
		*it = new_it;
		DiagonalBlock new_itp((*itp).jbar-1, d2, d3);
		*itp = new_itp;
	}
}

void Diagonal::checkConsistency(diag_iter it)
{
	if (!(*it).isReal() && isZero((*it).getBeta2())) {
		(*it).getBeta2() = 0.0; // put exact zero
		int jbar = (*it).getIndex();
		double* d2 = (*it).alpha.a2;
		(*it).alpha.a2 = (*it).alpha.a1;
		(*it).real = true;
		(*it).beta1 = 0;
		(*it).beta2 = 0;
		DiagonalBlock b(jbar+1, d2);
		blocks.insert((++it).iter(), b);
		num_real += 2;
		num_all++;
	}		
}

double Diagonal::getAverageSize(diag_iter start, diag_iter end)
{
	double res = 0;
	int num = 0;
	for (diag_iter run = start; run != end; ++run) {
		num++;
		res += (*run).getSize();
	}
	if (num > 0)
		res = res/num;
	return res;
}

Diagonal::diag_iter Diagonal::findClosestBlock(diag_iter start, diag_iter end, double a)
{
	diag_iter closest = start;
	double minim = 1.0e100;
	for (diag_iter run = start; run != end; ++run) {
		double dist = abs(a - (*run).getSize());
		if (dist < minim) {
			minim = dist;
			closest = run;
		}
	}
	return closest;
}

Diagonal::diag_iter Diagonal::findNextLargerBlock(diag_iter start, diag_iter end, double a)
{
	diag_iter closest = start;
	double minim = 1.0e100;
	for (diag_iter run = start; run != end; ++run) {
		double dist = (*run).getSize() - a;
		if ((0 <= dist) && (dist < minim)) {
			minim = dist;
			closest = run;
		}
	}
	return closest;
}	

void Diagonal::print() const
{
	printf("Num real: %d, num complex: %d\n",getNumReal(), getNumComplex());
	for (const_diag_iter it = begin(); it != end(); ++it) {
		if ((*it).isReal()) {
			printf("real: jbar=%d, alpha=%f\n", (*it).getIndex(), *((*it).getAlpha()));
		}
		else {
			printf("complex: jbar=%d, alpha=%f, beta1=%f, beta2=%f\n",
				   (*it).getIndex(), *((*it).getAlpha()), (*it).getBeta1(), (*it).getBeta2());
		}
	}
}

double Diagonal::EPS = 1.0e-300;

bool Diagonal::isZero(double p)
{
	return (abs(p)<EPS);
}


QuasiTriangular::const_col_iter
QuasiTriangular::col_begin(const DiagonalBlock& b) const
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return const_col_iter(&getData()[jbar*d_size], d_size, b.isReal(), 0);
}

QuasiTriangular::col_iter
QuasiTriangular::col_begin(const DiagonalBlock& b)
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return col_iter(&getData()[jbar*d_size], d_size, b.isReal(), 0);
}

QuasiTriangular::const_row_iter
QuasiTriangular::row_begin(const DiagonalBlock& b) const
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	int off = jbar*d_size+jbar+d_size;
	int col = jbar+1;
	if (!b.isReal()) {
		off = off + d_size;
		col++;
	}
	return const_row_iter(&getData()[off], d_size, b.isReal(), col);
}


QuasiTriangular::row_iter
QuasiTriangular::row_begin(const DiagonalBlock& b)
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	int off = jbar*d_size+jbar+d_size;
	int col = jbar+1;
	if (!b.isReal()) {
		off = off + d_size;
		col++;
	}
	return row_iter(&getData()[off], d_size, b.isReal(), col);
}

QuasiTriangular::const_col_iter
QuasiTriangular::col_end(const DiagonalBlock& b) const
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return const_col_iter(getData().base()+jbar*d_size+jbar, d_size, b.isReal(),
						  jbar);
}

QuasiTriangular::col_iter
QuasiTriangular::col_end(const DiagonalBlock& b)
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return col_iter(&getData()[jbar*d_size+jbar], d_size, b.isReal(), jbar);
}

QuasiTriangular::const_row_iter
QuasiTriangular::row_end(const DiagonalBlock& b) const
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return const_row_iter(&getData()[d_size*d_size+jbar], d_size, b.isReal(),
						  d_size);
}

QuasiTriangular::row_iter
QuasiTriangular::row_end(const DiagonalBlock& b)
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return row_iter(&getData()[d_size*d_size+jbar], d_size, b.isReal(), d_size);
}

QuasiTriangular::QuasiTriangular(double r, const QuasiTriangular& t)
	: SqSylvMatrix(t.numRows()), diagonal(getData().base(), t.diagonal)
{
	setMatrix(r, t);
}

QuasiTriangular::QuasiTriangular(double r, const QuasiTriangular& t,
								 double rr, const QuasiTriangular& tt)
	: SqSylvMatrix(t.numRows()), diagonal(getData().base(), t.diagonal)
{
	setMatrix(r, t);
	addMatrix(rr, tt);
}

QuasiTriangular::QuasiTriangular(const QuasiTriangular& t)
	: SqSylvMatrix(t), diagonal(getData().base(), t.diagonal)
{
}

QuasiTriangular::QuasiTriangular(const double* d, int d_size)
	: SqSylvMatrix(d, d_size), diagonal(getData().base(), d_size)
{}

QuasiTriangular::~QuasiTriangular()
{
}

QuasiTriangular::QuasiTriangular(int p, const QuasiTriangular& t)
	: SqSylvMatrix(t.numRows()), diagonal(getData().base(), t.diagonal)
{
	Vector aux(t.getData());
	blas_int d_size = diagonal.getSize();
	double alpha = 1.0;
	double beta = 0.0;
	BLAS_dgemm(const_cast<char*>("N"), const_cast<char*>("N"), &d_size, &d_size, &d_size, &alpha, aux.base(),
		   &d_size, const_cast<double*>(t.getData().base()), &d_size, &beta, getData().base(), &d_size);
}

QuasiTriangular::QuasiTriangular(const SchurDecomp& decomp)
	: SqSylvMatrix(decomp.getT()),
	  diagonal(getData().base(), decomp.getDim())
{
}

/* this pads matrix with intial columns with zeros */
QuasiTriangular::QuasiTriangular(const SchurDecompZero& decomp)
	: SqSylvMatrix(decomp.getDim())
{
	// nullify first decomp.getZeroCols() columns
	int zeros = decomp.getZeroCols()*decomp.getDim();
	Vector zv(getData(), 0, zeros);
	zv.zeros();
	// fill right upper part with decomp.getRU()
	for (int i = 0; i < decomp.getRU().numRows(); i++) {
		for (int j = 0; j < decomp.getRU().numCols(); j++) {
			getData()[(j+decomp.getZeroCols())*decomp.getDim()+i] = decomp.getRU().get(i,j);
		}
	}
	// fill right lower part with decomp.getT()
	for (int i = 0; i < decomp.getT().numRows(); i++) {
		for (int j = 0; j < decomp.getT().numCols(); j++) {
			getData()[(j+decomp.getZeroCols())*decomp.getDim()+decomp.getZeroCols()+i] = 
				decomp.getT().get(i,j);
		}
	}
	// construct diagonal
	Diagonal* const d = new Diagonal(getData().base(), decomp.getDim());
	diagonal = *d;
	delete d;
}

void QuasiTriangular::setMatrix(double r, const QuasiTriangular& t)
{
	getData().zeros();
	getData().add(r, t.getData());
}

void QuasiTriangular::setMatrixViaIter(double r, const QuasiTriangular& t)
{
	register double rr = r;
	diag_iter dil = diag_begin();
	const_diag_iter dir = t.diag_begin();
	for ( ; dil != diag_end(); ++dil, ++dir) {
		(*dil).getAlpha() = rr*(*(*dir).getAlpha());
		if (! (*dil).isReal()) {
			(*dil).getBeta1() = rr*(*dir).getBeta1();
			(*dil).getBeta2() = rr*(*dir).getBeta2();
		}
		col_iter cil = col_begin(*dil);
		const_col_iter cir = t.col_begin(*dir);
		for ( ; cil != col_end(*dil); ++cil, ++cir) {
			if ((*dil).isReal()) {
				*cil = rr*(*cir);
			} else {
				cil.a() = rr*cir.a();
				cil.b() = rr*cir.b();
			}
		}
	}
}

void QuasiTriangular::addMatrix(double r, const QuasiTriangular& t)
{
	getData().add(r, t.getData());
}

void QuasiTriangular::addMatrixViaIter(double r, const QuasiTriangular& t)
{
	register double rr = r;
	diag_iter dil = diag_begin();
	const_diag_iter dir = t.diag_begin();
	for ( ; dil != diag_end(); ++dil, ++dir) {
		(*dil).getAlpha() = (*(*dil).getAlpha()) + rr*(*(*dir).getAlpha());
		if (! (*dil).isReal()) {
			(*dil).getBeta1() += rr*(*dir).getBeta1();
			(*dil).getBeta2() += rr*(*dir).getBeta2();
		}
		col_iter cil = col_begin(*dil);
		const_col_iter cir = t.col_begin(*dir);
		for ( ; cil != col_end(*dil); ++cil, ++cir) {
			if ((*dil).isReal()) {
				*cil += rr*(*cir);
			} else {
				cil.a() += rr*cir.a();
				cil.b() += rr*cir.b();
			}
		}
	}
}

void QuasiTriangular::addUnit()
{
	for (diag_iter di = diag_begin(); di != diag_end(); ++di) {
		(*di).getAlpha() = *((*di).getAlpha()) + 1.0;
	}
}

void QuasiTriangular::solve(Vector& x, const ConstVector& b, double& eig_min)
{
	x = b;
	solvePre(x, eig_min);
}

void QuasiTriangular::solveTrans(Vector& x, const ConstVector& b, double& eig_min)
{
	x = b;
	solvePreTrans(x, eig_min);
}

void QuasiTriangular::solvePre(Vector& x, double& eig_min)
{
	addUnit();
	for (diag_iter di = diag_begin(); di != diag_end(); ++di) {
		double eig_size;
		if (!(*di).isReal()) {
			eig_size = (*di).getDeterminant();
			eliminateLeft((*di).getIndex()+1, (*di).getIndex(), x);
		} else {
			eig_size = *(*di).getAlpha()*(*(*di).getAlpha());
		}
		if (eig_size < eig_min)
			eig_min = eig_size;
	}

	blas_int nn = diagonal.getSize();
	blas_int lda = diagonal.getSize();
	blas_int incx = x.skip();
	BLAS_dtrsv(const_cast<char*>("U"), const_cast<char*>("N"), const_cast<char*>("N"), &nn, getData().base(), &lda, x.base(), &incx);
}

void QuasiTriangular::solvePreTrans(Vector& x, double& eig_min)
{
	addUnit();
	for (diag_iter di = diag_begin(); di != diag_end(); ++di) {
		double eig_size;
		if (!(*di).isReal()) {
			eig_size = (*di).getDeterminant();
			eliminateRight((*di).getIndex()+1, (*di).getIndex(), x);
		} else {
			eig_size = *(*di).getAlpha()*(*(*di).getAlpha());
		}
		if (eig_size < eig_min)
			eig_min = eig_size;
	}
	
	blas_int nn = diagonal.getSize();
	blas_int lda = diagonal.getSize();
	blas_int incx = x.skip();
	BLAS_dtrsv(const_cast<char*>("U"), const_cast<char*>("T"), const_cast<char*>("N"), &nn, getData().base(), &lda, x.base(), &incx);
}


/* calculates x = Tb */
void QuasiTriangular::multVec(Vector& x, const ConstVector& b) const
{
	x = b;
	blas_int nn = diagonal.getSize();
	blas_int lda = diagonal.getSize();
	blas_int incx = x.skip();
	BLAS_dtrmv(const_cast<char*>("U"), const_cast<char*>("N"), const_cast<char*>("N"), &nn, 
		   const_cast<double*>(getData().base()), &lda, x.base(), &incx);
	for (const_diag_iter di = diag_begin(); di != diag_end(); ++di) {
		if (!(*di).isReal()) {
			int jbar = (*di).getIndex();
			x[jbar+1] += (*di).getBeta2()*(b[jbar]);
		}
	}
}


void QuasiTriangular::multVecTrans(Vector& x, const ConstVector& b) const
{
	x = b;
	blas_int nn = diagonal.getSize();
	blas_int lda = diagonal.getSize();
	blas_int incx = x.skip();
	BLAS_dtrmv(const_cast<char*>("U"), const_cast<char*>("T"), const_cast<char*>("N"), &nn, 
		   const_cast<double*>(getData().base()), &lda, x.base(), &incx);
	for (const_diag_iter di = diag_begin(); di != diag_end(); ++di) {
		if (!(*di).isReal()) {
			int jbar = (*di).getIndex();
			x[jbar] += (*di).getBeta2()*b[jbar+1];
		}
	}	
}

void QuasiTriangular::multaVec(Vector& x, const ConstVector& b) const
{
	Vector tmp((const Vector&) x); // new copy
	multVec(x, b);
	x.add(1.0, tmp);
}

void QuasiTriangular::multaVecTrans(Vector& x, const ConstVector& b) const
{
	Vector tmp((const Vector&) x); // new copy
	multVecTrans(x, b);
	x.add(1.0, tmp);
}

/* calculates x=x+(T\otimes I)b, where size of I is given by b (KronVector) */
void QuasiTriangular::multaKron(KronVector& x, const ConstKronVector& b) const
{
	int id = b.getN()*power(b.getM(), b.getDepth()-1);
	ConstGeneralMatrix b_resh(b.base(), id, b.getM());
	GeneralMatrix x_resh(x.base(), id, b.getM());
	x_resh.multAndAdd(b_resh, ConstGeneralMatrix(*this), "trans");
}


/* calculates x=x+(T'\otimes I)b, where size of I is given by b (KronVector) */
void
QuasiTriangular::multaKronTrans(KronVector& x, const ConstKronVector& b) const
{
	int id = b.getN()*power(b.getM(), b.getDepth()-1);
	ConstGeneralMatrix b_resh(b.base(), id, b.getM());
	GeneralMatrix x_resh(x.base(), id, b.getM());
	x_resh.multAndAdd(b_resh, ConstGeneralMatrix(*this));
}


void QuasiTriangular::multKron(KronVector& x) const
{
	KronVector b((const KronVector&)x); // make copy
	x.zeros();
	multaKron(x, b);
}

void
QuasiTriangular::multKronTrans(KronVector& x) const
{
	KronVector b((const KronVector&)x); // make copy
	x.zeros();
	multaKronTrans(x, b);
}

void QuasiTriangular::multLeftOther(GeneralMatrix& a) const
{
	a.multLeft(*this);
}

void QuasiTriangular::multLeftOtherTrans(GeneralMatrix& a) const
{
	a.multLeftTrans(*this);
}

void QuasiTriangular::swapDiagLogically(diag_iter it)
{
	diagonal.swapLogically(it);
}

void QuasiTriangular::checkDiagConsistency(diag_iter it)
{
	diagonal.checkConsistency(it);
}

double QuasiTriangular::getAverageDiagSize(diag_iter start, diag_iter end)
{
	return diagonal.getAverageSize(start, end);
}

QuasiTriangular::diag_iter
QuasiTriangular::findClosestDiagBlock(diag_iter start, diag_iter end, double a)
{
	return diagonal.findClosestBlock(start, end, a);
}

QuasiTriangular::diag_iter
QuasiTriangular::findNextLargerBlock(diag_iter start, diag_iter end, double a)
{
	return diagonal.findNextLargerBlock(start, end, a);
}

int QuasiTriangular::getNumOffdiagonal() const
{
	return 	diagonal.getSize()*(diagonal.getSize()-1)/2 - diagonal.getNumComplex();
}
