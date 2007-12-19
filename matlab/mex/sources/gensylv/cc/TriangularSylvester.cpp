/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/TriangularSylvester.cpp,v 1.1.1.1 2004/06/04 13:00:59 kamenik Exp $ */

/* Tag $Name:  $ */

#include "TriangularSylvester.h"
#include "QuasiTriangularZero.h"
#include "KronUtils.h"
#include "BlockDiagonal.h"

#include <stdio.h>
#include <cmath>

double TriangularSylvester::diag_zero = 1.e-15;
double TriangularSylvester::diag_zero_sq = 1.e-30;

TriangularSylvester::TriangularSylvester(const QuasiTriangular& k,
										 const QuasiTriangular& f)
	: SylvesterSolver(k, f),
	  matrixKK(matrixK->clone(2, *matrixK)),
	  matrixFF(new QuasiTriangular(2, *matrixF))
{
}

TriangularSylvester::TriangularSylvester(const SchurDecompZero& kdecomp,
										 const SchurDecomp& fdecomp)
	: SylvesterSolver(kdecomp, fdecomp),
	  matrixKK(matrixK->clone(2, *matrixK)),
	  matrixFF(new QuasiTriangular(2, *matrixF))
{
}

TriangularSylvester::TriangularSylvester(const SchurDecompZero& kdecomp,
										 const SimilarityDecomp& fdecomp)
	: SylvesterSolver(kdecomp, fdecomp),
	  matrixKK(matrixK->clone(2, *matrixK)),
	  matrixFF(new BlockDiagonal(2, *matrixF))
{
}

TriangularSylvester::~TriangularSylvester()
{
	delete matrixKK;
	delete matrixFF;
}

void TriangularSylvester::print() const
{
	printf("matrix K (%d):\n",matrixK->getDiagonal().getSize());
	matrixK->print();
	printf("matrix F (%d):\n",matrixF->getDiagonal().getSize());
	matrixF->print();
}

void TriangularSylvester::solve(SylvParams& pars, KronVector& d) const
{
	double eig_min = 1e30;
	solvi(1., d, eig_min);
	pars.eig_min = sqrt(eig_min);
}

void TriangularSylvester::solvi(double r, KronVector& d, double& eig_min) const
{
	if (d.getDepth() == 0) {
		QuasiTriangular* t = matrixK->clone(r);
		t->solvePre(d, eig_min);
		delete t;
	} else {
		for (const_diag_iter di = matrixF->diag_begin();
			 di != matrixF->diag_end();
			 ++di) {
			if ((*di).isReal()) {
				solviRealAndEliminate(r, di, d, eig_min);
			} else {
				solviComplexAndEliminate(r, di, d, eig_min);
			}
		}
	}
}


void TriangularSylvester::solvii(double alpha, double beta1, double beta2,
								 KronVector& d1, KronVector& d2,
								 double& eig_min) const
{
	KronVector d1tmp(d1);
	KronVector d2tmp(d2);
	linEval(alpha, beta1, beta2, d1, d2, d1tmp, d2tmp);
	solviip(alpha, beta1*beta2, d1, eig_min);
	solviip(alpha, beta1*beta2, d2, eig_min);
}


void TriangularSylvester::solviip(double alpha, double betas,
								  KronVector& d, double& eig_min) const
{
	// quick exit to solvi if betas is small
	if (betas < diag_zero_sq) {
		solvi(alpha, d, eig_min);
		solvi(alpha, d, eig_min);
		return;
	}

	if (d.getDepth() == 0) {
		double aspbs = alpha*alpha+betas;
		QuasiTriangular* t= matrixK->clone(2*alpha, aspbs, *matrixKK);
		t->solvePre(d, eig_min);
		delete t;
	} else {
		const_diag_iter di = matrixF->diag_begin();
		const_diag_iter dsi = matrixFF->diag_begin();
		for (; di != matrixF->diag_end(); ++di, ++dsi) {
			if ((*di).isReal()) {
				solviipRealAndEliminate(alpha, betas, di, dsi, d, eig_min);
			} else {
				solviipComplexAndEliminate(alpha, betas, di, dsi, d, eig_min);
			}
		}
	}
}


void TriangularSylvester::solviRealAndEliminate(double r, const_diag_iter di,
												KronVector& d, double& eig_min) const
{
	// di is real
	int jbar = (*di).getIndex();
	double f = *((*di).getAlpha());
	KronVector dj(d, jbar);
	// solve system
	if (abs(r*f) > diag_zero) { 
		solvi(r*f, dj, eig_min);
	}
	// calculate y
	KronVector y((const KronVector&)dj);
	KronUtils::multKron(*matrixF, *matrixK, y);
	y.mult(r);
	double divisor = 1.0;
	solviEliminateReal(di, d, y, divisor);
}

void TriangularSylvester::solviEliminateReal(const_diag_iter di, KronVector& d,
											 const KronVector& y, double divisor) const
{
	for (const_row_iter ri = matrixF->row_begin(*di);
		 ri != matrixF->row_end(*di);
		 ++ri) {
		KronVector dk(d, ri.getCol());
		dk.add(-(*ri)/divisor, y);
	}
}

void TriangularSylvester::solviComplexAndEliminate(double r, const_diag_iter di,
												   KronVector& d, double& eig_min) const
{
	// di is complex
	int jbar = (*di).getIndex();
	// pick data
	double alpha = *(*di).getAlpha();
	double beta1 = (*di).getBeta2();
	double beta2 = -(*di).getBeta1();
	double aspbs = (*di).getDeterminant();
	KronVector dj(d, jbar);
	KronVector djj(d, jbar+1);
	// solve
	if (r*r*aspbs > diag_zero_sq) { 
		solvii(r*alpha, r*beta1, r*beta2, dj, djj, eig_min);
	}
	KronVector y1(dj);
	KronVector y2(djj);
	KronUtils::multKron(*matrixF, *matrixK, y1);
	KronUtils::multKron(*matrixF, *matrixK, y2);
	y1.mult(r);
	y2.mult(r);
	double divisor = 1.0;
	solviEliminateComplex(di, d, y1, y2, divisor);
}

void TriangularSylvester::solviEliminateComplex(const_diag_iter di, KronVector& d,
												const KronVector& y1, const KronVector& y2,
												double divisor) const
{
	for (const_row_iter ri = matrixF->row_begin(*di);
		 ri != matrixF->row_end(*di);
		 ++ri) {
		KronVector dk(d, ri.getCol());
		dk.add(-ri.a()/divisor, y1);
		dk.add(-ri.b()/divisor, y2);
	}
}

void TriangularSylvester::solviipRealAndEliminate(double alpha, double betas,
												  const_diag_iter di, const_diag_iter dsi,
												  KronVector& d, double& eig_min) const
{
	// di, and dsi are real		
	int jbar = (*di).getIndex();
	double aspbs = alpha*alpha+betas;
	// pick data
	double f = *((*di).getAlpha());
	double fs = f*f;
	KronVector dj(d, jbar);
	// solve
	if (fs*aspbs > diag_zero_sq) {
		solviip(f*alpha, fs*betas, dj, eig_min);
	}
	KronVector y1((const KronVector&)dj);
	KronVector y2((const KronVector&)dj);
	KronUtils::multKron(*matrixF, *matrixK, y1);
	y1.mult(2*alpha);
	KronUtils::multKron(*matrixFF, *matrixKK, y2);
	y2.mult(aspbs);
	double divisor = 1.0;
	double divisor2 = 1.0;
	solviipEliminateReal(di, dsi, d, y1, y2, divisor, divisor2);
}

void TriangularSylvester::solviipEliminateReal(const_diag_iter di, const_diag_iter dsi,
											   KronVector& d,
											   const KronVector& y1, const KronVector& y2,
											   double divisor, double divisor2) const
{
	const_row_iter ri = matrixF->row_begin(*di);
	const_row_iter rsi = matrixFF->row_begin(*dsi);
	for (; ri != matrixF->row_end(*di); ++ri, ++rsi) {
		KronVector dk(d, ri.getCol());
		dk.add(-(*ri)/divisor, y1);
		dk.add(-(*rsi)/divisor2, y2);
	}
} 

void TriangularSylvester::solviipComplexAndEliminate(double alpha, double betas,
													 const_diag_iter di, const_diag_iter dsi,
													 KronVector& d, double& eig_min) const
{
	// di, and dsi are complex
	int jbar = (*di).getIndex();
	double aspbs = alpha*alpha+betas;
	// pick data
	double gamma = *((*di).getAlpha());
	double delta1 = (*di).getBeta2(); // swap because of transpose
	double delta2 = -(*di).getBeta1();
	double gspds = (*di).getDeterminant();
	KronVector dj(d, jbar);
	KronVector djj(d, jbar+1);
	if (gspds*aspbs > diag_zero_sq) {
		solviipComplex(alpha, betas, gamma, delta1, delta2, dj, djj, eig_min);
	}
	// here dj, djj is solution, set y1, y2, y11, y22
	// y1
	KronVector y1((const KronVector&) dj);
	KronUtils::multKron(*matrixF, *matrixK, y1);
	y1.mult(2*alpha);
	// y11
	KronVector y11((const KronVector&) djj);
	KronUtils::multKron(*matrixF, *matrixK, y11);
	y11.mult(2*alpha);
	// y2
	KronVector y2((const KronVector&) dj);
	KronUtils::multKron(*matrixFF, *matrixKK, y2);
	y2.mult(aspbs);
	// y22
	KronVector y22((const KronVector&) djj);
	KronUtils::multKron(*matrixFF, *matrixKK, y22);
	y22.mult(aspbs);

	double divisor = 1.0;
	solviipEliminateComplex(di, dsi, d, y1, y11, y2, y22, divisor);
}


void TriangularSylvester::solviipComplex(double alpha, double betas, double gamma,
										 double delta1, double delta2,
										 KronVector& d1, KronVector& d2,
										 double& eig_min) const
{
	KronVector d1tmp(d1);
	KronVector d2tmp(d2);
	quaEval(alpha, betas, gamma, delta1, delta2,
			d1, d2, d1tmp, d2tmp);
	double delta = sqrt(delta1*delta2);
	double beta = sqrt(betas);
	double a1 = alpha*gamma - beta*delta;
	double b1 = alpha*delta + gamma*beta;
	double a2 = alpha*gamma + beta*delta;
	double b2 = alpha*delta - gamma*beta;
	solviip(a2, b2*b2, d1, eig_min);
	solviip(a1, b1*b1, d1, eig_min);
	solviip(a2, b2*b2, d2, eig_min);
	solviip(a1, b1*b1, d2, eig_min);
}

void TriangularSylvester::solviipEliminateComplex(const_diag_iter di, const_diag_iter dsi,
												  KronVector& d,
												  const KronVector& y1, const KronVector& y11,
												  const KronVector& y2, const KronVector& y22,
												  double divisor) const
{
	const_row_iter ri = matrixF->row_begin(*di);
	const_row_iter rsi = matrixFF->row_begin(*dsi);
	for (; ri != matrixF->row_end(*di); ++ri, ++rsi) {
		KronVector dk(d, ri.getCol());
		dk.add(-ri.a()/divisor, y1);
		dk.add(-ri.b()/divisor, y11);
		dk.add(-rsi.a()/divisor, y2);
		dk.add(-rsi.b()/divisor, y22);
	}
}

void TriangularSylvester::linEval(double alpha, double beta1, double beta2,
								  KronVector& x1, KronVector& x2,
								  const ConstKronVector& d1, const ConstKronVector& d2) const
{
	KronVector d1tmp(d1); // make copy
	KronVector d2tmp(d2); // make copy
	KronUtils::multKron(*matrixF, *matrixK, d1tmp);
	KronUtils::multKron(*matrixF, *matrixK, d2tmp);
	x1 = d1;
	x2 = d2;
	Vector::mult2a(alpha, beta1, -beta2, x1, x2, d1tmp, d2tmp);
}

void TriangularSylvester::quaEval(double alpha, double betas,
								  double gamma, double delta1, double delta2,
								  KronVector& x1, KronVector& x2,
								  const ConstKronVector& d1, const ConstKronVector& d2) const
{
	KronVector d1tmp(d1); // make copy
	KronVector d2tmp(d2); // make copy
	KronUtils::multKron(*matrixF, *matrixK, d1tmp);
	KronUtils::multKron(*matrixF, *matrixK, d2tmp);
	x1 = d1;
	x2 = d2;
	Vector::mult2a(2*alpha*gamma, 2*alpha*delta1, -2*alpha*delta2,
				   x1, x2, d1tmp, d2tmp);
	d1tmp = d1; // restore to d1
	d2tmp = d2; // restore to d2
	KronUtils::multKron(*matrixFF, *matrixKK, d1tmp);
	KronUtils::multKron(*matrixFF, *matrixKK, d2tmp);
	double aspbs = alpha*alpha + betas;
	double gspds = gamma*gamma - delta1*delta2;
	Vector::mult2a(aspbs*gspds, 2*aspbs*gamma*delta1, -2*aspbs*gamma*delta2,
				   x1, x2, d1tmp, d2tmp);
}


double TriangularSylvester::getEigSep(int depth) const
{
	int f_size = matrixF->getDiagonal().getSize();
	Vector feig(2*f_size);
	matrixF->getDiagonal().getEigenValues(feig);
	int k_size = matrixK->getDiagonal().getSize();
	Vector keig(2*k_size);
	matrixK->getDiagonal().getEigenValues(keig);
	
	KronVector eig(f_size, 2*k_size, depth);
	multEigVector(eig, feig, keig);

	double min = 1.0e20;
	for (int i = 0; i < eig.length()/2; i++) {
		double alpha = eig[2*i];
		double beta = eig[2*i+1];
		double ss = (alpha+1)*(alpha+1)+beta*beta;
		if (min > ss)
			min = ss;
	}

	return min;
}

void TriangularSylvester::multEigVector(KronVector& eig, const Vector& feig,
										const Vector& keig)
{
	int depth = eig.getDepth();
	int m = eig.getM();
	int n = eig.getN();

	if (depth == 0) {
		eig = keig;
	} else {
		KronVector aux(m, n, depth-1);
		multEigVector(aux, feig, keig);
		for (int i = 0; i < m; i++) {
			KronVector eigi(eig, i);
			eigi.zeros();
			eigi.add(&feig[2*i], aux);
		}
	}
}

// Local Variables:
// mode:C++
// End:
