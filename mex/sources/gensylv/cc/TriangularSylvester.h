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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/TriangularSylvester.h,v 1.1.1.1 2004/06/04 13:01:03 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef TRIANGULAR_SYLVESTER_H
#define TRIANGULAR_SYLVESTER_H

#include "SylvesterSolver.h"
#include "KronVector.h"
#include "QuasiTriangular.h"
#include "QuasiTriangularZero.h"
#include "SimilarityDecomp.h"

class TriangularSylvester : public SylvesterSolver {
	const QuasiTriangular* const matrixKK;
	const QuasiTriangular* const matrixFF;
public:
	TriangularSylvester(const QuasiTriangular& k, const QuasiTriangular& f);
	TriangularSylvester(const SchurDecompZero& kdecomp, const SchurDecomp& fdecomp);
	TriangularSylvester(const SchurDecompZero& kdecomp, const SimilarityDecomp& fdecomp);
	virtual ~TriangularSylvester();
	void print() const;
	void solve(SylvParams& pars, KronVector& d) const;

	void solvi(double r, KronVector& d, double& eig_min) const;
	void solvii(double alpha, double beta1, double beta2,
				KronVector& d1, KronVector& d2,
				double& eig_min) const;
	void solviip(double alpha, double betas,
				 KronVector& d, double& eig_min) const;
	/* evaluates:
	   |x1|   |d1| |alpha -beta1|                       |d1|
	   |  | = |  |+|            |\otimes F'...\otimes K |  |
	   |x2|   |d2| |beta2  alpha|                       |d2|
	*/
	void linEval(double alpha, double beta1, double beta2,
				 KronVector& x1, KronVector& x2,
				 const ConstKronVector& d1, const ConstKronVector& d2) const;
	void linEval(double alpha, double beta1, double beta2,
				 KronVector& x1, KronVector& x2,
				 const KronVector& d1, const KronVector& d2) const
		{linEval(alpha, beta1, beta2, x1, x2,
				 ConstKronVector(d1), ConstKronVector(d2));}

	/* evaluates:
	   |x1|   |d1|          |gamma -delta1|                       |d1|
	   |  | = |  | + 2alpha*|             |\otimes F'...\otimes K |  | +
	   |x2|   |d2|          |delta2  gamma|                       |d2|

                                     |gamma -delta1|^2                       |d1|
                   + (alpha^2+betas)*|             |\otimes F'2...\otimes K2 |  |
                                     |delta2  gamma|                         |d2|
	*/
	void quaEval(double alpha, double betas,
				 double gamma, double delta1, double delta2,
				 KronVector& x1, KronVector& x2,
				 const ConstKronVector& d1, const ConstKronVector& d2) const;
	void quaEval(double alpha, double betas,
				 double gamma, double delta1, double delta2,
				 KronVector& x1, KronVector& x2,
				 const KronVector& d1, const KronVector& d2) const
		{quaEval(alpha, betas, gamma, delta1, delta2, x1, x2,
				 ConstKronVector(d1), ConstKronVector(d2));}
private:
	/* returns square of size of minimal eigenvalue of the system solved,
	   now obsolete */ 
	double getEigSep(int depth) const;
	/* recursivelly calculates kronecker product of complex vectors (used in getEigSep) */
	static void multEigVector(KronVector& eig, const Vector& feig, const Vector& keig);
	/* auxiliary typedefs */
	typedef QuasiTriangular::const_diag_iter const_diag_iter;
	typedef QuasiTriangular::const_row_iter const_row_iter;
	/* called from solvi */
	void solviRealAndEliminate(double r, const_diag_iter di,
							   KronVector& d, double& eig_min) const;
	void solviComplexAndEliminate(double r, const_diag_iter di,
								  KronVector& d, double& eig_min) const;
	/* called from solviip */
	void solviipRealAndEliminate(double alpha, double betas,
								 const_diag_iter di, const_diag_iter dsi,
								 KronVector& d, double& eig_min) const;
	void solviipComplexAndEliminate(double alpha, double betas,
									const_diag_iter di, const_diag_iter dsi,
									KronVector& d, double& eig_min) const;
	/* eliminations */
	void solviEliminateReal(const_diag_iter di, KronVector& d,
							const KronVector& y, double divisor) const;
	void solviEliminateComplex(const_diag_iter di, KronVector& d,
							   const KronVector& y1, const KronVector& y2,
							   double divisor) const;
	void solviipEliminateReal(const_diag_iter di, const_diag_iter dsi,
							  KronVector& d,
							  const KronVector& y1, const KronVector& y2,
							  double divisor, double divisor2) const;
	void solviipEliminateComplex(const_diag_iter di, const_diag_iter dsi,
								 KronVector& d,
								 const KronVector& y1, const KronVector& y11,
								 const KronVector& y2, const KronVector& y22,
								 double divisor) const;
	/* Lemma 2 */
	void solviipComplex(double alpha, double betas, double gamma,
						double delta1, double delta2,
						KronVector& d1, KronVector& d2,
						double& eig_min) const;
	/* norms for what we consider zero on diagonal of F */
	static double diag_zero;
	static double diag_zero_sq; // square of diag_zero
};

#endif /* TRIANGULAR_SYLVESTER_H */


// Local Variables:
// mode:C++
// End:
