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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvMatrix.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLV_MATRIX_H
#define SYLV_MATRIX_H

#include "GeneralMatrix.h"
#include "KronVector.h"

class SqSylvMatrix;

class SylvMatrix : public GeneralMatrix {
public:
	SylvMatrix(int m, int n)
		: GeneralMatrix(m,n) {}
	SylvMatrix(const double* d, int m, int n)
		: GeneralMatrix(d, m, n) {}
	SylvMatrix(double* d, int m, int n)
		: GeneralMatrix(d, m, n) {}
	SylvMatrix(const GeneralMatrix& m)
		: GeneralMatrix(m) {}
	SylvMatrix(const GeneralMatrix& m, int i, int j, int nrows, int ncols)
		: GeneralMatrix(m, i, j, nrows, ncols) {}
	SylvMatrix(GeneralMatrix& m, int i, int j, int nrows, int ncols)
		: GeneralMatrix(m, i, j, nrows, ncols) {}
	SylvMatrix(const GeneralMatrix& a, const GeneralMatrix& b)
		: GeneralMatrix(a, b) {}

	/* this = |I 0|* this
	          |0 m|         */
	void multLeftI(const SqSylvMatrix& m);
	/* this = |I  0|* this
	          |0 m'|         */
	void multLeftITrans(const SqSylvMatrix& m);
	/* this = |0 a|*b, so that |0 a| is square */
	void multLeft(int zero_cols, const GeneralMatrix& a, const GeneralMatrix& b);
	/* this = this * (m\otimes m..\otimes m) */
	void multRightKron(const SqSylvMatrix& m, int order);
	/* this = this * (m'\otimes m'..\otimes m') */
	void multRightKronTrans(const SqSylvMatrix& m, int order);
	/* this = P*this, x = P*x, where P is gauss transformation setting
	 * a given element to zero */
	void eliminateLeft(int row, int col, Vector& x);
	/* this = this*P, x = P'*x, where P is gauss transformation setting
	 * a given element to zero */
	void eliminateRight(int row, int col, Vector& x);
};


class SqSylvMatrix : public SylvMatrix {
public:
	SqSylvMatrix(int m) : SylvMatrix(m, m) {}
	SqSylvMatrix(const double* d, int m) : SylvMatrix(d, m, m) {}
	SqSylvMatrix(double* d, int m) : SylvMatrix(d, m, m) {}
	SqSylvMatrix(const SqSylvMatrix& m) : SylvMatrix(m) {}
	SqSylvMatrix(const GeneralMatrix& m, int i, int j, int nrows)
		: SylvMatrix(m, i, j, nrows, nrows) {} 
	SqSylvMatrix(GeneralMatrix& m, int i, int j, int nrows)
		: SylvMatrix(m, i, j, nrows, nrows) {} 
	SqSylvMatrix(const GeneralMatrix& a, const GeneralMatrix& b);
	const SqSylvMatrix& operator=(const SqSylvMatrix& m)
		{GeneralMatrix::operator=(m); return *this;}
	/* x = (this \otimes this..\otimes this)*d */
	void multVecKron(KronVector& x, const KronVector& d) const;
	/* x = (this' \otimes this'..\otimes this')*d */
	void multVecKronTrans(KronVector& x, const KronVector& d) const;
	/* a = inv(this)*a, b=inv(this)*b */
	void multInvLeft2(GeneralMatrix& a, GeneralMatrix& b,
					  double& rcond1, double& rcondinf) const;
	/* this = I */
	void setUnit();
};


#endif /* SYLV_MATRIX_H */


// Local Variables:
// mode:C++
// End:
