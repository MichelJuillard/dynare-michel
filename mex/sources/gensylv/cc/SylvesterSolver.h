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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvesterSolver.h,v 1.1.1.1 2004/06/04 13:00:54 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLVESTER_SOLVER_H
#define SYLVESTER_SOLVER_H

#include "KronVector.h"
#include "QuasiTriangular.h"
#include "QuasiTriangularZero.h"
#include "SimilarityDecomp.h"
#include "SylvParams.h"
#include "SchurDecomp.h"

class SylvesterSolver {
protected:
	const QuasiTriangular* const matrixK;
	const QuasiTriangular* const matrixF;
private:
	/* return true when it is more efficient to use QuasiTriangular
	 * than QuasiTriangularZero */
	static bool zeroPad(const SchurDecompZero& kdecomp) {
		return ((kdecomp.getZeroCols()*3 < kdecomp.getDim()*2) ||
				(kdecomp.getZeroCols() < 10));
	}
public:
	SylvesterSolver(const QuasiTriangular& k, const QuasiTriangular& f)
		: matrixK(new QuasiTriangular(k)),
		  matrixF(new QuasiTriangular(f))
		{}
	SylvesterSolver(const SchurDecompZero& kdecomp, const SchurDecomp& fdecomp)
		: matrixK((zeroPad(kdecomp)) ?
				  new QuasiTriangular(kdecomp) : new QuasiTriangularZero(kdecomp)),
		  matrixF(new QuasiTriangular(fdecomp))
		{}
	SylvesterSolver(const SchurDecompZero& kdecomp, const SimilarityDecomp& fdecomp)
		: matrixK((zeroPad(kdecomp)) ?
				  new QuasiTriangular(kdecomp) : new QuasiTriangularZero(kdecomp)),
		  matrixF(new BlockDiagonal(fdecomp.getB()))
		{}
	virtual ~SylvesterSolver()
		{delete matrixK; delete matrixF;}
	virtual void solve(SylvParams& pars, KronVector& x) const = 0;
};

#endif /* SYLVESTER_SOLVER_H */


// Local Variables:
// mode:C++
// End:
