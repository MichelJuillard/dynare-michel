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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/IterativeSylvester.cpp,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#include "IterativeSylvester.h"
#include "KronUtils.h"

void IterativeSylvester::solve(SylvParams& pars, KronVector& x) const
{
	int max_steps = *(pars.max_num_iter);
	int steps = 1;
	double max_norm = *(pars.convergence_tol);
	double norm = performFirstStep(x);

	QuasiTriangular* kpow = matrixK->clone();
	QuasiTriangular* fpow = matrixF->clone();
	while (steps < max_steps && norm > max_norm) {
		kpow->multRight(SqSylvMatrix(*kpow)); // be careful to make copy
		fpow->multRight(SqSylvMatrix(*fpow)); // also here
		norm = performStep(*kpow, *fpow, x);
		steps++;
	}

	delete fpow;
	delete kpow;

	pars.converged = (norm <= max_norm);
	pars.iter_last_norm = norm;
	pars.num_iter = steps;
}

double IterativeSylvester::performFirstStep(KronVector& x) const
{
	KronVector xtmp((const KronVector&)x);
	KronUtils::multKron(*matrixF, *matrixK, xtmp);
	x.add(-1., xtmp);
	double norm = xtmp.getMax();
	return norm;
}

double IterativeSylvester::performStep(const QuasiTriangular& k, const QuasiTriangular& f,
									   KronVector& x)
{
	KronVector xtmp((const KronVector&)x);
	KronUtils::multKron(f, k, xtmp);
	x.add(1.0, xtmp);
	double norm = xtmp.getMax();
	return norm;
}

// Local Variables:
// mode:C++
// End:
