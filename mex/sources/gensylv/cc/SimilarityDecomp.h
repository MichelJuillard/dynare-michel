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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SimilarityDecomp.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SIMILARITY_DECOMP_H
#define SIMILARITY_DECOMP_H

#include "SylvMatrix.h"
#include "BlockDiagonal.h"
#include "SylvParams.h"

class SimilarityDecomp {
	SqSylvMatrix* q;
	BlockDiagonal* b;
	SqSylvMatrix* invq;
	typedef BlockDiagonal::diag_iter diag_iter;
public:
	SimilarityDecomp(const double* d, int d_size, double log10norm = 3.0);
	virtual ~SimilarityDecomp();
	const SqSylvMatrix& getQ() const
		{return *q;}
	const SqSylvMatrix& getInvQ() const
		{return *invq;}
	const BlockDiagonal& getB() const
		{return *b;}
	void check(SylvParams& pars, const GeneralMatrix& m) const;
	void infoToPars(SylvParams& pars) const;
protected:
	void getXDim(diag_iter start, diag_iter end, int& rows, int& cols) const;
	bool solveX(diag_iter start, diag_iter end, GeneralMatrix& X, double norm) const;
	void updateTransform(diag_iter start, diag_iter end, GeneralMatrix& X);
	void bringGuiltyBlock(diag_iter start, diag_iter& end);
	void diagonalize(double norm);
};

#endif /* SIMILARITY_DECOMP_H */


// Local Variables:
// mode:C++
// End:
