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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/BlockDiagonal.h,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef BLOCK_DIAGONAL_H
#define BLOCK_DIAGONAL_H

#include "QuasiTriangular.h"


class BlockDiagonal : public QuasiTriangular {
	int* const row_len;
	int* const col_len;
public:
	BlockDiagonal(const double* d, int d_size);
	BlockDiagonal(int p, const BlockDiagonal& b);
	BlockDiagonal(const BlockDiagonal& b);
	BlockDiagonal(const QuasiTriangular& t);
	const BlockDiagonal& operator=(const QuasiTriangular& t)
		{GeneralMatrix::operator=(t); return *this;}
	const BlockDiagonal& operator=(const BlockDiagonal& b);
	~BlockDiagonal() {delete [] row_len; delete [] col_len;}
	void setZeroBlockEdge(diag_iter edge);
	int getNumZeros() const;
	int getNumBlocks() const;
	int getLargestBlock() const;
	void printInfo() const;

	void multKron(KronVector& x) const;
	void multKronTrans(KronVector& x) const;

	const_col_iter col_begin(const DiagonalBlock& b) const;
	col_iter col_begin(const DiagonalBlock& b);
	const_row_iter row_end(const DiagonalBlock& b) const;
	row_iter row_end(const DiagonalBlock& b);
	QuasiTriangular* clone() const
		{return new BlockDiagonal(*this);}
private:
	void setZerosToRU(diag_iter edge);
	const_diag_iter findBlockStart(const_diag_iter from) const;
	static void savePartOfX(int si, int ei, const KronVector& x, Vector& work);
	void multKronBlock(const_diag_iter start, const_diag_iter end,
					   KronVector& x, Vector& work) const;
	void multKronBlockTrans(const_diag_iter start, const_diag_iter end,
							KronVector& x, Vector& work) const;
};

#endif /* BLOCK_DIAGONAL_H */


// Local Variables:
// mode:C++
// End:
