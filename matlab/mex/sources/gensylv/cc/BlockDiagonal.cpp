/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/BlockDiagonal.cpp,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#include "BlockDiagonal.h"

#include <stdio.h>
#include <string.h>

BlockDiagonal::BlockDiagonal(const double* d, int d_size)
	: QuasiTriangular(d, d_size),
	  row_len(new int[d_size]), col_len(new int[d_size])
{
	for (int i = 0; i < d_size; i++) {
		row_len[i] = d_size;
		col_len[i] = 0;
	}
}

BlockDiagonal::BlockDiagonal(const QuasiTriangular& t)
	: QuasiTriangular(t),
	  row_len(new int[t.numRows()]), col_len(new int[t.numRows()])
{
	for (int i = 0; i < t.numRows(); i++) {
		row_len[i] = t.numRows();
		col_len[i] = 0;
	}
}

BlockDiagonal::BlockDiagonal(int p, const BlockDiagonal& b)
	: QuasiTriangular(p, b),
	  row_len(new int[b.numRows()]), col_len(new int[b.numRows()])
{
	memcpy(row_len, b.row_len, b.numRows()*sizeof(int));
	memcpy(col_len, b.col_len, b.numRows()*sizeof(int));
}

BlockDiagonal::BlockDiagonal(const BlockDiagonal& b)
	: QuasiTriangular(b),
	  row_len(new int[b.numRows()]), col_len(new int[b.numRows()])
{
	memcpy(row_len, b.row_len, b.numRows()*sizeof(int));
	memcpy(col_len, b.col_len, b.numRows()*sizeof(int));
}

/* put zeroes to right upper submatrix whose first column is defined
 * by 'edge' */
void BlockDiagonal::setZerosToRU(diag_iter edge)
{
	int iedge = (*edge).getIndex();
	for (int i = 0; i < iedge; i++)
		for (int j = iedge; j < numCols(); j++)
			get(i,j) = 0.0;
}

/* Updates row_len and col_len so that there are zeroes in upper right part, this
 * |T1 0 |
 * |0  T2|. The first column of T2 is given by diagonal iterator 'edge'.

 * Note the semantics of row_len and col_len. row_len[i] is distance
 * of the right-most non-zero element of i-th row from the left, and
 * col_len[j] is distance of top-most non-zero element of j-th column
 * to the top. (First element has distance 1).
 */
void BlockDiagonal::setZeroBlockEdge(diag_iter edge)
{
	setZerosToRU(edge);

	int iedge = (*edge).getIndex();
	for (diag_iter run = diag_begin(); run != edge; ++run) {
		int ind = (*run).getIndex();
		if (row_len[ind] > iedge) {
			row_len[ind] = iedge;
			if (!(*run).isReal())
				row_len[ind+1] = iedge;
		} 
	}
	for (diag_iter run = edge; run != diag_end(); ++run) {
		int ind = (*run).getIndex();
		if (col_len[ind] < iedge) {
			col_len[ind] = iedge;
			if (!(*run).isReal())
				col_len[ind+1] = iedge;
		}
	}
}

BlockDiagonal::const_col_iter
BlockDiagonal::col_begin(const DiagonalBlock& b) const
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return const_col_iter(&getData()[jbar*d_size + col_len[jbar]], d_size,
						  b.isReal(), col_len[jbar]);
}

BlockDiagonal::col_iter
BlockDiagonal::col_begin(const DiagonalBlock& b)
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return col_iter(&getData()[jbar*d_size + col_len[jbar]], d_size,
					b.isReal(), col_len[jbar]);
}


BlockDiagonal::const_row_iter
BlockDiagonal::row_end(const DiagonalBlock& b) const
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return const_row_iter(&getData()[d_size*row_len[jbar]+jbar], d_size,
						  b.isReal(), row_len[jbar]);
}

BlockDiagonal::row_iter
BlockDiagonal::row_end(const DiagonalBlock& b)
{
	int jbar = b.getIndex();
	int d_size = diagonal.getSize();
	return row_iter(&getData()[d_size*row_len[jbar]+jbar], d_size,
					b.isReal(), row_len[jbar]);
}

int BlockDiagonal::getNumZeros() const
{
	int sum = 0;
	for (int i = 0; i < diagonal.getSize(); i++) {
		sum += diagonal.getSize() - row_len[i];
	}
	return sum;
}

QuasiTriangular::const_diag_iter
BlockDiagonal::findBlockStart(const_diag_iter from) const
{
	if (from != diag_end()) {
		++from;
		while (from != diag_end() &&
			   col_len[(*from).getIndex()] != (*from).getIndex())
			++from;
	}
	return from;
}

int BlockDiagonal::getLargestBlock() const
{
	int largest = 0;
	const_diag_iter start = diag_begin();
	const_diag_iter end = findBlockStart(start);
	while (start != diag_end()) {
		int si = (*start).getIndex();
		int ei = diagonal.getSize();
		if (end != diag_end())
			ei = (*end).getIndex();
		if (largest < ei-si)
			largest = ei-si;
		start = end;
		end = findBlockStart(start);
	}
	return largest;
}


void BlockDiagonal::savePartOfX(int si, int ei, const KronVector& x, Vector& work)
{
	for (int i = si; i < ei; i++) {
		ConstKronVector xi(x, i);
		Vector target(work, (i-si)*xi.length(), xi.length());
		target = xi;
	}
}

void BlockDiagonal::multKronBlock(const_diag_iter start, const_diag_iter end,
								  KronVector& x, Vector& work) const
{
	int si = (*start).getIndex();
	int ei = diagonal.getSize();
	if (end != diag_end())
		ei = (*end).getIndex();
	savePartOfX(si, ei, x, work);

	for (const_diag_iter di = start; di != end; ++di) {
		int jbar = (*di).getIndex();
		if ((*di).isReal()) {
			KronVector xi(x, jbar);
			xi.zeros();
			Vector wi(work, (jbar-si)*xi.length(), xi.length());
			xi.add(*((*di).getAlpha()), wi);
			for (const_row_iter ri = row_begin(*di); ri != row_end(*di); ++ri) {
				int col = ri.getCol();
				Vector wj(work, (col-si)*xi.length(), xi.length());
				xi.add(*ri, wj);
			}
		} else {
			KronVector xi(x, jbar);
			KronVector xii(x, jbar+1);
			xi.zeros();
			xii.zeros();
			Vector wi(work, (jbar-si)*xi.length(), xi.length());
			Vector wii(work, (jbar+1-si)*xi.length(), xi.length());
			xi.add(*((*di).getAlpha()), wi);
			xi.add((*di).getBeta1(), wii);
			xii.add((*di).getBeta2(), wi);
			xii.add(*((*di).getAlpha()), wii);
			for (const_row_iter ri = row_begin(*di); ri != row_end(*di); ++ri) {
				int col = ri.getCol();
				Vector wj(work, (col-si)*xi.length(), xi.length());
				xi.add(ri.a(), wj);
				xii.add(ri.b(), wj);
			}
		}
	}
}

void BlockDiagonal::multKronBlockTrans(const_diag_iter start, const_diag_iter end,
									   KronVector& x, Vector& work) const
{
	int si = (*start).getIndex();
	int ei = diagonal.getSize();
	if (end != diag_end())
		ei = (*end).getIndex();
	savePartOfX(si, ei, x, work);

	for (const_diag_iter di = start; di != end; ++di) {
		int jbar = (*di).getIndex();
		if ((*di).isReal()) {
			KronVector xi(x, jbar);
			xi.zeros();
			Vector wi(work, (jbar-si)*xi.length(), xi.length());
			xi.add(*((*di).getAlpha()), wi);
			for (const_col_iter ci = col_begin(*di); ci != col_end(*di); ++ci) {
				int row = ci.getRow();
				Vector wj(work, (row-si)*xi.length(), xi.length());
				xi.add(*ci, wj);
			}
		} else {
			KronVector xi(x, jbar);
			KronVector xii(x, jbar+1);
			xi.zeros();
			xii.zeros();
			Vector wi(work, (jbar-si)*xi.length(), xi.length());
			Vector wii(work, (jbar+1-si)*xi.length(), xi.length());
			xi.add(*((*di).getAlpha()), wi);
			xi.add((*di).getBeta2(), wii);
			xii.add((*di).getBeta1(), wi);
			xii.add(*((*di).getAlpha()), wii);
			for (const_col_iter ci = col_begin(*di); ci != col_end(*di); ++ci) {
				int row = ci.getRow();
				Vector wj(work, (row-si)*xi.length(), xi.length());
				xi.add(ci.a(), wj);
				xii.add(ci.b(), wj);
			}
		}
	}
}

void BlockDiagonal::multKron(KronVector& x) const
{
	int largest = getLargestBlock();
	Vector work(largest*x.getN()*power(x.getM(),x.getDepth()-1));
	const_diag_iter start = diag_begin();
	const_diag_iter end = findBlockStart(start);
	while (start != diag_end()) {
		multKronBlock(start, end, x, work);
		start = end;
		end = findBlockStart(start);
	}
}


void BlockDiagonal::multKronTrans(KronVector& x) const
{
	int largest = getLargestBlock();
	Vector work(largest*x.getN()*power(x.getM(),x.getDepth()-1));
	const_diag_iter start = diag_begin();
	const_diag_iter end = findBlockStart(start);
	while (start != diag_end()) {
		multKronBlockTrans(start, end, x, work);
		start = end;
		end = findBlockStart(start);
	}
}

void BlockDiagonal::printInfo() const
{
	printf("Block sizes:");
	int num_blocks = 0;
	const_diag_iter start = diag_begin();
	const_diag_iter end = findBlockStart(start);
	while (start != diag_end()) {
		int si = (*start).getIndex();
		int ei = diagonal.getSize();
		if (end != diag_end())
			ei = (*end).getIndex();
		printf(" %d", ei-si);
		num_blocks++;
		start = end;
		end = findBlockStart(start);
	}
	printf("\nNum blocks: %d\n", num_blocks);
	printf("There are %d zeros out of %d\n",
		   getNumZeros(), getNumOffdiagonal());
}

int BlockDiagonal::getNumBlocks() const
{
	int num_blocks = 0;
	const_diag_iter start = diag_begin();
	const_diag_iter end = findBlockStart(start);
	while (start != diag_end()) {
		num_blocks++;
		start = end;
		end = findBlockStart(start);
	}
	return num_blocks;
}

// Local Variables:
// mode:C++
// End:
