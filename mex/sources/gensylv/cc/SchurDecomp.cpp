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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SchurDecomp.cpp,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SchurDecomp.h"

#include "cpplapack.h"

SchurDecomp::SchurDecomp(const SqSylvMatrix& m)
	: q_destroy(true), t_destroy(true)
{
	int rows = m.numRows();
	q = new SqSylvMatrix(rows);
	SqSylvMatrix auxt(m);
	int sdim;
	double* const wr = new double[rows];
	double* const wi = new double[rows];
	int lwork = 6*rows;
	double* const work = new double[lwork];
	int info;
	LAPACK_dgees("V", "N", 0, &rows, auxt.base(), &rows, &sdim,
				 wr, wi, q->base(), &rows,
				 work, &lwork, 0, &info);
	delete [] work;
	delete [] wi;
	delete [] wr;
	t = new QuasiTriangular(auxt.base(), rows);
}

SchurDecomp::SchurDecomp(const QuasiTriangular& tr)
	: q_destroy(true), t_destroy(true)
{
	q = new SqSylvMatrix(tr.numRows());
	q->setUnit();
	t = new QuasiTriangular(tr);
}

SchurDecomp::SchurDecomp(QuasiTriangular& tr)
	: q_destroy(true), t_destroy(false)
{
	q = new SqSylvMatrix(tr.numRows());
	q->setUnit();
	t = &tr;
}

SchurDecomp::~SchurDecomp()
{
	if (t_destroy)
		delete t;
	if (q_destroy)
		delete q;
}

int SchurDecomp::getDim() const
{
	return t->numRows();
}

SchurDecompZero::SchurDecompZero(const GeneralMatrix& m)
	: SchurDecomp(SqSylvMatrix(m, m.numRows()-m.numCols(), 0, m.numCols())),
	  ru(m, 0, 0, m.numRows()-m.numCols(), m.numCols())
{
	ru.multRight(getQ());
}

int SchurDecompZero::getDim() const
{
	return getT().numRows()+ru.numRows();
}


