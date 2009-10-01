/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SchurDecomp.cpp,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SchurDecomp.h"

#include <dynlapack.h>

SchurDecomp::SchurDecomp(const SqSylvMatrix& m)
	: q_destroy(true), t_destroy(true)
{
	lapack_int rows = m.numRows();
	q = new SqSylvMatrix(rows);
	SqSylvMatrix auxt(m);
	lapack_int sdim;
	double* const wr = new double[rows];
	double* const wi = new double[rows];
	lapack_int lwork = 6*rows;
	double* const work = new double[lwork];
	lapack_int info;
	dgees("V", "N", 0, &rows, auxt.base(), &rows, &sdim,
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


