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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SchurDecomp.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SCHUR_DECOMP_H
#define SCHUR_DECOMP_H

#include "SylvMatrix.h"
#include "QuasiTriangular.h"

class QuasiTriangular;
class SchurDecomp {
	bool q_destroy;
	SqSylvMatrix* q;
	bool t_destroy;
	QuasiTriangular* t;
public:
	SchurDecomp(const SqSylvMatrix& m);
	SchurDecomp(const QuasiTriangular& tr);
	SchurDecomp(QuasiTriangular& tr);
	const SqSylvMatrix& getQ() const {return *q;}
	const QuasiTriangular& getT() const {return *t;}
	SqSylvMatrix& getQ() {return *q;}
	QuasiTriangular& getT() {return *t;}
	virtual int getDim() const;
	virtual ~SchurDecomp();
};

class SchurDecompZero : public SchurDecomp {
	GeneralMatrix ru; /* right upper matrix */
public:
	SchurDecompZero(const GeneralMatrix& m);
	const GeneralMatrix& getRU() const {return ru;}
	int getDim() const;
	int getZeroCols() const {return ru.numRows();}
};

#endif /* SCHUR_DECOMP_H */


// Local Variables:
// mode:C++
// End:
