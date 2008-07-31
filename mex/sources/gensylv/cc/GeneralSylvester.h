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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/GeneralSylvester.h,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef GENERAL_SYLVESTER_H
#define GENERAL_SYLVESTER_H

#include "SylvMatrix.h"
#include "SylvMemory.h"
#include "SimilarityDecomp.h"
#include "SylvesterSolver.h"

class GeneralSylvester {
	SylvParams pars;
	SylvMemoryDriver mem_driver;
	int order;
	const SqSylvMatrix a;
	const SylvMatrix b;
	const SqSylvMatrix c;
	SylvMatrix d;
	bool solved;
	SchurDecompZero* bdecomp;
	SimilarityDecomp* cdecomp;
	SylvesterSolver* sylv;
public:
	/* construct with my copy of d*/
	GeneralSylvester(int ord, int n, int m, int zero_cols,
					 const double* da, const double* db,
					 const double* dc, const double* dd,
					 const SylvParams& ps);
	GeneralSylvester(int ord, int n, int m, int zero_cols,
					 const double* da, const double* db,
					 const double* dc, const double* dd,
					 bool alloc_for_check = false);
	/* construct with provided storage for d */
	GeneralSylvester(int ord, int n, int m, int zero_cols,
					 const double* da, const double* db,
					 const double* dc, double* dd,
					 bool alloc_for_check = false);
	GeneralSylvester(int ord, int n, int m, int zero_cols,
					 const double* da, const double* db,
					 const double* dc, double* dd,
					 const SylvParams& ps);
	virtual ~GeneralSylvester();
	int getM() const {return c.numRows();}
	int getN() const {return a.numRows();}
	const double* getResult() const {return d.base();}
	const SylvParams& getParams() const {return pars;}
	SylvParams& getParams() {return pars;}
	void solve();
	void check(const double* ds);
private:
	void init();
};

#endif /* GENERAL_SYLVESTER_H */


// Local Variables:
// mode:C++
// End:
