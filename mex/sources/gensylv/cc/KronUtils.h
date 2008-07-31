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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/KronUtils.h,v 1.1.1.1 2004/06/04 13:00:31 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef KRON_UTILS_H
#define KRON_UTILS_H

#include "KronVector.h"
#include "QuasiTriangular.h"

class KronUtils {
public:
	/* multiplies I_m\otimes..\I_m\otimes T\otimes I_m...I_m\otimes I_n
	   with given b and returns x. T must be (m,m), number of
	   \otimes is b.getDepth(), level is a number of I_m's between T
	   and I_n plus 1. If level=0, then we multiply
       \I_m\otimes ..\otimes I_m\otimes T, T is (n,n) */
	static void multAtLevel(int level, const QuasiTriangular& t,
							KronVector& x);
	static void multAtLevelTrans(int level, const QuasiTriangular& t,
								 KronVector& x);

	/* multiplies x=(F'\otimes F'\otimes..\otimes K)x */
	static void multKron(const QuasiTriangular& f, const QuasiTriangular& k,
						 KronVector& x);
};

#endif /* KRON_UTILS_H */

// Local Variables:
// mode:C++
// End:
