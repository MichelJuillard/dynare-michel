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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/KronUtils.cpp,v 1.1.1.1 2004/06/04 13:00:31 kamenik Exp $ */

/* Tag $Name:  $ */

#include "KronUtils.h"

void KronUtils::multAtLevel(int level, const QuasiTriangular& t,
							KronVector& x)
{
	if (0 < level && level < x.getDepth()) {
		for (int i = 0; i < x.getM(); i++) {
			KronVector xi(x, i);
			multAtLevel(level, t, xi);
		}
	} else if (0 == level && 0 < x.getDepth()) {
		GeneralMatrix tmp(x.base(), x.getN(), power(x.getM(),x.getDepth()));
		t.multLeftOther(tmp);
	} else if (0 == level && 0 == x.getDepth()) {
		Vector b((const Vector&)x);
		t.multVec(x,b);
	} else { // 0 < level == depth
		t.multKron(x);
	}
}

void KronUtils::multAtLevelTrans(int level, const QuasiTriangular& t,
								 KronVector& x)
{
	if (0 < level && level < x.getDepth()) {
		for (int i = 0; i < x.getM(); i++) {
			KronVector xi(x, i);
			multAtLevelTrans(level, t, xi);
		}
	} else if (0 == level && 0 < x.getDepth()) {
		GeneralMatrix tmp(x.base(), x.getN(), power(x.getM(),x.getDepth()));
		t.multLeftOtherTrans(tmp);
	} else if (level == 0 && 0 == x.getDepth()) {
		Vector b((const Vector&)x);
		t.multVecTrans(x,b);
	} else { // 0 < level == depth
		t.multKronTrans(x);
	}
}

void KronUtils::multKron(const QuasiTriangular& f, const QuasiTriangular& k,
						 KronVector& x)
{
	multAtLevel(0, k, x);
	if (x.getDepth() > 0) {
		for (int level = 1; level <= x.getDepth(); level++)
			multAtLevelTrans(level, f, x);
	}
}
