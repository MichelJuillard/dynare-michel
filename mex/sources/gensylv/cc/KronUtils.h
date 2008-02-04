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
