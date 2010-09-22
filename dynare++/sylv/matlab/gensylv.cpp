/* $Header: /var/lib/cvs/dynare_cpp/sylv/matlab/gensylv.cpp,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#include "dynmex.h"
#include "mex.h"

#include "GeneralSylvester.h"
#include "SylvException.h"


void gen_sylv_solve(int order, int n, int m, int zero_cols,
					const double* A, const double* B,
					const double* C, double* X)
{
		GeneralSylvester s(order, n, m, zero_cols, A, B, C, X, false);
		s.solve();
}

void gen_sylv_solve_and_check(int order, int n, int m, int zero_cols,
							  const double* A, const double* B,
							  const double* C, const double* D,
							  double* X, mxArray*& p)
{
		GeneralSylvester s(order, n, m, zero_cols, A, B, C, X, true);
		s.solve();
		s.check(D);
		p = s.getParams().createStructArray();
}

extern "C" {
	void mexFunction(int nlhs, mxArray* plhs[],
					 int nhrs, const mxArray* prhs[])
	{
		if (nhrs != 5 || nlhs > 3 || nlhs < 2 )
                  DYN_MEX_FUNC_ERR_MSG_TXT("Gensylv: Must have exactly 5 input args and either 2 or 3 output args.");

		int order = (int)mxGetScalar(prhs[0]);
		const mxArray* const A = prhs[1];
		const mxArray* const B = prhs[2];
		const mxArray* const C = prhs[3];
		const mxArray* const D = prhs[4];
		const mwSize* const Adims = mxGetDimensions(A);
		const mwSize* const Bdims = mxGetDimensions(B);
		const mwSize* const Cdims = mxGetDimensions(C);
		const mwSize* const Ddims = mxGetDimensions(D);

                if (Adims[0] != Adims[1])
                  DYN_MEX_FUNC_ERR_MSG_TXT("Matrix A must be a square matrix.");
                if (Adims[0] != Bdims[0])
                  DYN_MEX_FUNC_ERR_MSG_TXT("Matrix A and matrix B must have the same number of rows.");
                if (Adims[0] != Ddims[0])
                  DYN_MEX_FUNC_ERR_MSG_TXT("Matrix A and matrix B must have the same number of rows.");
                if (Cdims[0] != Cdims[1])
                  DYN_MEX_FUNC_ERR_MSG_TXT("Matrix C must be square.");
                if (Bdims[0] < Bdims[1])
                  DYN_MEX_FUNC_ERR_MSG_TXT("Matrix B must not have more columns than rows.");
                if (Ddims[1] != power(Cdims[0], order))
                  DYN_MEX_FUNC_ERR_MSG_TXT("Matrix D has wrong number of columns.");

		int n = Adims[0];
		int m = Cdims[0];
		int zero_cols = Bdims[0] - Bdims[1];
		mxArray* X = mxCreateDoubleMatrix(Ddims[0], Ddims[1], mxREAL);
		// copy D to X
		Vector Xvec((double*)mxGetPr(X), power(m, order)*n);
		ConstVector Dvec((double*)mxGetPr(D), power(m, order)*n);
		Xvec = Dvec;
		// solve (or solve and check)
                try
                  {
                    if (nlhs == 2) {
                      gen_sylv_solve(order, n, m, zero_cols,
						   mxGetPr(A), mxGetPr(B), mxGetPr(C),
						   mxGetPr(X));
                    } else if (nlhs == 3) {
                      gen_sylv_solve_and_check(order, n, m, zero_cols,
									 mxGetPr(A), mxGetPr(B), mxGetPr(C),
									 mxGetPr(D), mxGetPr(X), plhs[2]);
                    }
                  }
                catch (const SylvException& e)
                  {
                    char mes[1000];
                    e.printMessage(mes, 999);
                    DYN_MEX_FUNC_ERR_MSG_TXT(mes);
                  }
                plhs[1] = X;
                plhs[0] = mxCreateDoubleScalar(0);
	}
};
