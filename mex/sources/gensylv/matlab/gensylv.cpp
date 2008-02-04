/* $Header: /var/lib/cvs/dynare_cpp/sylv/matlab/gensylv.cpp,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */


#include "mex.h"

#include "GeneralSylvester.h"
#include "SylvException.h"


void gen_sylv_solve(int order, int n, int m, int zero_cols,
					const double* A, const double* B,
					const double* C, double* X)
{
	try {
		GeneralSylvester s(order, n, m, zero_cols, A, B, C, X, false);
		s.solve();
	} catch (const SylvException& e) {
		char mes[1000];
		e.printMessage(mes, 999);
		mexErrMsgTxt(mes);
	}
}

void gen_sylv_solve_and_check(int order, int n, int m, int zero_cols,
							  const double* A, const double* B,
							  const double* C, const double* D,
							  double* X, mxArray*& p)
{
	try {
		GeneralSylvester s(order, n, m, zero_cols, A, B, C, X, true);
		s.solve();
		s.check(D);
		p = s.getParams().createStructArray();
	} catch (const SylvException& e) {
		char mes[1000];
		e.printMessage(mes, 999);
		mexErrMsgTxt(mes);
	}
}

void checkDimensions(const mwSize* const Adims, const mwSize* const Bdims,
					 const mwSize* const Cdims, const mwSize* const Ddims,
					 int order)
{
	if (Adims[0] != Adims[1])
		mexErrMsgTxt("Matrix A must be a square matrix.");
	if (Adims[0] != Bdims[0])
		mexErrMsgTxt("Matrix A and matrix B must have the same number of rows.");
	if (Adims[0] != Ddims[0])
		mexErrMsgTxt("Matrix A and matrix B must have the same number of rows.");
	if (Cdims[0] != Cdims[1])
		mexErrMsgTxt("Matrix C must be square.");
	if (Bdims[0] < Bdims[1])
		mexErrMsgTxt("Matrix B must not have more columns than rows.");
	if (Ddims[1] != power(Cdims[0], order))
		mexErrMsgTxt("Matrix D has wrong number of columns.");
}

extern "C" {
	void mexFunction(int nhls, mxArray* plhs[],
					 int nhrs, const mxArray* prhs[])
	{
		if (nhrs != 5)
			mexErrMsgTxt("Must have exactly 5 input args.");
		if (nhls !=1 && nhls != 2)
			mexErrMsgTxt("Must have 1 or 2 output args.");

		int order = (int)mxGetScalar(prhs[0]);
		const mxArray* const A = prhs[1];
		const mxArray* const B = prhs[2];
		const mxArray* const C = prhs[3];
		const mxArray* const D = prhs[4];
		const mwSize* const Adims = mxGetDimensions(A);
		const mwSize* const Bdims = mxGetDimensions(B);
		const mwSize* const Cdims = mxGetDimensions(C);
		const mwSize* const Ddims = mxGetDimensions(D);
		checkDimensions(Adims, Bdims, Cdims, Ddims, order);
		int n = Adims[0];
		int m = Cdims[0];
		int zero_cols = Bdims[0] - Bdims[1];
		mxArray* X = mxCreateDoubleMatrix(Ddims[0], Ddims[1], mxREAL);
		// copy D to X
		Vector Xvec((double*)mxGetPr(X), power(m, order)*n);
		ConstVector Dvec((double*)mxGetPr(D), power(m, order)*n);
		Xvec = Dvec;
		// solve (or solve and check)
		if (nhls == 1) {
			gen_sylv_solve(order, n, m, zero_cols,
						   mxGetPr(A), mxGetPr(B), mxGetPr(C),
						   mxGetPr(X));
		} else if (nhls == 2) {
			gen_sylv_solve_and_check(order, n, m, zero_cols,
									 mxGetPr(A), mxGetPr(B), mxGetPr(C),
									 mxGetPr(D), mxGetPr(X), plhs[1]);
		}
		plhs[0] = X;
	}
};
