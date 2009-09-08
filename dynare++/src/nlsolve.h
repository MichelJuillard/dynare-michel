// Copyright (C) 2006, Ondra Kamenik

// $Id: nlsolve.h 762 2006-05-22 13:00:07Z kamenik $

#ifndef OGU_NLSOLVE_H
#define OGU_NLSOLVE_H

#include "twod_matrix.h"
#include "journal.h"

namespace ogu {

	class OneDFunction {
	public:
		virtual ~OneDFunction() {}
		virtual double eval(double) = 0;
	};

	class GoldenSectionSearch {
	protected:
		static double tol;
		static double golden;
	public:
		static double search(OneDFunction& f, double x1, double x2);
	protected:
		/** This initializes a bracket by moving x2 and b (as a golden
		 * section of x1,x2) so that f(x1)>f(b) && f(b)<f(x2). The point
		 * x1 is not moved, since it is considered as reliable and f(x1)
		 * is supposed to be finite. If initialization of a bracket
		 * succeeded, then [x1, b, x2] is the bracket and true is
		 * returned. Otherwise, b is the minimum found and false is
		 * returned. */
		static bool init_bracket(OneDFunction& f, double x1, double& x2, double& b);
		/** This supposes that f(x1) is finite and it moves x2 toward x1
		 * until x2 and b (as a golden section of x1,x2) are finite. If
		 * succeeded, the routine returns true and x2, and b. Otherwise,
		 * it returns false. */
		static bool search_for_finite(OneDFunction& f, double x1, double& x2, double& b);
	};

	class VectorFunction {
	public:
		VectorFunction() {}
		virtual ~VectorFunction() {}
		virtual int inDim() const = 0;
		virtual int outDim() const = 0;
		/** Check dimensions of eval parameters. */
		void check_for_eval(const ConstVector& in, Vector& out) const;
		/** Evaluate the vector function. */
		virtual void eval(const ConstVector& in, Vector& out) = 0;
	};

	class Jacobian : public TwoDMatrix {
	public:
		Jacobian(int n)
			: TwoDMatrix(n,n) {}
		virtual ~Jacobian() {}
		virtual void eval(const Vector& in) = 0;
	};

	class NLSolver : public OneDFunction {
	protected:
		Journal& journal;
		VectorFunction& func;
		Jacobian& jacob;
		const int max_iter;
		const double tol;
	private:
		Vector xnewton;
		Vector xcauchy;
		Vector x;
	public:
		NLSolver(VectorFunction& f, Jacobian& j, int maxit, double tl, Journal& jr)
			: journal(jr), func(f), jacob(j), max_iter(maxit), tol(tl),
			  xnewton(f.inDim()), xcauchy(f.inDim()), x(f.inDim())
			{xnewton.zeros(); xcauchy.zeros(); x.zeros();}
		virtual ~NLSolver() {}
		/** Returns true if the problem has converged. xx as input is the
		 * starting value, as output it is a solution. */
		bool solve(Vector& xx, int& iter);
		/** To implement OneDFunction interface. It returns
		 * func(xx)^T*func(xx), where
		 * xx=x+lambda*xcauchy+(1-lambda)*xnewton. It is non-const only
		 * because it calls func, x, xnewton, xcauchy is not changed. */
		double eval(double lambda);
	};

};

#endif

// Local Variables:
// mode:C++
// End:
