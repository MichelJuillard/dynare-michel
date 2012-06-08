// Copyright (C) 2006, Ondra Kamenik

// $Id: nlsolve.cpp 762 2006-05-22 13:00:07Z kamenik $

#include "nlsolve.h"
#include "dynare_exception.h"

#include <cmath>

using namespace ogu;

/** This should not be greater than DBL_EPSILON^(1/2). */
double GoldenSectionSearch::tol = 1.e-4;

/** This is equal to the golden section ratio. */
double GoldenSectionSearch::golden = (3.-std::sqrt(5.))/2;

double GoldenSectionSearch::search(OneDFunction& f, double x1, double x2)
{
	double b;
	if (init_bracket(f, x1, x2, b)) {
		double fb = f.eval(b);
		double f1 = f.eval(x1);
		f.eval(x2);
		double dx;
		do {
			double w = (b-x1)/(x2-x1);
			dx = std::abs((1-2*w)*(x2-x1));
			double x;
			if (b-x1 > x2-b)
				x = b - dx;
			else
				x = b + dx;
			double fx = f.eval(x);
			if (! std::isfinite(fx))
				return x1;
			if (b-x1 > x2-b) {
				// x is on the left from b
				if (f1 > fx && fx < fb) {
					// pickup bracket [f1,fx,fb]
					x2 = b;
					fb = fx;
					b = x;
				} else {
					// pickup bracket [fx,fb,fx2]
					f1 = fx;
					x1 = x;
				}
			} else {
				// x is on the right from b
				if (f1 > fb && fb < fx) {
					// pickup bracket [f1,fb,fx]
					x2 = x;
				} else {
					// pickup bracket [fb,fx,f2]
					f1 = fb;
					x1 = b;
					fb = fx;
					b = x;
				}
			}
		} while(dx > tol);
	}
	return b;
}

bool GoldenSectionSearch::init_bracket(OneDFunction& f, double x1, double& x2, double& b)
{
	double f1 = f.eval(x1);
	if (! std::isfinite(f1))
		throw DynareException(__FILE__, __LINE__,
							  "Safer point not finite in GoldenSectionSearch::init_bracket");

	int cnt = 0;
	bool bracket_found = false;
	do {
		bool finite_found = search_for_finite(f, x1, x2, b);
		if (! finite_found) {
			b = x1;
			return false;
		}
		double f2 = f.eval(x2);
		double fb = f.eval(b);
		double bsym = 2*x2 - b;
		double fbsym = f.eval(bsym);
		// now we know that f1, f2, and fb are finite
		if (std::isfinite(fbsym)) {
			// we have four numbers f1, fb, f2, fbsym, we test for the
			// following combinations to find the bracket:
			// [f1,f2,fbsym], [f1,fb,fbsym] and [f1,fb,fbsym]
			if (f1 > f2 && f2 < fbsym) {
				bracket_found = true;
				b = x2;
				x2 = bsym;
			} else if (f1 > fb && fb < fbsym) {
				bracket_found = true;
				x2 = bsym;
			} else if (f1 > fb && fb < f2) {
				bracket_found = true;
			} else {
				double newx2 = b;
				// choose the smallest value in case we end
				if (f1 > fbsym) {
					// the smallest value is on the other end, we do
					// not want to continue
					b = bsym;
					return false;
				} else
					b = x1;
					// move x2 to b in case we continue 
					x2 = newx2;
			}
		} else {
			// we have only three numbers, we test for the bracket,
			// and if not found, we set b as potential result and
			// shorten x2 as potential init value for next cycle
			if (f1 > fb && fb < f2)
				bracket_found = true;
			else {
				double newx2 = b;
				// choose the smaller value in case we end
				if (f1 > f2)
					b = x2;
				else
					b = x1;
				// move x2 to b in case we continue
				x2 = newx2;
			}
		}
		cnt++;
	} while (! bracket_found && cnt < 5);
	
	return bracket_found;
}

/** This moves x2 toward to x1 until the function at x2 is finite and
 * b as a golden section between x1 and x2 yields also finite f. */
bool GoldenSectionSearch::search_for_finite(OneDFunction& f, double x1, double& x2, double&b)
{
	int cnt = 0;
	bool found = false;
	do {
		double f2 = f.eval(x2);
		b = (1-golden)*x1 + golden*x2;
		double fb = f.eval(b);
		found = std::isfinite(f2) && std::isfinite(fb);
		if (! found)
			x2 = b;
		cnt++;
	} while (! found && cnt < 5);

	return found;
}

void VectorFunction::check_for_eval(const ConstVector& in, Vector& out) const
{
	if (inDim() != in.length() || outDim() != out.length())
		throw DynareException(__FILE__, __LINE__,
							  "Wrong dimensions in VectorFunction::check_for_eval");
}

double NLSolver::eval(double lambda)
{
	Vector xx((const Vector&)x);
	xx.add(1-lambda, xcauchy);
	xx.add(lambda, xnewton);
	Vector ff(func.outDim());
	func.eval(xx, ff);
	return ff.dot(ff);
}

bool NLSolver::solve(Vector& xx, int& iter)
{
	JournalRecord rec(journal);
	rec << "Iter   lambda      residual" << endrec;
	JournalRecord rec1(journal);
	rec1 << "---------------------------" << endrec;
	char tmpbuf[14];

	x = (const Vector&)xx;
	iter = 0;
	// setup fx
	Vector fx(func.outDim());
	func.eval(x, fx);
	if (!fx.isFinite())
		throw DynareException(__FILE__,__LINE__,
							  "Initial guess does not yield finite residual in NLSolver::solve");
	bool converged = fx.getMax() < tol;
	JournalRecord rec2(journal);
	sprintf(tmpbuf, "%10.6g", fx.getMax());
	rec2 << iter << "         N/A   " << tmpbuf << endrec;
	while (! converged && iter < max_iter) {
		// setup Jacobian
		jacob.eval(x);
		// calculate cauchy step
		Vector g(func.inDim());
		g.zeros();
		ConstTwoDMatrix(jacob).multaVecTrans(g, fx);
		Vector Jg(func.inDim());
		Jg.zeros();
		ConstTwoDMatrix(jacob).multaVec(Jg, g);
		double m = -g.dot(g)/Jg.dot(Jg);
		xcauchy = (const Vector&) g;
		xcauchy.mult(m);
		// calculate newton step
		xnewton = (const Vector&) fx;
		ConstTwoDMatrix(jacob).multInvLeft(xnewton);
		xnewton.mult(-1);

		// line search
		double lambda = GoldenSectionSearch::search(*this, 0, 1);
		x.add(1-lambda, xcauchy);
		x.add(lambda, xnewton);
		// evaluate func
		func.eval(x, fx);
		converged = fx.getMax() < tol;

		// iter
		iter++;

		JournalRecord rec3(journal);
		sprintf(tmpbuf, "%10.6g", fx.getMax());
		rec3 << iter << "    " << lambda << "   " << tmpbuf << endrec;
	}
	xx = (const Vector&)x;

	return converged;
}
