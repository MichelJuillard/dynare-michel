/* $Id: tests.cpp 431 2005-08-16 15:41:01Z kamenik $ */
/* Copyright 2005, Ondra Kamenik */

#include "GeneralMatrix.h"
#include <dynlapack.h>
#include "SylvException.h"

#include "rfs_tensor.h"
#include "normal_moments.h"

#include "vector_function.h"
#include "quadrature.h"
#include "smolyak.h"
#include "product.h"
#include "quasi_mcarlo.h"

#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include <cmath>

const int num_threads = 2; // does nothing if DEBUG defined

// evaluates unfolded (Dx)^k power, where x is a vector, D is a
// Cholesky factor (lower triangular)
class MomentFunction : public VectorFunction {
	GeneralMatrix D;
	int k;
public:
	MomentFunction(const GeneralMatrix& inD, int kk)
		: VectorFunction(inD.numRows(), UFSTensor::calcMaxOffset(inD.numRows(), kk)),
		  D(inD), k(kk) {}
	MomentFunction(const MomentFunction& func)
		: VectorFunction(func), D(func.D), k(func.k) {}
	VectorFunction* clone() const
		{return new MomentFunction(*this);}
	void eval(const Vector& point, const ParameterSignal& sig, Vector& out);
};

void MomentFunction::eval(const Vector& point, const ParameterSignal& sig, Vector& out)
{
	if (point.length() != indim() || out.length() != outdim()) {
		printf("Wrong length of vectors in MomentFunction::eval\n");
		exit(1);
	}
	Vector y(point);
	y.zeros();
	D.multaVec(y, point);
	URSingleTensor ypow(y, k);
	out.zeros();
	out.add(1.0, ypow.getData());
}

class TensorPower : public VectorFunction {
	int k;
public:
	TensorPower(int nvar, int kk)
		: VectorFunction(nvar, UFSTensor::calcMaxOffset(nvar, kk)), k(kk) {}
	TensorPower(const TensorPower& func)
		: VectorFunction(func), k(func.k) {}
	VectorFunction* clone() const
		{return new TensorPower(*this);}
	void eval(const Vector& point, const ParameterSignal& sig, Vector& out);
};

void TensorPower::eval(const Vector& point, const ParameterSignal& sig, Vector& out)
{
	if (point.length() != indim() || out.length() != outdim()) {
		printf("Wrong length of vectors in TensorPower::eval\n");
		exit(1);
	}
	URSingleTensor ypow(point, k);
	out.zeros();
	out.add(1.0, ypow.getData());
}


// evaluates (1+1/d)^d*(x_1*...*x_d)^(1/d), its integral over <0,1>^d
// is 1.0, and its variation grows exponetially
// d = dim
class Function1 : public VectorFunction {
	int dim;
public:
	Function1(int d)
		: VectorFunction(d, 1), dim(d) {}
	Function1(const Function1& f)
		: VectorFunction(f.indim(), f.outdim()), dim(f.dim) {}
	VectorFunction* clone() const
		{return new Function1(*this);}
	virtual void eval(const Vector& point, const ParameterSignal& sig, Vector& out);
};

void Function1::eval(const Vector& point, const ParameterSignal& sig, Vector& out)
{
	if (point.length() != dim || out.length() != 1) {
		printf("Wrong length of vectors in Function1::eval\n");
		exit(1);
	}
	double r = 1;
	for (int i = 0; i < dim; i++)
		r *= point[i];
	r = pow(r, 1.0/dim);
	r *= pow(1.0 + 1.0/dim, (double)dim);
	out[0] = r;
}

// evaluates Function1 but with transformation x_i=0.5(y_i+1)
// this makes the new function integrate over <-1,1>^d to 1.0
class Function1Trans : public Function1 {
public:
	Function1Trans(int d)
		: Function1(d) {}
	Function1Trans(const Function1Trans& func)
		: Function1(func) {}
	VectorFunction* clone() const
		{return new Function1Trans(*this);}
	virtual void eval(const Vector& point, const ParameterSignal& sig, Vector& out);
};

void Function1Trans::eval(const Vector& point, const ParameterSignal& sig, Vector& out)
{
	Vector p(point.length());
	for (int i = 0; i < p.length(); i++)
		p[i] = 0.5*(point[i]+1);
	Function1::eval(p, sig, out);
	out.mult(pow(0.5,indim()));
}


// WallTimer class. Constructor saves the wall time, destructor
// cancels the current time from the saved, and prints the message
// with time information
class WallTimer {
	char mes[100];
	struct timeval start;
	bool new_line;
public:
	WallTimer(const char* m, bool nl = true)
		{strcpy(mes, m);new_line = nl; gettimeofday(&start, NULL);}
	~WallTimer()
		{
			struct timeval end;
			gettimeofday(&end, NULL);
			printf("%s%8.4g", mes,
				   end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)*1.0e-6);
			if (new_line)
				printf("\n");
		}
};

/****************************************************/
/*     declaration of TestRunnable class            */
/****************************************************/
class TestRunnable {
	char name[100];
public:
	int dim; // dimension of the solved problem
	int nvar; // number of variable of the solved problem
	TestRunnable(const char* n, int d, int nv)
		: dim(d), nvar(nv)
		{strncpy(name, n, 100);}
	bool test() const;
	virtual bool run() const =0;
	const char* getName() const
		{return name;}
protected:
	static bool smolyak_normal_moments(const GeneralMatrix& m, int imom, int level);
	static bool product_normal_moments(const GeneralMatrix& m, int imom, int level);
	static bool qmc_normal_moments(const GeneralMatrix& m, int imom, int level);
	static bool smolyak_product_cube(const VectorFunction& func, const Vector& res,
									 double tol, int level);
	static bool qmc_cube(const VectorFunction& func, double res, double tol, int level);
};

bool TestRunnable::test() const
{
	printf("Running test <%s>\n",name);
	bool passed;
	{
		WallTimer tim("Wall clock time ", false);
		passed = run();
	}
	if (passed) {
		printf("............................ passed\n\n");
		return passed;
	} else {
		printf("............................ FAILED\n\n");
		return passed;
	}
}


/****************************************************/
/*     definition of TestRunnable static methods    */
/****************************************************/
bool TestRunnable::smolyak_normal_moments(const GeneralMatrix& m, int imom, int level)
{
	// first make m*m' and then Cholesky factor
	GeneralMatrix mtr(m, "transpose");
	GeneralMatrix msq(m, mtr);

	// make vector function
	int dim = m.numRows();
	TensorPower tp(dim, imom);
	GaussConverterFunction func(tp, msq);

	// smolyak quadrature
	Vector smol_out(UFSTensor::calcMaxOffset(dim, imom));
	{
		WallTimer tim("\tSmolyak quadrature time:         ");
		GaussHermite gs;
		SmolyakQuadrature quad(dim, level, gs);
		quad.integrate(func, level, num_threads, smol_out);
		printf("\tNumber of Smolyak evaluations:    %d\n", quad.numEvals(level));
	}

	// check against theoretical moments
	UNormalMoments moments(imom, msq);
	smol_out.add(-1.0, (moments.get(Symmetry(imom)))->getData());
	printf("\tError:                         %16.12g\n", smol_out.getMax());
	return smol_out.getMax() < 1.e-7;
}

bool TestRunnable::product_normal_moments(const GeneralMatrix& m, int imom, int level)
{
	// first make m*m' and then Cholesky factor
	GeneralMatrix mtr(m, "transpose");
	GeneralMatrix msq(m, mtr);

	// make vector function
	int dim = m.numRows();
	TensorPower tp(dim, imom);
	GaussConverterFunction func(tp, msq);

	// product quadrature
	Vector prod_out(UFSTensor::calcMaxOffset(dim, imom));
	{
		WallTimer tim("\tProduct quadrature time:         ");
		GaussHermite gs;
		ProductQuadrature quad(dim, gs);
		quad.integrate(func, level, num_threads, prod_out);
		printf("\tNumber of product evaluations:    %d\n", quad.numEvals(level));
	}

	// check against theoretical moments
	UNormalMoments moments(imom, msq);
	prod_out.add(-1.0, (moments.get(Symmetry(imom)))->getData());
	printf("\tError:                         %16.12g\n", prod_out.getMax());
	return prod_out.getMax() < 1.e-7;
}

bool TestRunnable::qmc_normal_moments(const GeneralMatrix& m, int imom, int level)
{
	// first make m*m' and then Cholesky factor
	GeneralMatrix mtr(m, "transpose");
	GeneralMatrix msq(m, mtr);
	GeneralMatrix mchol(msq);
	int rows = mchol.numRows();
	for (int i = 0; i < rows; i++)
		for (int j = i+1; j < rows; j++)
			mchol.get(i,j) = 0.0;
	int info;
	dpotrf("L", &rows, mchol.base(), &rows, &info);

	// make vector function
	MomentFunction func(mchol, imom);

	// permutation schemes
	WarnockPerScheme wps;
	ReversePerScheme rps;
	IdentityPerScheme ips;
	PermutationScheme* scheme[] = {&wps, &rps, &ips};
	const char* labs[] = {"Warnock", "Reverse", "Identity"};

	// theoretical result
	int dim = mchol.numRows();
	UNormalMoments moments(imom, msq);
	Vector res((const Vector&)((moments.get(Symmetry(imom)))->getData()));

	// quasi monte carlo normal quadrature
	double max_error = 0.0;
	Vector qmc_out(UFSTensor::calcMaxOffset(dim, imom));
	for (int i = 0; i < 3; i++) {
		{
			char mes[100];
			sprintf(mes, "\tQMC normal quadrature time %8s:         ", labs[i]);
			WallTimer tim(mes);
			QMCarloNormalQuadrature quad(dim, level, *(scheme[i]));
			quad.integrate(func, level, num_threads, qmc_out);
		}
		qmc_out.add(-1.0, res);
		printf("\tError %8s:                         %16.12g\n", labs[i], qmc_out.getMax());
		if (qmc_out.getMax() > max_error) {
			max_error = qmc_out.getMax();
		}
	}

	return max_error < 1.e-7;
}


bool TestRunnable::smolyak_product_cube(const VectorFunction& func, const Vector& res,
										double tol, int level)
{
	if (res.length() != func.outdim()) {
		fprintf(stderr, "Incompatible dimensions of check value and function.\n");
		exit(1);
	}

	GaussLegendre glq;
	Vector out(func.outdim());
	double smol_error;
	double prod_error;
	{
		WallTimer tim("\tSmolyak quadrature time:         ");
		SmolyakQuadrature quad(func.indim(), level, glq);
		quad.integrate(func, level, num_threads, out);
		out.add(-1.0, res);
		smol_error = out.getMax();
		printf("\tNumber of Smolyak evaluations:    %d\n", quad.numEvals(level));
		printf("\tError:                            %16.12g\n", smol_error);
	}
	{
		WallTimer tim("\tProduct quadrature time:         ");
		ProductQuadrature quad(func.indim(), glq);
		quad.integrate(func, level, num_threads, out);
		out.add(-1.0, res);
		prod_error = out.getMax();
		printf("\tNumber of product evaluations:    %d\n", quad.numEvals(level));
		printf("\tError:                            %16.12g\n", prod_error);
	}

	return smol_error < tol && prod_error < tol;
}

bool TestRunnable::qmc_cube(const VectorFunction& func, double res, double tol, int level)
{
	Vector r(1);
	double error1;
	{
		WallTimer tim("\tQuasi-Monte Carlo (Warnock scrambling) time:  ");
		WarnockPerScheme wps;
		QMCarloCubeQuadrature qmc(func.indim(), level, wps);
//		qmc.savePoints("warnock.txt", level);
		qmc.integrate(func, level, num_threads, r);
		error1 = std::max(res - r[0], r[0] - res);
		printf("\tQuasi-Monte Carlo (Warnock scrambling) error: %16.12g\n",
			   error1);
	}
	double error2;
	{
		WallTimer tim("\tQuasi-Monte Carlo (reverse scrambling) time:  ");
		ReversePerScheme rps;
		QMCarloCubeQuadrature qmc(func.indim(), level, rps);
//		qmc.savePoints("reverse.txt", level);
		qmc.integrate(func, level, num_threads, r);
		error2 = std::max(res - r[0], r[0] - res);
		printf("\tQuasi-Monte Carlo (reverse scrambling) error: %16.12g\n",
			   error2);
	}
	double error3;
	{
		WallTimer tim("\tQuasi-Monte Carlo (no scrambling) time:       ");
		IdentityPerScheme ips;
		QMCarloCubeQuadrature qmc(func.indim(), level, ips);
//		qmc.savePoints("identity.txt", level);
		qmc.integrate(func, level, num_threads, r);
		error3 = std::max(res - r[0], r[0] - res);
		printf("\tQuasi-Monte Carlo (no scrambling) error:      %16.12g\n",
			   error3);
	}


	return error1 < tol && error2 < tol && error3 < tol;
}

/****************************************************/
/*     definition of TestRunnable subclasses        */
/****************************************************/
class SmolyakNormalMom1 : public TestRunnable {
public:
	SmolyakNormalMom1()
		: TestRunnable("Smolyak normal moments (dim=2, level=4, order=4)", 4, 2) {}

	bool run() const
		{
			GeneralMatrix m(2,2);
			m.zeros(); m.get(0,0)=1; m.get(1,1)=1;
			return smolyak_normal_moments(m, 4, 4);
		}
};

class SmolyakNormalMom2 : public TestRunnable {
public:
	SmolyakNormalMom2()
		: TestRunnable("Smolyak normal moments (dim=3, level=8, order=8)", 8, 3) {}

	bool run() const
		{
			GeneralMatrix m(3,3);
			m.zeros();
			m.get(0,0)=1; m.get(0,2)=0.5; m.get(1,1)=1;
			m.get(1,0)=0.5;m.get(2,2)=2;m.get(2,1)=4;
			return smolyak_normal_moments(m, 8, 8);
		}
};

class ProductNormalMom1 : public TestRunnable {
public:
	ProductNormalMom1()
		: TestRunnable("Product normal moments (dim=2, level=4, order=4)", 4, 2) {}

	bool run() const
		{
			GeneralMatrix m(2,2);
			m.zeros(); m.get(0,0)=1; m.get(1,1)=1;
			return product_normal_moments(m, 4, 4);
		}
};

class ProductNormalMom2 : public TestRunnable {
public:
	ProductNormalMom2()
		: TestRunnable("Product normal moments (dim=3, level=8, order=8)", 8, 3) {}

	bool run() const
		{
			GeneralMatrix m(3,3);
			m.zeros();
			m.get(0,0)=1; m.get(0,2)=0.5; m.get(1,1)=1;
			m.get(1,0)=0.5;m.get(2,2)=2;m.get(2,1)=4;
			return product_normal_moments(m, 8, 8);
		}
};

class QMCNormalMom1 : public TestRunnable {
public:
	QMCNormalMom1()
		: TestRunnable("QMC normal moments (dim=2, level=1000, order=4)", 4, 2) {}

	bool run() const
		{
			GeneralMatrix m(2,2);
			m.zeros(); m.get(0,0)=1; m.get(1,1)=1;
			return qmc_normal_moments(m, 4, 1000);
		}
};

class QMCNormalMom2 : public TestRunnable {
public:
	QMCNormalMom2()
		: TestRunnable("QMC normal moments (dim=3, level=10000, order=8)", 8, 3) {}

	bool run() const
		{
			GeneralMatrix m(3,3);
			m.zeros();
			m.get(0,0)=1; m.get(0,2)=0.5; m.get(1,1)=1;
			m.get(1,0)=0.5;m.get(2,2)=2;m.get(2,1)=4;
			return qmc_normal_moments(m, 8, 10000);
		}
};



// note that here we pass 1,1 to tls since smolyak has its own PascalTriangle
class F1GaussLegendre : public TestRunnable {
public:
	F1GaussLegendre()
		: TestRunnable("Function1 Gauss-Legendre (dim=6, level=13", 1, 1) {}

	bool run() const
		{
			Function1Trans f1(6);
			Vector res(1); res[0] = 1.0;
			return smolyak_product_cube(f1, res, 1e-2, 13);
		}
};


class F1QuasiMCarlo : public TestRunnable {
public:
	F1QuasiMCarlo()
		: TestRunnable("Function1 Quasi-Monte Carlo (dim=6, level=1000000)", 1, 1) {}

	bool run() const
		{
			Function1 f1(6);
			return qmc_cube(f1, 1.0, 1.e-4, 1000000);
		}
};

int main()
{
	TestRunnable* all_tests[50];
	// fill in vector of all tests
	int num_tests = 0;
	all_tests[num_tests++] = new SmolyakNormalMom1();
	all_tests[num_tests++] = new SmolyakNormalMom2();
	all_tests[num_tests++] = new ProductNormalMom1();
	all_tests[num_tests++] = new ProductNormalMom2();
	all_tests[num_tests++] = new QMCNormalMom1();
	all_tests[num_tests++] = new QMCNormalMom2();
/*
	all_tests[num_tests++] = new F1GaussLegendre();
	all_tests[num_tests++] = new F1QuasiMCarlo();
*/
	// find maximum dimension and maximum nvar
	int dmax=0;
	int nvmax = 0;
	for (int i = 0; i < num_tests; i++) {
		if (dmax < all_tests[i]->dim)
			dmax = all_tests[i]->dim;
		if (nvmax < all_tests[i]->nvar)
			nvmax = all_tests[i]->nvar;
	}
	tls.init(dmax, nvmax); // initialize library
	THREAD_GROUP::max_parallel_threads = num_threads;

	// launch the tests
	int success = 0;
	for (int i = 0; i < num_tests; i++) {
		try {
			if (all_tests[i]->test())
				success++;
		} catch (const TLException& e) {
			printf("Caugth TL exception in <%s>:\n", all_tests[i]->getName());
			e.print();
		} catch (SylvException& e) {
			printf("Caught Sylv exception in <%s>:\n", all_tests[i]->getName());
			e.printMessage();
		}
	}

	printf("There were %d tests that failed out of %d tests run.\n",
		   num_tests - success, num_tests);

	// destroy
	for (int i = 0; i < num_tests; i++) {
		delete all_tests[i];
	}

	return 0;
}
