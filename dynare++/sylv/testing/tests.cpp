/* $Header: /var/lib/cvs/dynare_cpp/sylv/testing/tests.cpp,v 1.2 2004/07/05 19:55:48 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvException.h"
#include "QuasiTriangular.h"
#include "QuasiTriangularZero.h"
#include "Vector.h"
#include "KronVector.h"
#include "KronUtils.h"
#include "TriangularSylvester.h"
#include "GeneralSylvester.h"
#include "SylvMemory.h"
#include "SchurDecompEig.h"
#include "SimilarityDecomp.h"
#include "IterativeSylvester.h"
#include "SylvMatrix.h"

#include "MMMatrix.h"

#include <cstdio>
#include <cstring>
#include <ctime>

#include <cmath>

class TestRunnable : public MallocAllocator {
	char name[100];
	static double eps_norm;
public:
	TestRunnable(const char* n){strncpy(name, n, 100);}
	bool test() const;
	virtual bool run() const =0;
	const char* getName() const {return name;}
protected:
	// declaration of auxiliary static methods
	static bool quasi_solve(bool trans, const char* mname, const char* vname);
	static bool mult_kron(bool trans, const char* mname, const char* vname,
						  const char* cname, int m, int n, int depth);
	static bool level_kron(bool trans, const char* mname, const char* vname,
						   const char* cname, int level, int m, int n, int depth);
	static bool kron_power(const char* m1name, const char* m2name, const char* vname,
						   const char* cname, int m, int n, int depth);
	static bool lin_eval(const char* m1name, const char* m2name, const char* vname,
						 const char* cname, int m, int n, int depth,
						 double alpha, double beta1, double beta2);
	static bool qua_eval(const char* m1name, const char* m2name, const char* vname,
						 const char* cname, int m, int n, int depth,
						 double alpha, double betas, double gamma,
						 double delta1, double delta2);
	static bool tri_sylv(const char* m1name, const char* m2name, const char* vname,
						 int m, int n, int depth);
	static bool gen_sylv(const char* aname, const char* bname, const char* cname,
						 const char* dname, int m, int n, int order);
	static bool eig_bubble(const char* aname, int from, int to);
	static bool block_diag(const char* aname, double log10norm = 3.0);
	static bool iter_sylv(const char* m1name, const char* m2name, const char* vname,
						  int m, int n, int depth);
};

double TestRunnable::eps_norm = 1.0e-10;

bool TestRunnable::test() const
{
	printf("Running test <%s>\n",name);
	clock_t start = clock();
	bool passed = run();
	clock_t end = clock();
	printf("CPU time %8.4g (CPU seconds)..................",
		   ((double)(end-start))/CLOCKS_PER_SEC);
	if (passed) {
		printf("passed\n\n");
		return passed;
	} else {
		printf("FAILED\n\n");
		return passed;
	}
}

/**********************************************************/
/*   auxiliary methods                                    */
/**********************************************************/

bool TestRunnable::quasi_solve(bool trans, const char* mname, const char* vname)
{
	MMMatrixIn mmt(mname);
	MMMatrixIn mmv(vname);

	SylvMemoryDriver memdriver(1, mmt.row(), mmt.row(), 1);
	QuasiTriangular* t;
	QuasiTriangular* tsave;
	if (mmt.row()==mmt.col()) {
		t = new QuasiTriangular(mmt.getData(), mmt.row());
		tsave = new QuasiTriangular(*t);
	} else if (mmt.row()>mmt.col()) {
		t = new QuasiTriangularZero(mmt.row()-mmt.col(), mmt.getData(), mmt.col());
		tsave = new QuasiTriangularZero((const QuasiTriangularZero&)*t);
	} else {
		printf("  Wrong quasi triangular dimensions, rows must be >= cols.\n");
		return false;
	}
	ConstVector v(mmv.getData(), mmv.row());
	Vector x(v.length());
	double eig_min = 1.0e20;
	if (trans)
		t->solveTrans(x, v, eig_min);
	else
		t->solve(x, v, eig_min);
	printf("eig_min = %8.4g\n", eig_min);
	Vector xx(v.length());
	if (trans)
		tsave->multVecTrans(xx, ConstVector(x));
	else
		tsave->multVec(xx, ConstVector(x));
	delete tsave;
	delete t;
	xx.add(-1.0, v);
	xx.add(1.0, x);
	double norm = xx.getNorm();
	printf("\terror norm = %8.4g\n",norm);
	return (norm < eps_norm);
}

bool TestRunnable::mult_kron(bool trans, const char* mname, const char* vname,
							 const char* cname, int m, int n, int depth)
{
	MMMatrixIn mmt(mname);
	MMMatrixIn mmv(vname);
	MMMatrixIn mmc(cname);

	int length = power(m,depth)*n;
	if (mmt.row() != m ||
		mmv.row() != length ||
		mmc.row() != length) {
		printf("  Incompatible sizes for krom mult action, len=%d, matrow=%d, m=%d, vrow=%d, crow=%d \n",length,mmt.row(), m, mmv.row(), mmc.row());
		return false;
	}

	SylvMemoryDriver memdriver(1, m, n, depth);
	QuasiTriangular t(mmt.getData(), mmt.row());
	Vector vraw(mmv.getData(), mmv.row());
	KronVector v(vraw, m, n, depth);
	Vector craw(mmc.getData(), mmc.row());
	KronVector c(craw, m, n, depth);
	if (trans)
		t.multKronTrans(v);
	else
		t.multKron(v);
	c.add(-1.0, v);
	double norm = c.getNorm();
	printf("\terror norm = %8.4g\n",norm);
	return (norm < eps_norm);
}

bool TestRunnable::level_kron(bool trans, const char* mname, const char* vname,
							  const char* cname, int level, int m, int n, int depth)
{
	MMMatrixIn mmt(mname);
	MMMatrixIn mmv(vname);
	MMMatrixIn mmc(cname);

	int length = power(m,depth)*n;
	if (level > 0 && mmt.row() != m ||
		level == 0 && mmt.row() != n ||
		mmv.row() != length ||
		mmc.row() != length) {
		printf("  Incompatible sizes for krom mult action, len=%d, matrow=%d, m=%d, n=%d, vrow=%d, crow=%d \n",length, mmt.row(), m, n, mmv.row(), mmc.row());
		return false;
	}

	SylvMemoryDriver memdriver(1, m, n, depth);
	QuasiTriangular t(mmt.getData(), mmt.row());
	Vector vraw(mmv.getData(), mmv.row());
	ConstKronVector v(vraw, m, n, depth);
	Vector craw(mmc.getData(), mmc.row());
	KronVector c(craw, m, n, depth);
	KronVector x(v);
	if (trans)
		KronUtils::multAtLevelTrans(level, t, x);
	else
		KronUtils::multAtLevel(level, t, x);
	x.add(-1, c);
	double norm = x.getNorm();
	printf("\terror norm = %8.4g\n",norm);
	return (norm < eps_norm);
}

bool TestRunnable::kron_power(const char* m1name, const char* m2name, const char* vname,
							  const char* cname, int m, int n, int depth)
{
	MMMatrixIn mmt1(m1name);
	MMMatrixIn mmt2(m2name);
	MMMatrixIn mmv(vname);
	MMMatrixIn mmc(cname);

	int length = power(m,depth)*n;
	if (mmt1.row() != m ||
		mmt2.row() != n ||
		mmv.row() != length ||
		mmc.row() != length) {
		printf("  Incompatible sizes for krom power mult action, len=%d, row1=%d, row2=%d, m=%d, n=%d, vrow=%d, crow=%d \n",length,mmt1.row(), mmt2.row(), m, n, mmv.row(), mmc.row());
		return false;
	}

	SylvMemoryDriver memdriver(2, m, n, depth);
	QuasiTriangular t1(mmt1.getData(), mmt1.row());
	QuasiTriangular t2(mmt2.getData(), mmt2.row());
	Vector vraw(mmv.getData(), mmv.row());
	ConstKronVector v(vraw, m, n, depth);
	Vector craw(mmc.getData(), mmc.row());
	KronVector c(craw, m, n, depth);
	KronVector x(v);
	memdriver.setStackMode(true);
	KronUtils::multKron(t1, t2, x);
	memdriver.setStackMode(false);
	x.add(-1, c);
	double norm = x.getNorm();
	printf("\terror norm = %8.4g\n",norm);
	return (norm < eps_norm);
}

bool TestRunnable::lin_eval(const char* m1name, const char* m2name, const char* vname,
							const char* cname, int m, int n, int depth,
							double alpha, double beta1, double beta2)
{
	MMMatrixIn mmt1(m1name);
	MMMatrixIn mmt2(m2name);
	MMMatrixIn mmv(vname);
	MMMatrixIn mmc(cname);

	int length = power(m,depth)*n;
	if (mmt1.row() != m ||
		mmt2.row() != n ||
		mmv.row() != 2*length ||
		mmc.row() != 2*length) {
		printf("  Incompatible sizes for lin eval action, len=%d, row1=%d, row2=%d, m=%d, n=%d, vrow=%d, crow=%d \n",length,mmt1.row(), mmt2.row(), m, n, mmv.row(), mmc.row());
		return false;
	}

	SylvMemoryDriver memdriver(1, m, n, depth);
	QuasiTriangular t1(mmt1.getData(), mmt1.row());
	QuasiTriangular t2(mmt2.getData(), mmt2.row());
	TriangularSylvester ts(t2, t1);
	Vector vraw1(mmv.getData(), length);
	ConstKronVector v1(vraw1, m, n, depth);
	Vector vraw2(mmv.getData()+length, length);
	ConstKronVector v2(vraw2, m, n, depth);
	Vector craw1(mmc.getData(), length);
	KronVector c1(craw1, m, n, depth);
	Vector craw2(mmc.getData()+length, length);
	KronVector c2(craw2, m, n, depth);
	KronVector x1(m, n, depth);
	KronVector x2(m, n, depth);
	memdriver.setStackMode(true);
	ts.linEval(alpha, beta1, beta2, x1, x2, v1, v2);
	memdriver.setStackMode(false);
	x1.add(-1, c1);
	x2.add(-1, c2);
	double norm1 = x1.getNorm();
	double norm2 = x2.getNorm();
	printf("\terror norm1 = %8.4g\n\terror norm2 = %8.4g\n",norm1,norm2);
	return (norm1*norm1+norm2*norm2 < eps_norm*eps_norm);
}


bool TestRunnable::qua_eval(const char* m1name, const char* m2name, const char* vname,
							const char* cname, int m, int n, int depth,
							double alpha, double betas, double gamma,
							double delta1, double delta2)
{
	MMMatrixIn mmt1(m1name);
	MMMatrixIn mmt2(m2name);
	MMMatrixIn mmv(vname);
	MMMatrixIn mmc(cname);

	int length = power(m,depth)*n;
	if (mmt1.row() != m ||
		mmt2.row() != n ||
		mmv.row() != 2*length ||
		mmc.row() != 2*length) {
		printf("  Incompatible sizes for qua eval action, len=%d, row1=%d, row2=%d, m=%d, n=%d, vrow=%d, crow=%d \n",length,mmt1.row(), mmt2.row(), m, n, mmv.row(), mmc.row());
		return false;
	}

	SylvMemoryDriver memdriver(3, m, n, depth);
	QuasiTriangular t1(mmt1.getData(), mmt1.row());
	QuasiTriangular t2(mmt2.getData(), mmt2.row());
	TriangularSylvester ts(t2, t1);
	Vector vraw1(mmv.getData(), length);
	ConstKronVector v1(vraw1, m, n, depth);
	Vector vraw2(mmv.getData()+length, length);
	ConstKronVector v2(vraw2, m, n, depth);
	Vector craw1(mmc.getData(), length);
	KronVector c1(craw1, m, n, depth);
	Vector craw2(mmc.getData()+length, length);
	KronVector c2(craw2, m, n, depth);
	KronVector x1(m, n, depth);
	KronVector x2(m, n, depth);
	memdriver.setStackMode(true);
	ts.quaEval(alpha, betas, gamma, delta1, delta2, x1, x2, v1, v2);
	memdriver.setStackMode(false);
	x1.add(-1, c1);
	x2.add(-1, c2);
	double norm1 = x1.getNorm();
	double norm2 = x2.getNorm();
	printf("\terror norm1 = %8.4g\n\terror norm2 = %8.4g\n",norm1,norm2);
	return (norm1*norm1+norm2*norm2 < 100*eps_norm*eps_norm); // relax norm
}

bool TestRunnable::tri_sylv(const char* m1name, const char* m2name, const char* vname,
							int m, int n, int depth)
{
	MMMatrixIn mmt1(m1name);
	MMMatrixIn mmt2(m2name);
	MMMatrixIn mmv(vname);

	int length = power(m,depth)*n;
	if (mmt1.row() != m ||
		mmt2.row() != n ||
		mmv.row() != length) {
		printf("  Incompatible sizes for triangular sylvester action, len=%d, row1=%d, row2=%d, m=%d, n=%d, vrow=%d\n",length,mmt1.row(), mmt2.row(), m, n, mmv.row());
		return false;
	}

	SylvMemoryDriver memdriver(4, m, n, depth); // need extra 2 for checks done via KronUtils::multKron
	memdriver.setStackMode(true);
	QuasiTriangular t1(mmt1.getData(), mmt1.row());
	QuasiTriangular t2(mmt2.getData(), mmt2.row());
	TriangularSylvester ts(t2, t1);
	Vector vraw(mmv.getData(), length);
	ConstKronVector v(vraw, m, n, depth);
	KronVector d(v); // copy of v
	SylvParams pars;
	ts.solve(pars, d);
	pars.print("\t");
	KronVector dcheck((const KronVector&)d);
	KronUtils::multKron(t1, t2, dcheck);
	dcheck.add(1.0, d);
	dcheck.add(-1.0, v);
	double norm = dcheck.getNorm();
	double xnorm = v.getNorm();
	printf("\trel. error norm = %8.4g\n",norm/xnorm);
	double max = dcheck.getMax();
	double xmax = v.getMax();
	printf("\trel. error max = %8.4g\n", max/xmax);
	memdriver.setStackMode(false);
	return (norm < xnorm*eps_norm);
}

bool TestRunnable::gen_sylv(const char* aname, const char* bname, const char* cname,
							const char* dname, int m, int n, int order)
{
	MMMatrixIn mma(aname);
	MMMatrixIn mmb(bname);
	MMMatrixIn mmc(cname);
	MMMatrixIn mmd(dname);

	if (m != mmc.row() || m != mmc.col() ||
		n != mma.row() || n != mma.col() ||
		n != mmb.row() || n <  mmb.col() ||
		n != mmd.row() || power(m, order) != mmd.col()) {
		printf("  Incompatible sizes for gen_sylv.\n");
		return false;
	}

	SylvParams ps(true);
	GeneralSylvester gs(order, n, m, n-mmb.col(),
						mma.getData(), mmb.getData(),
						mmc.getData(), mmd.getData(),
						ps);
	gs.solve();
	gs.check(mmd.getData());
	const SylvParams& pars = gs.getParams();
	pars.print("\t");
	return (*(pars.mat_err1) < eps_norm && *(pars.mat_errI) < eps_norm &&
			*(pars.mat_errF) < eps_norm && *(pars.vec_err1) < eps_norm &&
			*(pars.vec_errI) < eps_norm);
}

bool TestRunnable::eig_bubble(const char* aname, int from, int to)
{
	MMMatrixIn mma(aname);

	if (mma.row() != mma.col()) {
		printf("  Matrix is not square\n");
		return false;
	}

	int n = mma.row();
	SylvMemoryDriver memdriver(3, n, n, 2);
	QuasiTriangular orig(mma.getData(), n);
	SchurDecompEig dec((const QuasiTriangular&)orig);
	QuasiTriangular::diag_iter itf = dec.getT().diag_begin();
	QuasiTriangular::diag_iter itt = dec.getT().diag_begin();
	for (int i = 0; i < from; i++)
		++itf;
	for (int i = 0; i < to; i++)
		++itt;
	itt = dec.bubbleEigen(itf, itt);
	SqSylvMatrix check(dec.getQ(), dec.getT());
	check.multRightTrans(dec.getQ());
	check.add(-1, orig);
	double norm1 = check.getNorm1();
	double normInf = check.getNormInf();
	double onorm1 = orig.getNorm1();
	double onormInf = orig.getNormInf();
	printf("\tabs. error1 = %8.4g\n", norm1);
	printf("\tabs. errorI = %8.4g\n", normInf);
	printf("\trel. error1 = %8.4g\n", norm1/onorm1);
	printf("\trel. errorI = %8.4g\n", normInf/onormInf);
	return (norm1 < eps_norm*onorm1 && normInf < eps_norm*onormInf);
}

bool TestRunnable::block_diag(const char* aname, double log10norm)
{
	MMMatrixIn mma(aname);

	if (mma.row() != mma.col()) {
		printf("  Matrix is not square\n");
		return false;
	}

	int n = mma.row();
	SylvMemoryDriver memdriver(3, n, n, 2);
	SqSylvMatrix orig(mma.getData(), n);
	SimilarityDecomp dec(orig.base(), orig.numRows(), log10norm);
	dec.getB().printInfo();
	SqSylvMatrix check(dec.getQ(), dec.getB());
	check.multRight(dec.getInvQ());
	check.add(-1, orig);
	double norm1 = check.getNorm1();
	double normInf = check.getNormInf();
	double onorm1 = orig.getNorm1();
	double onormInf = orig.getNormInf();
	printf("\terror Q*B*invQ:\n");
	printf("\tabs. error1 = %8.4g\n", norm1);
	printf("\tabs. errorI = %8.4g\n", normInf);
	printf("\trel. error1 = %8.4g\n", norm1/onorm1);
	printf("\trel. errorI = %8.4g\n", normInf/onormInf);
	SqSylvMatrix check2(dec.getQ(), dec.getInvQ());
	SqSylvMatrix in(n);
	in.setUnit();
	check2.add(-1, in);
	double nor1 = check2.getNorm1();
	double norInf = check2.getNormInf();
	printf("\terror Q*invQ:\n");
	printf("\tabs. error1 = %8.4g\n", nor1);
	printf("\tabs. errorI = %8.4g\n", norInf);
	return (norm1 < eps_norm*pow(10, log10norm)*onorm1);
}

bool TestRunnable::iter_sylv(const char* m1name, const char* m2name, const char* vname,
							 int m, int n, int depth)
{
	MMMatrixIn mmt1(m1name);
	MMMatrixIn mmt2(m2name);
	MMMatrixIn mmv(vname);

	int length = power(m,depth)*n;
	if (mmt1.row() != m ||
		mmt2.row() != n ||
		mmv.row() != length) {
		printf("  Incompatible sizes for triangular sylvester iteration, len=%d, row1=%d, row2=%d, m=%d, n=%d, vrow=%d\n",length,mmt1.row(), mmt2.row(), m, n, mmv.row());
		return false;
	}

	SylvMemoryDriver memdriver(4, m, n, depth); // need extra 2 for checks done via KronUtils::multKron
	memdriver.setStackMode(true);
	QuasiTriangular t1(mmt1.getData(), mmt1.row());
	QuasiTriangular t2(mmt2.getData(), mmt2.row());
	IterativeSylvester is(t2, t1);
	Vector vraw(mmv.getData(), length);
	ConstKronVector v(vraw, m, n, depth);
	KronVector d(v); // copy of v
	SylvParams pars;
	pars.method = SylvParams::iter;
	is.solve(pars, d);
	pars.print("\t");
	KronVector dcheck((const KronVector&)d);
	KronUtils::multKron(t1, t2, dcheck);
	dcheck.add(1.0, d);
	dcheck.add(-1.0, v);
	double cnorm = dcheck.getNorm();
	double xnorm = v.getNorm();
	printf("\trel. error norm = %8.4g\n",cnorm/xnorm);
	double max = dcheck.getMax();
	double xmax = v.getMax();
	printf("\trel. error max = %8.4g\n", max/xmax);
	memdriver.setStackMode(false);
	return (cnorm < xnorm*eps_norm);
}

/**********************************************************/
/*   sub classes declarations                             */
/**********************************************************/

class PureTriangTest : public TestRunnable {
public:
	PureTriangTest() : TestRunnable("pure triangular solve (5)") {}
	bool run() const;
};

class PureTriangTransTest : public TestRunnable {
public:
	PureTriangTransTest() : TestRunnable("pure triangular solve trans (5)") {}
	bool run() const;
};

class PureTrLargeTest : public TestRunnable {
public:
	PureTrLargeTest() : TestRunnable("pure triangular large solve (300)") {}
	bool run() const;
};

class PureTrLargeTransTest : public TestRunnable {
public:
	PureTrLargeTransTest() : TestRunnable("pure triangular large solve trans (300)") {}
	bool run() const;
};

class QuasiTriangTest : public TestRunnable {
public:
	QuasiTriangTest() : TestRunnable("quasi triangular solve (7)") {}
	bool run() const;
};

class QuasiTriangTransTest : public TestRunnable {
public:
	QuasiTriangTransTest() : TestRunnable("quasi triangular solve trans (7)") {}
	bool run() const;
};

class QuasiTrLargeTest : public TestRunnable {
public:
	QuasiTrLargeTest() : TestRunnable("quasi triangular solve large (250)") {}
	bool run() const;
};

class QuasiTrLargeTransTest : public TestRunnable {
public:
	QuasiTrLargeTransTest() : TestRunnable("quasi triangular solve large trans (250)") {}
	bool run() const;
};

class QuasiZeroSmallTest : public TestRunnable {
public:
	QuasiZeroSmallTest() : TestRunnable("quasi tr. zero small test (2x1)") {}
	bool run() const;
};

class MultKronSmallTest : public TestRunnable {
public:
	MultKronSmallTest() : TestRunnable("kronecker small mult (2=2x1)") {}
	bool run() const;
};

class MultKronTest : public TestRunnable {
public:
	MultKronTest() : TestRunnable("kronecker mult (245=7x7x5)") {}
	bool run() const;
};

class MultKronSmallTransTest : public TestRunnable {
public:
	MultKronSmallTransTest() : TestRunnable("kronecker small trans mult (2=2x1)") {}
	bool run() const;
};

class MultKronTransTest : public TestRunnable {
public:
	MultKronTransTest() : TestRunnable("kronecker trans mult (245=7x7x5)") {}
	bool run() const;
};

class LevelKronTest : public TestRunnable {
public:
	LevelKronTest() : TestRunnable("kronecker level mult (1715=7x[7]x7x5)") {}
	bool run() const;
};

class LevelKronTransTest : public TestRunnable {
public:
	LevelKronTransTest() : TestRunnable("kronecker level trans mult (1715=7x[7]x7x5)") {}
	bool run() const;
};

class LevelZeroKronTest : public TestRunnable {
public:
	LevelZeroKronTest() : TestRunnable("kronecker level mult (1715=7x7x7x[5])") {}
	bool run() const;
};

class LevelZeroKronTransTest : public TestRunnable {
public:
	LevelZeroKronTransTest() : TestRunnable("kronecker level trans mult (1715=7x7x7x[5])") {}
	bool run() const;
};

class KronPowerTest : public TestRunnable {
public:
	KronPowerTest() : TestRunnable("kronecker power mult (1715=7x7x7x5)") {}
	bool run() const;
};

class SmallLinEvalTest : public TestRunnable {
public:
	SmallLinEvalTest() : TestRunnable("lin eval (24=2 x 2x2x3)") {}
	bool run() const;
};

class LinEvalTest : public TestRunnable {
public:
	LinEvalTest() : TestRunnable("lin eval (490=2 x 7x7x5)") {}
	bool run() const;
};

class SmallQuaEvalTest : public TestRunnable {
public:
	SmallQuaEvalTest() : TestRunnable("qua eval (24=2 x 2x2x3)") {}
	bool run() const;
};

class QuaEvalTest : public TestRunnable {
public:
	QuaEvalTest() : TestRunnable("qua eval (490=2 x 7x7x5)") {}
	bool run() const;
};

class TriSylvSmallRealTest : public TestRunnable {
public:
	TriSylvSmallRealTest() : TestRunnable("triangular sylvester small real solve (12=2x2x3)") {}
	bool run() const;
};

class TriSylvSmallComplexTest : public TestRunnable {
public:
	TriSylvSmallComplexTest() : TestRunnable("triangular sylvester small complx solve (12=2x2x3)") {}
	bool run() const;
};

class TriSylvTest : public TestRunnable {
public:
	TriSylvTest() : TestRunnable("triangular sylvester solve (245=7x7x5)") {}
	bool run() const;
};

class TriSylvBigTest : public TestRunnable {
public:
	TriSylvBigTest() : TestRunnable("triangular sylvester big solve (48000=40x40x30)") {}
	bool run() const;
};

class TriSylvLargeTest : public TestRunnable {
public:
	TriSylvLargeTest() : TestRunnable("triangular sylvester large solve (1920000=40x40x40x30)") {}
	bool run() const;
};

class IterSylvTest : public TestRunnable {
public:
	IterSylvTest() : TestRunnable("iterative sylvester solve (245=7x7x5)") {}
	bool run() const;
};

class IterSylvLargeTest : public TestRunnable {
public:
	IterSylvLargeTest() : TestRunnable("iterative sylvester large solve (1920000=40x40x40x30)") {}
	bool run() const;
};

class GenSylvSmallTest : public TestRunnable {
public:
	GenSylvSmallTest() : TestRunnable("general sylvester small solve (18=3x3x2)") {}
	bool run() const;
};

class GenSylvTest : public TestRunnable {
public:
	GenSylvTest() : TestRunnable("general sylvester solve (12000=20x20x30)") {}
	bool run() const;
};

class GenSylvSingTest : public TestRunnable {
public:
	GenSylvSingTest() : TestRunnable("general sylvester solve for sing. C (2500000=50x50x50x20)") {}
	bool run() const;
};

class GenSylvLargeTest : public TestRunnable {
public:
	GenSylvLargeTest() : TestRunnable("general sylvester solve (2500000=50x50x50x20)") {}
	bool run() const;
};

class EigBubFrankTest : public TestRunnable {
public:
	EigBubFrankTest() : TestRunnable("eig. bubble frank test (12x12)") {}
	bool run() const;
};

class EigBubSplitTest : public TestRunnable {
// complex eigenvalue is split by swapping it with real
public:
	EigBubSplitTest() : TestRunnable("eig. bubble complex split test (3x3)") {}
	bool run() const;
};

class EigBubSameTest : public TestRunnable {
// complex eigenevalue bypasses the same complex eigenvalue
public:
	EigBubSameTest() : TestRunnable("eig. bubble same test (5x5)") {}
	bool run() const;
};

class BlockDiagSmallTest : public TestRunnable {
public:
	BlockDiagSmallTest() : TestRunnable("block diagonalization small test (7x7)") {}
	bool run() const;
};

class BlockDiagFrankTest : public TestRunnable {
public:
	BlockDiagFrankTest() : TestRunnable("block diagonalization of frank (12x12)") {}
	bool run() const;
};

class BlockDiagIllCondTest : public TestRunnable {
public:
	BlockDiagIllCondTest() : TestRunnable("block diagonalization of ill conditioned (15x15)") {}
	bool run() const;
};

class BlockDiagBigTest : public TestRunnable {
public:
	BlockDiagBigTest() : TestRunnable("block diagonalization big test (50x50)") {}
	bool run() const;
};

/**********************************************************/
/*   run methods of sub classes                           */
/**********************************************************/

bool PureTriangTest::run() const
{
	return quasi_solve(false, "tr5x5.mm", "v5.mm");
}

bool PureTriangTransTest::run() const
{
	return quasi_solve(true, "tr5x5.mm", "v5.mm");
}

bool PureTrLargeTest::run() const
{
	return quasi_solve(false, "tr300x300.mm", "v300.mm");
}

bool PureTrLargeTransTest::run() const
{
	return quasi_solve(true, "tr300x300.mm", "v300.mm");
}

bool QuasiTriangTest::run() const
{
	return quasi_solve(false, "qt7x7.mm", "v7.mm");
}

bool QuasiTriangTransTest::run() const
{
	return quasi_solve(true, "qt7x7.mm", "v7.mm");
}

bool QuasiTrLargeTest::run() const
{
	return quasi_solve(false, "qt250x250.mm", "v250.mm");
}

bool QuasiTrLargeTransTest::run() const
{
	return quasi_solve(true, "qt250x250.mm", "v250.mm");
}

bool QuasiZeroSmallTest::run() const
{
	return quasi_solve(false, "b2x1.mm", "v2.mm");
}

bool MultKronSmallTest::run() const
{
	return mult_kron(false, "tr2x2.mm", "v2.mm", "vcheck2.mm", 2, 1, 1);
}

bool MultKronTest::run() const
{
	return mult_kron(false, "qt7x7.mm", "v245.mm", "vcheck245.mm", 7, 5, 2);
}

bool MultKronSmallTransTest::run() const
{
	return mult_kron(true, "tr2x2.mm", "v2.mm", "vcheck2a.mm", 2, 1, 1);
}

bool MultKronTransTest::run() const
{
	return mult_kron(true, "qt7x7.mm", "v245.mm", "vcheck245a.mm", 7, 5, 2);
}

bool LevelKronTest::run() const
{
	return level_kron(false, "qt7x7.mm", "v1715.mm", "vcheck1715.mm", 2, 7, 5, 3);
}

bool LevelKronTransTest::run() const
{
	return level_kron(true, "qt7x7.mm", "v1715.mm", "vcheck1715a.mm", 2, 7, 5, 3);
}

bool LevelZeroKronTest::run() const
{
	return level_kron(false, "tr5x5.mm", "v1715.mm", "vcheck1715b.mm", 0, 7, 5, 3);
}

bool LevelZeroKronTransTest::run() const
{
	return level_kron(true, "tr5x5.mm", "v1715.mm", "vcheck1715c.mm", 0, 7, 5, 3);
}

bool KronPowerTest::run() const
{
	return kron_power("qt7x7.mm", "tr5x5.mm", "v1715.mm", "vcheck1715d.mm", 7, 5, 3);
}

bool SmallLinEvalTest::run() const
{
	return lin_eval("qt2x2.mm", "qt3x3.mm", "v24.mm", "vcheck24.mm", 2, 3, 2,
					 2, 1, 3);
}

bool LinEvalTest::run() const
{
	return lin_eval("qt7x7.mm", "tr5x5.mm", "v490.mm", "vcheck490.mm", 7, 5, 2,
					 2, 1, 3);
}

bool SmallQuaEvalTest::run() const
{
	return qua_eval("qt2x2.mm", "qt3x3.mm", "v24.mm", "vcheck24q.mm", 2, 3, 2,
					 -0.5, 3, 2, 1, 3);
}

bool QuaEvalTest::run() const
{
	return qua_eval("qt7x7.mm", "tr5x5.mm", "v490.mm", "vcheck490q.mm", 7, 5, 2,
					 -0.5, 3, 2, 1, 3);
}

bool TriSylvSmallRealTest::run() const
{
	return tri_sylv("tr2x2.mm", "qt3x3.mm", "v12r.mm", 2, 3, 2);
}

bool TriSylvSmallComplexTest::run() const
{
	return tri_sylv("qt2x2.mm", "qt3x3.mm", "v12r.mm", 2, 3, 2);
}

bool TriSylvTest::run() const
{
	return tri_sylv("qt7x7eig06-09.mm", "tr5x5.mm", "v245r.mm", 7, 5, 2);
}

bool TriSylvBigTest::run() const
{
	return tri_sylv("qt40x40.mm", "qt30x30eig011-095.mm", "v48000.mm", 40, 30, 2);
}

bool TriSylvLargeTest::run() const
{
	return tri_sylv("qt40x40.mm", "qt30x30eig011-095.mm", "v1920000.mm", 40, 30, 3);
}

bool IterSylvTest::run() const
{
	return iter_sylv("qt7x7eig06-09.mm", "qt5x5.mm", "v245r.mm", 7, 5, 2);
}

bool IterSylvLargeTest::run() const
{
	return iter_sylv("qt40x40.mm", "qt30x30eig011-095.mm", "v1920000.mm", 40, 30, 3);
}

bool GenSylvSmallTest::run() const
{
	return gen_sylv("a2x2.mm", "b2x1.mm", "c3x3.mm", "d2x9.mm", 3, 2, 2);
}

bool GenSylvTest::run() const
{
	return gen_sylv("a30x30.mm", "b30x25.mm", "c20x20.mm", "d30x400.mm", 20, 30, 2);
}

bool GenSylvSingTest::run() const
{
	return gen_sylv("a20x20.mm", "b20x4.mm", "c50x50sing.mm", "d20x125000.mm", 50, 20, 3);
}

bool GenSylvLargeTest::run() const
{
	return gen_sylv("a20x20.mm", "b20x15.mm", "c50x50.mm", "d20x125000.mm", 50, 20, 3);
}

bool EigBubFrankTest::run() const
{
	return eig_bubble("qt_frank12x12.mm", 8, 0);
}

bool EigBubSplitTest::run() const
{
	return eig_bubble("qt_eps3x3.mm",1,0);
}

bool EigBubSameTest::run() const
{
	return eig_bubble("qt5x5.mm",2,0);
}

bool BlockDiagSmallTest::run() const
{
	return block_diag("qt7x7.mm", 0.1);
}

bool BlockDiagFrankTest::run() const
{
	return block_diag("qt_frank12x12.mm", 5);
}

bool BlockDiagIllCondTest::run() const
{
	return block_diag("ill_cond15x15.mm", 4.14);
}

bool BlockDiagBigTest::run() const
{
	return block_diag("c50x50.mm", 1.3);
}

/**********************************************************/
/*   main                                                 */
/**********************************************************/

int main()
{
	TestRunnable* all_tests[50];
	// fill in vector of all tests
	int num_tests = 0;
	all_tests[num_tests++] = new PureTriangTest();
	all_tests[num_tests++] = new PureTriangTransTest();
	all_tests[num_tests++] = new PureTrLargeTest();
	all_tests[num_tests++] = new PureTrLargeTransTest();
	all_tests[num_tests++] = new QuasiTriangTest();
	all_tests[num_tests++] = new QuasiTriangTransTest();
	all_tests[num_tests++] = new QuasiTrLargeTest();
	all_tests[num_tests++] = new QuasiTrLargeTransTest();
	all_tests[num_tests++] = new QuasiZeroSmallTest();
	all_tests[num_tests++] = new MultKronSmallTest();
	all_tests[num_tests++] = new MultKronTest();
	all_tests[num_tests++] = new MultKronSmallTransTest();
	all_tests[num_tests++] = new MultKronTransTest();
	all_tests[num_tests++] = new LevelKronTest();
	all_tests[num_tests++] = new LevelKronTransTest();
	all_tests[num_tests++] = new LevelZeroKronTest();
	all_tests[num_tests++] = new LevelZeroKronTransTest();
	all_tests[num_tests++] = new KronPowerTest();
	all_tests[num_tests++] = new SmallLinEvalTest();
	all_tests[num_tests++] = new LinEvalTest();
	all_tests[num_tests++] = new SmallQuaEvalTest();
	all_tests[num_tests++] = new QuaEvalTest();
	all_tests[num_tests++] = new EigBubFrankTest();
	all_tests[num_tests++] = new EigBubSplitTest();
	all_tests[num_tests++] = new EigBubSameTest();
	all_tests[num_tests++] = new BlockDiagSmallTest();
	all_tests[num_tests++] = new BlockDiagFrankTest();
	all_tests[num_tests++] = new BlockDiagIllCondTest();
	all_tests[num_tests++] = new BlockDiagBigTest();
	all_tests[num_tests++] = new TriSylvSmallRealTest();
	all_tests[num_tests++] = new TriSylvSmallComplexTest();
	all_tests[num_tests++] = new TriSylvTest();
	all_tests[num_tests++] = new TriSylvBigTest();
	all_tests[num_tests++] = new TriSylvLargeTest();
	all_tests[num_tests++] = new IterSylvTest();
	all_tests[num_tests++] = new IterSylvLargeTest();
	all_tests[num_tests++] = new GenSylvSmallTest();
	all_tests[num_tests++] = new GenSylvTest();
	all_tests[num_tests++] = new GenSylvSingTest();
	all_tests[num_tests++] = new GenSylvLargeTest();

	// launch the tests
	int success = 0;
	for (int i = 0; i < num_tests; i++) {
		try {
			if (all_tests[i]->test())
				success++;
		} catch (const MMException& e) {
			printf("Caugth MM exception in <%s>:\n%s", all_tests[i]->getName(),
				   e.getMessage());
		} catch (SylvException& e) {
			printf("Caught Sylv exception in %s:\n", all_tests[i]->getName());
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

