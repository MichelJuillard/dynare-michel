/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/GeneralSylvester.cpp,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#include "GeneralSylvester.h"
#include "SchurDecomp.h"
#include "SylvException.h"
#include "TriangularSylvester.h"
#include "IterativeSylvester.h"

#include <ctime>

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
								   const double* da, const double* db,
								   const double* dc, const double* dd,
								   const SylvParams& ps)
	: pars(ps), 
	  mem_driver(pars, 1, m, n, ord), order(ord), a(da, n),
	  b(db, n, n-zero_cols), c(dc, m), d(dd, n, power(m, order)),
	  solved(false)
{
	init();
}

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
								   const double* da, const double* db,
								   const double* dc, double* dd,
								   const SylvParams& ps)
	: pars(ps),
	  mem_driver(pars, 0, m, n, ord), order(ord), a(da, n),
	  b(db, n, n-zero_cols), c(dc, m), d(dd, n, power(m, order)),
	  solved(false)
{
	init();
}

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
								   const double* da, const double* db,
								   const double* dc, const double* dd,
								   bool alloc_for_check)
	: pars(alloc_for_check), 
	  mem_driver(pars, 1, m, n, ord), order(ord), a(da, n),
	  b(db, n, n-zero_cols), c(dc, m), d(dd, n, power(m, order)),
	  solved(false)
{
	init();
}

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
								   const double* da, const double* db,
								   const double* dc, double* dd,
								   bool alloc_for_check)
	: pars(alloc_for_check),
	  mem_driver(pars, 0, m, n, ord), order(ord), a(da, n),
	  b(db, n, n-zero_cols), c(dc, m), d(dd, n, power(m, order)),
	  solved(false)
{
	init();
}

void GeneralSylvester::init()
{
	GeneralMatrix ainvb(b);
	double rcond1;
	double rcondinf;
	a.multInvLeft2(ainvb, d, rcond1, rcondinf);
	pars.rcondA1 = rcond1;
	pars.rcondAI = rcondinf;
	bdecomp = new SchurDecompZero(ainvb);
	cdecomp = new SimilarityDecomp(c.getData().base(), c.numRows(), *(pars.bs_norm));
	cdecomp->check(pars, c);
	cdecomp->infoToPars(pars);
	if (*(pars.method) == SylvParams::recurse)
		sylv = new TriangularSylvester(*bdecomp, *cdecomp);
	else
		sylv = new IterativeSylvester(*bdecomp, *cdecomp);
}

void GeneralSylvester::solve()
{
	if (solved)
		throw SYLV_MES_EXCEPTION("Attempt to run solve() more than once.");

	mem_driver.setStackMode(true);

	clock_t start = clock();
	// multiply d
	d.multLeftITrans(bdecomp->getQ());
	d.multRightKron(cdecomp->getQ(), order);
	// convert to KronVector
	KronVector dkron(d.getData(), getM(), getN(), order);
	// solve
	sylv->solve(pars, dkron);
	// multiply d back
	d.multLeftI(bdecomp->getQ());
	d.multRightKron(cdecomp->getInvQ(), order);
	clock_t end = clock();
	pars.cpu_time = ((double)(end-start))/CLOCKS_PER_SEC;

	mem_driver.setStackMode(false);

	solved = true;
}

void GeneralSylvester::check(const double* ds)
{
	if (!solved)
		throw SYLV_MES_EXCEPTION("Cannot run check on system, which is not solved yet.");

	mem_driver.setStackMode(true);

	// calculate xcheck = AX+BXC^i-D
	SylvMatrix dcheck(d.numRows(), d.numCols());
	dcheck.multLeft(b.numRows()-b.numCols(), b, d);
	dcheck.multRightKron(c, order);
	dcheck.multAndAdd(a,d);
	ConstVector dv(ds, d.numRows()*d.numCols());
	dcheck.getData().add(-1.0, dv);
	// calculate relative norms
	pars.mat_err1 = dcheck.getNorm1()/d.getNorm1();
	pars.mat_errI = dcheck.getNormInf()/d.getNormInf();
	pars.mat_errF = dcheck.getData().getNorm()/d.getData().getNorm();
	pars.vec_err1 = dcheck.getData().getNorm1()/d.getData().getNorm1();
	pars.vec_errI = dcheck.getData().getMax()/d.getData().getMax();

	mem_driver.setStackMode(false);
}

GeneralSylvester::~GeneralSylvester()
{
	delete bdecomp;
	delete cdecomp;
	delete sylv;
}

// Local Variables:
// mode:C++
// End:
