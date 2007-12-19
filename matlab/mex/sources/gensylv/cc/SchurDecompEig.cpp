#include "SchurDecompEig.h"
#include "SylvException.h"
#include "cpplapack.h"

/* bubble diagonal 1-1, or 2-2 block from position 'from' to position
 * 'to'. If an eigenvalue cannot be swapped with its neighbour, the
 * neighbour is bubbled also in front. The method returns a new
 * position 'to', where the original block pointed by 'to' happens to
 * appear at the end. 'from' must be greater than 'to'.
 */
SchurDecompEig::diag_iter
SchurDecompEig::bubbleEigen(diag_iter from, diag_iter to)
{
	diag_iter run = from;
	while (run != to) {
		diag_iter runm = run;
		if (!tryToSwap(run, runm) && runm == to) {
			++to;
		} else {
			// bubble all eigenvalues from runm(incl.) to run(excl.),
			// this includes either bubbling generated eigenvalues due
			// to split, or an eigenvalue which couldn't be swapped
			while (runm != run) {
				to = bubbleEigen(runm, to);
				++runm;
			}
		}
	}
	return to;
}

/* this tries to swap two neighbouring eigenvalues, 'it' and '--it',
 * and returns 'itadd'. If the blocks can be swapped, new eigenvalues
 * can emerge due to possible 2-2 block splits. 'it' then points to
 * the last eigenvalue coming from block pointed by 'it' at the
 * begining, and 'itadd' points to the first. On swap failure, 'it' is
 * not changed, and 'itadd' points to previous eignevalue (which must
 * be moved backwards before). In either case, it is necessary to
 * resolve eigenvalues from 'itadd' to 'it', before the 'it' can be
 * resolved.
 * The success is signaled by returned true.
 */
bool SchurDecompEig::tryToSwap(diag_iter& it, diag_iter& itadd)
{
	itadd = it;
	--itadd;

	int n = getDim();
	int ifst = (*it).getIndex() + 1;
	int ilst = (*itadd).getIndex() + 1;
	double* work = new double[n];
	int info;
	LAPACK_dtrexc("V", &n, getT().base(), &n, getQ().base(), &n, &ifst, &ilst, work,
				  &info);
	delete [] work;
	if (info < 0) {
		throw SYLV_MES_EXCEPTION("Wrong argument to LAPACK_dtrexc.");
	}

	if (info == 0) {
		// swap successful
		getT().swapDiagLogically(itadd);
		//check for 2-2 block splits
		getT().checkDiagConsistency(it);
		getT().checkDiagConsistency(itadd);
		// and go back by 'it' in NEW eigenvalue set
		--it;
		return true;
	}
	return false;
}


void SchurDecompEig::orderEigen()
{
	diag_iter run = getT().diag_begin();
	diag_iter runp = run;
	++runp;
	double last_size = 0.0;
	while (runp != getT().diag_end()) {
		diag_iter least = getT().findNextLargerBlock(run, getT().diag_end(),
													 last_size);
		last_size = (*least).getSize();
		if (run == least)
			++run;
		else
			run = bubbleEigen(least, run);
		runp = run;
		++runp;
	}
}
