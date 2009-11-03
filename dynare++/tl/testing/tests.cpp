/* $Id: tests.cpp 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#include "SylvException.h"
#include "tl_exception.h"
#include "gs_tensor.h"
#include "factory.h"
#include "monoms.h"
#include "t_container.h"
#include "stack_container.h"
#include "t_polynomial.h"
#include "rfs_tensor.h"
#include "ps_tensor.h"
#include "tl_static.h"

#include <cstdio>
#include <cstring>
#include <ctime>


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
	template<class _Ttype>
	static bool index_forward(const Symmetry& s, const IntSequence& nvs);

	template <class _Ttype>
	static bool index_backward(const Symmetry& s, const IntSequence& nvs);

	template <class _Ttype>
	static bool index_offset(const Symmetry& s, const IntSequence& nvs);

	static bool fold_unfold(const FTensor* folded);
	static bool fs_fold_unfold(int r, int nv, int dim)
		{
			Factory f;
			FTensor* folded = f.make<FFSTensor>(r, nv, dim);
			return fold_unfold(folded); // folded deallocated in fold_unfold
		}
	static bool r_fold_unfold(int r, int nv, int dim)
		{
			Factory f;
			FTensor* folded = f.make<FRTensor>(r, nv, dim);
			return fold_unfold(folded); // folded deallocated in fold_unfold
		}
	static bool gs_fold_unfold(int r, const Symmetry& s, const IntSequence& nvs)
		{
			Factory f;
			FTensor* folded = f.make<FGSTensor>(r, s, nvs);
			return fold_unfold(folded); // folded deallocated in fold_unfold
		}

	static bool dense_prod(const Symmetry& bsym, const IntSequence& bnvs,
						   int hdim, int hnv, int rows);

	static bool folded_monomial(int ng, int nx, int ny, int nu, int dim);

	static bool unfolded_monomial(int ng, int nx, int ny, int nu, int dim);

	static bool fold_zcont(int nf, int ny, int nu, int nup, int nbigg,
						   int ng, int dim);

	static bool unfold_zcont(int nf, int ny, int nu, int nup, int nbigg,
							 int ng, int dim);

	static bool folded_contraction(int r, int nv, int dim);

	static bool unfolded_contraction(int r, int nv, int dim);

	static bool poly_eval(int r, int nv, int maxdim);


};

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


/****************************************************/
/*     definition of TestRunnable static methods    */
/****************************************************/
template <class _Ttype>
bool TestRunnable::index_forward(const Symmetry& s, const IntSequence& nvs)
{
	int fails = 0;
	int ndecr = 0;
	int nincr = 0;
	_Ttype dummy(0, TensorDimens(s, nvs));
	typename _Ttype::index run = dummy.end();
	do {
		--run;
		ndecr++;
		typename _Ttype::index run2 = dummy.begin();
		for (int i = 0; i < *run; i++) {
			++run2;
			nincr++;
		}
		if (! (run == run2))
			fails++;
	} while (run != dummy.begin());

	printf("\tnumber of columns    = %d\n",dummy.ncols());
	printf("\tnumber of increments = %d\n",nincr);
	printf("\tnumber of decrements = %d\n",ndecr);
	printf("\tnumber of failures   = %d\n",fails);

	return fails == 0;
}

template <class _Ttype>
bool TestRunnable::index_backward(const Symmetry& s, const IntSequence& nvs)
{
	int fails = 0;
	int ndecr = 0;
	int nincr = 0;
	_Ttype dummy(0, TensorDimens(s, nvs));
	typename _Ttype::index run = dummy.begin();
	while (run != dummy.end()) {
		typename _Ttype::index run2 = dummy.end();
		for (int i = 0; i < dummy.ncols() - *run; i++) {
			--run2;
			ndecr++;
		}
		if (! (run == run2))
			fails++;
		++run;
		nincr++;
	}

	printf("\tnumber of columns    = %d\n",dummy.ncols());
	printf("\tnumber of increments = %d\n",nincr);
	printf("\tnumber of decrements = %d\n",ndecr);
	printf("\tnumber of failures   = %d\n",fails);

	return fails == 0;
}

template <class _Ttype>
bool TestRunnable::index_offset(const Symmetry& s, const IntSequence& nvs)
{
	int fails = 0;
	int nincr = 0;
	_Ttype dummy(0, TensorDimens(s, nvs));
	for (typename _Ttype::index run = dummy.begin();
		 run != dummy.end(); ++run, nincr++) {
		typename _Ttype::index run2(&dummy, run.getCoor());
		if (! (run == run2))
			fails++;
	}

	printf("\tnumber of columns    = %d\n",dummy.ncols());
	printf("\tnumber of increments = %d\n",nincr);
	printf("\tnumber of failures   = %d\n",fails);

	return fails == 0;
}

bool TestRunnable::fold_unfold(const FTensor* folded)
{
	UTensor* unfolded = &(folded->unfold());
	FTensor* folded2 = &(unfolded->fold());
	folded2->add(-1.0, *folded);
	double normInf = folded2->getNormInf();
	double norm1 = folded2->getNorm1();
	printf("\tfolded size:       (%d, %d)\n",folded->nrows(), folded->ncols());
	printf("\tunfolded size:     (%d, %d)\n",unfolded->nrows(), unfolded->ncols());
	printf("\tdifference normInf: %8.4g\n", normInf);
	printf("\tdifference norm1:   %8.4g\n", norm1);

	delete folded;
	delete unfolded;
	delete folded2;

	return normInf < 1.0e-15;
}

bool TestRunnable::dense_prod(const Symmetry& bsym, const IntSequence& bnvs,
							  int hdim, int hnv, int rows)
{
	Factory f;
	FGSContainer* cont =
		f.makeCont<FGSTensor,FGSContainer>(hnv, bnvs, bsym.dimen()-hdim+1);
	FGSTensor* fh =
		f.make<FGSTensor>(rows, Symmetry(hdim), IntSequence(1, hnv));
	UGSTensor uh(*fh);
	FGSTensor fb(rows, TensorDimens(bsym, bnvs));
	fb.getData().zeros();
	clock_t s1 = clock();
	cont->multAndAdd(uh, fb);
	clock_t s2 = clock();
	UGSContainer ucont(*cont);
	clock_t s3 = clock();
	UGSTensor ub(rows, fb.getDims());
	ub.getData().zeros();
	clock_t s4 = clock();
	ucont.multAndAdd(uh, ub);
	clock_t s5 = clock();

	UGSTensor btmp(fb);
	btmp.add(-1, ub);
	double norm = btmp.getData().getMax();
	double norm1 = btmp.getNorm1();
	double normInf = btmp.getNormInf();

	printf("\ttime for folded product:     %8.4g\n",
		   ((double)(s2-s1))/CLOCKS_PER_SEC);
	printf("\ttime for unfolded product:   %8.4g\n",
		   ((double)(s5-s4))/CLOCKS_PER_SEC);
	printf("\ttime for container convert:  %8.4g\n",
		   ((double)(s3-s2))/CLOCKS_PER_SEC);
	printf("\tunfolded difference normMax: %10.6g\n", norm);
	printf("\tunfolded difference norm1:   %10.6g\n", norm1);
	printf("\tunfolded difference normInf: %10.6g\n", normInf);

	delete cont;
	delete fh;

	return norm < 1.e-13;
}

bool TestRunnable::folded_monomial(int ng, int nx, int ny, int nu, int dim)
{
	clock_t gen_time = clock();
	DenseDerivGenerator gen(ng, nx, ny, nu, 5, 0.3, dim);
	gen_time = clock()-gen_time;
	printf("\ttime for monom generation: %8.4g\n",
		   ((double)gen_time)/CLOCKS_PER_SEC);
	IntSequence nvs(2); nvs[0] = ny; nvs[1] = nu;
	double maxnorm = 0;
	for (int ydim = 0; ydim <= dim; ydim++) {
		Symmetry s(ydim, dim-ydim);
		printf("\tSymmetry: ");s.print();
		FGSTensor res(ng, TensorDimens(s, nvs));
		res.getData().zeros();
		clock_t stime = clock();
		for (int d = 1; d <= dim; d++) {
			gen.xcont->multAndAdd(*(gen.ts[d-1]), res);
		}
		stime = clock() - stime;
		printf("\t\ttime for symmetry: %8.4g\n",
			   ((double)stime)/CLOCKS_PER_SEC);
		const FGSTensor* mres = gen.rcont->get(s);
		res.add(-1.0, *mres);
		double normtmp = res.getData().getMax();
		printf("\t\terror normMax:     %10.6g\n", normtmp);
		if (normtmp > maxnorm)
			maxnorm = normtmp;
	}
	return maxnorm < 1.0e-10;
}

bool TestRunnable::unfolded_monomial(int ng, int nx, int ny, int nu, int dim)
{
	clock_t gen_time = clock();
	DenseDerivGenerator gen(ng, nx, ny, nu, 5, 0.3, dim);
	gen_time = clock()-gen_time;
	printf("\ttime for monom generation: %8.4g\n",
		   ((double)gen_time)/CLOCKS_PER_SEC);
	clock_t u_time = clock();
	gen.unfold();
	u_time = clock() - u_time;
	printf("\ttime for monom unfolding:  %8.4g\n",
		   ((double)u_time)/CLOCKS_PER_SEC);
	IntSequence nvs(2); nvs[0] = ny; nvs[1] = nu;
	double maxnorm = 0;
	for (int ydim = 0; ydim <= dim; ydim++) {
		Symmetry s(ydim, dim-ydim);
		printf("\tSymmetry: ");s.print();
		UGSTensor res(ng, TensorDimens(s, nvs));
		res.getData().zeros();
		clock_t stime = clock();
		for (int d = 1; d <= dim; d++) {
			gen.uxcont->multAndAdd(*(gen.uts[d-1]), res);
		}
		stime = clock() - stime;
		printf("\t\ttime for symmetry: %8.4g\n",
			   ((double)stime)/CLOCKS_PER_SEC);
		const FGSTensor* mres = gen.rcont->get(s);
		FGSTensor foldres(res);
		foldres.add(-1.0, *mres);
		double normtmp = foldres.getData().getMax();
		printf("\t\terror normMax:     %10.6g\n", normtmp);
		if (normtmp > maxnorm)
			maxnorm = normtmp;
	}
	return maxnorm < 1.0e-10;
}

bool TestRunnable::fold_zcont(int nf, int ny, int nu, int nup, int nbigg,
							  int ng, int dim)
{
	clock_t gen_time = clock();
	SparseDerivGenerator dg(nf, ny, nu, nup, nbigg, ng,
							5, 0.55, dim);
	gen_time = clock()-gen_time;
	for (int d = 1; d <= dim; d++) {
		printf("\tfill of dim=%d tensor:     %3.2f %%\n",
			   d, 100*dg.ts[d-1]->getFillFactor());
	}
	printf("\ttime for monom generation: %8.4g\n",
		   ((double)gen_time)/CLOCKS_PER_SEC);

	IntSequence nvs(4);
	nvs[0] = ny; nvs[1] = nu; nvs[2] = nup; nvs[3] = 1;
	double maxnorm = 0.0;

	// form ZContainer
	FoldedZContainer zc(dg.bigg, nbigg, dg.g, ng, ny, nu);

	for (int d = 2; d <= dim; d++) {
		SymmetrySet ss(d, 4);
		for (symiterator si(ss); !si.isEnd(); ++si) {
			printf("\tSymmetry: ");(*si).print();
			FGSTensor res(nf, TensorDimens(*si, nvs));
			res.getData().zeros();
			clock_t stime = clock();
			for (int l = 1; l <= (*si).dimen(); l++) {
				zc.multAndAdd(*(dg.ts[l-1]), res);
			}
			stime = clock() - stime;
			printf("\t\ttime for symmetry: %8.4g\n",
				   ((double)stime)/CLOCKS_PER_SEC);
			const FGSTensor* mres = dg.rcont->get(*si);
			res.add(-1.0, *mres);
			double normtmp = res.getData().getMax();
			printf("\t\terror normMax:     %10.6g\n", normtmp);
			if (normtmp > maxnorm)
				maxnorm = normtmp;
		}
	}
	return maxnorm < 1.0e-10;
}

bool TestRunnable::unfold_zcont(int nf, int ny, int nu, int nup, int nbigg,
								int ng, int dim)
{
	clock_t gen_time = clock();
	SparseDerivGenerator dg(nf, ny, nu, nup, nbigg, ng,
							5, 0.55, dim);
	gen_time = clock()-gen_time;
	for (int d = 1; d <= dim; d++) {
		printf("\tfill of dim=%d tensor:     %3.2f %%\n",
			   d, 100*dg.ts[d-1]->getFillFactor());
	}
	printf("\ttime for monom generation: %8.4g\n",
		   ((double)gen_time)/CLOCKS_PER_SEC);

	clock_t con_time = clock();
	UGSContainer uG_cont(*(dg.bigg));
	UGSContainer ug_cont(*(dg.g));
	con_time = clock()-con_time;
	printf("\ttime for container unfold: %8.4g\n",
		   ((double)con_time)/CLOCKS_PER_SEC);

	IntSequence nvs(4);
	nvs[0] = ny; nvs[1] = nu; nvs[2] = nup; nvs[3] = 1;
	double maxnorm = 0.0;

	// form ZContainer
	UnfoldedZContainer zc(&uG_cont, nbigg, &ug_cont, ng, ny, nu);

	for (int d = 2; d <= dim; d++) {
		SymmetrySet ss(d, 4);
		for (symiterator si(ss); !si.isEnd(); ++si) {
			printf("\tSymmetry: ");(*si).print();
			UGSTensor res(nf, TensorDimens(*si, nvs));
			res.getData().zeros();
			clock_t stime = clock();
			for (int l = 1; l <= (*si).dimen(); l++) {
				zc.multAndAdd(*(dg.ts[l-1]), res);
			}
			stime = clock() - stime;
			printf("\t\ttime for symmetry: %8.4g\n",
				   ((double)stime)/CLOCKS_PER_SEC);
			FGSTensor fold_res(res);
			const FGSTensor* mres = dg.rcont->get(*si);
			fold_res.add(-1.0, *mres);
			double normtmp = fold_res.getData().getMax();
			printf("\t\terror normMax:     %10.6g\n", normtmp);
			if (normtmp > maxnorm)
				maxnorm = normtmp;
		}
	}
	return maxnorm < 1.0e-10;
}

bool TestRunnable::folded_contraction(int r, int nv, int dim)
{
	Factory fact;
	Vector* x = fact.makeVector(nv);

	FFSTensor* forig = fact.make<FFSTensor>(r, nv, dim);
	FFSTensor* f = new FFSTensor(*forig);
	clock_t ctime = clock();
	for (int d = dim-1; d > 0; d--) {
		FFSTensor* fnew = new FFSTensor(*f, ConstVector(*x));
		delete f;
		f = fnew;
	}
	ctime = clock() - ctime;
	Vector res(forig->nrows());
	res.zeros();
	f->multaVec(res, *x);

	UFSTensor u(*forig);
	clock_t utime = clock();
	URSingleTensor ux(*x, dim);
	Vector v(u.nrows());
	v.zeros();
	u.multaVec(v, ux.getData());
	utime = clock() - utime;

	v.add(-1.0, res);
	printf("\ttime for folded contraction: %8.4g\n",
		   ((double)ctime)/CLOCKS_PER_SEC);
	printf("\ttime for unfolded power:     %8.4g\n",
		   ((double)utime)/CLOCKS_PER_SEC);
	printf("\terror normMax:     %10.6g\n", v.getMax());
	printf("\terror norm1:       %10.6g\n", v.getNorm1());

	delete f;
	delete x;

	return (v.getMax() < 1.e-10);
}

bool TestRunnable::unfolded_contraction(int r, int nv, int dim)
{
	Factory fact;
	Vector* x = fact.makeVector(nv);

	FFSTensor* forig = fact.make<FFSTensor>(r, nv, dim);
	UFSTensor uorig(*forig);
	delete forig;
	UFSTensor* u = new UFSTensor(uorig);
	clock_t ctime = clock();
	for (int d = dim-1; d > 0; d--) {
		UFSTensor* unew = new UFSTensor(*u, ConstVector(*x));
		delete u;
		u = unew;
	}
	ctime = clock() - ctime;
	Vector res(uorig.nrows());
	res.zeros();
	u->multaVec(res, *x);

	clock_t utime = clock();
	URSingleTensor ux(*x, dim);
	Vector v(uorig.nrows());
	v.zeros();
	uorig.multaVec(v, ux.getData());
	utime = clock() - utime;

	v.add(-1.0, res);
	printf("\ttime for unfolded contraction: %8.4g\n",
		   ((double)ctime)/CLOCKS_PER_SEC);
	printf("\ttime for unfolded power:       %8.4g\n",
		   ((double)utime)/CLOCKS_PER_SEC);
	printf("\terror normMax:     %10.6g\n", v.getMax());
	printf("\terror norm1:       %10.6g\n", v.getNorm1());

	delete u;
	delete x;

	return (v.getMax() < 1.e-10);
}

bool TestRunnable::poly_eval(int r, int nv, int maxdim)
{
	Factory fact;
	Vector* x = fact.makeVector(nv);

	Vector out_ft(r); out_ft.zeros();
	Vector out_fh(r); out_fh.zeros();
	Vector out_ut(r); out_ut.zeros();
	Vector out_uh(r); out_uh.zeros();

	UTensorPolynomial* up;
	{
		FTensorPolynomial* fp = fact.makePoly<FFSTensor, FTensorPolynomial>(r, nv, maxdim);

		clock_t ft_cl = clock();
		fp->evalTrad(out_ft, *x);
		ft_cl = clock() - ft_cl;
		printf("\ttime for folded power eval:    %8.4g\n",
			   ((double)ft_cl)/CLOCKS_PER_SEC);
		
		clock_t fh_cl = clock();
		fp->evalHorner(out_fh, *x);
		fh_cl = clock() - fh_cl;
		printf("\ttime for folded horner eval:   %8.4g\n",
			   ((double)fh_cl)/CLOCKS_PER_SEC);

		up = new UTensorPolynomial(*fp);
		delete fp;
	}

	clock_t ut_cl = clock();
	up->evalTrad(out_ut, *x);
	ut_cl = clock() - ut_cl;
	printf("\ttime for unfolded power eval:  %8.4g\n",
		   ((double)ut_cl)/CLOCKS_PER_SEC);

	clock_t uh_cl = clock();
	up->evalHorner(out_uh, *x);
	uh_cl = clock() - uh_cl;
	printf("\ttime for unfolded horner eval: %8.4g\n",
		   ((double)uh_cl)/CLOCKS_PER_SEC);

	out_ft.add(-1.0, out_ut);
	double max_ft = out_ft.getMax();
	out_fh.add(-1.0, out_ut);
	double max_fh = out_fh.getMax();
	out_uh.add(-1.0, out_ut);
	double max_uh = out_uh.getMax();

	printf("\tfolded power error norm max:     %10.6g\n", max_ft);
	printf("\tfolded horner error norm max:    %10.6g\n", max_fh);
	printf("\tunfolded horner error norm max:  %10.6g\n", max_uh);

	delete up;
	delete x;
	return (max_ft+max_fh+max_uh < 1.0e-10);
}


/****************************************************/
/*     definition of TestRunnable subclasses        */
/****************************************************/
class SmallIndexForwardFold : public TestRunnable {
public:
	SmallIndexForwardFold()
		: TestRunnable("small index forward for fold (44)(222)", 5, 4) {}
	bool run() const
		{
			Symmetry s(2,3);
			IntSequence nvs(2); nvs[0] = 4; nvs[1] = 2;
			return index_forward<FGSTensor>(s, nvs);
		}
};

class SmallIndexForwardUnfold : public TestRunnable {
public:
	SmallIndexForwardUnfold()
		: TestRunnable("small index forward for unfold (44)(222)", 5, 4) {}
	bool run() const
		{
			Symmetry s(2,3);
			IntSequence nvs(2); nvs[0] = 4; nvs[1] = 2;
			return index_forward<UGSTensor>(s, nvs);
		}
};

class IndexForwardFold : public TestRunnable {
public:
	IndexForwardFold()
		: TestRunnable("index forward for fold (55)(222)(22)", 7, 5) {}
	bool run() const
		{
			Symmetry s(2,3,2);
			IntSequence nvs(3); nvs[0] = 5; nvs[1] = 2; nvs[2] = 2;
			return index_forward<FGSTensor>(s, nvs);
		}
};

class IndexForwardUnfold : public TestRunnable {
public:
	IndexForwardUnfold()
		: TestRunnable("index forward for unfold (55)(222)(22)", 7, 5) {}
	bool run() const
		{
			Symmetry s(2,3,2);
			IntSequence nvs(3); nvs[0] = 5; nvs[1] = 2; nvs[2] = 2;
			return index_forward<UGSTensor>(s, nvs);
		}
};

class SmallIndexBackwardFold : public TestRunnable {
public:
	SmallIndexBackwardFold()
		: TestRunnable("small index backward for fold (3)(3)(222)", 5, 3) {}
	bool run() const
		{
			Symmetry s(1,1,3);
			IntSequence nvs(3); nvs[0] = 3; nvs[1] = 3; nvs[2] = 2;
			return index_backward<FGSTensor>(s, nvs);
		}
};

class IndexBackwardFold : public TestRunnable {
public:
	IndexBackwardFold()
		: TestRunnable("index backward for fold (44)(222)(44)", 7, 4) {}
	bool run() const
		{
			Symmetry s(2,3,2);
			IntSequence nvs(3); nvs[0] = 4; nvs[1] = 2; nvs[2] = 4;
			return index_backward<FGSTensor>(s, nvs);
		}
};

class SmallIndexBackwardUnfold : public TestRunnable {
public:
	SmallIndexBackwardUnfold()
		: TestRunnable("small index backward for unfold (3)(3)(222)", 5, 3) {}
	bool run() const
		{
			Symmetry s(1,1,3);
			IntSequence nvs(3); nvs[0] = 3; nvs[1] = 3; nvs[2] = 2;
			return index_backward<UGSTensor>(s, nvs);
		}
};

class IndexBackwardUnfold : public TestRunnable {
public:
	IndexBackwardUnfold()
		: TestRunnable("index backward for unfold (44)(222)(44)", 7, 4) {}
	bool run() const
		{
			Symmetry s(2,3,2);
			IntSequence nvs(3); nvs[0] = 4; nvs[1] = 2; nvs[2] = 4;
			return index_backward<UGSTensor>(s, nvs);
		}
};

class SmallIndexOffsetFold : public TestRunnable {
public:
	SmallIndexOffsetFold()
		: TestRunnable("small index offset for fold (44)(222)", 5, 4) {}
	bool run() const
		{
			Symmetry s(2,3);
			IntSequence nvs(2); nvs[0] = 4; nvs[1] = 2;
			return index_offset<FGSTensor>(s, nvs);
		}
};

class SmallIndexOffsetUnfold : public TestRunnable {
public:
	SmallIndexOffsetUnfold()
		: TestRunnable("small index offset for unfold (44)(222)", 5, 4) {}
	bool run() const
		{
			Symmetry s(2,3);
			IntSequence nvs(2); nvs[0] = 4; nvs[1] = 2;
			return index_offset<UGSTensor>(s, nvs);
		}
};

class IndexOffsetFold : public TestRunnable {
public:
	IndexOffsetFold()
		: TestRunnable("index offset for fold (55)(222)(22)", 5, 5) {}
	bool run() const
		{
			Symmetry s(2,3,2);
			IntSequence nvs(3); nvs[0] = 5; nvs[1] = 2; nvs[2] = 2;
			return index_offset<FGSTensor>(s, nvs);
		}
};

class IndexOffsetUnfold : public TestRunnable {
public:
	IndexOffsetUnfold()
		: TestRunnable("index offset for unfold (55)(222)(22)", 7, 5) {}
	bool run() const
		{
			Symmetry s(2,3,2);
			IntSequence nvs(3); nvs[0] = 5; nvs[1] = 2; nvs[2] = 2;
			return index_offset<UGSTensor>(s, nvs);
		}
};

class SmallFoldUnfoldFS : public TestRunnable {
public:
	SmallFoldUnfoldFS()
		: TestRunnable("small fold-unfold for full symmetry (444)", 3, 4) {}
	bool run() const
		{
			return fs_fold_unfold(5, 4, 3);
		}
};


class SmallFoldUnfoldGS : public TestRunnable {
public:
	SmallFoldUnfoldGS()
		: TestRunnable("small fold-unfold for gen symmetry (3)(33)(22)", 5, 3) {}
	bool run() const
		{
			Symmetry s(1,2,2);
			IntSequence nvs(3); nvs[0] = 3; nvs[1] = 3; nvs[2] = 2;
			return gs_fold_unfold(5, s, nvs);
		}
};

class FoldUnfoldFS : public TestRunnable {
public:
	FoldUnfoldFS()
		: TestRunnable("fold-unfold for full symmetry (9999)", 4, 9) {}
	bool run() const
		{
			return fs_fold_unfold(5, 9, 4);
		}
};


class FoldUnfoldGS : public TestRunnable {
public:
	FoldUnfoldGS()
		: TestRunnable("fold-unfold for gen symmetry (66)(2)(66)", 5, 6) {}
	bool run() const
		{
			Symmetry s(2,1,2);
			IntSequence nvs(3); nvs[0] = 6; nvs[1] = 2; nvs[2] = 6;
			return gs_fold_unfold(5, s, nvs);
		}
};

class SmallFoldUnfoldR : public TestRunnable {
public:
	SmallFoldUnfoldR()
		: TestRunnable("small fold-unfold for row full symmetry (333)", 3, 3) {}
	bool run() const
		{
			return r_fold_unfold(5, 3, 3);
		}
};

class FoldUnfoldR : public TestRunnable {
public:
	FoldUnfoldR()
		: TestRunnable("fold-unfold for row full symmetry (66666)", 5, 6) {}
	bool run() const
		{
			return r_fold_unfold(5, 6, 5);
		}
};

class SmallDenseProd : public TestRunnable {
public:
	SmallDenseProd()
		: TestRunnable("small dense prod bsym=1-2,nvs=3-2,h=2-3,r=2",3,3) {}
	bool run() const
		{
			IntSequence bnvs(2); bnvs[0]=3; bnvs[1]=2;
			return dense_prod(Symmetry(1,2), bnvs, 2, 3, 2);
		}
};

class DenseProd : public TestRunnable {
public:
	DenseProd()
		: TestRunnable("dense prod bsym=2-3,nvs=10-7,h=3-15,r=10",5,15) {}
	bool run() const
		{
			IntSequence bnvs(2); bnvs[0]=10; bnvs[1]=7;
			return dense_prod(Symmetry(2,3), bnvs, 3, 15, 10);
		}
};

class BigDenseProd : public TestRunnable {
public:
	BigDenseProd()
		: TestRunnable("dense prod bsym=3-2,nvs=13-11,h=3-20,r=20",6,20) {}
	bool run() const
		{
			IntSequence bnvs(2); bnvs[0]=13; bnvs[1]=11;
			return dense_prod(Symmetry(3,2), bnvs, 3, 20, 20);
		}
};

class SmallFoldedMonomial : public TestRunnable {
public:
	SmallFoldedMonomial()
		: TestRunnable("folded vrs. monoms (g,x,y,u)=(10,4,5,3), dim=4", 4, 8) {}
	bool run() const
		{
			return folded_monomial(10, 4, 5, 3, 4);
		}
};

class FoldedMonomial : public TestRunnable {
public:
	FoldedMonomial()
		: TestRunnable("folded vrs. monoms (g,x,y,u)=(20,12,10,5), dim=4", 4, 15) {}
	bool run() const
		{
			return folded_monomial(20, 12, 10, 5, 4);
		}
};

class SmallUnfoldedMonomial : public TestRunnable {
public:
	SmallUnfoldedMonomial()
		: TestRunnable("unfolded vrs. monoms (g,x,y,u)=(10,4,5,3), dim=4", 4, 8) {}
	bool run() const
		{
			return unfolded_monomial(10, 4, 5, 3, 4);
		}
};

class UnfoldedMonomial : public TestRunnable {
public:
	UnfoldedMonomial()
		: TestRunnable("unfolded vrs. monoms (g,x,y,u)=(20,12,10,5), dim=4", 4, 15) {}
	bool run() const
		{
			return unfolded_monomial(20, 12, 10, 5, 4);
		}
};

class FoldedContractionSmall : public TestRunnable {
public:
	FoldedContractionSmall()
		: TestRunnable("folded contraction small (r=5, nv=4, dim=3)", 3, 4) {}
	bool run() const
		{
			return folded_contraction(5, 4, 3);
		}
};

class FoldedContractionBig : public TestRunnable {
public:
	FoldedContractionBig()
		: TestRunnable("folded contraction big (r=20, nv=12, dim=5)", 5, 12) {}
	bool run() const
		{
			return folded_contraction(20, 12, 5);
		}
};

class UnfoldedContractionSmall : public TestRunnable {
public:
	UnfoldedContractionSmall()
		: TestRunnable("unfolded contraction small (r=5, nv=4, dim=3)", 3, 4) {}
	bool run() const
		{
			return unfolded_contraction(5, 4, 3);
		}
};

class UnfoldedContractionBig : public TestRunnable {
public:
	UnfoldedContractionBig()
		: TestRunnable("unfolded contraction big (r=20, nv=12, dim=5)", 5, 12) {}
	bool run() const
		{
			return unfolded_contraction(20, 12, 5);
		}
};

class PolyEvalSmall : public TestRunnable {
public:
	PolyEvalSmall()
		: TestRunnable("polynomial evaluation small (r=4, nv=5, maxdim=4)", 4, 5) {}
	bool run() const
		{
			return poly_eval(4, 5, 4);
		}
};

class PolyEvalBig : public TestRunnable {
public:
	PolyEvalBig()
		: TestRunnable("polynomial evaluation big (r=244, nv=97, maxdim=2)", 2, 97) {}
	bool run() const
		{
			return poly_eval(244, 97, 2);
		}
};

class FoldZContSmall : public TestRunnable {
public:
	FoldZContSmall()
		: TestRunnable("folded Z container (r=3,ny=2,nu=2,nup=1,G=2,g=2,dim=3)",
					   3, 8) {}
	bool run() const
		{
			return fold_zcont(3, 2, 2, 1, 2, 2, 3);
		}
};

class FoldZCont : public TestRunnable {
public:
	FoldZCont()
		: TestRunnable("folded Z container (r=13,ny=5,nu=7,nup=4,G=6,g=7,dim=4)",
					   4, 25) {}
	bool run() const
		{
			return fold_zcont(13, 5, 7, 4, 6, 7, 4);
		}
};

class UnfoldZContSmall : public TestRunnable {
public:
	UnfoldZContSmall()
		: TestRunnable("unfolded Z container (r=3,ny=2,nu=2,nup=1,G=2,g=2,dim=3)",
					   3, 8) {}
	bool run() const
		{
			return unfold_zcont(3, 2, 2, 1, 2, 2, 3);
		}
};

class UnfoldZCont : public TestRunnable {
public:
	UnfoldZCont()
		: TestRunnable("unfolded Z container (r=13,ny=5,nu=7,nup=4,G=6,g=7,dim=4",
					   4, 25) {}
	bool run() const
		{
			return unfold_zcont(13, 5, 7, 4, 6, 7, 4);
		}
};



int main()
{
	TestRunnable* all_tests[50];
	// fill in vector of all tests
	int num_tests = 0;
	all_tests[num_tests++] = new SmallIndexForwardFold();
	all_tests[num_tests++] = new SmallIndexForwardUnfold();
	all_tests[num_tests++] = new IndexForwardFold();
	all_tests[num_tests++] = new IndexForwardUnfold();
	all_tests[num_tests++] = new SmallIndexBackwardFold();
	all_tests[num_tests++] = new IndexBackwardFold();
	all_tests[num_tests++] = new SmallIndexBackwardUnfold();
	all_tests[num_tests++] = new IndexBackwardUnfold();
	all_tests[num_tests++] = new SmallIndexOffsetFold();
	all_tests[num_tests++] = new SmallIndexOffsetUnfold();
	all_tests[num_tests++] = new IndexOffsetFold();
	all_tests[num_tests++] = new IndexOffsetUnfold();
	all_tests[num_tests++] = new SmallFoldUnfoldFS();
	all_tests[num_tests++] = new SmallFoldUnfoldGS();
	all_tests[num_tests++] = new FoldUnfoldFS();
	all_tests[num_tests++] = new FoldUnfoldGS();
	all_tests[num_tests++] = new SmallFoldUnfoldR();
	all_tests[num_tests++] = new FoldUnfoldR();
	all_tests[num_tests++] = new SmallDenseProd();
	all_tests[num_tests++] = new DenseProd();
	all_tests[num_tests++] = new BigDenseProd();
	all_tests[num_tests++] = new SmallFoldedMonomial();
	all_tests[num_tests++] = new FoldedMonomial();
	all_tests[num_tests++] = new SmallUnfoldedMonomial();
	all_tests[num_tests++] = new UnfoldedMonomial();
	all_tests[num_tests++] = new FoldedContractionSmall();
	all_tests[num_tests++] = new FoldedContractionBig();
	all_tests[num_tests++] = new UnfoldedContractionSmall();
	all_tests[num_tests++] = new UnfoldedContractionBig();
	all_tests[num_tests++] = new PolyEvalSmall();
	all_tests[num_tests++] = new PolyEvalBig();
	all_tests[num_tests++] = new FoldZContSmall();
	all_tests[num_tests++] = new FoldZCont();
	all_tests[num_tests++] = new UnfoldZContSmall();
	all_tests[num_tests++] = new UnfoldZCont();

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
