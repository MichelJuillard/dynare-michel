/* $Id: monoms.cpp 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#include "monoms.h"
#include "tl_exception.h"
#include "fs_tensor.h"

#include <math.h>
#include <stdio.h>

IntGenerator intgen;

void IntGenerator::init(int nf, int ny, int nv, int nw, int nu,
						int mx, double prob)
{
	maxim = mx;
	probab = prob;
	long int seed = nf;
	seed = 256*seed + ny;
	seed = 256*seed + nv;
	seed = 256*seed + nw;
	seed = 256*seed + nu;
	srand48(seed);
}

int IntGenerator::get() const
{
	double d = drand48();
	int num_inter = (int)( ((double)2*maxim)/(1.0-probab));
	int num_zero_inter = num_inter - 2*maxim;
	if (d < ((double)num_zero_inter)/num_inter)
		return 0;
	return (int)(d*num_inter)-num_zero_inter-maxim;
}

Monom::Monom(int len)
	: IntSequence(len)
{
	for (int i = 0; i < len; i++)
		operator[](i) = intgen.get();
}

Monom::Monom(int len, int item)
	: IntSequence(len, item)
{
}

double Monom::deriv(const IntSequence& vars) const
{
	double res = 1.0;
	int first_same_i = 0;
	for (int i = 0; i < vars.size(); i++) {
		TL_RAISE_IF(vars[i] < 0 || vars[i] >= size(),
					"Wrong variable index in Monom::deriv");
		if (vars[i] != vars[first_same_i])
			first_same_i = i;
		int mult = operator[](vars[i]) - (i-first_same_i);
		if (mult == 0)
			return 0;
		res *= mult;
	}
	return res;
}

void Monom::multiplyWith(int ex, const Monom& m)
{
	TL_RAISE_IF(size() != m.size(),
				"Wrong sizes of monoms in Monom::multiplyWith");
	if (ex == 0)
		return;
	for (int i = 0; i < size(); i++)
		operator[](i) += m[i]*ex;
}

void Monom::print() const
{
	printf("[");
	for (int i = 0; i < size(); i++)
		printf("%3d", operator[](i));
	printf("]");
}

Monom1Vector::Monom1Vector(int nxx, int l)
	: nx(nxx), len(l), x(new (Monom*)[len])
{
	for (int i = 0; i < len; i++) {
		x[i] = new Monom(nx);
	}
}

Monom1Vector::~Monom1Vector()
{
	for (int i = 0; i < len; i++) {
		delete x[i];
	}
	delete [] x;
}

void Monom1Vector::deriv(const IntSequence& c, Vector& out) const
{
	TL_RAISE_IF(out.length() != len,
				"Wrong length of output vector in Monom1Vector::deriv");

	for (int i = 0; i < len; i++) {
		out[i] = x[i]->deriv(c);
	}
}


FGSTensor* Monom1Vector::deriv(int dim) const
{
	FGSTensor* res =
		new FGSTensor(len, TensorDimens(Symmetry(dim), IntSequence(1, nx)));
	for (Tensor::index it = res->begin(); it != res->end(); ++it) {
		Vector outcol(*res, *it);
		deriv(it.getCoor(), outcol);
	}
	return res;
}

void Monom1Vector::print() const
{
	printf("Variables: x(%d)\n", nx);
	printf("Rows: %d\n", len);
	for (int i = 0; i < len; i++) {
		printf("%2d: ", i);
		x[i]->print();
		printf("\n");
	}
}

Monom2Vector::Monom2Vector(int nyy, int nuu, int l)
	: ny(nyy), nu(nuu), len(l), y(new (Monom*)[len]), u(new (Monom*)[len])
{
	for (int i = 0; i < len; i++) {
		y[i] = new Monom(ny);
		u[i] = new Monom(nu);
	}
}

Monom2Vector::Monom2Vector(const Monom1Vector& g, const Monom2Vector& xmon)
	: ny(xmon.ny), nu(xmon.nu), len(g.len),
	  y(new (Monom*)[len]), u(new (Monom*)[len])	
{
	TL_RAISE_IF(xmon.len != g.nx,
				"Wrong number of x's in Monom2Vector constructor");

	for (int i = 0; i < len; i++) {
		y[i] = new Monom(ny, 0);
		u[i] = new Monom(nu, 0);
	}

	for (int i = 0; i < len; i++) {
		// multiply from xmon
		for (int j = 0; j < g.nx; j++) {
			int ex = g.x[i]->operator[](j);
			y[i]->multiplyWith(ex, *(xmon.y[j]));
			u[i]->multiplyWith(ex, *(xmon.u[j]));
		}
	}
}

Monom2Vector::~Monom2Vector()
{
	for (int i = 0; i < len; i++) {
		delete y[i];
		delete u[i];
	}
	delete [] y;
	delete [] u;
}

void Monom2Vector::deriv(const Symmetry& s, const IntSequence& c,
						 Vector& out) const
{
	TL_RAISE_IF(out.length() != len,
				"Wrong length of output vector in Monom2Vector::deriv");
	TL_RAISE_IF(s.num() != 2,
				"Wrong symmetry for Monom2Vector::deriv");
	TL_RAISE_IF(s.dimen() != c.size(),
				"Incompatible symmetry and coordinates in Monom2Vector::deriv");
	IntSequence cy(c, 0, s[0]);
	IntSequence cu(c, s[0], s.dimen());
	for (int i = 0; i < len; i++) {
		out[i] = y[i]->deriv(cy) * u[i]->deriv(cu);
	}
}

FGSTensor* Monom2Vector::deriv(const Symmetry& s) const
{
	IntSequence nvs(2); nvs[0] = ny; nvs[1] = nu;
	FGSTensor* t = new FGSTensor(len, TensorDimens(s, nvs));
	for (Tensor::index it = t->begin(); it != t->end(); ++it) {
		Vector col(*t, *it);
		deriv(s, it.getCoor(), col);
	}
	return t;
}

FGSContainer* Monom2Vector::deriv(int maxdim) const
{
	FGSContainer* res = new FGSContainer(2);
	for (int dim = 1; dim <= maxdim; dim++) {
		for (int ydim = 0; ydim <= dim; ydim++) {
			int udim = dim - ydim;
			Symmetry s(ydim, udim);
			res->insert(deriv(s));
		}
	}
	return res;
}

void Monom2Vector::print() const
{
	printf("Variables: y(%d) u(%d)\n", ny, nu);
	printf("Rows: %d\n", len);
	for (int i = 0; i < len; i++) {
		printf("%2d: ", i);
		y[i]->print();
		printf("    ");
		u[i]->print();
		printf("\n");
	}
}

Monom4Vector::~Monom4Vector()
{
	for (int i = 0; i < len; i++) {
		delete x1[i];
		delete x2[i];
		delete x3[i];
		delete x4[i];
	}
	delete [] x1;
	delete [] x2;
	delete [] x3;
	delete [] x4;
}

void Monom4Vector::init_random()
{
	for (int i = 0; i < len; i++) {
		x1[i] = new Monom(nx1);
		x2[i] = new Monom(nx2);
		x3[i] = new Monom(nx3);
		x4[i] = new Monom(nx4);
	}
}

Monom4Vector::Monom4Vector(int l, int ny, int nu)
	: len(l), nx1(ny), nx2(nu), nx3(0), nx4(1),
	  x1(new (Monom*)[len]),
	  x2(new (Monom*)[len]),
	  x3(new (Monom*)[len]),
	  x4(new (Monom*)[len])	  
{
	init_random();
}

Monom4Vector::Monom4Vector(int l, int ny, int nu, int nup)
	: len(l), nx1(ny), nx2(nu), nx3(nup), nx4(1),
	  x1(new (Monom*)[len]),
	  x2(new (Monom*)[len]),
	  x3(new (Monom*)[len]),
	  x4(new (Monom*)[len])
{
	init_random();
}

Monom4Vector::Monom4Vector(int l, int nbigg, int ng, int ny, int nu)
	: len(l), nx1(nbigg), nx2(ng), nx3(ny), nx4(nu),
	  x1(new (Monom*)[len]),
	  x2(new (Monom*)[len]),
	  x3(new (Monom*)[len]),
	  x4(new (Monom*)[len])
{
	init_random();
}

Monom4Vector::Monom4Vector(const Monom4Vector& f, const Monom4Vector& bigg,
						   const Monom4Vector& g)
	: len(f.len), nx1(bigg.nx1), nx2(bigg.nx2), nx3(bigg.nx3), nx4(1),
	  x1(new (Monom*)[len]),
	  x2(new (Monom*)[len]),
	  x3(new (Monom*)[len]),
	  x4(new (Monom*)[len])
{
	TL_RAISE_IF(!(bigg.nx1 == g.nx1 && bigg.nx2 == g.nx2 && g.nx3 == 0 &&
				  bigg.nx4 == 1 && g.nx4 == 1),
				"Incompatible g with G");
	TL_RAISE_IF(!(bigg.len == f.nx1 && g.len == f.nx2 &&
				  bigg.nx1 == f.nx3 && bigg.nx2 == f.nx4),
				"Incompatible g or G with f");

	for (int i = 0; i < len; i++) {
		x1[i] = new Monom(nx1, 0);
		x2[i] = new Monom(nx2, 0);
		x3[i] = new Monom(nx3, 0);
		x4[i] = new Monom(nx4, 0);
	}

	for (int i = 0; i < len; i++) {
		// multiply from G (first argument)
		for (int j = 0; j < f.nx1; j++) {
			int ex = f.x1[i]->operator[](j);
			x1[i]->multiplyWith(ex, *(bigg.x1[j]));
			x2[i]->multiplyWith(ex, *(bigg.x2[j]));
			x3[i]->multiplyWith(ex, *(bigg.x3[j]));
			x4[i]->multiplyWith(ex, *(bigg.x4[j]));
		}
		// multiply from g (second argument)
		for (int j = 0; j < f.nx2; j++) {
			int ex = f.x2[i]->operator[](j);
			x1[i]->multiplyWith(ex, *(g.x1[j]));
			x2[i]->multiplyWith(ex, *(g.x2[j]));
			x4[i]->multiplyWith(ex, *(g.x4[j]));
		}
		// add y as third argument of f
		x1[i]->add(1, *(f.x3[i]));
		// add u as fourth argument of f
		x2[i]->add(1, *(f.x4[i]));
	}
}

void Monom4Vector::deriv(const Symmetry& s, const IntSequence& coor,
						 Vector& out) const
{
	TL_RAISE_IF(out.length() != len,
				"Wrong length of output vector in Monom4Vector::deriv");
	TL_RAISE_IF(s.num() != 4,
				"Wrong symmetry for Monom4Vector::deriv");
	TL_RAISE_IF(s.dimen() != coor.size(),
				"Incompatible symmetry and coordinates in Monom4Vector::deriv");

	for (int i = 0; i < len; i++) {
		out[i] = 1;
		int off = 0;
		out[i] *= x1[i]->deriv(IntSequence(coor, off, off+s[0]));
		off += s[0];
		out[i] *= x2[i]->deriv(IntSequence(coor, off, off+s[1]));
		off += s[1];
		out[i] *= x3[i]->deriv(IntSequence(coor, off, off+s[2]));
		off += s[2];
		out[i] *= x4[i]->deriv(IntSequence(coor, off, off+s[3]));
	}
}

FGSTensor* Monom4Vector::deriv(const Symmetry& s) const
{
	IntSequence nvs(4);
	nvs[0] = nx1; nvs[1] = nx2; 
	nvs[2] = nx3; nvs[3] = nx4;

	FGSTensor* res = new FGSTensor(len, TensorDimens(s, nvs));
	for (Tensor::index run = res->begin(); run != res->end(); ++run) {
		Vector col(*res, *run);
		deriv(s, run.getCoor(), col);
	}
	return res;
}

FSSparseTensor* Monom4Vector::deriv(int dim) const
{
	IntSequence cum(4);
	cum[0] = 0; cum[1] = nx1; cum[2] = nx1+nx2; cum[3] = nx1+nx2+nx3;

	FSSparseTensor* res = new FSSparseTensor(dim, nx1+nx2+nx3+nx4, len);

	FFSTensor dummy(0, nx1+nx2+nx3+nx4, dim);
	for (Tensor::index run = dummy.begin(); run != dummy.end(); ++run) {
		Symmetry ind_sym(0,0,0,0);
		IntSequence ind(run.getCoor());
		for (int i = 0; i < ind.size(); i++) {
			int j = 3;
			while (j >= 0 && ind[i] < cum[j]) j--;
			ind_sym[j]++;
			ind[i] -= cum[j];
		}

		Vector col(len);
		deriv(ind_sym, ind, col);
		for (int i = 0; i < len; i++) {
			if (col[i] != 0.0) {
				res->insert(run.getCoor(), i, col[i]);
			}
		}
	}

	return res;
}

void Monom4Vector::print() const
{
	printf("Variables: x1(%d) x2(%d) x3(%d) x4(%d)\n",
		   nx1, nx2, nx3, nx4);
	printf("Rows: %d\n", len);
	for (int i = 0; i < len; i++) {
		printf("%2d: ", i);
		x1[i]->print();
		printf("    ");
		x2[i]->print();
		printf("    ");
		x3[i]->print();
		printf("    ");
		x4[i]->print();
		printf("\n");
	}
}

SparseDerivGenerator::SparseDerivGenerator(
	int nf, int ny, int nu, int nup, int nbigg, int ng,
	int mx, double prob, int maxdim)
	: maxdimen(maxdim), ts(new (FSSparseTensor*)[maxdimen])
{
	intgen.init(nf, ny, nu, nup, nbigg, mx, prob);

	Monom4Vector bigg_m(nbigg, ny, nu, nup);
	Monom4Vector g_m(ng, ny, nu);
	Monom4Vector f(nf, nbigg, ng, ny, nu);
	Monom4Vector r(f, bigg_m, g_m);
	bigg = new FGSContainer(4);
	g = new FGSContainer(4);
	rcont = new FGSContainer(4);

	for (int dim = 1; dim <= maxdimen; dim++) {
		SymmetrySet ss(dim, 4);
		for (symiterator si(ss); !si.isEnd(); ++si) {
			bigg->insert(bigg_m.deriv(*si));
			rcont->insert(r.deriv(*si));
			if ((*si)[2] == 0)
				g->insert(g_m.deriv(*si));
		}

		ts[dim-1] = f.deriv(dim);
	}
}

SparseDerivGenerator::~SparseDerivGenerator()
{
	delete bigg;
	delete g;
	delete rcont;
	for(int i = 0; i < maxdimen; i++)
		delete ts[i];
	delete [] ts;
}


DenseDerivGenerator::DenseDerivGenerator(int ng, int nx, int ny, int nu,
										 int mx, double prob, int maxdim)
	: maxdimen(maxdim), ts(new (FGSTensor*)[maxdimen]),
	  uts(new (UGSTensor*)[maxdimen])
{
	intgen.init(ng, nx, ny, nu, nu, mx, prob);
	Monom1Vector g(nx, ng);
	Monom2Vector x(ny, nu, nx);
	Monom2Vector r(g, x);
	xcont = x.deriv(maxdimen);
	rcont = r.deriv(maxdimen);
	uxcont = NULL;
	for (int d = 1; d <= maxdimen; d++) {
		ts[d-1] = g.deriv(d);
		uts[d-1] = NULL;
	}
}

void DenseDerivGenerator::unfold()
{
	uxcont = new UGSContainer(*xcont);
	for (int i = 0; i < maxdimen; i++) {
		uts[i] = new UGSTensor(*(ts[i]));
	}
}

DenseDerivGenerator::~DenseDerivGenerator()
{
	delete xcont;
	delete rcont;
	for (int i = 0; i < maxdimen; i++) {
		delete ts[i];
		if (uts[i])
			delete uts[i];
	}
	delete [] ts;	
	delete [] uts;	
}

// Local Variables:
// mode:C++
// End:
