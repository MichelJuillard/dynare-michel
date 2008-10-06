/* $Id: monoms.h 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#ifndef MONOMS_H
#define MONOMS_H

#include "int_sequence.h"
#include "gs_tensor.h"
#include "t_container.h"
#include "sparse_tensor.h"
#include "Vector.h"

class IntGenerator {
	int maxim;
	double probab;
public:
	IntGenerator()
		: maxim(5), probab(0.3) {}
	void init(int nf, int ny, int nv, int nw, int nu, int mx, double prob);
	int get() const;
};

extern IntGenerator intgen;


class Monom : public IntSequence {
public:
	Monom(int len); // generate a random monom
	Monom(int len, int item); // generate monom whose items are the given item
	double deriv(const IntSequence& vars) const;
	// this = this*m^ex (in monomial sense)
	void multiplyWith(int ex, const Monom& m);
	void print() const;
};

class Monom2Vector;
class Monom1Vector {
	friend class Monom2Vector;
	int nx;
	int len;
	Monom** const x;
public:
	Monom1Vector(int nxx, int l);
	~Monom1Vector();
	void deriv(const IntSequence& c, Vector& out) const;
	FGSTensor* deriv(int dim) const;
	void print() const;
};

//class Monom3Vector;
class Monom2Vector {
	int ny;
	int nu;
	int len;
	Monom** const y;
	Monom** const u;
public:
	// generate random vector of monom two vector
	Monom2Vector(int nyy, int nuu, int l);
	// calculate g(x(y,u))
	Monom2Vector(const Monom1Vector& g, const Monom2Vector& xmon);
	~Monom2Vector();
	void deriv(const Symmetry& s, const IntSequence& c, Vector& out) const;
	FGSTensor* deriv(const Symmetry& s) const;
	FGSContainer* deriv(int maxdim) const;
	void print() const;
};

class Monom4Vector {
	int len;
	int nx1;
	int nx2;
	int nx3;
	int nx4;
	Monom** const x1;
	Monom** const x2;
	Monom** const x3;
	Monom** const x4;
public:
    /* random for g(y,u,sigma) */
	Monom4Vector(int l, int ny, int nu);
	/* random for G(y,u,u',sigma) */
	Monom4Vector(int l, int ny, int nu, int nup);
	/* random for f(y+,y,y-,u) */
	Monom4Vector(int l, int nbigg, int ng, int ny, int nu);
	/* substitution f(G(y,u,u',sigma),g(y,u,sigma),y,u) */
	Monom4Vector(const Monom4Vector& f, const Monom4Vector& bigg,
				 const Monom4Vector& g);
	~Monom4Vector();
	FSSparseTensor* deriv(int dim) const;
	FGSTensor* deriv(const Symmetry& s) const;
	void deriv(const Symmetry& s, const IntSequence& coor, Vector& out) const;
	void print() const;
protected:
	void init_random();
};


struct SparseDerivGenerator {
	int maxdimen;
	FGSContainer* bigg;
	FGSContainer* g;
	FGSContainer* rcont;
	FSSparseTensor** const ts;
	SparseDerivGenerator(int nf, int ny, int nu, int nup, int nbigg, int ng,
						 int mx, double prob, int maxdim);
	~SparseDerivGenerator();
};


struct DenseDerivGenerator {
	int maxdimen;
	FGSContainer* xcont;
	FGSContainer* rcont;
	FGSTensor** const ts;
	UGSContainer* uxcont;
	UGSTensor** const uts;
	DenseDerivGenerator(int ng, int nx, int ny, int nu,
						int mx, double prob, int maxdim);
	void unfold();
	~DenseDerivGenerator();
};

#endif

// Local Variables:
// mode:C++
// End:
