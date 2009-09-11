/* $Id: factory.cpp 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#include "factory.h"

#include <cstdlib>
#include <math.h>

void Factory::init(const Symmetry& s, const IntSequence& nvs)
{
	IntSequence sym(s);
	long int seed = sym[0];
	seed = 256*seed + nvs[0];
	if (sym.size() > 1)
		seed = 256*seed + sym[1];
	if (nvs.size() > 1)
		seed = 256*seed + nvs[0];
	srand48(seed);
}

void Factory::init(int dim, int nv)
{
	long int seed = dim;
	seed = 256*seed + nv;
	srand48(seed);
}

double Factory::get() const
{
	return 1.0*(drand48()-0.5);
}

void Factory::fillMatrix(TwoDMatrix& m) const
{
	Vector& d = m.getData();
	for (int i = 0; i < d.length(); i++)
		d[i] = get();
}

Vector* Factory::makeVector(int n)
{
	init(n, n*n);

	Vector* v = new Vector(n);
	for (int i = 0; i < n; i++)
		(*v)[i] = get();

	return v;
}
