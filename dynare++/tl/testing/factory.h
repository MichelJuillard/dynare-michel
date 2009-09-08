/* $Id: factory.h 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#ifndef FACTORY_H
#define FACTORY_H 

#include "symmetry.h"
#include "int_sequence.h"
#include "twod_matrix.h"
#include "equivalence.h"
#include "rfs_tensor.h"
#include "t_container.h"

class Factory {
	void init(const Symmetry& s, const IntSequence& nvs);
	void init(int dim, int nv);
	void fillMatrix(TwoDMatrix& m) const;
public:
	double get() const;
	// this can be used with UGSTensor, FGSTensor
	template <class _Ttype>
	_Ttype* make(int r, const Symmetry& s, const IntSequence& nvs)
		{
			_Ttype* res = new _Ttype(r, TensorDimens(s, nvs));
			init(s, nvs);
			fillMatrix(*res);
			return res;
		}

	// this can be used with FFSTensor, UFSTensor, FRTensor, URTensor
	template <class _Ttype>
	_Ttype* make(int r, int nv, int dim)
		{
			_Ttype* res = new _Ttype(r, nv, dim);
			init(dim, nv);
			fillMatrix(*res);
			return res;
		}

	template <class _Ttype, class _Ctype>
	_Ctype* makeCont(int r, const IntSequence& nvs, int maxdim)
		{
			int symnum = nvs.size();
			_Ctype* res = new _Ctype(symnum);
			for (int dim = 1; dim <= maxdim; dim++) {
				if (symnum == 1) {
					// full symmetry
					Symmetry sym(dim);
					_Ttype* t = make<_Ttype>(r, sym, nvs);
					res->insert(t);
				} else {
					// general symmetry
					for (int i = 0; i <= dim; i++) {
						Symmetry sym(i, dim-i);
						_Ttype* t = make<_Ttype>(r, sym, nvs);
						res->insert(t);
					}
				}
			}
			return res;
		}

	template <class _Ttype, class _Ptype>
	_Ptype* makePoly(int r, int nv, int maxdim)
		{
			_Ptype* p = new _Ptype(r, nv);
			for (int d = 1; d <= maxdim; d++) {
				_Ttype* t = make<_Ttype>(r, nv, d);
				p->insert(t);
			}
			return p;
		}

	Vector* makeVector(int n);
};

#endif

// Local Variables:
// mode:C++
// End:
