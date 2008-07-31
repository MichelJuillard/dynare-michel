/*
 * Copyright (C) 2003-2005 Ondra Kamenik
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/QuasiTriangular.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef QUASI_TRIANGULAR_H
#define QUASI_TRIANGULAR_H

#include "Vector.h"
#include "KronVector.h"
#include "SylvMatrix.h"

#include <list>

using namespace std;

class DiagonalBlock;
class Diagonal;
class DiagPair {
private:
	double* a1;
	double* a2;
public:
	DiagPair() {}
	DiagPair(double* aa1, double* aa2) {a1 = aa1; a2 = aa2;}
	DiagPair(const DiagPair& p) {a1 = p.a1; a2 = p.a2;}
	const DiagPair& operator=(const DiagPair& p) {a1 = p.a1; a2 = p.a2; return *this;}
	const DiagPair& operator=(double v) {*a1 = v; *a2 = v; return *this;}
	const double& operator*() const {return *a1;}
	/** here we must not define double& operator*(), since it wouldn't 
	 rewrite both values, we use operator= for this */ 
	friend class Diagonal;
	friend class DiagonalBlock;
};

class DiagonalBlock {
private:
	int jbar;
	bool real;
	DiagPair alpha;
	double* beta1;
	double* beta2;

	void copy(const DiagonalBlock& b) {
		jbar = b.jbar;
		real = b.real;
		alpha = b.alpha;
		beta1 = b.beta1;
		beta2 = b.beta2;
	}

public:
	DiagonalBlock() {}
	DiagonalBlock(int jb, bool r, double* a1, double* a2,
				  double* b1, double* b2)
		: alpha(a1, a2)
		{
			jbar = jb;
			real = r;
			beta1 = b1;
			beta2 = b2;
		}
	// construct complex block
	DiagonalBlock(int jb, double* a1, double* a2)
		: alpha(a1, a2)
		{
			jbar = jb;
			real = false;
			beta1 = a2 - 1;
			beta2 = a1 + 1;
		}
	// construct real block
	DiagonalBlock(int jb, double* a1)
		: alpha(a1, a1)
		{
			jbar = jb;
			real = true;
			beta1 = 0;
			beta2 = 0;
		}
	DiagonalBlock(const DiagonalBlock& b)
		{copy(b);}
	const DiagonalBlock& operator=(const DiagonalBlock& b)
		{copy(b); return *this;}
	int getIndex() const
		{return jbar;}
	bool isReal() const
		{return real;}
	const DiagPair& getAlpha() const
		{return alpha;}
	DiagPair& getAlpha()
		{return alpha;}
	double& getBeta1() const
		{return *beta1;}
	double& getBeta2() const
		{return *beta2;}
	double getDeterminant() const;
	double getSBeta() const;
	double getSize() const;
	void setReal();
	// for debugging
	void checkBlock(const double* d, int d_size);
	friend class Diagonal;
};

template <class _Tdiag, class _Tblock, class _Titer>
struct _diag_iter {
	typedef _diag_iter<_Tdiag, _Tblock, _Titer> _Self;
	_Tdiag diag;
	_Titer it;
public:
	_diag_iter(_Tdiag d, _Titer iter) : diag(d), it(iter) {}
	_Tblock operator*() const {return *it;}
	_Self& operator++() {++it; return *this;}
	_Self& operator--() {--it; return *this;}
	bool operator==(const _Self& x) const {return x.it == it;}
	bool operator!=(const _Self& x) const {return x.it != it;}
	const _Self& operator=(const _Self& x) {it = x.it; return *this;}
	_Titer iter() const {return it;}
};

class Diagonal {
public:
	typedef _diag_iter<const Diagonal&, const DiagonalBlock&, list<DiagonalBlock>::const_iterator> const_diag_iter;
	typedef _diag_iter<Diagonal&, DiagonalBlock&, list<DiagonalBlock>::iterator> diag_iter;
private:
	int num_all;
	list<DiagonalBlock> blocks;
	int num_real;
	void copy(const Diagonal&);
public:
	Diagonal() : num_all(0), num_real(0) {}
	Diagonal(double* data, int d_size);
	Diagonal(double* data, const Diagonal& d);
	Diagonal(const Diagonal& d) {copy(d);}
	const Diagonal& operator =(const Diagonal& d) {copy(d); return *this;}
	virtual ~Diagonal() {}

	int getNumComplex() const {return num_all - num_real;}
	int getNumReal() const {return num_real;}
	int getSize() const {return getNumReal() + 2*getNumComplex();}
	int getNumBlocks() const {return num_all;}
	void getEigenValues(Vector& eig) const;
	void swapLogically(diag_iter it);
	void checkConsistency(diag_iter it);
	double getAverageSize(diag_iter start, diag_iter end);
	diag_iter findClosestBlock(diag_iter start, diag_iter end, double a);
	diag_iter findNextLargerBlock(diag_iter start, diag_iter end, double a);
	void print() const;

	diag_iter begin()
		{return diag_iter(*this, blocks.begin());}
	const_diag_iter begin() const
		{return const_diag_iter(*this, blocks.begin());}
	diag_iter end()
		{return diag_iter(*this, blocks.end());}
	const_diag_iter end() const
		{return const_diag_iter(*this, blocks.end());}

	/* redefine pointers as data start at p */
	void changeBase(double* p);
private:
	static double EPS;
	static int getNumComplex(const double* data, int d_size);
	static bool isZero(double p);
};

template <class _TRef, class _TPtr>
struct _matrix_iter {
	typedef _matrix_iter<_TRef, _TPtr> _Self;
	int d_size;
	bool real;
	_TPtr ptr;
public:
	_matrix_iter(_TPtr base, int ds, bool r)
		{ptr = base; d_size = ds; real = r;}
	virtual ~_matrix_iter() {}
	const _Self& operator=(const _Self& it)
		{ptr = it.ptr; d_size = it.d_size; real = it.real; return *this;}
	bool operator==(const _Self& it) const
		{return ptr == it.ptr;}
	bool operator!=(const _Self& it) const
		{return ptr != it.ptr;}
	_TRef operator*() const
		{return *ptr;}
	_TRef a() const
		{return *ptr;}
	virtual _Self& operator++() =0;
};

template <class _TRef, class _TPtr>
class _column_iter : public _matrix_iter<_TRef, _TPtr> {
	typedef _matrix_iter<_TRef, _TPtr> _Tparent; 
	typedef _column_iter<_TRef, _TPtr> _Self;
	int row;
public:
	_column_iter(_TPtr base, int ds, bool r, int rw)
		: _matrix_iter<_TRef, _TPtr>(base, ds, r), row(rw) {};
	_Self& operator++()
		{_Tparent::ptr++; row++; return *this;}
	_TRef b() const
		{
			if (_Tparent::real) {
				return *(_Tparent::ptr);
			} else {
				return *(_Tparent::ptr+_Tparent::d_size);
			}
		}
	int getRow() const {return row;}
};

template <class _TRef, class _TPtr>
class _row_iter : public _matrix_iter<_TRef, _TPtr> {
	typedef _matrix_iter<_TRef, _TPtr> _Tparent; 
	typedef _row_iter<_TRef, _TPtr> _Self;
	int col;
public:
	_row_iter(_TPtr base, int ds, bool r, int cl)
		: _matrix_iter<_TRef, _TPtr>(base, ds, r), col(cl) {};
	_Self& operator++()
		{_Tparent::ptr += _Tparent::d_size; col++; return *this;}
	virtual _TRef b() const
		{
			if (_Tparent::real) {
				return *(_Tparent::ptr);
			}else {
				return *(_Tparent::ptr+1);
			}
		}
	int getCol() const {return col;}
};

class SchurDecomp;
class SchurDecompZero;

class QuasiTriangular : public SqSylvMatrix {
public:
	typedef _column_iter<const double&, const double*> const_col_iter;
	typedef _column_iter<double&, double*> col_iter;
	typedef _row_iter<const double&, const double*> const_row_iter;
	typedef _row_iter<double&, double*> row_iter;	
	typedef Diagonal::const_diag_iter const_diag_iter;
	typedef Diagonal::diag_iter diag_iter;
protected:
	Diagonal diagonal;
public:
	QuasiTriangular(const double* d, int d_size);
	QuasiTriangular(double r, const QuasiTriangular& t);
	QuasiTriangular(double r, const QuasiTriangular& t,
					double rr, const QuasiTriangular& tt);
	QuasiTriangular(int p, const QuasiTriangular& t);
	QuasiTriangular(const SchurDecomp& decomp);
	QuasiTriangular(const SchurDecompZero& decomp);
	QuasiTriangular(const QuasiTriangular& t);
	virtual ~QuasiTriangular();
	const Diagonal& getDiagonal() const {return diagonal;}
	int getNumOffdiagonal() const;
	void swapDiagLogically(diag_iter it);
	void checkDiagConsistency(diag_iter it);
	double getAverageDiagSize(diag_iter start, diag_iter end);
	diag_iter findClosestDiagBlock(diag_iter start, diag_iter end, double a);
	diag_iter findNextLargerBlock(diag_iter start, diag_iter end, double a);


	/* (I+T)y = x, y-->x  */
	virtual void solvePre(Vector& x, double& eig_min);
	/* (I+T')y = x, y-->x */
	virtual void solvePreTrans(Vector& x, double& eig_min);
	/* (I+T)x = b */
	virtual void solve(Vector& x, const ConstVector& b, double& eig_min);
	/* (I+T')x = b */
	virtual void solveTrans(Vector& x, const ConstVector& b, double& eig_min);
	/* x = Tb */
	virtual void multVec(Vector& x, const ConstVector& b) const;
	/* x = T'b */
	virtual void multVecTrans(Vector& x, const ConstVector& b) const;
	/* x = x + Tb */
	virtual void multaVec(Vector& x, const ConstVector& b) const;
	/* x = x + T'b */
	virtual void multaVecTrans(Vector& x, const ConstVector& b) const;
	/* x = (T\otimes I)x */
	virtual void multKron(KronVector& x) const;
	/* x = (T'\otimes I)x */
	virtual void multKronTrans(KronVector& x) const;
	/* A = T*A */
	virtual void multLeftOther(GeneralMatrix& a) const;
	/* A = T'*A */
	virtual void multLeftOtherTrans(GeneralMatrix& a) const;

	const_diag_iter diag_begin() const
		{return diagonal.begin();}
	diag_iter diag_begin()
		{return diagonal.begin();}
	const_diag_iter diag_end() const 
		{return diagonal.end();}
	diag_iter diag_end()
		{return diagonal.end();}

	/* iterators for off diagonal elements */
	virtual const_col_iter col_begin(const DiagonalBlock& b) const;
	virtual col_iter col_begin(const DiagonalBlock& b);
	virtual const_row_iter row_begin(const DiagonalBlock& b) const;
	virtual row_iter row_begin(const DiagonalBlock& b);
	virtual const_col_iter col_end(const DiagonalBlock& b) const;
	virtual col_iter col_end(const DiagonalBlock& b);
	virtual const_row_iter row_end(const DiagonalBlock& b) const;
	virtual row_iter row_end(const DiagonalBlock& b);

	/* clone */
	virtual QuasiTriangular* clone() const
		{return new QuasiTriangular(*this);}
	virtual QuasiTriangular* clone(int p, const QuasiTriangular& t) const
		{return new QuasiTriangular(p, t);}
	virtual QuasiTriangular* clone(double r) const
		{return new QuasiTriangular(r, *this);}
	virtual QuasiTriangular* clone(double r, double rr, const QuasiTriangular& tt) const
		{return new QuasiTriangular(r, *this, rr, tt);}
protected:
	void setMatrix(double r, const QuasiTriangular& t);
	void addMatrix(double r, const QuasiTriangular& t);
private:
	void addUnit();
	/* x = x + (T\otimes I)b */
	void multaKron(KronVector& x, const ConstKronVector& b) const;
	/* x = x + (T'\otimes I)b */
	void multaKronTrans(KronVector& x, const ConstKronVector& b) const;
	/* implementation via iterators, useful for large matrices */
	void setMatrixViaIter(double r, const QuasiTriangular& t);
	void addMatrixViaIter(double r, const QuasiTriangular& t);	
	/* hide noneffective implementations of parents */
	void multsVec(Vector& x, const ConstVector& d) const;
	void multsVecTrans(Vector& x, const ConstVector& d) const;
};

#endif /* QUASI_TRIANGULAR_H */


// Local Variables:
// mode:C++
// End:
