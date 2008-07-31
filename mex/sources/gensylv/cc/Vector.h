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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/Vector.h,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef VECTOR_H
#define VECTOR_H

/* NOTE! Vector and ConstVector have not common super class in order
 * to avoid running virtual method invokation mechanism. Some
 * members, and methods are thus duplicated */ 

class GeneralMatrix;

class BaseVector {
protected:
	int len;
	int s;
	double* data;
public:
	BaseVector(int l, int si, double* d) : len(l), s(si), data(d) {}
	BaseVector(const BaseVector& v) : len(v.len), s(v.s), data(v.data) {}
	const BaseVector& operator=(const BaseVector& v)
		{len = v.len; s = v.s; data = v.data; return *this;}
	double& operator[](int i)
		{return data[s*i];}
	const double& operator[](int i) const
		{return data[s*i];}
	const double* base() const
		{return data;}
	double* base()
		{return data;}
	int length() const
		{return len;}
	int skip() const
		{return s;}
};

class ConstVector;

class Vector : public BaseVector {
protected:
	bool destroy;
public:
	Vector() : BaseVector(0, 1, (double*)0), destroy(false) {} 
	Vector(int l) : BaseVector(l, 1, new double[l]), destroy(true) {}
	Vector(Vector& v) : BaseVector(v.length(), v.skip(), v.base()), destroy(false) {}
	Vector(const Vector& v);
	Vector(const ConstVector& v);
	Vector(const double* d, int l);
	Vector(double* d, int l);
	Vector(double* d, int skip, int l); 
	Vector(Vector& v, int off, int l);
	Vector(const Vector& v, int off, int l);
	Vector(GeneralMatrix& m, int col);
	Vector(int row, GeneralMatrix& m);
	const Vector& operator=(Vector& v);
	const Vector& operator=(const Vector& v);
	const Vector& operator=(const ConstVector& v);
	virtual ~Vector();
	void zeros();
	void nans();
	void infs();
	bool toBeDestroyed() const {return destroy;}
	void rotatePair(double alpha, double beta1, double beta2, int i);
	void add(double r, const Vector& v);
	void add(double r, const ConstVector& v);
	void add(const double* z, const Vector& v);
	void add(const double* z, const ConstVector& v);
	void mult(double r);
	double getNorm() const;
	double getMax() const;
	double getNorm1() const;
	double dot(const Vector& y) const;
	bool isFinite() const;
	void print() const;

	/* multiplies | alpha -beta1|           |b1|   |x1|
                  |             |\otimes I .|  | = |  |
                  | -beta2 alpha|           |b2|   |x2|
	*/
	static void mult2(double alpha, double beta1, double beta2,
					  Vector& x1, Vector& x2,
					  const Vector& b1, const Vector& b2);
	/* same, but adds instead of set */
	static void mult2a(double alpha, double beta1, double beta2,
					   Vector& x1, Vector& x2,
					   const Vector& b1, const Vector& b2);
	/* same, but subtracts instead of add */
	static void mult2s(double alpha, double beta1, double beta2,
					   Vector& x1, Vector& x2,
					   const Vector& b1, const Vector& b2)
		{mult2a(-alpha, -beta1, -beta2, x1, x2, b1, b2);}
private:
	void copy(const double* d, int inc);
	const Vector& operator=(int); // must not be used (not implemented)
	const Vector& operator=(double); // must not be used (not implemented)
};


class BaseConstVector {
protected:
	int len;
	int s;
	const double* data;
public:
	BaseConstVector(int l, int si, const double* d) : len(l), s(si), data(d) {}
	BaseConstVector(const BaseConstVector& v) : len(v.len), s(v.s), data(v.data) {}
	const BaseConstVector& operator=(const BaseConstVector& v)
		{len = v.len; s = v.s; data = v.data; return *this;}
	const double& operator[](int i) const
		{return data[s*i];}
	const double* base() const
		{return data;}
	int length() const
		{return len;}
	int skip() const
		{return s;}
};

class ConstGeneralMatrix;

class ConstVector : public BaseConstVector {
public:
	ConstVector(const Vector& v) : BaseConstVector(v.length(), v.skip(), v.base()) {}
	ConstVector(const ConstVector& v) : BaseConstVector(v) {}
	ConstVector(const double* d, int l) : BaseConstVector(l, 1, d) {}
	ConstVector(const Vector& v, int off, int l);
	ConstVector(const ConstVector& v, int off, int l);
	ConstVector(const double* d, int skip, int l);
	ConstVector(const ConstGeneralMatrix& m, int col);
	ConstVector(int row, const ConstGeneralMatrix& m);
	
	virtual ~ConstVector() {}
	double getNorm() const;
	double getMax() const;
	double getNorm1() const;
	double dot(const ConstVector& y) const;
	bool isFinite() const;
	void print() const;
};

class ZeroPad {
public:
	static const int length = 16;
private:
	double pad[16];
public:
	ZeroPad();
	const double* getBase() const {return pad;} 
};

#endif /* VECTOR_H */


// Local Variables:
// mode:C++
// End:
