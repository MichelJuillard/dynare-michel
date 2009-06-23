/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/Vector.h,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef VECTOR_H
#define VECTOR_H

/* NOTE! Vector and ConstVector have not common super class in order
 * to avoid running virtual method invokation mechanism. Some
 * members, and methods are thus duplicated */ 


#ifdef MATLAB
#include "mex.h"
#define printf mexPrintf
#endif

#include <stdio.h>

class GeneralMatrix;
class ConstVector;

class Vector {
protected:
	int len;
	int s;
	double* data;
	bool destroy;
public:
	Vector() : len(0), s(1), data(0), destroy(false) {}
	Vector(int l) : len(l), s(1), data(new double[l]), destroy(true) {}
	Vector(Vector& v) : len(v.length()), s(v.skip()), data(v.base()), destroy(false) {}
	Vector(const Vector& v)
		: len(v.length()), s(1), data(new double[len]), destroy(true)
		{copy(v.base(), v.skip());}
	Vector(const ConstVector& v);
	Vector(const double* d, int l)
		: len(l), s(1), data(new double[len]), destroy(true)
		{copy(d, 1);}
	Vector(double* d, int l)
		: len(l), s(1), data(d), destroy(false) {}
	Vector(double* d, int skip, int l)
		: len(l), s(skip), data(d), destroy(false) {}
	Vector(Vector& v, int off, int l);
	Vector(const Vector& v, int off, int l);
	Vector(GeneralMatrix& m, int col);
	Vector(int row, GeneralMatrix& m);
	const Vector& operator=(const Vector& v);
	const Vector& operator=(const ConstVector& v);
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

	/** Exact equality. */
	bool operator==(const Vector& y) const;
	bool operator!=(const Vector& y) const;
	/** Lexicographic ordering. */
	bool operator<(const Vector& y) const;
	bool operator<=(const Vector& y) const;
	bool operator>(const Vector& y) const;
	bool operator>=(const Vector& y) const;

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
	/** Exact equality. */
	bool operator==(const ConstVector& y) const;
	bool operator!=(const ConstVector& y) const
		{return ! operator==(y);}
	/** Lexicographic ordering. */
	bool operator<(const ConstVector& y) const;
	bool operator<=(const ConstVector& y) const
		{return operator<(y) || operator==(y);}
	bool operator>(const ConstVector& y) const
		{return ! operator<=(y);}
	bool operator>=(const ConstVector& y) const
		{return ! operator<(y);}

	double getNorm() const;
	double getMax() const;
	double getNorm1() const;
	double dot(const ConstVector& y) const;
	bool isFinite() const;
	void print() const;
};

class ZeroPad {
public:
	//static const int length = 16;
	enum { length = 16};
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
