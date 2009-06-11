/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/Vector.cpp,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */


#include "Vector.h"
#include "GeneralMatrix.h"
#include "SylvException.h"
#include "cppblas.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace std;

ZeroPad zero_pad;

Vector::Vector(const ConstVector& v)
: len(v.length()), s(1), data(new double[len]), destroy(true)
  {
  copy(v.base(), v.skip());
  }

const Vector& Vector::operator=(const Vector& v)
  {
  if (this == &v)
    return *this;
  
  if (v.length() != length()) {
    throw SYLV_MES_EXCEPTION("Attempt to assign vectors with different lengths.");
    }
  if (s == v.s &&
    (data <= v.data && v.data < data+len*s ||
    v.data <= data && data < v.data+v.len*v.s) &&
    (data-v.data) % s == 0) {
    printf("this destroy=%d, v destroy=%d, data-v.data=%d, len=%d\n", destroy, v.destroy, data-v.data, len);
    throw SYLV_MES_EXCEPTION("Attempt to assign overlapping vectors.");
    }
  copy(v.base(), v.skip());
  return *this;
  }

const Vector& Vector::operator=(const ConstVector& v)
  {
  if (v.length() != length()) {
    throw SYLV_MES_EXCEPTION("Attempt to assign vectors with different lengths.");
    }
  if (v.skip() == 1 && skip() == 1 && (
    (base() < v.base() + v.length() && base() >= v.base()) ||
    (base() + length() < v.base() + v.length() &&
			 base() + length() > v.base()))) {
    throw SYLV_MES_EXCEPTION("Attempt to assign overlapping vectors.");
    }
  copy(v.base(), v.skip());
  return *this;
  }

void Vector::copy(const double* d, int inc)
  {
  int n = length();
  int incy = skip();
  BLAS_dcopy(&n, d, &inc, base(), &incy);
  }

Vector::Vector(Vector& v, int off, int l)
: len(l), s(v.skip()), data(v.base()+off*v.skip()), destroy(false)
  {
  if (off < 0 || off + length() > v.length())
    throw SYLV_MES_EXCEPTION("Subvector not contained in supvector.");
  }

Vector::Vector(const Vector& v, int off, int l)
: len(l), s(1), data(new double[len]), destroy(true)
  {
  if (off < 0 || off + length() > v.length())
    throw SYLV_MES_EXCEPTION("Subvector not contained in supvector.");
  copy(v.base()+off*v.skip(), v.skip());
  }

Vector::Vector(GeneralMatrix& m, int col)
: len(m.numRows()), s(1), data(&(m.get(0, col))), destroy(false)
  {
  }

Vector::Vector(int row, GeneralMatrix& m)
: len(m.numCols()), s(m.getLD()), data(&(m.get(row, 0))), destroy(false)
  {
  }

bool Vector::operator==(const Vector& y) const
  {
  return ConstVector(*this) == y;
  }

bool Vector::operator!=(const Vector& y) const
  {
  return ConstVector(*this) != y;
  }

bool Vector::operator<(const Vector& y) const
  {
  return ConstVector(*this) < y;
  }

bool Vector::operator<=(const Vector& y) const
  {
  return ConstVector(*this) <= y;
  }

bool Vector::operator>(const Vector& y) const
  {
  return ConstVector(*this) > y;
  }

bool Vector::operator>=(const Vector& y) const
  {
  return ConstVector(*this) >= y;
  }

void Vector::zeros()
  {
  if (skip() == 1) {
    double* p = base();
    for (int i = 0; i < length()/ZeroPad::length;
			 i++, p += ZeroPad::length)
         memcpy(p, zero_pad.getBase(), sizeof(double)*ZeroPad::length);
       for ( ; p < base()+length(); p++)
         *p = 0.0;
    } else {
    for (int i = 0; i < length(); i++)
      operator[](i) = 0.0;
      }
  }

void Vector::nans()
  {
  for (int i = 0; i < length(); i++)
    operator[](i) = std::numeric_limits<double>::quiet_NaN();
  }

void Vector::infs()
  {
  for (int i = 0; i < length(); i++)
    operator[](i) = std::numeric_limits<double>::infinity();
  }

Vector::~Vector()
  {
  if (destroy) {
    delete [] data;
    }
  }

void Vector::rotatePair(double alpha, double beta1, double beta2, int i)
  {
  double tmp = alpha*operator[](i) - beta1*operator[](i+1);
  operator[](i+1) = alpha*operator[](i+1) - beta2*operator[](i);
  operator[](i) = tmp;
  }

void Vector::add(double r, const Vector& v)
  {
  add(r, ConstVector(v));
  }

void Vector::add(double r, const ConstVector& v)
  {
  int n = length();
  int incx = v.skip();
  int incy = skip();
  BLAS_daxpy(&n, &r, v.base(), &incx, base(), &incy);
  }

void Vector::add(const double* z, const Vector& v)
  {
  add(z, ConstVector(v));
  }

void Vector::add(const double* z, const ConstVector& v)
  {
  int n = length()/2;
  int incx = v.skip();
  int incy = skip();
  BLAS_zaxpy(&n, z, v.base(), &incx, base(), &incy);
  }

void Vector::mult(double r)
  {
  int n = length();
  int incx = skip();
  BLAS_dscal(&n, &r, base(), &incx);
  }

void Vector::mult2(double alpha, double beta1, double beta2,
                   Vector& x1, Vector& x2,
                   const Vector& b1, const Vector& b2)
  {
  x1.zeros();
  x2.zeros();
  mult2a(alpha, beta1, beta2, x1, x2, b1, b2);
  }

void Vector::mult2a(double alpha, double beta1, double beta2,
                    Vector& x1, Vector& x2,
                    const Vector& b1, const Vector& b2)
  {
  x1.add(alpha, b1);
  x1.add(-beta1, b2);
  x2.add(alpha, b2);
  x2.add(-beta2, b1);
  }

double Vector::getNorm() const
  {
  ConstVector v(*this);
  return v.getNorm();
  }

double Vector::getMax() const
  {
  ConstVector v(*this);
  return v.getMax();
  }

double Vector::getNorm1() const
  {
  ConstVector v(*this);
  return v.getNorm1();
  }

double Vector::dot(const Vector& y) const
  {
  return ConstVector(*this).dot(ConstVector(y));
  }

bool Vector::isFinite() const
  {
  return (ConstVector(*this)).isFinite();
  }

void Vector::print() const
  {
  for (int i = 0; i < length(); i++) {
    printf("%d\t%8.4g\n", i, operator[](i));
    }
  }


ConstVector::ConstVector(const Vector& v, int off, int l) 
: BaseConstVector(l, v.skip(), v.base() + v.skip()*off)
  {
  if (off < 0 || off + length() > v.length()) {
    throw SYLV_MES_EXCEPTION("Subvector not contained in supvector.");
    }
  }

ConstVector::ConstVector(const ConstVector& v, int off, int l) 
: BaseConstVector(l, v.skip(), v.base() + v.skip()*off)
  {
  if (off < 0 || off + length() > v.length()) {
    throw SYLV_MES_EXCEPTION("Subvector not contained in supvector.");
    }
  }

ConstVector::ConstVector(const double* d, int skip, int l)
: BaseConstVector(l, skip, d)
  {
  }

ConstVector::ConstVector(const ConstGeneralMatrix& m, int col)
: BaseConstVector(m.numRows(), 1, &(m.get(0, col)))
  {
  }

ConstVector::ConstVector(int row, const ConstGeneralMatrix& m)
: BaseConstVector(m.numCols(), m.getLD(), &(m.get(row, 0)))
  {
  }

bool ConstVector::operator==(const ConstVector& y) const
  {
  if (length() != y.length())
    return false;
  if (length() == 0)
    return true;
  int i = 0;
  while (i < length() && operator[](i) == y[i])
    i++;
  return i == length();
  }

bool ConstVector::operator<(const ConstVector& y) const
  {
  int i = std::min(length(), y.length());
  int ii = 0;
  while (ii < i && operator[](ii) == y[ii])
    ii++;
  if (ii < i)
    return operator[](ii) < y[ii];
  else
    return length() < y.length();
  }

double ConstVector::getNorm() const
  {
  double s = 0;
  for (int i = 0; i < length(); i++) {
    s+=operator[](i)*operator[](i);
    }
  return sqrt(s);
  }

double ConstVector::getMax() const
  {
  double r = 0;
  for (int i = 0; i < length(); i++) {
    if (abs(operator[](i))>r)
      r = abs(operator[](i));
    }
  return r;
  }

double ConstVector::getNorm1() const
  {
  double norm = 0.0;
  for (int i = 0; i < length(); i++) {
    norm += abs(operator[](i));
    }
  return norm;
  }

double ConstVector::dot(const ConstVector& y) const
  {
  if (length() != y.length())
    throw SYLV_MES_EXCEPTION("Vector has different length in ConstVector::dot.");
  int n = length();
  int incx = skip();
  int incy = y.skip();
  return BLAS_ddot(&n, base(), &incx, y.base(), &incy);
  }

bool ConstVector::isFinite() const
  {
  int i = 0;
  while (i < length() && isfinite(operator[](i)))
    i++;
  return i == length();
  }

void ConstVector::print() const
  {
  for (int i = 0; i < length(); i++) {
    printf("%d\t%8.4g\n", i, operator[](i));
    }
  }


ZeroPad::ZeroPad()
  {
  for (int i = 0; i < length; i++)
    pad[i] = 0.0;
  }
