/*
 * Copyright (C) 2010-2011 Dynare Team
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

#include "Vector.hh"

#include <cstring>
#include <cassert>

Vector::Vector(size_t size_arg) : size(size_arg)
{
  data = new double[size];
}

Vector::Vector(const Vector &arg) : size(arg.size)
{
  data = new double[size];
  memcpy(data, arg.data, size*sizeof(double));
}

Vector::~Vector()
{
  delete[] data;
}

Vector &
Vector::operator=(const Vector &arg)
{
  assert(size == arg.size);
  memcpy(data, arg.data, size*sizeof(double));
  return *this;
}

std::ostream &
operator<<(std::ostream &out, const Vector &V)
{
  vec::print(out, V);
  return out;
}

VectorView::VectorView(double *data_arg, size_t size_arg, size_t stride_arg) :
  data(data_arg), size(size_arg), stride(stride_arg)
{
}

VectorView::VectorView(Vector &arg, size_t offset, size_t size_arg) :
  data(arg.getData() + offset*arg.getStride()), size(size_arg), stride(arg.getStride())
{
  assert(offset < arg.getSize() && offset + size <= arg.getSize());
}

VectorView::VectorView(VectorView &arg, size_t offset, size_t size_arg) :
  data(arg.getData() + offset*arg.getStride()), size(size_arg), stride(arg.getStride())
{
  assert(offset < arg.getSize() && offset + size <= arg.getSize());
}

VectorView &
VectorView::operator=(const VectorView &arg)
{
  assert(size == arg.getSize());
  double *p1;
  const double *p2;
  for (p1 = data, p2 = arg.getData(); p1 < data + size*stride; p1 += stride, p2 += arg.getStride())
    *p1 = *p2;
  return *this;
}

VectorConstView::VectorConstView(const double *data_arg, size_t size_arg, size_t stride_arg) :
  data(data_arg), size(size_arg), stride(stride_arg)
{
}

VectorConstView::VectorConstView(const Vector &arg, size_t offset, size_t size_arg) :
  data(arg.getData() + offset*arg.getStride()), size(size_arg), stride(arg.getStride())
{
  assert(offset < arg.getSize() && offset + size <= arg.getSize());
}

VectorConstView::VectorConstView(const VectorView &arg, size_t offset, size_t size_arg) :
  data(arg.getData() + offset*arg.getStride()), size(size_arg), stride(arg.getStride())
{
  assert(offset < arg.getSize() && offset + size <= arg.getSize());
}

VectorConstView::VectorConstView(const VectorConstView &arg, size_t offset, size_t size_arg) :
  data(arg.getData() + offset*arg.getStride()), size(size_arg), stride(arg.getStride())
{
  assert(offset < arg.getSize() && offset + size <= arg.getSize());
}

std::ostream &
operator<<(std::ostream &out, const VectorView &V)
{
  vec::print(out, V);
  return out;
}

std::ostream &
operator<<(std::ostream &out, const VectorConstView &V)
{
  vec::print(out, V);
  return out;
}
