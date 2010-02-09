/*
 * Copyright (C) 2010 Dynare Team
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

#include "Matrix.hh"

#include <cstring>
#include <cassert>

Matrix::Matrix(size_t rows_arg, size_t cols_arg) : rows(rows_arg), cols(cols_arg)
{
  data = new double[rows*cols];
}

Matrix::Matrix(size_t size_arg) : rows(size_arg), cols(size_arg)
{
  data = new double[rows*cols];
}

Matrix::Matrix(const Matrix &arg) : rows(arg.rows), cols(arg.cols)
{
  data = new double[rows*cols];
  memcpy(data, arg.data, rows*cols*sizeof(double));
}

Matrix::~Matrix()
{
  delete[] data;
}

Matrix &
Matrix::operator= (const Matrix &arg)
{
  assert(rows == arg.rows && cols == arg.cols);
  memcpy(data, arg.data, rows*cols*sizeof(double));
  return *this;
}

std::ostream &
operator<<(std::ostream &out, const Matrix &M)
{
  mat::print(out, M);
  return out;
}

MatrixView::MatrixView(double *data_arg, size_t rows_arg, size_t cols_arg, size_t ld_arg)
  : data(data_arg), rows(rows_arg), cols(cols_arg), ld(ld_arg)
{
}

MatrixView::MatrixView(Matrix &arg, size_t row_offset, size_t col_offset,
                       size_t rows_arg, size_t cols_arg) :
  data(arg.getData() + row_offset + col_offset*arg.getLd()), rows(rows_arg), cols(cols_arg), ld(arg.getLd())
{
  assert(row_offset < arg.getRows()
         && row_offset + rows_arg <= arg.getRows()
         && col_offset < arg.getCols()
         && col_offset + cols_arg <= arg.getCols());
}

std::ostream &
operator<<(std::ostream &out, const MatrixView &M)
{
  mat::print(out, M);
  return out;
}

MatrixView &
MatrixView::operator= (const MatrixView &arg)
{
  assert(rows == arg.getRows() && cols == arg.getCols());
  for (size_t j = 0; j < cols; j++)
    memcpy(data + j*ld, arg.getData() + j*arg.getLd(), rows*sizeof(double));
  return *this;
}

MatrixConstView::MatrixConstView(const double *data_arg, size_t rows_arg, size_t cols_arg, size_t ld_arg)
  : data(data_arg), rows(rows_arg), cols(cols_arg), ld(ld_arg)
{
}

MatrixConstView::MatrixConstView(const Matrix &arg, size_t row_offset, size_t col_offset,
                                 size_t rows_arg, size_t cols_arg) :
  data(arg.getData() + row_offset + col_offset*arg.getLd()), rows(rows_arg), cols(cols_arg), ld(arg.getLd())
{
  assert(row_offset < arg.getRows()
         && row_offset + rows_arg <= arg.getRows()
         && col_offset < arg.getCols()
         && col_offset + cols_arg <= arg.getCols());
}

std::ostream &
operator<<(std::ostream &out, const MatrixConstView &M)
{
  mat::print(out, M);
  return out;
}
