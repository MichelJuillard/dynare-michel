/*
 * Copyright (C) 2009-2010 Dynare Team
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

///////////////////////////////////////////////////////////
//  Proposal.cpp
//  Implementation of the Class Proposal
//  Created on:      15-Dec-2010 12:43:49
///////////////////////////////////////////////////////////

#include "Proposal.hh"

Proposal::Proposal(const VectorConstView &vJscale, const MatrixConstView &covariance) :
  len(covariance.getCols()),
  covarianceCholeskyDecomposition(len),
  newDraw(len),
  uniform_rng_type(0, 1), // uniform random number generator distribution type
  uniformVrng(base_rng, uniform_rng_type), // uniform random variate_generator
  normal_rng_type(0, 1), // normal random number generator distribution type (mean, standard)
  normalVrng(base_rng, normal_rng_type) // normal random variate_generator
{
  Matrix Jscale(len);
  Matrix DD(len);
  DD = covariance;

  lapack::choleskyDecomp(DD, "U");
  Jscale.setAll(0.0);
  for (size_t i = 0; i < len; i++)
    Jscale(i, i) = vJscale(i);
  blas::gemm("N", "N", 1.0, DD, Jscale, 0.0, covarianceCholeskyDecomposition);
}

void
Proposal::draw(Vector &mean, Vector &draw)
{
  assert(len == draw.getSize());
  assert(len == mean.getSize());

  draw = mean;
  for (size_t i = 0; i < len; ++i)
    newDraw(i) =  normalVrng();
  blas::gemv("T", 1.0, covarianceCholeskyDecomposition, newDraw, 1.0, draw);

}

Matrix &
Proposal::getVar()
{
  return covarianceCholeskyDecomposition;
}

/**
 * set or get if null arguments
 */
int
Proposal::seed()
{
  return curSeed;
}

void
Proposal::seed(int newSeed)
{
  curSeed = newSeed;
  base_rng.seed(curSeed);
}

/**
 * currently returns uniform for MH sampler,
 */
double
Proposal::selectionTestDraw()
{
  return uniformVrng();
}

