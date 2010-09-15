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

///////////////////////////////////////////////////////////
//  RandSampler.cpp
//  Implementation of the Class RandSampler
//  Created on:      09-Sep-2010 23:25:30
///////////////////////////////////////////////////////////

#include "RandSampler.hh"

/**
 * draw = Mean + randn(1,n) * Sigma_upper_chol;
 */
void
RandSampler::randMultiVar(Prior &distribution, Vector &draw, const Vector &mean, const Matrix &Scale, const size_t n)
{
  assert(n == draw.getSize());
  assert(n == mean.getSize() || 1 == mean.getSize());
  assert(n == Scale.getRows() || 1 == Scale.getRows());
  assert(n == Scale.getCols() || 1 == Scale.getCols());
  for (size_t i = 0; i < n; ++i)
    draw(i) = mean(1 == mean.getSize() ? 0 : i)
              + distribution.drand()
              *Scale((1 == Scale.getRows() ? 0 : i), (1 == Scale.getCols() ? 0 : i));
}

