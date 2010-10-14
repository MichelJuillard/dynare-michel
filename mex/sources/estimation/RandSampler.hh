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
//  RandSampler.hh
//  Implementation of the Class RandSampler
//  Created on:      09-Sep-2010 23:25:29
///////////////////////////////////////////////////////////

#if !defined(RS5D7E5E52_2A4F_4f98_9A5C_A7FD8C278E0A__INCLUDED_)
#define RS5D7E5E52_2A4F_4f98_9A5C_A7FD8C278E0A__INCLUDED_

#include "LogPosteriorDensity.hh"

class RandSampler {

public:
  virtual ~RandSampler(){};
  virtual double compute(VectorView &mhLogPostDens, MatrixView &mhParams, Matrix &steadyState,
                         Vector &estParams, Vector &deepParams, const MatrixConstView &data, Matrix &Q, Matrix &H,
                         const size_t presampleStart, int &info, const size_t startDraw, size_t nMHruns, const Matrix &Jscale,
                         LogPosteriorDensity &logPosteriorDensity, Prior &drawDistribution,
                         EstimatedParametersDescription &epd) = 0;
  /**
   * draw = Mean + randn(1,n) * Sigma_upper_chol;
   *
   */
  virtual void randMultiVar(Prior &distribution, Vector &draw, const Vector &mean, const Matrix &Scale, const size_t n);

};

#endif // !defined(5D7E5E52_2A4F_4f98_9A5C_A7FD8C278E0A__INCLUDED_)
