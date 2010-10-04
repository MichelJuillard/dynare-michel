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
//  RandomWalkMetropolisHastings.hh
//  Implementation of the Class RandomWalkMetropolisHastings
//  Created on:      07-Sep-2010 15:21:40
///////////////////////////////////////////////////////////

#if !defined(A6BBC5E0_598E_4863_B7FF_E87320056B80__INCLUDED_)
#define A6BBC5E0_598E_4863_B7FF_E87320056B80__INCLUDED_

#include "RandSampler.hh"

class RandomWalkMetropolisHastings : public RandSampler
{

private:
  UniformPrior uniform;
  Vector parDraw, newParDraw;

public:
  RandomWalkMetropolisHastings(size_t size) :
    uniform(0.0, 0.0, 0.0, 1.0, 0.0, 1.0),
    parDraw(size), newParDraw(size)
  {
  };
  virtual ~RandomWalkMetropolisHastings(){};
  virtual double compute(VectorView &mhLogPostDens, MatrixView &mhParams, Matrix &steadyState,
                         Vector &estParams, Vector &deepParams, const MatrixConstView &data, Matrix &Q, Matrix &H,
                         const size_t presampleStart, int &info, const size_t nMHruns, const Matrix &Jscale,
                         LogPosteriorDensity &logPosteriorDensity, Prior &drawDistribution,
                         EstimatedParametersDescription &epd);
};

#endif // !defined(A6BBC5E0_598E_4863_B7FF_E87320056B80__INCLUDED_)
