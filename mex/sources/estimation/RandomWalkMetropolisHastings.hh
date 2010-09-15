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

class RandomWalkMetropolisHastings : public RandSampler {

private:
  UniformPrior uniform;
  Matrix Dscale;
  Vector parDraw, newParDraw;

public:
  RandomWalkMetropolisHastings(size_t size) :
    uniform(0.0, 0.0, 0.0, 1.0, 0.0, 1.0),
    Dscale(size), parDraw(size), newParDraw(size)
  {
  };
  virtual ~RandomWalkMetropolisHastings(){};
  double compute(Vector &mhLogPostDens, MatrixView &mhParams, Matrix &steadyState,
                 Vector &estParams2, Vector &deepParams, const MatrixConstView &data, Matrix &Q, Matrix &H,
                 const size_t presampleStart, int &info, const size_t nMHruns, const Matrix &Jscale,
                 const Matrix &D, LogPosteriorDensity &logPosteriorDensity, Prior &drawDistribution);
//	void randMultiVar(const Prior &distribution, Vector &draw, const Matrix &Scale, const Vector &mean, const size_t n = 0);
  void saveDraws(const std::string &modName, const std::string &suffix, const MatrixView &Draws, const size_t block);

};

#endif // !defined(A6BBC5E0_598E_4863_B7FF_E87320056B80__INCLUDED_)
