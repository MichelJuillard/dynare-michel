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

#include "RandomWalkMetropolisHastings.hh"

double
RandomWalkMetropolisHastings::compute(Vector &mhLogPostDens, MatrixView &mhParams, Matrix &steadyState,
                                      Vector &estParams2, Vector &deepParams, const MatrixConstView &data, Matrix &Q, Matrix &H,
                                      const size_t presampleStart, int &info, const size_t nMHruns, const Matrix &Jscale,
                                      const Matrix &D, LogPosteriorDensity &lpd, Prior &drawDistribution)
{
  double logpost, newLogpost;
  size_t accepted = 0;
  parDraw = estParams2;
  blas::gemm("N", "N", 1.0, D, Jscale, 1.0, Dscale);
  for (size_t run = 0; run < nMHruns; ++run)
    {
      randMultiVar(drawDistribution, newParDraw, parDraw, Dscale, parDraw.getSize());
      try
        {
          newLogpost = lpd.compute(steadyState, newParDraw, deepParams, data, Q, H, presampleStart, info);
        }
      catch(...)
        {
          newLogpost = -INFINITY;
        }
      if ((newLogpost > -INFINITY) && log(uniform.drand()) < newLogpost-logpost)
        {
          mat::get_col(mhParams, run) = newParDraw;
          parDraw = newParDraw;
          mhLogPostDens(run) = newLogpost;
          logpost = newLogpost;
          accepted++;
        }
      else
        {
          mat::get_col(mhParams, run) = parDraw;
          mhLogPostDens(run) = logpost;
        }
    }
  return accepted/nMHruns;
}

void
RandomWalkMetropolisHastings::saveDraws(const std::string &modName, const std::string &suffix, const MatrixView &Draws, const size_t block)
{

}

