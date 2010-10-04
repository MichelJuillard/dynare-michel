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
RandomWalkMetropolisHastings::compute(VectorView &mhLogPostDens, MatrixView &mhParams, Matrix &steadyState,
                                      Vector &estParams, Vector &deepParams, const MatrixConstView &data, Matrix &Q, Matrix &H,
                                      const size_t presampleStart, int &info, const size_t nMHruns, const Matrix &Dscale,
                                      LogPosteriorDensity &lpd, Prior &drawDistribution, EstimatedParametersDescription &epd)
{
  bool overbound;
  double newLogpost, logpost;
  size_t count, accepted = 0;
  parDraw = estParams;
  logpost =  - lpd.compute(steadyState, estParams, deepParams, data, Q, H, presampleStart, info);
  for (size_t run = 0; run < nMHruns; ++run)
    {
      overbound=false;
      randMultiVar(drawDistribution, newParDraw, parDraw, Dscale, parDraw.getSize());
      for (count=0;count<parDraw.getSize();++count)
        {
          overbound=(newParDraw(count) <  epd.estParams[count].lower_bound || newParDraw(count) > epd.estParams[count].upper_bound );
          if (overbound) 
            {
              newLogpost = -INFINITY;
              break;
            }
        }
      if (!overbound)
        {
          try
            {
              newLogpost = - lpd.compute(steadyState, newParDraw, deepParams, data, Q, H, presampleStart, info);
            }
          catch(...)
            {
              newLogpost = -INFINITY;
            }
        }
      if ((newLogpost > -INFINITY) && log(uniform.drand()) < newLogpost-logpost)
        {
          mat::get_row(mhParams, run) = newParDraw;
          parDraw = newParDraw;
          mhLogPostDens(run) = newLogpost;
          logpost = newLogpost;
          accepted++;
        }
      else
        {
          mat::get_row(mhParams, run) = parDraw;
          mhLogPostDens(run) = logpost;
        }
    }
  return (double) accepted/nMHruns;
}

