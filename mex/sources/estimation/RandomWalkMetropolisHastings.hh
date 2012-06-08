/*
 * Copyright (C) 2010-2012 Dynare Team
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

#include <iostream>
#include <fstream>
#include "LogPosteriorDensity.hh"
#include "Proposal.hh"

class RandomWalkMetropolisHastings
{

private:
  Vector parDraw, newParDraw;

public:
  RandomWalkMetropolisHastings(size_t size) :
    parDraw(size), newParDraw(size)
  {
  };
  virtual ~RandomWalkMetropolisHastings() {};

  template<class VEC1>
  double compute(VectorView &mhLogPostDens, MatrixView &mhParams, VEC1 &steadyState,
		 Vector &estParams, VectorView &deepParams, const MatrixConstView &data, MatrixView &Q, Matrix &H,
		 const size_t presampleStart, int &info, const size_t startDraw, size_t nMHruns,
		 LogPosteriorDensity &lpd, Proposal &pDD, EstimatedParametersDescription &epd)
  {
    //streambuf *likbuf, *drawbuf *backup;
    std::ofstream urandfilestr, drawfilestr;
    urandfilestr.open("urand.csv");
    drawfilestr.open("paramdraws.csv");

    bool overbound;
    double newLogpost, logpost, urand;
    size_t count, accepted = 0;
    parDraw = estParams;

    logpost = -lpd.compute(steadyState, estParams, deepParams, data, Q, H, presampleStart, info);

    for (size_t run = startDraw - 1; run < nMHruns; ++run)
      {
	overbound = false;
	pDD.draw(parDraw, newParDraw);
	for (count = 0; count < parDraw.getSize(); ++count)
	  {
	    overbound = (newParDraw(count) < epd.estParams[count].lower_bound || newParDraw(count) > epd.estParams[count].upper_bound);
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
		newLogpost = -lpd.compute(steadyState, newParDraw, deepParams, data, Q, H, presampleStart, info);
	      }
	    catch (const std::exception &e)
	      {
		throw; // for now handle the system and other errors higher-up
	      }
	    catch (...)
	      {
		newLogpost = -INFINITY;
	      }
	  }
	urand = pDD.selectionTestDraw();
	if ((newLogpost > -INFINITY) && log(urand) < newLogpost-logpost)
	  {
	    parDraw = newParDraw;
	    logpost = newLogpost;
	    accepted++;
	  }
	mat::get_row(mhParams, run) = parDraw;
	mhLogPostDens(run) = logpost;

	urandfilestr << urand << "\n"; //","
	for (size_t c = 0; c < newParDraw.getSize()-1; ++c)
	  drawfilestr << newParDraw(c) << ",";
	drawfilestr <<  newParDraw(newParDraw.getSize()-1) << "\n";
      }

    urandfilestr.close();
    drawfilestr.close();

    return (double) accepted/(nMHruns-startDraw+1);
  };

};

#endif // !defined(A6BBC5E0_598E_4863_B7FF_E87320056B80__INCLUDED_)
