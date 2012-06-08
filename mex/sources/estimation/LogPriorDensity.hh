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
//  LogPriorDensity.hh
//  Implementation of the Class LogPriorDensity
//  Created on:      10-Feb-2010 20:56:08
///////////////////////////////////////////////////////////

#if !defined(LPD_011FD4CF_17CE_4805_882B_046AA07CF443__INCLUDED_)
#define LPD_011FD4CF_17CE_4805_882B_046AA07CF443__INCLUDED_

//#include <boost/random/variate_generator.hpp>
#include "Vector.hh"
#include "EstimatedParametersDescription.hh"

class LogPriorDensity
{

public:
  LogPriorDensity(EstimatedParametersDescription &estParsDesc);
  virtual ~LogPriorDensity();

  template<class VEC>
  double compute(VEC &ep)
  {
    assert(estParsDesc.estParams.size() == ep.getSize());
    double logPriorDensity = 0;
    for (size_t i = 0; i <  ep.getSize(); ++i)
      {
	logPriorDensity += log(((*(estParsDesc.estParams[i]).prior)).pdf(ep(i)));
	if (std::isinf(fabs(logPriorDensity)))
	  return logPriorDensity;
      }
    return logPriorDensity;
  };

  void computeNewParams(Vector &newParams);

private:
  const EstimatedParametersDescription &estParsDesc;

};

#endif // !defined(011FD4CF_17CE_4805_882B_046AA07CF443__INCLUDED_)
