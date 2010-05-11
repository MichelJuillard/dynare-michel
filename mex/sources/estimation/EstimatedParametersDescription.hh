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
//  EstimatedParametersDescription.h
//  Implementation of the Class EstimatedParametersDescription
//  Created on:      02-Feb-2010 13:06:46
///////////////////////////////////////////////////////////

#if !defined(E8F2C096_A301_42e8_80FF_A643291BF995__INCLUDED_)
#define E8F2C096_A301_42e8_80FF_A643291BF995__INCLUDED_

#include "EstimationSubsample.hh"
#include "EstimatedParameter.hh"
#include "Prior.hh"

/**
 * all estimation periods held in the vectors of EstPeriod-s are held in the
 * EstimatedParametersDescription together with the time-invariant attributes (e.g.
 * IDs..) for integration..
 *
 * The EstimatedParametersDescription. structure (or class with protected
 * elements) integrates all parametr and time period information and makes it
 * available to (friend ?) class LogPosteriorDensity on as needed basis and
 * updates all parameters (inc. H and Q) when time slot xparam1 is supplied.
 *
 *
 *
 */
class EstimatedParametersDescription
{
public:
  virtual ~EstimatedParametersDescription();
  EstimatedParametersDescription(std::vector<EstimationSubsample> &estSubsamples, std::vector<EstimatedParameter> &estParams);
  std::vector<EstimationSubsample> estSubsamples;
  std::vector<EstimatedParameter> estParams;
};

#endif // !defined(E8F2C096_A301_42e8_80FF_A643291BF995__INCLUDED_)
