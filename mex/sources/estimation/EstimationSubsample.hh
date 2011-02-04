/*
 * Copyright (C) 2009-2011 Dynare Team
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
//  EstimationSubsample.h
//  Implementation of the Class EstimationSubsample
//  Created on:      02-Feb-2010 13:06:01
///////////////////////////////////////////////////////////

#if !defined(EstSub34FB1351_4066_4993_B384_56AA65BC32C4__INCLUDED_)
#define EstSub34FB1351_4066_4993_B384_56AA65BC32C4__INCLUDED_

#include <cstdlib>
#include <vector>

/**
 * Contains longest common periods for different parameters with different start
 * and end periods:
 *
 * It defines start and end of that period.which has start higher/equal than any
 * individual parameter start and its end lower|equal than any individual
 * parameter end in that perod window.) integrated xparam1 and priors for that
 * period
 *
 * All those EstPeriod-s need to be constructed at some point once at a start of
 * estimation  by EstPeriod costructor  from the information such as invidual
 * parameter start-end periods and the relevant associated priors for those
 * periods. That initial info-set is held in the vector of vectors of
 * EstParamSubSampleDescription
 *
 * it constructs and contains indices of all estimated xparam1 parameters valid
 * for this period in the larger, extended xparam1x vector passed  in.  from the
 * gradient function.
 *
 * The EstimatedParametersDescription. structure (or class with protected
 * elements) integrates all that information and makes it available to (friend ?)
 * class LogPosteriorDensity on as needed basis and updates all parameters (inc. H
 * and Q) when time slot xparam1 is supplied.
 *
 * Time indices follow C convention: first period has index 0.
 */
class EstimationSubsample
{
public:
  EstimationSubsample(size_t startPeriod, size_t endPeriod);
  virtual ~EstimationSubsample();

  size_t startPeriod;
  size_t endPeriod;
  /**
   * indices of all estimated xparam1 parameters valid for this period in the larger,
   * extended xparam1x vector passed  in.  from the gradient function will be added at a later stage of development
   */
  // const std::vector<size_t> subSampleParamIDs;
};

#endif // !defined(34FB1351_4066_4993_B384_56AA65BC32C4__INCLUDED_)
