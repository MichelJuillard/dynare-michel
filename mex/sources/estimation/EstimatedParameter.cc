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
//  EstimatedParameter.cpp
//  Implementation of the Class EstimatedParameter
//  Created on:      02-Feb-2010 13:06:35
///////////////////////////////////////////////////////////

#include "EstimatedParameter.hh"

EstimatedParameter::EstimatedParameter(const EstimatedParameter::pType type_arg,
                                       size_t ID1_arg, size_t ID2_arg, const std::vector<size_t> &subSampleIDs_arg,
                                       double lower_bound_arg, double upper_bound_arg, Prior *prior_arg) :
  ptype(type_arg), ID1(ID1_arg), ID2(ID2_arg),
  lower_bound(lower_bound_arg), upper_bound(upper_bound_arg), prior(prior_arg),
  subSampleIDs(subSampleIDs_arg)
{
}

EstimatedParameter::~EstimatedParameter()
{
}

