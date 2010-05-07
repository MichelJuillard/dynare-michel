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
//  EstimatedParameter.h
//  Implementation of the Class EstimatedParameter
//  Created on:      02-Feb-2010 13:06:35
///////////////////////////////////////////////////////////

#if !defined(D879C8AE_5B69_4fc3_83BD_FA5A99030ECF__INCLUDED_)
#define D879C8AE_5B69_4fc3_83BD_FA5A99030ECF__INCLUDED_

#include "Prior.hh"
#include <cstdlib>
#include <vector>
struct EstimatedParameter
{
public:
  // parameter types
  enum pType
  {
    shock_SD = 1, // standard deviation of a structural shock
    measureErr_SD = 2, // standard deviation of a measurement error
    shock_Corr = 3, // correlation betwwen two structural shocks
    measureErr_Corr = 4, // correlation between two measurement errors
    deepPar = 5 // deep parameter
  };

  EstimatedParameter(const EstimatedParameter::pType type,
                     size_t ID1, size_t ID2, const std::vector<size_t> &subSampleIDs,
                     double lower_bound, double upper_bound, Prior prior
                     );
  virtual ~EstimatedParameter();

  const enum pType ptype;
  const size_t ID1;
  const size_t ID2;
  const double lower_bound;
  const double upper_bound;
  const Prior prior;
  const std::vector<size_t> subSampleIDs;
};

#endif // !defined(D879C8AE_5B69_4fc3_83BD_FA5A99030ECF__INCLUDED_)
