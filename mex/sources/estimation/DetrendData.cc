/*
 * Copyright (C) 2009-2013 Dynare Team
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
//  DetrendData.cpp
//  Implementation of the Class DetrendData
//  Created on:      02-Feb-2010 13:01:15
///////////////////////////////////////////////////////////

#include "DetrendData.hh"

DetrendData::DetrendData(const std::vector<size_t> &varobs_arg,
                         bool noconstant_arg)
  : varobs(varobs_arg), noconstant(noconstant_arg)
{
};

void
DetrendData::detrend(const VectorView &SteadyState, const MatrixConstView &dataView,
                     MatrixView &detrendedDataView)
{
  if (noconstant)
    detrendedDataView = dataView;
  else
    {
      for (size_t i = 0; i < varobs.size(); i++)
        for (size_t j = 0; j < dataView.getCols(); j++)
          detrendedDataView(i, j) = dataView(i, j) - SteadyState(varobs[i]);
    }
};

