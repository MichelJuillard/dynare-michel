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
//  DetrendData.h
//  Implementation of the Class DetrendData
//  Created on:      02-Feb-2010 13:01:15
///////////////////////////////////////////////////////////

#if !defined(DetrendData_312823A1_6248_4af0_B204_DB22F1237E9B__INCLUDED_)
#define DetrendData_312823A1_6248_4af0_B204_DB22F1237E9B__INCLUDED_

#include "Matrix.hh"

class DetrendData
{

public:
  virtual ~DetrendData(){};
  DetrendData(const bool logLinear); // add later Vector& trendCoeff);
  void detrend(const VectorView &SteadyState, const MatrixConstView &dataView, MatrixView &detrendedDataView);

private:
  const bool logLinear;
  //Vector trendCoeff;

};

#endif // !defined(312823A1_6248_4af0_B204_DB22F1237E9B__INCLUDED_)
