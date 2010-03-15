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
//  DetrendData.cpp
//  Implementation of the Class DetrendData
//  Created on:      02-Feb-2010 13:01:15
///////////////////////////////////////////////////////////

#include "DetrendData.hh"

DetrendData::DetrendData(const bool INlogLinear) //, Vector& INtrendCoeff)
  : logLinear(INlogLinear) //trendCoeff(INtrendCoeff)
{
};

void
DetrendData::detrend(const Vector &SteadyState, const MatrixView &dataView, Matrix &Yiew)
{
};

