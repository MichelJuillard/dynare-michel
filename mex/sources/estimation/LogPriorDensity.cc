/*
 * Copyright (C) 2009-2012 Dynare Team
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
//  LogPriorDensity.cpp
//  Implementation of the Class LogPriorDensity
//  Created on:      10-Feb-2010 20:56:08
///////////////////////////////////////////////////////////

#include "LogPriorDensity.hh"
LogPriorDensity::~LogPriorDensity()
{
};

LogPriorDensity::LogPriorDensity(EstimatedParametersDescription &estParsDesc_arg) :
  estParsDesc(estParsDesc_arg)
{
};

/**
 * Return random number for prior fromits distribution
 */
void
LogPriorDensity::computeNewParams(Vector &ep)
{

}
