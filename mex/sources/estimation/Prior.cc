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
//  Prior.cpp
//  Implementation of the Class Prior
//  Created on:      02-Feb-2010 13:06:20
///////////////////////////////////////////////////////////

#include "Prior.hh"

Prior::Prior(double mean_arg, double standard_arg, double lower_bound_arg,      double upper_bound_arg, double fhp_arg, double shp_arg) :
  mean(mean_arg), standard(standard_arg), lower_bound(lower_bound_arg),       upper_bound(upper_bound_arg),   fhp(fhp_arg), shp(shp_arg)
{

}

Prior::~Prior()
{

}

Prior *
Prior::constructPrior(pShape shape, double mean, double standard, double lower_bound, double upper_bound, double fhp, double shp)
{
  switch (shape)
    {
    case Beta:
      return new BetaPrior(mean, standard, lower_bound, upper_bound, fhp, shp);
    case Gamma:
      return new GammaPrior(mean, standard, lower_bound, upper_bound, fhp, shp);
    case Gaussian:
      return new GaussianPrior(mean, standard, lower_bound, upper_bound, fhp, shp);
    case Inv_gamma_1:
      return new InvGamma1_Prior(mean, standard, lower_bound, upper_bound, fhp, shp);
    case Uniform:
      return new UniformPrior(mean, standard, lower_bound, upper_bound, fhp, shp);
    case Inv_gamma_2:
      return new InvGamma2_Prior(mean, standard, lower_bound, upper_bound, fhp, shp);
    }
  std::cerr << "Unreacheable point" << std::endl;
  exit(EXIT_FAILURE);
}
