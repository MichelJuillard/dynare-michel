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
//  Prior.h
//  Implementation of the Class Prior
//  Created on:      02-Feb-2010 13:06:20
///////////////////////////////////////////////////////////

#if !defined(Prior_8D5F562F_C831_43f3_B390_5C4EF4433756__INCLUDED_)
#define Prior_8D5F562F_C831_43f3_B390_5C4EF4433756__INCLUDED_

struct Prior
{

public:
  Prior(int shape, double mean, double mode, double lower_boound,       double upper_boound,    double fhp,     double shp);
  virtual ~Prior();

  //! probablity density functions
  enum pShape
  {
    Beta = 1,
    Gamma = 2,
    Gaussian = 3, // i.e. Normal density
    Inv_gamma_1 = 4, // Inverse gamma (type 1) density
    Uniform = 5,
    Inv_gamma_2 = 6 //Inverse gamma (type 2) density
  };
  const pShape shape;
  const double mean;
  const double mode;
  const double lower_boound;
  const double upper_boound;
  /**
   * first  shape parameter
   */
  const double fhp;
  /**
   * second shape parameter
   */
  const double shp;

};

#endif // !defined(8D5F562F_C831_43f3_B390_5C4EF4433756__INCLUDED_)
