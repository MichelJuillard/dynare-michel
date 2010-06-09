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

// test for Prior pdfs and basic random number generators

#include "Prior.hh"

int
main(int argc, char **argv)
{
//BetaPrior bp( Prior::Beta, 0.1, 0.1, 0.1, 0.0, 1.0, 5.0, 1.0);
  BetaPrior bp(0.1, 0.1, 0.1, 0.0, 1.0, 5.0, 1.0);
  double pdf = bp.pdf(0.5);
  std::cout << "beta pdf of 5,1, 0.5: "  << std::setw(13) << pdf << std::endl;
  BetaPrior *bpp = new BetaPrior( //Prior::Beta,
    0.1, 0.1, 0.1, 0.0, 1.0, 1.0, 5.0);
  Prior *pp = bpp;
  pdf = (*pp).pdf(0.1);
  std::cout << "Parent (Beta) pdf of 1,5, 0.1: "  << std::setw(13) << pdf << std::endl;

  GammaPrior *gpp = new GammaPrior( //Prior::Beta,
    0.1, 0.1, 0.1, 0.0, 1.0, 1.0, 5.0);
  pp = gpp;
  pdf = (*pp).pdf(0.1);
  std::cout << "Parent (Gamma) pdf of 1,5, 0.1: "  << std::setw(13) << pdf << std::endl;

  UniformPrior up(5, 5, 2, 1, 10, 20, 100);
  std::cout << std::endl << "Uniform prior (5,5,2, 1,10,20,100): "  << std::endl;
  double ur, updf;
  for (int i = 0; i < 10; i++)
    {
      ur = up.drand();
      updf = up.pdf(ur);
      std::cout << "Uniform pdf of : "  << ur << " = " << updf << std::endl;
    }

  GaussianPrior gp(5, 5, 2, 1, 10, 20, 100);
  std::cout << std::endl << "Gaussian prior (5,5,2, 1,10,20,100): "  << std::endl;
  for (int i = 0; i < 10; i++)
    {
      ur = gp.drand();
      updf = gp.pdf(ur);
      std::cout << "Gaussian pdf of : "  << ur << " = " << updf << std::endl;
    }

  UniformPrior up2(5, 5, 2, 1, 100, 20, 100);
  std::cout << std::endl << "Uniform prior (5,5,2, 1,100,20,100): "  << std::endl;
  for (int i = 0; i < 10; i++)
    {
      ur = up2.drand();
      updf = up2.pdf(ur);
      std::cout << "Uniform pdf of : "  << ur << " = " << updf << std::endl;
    }

  GaussianPrior gp2(5, 5, 2, 1, 100, 5, 2);
  std::cout << std::endl << "Gaussian prior (5,5,2, 1,100,5,2): "  << std::endl;
  for (int i = 0; i < 10; i++)
    {
      ur = gp2.drand();
      updf = gp2.pdf(ur);
      std::cout << "Gaussian pdf of : "  << ur << " = " << updf << std::endl;
    }
};
