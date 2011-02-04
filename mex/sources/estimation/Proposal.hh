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
//  Proposal.hh
//  Implementation of the Class Proposal
//  Created on:      15-Dec-2010 12:43:49
///////////////////////////////////////////////////////////

#if !defined(CABFAB46_2FA5_4178_8D1C_05DFEAFB4570__INCLUDED_)
#define CABFAB46_2FA5_4178_8D1C_05DFEAFB4570__INCLUDED_

/**
 * Proposal class will then have the common, seed initialised base generator
 * (currently congruental rand48) and member functions such as seed, reset (to
 * initial,default value) and have all rand generators we need that generate from
 * that common base: single uniform and normal (either single or multivariate
 * determined by  the size of the variance matrix)  for now.
 *
 * The rng interfaces will be - a single, selectionTest() currently uniform for MH
 * sampler, draw(mean, newDraw) - single or multivariate parameters (currently
 * from Normal(mean,Sigma))  and poss. also the seed state (set and get), the
 * normal Chol_variance_decomp (constructed using Chol Decomp of init Variance
 * matrix) , as core members  .This class will handle the main boost random rng
 * intricacies. See enclosed updated diagram (if ok I will upload it)
 *
 */

#include "Matrix.hh"
#include "BlasBindings.hh"
#include "LapackBindings.hh"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

// typedef for base_uniform_generator_type
// rend48 seems better than the basic minstd_rand  but
// one my later try ecuyer1988 or taus88 which are better but not yet supported in the Boost versionwe use
typedef boost::rand48 base_uniform_generator_type;

class Proposal
{

public:
  Proposal();

public:
  Proposal(const VectorConstView &vJscale, const MatrixConstView &covariance);
  virtual ~Proposal() {};
  virtual void draw(Vector &mean, Vector &draw);
  virtual Matrix&getVar();
  virtual int seed();
  virtual void seed(int seedInit);
  virtual double selectionTestDraw();

private:
  size_t len;
  Matrix covarianceCholeskyDecomposition;
  /**
   * Vector of new draws
   *
   */
  Vector newDraw;

  base_uniform_generator_type base_rng;
  boost::uniform_real<> uniform_rng_type;
  boost::variate_generator<base_uniform_generator_type &,  boost::uniform_real<> > uniformVrng;
  boost::normal_distribution<double> normal_rng_type;
  boost::variate_generator<base_uniform_generator_type &, boost::normal_distribution<double> > normalVrng;

  int curSeed;

};

#endif // !defined(CABFAB46_2FA5_4178_8D1C_05DFEAFB4570__INCLUDED_)
