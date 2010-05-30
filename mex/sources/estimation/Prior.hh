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

#include <boost/math/distributions/beta.hpp> // for beta_distribution.
#include <boost/math/distributions/gamma.hpp> // for gamma_distribution.
#include <boost/math/distributions/normal.hpp> // for normal_distribution.
#include <boost/math/distributions/uniform.hpp> // for uniform_distribution.


struct Prior
  {

  public:
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

    Prior(double mean, double mode, double standard, double lower_bound, double upper_bound, double fhp, double shp);
    virtual ~Prior();

    const double mean;
    const double mode;
    const double standard;
    const double lower_bound;
    const double upper_bound;
    /**
    * first  shape parameter
    */
    const double fhp;
    /**
    * second shape parameter
    */
    const double shp;
    virtual pShape getShape()=0; // e.g. = Beta for  beta shape?

    virtual double pdf(double x) // probability density function value for x
      {
      std::cout << "Parent pdf undefined at parent level" << std::endl ;
      return 0.0;
      };

    virtual double drand(double x) // rand for density 
      {
      std::cout << "Parent rand undefined at parent level" << std::endl ;
      return 0.0;
      };

  };

struct BetaPrior: public Prior
  {
  public:
    boost::math::beta_distribution<double> distribution;

    BetaPrior( double mean, double mode, double standard, double lower_bound, double upper_bound, double fhp, double shp)
      : Prior( mean, mode, standard, lower_bound, upper_bound, fhp, shp),
      distribution(fhp, shp)
      {};
    virtual ~BetaPrior(){};
    virtual pShape getShape(){return Prior::Beta;}; // e.g. = Beta for  beta shape?

    virtual double pdf(double x)
      {
      double scalled=x;
      if (lower_bound ||1.0-upper_bound)
        scalled=(x- lower_bound)/(upper_bound- lower_bound);    
      return boost::math::pdf(distribution,scalled);
      };

    virtual double drand(double x) // rand for density 
      {
      return 0.0;
      };
  };

struct GammaPrior: public Prior
  {
  public:
    boost::math::gamma_distribution<double> distribution;

    GammaPrior( double mean, double mode, double standard, 
      double lower_bound, double upper_bound, double fhp, double shp)
      : Prior( mean, mode, standard, lower_bound, upper_bound, fhp, shp),
      distribution(fhp, shp)
      {};
    virtual ~GammaPrior(){};
    virtual pShape getShape(){return Prior::Gamma;}; // e.g. = Beta for  beta shape?
    virtual   double pdf(double x)
      {
      return boost::math::pdf(distribution,x- lower_bound);
      };
    virtual double drand(double x) // rand for density 
      {
      return 0.0;
      };
  };

//  X ~ IG1(s,nu) if X = sqrt(Y) where Y ~ IG2(s,nu) and Y = inv(Z) with Z ~ G(nu/2,2/s) (Gamma distribution) 
// lpdfig1(x,s,n)= lpdfgam(1/(x*x),n/2,2/s)-2*log(x*x)
struct InvGamma1_Prior: public Prior
  {
  public:
    boost::math::gamma_distribution<double> distribution;
    InvGamma1_Prior( double mean, double mode, double standard, 
      double lower_bound, double upper_bound, double fhp, double shp)
      : Prior( mean, mode, standard, lower_bound, upper_bound, fhp, shp),
      distribution(shp/2, 2/fhp)
      {};
    virtual ~InvGamma1_Prior(){};
    virtual pShape getShape(){return Prior::Inv_gamma_1;}; // e.g. = Beta for  beta shape?
    virtual  double pdf(double x)
      {
      double scalled=((x- lower_bound)*(x-lower_bound));
      if (x>lower_bound)
        return boost::math::pdf(distribution,1/scalled) / (scalled*scalled);
      else
        return 0;
      };
    virtual double drand(double x) // rand for density 
      {
      return 0.0;
      };
  };

// If x~InvGamma(a,b) , then  1/x ~Gamma(a,1/b) distribution
struct InvGamma2_Prior: public Prior
  {
  public:
    boost::math::gamma_distribution<double> distribution;

    InvGamma2_Prior( double mean, double mode, double standard, 
      double lower_bound, double upper_bound, double fhp, double shp)
      : Prior( mean, mode, standard, lower_bound, upper_bound, fhp, shp),
      distribution(shp/2, 2/fhp)
      {};
    virtual ~InvGamma2_Prior(){};
    virtual pShape getShape(){return Prior::Inv_gamma_2;}; // e.g. = Beta for  beta shape?

    virtual   double pdf(double x)
      {
      double scalled= x - lower_bound;
      if (scalled>0)
        return boost::math::pdf(distribution,1/scalled) / (scalled*scalled);
      else
        return 0;
      };
    virtual double drand(double x) // rand for density 
      {
      return 0.0;
      };
  };

struct GaussianPrior: public Prior
  {
  public:
    boost::math::normal_distribution<double> distribution;

    GaussianPrior( double mean, double mode, double standard, double lower_bound, double upper_bound, double fhp, double shp)
      : Prior( mean, mode, standard, lower_bound, upper_bound, fhp, shp),
      distribution(fhp, shp) //distribution(mean, standard)
      {};
    virtual ~GaussianPrior(){};
    virtual pShape getShape(){return Prior::Gaussian;}; // e.g. = Beta for  beta shape?
    virtual  double pdf(double x)
      {
      if (x>lower_bound && x<upper_bound )
        return boost::math::pdf(distribution,x);
      else
        return 0;
      };
    virtual double drand(double x) // rand for density 
      {
      return 0.0;
      };
  };

struct UniformPrior: public Prior
  {
  public:
    boost::math::uniform_distribution<double> distribution;

    UniformPrior( double mean, double mode, double standard, double lower_bound, double upper_bound, double fhp, double shp)
      : Prior( mean, mode, standard, lower_bound, upper_bound, fhp, shp),
      distribution(fhp, shp) //distribution(lower_bound, upper_bound)
      {};
    virtual ~UniformPrior(){};
    virtual pShape getShape(){return Prior::Uniform;}; // e.g. = Beta for  beta shape?
    virtual   double pdf(double x)
      {
      if (x>lower_bound && x<upper_bound )
        return boost::math::pdf(distribution,x);
      else
        return 0;
      };
    virtual double drand(double x) // rand for density 
      {
      return 0.0;
      };
  };

#endif // !defined(8D5F562F_C831_43f3_B390_5C4EF4433756__INCLUDED_)
