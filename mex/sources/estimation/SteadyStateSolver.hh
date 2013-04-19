/*
 * Copyright (C) 2013 Dynare Team
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

#include <string>

#include "Vector.hh"
#include "static_dll.hh"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

class SteadyStateSolver
{
private:
  StaticModelDLL static_dll;
  size_t n_endo;
  Vector residual; // Will be discarded, only used by df()
  Matrix g1; // Temporary buffer for computing transpose

  struct params
  {
    StaticModelDLL *static_dll;
    const double *deepParams;
    size_t n_params;
    const double *x;
    size_t n_exo;
    Vector *residual;
    Matrix *g1;
  };

  static int static_f(const gsl_vector *yy, void *p, gsl_vector *F);
  static int static_df(const gsl_vector *yy, void *p, gsl_matrix *J);
  static int static_fdf(const gsl_vector *yy, void *p, gsl_vector *F, gsl_matrix *J);

  const static double tolerance = 1e-7;
  const static size_t max_iterations = 1000;
public:
  class SteadyStateException
  {
  public:
    std::string message;
    SteadyStateException(const std::string &message_arg) : message(message_arg) 
    {
    }
  };

  SteadyStateSolver(const std::string &basename, size_t n_endo_arg);

  template <class Vec1, class Mat, class Vec2>
  void compute(Vec1 &steadyState, const Mat &Mx, const Vec2 &deepParams) throw (SteadyStateException)
  {
    assert(steadyState.getStride() == 1);
    assert(deepParams.getStride() == 1);
    assert(Mx.getLd() == Mx.getRows());

    assert(steadyState.getSize() == n_endo);

    gsl_vector_view ss = gsl_vector_view_array(steadyState.getData(), n_endo);

    params p = { &static_dll, deepParams.getData(), deepParams.getSize(), Mx.getData(), Mx.getCols(), &residual, &g1 };

    gsl_multiroot_function_fdf f = {&static_f, &static_df, &static_fdf,
                                    n_endo, &p};

    const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj;
    gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, n_endo);

    gsl_multiroot_fdfsolver_set(s, &f, &ss.vector);

    int status;
    size_t iter = 0;

    do
      {
        iter++;

        status = gsl_multiroot_fdfsolver_iterate(s);

        if (status)
          break;

        status = gsl_multiroot_test_residual(s->f, tolerance);
      }
    while(status == GSL_CONTINUE && iter < max_iterations);

    if (status != GSL_SUCCESS)
      throw SteadyStateException(std::string(gsl_strerror(status)));

    gsl_vector_memcpy(&ss.vector, gsl_multiroot_fdfsolver_root(s));

    gsl_multiroot_fdfsolver_free(s);
  }
};


