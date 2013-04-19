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

#include "SteadyStateSolver.hh"

SteadyStateSolver::SteadyStateSolver(const std::string &basename, size_t n_endo_arg)
  : static_dll(basename), n_endo(n_endo_arg), residual(n_endo), g1(n_endo)
{
  g1.setAll(0.0); // The static file does not initialize zero elements
}

int
SteadyStateSolver::static_f(const gsl_vector *yy, void *p, gsl_vector *F)
{
  params *pp = (params *) p;
  VectorConstView deepParams(pp->deepParams, pp->n_params, 1);
  MatrixConstView x(pp->x, 1, pp->n_exo, 1);
  
  VectorView y(yy->data, yy->size, yy->stride);
  VectorView residual(F->data, F->size, F->stride);

  pp->static_dll->eval(y, x, deepParams, residual, NULL, NULL);
  
  return GSL_SUCCESS;
}

int
SteadyStateSolver::static_df(const gsl_vector *yy, void *p, gsl_matrix *J)
{
  params *pp = (params *) p;
  VectorConstView deepParams(pp->deepParams, pp->n_params, 1);
  MatrixConstView x(pp->x, 1, pp->n_exo, 1);
  
  VectorView y(yy->data, yy->size, yy->stride);

  pp->static_dll->eval(y, x, deepParams, *pp->residual, pp->g1, NULL);
  
  assert(J->size1 == J->size2 && J->size1 == J->tda);
  MatrixView g1t(J->data, J->size1, J->size2, J->tda);
  mat::transpose(g1t, *pp->g1); // GSL wants row-major order

  return GSL_SUCCESS;
}

int
SteadyStateSolver::static_fdf(const gsl_vector *yy, void *p, gsl_vector *F, gsl_matrix *J)
{
  params *pp = (params *) p;
  VectorConstView deepParams(pp->deepParams, pp->n_params, 1);
  MatrixConstView x(pp->x, 1, pp->n_exo, 1);
  
  VectorView y(yy->data, yy->size, yy->stride);
  VectorView residual(F->data, F->size, F->stride);

  pp->static_dll->eval(y, x, deepParams, residual, pp->g1, NULL);
  
  assert(J->size1 == J->size2 && J->size1 == J->tda);
  MatrixView g1t(J->data, J->size1, J->size2, J->tda);
  mat::transpose(g1t, *pp->g1); // GSL wants row-major order

  return GSL_SUCCESS;
}

