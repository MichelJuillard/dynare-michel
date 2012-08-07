/*
 * Oct-file for bringing MATLAB's linsolve function to Octave.
 *
 * The implementation is incomplete:
 * - it only knows about the TRANSA, LT, UT and SYM options
 * - it only works with square matrices
 * - it only works on double matrices (no single precision or complex)
 *
 * Written by Sebastien Villemot <sebastien.villemot@ens.fr>.
 */

/*
 * Copyright (C) 2012 Dynare Team
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


#include <octave/oct.h>
#include <octave/ov-struct.h>

DEFUN_DLD(linsolve, args, nargout, "-*- texinfo -*-\n\
@deftypefn {Loadable Function} @var{x} = linsolve (@var{a}, @var{b})\n\
@deftypefnx {Loadable Function} [ @var{x}, @var{r} ] = linsolve (@var{a}, @var{b})\n\
@deftypefnx {Loadable Function} @var{x} = linsolve (@var{a}, @var{b}, @var{options})\n\
@deftypefnx {Loadable Function} [ @var{x}, @var{r} ] = linsolve (@var{a}, @var{b}, @var{options})\n\
\n\
Solves the linear system @math{A*X = B} and returns @var{X}.\n\
\n\
Alternatively, if @var{options} is provided and has a field @code{TRANSA} equal \
to @code{true}, then it solves the system @math{A'*X = B}.\n\
\n\
Also, the @code{LT} field of @var{options} (resp. the @code{UT} field) can be set \
to @code{true} to indicate that the matrix @var{a} is lower (resp. upper) \
triangular; similarly, the @code{SYM} field can be set to @code{true} to \
indicate that the matrix is symmetric.\n\
\n\
If requested, @var{r} will contain the reciprocal condition number.\n\
@end deftypefn\n\
")
{
  int nargin = args.length();
  octave_value_list retval;

  if (nargin > 3 || nargin < 2 || nargout > 2)
    {
      print_usage();
      return retval;
    }

  Matrix A = args(0).matrix_value();
  Matrix B = args(1).matrix_value();
  if (error_state)
    return retval;

  dim_vector dimA = A.dims();
  dim_vector dimB = B.dims();

  if (dimA(0) != dimB(0))
    {
      error("linsolve: must have same number of lines in A and B");
      return retval;
    }

  if (dimA(0) != dimA(1))
    {
      error("linsolve: rectangular A not yet supported");
      return retval;
    }

  
  MatrixType typA;
  typA.mark_as_full();
  
  bool transa = false;
  if (nargin == 3)
    {
      octave_scalar_map opts = args(2).scalar_map_value();
      if (error_state)
        return retval;

      octave_value tmp = opts.contents("TRANSA");
      transa = tmp.is_defined() && tmp.bool_matrix_value().elem(0);

      tmp = opts.contents("UT");
      if (tmp.is_defined() && tmp.bool_matrix_value().elem(0))
        typA.mark_as_upper_triangular();

      tmp = opts.contents("LT");
      if (tmp.is_defined() && tmp.bool_matrix_value().elem(0))
        typA.mark_as_lower_triangular();

      tmp = opts.contents("SYM");
      if (tmp.is_defined() && tmp.bool_matrix_value().elem(0))
        typA.mark_as_symmetric();
    }

  double rcond;
  octave_idx_type info;

  retval(0) = A.solve(typA, B, info, rcond, NULL, true, transa ? blas_trans : blas_no_trans);

  if (nargout == 2)
    retval(1) = rcond;

  return retval;
}
