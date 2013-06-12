/*
 * Oct-file for bringing MATLAB's ordschur function to Octave.
 * Simple wrapper around LAPACK's dtrsen.
 * Only supports real (double precision) decomposition.
 * Only selection of eigenvalues with a boolean vector is supported.
 *
 * Written by Sébastien Villemot <sebastien@dynare.org>.
 */

/*
 * Copyright (C) 2010-2012 Dynare Team
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
#include <octave/f77-fcn.h>

extern "C"
{
  F77_RET_T
  F77_FUNC(dtrsen, DTRSEN) (F77_CONST_CHAR_ARG_DECL, F77_CONST_CHAR_ARG_DECL,
                            const octave_idx_type *, const octave_idx_type &,
                            double *, const octave_idx_type &, double *, const octave_idx_type &,
                            double *, double *, octave_idx_type &, double &, double &, double *,
                            const octave_idx_type &, octave_idx_type *,
                            const octave_idx_type &, octave_idx_type &);
}

DEFUN_DLD(ordschur, args, nargout, "-*- texinfo -*-\n\
@deftypefn {Loadable Function} [ @var{us}, @var{ts} ] = ordschur (@var{u}, @var{t}, @var{select})\n\
\n\
Reorders the real Schur factorization @math{X = U*T*U'} so that selected\n\
eigenvalues appear in the upper left diagonal blocks of the quasi triangular\n\
Schur matrix @math{T}. The logical vector @var{select} specifies the selected\n\
eigenvalues as they appear along @math{T}'s diagonal.\n\
@end deftypefn\n\
")
{
  int nargin = args.length();
  octave_value_list retval;

  if (nargin != 3 || nargout != 2)
    {
      print_usage();
      return retval;
    }

  Matrix U = args(0).matrix_value();
  Matrix T = args(1).matrix_value();
  boolNDArray S = args(2).bool_array_value();
  if (error_state)
    return retval;

  dim_vector dimU = U.dims();
  dim_vector dimT = T.dims();
  octave_idx_type n = dimU(0);
  if (n != dimU(1) || n != dimT(0) || n != dimT(1))
    {
      error("ordschur: input matrices must be square and of same size");
      return retval;
    }
  if (S.nelem() != n)
    {
      error("ordschur: selection vector has wrong size");
      return retval;
    }

  octave_idx_type lwork = n, liwork = n;
  OCTAVE_LOCAL_BUFFER(double, wr, n);
  OCTAVE_LOCAL_BUFFER(double, wi, n);
  OCTAVE_LOCAL_BUFFER(double, work, lwork);
  OCTAVE_LOCAL_BUFFER(octave_idx_type, iwork, liwork);
  octave_idx_type m, info;
  double cond1, cond2;
  OCTAVE_LOCAL_BUFFER(octave_idx_type, S2, n);
  for (int i = 0; i < n; i++)
    S2[i] = S(i);

  F77_XFCN(dtrsen, dtrsen, (F77_CONST_CHAR_ARG("N"), F77_CONST_CHAR_ARG("V"),
                            S2, n, T.fortran_vec(), n, U.fortran_vec(), n,
                            wr, wi, m, cond1, cond2, work, lwork,
                            iwork, liwork, info));

  if (info != 0)
    {
      error("ordschur: dtrsen failed");
      return retval;
    }

  retval(0) = octave_value(U);
  retval(1) = octave_value(T);
  return retval;
}

/*

  %!test
  %! A = [1 2 3 -2; 4 5 6 -5 ; 7 8 9 -5; 10 11 12 4 ];
  %! [U, T] = schur(A);
  %! [US, TS] = ordschur(U, T, [ 0 0 1 1 ]);
  %! assert(US*TS*US', A, sqrt(eps))

*/
