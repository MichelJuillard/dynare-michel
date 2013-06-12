dnl Detect GSL.
dnl We don't use the official M4 macro since it relies on the script gsl-config,
dnl which does not work when cross-compiling.
dnl
dnl Copyright (C) 2010-2012 Dynare Team
dnl
dnl This file is part of Dynare.
dnl
dnl Dynare is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl Dynare is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

AC_DEFUN([AX_GSL],
[
AC_ARG_WITH(gsl, AC_HELP_STRING([--with-gsl=DIR], [prefix to GSL installation]),
            gsl_prefix="$withval", gsl_prefix="")

  has_gsl=yes

  if test "x$gsl_prefix" != "x"; then
    GSL_CPPFLAGS="-I$withval/include"
    GSL_LDFLAGS="-L$withval/lib"
  else
    GSL_CPPFLAGS=""
    GSL_LDFLAGS=""
  fi

  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_LIBS="$LIBS"

  LIBS=""
  CPPFLAGS="$CPPFLAGS $GSL_CPPFLAGS"
  LDFLAGS="$LDFLAGS $GSL_LDFLAGS"

	AC_LANG_PUSH(C)
  AC_CHECK_HEADER([gsl/gsl_cdf.h], [], [has_gsl=no])
	AC_LANG_POP(C)

  AC_CHECK_LIB([m], [cos])
  AC_CHECK_LIB([gslcblas], [cblas_dgemm], [LIBS="-lgslcblas $LIBS"], [has_gsl=no])
  AC_CHECK_LIB([gsl], [gsl_cdf_ugaussian_P], [LIBS="-lgsl $LIBS"], [has_gsl=no])

  if test "x$has_gsl" = "xyes"; then
    GSL_LIBS="$LIBS"
  else
    GSL_LIBS=""
  fi

  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"

  AC_SUBST(GSL_CPPFLAGS)
  AC_SUBST(GSL_LDFLAGS)
  AC_SUBST(GSL_LIBS)
])
