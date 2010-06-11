dnl Process this file with autoconf to produce a configure script.

dnl Copyright (C) 2009 Dynare Team
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

AC_PREREQ([2.61])
AC_INIT([dynare], [4.2-unstable])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CHECK_PROG([MKOCTFILE], [mkoctfile], [mkoctfile])

if test "x$MKOCTFILE" != "x"; then
  CC=`$MKOCTFILE -p CC`
  CXX=`$MKOCTFILE -p CXX`
  CFLAGS=`$MKOCTFILE -p CFLAGS`
  CXXFLAGS=`$MKOCTFILE -p CXXFLAGS`
fi

CFLAGS="$CFLAGS -Wall"
CXXFLAGS="$CXXFLAGS -Wall"

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AX_PROG_LN_S
AC_PROG_MKDIR_P
AX_PTHREAD
# Check for dlopen(), needed by k_order_perturbation DLL
AC_CHECK_LIB([dl], [dlopen], [LIBADD_DLOPEN="-ldl"], [])
AC_SUBST([LIBADD_DLOPEN])
# Check for GSL, needed by SWZ
AX_PATH_GSL
AM_CONDITIONAL([HAVE_GSL], [test "x$no_gsl" = "x"])

AM_CONDITIONAL([DO_SOMETHING], [test "x$MKOCTFILE" != "x"])

if test "x$MKOCTFILE" != "x"; then
  if test x"$ax_pthread_ok" = "xyes"; then
     BUILD_MEX_OCTAVE="yes"
  else
     BUILD_MEX_OCTAVE="yes (without POSIX threads)"
  fi
else
  BUILD_MEX_OCTAVE="no (missing mkoctfile)"
fi

if test "x$MKOCTFILE" != "x" -a "x$no_gsl" = "x"; then
   BUILD_SWZ_MEX_OCTAVE="yes"
else
   BUILD_SWZ_MEX_OCTAVE="no (missing GSL)"
fi


AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 MEX files for Octave (except SWZ):     $BUILD_MEX_OCTAVE
 SWZ MEX files for Octave:              $BUILD_SWZ_MEX_OCTAVE

])

AC_CONFIG_FILES([Makefile
                 mjdgges/Makefile
                 kronecker/Makefile
                 bytecode/Makefile
                 libdynare++/Makefile
                 gensylv/Makefile
                 k_order_perturbation/Makefile
                 dynare_simul_/Makefile
                 swz/Makefile])

AC_OUTPUT
