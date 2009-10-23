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
AC_INIT([dynare], [4.1-unstable])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AX_MATLAB

if test "x$ax_enable_matlab" = "xyes"; then
  AX_DOT_MEXEXT
  AX_MATLAB_ARCH
  AX_MATLAB_VERSION

  if test "x$ax_matlab_version_ok" = "xyes"; then
    AX_MEXOPTS

    # Determine subdirectory of mex/ to be used for storing MATLAB MEX binaries
    AX_COMPARE_VERSION([$MATLAB_VERSION], [le], [7.4], [DYNARE_MEXDIR="2007a"], [
      DYNARE_MEXDIR="2007b"
      if test ${MATLAB_ARCH%64} != ${MATLAB_ARCH}; then
         AX_COMPARE_VERSION([$MATLAB_VERSION], [ge], [7.8], [DYNARE_MEXDIR="2009a-64bit"])
      fi])
    AC_SUBST([DYNARE_MEXDIR])
  fi

  CFLAGS="$MATLAB_CFLAGS"
  CXXFLAGS="$MATLAB_CXXFLAGS"
fi

AC_PROG_CC
AC_PROG_CXX
AX_PROG_LN_S
AX_PTHREAD
# Check for dlopen(), needed by k_order_perturbation DLL
AC_CHECK_LIB([dl], [dlopen], [LIBADD_DLOPEN="-ldl"], [])
AC_SUBST([LIBADD_DLOPEN])

AM_CONDITIONAL([DO_SOMETHING], [test "x$ax_enable_matlab" = "xyes" -a "x$ax_matlab_version_ok" = "xyes" -a "x$ax_mexopts_ok" = "xyes"])

AC_CONFIG_FILES([Makefile
                 mjdgges/Makefile
                 kronecker/Makefile
                 gensylv/Makefile
                 bytecode/Makefile
                 k_order_perturbation/Makefile])

if test "x$ax_enable_matlab" = "xyes" -a "x$ax_matlab_version_ok" = "xyes" -a "x$ax_mexopts_ok" = "xyes"; then
  if test x"$ax_pthread_ok" = "xyes"; then
     BUILD_MEX_MATLAB="yes"
  else
     BUILD_MEX_MATLAB="yes (without POSIX threads)"
  fi
else
  BUILD_MEX_MATLAB="no (missing MATLAB, or unknown version, or unknown architecture)"
fi

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 MEX files for MATLAB:                  $BUILD_MEX_MATLAB

])

AC_OUTPUT
