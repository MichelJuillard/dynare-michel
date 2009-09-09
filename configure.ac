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
AC_CONFIG_SRCDIR([preprocessor/DynareMain.cc])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CC
AC_PROG_CXX

# Use C++ for testing headers
AC_LANG([C++])

CXXFLAGS="$CXXFLAGS -Wall"

AC_CANONICAL_HOST

case $host_os in
  cygwin*)
    # Create standalone binaries on Cygwin and eliminate spurious warnings
    CXXFLAGS="$CXXFLAGS -mno-cygwin -Wno-uninitialized"
    ;;
  mingw32*)
    # Suppress spurious warnings with old versions of MinGW
    CXXFLAGS="$CXXFLAGS -Wno-uninitialized"
    ;;
esac

AC_PROG_RANLIB

AM_PROG_LEX
AX_PROG_FLEX([], AC_MSG_ERROR([No suitable flex found.]))

AC_CHECK_HEADERS([FlexLexer.h], [], [AC_MSG_ERROR([Can't find FlexLexer.h])])

AC_PROG_YACC
AX_PROG_BISON([], AC_MSG_ERROR([No suitable bison found.]))

AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp], [], [AC_MSG_ERROR([Can't find Boost Graph Library])])
AC_CHECK_HEADERS([boost/graph/max_cardinality_matching.hpp], [], [AC_MSG_ERROR([Your Boost version is too old. Please install a newer one.])])

# Don't use deprecated hash structures
AC_DEFINE([BOOST_NO_HASH])

AC_CHECK_PROG([PDFLATEX], [pdflatex], [pdflatex])
AC_CHECK_PROG([BIBTEX], [bibtex], [bibtex])
AC_CHECK_PROG([DBLATEX], [dblatex], [dblatex])
AC_CHECK_PROG([XSLTPROC], [xsltproc], [xsltproc])

if test -n "${XSLTPROC}" -a -f /etc/xml/catalog
then
  AC_PROG_GREP
  AC_MSG_CHECKING([for DocBook XSL files in XML catalog])
  cat /etc/xml/catalog | ${GREP} -q http://docbook.sourceforge.net/release/xsl/
  if test $? -eq 0
  then
      AC_MSG_RESULT([yes])
  else
      AC_MSG_RESULT([no])
      AC_MSG_WARN([The creation of the HTML version will be very slow! You should install the docbook-xsl package.])
  fi
fi

AC_CHECK_PROG([OCTAVE], [octave], [octave])

AX_MATLAB

if test "$ax_enable_matlab" = "yes"
then
  AX_PATH_MEX
  AX_MEXEXT
  AX_MATLAB_ARCH
fi

AC_CONFIG_FILES([Makefile
                 preprocessor/macro/Makefile
                 preprocessor/Makefile
                 doc/Makefile
                 doc/preprocessor/Makefile
                 doc/macroprocessor/Makefile
                 doc/userguide/Makefile
                 tests/Makefile
                 matlab/dynare_version.m
                 windows/dynare-version.nsi
])

AC_OUTPUT
