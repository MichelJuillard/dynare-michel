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

AC_PROG_CXX

# Use C++ for testing headers
AC_LANG([C++])

CXXFLAGS="$CXXFLAGS -Wall"

AC_PROG_RANLIB
AC_PROG_LN_S
AC_PROG_MKDIR_P

AM_PROG_LEX
AX_PROG_FLEX([], AC_MSG_ERROR([No suitable flex found.]))

AC_CHECK_HEADERS([FlexLexer.h], [], [AC_MSG_ERROR([Can't find FlexLexer.h])])

AC_PROG_YACC
AX_PROG_BISON([], AC_MSG_ERROR([No suitable bison found.]))

AX_BOOST_BASE([1.34])

CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp], [], [AC_MSG_ERROR([Can't find Boost Graph Library])])
CPPFLAGS="$CPPFLAGS_SAVED"

# Don't use deprecated hash structures
AC_DEFINE([BOOST_NO_HASH])

AC_CHECK_PROG([PDFTEX], [pdftex], [pdftex])
AM_CONDITIONAL([HAVE_PDFTEX], [test "x$PDFTEX" != "x"])

AC_CHECK_PROG([PDFLATEX], [pdflatex], [pdflatex])
AM_CONDITIONAL([HAVE_PDFLATEX], [test "x$PDFLATEX" != "x"])

AC_CHECK_PROG([BIBTEX], [bibtex], [bibtex])
AM_CONDITIONAL([HAVE_BIBTEX], [test "x$BIBTEX" != "x"])

AC_CHECK_PROG([DBLATEX], [dblatex], [dblatex])
AM_CONDITIONAL([HAVE_DBLATEX], [test "x$DBLATEX" != "x"])

AX_PROG_XSLTPROC
AX_CHECK_DOCBOOK_XSLT_MIN(1.69.0)
AM_CONDITIONAL([HAVE_XSLT], [test "x$XSLTPROC" != "x" -a "x$DOCBOOK_XSLT_CURRENT_VERSION" != "xno"])

AC_CHECK_PROG([DOXYGEN], [doxygen], [doxygen])
AM_CONDITIONAL([HAVE_DOXYGEN], [test "x$DOXYGEN" != "x"])

AC_CHECK_PROG([OCTAVE], [octave], [octave])
AM_CONDITIONAL([HAVE_OCTAVE], [test "x$OCTAVE" != "x"])

AC_CHECK_PROG([MKOCTFILE], [mkoctfile], [mkoctfile])
AM_CONDITIONAL([HAVE_MKOCTFILE], [test "x$MKOCTFILE" != "x"])

AC_CHECK_PROG([CTANGLE], [ctangle], [ctangle])
AM_CONDITIONAL([HAVE_CTANGLE], [test "x$CTANGLE" != "x"])

AC_CHECK_PROG([CWEAVE], [cweave], [cweave])
AM_CONDITIONAL([HAVE_CWEAVE], [test "x$CWEAVE" != "x"])

AX_BLAS
AM_CONDITIONAL([HAVE_BLAS], [test x"$ax_blas_ok" = "xyes"])

AX_LAPACK
AM_CONDITIONAL([HAVE_LAPACK], [test x"$ax_lapack_ok" = "xyes"])

AX_PTHREAD
AM_CONDITIONAL([HAVE_PTHREAD], [test x"$ax_pthread_ok" = "xyes"])

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
fi
AM_CONDITIONAL([HAVE_MEXOPTS], [test "x$ax_enable_matlab" = "xyes" -a "x$ax_matlab_version_ok" = "xyes" -a "x$ax_mexopts_ok" = "xyes" ])


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
                 dynare++/Makefile
                 dynare++/parser/cc/Makefile
                 dynare++/sylv/Makefile
                 dynare++/sylv/cc/Makefile
                 dynare++/sylv/testing/Makefile
                 dynare++/tl/Makefile
                 dynare++/tl/cc/Makefile
                 dynare++/tl/testing/Makefile
                 dynare++/doc/Makefile
                 dynare++/utils/cc/Makefile
                 dynare++/integ/Makefile
                 dynare++/integ/cc/Makefile
                 dynare++/integ/src/Makefile
                 dynare++/integ/testing/Makefile
                 dynare++/kord/Makefile
                 dynare++/src/Makefile
                 mex/build/Makefile
                 mex/build/mjdgges/octave/Makefile
                 mex/build/mjdgges/matlab/Makefile
                 mex/build/kronecker/octave/Makefile
                 mex/build/kronecker/matlab/Makefile
                 mex/build/gensylv/octave/Makefile
                 mex/build/gensylv/matlab/Makefile
                 mex/build/bytecode/octave/Makefile
                 mex/build/bytecode/matlab/Makefile
])

# Construct final output message

BUILD_PREPROCESSOR="yes"

if test "x$CTANGLE" != "x" -a x"$ax_blas_ok" = "xyes" -a x"$ax_lapack_ok" = "xyes" -a x"$ax_pthread_ok" = "xyes"; then
  BUILD_DYNAREPLUSPLUS="yes"
else
  BUILD_DYNAREPLUSPLUS="no (missing one of: ctangle, BLAS, LAPACK, POSIX Threads)"
fi

if test "x$CWEAVE" != "x" -a x"$PDFTEX" != "x"; then
  BUILD_DYNAREPLUSPLUS_SRCDOC="yes"
else
  BUILD_DYNAREPLUSPLUS_SRCDOC="no (missing one of: cweave, pdftex)"
fi

if test "x$XSLTPROC" != "x" -a "x$DOCBOOK_XSLT_CURRENT_VERSION" != "xno"; then
  BUILD_DYNARE_HTML_MANUAL="yes"
else
  BUILD_DYNARE_HTML_MANUAL="no (missing one of: xsltproc, DocBook stylesheets)"
fi

if test "x$PDFLATEX" != "x" -a "x$DBLATEX" != "x"; then
  BUILD_DYNARE_PDF_MANUAL="yes"
else
  BUILD_DYNARE_PDF_MANUAL="no (missing one of: pdflatex, dblatex)"
fi

if test "x$PDFLATEX" != "x" -a "x$BIBTEX" != "x"; then
  BUILD_DYNARE_USERGUIDE="yes"
else
  BUILD_DYNARE_USERGUIDE="no (missing one of: pdflatex, bibtex)"
fi

if test "x$PDFLATEX" != "x"; then
  BUILD_OTHER_PDF_DOC="yes"
else
  BUILD_OTHER_PDF_DOC="no (missing pdflatex)"
fi

if test "x$DOXYGEN" != "x"; then
  BUILD_DYNARE_PREPROC_DOC="yes"
else
  BUILD_DYNARE_PREPROC_DOC="no (missing doxygen)"
fi

if test "x$MKOCTFILE" != "x"; then
  BUILD_MEX_OCTAVE="yes"
else
  BUILD_MEX_OCTAVE="no (missing mkoctfile)"
fi

if test "x$ax_enable_matlab" = "xyes" -a "x$ax_matlab_version_ok" = "xyes" -a "x$ax_mexopts_ok" = "xyes"; then
  BUILD_MEX_MATLAB="yes"
else
  BUILD_MEX_MATLAB="no (missing MATLAB, or unknown version, or unknown architecture)"
fi

if test "x$OCTAVE" != "x"; then
  TESTSUITE_OCTAVE="yes"
else
  TESTSUITE_OCTAVE="no"
fi

if test "$BUILD_DYNAREPLUSPLUS" = "yes"; then
  TESTSUITE_DYNAREPLUSPLUS="yes"
else
  TESTSUITE_DYNAREPLUSPLUS="no"
fi

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 Dynare preprocessor:                   $BUILD_PREPROCESSOR
 MEX files for Octave:                  $BUILD_MEX_OCTAVE
 MEX files for MATLAB:                  $BUILD_MEX_MATLAB
 Dynare++:                              $BUILD_DYNAREPLUSPLUS

PDF documentation (with "make pdf"):
 Dynare reference manual:               $BUILD_DYNARE_PDF_MANUAL
 Dynare user guide:                     $BUILD_DYNARE_USERGUIDE
 Dynare++ developer documentation:      $BUILD_DYNAREPLUSPLUS_SRCDOC
 Various other documents:               $BUILD_OTHER_PDF_DOC

HTML documentation (with "make html"):
 Dynare reference manual:               $BUILD_DYNARE_HTML_MANUAL
 Dynare preprocessor developer doc:     $BUILD_DYNARE_PREPROC_DOC

Testsuites (run with "make check"):
 Dynare for Octave:                     $TESTSUITE_OCTAVE
 Dynare++:                              $TESTSUITE_DYNAREPLUSPLUS
])

AC_OUTPUT
