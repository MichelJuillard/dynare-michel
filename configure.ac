# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([dynare], [4.1-Prerelease])
AC_CONFIG_SRCDIR([preprocessor/DynareMain.cc])

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

AC_CHECK_PROGS([FLEX], [flex])
if test -z "$FLEX"; then
  AC_MSG_ERROR([No suitable flex found.])
fi

AC_CHECK_HEADERS([FlexLexer.h], [], [AC_MSG_ERROR([Can't find FlexLexer.h])])

AC_CHECK_PROGS([BISON], [bison])
if test -z "$BISON"; then
  AC_MSG_ERROR([No suitable bison found.])
fi

AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp], [], [AC_MSG_ERROR([Can't find Boost Graph Library])])

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

AC_CONFIG_FILES([preprocessor/Makefile doc/Makefile doc/preprocessor/Makefile doc/macroprocessor/Makefile doc/userguide/Makefile])

AC_OUTPUT
