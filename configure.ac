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

AC_CONFIG_FILES([preprocessor/Makefile])

AC_OUTPUT
