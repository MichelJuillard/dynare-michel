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

AC_DEFUN([AX_MEXOPTS],
[dnl
AC_REQUIRE([AX_MEXEXT])
AC_REQUIRE([AX_MATLAB_ARCH])
AC_REQUIRE([AX_MATLAB_VERSION])

AC_MSG_CHECKING([for options to compile MEX for MATLAB])

MATLAB_CPPFLAGS="-I$MATLAB/extern/include"

case ${MATLAB_ARCH} in
  glnx86 | glnxa64)
    MATLAB_DEFS="-D_GNU_SOURCE -DNDEBUG"
    MATLAB_CC="gcc"
    MATLAB_CFLAGS="-ansi -fexceptions -fPIC -pthread -g -O2"
    MATLAB_CXX="g++"
    MATLAB_CXXFLAGS="-ansi -fPIC -pthread -g -O2"
    MATLAB_LDFLAGS="-shared -Wl,--version-script,$MATLAB/extern/lib/${MATLAB_ARCH}/mexFunction.map -Wl,--no-undefined -Wl,-rpath-link,$MATLAB/bin/${MATLAB_ARCH} -L$MATLAB/bin/${MATLAB_ARCH}"
    MATLAB_LIBS="-lmx -lmex -lmat -lm -lstdc++ -lmwlapack"
    # Starting from MATLAB 7.5, BLAS and LAPACK are in distinct libraries
    AX_COMPARE_VERSION([$MATLAB_VERSION], [ge], [7.5], [MATLAB_LIBS="${MATLAB_LIBS} -lmwblas"])
    ax_mexopts_ok="yes"
    ;;
  *)
    ax_mexopts_ok="no"
    ;;
esac

case ${MATLAB_ARCH} in
  glnx86)
    MATLAB_DEFS="$MATLAB_DEFS -D_FILE_OFFSET_BITS=64"
    MATLAB_CFLAGS="$MATLAB_CFLAGS -m32"
    MATLAB_CXXFLAGS="$MATLAB_CXXFLAGS -m32"
    ;;
  glnxa64)
    MATLAB_CFLAGS="$MATLAB_CFLAGS -fno-omit-frame-pointer"
    MATLAB_CXXFLAGS="$MATLAB_CXXFLAGS -fno-omit-frame-pointer"
    ;;
  *)
    ;;
esac

# mwSize, mwIndex and mwSignedIndex appeared in MATLAB 7.3
AX_COMPARE_VERSION([$MATLAB_VERSION], [lt], [7.3], [MATLAB_DEFS="$MATLAB_DEFS -DMWTYPES_NOT_DEFINED"])

# MATLAB Lapack expects mwSignedIndex arguments only starting with MATLAB 7.8
AX_COMPARE_VERSION([$MATLAB_VERSION], [ge], [7.8], [MATLAB_DEFS="$MATLAB_DEFS -DLAPACK_USE_MWSIGNEDINDEX"])

# blas.h and lapack.h appeared in MATLAB 7.5
AX_COMPARE_VERSION([$MATLAB_VERSION], [lt], [7.5], [MATLAB_DEFS="$MATLAB_DEFS -DNO_BLAS_H -DNO_LAPACK_H"])

if test "$ax_mexopts_ok" = "yes"; then
  AC_MSG_RESULT([ok])
else
  AC_MSG_RESULT([unknown])
fi

AC_SUBST([MATLAB_CPPFLAGS])
AC_SUBST([MATLAB_DEFS])
AC_SUBST([MATLAB_CC])
AC_SUBST([MATLAB_CFLAGS])
AC_SUBST([MATLAB_CXX])
AC_SUBST([MATLAB_CXXFLAGS])
AC_SUBST([MATLAB_LDFLAGS])
AC_SUBST([MATLAB_LIBS])
])
