dnl Detect the SLICOT Library.
dnl Called with an argument of either 'matlab' or 'octave', depending
dnl on the configure script from which we're calling it
dnl
dnl AX_SLICOT([matlab])
dnl AX_SLICOT([octave])
dnl
dnl Copyright (C) 2012 Dynare Team
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

AC_DEFUN([AX_SLICOT],
[
  if test "x$1" != "xmatlab" && test "x$1" != "xoctave"; then
    AC_MSG_ERROR([Argument to autoconf slicot macro must be either 'matlab' or 'octave'])
  fi

  AC_ARG_WITH(slicot, AC_HELP_STRING([--with-slicot=DIR], [prefix to SLICOT installation]),
              slicot_prefix="$withval", slicot_prefix="")
  has_slicot=yes

  if test "x$slicot_prefix" != "x"; then
    LDFLAGS_SLICOT="-L$withval/lib"
  else
    LDFLAGS_SLICOT=""
  fi
  ac_save_LDFLAGS="$LDFLAGS"
  LDFLAGS_SAVED="$LDFLAGS"

  AC_F77_FUNC(sb02od)

  if test "x$1" = "xmatlab"; then
    LDFLAGS="$MATLAB_LDFLAGS_NOMAP $LDFLAGS_SLICOT"

    case ${MATLAB_ARCH} in
       glnxa64 | win64 | maci64)
         AX_COMPARE_VERSION([$MATLAB_VERSION], [ge], [7.8], [use_64_bit_indexing=yes], [use_64_bit_indexing=no])
         ;;
       *)
         use_64_bit_indexing=no
         ;;
    esac

    if test "$use_64_bit_indexing" = "yes"; then
       AC_CHECK_LIB([slicot64_pic], [$sb02od], [LIBADD_SLICOT="-lslicot64_pic"], [has_slicot=no], [$MATLAB_LIBS])
    else
       AC_CHECK_LIB([slicot_pic], [$sb02od], [LIBADD_SLICOT="-lslicot_pic"], [has_slicot=no], [$MATLAB_LIBS])
    fi
  else
    LDFLAGS="$LDFLAGS $LDFLAGS_SLICOT"
    AC_CHECK_LIB([slicot], [$sb02od], [LIBADD_SLICOT="-lslicot"],
             [
               AC_CHECK_LIB([slicot_pic], [$sb02od], [LIBADD_SLICOT="-lslicot_pic"], [has_slicot=no], [`$MKOCTFILE -p BLAS_LIBS` `$MKOCTFILE -p LAPACK_LIBS`])
             ], # Fallback on libslicot_pic if dynamic libslicot not found
             [`$MKOCTFILE -p BLAS_LIBS` `$MKOCTFILE -p LAPACK_LIBS`])
  fi

  LDFLAGS="$ac_save_LDFLAGS"
  AC_SUBST(LDFLAGS_SLICOT)
  AC_SUBST(LIBADD_SLICOT)
])
