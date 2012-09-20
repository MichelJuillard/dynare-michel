dnl Detect the MATIO Library.
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

AC_DEFUN([AX_MATIO],
[
AC_ARG_WITH(matio, AC_HELP_STRING([--with-matio=DIR], [prefix to MATIO installation]),
            matio_prefix="$withval", matio_prefix="")

  has_matio=yes

  if test "x$matio_prefix" != "x"; then
    CPPFLAGS_MATIO="-I$withval/include"
    LDFLAGS_MATIO="-L$withval/lib"
  else
    CPPFLAGS_MATIO=""
    LDFLAGS_MATIO=""
  fi

  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LDFLAGS="$LDFLAGS"
  ac_save_LIBS="$LIBS"

  LIBADD_MATIO=""
  CPPFLAGS="$CPPFLAGS_MATIO $CPPFLAGS"
  LDFLAGS="$LDFLAGS_MATIO $LDFLAGS"

  AC_CHECK_HEADER([matio.h], [], [has_matio=no])
  AC_CHECK_LIB([matio], [Mat_Open], [LIBADD_MATIO="-lmatio"], [has_matio=no])

  CPPFLAGS="$ac_save_CPPFLAGS"
  LDFLAGS="$ac_save_LDFLAGS"
  LIBS="$ac_save_LIBS"

  AC_SUBST(CPPFLAGS_MATIO)
  AC_SUBST(LIBADD_MATIO)
  AC_SUBST(LDFLAGS_MATIO)
])
