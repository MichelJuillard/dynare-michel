dnl dynarepp.m4 --- check for Dynare++
dnl
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
dnl
dnl Code:

# AX_DYNAREPP
# ---------
# Check for Dynare++.
AC_DEFUN([AX_DYNAREPP],
[dnl
AC_PREREQ([2.63])
ax_enable_dynarepp=
AC_ARG_WITH([dynarepp], AC_HELP_STRING([--with-dynarepp=ARG], [path to Dynare++]),
[case $withval in
  no)
    # Explicitly enable or disable Dynare++ but determine
    # Matlab prefix automatically.
    ax_enable_dynarepp=no
    ;;
  *)
    # Enable Dynare++ and use ARG as the Dynare++ prefix.
    # ARG must be an existing directory.
    ax_enable_dynarepp=yes
    DYNAREPP=`cd "${withval-/}" > /dev/null 2>&1 && pwd`
    if test -z "$DYNAREPP" ; then
	AC_MSG_ERROR([invalid value '$withval' for --with-dynarepp])
    fi
    ;;
esac])
AC_MSG_CHECKING([for Dynare++])
if test x$ax_enable_dynarepp != xno ; then
    if test "${DYNAREPP+set}" = set && test -d "$DYNAREPP/kord" ; then
	ax_enable_dynarepp=yes
    elif test x$ax_enable_dynarepp = x ; then
	ax_enable_dynarepp=no
    else
	# Fail if Matlab was explicitly enabled.
	AC_MSG_RESULT([failure])
	AC_MSG_ERROR([check your Dynare++ setup])
    fi
fi
AC_MSG_RESULT([$ax_enable_dynarepp])
AC_SUBST([DYNAREPP])
])

dnl dynarepp.m4 ends here

dnl Local variables:
dnl tab-width: 8
dnl End:
