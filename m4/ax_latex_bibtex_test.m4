# ===========================================================================
#          http://www.nongnu.org/autoconf-archive/ax_latex_test.html
# ===========================================================================
#
# OBSOLETE MACRO
#
#   Deprecated because of licensing issues. The Lesser GPL imposes licensing
#   restrictions on the generated configure script unless it is augmented
#   with an Autoconf Exception clause.
#
# SYNOPSIS
#
#   AX_LATEX_BIBTEX_TEST(FILEDATA,BIBDATA,VARIABLETOSET,[NOCLEAN])
#
# DESCRIPTION
#
#   This macros creates a bib file called contest.bib with BIBDATA,
#   executes the latex application with FILEDATA as input, then runs
#   bibtex on the resulting aux file, and finally sets VARIABLETOSET
#   to yes or no depending on the result. If NOCLEAN is set, the folder
#   used for the test is not deleted after testing.
#
#   The macro assumes that the variables PDFLATEX and BIBTEX are set.
#
#   Adapted from the macro AX_LATEX_TEST by Sébastien Villemot.
#
# LICENSE
#
#   Copyright (c) 2008 Boretti Mathieu <boretti@eig.unige.ch>
#   Copyright (c) 2009 Dynare Team
#
#   This library is free software; you can redistribute it and/or modify it
#   under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation; either version 2.1 of the License, or (at
#   your option) any later version.
#
#   This library is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
#   General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with this library. If not, see <http://www.gnu.org/licenses/>.

AC_DEFUN([AX_LATEX_BIBTEX_TEST],[
rm -rf conftest.dir/.acltx
AS_MKDIR_P([conftest.dir/.acltx])
cd conftest.dir/.acltx
m4_ifval([$3],[$3="no"; export $3;])
cat > conftest.tex << ACLEOF
$1
ACLEOF
cat > conftest.bib << ACLEOF
$2
ACLEOF
$PDFLATEX conftest 2>&1 1>output
$BIBTEX conftest 2>&1 1>output2 m4_ifval([$3],[&& $3=yes])
cd ..
cd ..
sed 's/^/| /' conftest.dir/.acltx/conftest.tex >&5
echo "$as_me:$LINENO: executing $PDFLATEX conftest" >&5
sed 's/^/| /' conftest.dir/.acltx/output >&5
echo "$as_me:$LINENO: executing $BIBTEX conftest" >&5
sed 's/^/| /' conftest.dir/.acltx/output2 >&5
m4_ifval([$4],,[rm -rf conftest.dir/.acltx])
])
