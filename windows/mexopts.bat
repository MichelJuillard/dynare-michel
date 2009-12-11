@echo off
rem
rem    Compile and link options used for building MEX-files
rem    using Cygwin.
rem
rem    This file makes the assumption that you installed Cygwin in
rem    C:\CYGWIN, and that you install gcc-mingw package.
rem
rem    This file should be copied to:
rem    C:\Documents and Settings\<Username>\Application Data\MathWorks\MATLAB\<MATLAB version>\
rem
rem    Initial version by Michel Juillard, revised by Sebastien Villemot
rem

rem ********************************************************************
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%
set PATH=%PATH%;c:\cygwin\bin

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=gcc-3
set COMPFLAGS=-c -mno-cygwin -I"%MATLAB%\extern\include"
set OPTIMFLAGS=-O3
set DEBUGFLAGS=-g -Wall
set NAME_OBJECT=-o
set MW_TARGET_ARCH=win32

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC="%MATLAB%"\extern\lib\win32\microsoft\\
set PRELINK_CMDS1=echo EXPORTS > mex.def & echo mexFunction >> mex.def
set LINKER=gcc-3
set LINKFLAGS= -mno-cygwin -shared mex.def
set LINKFLAGSPOST= %LIBLOC%libmex.lib %LIBLOC%libmx.lib %LIBLOC%libmwlapack.lib %LIBLOC%libmwblas.lib -lstdc++
set LINKOPTIMFLAGS=-O3
set LINKDEBUGFLAGS= -g -Wall
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=-o "%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@
set POSTLINK_CMDS1=del mex.def
