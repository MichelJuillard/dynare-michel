% Build file for Dynare MEX Librairies under Matlab

% Copyright (C) 2007-2008 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

addpath '../../matlab'; % For matlab_ver_less_than

MATLAB_PATH = matlabroot;

COMPILE_OPTIONS = '';

if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer) ...
        || strcmpi('MACI', computer) || strcmpi('MAC', computer)
    % GNU/Linux (x86-32 or x86-64) or MacOS (Intel or PPC)
    LAPACK_PATH = '-lmwlapack';
    if matlab_ver_less_than('7.5')
        BLAS_PATH = LAPACK_PATH; % On <= 7.4, BLAS in included in LAPACK
    else
        BLAS_PATH = '-lmwblas';
    end
elseif strcmpi('PCWIN', computer) || strcmpi('PCWIN64', computer)
    % Windows (x86-32 or x86-64) with Microsoft or gcc compiler
    if strcmpi('PCWIN', computer)
      LIBRARY_PATH = [MATLAB_PATH '/extern/lib/win32/microsoft/'];
    else
      LIBRARY_PATH = [MATLAB_PATH '/extern/lib/win64/microsoft/'];
    end
    LAPACK_PATH = ['"' LIBRARY_PATH 'libmwlapack.lib"'];
    if matlab_ver_less_than('7.5')
        BLAS_PATH = LAPACK_PATH; % On <= 7.4, BLAS in included in LAPACK
    else
        BLAS_PATH = ['"' LIBRARY_PATH 'libmwblas.lib"'];
    end
else
    error('Unsupported platform')
end

% mwSize and mwIndex appeared in Matlab 7.3
if matlab_ver_less_than('7.3')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMWTYPES_NOT_DEFINED' ];
end

% Large array dims for 64 bits platforms appeared in Matlab 7.3
if (strcmpi('GLNXA64', computer) || strcmpi('PCWIN64', computer)) ...
      && ~matlab_ver_less_than('7.3')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
end

% blas.h and lapack.h appeared in Matlab 7.5
if matlab_ver_less_than('7.5')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DNO_BLAS_H -DNO_LAPACK_H' ];
end

if matlab_ver_less_than('7.5')
    OUTPUT_DIR = '../2007a';
else
    OUTPUT_DIR = '../2007b';
end

% Comment next line to suppress compilation debugging info
% COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -v' ];

COMPILE_COMMAND = [ 'mex ' COMPILE_OPTIONS ' -outdir ' OUTPUT_DIR ];

CFLAGS   = ' CFLAGS="\$CFLAGS -fopenmp" ';
CXXFLAGS = ' CXXFLAGS="\$CXXFLAGS -fopenmp" ';
LDFLAGS  = ' LDFLAGS="\$LDFLAGS -fopenmp" ';

disp('Compiling isopenmp...')
try
    eval([ 'mex ' COMPILE_OPTIONS  CFLAGS CXXFLAGS LDFLAGS  ' -outdir ' OUTPUT_DIR ' threads/isopenmp.cc ' ]);
    disp(' ')
    disp('|------------------------------------------------|')
    disp('|  OpenMp is used (multithreaded mex files) for: |')
    disp('|   * sparse_hessian_times_B_kronecker_C.cc      |')
    disp('|   * A_times_B_kronecker_C.cc                   |')
    disp('|   * simulate (SparseMatrix.cc)                 |')
    disp('|------------------------------------------------|')
    disp(' ')
catch
    eval([ 'mex ' COMPILE_OPTIONS ' -outdir ' OUTPUT_DIR ' threads/isopenmp.cc ' ]);
    disp(' ')
    disp('|------------------------------------------------|')
    disp('|  OpenMp is not available on this platform!     |')
    disp('|------------------------------------------------|')
    disp(' ')
    CFLAGS = [];
    CXXFLAGS = [];
    LDFLAGS = [];
end

COMPILE_COMMAND_OMP = [ 'mex ' COMPILE_OPTIONS  CFLAGS CXXFLAGS LDFLAGS  ' -outdir ' OUTPUT_DIR ];

system(['cp ' OUTPUT_DIR '/isopenmp.' mexext ' ./isopenmp.' mexext]);

disp('Compiling mjdgges...')
eval([ COMPILE_COMMAND ' mjdgges/mjdgges.c ' LAPACK_PATH ]);

disp('Compiling sparse_hessian_times_B_kronecker_C...')
eval([ COMPILE_COMMAND_OMP ' kronecker/sparse_hessian_times_B_kronecker_C.cc' ]);

disp('Compiling A_times_B_kronecker_C...')
if isopenmp
    eval([ COMPILE_COMMAND_OMP ' kronecker/A_times_B_kronecker_C.cc ']);    
else
    eval([ COMPILE_COMMAND ' kronecker/A_times_B_kronecker_C.cc ' BLAS_PATH]);
end

disp('Compiling gensylv...')
eval([ COMPILE_COMMAND ' -DMATLAB -Igensylv/cc ' ...
       'gensylv/matlab/gensylv.cpp ' ...
       'gensylv/cc/BlockDiagonal.cpp ' ... 
       'gensylv/cc/GeneralMatrix.cpp ' ...
       'gensylv/cc/GeneralSylvester.cpp ' ...
       'gensylv/cc/IterativeSylvester.cpp ' ...
       'gensylv/cc/KronUtils.cpp ' ...
       'gensylv/cc/KronVector.cpp ' ...
       'gensylv/cc/QuasiTriangular.cpp ' ...
       'gensylv/cc/QuasiTriangularZero.cpp ' ...
       'gensylv/cc/SchurDecomp.cpp ' ...
       'gensylv/cc/SchurDecompEig.cpp ' ...
       'gensylv/cc/SimilarityDecomp.cpp ' ...
       'gensylv/cc/SylvException.cpp ' ...
       'gensylv/cc/SylvMatrix.cpp ' ...
       'gensylv/cc/SylvMemory.cpp ' ...
       'gensylv/cc/SylvParams.cpp ' ...
       'gensylv/cc/TriangularSylvester.cpp ' ...
       'gensylv/cc/Vector.cpp ' ...
       BLAS_PATH ' ' LAPACK_PATH ]);

disp('Compiling simulate...')
eval([ COMPILE_COMMAND_OMP ' -Isimulate -I../../preprocessor simulate/simulate.cc simulate/Interpreter.cc simulate/Mem_Mngr.cc simulate/SparseMatrix.cc']);

% Clean up
system(['rm ./isopenmp.' mexext]);