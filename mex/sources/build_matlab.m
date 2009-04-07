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

% mwSize, mwIndex and mwSignedIndex appeared in Matlab 7.3
if matlab_ver_less_than('7.3')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMWTYPES_NOT_DEFINED' ];
end

% Matlab Lapack expects mwSignedIndex arguments only starting with Matlab 7.8
if ~matlab_ver_less_than('7.8')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DLAPACK_USE_MWSIGNEDINDEX' ];
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

disp(' ')
if exist(OUTPUT_DIR,'dir')
    disp('Delete old mex files.')
    delete([OUTPUT_DIR '/*.' mexext]);
else
    whereami = pwd;
    disp(['Create directory ' whereami(1:end-7) OUTPUT_DIR(4:end) '.'])
    mkdir(OUTPUT_DIR);
end
disp(' ')

COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DNO_OPENMP' ]; % Single thread version. To compile multithreaded versions of the mex file 
                                                      % use build_matlab_multithread.m instead.
                                                      % Read http://www.dynare.org/DynareWiki/UsingMultithreadedDlls. 
% Set Optimization and Debug flags
CXXOPTIMFLAGS = ' CXXOPTIMFLAGS=-O3 ';
COPTIMFLAGS = ' COPTIMFLAGS=-O3 ';
CXXDEBUGFLAGS = ' CXXDEBUGFLAGS= ';
CDEBUGFLAGS = ' CDEBUGFLAGS= ';
LDOPTIMFLAGS = ' LDOPTIMFLAGS=-O3 ';
LDDEBUGFLAGS = ' LDDEBUGFLAGS= ';
COMPILE_OPTIONS = [ COMPILE_OPTIONS CDEBUGFLAGS COPTIMFLAGS CXXDEBUGFLAGS CXXOPTIMFLAGS LDDEBUGFLAGS LDOPTIMFLAGS];

% Comment next line to suppress compilation debugging info
% COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -v' ];

COMPILE_COMMAND = [ 'mex ' COMPILE_OPTIONS ' -outdir ' OUTPUT_DIR ];

disp('Compiling mjdgges...')
eval([ COMPILE_COMMAND ' mjdgges/mjdgges.c ' LAPACK_PATH ]);

disp('Compiling sparse_hessian_times_B_kronecker_C...')
eval([ COMPILE_COMMAND ' kronecker/sparse_hessian_times_B_kronecker_C.cc' ]);

disp('Compiling A_times_B_kronecker_C...')
eval([ COMPILE_COMMAND ' kronecker/A_times_B_kronecker_C.cc ' BLAS_PATH]);

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
eval([ COMPILE_COMMAND ' -Isimulate -I../../preprocessor simulate/simulate.cc simulate/Interpreter.cc simulate/Mem_Mngr.cc simulate/SparseMatrix.cc']);