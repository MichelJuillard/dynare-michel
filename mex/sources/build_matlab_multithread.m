% Build file for Dynare MEX Librairies under Matlab with multithreading
% Read http://www.dynare.org/DynareWiki/UsingMultithreadedDlls

% Copyright (C) 2007-2009 Dynare Team
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

% Pass MATLAB_VERSION to C preprocessor in hexadecimal form
verstruct = ver('matlab');
matver = sscanf(verstruct.Version, '%d.%d.%d')';

COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMATLAB_MEX_FILE -DMATLAB_VERSION=0x' sprintf('%02d%02d', matver(1), matver(2)) ];

% Large array dims for 64 bits platforms appeared in Matlab 7.3
if (strcmpi('GLNXA64', computer) || strcmpi('PCWIN64', computer)) ...
      && ~matlab_ver_less_than('7.3')
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
end

OUTPUT_DIR = '../matlab';

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

CFLAGS   = ' CFLAGS="\$CFLAGS -fopenmp" ';
CXXFLAGS = ' CXXFLAGS="\$CXXFLAGS -fopenmp" ';
LDFLAGS  = ' LDFLAGS="\$LDFLAGS -fopenmp" ';

disp('Compiling isopenmp...')
try
    eval([ 'mex ' COMPILE_OPTIONS ' -DUSE_OMP ' CFLAGS CXXFLAGS LDFLAGS  ' -outdir ' OUTPUT_DIR ' threads/isopenmp.cc ' ]);
    disp(' ')
    disp('|------------------------------------------------|')
    disp('|  OpenMp is used (multithreaded mex files) for: |')
    disp('|   * sparse_hessian_times_B_kronecker_C.cc      |')
    disp('|   * A_times_B_kronecker_C.cc                   |')
    disp('|------------------------------------------------|')
    disp(' ')
    COMPILE_OPTIONS_OMP = [ COMPILE_OPTIONS ' -DUSE_OMP ' CFLAGS CXXFLAGS LDFLAGS ];
catch
    disp(' ')
    disp('|------------------------------------------------|')
    disp('|  OpenMp is not available on this platform!     |')
    disp('|------------------------------------------------|')
    disp(' ')
    CFLAGS = []; CXXFLAGS = []; LDFLAGS = [];
    COMPILE_OPTIONS_OMP = COMPILE_OPTIONS;
end

COMPILE_COMMAND_OMP = [ 'mex ' COMPILE_OPTIONS_OMP ' -outdir ' OUTPUT_DIR ];

disp('Compiling mjdgges...')
eval([ COMPILE_COMMAND ' -I. mjdgges/mjdgges.c ' LAPACK_PATH ]);

disp('Compiling sparse_hessian_times_B_kronecker_C...')
eval([ COMPILE_COMMAND_OMP ' -I. kronecker/sparse_hessian_times_B_kronecker_C.cc' ]);
    
disp('Compiling A_times_B_kronecker_C...')
eval([ COMPILE_COMMAND_OMP ' -I. kronecker/A_times_B_kronecker_C.cc ']);    

disp('Compiling gensylv...')
eval([ COMPILE_COMMAND ' -I. -I../../dynare++/sylv/cc ' ...
       '../../dynare++/sylv/matlab/gensylv.cpp ' ...
       '../../dynare++/sylv/cc/BlockDiagonal.cpp ' ... 
       '../../dynare++/sylv/cc/GeneralMatrix.cpp ' ...
       '../../dynare++/sylv/cc/GeneralSylvester.cpp ' ...
       '../../dynare++/sylv/cc/IterativeSylvester.cpp ' ...
       '../../dynare++/sylv/cc/KronUtils.cpp ' ...
       '../../dynare++/sylv/cc/KronVector.cpp ' ...
       '../../dynare++/sylv/cc/QuasiTriangular.cpp ' ...
       '../../dynare++/sylv/cc/QuasiTriangularZero.cpp ' ...
       '../../dynare++/sylv/cc/SchurDecomp.cpp ' ...
       '../../dynare++/sylv/cc/SchurDecompEig.cpp ' ...
       '../../dynare++/sylv/cc/SimilarityDecomp.cpp ' ...
       '../../dynare++/sylv/cc/SylvException.cpp ' ...
       '../../dynare++/sylv/cc/SylvMatrix.cpp ' ...
       '../../dynare++/sylv/cc/SylvMemory.cpp ' ...
       '../../dynare++/sylv/cc/SylvParams.cpp ' ...
       '../../dynare++/sylv/cc/TriangularSylvester.cpp ' ...
       '../../dynare++/sylv/cc/Vector.cpp ' ...
       BLAS_PATH ' ' LAPACK_PATH ]);

disp('Compiling bytecode...')
eval([ COMPILE_COMMAND_OMP ' -Ibytecode -I../../preprocessor bytecode/bytecode.cc bytecode/Interpreter.cc bytecode/Mem_Mngr.cc bytecode/SparseMatrix.cc']);