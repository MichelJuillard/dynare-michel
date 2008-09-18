% Build file for Dynare MEX Librairies under Matlab
% Copyright Dynare Team (2007-2008)
% GNU Public License

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
elseif strcmpi('PCWIN', computer)
    % Windows (x86-32) with Microsoft or gcc compiler
    LIBRARY_PATH = [MATLAB_PATH '/extern/lib/win32/microsoft/'];
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
if strcmpi('GLNXA64', computer) && ~matlab_ver_less_than('7.3')
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
COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -v' ];

COMPILE_COMMAND = [ 'mex ' COMPILE_OPTIONS ' -outdir ' OUTPUT_DIR ];

disp('Compiling mjdgges...')
eval([ COMPILE_COMMAND ' mjdgges/mjdgges.c ' LAPACK_PATH ]);
disp('Compiling sparse_hessian_times_B_kronecker_C...')
eval([ COMPILE_COMMAND ' kronecker/sparse_hessian_times_B_kronecker_C.cc ' BLAS_PATH ]);
disp('Compiling A_times_B_kronecker_C...')
eval([ COMPILE_COMMAND ' kronecker/A_times_B_kronecker_C.cc ' BLAS_PATH ]);
disp('Compiling gensylv...')
eval([ COMPILE_COMMAND ' -DMATLAB -Igensylv/cc gensylv/matlab/gensylv.cpp' ...
       ' gensylv/cc/BlockDiagonal.cpp ' ... 
       ' gensylv/cc/GeneralMatrix.cpp ' ...
       ' gensylv/cc/GeneralSylvester.cpp ' ...
       ' gensylv/cc/IterativeSylvester.cpp ' ...
       ' gensylv/cc/KronUtils.cpp ' ...
       ' gensylv/cc/KronVector.cpp ' ...
       ' gensylv/cc/QuasiTriangular.cpp ' ...
       ' gensylv/cc/QuasiTriangularZero.cpp ' ...
       ' gensylv/cc/SchurDecomp.cpp ' ...
       ' gensylv/cc/SchurDecompEig.cpp ' ...
       ' gensylv/cc/SimilarityDecomp.cpp ' ...
       ' gensylv/cc/SylvException.cpp ' ...
       ' gensylv/cc/SylvMatrix.cpp ' ...
       ' gensylv/cc/SylvMemory.cpp ' ...
       ' gensylv/cc/SylvParams.cpp ' ...
       ' gensylv/cc/TriangularSylvester.cpp ' ...
       ' gensylv/cc/Vector.cpp ' ...
       BLAS_PATH ' ' LAPACK_PATH ]);
disp('Compiling simulate...')
eval([ COMPILE_COMMAND ' -Isimulate -I../../preprocessor/include simulate/simulate.cc simulate/Interpreter.cc simulate/Mem_Mngr.cc simulate/SparseMatrix.cc simulate/linbcg.cc' ]);

% Forcing exit is necessary when autobuilding MEX files for Debian packages.
% In interactive mode, if anything fails during the compilation process,
% this line will not be reached and the user will see the error message.
exit
