% Build file for Dynare MEX Librairies
% Copyright Dynare Team (2007)
% GNU Public License

MATLAB  = ver('matlab');

% FIXME:
% It's not satisfactory to convert string versions into numbers, and to
% compare these numbers:
% - conversion will fail if version = 1.2.3
% - it will give 7.10 < 7.9
VERSION = str2num(MATLAB.Version);

MATLAB_PATH = matlabroot;

COMPILE_OPTIONS = '';

if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer) % GNU/Linux (i386)
    if strcmpi('GLNXA64', computer) % 64 bits ?
        COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
    end
    LAPACK_PATH = '-lmwlapack';
    if VERSION <= 7.4
        BLAS_PATH = LAPACK_PATH; % On <= 7.4, BLAS in included in LAPACK
        COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DNO_BLAS_H' ];
    else
        BLAS_PATH = '-lmwblas';
    end
elseif strcmpi('PCWIN', computer) % Windows (32 bits)
    LAPACK_PATH = 'libmwlapack.lib';
    if VERSION <= 7.4
        BLAS_PATH = LAPACK_PATH; % On <= 7.4, BLAS in included in LAPACK
        COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMWTYPES_NOT_DEFINED -DNO_BLAS_H' ];
    else
        BLAS_PATH = 'libmwblas.lib';
    end
else
    error('Unsupported platform')
end

if VERSION <= 7.4
    OUTPUT_DIR = '../2007a';
else
    OUTPUT_DIR = '../2007b';
end

% Comment next line to suppress compilation debugging info
COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -v' ]

COMPILE_COMMAND = [ 'mex ' COMPILE_OPTIONS ' -outdir ' OUTPUT_DIR ' ' ];

disp('Compiling mjdgges...')
system([ COMPILE_COMMAND ' mjdgges/mjdgges.c ' LAPACK_PATH ]);

disp('Compiling sparse_hessian_times_B_kronecker_C...')
system([ COMPILE_COMMAND ' kronecker/sparse_hessian_times_B_kronecker_C.cc ' BLAS_PATH ]);

disp('Compiling A_times_B_kronecker_C...')
system([ COMPILE_COMMAND ' kronecker/A_times_B_kronecker_C.cc ' BLAS_PATH ]);

disp('Compiling gensylv...')
system([ COMPILE_COMMAND ' -DMATLAB -Igensylv/cc gensylv/matlab/gensylv.cpp gensylv/cc/*.cpp ' BLAS_PATH ' ' LAPACK_PATH ]);
