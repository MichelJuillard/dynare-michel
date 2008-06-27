% Build file for Dynare MEX Librairies for Octave
% Copyright Dynare Team (2008)
% GNU Public License

COMPILE_OPTIONS = '-DNO_BLAS_H -DNO_LAPACK_H -DOCTAVE';

% Comment next line to suppress compilation debugging info
COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -v' ];

COMPILE_COMMAND = [ 'mex ' COMPILE_OPTIONS ];

disp('Compiling mjdgges...')
eval([ COMPILE_COMMAND ' mjdgges/mjdgges.c -o ../octave/mjdgges.mex']);
disp('Compiling sparse_hessian_times_B_kronecker_C...')
eval([ COMPILE_COMMAND ' kronecker/sparse_hessian_times_B_kronecker_C.cc -o ../octave/sparse_hessian_times_B_kronecker_C.mex']);
disp('Compiling A_times_B_kronecker_C...')
eval([ COMPILE_COMMAND ' kronecker/A_times_B_kronecker_C.cc -o ../octave/A_times_B_kronecker_C.mex']);
disp('Compiling gensylv...')
% With MS Visual C++ it is necessary to list all the source file names (a wildcard will not make it)
eval([ COMPILE_COMMAND ' -DMATLAB -Igensylv/cc gensylv/matlab/gensylv.cpp gensylv/cc/BlockDiagonal.cpp gensylv/cc/GeneralMatrix.cpp gensylv/cc/GeneralSylvester.cpp gensylv/cc/IterativeSylvester.cpp gensylv/cc/KronUtils.cpp gensylv/cc/KronVector.cpp gensylv/cc/QuasiTriangular.cpp gensylv/cc/QuasiTriangularZero.cpp gensylv/cc/SchurDecomp.cpp gensylv/cc/SchurDecompEig.cpp gensylv/cc/SimilarityDecomp.cpp gensylv/cc/SylvException.cpp gensylv/cc/SylvMatrix.cpp gensylv/cc/SylvMemory.cpp gensylv/cc/SylvParams.cpp gensylv/cc/TriangularSylvester.cpp gensylv/cc/Vector.cpp -o ../octave/gensylv.mex']);
disp('Compiling simulate...')
eval([ COMPILE_COMMAND ' -Isimulate -I../../preprocessor/include simulate/simulate.cc simulate/Interpreter.cc simulate/Mem_Mngr.cc simulate/SparseMatrix.cc simulate/linbcg.cc -o ../octave/simulate.mex']);
