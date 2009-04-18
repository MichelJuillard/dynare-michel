#if !defined(MATLAB_VERSIONS_COMPATIBILITY_H)
#define MATLAB_VERSIONS_COMPATIBILITY_H


#if !defined(LAPACK_USE_MWSIGNEDINDEX) || defined(OCTAVE)
typedef int lapack_int;
typedef int blas_int;
#else
typedef mwSignedIndex lapack_int;
typedef mwSignedIndex blas_int;
#endif

#endif
