#ifndef _MEXMOD
#define _MEXMOD
void swz_exit(int status);
void swz_fprintf_err(const char * str, ...);
extern int constant_seed;
#endif



#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)


#include <dynmex.h>
#include <dynblas.h>
#include <dynlapack.h>

#define swzMalloc mxMalloc
#define swzCalloc mxCalloc
#define swzRealloc mxRealloc
#define swzFree mxFree


#else


#define swz_fprintf_stdout printf
#define swzMalloc malloc
#define swzCalloc calloc
#define swzRealloc realloc
#define swzFree free


#endif
