#ifndef _MEXMOD
#define _MEXMOD
void swz_exit(int status);
void swz_fprintf_err(const char * str, ...);
#endif

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#define printf mexPrintf
#define exit swz_exit
#endif
