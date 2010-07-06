
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdarg.h>
#include <string.h>

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#include <dynmex.h>
#endif

int constant_seed;

void
swz_fprintf_err(char *str, ...)
{
  va_list ap;
  va_start(ap, str);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexPrintf(str, ap);
#else
  vfprintf(stderr, str, ap);
#endif

  va_end(ap);
}

void
swzExit(int status)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexErrMsgTxt("Error in mexfile.\n");
#else
  exit(status);
#endif
}
