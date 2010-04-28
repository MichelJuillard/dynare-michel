#include <stdlib.h>
#include "mex.h"

void
swz_exit(int status)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexErrMsgTxt("Error in mexfile.\n");
#else
  exit(status);
#endif
}
