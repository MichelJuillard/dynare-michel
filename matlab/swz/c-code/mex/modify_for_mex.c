#include "mex.h"

#include <stdlib.h>
#include <stdarg.h>
#include <stdarg.h>
#include <string.h>

void
swz_fprintf_err(FILE *cad, const char * str, ...)
{
  char *whole_str=(char*)NULL;
  char *msg_truncated = ".....MSG TRUNCATED\n";
  int num_args = 0;
  va_list ap;

  va_start(ap, str);

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  if (!(whole_str = (char *)malloc(sizeof(str)*sizeof(char))))
    {
      printf("Could not allocate memory\n");
      exit(0);
    }

  strcpy(whole_str, str);
  while (strtok_r(whole_str, "%", &whole_str) != NULL)
    num_args++;

  num_args = sizeof(str)*sizeof(char) + num_args*sizeof(long double);
  if (!(whole_str = (char *)realloc(whole_str, num_args + strlen(msg_truncated) + 1)))
    {
      printf("Could not allocate memory\n");
      exit(0);
    }

  vsnprintf(whole_str, num_args, str, ap);
  if(strlen(whole_str) + 1 == num_args)
    strcat(whole_str, msg_truncated);

  printf("%s", whole_str);

  free(whole_str);
#else
  vfprintf(stderr, str, ap);
#endif

  va_end(ap);
}

void
swz_exit(int status)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mexErrMsgTxt("Error in mexfile.\n");
#else
  exit(status);
#endif
}
