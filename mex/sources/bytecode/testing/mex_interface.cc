#include "mex_interface.hh"
#include <cstring>
#include <sstream>

using namespace std;

int
mexPrintf(const char *str, ...)
{
  va_list args;
  int retval;

  va_start (args, str);
  retval = vprintf (str, args);
  va_end (args);

  return retval;
}

void
mexErrMsgTxt(const string str)
{
  perror(str.c_str());
  exit(EXIT_FAILURE);
}

void
mxFree(void* to_release)
{
  free(to_release);
}

void*
mxMalloc(int amount)
{
  return malloc(amount);
}

void*
mxRealloc(void* to_extend, int amount)
{
  return realloc(to_extend, amount);
}


void
mexEvalString(const string str)
{
}
