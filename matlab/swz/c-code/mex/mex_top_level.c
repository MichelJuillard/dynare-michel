#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dynmex.h>

#include "modify_for_mex.h"

int main(int nargs, char **args);

/* MATLAB interface */
void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  int nargs = 0;
  int n = 0;
  int maxnargs = 0;

  char *mainarg = "./a.out";
  char *argument = NULL;
  char *beginarg = NULL;
  char **args = NULL;

  /*
   * Check args
   */
  if (nrhs != 1 || !mxIsChar(prhs[0]))
    mexErrMsgTxt("This function takes only one string argument.");

  if (nlhs != 0)
    mexWarnMsgTxt("This function has no return arguments.\n");

  /*
   * Allocate memory
   */
  maxnargs = (int)(mxGetN(prhs[0])/2+1);
  argument = (char *)swzCalloc(mxGetN(prhs[0])+1, sizeof(char));
  args = (char **)swzCalloc(maxnargs, sizeof(char *));
  if (argument==NULL || args==NULL)
    mexErrMsgTxt("In swz_mex: could not allocate memory. (1)");

  /*
   * Create argument string from prhs and parse to create args / nargs
   */
  if (!(args[nargs] = (char *)swzCalloc(strlen(mainarg)+1, sizeof(char))))
    mexErrMsgTxt("In swz_mex: could not allocate memory. (2)");
  strncpy(args[nargs++], mainarg, strlen(mainarg));

  if (mxGetString(prhs[0], argument, mxGetN(prhs[0])+1))
    mexErrMsgTxt("In swz_mex: error using mxGetString.\n");

  beginarg = &argument[0];
  while(n=strcspn(beginarg, " "))
    {
      if (!(args[nargs] = (char *)swzCalloc(n+1, sizeof(char))))
        mexErrMsgTxt("In swz_mex: could not allocate memory. (3)");
      strncpy(args[nargs++], beginarg, n);
      beginarg += (isspace(beginarg[n]) || isblank(beginarg[n]) ? ++n : n);
    }
  swzFree(argument);

  /*
   * Call top_level function (formerly main)
   */
  main(nargs, args);

  /*
   * free memory
   */
  for (n=0; n<nargs; n++)
    swzFree(args[n]);
  swzFree(args);
}
