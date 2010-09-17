/*
 * Copyright (C) 2010 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

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

  const char *mainarg = "./a.out";
  char *argument = NULL;
  char *beginarg = NULL;
  char **args = NULL;

  /*
   * Check args
   */
  if (nrhs != 1 || !mxIsChar(prhs[0]) || nlhs != 1)
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in swz_mex: this function takes 1 string input argument and returns 1 output argument.");

  /*
   * Allocate memory
   */
  maxnargs = (int)(mxGetN(prhs[0])/2+1);
  argument = (char *)swzCalloc(mxGetN(prhs[0])+1, sizeof(char));
  args = (char **)swzCalloc(maxnargs, sizeof(char *));
  if (argument==NULL || args==NULL)
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in swz_mex: could not allocate memory. (1)");

  /*
   * Create argument string from prhs and parse to create args / nargs
   */
  if (!(args[nargs] = (char *)swzCalloc(strlen(mainarg)+1, sizeof(char))))
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in swz_mex: could not allocate memory. (2)");

  strncpy(args[nargs++], mainarg, strlen(mainarg));

  if (mxGetString(prhs[0], argument, mxGetN(prhs[0])+1))
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in swz_mex: error using mxGetString.\n");

  beginarg = &argument[0];
  while((n=strcspn(beginarg, " ")))
    {
      if (!(args[nargs] = (char *)swzCalloc(n+1, sizeof(char))))
        DYN_MEX_FUNC_ERR_MSG_TXT("Error in swz_mex: could not allocate memory. (3)");

      strncpy(args[nargs++], beginarg, n);
      beginarg += (isspace(beginarg[n]) || isblank(beginarg[n]) ? ++n : n);
    }
  swzFree(argument);

  /*
   * Call top_level function (formerly main)
   */
  try
    {
      main(nargs, args);
    }
  catch (const char *str)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(str);
    }

  /*
   * free memory
   */
  for (n=0; n<nargs; n++)
    swzFree(args[n]);
  swzFree(args);

  plhs[0] = mxCreateDoubleScalar(0);
}
