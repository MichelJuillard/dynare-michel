
#include "switch.h"
#include "switchio.h"
#include "VARbase.h"
#include "VARio.h"
#include "VARio_matlab.h"
#include "dw_error.h"
#include "dw_ascii.h"

#include <stdlib.h>
#include <string.h>

#include "modify_for_mex.h"

/*
   Creates a standard initialization file from the matlab and specification file.  
*/
int main(int nargs, char **args)
{
  TStateModel *model;
  FILE *f_out;
  char *filename, *fmt="init_%s.dat", *header="Initial: ";
  
  dw_SetTerminalErrors(ALL_ERRORS);
  dw_SetVerboseErrors(ALL_ERRORS);

  if (nargs != 4)
    {
      fprintf(stderr,"Syntax:\n  create_init_file <matlab filename> <specs filename> <file tag>\n");
      exit(0);
    }

  model=Combine_matlab_standard(args[1],args[2]);
  ReadConstantParameters(args[1],model);
  sprintf(filename=(char*)malloc(strlen(fmt) + strlen(args[3]) - 1),fmt,args[3]);
  f_out=dw_CreateTextFile(filename);
  Write_VAR_Specification(f_out,(char*)NULL,model);
  WriteTransitionMatrices(f_out,(char*)NULL,header,model);
  Write_VAR_Parameters(f_out,(char*)NULL,header,model);
  fclose(f_out);
  FreeStateModel(model);
  
  return 0;
}
