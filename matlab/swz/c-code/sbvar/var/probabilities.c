
#include "switch.h"
#include "switchio.h"
#include "VARio.h"
#include "dw_parse_cmd.h"
#include "dw_ascii.h"

#include <stdlib.h>

#include "modify_for_mex.h"

/*
   Attempt to set up model from command line.  Command line options are the
   following

   -ft <filename tag>
      If this argument exists, then the following is attempted:
         specification file name = est_final_<tag>.dat
         output file name        = probabilites_<tag>.dat
         parameters file name    = est_final_<tag>.dat
         header                  = "Posterior mode: "

   -fs <filename>
      If this argument exists, then the specification file name is <filename>.
      The argument -fs takes precedence over -ft.

   -fo <filename>
      If this argument exists, then the output file name is <filename>.  The
      argument -fo takes precedence over -ft.  The default value is
      parameters.dat.

   -fp <filename>
      If this argument exists, then the parameters file name is <filename>.  The
      argument -fp takes precedence over -ft.  The default value is the filename
      associated with the argument -fs.

   -ph <header>
      If this argument exists, then the header for the parameters file is
      <header>.  The default value is "Posterior mode: ".

*/


int main(int nargs, char **args)
{
  char *spec=(char*)NULL, *parm=(char*)NULL, *head=(char*)NULL, *out=(char*)NULL, *buffer, *fmt;
  TStateModel *model;
  TVector *probabilities;
  int s, t;
  FILE *f_out;

/*    // specification filename   ansi-c*/
  if (buffer=dw_ParseString_String(nargs,args,"fs",(char*)NULL))
    strcpy(spec=(char*)malloc(strlen(buffer)+1),buffer);

/*    // output filename   ansi-c*/
  if (buffer=dw_ParseString_String(nargs,args,"fo",(char*)NULL))
    strcpy(out=(char*)malloc(strlen(buffer)+1),buffer);

/*    // parameter filename   ansi-c*/
  if (buffer=dw_ParseString_String(nargs,args,"fp",(char*)NULL))
    strcpy(parm=(char*)malloc(strlen(buffer)+1),buffer);

/*    // header   ansi-c*/
  if (buffer=dw_ParseString_String(nargs,args,"ph",(char*)NULL))
    strcpy(head=(char*)malloc(strlen(buffer)+1),buffer);

/*    // file tag   ansi-c*/
  if (buffer=dw_ParseString_String(nargs,args,"ft",(char*)NULL))
    {
      fmt="est_final_%s.dat";

/*        // specification filename   ansi-c*/
      if (!spec)
    sprintf(spec=(char*)malloc(strlen(fmt) + strlen(buffer) - 1),fmt,buffer);

/*        // parameter filename   ansi-c*/
      if (!parm)
    sprintf(parm=(char*)malloc(strlen(fmt) + strlen(buffer) - 1),fmt,buffer);

/*        // output filename   ansi-c*/
      if (!out)
    {
      fmt="probabilities_%s.dat";
      sprintf(out=(char*)malloc(strlen(fmt) + strlen(buffer) - 1),fmt,buffer);
    }
    }

  if (!spec)
    {
      swz_fprintf_err("No specification filename given\n");
      swz_fprintf_err("Command line syntax:\n"
                      "  -ft : file tag\n"
                      "  -fs : specification filename\n"
                      "  -fo : output filename (probablities.dat)\n"
                      "  -fp : parameters filename (specification filename)\n"
                      "  -fh : parameter header (Posterior mode: )\n"
                      );
      exit(1);
    }

  if (!parm)
    strcpy(parm=(char*)malloc(strlen(spec)+1),spec);

  if (!head)
    {
      buffer="Posterior mode: ";
      strcpy(head=(char*)malloc(strlen(buffer)+1),buffer);
    }

  if (!out)
    {
      buffer="probabilities.dat";
      strcpy(out=(char*)malloc(strlen(buffer)+1),buffer);
    }

  model=Read_VAR_Specification((FILE*)NULL,spec);
  ReadTransitionMatrices((FILE*)NULL,parm,head,model);
  Read_VAR_Parameters((FILE*)NULL,parm,head,model);

  probabilities=dw_CreateArray_vector(model->sv->nstates);
  for (s=model->sv->nstates-1; s >= 0; s--)
    probabilities[s]=ProbabilitiesState((TVector)NULL,s,model);

  f_out=dw_CreateTextFile(out);
  for (t=0; t <= model->sv->nobs; t++)
    {
      for (s=0; s < model->sv->nstates; s++)
    fprintf(f_out,"%lf ",ElementV(probabilities[s],t));
      fprintf(f_out,"\n");
    }

  free(spec);
  free(out);
  free(head);
  free(parm);

}
