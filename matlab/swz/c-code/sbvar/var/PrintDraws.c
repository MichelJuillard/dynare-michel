
#include "swzmatrix.h"
#include "dw_rand.h"
#include "dw_parse_cmd.h"
#include "dw_ascii.h"
#include "VARbase.h"
#include "VARio.h"
#include "switch.h"
#include "switchio.h"
#include "command_line_VAR.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "modify_for_mex.h"

int main(int nargs, char **args)
{
  TStateModel *model;
  T_VAR_Parameters *p;
  FILE *f_out;
  char *filename;
  int count, begin_time, end_time, tuning, burn_in, iterations, check, period=1000, seed, output, thinning,
      nd1;
  TVARCommandLine *cmd=(TVARCommandLine*)NULL;
  char *include_help[]={"-di","-do","-fs","-fp","-ph","-MLE",(char*)NULL},
       *additional_help[]={
    "-ft <tag>",
    "Tag for input file.  Input file name is est_final_<tag>.dat.",
    "-fto <tag>",
    "Tag for output file.  Output file names are draws_<tag>.dat and headers_<tag>.dat.  Default is -ft <tag>.",
    "-mh <integer>",
    "Tuning period for Metropolis-Hasting draws (default value = 30000)",
    "-b <integer>",
    "Burn-in period (default value = 0.1 * (number of iterations))",
    "-i <integer>",
    "Number of draw (default value = 1000), number saved is (-i)/(-t)",
    "-t <integer>",
    "Thinning factor.  Only one in t draws are written to file.",
    "-nd1",
    "Normalize diagonal of A0 to one (flat output only)",
    "-gs <integer>",
    "Seed value for generator - 0 gets seed from clock (default value = 0)",
    (char*)NULL,
    (char*)NULL};

/*    //=== Help Screen ===   ansi-c*/
  if (dw_FindArgument_String(nargs,args,"h") != -1)
    {
      printf("print_draws <options>\n");
      PrintHelpMessages(stdout,include_help,additional_help);
      return 0;
    }

/*    //=== Get seed, tuning peroid, burn-in period, number of iterations, and thinning factor   ansi-c*/
  seed=dw_ParseInteger_String(nargs,args,"gs",0);
  tuning=dw_ParseInteger_String(nargs,args,"mh",30000);
  iterations=dw_ParseInteger_String(nargs,args,"i",1000);
  burn_in=dw_ParseInteger_String(nargs,args,"b",iterations/10);
  thinning=dw_ParseInteger_String(nargs,args,"t",1);
  nd1=(dw_FindArgument_String(nargs,args,"nd1") >= 0) ? 1 : 0;

/*    //=== Initialize random number generator   ansi-c*/
  dw_initialize_generator(seed);

/*    //=== Setup model and initial parameters   ansi-c*/
  printf("Reading data...\n");
  if (!(model=CreateTStateModelFromEstimateFinal(nargs,args,&cmd)))
    {
      swz_fprintf_err("Unable to read model or parameters\n");
      exit(1);
    }
  p=(T_VAR_Parameters*)(model->theta);

/*    //=== Open header file and print headers   ansi-c*/
  filename=CreateFilenameFromTag("%sheader_%s.dat",cmd->out_tag,cmd->out_directory);
  f_out=fopen(filename,"wt");
  free(filename);
  WriteBaseTransitionMatricesFlat_Headers_SV(f_out,model->sv,"");
  Write_VAR_ParametersFlat_Headers(f_out,model);
  fprintf(f_out,"\n");
  fclose(f_out);

/*    //=== Open output file   ansi-c*/
  filename=CreateFilenameFromTag("%sdraws_%s.dat",cmd->out_tag,cmd->out_directory);
  f_out=fopen(filename,"wt");
  free(filename);

/*    // Burn-in period with calibration of jumping parameters   ansi-c*/
  printf("Calibrating jumping parameters - %d draws\n",tuning);
  begin_time=(int)time((time_t*)NULL);
  AdaptiveMetropolisScale(model,tuning,1000,1,(FILE*)NULL);       /*   tuning iterations - 1000 iterations before updating - verbose   ansi-c*/
  end_time=(int)time((time_t*)NULL);
  printf("Elapsed Time: %d seconds\n",end_time - begin_time);

/*    // Reset parametrers   ansi-c*/
  if (!ReadTransitionMatrices((FILE*)NULL,cmd->parameters_filename_actual,cmd->parameters_header_actual,model)
      || !Read_VAR_Parameters((FILE*)NULL,cmd->parameters_filename_actual,cmd->parameters_header_actual,model))
    printf("Unable to reset parameters after tuning\n");

/*    // Burn-in period   ansi-c*/
  printf("Burn-in period - %d draws\n",burn_in);
  for (check=period, count=1; count <= burn_in; count++)
    {
      DrawAll(model);

      if (count == check)
    {
      check+=period;
      printf("%d iterations completed out of %d\n",count,burn_in);
    }
    }
  end_time=(int)time((time_t*)NULL);
  printf("Elapsed Time: %d seconds\n",end_time - begin_time);
  ResetMetropolisInformation(p);

/*    // Simulation   ansi-c*/
  printf("Simulating - %d draws\n",iterations);
  for (check=period, output=thinning, count=1; count <= iterations; count++)
    {
      DrawAll(model);

      if (count == output)
        {
          WriteBaseTransitionMatricesFlat(f_out,model,"%lf ");
          if (nd1)
            Write_VAR_ParametersFlat_A0_Diagonal_One(f_out,model,"%lf ");
          else
            Write_VAR_ParametersFlat(f_out,model,"%lf ");
          fprintf(f_out,"\n");
          output+=thinning;
        }

      if (count == check)
    {
      check+=period;
      printf("%d(%d) iterations completed out of %d(%d)\n",count,thinning,iterations,thinning);
    }
    }
  end_time=(int)time((time_t*)NULL);
  printf("Elapsed Time: %d seconds\n",end_time - begin_time);

/*    // clean up   ansi-c*/
  fclose(f_out);
  FreeStateModel(model);
  Free_VARCommandLine(cmd);
}
