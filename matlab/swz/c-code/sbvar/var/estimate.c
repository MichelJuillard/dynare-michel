
#include "switch.h"
#include "switchio.h"
#include "switch_opt.h"
#include "VARbase.h"
#include "VARio.h"
#include "dw_error.h"
#include "dw_ascii.h"
#include "dw_parse_cmd.h"
#include "dw_rand.h"
#include "command_line_VAR.h"

#include "optpackage.h"
//#include "csminwel.h"
//#include "dw_csminwel.h"

#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "modify_for_mex.h"

#define FIND_POSTERIOR_MODE  1
#define FIND_LIKELIHOOD_MODE 2

typedef struct 
{
  int type;

  TVARCommandLine *cmd;

  char *csminwel_output_filename;
  char *intermediate_output_filename;

  PRECISION criterion_start;
  PRECISION criterion_end;
  PRECISION criterion_increment;

  int max_iterations_start;
  PRECISION max_iterations_increment;

  int max_block_iterations;

} TEstimateInfo;


void FindMode_VAR_csminwel(TStateModel *model, TEstimateInfo *estimate)
{
  int iteration, total_iteration, i, j, size_VAR, pos_VAR, size_Q, pos_Q;
  double objective, objective_last, likelihood, prior;
  int **block;
  FILE *f_out;
  char *header, *fmt="Iteration %d: ";

  // csminwel arguments 
  int itct, fcount, retcodeh, nit;
  double *x, fh, crit;
  TMatrix H;
  TVector g;

  f_out=dw_CreateTextFile(estimate->intermediate_output_filename);

  //==== Allocate memory  ===
  size_VAR=NumberFreeParametersTheta(model);
  size_Q=NumberFreeParametersQ(model);
  pos_VAR=0;
  pos_Q=size_VAR;
  x=(double*)malloc((size_VAR + size_Q)*sizeof(double));

  //=== Set starting value ===
  ConvertQToFreeParameters(model,x+pos_Q);
  ConvertThetaToFreeParameters(model,x+pos_VAR);

  //=== Set csminwel output file ===
  csminwel_SetPrintFile(estimate->csminwel_output_filename);

  //=== Print Initial Values ===
  fprintf(f_out,"\n//=== Initial Values ===//\n");
  fprintf(f_out,"Likelihood value:  %22.14le\n",objective=likelihood=LogLikelihood_StatesIntegratedOut(model));
  fprintf(f_out,"Prior value:  %22.14le\n",prior=LogPrior(model));
  fprintf(f_out,"Posterior value:  %22.14le\n\n",likelihood+prior);

  header="Initial: ";
  WriteTransitionMatrices(f_out,(char*)NULL,header,model);
  Write_VAR_Parameters(f_out,(char*)NULL,header,model);
  fflush(f_out);

  //=== Create blocking structure ===
  block=dw_CreateRectangularArray_int(2,2);
  block[0][0]=size_VAR;         block[0][1]=pos_VAR;
  block[1][0]=size_Q;           block[1][1]=pos_Q;

  //=== Objective ===
  if (estimate->type == FIND_POSTERIOR_MODE)
    objective=likelihood+prior;
  else
    objective=likelihood;

  for (total_iteration=1, crit=estimate->criterion_start, nit=estimate->max_iterations_start; 
       crit >= estimate->criterion_end; 
       crit*=estimate->criterion_increment, nit*=(int)estimate->max_iterations_increment)
    {
      for (iteration=1; iteration <= estimate->max_block_iterations; total_iteration++, iteration++)
	{
	  objective_last=objective;

	  fprintf(f_out,"\n\n//=== Iteration %d ===//\n",total_iteration);
	  fprintf(f_out,"Criterion/Max Iteration:  %le  %d\n",crit,nit);
	  fprintf(f_out,"Previous likelihood value:  %22.14le\n",likelihood);
	  fprintf(f_out,"Previous prior value:  %22.14le\n",prior);
	  fprintf(f_out,"Previous posterior value:  %22.14le\n\n",prior+likelihood);
	  fflush(f_out);

	  for (i=0; i < dw_DimA(block); i++)
	    if (block[i][0] > 0)
	      {
		g=CreateVector(block[i][0]);
		H=IdentityMatrix((TMatrix)NULL,block[i][0]);
		ProductMS(H,H,INI_H_CSMINWEL);

		SetupObjectiveFunction(model,x+block[i][1],x+pos_Q,x+pos_VAR);

		if (estimate->type == FIND_POSTERIOR_MODE)
		  csminwel(PosteriorObjectiveFunction_csminwel,x+block[i][1],block[i][0],pElementM(H),pElementV(g),NULL,
			   &fh,crit,&itct,nit,&fcount,&retcodeh,NULL,NULL);
		else
		  csminwel(MLEObjectiveFunction_csminwel,x+block[i][1],block[i][0],pElementM(H),pElementV(g),NULL,
			   &fh,crit,&itct,nit,&fcount,&retcodeh,NULL,NULL);

		ConvertFreeParametersToQ(model,x+pos_Q);
		ConvertFreeParametersToTheta(model,x+pos_VAR);

		FreeMatrix(H);
		FreeVector(g);

		fprintf(f_out,"Likelihood value after pass %d:  %22.14le\n",i,likelihood=LogLikelihood_StatesIntegratedOut(model));
		fprintf(f_out,"Prior value after pass %d:  %22.14le\n",i,prior=LogPrior(model));
		fprintf(f_out,"Posterior value after pass %d:  %22.14le\n",i,likelihood+prior);
		fprintf(f_out,"Csminwel return code: %d\n\n",retcodeh);
		fflush(f_out);
	      }

	  for (j=10, i=1; total_iteration >= j; j*=10, i++);
	  sprintf(header=(char*)malloc(strlen(fmt) + i - 1),fmt,total_iteration);
	  WriteTransitionMatrices(f_out,(char*)NULL,header,model);
	  Write_VAR_Parameters(f_out,(char*)NULL,header,model);
	  free(header);
	  fflush(f_out);

	  if (estimate->type == FIND_POSTERIOR_MODE)
	    objective=likelihood+prior;
	  else
	    objective=likelihood;

	  if (fabs(objective - objective_last) <= crit) break;
	}

      objective_last=objective;

      fprintf(f_out,"\n\n//=== Iteration %d ===//\n",++total_iteration);
      fprintf(f_out,"Criterion/Max Iteration:  %le  %d\n",crit,nit);
      fprintf(f_out,"Previous likelihood value:  %22.14le\n",likelihood);
      fprintf(f_out,"Previous prior value:  %22.14le\n",prior);
      fprintf(f_out,"Previous posterior value:  %22.14le\n\n",prior+likelihood);
      fflush(f_out);

      g=CreateVector(pos_Q+pos_VAR);
      H=IdentityMatrix((TMatrix)NULL,pos_Q+pos_VAR);
      ProductMS(H,H,INI_H_CSMINWEL);

      SetupObjectiveFunction(model,x,x+pos_Q,x+pos_VAR);

      if (estimate->type == FIND_POSTERIOR_MODE)
	csminwel(PosteriorObjectiveFunction_csminwel,x,pos_Q+pos_VAR,pElementM(H),pElementV(g),NULL,
		      &fh,crit,&itct,nit,&fcount,&retcodeh,NULL,NULL);
      else
	csminwel(MLEObjectiveFunction_csminwel,x,pos_Q+pos_VAR,pElementM(H),pElementV(g),NULL,
		      &fh,crit,&itct,nit,&fcount,&retcodeh,NULL,NULL);

      ConvertFreeParametersToQ(model,x+pos_Q);
      ConvertFreeParametersToTheta(model,x+pos_VAR);

      FreeMatrix(H);
      FreeVector(g);

      fprintf(f_out,"Likelihood value:  %22.14le\n",likelihood=LogLikelihood_StatesIntegratedOut(model));
      fprintf(f_out,"Prior value:  %22.14le\n",prior=LogPrior(model));
      fprintf(f_out,"Posterior value:  %22.14le\n",likelihood+prior);
      fprintf(f_out,"Csminwel return code: %d\n\n",retcodeh);
      fflush(f_out);

      for (j=10, i=1; total_iteration >= j; j*=10, i++);
      sprintf(header=(char*)malloc(strlen(fmt) + i - 1),fmt,total_iteration);
      WriteTransitionMatrices(f_out,(char*)NULL,header,model);
      Write_VAR_Parameters(f_out,(char*)NULL,header,model);
      free(header);
      fflush(f_out);

      if (estimate->type == FIND_POSTERIOR_MODE)
	objective=likelihood+prior;
      else
	objective=likelihood;
    }

  //=== Free memory ===
  free(x);
  dw_FreeArray(block);

  //=== Close File ===
  fclose(f_out);
}

/*
   filename -  
*
char** ReadInputFile(char *filename)
{
  char **args=(char**)NULL;
  char ***X;
  int i, j, n;
  X=dw_ReadDelimitedFile((FILE*)NULL,filename,' ',REMOVE_EMPTY_FIELDS | STRIP_WHITESPACE);
  if (X)
    {
      for (n=i=0; i < dw_DimA(X); i++)
	if (X[i]) n+=dw_DimA(X[i]);
      if (n > 0)
	{
	  args=dw_CreateArray_string(n);
	  for (n=i=0; i < dw_DimA(X); i++)
	    if (X[i])
	      for (j=0; j < dw_DimA(X[i]); j++)
		{
		  args[n]=X[i][j];
		  X[i][j]=(char*)NULL;
		}
	}
      dw_FreeArray(X);
    }
  return args;
}*/

/*
   Attempts to load parameters from given file.
*
int LoadParameters(FILE *f, char *filename, TStateModel *model)
{
  int i, terminal_errors, rtrn=0;
  char *header[5]={"","Initial: ","Current: ","Posterior mode: ","MLE: "};
  FILE *f_in=f ? f : dw_OpenTextFile(filename);

  terminal_errors=dw_SetTerminalErrors(ALL_ERRORS & (~USER_ERR));

  for (i=0; i < 5; i++)
    if (ReadTransitionMatrices(f_in,(char*)NULL,header[i],model) && Read_VAR_Parameters(f_in,(char*)NULL,header[i],model))
      {
	rtrn=1;
	break;
      }

  dw_SetTerminalErrors(terminal_errors);

  if (!f) fclose(f_in);
  return rtrn;
}*/

/*
   Attempts to get the parameters from the last iteration in the intermediate file.
*
int GetLastIteration(FILE *f_in, TStateModel *model, TEstimateInfo *estimate)
{
  char *id, *header, *fmt="//=== Iteration %d ===//";
  int terminal_errors, i, j, k=0;

  while (1)
    {
      for (j=10, i=1; k+1 >= j; j*=10, i++);
      sprintf(id=(char*)malloc(strlen(fmt) + i - 1),fmt,k+1);

      if (!dw_SetFilePosition(f_in,id))
	{
	  free(id);
	  terminal_errors=dw_SetTerminalErrors(ALL_ERRORS & (~USER_ERR));

	  fmt="Iteration %d: ";
          while (k > 0)
	    {
	      for (j=10, i=1; k >= j; j*=10, i++);
	      sprintf(header=(char*)malloc(strlen(fmt) + i - 1),fmt,k);
	      if (ReadTransitionMatrices(f_in,(char*)NULL,header,model) && Read_VAR_Parameters(f_in,(char*)NULL,header,model))
		{
		  printf("Using intermediate output - %s\n",header);
		  estimate->initialization_header=header;
		  dw_SetTerminalErrors(terminal_errors);
		  return 1;
		}
	      free(header);
	      k--;
	    }

	  dw_SetTerminalErrors(terminal_errors);
	  return 0;
	}

      free(id);
      k++;
    }
}

/*
   Attempt to set up model from command line.  Command line options are the following

   -di <directory>
      If this argument exists, then all input files are in specified directory.  

   -ft <filename tag>
      If this argument exists, then the following is attempted:

         specification file name:  est_final_<tag>.dat
         init/restart file name:   est_final_<tag>.dat with header="Posterior mode: "

         specification file name:  init_<tag>.dat
         init/restart file name:   est_intermediate_<tag>.dat with header="Iteration %d: "

         (not yet implemented)
         specification file name:  init_<tag>.dat
         init/restart file name:   est_csminwel_<tag>.dat  

         specification file name:  init_<tag>.dat
         init/restart file name:   init_<tag>.dat with header="Initial: "
    
      Failure to load both the specification and restart/init files causes the routine to exit.
   
   -fs <filename>
      If this argument exists, then the specification file name is <filename>.  The argument -ft
      takes precedence over -fs.

   -fr <filename>
      If this argument exists, then the init/restart file name is <filename>.  Must be used in 
      conjunction with the argument -fs.  The default value is the filename associated with the
      argument -fs.

   -rh <header>
      If this argument exists, then the header for the init/restart file is <header>.  Must be 
      used in conjuction with the arguments -fr or -fs.  The default value is "". 

   If no command line options are given, then attemps to use a default input file 
   with the name "default.ini".  Returns one valid pointer to a TStateModel upon
   success and null upon failure.
*
TStateModel* GetModelFromCommandLine(int nargs, char **args, TEstimateInfo *estimate)
{
  TStateModel *model;
  char *d1, *d2, *tag, *header, *filename, *fmt;
  FILE *f_in;

  d1=dw_ParseString_String(nargs,args,"di","");
  if (d1[0] && d1[strlen(d1)-1] != '/')
    {
      d2=(char*)malloc(strlen(d1)+2);
      strcat(strcpy(d2,d1),"/");
      d1=d2;
    }
  else
    d2=(char*)NULL;

  if (tag=dw_ParseString_String(nargs,args,"ft",(char*)NULL))
    {
      fmt="%sest_final_%s.dat";
      sprintf(filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
      if (f_in=fopen(filename,"rt"))
	{
	  model=Read_VAR_Specification(f_in,(char*)NULL);
	  header=dw_ParseString_String(nargs,args,"rh","Posterior mode: ");
	  ReadTransitionMatrices(f_in,(char*)NULL,header,model);
	  Read_VAR_Parameters(f_in,(char*)NULL,header,model);
	  fclose(f_in);
	  printf("Using final output\n");
	  estimate->specification_filename=filename;
	  estimate->initialization_filename=filename;
	  estimate->initialization_header=header;	  
	  if (d2) free(d2);
	  return model;
	}
      free(filename);

      fmt="%sinit_%s.dat";
      sprintf(filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
      if (f_in=fopen(filename,"rt"))
	{
	  model=Read_VAR_Specification(f_in,(char*)NULL);
	  estimate->specification_filename=filename;
	  fclose(f_in);

	  fmt="%sest_intermediate_%s.dat";
	  sprintf(filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
          if (f_in=fopen(filename,"rt"))
	    {
              if (GetLastIteration(f_in,model,estimate))
		{		  
		  fclose(f_in);
		  estimate->initialization_filename=filename;
		  if (d2) free(d2);
		  return model;
		}
	      fclose(f_in);
	    }
	  free(filename);

	  fmt="%sinit_%s.dat";
	  sprintf(filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
	  if (f_in=fopen(filename,"rt"))
	    {
	      header=dw_ParseString_String(nargs,args,"rh","Initial: ");
	      ReadTransitionMatrices(f_in,(char*)NULL,header,model);
	      Read_VAR_Parameters(f_in,(char*)NULL,header,model);
	      fclose(f_in);
	      printf("Using initial data\n");
	      estimate->initialization_filename=filename;
	      estimate->initialization_header=header;
	      if (d2) free(d2);	    
	      return model;
	    }

	  FreeStateModel(model);
	}
      free(filename);

      //if (d2) free(d2);
      //swz_fprintf_err("GetModelFromCommandLine():  Unable to create model.\n");
      goto ERROR;
    }

  if (tag=dw_ParseString_String(nargs,args,"fs",(char*)NULL))
    {
      sprintf(filename=(char*)malloc(strlen(d1) + strlen(tag) + 1),"%s%s",d1,tag);
      model=Read_VAR_Specification((FILE*)NULL,filename);
      estimate->specification_filename=filename;

      if (!(tag=dw_ParseString_String(nargs,args,"fr",(char*)NULL)))
	tag=dw_ParseString_String(nargs,args,"fs",(char*)NULL);
      sprintf(filename=(char*)malloc(strlen(d1) + strlen(tag) + 1),"%s%s",d1,tag);
      header=dw_ParseString_String(nargs,args,"rh","");
      ReadTransitionMatrices((FILE*)NULL,filename,header,model);
      Read_VAR_Parameters((FILE*)NULL,filename,header,model);
      estimate->initialization_filename=filename;
      estimate->initialization_header=header;

      if (d2) free(d2);
      return model;
    }

 ERROR:
  if (d2) free(d2);
  //swz_fprintf_err("GetModelFromCommandLine():  No specification file defined.\n");
  return (TStateModel*)NULL;
}

/*
   Attempt to set up model from command line.  Command line options are the following

   -do <directory>
      If this argument exists, then all output files are put in the specified directory.

   -fo <filename tag>
      If this argument exists, then the output files are

         est_csminwel_<tag>.dat
         est_intermediate_<tag>.dat
         est_final_<tag>.dat

      The default value is the filename tag associated with the argument -ft if it exists.  Otherwise
      it is "default". 

   //--- this is yet to be implemented
   -fa <filename>
      Aux output file.  The default value is est_aux_<filename tag>.dat.

   -MLE
      Find the maximum likelihood estimate

   -PM (default)
      Find the posterior mode

   -cb <floating point number> (default = 1.0e-3)
      Beginning csminwel exit criterion 

   -ce <floating point number> (default = 1.03-6)
      Ending csminwel exit criterion

   -ci <floating point number> (default = 0.1)
      csminwel exit criterion increment multiplier

   -ib <integer> (default = 50)
      Beginning csminwel maximum iteration value

   -ii <floating point number> (default = 2)
      csminwel maximum interation increment multiplier

   If no command line options are given, then attemps to use a default input file 
   with the name "default.ini".  Returns one valid pointer to a TStateModel upon
   success and null upon failure.
*
TEstimateInfo* GetEstimateInfoFromCommandLine(int nargs, char **args) //, TStateModel* model)
{
  TEstimateInfo *estimate;
  char *d1, *d2, *tag, *fmt;

  estimate=(TEstimateInfo*)malloc(sizeof(TEstimateInfo));

  // Output directory
  d1=dw_ParseString_String(nargs,args,"di","");
  if (d1[0] && d1[strlen(d1)-1] != '/')
    {
      d2=(char*)malloc(strlen(d1)+2);
      strcat(strcpy(d2,d1),"/");
      d1=d2;
    }
  else
    d2=(char*)NULL;

  // Output filenames
  if (!(tag=dw_ParseString_String(nargs,args,"fo",(char*)NULL)))
    tag=dw_ParseString_String(nargs,args,"ft","default");
  fmt="%sest_csminwel_%s.dat";
  sprintf(estimate->csminwel_output_filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
  fmt="%sest_intermediate_%s.dat";
  sprintf(estimate->intermediate_output_filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
  fmt="%sest_final_%s.dat";
  sprintf(estimate->final_output_filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
  fmt="%sest_aux_%s.dat";
  sprintf(estimate->aux_output_filename=(char*)malloc(strlen(d1) + strlen(fmt) + strlen(tag) - 3),fmt,d1,tag);
  if (d2) free(d2);

  // Posterior mode or MLE
  estimate->type=(dw_FindArgument_String(nargs,args,"MLE") >= 0) ? FIND_LIKELIHOOD_MODE : FIND_POSTERIOR_MODE;

  // Default values
  estimate->criterion_start=dw_ParseFloating_String(nargs,args,"cb",1.0e-3);
  estimate->criterion_end=dw_ParseFloating_String(nargs,args,"ce",1.0e-6);
  estimate->criterion_increment=dw_ParseFloating_String(nargs,args,"ci",0.1);
  estimate->max_iterations_start=dw_ParseInteger_String(nargs,args,"ib",50);
  estimate->max_iterations_increment=dw_ParseFloating_String(nargs,args,"ii",2.0);

  estimate->max_block_iterations=100;

  return estimate;
}

/*
   Creates TStateModel and reads parameters from command line.  Other estimate info
   is also obtained from command line.  
*/
static TStateModel* SetupFromCommandLine(int nargs, char **args, TEstimateInfo **p_info)
{ 
  TEstimateInfo *info;
  
  if (!(*p_info)) *p_info=(TEstimateInfo*)malloc(sizeof(TEstimateInfo));
  info=*p_info;

  info->cmd=Base_VARCommandLine(nargs,args,(TVARCommandLine*)NULL);

  // Posterior mode or MLE
  info->type=(dw_FindArgument_String(nargs,args,"MLE") >= 0) ? FIND_LIKELIHOOD_MODE : FIND_POSTERIOR_MODE;

  // Default values
  info->criterion_start=dw_ParseFloating_String(nargs,args,"cb",1.0e-3);
  info->criterion_end=dw_ParseFloating_String(nargs,args,"ce",1.0e-6);
  info->criterion_increment=dw_ParseFloating_String(nargs,args,"ci",0.1);
  info->max_iterations_start=dw_ParseInteger_String(nargs,args,"ib",50);
  info->max_iterations_increment=dw_ParseFloating_String(nargs,args,"ii",2.0);

  info->max_block_iterations=100;

   // Output filenames
  info->csminwel_output_filename=CreateFilenameFromTag("%sest_csminwel_%s.dat",info->cmd->out_tag,info->cmd->out_directory);
  info->intermediate_output_filename=CreateFilenameFromTag("%sest_intermediate_%s.dat",info->cmd->out_tag,info->cmd->out_directory);

  return CreateTStateModelForEstimate(nargs,args,&(info->cmd));
}

int main(int nargs, char **args)
{
  TStateModel *model;
  TEstimateInfo *estimate=(TEstimateInfo*)NULL;
  char *filename;
  FILE *f_out;
  time_t begin_time, end_time;
  int t, seed;
  TVector y;
  char *include_help[]={"-di","-do","-fs","-fp","-ph","-ft","-fto",(char*)NULL},
       *additional_help[]={
    "-MLE",
    "Finds the maximum likelihood estimate",
    "-PM",
    "Finds the posterior mode (default option)",
    "-cb <floating point number>",
    "Beginning csminwel exit criterion (default = 1.0e-3)",
    "-ce <floating point number>",
    "Ending csminwel exit criterion (default = 1.03-6)",
    "-ci <floating point number>",
    "csminwel exit criterion increment multiplier (default = 0.1)",
    "-ib <integer>",
    "Beginning csminwel maximum iteration value (default = 50)",
    "-ii <floating point number>",
    "csminwel maximum interation increment multiplier (default = 2)",
    "-nd1",
    "Normalize diagonal of A0 to one (flat output only)",
    "-gs <integer>",
    "Seed value for generator - 0 gets seed from clock (default value = 0)",
    (char*)NULL,
    (char*)NULL};

  //=== Help Screen ===
  if (dw_FindArgument_String(nargs,args,"h") != -1)
    {
      printf("print_draws <options>\n");
      PrintHelpMessages(stdout,include_help,additional_help);
      return 0;
    }

  // Generator seed
  seed=dw_ParseInteger_String(nargs,args,"gs",0);
  dw_initialize_generator(seed);

  printf("Reading initial data...\n");
  if (model=SetupFromCommandLine(nargs,args,&estimate))
    {
      // Estimation
      printf("Beginning estimation...\n");
      begin_time=time((time_t*)NULL);
      FindMode_VAR_csminwel(model,estimate);
      end_time=time((time_t*)NULL);

      // Write final output
      filename=CreateFilenameFromTag("%sest_final_%s.dat",estimate->cmd->out_tag,estimate->cmd->out_directory);
      if (f_out=fopen(filename,"wt"))
        {
          Write_VAR_Specification(f_out,(char*)NULL,model);
          fprintf(f_out,"Specification file: %s\n",estimate->cmd->specification_filename_actual);
          fprintf(f_out,"Initialization file: %s\n",estimate->cmd->parameters_filename_actual);
          fprintf(f_out,"Initialization header: \"%s\"\n",estimate->cmd->parameters_header_actual);

          fprintf(f_out,"Number free parameters in transition matrix: %d\n",NumberFreeParametersQ(model));
          fprintf(f_out,"Number free parameters in theta: %d\n",NumberFreeParametersTheta(model));

          fprintf(f_out,"Time stamp:  %s",ctime(&end_time));
          fprintf(f_out,"Elapsed time: %d seconds\n",(int)end_time-(int)begin_time);

          fprintf(f_out,"Likelihood Value: %g\n",LogLikelihood_StatesIntegratedOut(model));
          fprintf(f_out,"Prior Value: %g\n",LogPrior(model));
          fprintf(f_out,"Posterior Value: %g\n\n",LogPosterior_StatesIntegratedOut(model));

          WriteTransitionMatrices(f_out,(char*)NULL,estimate->cmd->out_header,model);
          Write_VAR_Parameters(f_out,(char*)NULL,estimate->cmd->out_header,model);

          fclose(f_out);
        }
      free(filename);

      // Write flat file
      filename=CreateFilenameFromTag("%sest_flat_header_%s.dat",estimate->cmd->out_tag,estimate->cmd->out_directory);
      if (f_out=fopen(filename,"wt"))
        {     
          WriteBaseTransitionMatricesFlat_Headers_SV(f_out,model->sv,"");
          Write_VAR_ParametersFlat_Headers(f_out,model);
          fprintf(f_out,"\n");
          fclose(f_out);
        }
      free(filename);
      filename=CreateFilenameFromTag("%sest_flat_%s.dat",estimate->cmd->out_tag,estimate->cmd->out_directory);
      if (f_out=fopen(filename,"wt"))
        {     
          WriteBaseTransitionMatricesFlat(f_out,model,"%lf ");
          if (dw_FindArgument_String(nargs,args,"nd1") >= 0)
            Write_VAR_ParametersFlat_A0_Diagonal_One(f_out,model,"%lf ");
          else
            Write_VAR_ParametersFlat(f_out,model,"%lf ");
          fprintf(f_out,"\n");
          fclose(f_out);
        }
       free(filename);

      // Write aux output
      filename=CreateFilenameFromTag("%sest_aux_%s.dat",estimate->cmd->out_tag,estimate->cmd->out_directory);
      if (f_out=fopen(filename,"wt"))
        {
          fprintf(f_out,"""ln(P(y[t]|Y[t-1],Z[t],theta,Q))"",""E[y[t]|Y[t-1],Z[t],theta,Q]""\n");

          y=CreateVector(((T_VAR_Parameters*)(model->theta))->nvars);
          for (t=1; t <= model->sv->nobs; t++)
            {
              fprintf(f_out,"%le,",LogConditionalLikelihood_StatesIntegratedOut(t,model));
              if (ExpectationSingleStep_StatesIntegratedOut(y,t,model))
                dw_PrintVector(f_out,y,"%le,");
              else
                fprintf(f_out,"\n");
            }

          FreeVector(y);
          fclose(f_out);  
        }
      free(filename);

      // Free memory
      FreeStateModel(model);
      Free_VARCommandLine(estimate->cmd);
    }
  else
    {
      // unable to create model
      if (estimate)  
        {
          if (estimate->cmd) Free_VARCommandLine(estimate->cmd);
          free(estimate);
        }
    }

  return 0;
}
