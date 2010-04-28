
#include "mhm_VAR.h"
#include "VARbase.h"
#include "VARio.h"
#include "switch.h"
#include "switchio.h"
#include "dw_rand.h"
#include "dw_error.h"
#include "dw_ascii.h"
#include "dw_parse_cmd.h"

#include <time.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "modify_for_mex.h"

static void ReadError_MHMio(char *id)
{
  char *errmsg, *fmt="Error after line identifier ""%s""";
  sprintf(errmsg=(char*)malloc(strlen(fmt) + strlen(id) - 1),fmt,id);
  dw_UserError(errmsg);
  free(errmsg);
}

/*
   Creates a copy of d and adds a trailing '/' if necessary.  The returned
   pointer, if not null, must be freed by calling routine.
*/
static char* AddSlash(char *d)
{
  char *d_out;
  int k=strlen(d);
  if (d[0] && d[k-1] != '/')
    {
      d_out=(char*)malloc(k+2);
      strcat(strcpy(d_out,d),"/");
    }
  else
    {
      d_out=(char*)malloc(k+2);
      strcpy(d_out,d);
    }
  return d_out;
}

T_MHM* RestartFromFinalFile(char *filename, T_MHM *mhm)
{
  FILE *f_in;
  char *id;
  TStateModel *model;
  if (f_in=fopen(filename,"rt"))
    {
      id="//== Specification after mhm draws ==//";
      if (dw_SetFilePosition(f_in,id))
	{
	  if (!mhm)
	    {
	      mhm=ReadMHM_Input(f_in,(char*)NULL,(T_MHM*)NULL);
	      mhm->mhm_filename=filename;
	    }
	  mhm->spec_filename=mhm->parameter_filename=filename;

	  model=Read_VAR_Specification(f_in,(char*)NULL);
	  mhm->parameter_header="Posterior mode: ";
	  ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
	  Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);
	  AddStateModel(model,mhm);
	  mhm->parameter_header="Final draw: ";
	  ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
	  Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);

	  mhm->n_burn1=-mhm->n_burn1;
          mhm->n_burn2=-mhm->n_burn2;
          mhm->n_mean_variance=-mhm->n_mean_variance;
	  ReadMeanVariance(f_in,mhm);
	  return mhm;
	}
    }
  return (T_MHM*)NULL;
}

/*

*/
T_MHM* RestartFromIntermediateFile(char *filename, T_MHM *mhm)
{
  FILE *f_in;
  char *id;
  TStateModel *model;
  if (f_in=fopen(filename,"rt"))
    {
      if (!mhm)
	{
	  mhm=ReadMHM_Input(f_in,(char*)NULL,(T_MHM*)NULL);
	  mhm->mhm_filename=filename;
	}
      mhm->spec_filename=mhm->parameter_filename=filename;
      id="//== Specification after mhm draws ==//";
      if (dw_SetFilePosition(f_in,id))
	{
	  model=Read_VAR_Specification(f_in,(char*)NULL);
	  mhm->parameter_header="Posterior mode: ";
	  ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
	  Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);
	  AddStateModel(model,mhm);
	  mhm->parameter_header="Final draw: ";
	  ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
	  Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);

	  mhm->n_burn1=-mhm->n_burn1;
          mhm->n_burn2=-mhm->n_burn2;
          mhm->n_mean_variance=-mhm->n_mean_variance;
	  ReadMeanVariance(f_in,mhm);
	}
      else
	{
	  id="//== Specification after mean-variance estimation ==//";
	  if (dw_SetFilePosition(f_in,id))
	    {
	      model=Read_VAR_Specification(f_in,(char*)NULL);
	      mhm->parameter_header="Posterior mode: ";
	      ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
	      Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);
	      AddStateModel(model,mhm);
	      mhm->parameter_header="Mean-variance: ";
	      ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
	      Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);

	      mhm->n_burn1=-mhm->n_burn1;
	      mhm->n_burn2=-mhm->n_burn2;
	      mhm->n_mean_variance=-mhm->n_mean_variance;
	      ReadMeanVariance(f_in,mhm);
	    }
	  else
	    {
	      id="//== Specification after second burn-in ==//";
	      if (dw_SetFilePosition(f_in,id))
		{
		  model=Read_VAR_Specification(f_in,(char*)NULL);
		  mhm->parameter_header="Posterior mode: ";
		  ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
		  Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);
		  AddStateModel(model,mhm);
		  mhm->parameter_header="Second burn-in: ";
		  ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
		  Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);

		  mhm->n_burn1=-mhm->n_burn1;
		  mhm->n_burn2=-mhm->n_burn2;
		}
	      else
		{
		  id="//== Specification after first burn-in ==//";
		  if (dw_SetFilePosition(f_in,id))
		    {
		      model=Read_VAR_Specification(f_in,(char*)NULL);
		      mhm->parameter_header="Posterior mode: ";
		      ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
		      Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);
		      AddStateModel(model,mhm);
		      mhm->parameter_header="First burn-in: ";
		      ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
		      Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);

		      mhm->n_burn1=-mhm->n_burn1;
		    }
		  else
		    {
		      model=Read_VAR_Specification(f_in,(char*)NULL);
		      mhm->parameter_header="Posterior mode: ";
		      ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
		      Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);
		      AddStateModel(model,mhm);
		    }
		}
	    }
	}
      fclose(f_in);
      return (mhm);
    }
  return (T_MHM*)NULL;
}

/*
   Attempt to set up model from command line.  Command line options are the 
   following

   -di <directory>
      If this argument exists, then all input files are in specified directory.  

   -do <directory>
      If this argument exists, then all output files are in specified directory.

   -ft <filename tag>
      If this argument exists, then the following is attempted:

         1) specification file name:  mhm_final_<tag>.dat
            mhm arguments file name:  mhm_final_<tag>.dat

         2) specification file name:  mhm_intermediate_<tag>.dat
            mhm arguments file name:  mhm_intermediate_<tag>.dat
         
         3) specification file name:  est_final_<tag>.dat
            mhm arguments file name:  -fi <filename>

   -fi <filename>
      If this argument exists, then additional mhm arguments are read from the 
      input file with the given filename.
   
   -fs <filename>
      If this argument exists, then the specification file name is <filename>.  
      The argument -ft takes precedence over -fs.

   -fp <filename>
      If this argument exists, then the posterior is read from <filename>.  Must 
      be used in conjunction with the argument -fs.  The default value is the 
      filename associated with the argument -fs.

   -ph <header>
      If this argument exists, then the header for the posterior file is 
      <header>.  Must be used in conjuction with the arguments -fp or -fs.  The 
      default value is "Posterior mode: ".

   -cm
      If this argument exists, then the mean of the posterior draws are used to
      center the quadratic form. 

   If no command line options are given, then attemps to use a default input file 
   with the name "default.ini".  Returns one valid pointer to a TStateModel upon
   success and null upon failure.
*/
#define LOG_TWO_PI 1.837877066409345
T_MHM* CreateMHM_CommandLine(int nargs, char **args)
{
  TStateModel *model;
  T_MHM *mhm=(T_MHM*)NULL, *rtrn=(T_MHM*)NULL;
  char *d_in, *d_out, *tag, *filename, *spec_filename, *mhm_filename, *id, *fmt;
  FILE *f_in;
  TVector alpha_scales;

  d_in=AddSlash(dw_ParseString_String(nargs,args,"di",""));
  d_out=AddSlash(dw_ParseString_String(nargs,args,"do",""));

  if (filename=dw_ParseString_String(nargs,args,"fi",(char*)NULL))
    {
      fmt="%s%s";
      sprintf(mhm_filename=(char*)malloc(strlen(d_in) + strlen(fmt) + strlen(filename) - 3),fmt,d_in,filename);
      mhm=ReadMHM_Input((FILE*)NULL,mhm_filename,(T_MHM*)NULL);
      mhm->mhm_filename=mhm_filename;
    }

  if (tag=dw_ParseString_String(nargs,args,"ft",(char*)NULL))
    {
      fmt="%smhm_final_%s.dat";
      sprintf(spec_filename=(char*)malloc(strlen(d_in) + strlen(fmt) + strlen(tag) - 3),fmt,d_in,tag);
      if (rtrn=RestartFromFinalFile(spec_filename,mhm))
	mhm=rtrn;
      else
	{
	  free(spec_filename);
	  fmt="%smhm_intermediate_%s.dat";
	  sprintf(spec_filename=(char*)malloc(strlen(d_in) + strlen(fmt) + strlen(tag) - 3),fmt,d_in,tag);
	  if (rtrn=RestartFromIntermediateFile(spec_filename,mhm))
	    mhm=rtrn;
	  else
	    {
	      free(spec_filename);
	      fmt="%sest_final_%s.dat";
	      sprintf(spec_filename=(char*)malloc(strlen(d_in) + strlen(fmt) + strlen(tag) - 3),fmt,d_in,tag);
	      if (!(f_in=fopen(spec_filename,"rt")))
		{
		  fprintf(stderr,"CreateMHM_CommandLine:  Unable to create model from %s tag.\n",tag);
		  if (mhm) FreeMHM(mhm);
		}
	      else
		if (mhm)
		  {		      
		    mhm->parameter_filename=mhm->spec_filename=spec_filename;		      
		    model=Read_VAR_Specification(f_in,(char*)NULL);		  
		    mhm->parameter_header="Posterior mode: ";
		    ReadTransitionMatrices(f_in,(char*)NULL,mhm->parameter_header,model);
		    Read_VAR_Parameters(f_in,(char*)NULL,mhm->parameter_header,model);
		    AddStateModel(model,mhm);
		    fclose(f_in);
		  }		  
	    }
	}
    }
  else 
    if (filename=dw_ParseString_String(nargs,args,"fs",(char*)NULL))
      {
	if (mhm)
	  {
	    fmt="%s%s";
	    sprintf(mhm->spec_filename=(char*)malloc(strlen(d_in) + strlen(fmt) + strlen(filename) - 3),fmt,d_in,filename);
	    model=Read_VAR_Specification((FILE*)NULL,mhm->spec_filename);
	    if (!(filename=dw_ParseString_String(nargs,args,"fp",(char*)NULL)))
	      filename=dw_ParseString_String(nargs,args,"fs",(char*)NULL);
	    sprintf(mhm->parameter_filename=(char*)malloc(strlen(d_in) + strlen(fmt) + strlen(filename) - 3),fmt,d_in,filename);
	    mhm->parameter_header=dw_ParseString_String(nargs,args,"ph","Posterior mode: ");
	    ReadTransitionMatrices((FILE*)NULL,mhm->parameter_filename,mhm->parameter_header,model);
	    Read_VAR_Parameters((FILE*)NULL,mhm->parameter_filename,mhm->parameter_header,model);
	    AddStateModel(model,mhm);
	  }
      }
    else
      {
	fprintf(stderr,"CreateMHM_CommandLine():  No specification file given.\n");
	if (mhm) FreeMHM(mhm);
	exit(0);
      }

  if (!mhm)
    {
      fprintf(stderr,"CreateMHM_CommandLine:  No mhm input data file specified.\n");
      exit(0);
    }

  // Output filenames
  if (!(tag=dw_ParseString_String(nargs,args,"fo",(char*)NULL)))
    tag=dw_ParseString_String(nargs,args,"ft","default");
  fmt="%smhm_intermediate_%s.dat";
  sprintf(mhm->intermediate_output_filename=(char*)malloc(strlen(d_out) + strlen(fmt) + strlen(tag) - 3),fmt,d_out,tag);
  fmt="%smhm_final_%s.dat";
  sprintf(mhm->final_output_filename=(char*)malloc(strlen(d_out) + strlen(fmt) + strlen(tag) - 3),fmt,d_out,tag);
  fmt="%smhm_intermediate_draws_%s.dat";
  sprintf(mhm->intermediate_draws_output_filename=(char*)malloc(strlen(d_out) + strlen(fmt) + strlen(tag) - 3),fmt,d_out,tag);
  fmt="%smhm_draws_%s.dat";
  sprintf(mhm->draws_output_filename=(char*)malloc(strlen(d_out) + strlen(fmt) + strlen(tag) - 3),fmt,d_out,tag);
  fmt="%smhm_regime_counts_%s.dat";
  sprintf(mhm->regime_counts_filename=(char*)malloc(strlen(d_out) + strlen(fmt) + strlen(tag) - 3),fmt,d_out,tag);
  //fmt="%smhm_draws_states_not_integrated_%s.dat";
  //sprintf(mhm->states_not_integrated_out_filename=(char*)malloc(strlen(d_out) + strlen(fmt) + strlen(tag) - 3),fmt,d_out,tag);

  free(d_in);
  free(d_out);
  return mhm;
}
#undef LOG_TWO_PI

int main(int nargs, char **args)
{
  T_MHM *mhm;
  char *header, *buffer[256];
  int initial_time, begin_time, end_time;
  FILE *f_out_intermediate, *f_out_final, *f_out_intermediate_draws;

  if (mhm=CreateMHM_CommandLine(nargs,args))
    {
      //=== Random seed ===//
      dw_initialize_generator(0);

      //=== Test new normalization code ===//
      /**
      TVector** A0;
      PRECISION x1, x2, x3, x4;
      A0=dw_CopyArray((TVector**)NULL,((T_VAR_Parameters*)mhm->model->theta)->A0);

      dw_initialize_generator(-1);

      // burn-in
      initial_time=begin_time=time((time_t*)NULL);
      BurnIn_AdaptiveMetropolisScale(mhm,0,1000);
      end_time=time((time_t*)NULL);

      // test
      while (1)
	{
	  Setup_WZ_Normalization((T_VAR_Parameters*)mhm->model->theta,A0);
	  printf("Likelihood = %lg (WZ normalization)\n",LogLikelihood_StatesIntegratedOut(mhm->model));
	  printf("Prior      = %lg (WZ normalization)\n",x2=LogPrior(mhm->model));
	  printf("Posterior  = %lg (WZ normalization)\n",x4=LogPosterior_StatesIntegratedOut(mhm->model));

	  Setup_No_Normalization((T_VAR_Parameters*)mhm->model->theta);
	  printf("Likelihood = %lg (no normalization)\n",LogLikelihood_StatesIntegratedOut(mhm->model));
	  printf("Prior      = %lg (no normalization)\n",x1=LogPrior(mhm->model));
	  printf("Posterior  = %lg (no normalization)\n",x3=LogPosterior_StatesIntegratedOut(mhm->model));

	  printf("Difference %lg %lg %lg\n\n",x2-x1,x4-x3,((T_VAR_Parameters*)mhm->model->theta)->nvars*log(2));

	  //Setup_WZ_Normalization((T_VAR_Parameters*)mhm->model->theta,A0);
	  DrawAll(mhm->model);

	  getchar();
	}
      /**/
      //=== End test new normalization code ===//

      // Use WZ normalization
      Setup_WZ_Normalization((T_VAR_Parameters*)mhm->model->theta,((T_VAR_Parameters*)mhm->model->theta)->A0);

      // Posterior mode - Initial specification
      f_out_intermediate=dw_AppendTextFile(mhm->intermediate_output_filename);
      fprintf(f_out_intermediate,"//== Initial Specification ==//\n\n");
      Write_VAR_Specification(f_out_intermediate,(char*)NULL,mhm->model);
      ((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies=0;
      Reset_VAR_Improper_Distribution_Counter();
      header="Posterior mode: ";
      WriteTransitionMatrices(f_out_intermediate,(char*)NULL,header,mhm->model);
      Write_VAR_Parameters(f_out_intermediate,(char*)NULL,header,mhm->model);
      fclose(f_out_intermediate);

/*       f_out_final=dw_CreateTextFile(mhm->final_output_filename); */
/*       header="Posterior mode: "; */
/*       WriteTransitionMatrices(f_out_final,(char*)NULL,header,mhm->model); */
/*       Write_VAR_Parameters(f_out_final,(char*)NULL,header,mhm->model); */
/*       fclose(f_out_final); */

      f_out_intermediate_draws=dw_CreateTextFile(mhm->intermediate_draws_output_filename);

      // First burn-in
      if (mhm->n_burn1 > 0)
	{
	  mhm->f_out=f_out_intermediate_draws;
	  initial_time=begin_time=time((time_t*)NULL);
	  BurnIn_AdaptiveMetropolisScale(mhm,mhm->n_burn1,1000);
	  end_time=time((time_t*)NULL);
	  printf("Elapsed Time: %d seconds\n",end_time - begin_time);
	}

      // After first burn-in
      f_out_intermediate=dw_AppendTextFile(mhm->intermediate_output_filename);
      fprintf(f_out_intermediate,"//== Specification after first burn-in ==//\n");
      fprintf(f_out_intermediate,"Number inconsistent normalizations: %d\n",((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies);
      fprintf(f_out_intermediate,"Number singular inverse variances: %d\n\n",Get_VAR_Improper_Distribution_Counter());
      ((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies=0;
      header="First burn-in: ";
      WriteTransitionMatrices(f_out_intermediate,(char*)NULL,header,mhm->model);
      Write_VAR_Parameters(f_out_intermediate,(char*)NULL,header,mhm->model);
      fclose(f_out_intermediate);

      // Second burn-in
      if (mhm->n_burn2 > 0)
	{
	  mhm->f_out=f_out_intermediate_draws;
	  initial_time=begin_time=time((time_t*)NULL);
	  BurnIn(mhm,mhm->n_burn2,1000);
	  end_time=time((time_t*)NULL);
	  printf("Elapsed Time: %d seconds\n",end_time - begin_time);

	  printf("Number inconsistent normalizations: %d\n",((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies);
	  printf("Number singular inverse variances: %d\n\n",Get_VAR_Improper_Distribution_Counter());
	}

      fclose(f_out_intermediate_draws);

      // After second burn-in
      f_out_intermediate=dw_AppendTextFile(mhm->intermediate_output_filename);
      fprintf(f_out_intermediate,"//== Specification after second burn-in ==//\n");
      fprintf(f_out_intermediate,"Number inconsistent normalizations: %d\n",((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies);
      fprintf(f_out_intermediate,"Number singular inverse variances: %d\n\n",Get_VAR_Improper_Distribution_Counter());
      ((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies=0;
      header="Second burn-in: ";
      WriteTransitionMatrices(f_out_intermediate,(char*)NULL,header,mhm->model);
      Write_VAR_Parameters(f_out_intermediate,(char*)NULL,header,mhm->model);
      fclose(f_out_intermediate);
  
      // Mean-variance estimation
      if (mhm->n_mean_variance > 0)
	{
	  begin_time=time((time_t*)NULL);
	  ComputeMeanVariance_MHM(mhm,mhm->n_mean_variance,10000);
	  end_time=time((time_t*)NULL);
	  printf("Elapsed Time: %d seconds\n",end_time - begin_time);

	  printf("Number inconsistent normalizations: %d\n",((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies);
	  printf("Number singular inverse variances: %d\n\n",Get_VAR_Improper_Distribution_Counter());
	}

      // Set center to mean if necessary
      if (dw_FindArgument_String(nargs,args,"cm") >= 0)
	{
	  printf("Using mean for center\n");
	  mhm->center=mhm->mean;
	}
      else
	printf("Using posterior mode for center\n");


      // After mean-variance estimation
      f_out_intermediate=dw_AppendTextFile(mhm->intermediate_output_filename);
      fprintf(f_out_intermediate,"//== Specification after mean-variance estimation ==//\n");
      fprintf(f_out_intermediate,"Number inconsistent normalizations: %d\n",((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies);
      fprintf(f_out_intermediate,"Number singular inverse variances: %d\n\n",Get_VAR_Improper_Distribution_Counter());
      ((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies=0;
      header="Mean-variance: ";
      WriteTransitionMatrices(f_out_intermediate,(char*)NULL,header,mhm->model);
      Write_VAR_Parameters(f_out_intermediate,(char*)NULL,header,mhm->model);
      WriteMeanVariance(f_out_intermediate,mhm);
      fclose(f_out_intermediate);

      // Open draw file and states file
      mhm->f_out=dw_CreateTextFile(mhm->draws_output_filename);
      WriteMHM_Input(mhm->f_out,mhm);
      WriteMeanVariance(mhm->f_out,mhm);
      //mhm->f_states_not_integrated_out=dw_CreateTextFile(mhm->states_not_integrated_out_filename);
      //WriteMHM_Input(mhm->f_states_not_integrated_out,mhm);
      //WriteMeanVariance(mhm->f_states_not_integrated_out,mhm);
      mhm->f_out_regime_counts=dw_CreateTextFile(mhm->regime_counts_filename);

      // Modified harmonic mean draws 
      fprintf(mhm->f_out,"\n//== Draws ==//\n");
      //fprintf(mhm->f_states_not_integrated_out,"\n//== Draws ==//\n");

      begin_time=time((time_t*)NULL);
      ComputeModifiedHarmonicMean(mhm,10000);
      end_time=time((time_t*)NULL);
      printf("Elapsed Time: %d seconds\n",end_time - begin_time);
      printf("Number inconsistent normalizations: %d\n",((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies);
      printf("Number singular inverse variances: %d\n\n",Get_VAR_Improper_Distribution_Counter());

      fclose(mhm->f_out);

      // After modified harmonic mean draws
      f_out_intermediate=dw_AppendTextFile(mhm->intermediate_output_filename);
      fprintf(f_out_intermediate,"//== Specification after mhm draws ==//\n");
      fprintf(f_out_intermediate,"Number inconsistent normalizations: %d\n",((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies);
      fprintf(f_out_intermediate,"Number singular inverse variances: %d\n\n",Get_VAR_Improper_Distribution_Counter());
      fprintf(f_out_intermediate,"//== RNG State ==//\n");
      dw_print_generator_state(f_out_intermediate);
      fprintf(f_out_intermediate,"\n");
      ((T_VAR_Parameters*)mhm->model->theta)->WZ_inconsistancies=0;
      header="Final draw: ";
      WriteTransitionMatrices(f_out_intermediate,(char*)NULL,header,mhm->model);
      Write_VAR_Parameters(f_out_intermediate,(char*)NULL,header,mhm->model);
      fclose(f_out_intermediate);
    }

  return 0;
}
