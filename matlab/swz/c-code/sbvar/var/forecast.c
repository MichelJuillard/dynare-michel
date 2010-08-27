
#include "switch.h"
#include "switchio.h"
#include "VARio.h"
#include "dw_parse_cmd.h"
#include "dw_ascii.h"
#include "dw_histogram.h"

#include <stdlib.h>
#include <string.h>

/*
   Assumes
    f_out : valid FILE pointer
    percentiles : vector of numbers between 0 and 1 inclusive
    draws : number of draws of shocks and regimes to make for each posterior draw
    posterior_file : FILE pointer to file containing posterior draws.  If null, current parameters are used.
    T : last observation to treat as data.  Usually equals model->nobs.
    h : non-negative integer
    model : point to valid TStateModel structure

   Results:
    Computes and prints to the file f_out the requested percentiles for forecasts 
    of the observables.

   Returns:
    One upon success and zero otherwise.

*/
int forecast_percentile(FILE *f_out, TVector percentiles, int draws, FILE *posterior_file, int T, int h, TStateModel *model)
{
  T_VAR_Parameters *p;
  int done=0, rtrn=0, *S, i=0, j, k, m, n=1000;
  TVector init_prob, prob, *shocks, initial;
  TMatrix forecast;
  TMatrixHistogram *histogram;

 // quick check of passed parameters
  if (!f_out || !percentiles || (draws <= 0) || (T < 0) || (h < 0) || !model) return 0;

  p=(T_VAR_Parameters*)(model->theta);

  if (T > p->nobs) return 0;

  // allocate memory
  S=(int*)malloc(h*sizeof(int));
  forecast=CreateMatrix(h,p->nvars);
  histogram=CreateMatrixHistogram(h,p->nvars,100,HISTOGRAM_VARIABLE);
  initial=CreateVector(p->npre);
  shocks=dw_CreateArray_vector(h);
  for (i=h-1; i >= 0; i--) shocks[i]=CreateVector(p->nvars);
  init_prob=CreateVector(p->nstates);
  prob=CreateVector(p->nstates);

  // Initial value
  EquateVector(initial,p->X[T]);

  i=0;
  while (!done)
    {
      // Read parameters and push them into model
      if (!posterior_file)
	done=1;
      else
	if (!ReadBaseTransitionMatricesFlat(posterior_file,model) || !Read_VAR_ParametersFlat(posterior_file,model))
	  {
	    done=2;
	    printf("total posterior draws processed - %d\n",i);
	  }
	else
	  if (i++ == n)
	    {
	      printf("%d posterior draws processed\n",i);
	      n+=1000;
	    }

	if (done != 2)
	  {
	    // Get filtered probability at time T
	    for (j=p->nstates-1; j >= 0; j--)
	      ElementV(init_prob,j)=ProbabilityStateConditionalCurrent(j,T,model);

	    for (k=draws; k > 0; k--)
	      {
		// Draw time T regime
		m=DrawDiscrete(init_prob);
              
		// Draw regimes from time T+1 through T+h inclusive
		for (j=0; j < h; j++)
		  {
		    ColumnVector(prob,model->sv->Q,m);
		    S[j]=m=DrawDiscrete(prob);
		  }

		// Draw shocks
		for (j=h-1; j >= 0; j--) dw_NormalVector(shocks[j]); // InitializeVector(shocks[i],0.0);

		// Compute forecast
		if (!forecast_base(forecast,h,initial,shocks,S,model))
		  goto ERROR_EXIT;

		// Accumulate impulse response
		AddMatrixObservation(forecast,histogram);
	      }
	  }
    }

  for (i=0; i < DimV(percentiles); i++)
    {
      MatrixPercentile(forecast,ElementV(percentiles,i),histogram);
      dw_PrintMatrix(f_out,forecast,"%lg ");
      fprintf(f_out,"\n");
    }

  rtrn=1;

ERROR_EXIT:
  FreeMatrixHistogram(histogram);
  FreeMatrix(forecast);
  free(S);
  FreeVector(initial);
  FreeVector(prob);
  FreeVector(init_prob);
  dw_FreeArray(shocks);

  return rtrn;
}

/*
   Assumes
    f_out : valid FILE pointer
    percentiles : vector of numbers between 0 and 1 inclusive
    draws : number of draws of shocks to make for each posterior draw
    posterior_file : FILE pointer to file containing posterior draws.  If null, current parameters are used.
    s : base state
    T : last observation to treat as data.  Usually equals model->nobs.
    h : non-negative integer
    model : point to valid TStateModel/T_MSStateSpace structure

   Results:
    Computes and prints to the file f_out the requested percentiles for forecasts 
    of the observables.

   Returns:
    One upon success and zero otherwise.

   Notes:
    The regime at time T is drawn from the filtered probabilities at time t, and
    is set to s there after. 

*/
int forecast_percentile_regime(FILE *f_out, TVector percentiles, int draws, 
			       FILE *posterior_file, int s, int T, int h, TStateModel *model)
{
  T_VAR_Parameters *p;
  int done=0, rtrn=0, *S, i=0, j, k, m, n=1000;
  TVector init_prob, prob, *shocks, initial;
  TMatrix forecast;
  TMatrixHistogram *histogram;

 // quick check of passed parameters
  if (!f_out || !percentiles || (draws <= 0) || (T < 0) || (h < 0) || !model) return 0;

  p=(T_VAR_Parameters*)(model->theta);

  if (T > p->nobs) return 0;

  // allocate memory
  S=(int*)malloc(h*sizeof(int));
  for (i=0; i < h; i++) S[i]=s;
  forecast=CreateMatrix(h,p->nvars);
  histogram=CreateMatrixHistogram(h,p->nvars,100,HISTOGRAM_VARIABLE);
  initial=CreateVector(p->npre);
  shocks=dw_CreateArray_vector(h);
  for (i=h-1; i >= 0; i--) shocks[i]=CreateVector(p->nvars);
  init_prob=CreateVector(p->nstates);
  prob=CreateVector(p->nstates);

  // Initial value
  EquateVector(initial,p->X[T]);

  i=0;
  while (!done)
    {
      // Read parameters and push them into model
      if (!posterior_file)
	done=1;
      else
	if (!ReadBaseTransitionMatricesFlat(posterior_file,model) || !Read_VAR_ParametersFlat(posterior_file,model))
	  {
	    done=2;
	    printf("total posterior draws processed - %d\n",i);
	  }
	else
	  if (i++ == n)
	    {
	      printf("%d posterior draws processed\n",i);
	      n+=1000;
	    }

      if (done != 2)
	{
	  for (k=draws; k > 0; k--)
	    {
	      // Draw shocks
	      for (j=h-1; j >= 0; j--) dw_NormalVector(shocks[j]); // InitializeVector(shocks[i],0.0);

	      // Compute forecast
	      if (!forecast_base(forecast,h,initial,shocks,S,model))
		goto ERROR_EXIT;

	      // Accumulate impulse response
	      AddMatrixObservation(forecast,histogram);
	    }
	}
    }

  for (i=0; i < DimV(percentiles); i++)
    {
      MatrixPercentile(forecast,ElementV(percentiles,i),histogram);
      dw_PrintMatrix(f_out,forecast,"%lg ");
      fprintf(f_out,"\n");
    }

  rtrn=1;

ERROR_EXIT:
  FreeMatrixHistogram(histogram);
  FreeMatrix(forecast);
  free(S);
  FreeVector(initial);
  FreeVector(prob);
  FreeVector(init_prob);
  dw_FreeArray(shocks);

  return rtrn;
}


/*
   Attempt to set up model from command line.  Command line options are the 
   following

   -ft <filename tag>
      If this argument exists, then the following is attempted:
         specification file name = est_final_<tag>.dat
         output file name        = ir_<tag>_regime_<k>.dat
         parameters file name    = est_final_<tag>.dat 
         header                  = "Posterior mode: "
   
   -fs <filename>
      If this argument exists, then the specification file name is <filename>.  
      The argument -fs takes precedence over -ft.

   -fp <filename>
      If this argument exists, then the parameters file name is <filename>.  The 
      argument -fp takes precedence over -ft.  The default value is the filename 
      associated with the argument -fs.

   -ph <header>
      If this argument exists, then the header for the parameters file is 
      <header>.  The default value is "Posterior mode: ".

   -horizon <integer>
      If this argument exists, then the horizon of the impulse responses is given
      by the passed integer.  The default value is 12.

   -error_bands 
      Output error bands.  (default = off - only median is computed)

   -percentiles n p_1 p_2 ... p_n
      Percentiles to compute. The first parameter after percentiles must be the 
      number of percentiles and the following values are the actual percentiles. 
      default = 3  0.16  0.50  0.84   if error_bands flag is set
              = 1  0.50               otherwise

   -parameter_uncertainty 
      Apply parameter uncertainty when computing error bands.

   -shocks_per_parameter <integer> 
      Number of shocks and regime paths to draw for each parameter draw.  The 
      default value is 1 if parameter_uncertainty is set and 10,000 otherwise.
     
   -thin 
      Thinning factor.  Only 1/thin of the draws in posterior draws file are 
      used. The default value is 1.

   -regimes 
      Produces forecasts as if each regime were permanent. (default = off)

   -regime <integer>
      Produces forecasts as if regime were permanent.

   -mean
      Produces mean forecast.  (default = off)

*/
int main(int nargs, char **args)
{
  char *spec=(char*)NULL, *parm=(char*)NULL, *head=(char*)NULL, *post=(char*)NULL, *out_filename, *tag, *buffer, *fmt;
  TStateModel *model;
  T_VAR_Parameters *p;
  TVector percentiles=(TVector)NULL;
  int s, horizon, thin, draws, i, j, n;
  FILE *f_out, *posterior_file;

  // specification filename
  if (buffer=dw_ParseString_String(nargs,args,"fs",(char*)NULL))
    strcpy(spec=(char*)malloc(strlen(buffer)+1),buffer);

  // parameter filename
  if (buffer=dw_ParseString_String(nargs,args,"fp",(char*)NULL))
    strcpy(parm=(char*)malloc(strlen(buffer)+1),buffer);

  // header
  if (buffer=dw_ParseString_String(nargs,args,"ph",(char*)NULL))
    strcpy(head=(char*)malloc(strlen(buffer)+1),buffer);

  // file tag
  if (tag=dw_ParseString_String(nargs,args,"ft",(char*)NULL))
    {
      fmt="est_final_%s.dat";

      // specification filename
      if (!spec)
	sprintf(spec=(char*)malloc(strlen(fmt) + strlen(tag) - 1),fmt,tag);

      // parameter filename
      if (!parm)
	sprintf(parm=(char*)malloc(strlen(fmt) + strlen(tag) - 1),fmt,tag);
    } 

  // horizon
  horizon=dw_ParseInteger_String(nargs,args,"horizon",12);

  if (!spec)
    {
      fprintf(stderr,"No specification filename given\n");
      fprintf(stderr,"Command line syntax:\n"
                     "  -ft : file tag\n"
	             "  -fs : specification filename (est_final_<tag>.dat)\n"
	             "  -fp : parameters filename (specification filename)\n"
                     "  -fh : parameter header (Posterior mode: )\n"
                     "  -horizon : horizon for the forecast (12)\n"
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

  model=Read_VAR_Specification((FILE*)NULL,spec);
  ReadTransitionMatrices((FILE*)NULL,parm,head,model);
  Read_VAR_Parameters((FILE*)NULL,parm,head,model);
  p=(T_VAR_Parameters*)(model->theta);

  free(spec);
  free(head);
  free(parm);

  //============================= Compute forecasts ============================= 

  // Mean forecast
  /* if (dw_FindArgument_String(nargs,args,"mean") != -1) */
  /*   { */
  /*     fmt="forecasts_mean_%s.prn"; */
  /*     sprintf(out_filename=(char*)malloc(strlen(fmt) + strlen(tag) - 1),fmt,tag); */
  /*     f_out=fopen(out_filename,"wt"); */
  /*     free(out_filename); */
  /*     printf("Constructing mean forecast\n"); */
  /*     if (F=dw_state_space_mean_unconditional_forecast((TVector*)NULL,h,statespace->nobs,model)) */
  /* 	for (i=0; i < h; i++) */
  /* 	  dw_PrintVector(f_out,F[i],"%le "); */
  /*     fclose(f_out); */
  /*     return; */
  /*   } */

  // Parameter uncertainty
  if (dw_FindArgument_String(nargs,args,"parameter_uncertainity") != -1)
    {
      // Open posterior draws file
      fmt="draws_%s.dat";
      sprintf(post=(char*)malloc(strlen(fmt) + strlen(tag) - 1),fmt,tag);
      if (!(posterior_file=fopen(post,"rt")))
	{
	  printf("Unable to open draws file: %s\n",post);
	  exit(0);
	}

      // Get thinning factor from command line
      thin=dw_ParseInteger_String(nargs,args,"thin",1);

      // Get shocks_per_parameter from command line
      draws=dw_ParseInteger_String(nargs,args,"shocks_per_parameter",1);
    }
  else
    {
      // Using posterior estimate
      posterior_file=(FILE*)NULL;

      // thinning factor not used
      thin=1;

      // Get shocks_per_parameter from command line
      draws=dw_ParseInteger_String(nargs,args,"shocks_per_parameter",10000);
    }

  // Setup percentiles
  if ((i=dw_FindArgument_String(nargs,args,"percentiles")) == -1)
    if (dw_FindArgument_String(nargs,args,"error_bands") == -1)
      {
	percentiles=CreateVector(1);
	ElementV(percentiles,0)=0.5;
      }
    else
      {
	percentiles=CreateVector(3);
	ElementV(percentiles,0)=0.16; ElementV(percentiles,1)=0.5; ElementV(percentiles,2)=0.84;
      }
  else
    if ((i+1 < nargs) && dw_IsInteger(args[i+1]) && ((n=atoi(args[i+1])) > 0) && (i+1+n < nargs))
      {
	percentiles=CreateVector(n);
	for (j=0; j < n; j++)
	  if (!dw_IsFloat(args[i+2+j])|| ((ElementV(percentiles,j)=atof(args[i+2+j])) <= 0.0)
	      || (ElementV(percentiles,j) >= 1.0)) break;
	if (j < n)
	  {
	    FreeVector(percentiles);
	    printf("forecasting command line: Error parsing percentiles\n");
	    exit(0);
	  }
      }
    else
      {
	printf("forecasting command line(): Error parsing percentiles\n");
	exit(0);
      }

  if (dw_FindArgument_String(nargs,args,"regimes") != -1)
    for (s=0; s < p->nstates; s++)
      {
	rewind(posterior_file);
	fmt="forecasts_percentiles_regime_%d_%s.prn";
	sprintf(out_filename=(char*)malloc(strlen(fmt) + strlen(tag) - 3),fmt,s,tag);
	f_out=fopen(out_filename,"wt");
	free(out_filename);
	printf("Constructing percentiles for forecasts - regime %d\n",s);
	forecast_percentile_regime(f_out,percentiles,draws,posterior_file,s,p->nobs,horizon,model);
	fclose(f_out);
      }
  else
    if (((s=dw_ParseInteger_String(nargs,args,"regime",-1)) >= 0) && (s < p->nstates))
      {
	fmt="forecasts_percentiles_regime_%d_%s.prn";
	sprintf(out_filename=(char*)malloc(strlen(fmt) + strlen(tag) - 3),fmt,s,tag);
	f_out=fopen(out_filename,"wt");
	free(out_filename);
	printf("Constructing percentiles for forecasts - regime %d\n",s);
	forecast_percentile_regime(f_out,percentiles,draws,posterior_file,s,p->nobs,horizon,model);
	fclose(f_out);
      }
    else
      {
	fmt="forecasts_percentiles_%s.prn";
	sprintf(out_filename=(char*)malloc(strlen(fmt) + strlen(tag) - 1),fmt,tag);
	f_out=fopen(out_filename,"wt");
	free(out_filename);
	printf("Constructing percentiles for forecasts - %d draws of shocks/regimes per posterior value\n",draws);
	forecast_percentile(f_out,percentiles,draws,posterior_file,p->nobs,horizon,model);
	fclose(f_out);
      }

  if (posterior_file) fclose(posterior_file);
  FreeVector(percentiles);
  //=============================================================================

  return 0;
}

