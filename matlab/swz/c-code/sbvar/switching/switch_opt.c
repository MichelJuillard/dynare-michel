
#include "switch_opt.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

//====== Static Global Variables ======
static TStateModel *Model=(TStateModel*)NULL;
static PRECISION *buffer=(PRECISION*)NULL;
static PRECISION *ModifiedFreeParameters=(PRECISION*)NULL;
static PRECISION *FreeParameters_Q=(PRECISION*)NULL;
static int NumberFreeParameters_Q=0;
static PRECISION *FreeParameters_Theta=(PRECISION*)NULL;
static int NumberFreeParameters_Theta=0;


void SetupObjectiveFunction(TStateModel *model, PRECISION *Modified, PRECISION *FreeQ, PRECISION *FreeTheta)
{
  if (buffer) free(buffer);
  Model=model;
  FreeParameters_Q=FreeQ;
  NumberFreeParameters_Q=NumberFreeParametersQ(model);
  FreeParameters_Theta=FreeTheta;
  NumberFreeParameters_Theta=model->routines->pNumberFreeParametersTheta(model);
  ModifiedFreeParameters=Modified;
}

void SetupObjectiveFunction_new(TStateModel *model, int FreeTheta_Idx, int FreeQ_Idx, int Modified_Idx)
{
  if (buffer) free(buffer);
  Model=model;
  NumberFreeParameters_Q=NumberFreeParametersQ(model);
  NumberFreeParameters_Theta=model->routines->pNumberFreeParametersTheta(model);
  buffer=(PRECISION*)malloc((NumberFreeParameters_Q + NumberFreeParameters_Theta)*sizeof(PRECISION));

  FreeParameters_Q=buffer+FreeQ_Idx;
  FreeParameters_Theta=buffer+FreeTheta_Idx;
  ModifiedFreeParameters=buffer+Modified_Idx;
}

PRECISION PosteriorObjectiveFunction(PRECISION *x, int n)
{
  if (x != ModifiedFreeParameters) memmove(ModifiedFreeParameters,x,n*sizeof(PRECISION));
  ConvertFreeParametersToQ(Model,FreeParameters_Q);
  ConvertFreeParametersToTheta(Model,FreeParameters_Theta);
  return -LogPosterior_StatesIntegratedOut(Model);

  //PRECISION lp_Q, lp_Theta, li;
  //FILE *f_out;
  //lp_Q=LogPrior_Q(Model);
  //lp_Theta=LogPrior_Theta(Model);
  //li=LogLikelihood_StatesIntegratedOut(Model);
  //if (isnan(lp_Q) || isnan(lp_Theta) || isnan(li))
  //  {
  //    f_out=fopen("tmp.tmp","wt");
  //    Write_VAR_Specification(f_out,(char*)NULL,Model);
  //    WriteTransitionMatrices(f_out,(char*)NULL,"Error: ",Model);
  //    Write_VAR_Parameters(f_out,(char*)NULL,"Error: ",Model);
  //    fprintf(f_out,"LogPrior_Theta(): %le\n",lp_Theta);
  //    fprintf(f_out,"LogPrior_Q(): %le\n",lp_Q);
  //    fprintf(f_out,"LogLikelihood_StatesIntegratedOut(): %le\n",li);
  //    fprintf(f_out,"Posterior: %le\n\n",lp_Q+lp_Theta+li);
  //    fclose(f_out);
  //    exit(0);
  //  }
  //return -(lp_Q+lp_Theta+li);
}

PRECISION PosteriorObjectiveFunction_csminwel(double *x, int n, double **args, int *dims)
{
  return PosteriorObjectiveFunction(x,n);
}

void PosteriorObjectiveFunction_npsol(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{
  *f=PosteriorObjectiveFunction(x,*n);
}

PRECISION MLEObjectiveFunction(PRECISION *x, int n)
{
  if (x != ModifiedFreeParameters) memmove(ModifiedFreeParameters,x,n*sizeof(PRECISION));
  ConvertFreeParametersToQ(Model,FreeParameters_Q);
  ConvertFreeParametersToTheta(Model,FreeParameters_Theta);
  return -LogLikelihood_StatesIntegratedOut(Model);
}

PRECISION MLEObjectiveFunction_csminwel(double *x, int n, double **args, int *dims)
{
  return MLEObjectiveFunction(x,n);
}

void MLEObjectiveFunction_npsol(int *mode, int *n, double *x, double *f, double *g, int *nstate)
{
  *f=MLEObjectiveFunction(x,*n);
}


PRECISION MLEObjectiveFunction_LogQ(PRECISION *x, int n)
{
  if (x != ModifiedFreeParameters) memmove(ModifiedFreeParameters,x,n*sizeof(PRECISION));
  ConvertLogFreeParametersToQ(Model,FreeParameters_Q);
  ConvertFreeParametersToTheta(Model,FreeParameters_Theta);
  return -LogLikelihood_StatesIntegratedOut(Model);
}

PRECISION MLEObjectiveFunction_LogQ_csminwel(double *x, int n, double **args, int *dims)
{
  return MLEObjectiveFunction_LogQ(x,n);
}





