
#include "mhm_VAR.h"
#include "VARio.h"
#include "switch.h"
#include "switchio.h"
#include "dw_ascii.h"
#include "dw_rand.h"
#include "dw_matrix_rand.h"
#include "dw_error.h"

#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "modify_for_mex.h"

  // Compute psudo-inverse of mhm->variance
static void PsudoInverse(TMatrix X, TMatrix Y)
{
  int i, j, k;
  TMatrix U, V;
  TVector d;
  PRECISION epsilon, tmp;
  k=RowM(Y);
  SVD(U=CreateMatrix(k,k),d=CreateVector(k),V=CreateMatrix(k,k),Y);
  for (epsilon=ElementV(d,0),i=k-1; i > 0; i--)
    if (ElementV(d,0) > epsilon) epsilon=ElementV(d,0);
  epsilon*=SQRT_MACHINE_EPSILON;
  for (j=k-1; j >= 0; j--)
    {
      tmp=(ElementV(d,j) > epsilon) ? 1.0/ElementV(d,j) : 0.0;
      for (i=k-1; i >= 0; i--)
    ElementM(V,i,j)*=tmp;
    }
  ProductTransposeMM(X,V,U);
  FreeMatrix(U);
  FreeVector(d);
  FreeMatrix(V);
}

void ResetMHM(T_MHM *mhm)
{
  mhm->N=0;

  mhm->sum=mhm->sum_square=0.0;

  mhm->max_log_posterior=mhm->max_log_likelihood=MINUS_INFINITY;
}

void FreeMHM(T_MHM *mhm)
{
  if (mhm)
    {
      FreeVector(mhm->mean);
      FreeVector(mhm->posterior_mode_VAR);
      FreeMatrix(mhm->variance);
      FreeMatrix(mhm->inverse_variance);

      FreeVector(mhm->free_parameters_VAR);

      dw_FreeArray(mhm->states);

      dw_FreeArray(mhm->BaseAlpha);
      dw_FreeArray(mhm->Alpha);

      free(mhm);
    }
}

T_MHM* AddStateModel(TStateModel *model, T_MHM *mhm)
{
  int i, nf_var=NumberFreeParametersTheta(model);

  if (!mhm) mhm=CreateMHM();

  // Allocate memory
  mhm->mean=CreateVector(nf_var);
  mhm->posterior_mode_VAR=CreateVector(nf_var);
  mhm->variance=CreateMatrix(nf_var,nf_var);
  mhm->inverse_variance=CreateMatrix(nf_var,nf_var);

  mhm->free_parameters_VAR=CreateVector(nf_var);

  mhm->states=dw_CreateArray_int(model->sv->nstates);

  mhm->BaseAlpha=dw_CreateArray_vector(dw_DimA(model->sv->ba));
  for (i=dw_DimA(model->sv->ba)-1; i >= 0; i--)
    mhm->BaseAlpha[i]=CreateVector(DimV(model->sv->ba[i]));

  // model information
  mhm->model=model;
  Setup_WZ_Normalization((T_VAR_Parameters*)mhm->model->theta,((T_VAR_Parameters*)mhm->model->theta)->A0);
  ConvertThetaToFreeParameters(model,pElementV(mhm->posterior_mode_VAR));
  mhm->log_likelihood_at_mode=LogLikelihood_StatesIntegratedOut(model);
  mhm->log_prior_at_mode=LogPrior(model);
  mhm->log_posterior_at_mode=mhm->log_likelihood_at_mode + mhm->log_prior_at_mode;

  // Center
  mhm->center=mhm->posterior_mode_VAR;

  return mhm;
}

T_MHM* AddDirichletScales(TVector alpha_scales, T_MHM *mhm)
{
  if (!mhm) mhm=CreateMHM();

  mhm->alpha_scales=EquateVector((TVector)NULL,alpha_scales);

  return mhm;
}

T_MHM* CreateMHM(void)
{
  int i, j;
  T_MHM* mhm;

  // Allocate structure
  mhm=(T_MHM*)malloc(sizeof(T_MHM));

  mhm->alpha_scales=(TVector)NULL;

  mhm->mean=(TVector)NULL;
  mhm->posterior_mode_VAR=(TVector)NULL;
  mhm->variance=(TMatrix)NULL;
  mhm->inverse_variance=(TMatrix)NULL;
  mhm->free_parameters_VAR=(TVector)NULL;

  mhm->BaseAlpha=(TVector*)NULL;
  mhm->Alpha=(TVector**)NULL;

  mhm->model=(TStateModel*)NULL;

  mhm->f_out=(FILE*)NULL;
  mhm->intermediate_output_filename=(char*)NULL;
  mhm->final_output_filename=(char*)NULL;
  mhm->intermediate_draws_output_filename=(char*)NULL;
  mhm->draws_output_filename=(char*)NULL;
  mhm->spec_filename=(char*)NULL;
  mhm->parameter_filename=(char*)NULL;
  mhm->parameter_header=(char*)NULL;
  mhm->mhm_filename=(char*)NULL;

  // Default values
  mhm->n_burn1=100000;
  mhm->n_burn2=0;
  mhm->n_mean_variance=200000;
  mhm->n_mhm=1000000;
  mhm->n_thin=1;

  ResetMHM(mhm);

  return mhm;
}

void BurnIn(T_MHM *mhm, int iterations, int period)
{
  int count, begin_time, check=period;
  printf("Beginning burn in -- %d iterations.\n",iterations);
  begin_time=time((time_t*)NULL);
  for (count=1; count <= iterations; count++)
    {
      DrawAll(mhm->model);

      if (count == check)
    {
      check+=period;
      if (mhm->f_out)
        {
          fprintf(mhm->f_out,"%d iterations completed out of %d\n",count,iterations);
          PrintJumps(mhm->f_out,(T_VAR_Parameters*)(mhm->model->theta));
          fflush(mhm->f_out);
        }

      printf("Total Elapsed Time: %d seconds\n",(int)time((time_t*)NULL) - begin_time);
      printf("%d iterations completed out of %d\n",count,iterations);
      PrintJumps(stdout,(T_VAR_Parameters*)(mhm->model->theta));
    }
    }
  ResetMetropolisInformation((T_VAR_Parameters*)(mhm->model->theta));
}

void BurnIn_AdaptiveMetropolisScale(T_MHM *mhm, int iterations, int period)
{
  int verbose=1;
  AdaptiveMetropolisScale(mhm->model,iterations,period,verbose,mhm->f_out);
}

/*
   Computes mean and variance and base alpha
*/
void ComputeMeanVariance_MHM(T_MHM *mhm, int iterations, int period)
{
  int i, j, begin_time, count, check=period;
  TVector *alpha;
  TMatrix S;
  PRECISION max, inc, tmp;

  dw_InitializeArray_vector(mhm->BaseAlpha,0.0);
  alpha=dw_CopyArray((void*)NULL,mhm->BaseAlpha);

  S=CreateMatrix(RowM(mhm->variance),ColM(mhm->variance));
  InitializeVector(mhm->mean,0.0);
  InitializeMatrix(mhm->variance,0.0);

  // loop and accumulate 1st and 2nd non-central moments
  printf("Beginning mean and variance estimation -- %d iterations.\n",iterations);
  begin_time=time((time_t*)NULL);
  for (count=1; count <= iterations; count++)
    {
      DrawAll(mhm->model);

      ConvertThetaToFreeParameters(mhm->model,pElementV(mhm->free_parameters_VAR));

      for (i=dw_DimA(alpha)-1; i >= 0; i--)
    {
      AddVV(mhm->BaseAlpha[i],mhm->BaseAlpha[i],mhm->model->sv->ba[i]);

      for (j=DimV(alpha[i])-1; j >= 0; j--)
        ElementV(alpha[i],j)+=ElementV(mhm->model->sv->ba[i],j)*ElementV(mhm->model->sv->ba[i],j);
    }

      AddVV(mhm->mean,mhm->mean,mhm->free_parameters_VAR);
      OuterProduct(S,mhm->free_parameters_VAR,mhm->free_parameters_VAR);
      AddMM(mhm->variance,mhm->variance,S);

      if (count == check)
    {
      check+=period;
      printf("Total Elapsed Time: %d seconds\n",(int)time((time_t*)NULL) - begin_time);
      printf("%d iterations completed out of %d\n",count,iterations);
      PrintJumps(stdout,(T_VAR_Parameters*)(mhm->model->theta));
    }
    }

  // compute 1st and 2nd central moments for normal terms
  ProductVS(mhm->mean,mhm->mean,1.0/(PRECISION)iterations);
  ProductMS(mhm->variance,mhm->variance,1.0/(PRECISION)iterations);
  OuterProduct(S,mhm->mean,mhm->mean);
  SubtractMM(mhm->variance,mhm->variance,S);

  // Psudo variance
  SubtractVV(mhm->free_parameters_VAR,mhm->mean,mhm->posterior_mode_VAR);
  OuterProduct(S,mhm->free_parameters_VAR,mhm->free_parameters_VAR);
  AddMM(mhm->variance,mhm->variance,S);

  // Compute psudo-inverse of mhm->variance
  PsudoInverse(mhm->inverse_variance,mhm->variance);

  FreeMatrix(S);

  // compute base alpha's for Dirichlet distribution
  for (i=dw_DimA(mhm->BaseAlpha)-1; i >= 0; i--)
    ProductVS(mhm->BaseAlpha[i],mhm->BaseAlpha[i],1.0/(PRECISION)iterations);
  for (i=dw_DimA(mhm->BaseAlpha)-1; i >= 0; i--)
    for (j=DimV(mhm->BaseAlpha[i])-1; j >= 0; j--)
      ElementV(alpha[i],j)=ElementV(alpha[i],j)/(PRECISION)iterations - ElementV(mhm->BaseAlpha[i],j)*ElementV(mhm->BaseAlpha[i],j);

  for (i=dw_DimA(mhm->BaseAlpha)-1; i >= 0; i--)
    {
      for (max=0.0, j=DimV(mhm->BaseAlpha[i])-1; j >= 0; j--)
    if ((tmp=ElementV(mhm->BaseAlpha[i],j)*(1.0-ElementV(mhm->BaseAlpha[i],j))/ElementV(alpha[i],j)) > max) max=tmp;
      max-=1.0;

      for (inc=0.0, j=DimV(mhm->BaseAlpha[i])-1; j >= 0; j--)
    if ((tmp=1.1 - max*ElementV(mhm->BaseAlpha[i],j)) > inc) inc=tmp;

      for (j=DimV(mhm->BaseAlpha[i])-1; j >= 0; j--)
    ElementV(mhm->BaseAlpha[i],j)=max*ElementV(mhm->BaseAlpha[i],j)+inc;
    }

  // Create Alpha's
  mhm->Alpha=dw_CreateArray_array(DimV(mhm->alpha_scales));
  for (i=dw_DimA(mhm->Alpha)-1; i >= 0; i--)
    {
      mhm->Alpha[i]=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
      for (j=dw_DimA(mhm->Alpha[i])-1; j >= 0; j--)
    mhm->Alpha[i][j]=ProductVS((TVector)NULL,mhm->BaseAlpha[j],ElementV(mhm->alpha_scales,i));
    }

  dw_FreeArray(alpha);
}

/* PRECISION UpdatePosteriorLikelihood(T_MHM *mhm) */
/* { */
/*   PRECISION log_likelihood, log_posterior, difference; */

/*   log_likelihood=LogLikelihood_StatesIntegratedOut(mhm->model); */
/*   log_posterior=log_likelihood + LogPrior(mhm->model); */
/*   if (mhm->N > 1) */
/*     { */
/*       mhm->sum+=(difference=log_posterior - mhm->old_log_posterior); */
/*       mhm->sum_square+=difference*difference; */
/*     } */
/*   mhm->old_log_posterior=log_posterior; */
/*   if (log_likelihood > mhm->max_log_likelihood) mhm->max_log_likelihood=log_likelihood; */
/*   if (log_posterior > mhm->max_log_posterior) mhm->max_log_posterior=log_posterior; */

/*   return log_posterior; */
/* } */

/* void UpdateModifiedHarmonicMean(T_MHM *mhm, int n_singular) */
/* { */
/*   int j, k; */
/*   PRECISION quadratic_form, log_posterior; */

/*   // Increment total number of observations */
/*   mhm->N++; */

/*    // Log posterior - log likelihood corrected for normalization (this is now done in VARbase.c) */
/*   log_posterior=UpdatePosteriorLikelihood(mhm); // + ((T_VAR_Parameters*)mhm->model->theta)->nvars*log(2); */

/*   // Compute quadratic form */
/*   ConvertThetaToFreeParameters(mhm->model,pElementV(mhm->free_parameters_VAR)); */
/*   SubtractVV(mhm->free_parameters_VAR,mhm->free_parameters_VAR,mhm->center); */
/*   quadratic_form=InnerProductSymmetric(mhm->free_parameters_VAR,mhm->inverse_variance); */

/*   // Print log posterior and quadratic form */
/*   fprintf(mhm->f_out,"%le %le",log_posterior,quadratic_form); */

/*   // Print Dirichlet PDF's */
/*   for (j=0; j < dw_DimA(mhm->Alpha); j++) */
/*     fprintf(mhm->f_out," %le",LogIndependentDirichlet_pdf(mhm->model->sv->ba,mhm->Alpha[j])); */

/*   // Print number of singular varinances */
/*   fprintf(mhm->f_out," %d\n",Get_VAR_Improper_Distribution_Counter()-n_singular); */

/*   // Print linefeed */
/*   //fprintf(mhm->f_out,"\n"); */

/*   // Tally states */
/*   for (j=mhm->model->sv->nstates-1; j >= 0; j--) mhm->states[j]=0; */
/*   for (j=mhm->model->sv->nobs; j > 1; j--) mhm->states[mhm->model->sv->S[j]]++; */
/*   for (j=mhm->model->sv->nstates-1; j >= 0; j--) fprintf(mhm->f_out_regime_counts,"%d ",mhm->states[j]); */
/*   fprintf(mhm->f_out_regime_counts,"\n"); */
/* } */

void UpdateModifiedHarmonicMean(T_MHM *mhm, int n_singular)
{
  int j, k;
  PRECISION quadratic_form, log_likelihood, log_likelihood_states_integrated_out,
    log_prior_theta, log_prior_Q, log_posterior, difference;

  // Increment total number of observations
  mhm->N++;

  // Compute likelihoods and priors
  log_likelihood=LogLikelihood(mhm->model);
  log_likelihood_states_integrated_out=LogLikelihood_StatesIntegratedOut(mhm->model);
  log_prior_theta=LogPrior_Theta(mhm->model);
  log_prior_Q=LogPrior_Q(mhm->model);
  log_posterior=log_likelihood_states_integrated_out + log_prior_theta + log_prior_Q;

  // Average change
  if (mhm->N > 1)
    {
      mhm->sum+=(difference=log_posterior - mhm->old_log_posterior);
      mhm->sum_square+=difference*difference;
    }
  mhm->old_log_posterior=log_posterior;

  // Maximum likelihoods and priors
  if (log_likelihood_states_integrated_out > mhm->max_log_likelihood) mhm->max_log_likelihood=log_likelihood_states_integrated_out;
  if (log_posterior > mhm->max_log_posterior) mhm->max_log_posterior=log_posterior;

  // Compute quadratic form
  ConvertThetaToFreeParameters(mhm->model,pElementV(mhm->free_parameters_VAR));
  SubtractVV(mhm->free_parameters_VAR,mhm->free_parameters_VAR,mhm->center);
   quadratic_form=InnerProductSymmetric(mhm->free_parameters_VAR,mhm->inverse_variance);

  /*** Standard output ***/
  // Print log posterior and quadratic form
  fprintf(mhm->f_out,"%le %le",log_posterior,quadratic_form);
  // Print Dirichlet PDF's
  for (j=0; j < dw_DimA(mhm->Alpha); j++)
    fprintf(mhm->f_out," %le",LogIndependentDirichlet_pdf(mhm->model->sv->ba,mhm->Alpha[j]));
  // Print number of singular varinances
  fprintf(mhm->f_out," %d\n",Get_VAR_Improper_Distribution_Counter()-n_singular);

  /*** States not integrated out output ***/
  //if (mhm->f_states_not_integrated_out)
  //  fprintf(mhm->f_states_not_integrated_out,"%le %le %le %le %le\n",quadratic_form,log_likelihood,log_prior_theta,log_prior_Q,log_posterior);

  // Tally states
  for (j=mhm->model->sv->nstates-1; j >= 0; j--) mhm->states[j]=0;
  for (j=mhm->model->sv->nobs; j > 1; j--) mhm->states[mhm->model->sv->S[j]]++;
  for (j=mhm->model->sv->nstates-1; j >= 0; j--) fprintf(mhm->f_out_regime_counts,"%d ",mhm->states[j]);
  fprintf(mhm->f_out_regime_counts,"\n");
}

void ComputeModifiedHarmonicMean(T_MHM *mhm, int period)
{
  FILE *f_tmp;
  char *header;

  int count, check=period, begin_time, i, n_singular;
  printf("Beginning modified harmonic mean calculation -- %d iterations.\n",mhm->n_mhm);
  begin_time=time((time_t*)NULL);
  for (count=1; count <= mhm->n_mhm; count++)
    {
      n_singular=Get_VAR_Improper_Distribution_Counter();
      for (i=mhm->n_thin; i > 0; i--)
    {
      DrawAll(mhm->model);
    }

      UpdateModifiedHarmonicMean(mhm,n_singular);
      if (count == check)
    {
      check+=period;
      printf("Total Elapsed Time: %d seconds\n",(int)time((time_t*)NULL) - begin_time);
      printf("%d iterations completed out of %d\n",count,mhm->n_mhm);
      PrintJumps(stdout,(T_VAR_Parameters*)(mhm->model->theta));
    }
    }
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/******************************** Input/Output *********************************/
/*******************************************************************************/
static int ReadError_MHMio(char *id)
{
  char *errmsg, *fmt="Error involving line identifier \"%s\"";
  sprintf(errmsg=(char*)malloc(strlen(fmt) + strlen(id) - 1),fmt,id);
  dw_UserError(errmsg);
  free(errmsg);
  return 0;
}

void PrintJumps(FILE *f, T_VAR_Parameters *p)
{
  if(f==stdout)
    printf("Jumping counts - Total: %d\n",p->Total_A0_Metropolis_Draws);
  else
    fprintf(f,"Jumping counts - Total: %d\n",p->Total_A0_Metropolis_Draws);

  dw_PrintArray(f,p->A0_Metropolis_Jumps,"%7d ");
}

void WriteMHM_Input(FILE *f_out, T_MHM *mhm)
{
  fprintf(f_out,"//== scale values for Dirichlet distribution ==//\n%d\n",DimV(mhm->alpha_scales));
  dw_PrintVector(f_out,mhm->alpha_scales,"%22.14le ");
  fprintf(f_out,"\n");

  fprintf(f_out,"//== number draws for first burn-in ==//\n%d\n\n",mhm->n_burn1);

  fprintf(f_out,"//== number draws for second burn-in ==//\n%d\n\n",mhm->n_burn2);

  fprintf(f_out,"//== number draws to estimate mean and variance ==//\n%d\n\n",mhm->n_mean_variance);

  fprintf(f_out,"//== number draws for modified harmonic mean process ==//\n%d\n\n",mhm->n_mhm);

  fprintf(f_out,"//== thinning factor for modified harmonic mean process ==//\n%d\n\n",mhm->n_thin);
}

T_MHM* ReadMHM_Input(FILE *f, char *filename, T_MHM *mhm)
{
  T_MHM *rtrn=mhm ? mhm : CreateMHM();
  FILE *f_in=f ? f : dw_OpenTextFile(filename);
  char *id;
  int m;

  id="//== scale values for Dirichlet distribution ==//";
  if (dw_SetFilePosition(f_in,id) && (fscanf(f_in," %d ",&m) == 1) && dw_ReadVector(f_in,rtrn->alpha_scales=CreateVector(m)))
    {
      id="//== number draws for first burn-in ==//";
      if (dw_SetFilePosition(f_in,id) && (fscanf(f_in," %d ",&(rtrn->n_burn1)) == 1))
    {
      id="//== number draws for second burn-in ==//";
      if (dw_SetFilePosition(f_in,id) && (fscanf(f_in," %d ",&(rtrn->n_burn2)) == 1))
        {
          id="//== number draws to estimate mean and variance ==//";
          if (dw_SetFilePosition(f_in,id) && (fscanf(f_in," %d ",&(rtrn->n_mean_variance)) == 1))
        {
          id="//== number draws for modified harmonic mean process ==//";
          if (dw_SetFilePosition(f_in,id) && (fscanf(f_in," %d ",&(rtrn->n_mhm)) == 1))
            {
              id="//== thinning factor for modified harmonic mean process ==//";
              if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %d ",&(rtrn->n_thin)) != 1))
            rtrn->n_thin=1;
              if (!f) fclose(f_in);
              return rtrn;
            }
        }
        }
    }
    }
  if (!mhm) FreeMHM(rtrn);
  if (!f) fclose(f_in);
  ReadError_MHMio(id);
  return (T_MHM*)NULL;
}

void WriteMeanVariance(FILE *f_out, T_MHM *mhm)
{
  fprintf(f_out,"//== Base Dirichlet parameters ==//\n");
  dw_PrintArray(f_out,mhm->BaseAlpha,"%22.14le ");

  fprintf(f_out,"//== Variance ==//\n");
  dw_PrintMatrix(f_out,mhm->variance,"%22.14le ");
  fprintf(f_out,"\n");

  fprintf(f_out,"//== Center ==//\n");
  dw_PrintVector(f_out,mhm->center,"%22.14le ");
  fprintf(f_out,"\n");

  fprintf(f_out,"//== Mean ==//\n");
  dw_PrintVector(f_out,mhm->mean,"%22.14le ");
  fprintf(f_out,"\n");

  fprintf(f_out,"//== Posterior mode VAR parameters ==//\n");
  dw_PrintVector(f_out,mhm->posterior_mode_VAR,"%22.14le ");
  fprintf(f_out,"\n");
}

int ReadMeanVariance(FILE *f_in, T_MHM *mhm)
{
  char *id;
  int i, j;

  id="//== Base Dirichlet parameters ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,mhm->BaseAlpha))
    return ReadError_MHMio(id);

  // Create alpha
  mhm->Alpha=dw_CreateArray_array(DimV(mhm->alpha_scales));
  for (i=dw_DimA(mhm->Alpha)-1; i >= 0; i--)
    {
      mhm->Alpha[i]=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
      for (j=dw_DimA(mhm->Alpha[i])-1; j >= 0; j--)
    mhm->Alpha[i][j]=ProductVS((TVector)NULL,mhm->BaseAlpha[j],ElementV(mhm->alpha_scales,i));
    }

  id="//== Variance ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,mhm->variance))
    return ReadError_MHMio(id);
  PsudoInverse(mhm->inverse_variance,mhm->variance);

  id="//== Center ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadVector(f_in,mhm->center))
    return ReadError_MHMio(id);

  id="//== Mean ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadVector(f_in,mhm->mean))
    return ReadError_MHMio(id);

  id="//== Posterior mode VAR parameters ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadVector(f_in,mhm->posterior_mode_VAR))
    return ReadError_MHMio(id);

  return 1;
}


/* T_MHM* Read_MHM(FILE *f, char *filename, TStateModel *model) */
/* { */
/*   char *id;  */

/*   id="//== p-values for gaussian truncation ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %d ",&n) != 1) || !dw_ReadVector(f_in,p_cuts=CreateVector(n))) dw_Error(PARSE_ERR); */

/*   id="//== zeta truncation values ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || !dw_ReadVector(f_in,zeta_cuts=CreateVector(2))) dw_Error(PARSE_ERR); */

/*   id="//== scale values for Dirichlet distribution ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %d ",&m) != 1) || !dw_ReadVector(f_in,dirichlet_scales=CreateVector(m))) dw_Error(PARSE_ERR); */

/*   id="//== initial Metropolis scale values ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || !dw_ReadVector(f_in,metropolis_scales=CreateVector(((T_VAR_Parameters*)(model->theta))->nvars))) dw_Error(PARSE_ERR); */
/*   SetupMetropolisInformation(metropolis_scales,(T_VAR_Parameters*)(model->theta)); */
/* } */

/* void PrintMHM(FILE *f, char *filename, TStateModel *model, T_MHM *mhm) */
/* { */
/*   int i, j; */
/*   FILE *f_out; */

/*   f_out=f ? f : dw_CreateTextFile(filename); */

/*   fprintf(f_out,"Log of marginal data density \n"); */
/*   for (i=0; i < RowM(mhm->log_sum); i++) */
/*     { */
/*       for (j=0; j < ColM(mhm->log_sum); j++) */
/*     fprintf(f_out,"%le ",log(mhm->N) - ElementM(mhm->log_sum,i,j)); */
/*       fprintf(f_out,"\n"); */
/*     } */
/*   fprintf(f_out,"\n"); */

/*   fprintf(f_out,"Total number of draws used to compute marginal data density\n%d\n\n",mhm->N); */

/*   fprintf(f_out,"For each p-value, percentage of non-zero terms in sum\n"); */
/*   for (i=0; i < mhm->n; i++) */
/*     if (mhm->N > 0) */
/*       fprintf(f_out,"%5.2lf ",(double)mhm->K[i]/(double)mhm->N); */
/*     else */
/*       fprintf(f_out,"0 "); */
/*   fprintf(f_out,"\n\n"); */

/*   fprintf(f_out,"Log values of sums of h( )/(Loglikelihood*Prior)\n"); */
/*   dw_PrintMatrix(f_out,mhm->log_sum,"%le "); */
/*   fprintf(f_out,"\n"); */

/*   fprintf(f_out,"Log value of the maximum of h( )/(Loglikelihood*Prior)\n"); */
/*   dw_PrintMatrix(f_out,mhm->log_max,"%le "); */
/*   fprintf(f_out,"\n"); */

/*   fprintf(f_out,"p-values for gaussian distribution\n"); */
/*   dw_PrintVector(f_out,mhm->p_values,"%5.3lf "); */
/*   fprintf(f_out,"\n"); */

/*   fprintf(f_out,"Cut points for the zeta\n"); */
/*   dw_PrintVector(f_out,mhm->zeta_cuts,"%lf "); */
/*   fprintf(f_out,"\n"); */

/*   fprintf(f_out,"Scaling factor for zeta truncation\n"); */
/*   dw_PrintVector(f_out,mhm->zeta_p_values,"%lf "); */
/*   fprintf(f_out,"\n"); */

/*   fprintf(f_out,"Log likelihood, posterior, and prior evauated at posterior peak\n"); */
/*   fprintf(f_out,"%lf  %lf  %lf\n\n",mhm->log_likelihood_at_mode,mhm->log_posterior_at_mode,mhm->log_prior_at_mode); */

/*   fprintf(f_out,"Maximum draw of log likelihood and posterior\n"); */
/*   fprintf(f_out,"%lf  %lf\n\n",mhm->max_log_likelihood,mhm->max_log_posterior); */

/*   fprintf(f_out,"Mean and standard deviation of the log ratio of the posterior kernel of successive draws\n%le  %lf\n\n", */
/*       mhm->sum/(double)mhm->N,sqrt((mhm->sum_square - mhm->sum*mhm->sum/(double)mhm->N)/(double)mhm->N)); */

/*   PrintJumps(f_out,(T_VAR_Parameters*)(model->theta)); */

/*   fprintf(f_out,"Total number of draws: %d\n",mhm->N); */
/*   fprintf(f_out,"Number of draws rejected because of zeta truncation: %d\n",mhm->zeta_truncations); */
/*   fprintf(f_out,"Number of draws rejected because of gaussian truncation.\n"); */
/*   dw_PrintArray(f_out,mhm->gaussian_truncations,"%d "); */

/*   if (!f) fclose(f_out); */
/* } */
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/



