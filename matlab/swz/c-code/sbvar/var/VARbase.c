
#include "VARbase.h"
#include "VARio.h"
#include "switch.h"
#include "switchio.h"
#include "dw_error.h"
#include "swzmatrix.h"
#include "bmatrix.h"
#include "dw_array.h"
#include "dw_matrix_array.h"
#include "dw_rand.h"
#include "dw_matrix_rand.h"
#include "dw_ascii.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "modify_for_mex.h"

/*  //=== Private counter for improper normal distribution ===   ansi-c*/
static int _VAR_IMPROPER_DISTRIBUTION_COUNTER = 0;
int Reset_VAR_Improper_Distribution_Counter(void)
{
  int rtrn=_VAR_IMPROPER_DISTRIBUTION_COUNTER;
  _VAR_IMPROPER_DISTRIBUTION_COUNTER=0;
  return rtrn;
}

int Get_VAR_Improper_Distribution_Counter(void)
{
  return _VAR_IMPROPER_DISTRIBUTION_COUNTER;
}

#define BPLUS_ERR  1
#define PSI_ERR    2
#define LAMBDA_ERR 3

static int _SINGULAR_ERROR = 0;

static int _VERBOSE_COUNT = 0;
FILE *V_FILE = (FILE*)NULL;

void Increment_Verbose(void)
{
  _VERBOSE_COUNT++;
}

void SetVerboseFile(FILE *f)
{
  V_FILE=f;
}

int Get_VAR_Sigular_Error(void)
{
  return _SINGULAR_ERROR;
}

TVector DrawNormal_InverseVariance_SVD(TVector x, TVector b, TMatrix S);
TVector SingularInverseVariance_RecoveryAttempt(TVector x, TVector b, TMatrix S, TMatrix InversePrior, TStateModel *model, int code);
/*  //=========================================================   ansi-c*/


/*  //=== Private Utility Functions ===   ansi-c*/
static int NumberStates(TMarkovStateVariable** sv, int n);
/*  //static int* CreateStateIndex(TMarkovStateVariable* top, TMarkovStateVariable** sv, int n);   ansi-c*/

extern FILE *fptr_debug;

/*******************************************************************************/
/************************** Constructors/Destructors ***************************/
/*******************************************************************************/
void FreeTheta_VAR(T_VAR_Parameters *p)
{
  int j;
  if (p)
    {
/*        // Free parameters   ansi-c*/
      dw_FreeArray(p->A0);
      dw_FreeArray(p->Aplus);
      dw_FreeArray(p->Zeta);

/*        // Free state variable translation   ansi-c*/
      dw_FreeArray(p->n_var_states);
      dw_FreeArray(p->var_states);
      dw_FreeArray(p->n_coef_states);
      dw_FreeArray(p->coef_states);
      dw_FreeArray(p->A0_states);
      dw_FreeArray(p->A0_column_states);

/*        // Free free parameters   ansi-c*/
      dw_FreeArray(p->dim_b0);
      dw_FreeArray(p->b0);
      dw_FreeArray(p->dim_bplus);
      dw_FreeArray(p->bplus);

/*        // Free Sims-Zha specification parameters and workspace   ansi-c*/
      dw_FreeArray(p->lambda);

/*        //--- Non-standard memory management ---   ansi-c*/
      if (p->constant)
    for (j=dw_DimA(p->constant)-1; j >= 0; j--)
      if (p->constant[j])
        pElementV(p->constant[j])=(PRECISION*)NULL;
      dw_FreeArray(p->constant);
/*        //--------------------------------------   ansi-c*/

      dw_FreeArray(p->psi);
      dw_FreeArray(p->inverse_psi_prior);

/*        // Free Priors   ansi-c*/
      FreeVector(p->zeta_a_prior);
      FreeVector(p->zeta_b_prior);
      dw_FreeArray(p->A0_prior);
      dw_FreeArray(p->Aplus_prior);

/*        // Free identifiying restrictions   ansi-c*/
      dw_FreeArray(p->U);
      dw_FreeArray(p->V);
      dw_FreeArray(p->W);
      dw_FreeArray(p->IsIdentity_V);

/*        // Free normalization   ansi-c*/
      dw_FreeArray(p->flipped);
      dw_FreeArray(p->Target);

/*        // Free workspace   ansi-c*/
      FreeVector(p->inverse_zeta_b_prior);
      dw_FreeArray(p->inverse_b0_prior);
      dw_FreeArray(p->inverse_bplus_prior);
      FreeVector(p->log_abs_det_A0);
      dw_FreeArray(p->A0_dot_products);
      dw_FreeArray(p->Aplus_dot_products);

/*        // Free state dependent fields   ansi-c*/
      dw_FreeArray(p->YY);
      dw_FreeArray(p->XY);
      dw_FreeArray(p->XX);
      dw_FreeArray(p->yy);
      dw_FreeArray(p->xy);
      dw_FreeArray(p->xx);
      dw_FreeArray(p->S);
      dw_FreeArray(p->T);

/*        // A0 Metropolis Info   ansi-c*/
      dw_FreeArray(p->A0_Metropolis_Scale);
      dw_FreeArray(p->A0_Metropolis_Jumps);

/*        // Free Data   ansi-c*/
      dw_FreeArray(p->Y);
      dw_FreeArray(p->X);

/*        // Free pointer   ansi-c*/
      swzFree(p);
    }
}

ThetaRoutines* CreateRoutines_VAR(void)
{
  ThetaRoutines *rtns=CreateThetaRoutines_empty();

  rtns->pLogConditionalLikelihood=LogConditionalProbability_VAR;
  rtns->pExpectationSingleStep=ExpectationSingleStep_VAR;
  rtns->pDestructor=(void (*)(void*))FreeTheta_VAR;
  rtns->pLogPrior=LogPrior_VAR;
  rtns->pNumberFreeParametersTheta=NumberFreeParametersVAR;
  rtns->pConvertFreeParametersToTheta=FreeParametersToVAR;
  rtns->pConvertThetaToFreeParameters=VARToFreeParameters;
  rtns->pDrawParameters=DrawParameters_VAR;
  rtns->pStatesChanged=StatesChanged_VAR;
  rtns->pThetaChanged=ThetaChanged_VAR;
  rtns->pInitializeForwardRecursion=InitializeForwardRecursion_VAR;

  return rtns;
}

T_VAR_Parameters* CreateTheta_VAR(int flag, int nvars, int nlags, int nexg, int nstates, int nobs,     /*   Specification and Sizes   ansi-c*/
                  int **coef_states, int **var_states,                                 /*   Translation Tables   ansi-c*/
                  TMatrix *U, TMatrix *V, TMatrix *W,                                  /*   Restrictions   ansi-c*/
                  TMatrix Y, TMatrix X)                                                /*   Data   ansi-c*/
{
  T_VAR_Parameters *p;
  int i, j, k, t, npre;
  TMatrix S;

  if ((nvars <= 0) || (nlags < 0) || (nexg < 0))
    {
      swz_fprintf_err("CreateTheta_VAR():  Invalid arguments passed.\n");
      swzExit(0);
    }

/*    //=== Allocate memory for T_VAR_Parameters ===   ansi-c*/
  if (!(p=(T_VAR_Parameters*)swzMalloc(sizeof(T_VAR_Parameters))))
    {
      swz_fprintf_err("Out of memory\n");
      swzExit(0);
    }

/*    //=== Model Specification ===   ansi-c*/
  if (flag & SPEC_SIMS_ZHA) flag|=SPEC_RANDOM_WALK;
  p->Specification=flag;

/*    //====== Flags ======   ansi-c*/
  p->valid_state_dependent_fields=0;
  p->valid_state_dependent_fields_previous=0;
  p->valid_log_abs_det_A0=0;
  p->valid_dot_products=0;
  p->valid_parameters=0;

/*    //=== Sizes ===//   ansi-c*/
  p->nobs=nobs;
  p->nstates=nstates;
  p->nvars=nvars;
  p->nlags=nlags;
  p->npre=npre=nvars*nlags+nexg;

/*    //====== Create n_coef_states, n_var_states ======   ansi-c*/
  p->n_coef_states=(int*)dw_CreateArray_int(nvars);
  p->n_var_states=(int*)dw_CreateArray_int(nvars);
  for (j=nvars-1; j >= 0; j--)
    {
      for (p->n_coef_states[j]=0, i=p->nstates-1; i >= 0; i--)
    if (coef_states[j][i] > p->n_coef_states[j])
      p->n_coef_states[j]=coef_states[j][i];
      p->n_coef_states[j]++;

      for (p->n_var_states[j]=0, i=nstates-1; i >= 0; i--)
    if (var_states[j][i] > p->n_var_states[j])
      p->n_var_states[j]=var_states[j][i];
      p->n_var_states[j]++;
    }

/*    //====== Create coef_states, var_states ======   ansi-c*/
  p->coef_states=(int**)dw_CopyArray((void*)NULL,coef_states);
  p->var_states=(int**)dw_CopyArray((void*)NULL,var_states);

/*    //====== Create A0, Aplus, Zeta ======   ansi-c*/
  p->A0=(TVector**)dw_CreateArray_array(nvars);
  p->Aplus=(TVector**)dw_CreateArray_array(nvars);
  p->Zeta=(PRECISION**)dw_CreateArray_array(nvars);
  for (j=nvars-1; j >= 0; j--)
    {
      p->A0[j]=dw_CreateArray_vector(p->n_coef_states[j]);
      for (k=p->n_coef_states[j]-1; k >= 0; k--)
        p->A0[j][k]=CreateVector(nvars);

      p->Aplus[j]=(TVector*)dw_CreateArray_vector(p->n_coef_states[j]);
      for (k=p->n_coef_states[j]-1; k >= 0; k--)
        p->Aplus[j][k]=CreateVector(npre);

      p->Zeta[j]=(PRECISION*)dw_CreateArray_scalar(p->n_var_states[j]);
    }

/*    //====== A0 Metropolis Info ======   ansi-c*/
  p->A0_Metropolis_Scale=dw_CreateArray_array(nvars);
  for (j=nvars-1; j >= 0; j--)
    p->A0_Metropolis_Scale[j]=dw_CreateArray_scalar(p->n_coef_states[j]);
  dw_InitializeArray_scalar(p->A0_Metropolis_Scale,1.0);

  p->A0_Metropolis_Jumps=dw_CreateArray_array(nvars);
  for (j=nvars-1; j >= 0; j--)
    p->A0_Metropolis_Jumps[j]=dw_CreateArray_int(p->n_coef_states[j]);
  dw_InitializeArray_int(p->A0_Metropolis_Jumps,0);

/*    //====== Create A0_states, A0_column_states, and log_det_abs_A0 ======   ansi-c*/
  p->A0_states=dw_CreateArray_int(nstates);
  for (p->n_A0_states=0, i=0; i < nstates; i++)
    {
      for (k=i-1; k >= 0; k--)
    {
      for (j=nvars-1; j >= 0; j--)
        if (coef_states[j][i] != coef_states[j][k]) break;
      if (j < 0) break;
    }
      p->A0_states[i]=(k < 0) ? p->n_A0_states++ : p->A0_states[k];
    }

  p->A0_column_states=dw_CreateRectangularArray_int(nvars,p->n_A0_states);
  for (i=0; i < p->n_A0_states; i++)
    for (k=0; k < nstates; k++)
      if (p->A0_states[k] == i)
    {
      for (j=nvars-1; j >= 0; j--)
        p->A0_column_states[j][i]=coef_states[j][k];
      break;
    }

  InitializeVector(p->log_abs_det_A0=CreateVector(p->n_A0_states),-1.0);

/*    //=== Set Restrictions ===   ansi-c*/
  p->U=dw_CopyArray(NULL,U);
  p->b0=(TVector**)dw_CreateArray_array(nvars);
  p->dim_b0=dw_CreateArray_int(nvars);
  for (j=nvars-1; j >= 0; j--)
    {
      p->b0[j]=(TVector*)dw_CreateArray_vector(p->n_coef_states[j]);
      for (k=p->n_coef_states[j]-1; k >= 0; k--)
    p->b0[j][k]=CreateVector(p->dim_b0[j]=ColM(U[j]));
    }

/*    //=== Normalization ===   ansi-c*/
  p->normalization_type=VAR_NORMALIZATION_NONE;
  p->Target=(TVector**)NULL;
  p->flipped=(int**)NULL;
  p->WZ_inconsistancies=0;

/*    //=== Specification ===   ansi-c*/
  if (flag & SPEC_RANDOM_WALK)
    {
      p->W=dw_CreateArray_matrix(nvars);
      InitializeMatrix(S=CreateMatrix(npre,nvars),0.0);
      for (j=nvars-1; j >= 0; j--) ElementM(S,j,j)=-1.0;
      for (j=nvars-1; j >= 0; j--) p->W[j]=EquateMatrix((TMatrix)NULL,S);
      FreeMatrix(S);
    }
  else
    p->W=dw_CopyArray((void*)NULL,W);

  if (flag & SPEC_SIMS_ZHA)
    {
      dw_InitializeArray_int(p->IsIdentity_V=dw_CreateArray_int(nvars),1);
      p->V=dw_CreateArray_matrix(nvars);
      for (j=p->nvars-1; j >= 0; j--) p->V[j]=IdentityMatrix((TMatrix)NULL,npre);

/*        // Setup psi and lambda parameters   ansi-c*/
      p->lambda=(TVector**)dw_CreateArray_array(nvars);
      p->psi=dw_CreateArray_vector(nvars);
      for (j=nvars-1; j >= 0; j--)
    {
      p->lambda[j]=dw_CreateArray_vector(p->n_coef_states[j]);
      for (k=dw_DimA(p->lambda[j])-1; k >= 0; k--)
        p->lambda[j][k]=CreateVector(nvars);

      p->psi[j]=CreateVector(npre - 1 + p->n_coef_states[j]);
    }

/*        //--- non-standard memory management ---   ansi-c*/
      p->constant=(TVector*)dw_CreateArray_vector(nvars);
      for (j=nvars-1; j >= 0; j--)
    {
      p->constant[j]=CreateVector(p->n_coef_states[j]);
      swzFree(pElementV(p->constant[j]));
      pElementV(p->constant[j])=pElementV(p->psi[j]) + npre - 1;
    }
/*        //--------------------------------------   ansi-c*/
    }
  else
    {
/*        //====== Sims-Zha Specification ======   ansi-c*/
      p->lambda=(TVector**)NULL;
      p->psi=(TVector*)NULL;
      p->constant=(TVector*)NULL;

/*        //====== If the number of columns in V[j] == npre then we may assume that V[j] is the identity. ======   ansi-c*/
      dw_InitializeArray_int(p->IsIdentity_V=dw_CreateArray_int(nvars),0);
      p->V=dw_CreateArray_matrix(nvars);
      for (j=nvars-1; j >= 0; j--)
    if (V[j])
      if (ColM(V[j]) < npre)
        p->V[j]=EquateMatrix((TMatrix)NULL,V[j]);
      else
        {
          p->V[j]=IdentityMatrix((TMatrix)NULL,npre);
          p->IsIdentity_V[j]=1;
        }
    }
  p->bplus=(TVector**)dw_CreateArray_array(nvars);
  p->dim_bplus=dw_CreateArray_int(nvars);
  for (j=nvars-1; j >= 0; j--)
    if (V[j])
      {
    p->bplus[j]=(TVector*)dw_CreateArray_vector(p->n_coef_states[j]);
    for (k=p->n_coef_states[j]-1; k >= 0; k--)
      p->bplus[j][k]=CreateVector(p->dim_bplus[j]=ColM(V[j]));
      }

/*    //====== Data ======   ansi-c*/
  p->Y=dw_CreateArray_vector(nobs+1);
  p->X=dw_CreateArray_vector(nobs+1);
  for (t=nobs; t > 0; t--)
    {
      p->Y[t]=CreateVector(nvars);
      for (i=nvars-1; i >= 0; i--) ElementV(p->Y[t],i)=ElementM(Y,t-1,i);

      p->X[t]=CreateVector(p->npre);
      for (i=p->npre-1; i >= 0; i--) ElementV(p->X[t],i)=ElementM(X,t-1,i);
    }

/*    //====== Workspace  ======   ansi-c*/
  p->minus_half_nvars_times_log2pi=-0.5*(double)nvars*log(2.0*3.141592653589793);

/*    // Dot products   ansi-c*/
  p->A0_dot_products=(PRECISION***)dw_CreateArray_array(nobs+1);
  for (t=0; t <= nobs; t++)
    {
      p->A0_dot_products[t]=(PRECISION**)dw_CreateArray_array(nvars);
      for (j=0; j < nvars; j++)
    p->A0_dot_products[t][j]=dw_CreateArray_scalar(p->n_coef_states[j]);
    }

  p->Aplus_dot_products=(PRECISION***)dw_CreateArray_array(nobs+1);
  for (t=0; t <= nobs; t++)
    {
      p->Aplus_dot_products[t]=(PRECISION**)dw_CreateArray_array(nvars);
      for (j=0; j < nvars; j++)
    p->Aplus_dot_products[t][j]=dw_CreateArray_scalar(p->n_coef_states[j]);
    }

/*    // State dependent data   ansi-c*/
  p->T=dw_CreateArray_int(nstates);
  p->YY=dw_CreateArray_matrix(nstates);
  p->XY=dw_CreateArray_matrix(nstates);
  p->XX=dw_CreateArray_matrix(nstates);
  for (k=nstates-1; k >= 0; k--)
    {
      p->YY[k]=CreateMatrix(nvars,nvars);
      p->XY[k]=CreateMatrix(p->npre,nvars);
      p->XX[k]=CreateMatrix(p->npre,p->npre);
    }
  p->yy=dw_CreateArray_matrix(nobs+1);
  p->xy=dw_CreateArray_matrix(nobs+1);
  p->xx=dw_CreateArray_matrix(nobs+1);
  p->S=dw_CreateArray_int(nobs+1);
  for (t=nobs; t > 0; t--)
    {
      p->yy[t]=OuterProduct((TMatrix)NULL,p->Y[t],p->Y[t]);
      p->xy[t]=OuterProduct((TMatrix)NULL,p->X[t],p->Y[t]);
      p->xx[t]=OuterProduct((TMatrix)NULL,p->X[t],p->X[t]);
    }

/*    //====== Set Priors to null ======   ansi-c*/
  p->A0_prior=(TMatrix*)NULL;
  p->Aplus_prior=(TMatrix*)NULL;
  p->zeta_a_prior=(TVector)NULL;
  p->zeta_b_prior=(TVector)NULL;
  p->lambda_prior=0.0;

  p->inverse_b0_prior=(TMatrix*)NULL;
  p->inverse_bplus_prior=(TMatrix*)NULL;
  p->inverse_zeta_b_prior=(TVector)NULL;
  p->inverse_lambda_prior=0.0;
  p->inverse_psi_prior=(TMatrix*)NULL;

/*    //====== Return ======   ansi-c*/
  return p;
}

void SetPriors_VAR(T_VAR_Parameters *theta, TMatrix* A0_prior, TMatrix* Aplus_prior, TVector zeta_a_prior, TVector zeta_b_prior)
{
  int j;
  TMatrix S;

/*    //====== Priors ======   ansi-c*/
  theta->A0_prior=dw_CopyArray(NULL,A0_prior);
  theta->Aplus_prior=dw_CopyArray(NULL,Aplus_prior);
  theta->zeta_a_prior=EquateVector((TVector)NULL,zeta_a_prior);
  theta->zeta_b_prior=EquateVector((TVector)NULL,zeta_b_prior);

/*    //====== Prior workspace ======   ansi-c*/
  theta->inverse_zeta_b_prior=CreateVector(theta->nvars);
  for (j=theta->nvars-1; j >= 0; j--)
    ElementV(theta->inverse_zeta_b_prior,j)=1.0/ElementV(zeta_b_prior,j);

  theta->inverse_b0_prior=dw_CreateArray_matrix(theta->nvars);
  for (j=theta->nvars-1; j >= 0; j--)
    {
      S=Inverse_LU((TMatrix)NULL,A0_prior[j]);
      theta->inverse_b0_prior[j]=MatrixInnerProductSymmetric((TMatrix)NULL,theta->U[j],S);
      ProductMS(theta->inverse_b0_prior[j],theta->inverse_b0_prior[j],theta->n_A0_states/theta->n_coef_states[j]);
      FreeMatrix(S);
    }

  theta->inverse_bplus_prior=dw_CreateArray_matrix(theta->nvars);
  for (j=theta->nvars-1; j >= 0; j--)
    if (theta->V[j])
      {
    S=Inverse_LU((TMatrix)NULL,Aplus_prior[j]);
    theta->inverse_bplus_prior[j]=MatrixInnerProductSymmetric((TMatrix)NULL,theta->V[j],S);
    ProductMS(theta->inverse_bplus_prior[j],theta->inverse_bplus_prior[j],theta->n_A0_states/theta->n_coef_states[j]);
    FreeMatrix(S);
      }

/*    //====== Set prior constant ======   ansi-c*/
  SetLogPriorConstant_VAR(theta);
}

void SetPriors_VAR_SimsZha(T_VAR_Parameters *theta, TMatrix* A0_prior, TMatrix* Aplus_prior, TVector zeta_a_prior,
               TVector zeta_b_prior, PRECISION lambda_prior)
{
  int j, k, n, m;
  TMatrix V, S, I, inverse_Aplus_prior;

  if (theta->Specification & SPEC_SIMS_ZHA)
    {
      theta->lambda_prior=lambda_prior;
      theta->inverse_lambda_prior=1.0/lambda_prior;

      theta->inverse_psi_prior=dw_CreateArray_matrix(theta->nvars);
      for (j=theta->nvars-1; j >= 0; j--)
    {
      inverse_Aplus_prior=Inverse_LU((TMatrix)NULL,Aplus_prior[j]);
      I=IdentityMatrix((TMatrix)NULL,theta->n_A0_states);
      S=KroneckerProduct((TMatrix)NULL,I,inverse_Aplus_prior);
      InitializeMatrix(V=CreateMatrix(theta->npre * theta->n_A0_states,theta->npre-1 + theta->n_coef_states[j]),0.0);
      for (k=theta->n_A0_states-1; k >= 0; k--)
        {
          for (n=theta->npre-2; n >= 0; n--)
        for (m=theta->npre-1; m >= 0; m--)
          ElementM(V,k*theta->npre + m,n)=ElementM(theta->V[j],m,n);

          n=theta->npre - 1 + theta->A0_column_states[j][k];
          for (m=theta->npre-1; m >= 0; m--)
        ElementM(V,k * theta->npre + m,n)=ElementM(theta->V[j],m,theta->npre - 1);
        }
      theta->inverse_psi_prior[j]=MatrixInnerProductSymmetric((TMatrix)NULL,V,S);
      FreeMatrix(V);
      FreeMatrix(S);
      FreeMatrix(I);
      FreeMatrix(inverse_Aplus_prior);
    }

      SetPriors_VAR(theta,A0_prior,Aplus_prior,zeta_a_prior,zeta_b_prior);
    }
  else
    {
      printf("Error SetPriors_VAR_SimsZha(): specification flag not set to SPEC_SIM_ZHA\n");
      swzExit(0);
    }
}

/*
   Assumes
    model: a properly initialized TStateModel structure.

   Returns
    A properly initialized TStateModel structure with the same number of first
    level state variables, but only one overall state.
*/
TStateModel* CreateConstantModel(TStateModel *model)
{
  TMarkovStateVariable *sv, **sv_array;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta), *theta;
  ThetaRoutines *routines;
  TMatrix X, Y;
  int **translation_table, i, j, t;

  if (model->sv->n_state_variables == 1)
    sv=CreateMarkovStateVariable_ConstantState(model->sv->nobs);
  else
    {
      sv_array=dw_CreateArray_pointer(model->sv->n_state_variables,(void (*)(void*))FreeMarkovStateVariable);
      for (i=model->sv->n_state_variables-1; i >= 0; i--)
    sv_array[i]=CreateMarkovStateVariable_ConstantState(model->sv->nobs);
      sv=CreateMarkovStateVariable_Multiple(model->sv->nobs,model->sv->n_state_variables,sv_array);
    }

  routines=CreateRoutines_VAR();

  Y=CreateMatrix(p->nobs,p->nvars);
  X=CreateMatrix(p->nobs,p->npre);
  for (t=p->nobs; t > 0; t--)
    {
      for (i=p->nvars-1; i >= 0; i--) ElementM(Y,t-1,i)=ElementV(p->Y[t],i);
      for (i=p->npre-1; i >= 0; i--) ElementM(X,t-1,i)=ElementV(p->X[t],i);
    }

  dw_InitializeArray_int(translation_table=dw_CreateRectangularArray_int(p->nvars,1),0);
  theta=CreateTheta_VAR(p->Specification,p->nvars,p->nlags,p->npre - p->nlags*p->nvars,1,p->nobs,translation_table,translation_table,p->U,p->V,p->W,Y,X);
  if (p->Specification & SPEC_SIMS_ZHA)
    SetPriors_VAR_SimsZha(theta,p->A0_prior,p->Aplus_prior,p->zeta_a_prior,p->zeta_b_prior,p->lambda_prior);
  else
    SetPriors_VAR(theta,p->A0_prior,p->Aplus_prior,p->zeta_a_prior,p->zeta_b_prior);

  FreeMatrix(X);
  FreeMatrix(Y);
  dw_FreeArray(translation_table);

/*    // Initialize parameters   ansi-c*/
  for (j=0; j < p->nvars; j++)
    {
      theta->Zeta[j][0]=p->Zeta[j][0];
      EquateVector(theta->A0[j][0],p->A0[j][0]);
      EquateVector(theta->Aplus[j][0],p->Aplus[j][0]);
    }
  Update_b0_bplus_from_A0_Aplus(theta);
  if ((theta->Specification & SPEC_SIMS_ZHA) == SPEC_SIMS_ZHA)
    Update_lambda_psi_from_bplus(theta);
  theta->valid_parameters=1;

/*    //InitializeParameters_VAR(theta);   ansi-c*/

  return CreateStateModel_new(sv,routines,theta);
}

/*
   Assumes
     table:  Properly initialized translation table.
     restricted table:  Properly initialize translation table.
*/
int **ExpandTranslationTable(int **table, TMarkovStateVariable *sv, TMarkovStateVariable *rsv, int s)
{
  int i, j, k, nstates;
  int **rtable, *idx, *master;

/*    // Compute size of new table.   ansi-c*/
  for (nstates=1, k=sv->n_state_variables-1; k >= 0; k--)
    nstates*=(k == s) ? rsv->state_variable[k]->nstates + 1 : rsv->state_variable[k]->nstates;

/*    // Create new table   ansi-c*/
  dw_InitializeArray_int(rtable=dw_CreateRectangularArray_int(dw_DimA(table),nstates),0);

/*    // Fill table   ansi-c*/
  idx=(int*)swzMalloc(sv->nstates*sizeof(int));
  master=(int*)swzMalloc(nstates*sizeof(int));
  for (k=i=0; k < sv->nstates; k++)
    {
      for (j=sv->n_state_variables-1; j >= 0; j--)
    if (sv->Index[k][j] > rsv->state_variable[j]->nstates)
      break;
    else
      if ((j != s) && (sv->Index[k][j] == rsv->state_variable[j]->nstates))
        break;
      if (j == -1) master[i++]=k;
    }
  for (j=dw_DimA(table)-1; j >= 0; j--)
    {
      for (k=sv->nstates-1; k >= 0; k--) idx[k]=0;
      for (k=nstates-1; k >= 0; k--)
    idx[table[j][master[k]]]=1;
      for (k=i=0; k < sv->nstates; k++) idx[k]=idx[k] ? i++ : -1;
      for (k=nstates-1; k >= 0; k--)
    rtable[j][k]=idx[table[j][master[k]]];
    }
  swzFree(master);
  swzFree(idx);

/*    // verbose   ansi-c*/
/*    //dw_PrintArray(stdout,rtable,(char*)NULL); printf("\n"); dw_PrintArray(stdout,table,(char*)NULL); getchar();   ansi-c*/
/*    //dw_PrintArray(stdout,sv->Index,(char*)NULL); getchar();   ansi-c*/

  return rtable;
}

/*
   Assumes
     model:
       Properly intialized TStateModel structure.  The Markov state variable
       structure must be flat with the same number state variables as
       restricted_model.  Each state variable must have as least as many states
       as the corresponding state variable in restricted_model.

     restricted_model:
       Properly intialized TStateModel structure.  The Markov state variable
       structure must be flat with the same number state variables as model. Each
       state variable must have no more states as the corresponding state
       variable in model.

     s:
       The state variable to expand.  It must be the case that s is between 0 and
       model->n_state_variables - 1, inclusive, and that the number of states in
       the sth state variable of model is strictly larger than the number of
       states in the sth state variable of restricted_model.

   Returns
     Properly initialized TStateModel structure with in which the sth state
     variable has one more state.

   Notes
     A Markov state variable structure is flat if state_variable[i] is a single
     Markov state variable for 0 <= i < n_state_variables.
*/
TStateModel* ExpandModel_VAR(TStateModel *model, TStateModel *restricted_model, int s)
{
  int i, j, k, m, q, t;
  TStateModel *expanded_model;
  T_VAR_Parameters *p=model->theta, *restricted_p=restricted_model->theta, *expanded_p;
  TMarkovStateVariable *sv=model->sv, *restricted_sv=restricted_model->sv, *expanded_sv, **sv_array;
  TMatrix X, Y;
  int **coef_table, **var_table;

/*    // Check sizes   ansi-c*/
  if ((p->nvars != restricted_p->nvars) || (p->nlags != restricted_p->nlags)
      || (p->npre != restricted_p->npre) || (p->nobs != restricted_p->nobs)
      || (sv->n_state_variables != restricted_sv->n_state_variables)
      || (s < 0) || (s > sv->n_state_variables))
    return (TStateModel*)NULL;

  for (k=sv->n_state_variables-1; k >= 0; k--)
    if (sv->state_variable[k]->nstates < restricted_sv->state_variable[k]->nstates) return (TStateModel*)NULL;

  if (sv->state_variable[s]->nstates == restricted_sv->state_variable[s]->nstates) return (TStateModel*)NULL;

/*    // Check VAR Restrictions   ansi-c*/

/*    // Check VAR Priors   ansi-c*/

/*    // Setup new Markov state variable   ansi-c*/
  sv_array=dw_CreateArray_pointer(sv->n_state_variables,(void (*)(void*))FreeMarkovStateVariable);
  for (k=model->sv->n_state_variables-1; k >= 0; k--)
    if (k != s)
      sv_array[k]=DuplicateMarkovStateVariable(restricted_sv->state_variable[k]);
    else
      sv_array[s]=RestrictMarkovStateVariable(sv->state_variable[k],restricted_sv->state_variable[k]->nstates+1);

/*    // Create multiple Markov state variable   ansi-c*/
  if (sv->n_state_variables == 1)
    {
      expanded_sv=sv_array[0];
      sv_array[0]=(TMarkovStateVariable*)NULL;
      dw_FreeArray(sv_array);
    }
  else
    expanded_sv=CreateMarkovStateVariable_Multiple(sv->nobs,sv->n_state_variables,sv_array);

/*    // Data   ansi-c*/
  Y=CreateMatrix(p->nobs,p->nvars);
  X=CreateMatrix(p->nobs,p->npre);
  for (t=p->nobs; t > 0; t--)
    {
      for (i=p->nvars-1; i >= 0; i--) ElementM(Y,t-1,i)=ElementV(p->Y[t],i);
      for (i=p->npre-1; i >= 0; i--) ElementM(X,t-1,i)=ElementV(p->X[t],i);
    }

/*    // Setup new translation tables   ansi-c*/
  coef_table=ExpandTranslationTable(p->coef_states,model->sv,restricted_model->sv,s);
  var_table=ExpandTranslationTable(p->var_states,model->sv,restricted_model->sv,s);

/*    // Setup new VAR   ansi-c*/
  expanded_p=CreateTheta_VAR(p->Specification,p->nvars,p->nlags,p->npre - p->nlags*p->nvars,expanded_sv->nstates,p->nobs,coef_table,var_table,p->U,p->V,p->W,Y,X);
  if (p->Specification & SPEC_SIMS_ZHA)
    SetPriors_VAR_SimsZha(expanded_p,p->A0_prior,p->Aplus_prior,p->zeta_a_prior,p->zeta_b_prior,p->lambda_prior);
  else
    SetPriors_VAR(expanded_p,p->A0_prior,p->Aplus_prior,p->zeta_a_prior,p->zeta_b_prior);

/*    // Create expanded model   ansi-c*/
  expanded_model=CreateStateModel_new(expanded_sv,CreateRoutines_VAR(),expanded_p);

/*    // Clean up   ansi-c*/
  FreeMatrix(X);
  FreeMatrix(Y);
  dw_FreeArray(coef_table);
  dw_FreeArray(var_table);

/*    // Set VAR parameters   ansi-c*/
  for (j=0; j < p->nvars; j++)
    {
      for (k=expanded_p->n_var_states[j]-1; k >= 0; k--)
    {
      for (i=expanded_p->nstates-1; i >= 0; i--)
        if (expanded_p->var_states[j][i] == k) break;

      for (q=m=0; q <= s; q++) m=m*expanded_sv->state_variable[q]->nstates + expanded_sv->Index[i][q];
      if (expanded_sv->Index[i][s] == expanded_sv->state_variable[s]->nstates - 1) m--;
      for ( ; q < expanded_sv->n_state_variables; q++) m=m*expanded_sv->state_variable[q]->nstates + expanded_sv->Index[i][q];

      expanded_p->Zeta[j][k]=restricted_p->Zeta[j][restricted_p->var_states[j][m]];
    }

      for (k=expanded_p->n_coef_states[j]-1; k >= 0; k--)
    {
      for (i=expanded_p->nstates-1; i >= 0; i--)
        if (expanded_p->coef_states[j][i] == k) break;

      for (q=m=0; q <= s; q++) m=m*expanded_sv->state_variable[q]->nstates + expanded_sv->Index[i][q];
      if (expanded_sv->Index[i][s] == expanded_sv->state_variable[s]->nstates - 1) m--;
      for ( ; q < expanded_sv->n_state_variables; q++) m=m*expanded_sv->state_variable[q]->nstates + expanded_sv->Index[i][q];

      EquateVector(expanded_p->A0[j][k],restricted_p->A0[j][restricted_p->coef_states[j][m]]);
      EquateVector(expanded_p->Aplus[j][k],restricted_p->Aplus[j][restricted_p->coef_states[j][m]]);
    }
    }
  Update_b0_bplus_from_A0_Aplus(expanded_p);
  if ((expanded_p->Specification & SPEC_SIMS_ZHA) == SPEC_SIMS_ZHA)
    Update_lambda_psi_from_bplus(expanded_p);
  expanded_p->valid_parameters=1;

/*    // Set transition matrices   ansi-c*/
  if (!restricted_model->sv->valid_transition_matrix)
    {
      FreeStateModel(expanded_model);
      return (TStateModel*)NULL;
    }
  for (k=expanded_sv->n_state_variables-1; k >= 0; k--)
    {
      if (k != s)
    {
/*        //EquateMatrix(expanded_sv->state_variable[k]->Q,restricted_sv->state_variable[k]->Q);   ansi-c*/
/*        //Update_B_from_Q_SV(expanded_sv->state_variable[k]);   ansi-c*/

      EquateMatrix(expanded_sv->state_variable[k]->Q,restricted_sv->state_variable[k]->Q);
      EquateVector(expanded_sv->state_variable[k]->B,restricted_sv->state_variable[k]->B);
    }
      else
    {
/*       j=restricted_sv->state_variable[k]->nstates; */
/*       for (i=0; i <= restricted_sv->state_variable[k]->nstates; i++) */
/*         for (j=0; j <= restricted_sv->state_variable[k]->nstates; j++) */
/*           ElementM(expanded_sv->state_variable[k]->Q,i,j)=0.5; */
/*  //=====================================================================================================================   ansi-c*/
/*        j=restricted_sv->state_variable[k]->nstates;  */
/*       InsertSubMatrix(expanded_sv->state_variable[k]->Q,restricted_sv->state_variable[k]->Q,0,0,0,0,j,j); */
/*       for (i=j; i >= 0; i--) */
/*         { */
/*           ElementM(expanded_sv->state_variable[k]->Q,i,j)=0.0; */
/*           ElementM(expanded_sv->state_variable[k]->Q,j,i)=0.0; */
/*         } */
/*       ElementM(expanded_sv->state_variable[k]->Q,j,j)= 0.5*ElementM(expanded_sv->state_variable[k]->Q,j-1,j-1); */
/*       ElementM(expanded_sv->state_variable[k]->Q,j-1,j)=0.5*ElementM(expanded_sv->state_variable[k]->Q,j-1,j-1); */
/*       for (i=0; i < j-1; i++) */
/*         ElementM(expanded_sv->state_variable[k]->Q,i,j)= ElementM(expanded_sv->state_variable[k]->Q,i,j-1); */
/*       Update_B_from_Q_SV(expanded_sv->state_variable[k]); */
/*  //=====================================================================================================================   ansi-c*/
/*        // Draw states for restricted model   ansi-c*/
      DrawStates(restricted_model);

/*        // Copy states for kth state variable   ansi-c*/
      dw_CopyArray(expanded_sv->state_variable[k]->S,restricted_sv->state_variable[k]->S);

/*        //dw_PrintArray(stdout,expanded_sv->state_variable[k]->S,(char*)NULL); getchar();   ansi-c*/

/*        // Draw transition matrix for kth state variable   ansi-c*/
      DrawTransitionMatrix_SV(expanded_sv->state_variable[k]);
/*  //=====================================================================================================================   ansi-c*/
/*        // Then new state for the sth state variable is made to be reflecting   ansi-c*/
/*       j=restricted_sv->state_variable[k]->nstates; */
/*       InsertSubMatrix(expanded_sv->state_variable[k]->Q,restricted_sv->state_variable[k]->Q,0,0,0,0,j,j); */
/*       for (i=j; i >= 0; i--) */
/*         { */
/*           ElementM(expanded_sv->state_variable[k]->Q,i,j)=0.0; */
/*           ElementM(expanded_sv->state_variable[k]->Q,j,i)=0.0; */
/*         } */
/*       ElementM(expanded_sv->state_variable[k]->Q,j-1,j)=1.0; */
/*       Update_B_from_Q_SV(expanded_sv->state_variable[k]); */
    }
    }
  PropagateTransitionMatrices_SV(expanded_sv);
/*    //  expanded_model->ValidTransitionMatrix=1;   ansi-c*/
  ValidateTransitionMatrices_SV(expanded_sv);
/*    // return   ansi-c*/
  return expanded_model;
}




/*
   The free transition matrix parameters are chosen to be equal.
*/
void FlatTransitionMatrix(TMarkovStateVariable *sv)
{
  int i, j;
  PRECISION p;
  if (sv->n_state_variables == 1)
    {
      for (i=dw_DimA(sv->b)-1; i >= 0; i--)
    for (p=1.0/(PRECISION)DimV(sv->b[i]), j=DimV(sv->b[i])-1; j >= 0; j--)
      ElementV(sv->b[i],j)=p;
      Update_Q_from_B_SV(sv);
    }
  else
    {
      for (i=sv->n_state_variables-1; i >= 0; i--)
    FlatTransitionMatrix(sv->state_variable[i]);
      MatrixTensor(sv->Q,sv->QA);
    }
}

/*

*/
int NestTransitionMatrices_SV(TMarkovStateVariable *sv, TMarkovStateVariable *restricted_sv)
{
  int i, j;
  PRECISION tmp;
  if (sv->n_state_variables < restricted_sv->n_state_variables) return 0;
  if (sv->n_state_variables == 1)
    if (sv->nstates < restricted_sv->nstates)
      return 0;
    else
      {
    for (i=0; i < restricted_sv->nstates-1; i++)
      {
        for (j=0; j < restricted_sv->nstates; j++)
          ElementM(sv->Q,i,j)=ElementM(restricted_sv->Q,i,j);
        for (tmp=ElementM(restricted_sv->Q,i,j-1); j < sv->nstates; j++)
          ElementM(sv->Q,i,j)=tmp;
      }

        for (j=0; j < restricted_sv->nstates; j++)
      {
        tmp=ElementM(restricted_sv->Q,restricted_sv->nstates-1,j)/(double)(sv->nstates - restricted_sv->nstates + 1);
/*          // tmp=0.0;   ansi-c*/
        for (i=restricted_sv->nstates-1; i < sv->nstates; i++)
          ElementM(sv->Q,i,j)=tmp;
      }
    for ( ; j < sv->nstates; j++)
      for (i=restricted_sv->nstates-1; i < sv->nstates; i++)
        ElementM(sv->Q,i,j)=tmp;

    return Update_B_from_Q_SV(sv);
      }
  else
    {
      for (i=sv->n_state_variables-1; i >= restricted_sv->n_state_variables; i--)
    FlatTransitionMatrix(sv->state_variable[i]);
      for ( ; i >= 0; i--)
    if (!NestTransitionMatrices_SV(sv->state_variable[i],restricted_sv->state_variable[i]))
      return 0;
      MatrixTensor(sv->Q,sv->QA);
    }
  return 1;
}

/*
   Sets the parameters of model so that it is equivalent to the parameters of
   restricted_model.  Currently, the routine only checks that the sizes of the
   two models are the same.  It should be the case that the restrictions are
   also identical.
*
int NestModel_VAR(TStateModel *model, TStateModel *restricted_model)
{
  int j, k;
  T_VAR_Parameters *p=model->theta, *restricted_p=restricted_model->theta;

  // Set transition matrices
  if (!restricted_model->ValidTransitionMatrix || !NestTransitionMatrices_SV(model->sv,restricted_model->sv))
    return 0;
  else
    model->ValidTransitionMatrix=1;

  // Check VAR sizes
  if ((p->nvars != restricted_p->nvars) || (p->nlags != restricted_p->nlags)
      || (p->npre != restricted_p->npre) || (p->nobs != restricted_p->nobs))
    return 0;

  // Set VAR parameters
  for (j=0; j < p->nvars; j++)
    {
      if (p->n_var_states[j] < restricted_p->n_var_states[j]) return 0;
      for (k=0; k < restricted_p->n_var_states[j]; k++)
    p->Zeta[j][k]=restricted_p->Zeta[j][k];
      for ( ; k < p->n_var_states[j]; k++)
    p->Zeta[j][k]=restricted_p->Zeta[j][restricted_p->n_var_states[j]-1];

      if (p->n_coef_states[j] < restricted_p->n_coef_states[j]) return 0;
      for (k=0; k < restricted_p->n_coef_states[j]; k++)
    {
      EquateVector(p->A0[j][k],restricted_p->A0[j][k]);
      EquateVector(p->Aplus[j][k],restricted_p->Aplus[j][k]);
    }
      for ( ; k < p->n_coef_states[j]; k++)
    {
      EquateVector(p->A0[j][k],restricted_p->A0[j][restricted_p->n_coef_states[j]-1]);
      EquateVector(p->Aplus[j][k],restricted_p->Aplus[j][restricted_p->n_coef_states[j]-1]);
    }
    }
  Update_b0_bplus_from_A0_Aplus(model->theta);
  return 1;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/*************************** Workspace computations ****************************/
/*******************************************************************************/
void ComputeDotProducts_All(T_VAR_Parameters *p)
{
  TVector *Y=p->Y, *X=p->X;
  int t, j, k;

  for (t=p->nobs; t > 0; t--)
    for (j=p->nvars-1; j >= 0; j--)
      for (k=p->n_coef_states[j]-1; k >= 0; k--)
    {
      p->A0_dot_products[t][j][k]=DotProduct(Y[t],p->A0[j][k]);
      p->Aplus_dot_products[t][j][k]=DotProduct(X[t],p->Aplus[j][k]);
    }

  p->valid_dot_products=1;

}

/*
   Assumes:
     p:  A valid T_VAR_Parameters structure and p->A0 has been initialized.

   Results:
     Fill the vector p->log_abs_det_A0 with the natural logarithm of the the
     abso
*/
void ComputeLogAbsDetA0_All(T_VAR_Parameters *p)
{
  TMatrix A0;
  int j, k=DimV(p->log_abs_det_A0)-1;

  A0=CreateMatrix(p->nvars,p->nvars);

/*    //=== Set initial A0 ===   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    memcpy(&ElementM(A0,0,j),pElementV(p->A0[j][p->A0_column_states[j][k]]),p->nvars*sizeof(PRECISION));

  ElementV(p->log_abs_det_A0,k)=LogAbsDeterminant_LU(A0);

  while (--k >= 0)
    {
/*        //=== Reset A0 ===   ansi-c*/
      for (j=p->nvars-1; j >= 0; j--)
    if (p->A0_column_states[j][k] != p->A0_column_states[j][k+1])
      memcpy(&ElementM(A0,0,j),pElementV(p->A0[j][p->A0_column_states[j][k]]),p->nvars*sizeof(PRECISION));

      ElementV(p->log_abs_det_A0,k)=LogAbsDeterminant_LU(A0);
    }

  FreeMatrix(A0);

  p->valid_log_abs_det_A0=1;
}

/*
   Computes the log of the absolute value of the determinant of A0(s) if
   A0_column_states[j][s] == k.
*/
void ComputeLogAbsDetA0(int j, int k, T_VAR_Parameters *p)
{
  TMatrix A0;
  int i, s;

  A0=CreateMatrix(p->nvars,p->nvars);

  for (s=DimV(p->log_abs_det_A0)-1; s >= 0; s--)
    if (p->A0_column_states[j][s] == k)
      {
    for (i=p->nvars-1; i >= 0; i--)
      memcpy(&ElementM(A0,0,i),pElementV(p->A0[i][p->A0_column_states[i][s]]),p->nvars*sizeof(PRECISION));

    ElementV(p->log_abs_det_A0,s)=LogAbsDeterminant_LU(A0);
      }

  FreeMatrix(A0);
}

/*
   If p->A_state_variable is non-negative, then it controls A0 and Aplus.  If
   p->A_state_variable is negative, then A0 and Aplus are constant across states.

   If Xi is not null and p->Xi_state_variable is non-negative, then it controls
   Xi.  If Xi is not null and p->Xi_state_variable is negative, then Xi is
   constant across states.  If Xi is null, then this parameter does not appear.

      -n*log(2*pi)/2 + log|A0[i]| + log|diag(Xi[j])|
            - (Y[t]'*A0[i] - X[t]'*Aplus[i])*diag(Xi[j])*diag(Xi[j])*(Y[t]'*A0[i] - X[t]'*Aplus[i])/2
*/
PRECISION LogConditionalProbability_VAR(int s, int t, TStateModel *model)
{
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);
  PRECISION  x, y, sum=0.0, logdet=0.0;
  int  j;

  if (!(p->valid_parameters)) return MINUS_INFINITY;

/*    //====== Computes log(abs(det(A0[i]))) and log(abs(det(Xi))) ======   ansi-c*/
  if (!p->valid_log_abs_det_A0) ComputeLogAbsDetA0_All(p);

/*    //====== Compute quadratic form ======   ansi-c*/
  if (p->valid_dot_products)
    for (j=p->nvars-1; j >= 0; j--)
      {
        y=p->Zeta[j][p->var_states[j][s]];
        x=p->A0_dot_products[t][j][p->coef_states[j][s]]
                       - p->Aplus_dot_products[t][j][p->coef_states[j][s]];
        if (y <= 0)
          return MINUS_INFINITY;
        else
          logdet+=log(y);
        sum+=y*x*x;
      }
   else
     for (j=p->nvars-1; j >= 0; j--)
      {
        y=p->Zeta[j][p->var_states[j][s]];
        x=DotProduct(p->Y[t],p->A0[j][p->coef_states[j][s]])
                           - DotProduct(p->X[t],p->Aplus[j][p->coef_states[j][s]]);
        if (y <= 0)
          return MINUS_INFINITY;
        else
          logdet+=log(y);
        sum+=y*x*x;
      }

/*    //====== Get log(det(A0)) ======   ansi-c*/
  logdet=0.5*logdet + ElementV(p->log_abs_det_A0,p->A0_states[s]);

  return p->minus_half_nvars_times_log2pi + logdet - 0.5*sum;
}

/*
    Since

      y(t)' * A0(s(t)) = x(t)' * Aplus(s(t)) + epsilon(t)' * Inverse(Xi(s(t)))

    The expectation is

      y(t)' = x(t)' * Aplus(s(t)) * Inverse(A0(s(t))
*/
TVector ExpectationSingleStep_VAR(TVector y, int s, int t, TStateModel *model)
{
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);
  TMatrix A0, Aplus;

  if ((t < 1) || (p->nobs < t))
    {
      dw_Error(SIZE_ERR);
      return (TVector)NULL;
    }

  A0=MakeA0((TMatrix)NULL,s,p);
  Aplus=MakeAplus((TMatrix)NULL,s,p);

  if (!y && !(y=CreateVector(p->nvars)))
    return (TVector)y;

  ProductVM(y,p->X[t],Aplus);
  ProductInverseVM(y,y,A0);

  FreeMatrix(Aplus);
  FreeMatrix(A0);

  return y;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/******************************** Notification *********************************/
void StatesChanged_VAR(TStateModel *model)
{
 ((T_VAR_Parameters*)(model->theta))->valid_state_dependent_fields=0;
}

void ThetaChanged_VAR(TStateModel *model)
{
  ((T_VAR_Parameters*)(model->theta))->valid_log_abs_det_A0
      =((T_VAR_Parameters*)(model->theta))->valid_dot_products=0;
}

void InitializeForwardRecursion_VAR(TStateModel *model)
{
  if (!((T_VAR_Parameters*)(model->theta))->valid_dot_products)
    ComputeDotProducts_All((T_VAR_Parameters*)(model->theta));
}

/*******************************************************************************/
/********************************* Simulations *********************************/
/*******************************************************************************/
void DrawParameters_VAR(TStateModel *model)
{
/*    // Draw unnormalized theta   ansi-c*/
  DrawZeta_DotProducts(model);
  DrawA0_Metropolis(model);
  DrawAplus(model);

/*    // Normalize   ansi-c*/
  Normalize_VAR((T_VAR_Parameters*)(model->theta));

/*    // Flags and notification that the VAR parameters have changed   ansi-c*/
  ((T_VAR_Parameters*)(model->theta))->valid_parameters=1;
  ThetaChanged(model);
}

/*
   Choose a random initial value for the VAR parameters.  The function
   ThetaChanged() cannot be called.
*/
void InitializeParameters_VAR(T_VAR_Parameters *p)
{
  int j, s;
  TMatrix X;

/*    // Initialize Zeta to one   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    for (s=p->n_var_states[j]-1; s >= 0; s--)
      p->Zeta[j][s]=1;

/*    // Draw b0 from prior - identical across states A0[j][s] = U[j]*b0[j][s].   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    {
      X=CholeskyUT((TMatrix)NULL,p->inverse_b0_prior[j]);
      Inverse_UT(X,X);
      dw_NormalVector(p->b0[j][0]);
      ProductMV(p->b0[j][0],X,p->b0[j][0]);
      ProductMV(p->A0[j][0],p->U[j],p->b0[j][0]);
      for (s=p->n_coef_states[j]-1; s > 0; s--)
    EquateVector(p->A0[j][s],p->A0[j][0]);
      FreeMatrix(X);
    }

/*    // Set Aplus[j][s] = W[j]*A0[j][s].   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    {
      if (!p->W[j])
    InitializeVector(p->Aplus[j][0],0.0);
      else
    ProductMV(p->Aplus[j][0],p->W[j],p->A0[j][0]);
      for (s=p->n_coef_states[j]-1; s > 0; s--)
    EquateVector(p->Aplus[j][s],p->Aplus[j][0]);
    }

/*    // Update b0, bplus, lambda, psi   ansi-c*/
  Update_b0_bplus_from_A0_Aplus(p);
  if ((p->Specification & SPEC_SIMS_ZHA) == SPEC_SIMS_ZHA) Update_lambda_psi_from_bplus(p);

/*    // Flags and notification that the VAR parameters have changed   ansi-c*/
  p->valid_parameters=1;
}

/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*  //--------------------------------- Draw Zeta ---------------------------------//   ansi-c*/
/*  //-----------------------------------------------------------------------------//   ansi-c*/
void DrawZeta_Aplus(TStateModel *model)
{
  int j, k, s, T;
  PRECISION v;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  if (!p->valid_state_dependent_fields) UpdateStateDependentFields(p,model->sv->S);

  for (j=p->nvars-1; j >= 0; j--)
    {
      for (s=p->n_var_states[j]-1; s > 0; s--)
    {
      v=0.0;
      T=0;
      for (k=p->nstates-1; k >= 0; k--)
        {
          if (p->var_states[j][k] == s)
        {
          v+=InnerProductSymmetric(p->A0[j][p->coef_states[j][k]],p->YY[k])
            - 2.0 * InnerProductNonSymmetric(p->Aplus[j][p->coef_states[j][k]],p->A0[j][p->coef_states[j][k]],p->XY[k])
            + InnerProductSymmetric(p->Aplus[j][p->coef_states[j][k]],p->XX[k]);
          T+=p->T[k];
        }
        }
      p->Zeta[j][s]=dw_gamma_rnd(0.5*(PRECISION)T + ElementV(p->zeta_a_prior,j))/(0.5*v + ElementV(p->inverse_zeta_b_prior,j));
    }

/*        //=== State 0 is normalized to one   ansi-c*/
      p->Zeta[j][0]=1.0;
    }
}

void DrawZeta_DotProducts(TStateModel *model)
{
  int j, s, t;
  PRECISION x;
  TVector v, T;
  int *S=model->sv->S;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);

  if (!p->valid_dot_products) ComputeDotProducts_All(p);

  for (j=p->nvars-1; j >= 0; j--)
    {
      v=InitializeVector(CreateVector(p->n_var_states[j]),0.0);
      T=InitializeVector(CreateVector(p->n_var_states[j]),0.0);
      for (t=p->nobs; t > 0; t--)
    {
      s=p->coef_states[j][S[t]];
      x=p->A0_dot_products[t][j][s] - p->Aplus_dot_products[t][j][s];
      s=p->var_states[j][S[t]];
      ElementV(v,s)+=x*x;
      ElementV(T,s)+=1.0;
    }
      for (s=p->n_var_states[j]-1; s > 0; s--)
        p->Zeta[j][s]=dw_gamma_rnd(0.5*ElementV(T,s) + ElementV(p->zeta_a_prior,j))/(0.5*ElementV(v,s) + ElementV(p->inverse_zeta_b_prior,j));
      FreeVector(v);
      FreeVector(T);

/*        //=== State 0 is normalized to one   ansi-c*/
      p->Zeta[j][0]=1.0;
    }
}

/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*  //-------------------------- Metropolis Draws of A0 ---------------------------//   ansi-c*/
/*  //-----------------------------------------------------------------------------//   ansi-c*/
#define MID          0.35
#define LOG_MID     -1.0498221244987
#define LOWER_BOUND  0.0052521875
#define UPPER_BOUND  0.81061308309895
void AdaptiveMetropolisScale(TStateModel *model, int iterations, int period, int verbose, FILE *f_posterior)
{
  struct TAdaptive
  {
    int begin_jump_ratio;
    int iterations;
    int end_iteration_count;
    PRECISION best_scale;
    PRECISION low_scale;
    PRECISION low_jump_ratio;
    PRECISION high_scale;
    PRECISION high_jump_ratio;
  } ***Adaptive;

  PRECISION new_scale, new_jump_ratio;
  int j, k, begin_time, count, check=period;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  Adaptive=(struct TAdaptive***)dw_CreateArray_array(p->nvars);
  for (j=p->nvars-1; j >= 0; j--)
    {
      Adaptive[j]=(struct TAdaptive**)dw_CreateArray_pointer(p->n_coef_states[j],swzFree);
      for (k=p->n_coef_states[j]-1; k >= 0; k--)
    {
      Adaptive[j][k]=((struct TAdaptive*)swzMalloc(sizeof(struct TAdaptive)));
      Adaptive[j][k]->begin_jump_ratio=p->A0_Metropolis_Jumps[j][k];
      Adaptive[j][k]->iterations=period;
      Adaptive[j][k]->end_iteration_count=period;
      Adaptive[j][k]->low_scale=Adaptive[j][k]->low_jump_ratio=Adaptive[j][k]->high_scale=Adaptive[j][k]->high_jump_ratio=-1.0;
      Adaptive[j][k]->best_scale=p->A0_Metropolis_Scale[j][k];
    }
    }

  ResetMetropolisInformation(p);

  if (verbose)
    {
      printf("Beginning adaptive burn in -- %d iterations.\n",iterations);
      begin_time=(int)time((time_t*)NULL);
    }

  for (count=1; count <= iterations; count++)
    {
      DrawAll(model);
      if (count == check)
    {
      if (f_posterior) fprintf(f_posterior,"%le\n",LogPosterior_StatesIntegratedOut(model));

          if (verbose)
        printf("%d iterations completed out of %d - elapsed time: %d seconds\n",count,iterations,(int)time((time_t*)NULL) - begin_time);

      for (j=p->nvars-1; j >= 0; j--)
        {
          for (k=p->n_coef_states[j]-1; k >= 0; k--)
        if (Adaptive[j][k]->end_iteration_count == count)
          {
/*              // Compute new jump ratio and get scale   ansi-c*/
            new_jump_ratio=(PRECISION)(p->A0_Metropolis_Jumps[j][k] - Adaptive[j][k]->begin_jump_ratio)
                                                                                 /(PRECISION)(Adaptive[j][k]->iterations);

/*              // Set new low or high bounds   ansi-c*/
            if (new_jump_ratio < MID)
              {
            Adaptive[j][k]->low_scale=p->A0_Metropolis_Scale[j][k];
            Adaptive[j][k]->low_jump_ratio=new_jump_ratio;
              }
            else
              {
            Adaptive[j][k]->high_scale=p->A0_Metropolis_Scale[j][k];
            Adaptive[j][k]->high_jump_ratio=new_jump_ratio;
              }

/*              // Compute new scale and best scale   ansi-c*/
            if (Adaptive[j][k]->low_jump_ratio < 0.0)
              {
            Adaptive[j][k]->best_scale=Adaptive[j][k]->high_scale;
            if (Adaptive[j][k]->low_scale < 0.0)
              new_scale=((new_jump_ratio > UPPER_BOUND) ? 5.0 : LOG_MID/log(new_jump_ratio))*Adaptive[j][k]->high_scale;
            else
              {
                new_scale=Adaptive[j][k]->low_scale;
                Adaptive[j][k]->low_scale=-1;
              }
              }
            else
              if (Adaptive[j][k]->high_jump_ratio < 0.0)
            {
              Adaptive[j][k]->best_scale=Adaptive[j][k]->low_scale;
              if (Adaptive[j][k]->high_scale < 0.0)
                new_scale=((new_jump_ratio < LOWER_BOUND) ? 0.2 : LOG_MID/log(new_jump_ratio))*Adaptive[j][k]->low_scale;
              else
                {
                  new_scale=Adaptive[j][k]->high_scale;
                  Adaptive[j][k]->high_scale=-1.0;
                }
            }
              else
            {
              new_scale=Adaptive[j][k]->best_scale=0.5*(Adaptive[j][k]->low_scale + Adaptive[j][k]->high_scale);
/*                //Adaptive[j][k]->iterations+=period;   ansi-c*/
              Adaptive[j][k]->iterations*=2;
              Adaptive[j][k]->low_jump_ratio=Adaptive[j][k]->high_jump_ratio=-1.0;
            }

/*              // Print data   ansi-c*/
            if (verbose)
              printf("col: %d  state: %d  (%d %lf %lf %lf)\n",j+1,k+1,p->A0_Metropolis_Jumps[j][k],
                                                                   new_jump_ratio,p->A0_Metropolis_Scale[j][k],new_scale);

/*              // Reset adaptive counts and A0_Metropolis_Scale   ansi-c*/
            Adaptive[j][k]->begin_jump_ratio=p->A0_Metropolis_Jumps[j][k];
            Adaptive[j][k]->end_iteration_count+=Adaptive[j][k]->iterations;
            p->A0_Metropolis_Scale[j][k]=new_scale;
          }
        else
          if (verbose)
            {
              new_jump_ratio=(PRECISION)(p->A0_Metropolis_Jumps[j][k] - Adaptive[j][k]->begin_jump_ratio)
                                             /(PRECISION)(Adaptive[j][k]->iterations - (Adaptive[j][k]->end_iteration_count - count));

              printf("col: %d  state: %d  (%d %lf %lf -)\n",j+1,k+1,p->A0_Metropolis_Jumps[j][k],
                                                                            new_jump_ratio,p->A0_Metropolis_Scale[j][k]);
            }
        }

      if (verbose) printf("\n");

      check+=period;
    }
    }

  for (j=p->nvars-1; j >= 0; j--)
    for (k=p->n_coef_states[j]-1; k >= 0; k--)
      p->A0_Metropolis_Scale[j][k]=Adaptive[j][k]->best_scale;

  ResetMetropolisInformation(p);

  dw_FreeArray(Adaptive);
}
#undef MID
#undef UPPER
#undef LOWER

void SetupMetropolisInformation(PRECISION **Scale, T_VAR_Parameters *p)
{
  dw_CopyArray(p->A0_Metropolis_Scale,Scale);
  ResetMetropolisInformation(p);
}

void ResetMetropolisInformation(T_VAR_Parameters *p)
{
  p->Total_A0_Metropolis_Draws=0;
  dw_InitializeArray_int(p->A0_Metropolis_Jumps,0);
}

static void GetProposedJump_A0(TVector b, int j, int k, T_VAR_Parameters *p)
{
  TMatrix YY, XX, XY, S, M0, M1;
  int s;
  PRECISION x;
  int terminal_errors;

/*    // Accumulate XX, XY, and YY   ansi-c*/
  InitializeMatrix(XX=CreateMatrix(p->npre,p->npre),0.0);
  InitializeMatrix(XY=CreateMatrix(p->npre,p->nvars),0.0);
  InitializeMatrix(YY=CreateMatrix(p->nvars,p->nvars),0.0);
  for (s=p->nstates-1; s >= 0; s--)
    if (p->coef_states[j][s] == k)
      {
    x=p->Zeta[j][p->var_states[j][s]];
    UpdateMS(XX,p->XX[s],x);
    UpdateMS(XY,p->XY[s],x);
    UpdateMS(YY,p->YY[s],x);
      }

/*    // S = inverse_b0_prior + U[j]'*(YY + W[j]'*XY + XY'*W[j] + W[j]'*XX*W[j])*U[j]   ansi-c*/
  if (p->W[j])
    {
/*        // M0 = W[j]'*(XX*W[j] + XY) + (XY'*W[j])   ansi-c*/
      M0=ProductMM((TMatrix)NULL,XX,p->W[j]);
      AddMM(M0,M0,XY);
      M1=TransposeProductMM((TMatrix)NULL,p->W[j],M0);
      AddMM(YY,YY,M1);
      TransposeProductMM(M1,XY,p->W[j]);
      AddMM(YY,YY,M1);
      FreeMatrix(M1);
      FreeMatrix(M0);
    }
  S=MatrixInnerProductSymmetric((TMatrix)NULL,p->U[j],YY);
/*    //dw_PrintMatrix(stdout,S,"%.17le "); printf("\n");   ansi-c*/
  AddMM(S,S,p->inverse_b0_prior[j]);

/*    //dw_PrintMatrix(stdout,S,"%.17le "); fgetc(stdin);   ansi-c*/

/*    // Simulate draw   ansi-c*/
  terminal_errors=dw_SetTerminalErrors(NO_ERR);
  dw_NormalVector(b);
  if (!InverseProductUV(b,CholeskyUT(S,S),b))
    {
      printf("Error in GetProposedJump_A0()\n");
      printf("j = %d, k = %d\n,Prior =\n",j,k);
      dw_PrintMatrix(stdout,p->inverse_b0_prior[j],"%lg ");
      printf("S =\n");
      dw_PrintMatrix(stdout,S,"%lg ");
      swzExit(1);
    }
  dw_SetTerminalErrors(terminal_errors);
/*   else */
/*     { */
/*       printf("GetProposedJump_A0()\n"); */
/*       printf("j = %d, k = %d\n,Prior =\n",j,k); */
/*       dw_PrintMatrix(stdout,p->inverse_b0_prior[j],"%lg "); */
/*       printf("S =\n"); */
/*       dw_PrintMatrix(stdout,S,"%lg "); */
/*       getchar(); */
/*     } */

/*    // Scale factor   ansi-c*/
  ProductVS(b,b,p->A0_Metropolis_Scale[j][k]);

/*    // Free memory   ansi-c*/
  FreeMatrix(S);
  FreeMatrix(YY);
  FreeMatrix(XY);
  FreeMatrix(XX);
}

/*
   Assumes
     j column
     k state A0[j], 0 <= k < dw_DimA(A[j])
     p pointer to valid T_VAR_Parameters structure

   Returns
     The log of the density of a[j][k] conditional on all other parameters.

   Notes
     Uses the following information in model->parametes = p and model->data = d.

        p->XX[], p->XY[], p->YY, p->nobs_by_state[], p->log_det_abs_A0[]

        p->A0[j][k], p->b0[j][k], p->Aplus[j][k], p->Xi[j][]

        p->coef_states[j][], p->var_states[j][], p->A0_states[],
        p->inverse_b0_prior[j]
*/
PRECISION LogKernel_A0(int j, int k, TStateModel *model)
{
  int s;
  PRECISION rtrn;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);

/*    //====== Prior ======   ansi-c*/
  rtrn=-0.5*InnerProductSymmetric(p->b0[j][k],p->inverse_b0_prior[j]);

  for (s=p->nstates-1; s >= 0; s--)
    if (p->coef_states[j][s] == k)
      rtrn+=ElementV(p->log_abs_det_A0,p->A0_states[s]) * p->T[s]
                   - 0.5 * p->Zeta[j][p->var_states[j][s]] * (InnerProductSymmetric(p->A0[j][k],p->YY[s])
                             - 2.0*InnerProductNonSymmetric(p->Aplus[j][k],p->A0[j][k],p->XY[s])
                          + InnerProductSymmetric(p->Aplus[j][k],p->XX[s]));
  return rtrn;
}

PRECISION LogKernel_A0_DotProduct(int j, int k, TStateModel *model)
{
  int s, t;
  int* S=model->sv->S;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);
  PRECISION  x, rtrn;

/*    //====== Prior ======   ansi-c*/
  rtrn=-0.5*InnerProductSymmetric(p->b0[j][k],p->inverse_b0_prior[j]);

  for (t=model->sv->nobs; t > 0; t--)
    if (p->coef_states[j][s=S[t]] == k)
      {
    x=DotProduct(p->Y[t],p->A0[j][k]) - DotProduct(p->X[t],p->Aplus[j][k]);
    rtrn+=ElementV(p->log_abs_det_A0,p->A0_states[s]) - 0.5*p->Zeta[j][p->var_states[j][s]]*x*x;
      }

  return rtrn;
}

/*
   Assumes:
     p - pointer to a valid T_VAR_Parameter structure
     Jumps - null pointer or pointer to a 2-dimensional integer array with the
             first dimensional at lease p->nvars and the second dimension at
             least the dimension of p->A0[j].

   Results:
     New values for p->b0 and p->A0 are obtained using the Metropolis algorithm.
     If Jumps is not null, then Jump[j][k] is updated for all the Metropolis
     jumps that were accepted.

   Notes:
     Calls GetA0MetropolisJumpSize() to get the variance for the normal jumping
     kernel used in this algorithm.  Calls LogConditionalDensity_A0() to compute
     the appropriate posterior density.
*/
void DrawA0_Metropolis(TStateModel *model)
{
  int j, k;
  PRECISION old_log_kernel, log_difference;
  TVector old_b0, old_a0, old_aplus, old_log_abs_det_A0;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  if (!p->valid_state_dependent_fields) UpdateStateDependentFields(p,model->sv->S);

  old_a0=CreateVector(p->nvars);
  old_aplus=CreateVector(p->npre);
  old_log_abs_det_A0=CreateVector(DimV(p->log_abs_det_A0));

  for (j=p->nvars-1; j >= 0; j--)
    {
      old_b0=CreateVector(DimV(p->b0[j][0]));

      for (k=dw_DimA(p->A0[j])-1; k >= 0; k--)
    {
/*        //=== Save old values ===   ansi-c*/
      EquateVector(old_b0,p->b0[j][k]);
      EquateVector(old_a0,p->A0[j][k]);
      EquateVector(old_aplus,p->Aplus[j][k]);
      EquateVector(old_log_abs_det_A0,p->log_abs_det_A0);
      old_log_kernel=LogKernel_A0(j,k,model);
/*        //old_log_kernel=LogKernel_A0_DotProduct(j,k,model);   ansi-c*/

/*        //=== Jump ===   ansi-c*/
      GetProposedJump_A0(p->b0[j][k],j,k,p);
      AddVV(p->b0[j][k],p->b0[j][k],old_b0);
      ProductMV(p->A0[j][k],p->U[j],p->b0[j][k]);
      Update_aplus_from_bplus_a0(j,k,p);
      ComputeLogAbsDetA0(j,k,p);

/*        //=== Accept Jump ===   ansi-c*/
      log_difference=LogKernel_A0(j,k,model) - old_log_kernel;
/*        //log_difference=LogKernel_A0_DotProduct(j,k,model) - old_log_kernel;   ansi-c*/
      if ((log_difference >= 0.0) || (dw_uniform_rnd() < exp(log_difference)))
        {
          p->A0_Metropolis_Jumps[j][k]++;
          p->valid_dot_products=0;
        }
      else
        {
          EquateVector(p->b0[j][k],old_b0);
          EquateVector(p->A0[j][k],old_a0);
          EquateVector(p->Aplus[j][k],old_aplus);
          EquateVector(p->log_abs_det_A0,old_log_abs_det_A0);
        }
    }

      FreeVector(old_b0);
    }

  FreeVector(old_log_abs_det_A0);
  FreeVector(old_aplus);
  FreeVector(old_a0);

  p->Total_A0_Metropolis_Draws++;
}

/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*  //-------------------------------- Draw Aplus ---------------------------------//   ansi-c*/
/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*
   The following matrices must be updated before calling this routine:

     p->A0  p->Xi  p->XX   p->XY

   The following matrices are used in this routine

     p->V  p->W  p->coef_states  p->var_states p->inverse_bplus_prior

   The following matrices are modified in this routine

     p->bplus  p->Aplus
*/
void DrawAplus(TStateModel *model)
{
  int j, k, s;
  TMatrix S, XX, XY, M;
  TVector v;
  PRECISION x;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);

  if (!p->valid_state_dependent_fields)
    UpdateStateDependentFields(p,model->sv->S);

  if (p->Specification & SPEC_SIMS_ZHA)
    {
      Draw_psi(model);
      Draw_lambda(model);
      Update_bplus_from_lambda_psi(p);
    }
  else
    {
      XX=CreateMatrix(p->npre,p->npre);
      XY=CreateMatrix(p->npre,p->nvars);
      v=CreateVector(p->npre);

      for (j=p->nvars-1; j >= 0; j--)
    if (p->bplus[j])
      {
        S=CreateMatrix(p->dim_bplus[j],p->dim_bplus[j]);

        for (k=p->n_coef_states[j]-1; k >= 0; k--)
          {
        InitializeMatrix(XX,0.0);
        InitializeMatrix(XY,0.0);
        for (s=p->nstates-1; s >= 0; s--)
          if (p->coef_states[j][s] == k)
            {
              x=p->Zeta[j][p->var_states[j][s]];
              UpdateMS(XX,p->XX[s],x);
              UpdateMS(XY,p->XY[s],x);
            }

/*          //=== Compute inverse variance ===   ansi-c*/
        if (!p->IsIdentity_V[j])
          MatrixInnerProductSymmetric(S,p->V[j],XX);
        else
          EquateMatrix(S,XX);
        AddMM(S,S,p->inverse_bplus_prior[j]);

/*          //=== Compute b ===   ansi-c*/
        if (p->W[j])
          if (p->Specification & SPEC_RANDOM_WALK)
            if (MajorForm(XY) && MajorForm(XX))
              bSubtract(pElementM(XY),pElementM(XY),pElementM(XX),RowM(XX)*p->nvars);
            else
              {
            M=SubMatrix((TMatrix)NULL,XX,0,0,RowM(XX),p->nvars);
            SubtractMM(XY,XY,M);
            FreeMatrix(M);
              }
          else
            {
              M=ProductMM((TMatrix)NULL,XX,p->W[j]);
              AddMM(XY,XY,M);
              FreeMatrix(M);
            }
        if (!p->IsIdentity_V[j])
          {
            ProductMV(v,XY,p->A0[j][k]);
            TransposeProductMV(p->bplus[j][k],p->V[j],v);
          }
        else
          ProductMV(p->bplus[j][k],XY,p->A0[j][k]);

/*          //=== Draw bplus ===   ansi-c*/
        if (!DrawNormal_InverseVariance(p->bplus[j][k],p->bplus[j][k],S))
          SingularInverseVariance_RecoveryAttempt(p->bplus[j][k],p->bplus[j][k],S,p->inverse_bplus_prior[j],model,BPLUS_ERR);
          }

        FreeMatrix(S);
      }

      FreeMatrix(XX);
      FreeMatrix(XY);
      FreeVector(v);
    }

  Update_Aplus_from_bplus_A0(p);
  ThetaChanged(model);
}

/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*  //--------------------------------- Draw psi ----------------------------------//   ansi-c*/
/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*
   Assumes:
     S:  (n*b + j) x (n*b + j) matrix with j > k.  S must be column major.
     XX:  (n*b + 1) x (n*b + 1) matrix.  XX must be column major.
     lambda:  n dimensional vector.

   Results:
     Adds

         PSI[k]'*diag(LAMBDA*lambda + e)*XX*diag(LAMBDA*lambda + e)*PSI[k]

     to S.  The matrices LAMBDA and PSI[k] are given by

                        -      -
                       |  I(n)  |
                       |   .    |                  -                  -
                       |   .    |                 |  I(n*b)   0(n*b,j) |
              LAMBDA = |   .    |        PSI[k] = |                    |
                       |        |                 | 0(1,n*b)  e(j,k+1)'|
                       |  I(n)  |                  -                  -
                       | 0(1,n) |
                        -      -

     where I(r) is the r x r identity matrix, 0(r,s) is the r x s zero matrix,
     and e(r,s) is the sth column of I(r).  e is the vector e(n*b+1,n*b+1).

   Notes:
     k is a zero based index.
*/
void update_psi_quadratic_form(TMatrix S, int n, int m, int k, TVector lambda, TMatrix XX)
{
  int i, j, p, u=n*m-1, v;
  PRECISION *x, *s, *z=pElementV(lambda), w;

  s=pElementM(S)+(n*m+k)*RowM(S)+u;
  x=pElementM(XX)+n*m*RowM(XX)+u;
  s[k+1]=x[1];
  for (v=m-1; v >= 0; v--)
    for (i=n-1; i >= 0; s--, x--, i--)
      (*s)+=z[i]*(*x);
  for (p=n-1, j=u; j >= 0; p--, j--)
    {
      if (p < 0) p=n-1;
      w=z[p];
      s=pElementM(S)+j*RowM(S)+u;
      x=pElementM(XX)+j*RowM(XX)+u;
      s[k+1]=x[1]*w;
      for (v=m-1; v >= 0; v--)
    for (i=n-1; i >= 0; s--, x--, i--)
      (*s)+=z[i]*(*x)*w;
    }
}

/*
   Assmes:
     model:  point to a valid TStateModel structure

   Results:
     A draw of psi conditional on A0,

   Notes:
     The matrices MUST be in column major format.  Basic matrix routines from
     bmatrix.c are called in this routine.
*/
void Draw_psi(TStateModel *model)
{
  int j, k, s, i, m;
  TVector b, v;
  TMatrix S, XX, XY;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);

/*    // Update state dependent matrices if necessary   ansi-c*/
  if (!p->valid_state_dependent_fields)
    UpdateStateDependentFields(p,model->sv->S);

/*    // Allocate memory   ansi-c*/
  XX=CreateMatrix(p->npre,p->npre);
  XY=CreateMatrix(p->npre,p->nvars);
  v=CreateVector(p->npre);

  if (_VERBOSE_COUNT)
    fprintf(V_FILE,"//=== Draw_psi() (count = %d) =====================================//\n",_VERBOSE_COUNT);

  for (j=p->nvars-1; j >= 0; j--)
    {
      InitializeVector(b=CreateVector(DimV(p->psi[j])),0.0);
      InitializeMatrix(S=CreateMatrix(DimV(p->psi[j]),DimV(p->psi[j])),0.0);

      for (k=p->n_coef_states[j]-1; k >= 0; k--)
    {
/*            // Accumulate XX and YY   ansi-c*/
      InitializeMatrix(XX,0.0);
      InitializeMatrix(XY,0.0);
      for (s=p->nstates-1; s >= 0; s--)
        if (p->coef_states[j][s] == k)
          {
        UpdateMS(XX,p->XX[s],p->Zeta[j][p->var_states[j][s]]);
        UpdateMS(XY,p->XY[s],p->Zeta[j][p->var_states[j][s]]);
          }

/*        // XY + XX*W[j]   ansi-c*/
      bSubtract(pElementM(XY),pElementM(XY),pElementM(XX),p->npre*p->nvars);

/*            // (XY + XX*W[j])*a0[j][k]   ansi-c*/
      ProductMV(v,XY,p->A0[j][k]);

/*        // b += PSI[k]'*diag(LAMBDA*lambda[j][k]+e)*(XY + XX*W[j])*a0[j][k]   ansi-c*/
      ElementV(b,p->npre-1+k)+=ElementV(v,i=p->npre-1);
      for (i--; i >= 0; )
        for (m=p->nvars-1; m >= 0; i--, m--)
          ElementV(b,i)+=ElementV(p->lambda[j][k],m)*ElementV(v,i);

/*        // S += PSI[k]'*diag(LAMBDA*lambda[j][k]+e)*XX*diag(LAMBDA*lambda[j][k]+e)*PSI[k]   ansi-c*/
      update_psi_quadratic_form(S,p->nvars,p->nlags,k,p->lambda[j][k],XX);

      if (_VERBOSE_COUNT)
        {
          fprintf(V_FILE,"//=== (j = %d   k = %d) ===//\n",j,k);
          fprintf(V_FILE,"XX =\n");
          dw_PrintMatrix(V_FILE,XX,"%lg ");

          fprintf(V_FILE,"lambda[%d][%d] =\n",j,k);
          dw_PrintVector(V_FILE,p->lambda[j][k],"%lg ");

          fprintf(V_FILE,"S =\n");
          dw_PrintMatrix(V_FILE,S,"%lg ");
        }
    }

/*        // Add inverse prior   ansi-c*/
      AddMM(S,S,p->inverse_psi_prior[j]);

/*       TMatrix U,V; */
/*       TVector d; */
/*       int size=RowM(p->inverse_psi_prior[j]); */
/*       printf("inverse psi prior (%d)\n",j); */
/*       SVD(U=CreateMatrix(size,size),d=CreateVector(size),V=CreateMatrix(size,size),p->inverse_psi_prior[j]); */
/*       dw_PrintVector(stdout,d,"%lg "); */
/*       FreeMatrix(U); FreeMatrix(V); FreeVector(d); */
/*       getchar(); */
      if (_VERBOSE_COUNT)
    {
      fprintf(V_FILE,"inverse prior=\n");
      dw_PrintMatrix(V_FILE,p->inverse_psi_prior[j],"%lg ");

      fprintf(V_FILE,"S =\n");
      dw_PrintMatrix(V_FILE,S,"%lg ");
      fprintf(V_FILE,"//=====================================================================//\n");
    }

/*        // Draw psi[j]   ansi-c*/
      if (!DrawNormal_InverseVariance(p->psi[j],b,S))
    SingularInverseVariance_RecoveryAttempt(p->psi[j],b,S,p->inverse_psi_prior[j],model,PSI_ERR);

      FreeMatrix(S);
      FreeVector(b);
    }

  FreeVector(v);
  FreeMatrix(XY);
  FreeMatrix(XX);
}

/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*  //-------------------------------- Draw lambda --------------------------------//   ansi-c*/
/*  //-----------------------------------------------------------------------------//   ansi-c*/
/*
   Assumes:
     S:  n x n matrix.
     XX: (n*b + 1) x (n*b + 1) symmetric matrix.  XX must be column major.
     psi:  (n*b+j) dimensional vector with j > 0.

   Results:
     Set S to

               LAMBDA'*diag(PSI[k]*psi)*XX*diag(PSI[k]*psi)*LAMBDA

     to S.  The matrices LAMBDA and PSI[k] are given by

                        -      -
                       |  I(n)  |
                       |   .    |                  -                  -
                       |   .    |                 |  I(n*b)   0(n*b,j) |
              LAMBDA = |   .    |        PSI[k] = |                    |
                       |        |                 | 0(1,n*b)   e(j,k)' |
                       |  I(n)  |                  -                  -
                       | 0(1,n) |
                        -      -

     where I(r) is the r x r identity matrix, 0(r,s) is the r x s zero matrix,
     and e(r,s) is the sth column of I(r).

   Notes:
     Even though PSI[k] depends on k, the term diag(PSI[k]*psi)*LAMBDA does not.
     For this reason, k is not passed.
*/
void lambda_quadratic_form(TMatrix S, int b, TVector psi, TMatrix XX)
{
  int bj, j, i, r, n=RowM(S);
  PRECISION *x=pElementM(XX)+(b*n-1)*RowM(XX), *y, *z=pElementV(psi), w;
  InitializeMatrix(S,0.0);
  for (bj=b-1; bj >= 0; bj--)
    for (y=pElementM(S)+(n-1)*n, j=n-1; j >= 0; x-=RowM(XX), y-=n, j--)
      {
    w=z[bj*n+j];
    for (r=b*n-1; r >= 0; )
      for (i=n-1; i >= 0; r--, i--)
        y[i]+=z[r]*x[r]*w;
      }
}

/*

   Notes:
     The matrices MUST be in column major format.  Basic matrix routines from
     bmatrix.c are called in this routine.
*/
void Draw_lambda(TStateModel *model)
{
  int j, k, s, i, m;
  TVector b, v;
  TMatrix S, XX, XY;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);

/*    // Update state dependent matrices if necessary   ansi-c*/
  if (!p->valid_state_dependent_fields)
    UpdateStateDependentFields(p,model->sv->S);

/*    // Allocate memory   ansi-c*/
  XX=CreateMatrix(p->npre,p->npre);
  XY=CreateMatrix(p->npre,p->nvars);
  v=CreateVector(p->npre);
  b=CreateVector(p->nvars);
  S=CreateMatrix(p->nvars,p->nvars);

  if (_VERBOSE_COUNT)
    fprintf(V_FILE,"//=== Draw_psi() (count = %d) =====================================//\n",_VERBOSE_COUNT);

  for (j=p->nvars-1; j >= 0; j--)
    {
      for (k=p->n_coef_states[j]-1; k > 0; k--)
    {
/*        // Accumulate XX and XY   ansi-c*/
      InitializeMatrix(XX,0.0);
      InitializeMatrix(XY,0.0);
      for (s=p->nstates-1; s >= 0; s--)
        if (p->coef_states[j][s] == k)
          {
        UpdateMS(XX,p->XX[s],p->Zeta[j][p->var_states[j][s]]);
        UpdateMS(XY,p->XY[s],p->Zeta[j][p->var_states[j][s]]);
          }

/*        // Compute mean   ansi-c*/
/*        // XY + XX*W[j]   ansi-c*/
      bSubtract(pElementM(XY),pElementM(XY),pElementM(XX),p->npre*p->nvars);

/*            // (XY + XX*W[j])*a0[j][k]   ansi-c*/
      ProductMV(v,XY,p->A0[j][k]);

/*        // (XY + XX*W[j])*a0[j][k] - XX*diag(PSI[k]*psi[j])*e   ansi-c*/
      bLinearUpdateScalar(pElementV(v),pElementM(XX)+p->npre*(p->npre-1),-ElementV(p->constant[j],k),p->npre);

/*        // b = LAMBDA'*diag(PSI[j][k]*psi[j])*((XY + XX*W[j])*a0[j][k] - XX*diag(PSI[k]*psi[j])*lambda_hat)   ansi-c*/
      InitializeVector(b,0.0);
      for (i=p->npre-2; i >= 0; )
        for (m=p->nvars-1; m >= 0; i--, m--)
          ElementV(b,m)+=ElementV(p->psi[j],i)*ElementV(v,i);

/*        // Compute inverse variance matrix   ansi-c*/
/*        // S = LAMBDA'*diag(PSI[j][k]*psi[j])*XX*diag(PSI[j][k]*psi[j])*LAMBDA   ansi-c*/
      lambda_quadratic_form(S,p->nlags,p->psi[j],XX);

      for (i=p->nvars-1; i >= 0; i--) ElementM(S,i,i)+=p->inverse_lambda_prior;

      if (_VERBOSE_COUNT)
        {
          fprintf(V_FILE,"//=== (j = %d   k = %d) ===//\n",j,k);
          fprintf(V_FILE,"XX =\n");
          dw_PrintMatrix(V_FILE,XX,"%lg ");

          fprintf(V_FILE,"psi[%d][%d] =\n",j,k);
          dw_PrintVector(V_FILE,p->psi[j],"%lg ");

          fprintf(V_FILE,"S =\n");
          dw_PrintMatrix(V_FILE,S,"%lg ");

          fprintf(V_FILE,"inverse prior = %lg\n",p->inverse_lambda_prior);
        }

/*        // Draw lambda[j][k]   ansi-c*/
      if (!DrawNormal_InverseVariance(p->lambda[j][k],b,S))
        {
          TMatrix InversePrior=IdentityMatrix((TMatrix)NULL,p->nvars);
          ProductMS(InversePrior,InversePrior,p->inverse_lambda_prior);
          SingularInverseVariance_RecoveryAttempt(p->lambda[j][k],b,S,InversePrior,model,LAMBDA_ERR);
          FreeMatrix(InversePrior);
        }
    }

/*        // State 0 normalized to one   ansi-c*/
      InitializeVector(p->lambda[j][k],1.0);
    }

  if (_VERBOSE_COUNT)
    fprintf(V_FILE,"//====================================================//\n");

/*    // Free memory   ansi-c*/
  FreeMatrix(S);
  FreeVector(b);
  FreeVector(v);
  FreeMatrix(XY);
  FreeMatrix(XX);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/********************** Updating T_VAR_Parameter Fields ***********************/
/******************************************************************************/
void UpdateStateDependentFields(T_VAR_Parameters *p, int *S)
{
  int i, s_prev, s_curr, t;

  dw_InitializeArray_int(p->T,0);

  if (!p->valid_state_dependent_fields_previous)
    {
/*        //=== Update Y'Y, X'Y, and X'X  ===   ansi-c*/
      for (i=dw_DimA(p->YY)-1; i >= 0; i--)
        {
          InitializeMatrix(p->YY[i],0.0);
          InitializeMatrix(p->XY[i],0.0);
          InitializeMatrix(p->XX[i],0.0);
        }
      for (t=p->nobs; t > 0; t--)
        {
          s_curr=S[t];
          AddMM(p->YY[s_curr],p->YY[s_curr],p->yy[t]);
          AddMM(p->XY[s_curr],p->XY[s_curr],p->xy[t]);
          AddMM(p->XX[s_curr],p->XX[s_curr],p->xx[t]);
          p->T[s_curr]++;
        }
/*        //p->valid_state_dependent_fields_previous=1;   ansi-c*/
    }
  else
    {
      for (t=p->nobs; t > 0; t--)
        {
          if ((s_curr=S[t]) != (s_prev=p->S[t]))
            {
              SubtractMM(p->YY[s_prev],p->YY[s_prev],p->yy[t]);
              SubtractMM(p->XY[s_prev],p->XY[s_prev],p->xy[t]);
              SubtractMM(p->XX[s_prev],p->XX[s_prev],p->xx[t]);
              AddMM(p->YY[s_curr],p->YY[s_curr],p->yy[t]);
              AddMM(p->XY[s_curr],p->XY[s_curr],p->xy[t]);
              AddMM(p->XX[s_curr],p->XX[s_curr],p->xx[t]);
            }
          p->T[s_curr]++;
        }
    }

/*    //memcpy(p->S,S,(p->nobs+1)*sizeof(int));   ansi-c*/

  p->valid_state_dependent_fields=1;
}

/*
   Assumes:
    p:  Pointer to valid T_VAR_Parameters structure with  b0 updated.

   Results:
    Updates A0

   Notes:
    Uses the relation

               A0[j][k] = U[j]*b0[j][k]

*/
void Update_A0_from_B0(T_VAR_Parameters *p)
{
  int j, k;

  for (j=p->nvars-1; j >= 0; j--)
    for (k=dw_DimA(p->A0[j])-1; k >= 0; k--)
      ProductMV(p->A0[j][k],p->U[j],p->b0[j][k]);
}

/*
    Sets

             Aplus[j][k] = V[j]*b0[j][k] - W[j]*A0[j][k]

    where 0 <= j < p->nvars and 0 <= k < p->n_coef_states[j].

    If (p->Specification & SPEC_RANDOMWALK) is set, uses the fact that
    W'[j] = [-I 0].  In this case a call to the base matrix function bAdd() is
    made.

    If p->IsIdentity_V[j] is set, uses the fact that V[j] = I.
*/
void Update_aplus_from_bplus_a0(int j, int k, T_VAR_Parameters *p)
{
  TVector v;

  if (p->IsIdentity_V[j])
    if (p->W[j])
      if (p->Specification & SPEC_RANDOM_WALK)
    {
      bAdd(pElementV(p->Aplus[j][k]),pElementV(p->bplus[j][k]),pElementV(p->A0[j][k]),p->nvars);
      memcpy(pElementV(p->Aplus[j][k]) + p->nvars,pElementV(p->bplus[j][k]) + p->nvars,(p->npre - p->nvars)*sizeof(PRECISION));
    }
      else
    {
      ProductMV(p->Aplus[j][k],p->W[j],p->A0[j][k]);
      SubtractVV(p->Aplus[j][k],p->bplus[j][k],p->Aplus[j][k]);
    }
    else
      EquateVector(p->Aplus[j][k],p->bplus[j][k]);
  else
    if (p->V[j])
      {
    ProductMV(p->Aplus[j][k],p->V[j],p->bplus[j][k]);
    if (p->W[j])
      if (p->Specification & SPEC_RANDOM_WALK)
        bAdd(pElementV(p->Aplus[j][k]),pElementV(p->Aplus[j][k]),pElementV(p->A0[j][k]),p->nvars);
      else
        {
          v=ProductMV((TVector)NULL,p->W[j],p->A0[j][k]);
          SubtractVV(p->Aplus[j][k],p->Aplus[j][k],v);
          FreeVector(v);
        }
      }
    else
      if (p->W[j])
    if (p->Specification & SPEC_RANDOM_WALK)
      {
        InitializeVector(p->Aplus[j][k],0.0);
        memcpy(pElementV(p->Aplus[j][k]),pElementV(p->A0[j][k]),p->nvars*sizeof(PRECISION));
      }
    else
      {
        ProductMV(p->Aplus[j][k],p->W[j],p->A0[j][k]);
        MinusV(p->Aplus[j][k],p->Aplus[j][k]);
      }
      else
    InitializeVector(p->Aplus[j][k],0.0);
}

/*
   Assumes:
    p:  Pointer to valid T_VAR_Parameters structure with A0 and bplus updated.

   Results:
    Updates Aplus

   Notes:
    Uses the relation

            Aplus[j][k] = V[j]*bplus[j][k] - W[j]*A0[j][k]

    If the Sims-Zha specification is used, special code is used to take advantage
    of the fact that V[j] is the identity and W[j] is diagonal with minus ones
    along the diagonal.

*/
void Update_Aplus_from_bplus_A0(T_VAR_Parameters *p)
{
  int j, k;
  for (j=p->nvars-1; j >= 0; j--)
    for (k=p->n_coef_states[j]-1; k >= 0; k--)
      Update_aplus_from_bplus_a0(j,k,p);
}

/*
   Assumes:
    p:  Pointer to valid TStateModel structure with psi and lambda updated.

   Results:
    Updates bplus and Aplus.

   Notes:
    Assumes Aplus[j][k] = V[j]*bplus[j][k] - W[j]*A0[j][k], V[j] is the identity,
    W[j] is a npre x nvar diagonal matrix with minus ones along the diagonal, and

                     | lambda[j][k][i % nvars]*psi[j][i] for 0 <= i < nvars*nlags
    bplus[j][k][i] = |
                     | psi[j][i]                         for nvar*nlags <= i

    This is the Sims-Zha specification.

*/
void Update_bplus_from_lambda_psi(T_VAR_Parameters *p)
{
  int j, k, i, m;
  PRECISION *p_bplus, *p_lambda, *p_psi;
  if (!(p->Specification & SPEC_SIMS_ZHA))
    {
      swz_fprintf_err("Update_bplus_from_lambda_psi() called without Sims-Zha specification\n");
      swzExit(0);
    }
  for (j=p->nvars-1; j >= 0; j--)
    {
      p_psi=pElementV(p->psi[j]);
      for (k=dw_DimA(p->bplus[j])-1; k >= 0; k--)
    {
      p_bplus=pElementV(p->bplus[j][k]);
      p_lambda=pElementV(p->lambda[j][k]);
          p_bplus[i=p->nlags * p->nvars]=ElementV(p->constant[j],k);
      for (i--; i >= 0; )
        for (m=p->nvars-1; m >= 0; i--, m--)
          p_bplus[i]=p_lambda[m]*p_psi[i];
    }
    }
}

/*
   Assumes:
    p:  Pointer to valid T_VAR_Parameters structure with  A0 and Aplus udated.

   Results:
    Updates b0 and bplus

   Notes:
    Uses the relations

               A0[j][k] = U[j]*b0[j][k]

            Aplus[j][k] = V[j]*bplus[j][k] - W[j]*A0[j][k]

             U'[j]*U[j] = identity

             V'[j]*V[j] = identity
*/
void Update_b0_bplus_from_A0_Aplus(T_VAR_Parameters *p)
{
  int i, j, k;
  TVector v;
  PRECISION *p_Aplus, *p_bplus, *p_A0;

/*    // A0   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    for (k=dw_DimA(p->A0[j])-1; k >= 0; k--)
      TransposeProductMV(p->b0[j][k],p->U[j],p->A0[j][k]);

/*    // Aplus   ansi-c*/
  if (p->Specification & SPEC_SIMS_ZHA)
    {
      for (j=p->nvars-1; j >= 0; j--)
    for (k=dw_DimA(p->Aplus[j])-1; k >= 0; k--)
      {
        p_A0=pElementV(p->A0[j][k]);
        p_Aplus=pElementV(p->Aplus[j][k]);
        p_bplus=pElementV(p->bplus[j][k]);
        i=p->nvars;
        memcpy(p_bplus+i,p_Aplus+i,(p->npre - p->nvars)*sizeof(PRECISION));
        for (i--; i >= 0; i--) p_bplus[i]=p_Aplus[i]-p_A0[i];
      }
    }
  else
    {
      v=CreateVector(p->npre);
      for (j=p->nvars-1; j >= 0; j--)
    if (p->V[j])
      for (k=dw_DimA(p->A0[j])-1; k >= 0; k--)
        if (p->W[j])
          {
        ProductMV(v,p->W[j],p->A0[j][k]);
        AddVV(v,v,p->Aplus[j][k]);
        TransposeProductMV(p->bplus[j][k],p->V[j],v);
          }
        else
          TransposeProductMV(p->bplus[j][k],p->V[j],p->Aplus[j][k]);
      FreeVector(v);
    }
}

/*
   Assumes:
    model:  pointer to valid TStateModel structure with A0 and bplus properly
            updated.

   Results:
    Updates lambda and psi.

   Notes:
    Assumes Aplus[j][k] = V[j]*bplus[j][k] - W[j]*A0[j][k], V[j] is the identity,
    W[j] is a npre x nvar matrix with minus ones along the diagonal and zeros
    elsewhere, and

                     | lambda[j][k][i % nvars]*psi[j][i] for 0 <= i < nvars*nlags
    bplus[j][k][i] = |
                     | psi[j][i]                        for nvar*nlags <= i

    The normalization lambda[j][0][i] = 1 is used.
*/
void Update_lambda_psi_from_bplus(T_VAR_Parameters *p)
{
  int j, k, i, m, n, dim=p->nlags*p->nvars;
  PRECISION *p_bplus, *p_lambda, *p_psi;
  for (j=p->nvars-1; j >= 0; j--)
    {
      p_psi=pElementV(p->psi[j]);
      p_bplus=pElementV(p->bplus[j][0]);
      memcpy(p_psi,p_bplus,dim*sizeof(PRECISION));
      ElementV(p->constant[j],0)=p_bplus[dim];
      InitializeVector(p->lambda[j][0],1.0);

      for (k=dw_DimA(p->bplus[j])-1; k > 0; k--)
    {
      p_bplus=pElementV(p->bplus[j][k]);
      p_lambda=pElementV(p->lambda[j][k]);
      ElementV(p->constant[j],k)=p_bplus[dim];
      for (m=p->nvars-1; m >= 0; m--)
        {
          for (n=dim+m, i=n-p->nvars; i >= 0; i-=p->nvars)
        if (fabs(p_psi[i]) > fabs(p_psi[n])) n=i;
          p_lambda[m]=(p_psi[n] != 0) ? p_bplus[n]/p_psi[n] : 0;
        }
    }
    }
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/********************************* Forecasts **********************************/
/******************************************************************************/
/*
   Assumes:
     forecast  : horizon x nvars matrix or null pointer
     horizon   : positive integer - forecast horizon
     initial   : initial value of predetermined variables.
     shocks    : array of length horizon of shocks or null pointer.  If null
                 pointer, then the shocks are all zero. Each vector is of length
                 nvar.
     S         : array of length horizon.  S[t] is the state at time T+1+t.
     model     : pointer to valid TStateModel structure.

   Results:
     Computes forecast

   Returns:
     The matrix forecast upon success and null upon failure.  If forecast is
     null, then it created.
*/
TMatrix forecast_base(TMatrix forecast, int horizon, TVector initial, TVector *shocks, int *S, TStateModel *model)
{
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);
  TMatrix *A0, *Aplus;
  TVector x, y;
  int i, t;

  // allocate forecast if necessary
  if (!forecast && !(forecast=CreateMatrix(horizon,p->nvars)))
    return (TMatrix)NULL;

  // allocate memory
  y=CreateVector(p->nvars);
  x=EquateVector((TVector)NULL,initial);
  A0=MakeA0_All((TMatrix*)NULL,p);
  Aplus=MakeAplus_All((TMatrix*)NULL,p);

  // forecast
  for (t=0; t < horizon; t++)
    {
      ProductVM(y,x,Aplus[S[t]]);
      if (shocks)
    {
      for (i=p->nvars-1; i >= 0; i--)
        ElementV(y,i)+=ElementV(shocks[t],i)/sqrt(p->Zeta[i][p->var_states[i][S[t]]]);
    }
      ProductInverseVM(y,y,A0[S[t]]);
      for (i=p->nvars-1; i >= 0; i--)
    ElementM(forecast,t,i)=ElementV(y,i);

      memmove(pElementV(x)+p->nvars,pElementV(x),(p->nlags-1)*p->nvars*sizeof(PRECISION));
      memcpy(pElementV(x),pElementV(y),p->nvars*sizeof(PRECISION));
    }

  // free memory
  dw_FreeArray(Aplus);
  dw_FreeArray(A0);
  FreeVector(x);
  FreeVector(y);

  return forecast;
}

/*
   For 1 <= k < h, y[k][i] is null if the ith coordinate of y(t0+1+k) is
   unrestricted and is its value otherwise.  In general, t0 is the last index for
   which we have full information.  It must be the case that t0 <= nobs.
*/
/* TVector* dw_state_space_mean_conditional_forecast(TVector *F, PRECISION ***y, int h, int t0, TStateModel *model) */
/* { */
/*   T_MSStateSpace *statespace=(T_MSStateSpace*)(model->theta); */
/*   TVector Pxi_i, *Pxi, *Pxi1, *SPxi, *Ez_i, **IEz, **Ez1, **Ez, **SEz, **ISEz, SPzeta, SPs, u, *z; */
/*   TMatrix Q, *Ezz_i, **IEzz, **Ezz1, **Ezz; */
/*   int i, k, s; */

/*   if ((t0 > statespace->t0) && !Filter(t0,model)) return (TVector*)NULL; */
/*   if ((t0 > model->t0) && !ForwardRecursion(t0,model)) return (TVector*)NULL; */
/*   if (((t0 < model->sv->t0) || (model->sv->t1 < t0)) && !sv_ComputeTransitionMatrix(t0,model->sv,model)) return (TVector*)NULL; */

/*   Pxi=dw_CreateArray_vector(h); */
/*   Pxi1=dw_CreateArray_vector(h); */
/*   SPxi=dw_CreateArray_vector(h); */
/*   IEz=dw_CreateRectangularArray_vector(h,statespace->zeta_modulus); */
/*   IEzz=dw_CreateRectangularArray_matrix(h,statespace->zeta_modulus); */
/*   Ez1=dw_CreateRectangularArray_vector(h,statespace->nstates); */
/*   Ezz1=dw_CreateRectangularArray_matrix(h,statespace->nstates); */
/*   Ez=dw_CreateRectangularArray_vector(h,statespace->nstates); */
/*   Ezz=dw_CreateRectangularArray_matrix(h,statespace->nstates); */
/*   SEz=dw_CreateRectangularArray_vector(h,statespace->nstates); */
/*   ISEz=dw_CreateRectangularArray_vector(h,statespace->zeta_modulus); */
/*   for (k=h-1; k >= 0; k--) */
/*     { */
/*       Pxi[k]=CreateVector(statespace->nstates); */
/*       Pxi1[k]=CreateVector(statespace->nstates); */
/*       SPxi[k]=CreateVector(statespace->nstates); */
/*       for (i=statespace->zeta_modulus-1; i >= 0; i--) */
/*     { */
/*       IEz[k][i]=CreateVector(statespace->nz); */
/*       IEzz[k][i]=CreateMatrix(statespace->nz,statespace->nz); */
/*       ISEz[k][i]=CreateVector(statespace->nz); */
/*     } */
/*       for (i=statespace->nstates-1; i >= 0; i--) */
/*     { */
/*       Ez1[k][i]=CreateVector(statespace->nz); */
/*       Ezz1[k][i]=CreateMatrix(statespace->nz,statespace->nz); */
/*       Ez[k][i]=CreateVector(statespace->nz); */
/*       Ezz[k][i]=CreateMatrix(statespace->nz,statespace->nz); */
/*       SEz[k][i]=CreateVector(statespace->nz); */
/*     } */
/*     } */

/*   ConditionalFilter(0,h-1,y,statespace->Ez[t0],statespace->Ezz[t0],model->V[t0],model->sv->Q, */
/*             Pxi,Pxi1,IEz,IEzz,Ez1,Ezz1,Ez,Ezz,statespace); */

/*   SmoothProbabilities_MSStateSpace(0,h-1,SPxi,Pxi,Pxi1,model->sv->Q); */

/*   SmoothMean_MSStateSpace(0,h-1,SEz,ISEz,Ez1,Ezz1,IEz,IEzz,SPxi,statespace); */

/*   SPzeta=CreateVector(statespace->zeta_modulus); */
/*   SPs=CreateVector(statespace->nbasestates); */
/*   u=CreateVector(statespace->ny); */
/*   z=dw_CreateArray_vector(statespace->nbasestates); */
/*   for (s=statespace->nbasestates-1; s >= 0; s--) z[s]=CreateVector(statespace->nz); */

/*   if (!F) */
/*     { */
/*       F=dw_CreateArray_vector(h); */
/*       for (k=h-1; k >= 0; k--) */
/*     F[k]=CreateVector(statespace->ny); */
/*     } */

/*   for (k=h-1; k >= 0; k--) */
/*     { */
/*       InitializeVector(F[k],0.0); */
/*       if (statespace->zeta_modulus > statespace->nbasestates) */
/*     { */
/*       IntegrateStatesSingle(SPzeta,SPxi[k],statespace->zeta_modulus,statespace->nbasestates,2); */
/*       IntegrateStatesSingleV(z,SPzeta,ISEz[k],statespace->nbasestates,statespace->zeta_modulus/statespace->nbasestates,2); */
/*       IntegrateStatesSingle(SPs,SPzeta,statespace->nbasestates,statespace->zeta_modulus/statespace->nbasestates,2); */
/*       for (s=statespace->nbasestates-1; s >= 0; s--) */
/*         { */
/*           ProductMV(u,statespace->H[s],z[s]); */
/*           AddVV(u,statespace->a[s],u); */
/*           LinearCombinationV(F[k],1.0,F[k],ElementV(SPs,s),u); */
/*         } */
/*     } */
/*       else */
/*     { */
/*       IntegrateStatesSingle(SPs,SPxi[k],statespace->nbasestates,statespace->zeta_modulus,2); */
/*       for (s=statespace->nbasestates-1; s >= 0; s--) */
/*         { */
/*           ProductMV(u,statespace->H[s],ISEz[k][s]); */
/*           AddVV(u,statespace->a[s],u); */
/*           LinearCombinationV(F[k],1.0,F[k],ElementV(SPs,s),u); */
/*         } */
/*     } */
/*     } */

/*   // Clean up */
/*   dw_FreeArray(z); */
/*   FreeVector(u); */
/*   FreeVector(SPs); */
/*   FreeVector(SPzeta); */
/*   dw_FreeArray(ISEz); */
/*   dw_FreeArray(SEz); */
/*   dw_FreeArray(Ezz); */
/*   dw_FreeArray(Ez); */
/*   dw_FreeArray(Ezz1); */
/*   dw_FreeArray(Ez1); */
/*   dw_FreeArray(IEzz); */
/*   dw_FreeArray(IEz); */
/*   dw_FreeArray(SPxi); */
/*   dw_FreeArray(Pxi1); */
/*   dw_FreeArray(Pxi); */

/*   return F; */
/* } */

/*

*/
/* TVector* dw_state_space_mean_unconditional_forecast(TVector *F, int h, int t0, TStateModel *model) */
/* { */
/*   T_MSStateSpace *statespace=(T_MSStateSpace*)(model->theta); */
/*   PRECISION ***y=(PRECISION***)dw_CreateMultidimensionalArrayList_scalar(3,h,statespace->ny,1); */
/*   int i, j; */
/*   for (i=h-1; i >= 0; i--) */
/*     for (j=statespace->ny-1; j >= 0; j--) */
/*       { dw_FreeArray(y[i][j]); y[i][j]=(PRECISION*)NULL; } */
/*   F=dw_state_space_mean_conditional_forecast(F,y,h,t0,model); */
/*   dw_FreeArray(y); */
/*   return F; */
/* } */


/******************************************************************************/
/** Impulse Response Routines                                                **/
/******************************************************************************/
/*
   Consider the model

     y'(t) = [y'(t-1) ... y'(t-p) z'(t)]*B + epsilon'(t) * Inverse(A0*Xi)


   The impulse response of variable j to shock i
   at horizon h is the element in row i and colum j of

             Inverse(A0*Xi) * J * S^(h-1) * J'

   where

                    B(1)   I ... 0
                     .     . .   .
               S =   .     .  .  .
                     .     .   . .
                   B(p-1)  0 ... I
                    B(p)   0 ... 0

               J = [I 0 ... 0]

               B'= [B'(1) ... B'(p) C']
*/
TMatrix ComputeImpulseResponseReducedForm(TMatrix R, int h, TMatrix A0_Xi_inv, TMatrix B, int nlags)
{
  TMatrix X=R, S, T, W;
  int n, t, i, j, m;

  dw_ClearError();

  if (!A0_Xi_inv || ((h > 1) && !B))
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  n=RowM(A0_Xi_inv);
  if (!X)
    {
      if (!(X=CreateMatrix(h,n*n))) return (TMatrix)NULL;
    }
  else
    if ((RowM(X) != h) || (ColM(X) != n*n))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }

  for (m=n*n-1, i=n-1; i >= 0; i--)
    for (j=n-1; j >= 0; m--, j--)
      ElementM(X,0,m)=ElementM(A0_Xi_inv,i,j);

  if (h > 1)
    {
      InitializeMatrix(S=CreateMatrix(n*nlags,n*nlags),0.0);
      for (i=n*(nlags - 1) - 1; i >= 0; i--) ElementM(S,i,i+n)=1.0;
      InsertSubMatrix(S,B,0,0,0,0,n*nlags,n);


      W=SubMatrix((TMatrix)NULL,S,0,0,n,n);
      ProductMM(W,A0_Xi_inv,W);
      for (m=n*n-1, i=n-1; i >= 0; i--)
    for (j=n-1; j >= 0; m--, j--)
      ElementM(X,1,m)=ElementM(W,i,j);

      if (h > 2)
    {
      T=EquateMatrix((TMatrix)NULL,S);

      for (t=2; t < h; t++)
        {
          ProductMM(T,T,S);
          SubMatrix(W,T,0,0,n,n);
          ProductMM(W,A0_Xi_inv,W);
          for (m=n*n-1, i=n-1; i >= 0; i--)
        for (j=n-1; j >= 0; m--, j--)
          ElementM(X,t,m)=ElementM(W,i,j);
        }

      FreeMatrix(T);
    }

      FreeMatrix(W);
      FreeMatrix(S);
    }

  if (dw_GetError() != NO_ERR)
    {
      if (X != R) FreeMatrix(X);
      return (TMatrix)NULL;
    }

  return X;
}

/*
   Consider the model

     y'(t)*A(0) = [y'(t-1) ... y'(t-p) z'(t)]*Aplus + epsilon'(t)*Inverse(Xi)


   The impulse response of variable j to shock i
   at horizon h is the element in row i and colum j of

             Inverse(A0*Xi) * J * S^(h-1) * J'

   where

                    B(1)   I ... 0
                     .     . .   .
               S =   .     .  .  .
                     .     .   . .
                   B(p-1)  0 ... I
                    B(p)   0 ... 0

               J = [I 0 ... 0]

               B = Aplus * Inverse(A0)

               B'= [B'(1) ... B'(p) C']

   Note that if Y'(t)=[y'(t) ... y'(t-p+1)], then

      Y'(t) = Y'(t-1)*S + J*z'(t)*C + J*epsilon'(t)*Inverse(A0*Xi)
*/
TMatrix ComputeImpulseResponseStructural(TMatrix R, int h, TMatrix A0, TMatrix Aplus, TVector Xi, int nlags)
{
  TMatrix X, B;
  int n=RowM(A0), i, j;
  PRECISION xi_inv;

  X=Inverse_LU((TMatrix)NULL,A0);

  B=(h > 1) ? ProductMM((TMatrix)NULL,Aplus,X) : (TMatrix)NULL;

  for (i=n-1; i >= 0; i--)
    for (xi_inv=1.0/ElementV(Xi,i), j=n-1; j >= 0; j--)
      ElementM(X,i,j)*=xi_inv;

  R=ComputeImpulseResponseReducedForm(R,h,X,B,nlags);

  FreeMatrix(X);
  FreeMatrix(B);

  return R;
}

/*
   Consider the model

     y(t)' * A(0)(s(t)) = Sum[y(t-i)' * A(i)(s(t)),i=1,...,p] +
                             + z(t)' * C + epsilon(t)' * Inverse(Xi(s(t)))


   Condition on state k occuring, the impulse response of variable j to shock i
   at horizon h is the element in row i and colum j of

             Inverse(Xi(k)) * Inverse(A(0)(k)) * J * S^(h-1) * J'

   where

               S =  A(1)(k)*Inverse(A(0)(k))   I ... 0
                             .                 . .   .
                             .                 .  .  .
                             .                 .   . .
                   A(p-1)(k)*Inverse(A(0)(k))  0 ... I
                    A(p)(k)*Inverse(A(0)(k))   0 ... 0
   and

                                   J = [I 0 ... 0].
*/
TMatrix ComputeImpulseResponse(TMatrix R, int h, int k, TStateModel *model)
{
  TMatrix X, Aplus, B;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);
  PRECISION xi_inv;
  int n=p->nvars, i, j;

  if ((k < 0) || (k >= p->nstates))
    {
      dw_Error(SIZE_ERR);
      return (TMatrix)NULL;
    }

  X=MakeA0((TMatrix)NULL,k,p);
  Inverse_LU(X,X);

  if (h > 1)
    {
      Aplus=MakeAplus((TMatrix)NULL,k,p);
      B=ProductMM((TMatrix)NULL,Aplus,X);
    }
  else
    Aplus=B=(TMatrix)NULL;

  for (i=n-1; i >= 0; i--)
    for (xi_inv=1.0/sqrt(p->Zeta[i][p->var_states[i][k]]), j=n-1; j >= 0; j--)
      ElementM(X,i,j)*=xi_inv;

  R=ComputeImpulseResponseReducedForm(R,h,X,B,p->nlags);

  FreeMatrix(B);
  FreeMatrix(Aplus);
  FreeMatrix(X);

  return R;
}

/*
   Computes the cummulative variance decomposition of the impulse responses (IR).
*/
TMatrix ComputeVarianceDecomposition(TMatrix X, TMatrix IR, int nvars)
{
  TMatrix Y=X;
  int t, i, j;
  PRECISION sum, tmp;

  dw_ClearError();

  if (!IR)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  else
    if (ColM(IR) != nvars*nvars)
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }
    else
      if (!Y)
    {
      if (!(Y=CreateMatrix(RowM(IR),ColM(IR)))) return (TMatrix)NULL;
    }
      else
    if ((RowM(Y) != RowM(IR)) || (ColM(Y) != ColM(IR)))
      {
        dw_Error(SIZE_ERR);
        return (TMatrix)NULL;
      }

/*    // Compute cummulative variation   ansi-c*/
  for (j=nvars*nvars-1; j >= 0; j--)
    {
      tmp=ElementM(IR,0,j);
      ElementM(Y,0,j)=tmp*tmp;
    }
  for (t=1; t < RowM(IR); t++)
    for (j=nvars*nvars-1; j >= 0; j--)
      {
    tmp=ElementM(IR,t,j);
    ElementM(Y,t,j)=tmp*tmp + ElementM(Y,t-1,j);
      }

/*    // Compute cummulative variance decomposition   ansi-c*/
  for (t=0; t < RowM(IR); t++)
    for (j=nvars-1; j >= 0; j--)
      {
    for (sum=0.0, i=nvars*(nvars-1)+j; i >= 0; i-=nvars) sum+=ElementM(Y,t,i);
    if (sum > 0)
      for (sum=1.0/sum, i=nvars*(nvars-1)+j; i >= 0; i-=nvars) ElementM(Y,t,i)*=sum;
      }

  return Y;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/



/*******************************************************************************/
/******************************** Optimization *********************************/
/*******************************************************************************/
/* PRECISION ComputeConstantSimsZha(TStateModel *model) */
/* { */
/*   int j; */
/*   PRECISION constant=0.0; */
/*   T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p); */
/*   for (j=p->nvars-1; j >= 0; j--) */
/*     { */
/*       // b0 */
/*       constant+=0.5*p->n_coef_states[j]*DimV(p->b0[j][0])*log(p->n_A0_states/p->n_coef_states[j]); */

/*       // bplus or psi and constant */
/*       if (p->bplus[j]) */
/*     if (p->Specification & SPEC_SIMS_ZHA) */
/*       constant+=0.5*((p->npre-1)*log(p->n_A0_states) + p->n_coef_states[j]*log(p->n_A0_states/p->n_coef_states[j])); */
/*     else */
/*       constant+=0.5*p->n_coef_states[j]*DimV(p->bplus[j][0])*log(p->n_A0_states/p->n_coef_states[j]); */
/*     } */
/*   return constant; */
/* } */

#define LN_TWO_PI 1.837877066409345
#define LN_TWO    0.6931471805599453
void SetLogPriorConstant_VAR(T_VAR_Parameters *p)
{
  int j;

  p->log_prior_constant=0.0;

  for (j=p->nvars-1; j >= 0; j--)
    {
/*        // Gamma prior on Zeta   ansi-c*/
      if (p->n_var_states[j] > 1)
     p->log_prior_constant+=(p->n_var_states[j] - 1)*(ElementV(p->zeta_a_prior,j)*log(ElementV(p->zeta_b_prior,j))
                            - dw_log_gamma(ElementV(p->zeta_a_prior,j)));

/*        // Normal prior on b0 (A0)   ansi-c*/
      p->log_prior_constant+=0.5*p->n_coef_states[j]*(-DimV(p->b0[j][0])*LN_TWO_PI + LogAbsDeterminant_LU(p->inverse_b0_prior[j]));

      if (p->Specification & SPEC_SIMS_ZHA)
    {
/*        // Independent normal prior each element of delta with variance equal to delta_prior   ansi-c*/
       p->log_prior_constant-=0.5 * (dw_DimA(p->lambda[j])-1) * p->nvars * (LN_TWO_PI + log(p->lambda_prior));

/*        // Normal prior on psi and constant   ansi-c*/
       p->log_prior_constant+=0.5*(-(p->npre-1+p->n_coef_states[j])*LN_TWO_PI + LogAbsDeterminant_LU(p->inverse_psi_prior[j]));
    }
      else
    {
/*        // Normal prior on bplus (Aplus)   ansi-c*/
      if (p->bplus[j])
        p->log_prior_constant+=0.5*p->n_coef_states[j]*(-DimV(p->bplus[j][0])*LN_TWO_PI + LogAbsDeterminant_LU(p->inverse_bplus_prior[j]));
    }
    }

/*    // Scale for normalization   ansi-c*/
  switch (p->normalization_type)
    {
    case VAR_NORMALIZATION_NONE:
      break;
    case VAR_NORMALIZATION_WZ:
      p->log_prior_constant+=p->nvars*log(2);
      break;
    default:
      printf("Unknown normalization type\n");
      swzExit(1);
    }
}
#undef LN_TWO_PI
#undef LN_TWO

PRECISION LogPrior_VAR(TStateModel *model)
{
  int j, k;
  PRECISION x, y, log_prior=0.0;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  if (!(p->valid_parameters)) return MINUS_INFINITY;

  for (j=p->nvars-1; j >= 0; j--)
    {
/*        // Gamma prior on Zeta   ansi-c*/
      if (p->n_var_states[j] > 1)
    {
      x=ElementV(p->zeta_b_prior,j);
      y=2*(ElementV(p->zeta_a_prior,j)-1);
      for (k=p->n_var_states[j]-1; k > 0; k--)
        log_prior+=y*sqrt(p->Zeta[j][k]) - x*p->Zeta[j][k];
    }

/*        // Normal prior on b0 (A0)   ansi-c*/
      for (k=p->n_coef_states[j]-1; k >= 0; k--)
     log_prior-=0.5*InnerProductSymmetric(p->b0[j][k],p->inverse_b0_prior[j]);

      if (p->Specification & SPEC_SIMS_ZHA)
    {
/*        // Independent normal prior each element of lambda with variance equal to lambda_prior   ansi-c*/
      if (dw_DimA(p->lambda[j]) > 1)
        {
          for (x=0.0, k=p->n_coef_states[j]-1; k > 0; k--)
        x+=DotProduct(p->lambda[j][k],p->lambda[j][k]);
          log_prior-=0.5 * x * p->inverse_lambda_prior;
        }

/*        // Normal prior on psi and constant   ansi-c*/
      log_prior-=0.5*InnerProductSymmetric(p->psi[j],p->inverse_psi_prior[j]);
    }
      else
    {
/*        // Normal prior on bplus (Aplus)   ansi-c*/
      for (k=p->n_coef_states[j]-1; k >= 0; k--)
        log_prior-=0.5*InnerProductSymmetric(p->bplus[j][k],p->inverse_bplus_prior[j]);
    }
    }

  return log_prior + p->log_prior_constant;
}

int NumberFreeParametersVAR(TStateModel *model)
{
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);
  int j,  k, size=0;

/*    // b0   ansi-c*/
  for (j=0; j < p->nvars; j++)
    for (k=0; k < dw_DimA(p->b0[j]); k++)
      size+=DimV(p->b0[j][k]);

  if (p->Specification & SPEC_SIMS_ZHA)
    {
/*        // lambda   ansi-c*/
      for (j=0; j < p->nvars; j++)
    for (k=1; k < dw_DimA(p->lambda[j]); k++)
      size+=DimV(p->lambda[j][k]);

/*        // psi   ansi-c*/
      for (j=0; j < p->nvars; j++)
    size+=DimV(p->psi[j]);
    }
  else
    {
/*        // bplus   ansi-c*/
      for (j=0; j < p->nvars; j++)
    if (p->bplus[j])
      for (k=0; k < dw_DimA(p->bplus[j]); k++)
        size+=DimV(p->bplus[j][k]);
    }

/*    // Zeta   ansi-c*/
  for (j=0; j < p->nvars; j++)
    size+=dw_DimA(p->Zeta[j])-1;

  return size;
}

void FreeParametersToVAR(TStateModel *model, PRECISION *f)
{
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);
  int k, j;

/*    // b0   ansi-c*/
  for (j=0; j < p->nvars; j++)
    for (k=0; k < dw_DimA(p->b0[j]); k++)
      {
    memcpy(pElementV(p->b0[j][k]),f,DimV(p->b0[j][k])*sizeof(PRECISION));
    f+=DimV(p->b0[j][k]);
      }

  if (p->Specification & SPEC_SIMS_ZHA)
    {
/*        // lambda   ansi-c*/
      for (j=0; j < p->nvars; j++)
    {
      InitializeVector(p->lambda[j][0],1.0);
      for (k=1; k < dw_DimA(p->lambda[j]); k++)
        {
          memcpy(pElementV(p->lambda[j][k]),f,DimV(p->lambda[j][k])*sizeof(PRECISION));
          f+=DimV(p->lambda[j][k]);
        }
    }

/*        // psi   ansi-c*/
      for (j=0; j < p->nvars; j++)
    {
      memcpy(pElementV(p->psi[j]),f,DimV(p->psi[j])*sizeof(PRECISION));
      f+=DimV(p->psi[j]);
    }
    }
  else
    {
/*        // bplus   ansi-c*/
      for (j=0; j < p->nvars; j++)
    if (p->bplus[j])
      for (k=0; k < dw_DimA(p->bplus[j]); k++)
        {
          memcpy(pElementV(p->bplus[j][k]),f,DimV(p->bplus[j][k])*sizeof(PRECISION));
          f+=DimV(p->bplus[j][k]);
        }
    }

/*    // Zeta   ansi-c*/
  for (j=0; j < p->nvars; j++)
    {
/*        // Zeta non-negative   ansi-c*/
      for (k=dw_DimA(p->Zeta[j])-2; k >= 0; k--)
    if (f[k] < 0.0)
      {
        p->valid_parameters=0;
        return;
      }

      p->Zeta[j][0]=1.0;
      memcpy(p->Zeta[j]+1,f,(dw_DimA(p->Zeta[j])-1)*sizeof(PRECISION));
      f+=dw_DimA(p->Zeta[j])-1;
    }

/*    // Valid normalization   ansi-c*/


/*    // Update A0 and Aplus   ansi-c*/
  if (p->Specification & SPEC_SIMS_ZHA) Update_bplus_from_lambda_psi(p);
  Update_A0_from_B0(p);
  Update_Aplus_from_bplus_A0(p);

/*    // Set flags   ansi-c*/
  p->valid_parameters=1;
  ThetaChanged(model);
}

void VARToFreeParameters(TStateModel *model, PRECISION *f)
{
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->p->p);
  int k, j;

/*    // b0   ansi-c*/
  for (j=0; j < p->nvars; j++)
    for (k=0; k < dw_DimA(p->b0[j]); k++)
      {
    memcpy(f,pElementV(p->b0[j][k]),DimV(p->b0[j][k])*sizeof(PRECISION));
    f+=DimV(p->b0[j][k]);
      }

  if (p->Specification & SPEC_SIMS_ZHA)
    {
/*        // lambda   ansi-c*/
      for (j=0; j < p->nvars; j++)
    for (k=1; k < dw_DimA(p->lambda[j]); k++)
      {
        memcpy(f,pElementV(p->lambda[j][k]),DimV(p->lambda[j][k])*sizeof(PRECISION));
        f+=DimV(p->lambda[j][k]);
      }

/*        // psi   ansi-c*/
      for (j=0; j < p->nvars; j++)
    {
      memcpy(f,pElementV(p->psi[j]),DimV(p->psi[j])*sizeof(PRECISION));
      f+=DimV(p->psi[j]);
    }
    }
  else
    {
/*        //bplus   ansi-c*/
      for (j=0; j < p->nvars; j++)
    if (p->bplus[j])
      for (k=0; k < dw_DimA(p->bplus[j]); k++)
        {
          memcpy(f,pElementV(p->bplus[j][k]),DimV(p->bplus[j][k])*sizeof(PRECISION));
          f+=DimV(p->bplus[j][k]);
        }
    }

/*    // Zeta   ansi-c*/
  for (j=0; j < p->nvars; j++)
    {
      memcpy(f,p->Zeta[j]+1,(dw_DimA(p->Zeta[j])-1)*sizeof(PRECISION));
      f+=dw_DimA(p->Zeta[j])-1;
    }
}

/*
   Assumes:
     p:  pointer to valid T_VAR_Parameters structure

   Returns:
     The starting position of the Zeta parameters in the array of free
     parameters.
*/
int ZetaIndex(T_VAR_Parameters *p)
{
  int j, k, index=0;

/*    // b0   ansi-c*/
  for (j=0; j < p->nvars; j++)
    for (k=0; k < dw_DimA(p->b0[j]); k++)
      index+=DimV(p->b0[j][k]);

  if (p->Specification & SPEC_SIMS_ZHA)
    {
/*        // lambda   ansi-c*/
      for (j=0; j < p->nvars; j++)
    for (k=1; k < dw_DimA(p->lambda[j]); k++)
      index+=DimV(p->lambda[j][k]);

/*        // psi   ansi-c*/
      for (j=0; j < p->nvars; j++)
    index+=DimV(p->psi[j]);
    }
  else
    {
/*        // bplus   ansi-c*/
      for (j=0; j < p->nvars; j++)
    if (p->bplus[j])
      for (k=0; k < dw_DimA(p->bplus[j]); k++)
        index+=DimV(p->bplus[j][k]);
    }

  return index;
}

/*
   Assumes:
     p:  pointer to valid T_VAR_Parameters structure

   Returns:
     The the number of Zeta parameters in the array of free parameters.
*/
int ZetaLength(T_VAR_Parameters *p)
{
  int j, length=0;
  for (j=0; j < p->nvars; j++) length+=dw_DimA(p->Zeta[j])-1;
  return length;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/******************************** Normalization ********************************/
/*******************************************************************************/
/*
   Assumes
     p : pointer to properly initialized T_VAR_Parameters structure.

   Returns
     1 : parameters normalized according to p->normalization_type and
         p->normalized.
     0 : parameters not normalized according to p->normalization_type and
         p->normalized.
*/
int IsNormalized_VAR(T_VAR_Parameters *p)
{
  return (p->normalized == p->normalization_type) ? 1 : 0;
}
/*
   Assumes
     p : pointer to properly initialized T_VAR_Parameters structure.

   Returns
     1 : successful normalization/at least one column changed
     0 : successful normalization/no column changed
    -1 : unsuccessful normalization (should not be returned)

*/
int Normalize_VAR(T_VAR_Parameters *p)
{
  switch(p->normalization_type)
    {
    case VAR_NORMALIZATION_WZ: return WZ_Normalize(p);
    case VAR_NORMALIZATION_NONE: return 0;
    default:
      printf("Unknown normalization type\n");
      swzExit(1);
    }
}


/*
   Changes the sign of the jth column for state k.  It must be the case that
   0 <= j < nvars and 0 <= k < n_coef_states[j]
*/
void ChangeSign(int j, int k, T_VAR_Parameters *p)
{
/*    // Change sign of A0[j][k]   ansi-c*/
  MinusV(p->A0[j][k],p->A0[j][k]);

/*    // Change sign of Aplus[j][k]   ansi-c*/
  MinusV(p->Aplus[j][k],p->Aplus[j][k]);

/*    // Change sign of b0[j][k]   ansi-c*/
  MinusV(p->b0[j][k],p->b0[j][k]);

/*    // Change sign of bplus[j][k]   ansi-c*/
  if (p->bplus[j]) MinusV(p->bplus[j][k],p->bplus[j][k]);

  if (p->Specification & SPEC_SIMS_ZHA)
    {
/*        // Change sign of constant   ansi-c*/
      ElementV(p->constant[j],k)=-ElementV(p->constant[j],k);

/*        // Change sign of lambda[j][p->A0_column_states[j][k]]   ansi-c*/
      MinusV(p->lambda[j][k],p->lambda[j][k]);
    }
}

void Setup_No_Normalization(T_VAR_Parameters *p)
{
  p->normalization_type=p->normalized=VAR_NORMALIZATION_NONE;

  if (p->flipped) dw_FreeArray(p->flipped);
  if (p->Target) dw_FreeArray(p->Target);

  p->flipped=(int**)NULL;
  p->Target=(TVector**)NULL;

  SetLogPriorConstant_VAR(p);
}

/*
   Sets up the WZ normalization.  The target parameters come from A0.
*/
void Setup_WZ_Normalization(T_VAR_Parameters *p, TVector **A0)
{
  int j;

  p->normalization_type=VAR_NORMALIZATION_WZ;

  if (p->flipped) dw_FreeArray(p->flipped);
  if (p->Target) dw_FreeArray(p->Target);

  p->flipped=dw_CreateArray_array(p->nvars);
  for (j=p->nvars-1; j >= 0; j--)
    p->flipped[j]=dw_CreateArray_int(p->n_coef_states[j]);

  p->Target=dw_CopyArray((TVector**)NULL,A0);

  p->WZ_inconsistancies=0;

  SetLogPriorConstant_VAR(p);

  WZ_Normalize(p);
}

/*
   Assumes
     p : pointer to properly initialized T_VAR_Parameters structure.

   Returns
     1 : successful normalization/at least one column changed
     0 : successful normalization/no column changed

*/
int WZ_Normalize(T_VAR_Parameters *p)
{
  int j, k, changed=0, inconsistent=0;
  TMatrix A0, M, Target;

/*    // zero flipped   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    for (k=p->n_coef_states[j]-1; k >= 0; k--)
      p->flipped[j][k]=0;

/*    // determine which columns to flip   ansi-c*/
  A0=CreateMatrix(p->nvars,p->nvars);
  Target=CreateMatrix(p->nvars,p->nvars);
  M=CreateMatrix(p->nvars,p->nvars);
  for (k=p->n_A0_states-1; k >= 0; k--)
    {
      for (j=p->nvars-1; j >= 0; j--)
    {
      memcpy(&ElementM(A0,0,j),pElementV(p->A0[j][p->A0_column_states[j][k]]),p->nvars*sizeof(PRECISION));
      memcpy(&ElementM(Target,0,j),pElementV(p->Target[j][p->A0_column_states[j][k]]),p->nvars*sizeof(PRECISION));
    }
      InverseProductMM(M,A0,Target);
      for (j=p->nvars-1; j >= 0; j--)
    if (ElementM(M,j,j) < 0.0)
      p->flipped[j][p->A0_column_states[j][k]]--;
    else
      p->flipped[j][p->A0_column_states[j][k]]++;
    }
  FreeMatrix(M);
  FreeMatrix(Target);
  FreeMatrix(A0);

/*    // flip columns and record inconsistencies   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    for (k=p->n_coef_states[j]-1; k >= 0; k--)
      if (p->flipped[j][k]*p->n_coef_states[j] == -p->n_A0_states)
    {
      changed=1;
      ChangeSign(j,k,p);
    }
      else
    if (p->flipped[j][k]*p->n_coef_states[j] != p->n_A0_states)
      {
        inconsistent=1;
        if ((p->flipped[j][k] < 0) || ((p->flipped[j][k] == 0) && (ElementV(p->A0[j][p->A0_column_states[j][k]],j) < 0.0)))
          {
        changed=1;
        ChangeSign(j,k,p);
          }
      }

  if (inconsistent) p->WZ_inconsistancies++;

  p->normalized=VAR_NORMALIZATION_WZ;

  return changed;
}

/* /\* */
/*    Assumes */
/*      A0:  p->nvars x p->nvars matrix */
/*      k:  0 <= k < p->n_A0_states */
/*      p:  valid pointer to T_VAR_Parameters stucture */

/*    Results */
/*      A0 is initialized to state k. */

/*        A0 = [A0[A0_column_state[0][k]], ..., A0[A0_column_state[nvars-1][k]]] */

/* *\/ */
/* void CreateA0_from_determinant_state(TMatrix A0, int k, T_VAR_Parameters *p) */
/* { */
/*   int j; */
/*   for (j=p->nvars-1; j >= 0; j--) */
/*     CopyColumnVector(A0,p->A0[j][p->A0_column_states[j][k]],j); */
/* } */

/* /\* */
/*    Assumes */
/*     p:  Valid pointer to T_VAR_Parameters structure. */
/*     Ref:  p->nvars x p->nvars matrix */

/*    Results */
/*     A0, Aplus, b0, bplus are normalized so that the diagonal of Inverse(A0)*Ref  */
/*     is positive. */

/*    Notes */
/*     The normalization described above is equivalent to requiring that the jth  */
/*     column of A0 and the jth column of Ref lie on the same side of the hyperplane  */
/*     spanned all the columns of A0 other than the jth. */
/* *\/ */
/* void Normalize_WZ(T_VAR_Parameters *p, TMatrix Ref) */
/* { */
/*   int j, k; */
/*   TMatrix A=CreateMatrix(p->nvars,p->nvars); */

/*   for (k=p->n_A0_states-1; k >= 0; k--) */
/*     { */
/*       CreateA0_from_determinant_state(A,k,p); */

/*       if (!InverseProductMM(A,A,Ref)) */
/*     { */
/*       printf("\nNormalize_WZ():  A0 not invertible\n"); */
/*       swzExit(0); */
/*     } */

/*       for (j=p->nvars-1; j >= 0; j--) */
/*     if (ElementM(A,j,j) < 0) */
/*       ChangeSign(j,k,p); */
/*     } */

/*   FreeMatrix(A); */
/* } */

/* /\* */
/*    Assumes */
/*     p:  Valid pointer to T_VAR_Parameters structure. */
/*     Ref:  Array of length p->nvars of vectors */

/*    Results */
/*     A0, Aplus, b0, and bplus are normalized so that the diagonal of A0'*Ref is  */
/*     positive. */

/*    Notes */
/*     If Ref has exactly one non-zero element in each column, this normalization */
/*     is equivalent to requiring that the corresponding elements of A0 are */
/*     positive. */
/* *\/ */
/* void Normalize_Traditional(T_VAR_Parameters *p, TVector *Ref) */
/* { */
/*   int j, k; */

/*   for (j=p->nvars-1; j >= 0; j--) */
/*     for (k=p->n_coef_states[j]-1; k >= 0; k--) */
/*       if (DotProduct(p->A0[j][k],Ref[j]) < 0) */
/*     { */
/*       // Change sign of A0[j][k] */
/*       MinusV(p->A0[j][k],p->A0[j][k]); */

/*       // Change sign of Aplus[j][k] */
/*       MinusV(p->Aplus[j][k],p->Aplus[j][k]); */

/*       // Change sign of b0[j][k] */
/*       MinusV(p->b0[j][k],p->b0[j][k]); */

/*       // Change sign of bplus[j][k] */
/*       if (p->bplus[j]) MinusV(p->bplus[j][k],p->bplus[j][k]); */
/*     } */
/* } */

/* /\* */
/*    Assumes */
/*     p:  Valid pointer to T_VAR_Parameters structure */

/*    Results */
/*     A0, Aplus, b0, and bplus are normalized so that the diagonal of A0 is  */
/*     positive. */
/* *\/ */
/* void Normalize_Diagonal(T_VAR_Parameters *p) */
/* { */
/*   int j, k; */

/*   for (j=p->nvars-1; j >= 0; j--) */
/*     for (k=p->n_coef_states[j]-1; k >= 0; k--) */
/*       if (ElementV(p->A0[j][k],j) < 0) */
/*     { */
/*       // Change sign of A0[j][k] */
/*       MinusV(p->A0[j][k],p->A0[j][k]); */

/*       // Change sign of Aplus[j][k] */
/*       MinusV(p->Aplus[j][k],p->Aplus[j][k]); */

/*       // Change sign of b0[j][k] */
/*       MinusV(p->b0[j][k],p->b0[j][k]); */

/*       // Change sign of bplus[j][k] */
/*       if (p->bplus[j]) MinusV(p->bplus[j][k],p->bplus[j][k]); */
/*     } */
/* } */
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/********************************** Utilities **********************************/
/*******************************************************************************/
/*
   Assumes:
     X : A m x n matrix or null pointer.
     Y : A array of pointers to vectors of length n.  For each i, Y[i] is an
         array of vectors of positive length.  For each i and j, Y[i][j] is
         a vector of length m.  Y must have be created via calls to the function
         CreateVectorMultidimensionArray() so that the macros DimA(Y) and
         DimA(Y[i]) are valid.

   Results:
     Creates X if X is null.  Sets X[i][j] to Y[j][k][i] if k is less than
     DimA(Y[j]) and to Y[j][0][i] otherwise.

   Returns:
     Returns the matrix X.

   Notes:
     The routine does not check to ensure that every Y[i] is non-null (and hence
     of positive length) nor does it check that all the vectors Y[i][j] are of
     the same length.
*/
TMatrix ConstructMatrixFromColumns(TMatrix X, TVector **Y, int k)
{
  int i, j, s;
  if (!Y)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  if (!X)
    {
      if (!(X=CreateMatrix(DimV(Y[0][0]),dw_DimA(Y))))
    return (TMatrix)NULL;
    }
  else
    if ((RowM(X) != DimV(Y[0][0])) || (ColM(X) != dw_DimA(Y)))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }
  if (MajorForm(X))
    for (i=RowM(X)*sizeof(PRECISION), j=ColM(X)-1; j >= 0; j--)
      memcpy(&ElementM(X,0,j),pElementV(Y[j][(k < dw_DimA(Y[j])) ? k : 0]),i);
  else
    for (j=ColM(X)-1; j >= 0; j--)
      {
    s=(k < dw_DimA(Y[j])) ? k : 0;
    for (i=RowM(X)-1; i >= 0; i--) ElementM(X,i,j)=ElementV(Y[j][s],i);
      }
  return X;
}

/*
   Assumes:
     A0 : p->nvars x p->nvars matrix or null pointer
     k  : 0 <= k < p->nstates
*/
TMatrix MakeA0(TMatrix A0, int s, T_VAR_Parameters *p)
{
  int j;
  if (!A0)
    {
      if (!(A0=CreateMatrix(p->nvars,p->nvars)))
    return (TMatrix)NULL;
    }
  else
    if ((RowM(A0) != p->nvars) || (ColM(A0) != p->nvars))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }
  for (j=0; j < p->nvars; j++)
    memcpy(&ElementM(A0,0,j),pElementV(p->A0[j][p->coef_states[j][s]]),p->nvars*sizeof(PRECISION));
  return A0;
}

/*
   Assumes:
     A0 : Matrix array of length n_states or null pointer.  A0[s] is either
          p->nvars x p->nvars matrix or null pointer
*/
TMatrix* MakeA0_All(TMatrix *A0, T_VAR_Parameters *p)
{
  int s;
  TMatrix *A0_in=A0;
  if (!A0)
    {
      if (!(A0=dw_CreateArray_matrix(p->nstates)))
    return (TMatrix*)NULL;
    }
  else
    if ((dw_DimA(A0) != p->nstates))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix*)NULL;
      }
  for (s=p->nstates-1; s >= 0; s--)
    if (!(A0[s]=MakeA0(A0[s],s,p)))
      {
    if (A0_in != A0) dw_FreeArray(A0);
    return (TMatrix*)NULL;
      }
  return A0;
}


/*
   Assumes:
     Aplus : p->npre x p->nvars matrix or null pointer
     k     : 0 <= k < p->nstates
*/
TMatrix MakeAplus(TMatrix Aplus, int k, T_VAR_Parameters *p)
{
  int j;
  if (!Aplus)
    {
      if (!(Aplus=CreateMatrix(p->npre,p->nvars)))
    return (TMatrix)NULL;
    }
  else
    if ((RowM(Aplus) != p->npre) || (ColM(Aplus) != p->nvars))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }
  for (j=0; j < p->nvars; j++)
    memcpy(&ElementM(Aplus,0,j),pElementV(p->Aplus[j][p->coef_states[j][k]]),p->npre*sizeof(PRECISION));
  return Aplus;
}

/*
   Assumes:
     Aplus : Matrix array of length n_states or null pointer.  Aplus[s] is either
             p->npre x p->nvars matrix or null pointer
*/
TMatrix* MakeAplus_All(TMatrix *Aplus, T_VAR_Parameters *p)
{
  int s;
  TMatrix *Aplus_in=Aplus;
  if (!Aplus)
    {
      if (!(Aplus=dw_CreateArray_matrix(p->nstates)))
    return (TMatrix*)NULL;
    }
  else
    if ((dw_DimA(Aplus) != p->nstates))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix*)NULL;
      }
  for (s=p->nstates-1; s >= 0; s--)
    if (!(Aplus[s]=MakeAplus(Aplus[s],s,p)))
      {
    if (Aplus_in != Aplus) dw_FreeArray(Aplus);
    return (TMatrix*)NULL;
      }
  return Aplus;
}


TMatrix MakeZeta(TMatrix Zeta, int k, T_VAR_Parameters *p)
{
  int j;
  if (!Zeta)
    {
      if (!(Zeta=CreateMatrix(p->nvars,p->nvars)))
    return (TMatrix)NULL;
    }
  else
    if ((RowM(Zeta) != p->nvars) || (ColM(Zeta) != p->nvars))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }
  InitializeMatrix(Zeta,0.0);
  for (j=0; j < p->nvars; j++)
    ElementM(Zeta,j,j)=p->Zeta[j][p->var_states[j][k]];
  return Zeta;
}

/*
   Assumes:
     Zeta : Matrix array of length n_states or null pointer.  Zeta[s] is either
             p->vars x p->nvars matrix or null pointer
*/
TMatrix* MakeZeta_All(TMatrix *Zeta, T_VAR_Parameters *p)
{
  int s;
  TMatrix *Zeta_in=Zeta;
  if (!Zeta)
    {
      if (!(Zeta=dw_CreateArray_matrix(p->nstates)))
    return (TMatrix*)NULL;
    }
  else
    if ((dw_DimA(Zeta) != p->nstates))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix*)NULL;
      }
  for (s=p->nstates-1; s >= 0; s--)
    if (!(Zeta[s]=MakeZeta(Zeta[s],s,p)))
      {
    if (Zeta_in != Zeta) dw_FreeArray(Zeta);
    return (TMatrix*)NULL;
      }
  return Zeta;
}

/*
   Assumes
     X:  n x n matrix or null pointer in column major format
     Y:  m x n matrix in column major format
     S:  m x m symmetric matrix in column major format

   Returns
     X = Y'*S*Y.  If X is null, it is created.

   Notes:
     The matrix X must be distinct from either Y or S.
*/
TMatrix MatrixInnerProductSymmetric(TMatrix X, TMatrix Y, TMatrix S)
{
  PRECISION *x, *y, *z, *s, w;
  int m, n, i, j, sj, k;
  if (!Y || !S)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  if (!X)
    {
      if (!(X=CreateMatrix(ColM(Y),ColM(Y)))) return (TMatrix)NULL;
    }
  else
    if ((RowM(S) != RowM(Y)) || (ColM(S) != RowM(Y)) || (RowM(X) != ColM(Y)) || (ColM(X) != ColM(Y)))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }

  InitializeMatrix(X,0.0);
  m=RowM(Y); n=ColM(Y);
  for (i=n-1, x=pElementM(X)+i, y=pElementM(Y)+i*m; i >= 0; x--, y-=m, i--)
    {
      for (s=pElementM(S)+(sj=m-1)*m; sj >= 0; z--, s-=m, sj--)
    {
      for (w=0.0, k=m-1; k >= 0; k--)
        w+=y[k]*s[k];
      z=pElementM(Y)+(n-1)*m+sj;
      for (j=n*(n-1); j >= 0; z-=m, j-=n)
        x[j]+=w*(*z);
    }
    }
  return X;
}

/*
   Assumes
     x : m-vector
     S : m x m symmetric matrix in column major format

   Results
     returns x'*S*x
*/
PRECISION InnerProductSymmetric(TVector x, TMatrix S)
{
  PRECISION *s, result=0.0, tmp;
  int i, j;
  if ((DimV(x) != RowM(S)) || (DimV(x) != ColM(S)))
    {
      dw_Error(SIZE_ERR);
      return 0.0;
    }
  for (s=pElementM(S), j=0; j < DimV(x); j++)
    {
      for (tmp=0.0, i=0; i < j; s++, i++) tmp+=ElementV(x,i)*(*s);
      result+=(2.0*tmp + ElementV(x,j)*(*s))*ElementV(x,j);
      s+=(DimV(x) - j);
    }
 return result;
}

/*
   Assumes
     x : m-vector
     y : n-vector
     S : m x n matrix in column major format

   Results
     returns x'*S*y
*/
PRECISION InnerProductNonSymmetric(TVector x, TVector y, TMatrix S)
{
  PRECISION *s, result=0.0, tmp;
  int i, j;
  if ((DimV(x) != RowM(S)) || (DimV(y) != ColM(S)))
    {
      dw_Error(SIZE_ERR);
      return 0.0;
    }
  for (s=pElementM(S)+DimV(x)*DimV(y)-1, j=DimV(y)-1; j >= 0; j--)
    {
      for (tmp=0.0, i=DimV(x)-1; i >= 0; s--, i--) tmp+=ElementV(x,i)*(*s);
      result+=tmp*ElementV(y,j);
    }
  return result;
}


/*
   Assumes
     x : n-vector or null pointer
     b : n-vector
     S : n x n symmetric and positive definite matrix

   Results
     The vector x is drawn from a multivariate normal distribution with

                         variance = Inverse(S)
     and
                             mean = Inverse(S)*b

     If x is null, it is created.  The matrix S is modified.

   Returns
     The vector x upon success and null upon failure.

   Notes
     Uses the Cholesky decomposition of S and the function
     DrawNormal_UpperTriangular().
*/
TVector DrawNormal_InverseVariance(TVector x, TVector b, TMatrix S)
{
  int terminal_error=dw_SetTerminalErrors(0);
  TMatrix U=CholeskyUT((TMatrix)NULL,S);
  dw_SetTerminalErrors(terminal_error);
  if (U)
    {
      x=DrawNormal_InverseUpperTriangular(x,b,U);
      FreeMatrix(U);
      return x;
    }
  else
    return (TVector)NULL;
/*      //return DrawNormal_InverseVariance_SVD(x,b,S);   ansi-c*/
}

/*
   Assumes
     x : n-vector or null pointer
     b : n-vector
     S : n x n symmetric and positive definite matrix

   Results
     The vector x is drawn from a multivariate normal distribution with

                         variance = Inverse(S)
     and
                             mean = Inverse(S)*b

     If x is null, it is created.  The matrix S is modified.

   Returns
     The vector x upon success and null upon failure.

   Notes
     Uses the Singular value decomposition of S to compute the square root of the
     inverse of S.  If S = A'*A, and c is drawn from a standard normal
     distribution, then

                          Inverse(A)*(c + Inverse(A')*b)

     is drawn from the required distribution.
*/
TVector DrawNormal_InverseVariance_SVD(TVector x, TVector b, TMatrix S)
{
  PRECISION tol, scale;
  TMatrix U, V;
  TVector d, rtrn;
  int i, j, n=DimV(b);

  _VAR_IMPROPER_DISTRIBUTION_COUNTER++;

  SVD(U=CreateMatrix(n,n),d=CreateVector(n),V=CreateMatrix(n,n),S);
  for (tol=ElementV(d,0), i=n-1; i > 0; i--)
    if (tol < ElementV(d,i)) tol=ElementV(d,i);
/*    //tol*=n*MACHINE_EPSILON;   ansi-c*/
  tol*=SQRT_MACHINE_EPSILON;
  for (j=n-1; j >= 0; j--)
    {
      scale=(ElementV(d,j) < tol) ? 1.0/sqrt(tol) : 1.0/sqrt(ElementV(d,j));
      for (i=n-1; i >= 0; i--)
    ElementM(V,i,j)*=scale;
    }
  rtrn=ProductTransposeVM(x,b,V);
  for (i=n-1; i >= 0; i--) ElementV(rtrn,i)+=dw_gaussian_rnd();
  ProductMV(rtrn,V,rtrn);
  FreeMatrix(U);
  FreeMatrix(V);
  FreeVector(d);
  return rtrn;
}

/*
   Assumes
     x : n-vector or null pointer
     b : n-vector
     U : n x n upper triangular matrix with non-zero diagonal

   Results
     The vector x is drawn from a multivariate normal distribution with

                         variance = Inverse(U'*U)
     and
                             mean = Inverse(U'*U)*b

    If x is null, it is created.

   Returns
     The vector x upon success and null upon failure.

   Notes
     If c is drawn from a standard normal distribution, then

                          Inverse(U)*(c + Inverse(U')*b)

     is drawn from the required distribution.  If S=U'U, then x is drawn from the
     multivariate normal distribution with mean Inverse(S)*b and variance
     Inverse(S).  The matrix U can be obtained from S by calling CholeskyUT(U,S).
*/
TVector DrawNormal_InverseUpperTriangular(TVector x, TVector b, TMatrix U)
{
  int i;
  TVector rtrn=ProductInverseVU(x,b,U);
  if (rtrn)
    {
      for (i=DimV(rtrn)-1; i >= 0; i--) ElementV(rtrn,i)+=dw_gaussian_rnd();
      if (!InverseProductUV(rtrn,U,rtrn))
    {
      if (!x) FreeVector(rtrn);
      return (TVector)NULL;
    }
    }
  return rtrn;
}


/*
   Attempts recovery from a singular inverse variance matrix.
*/
TVector SingularInverseVariance_RecoveryAttempt(TVector x, TVector b, TMatrix S, TMatrix InversePrior, TStateModel *model, int code)
{
  FILE *f_out;
  char *header;
  char filename[256];

/*    // print warning   ansi-c*/
  printf("Singular error\n");

/*    // Construct file name and open file   ansi-c*/
  if (!V_FILE)
    {
      sprintf(filename,"singular_error_%03d_test.dat",_VAR_IMPROPER_DISTRIBUTION_COUNTER);
      f_out=dw_CreateTextFile(filename);
    }
  else
    f_out=V_FILE;

  fprintf(f_out,"//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//\n");

/*    // Print error message, inverse variance matrix, and inverse prior   ansi-c*/
  switch (code)
    {
    case BPLUS_ERR:
      fprintf(f_out,"Error in DrawAplus(): Inverse of variance is singular\n");
      break;
    case PSI_ERR:
      fprintf(f_out,"Error in Draw_psi(): Inverse of variance is singular\n");
      break;
    case LAMBDA_ERR:
      fprintf(f_out,"Error in Draw_Lambda(): Inverse of variance is singular\n");
      break;
    default:
      fprintf(f_out,"Unknown routine code: Inverse of variance is singular\n");
      break;
    }
  fprintf(f_out,"Inverse variance =\n");
  dw_PrintMatrix(f_out,S,"%lg ");
  fprintf(f_out,"Inverse prior =\n");
  dw_PrintMatrix(f_out,InversePrior,"%lg ");
  fprintf(f_out,"\n");

/*    // Print generator state   ansi-c*/
  fprintf(f_out,"\\== Generator State ==\\\n");
  dw_print_generator_state(f_out);
  fprintf(f_out,"\n");

/*    // Attempt recovery   ansi-c*/
  x=DrawNormal_InverseVariance_SVD(x,b,S);

/*    // Print parameters   ansi-c*/
  header="Error draw: ";
  WriteStates(f_out,(char*)NULL,header,model);
  WriteTransitionMatrices(f_out,(char*)NULL,header,model);
  Write_VAR_Parameters(f_out,(char*)NULL,header,model);

  fprintf(f_out,"//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//\n");

/*    // Close file   ansi-c*/
  if (f_out != V_FILE) fclose(f_out);

  _SINGULAR_ERROR=1;

  return x;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
