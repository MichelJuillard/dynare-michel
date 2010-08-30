
#include "VARio.h"
#include "switchio.h"
#include "dw_error.h"
#include "dw_ascii.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "modify_for_mex.h"

static int strlen_int(int n)
{
  int i, j;
  for (i=1, j=10; n >= j; i++, j*=10);
  return i;
}

static void ReadError_VARio(char *id)
{
  char *errmsg, *fmt="Error after line identifier ""%s""";
  sprintf(errmsg=(char*)swzMalloc(strlen(fmt) + strlen(id) - 1),fmt,id);
  dw_UserError(errmsg);
  swzFree(errmsg);
}

static int ReadInteger_VARio(FILE *f_in, char *id)
{
  int i;
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %d ",&i) != 1)) ReadError_VARio(id);
  return i;
}

static PRECISION ReadScalar_VARio(FILE *f_in, char *id)
{
  double x;
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %lf ",&x) != 1)) ReadError_VARio(id);
  return (PRECISION)x;
}

static void ReadMatrix_VARio(FILE *f_in, char *id, TMatrix X)
{
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,X)) ReadError_VARio(id);
}

static void ReadVector_VARio(FILE *f_in, char *id, TVector X)
{
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadVector(f_in,X)) ReadError_VARio(id);
}

static void ReadArray_VARio(FILE *f_in, char *id, void *X)
{
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,X)) ReadError_VARio(id);
}

static FILE* OpenFile_VARio(FILE *f, char *filename)
{
  char *errmsg, *fmt="Unable to open %s";
  if (!f)
    if (!filename)
      dw_UserError("Filename pointer is null.");
    else
      if (!(f=fopen(filename,"rt")))
        {
          sprintf(errmsg=(char*)swzMalloc(strlen(fmt) + strlen(filename) - 1),fmt,filename);
          dw_UserError(errmsg);
          swzFree(errmsg);
        }
  return f;
}

/*
   Assumes:
    f:  valid file pointer or null
    filename:  pointer to null terminated string or null

   Returns:
    A pointer to a valid TStateModel upon success and null pointer upon failure.
    Upon failure, the routine prints an error message if USER_ERR is a verbose
    error and terminates if USER_ERR is a terminal error.  The terminal errors
    and verbose errors can be set with dw_SetTerminalErrors() and
    dw_SetVerboseErrors().

   Results:
    Upon success, a valid TStateModel is created and initialized.

   Notes:
    One of f and filename must not be null.
*/
TStateModel* Read_VAR_Specification(FILE *f, char *filename)
{
  TMarkovStateVariable *sv;
  T_VAR_Parameters *p;
  char *id, *fmt;
  int *IV;
  int j, spec, nvars, nlags, nexg, npre, nstates, nobs;
  PRECISION lambda_prior;
  TVector zeta_a_prior, zeta_b_prior;
  TMatrix *U, *V, *W, *A0_prior, *Aplus_prior, Y, X;
  int **coef_states, **var_states;
  PRECISION** A0_Metropolis_Scale=(PRECISION**)NULL;

/*    // Valid file   ansi-c*/
   FILE *f_in=OpenFile_VARio(f,filename);
  if (!f_in) return (TStateModel*)NULL;

/*    // Read Markov specifications   ansi-c*/
  sv=ReadMarkovSpecification(f_in,(char*)NULL);

/*    //=== Sizes ===//   ansi-c*/
  nvars=ReadInteger_VARio(f_in,"//== Number Variables ==//");
  nlags=ReadInteger_VARio(f_in,"//== Number Lags ==//");
  nexg=ReadInteger_VARio(f_in,"//== Exogenous Variables ==//");
  nstates=ReadInteger_VARio(f_in,"//== Number States ==//");
  nobs=ReadInteger_VARio(f_in,"//== Number Observations ==//");
  npre=nvars*nlags+nexg;
  if ((nobs != sv->nobs) || (nstates != sv->nstates))
    {
      dw_UserError("Read_VAR_Specification():  different values for nobs or nstates.");
      return (TStateModel*)NULL;
    }

/*    //=== Restrictions - U[j] ===//   ansi-c*/
  ReadArray_VARio(f_in,"//== Number of free parameters in each column of A0 ==//",IV=dw_CreateArray_int(nvars));
  U=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    {
      fmt="//== U[%d] ==//";
      sprintf(id=(char*)swzMalloc(strlen(fmt) + strlen_int(j+1) - 1),fmt,j+1);
      ReadMatrix_VARio(f_in,id,U[j]=CreateMatrix(nvars,IV[j]));
      swzFree(id);
    }
  dw_FreeArray(IV);

/*    //=== Restrictions - V[j] ===//   ansi-c*/
  ReadArray_VARio(f_in,"//== Number of free parameters in each column of Aplus ==//",IV=dw_CreateArray_int(nvars));
  V=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (IV[j] > 0)
      {
    fmt="//== V[%d] ==//";
    sprintf(id=(char*)swzMalloc(strlen(fmt) + strlen_int(j+1) - 1),fmt,j+1);
    ReadMatrix_VARio(f_in,id,V[j]=CreateMatrix(npre,IV[j]));
    swzFree(id);
      }
  dw_FreeArray(IV);

/*    //=== Restrictions - W[j] ===//   ansi-c*/
  ReadArray_VARio(f_in,"//== Non-zero W[j] ==//",IV=dw_CreateArray_int(nvars));
  W=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (IV[j])
      {
    fmt="//== W[%d] ==//";
    sprintf(id=(char*)swzMalloc(strlen(fmt) + strlen_int(j+1) - 1),fmt,j+1);
    ReadMatrix_VARio(f_in,id,W[j]=CreateMatrix(npre,nvars));
    swzFree(id);
      }
  dw_FreeArray(IV);

/*    //====== Priors ======   ansi-c*/
  ReadVector_VARio(f_in,"//== Gamma prior on zeta - a ==//",zeta_a_prior=CreateVector(nvars));
  ReadVector_VARio(f_in,"//== Gamma prior on zeta - b ==//",zeta_b_prior=CreateVector(nvars));

  A0_prior=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    {
      fmt="//== Variance of Gaussian prior on column %d of A0 ==//";
      sprintf(id=(char*)swzMalloc(strlen(fmt) + strlen_int(j+1) - 1),fmt,j+1);
      ReadMatrix_VARio(f_in,id,A0_prior[j]=CreateMatrix(nvars,nvars));
      swzFree(id);
    }

  Aplus_prior=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    {
      fmt="//== Variance of Gaussian prior on column %d of Aplus ==//";
      sprintf(id=(char*)swzMalloc(strlen(fmt) + strlen_int(j+1) - 1),fmt,j+1);
      ReadMatrix_VARio(f_in,id,Aplus_prior[j]=CreateMatrix(npre,npre));
      swzFree(id);
    }

/*    //=== Specification ===//   ansi-c*/
  spec=ReadInteger_VARio(f_in,"//== Specification (0=default  1=Sims-Zha  2=Random Walk) ==//");
  switch (spec)
    {
    case 0: spec=0; break;
    case 1: spec=SPEC_SIMS_ZHA | SPEC_RANDOM_WALK; break;
    case 2: spec=SPEC_RANDOM_WALK; break;
    default: ReadError_VARio("//== Specification (0=default  1=Sims-Zha  2=Random Walk) ==//"); swzExit(0);
    }
  if (spec & SPEC_SIMS_ZHA)
    lambda_prior=ReadScalar_VARio(f_in,"//== Variance of Gaussian prior on lambda ==//");

/*    //====== coefficient and variance state variables ======   ansi-c*/
  ReadArray_VARio(f_in,"//== Translation table for coefficient states ==//",coef_states=dw_CreateRectangularArray_int(nvars,nstates));
  ReadArray_VARio(f_in,"//== Translation table for variance states ==//",var_states=dw_CreateRectangularArray_int(nvars,nstates));

/*    //====== Metropolis jumping kernel info for A0 ======   ansi-c*/
  if (dw_SetFilePosition(f_in,"//== Metropolis kernel scales for A0 ==//"))
    {
      A0_Metropolis_Scale=dw_CreateArray_array(nvars);
      for (j=nvars-1; j >= 0; j--)
    A0_Metropolis_Scale[j]=dw_CreateArray_scalar(GetNumberStatesFromTranslationMatrix(j,coef_states));
      if (!dw_ReadArray(f_in,A0_Metropolis_Scale)) ReadError_VARio(id);
    }

/*    //=== Data  ===   ansi-c*/
  ReadMatrix_VARio(f_in,"//== Data Y (nobs x nvars) ==//",Y=CreateMatrix(nobs,nvars));
  ReadMatrix_VARio(f_in,"//== Data X (nobs x npre) ==//",X=CreateMatrix(nobs,npre));

/*    //=== Create T_VAR_Parameters structure ===   ansi-c*/
  p=CreateTheta_VAR(spec,nvars,nlags,nexg,nstates,nobs,coef_states,var_states,U,V,W,Y,X);
  if (spec & SPEC_SIMS_ZHA)
    SetPriors_VAR_SimsZha(p,A0_prior,Aplus_prior,zeta_a_prior,zeta_b_prior,lambda_prior);
  else
    SetPriors_VAR(p,A0_prior,Aplus_prior,zeta_a_prior,zeta_b_prior);

  if (A0_Metropolis_Scale) SetupMetropolisInformation(A0_Metropolis_Scale,p);

/*    //=== Close output file ===   ansi-c*/
  if (!f) fclose(f_in);

/*    //=== Free memory ===   ansi-c*/
  dw_FreeArray(U);
  dw_FreeArray(V);
  dw_FreeArray(W);
  FreeVector(zeta_a_prior);
  FreeVector(zeta_b_prior);
  dw_FreeArray(A0_prior);
  dw_FreeArray(Aplus_prior);
  FreeMatrix(X);
  FreeMatrix(Y);
  dw_FreeArray(coef_states);
  dw_FreeArray(var_states);
  dw_FreeArray(A0_Metropolis_Scale);

/*    //=== return TStateModel structure ===   ansi-c*/
  return CreateStateModel_new(sv,CreateRoutines_VAR(),p);
}

/*
   Writes the specification
*/
void Write_VAR_Specification(FILE *f, char *filename, TStateModel *model)
{
  int j, t;
  FILE *f_out=f ? f : dw_CreateTextFile(filename);
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

/*    // Write Markov specifications   ansi-c*/
  WriteMarkovSpecification(f_out,(char*)NULL,model);

/*    //=== Sizes ===//   ansi-c*/
  fprintf(f_out,"//== Number Variables ==//\n%d\n\n",p->nvars);
  fprintf(f_out,"//== Number Lags ==//\n%d\n\n",p->nlags);
  fprintf(f_out,"//== Exogenous Variables ==//\n%d\n\n",p->npre - p->nvars * p->nlags);
  fprintf(f_out,"//== Number States ==//\n%d\n\n",p->nstates);
  fprintf(f_out,"//== Number Observations ==//\n%d\n\n",p->nobs);

/*    //=== Restrictions - U[j] ===//   ansi-c*/
  fprintf(f_out,"//== Number of free parameters in each column of A0 ==//\n");
  for (j=0; j < p->nvars; j++)
    fprintf(f_out,"%d ",ColM(p->U[j]));
  fprintf(f_out,"\n\n");
  for (j=0; j < p->nvars; j++)
    {
      fprintf(f_out,"//== U[%d] ==//\n",j+1);
      dw_PrintMatrix(f_out,p->U[j],"%22.14le ");
      fprintf(f_out,"\n");
    }

/*    //=== Restrictions - V[j] ===//   ansi-c*/
  fprintf(f_out,"//== Number of free parameters in each column of Aplus ==//\n");
  for (j=0; j < p->nvars; j++)
    fprintf(f_out,"%d ",p->V[j] ? ColM(p->V[j]) : 0);
  fprintf(f_out,"\n\n");
  for (j=0; j < p->nvars; j++)
    if (p->V[j])
      {
    fprintf(f_out,"//== V[%d] ==//\n",j+1);
    dw_PrintMatrix(f_out,p->V[j],"%22.14le ");
    fprintf(f_out,"\n");
      }

/*    //=== Restrictions - W[j] ===//   ansi-c*/
  fprintf(f_out,"//== Non-zero W[j] ==//\n");
  for (j=0; j < p->nvars; j++)
    fprintf(f_out,"%d ",p->W[j] ? 1 : 0);
  fprintf(f_out,"\n\n");
  for (j=0; j < p->nvars; j++)
    if (p->W[j])
      {
    fprintf(f_out,"//== W[%d] ==//\n",j+1);
    dw_PrintMatrix(f_out,p->W[j],"%22.14le ");
    fprintf(f_out,"\n");
      }

/*    //====== Priors ======   ansi-c*/
  fprintf(f_out,"//== Gamma prior on zeta - a ==//\n");
  dw_PrintVector(f_out,p->zeta_a_prior,"%22.14le ");
  fprintf(f_out,"\n");
  fprintf(f_out,"//== Gamma prior on zeta - b ==//\n");
  dw_PrintVector(f_out,p->zeta_b_prior,"%22.14le ");
  fprintf(f_out,"\n");

  for (j=0; j < p->nvars; j++)
    {
      fprintf(f_out,"//== Variance of Gaussian prior on column %d of A0 ==//\n",j+1);
      dw_PrintMatrix(f_out,p->A0_prior[j],"%22.14le ");
      fprintf(f_out,"\n");
    }

  for (j=0; j < p->nvars; j++)
    {
      fprintf(f_out,"//== Variance of Gaussian prior on column %d of Aplus ==//\n",j+1);
      dw_PrintMatrix(f_out,p->Aplus_prior[j],"%22.14le ");
      fprintf(f_out,"\n");
    }

/*    //=== Model specification ===//   ansi-c*/
  fprintf(f_out,"//== Specification (0=default  1=Sims-Zha  2=Random Walk) ==//\n");
  if (p->Specification & SPEC_SIMS_ZHA)
    fprintf(f_out,"1\n\n");
  else
    if (p->Specification & SPEC_RANDOM_WALK)
      fprintf(f_out,"2\n\n");
    else
      fprintf(f_out,"0\n\n");
  if ((p->Specification & SPEC_SIMS_ZHA) == SPEC_SIMS_ZHA)
    fprintf(f_out,"//== Variance of Gaussian prior on lambda ==//\n%22.14le\n\n",p->lambda_prior);

/*    //====== coefficient and variance state variables ======   ansi-c*/
  fprintf(f_out,"//== Translation table for coefficient states ==//\n");
  dw_PrintArray(f_out,p->coef_states,"%4d ");

  fprintf(f_out,"//== Translation table for variance states ==//\n");
  dw_PrintArray(f_out,p->var_states,"%4d ");

/*    //====== Metropolis jumping kernel info for A0 ======   ansi-c*/
  fprintf(f_out,"//== Metropolis kernel scales for A0 ==//\n");
  dw_PrintArray(f_out,p->A0_Metropolis_Scale,"%22.14le ");

/*    //=== Data  ===   ansi-c*/
  fprintf(f_out,"//== Data Y (nobs x nvars) ==//\n");
  for (t=1; t <= p->nobs; t++)
    dw_PrintVector(f_out,p->Y[t],"%22.14le ");
  fprintf(f_out,"\n");

  fprintf(f_out,"//== Data X (nobs x npre) ==//\n");
  for (t=1; t <= p->nobs; t++)
    dw_PrintVector(f_out,p->X[t],"%22.14le ");
  fprintf(f_out,"\n");

/*    //=== Close output file ===   ansi-c*/
  if (!f) fclose(f_out);
}

/*
   Assumes:
    f:  valid file pointer or null
    filename:  pointer to null terminated string or null
    model:  pointer to valid TStateModel structure

   Returns:
    One upon success.  Upon failure, the routine prints an error message if
    USER_ERR is a verbose error, terminates if USER_ERR is a terminal error and
    returns zero if USER_ERR is not a terminal error.  The terminal errors and
    verbose errors can be set with dw_SetTerminalErrors() and
    dw_SetVerboseErrors().

   Results:
    Upon success, the following fields of p will be filled:

                         A0, Aplus, Zeta, b0, bplus.

    If the Sims-Zha specification is used, the following fields will also be
    filled

                               lambda, psi.

    The routine Thetahanged() will be called.

   Notes:
    One of f and filename must not be null.

    The file must contain line identifiers of the form

                           //== A0[s] ==//
                           //== Aplus[s] ==//
                           //== Zeta[s] ==//

    for 1 <= s <= p->nstates.

    Zeta is checked for non-negativity.  No checks are made to ensure that A0[s],
    Aplus[s], or Zeta[s] satisfy any restrictions.
*/
int Read_VAR_Parameters(FILE *f, char *filename, char *header, TStateModel *model)
{
  FILE *f_in;
  char *idbuffer, *fmt;
  TMatrix *A0, *Aplus, *Zeta;
  int i, j, s;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

/*    // Valid file   ansi-c*/
  f_in=OpenFile_VARio(f,filename);
  if (!f_in) return 0;

  if (!header) header="";

/*    // Allocate memory   ansi-c*/
  A0=dw_CreateArray_matrix(p->nstates);
  Aplus=dw_CreateArray_matrix(p->nstates);
  Zeta=dw_CreateArray_matrix(p->nstates);

/*    // Read File   ansi-c*/
  for (s=0; s < p->nstates; s++)
    {
      fmt="//== %sA0[%d] ==//";
      sprintf(idbuffer=(char*)swzMalloc(strlen(fmt)+strlen(header)+strlen_int(s+1)-3),fmt,header,s+1);
      if (!dw_SetFilePosition(f_in,idbuffer) || !dw_ReadMatrix(f_in,A0[s]=CreateMatrix(p->nvars,p->nvars)))
    {
      ReadError_VARio(idbuffer);
      swzFree(idbuffer);
      return 0;
    }
      swzFree(idbuffer);

      fmt="//== %sAplus[%d] ==//";
      sprintf(idbuffer=(char*)swzMalloc(strlen(fmt)+strlen(header)+strlen_int(s+1)-3),fmt,header,s+1);
      if (!dw_SetFilePosition(f_in,idbuffer) || !dw_ReadMatrix(f_in,Aplus[s]=CreateMatrix(p->npre,p->nvars)))
    {
      ReadError_VARio(idbuffer);
      swzFree(idbuffer);
      return 0;
    }
      swzFree(idbuffer);

      fmt="//== %sZeta[%d] ==//";
      sprintf(idbuffer=(char*)swzMalloc(strlen(fmt)+strlen(header)+strlen_int(s+1)-3),fmt,header,s+1);
      if (!dw_SetFilePosition(f_in,idbuffer) || !dw_ReadMatrix(f_in,Zeta[s]=CreateMatrix(p->nvars,p->nvars)))
    {
      ReadError_VARio(idbuffer);
      swzFree(idbuffer);
      return 0;
    }
      swzFree(idbuffer);
    }

/*    // Set A0, Aplus, and Zeta   ansi-c*/
  for (j=0; j < p->nvars; j++)
    for (s=0; s < p->nstates; s++)
      {
    for (i=0; i < p->nvars; i++)
      ElementV(p->A0[j][p->coef_states[j][s]],i)=ElementM(A0[s],i,j);

    for (i=0; i < p->npre; i++)
      ElementV(p->Aplus[j][p->coef_states[j][s]],i)=ElementM(Aplus[s],i,j);

    p->Zeta[j][p->var_states[j][s]]=ElementM(Zeta[s],j,j);
      }

/*    // Free memory   ansi-c*/
  dw_FreeArray(A0);
  dw_FreeArray(Aplus);
  dw_FreeArray(Zeta);

/*    // Check Zeta non-negative   ansi-c*/
  for (j=p->nvars-1; j >= 0; j--)
    for (s=p->n_var_states[j]-1; s >= 0; s--)
      if (p->Zeta[j][s] < 0.0)
    {
      dw_UserError("Zeta has negative value.");
      p->valid_parameters=0;
      ThetaChanged(model);
      return 0;
    }

/*    // Update b0, bplus, lambda, psi   ansi-c*/
  Update_b0_bplus_from_A0_Aplus(p);
  if ((p->Specification & SPEC_SIMS_ZHA) == SPEC_SIMS_ZHA) Update_lambda_psi_from_bplus(p);

/*    // Flags and notification that the VAR parameters have changed   ansi-c*/
  p->valid_parameters=1;
  ThetaChanged(model);

  return 1;
}

/*
   Writes the VAR parameters to a file.  The identifiers are

     //== A0[s] ==//
     //== Aplus[s] ==//
     //== Zeta[s] ==//

   for 1 <= s <= nstates
*/
int Write_VAR_Parameters(FILE *f, char *filename, char *header, TStateModel *model)
{
  TMatrix X;
  int s;
  FILE *f_out;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  f_out=f ? f :dw_CreateTextFile(filename);

  if (!header) header="";

  for (s=0; s < p->nstates; s++)
    {
      X=MakeA0((TMatrix)NULL,s,p);
      fprintf(f_out,"//== %sA0[%d] ==//\n",header,s+1);
      dw_PrintMatrix(f_out,X,"%22.14le ");
      fprintf(f_out,"\n");
      FreeMatrix(X);

      X=MakeAplus((TMatrix)NULL,s,p);
      fprintf(f_out,"//== %sAplus[%d] ==//\n",header,s+1);
      dw_PrintMatrix(f_out,X,"%22.14le ");
      fprintf(f_out,"\n");
      FreeMatrix(X);

      X=MakeZeta((TMatrix)NULL,s,p);
      fprintf(f_out,"//== %sZeta[%d] ==//\n",header,s+1);
      dw_PrintMatrix(f_out,X,"%22.14le ");
      fprintf(f_out,"\n");
      FreeMatrix(X);
    }

  if (!f) fclose(f_out);

  return 1;
}

/*
   Writes the headers for Write_VAR_ParametersFlat().  This routine can
   be used to give the ordering for Write_VAR_ParametersFlat().
*/
int Write_VAR_ParametersFlat_Headers(FILE *f_out, TStateModel *model)
{
  int i, j, s;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  if (!f_out) return 0;

  for (s=0; s < p->nstates; s++)
    {
      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->nvars; i++)
      fprintf(f_out,"A0[%d](%d,%d) ",s+1,i+1,j+1);

      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->npre; i++)
      fprintf(f_out,"Aplus[%d](%d,%d) ",s+1,i+1,j+1);

      for (j=0; j < p->nvars; j++)
    fprintf(f_out,"Zeta[%d](%d,%d) ",s+1,j+1,j+1);
    }

  return 1;
}

int Read_VAR_ParametersFlat(FILE *f_in, TStateModel *model)
{
  TMatrix *A0, *Aplus;
  TVector *Zeta;
  int i, j, s, rtrn=0;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  // Allocate memory
  A0=dw_CreateArray_matrix(p->nstates);
  Aplus=dw_CreateArray_matrix(p->nstates);
  Zeta=dw_CreateArray_vector(p->nstates);

  // Read File
  for (s=0; s < p->nstates; s++)
    {
      A0[s]=CreateMatrix(p->nvars,p->nvars);
      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->nvars; i++)
      if (fscanf(f_in," %lf ",&ElementM(A0[s],i,j)) != 1)
        goto ERROR;

      Aplus[s]=CreateMatrix(p->npre,p->nvars);
      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->npre; i++)
      if (fscanf(f_in," %lf ",&ElementM(Aplus[s],i,j)) != 1)
        goto ERROR;

      Zeta[s]=CreateVector(p->nvars);
      for (j=0; j < p->nvars; j++)
    if (fscanf(f_in," %lf ",&ElementV(Zeta[s],j)) != 1)
      goto ERROR;
    else
      if (ElementV(Zeta[s],j) < 0.0)
        goto ERROR;
    }

  // Set A0, Aplus, and Zeta
  for (j=0; j < p->nvars; j++)
    for (s=0; s < p->nstates; s++)
      {
    for (i=0; i < p->nvars; i++)
      ElementV(p->A0[j][p->coef_states[j][s]],i)=ElementM(A0[s],i,j);

    for (i=0; i < p->npre; i++)
      ElementV(p->Aplus[j][p->coef_states[j][s]],i)=ElementM(Aplus[s],i,j);

    p->Zeta[j][p->var_states[j][s]]=ElementV(Zeta[s],j);
      }

  // Update b0, bplus, lambda, psi
  Update_b0_bplus_from_A0_Aplus(p);
  if ((p->Specification & SPEC_SIMS_ZHA) == SPEC_SIMS_ZHA) Update_lambda_psi_from_bplus(p);

  // Flags and notification that the VAR parameters have changed
  p->valid_parameters=1;
  ThetaChanged(model);
  rtrn=1;

 ERROR:

  // Free memory
  dw_FreeArray(A0);
  dw_FreeArray(Aplus);
  dw_FreeArray(Zeta);

  return rtrn;
}

/*
   For each state the VAR parameters are printed as follows
    A0    (by columns)
    Aplus (by columns)
    Zeta  (diagonal)
*/
int Write_VAR_ParametersFlat(FILE *f, TStateModel *model, char *fmt)
{
  TMatrix A0, Aplus;
  int s, i, j;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  if (!f) return 0;

  if (!fmt) fmt="%lf ";

  A0=CreateMatrix(p->nvars,p->nvars);
  Aplus=CreateMatrix(p->npre,p->nvars);

  for (s=0; s < p->nstates; s++)
    {
      MakeA0(A0,s,p);
      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->nvars; i++)
      fprintf(f,fmt,ElementM(A0,i,j));


      MakeAplus(Aplus,s,p);
      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->npre; i++)
      fprintf(f,fmt,ElementM(Aplus,i,j));

      for (j=0; j < p->nvars; j++)
    fprintf(f,fmt,p->Zeta[j][p->var_states[j][s]]);
    }

  FreeMatrix(Aplus);
  FreeMatrix(A0);

  return 1;
}

/*
   For each state the VAR parameters are printed as follows
    A0    (by columns)
    Aplus (by columns)
    Zeta  (diagonal)
   The system is normalized so that the diagonal of A0 is one.
*/
int Write_VAR_ParametersFlat_A0_Diagonal_One(FILE *f, TStateModel *model, char *fmt)
{
  TMatrix A0, Aplus;
  TVector i_diagonal, s_diagonal;
  int s, i, j;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);
  PRECISION x;

  if (!f) return 0;

  if (!fmt) fmt="%lf ";

  A0=CreateMatrix(p->nvars,p->nvars);
  Aplus=CreateMatrix(p->npre,p->nvars);
  s_diagonal=CreateVector(p->nvars);
  i_diagonal=CreateVector(p->nvars);

  for (s=0; s < p->nstates; s++)
    {
      MakeA0(A0,s,p);
      for (i=p->nvars-1; i >= 0; i--)
    {
      ElementV(i_diagonal,i)=1.0/(x=ElementM(A0,i,i));
      ElementV(s_diagonal,i)=x*x;
    }

      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->nvars; i++)
      fprintf(f,fmt,ElementM(A0,i,j)*ElementV(i_diagonal,j));

      MakeAplus(Aplus,s,p);
      for (j=0; j < p->nvars; j++)
    for (i=0; i < p->npre; i++)
      fprintf(f,fmt,ElementM(Aplus,i,j)*ElementV(i_diagonal,j));

      for (j=0; j < p->nvars; j++)
    fprintf(f,fmt,p->Zeta[j][p->var_states[j][s]] * ElementV(s_diagonal,j));
    }

  FreeVector(i_diagonal);
  FreeVector(s_diagonal);
  FreeMatrix(Aplus);
  FreeMatrix(A0);

  return 1;
}

/*
   Attempts to read all parameters.  The identifiers are

     //== <id>States ==//
     //== <id>Transition matrix[] ==//
     //== <id>A0[s] ==//
     //== <id>Aplus[s] ==//
     //== <id>Zeta[s] ==//

   for 1 <= s <= nstates
*/
void ReadAllParameters(FILE *f, char *filename, char *id, TStateModel *model)
{
  char *buffer, *fmt="//== %sStates ==//";
  FILE *f_in=f ? f :dw_OpenTextFile(filename);

  if (!id) id="";

  sprintf(buffer=(char*)swzMalloc(strlen(fmt) + strlen(id) - 1),fmt,id);
  ReadArray_VARio(f_in,buffer,model->sv->S);
  swzFree(buffer);

  ReadTransitionMatrices(f_in,(char*)NULL,id,model);
  Read_VAR_Parameters(f_in,(char*)NULL,id,model);
}

/*
   Attempts to write all parameters using a format readable by the routine
   ReadAllParameters().
*/
void WriteAllParameters(FILE *f, char *filename, char *id, TStateModel *model)
{
  FILE *f_in=f ? f : dw_CreateTextFile(filename);

  if (!id) id="";

  fprintf(f_in,"//== %sStates ==//\n",id);
  dw_PrintArray(f_in,model->sv->S,(char*)NULL);
  fprintf(f_in,"\n");

  WriteTransitionMatrices(f_in,(char*)NULL,id,model);
  Write_VAR_Parameters(f_in,(char*)NULL,id,model);

  if(!f) fclose(f_in);
}

/*******************************************************************************/
/******************************** Input/Output *********************************/
/*******************************************************************************/
void Write_ReducedFormVAR_Parameters(FILE *f, char *filename, T_VAR_Parameters *p)
{
  TMatrix A0, Aplus, Zeta, C, Sigma;
  int k;
  FILE *f_out;

  f_out=f ? f :dw_CreateTextFile(filename);

  A0=CreateMatrix(p->nvars,p->nvars);
  Aplus=CreateMatrix(p->npre,p->nvars);
  Zeta=CreateMatrix(p->nvars,p->nvars);
  C=CreateMatrix(p->npre,p->nvars);
  Sigma=CreateMatrix(p->nvars,p->nvars);

  for (k=0; k < p->nstates; k++)
    {
      MakeA0(A0,k,p);
      MakeAplus(Aplus,k,p);
      MakeZeta(Zeta,k,p);

/*        //ProductInverseMM(C,Aplus,A0);   ansi-c*/
/*        //ProductMM(A0,A0,Xi);   ansi-c*/
/*        //ProductTransposeMM(Sigma,A0,A0);   ansi-c*/
/*        //Inverse_LU(Sigma,Sigma);   ansi-c*/

      fprintf(f_out,"//== Reduced Form[%d] ==//\n",k+1);
      dw_PrintMatrix(f_out,C,"%lf ");
      fprintf(f_out,"\n");

      fprintf(f_out,"//== Variance[%d] ==//\n",k+1);
      dw_PrintMatrix(f_out,Sigma,"%lf ");
      fprintf(f_out,"\n");
    }

  FreeMatrix(A0);
  FreeMatrix(Aplus);
  FreeMatrix(Zeta);
  FreeMatrix(C);
  FreeMatrix(Sigma);

  if (!f) fclose(f_out);
}

/*
   Create Model from data file.  Assumes that the state variables have a flat
   structure.
*/
void Write_VAR_Info(FILE *f, char *filename, T_VAR_Parameters *p)
{
  FILE *f_out;
  int j;

  if (!f)
    f_out=dw_CreateTextFile(filename);
  else
    f_out=f;

/*    //=== Write sizes ===//   ansi-c*/
  fprintf(f_out,"//== Number Observations ==//\n%d\n\n",p->nobs);
  fprintf(f_out,"//== Number Variables ==//\n%d\n\n",p->nvars);
  fprintf(f_out,"//== Number Lags ==//\n%d\n\n",p->nlags);
  fprintf(f_out,"//== Exogenous Variables ==//\n%d\n\n",p->npre - p->nvars * p->nlags);

/*    //=== Restrictions - U[j] ===//   ansi-c*/
  fprintf(f_out,"//== Number of free parameters in jth column of A0 ==//\n");
  for (j=0; j < p->nvars; j++)
    fprintf(f_out,"%d ",ColM(p->U[j]));
  fprintf(f_out,"\n\n");
  fprintf(f_out,"//== U[j] 0 <= j < nvars ==//\n");
  for (j=0; j < p->nvars; j++)
    {
      dw_PrintMatrix(f_out,p->U[j],"%lf ");
      fprintf(f_out,"\n");
    }

/*    //=== Restrictions - V[j] ===//   ansi-c*/
  fprintf(f_out,"//== Number of free parameters in jth column of Aplus ==//\n");
  for (j=0; j < p->nvars; j++)
    fprintf(f_out,"%d ",p->V[j] ? ColM(p->V[j]) : 0);
  fprintf(f_out,"\n\n");
  fprintf(f_out,"//== V[j] 0 <= j < nvars ==//\n");
  for (j=0; j < p->nvars; j++)
    if (p->V[j])
      {
    dw_PrintMatrix(f_out,p->V[j],"%lf ");
    fprintf(f_out,"\n");
      }

/*    //=== Restrictions - W[j] ===//   ansi-c*/
  fprintf(f_out,"//== Non-zero W[j] ==//\n");
  for (j=0; j < p->nvars; j++)
    fprintf(f_out,"%d ",p->W[j] ? 1 : 0);
  fprintf(f_out,"\n\n");
  fprintf(f_out,"//== W[j] 0 <= j < nvars ==//\n");
  for (j=0; j < p->nvars; j++)
    if (p->W[j])
      {
    dw_PrintMatrix(f_out,p->W[j],"%lf ");
    fprintf(f_out,"\n");
      }

/*    //====== Priors ======   ansi-c*/
  fprintf(f_out,"//== Gamma prior on Xi ==//\n");
  for (j=0; j < p->nvars; j++)
    fprintf(f_out,"%lf %lf\n",ElementV(p->zeta_a_prior,j),ElementV(p->zeta_b_prior,j));
  fprintf(f_out,"\n");

  fprintf(f_out,"//== Prior on jth column of A0 - Gaussian variance ==//\n");
  for (j=0; j < p->nvars; j++)
    {
      dw_PrintMatrix(f_out,p->A0_prior[j],"%lf ");
      fprintf(f_out,"\n");
    }

  fprintf(f_out,"//== Prior on jth column of Aplus - Gaussian variance ==//\n");
  for (j=0; j < p->nvars; j++)
    {
      dw_PrintMatrix(f_out,p->Aplus_prior[j],"%lf ");
      fprintf(f_out,"\n");
    }

/*   //====== coefficient/variance state variables ====== */
/*   CStates=dw_CreateRegularArrayList_int(2,p->nvars,sv->n_state_variables); */
/*   id="//== Controlling states variables for coefficients ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,CStates)) dw_Error(PARSE_ERR); */

/*   VStates=dw_CreateRegularArrayList_int(2,p->nvars,sv->n_state_variables); */
/*   id="//== Controlling states variables for variance ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,VStates)) dw_Error(PARSE_ERR); */

/*   //=== Read Data  === */
/*   if (!dw_SetFilePosition(f_in,"//== Data Y (T x nvars) ==//") */
/*       || !dw_SetFilePosition(f_in,"//== Data X (T x npre) ==//")) */
/*     p->X=p->Y=(TVector*)NULL; */
/*   else */
/*     { */
/*       // Initialize Y */
/*       id="//== Data Y (T x nvars) ==//"; */
/*       if (!dw_SetFilePosition(f_in,id)) dw_Error(PARSE_ERR); */
/*       p->Y=dw_CreateArray_vector(p->nobs+1); */
/*       for (t=1; t <= p->nobs; t++) */
/*     if (!dw_ReadVector(f_in,p->Y[t]=CreateVector(p->nvars))) dw_Error(PARSE_ERR); */

/*       // Initialize X */
/*       id="//== Data X (T x npre) ==//"; */
/*       if (!dw_SetFilePosition(f_in,id)) dw_Error(PARSE_ERR); */
/*       p->X=dw_CreateArray_vector(p->nobs+1); */
/*       for (t=1; t <= p->nobs; t++) */
/*     if (!dw_ReadVector(f_in,p->X[t]=CreateVector(p->npre))) */
/*       dw_Error(PARSE_ERR); */
/*     } */

/*    //=== Close output file ===   ansi-c*/
  if (!f) fclose(f_out);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*
   Create Model from data file.  Assumes that the state variables have a flat
   structure.
*/
/**
TStateModel* CreateStateModel_VAR_File(FILE *f, char *filename)
{
  TMarkovStateVariable *sv;
  T_VAR_Parameters *p;

  //=== Create Markov State Variable ===
  sv=CreateMarkovStateVariable_File(f,filename,0);

  //=== Create VAR Parameters
  p=Create_VAR_Parameters_File(f,filename,sv);

  //=== Create TStateModel ===
  return CreateStateModel_new(sv,CreateRoutines_VAR(),p);
}
/**/
