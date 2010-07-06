
#include "switch.h"
#include "dw_array.h"
#include "dw_matrix_array.h"
#include "dw_error.h"
#include "dw_rand.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "modify_for_mex.h"

/*******************************************************************************/
/**************************** TMarkovStateVariable *****************************/
/*******************************************************************************/

/********************** TMarkovStateVariable Destructors ***********************/
/*
   Assumes:
     sv: valid pointer to TMarkovStateVariable structure or null pointer

   Results:
     Frees all memory allocated to sv.
*/
void FreeMarkovStateVariable(TMarkovStateVariable *sv)
{
  int i;
  if (sv)
    {
      dw_FreeArray(sv->S);

      FreeMatrix(sv->Q);

      FreeVector(sv->B);

/*        //====== Non-standard memory managment ======   ansi-c*/
      if (sv->b)
    {
      for (i=dw_DimA(sv->b)-1; i >= 0; i--)
        if (sv->b[i]) pElementV(sv->b[i])=(PRECISION*)NULL;
      dw_FreeArray(sv->b);
    }
/*        //===========================================   ansi-c*/

      FreeMatrix(sv->Prior);

      FreeVector(sv->Prior_B);

/*        //====== Non-standard memory managment ======   ansi-c*/
      if (sv->Prior_b)
    {
      for (i=dw_DimA(sv->Prior_b)-1; i >= 0; i--)
        if (sv->Prior_b[i]) pElementV(sv->Prior_b[i])=(PRECISION*)NULL;
      dw_FreeArray(sv->Prior_b);
    }
/*        //===========================================   ansi-c*/

      dw_FreeArray(sv->FreeDim);

      dw_FreeArray(sv->NonZeroIndex);

      FreeMatrix(sv->MQ);

      if (sv->n_state_variables > 1)
    dw_FreeArray(sv->state_variable);

      dw_FreeArray(sv->Index);

      dw_FreeArray(sv->lag_index);

      if (sv->SA)
    {
      for (i=dw_DimA(sv->SA)-1; i >= 0; i--) sv->SA[i]=(int*)NULL;
      dw_FreeArray(sv->SA);
    }

      if (sv->QA)
    {
      for (i=dw_DimA(sv->QA)-1; i >= 0; i--) sv->QA[i]=(TMatrix)NULL;
      dw_FreeArray(sv->QA);
    }

      if (sv->ba)
    {
      for (i=dw_DimA(sv->ba)-1; i >= 0; i--) sv->ba[i]=(TVector)NULL;
      dw_FreeArray(sv->ba);
    }

      if (sv->Prior_ba)
    {
      for (i=dw_DimA(sv->Prior_ba)-1; i >= 0; i--) sv->Prior_ba[i]=(TVector)NULL;
      dw_FreeArray(sv->Prior_ba);
    }

      swzFree(sv);
    }
}
/*******************************************************************************/

/********************** TMarkovStateVariable Constructors **********************/
/*
   Assumes
    nstates      : positive integer
    nobs         : positive integer
    Prior        : nstates x nstates matrix
    FreeDim      : integer array
    NonZeroIndex : nstates x nstates integer matrix
    MQ           : nstates x nstates matrix

   Returns
    A valid pointer to a single TMarkovStateVariable structure.

   Notes
    This is the basic constructor for the TMarkovStateVariable structure.  Upon
    error, this procedure terminates the program.
*/
TMarkovStateVariable* CreateMarkovStateVariable_Single(int nstates, int nobs, TMatrix Prior, int* FreeDim, int** NonZeroIndex, TMatrix MQ)
{
  TMarkovStateVariable *sv;
  int i, j, k, q, total_free=0, terminal_errors;

  if (!CheckRestrictions(FreeDim,NonZeroIndex,MQ,nstates))
    {
      swz_fprintf_err("CreateMarkovStateVariable_Single(): Error in restrictions\n");
      swzExit(0);
    }

  if (!CheckPrior(Prior,nstates) || !CheckPriorOnFreeParameters(Prior,NonZeroIndex,nstates))
    {
      swz_fprintf_err("CreateMarkovStateVariable_Single(): Error in priors\n");
      swzExit(0);
    }

/*    //=== Compute total number of free transition matrix parameters   ansi-c*/
  for (k=dw_DimA(FreeDim)-1; k >= 0; k--) total_free+=FreeDim[k];

  if ((nstates <= 0) || (nobs <= 0))
    {
      swz_fprintf_err("CreateMarkovStateVariable(): improper argument values\n");
      swzExit(0);
    }

  if (!(sv=(TMarkovStateVariable*)swzMalloc(sizeof(TMarkovStateVariable))))
    {
      swz_fprintf_err("CreateMarkovStateVariable(): out of memory\n");
      swzExit(0);
    }

/*    //=== Set flags ===   ansi-c*/
  sv->valid_transition_matrix=0;

/*    //=== Set to terminate on memory error ===   ansi-c*/
  terminal_errors=dw_SetTerminalErrors(MEM_ERR);

/*    //=== Sizes ===   ansi-c*/
  sv->nstates=nstates;
  sv->nobs=nobs;

/*    //== State vector ===   ansi-c*/
  dw_InitializeArray_int(sv->S=dw_CreateArray_int(nobs+1),0);

/*    //=== Number of lagged values of base state variable to encode ===   ansi-c*/
  sv->nlags_encoded=0;
  sv->nbasestates=nstates;
  sv->lag_index=CreateLagIndex(sv->nbasestates,sv->nlags_encoded,sv->nstates);

/*    //=== Transition matrix ===   ansi-c*/
  sv->Q=CreateMatrix(nstates,nstates);

/*    //=== Free transition matrix parameters ===   ansi-c*/
  sv->B=CreateVector(total_free);

  sv->b=dw_CreateArray_vector(dw_DimA(FreeDim));
  for (q=k=0; k < dw_DimA(FreeDim); k++)
    {
      sv->b[k]=CreateVector(FreeDim[k]);
/*        // seting up non-standard memory management   ansi-c*/
      swzFree(pElementV(sv->b[k]));
      pElementV(sv->b[k])=pElementV(sv->B)+q;
      q+=FreeDim[k];
    }

/*    //=== Prior information ===   ansi-c*/
  sv->Prior=EquateMatrix((TMatrix)NULL,Prior);

  sv->Prior_B=CreateVector(total_free);
  InitializeVector(sv->Prior_B,1.0);
  for (j=nstates-1; j >= 0; j--)
    for (i=nstates-1; i >= 0; i--)
      if ((k=NonZeroIndex[i][j]) >= 0)
    ElementV(sv->Prior_B,k)+=ElementM(Prior,i,j)-1.0;

  sv->Prior_b=dw_CreateArray_vector(dw_DimA(FreeDim));
  for (q=k=0; k < dw_DimA(FreeDim); k++)
    {
      sv->Prior_b[k]=CreateVector(FreeDim[k]);
/*        // seting up non-standard memory management   ansi-c*/
      swzFree(pElementV(sv->Prior_b[k]));
      pElementV(sv->Prior_b[k])=pElementV(sv->Prior_B)+q;
      q+=FreeDim[k];
    }

/*    //=== Restrictions ===   ansi-c*/
  sv->FreeDim=(int*)dw_CopyArray(NULL,FreeDim);
  sv->NonZeroIndex=(int**)dw_CopyArray(NULL,NonZeroIndex);
  sv->MQ=EquateMatrix((TMatrix)NULL,MQ);

/*    //=== Multiple state variables ===   ansi-c*/
  sv->parent=sv;
  sv->n_state_variables=1;
  sv->state_variable=(TMarkovStateVariable**)dw_CreateArray_array(1);
  sv->state_variable[0]=sv;
  sv->Index=(int**)dw_CreateArray_array(sv->nstates);
  sv->SA=(int**)dw_CreateArray_array(sv->n_state_variables);
  sv->QA=dw_CreateArray_matrix(sv->n_state_variables);
  sv->ba=dw_CreateArray_vector(dw_DimA(sv->b));
  for (k=dw_DimA(sv->ba)-1; k >= 0; k--) sv->ba[k]=sv->b[k];
  sv->Prior_ba=dw_CreateArray_vector(dw_DimA(sv->Prior_b));
  for (k=dw_DimA(sv->Prior_ba)-1; k >= 0; k--) sv->Prior_ba[k]=sv->Prior_b[k];

/*    //=== Initialize Index ===   ansi-c*/
  sv->Index[i=sv->nstates-1]=dw_CreateArray_int(sv->n_state_variables);
  for (k=sv->n_state_variables-1; k >= 0; k--)
    sv->Index[i][k]=sv->state_variable[k]->nstates-1;
  for (i--; i >= 0; i--)
    {
      sv->Index[i]=(int*)dw_CopyArray(NULL,sv->Index[i+1]);
      for (k=sv->n_state_variables-1; k >= 0; k--)
    if (--sv->Index[i][k] >= 0)
      break;
    else
      sv->Index[i][k]=sv->state_variable[k]->nstates-1;
    }

/*    //=== Initialize SA and QA ===   ansi-c*/
  for (k=sv->n_state_variables-1; k >= 0; k--)
    {
      sv->SA[k]=sv->state_variable[k]->S;
      sv->QA[k]=sv->state_variable[k]->Q;
    }

/*    //=== Control variables ===   ansi-c*/
  sv->UseErgodic=0;

/*    //=== Set Constants ===   ansi-c*/
  SetLogPriorConstant_SV(sv);

/*    //=== Set transition matrix to mean of prior ===   ansi-c*/
  SetTransitionMatrixToPriorMean_SV(sv);

/*    //=== Reset terminal errors ===   ansi-c*/
  dw_SetTerminalErrors(terminal_errors);

  return sv;
}

TMarkovStateVariable* CreateMarkovStateVariable_Multiple(int nobs, int n_state_variables, TMarkovStateVariable **state_variable)
{
  int i, j, k, terminal_errors;
  TMarkovStateVariable *sv;

  if ((n_state_variables <= 1) || (nobs <= 0) || !state_variable || (dw_DimA(state_variable) != n_state_variables))
    {
      printf("CreateMarkovStateVariable_Multiple(): invalid arguments\n");
      swzExit(0);
    }

  if (!(sv=(TMarkovStateVariable*)swzMalloc(sizeof(TMarkovStateVariable))))
    {
      printf("CreateMarkovStateVariable_Multiple(): out of memory\n");
      swzExit(0);
    }

/*    //=== Set to terminate on memory error ===   ansi-c*/
  terminal_errors=dw_SetTerminalErrors(MEM_ERR);

/*    //=== Flags ===   ansi-c*/
  sv->valid_transition_matrix=0;

/*    //=== Sizes ===   ansi-c*/
  for (sv->nstates=1, k=n_state_variables-1; k >= 0; k--) sv->nstates*=state_variable[k]->nstates;
  sv->nobs=nobs;

/*    //== State vector ===   ansi-c*/
  dw_InitializeArray_int(sv->S=dw_CreateArray_int(nobs+1),0);

/*    //=== Transition matrix ===   ansi-c*/
  sv->Q=CreateMatrix(sv->nstates,sv->nstates);

/*   //=== Free transition matrix parameters ===   ansi-c*/
  sv->B=(TVector)NULL;
  sv->b=(TVector*)NULL;

/*    //=== Number of lagged values of base state variable to encode ===   ansi-c*/
  sv->nlags_encoded=0;
  sv->nbasestates=sv->nstates;
  sv->lag_index=CreateLagIndex(sv->nbasestates,sv->nlags_encoded,sv->nstates);

/*    //=== Prior information ===   ansi-c*/
  sv->Prior=(TMatrix)NULL;
  sv->Prior_B=(TVector)NULL;
  sv->Prior_b=(TVector*)NULL;

/*    //=== Restrictions ===   ansi-c*/
  sv->FreeDim=(int*)NULL;
  sv->NonZeroIndex=(int**)NULL;
  sv->MQ=(TMatrix)NULL;

/*    //=== Multiple state variables ===   ansi-c*/
  sv->parent=sv;
  sv->n_state_variables=n_state_variables;
  sv->state_variable=state_variable;
  for (k=n_state_variables-1; k >= 0; k--) state_variable[k]->parent=sv;

/*    //=== Initialize Index ===   ansi-c*/
  sv->Index=(int**)dw_CreateArray_array(sv->nstates);
  sv->Index[i=sv->nstates-1]=dw_CreateArray_int(n_state_variables);
  for (k=n_state_variables-1; k >= 0; k--)
    sv->Index[i][k]=state_variable[k]->nstates-1;
  for (i--; i >= 0; i--)
    {
      sv->Index[i]=dw_CopyArray((int*)NULL,sv->Index[i+1]);
      for (k=n_state_variables-1; k >= 0; k--)
    if (--(sv->Index[i][k]) >= 0)
      break;
    else
      sv->Index[i][k]=state_variable[k]->nstates-1;
    }

/*    //=== Initialize SA and QA ===   ansi-c*/
  sv->SA=(int**)dw_CreateArray_array(n_state_variables);
  sv->QA=dw_CreateArray_matrix(n_state_variables);
  for (k=0; k < n_state_variables; k++)
    {
      sv->SA[k]=state_variable[k]->S;
      sv->QA[k]=state_variable[k]->Q;
    }

/*    //=== Initialize ba and Prior_ba ===   ansi-c*/
  for (i=k=0; k < n_state_variables; k++) i+=dw_DimA(state_variable[k]->ba);
  sv->ba=dw_CreateArray_vector(i);
  for (i=k=0; k < n_state_variables; k++)
    for (j=0; j < dw_DimA(state_variable[k]->ba); i++, j++)
      sv->ba[i]=state_variable[k]->ba[j];

  for (i=k=0; k < n_state_variables; k++) i+=dw_DimA(state_variable[k]->Prior_ba);
  sv->Prior_ba=dw_CreateArray_vector(i);
  for (i=k=0; k < n_state_variables; k++)
    for (j=0; j < dw_DimA(state_variable[k]->Prior_ba); i++, j++)
      sv->Prior_ba[i]=state_variable[k]->Prior_ba[j];

/*    //=== Control variables ===   ansi-c*/
  sv->UseErgodic=0;

/*    //=== Set Constants ===   ansi-c*/
  SetLogPriorConstant_SV(sv);

/*    //=== Set transition matrix to mean of prior ===   ansi-c*/
  SetTransitionMatrixToPriorMean_SV(sv);

/*    //=== Reset terminal errors ===   ansi-c*/
  dw_SetTerminalErrors(terminal_errors);

  return sv;
}

TMarkovStateVariable* CreateMarkovStateVariable_Mixture(int nstates, int nobs, TMatrix Prior)
{
  int i, j;
  TMarkovStateVariable *sv;
  int* FreeDim;
  int** NonZeroIndex;
  TMatrix MQ;
  NonZeroIndex=dw_CreateRectangularArray_int(nstates,nstates);
  for (i=nstates-1; i >= 0; i--)
    for (j=nstates-1; j >= 0; j--)
      NonZeroIndex[i][j]=i;
  MQ=InitializeMatrix(CreateMatrix(nstates,nstates),1.0);
  FreeDim=dw_CreateArray_int(1);
  FreeDim[0]=nstates;
  sv=CreateMarkovStateVariable_Single(nstates,nobs,Prior,FreeDim,NonZeroIndex,MQ);
  dw_FreeArray(NonZeroIndex);
  FreeMatrix(MQ);
  dw_FreeArray(FreeDim);
  return sv;
}

TMarkovStateVariable* CreateMarkovStateVariable_NoRestrictions(int nstates, int nobs, TMatrix Prior)
{
  int i, j;
  TMarkovStateVariable *sv;
  int* FreeDim;
  int** NonZeroIndex;
  TMatrix MQ;
  NonZeroIndex=dw_CreateRectangularArray_int(nstates,nstates);
  for (i=nstates-1; i >= 0; i--)
    for (j=nstates-1; j >= 0; j--)
      NonZeroIndex[i][j]=i+nstates*j;
  InitializeMatrix(MQ=CreateMatrix(nstates,nstates),1.0);
  FreeDim=dw_CreateArray_int(nstates);
  for (i=nstates-1; i >= 0; i--) FreeDim[i]=nstates;
  sv=CreateMarkovStateVariable_Single(nstates,nobs,Prior,FreeDim,NonZeroIndex,MQ);
  dw_FreeArray(NonZeroIndex);
  FreeMatrix(MQ);
  dw_FreeArray(FreeDim);
  return sv;
}

TMarkovStateVariable* CreateMarkovStateVariable_Exclusion(int nstates, int nobs, TMatrix Prior, TMatrix Exclusion)
{
  int i, j, k, q;
  TMarkovStateVariable *sv;
  int* FreeDim;
  int** NonZeroIndex;
  TMatrix MQ;
  dw_InitializeArray_int(NonZeroIndex=dw_CreateRectangularArray_int(nstates,nstates),-1);
  MQ=InitializeMatrix(CreateMatrix(nstates,nstates),0.0);
  FreeDim=dw_CreateArray_int(nstates);
  for (k=j=0; j < nstates; j++)
    {
      for (q=i=0; i < nstates; i++)
    if (ElementM(Exclusion,i,j) > 0)
      {
        NonZeroIndex[i][j]=k++;
        ElementM(MQ,i,j)=1.0;
        q++;
      }
      FreeDim[j]=q;
    }
  sv=CreateMarkovStateVariable_Single(nstates,nobs,Prior,FreeDim,NonZeroIndex,MQ);
  dw_FreeArray(NonZeroIndex);
  FreeMatrix(MQ);
  dw_FreeArray(FreeDim);
  return sv;
}

TMarkovStateVariable* CreateMarkovStateVariable_SimpleRestrictions(int nstates, int nobs, TMatrix Prior, TMatrix* Restriction)
{
  int i, j, k, q;
  TMarkovStateVariable *sv;
  int* FreeDim;
  int** NonZeroIndex;
  TMatrix MQ;
  dw_InitializeArray_int(NonZeroIndex=dw_CreateRectangularArray_int(nstates,nstates),-1);
  InitializeMatrix(MQ=CreateMatrix(nstates,nstates),0.0);
  FreeDim=dw_CreateArray_int(nstates);
  for (q=k=0; k < nstates; k++)
    {
      FreeDim[k]=ColM(Restriction[k]);
      for (i=0; i < nstates; i++)
        {
          for (j=0; j < ColM(Restriction[k]); j++)
            if (ElementM(Restriction[k],i,j) > 0)
              {
                NonZeroIndex[i][k]=q+j;
                ElementM(MQ,i,k)=ElementM(Restriction[k],i,j);
                break;
              }
        }
      q+=FreeDim[k];
    }
  sv=CreateMarkovStateVariable_Single(nstates,nobs,Prior,FreeDim,NonZeroIndex,MQ);
  dw_FreeArray(NonZeroIndex);
  FreeMatrix(MQ);
  dw_FreeArray(FreeDim);
  return sv;
}

TMarkovStateVariable* CreateMarkovStateVariable_ConstantState(int nobs)
{
  TMarkovStateVariable *sv;
  int* FreeDim;
  int** NonZeroIndex;
  TMatrix MQ, Prior;
  dw_InitializeArray_int(NonZeroIndex=dw_CreateRectangularArray_int(1,1),0);
  InitializeMatrix(MQ=CreateMatrix(1,1),1.0);
  dw_InitializeArray_int(FreeDim=dw_CreateArray_int(1),1);
  InitializeMatrix(Prior=CreateMatrix(1,1),1.0);
  sv=CreateMarkovStateVariable_Single(1,nobs,Prior,FreeDim,NonZeroIndex,MQ);
  dw_FreeArray(NonZeroIndex);
  FreeMatrix(MQ);
  dw_FreeArray(FreeDim);
  FreeMatrix(Prior);
  return sv;
}

TMarkovStateVariable* DuplicateMarkovStateVariable(TMarkovStateVariable *sv)
{
  TMarkovStateVariable **sv_array, *dup;
  int k;
  if (sv->n_state_variables > 1)
    {
      sv_array=dw_CreateArray_pointer(sv->n_state_variables,NULL);
      for (k=sv->n_state_variables-1; k >= 0; k--)
    sv_array[k]=DuplicateMarkovStateVariable(sv->state_variable[k]);
      dup=CreateMarkovStateVariable_Multiple(sv->nobs,sv->n_state_variables,sv_array);
      EquateMatrix(dup->Q,sv->Q);
      memcpy(dup->S,sv->S,(sv->nobs+1)*sizeof(int));
      dup->valid_transition_matrix=sv->valid_transition_matrix;
    }
  else
    {
      dup=CreateMarkovStateVariable_Single(sv->nstates,sv->nobs,sv->Prior,sv->FreeDim,sv->NonZeroIndex,sv->MQ);
      EquateMatrix(dup->Q,sv->Q);
      EquateVector(dup->B,sv->B);
      memcpy(dup->S,sv->S,(sv->nobs+1)*sizeof(int));
      dup->valid_transition_matrix=sv->valid_transition_matrix;
    }
  return dup;
}

/*

*/
TMarkovStateVariable* RestrictMarkovStateVariable(TMarkovStateVariable *sv, int nstates)
{
  TMarkovStateVariable *rsv;
  int* free_translation;
  int** NonZeroIndex;
  int *FreeDim;
  TMatrix MQ, Prior;
  int i, j, k, m, n, q;
  PRECISION sum_in, sum_out, scale;

  if (nstates == 1) return CreateMarkovStateVariable_ConstantState(sv->nobs);

  if (nstates == sv->nstates) return DuplicateMarkovStateVariable(sv);

  free_translation=(int*)swzMalloc(DimV(sv->B)*sizeof(int));

/*    // free_translation[i] = 1 if B[i] is accessed   ansi-c*/
  for (i=DimV(sv->B)-1; i >= 0; i--) free_translation[i]=0;
  for (i=nstates-1; i >= 0; i--)
    for (j=nstates-1; j >= 0; j--)
      if (sv->NonZeroIndex[i][j] >= 0)
    free_translation[sv->NonZeroIndex[i][j]]=1;

/*    // free_translation[i] will be index into new quasi-free Dirichlet vector B   ansi-c*/
  for (i=j=0; j < DimV(sv->B); j++)
    free_translation[j]=free_translation[j] ? i++ : -1;

/*    // Get number of new quasi-free variables   ansi-c*/
  for (k=m=n=0; k < dw_DimA(sv->FreeDim); k++)
    for (q=n, n+=sv->FreeDim[k]; q < n; q++)
      if (free_translation[q] != -1)
    {
      m++;
      break;
    }
  dw_InitializeArray_int(FreeDim=dw_CreateArray_int(m),0);

/*    // Compute size of new quasi-free variables   ansi-c*/
  for (k=m=n=0; k < dw_DimA(sv->FreeDim); k++)
    for (q=n, n+=sv->FreeDim[k]; q < n; q++)
      if (free_translation[q] != -1)
    {
      for (FreeDim[m]=1, q++; q < n; q++)
        if (free_translation[q] != -1) FreeDim[m]++;
      m++;
      break;
    }

/*    // Create new NonZeroIndex   ansi-c*/
  NonZeroIndex=dw_CreateRectangularArray_int(nstates,nstates);
  for (i=nstates-1; i >= 0; i--)
    for (j=nstates-1; j >= 0; j--)
      NonZeroIndex[i][j]=free_translation[sv->NonZeroIndex[i][j]];

/*    // Create new MQ   ansi-c*/
  InitializeMatrix(MQ=CreateMatrix(nstates,nstates),0.0);
  for (j=0; j < nstates; j++)
    for (q=DimV(sv->B)-1; q >= 0; q--)
      if (free_translation[q] != -1)
    {
      sum_in=sum_out=0.0;
      for (i=0; i < nstates; i++)
        if (sv->NonZeroIndex[i][j] == q) sum_in+=ElementM(sv->MQ,i,j);
      for ( ; i < sv->nstates; i++)
        if (sv->NonZeroIndex[i][j] == q) sum_out+=ElementM(sv->MQ,i,j);
      scale=(sum_in > 0) ? (sum_in + sum_out)/sum_in : 1.0;
      for (i=0; i < nstates; i++)
        if (sv->NonZeroIndex[i][j] == q) ElementM(MQ,i,j)=scale*ElementM(sv->MQ,i,j);
    }

/*    // Prior   ansi-c*/
/*    // The new prior is chosen so that Prior[i][j] = scale*sv->Prior[i][j] where   ansi-c*/
/*    // scale = 1 if sv->Prior[i][j] <= 1 and scale = (nstates-1)/(sv->nstates-1)   ansi-c*/
/*    // otherwise.  This tends to keep the old and new prior mean roughly equal for   ansi-c*/
/*    // those elements with larger prior means and scales the prior mean equally   ansi-c*/
/*    // for those with smaller prior means.  This ie exactly true if sv->Prior[i][j]   ansi-c*/
/*    // were equal to one for all j except one.   ansi-c*/
  Prior=CreateMatrix(nstates,nstates);
  for (i=nstates-1; i >= 0; i--)
    for (j=nstates-1; j >= 0; j--)
      ElementM(Prior,i,j)=(ElementM(sv->Prior,i,j) > 1) ? (nstates-1)*ElementM(sv->Prior,i,j)/(sv->nstates-1) : ElementM(sv->Prior,i,j);

/*    // Attempt to make new Markov state variable.  MQ/NonZeroIndex may not be valid   ansi-c*/
  rsv=CreateMarkovStateVariable_Single(nstates,sv->nobs,Prior,FreeDim,NonZeroIndex,MQ);

/*    // Clean up   ansi-c*/
  FreeMatrix(Prior);
  FreeMatrix(MQ);
  dw_FreeArray(NonZeroIndex);
  dw_FreeArray(FreeDim);
  swzFree(free_translation);

  return rsv;
}

/*

*/
TMatrix ConvertBaseTransitionMatrix(TMatrix Q, TMatrix bQ, int nlags_encoded)
{
  int n, a, b, c, d, e, i, j, k, m;
  PRECISION p;

  if (nlags_encoded == 0)
    return EquateMatrix(Q,bQ);

  if (!bQ)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  n=RowM(bQ);
  if (n != ColM(bQ))
    {
      dw_Error(SIZE_ERR);
      return (TMatrix)NULL;
    }
  for (a=n, i=nlags_encoded; i > 0; i--) a*=n;
  if (!Q)
    {
      if (!(Q=CreateMatrix(a,a)))
    return (TMatrix)NULL;
    }
  else
    if ((a != RowM(Q)) || (a != ColM(Q)))
      {
    dw_Error(SIZE_ERR);
    return (TMatrix)NULL;
      }

  a/=n;
  b=a/n;
  InitializeMatrix(Q,0.0);
  for (j=0; j < n; j++)
    for (i=0; i < n; i++)
      for (p=ElementM(bQ,i,j), k=0; k < b; k++)
    for (c=b*j+k, d=i*a+c, e=n*c, m=0; m < n; m++)
      ElementM(Q,d,e+m)=p;

  return Q;
}

/*


*/
int** CreateLagIndex(int nbasestates, int nlags, int nstates)
{
  int **lag_index;
  int j, k;

  lag_index=dw_CreateRectangularArray_int(nstates,nlags+1);
  for (j=nlags; j >= 0; j--) lag_index[0][j]=0;
  for (k=1; k < nstates; k++)
    {
      for (j=nlags; j >= 0; j--)
    if (lag_index[k-1][j] < nbasestates-1)
      {
        lag_index[k][j]=lag_index[k-1][j]+1;
        break;
      }
    else
      lag_index[k][j]=0;
      for (--j; j >= 0; j--) lag_index[k][j]=lag_index[k-1][j];
    }

  return lag_index;
}

/*

*/
TMarkovStateVariable* CreateMarkovStateVariable_Lags(int nlags, TMarkovStateVariable *base)
{
  TMarkovStateVariable *sv;
  int** NonZeroIndex;
  TMatrix MQ;
  TMatrix Prior;
  int nstates, i, j, k, m, n, p, q;
  PRECISION scale;

  if (base->n_state_variables > 1)
    {
      swz_fprintf_err("CreateMarkovStateVariable_Lags():  multiple state variable for base.");
      swzExit(0);
    }

  if (nlags > 0)
    {
      nstates=base->nstates;
      for (i=nlags; i > 0; i--) nstates*=base->nstates;
      dw_InitializeArray_int(NonZeroIndex=dw_CreateRectangularArray_int(nstates,nstates),-1);
      InitializeMatrix(MQ=CreateMatrix(nstates,nstates),0.0);
      InitializeMatrix(Prior=CreateMatrix(nstates,nstates),1.0);

      for (m=base->nstates, n=1, i=nlags-1; i > 0; i--) n*=m;
      scale=pow(base->nstates,-nlags);
      for (p=i=0; i < base->nstates; i++)
    for (q=j=0; j < base->nstates; j++)
      for (k=0; k < n; p++, k++)
        for (m=0; m < base->nstates; q++, m++)
          {
        NonZeroIndex[p][q]=base->NonZeroIndex[i][j];
        ElementM(MQ,p,q)=ElementM(base->MQ,i,j);
        ElementM(Prior,p,q)=scale*(ElementM(base->Prior,i,j)-1.0)+1.0;
          }

      sv=CreateMarkovStateVariable_Single(nstates,base->nobs,Prior,base->FreeDim,NonZeroIndex,MQ);

      dw_FreeArray(sv->lag_index);
      sv->nlags_encoded=nlags;
      sv->nbasestates=base->nstates;
      sv->lag_index=CreateLagIndex(sv->nbasestates,sv->nlags_encoded,sv->nstates);

      dw_FreeArray(NonZeroIndex);
      FreeMatrix(MQ);
      FreeMatrix(Prior);
    }
  else
    sv=CreateMarkovStateVariable_Single(base->nstates,base->nobs,base->Prior,base->FreeDim,base->NonZeroIndex,base->MQ);
  return sv;
}
/*******************************************************************************/

/**************** TMarkovStateVariable Data Extraction Routines ****************/
/*
   Get the base transition matrix, which is the transition matrix with any
   encoded lags removed.
*/
TMatrix GetBaseTransitionMatrix_SV(TMatrix Q, TMarkovStateVariable *sv)
{
  int i, j, k, n;
  TMatrix *A;

  if (!sv->valid_transition_matrix) return (TMatrix)NULL;

  if (sv->n_state_variables > 1)
    {
      A=dw_CreateArray_matrix(sv->n_state_variables);
      for (i=sv->n_state_variables-1; i >= 0; i--)
    A[i]=GetBaseTransitionMatrix_SV((TMatrix)NULL,sv->state_variable[i]);
      Q=MatrixTensor(Q,A);
      dw_FreeArray(A);
      return Q;
    }
  else
    if (sv->nlags_encoded > 0)
      {
    n=sv->nbasestates;
    for (k=1, i=sv->nlags_encoded-1; i > 0; i--) k=k*n;
    if (!Q)
      {
        if (!(Q=CreateMatrix(n,n)))
          return (TMatrix)NULL;
      }
    else
      if ((RowM(Q) != n) || (ColM(Q) != n))
        {
          dw_Error(SIZE_ERR);
          return (TMatrix)NULL;
        }
    for (i=n-1; i >= 0; i--)
      for (j=n-1; j >= 0; j--)
        ElementM(Q,i,j)=ElementM(sv->Q,(i*n+j)*k,j*n*k);
    return Q;
      }
    else
      return EquateMatrix(Q,sv->Q);
}
/*******************************************************************************/

/***************** TMarkovStateVariable Normalization Routines *****************/
void PropagateSwap_SV(TMarkovStateVariable *sv)
{
  if (sv->parent != sv)
    {
      sv=sv->parent;
      sv->Q=MatrixTensor(sv->Q,sv->QA);
      PropagateSwap_SV(sv);
    }
}

void Swap_SV(TMarkovStateVariable *sv, int i, int j)
{
  TMatrix Q, Y;
  TPermutation X;
  if (sv->n_state_variables > 1)
    {
      swz_fprintf_err("Swap(): Can only swap indices for terminal state variables.\n");
      swzExit(0);
    }

  if ((i < 0) || (j < 0) || (i >= sv->nbasestates) || (j >= sv->nbasestates))
    {
      swz_fprintf_err("Swap(): Indicies out of range.\n");
      swzExit(0);
    }

  X=TranspositionPermutation((TPermutation)NULL,i,j,sv->nbasestates);
  if (sv->nlags_encoded)
    {
      Q=GetBaseTransitionMatrix_SV((TMatrix)NULL,sv);
      Y=ProductMP((TMatrix)NULL,Q,X);
      TransposeProductPM(Q,X,Y);
      ConvertBaseTransitionMatrix(sv->Q,Q,sv->nlags_encoded);
      FreeMatrix(Q);
    }
  else
    {
      Y=ProductMP((TMatrix)NULL,sv->Q,X);
      TransposeProductPM(sv->Q,X,Y);
    }
  FreeMatrix(Y);
  FreePermutation(X);
  if (!Update_B_from_Q_SV(sv))
    {
      swz_fprintf_err("Swap(): Restrictions violated.\n");
      swzExit(0);
    }
  PropagateSwap_SV(sv);
}

/*******************************************************************************/

/********************* TMarkovStateVariable Prior Routines *********************/
void SetLogPriorConstant_SV(TMarkovStateVariable *sv)
{
  PRECISION sum, x;
  int i, j;
  sv->LogPriorConstant=0.0;
  for (i=dw_DimA(sv->ba)-1; i >= 0; i--)
    {
      for (sum=0.0, j=DimV(sv->ba[i])-1; j >= 0; j--)
    {
      sum+=(x=ElementV(sv->Prior_ba[i],j));
      sv->LogPriorConstant-=dw_log_gamma(x);
    }
      sv->LogPriorConstant+=dw_log_gamma(sum);
    }
}

PRECISION LogPrior_SV(TMarkovStateVariable *sv)
{
  PRECISION log_prior=sv->LogPriorConstant, y;
  int i, j;
  for (i=dw_DimA(sv->ba)-1; i >= 0; i--)
    for (j=DimV(sv->ba[i])-1; j >= 0; j--)
      if ((y=ElementV(sv->ba[i],j)) > 0.0)
    log_prior+=(ElementV(sv->Prior_ba[i],j)-1.0)*log(y);
      else
    if (ElementV(sv->Prior_ba[i],j) > 1.0)
      return MINUS_INFINITY;
    else
      if (ElementV(sv->Prior_ba[i],j) < 1.0)
        return PLUS_INFINITY;
  return log_prior;
}
/******************************************************************************/

/********************** TMarkovStateVariables Simulation **********************/
/*
   Assumes:
     sv is a pointer to a properly initialized TMarkovStateVariable and all
     state vectors have been computed.

   Results:
    Draws the transition matrix sv->P from the posterior distribution,
    conditional on sv->S and the parameters other than sv->P.
*/
void DrawTransitionMatrix_SV(TMarkovStateVariable *sv)
{
  int i, j, k, t;

  if (sv->n_state_variables > 1)
    {
      for (k=sv->n_state_variables-1; k >= 0; k--)
    DrawTransitionMatrix_SV(sv->state_variable[k]);
      MatrixTensor(sv->Q,sv->QA);

/*        // Flags   ansi-c*/
      sv->valid_transition_matrix=1;
    }
  else
    {
/*        // Set Dirichlet parameters   ansi-c*/
      EquateVector(sv->B,sv->Prior_B);
      for (t=sv->nobs; t > 0; t--)
    if ((k=sv->NonZeroIndex[sv->S[t]][sv->S[t-1]]) >= 0)
      ElementV(sv->B,k)+=1.0;

/*        // Generate b[j]   ansi-c*/
      for (k=dw_DimA(sv->FreeDim)-1; k >= 0; k--)
    if (!DrawDirichletVector(sv->b[k],sv->b[k]))
      {
        swz_fprintf_err("Error drawing Dirichlet vector\n");
        swzExit(0);
      }

/*        // Compute Q   ansi-c*/
      for (j=sv->nstates-1; j >= 0; j--)
    for (i=sv->nstates-1; i >= 0; i--)
      ElementM(sv->Q,i,j)=((k=sv->NonZeroIndex[i][j]) >= 0) ? ElementM(sv->MQ,i,j)*ElementV(sv->B,k) : 0.0;

/*        // Flags   ansi-c*/
      sv->valid_transition_matrix=1;
    }
}

/*
   Assumes:
    sv is a pointer to a properly initialized TMarkovStateVariable.

   Results:
    Draws the transition matrix sv->P from the prior distribution.
*/
void DrawTransitionMatrixFromPrior_SV(TMarkovStateVariable *sv)
{
  int i, j, k;

  if (sv->n_state_variables > 1)
    {
      for (k=sv->n_state_variables-1; k >= 0; k--)
    DrawTransitionMatrixFromPrior_SV(sv->state_variable[k]);
      MatrixTensor(sv->Q,sv->QA);

/*        // Flags   ansi-c*/
      sv->valid_transition_matrix=1;
    }
  else
    {
/*        // Set Dirichlet parameters   ansi-c*/
      EquateVector(sv->B,sv->Prior_B);

/*        // Generate b[j]   ansi-c*/
      for (j=dw_DimA(sv->FreeDim)-1; j >= 0; j--)
    if (!DrawDirichletVector(sv->b[j],sv->b[j]))
      {
        swz_fprintf_err("Error drawing Dirichlet vector\n");
        swzExit(0);
      }

/*        // Compute P   ansi-c*/
      for (j=sv->nstates-1; j >= 0; j--)
    for (i=sv->nstates-1; i >= 0; i--)
      ElementM(sv->Q,i,j)=((k=sv->NonZeroIndex[i][j]) >= 0) ? ElementM(sv->MQ,i,j)*ElementV(sv->B,k) : 0.0;

/*        // Flags   ansi-c*/
      sv->valid_transition_matrix=1;
    }
}

void SetTransitionMatrixToPriorMean_SV(TMarkovStateVariable *sv)
{
  int i, j, k;
  PRECISION sum;
  if (sv->n_state_variables > 1)
    {
      for (k=sv->n_state_variables-1; k >= 0; k--)
    SetTransitionMatrixToPriorMean_SV(sv->state_variable[k]);
      MatrixTensor(sv->Q,sv->QA);
      sv->valid_transition_matrix=1;
    }
  else
    {
      for (j=dw_DimA(sv->Prior_b)-1; j >= 0; j--)
    {
      for (sum=ElementV(sv->Prior_b[j],0), i=DimV(sv->Prior_b[j])-1; i > 0; i--)
        sum+=ElementV(sv->Prior_b[j],i);
      for (i=DimV(sv->Prior_b[j])-1; i >= 0; i--)
        ElementV(sv->b[j],i)=ElementV(sv->Prior_b[j],i)/sum;
    }
      Update_Q_from_B_SV(sv);
      sv->valid_transition_matrix=1;
    }
}

/*
   Assumes:
     sv is a pointer to a properly initialized TMarkovStateVariable structure and
     the transition matrix has been.

   Results:
     Draws the vector of states S from the distribution defined by the transition
     matrix.

   Notes:
     If sv->Q was drawn from the prior distribution, then the vector of states
     can be considered as being drawing from the prior distribution.
*/
void DrawStatesFromTransitionMatrix_SV(TMarkovStateVariable *sv)
{
  int i, j, t;
  PRECISION u, s;

  if (!sv->valid_transition_matrix)
    {
      printf("DrawStatesFromTransitionMatrix_SV():  Invalid transition matrix.\n");
      swzExit(0);
    }

/*    //=== Draw initial state from ergodic or uniform distribution ===   ansi-c*/
  if (sv->UseErgodic)
    {
      TVector v;
      if (!(v=Ergodic((TVector)NULL,sv->Q)))
    {
      printf("DrawStatesFromTransitionMatrix_SV():  Ergodic distribution does not exists.\n");
      swzExit(0);
    }
      if ((u=dw_uniform_rnd()) >= (s=ElementV(v,i=sv->nstates-1)))
    while (--i > 0)
      if (u < (s+=ElementV(v,i))) break;
      FreeVector(v);
    }
  else
    if ((i=(int)floor(dw_uniform_rnd()*sv->nstates)) >= sv->nstates) i=sv->nstates-1;
  sv->S[0]=j=i;

/*    //=== Draw subsequent states from transition matrix ===   ansi-c*/
  for (t=1; t <= sv->nobs; t++)
    {
      if ((u=dw_uniform_rnd()) >= (s=ElementM(sv->Q,i=sv->nstates-1,j)))
    while (--i > 0)
      if (u < (s+=ElementM(sv->Q,i,j))) break;
      sv->S[t]=j=i;
    }

/*    //====== Propagate states ======   ansi-c*/
  PropagateStates_SV(sv);
}
/******************************************************************************/

/******************* TMarkovStateVariables Utility Routines *******************/
/*
   Set sv->valid_transition_matrix to zero for all state variables in sv.
*/
void InvalidateTransitionMatrices_SV(TMarkovStateVariable *sv)
{
  int k;
  if (sv->n_state_variables > 1)
    for (k=sv->n_state_variables-1; k >=0; k--)
      InvalidateTransitionMatrices_SV(sv->state_variable[k]);
  sv->valid_transition_matrix=0;
}

/*
   Set sv->valid_transition_matrix to one for all state variables in sv.
*/
void ValidateTransitionMatrices_SV(TMarkovStateVariable *sv)
{
  int k;
  if (sv->n_state_variables > 1)
    for (k=sv->n_state_variables-1; k >=0; k--)
      ValidateTransitionMatrices_SV(sv->state_variable[k]);
  sv->valid_transition_matrix=1;
}


/*
   Distributes the vector of states throughout the tree of states variables
*/
void PropagateStates_SV(TMarkovStateVariable *sv)
{
  int k, t;
  if (sv->n_state_variables > 1)
    {
      for (t=sv->nobs; t >= 0; t--)
    for (k=sv->n_state_variables-1; k >= 0; k--)
      sv->SA[k][t]=sv->Index[sv->S[t]][k];
      for (k=sv->n_state_variables-1; k >=0; k--)
    PropagateStates_SV(sv->state_variable[k]);
    }
}

/*
   If there are multiple state variables, computes the transition matrix as the
   tensor product of transition matrics.  The single Markov state variables
   must have a valid transition matrix.  Returns one upon success and zero upon
   failure.  Failure occurs if a single Markov state variable does not have a
   valid transition matrix.
*/
int PropagateTransitionMatrices_SV(TMarkovStateVariable *sv)
{
  int k;
  if (sv->n_state_variables > 1)
    {
      for (k=sv->n_state_variables-1; k >= 0; k--)
    if (!PropagateTransitionMatrices_SV(sv->state_variable[k]))
      return sv->valid_transition_matrix=0;
      sv->valid_transition_matrix=1;
      MatrixTensor(sv->Q,sv->QA);
    }
  return sv->valid_transition_matrix;
}

/*
   Assumes
     sv->Q has been updated and that sv->n_state_variables = 1.

   Returns
     One upon success and zero upon failure

   Results
     Updates sv->B (and hence sv->b and sv->ba) from sv->Q.  Then sv->Q is
     updated from the computed sv->B.

   Notes
     If sv->n_state_variables != 1, then the routine returns zero.

     If any of the elements of sv->b[i] is negative or their sum differs from one
     by more than DimV(sv->b[i])*MACHINE_EPSILON, then the routine returns zero.

     If any element of the updated sv->Q differs from the original sv->Q by more
     than SQRT_MACHINE_EPSILON, then the routine returns zero.

     Because of the checking done, care should be taken in calling this routine
     in instances where efficiency is important.
*/
int Update_B_from_Q_SV(TMarkovStateVariable *sv)
{
  int i, j, k;
  TVector scale;
  PRECISION x, sum;
  if (sv->n_state_variables > 1)
    return 0;

/*    // Compute B   ansi-c*/
  InitializeVector(scale=CreateVector(DimV(sv->B)),0.0);
  InitializeVector(sv->B,0.0);
  for (i=sv->nstates-1; i >= 0; i--)
    for (j=sv->nstates-1; j >= 0; j--)
      if ((k=sv->NonZeroIndex[i][j]) >= 0)
    {
      x=ElementM(sv->MQ,i,j);
      ElementV(sv->B,k)+=x*ElementM(sv->Q,i,j);
      ElementV(scale,k)+=x*x;
    }
  for (i=DimV(scale)-1; i >= 0; i--) ElementV(sv->B,i)/=ElementV(scale,i);
  FreeVector(scale);

/*    // Check computed B   ansi-c*/
  for (i=dw_DimA(sv->b)-1; i >= 0; i--)
    {
      for (sum=0.0, j=DimV(sv->b[i])-1; j >= 0; j--)
    if (ElementV(sv->b[i],j) < 0)
      return 0;
    else
      sum+=ElementV(sv->b[i],j);
      if (fabs(sum-1.0) > SQRT_MACHINE_EPSILON)
    return 0;
    }

/*    // Check computed Q   ansi-c*/
  for (i=sv->nstates-1; i >= 0; i--)
    for (j=sv->nstates-1; j >= 0; j--)
      {
    x=ElementM(sv->Q,i,j);
    if ((k=sv->NonZeroIndex[i][j]) >= 0)
      x-=ElementM(sv->MQ,i,j)*ElementV(sv->B,k);
    if (fabs(x) > SQRT_MACHINE_EPSILON)
      return 0;
    else
      ElementM(sv->Q,i,j)=(k >= 0) ? ElementM(sv->MQ,i,j)*ElementV(sv->B,k) : 0.0;
      }

/*    // Success   ansi-c*/
  return 1;
}

/*
   Assumes:
     sv->ba has been updated

   Results:
     Computes sv->Q

   Notes:
     Uses recursion to compute sv->Q
*/
void Update_Q_from_B_SV(TMarkovStateVariable *sv)
{
  int i, j, k;
  if (sv->n_state_variables > 1)
    {
      for (k=sv->n_state_variables-1; k >= 0; k--) Update_Q_from_B_SV(sv->state_variable[k]);
      MatrixTensor(sv->Q,sv->QA);
    }
  else
    {
      for (j=sv->nstates-1; j >= 0; j--)
    for (i=sv->nstates-1; i >= 0; i--)
      ElementM(sv->Q,i,j)=((k=sv->NonZeroIndex[i][j]) >= 0) ? ElementM(sv->MQ,i,j)*ElementV(sv->B,k) : 0.0;
    }
}

/*
    Assumes:
      sv - valid pointer to TMarkovStateVariable structure.

    Returns:
      total number of states
*/
int TotalNumberStateVariables_SV(TMarkovStateVariable *sv)
{
  int i, n;
  if (sv->n_state_variables == 1)
    return 1;
  else
    {
      for (n=0, i=sv->n_state_variables-1; i >= 0; i--)
    n+=TotalNumberStateVariables_SV(sv->state_variable[i]);
      return n;
    }
}

/******************************************************************************/
/** Translation Arrays for Collections of Independent Markov State Variables **/
/******************************************************************************/
/*
   Counts the total number of states in the array of Markov state variables
*
static int NumberStates(TMarkovStateVariable** sv, int n)
{
  int i, m=1;
  for (i=0; i < n; i++) m*=sv[i]->nstates;
  return m;
}

/*
   Assumes
     j:  An integer with 0 <= j < dw_DimA(states)
     states:  Rectangular array of integers
*/
int GetNumberStatesFromTranslationMatrix(int j, int **states)
{
  int i, n=0;
  for (i=dw_DimA(states[j])-1; i >= 0; i--)
    if (states[j][i] > n) n=states[j][i];
  return n+1;
}


/*
   Attempts to find the Markov state variable sv in the chain of Markov variables
   given by top.  If found, multiples the existing index by the number of states
   in sv and adds the appropriate index number.
*/
static int IncrementIndex(int idx, int i, TMarkovStateVariable *top, TMarkovStateVariable *sv)
{
  int j, k;
  if (top == sv) return i + idx * sv->nstates;
  if (top->n_state_variables > 1)
    for (j=0; j < top->n_state_variables; j++)
      if ((k=IncrementIndex(idx,top->Index[i][j],top->state_variable[j],sv)) >= 0)
    return k;
  return -1;
}

/*
   Assumes:
     sv    - Valid pointer to TMarkovStateVariable structure.
     list  - Array of valid pointers of TMarkovStateVariable structures.  Each
             structure in the array must be in the tree whose root is sv.

   Returns:
     Integer array of length sv->nstates.
*/
int* CreateStateIndex(TMarkovStateVariable* sv, TMarkovStateVariable** list, int n)
{
  int i, j, k;
  int* index;
  if (!(index=dw_CreateArray_int(sv->nstates)))
    {
      printf("CreateStateIndex():  Out of memory.\n");
      swzExit(0);
    }
  for (i=sv->nstates-1; i >= 0; i--)
    {
      for (k=j=0; j < n; j++)
    if ((k=IncrementIndex(k,i,sv,list[j])) == -1)
      {
        printf("CreateStateIndex():  Unable to find required state variable.\n");
        swzExit(0);
      }
     index[i]=k;
    }
  return index;
}

/*
   Assumes:
     list :
*/
int **CreateTranslationMatrix(TMarkovStateVariable ***list, TMarkovStateVariable *sv)
{
  int j, **states=(int**)dw_CreateArray_array(dw_DimA(list));
  for (j=dw_DimA(list)-1; j >= 0; j--)
    states[j]=CreateStateIndex(sv,list[j],list[j] ? dw_DimA(list[j]) : 0);
  return states;
}

/*
   Assumes
     states -  An m x sv->n_state_variables array of integers.  The integer
               states[j][k] will be non-zero of the kth state effects the jth
               group.
     sv     -  Valid pointer to TMarkovStateVariable structure.  It must be the
               case that the number of columns in states is equal to
               sv->n_state_variables.

   Returns
     An m x nstates array of integers that is a translation from the values of
     the state variable sv to the state variable that is a product of the state
     variables with non-zero entries in states.
*/
int** CreateTranslationMatrix_Flat(int **states, TMarkovStateVariable *sv)
{
  int **translation;
  int i, j, k;
  translation=dw_CreateRectangularArray_int(dw_DimA(states),sv->nstates);

  for (j=dw_DimA(states)-1; j >= 0; j--)
    for (i=sv->nstates-1; i >= 0; i--)
      for (translation[j][i]=k=0; k < sv->n_state_variables; k++)
    if (states[j][k])
      translation[j][i]=translation[j][i]*sv->state_variable[k]->nstates + sv->Index[i][k];

  return translation;
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/******************************** ThetaRoutines ********************************/
/*******************************************************************************/
ThetaRoutines* CreateThetaRoutines_empty(void)
{
  ThetaRoutines* pRoutines=(ThetaRoutines*)swzMalloc(sizeof(ThetaRoutines));

  if (pRoutines)
    {
      pRoutines->pLogConditionalLikelihood=NULL;
      pRoutines->pExpectationSingleStep=NULL;
      pRoutines->pDestructor=NULL;
      pRoutines->pDrawParameters=NULL;
      pRoutines->pLogPrior=NULL;
      pRoutines->pNumberFreeParametersTheta=NULL;
      pRoutines->pConvertFreeParametersToTheta=NULL;
      pRoutines->pConvertThetaToFreeParameters=NULL;
      pRoutines->pStatesChanged=NULL;
      pRoutines->pThetaChanged=NULL;
      pRoutines->pTransitionMatrixChanged=NULL;
      pRoutines->pValidTheta=NULL;
      pRoutines->pInitializeForwardRecursion=NULL;
      pRoutines->pPermuteTheta=NULL;
    }

  return pRoutines;
}
/*******************************************************************************/
/********************************* TStateModel *********************************/
/*******************************************************************************/

/*************************** TStateModel Destructors ***************************/
/*
   Frees TStateModel
*/
void FreeStateModel(TStateModel *model)
{
  if (model)
    {
      if (model->routines)
        {
          if (model->routines->pDestructor)
            model->routines->pDestructor(model->theta);
          swzFree(model->routines);
        }
      FreeMarkovStateVariable(model->sv);
      dw_FreeArray(model->V);
      dw_FreeArray(model->Z);
      dw_FreeArray(model->states_count);
      swzFree(model);
    }
}
/*******************************************************************************/

/************************** TStateModel Constructors ***************************/
/*
   Creates TStateModel structure.  The structures sv, routines, and theta become the
   property of the the TStateModel structure.  The structures sv and routines will
   be destroyed by the TStateModel destructor.  If the pDestructor field of the
   routines structure has been set, then the structure theta will be destroyed by
   TStateModel destructor.
*/
TStateModel* CreateStateModel_new(TMarkovStateVariable *sv, ThetaRoutines *routines, void *theta)
{
  int t;
  TStateModel *model;

  if (!(model=(TStateModel*)swzMalloc(sizeof(TStateModel))))
    {
      swz_fprintf_err("CreateStateModel():  Out of memory\n");
      swzExit(0);
    }

  model->sv=sv;
  model->theta=theta;
  model->routines=routines;

/*    //=== Common work space ===   ansi-c*/
  model->Z=dw_CreateArray_vector(sv->nobs+1);
  model->V=dw_CreateArray_vector(sv->nobs+1);
  for (t=sv->nobs; t >= 0; t--)
    {
      model->Z[t]=CreateVector(sv->nstates);
      model->V[t]=CreateVector(sv->nstates);
    }
  InitializeVector(model->V[0],1.0/(PRECISION)(sv->nstates));

  model->n_degenerate_draws=0;
  model->states_count=dw_CreateArray_int(sv->nstates);

/*    //== Set control variables   ansi-c*/
  model->ValidForwardRecursion=0;
  model->UseLogFreeParametersQ=0;

/*    //=== Obsolete fields ===   ansi-c*/
  model->parameters=theta;
  model->p=(TParameters*)swzMalloc(sizeof(TParameters));
  model->p->pLogConditionalLikelihood=model->routines->pLogConditionalLikelihood;
  model->p->pParameterDestructor=model->routines->pDestructor;
  model->p->pLogPrior=model->routines->pLogPrior;
  model->p->pNumberFreeParametersTheta=model->routines->pNumberFreeParametersTheta;
  model->p->pConvertFreeParametersToTheta=model->routines->pConvertFreeParametersToTheta;
  model->p->pConvertThetaToFreeParameters=model->routines->pConvertThetaToFreeParameters;
  model->p->pDrawParameters=model->routines->pDrawParameters;
  model->p->p=theta;

/*    //=== Set Transition matrix to prior mean ===   ansi-c*/
  SetTransitionMatrixToPriorMean(model);

  return model;
}
/*******************************************************************************/

/************************** TStateModel Notifications **************************/
void StatesChanged(TStateModel *model)
{
  if (model->routines->pStatesChanged) model->routines->pStatesChanged(model);
}

void TransitionMatricesChanged(TStateModel *model)
{
  model->ValidForwardRecursion=0;
}

void ThetaChanged(TStateModel *model)
{
  model->ValidForwardRecursion=0;
  if (model->routines->pThetaChanged) model->routines->pThetaChanged(model);
}
/*******************************************************************************/

/**************************** TStateModel Functions ****************************/
/*
   Assumes:
     model is a pointer to a properly initialized TStateModel.

   Results:
     Computes the following:

       V[t][i] = P(s[t]=i | Y[t],Z[t],theta,Q)  for 0 <= t <= T and 0 <= i < nstates
       Z[t][i] = P(s[t]=i | Y[t-1],Z[t-1],theta,Q)  for 0 < t <= T and 0 <= i < nstates
       L[t] = Sum(ln(Sum(P(y[t] | z[t],theta,s[t])*P(s[t] | Y[t-1],Z[t-1],theta,Q),0 <= s[t] < nstates)),0 < t <= T)
*/
void ForwardRecursion(TStateModel *model)
{
  int s, t;
  PRECISION scale, u;
  TMarkovStateVariable *sv=model->sv;

/*    //====== Initializes prior to forward recursion if necessary ======   ansi-c*/
  if (model->routines->pInitializeForwardRecursion)
    model->routines->pInitializeForwardRecursion(model);

/*    //====== Initialize L ======   ansi-c*/
  model->L=0;

/*    //====== Get ergodic distribution if necessary ======   ansi-c*/
  if (sv->UseErgodic)
    if (!Ergodic(model->V[0],sv->Q))
      {
    printf("ForwardRecursion():  Ergodic distribution does not exist.\n");
    swzExit(0);
      }

/*    //====== forward recursion ======   ansi-c*/
  for (t=1; t <= sv->nobs; t++)
    {
/*        //------ compute Z[t] ------   ansi-c*/
      ProductMV(model->Z[t],sv->Q,model->V[t-1]);

/*        //------ compute log conditional probabilities and scale ------   ansi-c*/
      scale=MINUS_INFINITY;
      for (s=sv->nstates-1; s >= 0; s--)
    {
      if (ElementV(model->Z[t],s) > 0.0)
        {
          u=LogConditionalLikelihood(s,t,model);
          scale=AddScaledLogs(1.0,scale,ElementV(model->Z[t],s),u);
          ElementV(model->V[t],s)=log(ElementV(model->Z[t],s)) + u;
        }
      else
        ElementV(model->V[t],s)=MINUS_INFINITY;
    }

/*        //------ update L ------   ansi-c*/
      model->L+=scale;

/*        //------ scale V[t] ------   ansi-c*/
      for (s=sv->nstates-1; s >= 0; s--)
    if (ElementV(model->V[t],s) != MINUS_INFINITY)
      ElementV(model->V[t],s)=exp(ElementV(model->V[t],s) - scale);
    else
      ElementV(model->V[t],s)=0.0;
    }

  model->ValidForwardRecursion=1;
}

/*
   Assumes:
     model is a pointer to a properly initialized TStateModel.

   Results:
     Draws the Markov state variables from the posterior distribution using
     backward recursion.  The state variable values are stored in S.

   Notes:
     If model->ValidForwardRecursion is not set, calls ForwardRecusion().
*/
void DrawStates(TStateModel *model)
{
  int i, j, t;
  PRECISION scale, u, s;
  TMarkovStateVariable *sv=model->sv;
  int *states_count=model->states_count;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

/*    //====== Backward recursion ======   ansi-c*/
  if ((u=dw_uniform_rnd()) >= (s=ElementV(model->V[t=sv->nobs],i=sv->nstates-1)))
    while (--i > 0)
      if (u < (s+=ElementV(model->V[t],i))) break;
  states_count[sv->S[t]=i]++;
  for (t--; t >= 0; t--)
    {
      scale=1.0/ElementV(model->Z[t+1],j=i);
      i=sv->nstates-1;
      if ((u=dw_uniform_rnd()) >= (s=ElementV(model->V[t],i)*ElementM(sv->Q,j,i)*scale))
    while (--i > 0)
      if (u < (s+=ElementV(model->V[t],i)*ElementM(sv->Q,j,i)*scale)) break;
      states_count[sv->S[t]=i]++;
    }


/*    //====== Backward recursion ======   ansi-c*/
/*   for (iteration=10000; iteration > 0; iteration--) */
/*     { */
/*       for (i=sv->nstates-1; i >= 0; i--) states_count[i]=0; */

/*       if ((u=dw_uniform_rnd()) >= (s=ElementV(model->V[t=sv->nobs],i=sv->nstates-1))) */
/*     while (--i > 0) */
/*       if (u < (s+=ElementV(model->V[t],i))) break; */
/*       states_count[sv->S[t]=i]++; */
/*       for (t--; t >= 0; t--) */
/*     { */
/*       scale=1.0/ElementV(model->Z[t+1],j=i); */
/*       i=sv->nstates-1; */
/*       if ((u=dw_uniform_rnd()) >= (s=ElementV(model->V[t],i)*ElementM(sv->Q,j,i)*scale)) */
/*         while (--i > 0) */
/*           if (u < (s+=ElementV(model->V[t],i)*ElementM(sv->Q,j,i)*scale)) break; */
/*       states_count[sv->S[t]=i]++; */
/*     } */

/*       states_count[sv->S[0]]--; */
/*       for (i=sv->nstates-1; i >= 0; i--) */
/*     if (!states_count[i]) break; */
/*       if (i < 0) break; */
/*       model->n_degenerate_draws++; */
/*     } */

/*    //====== Propagate states ======   ansi-c*/
  PropagateStates_SV(model->sv);

/*    //====== State change notification ======   ansi-c*/
  StatesChanged(model);
}

/*
   Assumes:
     model is a pointer to a properly initialized TStateModel.

   Results:
     Draws the vector of states S from the distribution defined by the transition
     matrix.

   Notes:
     The data or the values of the parameters other than sv->P play no role in
     the distribution of the vector of states.  If sv->P was drawn from the prior
     distribution, then the vector of states can be considered as being drawing
     from the prior distribution.
*/
void DrawStatesFromTransitionMatrix(TStateModel *model)
{
/*    //====== Draw states ======   ansi-c*/
  DrawStatesFromTransitionMatrix_SV(model->sv);

/*    //====== State change notification ======   ansi-c*/
  StatesChanged(model);
}

/*
   Assumes:
     model:  pointer to a properly initialized TStateModel structure

   Results:
    Draws the transition matrix from the posterior distribution, conditional on
    the states and theta.
*/
void DrawTransitionMatrix(TStateModel *model)
{
/*    //== Draw transition matrix using recursive call ===   ansi-c*/
  DrawTransitionMatrix_SV(model->sv);

/*    //====== Transition Matrix Change change notification ======   ansi-c*/
  TransitionMatricesChanged(model);
}

/*
   Assumes:
    model:  pointer to a properly initialized TStateModel structure

   Results:
    Draws the transition matrix sv->Q from the prior distribution.
*/
void DrawTransitionMatrixFromPrior(TStateModel *model)
{
/*    //====== Draw transition matrix using recursive call ======   ansi-c*/
  DrawTransitionMatrixFromPrior_SV(model->sv);

/*    //====== Transition matrix change notification ======   ansi-c*/
  TransitionMatricesChanged(model);
}

void SetTransitionMatrixToPriorMean(TStateModel *model)
{
/*    //====== Set transition matrix using recursive call ======   ansi-c*/
  SetTransitionMatrixToPriorMean_SV(model->sv);

/*    //====== Transition matrix change notification ======   ansi-c*/
  TransitionMatricesChanged(model);
}

void DrawTheta(TStateModel *model)
{
/*    //====== Draw theta ======   ansi-c*/
  model->routines->pDrawParameters(model);

/*    //====== Theta change notification ======   ansi-c*/
  ThetaChanged(model);
}

void DrawAll(TStateModel *model)
{
  DrawStates(model);
  DrawTransitionMatrix(model);
  DrawTheta(model);
}

/*******************************************************************************/
/******************************* Normalization  ********************************/
/*******************************************************************************/
int SetStateNormalizationMethod(int (*pGetNormalization)(int*, struct TStateModel_tag*),int (*pPermuteTheta)(int*, struct TStateModel_tag*),TStateModel *model)
{
  if (pGetNormalization && pPermuteTheta)
    {
      model->NormalizeStates=1;
      model->routines->pGetNormalization=pGetNormalization;
      model->routines->pPermuteTheta=pPermuteTheta;
    }
  else
    {
      model->NormalizeStates=0;
      model->routines->pGetNormalization=NULL;
      model->routines->pPermuteTheta=NULL;
    }
  return 1;
}

/*
   Returns 1 if states successfully normalized and zero otherwise.
*/
int NormalizeStates(TStateModel *model)
{
/*   int *p, rtrn=0; */

/*   if (!(model->NormalizeStates)) return 1; */

/*   if (p=(int*)swzMalloc(model->sv->nstates)) */
/*     { */
/*       if (model->routines->pGetNormalization(p,model)) */
/*     if (Permute_SV(p,model->sv)) */
/*       if (model->routiens->pPermuteTheta(p,model)) */
/*         rtrn=1; */
/*       swzFree(p); */
/*     } */

/*   return rtrn; */
  return 1;
}


/*******************************************************************************/
/**************************** Compute Probabilities ****************************/
/*******************************************************************************/

/*
   Assumes:
     t : 1 <= t <= nobs
     model : pointer to valid TStateModel

   Returns:
     If valid transition matrix, returns the natural logrithm of

                        P(y[t] | Y[t-1], Z[t], theta, Q).

     If the transition matrix is invalid, return minus infinity.  Note that minus
     infinity is also a valid return value.

   Notes:
     If not ValidForwardRecursion, then calls ForwardRecursion().  The time
     parameter t is one based.  Computation is based on the following:

      P(y[t] | Y[t-1], Z[t], theta, Q) =
         Sum(P(s[t]=s | Y[t-1], Z[t], theta, Q)
                * P(y[t] | Y[t-1], Z[t], theta, Q, s), 0 <= s < nstates)

      P(s[t]=s | Y[t-1], Z[t], theta, Q) = P(s[t]=s | Y[t-1], Z[t-1], theta, Q)
*/
PRECISION LogConditionalLikelihood_StatesIntegratedOut(int t, TStateModel *model)
{
  int s;
  PRECISION x;

/*    //====== Check if transition matrix is valid ======   ansi-c*/
  if (!(model->sv->valid_transition_matrix)) return MINUS_INFINITY;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

  x=LogConditionalLikelihood(0,t,model) + log(ProbabilityStateConditionalPrevious(0,t,model));

  for (s=model->sv->nstates-1; s > 0; s--)
    x=AddLogs(x,LogConditionalLikelihood(s,t,model) + log(ProbabilityStateConditionalPrevious(s,t,model)));

  return x;
}

/*
   Assumes:
     t : 1 <= t <= nobs
     model : pointer to valid TStateModel

   Returns:
    Returns

          E[y[t] | Y[t-1], Z[t], theta, Q]

     upon success and MINUS_INFINITY upon failure.  Check
     model->ValidTransitionMatrix to determine if a failure has occured.

   Notes:
     If not ValidForwardRecursion, then calls ForwardRecursion().  The time
     parameter t is one based.  Computation is based on the following:

      P(y[t] | Y[t-1], Z[t], theta, Q) =
         Sum(P(s[t]=s | Y[t-1], Z[t], theta, Q)
                * E(y[t] | Y[t-1], Z[t], theta, Q, s), 0 <= s < nstates)

      P(s[t]=s | Y[t-1], Z[t], theta, Q) = P(s[t]=s | Y[t-1], Z[t-1], theta, Q)
*/
TVector ExpectationSingleStep_StatesIntegratedOut(TVector y, int t, TStateModel *model)
{
  int s;
  TVector y_tmp;

/*    //====== Check if transition matrix is valid ======   ansi-c*/
  if (!(model->sv->valid_transition_matrix)) (TVector)NULL;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

  y=ExpectationSingleStep(y,0,t,model);
  ProductVS(y,y,ProbabilityStateConditionalPrevious(0,t,model));

  y_tmp=CreateVector(DimV(y));
  for (s=model->sv->nstates-1; s > 0; s--)
    {
      ExpectationSingleStep(y_tmp,s,t,model);
      ProductVS(y_tmp,y_tmp,ProbabilityStateConditionalPrevious(s,t,model));
      AddVV(y,y,y_tmp);
    }

  FreeVector (y_tmp);

  return y;
}

/*
   Assumes:
     s : 0 <= s < nstates
     t : 1 <= t <= nobs
     model : pointer to valid TStateModel

   Returns:
     Returns P(s[t] = s | Y[t], Z[t], theta, Q) if there is a valid transition
     matrix and -1 otherwise.

   Notes:
     If not ValidForwardRecursion, then calls ForwardRecursion().  The time
     parameter t is one based.
*/
PRECISION ProbabilityStateConditionalCurrent(int s, int t, TStateModel *model)
{
/*    //====== Check if transition matrix is valid ======   ansi-c*/
  if (!(model->sv->valid_transition_matrix)) return -1.0;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

  return ElementV(model->V[t],s);
}

/*
   Assumes:
     s : 0 <= s < nstates
     t : 1 <= t <= nobs
     model : pointer to valid TStateModel

   Returns:

          P(s[t] | Y[t-1], Z[t-1], theta, Q)

   Notes:
     If not ValidForwardRecursion, then calls ForwardRecursion().  The time
     parameter t is one based.
*/
PRECISION ProbabilityStateConditionalPrevious(int s, int t, TStateModel *model)
{
/*    //====== Check if transition matrix is valid ======   ansi-c*/
  if (!(model->sv->valid_transition_matrix)) return -1.0;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

  return ElementV(model->Z[t],s);
}

/*
   Assumes:
     S : vector of length at least nobs+1
     s : 0 <= s < nstates
     model : pointer to valid TStateModel

   Returns:
     Upon success, the vector P is returned and upon failure a null pointer is
     returned. If P is initially the null pointer, P is created.

   Results:

          P[t] = P(s[t] = s | Y[T], Z[T], theta, Q)

     for 0 <= t <= T.

   Notes:
     If not ValidForwardRecursion, then calls ForwardRecursion().  The time
     parameter t is one based.  The computations use the following facts:

      P(s[t]=k | Y[T], Z[T], theta, Q)
              = Sum(P(s[t]=k, s[t+1]=i | Y[T], Z[T], theta, Q), 0 <= i < nstates)

      P(s[t]=k, s[t+1]=i | Y[T], Z[T], theta, Q)
                    = P(s[t]=k | Y[T], Z[T], theta, Q, s[t+1]=i)
                                             * P(s[t+1]=i | Y[T], Z[T], theta, Q)

      P(s[t]=k | Y[t], Z[t], theta, Q, s[t+1]=i)
                                     = P(s[t]=k | Y[T], Z[T], theta, Q, s[t+1]=i)

      P(s[t]=k | Y[t], Z[t], theta, Q, s[t+1]=i) = Q[i][k]
          * P(s[t]=k | Y[t], Z[t], theta, Q) / P(s[t+1]=i | Y[t], Z[t], theta, Q)

     The third of these facts follows from

       P(s[t] | Y[t], Z[t], theta, Q, s[t+1])
                                       = P(s[t] | Y[T], Z[T], theta, Q, S(t+1,T))

     where S(t+1,T)={s[t+1],...,s[T]}.  This is proven in the related paper.
*/
TVector ProbabilitiesState(TVector P, int s, TStateModel *model)
{
  int nobs=model->sv->nobs, nstates=model->sv->nstates;
  int i, k, t;
  TVector p1, p2, p3;

/*    //====== Check if transition matrix is valid ======   ansi-c*/
  if (!(model->sv->valid_transition_matrix)) return (TVector)NULL;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

  if (!P)
    P=CreateVector(nobs+1);
  else
    if (DimV(P) != nobs+1) return (TVector)NULL;

  p1=CreateVector(nstates);
  p2=CreateVector(nstates);

/*    // set p1 and the S[nobs]   ansi-c*/
  EquateVector(p1,model->V[nobs]);
  ElementV(P,nobs)=ElementV(p1,s);

  for (t=nobs-1; t >= 0; t--)
    {
/*        // s[t] = k and s[t+1] = i   ansi-c*/
      for (k=nstates-1; k >= 0; k--)
    {
      for (ElementV(p2,k)=0.0, i=nstates-1; i >= 0; i--)
        if (ElementV(model->Z[t+1],i) > 0.0)
          ElementV(p2,k)+=ElementV(p1,i)*ElementM(model->sv->Q,i,k)/ElementV(model->Z[t+1],i);
      ElementV(p2,k)*=ElementV(model->V[t],k);
    }

/*        // interchange p1 and p2   ansi-c*/
      p3=p1;
      p1=p2;
      p2=p3;

      ElementV(P,t)=ElementV(p1,s);
    }

  FreeVector(p1);
  FreeVector(p2);

  return P;
}

/*
   Assumes:
     S     : integer array of length model->sv->nobs+1
     model : pointer to valid TStateModel

   Returns:
     If the transition matrix is valid, returns the natural logrithm of
     P(S[T] | Y[T], Z[T], theta, Q).  If the transition matrix is invalid, return
     minus infinity.  Note that minus infinity is also a valid return value.

   Notes:
     If not ValidForwardRecursion, then calls ForwardRecursion().  The time
     parameter t is one based.  The computations use the following facts:

      P(S[T] | Y[T], Z[T], theta, Q)
                  = Product(P(s[t] | Y[T], Z[T], theta, Q, S(t+1,T)),1 <= t <= T)

      P(s[t] | Y[t], Z[t], theta, Q, s[t+1])
                                       = P(s[t] | Y[T], Z[T], theta, Q, S(t+1,T))

      P(s[t] | Y[t], Z[t], theta, Q, s[t+1]) = Q[s(t+1)][s(t)]
              * P(s[t] | Y[t], Z[t], theta, Q) / P(s[t+1] | Y[t], Z[t], theta, Q)

     where S(t+1,T)={s[t+1],...,s[T]}.  The second of these facts is proven in
     the related paper.
*/
PRECISION LogConditionalProbabilityStates(int *S, TStateModel *model)
{
  PRECISION rtrn;
  TMarkovStateVariable *sv=model->sv;
  TMatrix Q=sv->Q;
  int t=sv->nobs;

/*    //====== Check if transition matrix is valid ======   ansi-c*/
  if (!(model->sv->valid_transition_matrix)) return MINUS_INFINITY;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

/*    //====== Log probability ======   ansi-c*/
  rtrn=log(ElementV(model->V[t],S[t]));
  while (--t >= 0)
    if (ElementV(model->Z[t+1],S[t+1]) > 0.0)
      rtrn+=log(ElementM(Q,S[t+1],S[t])*ElementV(model->V[t],S[t])/ElementV(model->Z[t+1],S[t+1]));
    else
      return MINUS_INFINITY;

/*    //====== Probability ======   ansi-c*/
/*   rtrn=ElementV(model->V[t],S[t]); */
/*   while (--t >= 0) */
/*     rtrn*=ElementM(Q,S[t+1],S[t])*ElementV(model->V[t],S[t])/ElementV(model->Z[t+1],S[t+1]); */

  return rtrn;
}

/*
   Assumes
     model:  pointer to valid TStateModel structure

   Returns
     The natural logrithm of

                   P(S[T] | Theta, Q)
*/
PRECISION LogConditionalPrior_S(TStateModel *model)
{
  int t, *S=model->sv->S;
  TMatrix Q=model->sv->Q;
  PRECISION p=-log(model->sv->nstates);
  for (t=model->sv->nobs; t > 0; t--)
    p+=log(ElementM(Q,S[t],S[t-1]));
  return p;
}

/*
   Assumes:
     model : pointer to valid TStateModel structure

   Returns:
     The natural logrithm of

        P(Y[T] | Z[T], Theta, Q, S[T])

*/
PRECISION LogLikelihood(TStateModel *model)
{
  TMarkovStateVariable *sv=model->sv;
  int t, *S=sv->S;
  PRECISION loglikelihood=0.0;
  for (t=sv->nobs; t > 0; t--)
    loglikelihood+=LogConditionalLikelihood(S[t],t,model);
  return loglikelihood;
}

/*
   Assumes:
     model : pointer to valid TStateModel structure

   Returns:
     If the transition matrix is valid, returns the natural logrithm of

                        P(Y[T] | Z[t], theta, Q).

     If the transition matrix is invalid, then minus infinity is returned.  Note
     that minus infinity is also a valid return value.

   Notes:
     If not ValidForwardRecursion, then calls ForwardRecursion().
*/
PRECISION LogLikelihood_StatesIntegratedOut(TStateModel *model)
{
/*    //====== Check if transition matrix is valid ======   ansi-c*/
  if (!(model->sv->valid_transition_matrix)) return MINUS_INFINITY;

/*    //====== Check if ForwardRecursion() has been called ======   ansi-c*/
  if (!model->ValidForwardRecursion) ForwardRecursion(model);

  return model->L;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/




/*******************************************************************************/
/************************** Free parameters routines ***************************/
/*******************************************************************************/
/*
    Returns the number of free parameters in sv->Q.  Each column of Q is a linear
    combination of random vectors with Dirichlet distributions.  A n-dimensional
    random vector with the Dirichlet distributions only has n-1 degrees of
    freedom.
*/
int NumberFreeParametersQ(TStateModel *model)
{
  int i, n=0;
  TVector* ba=model->sv->ba;
  for (i=dw_DimA(ba)-1; i >= 0; i--) n+=DimV(ba[i]);
  return n - dw_DimA(ba);
}

/*
    Assumes
      model:  valid pointer to a TStateModel structure.
      f:  array of PRECISION of length at least NumberFreeParametersQ(sv).

    Results
      Copies the free parameters of Q into the array f.

    Notes
      If model->ValidTransitionMatrix is not set, then routine prints error
      message and exits.
*/
void ConvertQToFreeParameters(TStateModel *model, PRECISION *f)
{
  int i, k;
  TVector* ba=model->sv->ba;
  if (!(model->sv->valid_transition_matrix))
    {
      swz_fprintf_err("ConvertQToFreeParameters():  Transition matrices not valid.\n");
      swzExit(0);
    }
  for (i=0; i < dw_DimA(ba); f+=k, i++)
    memcpy(f,pElementV(ba[i]),(k=DimV(ba[i])-1)*sizeof(PRECISION));
}

/*
    Assumes
      sv : valid pointer to a TMarkovStateVariable structure.
      f:  array of PRECISION of length at least NumberFreeParametersQ(sv).  The
          elements of f must be non-negative an  0 and 1 inclusive.

    Results
      Converts the array of free parameters into the matrix Q.

    Notes
      The routine TransitionMatricesChanged() is called.  If any of the elements
      of f are negative or if one minus their sum is less than the number of
      states times machine epsilon, then model->sv->valid_transition_matrix is
      set to zero.
*/
void ConvertFreeParametersToQ(TStateModel *model, PRECISION *f)
{
  int i, j, k;
  PRECISION scale;
  TVector* ba=model->sv->ba;
  for (i=0; i < dw_DimA(ba); f+=k, i++)
    {
      memcpy(pElementV(ba[i]),f,(k=DimV(ba[i])-1)*sizeof(PRECISION));
      for (scale=1.0, j=k-1; j >= 0; j--)
    if (f[j] < 0)
      {
        InvalidateTransitionMatrices_SV(model->sv);
        TransitionMatricesChanged(model);
        return;
      }
    else
      scale-=f[j];
      if (scale < 0)
    if (scale < -DimV(ba[i])*MACHINE_EPSILON)
      {
        InvalidateTransitionMatrices_SV(model->sv);
        TransitionMatricesChanged(model);
        return;
      }
    else
      ElementV(ba[i],k)=0.0;
      else
    ElementV(ba[i],k)=scale;
    }
  Update_Q_from_B_SV(model->sv);
  if (!(model->sv->valid_transition_matrix)) ValidateTransitionMatrices_SV(model->sv);
  TransitionMatricesChanged(model);
}

/*
    Assumes:
      sv : valid pointer to a TMarkovStateVariable structure.
      f : array of PRECISION of length at least NumberFreeParametersQ(model).

    Results:
      Copies the free parameters of Q into the array f.

    Notes:
      The natural logrithms of free parameters of Q are stored in f.
*/
void ConvertQToLogFreeParameters(TStateModel *model, PRECISION *f)
{
  int i, j, k;
  TVector* ba=model->sv->ba;

  if (!(model->sv->valid_transition_matrix))
    {
      swz_fprintf_err("ConvertQToFreeParameters():  Transition matrices not valid.\n");
      swzExit(0);
    }

  for (i=0; i < dw_DimA(ba); f+=k, i++)
    for (j=(k=DimV(ba[i])-1)-1; j >= 0; j--)
      f[j]=(ElementV(ba[i],j) > 0) ? log(ElementV(ba[i],j)) : MINUS_INFINITY;
}

/*
    Assumes:
      model : valid pointer to a TStateModel structure.
      f : array of PRECISION of length at least NumberFreeParametersQ(model). The
          elements of f must non-positive and the sum of the natural exponential
          of the elements must the less than or equal to 1.

    Results:
      Converts the array of free parameters into the matrix Q.
*/
void ConvertLogFreeParametersToQ(TStateModel *model, PRECISION *f)
{
  int i, j, k;
  PRECISION scale;
  TVector* ba=model->sv->ba;
  for (i=0; i < dw_DimA(ba); f+=k, i++)
    {
      for (scale=1.0, j=(k=DimV(ba[i])-1)-1; j >= 0; j--)
    if (f[j] > 0)
      {
        InvalidateTransitionMatrices_SV(model->sv);
        TransitionMatricesChanged(model);
        return;
      }
    else
      scale-=(ElementV(ba[i],j)=exp(f[j]));
      if (scale < 0)
    {
      InvalidateTransitionMatrices_SV(model->sv);
      TransitionMatricesChanged(model);
      return;
    }
      ElementV(ba[i],k)=scale;
    }
  Update_Q_from_B_SV(model->sv);
  if (!(model->sv->valid_transition_matrix)) ValidateTransitionMatrices_SV(model->sv);
  TransitionMatricesChanged(model);
}

void ConvertFreeParametersToTheta(TStateModel *model, PRECISION *f)
{
  model->routines->pConvertFreeParametersToTheta(model,f);
  ThetaChanged(model);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/****************************** Utility Routines *******************************/
/*******************************************************************************/
/*
   Checks that FreeDim, NonZeroIndex, and MQ satisfies the appropriate
   conditions.  Returns 1 if the conditions are satisfied and 0 otherwise.
*/
int CheckRestrictions(int* FreeDim, int** NonZeroIndex, TMatrix MQ, int nstates)
{
  int i, j, k, q, r, total_free;
  PRECISION sum, common_sum, total_sum;

/*    //====== Check for null pointer and sizes ======   ansi-c*/
  if (!FreeDim || !MQ || (RowM(MQ) != nstates) || (ColM(MQ) != nstates))
    return 0;
  if (!NonZeroIndex || (dw_DimA(NonZeroIndex) != nstates))
    return 0;
  else
    for (i=nstates-1; i >= 0; i--)
      if (dw_DimA(NonZeroIndex[i]) != nstates) return 0;

/*    // Checks FreeDim[i] > 0   ansi-c*/
/*    // Computes total_free = FreeDim[0] + ... + FreeDim[dw_DimA(FreeDim)-1]   ansi-c*/
  for (total_free=0, i=dw_DimA(FreeDim)-1; i >= 0; i--)
    if (FreeDim[i] <= 0)
      return 0;
    else
      total_free+=FreeDim[i];

/*    // Checks -1 <= NonZeroIndex[i][j] < total_free.   ansi-c*/
/*    // Checks NonZeroIndex[i][j] >= 0, ==> MQ[i][j] > 0.   ansi-c*/
/*    // Checks NonZeroIndex[i][j] = -1, ==> MQ[i][j] = 0.   ansi-c*/
  for (j=0; j < nstates; j++)
    for (i=0; i < nstates; i++)
      if ((NonZeroIndex[i][j] < -1) || (total_free <= NonZeroIndex[i][j]))
        return 0;
      else
        if (NonZeroIndex[i][j] >= 0)
          {
            if (ElementM(MQ,i,j) <= 0.0) return 0;
          }
        else
          {
            if (ElementM(MQ,i,j) != 0.0) return 0;
          }

/*    //====== Check that column sums are correct ======   ansi-c*/
  for (j=0; j < nstates; j++)
    {
      total_sum=0.0;
      for (q=r=k=0; k < total_free; k++)
        if (k == q)
          {
            common_sum=0.0;
            for (i=0; i < nstates; i++)
              if (NonZeroIndex[i][j] == k)
        common_sum+=ElementM(MQ,i,j);
            q+=FreeDim[r++];
            total_sum+=common_sum;
          }
        else
          {
            sum=0.0;
            for (i=0; i < nstates; i++)
          if (NonZeroIndex[i][j] == k)
        sum+=ElementM(MQ,i,j);
            if (fabs(sum - common_sum) > SQRT_MACHINE_EPSILON) return 0;
          }
      if (fabs(total_sum - 1.0) > SQRT_MACHINE_EPSILON) return 0;
    }

  return 1;
}

/*
   Checks that the Prior matrix is of the correct size and that all its elements
   are positive.
*/
int CheckPrior(TMatrix Prior, int nstates)
{
  int i, j;
  if (!Prior || (RowM(Prior) != nstates) || (ColM(Prior) != nstates)) return 0;
  for (i=0; i < nstates; i++)
    for (j=0; j < nstates; j++)
      if (ElementM(Prior,i,j) <= 0) return 0;
  return 1;
}

/*
   Checks that both Prior and NonZeroIndex are of the correct size and that Prior
   elements are large enough.
*/
int CheckPriorOnFreeParameters(TMatrix Prior, int** NonZeroIndex, int nstates)
{
  int i, j, q=0;
  PRECISION alpha;

/*    // non-null pointers and matrices of correct size   ansi-c*/
  if (!Prior || (RowM(Prior) != nstates) || (ColM(Prior) != nstates))
    return 0;
  if (!NonZeroIndex || (dw_DimA(NonZeroIndex) != nstates))
    return 0;
  else
    for (i=nstates-1; i >= 0; i--)
      if (dw_DimA(NonZeroIndex[i]) != nstates) return 0;
  for (j=0; j < nstates; j++)
    for (i=0; i < nstates; i++)
      if (q < NonZeroIndex[i][j]) q=NonZeroIndex[i][j];
  for ( ; q >= 0; q--)
    {
      alpha=1.0;
      for (j=0; j < nstates; j++)
        for (i=0; i < nstates; i++)
          if (NonZeroIndex[i][j] == q)
            alpha+=ElementM(Prior,i,j)-1.0;
      if (alpha <= 0) return 0;
    }
  return 1;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/****************************** Utility Routines *******************************/
/*******************************************************************************/
/*
   Assumes:
     v : An n vector or null pointer
     P : An n x n transition matrix.  It must be the case that P(i,j) >= 0 for
         every i and j and P(0,j) + ... + P(n-1,j) = 1 for every j.

   Results:
     Computes the ergodic distribution of the transition matrix P if it exists.
     The ergodic distribution will exist if there is a unique vector v such that
     the elements of v are non-negative, sum to one and (P - I)v = 0.

   Notes:
     If w is the n dimensional row vector of ones, then w(P - I) = 0.  This
     implies that there exists a non-zero v such that (P - I)v = 0.  It is easy
     to show that if (P - I)v = 0 and v1 is the positive component of v and v2 is
     the negative component of v, then both (P - I)v1 = 0 and (P - I)v2 = 0.  So
     there is always exists a v such that the elements of v are non-negative, sum
     to one, (P - I)v = 0.  However, a unique such element might not exist.

     This routine does not check for uniqueness.

     This version uses the LU decomposition, which is fast and fairly stable.
*/
TVector Ergodic(TVector v, TMatrix P)
{
  TVector rtrn;
  TMatrix IP;
  TPermutation S;
  PRECISION sum;
  int i, j, k;

  if (!P)
    {
      dw_Error(NULL_ERR);
      return (TVector)NULL;
    }
  if (RowM(P) != ColM(P))
    {
      dw_Error(SIZE_ERR);
      return (TVector)NULL;
    }
  if (!v)
    {
      if (!(rtrn=CreateVector(RowM(P))))
    return (TVector)NULL;
    }
  else
    if (RowM(P) != DimV(v))
      {
    dw_Error(SIZE_ERR);
    return (TVector)NULL;
      }
    else
      rtrn=v;

  IP=EquateMatrix((TMatrix)NULL,P);
  for (i=RowM(P)-1; i >= 0; i--)
    ElementM(IP,i,i)-=1.0;
  S=CreatePermutation(RowM(P));

  if (LU(S,IP,IP))
    {
      for (k=0; k < RowM(P); k++)
    if (fabs(ElementM(IP,k,k)) < SQRT_MACHINE_EPSILON) break;

      if (k < RowM(P))
    {
      for (i=DimV(rtrn)-1; i > k; i--)
        ElementV(rtrn,i)=0.0;
      ElementV(rtrn,i--)=1.0;
      for ( ; i >= 0; i--)
        {
          for (sum=-ElementM(IP,i,k), j=k-1; j > i; j--) sum-=ElementM(IP,i,j)*ElementV(rtrn,j);
          ElementV(rtrn,i)=sum/ElementM(IP,i,i);
        }

      for (sum=1.0, i=k-1; i >= 0; i--)
        if (ElementV(rtrn,i) < 0.0)
          ElementV(rtrn,i)=0.0;
        else
          sum+=ElementV(rtrn,i);

      ProductSV(rtrn,1.0/sum,rtrn);
      FreePermutation(S);
      FreeMatrix(IP);
      return rtrn;
    }
    }

  FreePermutation(S);
  FreeMatrix(IP);

  if (!v) FreeVector(rtrn);
  return (TVector) NULL;
}

/*
   Assumes:
     v : An n vector or null pointer
     P : An n x n transition matrix.  It must be the case that P(i,j) >= 0 for
         every i and j and P(0,j) + ... + P(n-1,j) = 1 for every j.

   Results:
     Computes the ergodic distribution of the transition matrix P if it exists.
     The ergodic distribution will exist if there is a unique vector v such that
     the elements of v are non-negative, sum to one and (P - I)v = 0.

   Notes:
     If w is the n dimensional row vector of ones, then w(P - I) = 0.  This
     implies that there exists a non-zero v such that (P - I)v = 0.  It is easy
     to show that if (P - I)v = 0 and v1 is the positive component of v and v2 is
     the negative component of v, then both (P - I)v1 = 0 and (P - I)v2 = 0.  So
     there is always exists a v such that the elements of v are non-negative, sum
     to one, (P - I)v = 0.  However, a unique such element might not exist.

     This version use the SVD and checks for uniqueness.
*/
TVector Ergodic_SVD(TVector v, TMatrix P)
{
  TVector rtrn;
  TMatrix IP, U, V;
  PRECISION sum;
  int i, j, k;

  if (!P)
    {
      dw_Error(NULL_ERR);
      return (TVector)NULL;
    }
  if (RowM(P) != ColM(P))
    {
      dw_Error(SIZE_ERR);
      return (TVector)NULL;
    }
  if (!v)
    {
      if (!(rtrn=CreateVector(RowM(P))))
    return (TVector)NULL;
    }
  else
    if (RowM(P) != DimV(v))
      {
    dw_Error(SIZE_ERR);
    return (TVector)NULL;
      }
    else
      rtrn=v;

  IP=EquateMatrix((TMatrix)NULL,P);
  for (i=RowM(P)-1; i >= 0; i--)
    ElementM(IP,i,i)-=1.0;
  U=CreateMatrix(RowM(P),RowM(P));
  V=CreateMatrix(RowM(P),RowM(P));

  if (SVD(U,rtrn,V,IP))
    {
      k=(ElementV(rtrn,0) < SQRT_MACHINE_EPSILON) ? 1 : 0;
      for (j=0, i=RowM(P)-1; i > 0; i--)
    if (ElementV(rtrn,i) < SQRT_MACHINE_EPSILON)
      {
        k++;
        if (ElementV(rtrn,i) < ElementV(rtrn,j)) j=i;
      }

      if (k == 1)
    {
      ColumnVector(rtrn,V,j);
      for (sum=ElementV(rtrn,0), i=DimV(rtrn)-1; i > 0; i--)
        sum+=ElementV(rtrn,i);
      if (sum > 0.0)
        {
          for (sum=0.0, i=DimV(rtrn)-1; i >= 0; i--)
        if (ElementV(rtrn,i) < 0.0)
          ElementV(rtrn,i)=0;
        else
          sum+=ElementV(rtrn,i);
        }
      else
        {
          for (sum=0.0, i=DimV(rtrn)-1; i >= 0; i--)
        if (ElementV(rtrn,i) > 0.0)
          ElementV(rtrn,i)=0;
        else
          sum+=ElementV(rtrn,i);
        }

      ProductSV(rtrn,1.0/sum,rtrn);

      FreeMatrix(V);
      FreeMatrix(U);
      FreeMatrix(IP);
      return rtrn;
    }
    }

  FreeMatrix(V);
  FreeMatrix(U);
  FreeMatrix(IP);

  if (!v) FreeVector(rtrn);
  return (TVector) NULL;
}

/*
   Assumes:
     v : An n vector or null pointer
     P : An n x n transition matrix.  It must be the case that P(i,j) >= 0 for
         every i and j and P(0,j) + ... + P(n-1,j) = 1 for every j.

   Results:
     Computes the ergodic distribution of the transition matrix P if it exists.
     The ergodic distribution will exist if there is a unique vector v such that
     the elements of v are non-negative, sum to one and (P - I)v = 0.

   Notes:
     If w is the n dimensional row vector of ones, then w(P - I) = 0.  This
     implies that there exists a non-zero v such that (P - I)v = 0.  It is easy
     to show that if (P - I)v = 0 and v1 is the positive component of v and v2 is
     the negative component of v, then both (P - I)v1 = 0 and (P - I)v2 = 0.  So
     there is always exists a v such that the elements of v are non-negative, sum
     to one, (P - I)v = 0.  However, a unique such element might not exist.

     This version use the SVD and checks for uniqueness.
*/
TVector* ErgodicAll_SVD(TMatrix P)
{
  TVector d, *rtrn;
  TMatrix IP, U, V;
  PRECISION sum;
  int i, j, k, p, n;

  if (!P)
    {
      dw_Error(NULL_ERR);
      return (TVector*)NULL;
    }
  if (RowM(P) != ColM(P))
    {
      dw_Error(SIZE_ERR);
      return (TVector*)NULL;
    }

  IP=EquateMatrix((TMatrix)NULL,P);
  for (i=RowM(P)-1; i >= 0; i--)
    ElementM(IP,i,i)-=1.0;
  U=CreateMatrix(RowM(P),RowM(P));
  V=CreateMatrix(RowM(P),RowM(P));
  d=CreateVector(RowM(P));

  if (SVD(U,d,V,IP))
    {
      for (k=0, i=RowM(P)-1; i >= 0; i--)
    if (ElementV(d,i) < SQRT_MACHINE_EPSILON)
      {
        for (p=n=0, j=RowM(P)-1; j >= 0; j--)
          if (ElementM(V,j,i) > SQRT_MACHINE_EPSILON)
        p=1;
          else
        if (ElementM(V,j,i) < -SQRT_MACHINE_EPSILON)
          n=1;
        k+=n+p;
      }

      if (k > 0)
    {
      rtrn=dw_CreateArray_vector(k);
      for (k=0, i=RowM(P)-1; i >= 0; i--)
        if (ElementV(d,i) < SQRT_MACHINE_EPSILON)
          {
        for (p=n=0, j=RowM(P)-1; j >= 0; j--)
          if (ElementM(V,j,i) > SQRT_MACHINE_EPSILON)
            p=1;
          else
            if (ElementM(V,j,i) < -SQRT_MACHINE_EPSILON)
              n=1;

        if (p > 0)
          {
            rtrn[k]=ColumnVector((TVector)NULL,V,i);
            for (sum=0.0, j=DimV(rtrn[k])-1; j >= 0; j--)
              if (ElementV(rtrn[k],j) < 0.0)
            ElementV(rtrn[k],j)=0.0;
              else
            sum+=ElementV(rtrn[k],j);
            if (sum < SQRT_MACHINE_EPSILON) { printf("Sum is not positive %le\n",sum); getchar(); }
            ProductSV(rtrn[k],1.0/sum,rtrn[k]);
            k++;
          }

        if (n > 0)
          {
            rtrn[k]=ColumnVector((TVector)NULL,V,i);
            for (sum=0.0, j=DimV(rtrn[k])-1; j >= 0; j--)
              if (ElementV(rtrn[k],j) > 0.0)
            ElementV(rtrn[k],j)=0.0;
              else
            sum+=ElementV(rtrn[k],j);
            if (sum > -SQRT_MACHINE_EPSILON) { printf("Sum is not negative %le\n",sum); getchar(); }
            ProductSV(rtrn[k],1.0/sum,rtrn[k]);
            k++;
          }
          }

      FreeVector(d);
      FreeMatrix(V);
      FreeMatrix(U);
      FreeMatrix(IP);
      return rtrn;
    }
    }

  FreeVector(d);
  FreeMatrix(V);
  FreeMatrix(U);
  FreeMatrix(IP);

  return (TVector*)NULL;
}

/*
   Assumes:
     Q : n-vector or null pointer
     A : n-vector with elements greater than -1

   Results:
     Fills the vector Q with a draw from the Dirichlet distribution with
     parameters given by A.  The density of the Dirichlet distribution is


        Product(Gamma(A[i]), 0 <= i < n)
       ---------------------------------- * Product(Q[i]^A[i], 0 <= i < n)
          Gamma(Sum(A[i], 0 <= i < n))

   Returns:
     Upon success, returns Q.  Upon failure, returns null.  If Q is null, then it
     is created.

   Notes:
     The arguments A and Q do not have to be distinct.

     Sometimes the Dirichlet parameters are restricted to be greater than -1, as
     is done here.  Other times the paramaters are restricted to be positive and
     the exponent in the formula for the density is A[i] - 1.  Care must be taken
     to ensure the the correct form is used.
*/
TVector DrawDirichletVector(TVector Q, TVector A)
{
  int i;
  PRECISION sum;
  TVector rtrn;
  if (!A)
    {
      dw_Error(NULL_ERR);
      return (TVector)NULL;
    }
  if (!Q)
    {
      if (!(rtrn=CreateVector(DimV(A))))
    return (TVector)NULL;
    }
  else
    if (DimV(Q) != DimV(A))
      {
    dw_Error(SIZE_ERR);
    return (TVector)NULL;
      }
    else
      rtrn=Q;

  for (sum=0.0, i=DimV(A)-1; i >= 0; i--)
    sum+=ElementV(rtrn,i)=dw_gamma_rnd(ElementV(A,i));

  if (sum > 0)
    for (sum=1.0/sum, i=DimV(A)-1; i >= 0; i--)
      ElementV(rtrn,i)*=sum;
  else
    {
      if (!Q) FreeVector(rtrn);
      return (TVector)NULL;
    }
  return rtrn;
}

/*
   Assumes:
     Q : array of vectors
     A : arrat if vectorsr with elements greater than -1

   Results:
     Fills the vector Q[i] with a draw from the Dirichlet distribution with
     parameters given by A[i].

   Returns:
     Upon success, returns Q.  Upon failure, returns null.  If Q is null, then it
     is created.

   Notes:
     The arguments A and Q do not have to be distinct.

     Sometimes the Dirichlet parameters are restricted to be greater than -1, as
     is done here.  Other times the paramaters are restricted to be positive and
     the exponent in the formula for the density is A[i] - 1.  Care must be taken
     to ensure the the correct form is used.
*/
TVector* DrawIndependentDirichletVector(TVector *Q, TVector *A)
{
  int i;
  TVector v, *rtrn;
  if (!A)
    {
      dw_Error(NULL_ERR);
      return (TVector*)NULL;
    }
  if (!Q)
    {
      if (!(rtrn=dw_CreateArray_vector(dw_DimA(A))))
    return (TVector*)NULL;
    }
  else
    if (dw_DimA(Q) != dw_DimA(A))
      {
    dw_Error(SIZE_ERR);
    return (TVector*)NULL;
      }
    else
      rtrn=Q;

  for (i=dw_DimA(Q) - 1; i >= 0; i--)
    if (v=DrawDirichletVector(rtrn[i],A[i]))
      rtrn[i]=v;
    else
      {
    if (!Q) dw_FreeArray(rtrn);
    return (TVector*)NULL;
      }

  return rtrn;
}

/*
   Assumes:
     Q : n dimensional vector with elements between 0 and 1.
     A : n dimensional vector with positive elements

   Returns:
     The log value of the Dirichlet density with parameters A evaluated at Q.

   Notes:
     The Dirichlet density is given by

         Gamma(Sum(A[i],0 <= i < n))
      --------------------------------- Product(Q[i]^(A[i] - 1), 0 <= i < n)
       Product(Gamma(A[i]), 0 <= i < n)
*/
PRECISION LogDirichlet_pdf(TVector Q, TVector A)
{
  PRECISION sum=0.0, log_pdf=0.0, x, y;
  int i;
  if (!Q || !A)
    {
      dw_Error(NULL_ERR);
      return MINUS_INFINITY;
    }
  if (DimV(Q) != DimV(A))
    {
      dw_Error(SIZE_ERR);
      return MINUS_INFINITY;
    }
  for (i=DimV(Q)-1; i >= 0; i--)
    {
      sum+=(x=ElementV(A,i));
      log_pdf-=dw_log_gamma(x);
      if ((y=ElementV(Q,i)) > 0.0)
    log_pdf+=(x-1.0)*log(y);
      else
    if (x > 1.0)
      return PLUS_INFINITY;
    else
      if (x < 1.0)
        return MINUS_INFINITY;
    }
  return log_pdf + dw_log_gamma(sum);
}

/*
   Assumes:
     Q : m dimensional vector array.  The elements of each vector in the array
         are between 0 and 1.
     A : m dimensional vector array.  The elements of each vector in the array
         are positive

   Returns:
     The sum of the log values of the Dirichlet densities with parameters A[i]
     evaluated at Q[i].
*/
PRECISION LogIndependentDirichlet_pdf(TVector *Q, TVector *A)
{
  int i;
  PRECISION log_pdf=0.0;
  for (i=dw_DimA(Q)-1; i >= 0; i--)
    log_pdf+=LogDirichlet_pdf(Q[i],A[i]);
  return log_pdf;
}

/*
   Returns ln(exp(a) + exp(b)) computed to avoid overflow.  If
   a = ln(c) and b = ln(d), as is usually the case, then the
   routine returns ln(c + d).

*/
PRECISION AddLogs(PRECISION a, PRECISION b)
{
  return (a > b) ? a + log(1.0 + exp(b-a)) : b + log(exp(a-b) + 1.0);
}

/*
   Returns ln(x*exp(a) + y*exp(b)) computed to avoid overflow.  If a = ln(c) and
   b = ln(d), as is usually the case, then the routine returns ln(x*c + y*d).
   Assumes that x*exp(a) + y*exp(b) is positive, but no checking is done.  This
   condition will always be satisfied if both x and y are positive.

*/
PRECISION AddScaledLogs(PRECISION x, PRECISION a, PRECISION y, PRECISION b)
{
  return (a > b) ? a + log(x + y*exp(b-a)) : b + log(x*exp(a-b) + y);
}

/*
   Assumes:
    There are DimV(p) states, denoted by 0, ... , Dim(p)-1.  The probability of
    state i occuring is given by p[i].  It must be the case that the sum of the
    p[i] is equal to one.

   Returns:
    A random draw of one of the states.
*/
int DrawDiscrete(TVector p)
{
  int i=DimV(p)-1;
  PRECISION u=dw_uniform_rnd(), s=ElementV(p,i);
  while (i > 0)
    if (u < s)
      return i;
    else
      s+=ElementV(p,--i);
  return 0;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/****************************** Obsolete Routines ******************************/
/*******************************************************************************
/*
   Creates TStateModel
*/
TStateModel* CreateStateModel(TMarkovStateVariable *sv, TParameters *p)
{
  ThetaRoutines *routines;
  TStateModel *model;

  routines=CreateThetaRoutines_empty();
  routines->pLogConditionalLikelihood=p->pLogConditionalLikelihood;
  routines->pDestructor=p->pParameterDestructor;
  routines->pLogPrior=p->pLogPrior;
  routines->pNumberFreeParametersTheta=p->pNumberFreeParametersTheta;
  routines->pConvertFreeParametersToTheta=p->pConvertFreeParametersToTheta;
  routines->pConvertThetaToFreeParameters=p->pConvertThetaToFreeParameters;
  routines->pDrawParameters=p->pDrawParameters;

  model=CreateStateModel_new(sv,routines,p->p);
  return model;
}

void FreeParameters(TParameters *p)
{
  if (p)
    {
      if (p->pParameterDestructor) p->pParameterDestructor(p->p);
      swzFree(p);
    }
}

/*
   Constructs a TParameters type using passed data.  If pParameterDestructor is
   null, then the last argument will not be freed by FreeParameters().  If
   pParameterDestructor is not null, then the last argument will be freed by
   FreeParameters().
*/
TParameters* CreateParameters(PRECISION (*pLogConditionalLikelihood)(int,int,struct TStateModel_tag*),
                              void (*pParameterDestructor)(void*),
                  PRECISION (*pLogPrior)(struct TStateModel_tag*),
                  int (*pNumberFreeParametersTheta)(struct TStateModel_tag*),
                  void (*pConvertFreeParametersToTheta)(struct TStateModel_tag*, PRECISION*),
                  void (*pConvertThetaToFreeParameters)(struct TStateModel_tag*, PRECISION*),
                              void (*pDrawParameters)(struct TStateModel_tag*),
                              void *parameters)
{
  TParameters* p=(TParameters*)swzMalloc(sizeof(TParameters));
  if (p)
    {
      p->pLogConditionalLikelihood=pLogConditionalLikelihood;
      p->pParameterDestructor=pParameterDestructor;
      p->pLogPrior=pLogPrior;
      p->pNumberFreeParametersTheta=pNumberFreeParametersTheta;
      p->pConvertFreeParametersToTheta=pConvertFreeParametersToTheta;
      p->pConvertThetaToFreeParameters=pConvertThetaToFreeParameters;
      p->pDrawParameters=pDrawParameters;
      p->p=parameters;
    }
  return p;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


