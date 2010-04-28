
#include "VARio_matlab.h"
#include "switchio.h"

#include "dw_error.h"
#include "dw_ascii.h"

#include <stdlib.h>
#include <string.h>

#include "modify_for_mex.h"

static int ReadError_VARio_matlab(char *id)
{
  char *errmsg, *fmt="Error after line identifier ""%s""";
  sprintf(errmsg=(char*)malloc(strlen(fmt) + strlen(id) - 1),fmt,id);
  dw_UserError(errmsg);
  free(errmsg);
  return 1;
}

TStateModel* Combine_matlab_standard(char *matlabfile, char *standardfile)
{
  FILE *f_in;
  TStateModel *model;
  TMarkovStateVariable *sv, ***coef_sv, ***var_sv;
  T_VAR_Parameters *p;
  char *id;
  int *IV, **States;
  int nlags, nvars, nexg, npre, nstates, nobs, i, j, n, SimsZha=1, RandomWalk=1, flag;
  PRECISION scalar_zeta_a_prior, scalar_zeta_b_prior, lambda_prior;
  int **coef_states, **var_states;
  TMatrix *U, *V, *W, *A0_prior, *Aplus_prior, X, Y, S;
  TVector zeta_a_prior, zeta_b_prior;

  //=== Open matlab input file
  f_in=dw_OpenTextFile(matlabfile);

  //=== Read sizes ===//
  id="//== lags, nvar, nStates, T ==//";
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %d %d %d %d ",&nlags,&nvars,&nstates,&nobs) != 4)) ReadError_VARio_matlab(id);

  //=== A single constant term ===//
  nexg=1;
  npre=nvars * nlags + nexg;

  //=== Restrictions - U[j] ===//
  IV=dw_CreateArray_int(nvars);
  id="//== n0const: nvar-by-1 ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,IV)) ReadError_VARio_matlab(id);
  id="//== Uiconst: cell(nvar,1) and nvar-by-n0const(i) for the ith cell (equation) ==//";
  if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id);
  U=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (!dw_ReadMatrix(f_in,U[j]=CreateMatrix(nvars,IV[j]))) ReadError_VARio_matlab(id);
  dw_FreeArray(IV);

  //=== Restrictions - V[j] (V[j] should be an npre x npre identity matrix) ===//
  IV=dw_CreateArray_int(nvars);
  id="//== npconst: nvar-by-1 ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,IV)) ReadError_VARio_matlab(id);
  for (j=nvars-1; j >= 0; j--) 
    if (IV[j] != npre) SimsZha=0; 
  V=dw_CreateArray_matrix(nvars);
  if (SimsZha)
    {
      for (j=nvars-1; j >= 0; j--)
	V[j]=IdentityMatrix((TMatrix)NULL,npre);
    }
  else
    {
      id="//== Viconst: cell(nvar,1) and ncoef-by-n0const(i) for the ith cell (equation) ==//";
      if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id);
      for (j=0; j < nvars; j++)
	if (!dw_ReadMatrix(f_in,V[j]=CreateMatrix(npre,IV[j]))) ReadError_VARio_matlab(id);
    }
  dw_FreeArray(IV);

  //=== Restrictions - W[j] (Random walk specification) ===//
  InitializeMatrix(S=CreateMatrix(npre,nvars),0.0);
  for (j=nvars-1; j >= 0; j--) ElementM(S,j,j)=-1.0;
  W=dw_CreateArray_matrix(nvars);
  for (j=nvars-1; j >= 0; j--)
    W[j]=EquateMatrix((TMatrix)NULL,S);
  FreeMatrix(S);

  //====== Priors ======
  id="//== gxia: alpha parameter for gamma prior of xi ==//";
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %lf ",&scalar_zeta_a_prior) != 1)) ReadError_VARio_matlab(id);
  id="//== gxib: beta parameter for gamma prior of xi ==//";
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %lf ",&scalar_zeta_b_prior) != 1)) ReadError_VARio_matlab(id);
  zeta_a_prior=CreateVector(nvars);
  zeta_b_prior=CreateVector(nvars);
  for (j=nvars-1; j >= 0; j--)
    {
      ElementV(zeta_a_prior,j)=scalar_zeta_a_prior;
      ElementV(zeta_b_prior,j)=scalar_zeta_b_prior;
    }

  id="//== H0barconstcell: cell(nvar,1) and n-by-n for the ith cell (equation) ==//";
  if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id);
  A0_prior=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (!dw_ReadMatrix(f_in,A0_prior[j]=CreateMatrix(nvars,nvars))) ReadError_VARio_matlab(id);

  id="//== Hpbarconstcell: cell(nvar,1) and ncoef-by-ncoef for the ith cell (equation) ==//";
  if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id);
  Aplus_prior=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (!dw_ReadMatrix(f_in,Aplus_prior[j]=CreateMatrix(npre,npre))) ReadError_VARio_matlab(id);

  // Initialize Y
  id="//== Yleft -- Y: T-by-nvar ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,Y=CreateMatrix(nobs,nvars))) ReadError_VARio_matlab(id);

  // Initialize X
  id="//== Xright -- X: T-by-ncoef ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,X=CreateMatrix(nobs,npre))) ReadError_VARio_matlab(id);

  //=== Sims-Zha specification ===
  id="//== glamdasig: sigma parameter for normal prior of lamda ==//";
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %lf ",&lambda_prior) != 1)) ReadError_VARio_matlab(id);
  lambda_prior*=lambda_prior;

  //=== Close matlab input file ===
  fclose(f_in);

  //=== Open standard input file
  f_in=dw_OpenTextFile(standardfile);

  //=== Create Markov state variable ===//
  sv=CreateMarkovStateVariable_File(f_in,(char*)NULL,nobs);
  
  //====== coefficient/variance state variables ======
  id="//== Controlling states variables for coefficients ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,States=dw_CreateRectangularArray_int(nvars,sv->n_state_variables))) 
    ReadError_VARio_matlab(id);
  coef_sv=(TMarkovStateVariable ***)dw_CreateArray_array(nvars);
  for (j=nvars-1; j >= 0; j--)
    {
      for (n=i=0; i < sv->n_state_variables; i++)
	if (States[j][i]) n++;
      if (n > 0)
	{
	  coef_sv[j]=(TMarkovStateVariable **)dw_CreateArray_pointer(n,NULL);
	  for (n=i=0; i < sv->n_state_variables; i++)
	    if (States[j][i]) coef_sv[j][n++]=sv->state_variable[i];
	}
    }
  coef_states=CreateTranslationMatrix(coef_sv,sv);
  dw_FreeArray(States);
  dw_FreeArray(coef_sv);

  id="//== Controlling states variables for variance ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,States=dw_CreateRectangularArray_int(nvars,sv->n_state_variables))) 
    ReadError_VARio_matlab(id);
  var_sv=(TMarkovStateVariable ***)dw_CreateArray_array(nvars);
  for (j=nvars-1; j >= 0; j--)
    {
      for (n=i=0; i < sv->n_state_variables; i++)
	if (States[j][i]) n++;
      if (n > 0)
	{
	  var_sv[j]=(TMarkovStateVariable **)dw_CreateArray_pointer(n,NULL);
	  for (n=i=0; i < sv->n_state_variables; i++)
	    if (States[j][i]) var_sv[j][n++]=sv->state_variable[i];
	}
    }
  var_states=CreateTranslationMatrix(var_sv,sv);
  dw_FreeArray(States);
  dw_FreeArray(var_sv);

  //=== Close standard input file ===
  fclose(f_in);

  //=== Create T_VAR_Parameters structure ===
  flag=SimsZha ? SPEC_SIMS_ZHA : 0;
  flag|=RandomWalk ? SPEC_RANDOM_WALK : 0;
  p=CreateTheta_VAR(flag,nvars,nlags,nexg,sv->nstates,sv->nobs,coef_states,var_states,U,V,W,Y,X);
  if (flag & SPEC_SIMS_ZHA)
    SetPriors_VAR_SimsZha(p,A0_prior,Aplus_prior,zeta_a_prior,zeta_b_prior,lambda_prior);
  else
    SetPriors_VAR(p,A0_prior,Aplus_prior,zeta_a_prior,zeta_b_prior);

  //p=Create_VAR_Parameters(nvars,nlags,nexg,sv->nstates,sv->nobs,U,V,W,Zeta_a_prior,Zeta_b_prior,A0_prior,Aplus_prior,Y,X,coef_states,var_states);
  //SetupSimsZhaSpecification(p,delta_prior*delta_prior);

  //=== Create TStateModel ===
  model=CreateStateModel_new(sv,CreateRoutines_VAR(),p);

  //=== Print Model specifications to file ===

  //=== Free memory ===
  FreeMatrix(X);
  FreeMatrix(Y);
  dw_FreeArray(Aplus_prior);
  dw_FreeArray(A0_prior);
  FreeVector(zeta_b_prior);
  FreeVector(zeta_a_prior);
  dw_FreeArray(W);
  dw_FreeArray(V);
  dw_FreeArray(U);

  return model;
}

/*
    This reads the constant parameters from filename, which was created
    from Matlab and then sets all the parameters to the constant parameters.
*/
void ReadConstantParameters(char *filename, TStateModel *model)
{
  char *id;
  int i, j, s;
  FILE *f_in;
  TMatrix A0, Aplus;
  T_VAR_Parameters *p=(T_VAR_Parameters*)(model->theta);

  if (!(f_in=fopen(filename,"rt")))
    {
      printf("Unable to read the input data file: %s\n", filename);
      exit(0);
    }

  // A0
  id="//== A0hat: nvar-by-nvar ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,A0=CreateMatrix(p->nvars,p->nvars))) ReadError_VARio_matlab(id);
  for (j=p->nvars-1; j >= 0; j--)
    for (s=p->n_coef_states[j]-1; s >= 0; s--)
      for (i=p->nvars-1; i >= 0; i--)
	ElementV(p->A0[j][s],i)=ElementM(A0,i,j);
  FreeMatrix(A0);

  // Aplus
  id="//== Aphat: ncoef(lags*nvar+1)-by-nvar ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,Aplus=CreateMatrix(p->npre,p->nvars))) ReadError_VARio_matlab(id);
  for (j=p->nvars-1; j >= 0; j--)
    for (s=p->n_coef_states[j]-1; s >= 0; s--)
      for (i=p->npre-1; i >= 0; i--)
	ElementV(p->Aplus[j][s],i)=ElementM(Aplus,i,j);
  FreeMatrix(Aplus);  

  // Zeta
  for (j=p->nvars-1; j >= 0; j--)
    for (s=p->n_var_states[j]-1; s >= 0; s--)
      p->Zeta[j][s]=1.0;

  // b0, bplus, lambda, and psi
  Update_b0_bplus_from_A0_Aplus(p);
  if (p->Specification & SPEC_SIMS_ZHA) Update_lambda_psi_from_bplus(p);

  // Flags
  p->valid_parameters=1;

  // Transition matrix
  SetTransitionMatrixToPriorMean(model);

  ThetaChanged(model);
}

/*
   Create Model from Matlab data file
*/
TStateModel* CreateStateModel_VAR_matlab(char *filename)
{
  T_VAR_Parameters *p;
  FILE *f_in;
  char *id;
  TMatrix PriorTransitionMatrix, S;
  int *IV, **IM;
  PRECISION scalar_Zeta_a_prior, scalar_Zeta_b_prior, lambda_prior;
  int i, j, nvars, nlags, nexg, npre, nobs, nstates;                                                        
  TMatrix *U, *V, *W;                                                 
  TVector Zeta_a_prior, Zeta_b_prior; 
  TMatrix *A0_prior, *Aplus_prior;
  TMatrix Y, X; 
  int **coef_states, **var_states;  
  TMarkovStateVariable *sv;                                                              

  //=== Open file ===
  if (!(f_in=fopen(filename,"rt")))
    {
      printf("Unable to read the input data file: %s\n", filename);
      exit(0);
    }

  //=== Read sizes ===//
  id="//== lags, nvar, nStates, T ==//";
  if (!dw_SetFilePosition(f_in,id)
      || (fscanf(f_in," %d %d %d %d ",&nlags,&nvars,&nstates,&nobs) != 4)) ReadError_VARio_matlab(id);

  //=== A single constant term ===//
  nexg=1;
  npre=nvars * nlags + nexg;

  //=== Restrictions - U[j] ===//
  IV=dw_CreateArray_int(nvars);
  id="//== n0const: nvar-by-1 ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,IV)) ReadError_VARio_matlab(id);
  id="//== Uiconst: cell(nvar,1) and nvar-by-n0const(i) for the ith cell (equation) ==//";
  if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id);
  U=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (!dw_ReadMatrix(f_in,U[j]=CreateMatrix(nvars,IV[j]))) ReadError_VARio_matlab(id);
  dw_FreeArray(IV);

  //=== Restrictions - V[j] (V[j] should be an npre x npre identity matrix) ===//
  IV=dw_CreateArray_int(nvars);
  id="//== npconst: nvar-by-1 ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,IV)) ReadError_VARio_matlab(id);
  for (j=nvars-1; j >= 0; j--) 
    if (IV[j] != npre) 
      {
	fprintf(stderr,"V[%d] not %d x %d\n",j,npre,npre);
	exit(0);
      }
  V=dw_CreateArray_matrix(nvars);
  for (j=nvars-1; j >= 0; j--)
    V[j]=IdentityMatrix((TMatrix)NULL,npre);
  dw_FreeArray(IV);

  //=== Restrictions - W[j] (Random walk specification) ===//
  InitializeMatrix(S=CreateMatrix(npre,nvars),0.0);
  for (j=nvars-1; j >= 0; j--) ElementM(S,j,j)=-1.0;
  W=dw_CreateArray_matrix(nvars);
  for (j=nvars-1; j >= 0; j--)
    W[j]=EquateMatrix((TMatrix)NULL,S);
  FreeMatrix(S);

  //=== Create TMarkovStateVariable  ===//
  PriorTransitionMatrix=CreateMatrix(nstates,nstates);
  id="//== Galpha: nStates-by-nStates ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,PriorTransitionMatrix)) ReadError_VARio_matlab(id);
  sv=CreateMarkovStateVariable_NoRestrictions(nstates,nobs,PriorTransitionMatrix);
  FreeMatrix(PriorTransitionMatrix);

  //====== regime/shock state variables ======
  coef_states=dw_CreateRectangularArray_int(nvars,nstates);
  var_states=dw_CreateRectangularArray_int(nvars,nstates);
  IM=dw_CreateRectangularArray_int(nvars,2);
  id="//== indxEqnTv_m: nvar-by-2 ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,IM)) ReadError_VARio_matlab(id);
  for (j=nvars-1; j >= 0; j--)
    switch (IM[j][1])
      {
      case 1:
	for (i=nstates-1; i >= 0; i--)
	  coef_states[j][i]=var_states[j][i]=0;
	break;
      case 2:
        for (i=nstates-1; i >= 0; i--)
          {
	    coef_states[j][i]=0;
	    var_states[j][i]=i;
	  }
	break;
      case 3:
        for (i=nstates-1; i >= 0; i--)
          {
	    coef_states[j][i]=i;
	    var_states[j][i]=0;
	  }
	break;
      case 4:
	fprintf(stderr,"Case %d not implimented.\n",4);
	exit(0);
      default:
	fprintf(stderr,"Unknown type.\n");
	exit(0);
      }
  dw_FreeArray(IM);

  //====== Priors ======
  id="//== gxia: alpha parameter for gamma prior of xi ==//";
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %lf ",&scalar_Zeta_a_prior) != 1)) ReadError_VARio_matlab(id);
  id="//== gxib: beta parameter for gamma prior of xi ==//";
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %lf ",&scalar_Zeta_b_prior) != 1)) ReadError_VARio_matlab(id);
  Zeta_a_prior=CreateVector(nvars);
  Zeta_b_prior=CreateVector(nvars);
  for (j=nvars-1; j >= 0; j--)
    {
      ElementV(Zeta_a_prior,j)=scalar_Zeta_a_prior;
      ElementV(Zeta_b_prior,j)=scalar_Zeta_b_prior;
    }

  id="//== H0barconstcell: cell(nvar,1) and n-by-n for the ith cell (equation) ==//";
  if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id);
  A0_prior=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (!dw_ReadMatrix(f_in,A0_prior[j]=CreateMatrix(nvars,nvars))) ReadError_VARio_matlab(id);

  id="//== Hpbarconstcell: cell(nvar,1) and ncoef-by-ncoef for the ith cell (equation) ==//";
  if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id);
  Aplus_prior=dw_CreateArray_matrix(nvars);
  for (j=0; j < nvars; j++)
    if (!dw_ReadMatrix(f_in,Aplus_prior[j]=CreateMatrix(npre,npre))) ReadError_VARio_matlab(id);

/*   //=========================== Checks */
/*   TMatrix *H0, *Ui; */
/*   int *n0s; */
/*   TMatrix Sigma, XX, YY, ZZ; */
/*   int i, ii, jj; */
/*   PRECISION max; */

/*   id="//== n0s: nvar-by-1 ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,n0s=dw_CreateArray_int(nvars))) ReadError_VARio_matlab(id); */

/*   id="//== H0tldcell_inv: cell(nvar,1) and n0s(i)-by-n0s(i) for the ith cell ==//"; */
/*   if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id); */
/*   H0=dw_CreateArray_matrix(nvars); */
/*   for (j=0; j < nvars; j++) */
/*     if (!dw_ReadMatrix(f_in,H0[j]=CreateMatrix(n0s[j],n0s[j]))) ReadError_VARio_matlab(id); */

/*   id="//== Ui: cell(nvar,1) and nvar*nStates-by-n0s(i) for the ith cell ==//"; */
/*   if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id); */
/*   Ui=dw_CreateArray_matrix(nvars); */
/*   for (j=0; j < nvars; j++) */
/*     if (!dw_ReadMatrix(f_in,Ui[j]=CreateMatrix(nvars*nstates,n0s[j]))) ReadError_VARio_matlab(id); */

/*   Sigma=CreateMatrix(nvars*nstates,nvars*nstates); */
/*   XX=CreateMatrix(nvars,nvars); */
/*   for (j=0; j < nvars; j++) */
/*     { */
/*       InitializeMatrix(Sigma,0.0); */
/*       Inverse_LU(XX,A0_prior[j]); */
/*       for (i=0; i < nstates; i++) */
/* 	for (ii=0; ii < nvars; ii++) */
/* 	  for (jj=0; jj < nvars; jj++) */
/* 	    ElementM(Sigma,i*nvars+ii,i*nvars+jj)=ElementM(XX,ii,jj); */
/*       YY=TransposeProductMM((TMatrix)NULL,Ui[j],Sigma); */
/*       ZZ=ProductMM((TMatrix)NULL,YY,Ui[j]); */

/*       fprintf(stdout,"Computed[%d]\n",j); dw_PrintMatrix(stdout,ZZ,"%le "); fprintf(stdout,"\n"); */
/*       fprintf(stdout,"File[%d]\n",j); dw_PrintMatrix(stdout,H0[j],"%le "); fprintf(stdout,"\n"); */
/*       max=0.0; */
/*       for (ii=0; ii < RowM(ZZ); ii++) */
/* 	  for (jj=0; jj < ColM(ZZ); jj++) */
/* 	    if (max < fabs(ElementM(H0[j],ii,jj) - ElementM(ZZ,ii,jj))) max=fabs(ElementM(H0[j],ii,jj) - ElementM(ZZ,ii,jj)); */
/*       fprintf(stdout,"H0: max[%d] = %le\n",j,max); */

/*       FreeMatrix(ZZ); */
/*       FreeMatrix(YY); */
/*       getc(stdin); */
/*     } */
/*   //exit(0); */

/*   id="//== nps: nvar-by-1 ==//"; */
/*   if (!dw_SetFilePosition(f_in,id) || !dw_ReadArray(f_in,n0s=dw_CreateArray_int(nvars))) ReadError_VARio_matlab(id); */

/*   id="//== Hptldcell_inv: cell(nvar,1) and nps(i)-by-nps(i) for the ith cell ==//"; */
/*   if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id); */
/*   H0=dw_CreateArray_matrix(nvars); */
/*   for (j=0; j < nvars; j++) */
/*     if (!dw_ReadMatrix(f_in,H0[j]=CreateMatrix(n0s[j],n0s[j]))) ReadError_VARio_matlab(id); */

/*   id="//== Vi: cell(nvar,1) and k*nStates-by-nps(i) for the ith cell ==//"; */
/*   if (!dw_SetFilePosition(f_in,id)) ReadError_VARio_matlab(id); */
/*   Ui=dw_CreateArray_matrix(nvars); */
/*   for (j=0; j < nvars; j++) */
/*     if (!dw_ReadMatrix(f_in,Ui[j]=CreateMatrix(npre*nstates,n0s[j]))) ReadError_VARio_matlab(id); */

/*   Sigma=CreateMatrix(npre*nstates,npre*nstates); */
/*   XX=CreateMatrix(npre,npre); */
/*   for (j=0; j < nvars; j++) */
/*     { */
/*       InitializeMatrix(Sigma,0.0); */
/*       Inverse_LU(XX,Aplus_prior[j]); */
/*       for (i=0; i < nstates; i++) */
/* 	for (ii=0; ii < npre; ii++) */
/* 	  for (jj=0; jj < npre; jj++) */
/* 	    ElementM(Sigma,i*npre+ii,i*npre+jj)=ElementM(XX,ii,jj); */
/*       YY=TransposeProductMM((TMatrix)NULL,Ui[j],Sigma); */
/*       ZZ=ProductMM((TMatrix)NULL,YY,Ui[j]); */

/*       max=0.0; */
/*       for (ii=0; ii < RowM(ZZ); ii++) */
/* 	  for (jj=0; jj < ColM(ZZ); jj++) */
/* 	    if (max < fabs(ElementM(H0[j],ii,jj) - ElementM(ZZ,ii,jj))) max=fabs(ElementM(H0[j],ii,jj) - ElementM(ZZ,ii,jj)); */
/*       fprintf(stdout,"max[%d] = %le\n",j,max); */

/*       FreeMatrix(ZZ); */
/*       FreeMatrix(YY); */
/*       getc(stdin); */
/*     } */
/*   exit(0); */
/*   //=========================== Checks */

  // Initialize Y
  id="//== Yleft -- Y: T-by-nvar ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,Y=CreateMatrix(nobs,nvars))) ReadError_VARio_matlab(id);

  // Initialize X
  id="//== Xright -- X: T-by-ncoef ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,X=CreateMatrix(nobs,npre))) ReadError_VARio_matlab(id);

  //=== Create T_VAR_Parameters structure ===
  p=CreateTheta_VAR(SPEC_SIMS_ZHA | SPEC_RANDOM_WALK,nvars,nlags,nexg,nstates,nobs,coef_states,var_states,U,V,W,Y,X);

  //=== Sims-Zha specification ===
  id="//== glamdasig: sigma parameter for normal prior of lamda ==//";
  if (!dw_SetFilePosition(f_in,id) || (fscanf(f_in," %lf ",&lambda_prior) != 1)) ReadError_VARio_matlab(id);
  SetPriors_VAR_SimsZha(p,A0_prior,Aplus_prior,Zeta_a_prior,Zeta_b_prior,lambda_prior*lambda_prior);

  //=== Close input file ===
  fclose(f_in);

  //=== Create TStateModel ===
  return CreateStateModel_new(sv,CreateRoutines_VAR(),p);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
