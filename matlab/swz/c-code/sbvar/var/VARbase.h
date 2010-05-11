
#ifndef __VAR_BASE_MODEL__
#define __VAR_BASE_MODEL__

#include "switch.h"
#include "swzmatrix.h"
#include "dw_matrix_array.h"

#define standard_ordering 1  // for future implementation
#define SPEC_RANDOM_WALK  0x00000001
#define SPEC_SIMS_ZHA     0x00000002

//=== Normalization types (must be mutually exclusive) ===//
#define VAR_NORMALIZATION_NONE    0x00000000
#define VAR_NORMALIZATION_WZ      0x00000001

typedef struct
{
  //====== Model specification ======
  int Specification;
  int *IsIdentity_V;

  //====== Free parameter specification (for future implementation) ======
  int FreeParameterType;

  //====== Sizes ======
  int nvars;
  int nlags;
  int npre;
  int nobs;
  int nstates;

  //====== State variable translation ======
  int*  n_var_states;          // nvars                 n_var_states[j] is the number of variance states for column j
  int** var_states;            // nvars x n_states      translation table for variance states
  int*  n_coef_states;         // nvars                 n_coef_states[j] is the number of coefficients states for column j
  int** coef_states;           // nvars x n_states      translation table for coefficient states
  int   n_A0_states;           //                       number of states for the matrix A0
  int*  A0_states;             // n_states              translation table for the matrix A0
  int** A0_column_states;      // nvars x n_A0_states   translation table from determinant of A0 states to coefficient states

  //====== Parameters ======
  PRECISION** Zeta;            // nvars x n_var_states[j]
  TVector** A0;                // nvars x n_coef_states[j] x nvars
  TVector** Aplus;             // nvars x n_coef_states[j] x npre

  //====== Free parameters ======
  int* dim_b0;
  TVector** b0;
  int* dim_bplus;
  TVector** bplus;

  //====== Priors ======
  TVector  zeta_a_prior;            //
  TVector  zeta_b_prior;            //
  TMatrix* A0_prior;                // A0_prior[j] = constant parameter variance of the normal prior on the jth column of A0
  TMatrix* Aplus_prior;             // Aplus_prior[j] = constant parameter variance of the normal prior on the jth column of Aplus

  //====== Identifying restrictions ======
  TMatrix *U;
  TMatrix *V;
  TMatrix *W;

  //====== Sims-Zha specification parameters and workspace ======
  TVector** lambda;                               // nvars x n_coef_states[j] array of nvars-dimensional vectors
  TVector*  constant;                             // nvars x n_coef_states[j] -- constant[j][k] == psi[j][npre - 1 + k]
  TVector*  psi;                                  // nvars x (npre - 1 + n_coef_states[j])
  PRECISION lambda_prior;                         // prior variance of each element of lambda
  PRECISION inverse_lambda_prior;
  TMatrix*  inverse_psi_prior;

  //====== Normalization ======
  int normalization_type;                         // type of normalization used
  int normalized;                                 // type of normalization actually used to normalize current draw
  TVector** Target;                               // nvar x n_coef_states[j] array of nvars-dimensional vectors
  int** flipped;                                  // nvar x n_coef_states[j] array of integers
  int WZ_inconsistancies;

  //====== Workspace  ======
  PRECISION log_prior_constant;                    // Constant of integrate for the log prior
  PRECISION minus_half_nvars_times_log2pi;         // Constant used in LogConditionalProbability functions
  TVector  inverse_zeta_b_prior;                   // inverse_zeta_b_prior = 1.0/zeta_b_prior
  TMatrix* inverse_b0_prior;                       // inverse_b0_prior = U'[j]*Inverse(A0_prior[j])*U[j]
  TMatrix* inverse_bplus_prior;                    // inverse_bplus_prior = V'[j]*Inverse(Aplus_prior[j])*V[j]

  TVector log_abs_det_A0;                          // log(abs(det(A0[k])))

  PRECISION*** A0_dot_products;                    // A0_dot_products[t][j][k] = Y'[t] * A0[j][k]
  PRECISION*** Aplus_dot_products;                 // Aplus_dot_products[t][j][k] = X'[t] * Aplus[j][k]

  // A0 Metropolis Info
  PRECISION** A0_Metropolis_Scale;
  int Total_A0_Metropolis_Draws;
  int** A0_Metropolis_Jumps;

  // State dependent fields
  TMatrix* YY;                                     // YY[k] = sum(Y[t]*Y'[t], 1 <= t <= nobs and S[t] == k)
  TMatrix* XY;                                     // YX[k] = sum(X[t]*Y'[t], 1 <= t <= nobs and S[t] == k)
  TMatrix* XX;                                     // XX[k] = sum(X[t]*X'[t], 1 <= t <= nobs and S[t] == k)
  int* T;                                          // T[k]  = number of t with 1 <= t <= nobs and S[t] == k
  int *S;                                          // S[t]  = state variable used to compute YY, XY, XX, and T
  TMatrix* yy;                                     // yy[t] = Y[t]*Y'[t]
  TMatrix* xy;                                     // xy[t] = X[t]*Y'[t]
  TMatrix* xx;                                     // xx[t] = X[t]*X'[t]

  // Flags for validity of workspace fields
  int valid_log_abs_det_A0;                        // Invalid after A0 changes
  int valid_dot_products;                          // Invalid after A0 or Aplus changes
  int valid_state_dependent_fields;                // Invalid after states change
  int valid_state_dependent_fields_previous;       // Initially invalid.
  int valid_parameters;                            // Initially invalid.  Valid after successful read or draw of parameters.
                                                   //   Parametes are invalid if Zeta is negative or if they do not satisfy
                                                   //   the normalization.

  //=== Data ===
  TVector* Y;         // Y[t] nvar vector of time t data for 1 <= t <= T
  TVector* X;         // X[t] npre vector of time t predetermined variables for 1 <= t <= T

} T_VAR_Parameters;

// Constructors-Destructors
void FreeTheta_VAR(T_VAR_Parameters *p);
ThetaRoutines* CreateRoutines_VAR(void);
T_VAR_Parameters* CreateTheta_VAR(int flag, int nvars, int nlags, int nexg, int nstates, int nobs,    // Specification and Sizes
                  int **coef_states, int **var_states,                                // Translation Tables
                  TMatrix *U, TMatrix *V, TMatrix *W,                                 // Restrictions
                  TMatrix Y, TMatrix X);                                              // Data
int** CreateTranslationMatrix_Flat(int **states, TMarkovStateVariable *sv);

void SetPriors_VAR(T_VAR_Parameters *theta, TMatrix* A0_prior, TMatrix* Aplus_prior, TVector zeta_a_prior, TVector zeta_b_prior);
void SetPriors_VAR_SimsZha(T_VAR_Parameters *theta, TMatrix* A0_prior, TMatrix* Aplus_prior, TVector zeta_a_prior,
               TVector zeta_b_prior, PRECISION lambda_prior);


TStateModel* CreateConstantModel(TStateModel *model);
TStateModel* ExpandModel_VAR(TStateModel *model, TStateModel *restricted_model, int s);


void SetupSimsZhaSpecification(T_VAR_Parameters *p, PRECISION lambda_prior);

PRECISION LogConditionalProbability_VAR(int i, int t, TStateModel *model);
TVector ExpectationSingleStep_VAR(TVector y, int s, int t, TStateModel *model);


void DrawParameters_VAR(TStateModel *model);
void InitializeParameters_VAR(T_VAR_Parameters *p);

// Priors
void SetLogPriorConstant_VAR(T_VAR_Parameters *p);
PRECISION LogPrior_VAR(TStateModel *model);

// Normalization
int IsNormalized_VAR(T_VAR_Parameters *p);
int Normalize_VAR(T_VAR_Parameters *p);
void Setup_No_Normalization(T_VAR_Parameters *p);
void Setup_WZ_Normalization(T_VAR_Parameters *p, TVector **A0);
int WZ_Normalize(T_VAR_Parameters *p);

// Notification
void StatesChanged_VAR(TStateModel *model);
void ThetaChanged_VAR(TStateModel *model);
void InitializeForwardRecursion_VAR(TStateModel *model);

// Utility Routines
int Reset_VAR_Improper_Distribution_Counter(void);
int Get_VAR_Improper_Distribution_Counter(void);
void Increment_Verbose(void);
void SetVerboseFile(FILE *f);

// Optimization
int NumberFreeParametersVAR(TStateModel *model);
void FreeParametersToVAR(TStateModel *model, PRECISION *f);
void VARToFreeParameters(TStateModel *model, PRECISION *f);
int ZetaIndex(T_VAR_Parameters *p);
int ZetaLength(T_VAR_Parameters *p);

//PRECISION ComputeConstantSimsZha(TStateModel *model);

//
void PsiDeltaToAplus(TStateModel *model);

// Impulse Response
TMatrix ComputeImpulseResponseReducedForm(TMatrix R, int h, TMatrix A0_Xi_inv, TMatrix B, int nlags);
TMatrix ComputeImpulseResponseStructural(TMatrix R, int h, TMatrix A0, TMatrix Aplus, TVector Xi, int nlags);
TMatrix ComputeImpulseResponse(TMatrix R, int h, int k, TStateModel *model);
TMatrix ComputeVarianceDecomposition(TMatrix X, TMatrix IR, int nvars);

// Simulation
void DrawZeta_Aplus(TStateModel *model);
void DrawZeta_DotProducts(TStateModel *model);
void AdaptiveMetropolisScale(TStateModel *model, int iterations, int period, int verbose, FILE *f_posterior);
void SetupMetropolisInformation(PRECISION **Scale, T_VAR_Parameters *p);
void ResetMetropolisInformation(T_VAR_Parameters *p);
PRECISION LogKernel_A0_DotProducts(int j, int k, TStateModel *model);
PRECISION LogKernel_A0(int j, int k, TStateModel *model);
void DrawA0_Metropolis(TStateModel *model);
void DrawAplus(TStateModel *model);
void Draw_psi(TStateModel *model);
void Draw_lambda(TStateModel *model);

/* Utilities */
void ComputeDotProducts_All(T_VAR_Parameters *p);
void ComputeLogAbsDetA0_All(T_VAR_Parameters *p);
void ComputeLogAbsDetA0(int j, int k, T_VAR_Parameters *p);

TMatrix MakeA0(TMatrix A0, int k, T_VAR_Parameters *p);
TMatrix MakeAplus(TMatrix Aplus, int k, T_VAR_Parameters *p);
TMatrix MakeZeta(TMatrix Zeta, int k, T_VAR_Parameters *p);
TMatrix ConstructMatrixFromColumns(TMatrix X, TVector **, int k);

void UpdateStateDependentFields(T_VAR_Parameters *p, int *S);
void Update_aplus_from_bplus_a0(int j, int k, T_VAR_Parameters *p);
void Update_A0_from_b0(T_VAR_Parameters *p);
void Update_Aplus_from_bplus_A0(T_VAR_Parameters *p);
void Update_bplus_from_lambda_psi(T_VAR_Parameters *p);
void Update_b0_bplus_from_A0_Aplus(T_VAR_Parameters *p);
void Update_lambda_psi_from_bplus(T_VAR_Parameters *p);

int GetNumberStatesFromTranslationMatrix(int j, int **states);
int **CreateTranslationMatrix(TMarkovStateVariable ***list, TMarkovStateVariable *sv);

//PRECISION InnerProductSymmetric(TVector x, TMatrix S);
//PRECISION InnerProductNonSymmetric(TVector x, TVector y, TMatrix S);

void update_psi_quadratic_form(TMatrix S, int n, int m, int k, TVector lambda, TMatrix XX);
TMatrix MatrixInnerProductSymmetric(TMatrix X, TMatrix Y, TMatrix S);
PRECISION InnerProductSymmetric(TVector x, TMatrix S);
PRECISION InnerProductNonSymmetric(TVector x, TVector y, TMatrix S);
TVector DrawNormal_InverseVariance(TVector x, TVector b, TMatrix S);
TVector DrawNormal_InverseVariance_SVD(TVector x, TVector b, TMatrix S);
TVector DrawNormal_InverseUpperTriangular(TVector x, TVector b, TMatrix T);


// Obsolete routines


#endif  // __VAR_BASE_MODEL__


/********************************************************************************
Notes:

The model:

    y(t)' * A0(s(t)) = x(t)' * Aplus(s(t)) + epsilon(t)' * Inverse(Xi(s(t)))

   where
           y(t) is nvars x 1
           x(t) is npre x 1
             x(t)=[y(t-1),...,y(t-p),z(t)], where z(t) is exogenous
           epsilon(t) is nvars x 1
           A0(k) is nvars x nvars
           Aplus(k) is npre x nvars
           Xi(k) is an nvars x nvars diagonal matrix
           s(t) is an integer with 0 <= s(t) < nstates

   Furthermore

           A0(j,k) = U(j) * b0(j,k)

           Aplus(j,k) = V(j) * bplus(j,k) - W(j) * A0(j,k)

   and

           Zeta(j,k) = Xi(j,k)*Xi(j,k)
   where

           A0(j,k) is the jth column of A0(k)
           Aplus(j,k) is the jth column of A0(k)
           Xi(j,k) is the jth diagonal element of Xi(k)
           b0(j,k) is q(j) x 1
           bplus(j,k) is r(j) x 1
           U(j) is nvars x q(j) with orthonormal columns
           V(j) is npre x r(j) with orthonormal columns
           W(j) is npre x nvar
           e(j) is the jth column of an identity matrix

Sims-Zha Specification:
   This specification imposes that r(j) == npre, V(j) is the identity matrix, and
   W(j) is equal to a npre x nvars diagonal matrix with minus ones along the
   diagonal.  Further restrictions are imposed of the form.

           bplus(j,k) = f(psi(j),lambda(j,k))

   where

           psi(j) is npre x 1
           lambda(j,k) is nvars x 1

   and f is the function

           f(a,b) = diag(vec(a))

Random Walk Specification:
   This specification imposes that W(j) is equal to a npre x nvars diagonal
   matrix with minus ones along the diagonal.  Though it is not imposed, we
   usually want Aplus(j,k) to satisfy linear restrictions implicit in the
   matrix V(j).  This means that W(j) must be in the span of V(j) and hence
   (I - V(j)*V'(j))*W(j) = 0.


Normalization:
   We normalize by requiring Xi(0) to be the identity matrix and
   delta(j,0) to be a vector of ones.

Prior:
   A0(j,k) - The prior on A0(j,k) is normal with mean zero and covariance matrix
   A0_Prior(j).  This implies that the prior on b0(j,k) is normal with mean zero
   and covariance matrix Inverse(U'[j]*Inverse(A0_prior[j])*U[j]).

   Aplus(j,k) - The prior on Aplus(j,k) conditional on A0(j,k) is normal with mean
   -W(j) * A0(j,k) and covariance Aplus_Prior(j).  This implies that the prior on
   bplus(j,k) is normal with mean zero and covariance matrix
   Inverse(V'[j]*Inverse(Aplus_prior[j])*V[j])

   Zeta(j,k) - The prior on Zeta(j,k) is Gamma(zeta_a_prior(j),zeta_b_prior(j)).

---------------------------------------------------------------------------------

TVector** A0
  The length of A0 is nvars.  The vector A0[j][k] is the jth column of A0 when
  the jth coefficient state variable is equal to k.  Note that when the Markov
  state variable is equal to s, the jth coefficient state variable is equal to
  coef_states[j][s].  The number of distinct values for the jth coefficient state
  variable is equal to the dimension of A0[j].  This field is created with
  dw_CreateArray_array() and freed with dw_FreeArray().

TVector** b0
  The length of b0 is nvars.  The vector b0[j][k] consists of the free parameters
  in the jth column of A0 when the jth coefficient state variable is equal to k.
  The dimension of b0[j][k] does not vary across k.  Note that when the Markov
  state variable is equal to s, the jth coefficient state variable is equal to
  coef_states[j][s].  The dimension of b0[j] is equal to the dimension of A0[j].
  This field is created with dw_CreateArray_array() and freed with
  dw_FreeArray().

TVector** Aplus
  The length of Aplus is nvars.  The vector Aplus[j][k] is the jth column of
  Aplus when the jth coefficient state variable is equal to k.  Note that when
  the Markov state variable is equal to s, the jth coefficient state variable is
  equal to coef_states[j][s].  The dimension of Aplus[j] is equal to the
  dimension of A0[j].  This field is created with dw_CreateArray_array() and
  freed with dw_FreeArray().

TVector** bplus
  The length of bplus is nvars.  The vector bplus[j][k] consists of the free
  parameters in the jth column of Aplus when the jth coefficient state variable
  is equal to k.  The dimension of bplus[j][k] does not vary across k.  Note that
  when the Markov state variable is equal to s, the jth coefficient state
  variable is equal to coef_states[j][s].  The dimension of bplus[j] is equal to
  the dimension of A0[j].  This field is created with dw_CreateArray_array() and
  freed with dw_FreeArray().

PRECISION** Zeta
  The length of Zeta is nvars.  The value of Zeta[j][k] is the square of the
  value of the jth diagonal element of Xi when the jth variance state variable is
  equal to k.  Note that the the Markov state variable is equal to s, the jth
  variance state variable is equal to var_states[j][s].  The number of distinct
  values for the jth variance state variable is equal to the dimension of
  Zeta[j].  This field is created with dw_CreateArray_array() and freed with
  dw_FreeArray().

TVector** delta
  The length of delta is nvars.  The vector bplus[j][k] is a non-linear function
  of delta[j][k] and psi[k].  The length of delta[j][k] is nvars.  This field is
  non-null only when using the Sims-Zha specification.

TVector* psi
  The length of psi is nvars.  The vector bplus[j][k] is a non-linear function
  of psi[k] and delta[j][k].  The length of psi[k] is npre.  This field is non-
  null only when using the Sims-Zha specification.


=============================== State Translation ===============================

int*  n_var_states
  An integer array of dimension nvars.  The value of n_var_states[j] is the
  number of variance states for column j.

int** var_states
  An integer array of dimension nvars by nstates.  The value of var_states[j][k]
  is the value of the variance state for column j when the overall Markov state
  variable is equal to k.  It is used as an index into Xi[j].  It must be the
  case that

                       0 <= var_states[j][k] < n_var_states[j].

int*  n_coef_states
  An integer arrary of dimension nvars.  The value of n_coef_states[j] is the
  number of coefficient states for column j.

int** coef_states
  An integer array of dimension nvar by nstates.  The value of coef_states[j][k]
  is the value of the coefficient state for column j when the overall Markov
  state variable is equal to k.  It is used as an index into A0[j], b0[j],
  Aplus[j] or bplus[j].  It must be the case that

                      0 <= coef_states[j][k] < n_coef_states[j].

int n_A0_states
  The number of distinct values for the matrix A0.

int* A0_states
  An integer array of dimension nstates.  The value of A0_states[k] is the value
  of the state variable controlling A0 when the value of the overall Markov state
  variable is k.  It is used as an index into the vector log_abs_det_A0. It must
  be the case that

                          0 <= A0_states[k] < n_A0_states.

int** A0_column_states
  An integer array of dimension nvars by n_A0_states.  The value of
  A0_column_states[j][k] is the value of the coefficient state for column j when
  value of the state variable controlling the matrix A0 is k.  It is used as an
  index into A0[j].  It must be the case that

                   0 <= A0_column_states[j][k] < n_coef_states[j].


================================= Normalization =================================
For 0 <= k < n_A0_states, the contemporaneous coefficient matrix A[k] is formed.
For 0 <= j < nvars and 0 <= k < n_A0_states, the number

                e[j]*Inverse(A[k])*Target[j][A0_column_states[j][k]]

is computed.  If this number is negative, then the sign of

                              A0[j][A0_column_states[j][k]]

is flipped.  If the sign of any element of A0[j][.] is flipped more than once,
this event is recorded.

********************************************************************************/
