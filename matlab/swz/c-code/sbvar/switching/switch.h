

#ifndef __MARKOV_SWITCHING__
#define __MARKOV_SWITCHING__

#define __SWITCHING_VER_100__

#include "swzmatrix.h"

//=== Declaring structures so pointers can be defined ===
struct TStateModel_tag;
struct TParameters_tag;

/*******************************************************************************/
/**************************** TMarkovStateVariable *****************************/
/*******************************************************************************/

typedef struct TMarkovStateVariable_tag
{
  //=== Flags ===
  int valid_transition_matrix;

  //=== Sizes ===
  int nobs;
  int nstates;

  //=== State vector ===
  int* S;

  //=== Transition matrix
  TMatrix Q;

  //=== Quasi-free transition matrix parameters ===
  TVector *b;                   // The elements of b[k] are non-negative and their sum equals one up to DimV(b[k])*MACHINE_EPSILON.
  TVector  B;                   // b stacked into single vector

  //=== Prior information ===
  TMatrix  Prior;               // Dirichlet prior on the columns of Q.  Must be nstates x nstates with positive elements.
  TVector *Prior_b;             // Dirichlet prior on the quasi-free parameters b
  TVector  Prior_B;             // Prior_b stacked into single vector

  //=== Lag information encoding ===
  int   nlags_encoded;          // Number of lags encoded in the restrictions
  int   nbasestates;            // Number of base states nbasestates^(nlags_encoded) = nstates
  int** lag_index;              // nstates x (nlags_encoded + 1) lag_index[i][j] is the value of the jth lag when the overall state is k

  //=== Restrictions ===
  int* FreeDim;                 //
  int** NonZeroIndex;           // nstates x nstates
  TMatrix MQ;                   // nstates x nstates

  //=== Parent Markov state variable ===
  struct TMarkovStateVariable_tag *parent;           // either parent state variable or pointer to itself

  //=== Multiple state variables ===
  int n_state_variables;
  struct TMarkovStateVariable_tag **state_variable;
  TMatrix *QA;
  TVector *ba;
  TVector *Prior_ba;
  int** SA;
  int** Index;

  //=== Control variables ===
  int UseErgodic;

  //=== Workspace ===
  PRECISION LogPriorConstant;

} TMarkovStateVariable;


//=== Destructors ===
void FreeMarkovStateVariable(TMarkovStateVariable *sv);

//=== Constructors ===
TMarkovStateVariable* CreateMarkovStateVariable_Single(int nstates, int nobs, TMatrix Prior, int* FreeDim, int** NonZeroIndex, TMatrix MQ);
TMarkovStateVariable* CreateMarkovStateVariable_Multiple(int nobs, int n_state_variables, TMarkovStateVariable **state_variable);

TMarkovStateVariable* CreateMarkovStateVariable_Mixture(int nstates, int nobs, TMatrix Prior);
TMarkovStateVariable* CreateMarkovStateVariable_NoRestrictions(int nstates, int nobs, TMatrix Prior);
TMarkovStateVariable* CreateMarkovStateVariable_Exclusion(int nstates, int nobs, TMatrix Prior, TMatrix Exclusion);
TMarkovStateVariable* CreateMarkovStateVariable_SimpleRestrictions(int nstates, int nobs, TMatrix Prior, TMatrix* Restriction);
TMarkovStateVariable* CreateMarkovStateVariable_ConstantState(int nobs);
TMarkovStateVariable* DuplicateMarkovStateVariable(TMarkovStateVariable *sv);
TMarkovStateVariable* RestrictMarkovStateVariable(TMarkovStateVariable *sv, int nstates);

//=== Encoding lagged states into Markov state variable ===
TMarkovStateVariable* CreateMarkovStateVariable_Lags(int nlags, TMarkovStateVariable *base);
int** CreateLagIndex(int nbasestates, int nlags, int nstates);
TMatrix ConvertBaseTransitionMatrix(TMatrix T, TMatrix bT, int nlags);

//=== Data extractions routines ===
TMatrix GetTransitionMatrix_SV(TMatrix Q, TMarkovStateVariable *sv);
TMatrix GetBaseTransitionMatrix_SV(TMatrix Q, TMarkovStateVariable *sv);
#define GetTransitionProbability_SV(sv,j,i) (ElementM((sv)->Q,N_E2I[i],N_E2I[j]))
#define DecomposeIndexInd_SV(sv,i,j)  ((sv)->state_variable[j]->N_I2E[(sv)->Index[N_E2I[i]][j]])
#define DecomposeIndexLag_SV(sv,i,lag) ((sv)->baseN_I2E[(sv)->lag_index[N_E2I[i]][lag]])

//=== Normalization ===
void PropagateSwap_SV(TMarkovStateVariable *sv);
void Swap_SV(TMarkovStateVariable *sv, int i, int j);

//=== Prior routines ===
void SetLogPriorConstant_SV(TMarkovStateVariable *sv);
PRECISION LogPrior_SV(TMarkovStateVariable *sv);

//=== Simulation ===
void DrawTransitionMatrix_SV(TMarkovStateVariable *sv);
void DrawTransitionMatrixFromPrior_SV(TMarkovStateVariable *sv);
void SetTransitionMatrixToPriorMean_SV(TMarkovStateVariable *sv);
void DrawStatesFromTransitionMatrix_SV(TMarkovStateVariable *sv);

//=== Utility routines ===
void InvalidateTransitionMatrices_SV(TMarkovStateVariable *sv);
void ValidateTransitionMatrices_SV(TMarkovStateVariable *sv);
void PropagateStates_SV(TMarkovStateVariable *sv);
int PropagateTransitionMatrices_SV(TMarkovStateVariable *sv);
void Update_Q_from_B_SV(TMarkovStateVariable *sv);
int Update_B_from_Q_SV(TMarkovStateVariable *sv);
int TotalNumberStateVariables_SV(TMarkovStateVariable *sv);

int*  CreateStateIndex(TMarkovStateVariable* sv, TMarkovStateVariable** list, int n);
int** CreateTranslationMatrix(TMarkovStateVariable ***list, TMarkovStateVariable *sv);
int** CreateTranslationMatrix_Flat(int **states, TMarkovStateVariable *sv);
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/******************************** ThetaRoutines ********************************/
/*******************************************************************************/
typedef struct
{
  //=== Computes ln(P(y[t] | Y[t-1], Z[t], theta, s[t] = s)) ===
  PRECISION (*pLogConditionalLikelihood)(int s, int t, struct TStateModel_tag *model);

  //=== Computes E[y[t] | Y[t-1], Z[t], theta, Q, s[t]]
  TVector (*pExpectationSingleStep)(TVector y, int s, int t, struct TStateModel_tag *model);

  //=== Destructs parameters ===
  void (*pDestructor)(void *);

  //=== Draws parameters conditional states and transition probability ===
  void (*pDrawParameters)(struct TStateModel_tag *);

  //=== Computes Log of the prior on the model specific parameters ===
  PRECISION (*pLogPrior)(struct TStateModel_tag *);

  //=== Converts between free parameters and model specific parameters ===
  int (*pNumberFreeParametersTheta)(struct TStateModel_tag*);
  void (*pConvertFreeParametersToTheta)(struct TStateModel_tag*, PRECISION*);
  void (*pConvertThetaToFreeParameters)(struct TStateModel_tag*, PRECISION*);

  //=== Notification routines ===
  void (*pStatesChanged)(struct TStateModel_tag*);
  void (*pThetaChanged)(struct TStateModel_tag*);
  void (*pTransitionMatrixChanged)(struct TStateModel_tag*);
  int (*pValidTheta)(struct TStateModel_tag*);

  //=== Allows for initialization of data structures before forward recursion ===
  void (*pInitializeForwardRecursion)(struct TStateModel_tag*);

  //=== Permutes the elements of Theta.
  int (*pGetNormalization)(int*, struct TStateModel_tag*);
  int (*pPermuteTheta)(int*, struct TStateModel_tag*);

} ThetaRoutines;

//=== Constructors ===
ThetaRoutines* CreateThetaRoutines_empty(void);

/*******************************************************************************/
/********************************* TStateModel *********************************/
/*******************************************************************************/
typedef struct TStateModel_tag
{
  TMarkovStateVariable *sv;
  ThetaRoutines *routines;
  void *theta;

  //=== Control variables ===
  int ValidForwardRecursion;
  int UseLogFreeParametersQ;
  int NormalizeStates;

  //=== Common work space ===
  TVector* V;       // V[t][i] = P(s[t] = i | Y[t], Z[t], theta, Q)  0 <= t <= T and 0 <= i < nstates
  TVector* Z;       // Z[t][i] = P(s[t] = i | Y[t-1], Z[t-1], theta, Q)  0 < t <= T and 0 <= i < nstates
  PRECISION L;      // L = Sum(ln(Sum(P(y[t] | Y[t-1], Z[t], theta, s[t]) * P(s[t] | Y[t-1], Z[t-1], theta, Q),0 <= s[t] < nstates)),0 < t <= T)

  //=== Simulation status fields
  int n_degenerate_draws;    // counter for number of degenerate draws
  int *states_count;         // integer array of length nstates to count the number of each of the states

  //=== Obsolete fields retained for backward compatibility ===
  void *parameters;
  struct TParameters_tag *p;

} TStateModel;

//=== Destructors ===
void FreeStateModel(TStateModel *model);

//=== Constructors ===
TStateModel* CreateStateModel_new(TMarkovStateVariable *sv, ThetaRoutines *routines, void *theta);

//=== Notification routines ===
void StatesChanged(TStateModel *model);
void TransitionMatricesChanged(TStateModel *model);
void ThetaChanged(TStateModel *model);
#define ValidTheta(model)             (((model)->routines->pValidTheta) ? (model)->routines->pValidTheta(model) : 1)
#define ValidTransitionMatrix(model)  ((model)->sv->valid_transition_matrix)

//=== Simulation routines ===
void ForwardRecursion(TStateModel *model);
void DrawStates(TStateModel *model);
void DrawStatesFromTransitionMatrix(TStateModel *model);
void DrawTransitionMatrix(TStateModel *model);
void DrawTransitionMatrixFromPrior(TStateModel *model);
void SetTransitionMatrixToPriorMean(TStateModel *model);
void DrawTheta(TStateModel *model);
void DrawAll(TStateModel *model);

//=== Normalization ===
int SetStateNormalizationMethod(int (*pGetNormalization)(int*, struct TStateModel_tag*),int (*pPermuteTheta)(int*, struct TStateModel_tag*),TStateModel *model);
int NormalizeStates(TStateModel *model);

//=== Probability routines ===
// ln(P(y[t] | Y[t-1], Z[t], theta, Q, s[t]))
#define LogConditionalLikelihood(s,t,model)  ((model)->routines->pLogConditionalLikelihood(s,t,model))

// ln(P(y[t] | Y[t-1], Z[t], theta, Q))
PRECISION LogConditionalLikelihood_StatesIntegratedOut(int t, TStateModel *model);

// E[y[t] | Y[t-1], Z[t], theta, Q, s[t]]
#define ExpectationSingleStep(y,s,t,model) ((model)->routines->pExpectationSingleStep(y,s,t,model))

// E[y[t] | Y[t-1], Z[t], theta, Q]
TVector ExpectationSingleStep_StatesIntegratedOut(TVector y, int t, TStateModel *model);

// ln(P(Q))
#define LogPrior_Q(model)  (LogPrior_SV((model)->sv))

// ln(P(theta))
#define LogPrior_Theta(model)  ((model)->routines->pLogPrior(model))

// ln(P(theta, Q))
#define LogPrior(model)  (LogPrior_Theta(model) + LogPrior_Q(model))

// ln(P(S[T] | theta, Q))
PRECISION LogConditionalPrior_S(TStateModel *model);

// ln(P(Y[T] | Z[T], theta, Q, S[T]))
PRECISION LogLikelihood(TStateModel *model);

// ln(P(Y[T] | Z[T], theta, Q, S[T]) * P(S[T] | Theta, Q) * P(Theta, Q))
#define LogPosterior(model)  (LogLikelihood(model) + LogConditionalPrior_S(model) + LogPrior(model))

// ln(P(Y[T] | Z[T], theta, Q))
PRECISION LogLikelihood_StatesIntegratedOut(TStateModel *model);

// ln(P(Y[T] | Z[T], theta, Q) * P(Theta, Q))
#define LogPosterior_StatesIntegratedOut(model)  (LogLikelihood_StatesIntegratedOut(model) + LogPrior(model))

// ln(P(S[T] | Y[T], Z[T], theta, Q))
PRECISION LogConditionalProbabilityStates(int *S, TStateModel *model);

// P(s[t] | Y[t], Z[t], theta, Q)
PRECISION ProbabilityStateConditionalCurrent(int s, int t, TStateModel *model);

// P(s[t] | Y[t-1], Z[t-1], theta, Q)
PRECISION ProbabilityStateConditionalPrevious(int s, int t, TStateModel *model);

// P[t] = P(s[t] = s | Y[T], Z[T], theta, Q)
TVector ProbabilitiesState(TVector P, int s, TStateModel *model);

//=== Free parameters routines ===
int NumberFreeParametersQ(TStateModel *model);
void ConvertQToFreeParameters(TStateModel *model, PRECISION *f);                                                // needs to be modified
void ConvertFreeParametersToQ(TStateModel *model, PRECISION *f);                                                // needs to be modified
void ConvertQToLogFreeParameters(TStateModel *model, PRECISION *f);                                             // needs to be modified
void ConvertLogFreeParametersToQ(TStateModel *model, PRECISION *f);                                             // needs to be modified
#define NumberFreeParametersTheta(model)  ((model)->routines->pNumberFreeParametersTheta(model))
void ConvertFreeParametersToTheta(TStateModel *model, PRECISION *f);
#define ConvertThetaToFreeParameters(model,f)  ((model)->routines->pConvertThetaToFreeParameters(model,f))

//=== Setup integrity routines ===
int CheckRestrictions(int* FreeDim, int** NonZeroIndex, TMatrix MP, int nstates);
int CheckPrior(TMatrix Prior, int nstates);
int CheckPriorOnFreeParameters(TMatrix Prior, int** NonZeroIndex, int nstates);

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

//=== Utility routines ===
TVector Ergodic(TVector v, TMatrix P);
TVector Ergodic_SVD(TVector v, TMatrix P);
TVector* ErgodicAll_SVD(TMatrix P);
TVector DrawDirichletVector(TVector Q, TVector Alpha);
TVector* DrawIndependentDirichletVector(TVector *Q, TVector *A);
PRECISION LogDirichlet_pdf(TVector Q, TVector Alpha);
PRECISION LogIndependentDirichlet_pdf(TVector *Q, TVector *Alpha);
int DrawDiscrete(TVector p);
PRECISION AddLogs(PRECISION a, PRECISION b);
PRECISION AddScaledLogs(PRECISION x, PRECISION a, PRECISION y, PRECISION b);

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/****** Obsolete names and structures retained for backward compatibility ******/
/*******************************************************************************/
typedef struct TParameters_tag
{
  //=== Computes ln(P(y[t] | Y[t-1], Z[t], theta, s[t] = s)) ===
  PRECISION (*pLogConditionalLikelihood)(int s, int t, struct TStateModel_tag *model);

  //=== Destructs parameters ===
  void (*pParameterDestructor)(void *parameters);

  //=== Draws parameters conditional states and transition probability ===
  void (*pDrawParameters)(struct TStateModel_tag *);

  //=== Computes Log of the prior on the model specific parameters ===
  PRECISION (*pLogPrior)(struct TStateModel_tag *);

  //=== Converts between free parameters and model specific parameters ===
  int (*pNumberFreeParametersTheta)(struct TStateModel_tag*);
  void (*pConvertFreeParametersToTheta)(struct TStateModel_tag*, PRECISION*);
  void (*pConvertThetaToFreeParameters)(struct TStateModel_tag*, PRECISION*);

  // Obsolete fields retained for backward compatibility
  void *p;

} TParameters;

//=== Constructors/Destructors ===
void FreeParameters(TParameters *p);

TParameters* CreateParameters(PRECISION (*)(int,int,struct TStateModel_tag*),    // pLogConditionalLikelihood
                              void (*)(void*),                                   // Destructor for parameters
                  PRECISION (*)(struct TStateModel_tag*),            // pLogPrior
                  int (*)(struct TStateModel_tag*),                  // pNumberFreeModelSpecificParameters
                  void (*)(struct TStateModel_tag*, PRECISION*),     // pConvertFreeParametersToModelSpecificParameters
                  void (*)(struct TStateModel_tag*, PRECISION*),     // pConvertModelSpecificParametersToFreeParameters
                              void (*)(struct TStateModel_tag*),                 // pDrawParameters
                              void *);                                           // pointer to user defined parameters

TStateModel* CreateStateModel(TMarkovStateVariable *sv, TParameters *p);

//=== Obsolete names ===
#define ProbabilityStateGivenCurrentData(s,t,model) ProbabilityStateConditionalCurrent(s,t,model)
#define ProbabilityStateGivenPreviousData(s,t,model) ProbabilityStateConditionalPrevious(s,t,model)
#define ProbabilityStateGivenAllData(P,s,model) ProbabilitiesState(P,s,model)
//PRECISION ProbabilityStatesGivenData(TStateModel *model);
#define LogLikelihoodGivenParameters(model)  LogLikelihood_StatesIntegratedOut(model);
//PRECISION LogMarginalPosterior(TStateModel *model);
//void DrawAllParameters(TMarkovStateVariable *sv);
//void DrawAllParametersAndNormalizeStates(TMarkovStateVariable *sv);
//PRECISION ComputeMarginalLogLikelihood(TMarkovStateVariable *sv);
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

#endif

/***************************** TMarkovStateVariable *****************************
  The TMarkovStateVariable type can represent either a single Markov state
  variable or a collection of independent Markov state variables.

  The transition matrix Q is generated for a single Markov state variable via
  the routines DrawTransitionMatrixFromPrior_SV() or DrawTransitionMatrix_SV().
  Calls to these functions by multiple independent Markov state variables result
  in recursive call to these functions.

  The vector of states S is generated only by a TStateModel type containing the
  TMarkovStateVariable type.  The state is only generated at the top level and
  can be propagated downward with a call to PropagateStates_SV().

  The following set of fields are set for both types.
  ===============================================================================
    int UseErgodic
      Uses the ergodic distribution if non-zero and use the uniform distribution
      otherwise.

    int nstates
      Number of states.  Always positive.

    int nobs
      Number of observations.  Always positive.

    int* S
      S[t] is the state at time t, for 0 <= t <= nobs.  S is created via a call
      to dw_CreateArray_int().  It is guaranteed that 0 <= S[t] < nstates.

    TMatrix Q
      Markov transition matrix.

    struct TMarkovStateVariable_tag *parent
      Parent of the Markov state variable.  If the Markov state variable has no
      parent, then parent is a pointer to the structure itself.

    int n_state_variables
      Number of state variables.  Will be equal to one for single Markov state
      variables and larger than one in the case of multiple independent Markov
      state variables

    struct TMarkovStateVariable_tag **state_variable
      An array of markov state variables of length n_state_variables.  If
      n_state_variables is equal to one, then state_variable[0] is a pointer to
      the structure itself.  Care must be taken to ensure that infinite loops do
      not result when transversing through state variables.  When creating a
      mulitple Markov state variable via a call to the routine
      CreateMarkovStateVariable_Multiple(), the last argument which is a pointer
      to a pointer to a TMarkovStateVariable must have been created with

        dw_CreateArray_pointer(n_state_variables,(void (*)(void*))FreeMarkovStateVariable);

      Furthermore, the structure receives ownership of this argument and is
      responsible for its freeing.

    int** Index
      This is a nstates x n_state_variables rectangular array of integers.  State
      s corresponds to state_variable[i] being equal to Index[s][i].

    int** SA
      An array of integer pointers of length n_state_variables.  The pointers SA[i]
      and state_variable[i]->S point to the same address.

    TMatrix* QA
      An array of matrices of length n_state_variables.  The matrix QA[i] is
      is the matrix state_variable[i]->Q.

    TVector* ba
      For single Markov state variables, ba[i] = b[i].  For multiple state
      variables ba[i] = state_variable[k]->ba[j] where

        i = j + dw_DimA(state_variable[0]->ba) + ... + dw_DimA(state_variable[k-1]->ba)

    TVector* Prior_ba
      For single Markov state variables, Prior_ba[i] = Prior_b[i].  For multiple
      state variables Prior_ba[i] = state_variable[k]->Prior_ba[j] where

        i = j + dw_DimA(state_variable[0]->Prior_ba) + ... + dw_DimA(state_variable[k-1]->Prior_ba)

  ===============================================================================
  The following fields are set only for single Markov state variables and are
  set to null for multiple independent Markov state variables.

    TVector B
      The vector B is the vector of quasi-free parameters.

    TVector *b
      Array of vectors of length DimA(FreeDim).  The element b[k] is of length
      FreeDim[k].  Non-standard memory management is used so that

               &(b[k][i])=&B[FreeDim[0] + ... + FreeDim[k-1] + i])

      The elements of b[k] are non-negative and their sum equals one up to
      DimV(b[k])*MACHINE_EPSILON.

    TMatrix Prior
       Prior Dirichlet parameters for Q.

    TVector *Prior_b
      The Dirichlet prior parametrs for b.  Array of vectors of length
      DimA(FreeDim).  The element Prior_b[k] is of length FreeDim[k].
      Non-standard memory management is used so that

            &(Prior_b[k][i])=&B[FreeDim[0] + ... + FreeDim[k-1] + i])

    TVector Prior_B
      The Dirichlet prior parameters for B.  This vector is created and
      initialized by CreateMarkovStateVariable().  The element B[k]-1 is the sum
      of Prior[i][j]-1 over all (i,j) such that NonZeroIndex[i][j] == k.

    int* FreeDim
      FreeDim[k] is the length of the kth free Dirichlet vector.  The length of B
      must be equal to FreeDim[0] + ... + FreeDim[dw_DimA(FreeDim)-1].

    int** NonZeroIndex
      Defines the relationship between Q and B.
                  --
                 | MQ[i][j]*B[NonZeroIndex[i][j]]  if NonZeroIndex[i][j] >=  0
       Q[i][j] = |
                 | 0.0                             if NonZeroIndex[i][j] == -1
                  --
    TMatrix MQ
      Coefficients for the elements of Q in terms of the free parameters B.


  int   nlags_encoded;          // Number of lags encoded in the restrictions
  int   nbasestates;            // Number of base states nbasestates^(nlags_encoded) = nstates
  int** lag_index;              // nstates x (nlags_encoded + 1) lag_index[i][j] is the value of the jth lag when the overall state is k


  ===============================================================================


  ===============================================================================
  Normalization:
   In general, a permutation of the states, together with the corresponding
   permutation of the rows and columns of the transition matrix and the model
   dependent parameters theta, does not change the value of the likelihood
   function and so presents a normalization problem.  However, some permutations
   of the states are not permissible in that they may violate the restrictions
   placed on the transition matrix or restrictions on the model dependent
   parameters.  Furthermore, even if a permutation did not cause a violation of
   any restrictions, a non-symmetric prior on the model dependent parameters
   could cause the value of the likelihood to change.

********************************************************************************/

/********************************************************************************
   Constrained optimization
    Because one of the tasks of this suite of programs is to find the maximum
    likelihood estimate, constrained optimization routines are often called.
    These routines almost always require a single or double precision array
    containing the parameters to be optimized over to be passed.  To facilitate
    this, the TMarkovStateVariable type allows the memory being allocated for
    parameters to be passed, which allows the user to control exactly how the
    memory is laid out.  Also, the TParameter type is defined which automates
    some of this process.  The goal is to allow different aspects of this suite
    to access the common memory area for parameters in a manner appropriate to
    the task at hand while minimizing the amount of memory copying.

   Memory Management
    Different parts of this suite of routines need to have the parameters
    available in certain formats.  However, copying parameters from one format to
    another should be avoided.  Our approach is as follows.  We define a
    structure TParameters which contains the following:

      int n_real_parameters;
      int n_integer_parameters;
      int n_miscellaneous_parameters;
      void* v;
      void* extra;
        The pointer v points to contiguous memory that contains all the
        parameters.  The first n_real_parameters*sizeof(PRECISION) bytes contains
        floating point parameters.  The next n_integer_parameters*sizeof(int)
        bytes contain integer parameters.  The next n_miscellaneous_parameters
        bytes contain whatever parameters do not fit into the first two
        catagories.

        The idea is that many maximization routines work on a vector of floating
        point numbers.  The integer parameters are included to record the values
        of Markov state variables, though they could certainly be used for other
        purposes.  The general purpose memory is to provide flexibility.

        The pointer extra will point to some structure needed for the
        implementation of model in each state.

   Markov state variable parameter pointers
    The Markov state variables require transition matrix parameters and the
    values of the Markov state variables.  For this two structures are defined.

      int n_state_variables;
      int T;
      TVectorArray BA;
      TMatrixArray QA;
      TIntVectorArray SA;
        All of pElementV(BA[i]), pElementM(QA[i]), and pElementIV(SA[i]) point to
        memory positions in the vector v.  Non-standard memory management
        techniques are used and the following must be true.

          All of pElementV(BA[i]), pElementM(QA[i]), and pElementIV(SA[i]) were
          allocated with malloc(), can be freed with free(), and none of
          FreeVector(BA[i]), FreeMatrix(QA[i]), or FreeIntMatrix(SA[i]) attempt
          to free pElemementV(BA[i]), pElementM(QA[i]), or pElementIV(SA[i]) if
          these are null pointers.

        The default behavior is for the memory allocation of v is:

                           Theta
                           BA[0]
                             .
                             .
                   BA[n_state_variables-1]
                           QA[0]
                             .
                             .
                   QA[n_state_variables-1]
                           SA[0]
                             .
                             .
                   SA[n_state_variables-1]

        The vector BA[i] stores the free parameters for the transition matrix
        of the ith Markov state variable.  The matrix QA[i] stores the full
        transition matrix for the ith Markov state variable.  The integer array
        SA[i] stores the values of the ith Markov state variable.  The number of
        Markov state variables is n_state_variables and there are T+1 values
        stored for each Markov state variable.  The total number of states, which
        is the product of the number of states for each state variables is
        nstates.  If there are no restrictions on the ith transition matrix other
        than its elements being non-negative and the sum of its columns being
        one, then BA[i] can be equal to QA[i].


*/

