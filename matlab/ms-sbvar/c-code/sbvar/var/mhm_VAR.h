
#ifndef __MODIFIED_HARMONIC_MEAN_
#define __MODIFIED_HARMONIC_MEAN_

#include "swzmatrix.h"
#include "switch.h"
#include "VARbase.h"

typedef struct
{
/*    // Sample sizes to use in computations   ansi-c*/
  int n_burn1;                /*   set to negative of actual value if first burn-in has been completed   ansi-c*/
  int n_burn2;                /*   set to negative of actual value if second burn-in has been completed   ansi-c*/
  int n_mean_variance;        /*   set to negative of actual value if mean and variance have been computed   ansi-c*/
  int n_mhm;                  /*   number of draws for computing the modified harmonic mean   ansi-c*/
  int n_thin;                 /*   thinning factor   ansi-c*/

/*    // Accumulation fields   ansi-c*/
  int N;                      /*   Total number observations   ansi-c*/

  PRECISION old_log_posterior;
  PRECISION sum;
  PRECISION sum_square;
  PRECISION max_log_posterior;
  PRECISION max_log_likelihood;

/*    // mhm info   ansi-c*/
  TVector mean;               /*   Gaussian mean   ansi-c*/
  TMatrix variance;           /*   Gaussian variance   ansi-c*/
  TMatrix inverse_variance;   /*   Inverse of Gaussian variance   ansi-c*/
  TVector center;             /*   Used to center gaussian.  Must be equal to posterior_mode_VAR or mean.   ansi-c*/

  TVector   alpha_scales;     /*   scaling values for base dirichlet pdf parameters   ansi-c*/
  TVector*  BaseAlpha;        /*   base dirichlet pdf parameters   ansi-c*/
  TVector** Alpha;            /*   base dirichlet pdf parameters times the scale factors   ansi-c*/

/*    // Model info   ansi-c*/
  TStateModel *model;
  TVector posterior_mode_VAR;
  PRECISION log_likelihood_at_mode;
  PRECISION log_posterior_at_mode;
  PRECISION log_prior_at_mode;

/*    // Workspace   ansi-c*/
  TVector free_parameters_VAR;       /*   workspace for free parameters for VAR   ansi-c*/

/*   // Workspace for states   ansi-c*/
  int *states;

/*    // Files   ansi-c*/
  FILE *f_out;
  FILE *f_out_regime_counts;
  char *regime_counts_filename;
  char *intermediate_output_filename;
  char *final_output_filename;
  char *intermediate_draws_output_filename;
  char *draws_output_filename;
  char *spec_filename;
  char *parameter_filename;
  char *parameter_header;
  char *mhm_filename;

} T_MHM;

/*  // Constructors   ansi-c*/
void FreeMHM(T_MHM *mhm);
T_MHM* CreateMHM(void);
T_MHM* AddDirichletScales(TVector alpha_scales, T_MHM *mhm);
T_MHM* AddStateModel(TStateModel *model, T_MHM *mhm);

void ResetMHM(T_MHM *mhm);
void BurnIn(T_MHM *mhm, int iterations, int period);
void BurnIn_AdaptiveMetropolisScale(T_MHM *mhm, int iterations, int period);
void ComputeMeanVariance_MHM(T_MHM *mhm, int iterations, int period);
int IsValidZeta(PRECISION* zeta, int n, PRECISION* gamma_cuts);
PRECISION UpdatePosteriorLikelihood(T_MHM *mhm);
void UpdateModifiedHarmonicMean(T_MHM *mhm, int n_singular);
void ComputeModifiedHarmonicMean(T_MHM *mhm, int period);

void WriteMHM_Input(FILE *f_out, T_MHM *mhm);
T_MHM* ReadMHM_Input(FILE *f_in, char *filename, T_MHM *mhm);
void WriteMeanVariance(FILE *f_out, T_MHM *mhm);
int ReadMeanVariance(FILE *f_in, T_MHM *mhm);

void PrintJumps(FILE *f, T_VAR_Parameters *p);
void PrintMHM(FILE *f, char *filename, TStateModel *model, T_MHM *mhm);

#endif
