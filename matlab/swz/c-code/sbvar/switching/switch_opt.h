

#ifndef __MARKOV_SWITCHING_OPTIMIZATION__
#define __MARKOV_SWITCHING_OPTIMIZATION__

#include "switch.h"

void SetupObjectiveFunction(TStateModel *model, PRECISION *MFPparms, PRECISION *FTMparms, PRECISION *FMSparms);

PRECISION PosteriorObjectiveFunction(PRECISION *x, int n);
PRECISION PosteriorObjectiveFunction_csminwel(double *x, int n, double **args, int *dims);
void PosteriorObjectiveFunction_npsol(int *mode, int *n, double *x, double *f, double *g, int *nstate);

PRECISION MLEObjectiveFunction(PRECISION *x, int n);
PRECISION MLEObjectiveFunction_csminwel(double *x, int n, double **args, int *dims);
void MLEObjectiveFunction_npsol(int *mode, int *n, double *x, double *f, double *g, int *nstate);

PRECISION MLEObjectiveFunction_LogQ(PRECISION *x, int n);
PRECISION MLEObjectiveFunction_LogQ_csminwel(double *x, int n, double **args, int *dims);

#endif
