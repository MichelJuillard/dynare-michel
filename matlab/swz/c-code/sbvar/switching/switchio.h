
#include "switch.h"

/*
   Base routines for reading/writing Markov state variables and transition
   matrices in native ascii format.
*/
TMarkovStateVariable* ReadMarkovSpecification_SV(FILE *f_in, char *idstring, int nobs);
int WriteMarkovSpecification_SV(FILE *f_out, TMarkovStateVariable *sv, char *idstring);
int ReadTransitionMatrices_SV(FILE *f_in, TMarkovStateVariable* sv, char *header, char *idstring);
int WriteTransitionMatrices_SV(FILE *f_out, TMarkovStateVariable* sv, char *header, char *idstring);
int ReadBaseTransitionMatrices_SV(FILE *f_out, TMarkovStateVariable *sv, char *header, char *idstring);
int WriteBaseTransitionMatrices_SV(FILE *f_out, TMarkovStateVariable *sv, char *header, char *idstring);

int WriteBaseTransitionMatricesFlat_SV(FILE *f_out, TMarkovStateVariable *sv, char *fmt);
void WriteBaseTransitionMatricesFlat_Headers_SV(FILE *f_out, TMarkovStateVariable* sv, char *idstring);

/*
   Routines for reading/writing Markov state variables and transition matrices
   from TStateModel.  Calls base routines.
*/
TMarkovStateVariable* ReadMarkovSpecification(FILE *f, char *filename);
int WriteMarkovSpecification(FILE *f, char *filename, TStateModel *model);
int ReadTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model);
int WriteTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model);
int ReadStates(FILE *f, char *filename, char *header, TStateModel *model);
int WriteStates(FILE *f, char *filename, char *header, TStateModel *model);
int ReadBaseTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model);
int WriteBaseTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model);

int WriteBaseTransitionMatricesFlat(FILE *f, TStateModel *model, char *fmt);

/*
   Read flat markov state variable specification from file.
*/
TMarkovStateVariable* CreateMarkovStateVariable_File(FILE *f, char *filename, int nobs);
