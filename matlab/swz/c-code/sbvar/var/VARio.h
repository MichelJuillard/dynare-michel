
#ifndef __VAR_INPUT_OUTPUT__
#define __VAR_INPUT_OUTPUT__

#include "switch.h"
#include "VARbase.h"

void Write_VAR_Specification(FILE *f, char *filename, TStateModel *model);
TStateModel* Read_VAR_Specification(FILE *f, char *filename);

int Write_VAR_Parameters(FILE *f, char *filename, char *id, TStateModel *model);
int Read_VAR_Parameters(FILE *f, char *filename, char *id, TStateModel *model);
int Write_VAR_ParametersFlat(FILE *f, TStateModel *model, char *fmt);
int Write_VAR_ParametersFlat_Headers(FILE *f_out, TStateModel *model);
int Write_VAR_ParametersFlat_A0_Diagonal_One(FILE *f, TStateModel *model, char *fmt);

void ReadAllParameters(FILE *f, char *filename, char *id, TStateModel *model);
void WriteAllParameters(FILE *f, char *filename, char *id, TStateModel *model);

//T_VAR_Parameters* Create_VAR_Parameters_File(FILE *f, char *filename, TMarkovStateVariable *sv);
//TStateModel* CreateStateModel_VAR_File(FILE *f, char *filename);

//void PrintParametersVAR(FILE *f_out, TStateModel *model);
//void Write_VAR_Info(FILE *f, char *filename, T_VAR_Parameters *p);

#endif
