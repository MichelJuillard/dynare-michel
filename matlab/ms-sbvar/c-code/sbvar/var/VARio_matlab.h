
#ifndef __VAR_INPUT_OUTPUT_MATLAB__
#define __VAR_INPUT_OUTPUT_MATLAB__

#include "switch.h"
#include "VARbase.h"

TStateModel* Combine_matlab_standard(char *inmatlab, char *instandard);
void ReadConstantParameters(char *filename, TStateModel *model);
TStateModel* CreateStateModel_VAR_matlab(char *filename);

#endif
