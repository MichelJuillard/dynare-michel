/*
 * Copyright (C) 2007-2009 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MEX_INTERFACE_HH_INCLUDED
#define MEX_INTERFACE_HH_INCLUDED
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <map>
#include <vector>
#include <algorithm>
using namespace std;


typedef unsigned int mwIndex;
typedef unsigned int mwSize;

enum mxData_type
  {
    mxREAL,
    mxCOMPLEX
  };

enum mxArray_type
  {
    mxDOUBLE_CLASS = 0,
    mxSINGLE_CLASS = 1,
    mxLOGICAL_CLASS = 2,
    mxCHAR_CLASS = 3,
    mxSPARSE_CLASS = 4,
    mxSTRUCT_CLASS = 5,
    mxCELL_CLASS = 6,
    mxFUNCTION_CLASS = 7,
    mxUINT8_CLASS = 8,
    mxINT16_CLASS = 9,
    mxUINT16_CLASS = 10,
    mxINT32_CLASS = 11,
    mxUINT32_CLASS = 12,
    mxINT64_CLASS = 13,
    mxUINT64_CLASS = 14
  };

class mxArray_tag {
  public:
  unsigned int size_1, size_2;
  mwIndex *Ir, *Jc;
  int Nzmax;
  double *data;
  mxArray_type type;
  vector<string> field_name;
  vector<vector<mxArray_tag*> > field_array;
  mxArray_tag();
};



typedef mxArray_tag mxArray;

void load_global(char* file_name);
/*Matlab clone function*/
int mexPrintf(const char *str, ...);
void mexErrMsgTxt(const string str);
void mexWarnMsgTxt(const string str);
void* mxMalloc(unsigned int amount);
void* mxRealloc(void *to_extend, unsigned int amount);
void* mxCalloc(unsigned int nb_elements, unsigned int amount_per_element);
void mxFree(void *to_release);
void mexEvalString(const string str);
double* mxGetPr(const mxArray *b_m);
inline mwIndex* mxGetIr(const mxArray *A_m) {return (mwIndex*)A_m->Ir;};
inline mwIndex* mxGetJc(const mxArray *A_m) {return (mwIndex*)A_m->Jc;};
inline int mxGetNzmax(const mxArray *A_m) {return A_m->Nzmax;};
inline int mxGetM(const mxArray *A_m) {return A_m->size_1;};
inline int mxGetN(const mxArray *A_m) {return A_m->size_2;};
inline bool mxIsChar(const mxArray *A) {return A->type == mxCHAR_CLASS;};
inline bool mxIsFloat(const mxArray *A) {return A->type == mxDOUBLE_CLASS || A->type == mxSINGLE_CLASS;};
inline bool mxIsStruct(const mxArray *A) {return A->type == mxSTRUCT_CLASS;};
mxArray* mxCreateDoubleMatrix(unsigned int rows,unsigned int cols, mxData_type mx_type);
mxArray* mxCreateSparse(unsigned int rows, unsigned int cols, unsigned int nz_max, mxData_type mx_type);
mxArray* mxCreateCharArray(unsigned int rows, unsigned int cols, mxData_type mx_type);
mxArray* mxCreateDoubleScalar(double value);
mxArray* mxCreateStructMatrix(unsigned int rows, unsigned int cols, unsigned int nfields, const string &fieldnames);
inline mxArray* mxCreateStructArray(unsigned int rows, mwSize* cols, int nfields, const string &fieldnames) {return mxCreateStructMatrix(rows, *cols, nfields, fieldnames);};
mxArray* mxCreatNULLMatrix();
void mexCallMATLAB(unsigned int n_lhs, mxArray* lhs[], unsigned int n_rhs, mxArray* rhs[], const char* function);
void mxDestroyArray(mxArray* A_m);
mxArray* read_struct(FILE *fid);
mxArray* read_Array(FILE *fid);
mxArray* read_double_array(FILE *fid);
mxArray* mexGetVariable(const char* space_name, const char* matrix_name);
int mxGetFieldNumber(const mxArray* Struct, const char* field_name);
mxArray* mxGetFieldByNumber(mxArray* Struct, unsigned int pos, unsigned int field_number);
void mxSetFieldByNumber(mxArray *Struct, mwIndex index, unsigned int field_number, mxArray *pvalue);
int mxAddField(mxArray* Struct, const char* field_name);
void mxSetFieldByNumber(mxArray *Struct, mwIndex index, unsigned int fieldnumber, mxArray *pvalue);
inline unsigned int mxGetNumberOfElements(const mxArray* A) {return A->size_1 * A->size_2;};
int mxGetString(const mxArray *array, char *buf, unsigned int buflen);
inline double mxGetScalar(const mxArray *array) {if (!mxIsFloat(array)) mexErrMsgTxt("not a float array\n");return array->data[0];};
mxArray* mxDuplicateArray(const mxArray *array);
void Free_simple_array(mxArray* array);
void Free_struct_array(mxArray* array);
void Free_Array(mxArray* array);
void Free_global();
#endif
