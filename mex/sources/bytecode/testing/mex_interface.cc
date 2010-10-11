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

#include "mex_interface.hh"
#include <cstring>
#include <sstream>

using namespace std;

map<string, mxArray*> mxglobal;

mxArray_tag::mxArray_tag()
{
  type = mxDOUBLE_CLASS;
  Ir = NULL;
  Jc = NULL;
  Nzmax = 0;
  field_name.clear();
  field_array.clear();
  size_1 = 0;
  size_2 = 0;
  data = NULL;
}

int
mexPrintf(const char *str, ...)
{
  va_list args;
  int retval;

  va_start(args, str);
  retval = vprintf(str, args);
  va_end(args);

  return retval;
}

void
mexErrMsgTxt(const string str)
{
  perror(str.c_str());
  exit(EXIT_FAILURE);
}

void
mexWarnMsgTxt(const string str)
{
  perror(str.c_str());
}

void
mxFree(void *to_release)
{
  free(to_release);
}

void *
mxMalloc(unsigned int amount)
{
  return malloc(amount);
}

void *
mxRealloc(void *to_extend, unsigned int amount)
{
  return realloc(to_extend, amount);
}

void*
mxCalloc(unsigned int nb_elements, unsigned int amount_per_element)
{
  return calloc(nb_elements, amount_per_element);
}

void
mexEvalString(const string str)
{
}

double*
mxGetPr(const mxArray *b_m)
{
  return b_m->data;
}

mxArray*
mxCreateDoubleMatrix(unsigned int rows,unsigned int cols, mxData_type mx_type)
{
  mxArray *Array = new mxArray;
  Array->type = mxDOUBLE_CLASS;
  Array->size_1 = rows;
  Array->size_2 = cols;
  Array->data = (double*)mxMalloc(rows*cols*sizeof(double));
  return(Array);
}

mxArray*
mxCreateCharArray(unsigned int rows,unsigned int cols, mxData_type mx_type)
{
  mxArray *Array = new mxArray;
  Array->type = mxCHAR_CLASS;
  Array->size_1 = rows;
  Array->size_2 = cols;
  Array->data = (double*)mxMalloc(2*rows*cols*sizeof(char));
  return(Array);
}

mxArray*
mxCreateSparse(unsigned int rows, unsigned int cols, unsigned int nz_max, mxData_type mx_type)
{
  mxArray *Array = new mxArray;
  Array->type = mxSPARSE_CLASS;
  Array->size_1 = rows;
  Array->size_2 = cols;
  Array->Nzmax = nz_max;
  Array->data = (double*)mxMalloc(nz_max*sizeof(double));
  Array->Ir = (mwIndex*)mxMalloc(nz_max*sizeof(mwIndex));
  Array->Jc = (mwIndex*)mxMalloc(cols*sizeof(mwIndex));
  return(Array);
}

mxArray*
mxCreateDoubleScalar(double value)
{
  mxArray *Array = new mxArray;
  Array->type = mxSINGLE_CLASS;
  Array->size_1 = 1;
  Array->size_2 = 1;
  Array->data = (double*)mxMalloc(sizeof(double));
  Array->data[0] = value;
  return(Array);
}

mxArray*
mxCreatNULLMatrix()
{
  return NULL;
}

mxArray*
mxCreateStructMatrix(unsigned int rows, unsigned int cols, unsigned int nfields, const vector<string> &fieldnames)
{
  mxArray *Array = new mxArray;
  Array->type = mxSTRUCT_CLASS;
  Array->size_1 = rows;
  Array->size_2 = cols;
  for (unsigned int i = 0; i < nfields; i++)
    Array->field_name.push_back(fieldnames[i]);
  return Array;
}

void
mxDestroyArray(mxArray* A_m)
{
  mxFree(A_m->data);
  if (A_m->Ir)
    mxFree(A_m->Ir);
  if (A_m->Jc)
    mxFree(A_m->Jc);
  mxFree(A_m);
}

void
mexCallMATLAB(unsigned int n_lhs, mxArray* matrix_lhs[], unsigned int n_rhs, mxArray* matrix_rhs[], const char* function)
{
  if (strncmp(function, "disp", 4) == 0)
    {
      int cols = mxGetN(matrix_rhs[0]);
      int rows = mxGetM(matrix_rhs[0]);
      if (matrix_rhs[0]->type == mxCHAR_CLASS)
        {
          char* p = (char*)matrix_rhs[0]->data;
          for (int i = 0; i < rows; i++)
            {
              for (int j = 0; j < cols; j++)
                mexPrintf("%c", p[2*(i+j*rows)]);
              mexPrintf("\n");
            }
        }
      else
        for (int i = 0 ; i < rows; i++)
          {
            for (int j = 0; j < cols; j++)
              mexPrintf("%8.4f ", matrix_rhs[0]->data[i+j*rows]);
            mexPrintf("\n");
          }
    }
}

mxArray*
mexGetVariable(const char* space_name, const char* matrix_name)
{
  if (strncmp(space_name, "global", 6) != 0)
    mexErrMsgTxt("space_name not handle in mexGetVariable\n");
  return mxglobal[matrix_name];
}

int
mxGetFieldNumber(const mxArray* Struct, const char* field_name)
{
  vector<string>::const_iterator it = find(Struct->field_name.begin(), Struct->field_name.end(), field_name);
  if (it == Struct->field_name.end())
    {
      stringstream tmp;
      tmp << "unknown field name '" << field_name << "' in mxGetFieldNumber\n";
      mexErrMsgTxt(tmp.str());
    }
  return it - Struct->field_name.begin();
}

mxArray*
mxGetFieldByNumber(mxArray* Struct, unsigned int pos, unsigned int field_number)
{
  if (pos > Struct->size_1 * Struct->size_2)
    mexErrMsgTxt("index out of range in mxGetFieldByNumber\n");
  if (field_number > Struct->field_name.size())
    mexErrMsgTxt("field_number out of range in mxGetFieldByNumber\n");
  return Struct->field_array[pos][field_number];
}

int
mxAddField(mxArray* Struct, const char* fieldname)
{
  if (!mxIsStruct(Struct))
    return -1;
  else
    {
      Struct->field_name.push_back(fieldname);
      return Struct->field_name.size();
    }
}

void
mxSetFieldByNumber(mxArray *Struct, mwIndex index, unsigned int field_number, mxArray *pvalue)
{
  if (index >= Struct->size_1 * Struct->size_2)
    mexErrMsgTxt("index out of range in mxSetFieldByNumber\n");
  unsigned int nfields = Struct->field_name.size();
  if (field_number >= nfields)
    mexErrMsgTxt("field_number out of range in mxSetFieldByNumber\n");
  while (Struct->field_array.size() <= index)
    Struct->field_array.push_back(vector<mxArray*>(nfields, NULL));
  Struct->field_array[index][field_number] = pvalue;
}

int
mxGetString(const mxArray *array, char *buf, unsigned int buflen)
{
  unsigned int size;
  if (!mxIsChar(array))
    return -1;
  else if (buflen <= (size = mxGetNumberOfElements(array)))
    return -1;
  else
    {
      char *pchar = (char*)array->data;
      for (unsigned int i = 0; i < size; i++)
        buf[i] = pchar[2*i];
      buf[size+1] = ' ';
    }
  return 0;
}

mxArray*
mxDuplicateArray(const mxArray *array)
{
  unsigned int Size;
  unsigned int i;
  vector<mxArray*>::const_iterator it_v_array;
  vector<string>::const_iterator it_v_string;
  mxArray *Array = (mxArray*)mxMalloc(sizeof(mxArray));
  Array->type = array->type;
  Array->size_1 = array->size_1;
  Array->size_2 = Array->size_2;
  switch(array->type)
    {
    case mxDOUBLE_CLASS:
    case mxSINGLE_CLASS:
      Size = array->size_1*array->size_2*sizeof(double);
      Array->data = (double*)mxMalloc(Size);
      memcpy(Array->data, array->data, Size);
      break;
    case mxSTRUCT_CLASS:
      for (i = 0; i < array->size_1 * array->size_2; i++)
        for (it_v_array = array->field_array[i].begin(); it_v_array != array->field_array[i].end(); it_v_array++)
          Array->field_array[i].push_back(mxDuplicateArray(*it_v_array));
      for (it_v_string = array->field_name.begin(); it_v_string != array->field_name.end(); it_v_string++)
        Array->field_name.push_back(*it_v_string);
      break;
    case mxSPARSE_CLASS:
      Array->Nzmax = array->Nzmax;
      Size = array->Nzmax * sizeof(double);
      Array->data = (double*)mxMalloc(Size);
      memcpy(Array->data, array->data, Size);
      Size = array->Nzmax * sizeof(mwIndex);
      Array->Ir = (mwIndex*)mxMalloc(Size);
      memcpy(Array->Ir, array->Ir, Size);
      Size = array->size_2*sizeof(mwIndex);
      Array->Jc = (mwIndex*)mxMalloc(Size);
      memcpy(Array->Jc, array->Jc, Size);
      break;
    case mxCHAR_CLASS:
      Size = 2*array->size_1*array->size_2*sizeof(char);
      Array->data = (double*)mxMalloc(Size);
      memcpy(Array->data , array->data, Size);
      break;
    default:
      ostringstream tmp;
      tmp << "Array type not handle: " << array->type << "\n";
      mexErrMsgTxt(tmp.str());
    }
  return(Array);
}

mxArray*
read_struct(FILE *fid)
{
  unsigned int nfields = 0;
  fscanf(fid, "%d", &nfields);
  unsigned int size_1, size_2;
  fscanf(fid, "%d", &size_1);
  fscanf(fid, "%d", &size_2);
  vector<string> fieldnames;
  vector<mxArray*> v_Array;
  vector<vector<mxArray* > > vv_Array;
  for (unsigned int j = 0; j < size_1 * size_2; j++)
    {
      for (unsigned int i = 0; i < nfields; i++)
        {
          char name[512];
          fscanf(fid, "%s", name);
          if (j == 0)
            fieldnames.push_back(name);
          v_Array.push_back(read_Array(fid));
        }
      vv_Array.push_back(v_Array);
      v_Array.clear();
    }
  mxArray* Struct;
  Struct = mxCreateStructMatrix(size_1, size_2, nfields, fieldnames);
  for (unsigned int j = 0; j < size_1 * size_2; j++)
    {
      v_Array = vv_Array[j];
      for (unsigned int i = 0; i < nfields; i++)
        {
          mxSetFieldByNumber(Struct, j, i, v_Array[i]);
        }
    }
  return Struct;
}

mxArray*
read_double_array(FILE *fid)
{
  int size_1, size_2;
  fscanf(fid, "%d", &size_1);
  fscanf(fid, "%d", &size_2);
  mxArray *A = mxCreateDoubleMatrix(size_1, size_2, mxREAL);
  double* data = mxGetPr(A);
  float f;
  for (int i = 0; i < size_1*size_2; i++)
    {
      fscanf(fid, "%f", &f);
      data[i] = f;
    }
  return A;
}

mxArray*
read_char_array(FILE *fid)
{
  unsigned int size_1, size_2;
  fscanf(fid, "%d", &size_1);
  fscanf(fid, "%d", &size_2);
  mxArray *A = mxCreateCharArray(size_1, size_2, mxREAL);
  char* data = (char*)mxGetPr(A);
  char pchar[size_2+1];
  for (unsigned int i = 0; i < size_1; i++)
    {
      fscanf(fid, "%s", pchar);
      for(unsigned int j = 0; j < strlen(pchar); j++)
        data[2*(i+j*size_1)] = pchar[j];
      for(unsigned int j = strlen(pchar); j < size_2; j++)
        data[2*(i+j*size_1)] = ' ';
    }
  return A;
}

mxArray*
read_Array(FILE *fid)
{
  mxArray *Array;
  mxArray_type array_type;
  unsigned int j;
  fscanf(fid, "%d", &j);
  array_type = static_cast<mxArray_type>(j);
  switch(array_type)
    {
    case mxDOUBLE_CLASS:
    case mxSINGLE_CLASS:
      Array = read_double_array(fid);
      break;
    case mxSTRUCT_CLASS:
      Array = read_struct(fid);
      break;
    case mxCHAR_CLASS:
      Array = read_char_array(fid);
      break;
    default:
      ostringstream tmp;
      tmp << "Array type not handle in read_Array: " << array_type << "\n";
      mexErrMsgTxt(tmp.str());
    }
  return(Array);
}

void
load_global(char* file_name)
{
  FILE *fid;
  ostringstream tmp_out("");
  tmp_out << file_name << "_options.txt";
  fid = fopen(tmp_out.str().c_str(), "r");
  if (fid == NULL)
    {
      string s = "Can't open file ";
      s += tmp_out.str().c_str();
      mexErrMsgTxt(s);
    }
  mxglobal["options_"] = read_struct(fid);
  fclose(fid);

  tmp_out.str("");
  tmp_out << file_name << "_M.txt";
  fid = fopen(tmp_out.str().c_str(), "r");
  if (fid == NULL)
    {
      string s = "Can't open file ";
      s += tmp_out.str().c_str();
      mexErrMsgTxt(s);
    }
  mxglobal["M_"] = read_struct(fid);
  fclose(fid);

  tmp_out.str("");
  tmp_out << file_name << "_oo.txt";
  fid = fopen(tmp_out.str().c_str(), "r");
  if (fid == NULL)
    {
      string s = "Can't open file ";
      s += tmp_out.str().c_str();
      mexErrMsgTxt(s);
    }
  mxglobal["oo_"] = read_struct(fid);
  fclose(fid);
}

void
Free_simple_array(mxArray* array)
{
  mxFree(array->data);
  delete array;
}

void
Free_struct_array(mxArray* array)
{
  for (unsigned int i = 0; i < array->size_1 * array->size_2; i++)
    for (unsigned int j = 0; j < array->field_name.size(); j++)
      Free_Array(array->field_array[i][j]);
  delete array;
}

void
Free_Array(mxArray* array)
{
  switch(array->type)
    {
    case mxDOUBLE_CLASS:
    case mxSINGLE_CLASS:
    case mxCHAR_CLASS:
      Free_simple_array(array);
      break;
    case mxSTRUCT_CLASS:
      Free_struct_array(array);
      break;
    default:
      mexErrMsgTxt("Array type not handle in read_Array\n");
    }
}

void
Free_global()
{
  for (map<string, mxArray*>::iterator it = mxglobal.begin(); it != mxglobal.end(); it++)
    Free_Array(it->second);
}
