/*
* Copyright (C) 2008-2009 Dynare Team
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

#include "mexutils.h"


MexStruct::MexStruct( const int numparstruct):
numParStruct(numparstruct), structName(string(""))
  { 
  // get Dynare mexSturcture pointers and store them locally
#ifdef DEBUG
    mexPrintf("MexStruct reserve=%d \n", numParStruct);
#endif
  parStruct.reserve(numParStruct);
  parStructBase.reserve(numParStruct);
  for (int i=0;i<numParStruct;++i)
    {
    parStruct[i]=mexGetVariablePtr("caller", DynareParamStructsNm[i]);
    parStructBase[i]=caller;
#ifdef DEBUG
    mexPrintf("MexStruct to insert i=%d parStructNm[i]=%s using base[i]=%d  %s\n", i, DynareParamStructsNm[i], parStructBase[i], mexBase[parStructBase[i]]);
#endif
    // get field names into the map:
    pair <map<string, int>::iterator, bool> ret;
    int j=0;
    const char* field;
    while ((field=mxGetFieldNameByNumber(parStruct[i],j))!=NULL)
      {
#ifdef DEBUG
        mexPrintf("MexStruct insert field= %s\n", field);
#endif
      ret=parNamStructMap.insert(make_pair(string(field),i));
      if (!ret.second)
            mexPrintf("MexStruct failed to insert field %s from struct %d %s using base %d %s\n"
            , field, i, DynareParamStructsNm[i], parStructBase[i], mexBase[parStructBase[i]]);
//        mexErrMsgTxt("MexStruct: Failed to insert param \n");
      j++;
      }
    }
#ifdef DEBUG
        mexPrintf("MexStruct insert finished \n");
#endif
  }

MexStructParam& 
MexStruct::getMexStructField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  return *(new MexStructParam(mxf, this, field));
  }

void * 
MexStruct::getField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  return mxGetData(mxf);
  }

GeneralMatrix & 
MexStruct::getMatrixField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  return *(new GeneralMatrix(mxGetPr(mxf), mxGetM(mxf),  mxGetN(mxf))) ;
  }

double& 
MexStruct::getDoubleField(const string& field)
  {
  const mxArray* mxf = getMxField( field);
  if (!mxIsDouble(mxf))
    mexErrMsgTxt("Input must be of type double.");
//  double ret=mxGetScalar(mxf);
//  return ret;
  double *ret=(new double);
  *ret=mxGetScalar(mxf);
  return *ret;
  }
 
string& 
MexStruct::getStringField(const string& field)
  {
  mxArray* mxf = getMxField( field);
  if (!mxIsChar(mxf))
    mexErrMsgTxt("Input must be of type char.");
  return *(new string(mxArrayToString(mxf)));
  }

char * 
MexStruct::getCharField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  if (!mxIsChar(mxf))
    mexErrMsgTxt("Input must be of type char.");
  return mxArrayToString(mxf);
  }

vector<string>& 
MexStruct::getStringVectorField(const string& field)
  {
  charArraySt * cap = getCharArrayField(field);
  vector <string>*sv=(new vector<string>(cap->len));
  for (int i= 0; i<cap->len;++i)
    (*sv)[i]=string(cap->charArrayPtr[i]);
  return *sv;
  }

vector<int>& 
MexStruct::getIntVectorField(const string& field)
  {
  mxArray* mxfp = getMxField( field);
  if (MIN(mxGetM(mxfp),mxGetN(mxfp))!=1)
    throw SYLV_MES_EXCEPTION("Int vector is a 2D Matrix .");
  double* dparams = (double *) mxGetData(mxfp);
  int npar = (int) MAX(mxGetM(mxfp), mxGetN(mxfp));
  vector<int> *vars = (new vector<int>(npar));
  for (int v = 0; v < npar; v++)
    {
    (*vars)[v] = (int)(*(dparams++)); //[v];
#ifdef DEBUG
    mexPrintf("%s[%d]=%d.\n", field.c_str() , v, (*vars)[v]);
#endif
    }
  
  return *vars;
  };


Vector& 
MexStruct::getDoubleVectorField(const string& field)
  {
  mxArray* mxfp = getMxField( field);
  if (MIN(mxGetM(mxfp), mxGetN(mxfp))!=1)
    throw SYLV_MES_EXCEPTION("Double Vector is a 2D Matrix .");
  double* dparams = (double *) mxGetData(mxfp);
  int npar = (int) MAX(mxGetM(mxfp), mxGetN(mxfp));
  Vector *vars = (new Vector(dparams,npar));
  return *vars;
  };



charArraySt * 
MexStruct::getCharArrayField( const string& field)
  {
  mxArray* mxfp = getMxField( field);
  const int len = (int) mxGetM(mxfp);
  const int width = (int) mxGetN(mxfp);
  if (!mxIsChar(mxfp))
    mexErrMsgTxt("Input must be of type char.");

  charArraySt * cap=new charArraySt;
  cap->charArrayPtr = (char**)MxArrayToStringArray(mxfp, len, width);
  cap->len=len;
  return cap;
  }

void 
MexStruct::ReorderCols(GeneralMatrix &tdx, const vector<int> &vOrder)
  {MexStructParam::ReorderCols(tdx,vOrder, "static");};


void
MexStruct::setField( const string& field,  const string& val)
  {
  mxArray *newVal=mxCreateString(val.c_str());
  UpdateMexField( field,  newVal);  
  }
  /*******
  {
  mxArray *fp = mxGetField(mexStruct, 0, field);
  mxDestroy(fp)
  mxSetField(mexStruct,0,field, newVal);
  }
  ***************/

void
MexStruct::setField( const string& field, const double val)
  {
  mxArray *newVal=mxCreateDoubleScalar(val);
  UpdateMexField( field,  newVal);
  }

void
MexStruct::setField( const string& field, const GeneralMatrix& gmval)
  {
  mxArray *newVal=mxCreateDoubleMatrix(gmval.numRows(),gmval.numCols(),mxREAL);
  memcpy(mxGetPr(newVal),gmval.base(),gmval.numRows()*gmval.numCols()*sizeof(double));
  UpdateMexField( field, newVal);
  }

void
MexStruct::UpdateMexField( const string& field, mxArray *newVal)
  {
  StructBaseNameMap* dpsm=getMxFieldStruct(field);
  mxArray *sp =mexGetVariable(dpsm->baseName, dpsm->structName);
  if (sp == NULL )
  {
    mexPrintf("Variable not found : base: %s, structName %s\n", dpsm->baseName, dpsm->structName);
    mexErrMsgTxt(" \n");
  }
  mxArray *fp = mxGetField(sp, 0, field.c_str());
  mxDestroyArray(fp);
  mxSetField(sp,0,field.c_str(), newVal);
  mexPutVariable(dpsm->baseName, dpsm->structName, sp);
  //mxDestroyArray(newVal);
  }

/***********************************************
* MexStructParam 
* holds single Matlab structure passed as parameter
***********************************************/

MexStructParam& 
MexStructParam::getMexStructField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  return *(new MexStructParam(mxf, this, field));
  }

void * 
MexStructParam::getField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  return mxGetData(mxf);
  }

double& 
MexStructParam::getDoubleField(const string& field)
  {
  const mxArray* mxf = getMxField( field);
  if (!mxIsDouble(mxf))
    mexErrMsgTxt("Input must be of type double.");
  double *ret=(new double);
  *ret=mxGetScalar(mxf);
  return *ret;
  }

GeneralMatrix&
MexStructParam::getMatrixField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  return *(new GeneralMatrix(mxGetPr(mxf), mxGetM(mxf),  mxGetN(mxf)) );
  }

string& 
MexStructParam::getStringField(const string& field)
  {
  mxArray* mxf = getMxField( field);
  if (!mxIsChar(mxf))
    mexErrMsgTxt("Input must be of type char.");
  return *(new string(mxArrayToString(mxf)));
  }

char * 
MexStructParam::getCharField( const string& field)
  {
  mxArray* mxf = getMxField( field);
  if (!mxIsChar(mxf))
    mexErrMsgTxt("Input must be of type char.");
  return mxArrayToString(mxf);
  }

vector<string>& 
MexStructParam::getStringVectorField(const string& field)
  {
  charArraySt * cap = getCharArrayField( field);
  vector <string>*sv= (new vector<string>(cap->len));
  for (int i= 0; i<cap->len;++i)
    (*sv)[i]=string(cap->charArrayPtr[i]);
  return *sv;
  }

vector<int>& 
MexStructParam::getIntVectorField(const string& field)
  {
  mxArray* mxfp = getMxField( field);
  if (MIN(mxGetM(mxfp), mxGetN(mxfp))!=1)
    throw SYLV_MES_EXCEPTION("Int vector is a 2D Matrix .");
  double* dparams = (double *) mxGetData(mxfp);
  int npar = (int) MAX(mxGetM(mxfp), mxGetN(mxfp));
  vector<int> *vars = (new vector<int>(npar));
  for (int v = 0; v < npar; v++)
    {
    (*vars)[v] = (int)(*(dparams++)); //[v];
#ifdef DEBUG
    mexPrintf("%s[%d]=%d.\n", field.c_str(), v, (*vars)[v]);
#endif
    }
  
  return *vars;
  };


Vector& 
MexStructParam::getDoubleVectorField(const string& field)
  {
  mxArray* mxfp = getMxField( field);
  if (MIN(mxGetM(mxfp), mxGetN(mxfp))!=1)
    throw SYLV_MES_EXCEPTION("Double Vector is a 2D Matrix .");
  double* dparams = (double *) mxGetData(mxfp);
  int npar = (int) MAX(mxGetM(mxfp), mxGetN(mxfp));
  Vector *vars = (new Vector(dparams,npar));
  return *vars;
  };




charArraySt * 
MexStructParam::getCharArrayField( const string& field)
  {
  mxArray* mxfp = getMxField( field);
  const int len = (int) mxGetM(mxfp);
  const int width = (int) mxGetN(mxfp);
  if (!mxIsChar(mxfp))
    mexErrMsgTxt("Input must be of type char.");

  charArraySt * cap=new charArraySt;
  cap->charArrayPtr = (char**)MxArrayToStringArray(mxfp, len, width);
  cap->len=len;
  return cap;
  }

void
MexStructParam::setField( const string& field,  const string& val)
  {
  mxArray *newVal=mxCreateString(val.c_str());
  mxSetField((mxArray*)parStruct,0,field.c_str(), newVal);
  //UpdateMexField( field,  newVal);  
  }
  /*******
  {
  mxArray *fp = mxGetField(mexStruct, 0, field);
  mxDestroy(fp)
  mxSetField(mexStruct,0,field, newVal);
  }
  ***************/

void
MexStructParam::setField( const string& field, const double val)
  {
  mxArray *newVal=mxCreateDoubleScalar(val);
  mxSetField((mxArray*)parStruct,0,field.c_str(), newVal);
  //UpdateMexField( field,  newVal);
  }

void
MexStructParam::setField( const string& field, const GeneralMatrix& gmval)
  {
  mxArray *newVal=mxCreateDoubleMatrix(gmval.numRows(),gmval.numCols(),mxREAL);
  memcpy(mxGetPr(newVal),gmval.base(),gmval.numRows()*gmval.numCols()*sizeof(double));
  mxSetField((mxArray*)parStruct,0,field.c_str(), newVal);
  //UpdateMexField( field, newVal);
  }

void
MexStructParam::UpdateMexField( const string& field, mxArray *newVal)
  {
 
  mxSetField((mxArray*)parStruct,0,field.c_str(), newVal);
 
 /************
  if (parStructParent!=NULL)
    parStructParent->setField(field,newVal);
  else
    {
  
    mxArray *sp =mexGetVariable("caller",structName);
    if (sp == NULL )
      {
        mexPrintf("Variable not found : base: %s, structName %s\n", dpsm->baseName, dpsm->structName);
        mexErrMsgTxt(" \n");
      }
    mxArray *fp = mxGetField(sp, 0, field.c_str());
    mxDestroyArray(fp);
    mxSetField(sp,0,field.c_str(), newVal);
    mexPutVariable(dpsm->baseName, dpsm->structName, sp);
    //mxDestroyArray(newVal);
    }
  ************/
  }




/************************************
* Reorder first  variables in a vector
* according to order given in  varsOrder 

************************************/

void
//MexStructParam::ReorderCols(GeneralMatrix &tdx, const vector<int> *vOrder)
MexStructParam::ReorderCols(GeneralMatrix &tdx, const vector<int> &vOrder, char* stat)
{

  if (tdx.numCols() > vOrder.size())
    {
      mexPrintf(" Error in ReorderColumns - size of order var is too small");
      return;
    }
//  GeneralMatrix tmp(*tdx); // temporary 2D matrix
  GeneralMatrix tmpR(tdx); // temporary 2D matrix
//  GeneralMatrix &tmpR = tmp;
  tdx.zeros(); // empty original matrix
  // reorder the columns
  try
    {
      for (int i = 0; i < tdx.numCols(); i++)
        tdx.copyColumns(tmpR, (vOrder)[i],vOrder[i], i);
//        tdx->copyColumn(tmpR, (*vOrder)[i], i);
    }
  catch (const SylvException &e)
    {
      printf("Caugth exception in ReorderColumns: ");
      e.printMessage();
      return; // 255;
    }
  catch (...)
    {
      mexPrintf(" Error in ReorderColumns - wrong index?");
    }
}
void
MexStructParam::ReorderCols(GeneralMatrix &tdx, const int *vOrder)
{

//  GeneralMatrix tmp(*tdx); // temporary 2D matrix
  GeneralMatrix tmpR(tdx); // temporary 2D matrix
//  GeneralMatrix &tmpR = tdztmp;
  tdx.zeros(); // empty original matrix
  // reorder the columns
  try
    {
      for (int i = 0; i < tdx.numCols(); i++)
        tdx.copyColumns(tmpR, vOrder[i],vOrder[i], i);
    }
  catch (const SylvException &e)
    {
      printf("Caugth  SYLV_EXCEPTION in ReorderColumns: ");
      e.printMessage();
      return; // 255;
    }
  catch (...)
    {
      mexPrintf(" Error in ReorderColumns - wrong index?");
    }
}

/**************
void
//MexStructParam::ReorderCols(GeneralMatrix *tdx, const vector<int> *vOrder)
MexStructParam::ReorderCols(GeneralMatrix *tdx, const vector<int> &vOrder)
{

  if (tdx->ncols() > vOrder->size())
    {
      mexPrintf(" Error in ReorderColumns - size of order var is too small");
      return;
    }
  GeneralMatrix tmp(*tdx); // temporary 2D matrix
//  GeneralMatrix& tmpR(tdx); // temporary 2D matrix
  GeneralMatrix &tmpR = tmp;
  tdx->zeros(); // empty original matrix
  // reorder the columns
  try
    {
      for (int i = 0; i < tdx->ncols(); i++)
        tdx->copyColumn(tmpR, (vOrder)[i], i);
//        tdx->copyColumn(tmpR, (*vOrder)[i], i);
    }
  catch (const TLException &e)
    {
      printf("Caugth TL exception in ReorderColumns: ");
      e.print();
      return; // 255;
    }
  catch (...)
    {
      mexPrintf(" Error in ReorderColumns - wrong index?");
    }
}
void
MexStructParam::ReorderCols(GeneralMatrix *tdx, const int *vOrder)
{

  GeneralMatrix tmp(*tdx); // temporary 2D matrix
//  GeneralMatrix tmpR(tdx); // temporary 2D matrix
  GeneralMatrix &tmpR = tdztmp;
  tdx->zeros(); // empty original matrix
  // reorder the columns
  try
    {
      for (int i = 0; i < tdx->ncols(); i++)
        tdx->copyColumn(tmpR, vOrder[i], i);
    }
  catch (const TLException &e)
    {
      printf("Caugth TL exception in ReorderColumns: ");
      e.print();
      return; // 255;
    }
  catch (...)
    {
      mexPrintf(" Error in ReorderColumns - wrong index?");
    }
}

*******************/
//////////////////////////////////////////////////////
// Convert Matlab string array to C type array of string pointers
// Poblem is that Matlab mx function returns a long string concatenated by columns rather than rows
// hence a rather low level approach is needed
///////////////////////////////////////////////////////
const char **
MxArrayToStringArray(const mxArray *mxFldp, const int len, const int width)
  {
  char *cNamesCharStr = mxArrayToString(mxFldp);
  const char **ret = MxArrayToStringArray(cNamesCharStr, len, width);
  return ret;
  }

const char **
MxArrayToStringArray(const char *cNamesCharStr, const int len=1, const int width=1)
  {
  char cNamesMX[len][width+1]; //
#ifdef DEBUG
  mexPrintf("loop MxArrayToStringArray cNamesCharStr = %s \n", cNamesCharStr);
#endif
  for (int i = 0; i < width; i++)
    {
    for (int j = 0; j < len; j++)
      {
      // Allow alphanumeric and underscores "_" only:
      if (isalnum(cNamesCharStr[j+i*len]) || ('_' == cNamesCharStr[j+i*len]))
        {
        cNamesMX[j][i] = cNamesCharStr[j+i*len];
        }
      else cNamesMX[j][i] = '\0';
      }
    }
  const char **ret = (const char **) mxCalloc(len, sizeof(char *));
  for (int j = 0; j < len; j++)
    {
    cNamesMX[j][width] = '\0';
#ifdef DEBUG
    //				mexPrintf("String [%d]= %s \n", j, cNamesMX[j]);
#endif
    char *token = (char *) mxCalloc(strlen(cNamesMX[j])+1, sizeof(char));
    strcpy(token, cNamesMX[j]);
    ret[j] = token;
#ifdef DEBUG
    mexPrintf("ret [%d]= %s \n", j, ret[j]);
#endif
    }
  return ret;
  }
