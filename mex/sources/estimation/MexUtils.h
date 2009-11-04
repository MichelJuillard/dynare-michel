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
#ifndef MEXUTILS_H
#define MEXUTILS_H
#include <utility>
#include <vector>
#include "ioutils.h"
#include "mex.h"
#include "matrix.h"
//#include "SylvException.h"
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


//#include "k_ord_dynare.h"
//using namespace std;
enum {base, caller, global};
enum  {M_,  oo_, options_,bayestopt_, estim_params_, dr};
extern const char *DynareParamStructsNm [];
extern const char* mexBase[];

const int numParStructs=5;
//#define numParStructsM 6


/**
struct DynareParamsStructMap
{ char * paramName;
int dynParStruct;
}

  struct DynareParamsStructBasePair
  { 
  int dynParStruct;
  int mexBase;
  }
**/



struct StructBaseNameMap
  { char* baseName;
char* structName;
  };

class MexStructParam;

class MexStruct :public  virtual GeneralParams 
  {
  vector <int> parStructBase;
  vector <const mxArray*> parStruct;  // array of struct pointers
  map <string, int> parNamStructMap; // which struct par belongs
  const int numParStruct;
  const string structName;
  
  mxArray* getMxField( const string& field) 
    {
    map<string, int>::iterator it=parNamStructMap.find(field);
    if (it==parNamStructMap.end())
      throw(SYLV_MES_EXCEPTION("no parameter with such name"));
    return mxGetField(parStruct[it->second], 0, field.c_str() );
    }
  
  StructBaseNameMap* getMxFieldStruct( const string& field) 
    {
    map<string, int>::iterator it=parNamStructMap.find(field);
    if (it==parNamStructMap.end())
      throw(SYLV_MES_EXCEPTION("no parameter with such name"));
    StructBaseNameMap* dpsm=new StructBaseNameMap;
    dpsm->baseName=(char*)mexBase[parStructBase[it->second]];
    dpsm->structName=(char*)DynareParamStructsNm[it->second];
    
    return dpsm;
    }

//  void ReorderCols(GeneralMatrix* tdx, const vector<int>*vOrder);
//    {KordpDynare::ReorderCols((TwoDMatrix*) tdx,  *vOrder)};
//  void ReorderCols(GeneralMatrix* tdx, const char*vOrder);
//    {KordpDynare::ReorderCols((TwoDMatrix*) tdx,  *vOrder)};
  public:
    MexStruct(int numParStr=1);
    /**
    MexStruct( char *sBases,  char * mexStructNames, int numParStr=1)
    : base(sBases), structName(mexStructNames), numParStruct(numParStr)
    { 
    
      mexStruct=mexGetVariable(base, structName);
      }
    **/
    ~MexStruct(){};
    void * 
      getMxArrayField(const string& field)
      {
      return getMxField( field);
      }
    string& name(){return *(new string(structName));};
    MexStructParam& getMexStructField( const string& field);
    GeneralParams& getStructField( const string& field)
      {return (GeneralParams&) getMexStructField( field);};
    void *  getField( const string& field);
    double& getDoubleField(const string& field);
    string& getStringField(const string& field);
    vector<string>& getStringVectorField(const string& field);
    vector<int>& getIntVectorField(const string& field);
    Vector& getDoubleVectorField(const string& field);
    char *  getCharField( const string& field);
    charArraySt * getCharArrayField( const string& field);
    GeneralMatrix & getMatrixField( const string& field);
    void ReorderCols(GeneralMatrix& tdx, const vector<int>&vOrder);
//      (MexStructParam::ReorderCols(tdx, vOrder););
//    static void ReorderCols(GeneralMatrix& tdx, const int*vOrder);
//      (MexStructParam::ReorderCols(tdx, vOrder););
    void setField(const string& field, const string& newVal);
    void setField(const string& field, const double val);
    void setField(const string& field, const GeneralMatrix& val);
    void UpdateMexField( const string& field, mxArray *newVal);
    //void putMexStruct();//{mexPutVariable(base, structName, mexStruct);};
  };

/***********************************************
* MexStructParam 
* holds single Matlab structure passed as parameter
***********************************************/

class MexStructParam :public virtual GeneralParams 
  {
  const mxArray* parStruct; // struct pointer
  const GeneralParams* parStructParent; // if any
  const string structName; // if any, param name of the structure in its parent.
  
  mxArray* getMxField( const string& field)
    {
    return mxGetField(parStruct, 0, field.c_str());
    }
  
//  void ReorderCols(GeneralMatrix* tdx, const vector<int>*vOrder);
//  void ReorderCols(GeneralMatrix* tdx, const char*vOrder);

  public:
    MexStructParam(const mxArray* paramStruct, const GeneralParams* parent, const string& name):
        parStruct(paramStruct), parStructParent(parent), structName(name) {};
    MexStructParam( const mxArray* paramStruct, const GeneralParams* parent):
        parStruct(paramStruct), parStructParent(parent), structName(string("")){};
    MexStructParam(const mxArray* paramStruct, const string& name):
        parStruct(paramStruct), parStructParent(NULL), structName(name){};
    ~MexStructParam(){};
    void * 
      getMxArrayField(const string& field)
      {
      return getMxField( field);
      }
    string& name(){return *(new string(structName));};
    MexStructParam& getMexStructField( const string& field);
    GeneralParams& getStructField( const string& field)
      {return (GeneralParams&) getMexStructField( field);};
    void *  getField( const string& field);
    double& getDoubleField(const string& field);
    string& getStringField(const string& field);
    vector<string>& getStringVectorField(const string& field);
    vector<int>& getIntVectorField(const string& field);
    Vector& getDoubleVectorField(const string& field);
    char *  getCharField( const string& field);
    charArraySt * getCharArrayField( const string& field);
    GeneralMatrix & getMatrixField( const string& field);
    void ReorderCols(GeneralMatrix& tdx, const vector<int>&vOrder)
      {ReorderCols(tdx, vOrder, "stat");};
    static void ReorderCols(GeneralMatrix& tdx, const vector<int>&vOrder, char*stat);
    static void ReorderCols(GeneralMatrix& tdx, const int*vOrder);
    void setField(const string& field, const string& newVal);
    void setField(const string& field, const double val);
    void setField(const string& field, const GeneralMatrix& val);
    void UpdateMexField( const string& field, mxArray *newVal);
    //void putMexStruct();//{mexPutVariable(base, structName, mexStruct);};
  };

  /***************
  class ConstMexStruct :public GeneralParams 
  {
  const mxArray* mexStruct;
  DynareParamsStructMap dynpsm;
  
    const char *base;  // one of: base, caller, global
    const char *structName;
    MexStruct( char *sBase,  char * mexStructName)
    : base(sBase), structName(mexStructName)
    { 
    mexStruct=mexGetVariablePtr(base, structName);
    }
    virtual ~MexStruct(){};
    mxArray * 
    getField(string& field)
    {
    return mxGetField(mexStruct, 0, field);
    }
    };
    
*******************/

//////////////////////////////////////////////////////
// Convert Matlab Dynare string array to C type array of string pointers
// Poblem is that Matlab mx function returns a long string concatenated by columns rather than rows
// hence a rather low level approach is needed
///////////////////////////////////////////////////////
const char **
MxArrayToStringArray(const mxArray *mxFldp, const int len, const int width);
const char **
MxArrayToStringArray(const char *cNamesCharStr, const int len, const int width);


#endif