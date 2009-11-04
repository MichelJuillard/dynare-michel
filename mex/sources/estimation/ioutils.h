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
#ifndef IOUTILS_H
#define  IOUTILS_H
#include "GeneralMatrix.h"
#include "SylvException.h"
#include <map>
#include <string>
//using namespace std;


struct charArraySt
  {
  char ** charArrayPtr;
  int len;
  };

class GeneralParams 
  {
//  map <string, int> params;
//  const char *structName;
public:
  GeneralParams(){};
  virtual ~GeneralParams(){};
  virtual string& name()=0;
  virtual void * 
    getField(const string& field)=0;
  virtual double&
    getDoubleField(const string& field)=0;
  virtual string &
    getStringField(const string& field)=0;
  virtual vector<string>&
    getStringVectorField(const string& field)=0;
  virtual vector<int>&
    getIntVectorField(const string& field)=0;
  virtual Vector&
    getDoubleVectorField(const string& field)=0;
  virtual GeneralParams&
     getStructField( const string& field)=0;
  virtual char * 
    getCharField(const string& field)=0;
  virtual charArraySt * 
    getCharArrayField( const string& field)=0;
  // uses General Matrix from sylv library
  virtual GeneralMatrix& 
    getMatrixField(const string& field)=0;
  virtual void ReorderCols(GeneralMatrix& tdx, const vector<int>&vOrder)=0;
//  virtual void ReorderCols(GeneralMatrix& tdx, const int*vOrder)=0;
  virtual void
 // overloaded pramater update "set" methods:
    setField(const string& field, const string&newVal)=0;
  virtual void
    setField(const string& field, const double val)=0;
  virtual void
    setField(const string& field, const GeneralMatrix& val)=0;
  };


#endif
