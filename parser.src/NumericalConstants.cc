/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the NumericalConstants class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
using namespace std;
//------------------------------------------------------------------------------
#include "NumericalConstants.hh"

NumericalConstants::NumericalConstants()
{
  AddConstant("0.0");
  AddConstant("1.0");
}

//------------------------------------------------------------------------------
NumericalConstants::~NumericalConstants()
{
  // Empty
}

//------------------------------------------------------------------------------
int NumericalConstants::AddConstant(string iConst)
{
  map<string, int, less<string> >::iterator iter = numConstantsIndex.find(iConst);

  if (iter != numConstantsIndex.end())
    return iter->second;

  int id = (int) mNumericalConstants.size();
  mNumericalConstants.push_back(iConst);
  numConstantsIndex[iConst] = id;
  return id;
}

//------------------------------------------------------------------------------
string NumericalConstants::get(int ID) const
{
  if (ID < (int)mNumericalConstants.size())
    {
      return mNumericalConstants[ID];
    }
  else
    {
      return "";
    }
}

//------------------------------------------------------------------------------
