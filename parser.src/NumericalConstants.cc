#include <iostream>

#include "NumericalConstants.hh"

NumericalConstants::NumericalConstants()
{
  AddConstant("0");
  AddConstant("1");
}

int
NumericalConstants::AddConstant(const string &iConst)
{
  map<string, int>::iterator iter = numConstantsIndex.find(iConst);
  //cout << "iConst=" << iConst << "\n" ;
  if (iter != numConstantsIndex.end())
    return iter->second;

  if (atof(iConst.c_str()) < 0)
    {
      cerr << "Can't handle a negative constant..!" << endl;
      exit(-1);
    }

  int id = (int) mNumericalConstants.size();
  mNumericalConstants.push_back(iConst);
  numConstantsIndex[iConst] = id;
  return id;
}

string
NumericalConstants::get(int ID) const
{
  if (ID < (int) mNumericalConstants.size())
    return mNumericalConstants[ID];
  else
    {
      cerr << "Unknown constant" << endl;
      exit(-1);
    }
}

double
NumericalConstants::getDouble(int iID) const
{
  return(atof(get(iID).c_str()));
}
