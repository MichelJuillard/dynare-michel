/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the NumericalConstants class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
using namespace std;
//------------------------------------------------------------------------------
#include "NumericalConstants.h"
//------------------------------------------------------------------------------
vector<string> NumericalConstants::mNumericalConstants = *(new vector<string>);	
//------------------------------------------------------------------------------
NumericalConstants::NumericalConstants()
{
	mNumericalConstants.push_back("0.0");
	mNumericalConstants.push_back("1.0");	

}
//------------------------------------------------------------------------------
NumericalConstants::~NumericalConstants()
{
	// Empty
}
//------------------------------------------------------------------------------
int NumericalConstants::AddConstant(string iConst)
{
	if (iConst == "0.0")
		return 0;
	else if (iConst == "1.0")
		return 1;
	mNumericalConstants.push_back(iConst);	
	return (int) mNumericalConstants.size()-1;		
}
//------------------------------------------------------------------------------
string NumericalConstants::get(int ID)
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
