/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the Shocks class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
using namespace std;
//------------------------------------------------------------------------------
#include "Shocks.h"
#include "ModelParameters.h"
#include "Interface.h"

static int mshock_flag = 0;
static int exo_det_length = 0;

//------------------------------------------------------------------------------
//ostringstream	Shocks::output;
//------------------------------------------------------------------------------
Shocks::Shocks()
{ 
	// Empty
}
//------------------------------------------------------------------------------
Shocks::~Shocks()
{
	// Empty
}
//------------------------------------------------------------------------------
void Shocks::setOutput(ostringstream* iOutput)
{
	output = iOutput;	
}
//------------------------------------------------------------------------------
void Shocks::BeginShocks(void)
{	
  mshock_flag = 0;
  // Writing a Matlab comment
  *output << interfaces::comment() << "\n" << interfaces::comment() << "SHOCKS instructions \n"
          << interfaces::comment() << "\n";
  //Writing intstruction that initialize a shocks 
  *output << "make_ex_;\n";
}							
//------------------------------------------------------------------------------
void Shocks::BeginMShocks(void)
{	
  mshock_flag = 1;
  // Writing a Matlab comment
  *output << interfaces::comment() << "\n" << interfaces::comment() << "MSHOCKS instructions \n"
          << interfaces::comment() << "\n";
  //Writing intstruction that initialize a shocks 
  *output << "make_ex_;\n";
}
void Shocks::EndShocks(void)
{
  *output << "M_.exo_det_length = " << exo_det_length << ";\n";
}							
//------------------------------------------------------------------------------
void Shocks::AddDetShockExo(int id1)
{
  if (mPeriod1.size() != mPeriod2.size() ||
      mPeriod1.size() != mValues.size())
    {
      string msg = "in shocks statement, number of periods and values dont agree";
      (* error) (msg.c_str());
    }
  for (unsigned int i = 0; i < mPeriod1.size(); i++)
    {
      string period1 = mPeriod1[i];
      string period2 = mPeriod2[i];
      if (period1 == period2)
	{
	  *output << "set_shocks(" << mshock_flag << "," << period1 
		  << ", " << id1+1 << ", " << mValues[i] 
		  << ");\n";
	}
      else
	{
	  *output << "set_shocks(" << mshock_flag << "," << period1
		  << ":" << period2 << ", " << id1+1 
		  << ", " << mValues[i] << ");\n";
	}
    }
  mPeriod1.clear();
  mPeriod2.clear();
  mValues.clear();
}
void Shocks::AddDetShockExoDet(int id1)
{
  if (mPeriod1.size() != mPeriod2.size() ||
      mPeriod1.size() != mValues.size())
    {
      string msg = "in shocks statement, number of periods and values dont agree";
      (* error) (msg.c_str());
    }
  for (unsigned int i = 0; i < mPeriod1.size(); i++)
    {
      string period1 = mPeriod1[i];
      string period2 = mPeriod2[i];
      if (period1 == period2)
	{
	  *output << "set_shocks(" << mshock_flag + 2 << "," << period1 
		  << ", " << id1+1 << ", " << mValues[i] 
		  << ");\n";
	}
      else
	{
	  *output << "set_shocks(" << mshock_flag+2 << "," << period1
		  << ":" << period2 << ", " << id1+1 
		  << ", " << mValues[i] << ");\n";
	}
      int p2_int = atoi(period2.c_str()); 
      if (p2_int > exo_det_length)
	{
	  exo_det_length = p2_int;
	}
    }
  mPeriod1.clear();
  mPeriod2.clear();
  mValues.clear();
}
void Shocks::AddSTDShock(int id1, std::string value)
{
  *output << "M_.Sigma_e(" << id1+1 << ", " << id1+1 << ") = " << value << "^2;\n";
  mPeriod1.clear();
  mPeriod2.clear();
  mValues.clear();
}
void Shocks::AddVARShock(int id1, std::string value)
{
  *output << "M_.Sigma_e(" << id1+1 << ", " << id1+1 << ") = " << value << ";\n";
  mPeriod1.clear();
  mPeriod2.clear();
  mValues.clear();
}
void Shocks::AddCOVAShock(int id1, int id2 , std::string value)
{
  *output << "M_.Sigma_e(" << id1+1 << ", " << id2+1 << ") = " <<	value <<
    "; M_.Sigma_e(" << id2+1 << ", " << id1+1 << ") = M_.Sigma_e(" << 
    id1+1 << ", " << id2+1 << ");\n";
  mPeriod1.clear();
  mPeriod2.clear();
  mValues.clear();
}
void Shocks::AddCORRShock(int id1, int id2 , std::string value)
{
  *output << "M_.Sigma_e(" << id1+1 << ", " << id2+1 << ") = " <<
    value << "*sqrt(M_.Sigma_e(" << id1 << ", " << id1+1 << ")*M_.Sigma_e(" <<
    id2 << ", " << id2+1 << "); M_.Sigma_e(" << id2+1 << ", " << 
    id1 << ") = M_.Sigma_e(" << id1+1 << ", " << id2+1 << ");\n";
  mPeriod1.clear();
  mPeriod2.clear();
  mValues.clear();
}
  /*
    sets a shock 
    for type == 0, set type, id1, shock_elem will be set by another method
    for type == 1 or 2 set type, id1 and value
    for type == 3 or 4 set type, id1, id2 and value
  */
//------------------------------------------------------------------------------
/*
string Shocks::get(void)
{
	return output.str();
}
*/
//------------------------------------------------------------------------------
void Shocks::AddPeriod(string p1, string p2)
{
	mPeriod1.push_back(p1);
	mPeriod2.push_back(p2);
}
//------------------------------------------------------------------------------
void Shocks::AddPeriod(string p1)
{
	mPeriod1.push_back(p1);
	mPeriod2.push_back(p1);

}
//------------------------------------------------------------------------------
void Shocks::AddValue(string value)
{
	mValues.push_back(value);
}
//------------------------------------------------------------------------------
