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
	// Writing a Matlab comment
	*output << "%\n% (M)SHOCKS instructions \n%\n";
	//Writing intstruction that initialize a shocks 
  	*output << "make_ex_;\n";
}							
//------------------------------------------------------------------------------
void Shocks::AddShock(shocktype type, int id1, int id2, string value)
{
	Shock shock;
	
	switch (type)
	{
		case eDeterministic :
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
				if  (ModelParameters::exo_det_nbr > 0)
					if (period1 == period2)
					{
						*output << "ex_det_(M_.maximum_lag + " << period1 <<
							", " << id1+1 << ") = " << mValues[i] << ";\n";
					}
					else
					{
						*output << "ex_det_(" << period1 << ":" << period2 <<
						", " << id1+1 << ") = repmat(" << mValues[i] << ", " <<
						period2 << "-" << period1 << "+1, 1);\n";
					}
				if (ModelParameters::exo_det_nbr == 0)
					if (period1 == period2)
					{
						*output << "oo_.exo_simul(M_.maximum_lag+" << period1 << ", " << id1+1 <<
						 ") = " <<  mValues[i] << ";\n";
					}
					else
					{
						*output << "oo_.exo_simul(" << period1 << ":" << period2 << ", " <<
							id1+1 << ") = repmat(" << mValues[i] << ", " << period2 << "-" << 
							period1 << "+1,1);\n";
					}
			}
			break;
		case eSTDerror		:
			*output << "M_.Sigma_e(" << id1+1 << ", " << id1+1 << ") = " << value << "^2;\n";
			break;
		case eVariance 		:
			*output << "M_.Sigma_e(" << id1+1 << ", " << id1+1 << ") = " << value << ";\n";
			break;
		case eCovariance	:
			*output << "M_.Sigma_e(" << id1+1 << ", " << id2+1 << ") = " <<	value <<
				"; M_.Sigma_e(" << id2+1 << ", " << id1+1 << ") = M_.Sigma_e(" << 
				id1+1 << ", " << id2+1 << ");\n";
			break;
		case eCorrelation	:
			*output << "M_.Sigma_e(" << id1+1 << ", " << id2+1 << ") = " <<
				value << "*sqrt(M_.Sigma_e(" << id1 << ", " << id1+1 << ")*M_.Sigma_e(" <<
				id2 << ", " << id2+1 << "); M_.Sigma_e(" << id2+1 << ", " << 
				id1 << ") = M_.Sigma_e(" << id1+1 << ", " << id2+1 << ");\n";
			break;
		default :
			;
	}		
	mPeriod1.clear();
	mPeriod2.clear();
	mValues.clear();

/*
	sets a shock 
for type == 0, set type, id1, shock_elem will be set by another method
for type == 1 or 2 set type, id1 and value
for type == 3 or 4 set type, id1, id2 and value
*/
}
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
