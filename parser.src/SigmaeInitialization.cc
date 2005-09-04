/*! \file 
 \version 1.0
 \date 04/13/2004
 \par This file defines the SigmaeInitialization class.
*/
//------------------------------------------------------------------------------
#include <iostream>

using namespace std;

#include "SigmaeInitialization.h"
//------------------------------------------------------------------------------
SigmaeInitialization::SigmaeInitialization()
{
	// Empty
}
//------------------------------------------------------------------------------
SigmaeInitialization::~SigmaeInitialization()
{
	// Empty
}
//------------------------------------------------------------------------------
void SigmaeInitialization::setOutput(ostringstream* iOutput)
{
	output = iOutput;	
}
//------------------------------------------------------------------------------
void SigmaeInitialization::AddExpression(string expression)
{
	row.push_back(expression);
}
//------------------------------------------------------------------------------
void SigmaeInitialization::EndOfRow()
{
	
	matrix.push_back(row);
	row.clear();
}
//------------------------------------------------------------------------------
void SigmaeInitialization::CheckMatrix(void)
{
	vector<vector<string> >::iterator ir;
	unsigned int nbe;
	int inc;
	// Checking if first or last row has one element. 
	if (matrix.front().size() == 1)
	{
		inc = 1;
		nbe = 2;
		type = eLower;		
	}
	else if (matrix.back().size() == 1)
	{
		inc = -1;
		nbe = matrix.front().size()-1;
		type = eUpper;
	}
	else
	{
		string msg = "sigma_e isn't in triangular form";
   		(* error) (msg.c_str());
	}	
	// Checking if matrix is triangular (upper or lower):
	// each row has one element more or less than the previous one
	// and first or last one has one element. 	
	for (ir = matrix.begin(), ir++; ir != matrix.end(); ir++, nbe += inc )
		if (ir->size() != nbe)
		{
			string msg = "sigma_e isn't in triangular form!";
	   		(* error) (msg.c_str());
		}
}
//------------------------------------------------------------------------------
void SigmaeInitialization::SetMatrix(void)
{
	unsigned int ic, ic1;
	unsigned int ir, ir1;
	
	*output << "M_.Sigma_e = [...\n";
	for (ir = 0; ir < matrix.size(); ir++)
	{
		*output << "\t";
		for (ic = 0; ic < matrix.size(); ic++)
		{
			if (ic >= ir && type == eUpper)
			{
				ic1 = ic-ir;
				ir1 = ir;
			}
			else if (ic < ir && type == eUpper)
			{
				ic1 = ir-ic;
				ir1 = ic;
			}
			else if (ic > ir && type == eLower)
			{
				ic1 = ir;
				ir1 = ic;
			}
			else if (ic <= ir && type == eLower)
			{
				ic1 = ic;
				ir1 = ir;
			}

			*output << matrix[ir1][ic1] << " ";
		}
		*output << ";...\n";
	}
	*output << "\t];...\n";
}
//------------------------------------------------------------------------------
void SigmaeInitialization::set(void)
{
	CheckMatrix();
	
	SetMatrix();		
	matrix.clear();
	row.clear();
}
//------------------------------------------------------------------------------
