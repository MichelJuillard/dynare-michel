/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the VariableTable class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
//------------------------------------------------------------------------------
#include "VariableTable.h"
//------------------------------------------------------------------------------
map<varKey,int> VariableTable::mVariableTable = *(new map<varKey,int>);
vector<varKey> VariableTable::mVariableIndex = *(new vector<varKey>);	
vector<int>	VariableTable::mSortedVariableID = *(new vector<int>);
vector<int>	VariableTable::mPrintFormatIndex = *(new vector<int>);
void (* VariableTable::error) (const char* ) = NULL;
//------------------------------------------------------------------------------
VariableTable::VariableTable()
{ 
	// Empty
}
//------------------------------------------------------------------------------
VariableTable::~VariableTable()
{
	// Empty
}
//------------------------------------------------------------------------------
int VariableTable::AddVariable(string iName, int iLag)
{
	int 		lVariableID;	
	//Variable 	lVariable;
	varKey			key;
	// Testing if symbol exists 
	if (!SymbolTable::Exist(iName))
    {
      		string msg = "unknown symbol: " + iName;
	   		(* error) (msg.c_str());
      		exit(-1);
    }
    // Testing if variable exists in VaiableTable
    lVariableID = getID(iName,iLag);
    if (lVariableID != -1)
    {
    	return lVariableID;
    }
    lVariableID = mVariableIndex.size();
    // Making variable struct and key
    //lVariable.type = SymbolTable::getType(iName);
    //lVariable.symbol_id = SymbolTable::getID(iName);
    //lVariable.variable_id = lVariableID;
    key = make_pair(iName,iLag);
    // Pushing variable on VariableTable 
  	//mVariableTable[key] = lVariable;
  	mVariableTable[key] = lVariableID;
  	mVariableIndex.push_back(key);
  	// Setting lags field in symbol table
  	SymbolTable::AddLag(iName, iLag);
  	// Setting variable numbers
  	Type type = getType(lVariableID);
  	if (type == eEndogenous)	ModelParameters::var_endo_nbr++;
  	if (type == eExogenous)	ModelParameters::var_exo_nbr++;
  	// Setting Maximum and minimum lags
	if (ModelParameters::max_lead < iLag && iLag >= 0)
			ModelParameters::max_lead = iLag;
	if (ModelParameters::max_lag < -iLag && iLag <= 0)
			ModelParameters::max_lag = -iLag;
  	return mVariableIndex.size()-1;
}
//------------------------------------------------------------------------------
void VariableTable::decSymbolID(string iName, int id, int iLag, Type iType)
{
	int 		lVariableID;	
	Variable 	lVariable;
	varKey			key;

    // Testing if variable exists in VaiableTable
    lVariableID = getID(iName,iLag);
    if (lVariableID == -1)
    {
    	return;
    }
    // Making variable struct and key to update
    lVariable.type = iType;
    lVariable.symbol_id = id-1;
    lVariable.variable_id = lVariableID;
    // Updating VariableTable with new variable
    key = make_pair(iName,iLag);
  	mVariableTable[key] = lVariableID;

}
//------------------------------------------------------------------------------
void	VariableTable::Sort(void)
{
	varKey 					key;
	int 					variable;
	vector<pair<unsigned long long int,int> >	VarToSort; // To store variable lags
	vector<int>				IDs;
	vector<int>				Lags;
	vector<Type>			Types;
	
	if (mVariableIndex.size() == 1)
	{
		mSortedVariableID.push_back(0);
		mPrintFormatIndex.push_back(0);
		return;
	}
	// First putting types into TypesToSort
	for (int id=0; id < mVariableIndex.size(); id++)
	{
		key = mVariableIndex[id];
		variable = mVariableTable[key];
		//IDs.push_back(variable.symbol_id);
		//Types.push_back(variable.type);
		IDs.push_back(getSymbolID(variable));
		Types.push_back(getType(variable));
		Lags.push_back(key.second);
		unsigned long long int lag = Lags[id]+ModelParameters::max_lag;
		lag = lag << (4*sizeof(int));
		unsigned long long int  type = Types[id];
		type = type <<  8*sizeof(int);
		unsigned long long int sort_pound = IDs[id]+lag+type;
		VarToSort.push_back(make_pair(sort_pound,id));	
	}
	// Uncomment this to debug
	/*
	cout << "Before sorting\n";
	cout << "S T L ID  pound \n";
	for (int id=0; id < VarToSort.size(); id++)
	{	
		Type type = Types[VarToSort[id].second];
		int lag = Lags[VarToSort[id].second];
		int ID = IDs[VarToSort[id].second];
		cout << SymbolTable::getNameByID(type, ID) << " " 
			 << type << " " 
			 << lag <<  " " 
			 << ID << " "
			 << VarToSort[id].first << "\n";
	}
	*/
	// Sorting variables
	sort(VarToSort.begin(), VarToSort.end());
	// Puting "sorted Ids" into mSortedVariableID
	mSortedVariableID.resize(VarToSort.size());
	mPrintFormatIndex.resize(VarToSort.size());
	Type type = Types[VarToSort[0].second];
	int index = 0;
	for (int id = 0; id < VarToSort.size(); id++)
	{	
		int id2 = VarToSort[id].second;
		mSortedVariableID[id2] = id;
		if (type == Types[id2])
		{
			mPrintFormatIndex[id2] = index;
			index++;
		}
		else
		{						
			mPrintFormatIndex[id2] = 0;
			type = Types[id2];
			index = 1;
		}
	}
	// Uncomment this to debug
	/*
	cout << "After sorting\n";
	cout << "S T L ID   SVID PIDX\n";
	for (int id=0; id < VarToSort.size(); id++)
	{	
		Type type = Types[VarToSort[id].second];
		int lag = Lags[VarToSort[id].second];
		int ID = IDs[VarToSort[id].second];
		cout << SymbolTable::getNameByID(Types[id], IDs[id]) << " " 
			 << Types[id] << " " 
			 << Lags[id] << " " 
			 << IDs[id] << " "
			 << mSortedVariableID[id] << " "
			 << mPrintFormatIndex[id] << "\n";
	}
	*/
}
//------------------------------------------------------------------------------
