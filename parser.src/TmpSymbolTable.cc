/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the TmpSymbolTable class methodes.
*/
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
#include "SymbolTable.h"
#include "TmpSymbolTable.h"
#include "Interface.h"
//------------------------------------------------------------------------------
TmpSymbolTable::TmpSymbolTable()
{
	// Empty
}
//------------------------------------------------------------------------------
TmpSymbolTable::TmpSymbolTable(const TmpSymbolTable &tst)
{
	tmpsymboltable = tst.tmpsymboltable;
	output.str(tst.output.str());
}
//------------------------------------------------------------------------------
TmpSymbolTable::~TmpSymbolTable()
{
	// Empty
}
//------------------------------------------------------------------------------
void TmpSymbolTable::AddTempSymbol(string symbol)
{
	if (SymbolTable::Exist(symbol))
		tmpsymboltable.push_back(symbol);
	else
	{
		string msg = "Unknown symbol : "+symbol;
		error(msg.c_str());
	}
}
//------------------------------------------------------------------------------
void 	TmpSymbolTable::AddTempSymbol(string symbol1, string symbol2)
{
	if (SymbolTable::Exist(symbol1))
		tmpsymboltable.push_back(symbol1);
	else
	{
		string msg = "Unknown symbol : "+symbol1;
		error(msg.c_str());
	}

	if (SymbolTable::Exist(symbol2))
		NameTable.push_back(symbol2);
	else
	{
		string msg = "Unknown symbol : "+symbol2;
		error(msg.c_str());
	}
}
//------------------------------------------------------------------------------
void TmpSymbolTable::set(string varname)				
{
	list<string>::iterator it;
	output.str("");
	output << "\n" << varname << "=[];\n";
	for (it = tmpsymboltable.begin(); it!= tmpsymboltable.end(); it++)
	{
		if (SymbolTable::isReferenced(*it) == eReferenced)
		{
		  output << varname << " = ";
		  output << interfaces::strvcat(varname,"'"+*it+"'")+";\n";			
		}
	}
}
//------------------------------------------------------------------------------
string TmpSymbolTable::get(void)
{
	tmpsymboltable.clear();
	NameTable.clear();
	return output.str();
}
//------------------------------------------------------------------------------
int TmpSymbolTable::size(void)
{
	return tmpsymboltable.size();
}
//------------------------------------------------------------------------------
