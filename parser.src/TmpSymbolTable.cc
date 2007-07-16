/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the TmpSymbolTable class methodes.
*/
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
#include "SymbolTable.hh"
#include "TmpSymbolTable.hh"
#include "Interface.hh"
//------------------------------------------------------------------------------
TmpSymbolTable::TmpSymbolTable(const SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg)
{
  // Empty
}

TmpSymbolTable::~TmpSymbolTable()
{
  // Empty
}

void
TmpSymbolTable::AddTempSymbol(const string &symbol)
{
  // FIXME: add check to verify that symbol exists in symbol_table
  // FIXME: add check to verify that symbol doesn't yet exist in the present table
  tmpsymboltable.push_back(symbol);
}

void
TmpSymbolTable::AddTempSymbol(const string &symbol1, const string &symbol2)
{
  // FIXME: add checks to verify that symbol1 and symbol2 exist in symbol_table
  // FIXME: add check to verify that symbol1 doesn't yet exist in the present table
  tmpsymboltable.push_back(symbol1);
  nameTable.push_back(symbol2);
}

void
TmpSymbolTable::writeOutput(const string &varname, ostream &output) const
{
  output << varname << "=[];" << endl;
  for (vector<string>::const_iterator it = tmpsymboltable.begin();
       it != tmpsymboltable.end(); it++)
    {
      output << varname << " = ";
      output << interfaces::strvcat(varname, "'" + *it + "'") << ";" << endl;
    }
}

void
TmpSymbolTable::clear()
{
  tmpsymboltable.clear();
  nameTable.clear();
}
