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
TmpSymbolTable::TmpSymbolTable()
{
  // Empty
}

TmpSymbolTable::~TmpSymbolTable()
{
  // Empty
}

void TmpSymbolTable::setGlobalSymbolTable(SymbolTable *symbol_table_arg)
{
  symbol_table = symbol_table_arg;
}

void TmpSymbolTable::AddTempSymbol(string symbol)
{
  if (symbol_table->Exist(symbol))
    tmpsymboltable.push_back(symbol);
  else
    {
      string msg = "Unknown symbol : "+symbol;
      error(msg.c_str());
    }
}

//------------------------------------------------------------------------------
void  TmpSymbolTable::AddTempSymbol(string symbol1, string symbol2)
{
  if (symbol_table->Exist(symbol1))
    tmpsymboltable.push_back(symbol1);
  else
    {
      string msg = "Unknown symbol : "+symbol1;
      error(msg.c_str());
    }

  if (symbol_table->Exist(symbol2))
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
      if (symbol_table->isReferenced(*it) == eReferenced)
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
