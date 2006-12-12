#ifndef _TMPSYMBOLTABLE_HH
#define _TMPSYMBOLTABLE_HH
//------------------------------------------------------------------------------
/** \file
 * \version 1.0
 * \date 04/26/2004
 * \par This file defines the TmpSymbolTable class.
 */
//------------------------------------------------------------------------------
#include <string>
#include <vector>
#include <ostream>

#include "SymbolTable.hh"

/*!
  \class  TmpSymbolTable
  \brief  Defines temparary symbol table used with computing tasks
*/
class TmpSymbolTable
{
private :
  /*! list of string TempSymbolTable */
  std::vector<std::string> tmpsymboltable;
  /*! List of symbol Values */
  std::vector<std::string> nameTable;
  //! A reference to enclosing symbol table
  const SymbolTable &symbol_table;
public :
  /*! Constrcutor */
  TmpSymbolTable(const SymbolTable &symbol_table_arg);
  /*! Destructor*/
  ~TmpSymbolTable();
  /*! Adds a temp symbol */
  void AddTempSymbol(const std::string &symbol);
  /*! Adds a temp symbol and its value */
  void AddTempSymbol(const std::string &symbol1, const std::string &symbol2);
  /*! Write TempSymbolTable to output string */
  void writeOutput(const std::string &varname, std::ostream &output) const;
  //! Clears all content
  void clear();
};
//------------------------------------------------------------------------------
#endif
