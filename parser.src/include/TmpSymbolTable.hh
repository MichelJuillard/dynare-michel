#ifndef _TMPSYMBOLTABLE_HH
#define _TMPSYMBOLTABLE_HH
//------------------------------------------------------------------------------
/** \file
 * \version 1.0
 * \date 04/26/2004
 * \par This file defines the TmpSymbolTable class.
 */
//------------------------------------------------------------------------------
#include <list>
#include <sstream>

#include "SymbolTable.hh"

/*!
  \class  TmpSymbolTable
  \brief  Defines temparary symbol table used with computing tasks
*/
class TmpSymbolTable
{
private :
  /*! list of string TempSymbolTable */
  std::list<std::string>  tmpsymboltable;
  /*! List of symbol Values */
  std::list<std::string>  NameTable;
  /*! Output of this class */
  std::ostringstream  output;
  //! Pointer to global symbol table
  SymbolTable *symbol_table;
public :
  /*! Constrcutor */
  TmpSymbolTable();
  /*! Destructor*/
  ~TmpSymbolTable();
  //! Sets global symbol table pointer
  void setGlobalSymbolTable(SymbolTable *symbol_table_arg);
  /*! Pointer to error function of parser class */
  void (* error) (const char* m);
  /*! Adds a temp symbol */
  void  AddTempSymbol(std::string symbol);
  /*! Adds a temp symbol and its value */
  void  AddTempSymbol(std::string symbol1, std::string symbol2);
  /*! Write TempSymbolTable to output string */
  void  set(std::string varname);
  /*! Gets size of TempSymbolTable */
  int   size(void);
  /*! Gets output of this class*/
  std::string  get(void);
};
//------------------------------------------------------------------------------
#endif
