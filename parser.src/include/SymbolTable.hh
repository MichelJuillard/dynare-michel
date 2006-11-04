#ifndef _SYMBOLTABLE_HH
#define _SYMBOLTABLE_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 12/01/2003
  \par This file defines the SymbolTable class .
*/
//------------------------------------------------------------------------------
#include <map>
#include <string>
#include <vector>
#include <sstream>
//------------------------------------------------------------------------------
#include "ModelParameters.hh"
#include "SymbolTableTypes.hh"
//------------------------------------------------------------------------------
/*!
  \class SymbolTable
  \brief This class keeps track of symbols
*/

class SymbolTable
{

private :
  static std::ostringstream output;
  /*! Adds symbol into symbol table
    \param name a string.
    \param type a Type struct.
    \par Description
    - warning if symbol is already set with same type \n
    - error if symbol is already set with different type\n
    - set Name and Type\n
    - increase corresponding counter in ModelParameters class\n
  */
  static int AddSymbol(std::string name,Type type, std::string tex_name);
protected :
  /*! Symbol table map */
  static std::map<std::string, Symbol, std::less<std::string> > symboltable;
  /*! Symbol name table indexed by type and ID */
  static std::vector< std::vector<std::string> >          name_table;
  static std::vector< std::vector<std::string> >          tex_name_table;
protected :
  /*! Changes type of a symbol */
  static  void  ResetType(std::string name,Type new_type);
public :
  /*! Constructor */
  SymbolTable();
  /*! Destructor*/
  ~SymbolTable();
  /*! Pointer to error function of parser class */
  static void (* error) (const char* m);
  /*! Adds a symbol apearing in declaration
    - warning if symbol is already set with same type
    - error if symbol is already set with different type
    - set name, type
    - increase corresponding counter in ModelParameters
  */
  static  int   AddSymbolDeclar(std::string name,Type type, std::string tex_name);
  /*! Adds symbol range */
  static  void  AddSymbolRange(std::string name,int nbr,Type type, std::string tex_name);
  /*! Adds a lag to field lags */
  static void AddLag(std::string name,int lag);
  /*! Sets a symbol as referenced */
  static  void  SetReferenced(std::string name);
  /*! Return eReferenced if symbol is referenced eNotReferenced otherwise*/
  static  Reference   isReferenced(std::string name);
  /*! Tests if symbol exists in symbol table
    \return true if exists, false outherwise
  */
  inline static   bool  Exist(std::string name);
  /*! Gets name by type and ID */
  inline static std::string getNameByID(Type type,int id);
  /*! Gets tex name by type and ID */
  inline static std::string getTexNameByID(Type type,int id);
  /*! Gets type by name */
  inline static Type  getType(std::string name);
  /*! Gets ID by name */
  inline static int   getID(std::string name);
  /*! Gets output string of this class */
  static std::string  get();
  /*! Checks if symbols are used in model equations, removes unused symbol */
  void clean();
  void erase_local_parameters();
};
inline bool SymbolTable::Exist(std::string name)
{
  std::map<std::string, Symbol, std::less<std::string> >::iterator iter;

  iter = symboltable.find(name);
  //Testing if symbol exists
  if (iter == symboltable.end()) return false;
  else return true;
}

//------------------------------------------------------------------------------
inline std::string  SymbolTable::getNameByID(Type type,int id)
{
  if (id >= 0 && (int)name_table[type].size() > id)
    return(name_table[type][id]);
  else return "";
}

//------------------------------------------------------------------------------
inline std::string  SymbolTable::getTexNameByID(Type type,int id)
{
  if (id >= 0 && (int)tex_name_table[type].size() > id)
    return(tex_name_table[type][id]);
  else return "";
}

//------------------------------------------------------------------------------
inline Type SymbolTable::getType(std::string name)
{
  if (Exist(name))
    return(symboltable[name].type);
  else
    return eUNDEF;
}

//------------------------------------------------------------------------------
inline int  SymbolTable::getID(std::string name)
{
  if (Exist(name))
    return(symboltable[name].id);
  else
    return -1;
}

//------------------------------------------------------------------------------
#endif
