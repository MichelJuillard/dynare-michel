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
#include <ostream>

#include "SymbolTableTypes.hh"

/*!
  \class SymbolTable
  \brief This class keeps track of symbols
*/
class SymbolTable
{
private:
  /*! Adds symbol into symbol table
    \param name a string.
    \param type a Type struct.
    \param tex_name a string for the TeX name.
    \par Description
    - warning if symbol is already set with same type
    - error if symbol is already set with different type
    - set Name and Type
    - increase corresponding counter in ModelParameters class
  */
  int AddSymbol(std::string name, Type type, std::string tex_name);
  /*! Symbol table map */
  std::map<std::string, Symbol, std::less<std::string> > symboltable;
  //! Typedef for const iterator on symboltable map
  typedef std::map<std::string, Symbol, std::less<std::string> >::const_iterator symboltable_const_iterator;
  /*! Symbol name table indexed by type and ID */
  std::vector< std::vector<std::string> > name_table;
  std::vector< std::vector<std::string> > tex_name_table;
  /*! Changes type of a symbol */
  void ResetType(std::string name, Type new_type);
public :
  /*! Constructor */
  SymbolTable();
  //! Number of declared endogenous variables
  int endo_nbr;
  //! Number of declared exogenous variables
  int exo_nbr;
  //! Number of declared deterministic exogenous variables
  int exo_det_nbr;
  //! Number of declared parameters
  int parameter_nbr;
  //! Number of declared local parameters
  int local_parameter_nbr;
  //! Number of declared recursive variables
  int recur_nbr;
  /*! Pointer to error function of parser class */
  void (* error) (const char* m);
  /*! Adds a symbol apearing in declaration
    - warning if symbol is already set with same type
    - error if symbol is already set with different type
    - set name, type
    - increase corresponding counter in ModelParameters
  */
  int AddSymbolDeclar(std::string name, Type type, std::string tex_name);
  /*! Adds symbol range */
  void AddSymbolRange(std::string name, int nbr, Type type, std::string tex_name);
  /*! Sets a symbol as referenced */
  void SetReferenced(std::string name);
  /*! Return eReferenced if symbol is referenced eNotReferenced otherwise*/
  Reference isReferenced(const std::string &name) const;
  /*! Tests if symbol exists in symbol table
    \return true if exists, false outherwise
  */
  inline bool Exist(const std::string &name) const;
  /*! Gets name by type and ID */
  inline std::string getNameByID(Type type, int id);
  /*! Gets tex name by type and ID */
  inline std::string getTexNameByID(Type type, int id);
  /*! Gets type by name */
  inline Type getType(const std::string &name) const;
  /*! Gets ID by name */
  inline int getID(const std::string &name) const;
  //! Write output of this class
  void writeOutput(std::ostream &output);
};

inline bool SymbolTable::Exist(const std::string &name) const
{
  symboltable_const_iterator iter = symboltable.find(name);
  return (iter != symboltable.end());
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
inline Type SymbolTable::getType(const std::string &name) const
{
  symboltable_const_iterator iter = symboltable.find(name);
  if (iter == symboltable.end())
    return eUNDEF;
  else
    return iter->second.type;
}

//------------------------------------------------------------------------------
inline int SymbolTable::getID(const std::string &name) const
{
  symboltable_const_iterator iter = symboltable.find(name);
  if (iter == symboltable.end())
    return -1;
  else
    return(iter->second.id);
}

//------------------------------------------------------------------------------
#endif
