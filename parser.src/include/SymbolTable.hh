#ifndef _SYMBOLTABLE_HH
#define _SYMBOLTABLE_HH

using namespace std;

#include <map>
#include <string>
#include <vector>
#include <ostream>
#include <iostream>

//! Enumeration of possible symbol types
/*! Be careful not to change the order of the enumeration, it matters for VariableTable (at least up to eParameter) */
enum Type
  {
    eEndogenous = 0,               //!< Endogenous
    eExogenous = 1,                //!< Exogenous
    eExogenousDet = 2,             //!< Exogenous deterministic
    eRecursiveVariable = 3,        //!< Recursive variable (reserved for future use)
    eParameter = 4,                //!< Parameter
    eModelLocalVariable = 10,      //!< Local variable whose scope is model (pound expression)
    eModFileLocalVariable = 11,    //!< Local variable whose scope is mod file (model excluded)
    eUnknownFunction = 12          //!< Function unknown to the preprocessor
  };

//! Stores the symbol table
/*!
    A symbol is given by its name, and is internally represented by a pair (type, id).

    There is a distinct sequence of ids for each type, so two symbol of different types can have the same id.

    Also manages a TeX name for each symbol, which by default is an empty string.
*/
class SymbolTable
{
private:
  //! A symbol is represented by a pair (type, id)
  typedef pair<Type, int> symbol_type;

  //! Type for map: symbol_name -> (type, id)
  typedef map<string, symbol_type> symbol_table_type;
  //! Maps strings to pairs (type,id)
  symbol_table_type symbol_table;

  //! Type for map: (type, id) -> symbol_name
  typedef map<symbol_type, string> inv_symbol_table_type;

  //! Maps pairs (type, id) to names
  inv_symbol_table_type name_table;
  //! Maps pairs (type, id) to TeX names
  inv_symbol_table_type tex_name_table;
public:
  SymbolTable();
  //! Thrown when trying to access an unknown symbol (by name)
  class UnknownSymbolNameException
  {
  public:
    //! Symbol name
    string name;
    UnknownSymbolNameException(const string &name_arg) : name(name_arg) {}
  };
  //! Thrown when trying to access an unknown symbol (by type+id pair)
  class UnknownSymbolIDException
  {
  public:
    //! Symbol type
    Type type;
    //! Symbol ID
    int id;
    UnknownSymbolIDException(Type type_arg, int id_arg) : type(type_arg), id(id_arg) {}
  };
  //! Thrown when trying to declare a symbol twice
  class AlreadyDeclaredException
  {
  public:
    //! Symbol name
    string name;
    //! Was the previous declaration done with the same symbol type ?
    bool same_type;
    AlreadyDeclaredException(const string &name_arg, bool same_type_arg) : name(name_arg), same_type(same_type_arg) {}
  };
  //! Number of declared endogenous variables
  int endo_nbr;
  //! Number of declared exogenous variables
  int exo_nbr;
  //! Number of declared deterministic exogenous variables
  int exo_det_nbr;
  //! Number of declared recursive variables
  int recur_nbr;
  //! Number of declared parameters
  int parameter_nbr;
  //! Number of model local variables
  int model_local_variable_nbr;
  //! Number of modfile local variables
  int modfile_local_variable_nbr;
  //! Number of unknown functions
  int unknown_function_nbr;
  //! Add a symbol
  void addSymbol(const string &name, Type type, const string &tex_name = "") throw (AlreadyDeclaredException);
  //! Tests if symbol already exists
  inline bool exists(const string &name) const;
  //! Get symbol name by type and ID
  inline string getNameByID(Type type, int id) const throw (UnknownSymbolIDException);
  //! Get TeX name by type and ID
  inline string getTeXNameByID(Type type, int id) const throw (UnknownSymbolIDException);
  //! Get type by name
  inline Type getType(const string &name) const throw (UnknownSymbolNameException);
  //! Get ID by name
  inline int getID(const string &name) const throw (UnknownSymbolNameException);
  //! Write output of this class
  void writeOutput(ostream &output) const;

};

inline bool
SymbolTable::exists(const string &name) const
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  return (iter != symbol_table.end());
}

inline string
SymbolTable::getNameByID(Type type, int id) const throw (UnknownSymbolIDException)
{
  inv_symbol_table_type::const_iterator iter = name_table.find(make_pair(type, id));
  if (iter != name_table.end())
    return iter->second;
  else
    throw UnknownSymbolIDException(type, id);
}

inline string
SymbolTable::getTeXNameByID(Type type, int id) const throw (UnknownSymbolIDException)
{
  inv_symbol_table_type::const_iterator iter = tex_name_table.find(make_pair(type, id));
  if (iter != tex_name_table.end())
    return iter->second;
  else
    throw UnknownSymbolIDException(type, id);
}

inline Type
SymbolTable::getType(const string &name) const throw (UnknownSymbolNameException)
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  if (iter != symbol_table.end())
    return iter->second.first;
  else
    throw UnknownSymbolNameException(name);
}

inline int
SymbolTable::getID(const string &name) const throw (UnknownSymbolNameException)
{
  symbol_table_type::const_iterator iter = symbol_table.find(name);
  if (iter != symbol_table.end())
    return iter->second.second;
  else
    throw UnknownSymbolNameException(name);
}

#endif
