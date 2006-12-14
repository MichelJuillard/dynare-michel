#ifndef _VARIABLETABLE_HH
#define _VARIABLETABLE_HH

using namespace std;

#include <map>
#include <string>
#include <vector>

#include "SymbolTable.hh"

//! This class is used to store variables as they appear in the model (with their lead or lag)
/*! \todo Raise exceptions when requesting ordered IDs before calling to Sort() */
class VariableTable
{
private:
  //! A reference to the symbol table
  const SymbolTable &symbol_table;
  //! A reference to model parameters
  ModelParameters &mod_param;
  //! Variable key type to acced variable table elements
  typedef pair<string, int> varKey;
  //! Maps a pair (symbol, lag) to an ID
  map<varKey, int> mVariableTable;
  //! Maps an ID to a pair (symbol, lag)
  /*! It is the reverse map of mVariableTable */
  vector<varKey> mVariableIndex;
  //! Variable IDs of sorted variable table
  vector<int> mSortedVariableID;
  //! For each variable, gives its index number among variables of the same type
  /*! It is the index used in the output file:
      - in the lead/lag matrix
      - in the right hand side of equations (such as y(index))
  */
  vector<int> mPrintFormatIndex;
public:
  VariableTable(const SymbolTable &symbol_table_arg, ModelParameters &mod_param_arg);
  //! Adds a variable in the table, and returns its (newly allocated) varID
  /*! Also works if the variable already exists */
  int AddVariable(const string &iName, int iLag);
  //! Return variable ID
  inline int getID(const string &iName, int iLag) const;
  //! Return lag of variable
  inline int getLag(int ivarID) const;
  //! Return symbol ID of variable
  inline int getSymbolID(int ivarID) const;
  //! Get variable type
  inline Type getType(int ivarID) const;
  //! Get number of variables in mVariableTable
  inline int size() const;
  //! Get variable ID of sorted variable table
  inline int getSortID(int iVarID) const;
  //! Return variable index to print in format : y(index) or oo_.y_simul(index) ...
  inline int getPrintIndex(int iVarID) const;
  //! Sorts variable table
  /*! The order used is a lexicographic order over the tuple (type, lag, symbolID) */
  void Sort();
};

inline int
VariableTable::getSortID(int iVarID) const
{
  return mSortedVariableID[iVarID];
}

inline int
VariableTable::getPrintIndex(int iVarID) const
{
  return mPrintFormatIndex[iVarID];
}

inline int
VariableTable::getID(const string &iName, int iLag) const
{
  map<varKey, int>::const_iterator it = mVariableTable.find(make_pair(iName, iLag));
  if (it == mVariableTable.end())
    return -1;
  else
    return it->second;
}

inline Type
VariableTable::getType(int ivarID) const
{
  varKey key = mVariableIndex[ivarID];
  return symbol_table.getType(key.first);
}

inline int
VariableTable::getSymbolID(int ivarID) const
{
  varKey key = mVariableIndex[ivarID];
  return symbol_table.getID(key.first);
}

inline int
VariableTable::getLag(int ivarID) const
{
  return mVariableIndex[ivarID].second;
}

inline int
VariableTable::size() const
{
  return mVariableTable.size();
}

#endif
