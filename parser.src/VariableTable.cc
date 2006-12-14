#include <iostream>
#include <algorithm>

#include "VariableTable.hh"

VariableTable::VariableTable(const SymbolTable &symbol_table_arg,
                             ModelParameters &mod_param_arg) :
  symbol_table(symbol_table_arg),
  mod_param(mod_param_arg)
{
}

int
VariableTable::AddVariable(const string &iName, int iLag)
{
  // Testing if symbol exists
  if (!symbol_table.Exist(iName))
    {
      cerr << "Unknown symbol: " << iName << endl;
      exit(-1);
    }
  // Testing if variable exists in VariableTable
  int lVariableID = getID(iName,iLag);
  if (lVariableID != -1)
    return lVariableID;

  lVariableID = mVariableIndex.size();
  // Making variable struct and key
  varKey key = make_pair(iName, iLag);
  // Pushing variable on VariableTable
  mVariableTable[key] = lVariableID;
  mVariableIndex.push_back(key);

  // Setting dynamic variables numbers
  Type type = getType(lVariableID);
  if (type == eEndogenous)
    mod_param.var_endo_nbr++;
  if (type == eExogenous)
    mod_param.var_exo_nbr++;
  if (type == eExogenousDet)
    mod_param.var_exo_det_nbr++;

  // Setting maximum and minimum lags
  if (mod_param.max_lead < iLag)
    mod_param.max_lead = iLag;
  else if (-mod_param.max_lag > iLag)
    mod_param.max_lag = -iLag;

  switch(type)
    {
    case eEndogenous:
      if (mod_param.max_endo_lead < iLag)
        mod_param.max_endo_lead = iLag;
      else if (-mod_param.max_endo_lag > iLag)
        mod_param.max_endo_lag = -iLag;
      break;
    case eExogenous:
      if (mod_param.max_exo_lead < iLag)
        mod_param.max_exo_lead = iLag;
      else if (-mod_param.max_exo_lag > iLag)
        mod_param.max_exo_lag = -iLag;
      break;
    case eExogenousDet:
      if (mod_param.max_exo_det_lead < iLag)
        mod_param.max_exo_det_lead = iLag;
      else if (-mod_param.max_exo_det_lag > iLag)
        mod_param.max_exo_det_lag = -iLag;
      break;
    case eRecursiveVariable:
      if (mod_param.max_recur_lead < iLag)
        mod_param.max_recur_lead = iLag;
      else if (-mod_param.max_recur_lag > iLag)
        mod_param.max_recur_lag = -iLag;
      break;
    default:
      ;
    }
  return lVariableID;
}

void
VariableTable::Sort()
{
  // Trivial case where no ordering is necessary
  if (mVariableIndex.size() == 1)
    {
      mSortedVariableID.push_back(0);
      mPrintFormatIndex.push_back(0);
      return;
    }

  /* The type of key for lexicographic ordering over variables:
     the key is equal to (type, lag, symbol_id) */
  typedef pair<Type, pair<int, int> > lexicographic_key_type;

  // Construct the vector matching keys to their varIDs
  vector<pair<lexicographic_key_type, int> > VarToSort;
  for (unsigned int varID = 0; varID < mVariableIndex.size(); varID++)
    {
      Type type = getType(varID);
      int symbolID = getSymbolID(varID);
      int lag = mVariableIndex[varID].second;
      VarToSort.push_back(make_pair(make_pair(type, make_pair(lag, symbolID)), varID));
    }

  // Sort variables using the lexicographic ordering
  sort(VarToSort.begin(), VarToSort.end());

  // Fill mSortedVariableID and mPrintFormatIndex
  mSortedVariableID.resize(VarToSort.size());
  mPrintFormatIndex.resize(VarToSort.size());
  Type type = getType(VarToSort[0].second);
  int index = 0;
  for (unsigned int sortedID = 0; sortedID < VarToSort.size(); sortedID++)
    {
      int varID = VarToSort[sortedID].second;
      mSortedVariableID[varID] = sortedID;
      if (type == getType(varID))
        {
          mPrintFormatIndex[varID] = index;
          index++;
        }
      else
        {
          mPrintFormatIndex[varID] = 0;
          type = getType(varID);
          index = 1;
        }
    }
}
