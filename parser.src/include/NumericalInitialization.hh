#ifndef _NUMERICALINITIALIZATION_HH
#define _NUMERICALINITIALIZATION_HH

using namespace std;

#include <string>
#include <map>

#include "SymbolTable.hh"
#include "Statement.hh"

class InitParamStatement : public Statement
{
private:
  const string param_name;
  const string param_value;
  const SymbolTable &symbol_table;
public:
  InitParamStatement(const string &param_name_arg, const string &param_value_arg,
                     const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class InitOrEndValStatement : public Statement
{
public:
  /*!
    We use a vector instead of a map, since the order of declaration matters:
    an initialization can depend on a previously initialized variable inside the block
  */
  typedef vector<pair<string, string> > init_values_type;
protected:
  const init_values_type init_values;
  const SymbolTable &symbol_table;
public:
  InitOrEndValStatement(const init_values_type &init_values_arg,
                        const SymbolTable &symbol_table_arg);
protected:
  void writeInitValues(ostream &output) const;
};

class InitValStatement : public InitOrEndValStatement
{
private:
  const ModelParameters &mod_param;
public:
  InitValStatement(const init_values_type &init_values_arg,
                   const SymbolTable &symbol_table_arg,
                   const ModelParameters &mod_param_arg);
  virtual void writeOutput(ostream &output) const;
};

class EndValStatement : public InitOrEndValStatement
{
public:
  EndValStatement(const init_values_type &init_values_arg,
                  const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class HistValStatement : public Statement
{
public:
  /*!
    Contrary to Initval and Endval, we use a map, since it is impossible to reuse
    a given initialization value in a second initialization inside the block.
  */
  typedef map<pair<string, int>, string> hist_values_type;
private:
  const hist_values_type hist_values;
  const SymbolTable &symbol_table;
public:
  HistValStatement(const hist_values_type &hist_values_arg,
                   const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

#endif
