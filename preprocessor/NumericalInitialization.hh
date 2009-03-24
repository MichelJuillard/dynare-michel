/*
 * Copyright (C) 2003-2009 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NUMERICALINITIALIZATION_HH
#define _NUMERICALINITIALIZATION_HH

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "SymbolTable.hh"
#include "ExprNode.hh"
#include "Statement.hh"

class InitParamStatement : public Statement
{
private:
  const int symb_id;
  const NodeID param_value;
  const SymbolTable &symbol_table;
public:
  InitParamStatement(int symb_id_arg, const NodeID param_value_arg,
                     const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
  //! Fill eval context with parameter value
  void fillEvalContext(eval_context_type &eval_context) const;
};

class InitOrEndValStatement : public Statement
{
public:
  /*!
    We use a vector instead of a map, since the order of declaration matters:
    an initialization can depend on a previously initialized variable inside the block
  */
  typedef vector<pair<int, NodeID> > init_values_type;
protected:
  const init_values_type init_values;
  const SymbolTable &symbol_table;
public:
  InitOrEndValStatement(const init_values_type &init_values_arg,
                        const SymbolTable &symbol_table_arg);
  //! Fill eval context with variables values
  void fillEvalContext(eval_context_type &eval_context) const;
protected:
  void writeInitValues(ostream &output) const;
};

class InitValStatement : public InitOrEndValStatement
{
public:
  InitValStatement(const init_values_type &init_values_arg,
                   const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class EndValStatement : public InitOrEndValStatement
{
public:
  EndValStatement(const init_values_type &init_values_arg,
                  const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class HistValStatement : public Statement
{
public:
  /*!
    Contrary to Initval and Endval, we use a map, since it is impossible to reuse
    a given initialization value in a second initialization inside the block.
    Maps pairs (symbol_id, lag) to NodeID
  */
  typedef map<pair<int, int>, NodeID> hist_values_type;
private:
  const hist_values_type hist_values;
  const SymbolTable &symbol_table;
public:
  HistValStatement(const hist_values_type &hist_values_arg,
                   const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class InitvalFileStatement : public Statement
{
private:
  const string filename;
public:
  InitvalFileStatement(const string &filename_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class HomotopyStatement : public Statement
{
public:
  //! Stores the declarations of homotopy_setup
  /*! Order matter so we use a vector. First NodeID can be NULL if no initial value given. */
  typedef vector<pair<int, pair<NodeID, NodeID> > > homotopy_values_type;
private:
  const homotopy_values_type homotopy_values;
  const SymbolTable &symbol_table;
public:
  HomotopyStatement(const homotopy_values_type &homotopy_values_arg,
                    const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class SaveParamsAndSteadyStateStatement : public Statement
{
private:
  const string filename;
public:
  SaveParamsAndSteadyStateStatement(const string &filename_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class LoadParamsAndSteadyStateStatement : public Statement
{
private:
  const SymbolTable &symbol_table;
  //! Content of the file
  /*! Maps symbol ID to numeric value (stored as string) */
  map<int, string> content;
public:
  LoadParamsAndSteadyStateStatement(const string &filename,
                                    const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
  //! Fill eval context with parameters/variables values
  void fillEvalContext(eval_context_type &eval_context) const;
};

#endif