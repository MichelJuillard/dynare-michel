#ifndef _MOD_FILE_HH
#define _MOD_FILE_HH

using namespace std;

#include <ostream>

#include "ModelParameters.hh"
#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "ModelTree.hh"
#include "VariableTable.hh"
#include "Statement.hh"

//! The abstract representation of a "mod" file
class ModFile
{
public:
  ModFile();
  ~ModFile();
  //! Model parameters
  ModelParameters model_parameters;
  //! Symbol table
  SymbolTable symbol_table;
  //! Numerical constants table
  NumericalConstants num_constants;
  //! Model equations and their derivatives
  ModelTree model_tree;
  //! Option order
  int order;
  //! Option linear
  int linear;
private:
  //! List of statements
  vector<Statement *> statements;

public:
  //! Add a statement
  void addStatement(Statement *st);
  //! Writes Matlab/Scilab output files
  /*!
    \param basename The base name used for writing output files. Should be the name of the mod file without its extension
    \param clear_all Should a "clear all" instruction be written to output ?
  */
  void writeOutputFiles(const string &basename, bool clear_all);
};

#endif // ! MOD_FILE_HH
