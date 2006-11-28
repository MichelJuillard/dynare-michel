#ifndef _MOD_FILE_HH
#define _MOD_FILE_HH

#include "ModelParameters.hh"
#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "NumericalInitialization.hh"
#include "Shocks.hh"
#include "SigmaeInitialization.hh"
#include "ComputingTasks.hh"
#include "ModelTree.hh"
#include "VariableTable.hh"

//! The abstract representation of a "mod" file
class ModFile
{
public:
  //! Constructor
  ModFile();
  //! Model parameters
  ModelParameters model_parameters;
  //! Symbol table
  SymbolTable symbol_table;
  //! Variable table
  VariableTable variable_table;
  //! Numerical constants
  NumericalConstants num_constants;
  //! Numerical initalisations
  NumericalInitialization numerical_initialization;
  //! Shocks instructions
  Shocks shocks;
  //! Sigma_e instructions
  SigmaeInitialization sigmae;
  //! Computing tasks instructions
  ComputingTasks computing_tasks;
  //! Model equations and their derivatives
  ModelTree model_tree;
  //! Option order
  int order;
  //! Option linear
  int linear;
};

#endif // ! MOD_FILE_HH
