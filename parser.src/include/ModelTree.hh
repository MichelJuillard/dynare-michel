#ifndef _MODELTREE_HH
#define _MODELTREE_HH

using namespace std;

#include <string>
#include <vector>
#include <map>
#include <ostream>

#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "DataTree.hh"

//! Stores a model's equations and derivatives
class ModelTree : public DataTree
{
private:
  //! Stores declared equations
  vector<BinaryOpNode *> equations;

  typedef map<pair<int, int>, NodeID> first_derivatives_type;
  //! First order derivatives
  /*! First index is equation number, second is variable w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Variable indexes used are those of the variable_table, before sorting.
  */
  first_derivatives_type first_derivatives;

  typedef map<pair<int, pair<int, int> >, NodeID> second_derivatives_type;
  //! Second order derivatives
  /*! First index is equation number, second and third are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only second order derivatives where var1 >= var2 (for obvious symmetry reasons).
    Variable indexes used are those of the variable_table, before sorting.
  */
  second_derivatives_type second_derivatives;

  //! Temporary terms (those which will be noted Txxxx)
  temporary_terms_type temporary_terms;

  //! Computes derivatives of ModelTree
  void derive(int order);
  //! Computes temporary terms
  void computeTemporaryTerms(int order);

  //! Writes temporary terms
  void writeTemporaryTerms(ostream &output, bool is_dynamic) const;
  //! Writes local parameters
  void writeLocalParameters(ostream &output, bool is_dynamic) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, bool is_dynamic) const;
  //! Writes the static model equations and its derivatives
  /*! \todo handle hessian in C output */
  void writeStaticModel(ostream &StaticOutput) const;
  //! Writes the dynamic model equations and its derivatives
  void writeDynamicModel(ostream &DynamicOutput) const;
  //! Writes static model file (Matlab version)
  void writeStaticMFile(const string &static_basename) const;
  //! Writes static model file (C version)
  void writeStaticCFile(const string &static_basename) const;
  //! Writes dynamic model file (Matlab version)
  void writeDynamicMFile(const string &dynamic_basename) const;
  //! Writes dynamic model file (C version)
  void writeDynamicCFile(const string &dynamic_basename) const;

public:
  ModelTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Declare a node as an equation of the model
  void addEquation(NodeID eq);
  //! Do some checking
  void checkPass() const;
  //! Whether dynamic Jacobian (w.r. to endogenous) should be written
  bool computeJacobian;
  //! Whether dynamic Jacobian (w.r. to endogenous and exogenous) should be written
  bool computeJacobianExo;
  //! Whether dynamic Hessian (w.r. to endogenous and exogenous) should be written
  bool computeHessian;
  //! Whether static Hessian (w.r. to endogenous only) should be written
  bool computeStaticHessian;
  //! Execute computations (variable sorting + derivation)
  /*! You must set computeJacobian, computeJacobianExo and computeHessian to correct values before
    calling this function */
  void computingPass();
  //! Writes model initialization and lead/lag incidence matrix to output
  void writeOutput(ostream &output) const;
  //! Writes static model file
  void writeStaticFile(const string &basename) const;
  //! Writes dynamic model file
  void writeDynamicFile(const string &basename) const;
};

#endif
