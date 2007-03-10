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
#include "BlockTriangular.hh"

#define LCC_COMPILE 0
#define GCC_COMPILE 1

//! The three in which ModelTree can work
enum ModelTreeMode
  {
    eStandardMode, //!< Standard mode (static and dynamic files in Matlab)
    eDLLMode,      //!< DLL mode (static and dynamic files in C)
    eSparseDLLMode //!< Sparse DLL mode (static file in Matlab, dynamic file in C with block decomposition plus a binary file)
  };

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

  typedef map<pair<int, pair<int, pair<int, int> > >, NodeID> third_derivatives_type;
  //! Third order derivatives
  /*! First index is equation number, second, third and fourth are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only third order derivatives where var1 >= var2 >= var3 (for obvious symmetry reasons).
    Variable indexes used are those of the variable_table, before sorting.
  */
  third_derivatives_type third_derivatives;

  //! Temporary terms (those which will be noted Txxxx)
  temporary_terms_type temporary_terms;

  //! Computes derivatives of ModelTree
  void derive(int order);
  //! Write derivative of an equation w.r. to a variable
  void writeDerivative(ostream &output, int eq, int var, int lag, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  //! Computes temporary terms
  void computeTemporaryTerms(int order);
  void computeTemporaryTermsOrdered(int order, Model_Block *ModelBlock);
  //! Writes temporary terms
  void writeTemporaryTerms(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes local parameters
  /*! No temporary term is used in the output, so that local parameters declarations can be safely put before temporary terms declaration in the output files */
  void writeLocalParameters(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type) const;
  //! Writes the static model equations and its derivatives
  /*! \todo handle hessian in C output */
  void writeStaticModel(ostream &StaticOutput) const;
  //! Writes the dynamic model equations and its derivatives
  /*! \todo add third derivatives handling in C output */
  void writeDynamicModel(ostream &DynamicOutput) const;
  void writeModelEquationsOrdered(ostream &output, Model_Block *ModelBlock) const;
  //! Writes static model file (Matlab version)
  void writeStaticMFile(const string &static_basename) const;
  //! Writes static model file (C version)
  void writeStaticCFile(const string &static_basename) const;
  //! Writes dynamic model file (Matlab version)
  void writeDynamicMFile(const string &dynamic_basename) const;
  //! Writes dynamic model file (C version)
  /*! \todo add third derivatives handling */
  void writeDynamicCFile(const string &dynamic_basename) const;
  //! Writes dynamic model header file when SparseDLL option is on
  void writeSparseDLLDynamicHFile(const string &dynamic_basename) const;
  //! Writes dynamic model file when SparseDLL option is on
  void writeSparseDLLDynamicCFileAndBinFile(const string &dynamic_basename, const string &bin_basename) const;
  void evaluateJacobian(const eval_context_type &eval_context);
  void BlockLinear(Model_Block *ModelBlock);
  string reform(string name) const;

public:
  ModelTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants);
  //! Mode in which the ModelTree is supposed to work (Matlab, DLL or SparseDLL)
  ModelTreeMode mode;
  //! Type of compiler used in matlab for SPARSE_DLL option: 0 = LCC or 1 = GCC
  int compiler;
  //! Absolute value under which a number is considered to be zero
  double cutoff;

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
  //! Whether dynamic third order derivatives (w.r. to endogenous and exogenous) should be written
  bool computeThirdDerivatives;
  //! Execute computations (variable sorting + derivation)
  /*! You must set computeJacobian, computeJacobianExo, computeHessian, computeStaticHessian and computeThirdDerivatives to correct values before calling this function */
  void computingPass(const eval_context_type &eval_context);
  //! Writes model initialization and lead/lag incidence matrix to output
  void writeOutput(ostream &output) const;
  //! Writes static model file
  void writeStaticFile(const string &basename) const;
  //! Writes dynamic model file
  void writeDynamicFile(const string &basename) const;
  //! Complete set to block decompose the model
  BlockTriangular block_triangular;
  int equation_number() const;
};

#endif
