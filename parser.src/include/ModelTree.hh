#ifndef _MODELTREE_HH
#define _MODELTREE_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/13/2003
  \par This file defines the ModelTree class.
*/
//------------------------------------------------------------------------------
#include <string>
#include <vector>
#include <list>
#include <stack>
#include <sstream>
#include <fstream>
#include <ostream>
//------------------------------------------------------------------------------
#include "SymbolTable.hh"
#include "OperatorTable.hh"
#include "NumericalConstants.hh"
#include "ModelTypes.hh"
#include "DataTree.hh"
//------------------------------------------------------------------------------
/*!
  \class  ModelTree
  \brief  Stores a model equations and derivatives.
*/
class ModelTree : public DataTree
{
private :
  /*! Stores ID of equations and their derivatives */
  std::vector<std::vector<DerivativeIndex> > mDerivativeIndex;
  /*!
    A token is writen as a temporary expression
    if its cost is greater than min_cost
  */
  int           min_cost;
  /*! left and right parentheses can be (,[ or ),] */
  char          lpar, rpar;
  //! Reference to numerical constants table
  const NumericalConstants &num_constants;

  /*! Computes argument derivative */
  inline NodeID     DeriveArgument(NodeID iArg, Type iType, int iVarID);
  /*! Gets output argument of terminal token */
  inline std::string    getArgument(NodeID id, Type type, EquationType  iEquationType);
  /*! Gets expression of part of model tree */
  inline std::string    getExpression(NodeID StartID, EquationType  iEquationType, int iEquationID = -1);
  inline int optimize(NodeID id);
  /*! Computes derivatives of ModelTree */
  void    derive(int iOrder);
  //! Writes the static model equations and its derivatives
  void writeStaticModel(std::ostream &StaticOutput);
  //! Writes the dynamic model equations and its derivatives
  void writeDynamicModel(std::ostream &DynamicOutput);
  //! Writes static model file (Matlab version)
  void writeStaticMFile(const std::string &static_basename);
  //! Writes static model file (C version)
  void writeStaticCFile(const std::string &static_basename);
  //! Writes dynamic model file (Matlab version)
  void writeDynamicMFile(const std::string &dynamic_basename);
  //! Writes dynamic model file (C version)
  void writeDynamicCFile(const std::string &dynamic_basename);

public:
  //! Constructor
  ModelTree(SymbolTable &symbol_table_arg, const NumericalConstants &num_constants);
  //! Number of equations contained in this model tree
  int eq_nbr;
  //! Do some checking
  void checkPass() const;
  //! Whether Jacobian (vs endogenous) should be written
  bool computeJacobian;
  //! Whether Jacobian (vs endogenous and exogenous) should be written
  bool computeJacobianExo;
  //! Whether Hessian (vs endogenous and exogenous) should be written
  bool computeHessian;
  //! Execute computations (variable sorting + derivation)
  /*! You must set computeJacobian, computeJacobianExo and computeHessian to correct values before
    calling this function */
  void computingPass();
  //! Writes model initialization and lead/lag incidence matrix to output
  void writeOutput(std::ostream &output) const;
  //! Writes static model file
  /*! \todo make this method const */
  void writeStaticFile(const std::string &basename);
  //! Writes dynamic model file
  /*! \todo make this method const */
  void writeDynamicFile(const std::string &basename);
};
//------------------------------------------------------------------------------
#endif
