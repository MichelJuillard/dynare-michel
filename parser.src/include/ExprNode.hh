#ifndef _EXPR_NODE_HH
#define _EXPR_NODE_HH

using namespace std;

#include <set>
#include <map>

#include "SymbolTableTypes.hh"

class DataTree;

typedef class ExprNode *NodeID;

struct ExprNodeLess;

//! Type for set of temporary terms
/*! They are ordered by index number thanks to ExprNodeLess */
typedef set<NodeID, ExprNodeLess> temporary_terms_type;

//! Base class for expression nodes
class ExprNode
{
  friend class DataTree;
  friend class ModelTree;
  friend class ExprNodeLess;
  friend class NumConstNode;
  friend class VariableNode;
  friend class UnaryOpNode;
  friend class BinaryOpNode;

private:
  //! Computes derivative w.r. to variable varID (but doesn't store it in derivatives map)
  /*! You shoud use getDerivative() to get the benefit of symbolic a priori and of caching */
  virtual NodeID computeDerivative(int varID) = 0;

protected:
  //! Reference to the enclosing DataTree
  DataTree &datatree;

  //! Index number
  int idx;

  //! Set of variable IDs with respect to which the derivative is potentially non-null
  set<int> non_null_derivatives;

  //! Used for caching of first order derivatives (when non-null)
  map<int, NodeID> derivatives;

  //! Cost of computing current node
  /*! Nodes included in temporary_terms are considered having a null cost */
  virtual int cost(temporary_terms_type &temporary_terms) const;

public:
  ExprNode(DataTree &datatree_arg);
  virtual ~ExprNode();

  //! Returns derivative w.r. to variable varID
  /*! Uses a symbolic a priori to pre-detect null derivatives, and caches the result for other derivatives (to avoid computing it several times)
   For an equal node, returns the derivative of lhs minus rhs */
  NodeID getDerivative(int varID);

  //! Returns precedence of node
  /*! Equals 100 for constants, variables, unary ops, and temporary terms */
  virtual int precedence(const temporary_terms_type &temporary_terms) const;

  //! Fills temporary_terms set, using reference counts
  /*! A node will be marked as a temporary term if it is referenced at least two times (i.e. has at least two parents), and has a computing cost (multiplied by reference count) greater to datatree.min_cost */
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms) const;

  //! Writes output of node, using a Txxx notation for nodes in temporary_terms
  virtual void writeOutput(ostream &output, bool is_dynamic, const temporary_terms_type &temporary_terms) const = 0;
};

//! Object used to compare two nodes (using their indexes)
struct ExprNodeLess
{
  bool operator()(NodeID arg1, NodeID arg2) const
  {
    return arg1->idx < arg2->idx;
  }
};

//! Numerical constant node
class NumConstNode : public ExprNode
{
private:
  //! Id from numerical constants table
  const int id;
  virtual NodeID computeDerivative(int varID);
public:
  NumConstNode(DataTree &datatree_arg, int id_arg);
  virtual void writeOutput(ostream &output, bool is_dynamic, const temporary_terms_type &temporary_terms) const;
};

//! Symbol or variable node
class VariableNode : public ExprNode
{
private:
  //! Id of the symbol/variable
  /*! For an endogenous, exogenous or recursive variable, the id is taken from the variable table (before sorting). For other types of symbols, it's the id from the symbol table. */
  const int id;
  const Type type;
  virtual NodeID computeDerivative(int varID);
public:
  VariableNode(DataTree &datatree_arg, int id_arg, Type type_arg);
  virtual void writeOutput(ostream &output, bool is_dynamic, const temporary_terms_type &temporary_terms) const;
};

enum UnaryOpcode
  {
    oUminus,
    oExp,
    oLog,
    oLog10,
    oCos,
    oSin,
    oTan,
    oAcos,
    oAsin,
    oAtan,
    oCosh,
    oSinh,
    oTanh,
    oAcosh,
    oAsinh,
    oAtanh,
    oSqrt
  };

//! Unary operator node
class UnaryOpNode : public ExprNode
{
  friend class DataTree;
private:
  const NodeID arg;
  const UnaryOpcode op_code;
  virtual NodeID computeDerivative(int varID);

  int cost(temporary_terms_type &temporary_terms) const;
public:
  UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const NodeID arg_arg);
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms) const;
  virtual void writeOutput(ostream &output, bool is_dynamic, const temporary_terms_type &temporary_terms) const;
};

enum BinaryOpcode
  {
    oPlus,
    oMinus,
    oTimes,
    oDivide,
    oPower,
    oEqual
  };

//! Binary operator node
class BinaryOpNode : public ExprNode
{
  friend class ModelTree;
private:
  const NodeID arg1, arg2;
  const BinaryOpcode op_code;
  virtual NodeID computeDerivative(int varID);
  int cost(temporary_terms_type &temporary_terms) const;
public:
  BinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
               BinaryOpcode op_code_arg, const NodeID arg2_arg);
  virtual int precedence(const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms) const;
  virtual void writeOutput(ostream &output, bool is_dynamic, const temporary_terms_type &temporary_terms) const;
};

#endif
