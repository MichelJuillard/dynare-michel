#ifndef _DATATREE_HH
#define _DATATREE_HH

using namespace std;

#include <string>
#include <map>
#include <list>
#include <sstream>

#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "VariableTable.hh"
#include "ExprNode.hh"

class DataTree
{
  friend class ExprNode;
  friend class NumConstNode;
  friend class VariableNode;
  friend class UnaryOpNode;
  friend class BinaryOpNode;
protected:
  //! A reference to the symbol table
  SymbolTable &symbol_table;
  //! Reference to numerical constants table
  NumericalConstants &num_constants;
  
  typedef list<NodeID> node_list_type;
  //! The list of nodes
  node_list_type node_list;
  //! A counter for filling ExprNode's idx field
  int node_counter;

  //! Stores local variables value
  map<int, NodeID> local_variables_table;

  typedef map<int, NodeID> num_const_node_map_type;
  num_const_node_map_type num_const_node_map;
  //! Type (symbol_id, type, lag) used as key
  typedef map<pair<pair<int, Type>, int>, NodeID> variable_node_map_type;
  variable_node_map_type variable_node_map;
  typedef map<pair<NodeID, int>, NodeID> unary_op_node_map_type;
  unary_op_node_map_type unary_op_node_map;
  typedef map<pair<pair<NodeID, NodeID>, int>, NodeID> binary_op_node_map_type;
  binary_op_node_map_type binary_op_node_map;

  inline NodeID AddPossiblyNegativeConstant(double val);
  inline NodeID AddUnaryOp(UnaryOpcode op_code, NodeID arg);
  inline NodeID AddBinaryOp(NodeID arg1, BinaryOpcode op_code, NodeID arg2);
public:
  DataTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants_arg);
  virtual ~DataTree();
  //! The variable table
  VariableTable variable_table;
  NodeID Zero, One, MinusOne;

  //! Raised when a local parameter is declared twice
  class LocalParameterException
  {
  public:
    string name;
    LocalParameterException(const string &name_arg) : name(name_arg) {}
  };

  NodeID AddNumConstant(const string &value);
  NodeID AddVariable(const string &name, int lag = 0);
  //! Adds "arg1+arg2" to model tree
  NodeID AddPlus(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1-arg2" to model tree
  NodeID AddMinus(NodeID iArg1, NodeID iArg2);
  //! Adds "-arg" to model tree
  NodeID AddUMinus(NodeID iArg1);
  //! Adds "arg1*arg2" to model tree
  NodeID AddTimes(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1/arg2" to model tree
  NodeID AddDivide(NodeID iArg1, NodeID iArg2);
  //! Adds "arg1^arg2" to model tree
  NodeID AddPower(NodeID iArg1, NodeID iArg2);
  //! Adds "exp(arg)" to model tree
  NodeID AddExp(NodeID iArg1);
  //! Adds "log(arg)" to model tree
  NodeID AddLog(NodeID iArg1);
  //! Adds "log10(arg)" to model tree
  NodeID AddLog10(NodeID iArg1);
  //! Adds "cos(arg)" to model tree
  NodeID AddCos(NodeID iArg1);
  //! Adds "sin(arg)" to model tree
  NodeID AddSin(NodeID iArg1);
  //! Adds "tan(arg)" to model tree
  NodeID AddTan(NodeID iArg1);
  //! Adds "acos(arg)" to model tree
  NodeID AddACos(NodeID iArg1);
  //! Adds "asin(arg)" to model tree
  NodeID AddASin(NodeID iArg1);
  //! Adds "atan(arg)" to model tree
  NodeID AddATan(NodeID iArg1);
  //! Adds "cosh(arg)" to model tree
  NodeID AddCosH(NodeID iArg1);
  //! Adds "sinh(arg)" to model tree
  NodeID AddSinH(NodeID iArg1);
  //! Adds "tanh(arg)" to model tree
  NodeID AddTanH(NodeID iArg1);
  //! Adds "acosh(arg)" to model tree
  NodeID AddACosH(NodeID iArg1);
  //! Adds "asinh(arg)" to model tree
  NodeID AddASinH(NodeID iArg1);
  //! Adds "atanh(args)" to model tree
  NodeID AddATanH(NodeID iArg1);
  //! Adds "sqrt(arg)" to model tree
  NodeID AddSqRt(NodeID iArg1);
  //! Adds "arg1=arg2" to model tree
  NodeID AddEqual(NodeID iArg1, NodeID iArg2);
  void AddLocalParameter(const string &name, NodeID value) throw (LocalParameterException);
  //! Adds an unknown function node
  /*! \todo Use a map to share identical nodes */
  NodeID AddUnknownFunction(const string &function_name, const vector<NodeID> &arguments);
};

inline NodeID
DataTree::AddPossiblyNegativeConstant(double v)
{
  bool neg = false;
  if (v < 0)
    {
      v = -v;
      neg = true;
    }
  ostringstream ost;
  ost << v;
  
  NodeID cnode = AddNumConstant(ost.str());

  if (neg)
    return AddUMinus(cnode);
  else
    return cnode;
}

inline NodeID
DataTree::AddUnaryOp(UnaryOpcode op_code, NodeID arg)
{
  // If the node already exists in tree, share it
  unary_op_node_map_type::iterator it = unary_op_node_map.find(make_pair(arg, op_code));
  if (it != unary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  // Case where arg is a constant and op_code == oUminus (i.e. we're adding a negative constant) is skipped
  NumConstNode *carg = dynamic_cast<NumConstNode *>(arg);
  if (op_code != oUminus || carg == NULL)
    {
      try
        {
          double argval = arg->eval(eval_context_type());
          double val = UnaryOpNode::eval_opcode(op_code, argval);
          return AddPossiblyNegativeConstant(val);
        }
      catch(ExprNode::EvalException &e)
        {
        }
    }

  return new UnaryOpNode(*this, op_code, arg);
}

inline NodeID
DataTree::AddBinaryOp(NodeID arg1, BinaryOpcode op_code, NodeID arg2)
{
  binary_op_node_map_type::iterator it = binary_op_node_map.find(make_pair(make_pair(arg1, arg2), op_code));
  if (it != binary_op_node_map.end())
    return it->second;

  // Try to reduce to a constant
  try
    {
      double argval1 = arg1->eval(eval_context_type());
      double argval2 = arg2->eval(eval_context_type());
      double val = BinaryOpNode::eval_opcode(argval1, op_code, argval2);
      return AddPossiblyNegativeConstant(val);
    }
  catch(ExprNode::EvalException &e)
    {
    }

  return new BinaryOpNode(*this, arg1, op_code, arg2);
}

#endif
