#ifndef _EXPR_NODE_HH
#define _EXPR_NODE_HH

using namespace std;

#include <set>
#include <map>

#include "SymbolTableTypes.hh"

class DataTree;

typedef class ExprNode *NodeID;

typedef struct Model_Block;

struct ExprNodeLess;

//! Type for set of temporary terms
/*! They are ordered by index number thanks to ExprNodeLess */
typedef set<NodeID, ExprNodeLess> temporary_terms_type;

//! Possible types of output when writing ExprNode(s)
enum ExprNodeOutputType
  {
    oMatlabStaticModel,       //!< Matlab code, static model declarations
    oMatlabDynamicModel,      //!< Matlab code, dynamic model declarations
    oCStaticModel,            //!< C code, static model declarations
    oCDynamicModel,           //!< C code, dynamic model declarations
    oCDynamicModelSparseDLL,  //!< C code, dynamic model declarations in SparseDLL module
    oMatlabOutsideModel       //!< Matlab code, outside model block (for example in initval)
  };

/* Equal to 1 for Matlab langage, or to 0 for C language
   In Matlab, array indexes begin at 1, while they begin at 0 in C */ 
#define OFFSET(output_type) ((output_type == oMatlabStaticModel)      \
                             || (output_type == oMatlabDynamicModel)  \
                             || (output_type == oMatlabOutsideModel))

// Left parenthesis: '(' for Matlab, '[' for C
#define LPAR(output_type) (OFFSET(output_type) ? '(' : '[')

// Right parenthesis: ')' for Matlab, ']' for C
#define RPAR(output_type) (OFFSET(output_type) ? ')' : ']')

// Computing cost above which a node can be declared a temporary term
#define MIN_COST_MATLAB (40*90)
#define MIN_COST_C (40*4)
#define MIN_COST(is_matlab) (is_matlab ? MIN_COST_MATLAB : MIN_COST_C)

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
  virtual int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;

  //! set of endogenous variables in the current expression
  //! <symbolID, lag>
  set< pair<int,int> > present_endogenous;

public:
  ExprNode(DataTree &datatree_arg);
  virtual ~ExprNode();

  //! Returns derivative w.r. to variable varID
  /*! Uses a symbolic a priori to pre-detect null derivatives, and caches the result for other derivatives (to avoid computing it several times)
   For an equal node, returns the derivative of lhs minus rhs */
  NodeID getDerivative(int varID);

  //! Returns precedence of node
  /*! Equals 100 for constants, variables, unary ops, and temporary terms */
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;

  //! Fills temporary_terms set, using reference counts
  /*! A node will be marked as a temporary term if it is referenced at least two times (i.e. has at least two parents), and has a computing cost (multiplied by reference count) greater to datatree.min_cost */
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;

  //! Writes output of node, using a Txxx notation for nodes in temporary_terms
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const = 0;

  //! Collects the Endogenous in a expression
  virtual void collectEndogenous(NodeID &Id) = 0;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock) const;
  int present_endogenous_size() const;
  int present_endogenous_find(int var, int lag) const;
  virtual void Evaluate() const = 0;
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
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void collectEndogenous(NodeID &Id);
  virtual void Evaluate() const;
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
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms = temporary_terms_type()) const;
  virtual void collectEndogenous(NodeID &Id);
  virtual void Evaluate() const;
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

  int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const NodeID arg_arg);
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock) const;
  virtual void collectEndogenous(NodeID &Id);
  virtual void Evaluate() const;
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
  int cost(const temporary_terms_type &temporary_terms, bool is_matlab) const;
public:
  BinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
               BinaryOpcode op_code_arg, const NodeID arg2_arg);
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count, temporary_terms_type &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const;
  virtual void computeTemporaryTerms(map<NodeID, int> &reference_count,
                                     temporary_terms_type &temporary_terms,
                                     map<NodeID, int> &first_occurence,
                                     int Curr_block,
                                     Model_Block *ModelBlock) const;
  virtual void collectEndogenous(NodeID &Id);
  virtual void Evaluate() const;
};

typedef struct IM_compact
{
  int size, u_init, u_finish, nb_endo;
  int *u, *us, *Var, *Equ, *Var_Index, *Equ_Index, *Var_dyn_Index;
};

typedef struct Block
{
  int Size, Sized, Type, Simulation_Type, Max_Lead, Max_Lag, Nb_Lead_Lag_Endo;
  bool is_linear;
  int* Equation;
  int *Variable, *Variable_Sorted, *dVariable;
  int *variable_dyn_index, *variable_dyn_leadlag;
  temporary_terms_type *Temporary_terms;
  IM_compact *IM_lead_lag;
};

typedef struct Model_Block
{
  int Size, Periods;
  Block* Block_List;
  int *in_Block_Equ, *in_Block_Var, *in_Equ_of_Block, *in_Var_of_Block;
};

#endif
