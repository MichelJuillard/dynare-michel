#include <iostream>
#include <iterator>
#include <algorithm>

#include "ExprNode.hh"
#include "DataTree.hh"

ExprNode::ExprNode(DataTree &datatree_arg) : datatree(datatree_arg)
{
  // Add myself to datatree
  datatree.node_list.push_back(this);

  // Set my index and increment counter
  idx = datatree.node_counter++;
}

ExprNode::~ExprNode()
{
}

NodeID
ExprNode::getDerivative(int varID)
{
  // Return zero if derivative is necessarily null (using symbolic a priori)
  set<int>::const_iterator it = non_null_derivatives.find(varID);
  if (it == non_null_derivatives.end())
    return datatree.Zero;

  // If derivative is stored in cache, use the cached value, otherwise compute it (and cache it)
  map<int, NodeID>::const_iterator it2 = derivatives.find(varID);
  if (it2 != derivatives.end())
    return it2->second;
  else
    {
      NodeID d = computeDerivative(varID);
      derivatives[varID] = d;
      return d;
    }
}

int
ExprNode::precedence(const temporary_terms_type &temporary_terms) const
{
  // For a constant, a variable, or a unary op, the precedence is maximal
  return 100;
}

int
ExprNode::cost(temporary_terms_type &temporary_terms) const
{
  // For a terminal node, the cost is null
  return 0;
}

void
ExprNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                temporary_terms_type &temporary_terms) const
{
  // Nothing to do for a terminal node
}

NumConstNode::NumConstNode(DataTree &datatree_arg, int id_arg) :
  ExprNode(datatree_arg),
  id(id_arg)
{
  // Add myself to the num const map
  datatree.num_const_node_map[id] = this;

  // All derivatives are null, so non_null_derivatives is left empty
}

NodeID
NumConstNode::computeDerivative(int varID)
{
  return datatree.Zero;
}

void
NumConstNode::writeOutput(ostream &output, bool is_dynamic,
                          const temporary_terms_type &temporary_terms) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
  if (it != temporary_terms.end())
    output << "T" << idx;
  else
    output << datatree.num_constants.get(id);
}

VariableNode::VariableNode(DataTree &datatree_arg, int id_arg, Type type_arg) :
  ExprNode(datatree_arg),
  id(id_arg),
  type(type_arg)
{
  // Add myself to the variable map
  datatree.variable_node_map[make_pair(id, type)] = this;

  // Fill in non_null_derivatives
  switch(type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eRecursiveVariable:
      // For a variable, the only non-null derivative is with respect to itself
      non_null_derivatives.insert(id);
      break;
    case eParameter:
      // All derivatives are null, do nothing
      break;
    case eLocalParameter:
      // Non null derivatives are those of the value of the local parameter
      non_null_derivatives = datatree.local_parameters_table[id]->non_null_derivatives;
      break;
    case eNumericalConstant:
    case eUNDEF:
    case eTempResult:
      // Impossible cases
      cerr << "Incorrect symbol type used in VariableNode" << endl;
      exit(-1);
    }
}

NodeID
VariableNode::computeDerivative(int varID)
{
  switch(type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eRecursiveVariable:
      if (varID == id)
        return datatree.One;
      else
        return datatree.Zero;
    case eParameter:
      return datatree.Zero;
    case eLocalParameter:
      return datatree.local_parameters_table[id]->getDerivative(varID);
    case eNumericalConstant:
    case eUNDEF:
    case eTempResult:
      // Impossible cases
      cerr << "Incorrect symbol type used in VariableNode" << endl;
      exit(-1);
    }
  cerr << "Impossible case!" << endl;
  exit(-1);
}

void
VariableNode::writeOutput(ostream &output, bool is_dynamic,
                          const temporary_terms_type &temporary_terms) const
{
  // If node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  char &lpar = datatree.lpar;
  char &rpar = datatree.rpar;

  int idx;
  switch(type)
    {
    case eParameter:
      output << "params" << lpar << id + datatree.offset << rpar;
      break;
    case eLocalParameter:
      output << datatree.symbol_table.getNameByID(eLocalParameter, id);
      break;
    case eEndogenous:
      if (is_dynamic)
        {
          idx = datatree.variable_table.getPrintIndex(id) + datatree.offset;
          output <<  "y" << lpar << idx << rpar;
        }
      else
        {
          idx = datatree.variable_table.getSymbolID(id) + datatree.offset;
          output <<  "y" << lpar << idx << rpar;
        }
      break;
    case eExogenous:
      idx = datatree.variable_table.getSymbolID(id) + datatree.offset;
      if (is_dynamic)
          {
            int lag = datatree.variable_table.getLag(id);
            if (datatree.offset == 1)
              if (lag != 0)
                output <<  "x" << lpar << "it_ + " << lag << ", " << idx << rpar;
              else
                output <<  "x" << lpar << "it_, " << idx << rpar;
            else
              if (lag != 0)
                output <<  "x" << lpar << "it_+" << lag << "+" << idx << "*nb_row_x" << rpar;
              else
                output <<  "x" << lpar << "it_+" << idx << "*nb_row_x" << rpar;
          }
      else
        output << "x" << lpar << idx << rpar;
      break;
    case eExogenousDet:
      idx = datatree.variable_table.getSymbolID(id) + datatree.symbol_table.exo_nbr
        + datatree.offset;
      if (is_dynamic)
          {
            int lag = datatree.variable_table.getLag(id);
            if (datatree.offset == 1)
              if (lag != 0)
                output <<  "x" << lpar << "it_ + " << lag << ", " << idx << rpar;
              else
                output <<  "x" << lpar << "it_, " << idx << rpar;
            else
              if (lag != 0)
                output <<  "x" << lpar << "it_ + " << lag << "+" << idx <<  "*nb_row_xd" << rpar;
              else
                output <<  "x" << lpar << "it_+" << idx << "*nb_row_xd" <<  rpar;
          }
      else
        output <<  "x" << lpar << idx << rpar;
      break;
    case eRecursiveVariable:
      cerr << "Recursive variable not implemented" << endl;
      exit(-1);
    case eTempResult:
    case eNumericalConstant:
    case eUNDEF:
      // Impossible cases
      cerr << "Incorrect symbol type used in VariableNode" << endl;
      exit(-1);
    }
}

UnaryOpNode::UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const NodeID arg_arg) :
  ExprNode(datatree_arg),
  arg(arg_arg),
  op_code(op_code_arg)
{
  // Add myself to the unary op map
  datatree.unary_op_node_map[make_pair(arg, op_code)] = this;

  // Non-null derivatives are those of the argument
  non_null_derivatives = arg->non_null_derivatives;
}

NodeID
UnaryOpNode::computeDerivative(int varID)
{
  NodeID darg = arg->getDerivative(varID);

  NodeID t11, t12, t13;

  switch(op_code)
    {
    case oUminus:
      return datatree.AddUMinus(darg);
    case oExp:
      return datatree.AddTimes(darg, this);
    case oLog:
      return datatree.AddDivide(darg, arg);
    case oLog10:
      t11 = datatree.AddExp(datatree.One);
      t12 = datatree.AddLog10(t11);
      t13 = datatree.AddDivide(darg, arg);
      return datatree.AddTimes(t12, t13);
    case oCos:
      t11 = datatree.AddSin(arg);
      t12 = datatree.AddUMinus(t11);
      return datatree.AddTimes(darg, t12);
    case oSin:
      t11 = datatree.AddCos(arg);
      return datatree.AddTimes(darg, t11);
    case oTan:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddPlus(t11, datatree.One);
      return datatree.AddTimes(darg, t12);
    case oAcos:
      t11 = datatree.AddSin(this);
      t12 = datatree.AddDivide(darg, t11);
      return datatree.AddUMinus(t12);
    case oAsin:
      t11 = datatree.AddCos(this);
      return datatree.AddDivide(darg, t11);
    case oAtan:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddPlus(datatree.One, t11);
      return datatree.AddDivide(darg, t12);
    case oCosh:
      t11 = datatree.AddSinH(arg);
      return datatree.AddTimes(darg, t11);
    case oSinh:
      t11 = datatree.AddCosH(arg);
      return datatree.AddTimes(darg, t11);
    case oTanh:
      t11 = datatree.AddTimes(this, this);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oAcosh:
      t11 = datatree.AddSinH(this);
      return datatree.AddDivide(darg, t11);
    case oAsinh:
      t11 = datatree.AddCosH(this);
      return datatree.AddDivide(darg, t11);
    case oAtanh:
      t11 = datatree.AddTimes(arg, arg);
      t12 = datatree.AddMinus(datatree.One, t11);
      return datatree.AddTimes(darg, t12);
    case oSqrt:
      t11 = datatree.AddPlus(this, this);
      return datatree.AddDivide(darg, t11);
    }
  cerr << "Impossible case!" << endl;
  exit(-1);
}

int
UnaryOpNode::cost(temporary_terms_type &temporary_terms) const
{
  // For a temporary term, the cost is null
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    return 0;

  int cost = arg->cost(temporary_terms);

  if (datatree.offset)
    // Cost for Matlab files
    switch(op_code)
      {
      case oUminus:
        return cost + 70;
      case oExp:
        return cost + 160;
      case oLog:
        return cost + 300;
      case oLog10:
        return cost + 16000;
      case oCos:
      case oSin:
      case oCosh:
        return cost + 210;
      case oTan:
        return cost + 230;
      case oAcos:
        return cost + 300;
      case oAsin:
        return cost + 310;
      case oAtan:
        return cost + 140;
      case oSinh:
        return cost + 240;
      case oTanh:
        return cost + 190;
      case oAcosh:
        return cost + 770;
      case oAsinh:
        return cost + 460;
      case oAtanh:
        return cost + 350;
      case oSqrt:
        return cost + 570;
      }
  else
    // Cost for C files
    switch(op_code)
      {
      case oUminus:
        return cost + 3;
      case oExp:
      case oAcosh:
        return cost + 210;
      case oLog:
        return cost + 137;
      case oLog10:
        return cost + 139;
      case oCos:
      case oSin:
        return cost + 160;
      case oTan:
        return cost + 170;
      case oAcos:
      case oAtan:
        return cost + 190;
      case oAsin:
        return cost + 180;
      case oCosh:
      case oSinh:
      case oTanh:
        return cost + 240;
      case oAsinh:
        return cost + 220;
      case oAtanh:
        return cost + 150;
      case oSqrt:
        return cost + 90;
      }
  cerr << "Impossible case!" << endl;
  exit(-1);
}       

void
UnaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                   temporary_terms_type &temporary_terms) const
{
  NodeID this2 = const_cast<UnaryOpNode *>(this);

  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      arg->computeTemporaryTerms(reference_count, temporary_terms);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms) > datatree.min_cost)
        temporary_terms.insert(this2);
    }
}

void
UnaryOpNode::writeOutput(ostream &output, bool is_dynamic,
                         const temporary_terms_type &temporary_terms) const
{
  // If node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  // Always put parenthesis around uminus nodes
  if (op_code == oUminus)
    output << "(";

  switch(op_code)
    {
    case oUminus:
      output << "-";
      break;
    case oExp:
      output << "exp";
      break;
    case oLog:
      output << "log";
      break;
    case oLog10:
      output << "log10";
      break;
    case oCos:
      output << "cos";
      break;
    case oSin:
      output << "sin";
      break;
    case oTan:
      output << "tan";
      break;
    case oAcos:
      output << "acos";
      break;
    case oAsin:
      output << "asin";
      break;
    case oAtan:
      output << "atan";
      break;
    case oCosh:
      output << "cosh";
      break;
    case oSinh:
      output << "sinh";
      break;
    case oTanh:
      output << "tanh";
      break;
    case oAcosh:
      output << "acosh";
      break;
    case oAsinh:
      output << "asinh";
      break;
    case oAtanh:
      output << "atanh";
      break;
    case oSqrt:
      output << "sqrt";
      break;
    }

  bool close_parenthesis = false;

  /* Enclose argument with parentheses if:
     - current opcode is not uminus, or
     - current opcode is uminus and argument has lowest precedence
  */
  if (op_code != oUminus
      || (op_code == oUminus && arg->precedence(temporary_terms) < precedence(temporary_terms)))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write argument
  arg->writeOutput(output, is_dynamic, temporary_terms);

  if (close_parenthesis)
    output << ")";

  // Close parenthesis for uminus
  if (op_code == oUminus)
    output << ")";
}

BinaryOpNode::BinaryOpNode(DataTree &datatree_arg, const NodeID arg1_arg,
                           BinaryOpcode op_code_arg, const NodeID arg2_arg) :
  ExprNode(datatree_arg),
  arg1(arg1_arg),
  arg2(arg2_arg),
  op_code(op_code_arg)
{
  datatree.binary_op_node_map[make_pair(make_pair(arg1, arg2), op_code)] = this;

  // Non-null derivatives are the union of those of the arguments
  // Compute set union of arg1->non_null_derivatives and arg2->non_null_derivatives
  set_union(arg1->non_null_derivatives.begin(),
            arg1->non_null_derivatives.end(),
            arg2->non_null_derivatives.begin(),
            arg2->non_null_derivatives.end(),
            inserter(non_null_derivatives, non_null_derivatives.begin()));
}

NodeID
BinaryOpNode::computeDerivative(int varID)
{
  NodeID darg1 = arg1->getDerivative(varID);
  NodeID darg2 = arg2->getDerivative(varID);

  NodeID t11, t12, t13, t14, t15;

  switch(op_code)
    {
    case oPlus:
      return datatree.AddPlus(darg1, darg2);
    case oMinus:
      return datatree.AddMinus(darg1, darg2);
    case oTimes:
      t11 = datatree.AddTimes(darg1, arg2);
      t12 = datatree.AddTimes(darg2, arg1);
      return datatree.AddPlus(t11, t12);
    case oDivide:
      t11 = datatree.AddTimes(darg1, arg2);
      t12 = datatree.AddTimes(darg2, arg1);
      t13 = datatree.AddMinus(t11, t12);
      t14 = datatree.AddTimes(arg2, arg2);
      return datatree.AddDivide(t13, t14);
    case oPower:
      if (darg2 == datatree.Zero)
        {
          if (darg1 == datatree.Zero)
            return datatree.Zero;
          else
            {
              t11 = datatree.AddMinus(arg2, datatree.One);
              t12 = datatree.AddPower(arg1, t11);
              t13 = datatree.AddTimes(arg2, t12);
              return datatree.AddTimes(darg1, t13);
            }
        }
      else
        {
          t11 = datatree.AddLog(arg1);
          t12 = datatree.AddTimes(darg2, t11);
          t13 = datatree.AddTimes(darg1, arg2);
          t14 = datatree.AddDivide(t13, arg1);
          t15 = datatree.AddPlus(t12, t14);
          return datatree.AddTimes(t15, this);
        }
    case oEqual:
      return datatree.AddMinus(darg1, darg2);
    }
  cerr << "Impossible case!" << endl;
  exit(-1);
}

int
BinaryOpNode::precedence(const temporary_terms_type &temporary_terms) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  // A temporary term behaves as a variable
  if (it != temporary_terms.end())
    return 100;

  switch(op_code)
    {
    case oEqual:
    case oPlus:
    case oMinus:
      return 0;
    case oTimes:
    case oDivide:
      return 1;
    case oPower:
      if (datatree.offset)
        // In C, power operator is of the form pow(a, b)
        return 100;
      else
        return 3;
    }
  cerr << "Impossible case!" << endl;
  exit(-1);
}

int
BinaryOpNode::cost(temporary_terms_type &temporary_terms) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  // For a temporary term, the cost is null
  if (it != temporary_terms.end())
    return 0;

  int cost = arg1->cost(temporary_terms);
  cost += arg2->cost(temporary_terms);

  if (datatree.offset)
    // Cost for Matlab files
    switch(op_code)
      {
      case oPlus:
      case oMinus:
      case oTimes:
        return cost + 90;
      case oDivide:
        return cost + 990;
      case oPower:
        return cost + 1160;
      case oEqual:
        return cost;
      }
  else
    // Cost for C files
    switch(op_code)
      {
      case oPlus:
      case oMinus:
      case oTimes:
        return cost + 4;
      case oDivide:
        return cost + 15;
      case oPower:
        return cost + 520;
      case oEqual:
        return cost;
      }
  cerr << "Impossible case!" << endl;
  exit(-1);
}

void
BinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms) const
{
  NodeID this2 = const_cast<BinaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = 1;
      arg1->computeTemporaryTerms(reference_count, temporary_terms);
      arg2->computeTemporaryTerms(reference_count, temporary_terms);
    }
  else
    {
      // If the node has already been encountered, increment its ref count
      //  and declare it as a temporary term if it is too costly
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms) > datatree.min_cost)
        temporary_terms.insert(this2);
    }
}

void
BinaryOpNode::writeOutput(ostream &output, bool is_dynamic,
                          const temporary_terms_type &temporary_terms) const
{
  // If current node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      output << "T" << idx;
      return;
    }

  // Treat special case of power operator in C
  if (op_code == oPower && datatree.offset == 0)
    {
      output << "pow(";
      arg1->writeOutput(output, is_dynamic, temporary_terms);
      output << ",";
      arg2->writeOutput(output, is_dynamic, temporary_terms);
      output << ")";
      return;
    }

  int prec = precedence(temporary_terms);

  bool close_parenthesis = false;

  // If left argument has a lower precedence, or if current and left argument are both power operators, add parenthesis around left argument
  BinaryOpNode *barg1 = dynamic_cast<BinaryOpNode *>(arg1);
  if (arg1->precedence(temporary_terms) < prec
      || (op_code == oPower && barg1 != NULL && barg1->op_code == oPower))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write left argument
  arg1->writeOutput(output, is_dynamic, temporary_terms);

  if (close_parenthesis)
    output << ")";

  // Write current operator symbol
  switch(op_code)
    {
    case oPlus:
      output << "+";
      break;
    case oMinus:
      output << "-";
      break;
    case oTimes:
      output << "*";
      break;
    case oDivide:
      output << "/";
      break;
    case oPower:
      output << "^";
      break;
    case oEqual:
      output << "=";
      break;
    }

  close_parenthesis = false;

  /* Add parenthesis around right argument if:
     - its precedence is lower than those of the current node
     - it is a power operator and current operator is also a power operator
     - it is a minus operator with same precedence than current operator
     - it is a divide operator with same precedence than current operator */
  BinaryOpNode *barg2 = dynamic_cast<BinaryOpNode *>(arg2);
  int arg2_prec = arg2->precedence(temporary_terms);
  if (arg2_prec < prec
      || (op_code == oPower && barg2 != NULL && barg2->op_code == oPower)
      || (op_code == oMinus && arg2_prec == prec)
      || (op_code == oDivide && arg2_prec == prec))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write right argument
  arg2->writeOutput(output, is_dynamic, temporary_terms);

  if (close_parenthesis)
    output << ")";
}
