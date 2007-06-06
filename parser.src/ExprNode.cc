#include <iostream>
#include <iterator>
#include <algorithm>

#include <math.h>

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
ExprNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
{
  // For a constant, a variable, or a unary op, the precedence is maximal
  return 100;
}

int
ExprNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
{
  // For a terminal node, the cost is null
  return 0;
}

int
ExprNode::present_endogenous_size() const
{
  return(present_endogenous.size());
}

int
ExprNode::present_endogenous_find(int var, int lag) const
{
  return(present_endogenous.find(make_pair(var,lag))!=present_endogenous.end());
}

void
ExprNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                temporary_terms_type &temporary_terms,
                                bool is_matlab) const
{
  // Nothing to do for a terminal node
}

void
ExprNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                temporary_terms_type &temporary_terms,
                                map<NodeID, int> &first_occurence,
                                int Curr_block,
                                Model_Block *ModelBlock) const
{
  // Nothing to do for a terminal node
}

void
ExprNode::writeOutput(ostream &output)
{
  writeOutput(output, oMatlabOutsideModel, temporary_terms_type());
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
NumConstNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<NumConstNode *>(this));
  if (it != temporary_terms.end())
    if (output_type != oCDynamicModelSparseDLL)
      output << "T" << idx;
    else
      output << "T" << idx << "[it_]";
  else
    output << datatree.num_constants.get(id);
}

double
NumConstNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  return(datatree.num_constants.getDouble(id));
}

void
NumConstNode::collectEndogenous(NodeID &Id)
{
}

VariableNode::VariableNode(DataTree &datatree_arg, int symb_id_arg, Type type_arg, int lag_arg) :
  ExprNode(datatree_arg),
  symb_id(symb_id_arg),
  type(type_arg),
  lag(lag_arg)
{
  // Add myself to the variable map
  datatree.variable_node_map[make_pair(make_pair(symb_id, type), lag)] = this;

  // Add myself to the variable table if necessary and initialize var_id
  if (type == eEndogenous
      || type == eExogenousDet
      || type == eExogenous
      || type == eRecursiveVariable)
    var_id = datatree.variable_table.AddVariable(datatree.symbol_table.getNameByID(type, symb_id), lag);
  else
    var_id = -1;

  // Fill in non_null_derivatives
  switch(type)
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
    case eRecursiveVariable:
      // For a variable, the only non-null derivative is with respect to itself
      non_null_derivatives.insert(var_id);
      break;
    case eParameter:
      // All derivatives are null, do nothing
      break;
    case eModelLocalVariable:
      // Non null derivatives are those of the value of the local parameter
      non_null_derivatives = datatree.local_variables_table[symb_id]->non_null_derivatives;
      break;
    case eModFileLocalVariable:
      // Such a variable is never derived
      break;
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
      if (varID == var_id)
        return datatree.One;
      else
        return datatree.Zero;
    case eParameter:
      return datatree.Zero;
    case eModelLocalVariable:
      return datatree.local_variables_table[symb_id]->getDerivative(varID);
    case eModFileLocalVariable:
      cerr << "ModFileLocalVariable is not derivable" << endl;
      exit(-1);
    }
  cerr << "Impossible case!" << endl;
  exit(-1);
}

void
VariableNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
{
  // If node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<VariableNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type != oCDynamicModelSparseDLL)
        output << "T" << idx;
      else
        output << "T" << idx << "[it_]";
      return;
    }

  int i;
  switch(type)
    {
    case eParameter:
      if (output_type == oMatlabOutsideModel)
        output << "M_.params" << "(" << symb_id + 1 << ")";
      else
        output << "params" << LPAR(output_type) << symb_id + OFFSET(output_type) << RPAR(output_type);
      break;

    case eModelLocalVariable:
    case eModFileLocalVariable:
      output << datatree.symbol_table.getNameByID(type, symb_id);
      break;

    case eEndogenous:
      switch(output_type)
        {
        case oMatlabDynamicModel:
        case oCDynamicModel:
          i = datatree.variable_table.getPrintIndex(var_id) + OFFSET(output_type);
          output <<  "y" << LPAR(output_type) << i << RPAR(output_type);
          break;
        case oMatlabStaticModel:
        case oCStaticModel:
          i = symb_id + OFFSET(output_type);
          output <<  "y" << LPAR(output_type) << i << RPAR(output_type);
          break;
        case oCDynamicModelSparseDLL:
          if (lag > 0)
            output << "y" << LPAR(output_type) << "(it_+" << lag << ")*y_size+" << symb_id << RPAR(output_type);
          else if (lag < 0)
            output << "y" << LPAR(output_type) << "(it_" << lag << ")*y_size+" << symb_id << RPAR(output_type);
          else
            output << "y" << LPAR(output_type) << "Per_y_+" << symb_id << RPAR(output_type);
          break;
        case oMatlabOutsideModel:
          output << "oo_.steady_state" << "(" << symb_id + 1 << ")";
          break;
        }
      break;

    case eExogenous:
      i = symb_id + OFFSET(output_type);
      switch(output_type)
        {
        case oMatlabDynamicModel:
          if (lag > 0)
            output <<  "x(it_+" << lag << ", " << i << ")";
          else if (lag < 0)
            output <<  "x(it_" << lag << ", " << i << ")";
          else
            output <<  "x(it_, " << i << ")";
          break;
        case oCDynamicModel:
        case oCDynamicModelSparseDLL:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_x]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_x]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_x]";
          break;
        case oMatlabStaticModel:
        case oCStaticModel:
          output << "x" << LPAR(output_type) << i << RPAR(output_type);
          break;
        case oMatlabOutsideModel:
          if (lag != 0)
            {
              cerr << "VariableNode::writeOutput: lag != 0 for exogenous variable outside model scope!" << endl;
              exit(-1);
            }
          output <<  "oo_.exo_steady_state" << "(" << i << ")";
          break;
        }
      break;

    case eExogenousDet:
      i = symb_id + datatree.symbol_table.exo_nbr + OFFSET(output_type);
      switch(output_type)
        {
        case oMatlabDynamicModel:
          if (lag > 0)
            output <<  "x(it_+" << lag << ", " << i << ")";
          else if (lag < 0)
            output <<  "x(it_" << lag << ", " << i << ")";
          else
            output <<  "x(it_, " << i << ")";
          break;
        case oCDynamicModel:
        case oCDynamicModelSparseDLL:
          if (lag == 0)
            output <<  "x[it_+" << i << "*nb_row_xd]";
          else if (lag > 0)
            output <<  "x[it_+" << lag << "+" << i << "*nb_row_xd]";
          else
            output <<  "x[it_" << lag << "+" << i << "*nb_row_xd]";
          break;
        case oMatlabStaticModel:
        case oCStaticModel:
          output << "x" << LPAR(output_type) << i << RPAR(output_type);
          break;
        case oMatlabOutsideModel:
          if (lag != 0)
            {
              cerr << "VariableNode::writeOutput: lag != 0 for exogenous determistic variable outside model scope!" << endl;
              exit(-1);
            }
          output <<  "oo_.exo_det_steady_state" << "(" << symb_id + 1 << ")";
          break;
        }
      break;

    case eRecursiveVariable:
      cerr << "Recursive variable not implemented" << endl;
      exit(-1);
    }
}

double
VariableNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  /*if (lag != 0)
    throw EvalException();*/
  eval_context_type::const_iterator it = eval_context.find(make_pair(symb_id, type));
  if (it == eval_context.end())
    {
      if (eval_context.size()>0)
        {
          cout << "Error: the variable or parameter (" << datatree.symbol_table.getNameByID( type, symb_id) << ") has not been initialized (in derivatives evaluation)\n";
          cout.flush();
        }
      throw EvalException();
    }
  return it->second;
}

void
VariableNode::collectEndogenous(NodeID &Id)
{
  if (type == eEndogenous)
    Id->present_endogenous.insert(make_pair(symb_id, lag));
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
UnaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
{
  // For a temporary term, the cost is null
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    return 0;

  int cost = arg->cost(temporary_terms, is_matlab);

  if (is_matlab)
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
                                   temporary_terms_type &temporary_terms,
                                   bool is_matlab) const
{
  NodeID this2 = const_cast<UnaryOpNode *>(this);

  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      arg->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab))
        temporary_terms.insert(this2);
    }
}

void
UnaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                   temporary_terms_type &temporary_terms,
                                   map<NodeID, int> &first_occurence,
                                   int Curr_block,
                                   Model_Block *ModelBlock) const
{
  NodeID this2 = const_cast<UnaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = Curr_block;
      arg->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
        {
          temporary_terms.insert(this2);
          ModelBlock->Block_List[first_occurence[this2]].Temporary_terms->insert(this2);
        }
    }
}

void
UnaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                         const temporary_terms_type &temporary_terms) const
{
  // If node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<UnaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type != oCDynamicModelSparseDLL)
        output << "T" << idx;
      else
        output << "T" << idx << "[it_]";
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
      || (op_code == oUminus
          && arg->precedence(output_type, temporary_terms) < precedence(output_type, temporary_terms)))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write argument
  arg->writeOutput(output, output_type, temporary_terms);

  if (close_parenthesis)
    output << ")";

  // Close parenthesis for uminus
  if (op_code == oUminus)
    output << ")";
}

double
UnaryOpNode::eval_opcode(UnaryOpcode op_code, double v) throw (EvalException)
{
  switch(op_code)
    {
    case oUminus:
      return(-v);
    case oExp:
      return(exp(v));
    case oLog:
      return(log(v));
    case oLog10:
      return(log10(v));
    case oCos:
      return(cos(v));
    case oSin:
      return(sin(v));
    case oTan:
      return(tan(v));
    case oAcos:
      return(acos(v));
    case oAsin:
      return(asin(v));
    case oAtan:
      return(atan(v));
    case oCosh:
      return(cosh(v));
    case oSinh:
      return(sinh(v));
    case oTanh:
      return(tanh(v));
    case oAcosh:
      return(acosh(v));
    case oAsinh:
      return(asinh(v));
    case oAtanh:
      return(atanh(v));
    case oSqrt:
      return(sqrt(v));
    }
  // Impossible
  throw EvalException();
}

double
UnaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v = arg->eval(eval_context);

  return eval_opcode(op_code, v);
}

void
UnaryOpNode::collectEndogenous(NodeID &Id)
{
  arg->collectEndogenous(Id);
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
BinaryOpNode::precedence(ExprNodeOutputType output_type, const temporary_terms_type &temporary_terms) const
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
      if (!OFFSET(output_type))
        // In C, power operator is of the form pow(a, b)
        return 100;
      else
        return 3;
    }
  cerr << "Impossible case!" << endl;
  exit(-1);
}

int
BinaryOpNode::cost(const temporary_terms_type &temporary_terms, bool is_matlab) const
{
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  // For a temporary term, the cost is null
  if (it != temporary_terms.end())
    return 0;

  int cost = arg1->cost(temporary_terms, is_matlab);
  cost += arg2->cost(temporary_terms, is_matlab);

  if (is_matlab)
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
                                    temporary_terms_type &temporary_terms,
                                    bool is_matlab) const
{
  NodeID this2 = const_cast<BinaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      // If this node has never been encountered, set its ref count to one,
      //  and travel through its children
      reference_count[this2] = 1;
      arg1->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, is_matlab);
    }
  else
    {
      // If the node has already been encountered, increment its ref count
      //  and declare it as a temporary term if it is too costly
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, is_matlab) > MIN_COST(is_matlab))
        temporary_terms.insert(this2);
    }
}

void
BinaryOpNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                    temporary_terms_type &temporary_terms,
                                    map<NodeID, int> &first_occurence,
                                    int Curr_block,
                                    Model_Block *ModelBlock) const
{
  NodeID this2 = const_cast<BinaryOpNode *>(this);
  map<NodeID, int>::iterator it = reference_count.find(this2);
  if (it == reference_count.end())
    {
      reference_count[this2] = 1;
      first_occurence[this2] = Curr_block;
      arg1->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock);
      arg2->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, Curr_block, ModelBlock);
    }
  else
    {
      reference_count[this2]++;
      if (reference_count[this2] * cost(temporary_terms, false) > MIN_COST_C)
        {
          temporary_terms.insert(this2);
          ModelBlock->Block_List[first_occurence[this2]].Temporary_terms->insert(this2);
        }
    }
}

double
BinaryOpNode::eval_opcode(double v1, BinaryOpcode op_code, double v2) throw (EvalException)
{
  switch(op_code)
    {
    case oPlus:
      return(v1 + v2);
    case oMinus:
      return(v1 - v2);
    case oTimes:
      return(v1 * v2);
    case oDivide:
      return(v1 / v2);
    case oPower:
      return(pow(v1, v2));
    case oEqual:
    default:
      throw EvalException();
    }
}

double
BinaryOpNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  double v1 = arg1->eval(eval_context);
  double v2 = arg2->eval(eval_context);

  return eval_opcode(v1, op_code, v2);
}

void
BinaryOpNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                          const temporary_terms_type &temporary_terms) const
{
  // If current node is a temporary term
  temporary_terms_type::const_iterator it = temporary_terms.find(const_cast<BinaryOpNode *>(this));
  if (it != temporary_terms.end())
    {
      if (output_type != oCDynamicModelSparseDLL)
        output << "T" << idx;
      else
        output << "T" << idx << "[it_]";
      return;
    }

  // Treat special case of power operator in C
  if (op_code == oPower && (!OFFSET(output_type)))
    {
      output << "pow(";
      arg1->writeOutput(output, output_type, temporary_terms);
      output << ",";
      arg2->writeOutput(output, output_type, temporary_terms);
      output << ")";
      return;
    }

  int prec = precedence(output_type, temporary_terms);

  bool close_parenthesis = false;

  // If left argument has a lower precedence, or if current and left argument are both power operators, add parenthesis around left argument
  BinaryOpNode *barg1 = dynamic_cast<BinaryOpNode *>(arg1);
  if (arg1->precedence(output_type, temporary_terms) < prec
      || (op_code == oPower && barg1 != NULL && barg1->op_code == oPower))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write left argument
  arg1->writeOutput(output, output_type, temporary_terms);

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
  int arg2_prec = arg2->precedence(output_type, temporary_terms);
  if (arg2_prec < prec
      || (op_code == oPower && barg2 != NULL && barg2->op_code == oPower)
      || (op_code == oMinus && arg2_prec == prec)
      || (op_code == oDivide && arg2_prec == prec))
    {
      output << "(";
      close_parenthesis = true;
    }

  // Write right argument
  arg2->writeOutput(output, output_type, temporary_terms);

  if (close_parenthesis)
    output << ")";
}

void
BinaryOpNode::collectEndogenous(NodeID &Id)
{
  arg1->collectEndogenous(Id);
  arg2->collectEndogenous(Id);
}

UnknownFunctionNode::UnknownFunctionNode(DataTree &datatree_arg,
                                         const string &function_name_arg,
                                         const vector<NodeID> &arguments_arg) :
  ExprNode(datatree_arg),
  function_name(function_name_arg),
  arguments(arguments_arg)
{
}

NodeID
UnknownFunctionNode::computeDerivative(int varID)
{
  cerr << "UnknownFunctionNode::computeDerivative: operation impossible!" << endl;
  exit(-1);
}

void
UnknownFunctionNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                           temporary_terms_type &temporary_terms,
                                           bool is_matlab) const
{
  cerr << "UnknownFunctionNode::computeTemporaryTerms: operation impossible!" << endl;
  exit(-1);
}

void UnknownFunctionNode::writeOutput(ostream &output, ExprNodeOutputType output_type,
                                      const temporary_terms_type &temporary_terms) const
{
  output << function_name << "(";
  for(vector<NodeID>::const_iterator it = arguments.begin();
      it != arguments.end(); it++)
    {
      if (it != arguments.begin())
        output << ",";

      (*it)->writeOutput(output, output_type, temporary_terms);
    }
  output << ")";
}

void
UnknownFunctionNode::computeTemporaryTerms(map<NodeID, int> &reference_count,
                                           temporary_terms_type &temporary_terms,
                                           map<NodeID, int> &first_occurence,
                                           int Curr_block,
                                           Model_Block *ModelBlock) const
{
  cerr << "UnknownFunctionNode::computeTemporaryTerms: not implemented" << endl;
  exit(-1);
}

void
UnknownFunctionNode::collectEndogenous(NodeID &Id)
{
  cerr << "UnknownFunctionNode::collectEndogenous: not implemented" << endl;
  exit(-1);
}

double
UnknownFunctionNode::eval(const eval_context_type &eval_context) const throw (EvalException)
{
  cout << "Unknown function\n";
  cout.flush();
  throw EvalException();
}
