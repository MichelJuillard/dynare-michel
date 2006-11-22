/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the Expression class methodes.
*/
//------------------------------------------------------------------------------
#include <stack>
using namespace std;
//------------------------------------------------------------------------------
#include "Expression.hh"
//------------------------------------------------------------------------------
ostringstream Expression::output;
//------------------------------------------------------------------------------
Expression::Expression()
{
  // Empty
}

//------------------------------------------------------------------------------
Expression::~Expression()
{
  // Empty
}

//------------------------------------------------------------------------------
int Expression::AddToken(int id1,Type type1,int id2,Type type2,int op_code)
{
  Token token;
  // Making token structure
  token.id1 = id1;
  token.type1 = type1;
  token.id2 = id2;
  token.type2 = type2;
  token.op_code = op_code;
  token.op_name = operator_table.str(op_code);
  //  Inserting token into expression_list
  expression_list.push_back(token);
  return expression_list.size() -1;
}

//------------------------------------------------------------------------------
int Expression::AddToken(int id1,Type type1,int op_code)
{
  Token token;
  // Making token structure
  token.id1 = id1;
  token.type1 = type1;
  token.id2 = -1;
  token.type2 = eUNDEF;
  token.op_code = op_code;
  token.op_name = operator_table.str(op_code);;
  //  Inserting token into expression_list
  expression_list.push_back(token);
  return expression_list.size() -1;
}

//------------------------------------------------------------------------------
int Expression::AddToken(int id1,Type type1, string ufunction)
{
  Token token;
  // Making token structure
  token.id1 = id1;
  token.type1 = type1;
  token.id2 = -1;
  token.type2 = eUNDEF;
  token.op_code = token::NAME;
  token.op_name = ufunction;
  //  Inserting token into expression_list
  expression_list.push_back(token);
  return expression_list.size() -1;
}

//------------------------------------------------------------------------------
void Expression::set(void)
{
  // Stack of temporary tokens
  stack <int, vector<Token> > stack_token;
  // Dtack of temporary expressions
  stack <int, vector<string> > stack_expression;
  // Temporary output
  ostringstream exp;
  // temporary variables for saving arguments and name oparator
  string argument1, argument2, op_name;
  // Define type for type operator (binary or unary)
  enum OperatorType
  {
    unary,
    binary
  };
  OperatorType op_type;
  int current_op, last_op;

  // Clearing output string
  output.str("");
  // Starting from the end of list
  stack_token.push(expression_list.back());
  // Main loop :
  // Repeat for last token from the stack
  // (1) 	if argument is temporary result, and not yet followed,
  //			set it as followed (flag) and push corresponding token
  //			on the token stack
  // (2) argument followed, or final argument
  //		(2.1) if argument is followed
  //			- set argument1 (or argument2) by last expression on
  //			expression tack
  //			- pop last expression from expression stack
  //		(2.2) if final argument
  //			  set argument1 (or argument2) by final argument
  // (3) set op_name by last token from the token stack
  // (3) pop last token from the token stack
  // (4) write temporary expression (using argument1, argument2
  //		and op_name) and push it on the expression stack
  // (5)

  while (stack_token.size() > 0)
    {
      // First argument is a temporary result,
      // pushing token on token stack and setting that argument to be followed
      if ((stack_token.top().type1 == eTempResult) &&
          (stack_token.top().followed1 == false))
        {
          stack_token.top().followed1 = true;
          stack_token.push(expression_list[stack_token.top().id1]);
        }
      // Second argument is a temporary result,
      // pushing token on stack and setting that argument to be followed
      else if ((stack_token.top().type2 == eTempResult) &&
               (stack_token.top().followed2 == false))
        {
          stack_token.top().followed2 = true;
          stack_token.push(expression_list[stack_token.top().id2]);
        }
      // Writing expression
      else
        {
          // Final token, no argment followed
          if ((stack_token.top().followed1 == false) &&
              (stack_token.top().followed2 == false))
            {
              argument1 = getArgument(stack_token.top().type1,stack_token.top().id1);
              current_op = stack_token.top().op_code;
              // Testing if unary or binary token
              if (stack_token.top().id2 >= 0)
                {
                  argument2 = getArgument(stack_token.top().type2,stack_token.top().id2);
                  op_type = binary;
                }
              else
                {
                  op_type = unary;
                }
            }
          // Both arguments are followed, writing stacked expressions
          else if ((stack_token.top().followed1 == true) &&
                   (stack_token.top().followed2 == true))
            {
              // Testing if unary or binary token
              if (stack_token.top().id2 >= 0)
                {
                  argument2 = stack_expression.top();
                  stack_expression.pop();
                  op_type = binary;

                }
              else
                {
                  op_type = unary;
                }
              argument1 = stack_expression.top();
              current_op = stack_token.top().op_code;
              stack_expression.pop();

            }
          // Only argument 1 is followed, combining expressions
          else if (stack_token.top().followed1 == true)
            {
              argument1 = stack_expression.top();
              current_op = stack_token.top().op_code;
              stack_expression.pop();
              // Testing if unary or binary token
              if (stack_token.top().id2 >= 0)
                {
                  argument2 = getArgument(stack_token.top().type2,stack_token.top().id2);
                  op_type = binary;
                }
              else
                {
                  op_type = unary;
                }
            }
          // Only argument 2 is followed, combining experssions
          else if (stack_token.top().followed2 == true)
            {
              argument1 = getArgument(stack_token.top().type1,stack_token.top().id1);
              argument2 = stack_expression.top();
              stack_expression.pop();
              current_op = stack_token.top().op_code;
              op_type = binary;
            }
          op_name = stack_token.top().op_name;
          exp.str("");
          stack_token.pop();
          // Saving last operator for the followed argument
          if (stack_token.size() > 0)
            {
              last_op = stack_token.top().op_code;
            }
          else
            {
              last_op = current_op;
            }
          if (op_type == binary)
            {
              // Comma operator, no parentheses
              // parentheses are writing with function operator
              if (current_op == token::COMMA)
                {
                  exp <<  argument1 << op_name << argument2 ;
                }
              else
                {
                  exp <<  "(" << argument1 << op_name << argument2 << ")";
                }
            }
          else
            {
              // Case of functions
              if (operator_table.isfunction(current_op) == true)
                exp <<  op_name << "(" << argument1 << ")";
              else
                exp <<  "(" << op_name << argument1 << ")";
            }
          stack_expression.push(exp.str());

        }
    }
  output << stack_expression.top();
  expression_list.clear();
}

//------------------------------------------------------------------------------
string Expression::getArgument(Type type,int id)
{
  ostringstream argument;

  if (type == eExogenous)
    {
      argument << "oo_.exo_steady_state"<< "(" << id+1 << ")";
    }
  else if (type == eExogenousDet)
    {
      argument <<  "oo_.exo_det_steady_state" << "(" << id+1 << ")";
    }
  else if (type == eEndogenous)
    {
      argument <<  "oo_.steady_state" << "(" << id+1 << ")";
    }
  else if (type == eParameter)
    {
      argument << "M_.params" << "(" << id+1 << ")";
    }
  else if (type == eNumericalConstant)
    {
      argument << NumericalConstants::get(id);
    }
  return argument.str();
}

//------------------------------------------------------------------------------
string Expression::get(void)
{
  return output.str();
}

//------------------------------------------------------------------------------
void Expression::clear(void)
{
  expression_list.clear();
}

//------------------------------------------------------------------------------
