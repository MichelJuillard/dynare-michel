#ifndef _EXPRESSION_HH
#define _EXPRESSION_HH
//------------------------------------------------------------------------------
/** \file
 * \version 1.0
 * \date 04/13/2004
 * \par This file defines the Expression class .
 */
//------------------------------------------------------------------------------
#include <map>
#include <vector>
#include <sstream>

#include "OperatorTable.hh"
#include "NumericalConstants.hh"

struct Token
{
  /*! ID of first operand */
  int id1;
  /*! Type of first operand  */
  Type type1;
  /*! Flag : operand 1 is followed or not */
  bool followed1;
  /*! ID of second operand  */
  int id2;
  /*! Type of second operand  */
  Type type2;
  /*! Flag : operand 2 is followed or not */
  bool followed2;
  /*! Operator code */
  int op_code;
  /*! Operator name */
  std::string op_name;
  /*! costructor */
  Token(){followed1=followed2=false;};
};
//------------------------------------------------------------------------------
/*!
  \class Expression
  \brief Handles expressions appearing
  in initialization statements. These expressions aren't meant to be derived
*/
class Expression
{
private :
  /*! Vector tokens */
  std::vector<Token>   expression_list;
  /*! Operator table : names and precedences */
  OperatorTable operator_table;
  /*! Output string of the class */
  std::ostringstream output;
  //! Pointer to numerical constants table
  NumericalConstants *num_constants;

public :
  /*! Constructor */
  Expression();
  /*! Destructor */
  ~Expression();
  //! Set numerical constants pointer
  void setNumericalConstants(NumericalConstants *num_constants_arg);
  /*! Adds binary token to expression list */
  int   AddToken(int id1,Type type1, int id2,Type type2,int op_code);
  /*! Adds unary token to expression list */
  int   AddToken(int id1,Type type1, int op_code);
  /*! Adds unkown function to expression list */
  int   AddToken(int id1,Type type1, std::string ufunction);
  /*! Returns output string */
  std::string get();
  /*! Clear expression list */
  void  clear(void);
  /*! Print expression to output string */
  void  set(void);
  /*! Gets output argument name */
  std::string getArgument(Type type, int id);
};
//------------------------------------------------------------------------------
#endif
