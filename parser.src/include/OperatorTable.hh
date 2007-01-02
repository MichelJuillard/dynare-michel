#ifndef _OPERATORTABLE_HH
#define _OPERATORTABLE_HH

using namespace std;

#include <map>
#include <string>

#include "DynareBison.hh"

//! Shortcut to access tokens defined by Bison
typedef yy::parser::token token;

//! Stores informations about operators
class OperatorTable
{
private:
  //! Stores informations about a given operator
  struct Operator
  {
    //! Operator name for both Matlab and C
    string str;
    //! Operator precedence
    int precedence;
    //! True if operator is a function
    bool isfunction;
    //! Time computation cost of operator in C
    int c_cost;
    //! Time computation cost of operator in Matlab
    int m_cost;
  };
  //! Type of operators map indexed by code
  typedef map<int, Operator> operator_map;
  //! Operator table
  static operator_map operator_table;
  //! Has the table been initialized ?
  static bool initialized;
  //! Initialize the table (or does nothing if it has already be done)
  static void init();
public:
  //! Get operator name
  static string str(int op_code);
  //! Get operator precedence
  static int precedence(int op_code);
  //! Is operator a function ?
  static bool isfunction(int op_code);
  //! Get cost of operator
  static int cost(int op_code, int offset);
};
#endif
