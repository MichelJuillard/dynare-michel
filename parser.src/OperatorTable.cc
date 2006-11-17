/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the OperatorTable class methodes.
*/
//------------------------------------------------------------------------------
#include <map>
#include <string>
using namespace std;
//------------------------------------------------------------------------------
#include "OperatorTable.hh"
//------------------------------------------------------------------------------
OperatorTable::OperatorTable()
{
  operator_table[token::COMMA].str = ",";
  operator_table[token::EQUAL].str   = "=";
  operator_table[token::PLUS].str  = "+";
  operator_table[token::MINUS].str   = "-";
  operator_table[token::TIMES].str   = "*";
  operator_table[token::DIVIDE].str  = "/";
  operator_table[token::UMINUS].str  = "-";
  operator_table[token::POWER].str   = "^";
  operator_table[token::EXP].str   = "exp";
  operator_table[token::LOG].str   = "log";
  operator_table[token::LOG10].str   = "log10";
  operator_table[token::COS].str   = "cos";
  operator_table[token::SIN].str   = "sin";
  operator_table[token::TAN].str   = "tan";
  operator_table[token::ACOS].str  = "acos";
  operator_table[token::ASIN].str  = "asin";
  operator_table[token::ATAN].str  = "atan";
  operator_table[token::COSH].str  = "cosh";
  operator_table[token::SINH].str  = "sinh";
  operator_table[token::TANH].str  = "tanh";
  operator_table[token::ACOSH].str   = "acosh";
  operator_table[token::ASINH].str   = "asinh";
  operator_table[token::ATANH].str   = "atanh";
  operator_table[token::SQRT].str  = "sqrt";
  operator_table[token::NAME].str  = "";

  operator_table[token::COMMA].precedence  = -1;
  operator_table[token::EQUAL].precedence  = 0;
  operator_table[token::PLUS].precedence   = 0;
  operator_table[token::MINUS].precedence  = 1;
  operator_table[token::TIMES].precedence  = 2;
  operator_table[token::DIVIDE].precedence   = 3;
  operator_table[token::UMINUS].precedence   = 4;
  operator_table[token::POWER].precedence  = 5;
  operator_table[token::EXP].precedence    =
    operator_table[token::LOG].precedence    =
    operator_table[token::LOG10].precedence  =
    operator_table[token::COS].precedence    =
    operator_table[token::SIN].precedence    =
    operator_table[token::TAN].precedence    =
    operator_table[token::ACOS].precedence   =
    operator_table[token::ASIN].precedence   =
    operator_table[token::ATAN].precedence   =
    operator_table[token::COSH].precedence   =
    operator_table[token::SINH].precedence   =
    operator_table[token::TANH].precedence   =
    operator_table[token::ACOSH].precedence  =
    operator_table[token::ASINH].precedence  =
    operator_table[token::ATANH].precedence  =
    operator_table[token::SQRT].precedence   =
    operator_table[token::NAME].precedence   = 6;

  // Operator costs for M files
  operator_table[token::COMMA].cost[1]   = 0;
  operator_table[token::EQUAL].cost[1]   = 0;
  operator_table[token::PLUS].cost[1]  = 90;
  operator_table[token::MINUS].cost[1]   = 90;
  operator_table[token::TIMES].cost[1]   = 90;
  operator_table[token::DIVIDE].cost[1]  = 990;
  operator_table[token::UMINUS].cost[1]  = 70;
  operator_table[token::POWER].cost[1]   = 1160;
  operator_table[token::EXP].cost[1]   = 160;
  operator_table[token::LOG].cost[1]   = 300;
  operator_table[token::LOG10].cost[1]   = 16000;
  operator_table[token::COS].cost[1]   = 210;
  operator_table[token::SIN].cost[1]   = 210;
  operator_table[token::TAN].cost[1]   = 230;
  operator_table[token::ACOS].cost[1]  = 300;
  operator_table[token::ASIN].cost[1]  = 310;
  operator_table[token::ATAN].cost[1]  = 140;
  operator_table[token::COSH].cost[1]  = 210;
  operator_table[token::SINH].cost[1]  = 240;
  operator_table[token::TANH].cost[1]  = 190;
  operator_table[token::ACOSH].cost[1]   = 770;
  operator_table[token::ASINH].cost[1]   = 460;
  operator_table[token::ATANH].cost[1]   = 350;
  operator_table[token::SQRT].cost[1]  = 570;
  operator_table[token::NAME].cost[1]  = 0;

  // Operator costs for C files
  operator_table[token::COMMA].cost[0]   = 0;
  operator_table[token::EQUAL].cost[0]   = 0;
  operator_table[token::PLUS].cost[0]  = 4;
  operator_table[token::MINUS].cost[0]   = 4;
  operator_table[token::TIMES].cost[0]   = 4;
  operator_table[token::DIVIDE].cost[0]  = 15;
  operator_table[token::UMINUS].cost[0]  = 3;
  operator_table[token::POWER].cost[0]   = 520;
  operator_table[token::EXP].cost[0]   = 210;
  operator_table[token::LOG].cost[0]   = 137;
  operator_table[token::LOG10].cost[0]   = 139;
  operator_table[token::COS].cost[0]   = 160;
  operator_table[token::SIN].cost[0]   = 160;
  operator_table[token::TAN].cost[0]   = 170;
  operator_table[token::ACOS].cost[0]  = 190;
  operator_table[token::ASIN].cost[0]  = 180;
  operator_table[token::ATAN].cost[0]  = 190;
  operator_table[token::COSH].cost[0]  = 240;
  operator_table[token::SINH].cost[0]  = 240;
  operator_table[token::TANH].cost[0]  = 240;
  operator_table[token::ACOSH].cost[0]   = 210;
  operator_table[token::ASINH].cost[0]   = 220;
  operator_table[token::ATANH].cost[0]   = 150;
  operator_table[token::SQRT].cost[0]  = 90;
  operator_table[token::NAME].cost[0]  = 0;

  operator_table[token::COMMA].isfunction  = false;
  operator_table[token::EQUAL].isfunction  = false;
  operator_table[token::PLUS].isfunction   = false;
  operator_table[token::MINUS].isfunction  = false;
  operator_table[token::TIMES].isfunction  = false;
  operator_table[token::DIVIDE].isfunction   = false;
  operator_table[token::UMINUS].isfunction   = false;
  operator_table[token::POWER].isfunction  = false;
  operator_table[token::EXP].isfunction    = true;
  operator_table[token::LOG].isfunction    = true;
  operator_table[token::LOG10].isfunction  = true;
  operator_table[token::COS].isfunction    = true;
  operator_table[token::SIN].isfunction    = true;
  operator_table[token::TAN].isfunction    = true;
  operator_table[token::ACOS].isfunction   = true;
  operator_table[token::ASIN].isfunction   = true;
  operator_table[token::ATAN].isfunction   = true;
  operator_table[token::COSH].isfunction   = true;
  operator_table[token::SINH].isfunction   = true;
  operator_table[token::TANH].isfunction   = true;
  operator_table[token::ACOSH].isfunction  = true;
  operator_table[token::ASINH].isfunction  = true;
  operator_table[token::ATANH].isfunction  = true;
  operator_table[token::SQRT].isfunction   = true;
  operator_table[token::NAME].isfunction   = true;
}

//------------------------------------------------------------------------------
OperatorTable::~OperatorTable()
{
  // Empty
}

//------------------------------------------------------------------------------
