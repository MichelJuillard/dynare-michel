#include "OperatorTable.hh"

OperatorTable::operator_map OperatorTable::operator_table;

bool OperatorTable::initialized = false;

void
OperatorTable::init()
{
  if (initialized)
    return;
  initialized = true;

  operator_table[token::COMMA].str = ",";
  operator_table[token::EQUAL].str = "=";
  operator_table[token::PLUS].str = "+";
  operator_table[token::MINUS].str = "-";
  operator_table[token::TIMES].str = "*";
  operator_table[token::DIVIDE].str = "/";
  operator_table[token::UMINUS].str = "-";
  operator_table[token::POWER].str = "^";
  operator_table[token::EXP].str = "exp";
  operator_table[token::LOG].str = "log";
  operator_table[token::LOG10].str = "log10";
  operator_table[token::COS].str = "cos";
  operator_table[token::SIN].str = "sin";
  operator_table[token::TAN].str = "tan";
  operator_table[token::ACOS].str = "acos";
  operator_table[token::ASIN].str = "asin";
  operator_table[token::ATAN].str = "atan";
  operator_table[token::COSH].str = "cosh";
  operator_table[token::SINH].str = "sinh";
  operator_table[token::TANH].str = "tanh";
  operator_table[token::ACOSH].str = "acosh";
  operator_table[token::ASINH].str = "asinh";
  operator_table[token::ATANH].str = "atanh";
  operator_table[token::SQRT].str = "sqrt";
  operator_table[token::NAME].str = "";

  operator_table[token::COMMA].precedence = -1;
  operator_table[token::EQUAL].precedence = 0;
  operator_table[token::PLUS].precedence = 0;
  operator_table[token::MINUS].precedence = 1;
  operator_table[token::TIMES].precedence = 2;
  operator_table[token::DIVIDE].precedence = 3;
  operator_table[token::UMINUS].precedence = 4;
  operator_table[token::POWER].precedence = 5;
  operator_table[token::EXP].precedence =
    operator_table[token::LOG].precedence =
    operator_table[token::LOG10].precedence =
    operator_table[token::COS].precedence =
    operator_table[token::SIN].precedence =
    operator_table[token::TAN].precedence =
    operator_table[token::ACOS].precedence =
    operator_table[token::ASIN].precedence =
    operator_table[token::ATAN].precedence =
    operator_table[token::COSH].precedence =
    operator_table[token::SINH].precedence =
    operator_table[token::TANH].precedence =
    operator_table[token::ACOSH].precedence =
    operator_table[token::ASINH].precedence =
    operator_table[token::ATANH].precedence =
    operator_table[token::SQRT].precedence =
    operator_table[token::NAME].precedence = 6;

  // Operator costs for M files
  operator_table[token::COMMA].m_cost = 0;
  operator_table[token::EQUAL].m_cost = 0;
  operator_table[token::PLUS].m_cost = 90;
  operator_table[token::MINUS].m_cost = 90;
  operator_table[token::TIMES].m_cost = 90;
  operator_table[token::DIVIDE].m_cost = 990;
  operator_table[token::UMINUS].m_cost = 70;
  operator_table[token::POWER].m_cost = 1160;
  operator_table[token::EXP].m_cost = 160;
  operator_table[token::LOG].m_cost = 300;
  operator_table[token::LOG10].m_cost = 16000;
  operator_table[token::COS].m_cost = 210;
  operator_table[token::SIN].m_cost = 210;
  operator_table[token::TAN].m_cost = 230;
  operator_table[token::ACOS].m_cost = 300;
  operator_table[token::ASIN].m_cost = 310;
  operator_table[token::ATAN].m_cost = 140;
  operator_table[token::COSH].m_cost = 210;
  operator_table[token::SINH].m_cost = 240;
  operator_table[token::TANH].m_cost = 190;
  operator_table[token::ACOSH].m_cost = 770;
  operator_table[token::ASINH].m_cost = 460;
  operator_table[token::ATANH].m_cost = 350;
  operator_table[token::SQRT].m_cost = 570;
  operator_table[token::NAME].m_cost = 0;

  // Operator costs for C files
  operator_table[token::COMMA].c_cost = 0;
  operator_table[token::EQUAL].c_cost = 0;
  operator_table[token::PLUS].c_cost = 4;
  operator_table[token::MINUS].c_cost = 4;
  operator_table[token::TIMES].c_cost = 4;
  operator_table[token::DIVIDE].c_cost = 15;
  operator_table[token::UMINUS].c_cost = 3;
  operator_table[token::POWER].c_cost = 520;
  operator_table[token::EXP].c_cost = 210;
  operator_table[token::LOG].c_cost = 137;
  operator_table[token::LOG10].c_cost = 139;
  operator_table[token::COS].c_cost = 160;
  operator_table[token::SIN].c_cost = 160;
  operator_table[token::TAN].c_cost = 170;
  operator_table[token::ACOS].c_cost = 190;
  operator_table[token::ASIN].c_cost = 180;
  operator_table[token::ATAN].c_cost = 190;
  operator_table[token::COSH].c_cost = 240;
  operator_table[token::SINH].c_cost = 240;
  operator_table[token::TANH].c_cost = 240;
  operator_table[token::ACOSH].c_cost = 210;
  operator_table[token::ASINH].c_cost = 220;
  operator_table[token::ATANH].c_cost = 150;
  operator_table[token::SQRT].c_cost = 90;
  operator_table[token::NAME].c_cost = 0;

  operator_table[token::COMMA].isfunction = false;
  operator_table[token::EQUAL].isfunction = false;
  operator_table[token::PLUS].isfunction = false;
  operator_table[token::MINUS].isfunction = false;
  operator_table[token::TIMES].isfunction = false;
  operator_table[token::DIVIDE].isfunction = false;
  operator_table[token::UMINUS].isfunction = false;
  operator_table[token::POWER].isfunction = false;
  operator_table[token::EXP].isfunction = true;
  operator_table[token::LOG].isfunction = true;
  operator_table[token::LOG10].isfunction = true;
  operator_table[token::COS].isfunction = true;
  operator_table[token::SIN].isfunction = true;
  operator_table[token::TAN].isfunction = true;
  operator_table[token::ACOS].isfunction = true;
  operator_table[token::ASIN].isfunction = true;
  operator_table[token::ATAN].isfunction = true;
  operator_table[token::COSH].isfunction = true;
  operator_table[token::SINH].isfunction = true;
  operator_table[token::TANH].isfunction = true;
  operator_table[token::ACOSH].isfunction = true;
  operator_table[token::ASINH].isfunction = true;
  operator_table[token::ATANH].isfunction = true;
  operator_table[token::SQRT].isfunction = true;
  operator_table[token::NAME].isfunction = true;
}

string
OperatorTable::str(int op_code)
{
  init();
  return operator_table[op_code].str;
}

int
OperatorTable::precedence(int op_code)
{
  init();
  return operator_table[op_code].precedence;
}

bool
OperatorTable::isfunction(int op_code)
{
  init();
  return operator_table[op_code].isfunction;
}

int
OperatorTable::cost(int op_code, int offset)
{
  init();
  if (offset == 0)
    return operator_table[op_code].c_cost;
  else
    return operator_table[op_code].m_cost;
}
