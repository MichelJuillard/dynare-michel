#ifndef _INTERPRETE_HH
#define _INTERPRETE_HH
#include <stack>
#include <map>
#include <string>
#include <iostream>
#include "SymbolTable.hh"
#include "VariableTable.hh"
#include "ExprNode.hh"
//#define PRINT_IT


using namespace std;
typedef struct t_val_index
{
  int index,type,indexed;
  double value;
};

typedef struct t_val1_index
{
  string name;
  double value;
};


typedef stack<double> STACK_SIMPLE00;
typedef map<string, t_val_index, less<string > > t_map_to_val_index;
typedef map<int, t_val1_index > t_map_int;

class interprete
{
public:
  double cutoff;
  bool eval;
  double u1, u2;
  STACK_SIMPLE00 Stack;
  interprete();
  double S_to_Val(string *str);
  double get_value(string *str, double ll);
  void put_value(string *str, int num,int Type, double val);
  void print_all();
  void create_id_map(int* Table, int Size, int HSize);
  double GetDataValue(/*NodeID*/int id, Type type);
  void set_cutoff(double r);
private:
  string tmp_put_value_name;
  t_map_to_val_index variable;
  t_map_int i_endo_variable, i_exo_variable, i_param_variable;
};
#endif
