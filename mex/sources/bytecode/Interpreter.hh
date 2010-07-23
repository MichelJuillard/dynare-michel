/*
 * Copyright (C) 2007-2009 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INTERPRETER_HH_INCLUDED
#define INTERPRETER_HH_INCLUDED

#include <stack>
#include <vector>
#include <string>
#include <cmath>
#define BYTE_CODE
#include "CodeInterpreter.hh"
#include "SparseMatrix.hh"
#ifdef LINBCG
# include "linbcg.hh"
#endif
#ifndef DEBUG_EX
# include "mex.h"
#else
# include "mex_interface.hh"
#endif

//#define DEBUGC

using namespace std;

#define pow_ pow
typedef vector<pair<Tags, void * > > code_liste_type;
typedef code_liste_type::const_iterator it_code_type;

class Interpreter : SparseMatrix
{
private:
  #ifndef DEBUG_EX
  vector<mxArray*> jacobian_block, jacobian_other_endo_block, jacobian_exo_block, jacobian_det_exo_block;
  #endif
  ExpressionType EQN_type;
  unsigned int EQN_equation, EQN_block, EQN_block_number;
  unsigned int EQN_dvar1, EQN_dvar2, EQN_dvar3;
  int EQN_lag1, EQN_lag2, EQN_lag3;
  it_code_type it_code_expr;
  protected:
  double pow1(double a, double b, bool evaluate);
  double log1(double a, bool evaluate);
  double log10_1(double a, bool evaluate);
  string remove_white(string str);
  string add_underscore_to_fpe(string str);
  string get_variable(SymbolType variable_type, unsigned int variable_num);
  string error_location(bool evaluate);
  void compute_block_time(int Per_u_, bool evaluate, int block_num, int size, bool steady_state);
  string print_expression(it_code_type it_code, bool evaluate);
  void evaluate_a_block(const int size, const int type, string bin_basename, bool steady_state, int block_num,
                        const bool is_linear = false, const int symbol_table_endo_nbr = 0, const int Block_List_Max_Lag = 0, const int Block_List_Max_Lead = 0, const int u_count_int = 0);
  bool simulate_a_block(const int size, const int type, string file_name, string bin_basename, bool Gaussian_Elimination, bool steady_state, int block_num,
                        const bool is_linear = false, const int symbol_table_endo_nbr = 0, const int Block_List_Max_Lag = 0, const int Block_List_Max_Lead = 0, const int u_count_int = 0);
  double *T;
  vector<Block_contain_type> Block_Contain;
  code_liste_type code_liste;
  it_code_type it_code;
  stack<double> Stack;
  int Block_Count, Per_u_, Per_y_;
  int it_, nb_row_x, nb_row_xd, maxit_, size_of_direction;
  double *g1, *r;
  double solve_tolf;
  bool GaussSeidel;
  double *x, *params;
  double *steady_y, *steady_x;
  map<pair<pair<int, int>, int>, int> IM_i;
  int equation, derivative_equation, derivative_variable;
  string filename;
  int minimal_solving_periods;
  int stack_solve_algo, solve_algo;
public:
  ~Interpreter();
  Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg, double *steady_x_arg,
              double *direction_arg, int y_size_arg, int nb_row_x_arg,
              int nb_row_xd_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg, int maxit_arg_, double solve_tolf_arg, int size_o_direction_arg,
              double slowc_arg, int y_decal_arg, double markowitz_c_arg, string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg);
  bool compute_blocks(string file_name, string bin_basename, bool steady_state, bool evaluate, bool block, int &nb_blocks);
  #ifndef DEBUG_EX
  inline mxArray* get_jacob(int block_num) {return jacobian_block[block_num];};
  inline mxArray* get_jacob_exo(int block_num) {return jacobian_exo_block[block_num];};
  inline mxArray* get_jacob_exo_det(int block_num) {return jacobian_det_exo_block[block_num];};
  inline mxArray* get_jacob_other_endo(int block_num) {return jacobian_other_endo_block[block_num];};
  #endif
};

#endif
