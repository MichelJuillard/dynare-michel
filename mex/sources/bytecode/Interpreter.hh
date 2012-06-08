/*
 * Copyright (C) 2007-2012 Dynare Team
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
# include <dynmex.h>
#else
# include "mex_interface.hh"
#endif

//#define DEBUGC

using namespace std;

#define pow_ pow

class Interpreter : public SparseMatrix
{
private:
  unsigned int EQN_dvar1, EQN_dvar2, EQN_dvar3;
  int EQN_lag1, EQN_lag2, EQN_lag3;
  mxArray *GlobalTemporaryTerms;
protected:
  double pow1(double a, double b);
  double divide(double a, double b);
  double log1(double a);
  double log10_1(double a);
  void compute_block_time(int Per_u_, bool evaluate, int block_num, int size, bool steady_state);
  void evaluate_a_block(const int size, const int type, string bin_basename, bool steady_state, int block_num,
                        const bool is_linear = false, const int symbol_table_endo_nbr = 0, const int Block_List_Max_Lag = 0, const int Block_List_Max_Lead = 0, const int u_count_int = 0, int block = -1);
  int simulate_a_block(const int size, const int type, string file_name, string bin_basename, bool Gaussian_Elimination, bool steady_state, bool print_it, int block_num,
                       const bool is_linear = false, const int symbol_table_endo_nbr = 0, const int Block_List_Max_Lag = 0, const int Block_List_Max_Lead = 0, const int u_count_int = 0);
  void print_a_block(const int size, const int type, string bin_basename, bool steady_state, int block_num,
                     const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag,
                     const int Block_List_Max_Lead, const int u_count_int, int block);
  void SingularDisplay(int Per_u_, bool evaluate, int Block_Count, int size, bool steady_state, it_code_type begining);
  vector<Block_contain_type> Block_Contain;
  code_liste_type code_liste;
  it_code_type it_code;
  int Block_Count, Per_u_, Per_y_;
  int it_, maxit_, size_of_direction;
  double solve_tolf;
  bool GaussSeidel;
  map<pair<pair<int, int>, int>, int> IM_i;
  int equation, derivative_equation, derivative_variable;
  string filename;
  int minimal_solving_periods;
  int stack_solve_algo, solve_algo;
  bool global_temporary_terms;
  bool print, print_error;
public:
  ~Interpreter();
  Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg, double *steady_x_arg,
              double *direction_arg, int y_size_arg, int nb_row_x_arg,
              int nb_row_xd_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg, int maxit_arg_, double solve_tolf_arg, int size_o_direction_arg,
              double slowc_arg, int y_decal_arg, double markowitz_c_arg, string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
              bool global_temporary_terms_arg, bool print_arg, bool print_error_arg, mxArray *GlobalTemporaryTerms_arg);
  bool compute_blocks(string file_name, string bin_basename, bool steady_state, bool evaluate, int block, int &nb_blocks, bool print_it);
  inline mxArray *
  get_jacob(int block_num)
  {
    return jacobian_block[block_num];
  };
  inline mxArray *
  get_jacob_exo(int block_num)
  {
    return jacobian_exo_block[block_num];
  };
  inline mxArray *
  get_jacob_exo_det(int block_num)
  {
    return jacobian_det_exo_block[block_num];
  };
  inline mxArray *
  get_jacob_other_endo(int block_num)
  {
    return jacobian_other_endo_block[block_num];
  };
  inline vector<double>
  get_residual()
  {
    return residual;
  };
  inline mxArray *
  get_Temporary_Terms()
  {
    return GlobalTemporaryTerms;
  };
};

#endif
