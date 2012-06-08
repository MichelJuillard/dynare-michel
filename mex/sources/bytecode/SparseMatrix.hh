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

#ifndef SPARSEMATRIX_HH_INCLUDED
#define SPARSEMATRIX_HH_INCLUDED

#include <fstream>
#include <stack>
#include <cmath>
#include <map>
#include <ctime>

#ifdef OCTAVE_MEX_FILE
# define CHAR_LENGTH 1
#else
# define CHAR_LENGTH 2
#endif

#include "Mem_Mngr.hh"
#include "ErrorHandling.hh"
#define NEW_ALLOC
#define MARKOVITZ

using namespace std;

struct t_save_op_s
{
  short int lag, operat;
  int first, second;
};

const int IFLD  = 0;
const int IFDIV = 1;
const int IFLESS = 2;
const int IFSUB = 3;
const int IFLDZ = 4;
const int IFMUL = 5;
const int IFSTP = 6;
const int IFADD = 7;
const double eps = 1e-10;
const double very_big = 1e24;
const int alt_symbolic_count_max = 1;
const double mem_increasing_factor = 1.1;

class SparseMatrix : public ErrorMsg
{
public:
  SparseMatrix();
  void Simulate_Newton_Two_Boundaries(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, bool cvg, int &iter, int minimal_solving_periods, int stack_solve_algo, unsigned int endo_name_length, char *P_endo_names) /*throw(ErrorHandlingException)*/;
  bool Simulate_Newton_One_Boundary(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, bool print_it, bool cvg, int &iter, bool steady_state, int stack_solve_algo, int solve_algo);
  void Direct_Simulate(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, int iter);
  void fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1);
  void Read_SparseMatrix(string file_name, const int Size, int periods, int y_kmin, int y_kmax, bool steady_state, bool two_boundaries, int stack_solve_algo, int solve_algo);
  void Read_file(string file_name, int periods, int u_size1, int y_size, int y_kmin, int y_kmax, int &nb_endo, int &u_count, int &u_count_init, double *u);
  void Singular_display(int block, int Size, bool steady_state, it_code_type it_code);
  double g0, gp0, glambda2, try_at_iteration;

private:
  void Init_GE(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM);
  void Init_Matlab_Sparse(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM, mxArray *A_m, mxArray *b_m, mxArray *x0_m);
  void Init_Matlab_Sparse_Simple(int Size, map<pair<pair<int, int>, int>, int> &IM, mxArray *A_m, mxArray *b_m, bool &zero_solution, mxArray *x0_m);
  void Simple_Init(int Size, std::map<std::pair<std::pair<int, int>, int>, int> &IM, bool &zero_solution);
  void End_GE(int Size);
  void Solve_ByteCode_Symbolic_Sparse_GaussianElimination(int Size, bool symbolic, int Block_number);
  bool Solve_ByteCode_Sparse_GaussianElimination(int Size, int blck, bool steady_state, int it_);
  void Solve_Matlab_Relaxation(mxArray *A_m, mxArray *b_m, unsigned int Size, double slowc_l, bool is_two_boundaries, int  it_);
  void Solve_Matlab_LU_UMFPack(mxArray *A_m, mxArray *b_m, int Size, double slowc_l, bool is_two_boundaries, int it_);
  void Solve_Matlab_GMRES(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, bool steady_state, mxArray *x0_m);
  void Solve_Matlab_BiCGStab(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m, bool steady_state);
  bool compare(int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long int nop4,  int Size
#ifdef PROFILER
               , long int *ndiv, long int *nsub
#endif
               );
  void Insert(const int r, const int c, const int u_index, const int lag_index);
  void Delete(const int r, const int c);
  int At_Row(int r, NonZeroElem **first);
  int At_Pos(int r, int c, NonZeroElem **first);
  int At_Col(int c, NonZeroElem **first);
  int At_Col(int c, int lag, NonZeroElem **first);
  int NRow(int r);
  int NCol(int c);
  int Union_Row(int row1, int row2);
  void Print(int Size, int *b);
  int Get_u();
  void Delete_u(int pos);
  void Clear_u();
  void Print_u();
  void CheckIt(int y_size, int y_kmin, int y_kmax, int Size, int periods, int iter);
  void Check_the_Solution(int periods, int y_kmin, int y_kmax, int Size, double *u, int *pivot, int *b);
  int complete(int beg_t, int Size, int periods, int *b);
  void bksub(int tbreak, int last_period, int Size, double slowc_l
#ifdef PROFILER
               , long int *nmul
#endif
               );
  void simple_bksub(int it_, int Size, double slowc_l);
  mxArray *Sparse_transpose(mxArray *A_m);
  mxArray *Sparse_mult_SAT_SB(mxArray *A_m, mxArray *B_m);
  mxArray *Sparse_mult_SAT_B(mxArray *A_m, mxArray *B_m);
  mxArray *mult_SAT_B(mxArray *A_m, mxArray *B_m);
  mxArray *Sparse_substract_SA_SB(mxArray *A_m, mxArray *B_m);
  mxArray *Sparse_substract_A_SB(mxArray *A_m, mxArray *B_m);
  mxArray *substract_A_B(mxArray *A_m, mxArray *B_m);

  stack<double> Stack;
  int nb_prologue_table_u, nb_first_table_u, nb_middle_table_u, nb_last_table_u;
  int nb_prologue_table_y, nb_first_table_y, nb_middle_table_y, nb_last_table_y;
  int middle_count_loop;
  char type;
  fstream SaveCode;
  string filename;
  int max_u, min_u;
  clock_t time00;

  Mem_Mngr mem_mngr;
  vector<int> u_liste;
  map<pair<int, int>, NonZeroElem *> Mapped_Array;
  int *NbNZRow, *NbNZCol;
  NonZeroElem **FNZE_R, **FNZE_C;
  int u_count_init;

  int *pivot, *pivotk, *pivot_save;
  double *pivotv, *pivotva;
  int *b;
  bool *line_done;
  bool symbolic, alt_symbolic;
  int alt_symbolic_count;
  int *g_save_op;
  int first_count_loop;
  int g_nop_all;
  double markowitz_c_s;
  double res1a;
  long int nop_all, nop1, nop2;
  map<pair<pair<int, int>, int>, int> IM_i;
protected:
  vector<double> residual;
  int u_count_alloc, u_count_alloc_save;
  vector<double *> jac;
  double *jcb;
  double res1, res2, max_res;
  int max_res_idx;
  double slowc, slowc_save, prev_slowc_save, markowitz_c;
  int y_decal;
  int  *index_vara, *index_equa;
  int u_count, tbreak_g;
  int iter;
  double *direction;
  int start_compare;
  int restart;
  bool error_not_printed;
  double g_lambda1, g_lambda2, gp_0;
  double lu_inc_tol;
};

#endif
