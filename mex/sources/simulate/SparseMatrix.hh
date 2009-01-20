/*
 * Copyright (C) 2007-2008 Dynare Team
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
#include "Mem_Mngr.hh"
//! Openmp is available in GCC since version 4.3.2
//! Test if GCC version is greater then 4.3.2 order to avoid a copilation error
#define GNUVER 100*__GNUC__+10*__GNUC_MINOR__+__GNUC_PATCHLEVEL__
#if GNUVER >= 432
  #include <omp.h>
#endif
#define NEW_ALLOC
#define MARKOVITZ
//#define MEMORY_LEAKS

using namespace std;

struct t_save_op_s
{
  short int lag, operat;
  int first, second;
};

const int IFLD  =0;
const int IFDIV =1;
const int IFLESS=2;
const int IFSUB =3;
const int IFLDZ =4;
const int IFMUL =5;
const int IFSTP =6;
const int IFADD =7;
const double eps=1e-7;
const double very_big=1e24;
const int alt_symbolic_count_max=1;

struct t_table_y
  {
    int index, nb;
    int *u_index, *y_index;
  };

struct t_table_u
  {
    t_table_u* pNext;
    unsigned char type;
    int index;
    int op1, op2;
  };


class SparseMatrix
  {
  public:
    SparseMatrix();
    int simulate_NG1(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, bool cvg, int &iter);
    int simulate_NG(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, bool print_it, bool cvg, int &iter);
    void Direct_Simulate(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, int iter);
    void fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1);
    //void initialize(int periods_arg, int nb_endo_arg, int y_kmin_arg, int y_kmax_arg, int y_size_arg, int u_count_arg, int u_count_init_arg, double *u_arg, double *y_arg, double *ya_arg, double slowc_arg, int y_decal_arg, double markowitz_c_arg, double res1_arg, double res2_arg, double max_res_arg);
    void Read_SparseMatrix(string file_name, int Size, int periods, int y_kmin, int y_kmax);
    void Read_file(string file_name, int periods, int u_size1, int y_size, int y_kmin, int y_kmax, int &nb_endo, int &u_count, int &u_count_init, double* u);

 private:
    void Init(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int> ,int>, int> IM);
    void ShortInit(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int> ,int>, int> IM);
    void Simple_Init(int it_, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM);
    void End(int Size);
    bool compare( int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long int nop4,  int Size
#ifdef PROFILER
                 , long int *ndiv, long int *nsub
#endif
                );
    void Insert(int r, int c, int u_index, int lag_index);
    void Delete(int r,int c, int Size, int *b);
    int At_Row(int r, NonZeroElem **first);
    int At_Pos(int r, int c, NonZeroElem **first);
    int At_Col(int c, NonZeroElem **first);
    int At_Col(int c, int lag, NonZeroElem **first);
    int NRow(int r);
    int NCol(int c);
    int Union_Row(int row1, int row2);
    void Print(int Size,int *b);
    int Get_u();
    void Delete_u(int pos);
    void Clear_u();
    void Print_u();
    int complete(int beg_t, int Size, int periods, int *b);
    double bksub( int tbreak, int last_period, int Size, double slowc_l
#ifdef PROFILER
    , long int *nmul
#endif
    );
    double simple_bksub(int it_, int Size, double slowc_l);
    void run_triangular(int nop_all,int *op_all);
    void run_it(int nop_all,int *op_all);
    void run_u_period1(int periods);
    void close_swp_file();

    void read_file_table_u(t_table_u **table_u, t_table_u **F_table_u, t_table_u **i_table_u, t_table_u **F_i_table_u, int *nb_table_u, bool i_to_do, bool shifting, int *nb_add_u_count, int y_kmin, int y_kmax, int u_size);
    void read_file_table_y(t_table_y **table_y, t_table_y **i_table_y, int *nb_table_y, bool i_to_do, bool shifting, int y_kmin, int y_kmax, int u_size, int y_size);
    void close_SaveCode();

    stack<double> Stack;
    int nb_prologue_table_u, nb_first_table_u, nb_middle_table_u, nb_last_table_u;
    int nb_prologue_table_y, nb_first_table_y, nb_middle_table_y, nb_last_table_y;
    int middle_count_loop;
    char type;
    t_table_u *prologue_table_u, *first_table_u, *first_i_table_u, *middle_table_u, *middle_i_table_u, *last_table_u;
    t_table_y *prologue_table_y, *first_table_y, *middle_table_y, *middle_i_table_y, *last_table_y;
    t_table_u *F_prologue_table_u, *F_first_table_u, *F_first_i_table_u, *F_middle_table_u, *F_middle_i_table_u, *F_last_table_u;
    fstream SaveCode;
    string filename;
    int max_u, min_u;
    clock_t time00;

    Mem_Mngr mem_mngr;
    vector<int> u_liste;
    int *NbNZRow, *NbNZCol;
    NonZeroElem **FNZE_R, **FNZE_C;
    int nb_endo, u_count_init;

    int *pivot, *pivotk;
    double *pivotv, *pivotva;
    int *b;
    bool *line_done;
    bool symbolic, alt_symbolic;
    int alt_symbolic_count;
    int *g_save_op;
    int first_count_loop;
    int g_nop_all;
    int u_count_alloc, u_count_alloc_save;
    double markowitz_c_s;
    double res1a;
    long int nop_all, /*nopa_all,*/ nop1, nop2;
    map<pair<pair<int, int> ,int>, int> IM_i;
protected:
    double *u, *y, *ya;
    double res1, res2, max_res;
    double slowc, slowc_save, markowitz_c;
    int y_kmin, y_kmax, y_size, periods, y_decal;
//public:
    int  *index_vara, *index_equa;
    int u_count, tbreak_g;
    int iter;
    double *direction;

  };





#endif
