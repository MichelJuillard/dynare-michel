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

#ifndef SPARSEMATRIX_HH_INCLUDED
#define SPARSEMATRIX_HH_INCLUDED

#include <fstream>
#include <stack>
#include <cmath>
#include <map>
#include <ctime>
#include "Mem_Mngr.hh"
#ifdef _MSC_VER
  #include <limits>
#endif
#define NEW_ALLOC
#define MARKOVITZ

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
const double eps=1e-10;
const double very_big=1e24;
const int alt_symbolic_count_max=1;


class SparseMatrix
  {
  public:
    SparseMatrix();
    int simulate_NG1(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, bool cvg, int &iter, int minimal_solving_periods);
    bool simulate_NG(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, bool print_it, bool cvg, int &iter, bool steady_state);
    void Direct_Simulate(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, int iter);
    void fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1);
    void Read_SparseMatrix(string file_name, const int Size, int periods, int y_kmin, int y_kmax, bool steady_state, bool two_boundaries);
    void Read_file(string file_name, int periods, int u_size1, int y_size, int y_kmin, int y_kmax, int &nb_endo, int &u_count, int &u_count_init, double* u);

#ifdef _MSC_VER
    unsigned long nan__[2];
    double NAN;

    inline bool isnan(double value)
     {
       return value != value;
     }

    inline bool isinf(double value)
      {
        return (std::numeric_limits<double>::has_infinity &&
        value == std::numeric_limits<double>::infinity());
      }


    inline double asinh(double x)
     {
       return log(x+sqrt(x*x+1));
     }

    template<typename T>
    inline T acosh(T x)
      {
        if(!(x>=1.0)) return sqrt(-1.0);
        return log(x+sqrt(x*x-1.0));
      }

    template<typename T>
    inline T atanh(T x)
      {
        if(!(x>-1.0 && x<1.0)) return sqrt(-1.0);
        return log((1.0+x)/(1.0-x))/2.0;
      }

#endif

 private:
    void Init(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int> ,int>, int> &IM);
    void ShortInit(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int> ,int>, int> &IM);
    void Simple_Init(int it_, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> &IM);
    void End(int Size);
    bool compare( int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long int nop4,  int Size
#ifdef PROFILER
                 , long int *ndiv, long int *nsub
#endif
                );
    void Insert(const int r, const int c, const int u_index, const int lag_index);
    void Delete(const int r,const int c);
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
    void CheckIt(int y_size, int y_kmin, int y_kmax, int Size, int periods, int iter);
    void Check_the_Solution(int periods, int y_kmin, int y_kmax, int Size, double *u, int* pivot, int* b);
    int complete(int beg_t, int Size, int periods, int *b);
    double bksub( int tbreak, int last_period, int Size, double slowc_l
#ifdef PROFILER
    , long int *nmul
#endif
    );
    double simple_bksub(int it_, int Size, double slowc_l);
    /*void close_swp_file();*/
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
    map<pair<int, int>,NonZeroElem*> Mapped_Array;
    int *NbNZRow, *NbNZCol;
    NonZeroElem **FNZE_R, **FNZE_C;
    int nb_endo, u_count_init;

    int *pivot, *pivotk, *pivot_save;
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
    long int nop_all, nop1, nop2;
    map<pair<pair<int, int> ,int>, int> IM_i;
protected:
    double *u, *y, *ya;
    double res1, res2, max_res, max_res_idx;
    double slowc, slowc_save, markowitz_c;
    int y_kmin, y_kmax, y_size, periods, y_decal;
//public:
    int  *index_vara, *index_equa;
    int u_count, tbreak_g;
    int iter;
    double *direction;
    int start_compare;
    int restart;
    bool error_not_printed;
  };



#endif
