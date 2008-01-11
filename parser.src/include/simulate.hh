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

#ifndef SIMULATE_HH_INCLUDED
#define SIMULATE_HH_INCLUDED

/*#include <stack>
#include <set>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>*/
//#include <time.h>
//#include <string>
/*#include <map>
#include <algorithm>
#include "CodeInterpreter.hh"
#include "SymbolTableTypes.hh"*/
#include "Interpreter.hh"
//#include "Mem_Mngr.hh"
/*#include "LinBCG.hh"*/
#include "mex.h"
/*#include "ExprNode.hh"*/



//#define pow pow1

/*typedef struct Variable_l
{
  int* Index;
};
typedef struct tBlock
{
    int Size, Sized, Type, Max_Lead, Max_Lag, Simulation_Type, Nb_Lead_Lag_Endo;
    int *Variable, *dVariable, *Equation;
    int *variable_dyn_index, *variable_dyn_leadlag;
    IM_compact *IM_lead_lag;
};

typedef struct tModel_Block
{
    int Size;
    tBlock * List;
};
#define MARKOVITZ
#define PRINT_OUT_p
//#define RECORD_ALL
//#define DEBUGC
//#define PRINT_OUT
//#define PRINT_OUT_y1
//#define PRINT_u
//#define PRINT_OUT_b
//#define PRINT_OUT_y
//#define DEBUG
//#define EXTENDED
//#define FLOAT
//#define WRITE_u
//#define MEMORY_LEAKS
//#define N_MX_ALLOC
//#define MEM_ALLOC_CHK
#define NEW_ALLOC
//#define PROFILER
//#ifdef EXTENDED
//typedef long double double;
//#else
//typedef double double;
//#endif
*/

using namespace std;




/*std::multimap<std::pair<int,int>, int> var_in_equ_and_lag_i, equ_in_var_and_lag_i;

int *pivota=NULL, *save_op_all=NULL;



#ifdef RECORD_ALL
bool record_all=false;
#endif
long int nopa_all;
*/
int /*Per_y_, Per_u_, it_, */nb_row_x, nb_row_xd, u_size, y_size, x_size, y_kmin, y_kmax, y_decal;
int periods, maxit_;
double *params, markowitz_c, slowc, slowc_save;
double  *u, *y, *x, *r, *g1, *g2, *ya;
double solve_tolf;
//pctimer_t t0, t1;
clock_t t0, t1;
int size_of_direction;
int i, j, k;


/*double err;






double res1, res2;

double max_res;
bool cvg;

*/


#endif // SIMULATE_HH_INCLUDED
