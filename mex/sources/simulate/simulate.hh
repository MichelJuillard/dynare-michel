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


#include "Interpreter.hh"
#ifndef DEBUG_EX
  #include "mex.h"
#else
  #include "mex_interface.hh"
#endif

using namespace std;




int nb_row_x, nb_row_xd, u_size, y_size, x_size, y_kmin, y_kmax, y_decal;
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
