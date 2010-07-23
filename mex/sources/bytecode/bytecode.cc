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

////////////////////////////////////////////////////////////////////////
//                           simulate.cc                              //
//              simulate file designed for GNU GCC C++ compiler       //
////////////////////////////////////////////////////////////////////////
//#define _GLIBCXX_USE_C99_FENV_TR1 1
//#include <cfenv>


#include <cstring>
#include "Interpreter.hh"
#ifndef DEBUG_EX
# include "mex.h"
#else
# include "mex_interface.hh"
#endif

#include "Mem_Mngr.hh"


#ifdef DEBUG_EX

using namespace std;
# include <sstream>

string
Get_Argument(const char *argv)
{
  string f(argv);
  return f;
}

int
main(int argc, const char *argv[])
{
  FILE *fid;
  bool steady_state = false;
  bool evaluate = false;
  bool block = false;

  if (argc < 2)
    {
      mexPrintf("model filename expected\n");
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("Exit from Dynare");
    }
  float f_tmp;
  ostringstream tmp_out("");
  tmp_out << argv[1] << "_options.txt";

  int nb_params;
  int i, row_y, col_y, row_x, col_x;
  double *yd, *xd;
  double *direction;
  int minimal_solving_periods;

  string file_name(argv[1]);

  int count_array_argument = 0;
  for (i = 2; i < argc; i++)
    {
      if (Get_Argument(argv[i]) == "static")
        steady_state = true;
      else if (Get_Argument(argv[i]) == "dynamic")
        steady_state = false;
      else if (Get_Argument(argv[i]) == "evaluate")
        evaluate = true;
      else if (Get_Argument(argv[i]) == "block")
        block = true;
      else
        {
          mexPrintf("Unknown argument : ");
          mexEvalString("st=fclose('all');clear all;");
          string f;
          f = Get_Argument(argv[i]);
          f.append("\n");
          mexErrMsgTxt(f.c_str());
        }
    }

  fid = fopen(tmp_out.str().c_str(), "r");
  int periods = 1;
  if (!steady_state)
    fscanf(fid, "%d", &periods);
  int maxit_;
  fscanf(fid, "%d", &maxit_);
  fscanf(fid, "%f", &f_tmp);
  double slowc = f_tmp;
  fscanf(fid, "%f", &f_tmp);
  double markowitz_c = f_tmp;
  fscanf(fid, "%f", &f_tmp);
  double solve_tolf = f_tmp;
  fscanf(fid, "%d", &minimal_solving_periods);
  fclose(fid);

  tmp_out.str("");
  tmp_out << argv[1] << "_M.txt";
  fid = fopen(tmp_out.str().c_str(), "r");
  int y_kmin;
  fscanf(fid, "%d", &y_kmin);
  int y_kmax;
  fscanf(fid, "%d", &y_kmax);
  int y_decal;
  fscanf(fid, "%d", &y_decal);
  fscanf(fid, "%d", &nb_params);
  fscanf(fid, "%d", &row_x);
  fscanf(fid, "%d", &col_x);
  fscanf(fid, "%d", &row_y);
  fscanf(fid, "%d", &col_y);
  int steady_row_y, steady_col_y;
  int steady_row_x, steady_col_x;
  fscanf(fid, "%d", &steady_row_y);
  fscanf(fid, "%d", &steady_col_y);
  fscanf(fid, "%d", &steady_row_x);
  fscanf(fid, "%d", &steady_col_x);
  int nb_row_xd;
  fscanf(fid, "%d", &nb_row_xd);
  double *params = (double *) malloc(nb_params*sizeof(params[0]));
  for (i = 0; i < nb_params; i++)
    {
      fscanf(fid, "%f", &f_tmp);
      params[i] = f_tmp;
    }
  fclose(fid);

  yd = (double *) malloc(row_y*col_y*sizeof(yd[0]));
  xd = (double *) malloc(row_x*col_x*sizeof(xd[0]));
  tmp_out.str("");
  tmp_out << argv[1] << "_oo.txt";

  fid = fopen(tmp_out.str().c_str(), "r");
  for (i = 0; i < col_y*row_y; i++)
    {
      fscanf(fid, "%f", &f_tmp);
      yd[i] = f_tmp;
    }
  for (i = 0; i < col_x*row_x; i++)
    {
      fscanf(fid, "%f", &f_tmp);
      xd[i] = f_tmp;
    }
  double *steady_yd, *steady_xd;
  steady_yd = (double *) malloc(steady_row_y*steady_col_y*sizeof(steady_yd[0]));
  steady_xd = (double *) malloc(steady_row_x*steady_col_x*sizeof(steady_xd[0]));
  for (i = 0; i < steady_row_y*steady_col_y; i++)
    {
      fscanf(fid, "%f", &f_tmp);
      steady_yd[i] = f_tmp;
    }
  for (i = 0; i < steady_row_x*steady_col_x; i++)
    {
      fscanf(fid, "%f", &f_tmp);
      steady_xd[i] = f_tmp;
    }
  fclose(fid);

  int size_of_direction = col_y*row_y*sizeof(double);
  double *y = (double *) mxMalloc(size_of_direction);
  double *ya = (double *) mxMalloc(size_of_direction);
  direction = (double *) mxMalloc(size_of_direction);
  memset(direction, 0, size_of_direction);
  double *x = (double *) mxMalloc(col_x*row_x*sizeof(double));
  for (i = 0; i < row_x*col_x; i++)
    x[i] = double (xd[i]);
  for (i = 0; i < row_y*col_y; i++)
    y[i] = double (yd[i]);
  free(yd);
  free(xd);

  int y_size = row_y;
  int nb_row_x = row_x;
  clock_t t0 = clock();

  Interpreter interprete(params, y, ya, x, steady_yd, steady_xd, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal, markowitz_c, file_name, minimal_solving_periods);

  string f(file_name);
  int nb_blocks = 0;
  interprete.compute_blocks(f, f, steady_state, evaluate, block, nb_blocks);
  clock_t t1 = clock();
  if (!evaluate)
    mexPrintf("Simulation Time=%f milliseconds\n", 1000.0*(double (t1)-double (t0))/double (CLOCKS_PER_SEC));
  if (x)
    mxFree(x);
  if (y)
    mxFree(y);
  if (ya)
    mxFree(ya);
  if (direction)
    mxFree(direction);
  if (steady_yd)
    mxFree(steady_yd);
  if (steady_xd)
    mxFree(steady_xd);
  free(params);
}

#else

string
Get_Argument(const mxArray *prhs)
{
  const mxArray *mxa = prhs;
  int buflen = mxGetM(mxa) * mxGetN(mxa) + 1;
  char *first_argument;
  first_argument = (char *) mxCalloc(buflen, sizeof(char));
  int status = mxGetString(mxa, first_argument, buflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space. The first argument is truncated.");
  string f(first_argument);
  mxFree(first_argument);
  return f;
}


/* The gateway routine */
void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *M_, *oo_, *options_;
  mxArray *block_structur = NULL;
  int i, row_y = 0, col_y = 0, row_x = 0, col_x = 0, nb_row_xd = 0;
  int steady_row_y, steady_col_y, steady_row_x, steady_col_x, steady_nb_row_xd;
  int y_kmin = 0, y_kmax = 0, y_decal = 0, periods = 1;
  double *pind;
  double *direction;
  bool steady_state = false;
  bool evaluate = false;
  bool block = false;
  double *params = NULL;
  double *yd = NULL, *xd = NULL;
  int count_array_argument = 0;
  for (i = 0; i < nrhs; i++)
    {
      if (!mxIsChar(prhs[i]))
        {
          switch (count_array_argument)
            {
            case 0:
              yd = mxGetPr(prhs[i]);
              row_y = mxGetM(prhs[i]);
              col_y = mxGetN(prhs[i]);
              break;
            case 1:
              xd =  mxGetPr(prhs[i]);
              row_x = mxGetM(prhs[i]);
              col_x = mxGetN(prhs[i]);
              break;
            case 2:
              params = mxGetPr(prhs[i]);
              break;
            case 3:
              periods = mxGetScalar(prhs[i]);
              break;
            case 4:
              block_structur = mxDuplicateArray(prhs[i]);
              break;
            default:
              mexPrintf("Unknown argument\n");
            }
          count_array_argument++;
        }
      else if (Get_Argument(prhs[i]) == "static")
        steady_state = true;
      else if (Get_Argument(prhs[i]) == "dynamic")
        steady_state = false;
      else if (Get_Argument(prhs[i]) == "evaluate")
        evaluate = true;
      else if (Get_Argument(prhs[i]) == "block")
        block = true;
      else
        {
          mexPrintf("Unknown argument : ");
          string f;
          f = Get_Argument(prhs[i]);
          f.append("\n");
          mexErrMsgTxt(f.c_str());
        }
    }
  if (count_array_argument > 0 && count_array_argument < 4)
    {
      mexPrintf("Missing arguments. All the following arguments have to be indicated y, x, params, it_\n");
      mexErrMsgTxt("end of bytecode\n");
    }

  M_ = mexGetVariable("global", "M_");
  if (M_ == NULL)
    {
      mexPrintf("Global variable not found : ");
      mexErrMsgTxt("M_ \n");
    }
  /* Gets variables and parameters from global workspace of Matlab */
  oo_ = mexGetVariable("global", "oo_");
  if (oo_ == NULL)
    {
      mexPrintf("Global variable not found : ");
      mexErrMsgTxt("oo_ \n");
    }
  options_ = mexGetVariable("global", "options_");
  if (options_ == NULL)
    {
      mexPrintf("Global variable not found : ");
      mexErrMsgTxt("options_ \n");
    }
  if (!count_array_argument)
    params = mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "params")));
  double *steady_yd = NULL, *steady_xd = NULL;
  if (!steady_state)
    {
      if (!count_array_argument)
        {
          yd = mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "endo_simul")));
          row_y = mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "endo_simul")));
          col_y = mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "endo_simul")));;
        }
      if (!count_array_argument)
        {
          xd = mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_simul")));
          row_x = mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_simul")));
          col_x = mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_simul")));
        }
      nb_row_xd = int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_det_nbr"))))));

      y_kmin = int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "maximum_lag"))))));
      y_kmax = int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "maximum_lead"))))));
      y_decal = max(0, y_kmin-int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "maximum_endo_lag")))))));
      if (!count_array_argument)
        periods = int (floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "periods"))))));

      steady_yd = mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "steady_state")));
      steady_row_y = mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "steady_state")));
      steady_col_y = mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "steady_state")));;
      steady_xd = mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_steady_state")));
      steady_row_x = mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_steady_state")));
      steady_col_x = mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_steady_state")));
      steady_nb_row_xd = int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_det_nbr"))))));
    }
  else
    {
      if (!count_array_argument)
        {
          yd = mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "steady_state")));
          row_y = mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "steady_state")));
          col_y = mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "steady_state")));;
        }
      if (!count_array_argument)
        {
          xd = mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_steady_state")));
          row_x = mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_steady_state")));
          col_x = mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_, "exo_steady_state")));
        }
      nb_row_xd = int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_det_nbr"))))));
    }
  int maxit_ = int (floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "maxit_"))))));
  double slowc = double (*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "slowc")))));
  double markowitz_c = double (*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "markowitz")))));
  int minimal_solving_periods = int (*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "minimal_solving_periods")))));
  int stack_solve_algo = int (*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "stack_solve_algo")))));
  int solve_algo = int (*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "solve_algo")))));
  double solve_tolf;
  if (steady_state)
    solve_tolf = *(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "solve_tolf"))));
  else
    solve_tolf = *(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_, "dynatol"))));
  mxArray *mxa = mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "fname"));
  int buflen = mxGetM(mxa) * mxGetN(mxa) + 1;
  char *fname;
  fname = (char *) mxCalloc(buflen, sizeof(char));
  string file_name = fname;
  int status = mxGetString(mxa, fname, buflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space. Filename is truncated.");

  int size_of_direction = col_y*row_y*sizeof(double);
  double *y = (double *) mxMalloc(size_of_direction);
  double *ya = (double *) mxMalloc(size_of_direction);
  direction = (double *) mxMalloc(size_of_direction);
  memset(direction, 0, size_of_direction);
  double *x = (double *) mxMalloc(col_x*row_x*sizeof(double));
  for (i = 0; i < row_x*col_x; i++)
    x[i] = double (xd[i]);
  for (i = 0; i < row_y*col_y; i++)
    {
      y[i]  = double (yd[i]);
      ya[i] = double (yd[i]);
    }
  int y_size = row_y;
  int nb_row_x = row_x;

  clock_t t0 = clock();
  Interpreter interprete(params, y, ya, x, steady_yd, steady_xd, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal, markowitz_c, file_name, minimal_solving_periods, stack_solve_algo, solve_algo);
  string f(fname);
  int nb_blocks = 0;
  bool result = interprete.compute_blocks(f, f, steady_state, evaluate, block, nb_blocks);
  clock_t t1 = clock();
  if (!steady_state && !evaluate)
    mexPrintf("Simulation Time=%f milliseconds\n", 1000.0*(double (t1)-double (t0))/double (CLOCKS_PER_SEC));
  if (nlhs > 0)
    {
      plhs[0] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);
      pind = mxGetPr(plhs[0]);
      if (evaluate)
        for (i = 0; i < row_y*col_y; i++)
          pind[i] = y[i]-ya[i];
      else
        for (i = 0; i < row_y*col_y; i++)
          pind[i] = y[i];
      if (nlhs > 1)
        {
          plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
          pind = mxGetPr(plhs[1]);
          if (result)
            pind[0] = 0;
          else
            pind[0] = 1;
          if (nlhs > 2)
            {
              int jacob_field_number = 0, jacob_exo_field_number = 0, jacob_exo_det_field_number = 0, jacob_other_endo_field_number = 0;
              if (!block_structur)
                {
                  const char *field_names[] = {"jacob","jacob_exo","jacob_exo_det","jacob_other_endo"};
                  jacob_field_number=0;
                  jacob_exo_field_number=1;
                  jacob_exo_det_field_number=2;
                  jacob_other_endo_field_number=2;
                  mwSize dims[1] = {nb_blocks };
                  plhs[2] = mxCreateStructArray(1, dims, 4, field_names);
                  mexPrintf("the structure has been created\n");
                }
              else if (!mxIsStruct(block_structur))
                {
                  mexPrintf("The third output argument must be a structure\n");
                  mexErrMsgTxt("end of bytecode\n");
                }
              else
                {
                  mexPrintf("Adding Fields\n");
                  jacob_field_number = mxAddField(block_structur, "jacob");
                  if (jacob_field_number == -1)
                    {
                      mexPrintf("Cannot add extra field to the structArray\n");
                      mexErrMsgTxt("end of bytecode\n");
                    }
                  jacob_exo_field_number = mxAddField(block_structur, "jacob_exo");
                  if (jacob_exo_field_number == -1)
                    {
                      mexPrintf("Cannot add extra field to the structArray\n");
                      mexErrMsgTxt("end of bytecode\n");
                    }
                  jacob_exo_det_field_number = mxAddField(block_structur, "jacob_exo_det");
                  if (jacob_exo_det_field_number == -1)
                    {
                      mexPrintf("Cannot add extra field to the structArray\n");
                      mexErrMsgTxt("end of bytecode\n");
                    }
                  jacob_other_endo_field_number = mxAddField(block_structur, "jacob_other_endo");
                  if (jacob_other_endo_field_number == -1)
                    {
                      mexPrintf("Cannot add extra field to the structArray\n");
                      mexErrMsgTxt("end of bytecode\n");
                    }
                }
              for (int i = 0; i < nb_blocks; i++)
                {
                  mxSetFieldByNumber(block_structur,i,jacob_field_number,interprete.get_jacob(i));
                  mxSetFieldByNumber(block_structur,i,jacob_exo_field_number,interprete.get_jacob_exo(i));
                  mxSetFieldByNumber(block_structur,i,jacob_exo_det_field_number,interprete.get_jacob_exo_det(i));
                  mxSetFieldByNumber(block_structur,i,jacob_other_endo_field_number,interprete.get_jacob_other_endo(i));
                }
              plhs[2] = block_structur;
            }
        }
    }
  if (x)
    mxFree(x);
  if (y)
    mxFree(y);
  if (ya)
    mxFree(ya);
  if (direction)
    mxFree(direction);
}
#endif
