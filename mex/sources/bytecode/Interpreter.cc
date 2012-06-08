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

#include <cstring>
#include <sstream>
#include "Interpreter.hh"
#define BIG 1.0e+8;
#define SMALL 1.0e-5;
//#define DEBUG

Interpreter::~Interpreter()
{
}

Interpreter::Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg, double *steady_x_arg,
                         double *direction_arg, int y_size_arg,
                         int nb_row_x_arg, int nb_row_xd_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
                         int maxit_arg_, double solve_tolf_arg, int size_of_direction_arg, double slowc_arg, int y_decal_arg, double markowitz_c_arg,
                         string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
                         bool global_temporary_terms_arg, bool print_arg, bool print_error_arg, mxArray *GlobalTemporaryTerms_arg)
{
  params = params_arg;
  y = y_arg;
  ya = ya_arg;
  x = x_arg;
  steady_y = steady_y_arg;
  steady_x = steady_x_arg;
  direction = direction_arg;
  y_size = y_size_arg;
  nb_row_x = nb_row_x_arg;
  nb_row_xd = nb_row_xd_arg;
  periods = periods_arg;
  y_kmax = y_kmax_arg;
  y_kmin = y_kmin_arg;
  maxit_ = maxit_arg_;
  solve_tolf = solve_tolf_arg;
  size_of_direction = size_of_direction_arg;
  slowc = slowc_arg;
  slowc_save = slowc;
  y_decal = y_decal_arg;
  markowitz_c = markowitz_c_arg;
  filename = filename_arg;
  T = NULL;
  error_not_printed = true;
  minimal_solving_periods = minimal_solving_periods_arg;
  stack_solve_algo = stack_solve_algo_arg;
  solve_algo = solve_algo_arg;
  global_temporary_terms = global_temporary_terms_arg;
  print = print_arg;
  GlobalTemporaryTerms = GlobalTemporaryTerms_arg;
  print_error = print_error_arg;
}

double
Interpreter::pow1(double a, double b)
{
  double r = pow_(a, b);
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = 0.0000000000000000000000001;
      if (print_error)
        throw PowExceptionHandling(a, b);
    }
  return r;
}

double
Interpreter::divide(double a, double b)
{
  double r = a / b;
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = 1e70;
      if (print_error)
        throw DivideExceptionHandling(a, b);
    }
  return r;
}

double
Interpreter::log1(double a)
{
  double r = log(a);
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = -1e70;
      if (print_error)
        throw LogExceptionHandling(a);
    }
  return r;
}

double
Interpreter::log10_1(double a)
{
  double r = log(a);
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = -1e70;
      if (print_error)
        throw Log10ExceptionHandling(a);
    }
  return r;
}

void
Interpreter::compute_block_time(int Per_u_, bool evaluate, int block_num, int size, bool steady_state)
{
  int var = 0, lag = 0, op;
  unsigned int eq, pos_col;
  ostringstream tmp_out;
  double v1, v2, v3;
  bool go_on = true;
  double ll;
  double rr;
  double *jacob = NULL, *jacob_other_endo = NULL, *jacob_exo = NULL, *jacob_exo_det = NULL;
  EQN_block = block_num;
  stack<double> Stack;
  external_function_type function_type = ExternalFunctionWithoutDerivative;

#ifdef DEBUG
  mexPrintf("compute_block_time\n");
#endif
  if (evaluate /*&& !steady_state*/)
    {
      jacob = mxGetPr(jacobian_block[block_num]);
      if (!steady_state)
        {
          jacob_other_endo = mxGetPr(jacobian_other_endo_block[block_num]);
          jacob_exo = mxGetPr(jacobian_exo_block[block_num]);
          jacob_exo_det = mxGetPr(jacobian_det_exo_block[block_num]);
        }
    }

  while (go_on)
    {
      switch (it_code->first)
        {
        case FNUMEXPR:
#ifdef DEBUG
          mexPrintf("FNUMEXPR\n");
#endif
          it_code_expr = it_code;
          switch (((FNUMEXPR_ *) it_code->second)->get_expression_type())
            {
            case TemporaryTerm:
#ifdef DEBUG
              mexPrintf("TemporaryTerm\n");
#endif
              EQN_type = TemporaryTerm;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
#ifdef DEBUG
              mexPrintf("EQN_equation=%d\n", EQN_equation); mexEvalString("drawnow;");
#endif
              break;
            case ModelEquation:
#ifdef DEBUG
              mexPrintf("ModelEquation\n");
#endif
              EQN_type = ModelEquation;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              break;
            case FirstEndoDerivative:
#ifdef DEBUG
              mexPrintf("FirstEndoDerivative\n");
#endif
              EQN_type = FirstEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstOtherEndoDerivative:
#ifdef DEBUG
              mexPrintf("FirstOtherEndoDerivative\n");
#endif
              EQN_type = FirstOtherEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExoDerivative:
#ifdef DEBUG
              mexPrintf("FirstExoDerivative\n");
#endif
              EQN_type = FirstExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExodetDerivative:
#ifdef DEBUG
              mexPrintf("FirstExodetDerivative\n");
#endif
              EQN_type = FirstExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstParamDerivative:
#ifdef DEBUG
              mexPrintf("FirstParamDerivative\n");
#endif
              EQN_type = FirstParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              break;
            case SecondEndoDerivative:
#ifdef DEBUG
              mexPrintf("SecondEndoDerivative\n");
#endif
              EQN_type = SecondEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExoDerivative:
#ifdef DEBUG
              mexPrintf("SecondExoDerivative\n");
#endif
              EQN_type = SecondExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExodetDerivative:
#ifdef DEBUG
              mexPrintf("SecondExodetDerivative\n");
#endif
              EQN_type = SecondExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondParamDerivative:
#ifdef DEBUG
              mexPrintf("SecondParamDerivative\n");
#endif
              EQN_type = SecondParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              break;
            case ThirdEndoDerivative:
#ifdef DEBUG
              mexPrintf("ThirdEndoDerivative\n");
#endif
              EQN_type = ThirdEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              EQN_lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdExoDerivative:
#ifdef DEBUG
              mexPrintf("ThirdExoDerivative\n");
#endif
              EQN_type = ThirdExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              EQN_lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdExodetDerivative:
#ifdef DEBUG
              mexPrintf("ThirdExodetDerivative\n");
#endif
              EQN_type = ThirdExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              EQN_lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdParamDerivative:
#ifdef DEBUG
              mexPrintf("ThirdParamDerivative\n");
#endif
              EQN_type = ThirdParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              break;
            }
          break;
        case FLDV:
          //load a variable in the processor
          switch (((FLDV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDV Param[var=%d]\n", var);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV y[var=%d, lag=%d, it_=%d], y_size=%d evaluate=%d\n", var, lag, it_, y_size, evaluate);
#endif
              if (evaluate)
                Stack.push(ya[(it_+lag)*y_size+var]);
              else
                Stack.push(y[(it_+lag)*y_size+var]);
#ifdef DEBUG
              tmp_out << " y[" << it_+lag << ", " << var << "](" << y[(it_+lag)*y_size+var] << ")";
#endif
              break;
            case eExogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV x[var=%d, lag=%d, it_=%d], nb_row_x=%d evaluate=%d\n", var, lag, it_, nb_row_x, evaluate);
              tmp_out << " x[" << it_+lag << ", " << var << "](" << x[it_+lag+var*nb_row_x] << ")";
#endif
              Stack.push(x[it_+lag+var*nb_row_x]);
              break;
            case eExogenousDet:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
              Stack.push(x[it_+lag+var*nb_row_xd]);
              break;
            case eModelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDV a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDV: Unknown variable type\n");
            }
          break;
        case FLDSV:
          //load a variable in the processor
          switch (((FLDSV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV Param[var=%d]=%f\n", var, params[var]);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV y[var=%d]=%f\n", var, ya[var]);
              tmp_out << " y[" << var << "](" << y[var] << ")";
#endif
              if (evaluate)
                Stack.push(ya[var]);
              else
                Stack.push(y[var]);
              break;
            case eExogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV x[var=%d]\n", var);
              tmp_out << " x[" << var << "](" << x[var] << ")";
#endif
              Stack.push(x[var]);
              break;
            case eExogenousDet:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV xd[var=%d]\n", var);
#endif
              Stack.push(x[var]);
              break;
            case eModelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDSV a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDSV: Unknown variable type\n");
            }
          break;
        case FLDVS:
          //load a variable in the processor
          switch (((FLDVS_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("params[%d]\n", var);
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS steady_y[%d]\n", var);
#endif
              Stack.push(steady_y[var]);
              break;
            case eExogenous:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS x[%d] \n", var);
#endif
              Stack.push(x[var]);
              break;
            case eExogenousDet:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS xd[%d]\n", var);
#endif
              Stack.push(x[var]);
              break;
            case eModelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDVS a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDVS: Unknown variable type\n");
            }
          break;
        case FLDT:
          //load a temporary variable in the processor
          var = ((FLDT_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("T[it_=%d var=%d, y_kmin=%d, y_kmax=%d == %d]=>%f\n", it_, var, y_kmin, y_kmax, var*(periods+y_kmin+y_kmax)+it_, var);
          tmp_out << " T[" << it_ << ", " << var << "](" << T[var*(periods+y_kmin+y_kmax)+it_] << ")";
#endif
          Stack.push(T[var*(periods+y_kmin+y_kmax)+it_]);
          break;
        case FLDST:
          //load a temporary variable in the processor
          var = ((FLDST_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDST T[%d]", var);
#endif
          Stack.push(T[var]);
#ifdef DEBUG
          mexPrintf("=%f\n", T[var]);
          tmp_out << " T[" << var << "](" << T[var] << ")";
#endif
          break;
        case FLDU:
          //load u variable in the processor
          var = ((FLDU_ *) it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FLDU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDSU:
          //load u variable in the processor
          var = ((FLDSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDSU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDR:
          //load u variable in the processor
          var = ((FLDR_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDR r[%d]\n", var);
#endif
          Stack.push(r[var]);
          break;
        case FLDZ:
          //load 0 in the processor
#ifdef DEBUG
          mexPrintf("FLDZ\n");
#endif
          Stack.push(0.0);
#ifdef DEBUG
          tmp_out << " 0";
#endif
          break;
        case FLDC:
          //load a numerical constant in the processor
          ll = ((FLDC_ *) it_code->second)->get_value();
#ifdef DEBUG
          mexPrintf("FLDC = %f\n", ll);
          tmp_out << " " << ll;
#endif

          Stack.push(ll);
          break;
        case FSTPV:
          //load a variable in the processor
          switch (((FSTPV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FSTPV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FSTPV params[%d]\n", var);
#endif
              params[var] = Stack.top();
              Stack.pop();
              break;
            case eEndogenous:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              y[(it_+lag)*y_size+var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d, %d](%f)=%s\n", it_+lag, var, y[(it_+lag)*y_size+var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case eExogenous:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              x[it_+lag+var*nb_row_x]  = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[it_+lag+var*nb_row_x], tmp_out.str().c_str());
              tmp_out.str("");
#endif

              Stack.pop();
              break;
            case eExogenousDet:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              x[it_+lag+var*nb_row_xd] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[it_+lag+var*nb_row_xd], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            default:
              mexPrintf("FSTPV: Unknown variable type\n");
            }
          break;
        case FSTPSV:
          //load a variable in the processor
          switch (((FSTPSV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              params[var] = Stack.top();
              Stack.pop();
              break;
            case eEndogenous:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              y[var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d](%f)=%s\n", var, y[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case eExogenous:
            case eExogenousDet:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              x[var]  = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            default:
              mexPrintf("FSTPSV: Unknown variable type\n");
            }
          break;
        case FSTPT:
          //store in a temporary variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPT\n");
#endif
          var = ((FSTPT_ *) it_code->second)->get_pos();
          T[var*(periods+y_kmin+y_kmax)+it_] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" T[%d, %d](%f)=%s\n", it_, var, T[var*(periods+y_kmin+y_kmax)+it_], tmp_out.str().c_str());
          tmp_out.str("");
#endif

          Stack.pop();
          break;
        case FSTPST:
          //store in a temporary variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPST\n");
#endif
          var = ((FSTPST_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("var=%d\n", var);
#endif
          T[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" T[%d](%f)=%s\n", var, T[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPU:
          //store in u variable from the processor
          var = ((FSTPU_ *) it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FSTPU\n");
          mexPrintf("var=%d\n", var);
#endif
          u[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPSU:
          //store in u variable from the processor
          var = ((FSTPSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
          if (var >= u_count_alloc || var < 0)
            mexPrintf("Erreur var=%d\n", var);
#endif
          u[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPR:
          //store in residual variable from the processor
          var = ((FSTPR_ *) it_code->second)->get_pos();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf("FSTPR r[%d]", var);
          tmp_out.str("");
#endif
          r[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf("(%f)=%s\n", r[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPG:
          //store in derivative (g) variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPG\n");
          mexEvalString("drawnow;");
#endif
          var = ((FSTPG_ *) it_code->second)->get_pos();
          g1[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;

        case FSTPG2:
          //store in the jacobian matrix
          rr = Stack.top();
          if (EQN_type != FirstEndoDerivative)
            {
              ostringstream tmp;
              tmp << " in compute_block_time, impossible case " << EQN_type << " not implement in static jacobian\n";
              throw FatalExceptionHandling(tmp.str());
            }
          eq = ((FSTPG2_ *) it_code->second)->get_row();
          var = ((FSTPG2_ *) it_code->second)->get_col();
#ifdef DEBUG
          mexPrintf("FSTPG2 eq=%d, var=%d\n", eq, var);
          mexEvalString("drawnow;");
#endif
          jacob[eq + size*var] = rr;
          break;
        case FSTPG3:
          //store in derivative (g) variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPG3\n");
          mexEvalString("drawnow;");
#endif
          rr = Stack.top();
          switch (EQN_type)
            {
            case FirstEndoDerivative:
              eq = ((FSTPG3_ *) it_code->second)->get_row();
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Endo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
#endif
              jacob[eq + size*pos_col] = rr;
              break;
            case FirstOtherEndoDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
              jacob_other_endo[eq + size*pos_col] = rr;
              break;
            case FirstExoDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Exo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
#endif
              jacob_exo[eq + size*pos_col] = rr;
              break;
            case FirstExodetDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
              jacob_exo_det[eq + size*pos_col] = rr;
              break;
            default:
              ostringstream tmp;
              tmp << " in compute_block_time, variable " << EQN_type << " not used yet\n";
              throw FatalExceptionHandling(tmp.str());
            }
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;

        case FBINARY:
          op = ((FBINARY_ *) it_code->second)->get_op_type();
#ifdef DEBUG
          mexPrintf("FBINARY, op=%d\n", op);
#endif
          v2 = Stack.top();
          Stack.pop();
          v1 = Stack.top();
          Stack.pop();
          switch (op)
            {
            case oPlus:
              Stack.push(v1 + v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "+" << v2 << "|";
#endif
              break;
            case oMinus:
              Stack.push(v1 - v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "-" << v2 << "|";
#endif
              break;
            case oTimes:
              Stack.push(v1 * v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "*" << v2 << "|";
#endif
              break;
            case oDivide:
              double tmp;
#ifdef DEBUG
              mexPrintf("v1=%f / v2=%f\n", v1, v2);
#endif
              try
                {
                  tmp = divide(v1, v2);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
#ifdef DEBUG
              tmp_out << " |" << v1 << "/" << v2 << "|";
#endif
              break;
            case oLess:
              Stack.push(double (v1 < v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "<" << v2 << "|";
#endif
              break;
            case oGreater:
              Stack.push(double (v1 > v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << ">" << v2 << "|";
#endif
              break;
            case oLessEqual:
              Stack.push(double (v1 <= v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "<=" << v2 << "|";
#endif
              break;
            case oGreaterEqual:
              Stack.push(double (v1 >= v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << ">=" << v2 << "|";
#endif
              break;
            case oEqualEqual:
              Stack.push(double (v1 == v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "==" << v2 << "|";
#endif
              break;
            case oDifferent:
              Stack.push(double (v1 != v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "!=" << v2 << "|";
#endif
              break;
            case oPower:
#ifdef DEBUG
              mexPrintf("pow\n");
#endif
              try
                {
                  tmp = pow1(v1, v2);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);

#ifdef DEBUG
              tmp_out << " |" << v1 << "^" << v2 << "|";
#endif
              break;
            case oPowerDeriv:
              {
                int derivOrder = nearbyint(Stack.top());
                Stack.pop();
                try
                  {
                    if (fabs(v1) < NEAR_ZERO && v2 > 0
                        && derivOrder > v2
                        && fabs(v2-nearbyint(v2)) < NEAR_ZERO)
                      Stack.push(0.0);
                    else
                      {
                        double dxp = pow1(v1, v2-derivOrder);
                        for (int i = 0; i < derivOrder; i++)
                          dxp *= v2--;
                        Stack.push(dxp);
                      }
                  }
                catch (FloatingPointExceptionHandling &fpeh)
                  {
                    mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                    go_on = false;
                  }
              }

#ifdef DEBUG
              tmp_out << " |PowerDeriv(" << v1 << ", " << v2 << ")|";
#endif
              break;
            case oMax:
              Stack.push(max(v1, v2));
#ifdef DEBUG
              tmp_out << " |max(" << v1 << "," << v2 << ")|";
#endif
              break;
            case oMin:
              Stack.push(min(v1, v2));
#ifdef DEBUG
              tmp_out << " |min(" << v1 << "," << v2 << ")|";
#endif
              break;
            case oEqual:
              // Nothing to do
              break;
            default:
              {
                mexPrintf("Error\n");
                ostringstream tmp;
                tmp << " in compute_block_time, unknown binary operator " << op << "\n";
                throw FatalExceptionHandling(tmp.str());
              }
            }
          break;
        case FUNARY:
          op = ((FUNARY_ *) it_code->second)->get_op_type();
          v1 = Stack.top();
          Stack.pop();
#ifdef DEBUG
          mexPrintf("FUNARY, op=%d\n", op);
#endif
          switch (op)
            {
            case oUminus:
              Stack.push(-v1);
#ifdef DEBUG
              tmp_out << " |-(" << v1 << ")|";
#endif

              break;
            case oExp:
              Stack.push(exp(v1));
#ifdef DEBUG
              tmp_out << " |exp(" << v1 << ")|";
#endif
              break;
            case oLog:
              double tmp;
              try
                {
                  tmp = log1(v1);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
              //if (isnan(res1))

#ifdef DEBUG
              tmp_out << " |log(" << v1 << ")|";
#endif
              break;
            case oLog10:
              try
                {
                  tmp = log10_1(v1);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
#ifdef DEBUG
              tmp_out << " |log10(" << v1 << ")|";
#endif
              break;
            case oCos:
              Stack.push(cos(v1));
#ifdef DEBUG
              tmp_out << " |cos(" << v1 << ")|";
#endif
              break;
            case oSin:
              Stack.push(sin(v1));
#ifdef DEBUG
              tmp_out << " |sin(" << v1 << ")|";
#endif
              break;
            case oTan:
              Stack.push(tan(v1));
#ifdef DEBUG
              tmp_out << " |tan(" << v1 << ")|";
#endif
              break;
            case oAcos:
              Stack.push(acos(v1));
#ifdef DEBUG
              tmp_out << " |acos(" << v1 << ")|";
#endif
              break;
            case oAsin:
              Stack.push(asin(v1));
#ifdef DEBUG
              tmp_out << " |asin(" << v1 << ")|";
#endif
              break;
            case oAtan:
              Stack.push(atan(v1));
#ifdef DEBUG
              tmp_out << " |atan(" << v1 << ")|";
#endif
              break;
            case oCosh:
              Stack.push(cosh(v1));
#ifdef DEBUG
              tmp_out << " |cosh(" << v1 << ")|";
#endif
              break;
            case oSinh:
              Stack.push(sinh(v1));
#ifdef DEBUG
              tmp_out << " |sinh(" << v1 << ")|";
#endif
              break;
            case oTanh:
              Stack.push(tanh(v1));
#ifdef DEBUG
              tmp_out << " |tanh(" << v1 << ")|";
#endif
              break;
            case oAcosh:
              Stack.push(acosh(v1));
#ifdef DEBUG
              tmp_out << " |acosh(" << v1 << ")|";
#endif
              break;
            case oAsinh:
              Stack.push(asinh(v1));
#ifdef DEBUG
              tmp_out << " |asinh(" << v1 << ")|";
#endif
              break;
            case oAtanh:
              Stack.push(atanh(v1));
#ifdef DEBUG
              tmp_out << " |atanh(" << v1 << ")|";
#endif
              break;
            case oSqrt:
              Stack.push(sqrt(v1));
#ifdef DEBUG
              tmp_out << " |sqrt(" << v1 << ")|";
#endif
              break;
            case oErf:
              Stack.push(erf(v1));
#ifdef DEBUG
              tmp_out << " |erf(" << v1 << ")|";

#endif
              break;
            default:
              {
                mexPrintf("Error\n");
                ostringstream tmp;
                tmp << " in compute_block_time, unknown unary operator " << op << "\n";
                throw FatalExceptionHandling(tmp.str());
              }
            }
          break;
        case FTRINARY:
          op = ((FTRINARY_ *) it_code->second)->get_op_type();
          v3 = Stack.top();
          Stack.pop();
          v2 = Stack.top();
          Stack.pop();
          v1 = Stack.top();
          Stack.pop();
          switch (op)
            {
            case oNormcdf:
              Stack.push(0.5*(1+erf((v1-v2)/v3/M_SQRT2)));
#ifdef DEBUG
              tmp_out << " |normcdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            case oNormpdf:
              Stack.push(1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3, 2)/2)));
#ifdef DEBUG
              tmp_out << " |normpdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            default:
              {
                mexPrintf("Error\n");
                ostringstream tmp;
                tmp << " in compute_block_time, unknown trinary operator " << op << "\n";
                throw FatalExceptionHandling(tmp.str());
              }
            }
          break;

        case FPUSH:
          break;

        case FCALL:
          {
#ifdef DEBUG
            mexPrintf("------------------------------\n");
            mexPrintf("CALL "); mexEvalString("drawnow;");
#endif
            FCALL_ *fc = (FCALL_ *) it_code->second;
            string function_name = fc->get_function_name();
#ifdef DEBUG
            mexPrintf("function_name=%s ", function_name.c_str()); mexEvalString("drawnow;");
#endif
            unsigned int nb_input_arguments = fc->get_nb_input_arguments();
#ifdef DEBUG
            mexPrintf("nb_input_arguments=%d ", nb_input_arguments); mexEvalString("drawnow;");
#endif
            unsigned int nb_output_arguments = fc->get_nb_output_arguments();
#ifdef DEBUG
            mexPrintf("nb_output_arguments=%d\n", nb_output_arguments); mexEvalString("drawnow;");
#endif

            mxArray *output_arguments[3];
            string arg_func_name = fc->get_arg_func_name();
#ifdef DEBUG
            mexPrintf("arg_func_name.length() = %d\n", arg_func_name.length());
            mexPrintf("arg_func_name.c_str() = %s\n", arg_func_name.c_str());
#endif
            unsigned int nb_add_input_arguments = fc->get_nb_add_input_arguments();
            function_type = fc->get_function_type();
#ifdef DEBUG
            mexPrintf("function_type=%d ExternalFunctionWithoutDerivative=%d\n", function_type, ExternalFunctionWithoutDerivative);
            mexEvalString("drawnow;");
#endif
            mxArray **input_arguments;
            switch (function_type)
              {
              case ExternalFunctionWithoutDerivative:
              case ExternalFunctionWithFirstDerivative:
              case ExternalFunctionWithFirstandSecondDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
#ifdef DEBUG
                  mexPrintf("Stack.size()=%d\n", Stack.size());
                  mexEvalString("drawnow;");
#endif
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[nb_input_arguments - i - 1] = vv;
                      Stack.pop();
                    }
                  mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                  if (function_type == ExternalFunctionWithFirstDerivative || function_type == ExternalFunctionWithFirstandSecondDerivative)
                    {
                      unsigned int indx = fc->get_indx();
                      double *FD1 = mxGetPr(output_arguments[1]);
                      unsigned int rows = mxGetN(output_arguments[1]);
                      for (unsigned int i = 0; i < rows; i++)
                        TEFD[make_pair(indx, i)] = FD1[i];
                    }
                  if (function_type == ExternalFunctionWithFirstandSecondDerivative)
                    {
                      unsigned int indx = fc->get_indx();
                      double *FD2 = mxGetPr(output_arguments[2]);
                      unsigned int rows = mxGetM(output_arguments[2]);
                      unsigned int cols = mxGetN(output_arguments[2]);
                      unsigned int k = 0;
                      for (unsigned int j = 0; j < cols; j++)
                        for (unsigned int i = 0; i < rows; i++)
                          TEFDD[make_pair(indx, make_pair(i, j))] = FD2[k++];
                    }
                }
                break;
              case ExternalFunctionNumericalFirstDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                  mxArray *vv = mxCreateString(arg_func_name.c_str());
                  input_arguments[0] = vv;
                  vv = mxCreateDoubleScalar(fc->get_row());
                  input_arguments[1] = vv;
                  vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                  for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                    {
                      double rr = Stack.top();
#ifdef DEBUG
                      mexPrintf("i=%d rr = %f Stack.size()=%d\n", i, rr, Stack.size());
#endif
                      mxSetCell(vv, nb_add_input_arguments - (i+1), mxCreateDoubleScalar(rr));
                      Stack.pop();
                    }
                  input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
#ifdef DEBUG
                  mexCallMATLAB(0, NULL, 1, &input_arguments[0], "disp");
                  mexCallMATLAB(0, NULL, 1, &input_arguments[1], "disp");
                  mexCallMATLAB(0, NULL, 1, &input_arguments[2], "celldisp");
                  mexPrintf("OK\n");
                  mexEvalString("drawnow;");
#endif
                  nb_input_arguments = 3;
                  mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                  double *rr = mxGetPr(output_arguments[0]);
#ifdef DEBUG
                  mexPrintf("*rr=%f\n", *rr);
#endif
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionFirstDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[(nb_input_arguments - 1) - i] = vv;
                      Stack.pop();
                    }
                  mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                  unsigned int indx = fc->get_indx();
                  double *FD1 = mxGetPr(output_arguments[0]);
                  //mexPrint
                  unsigned int rows = mxGetN(output_arguments[0]);
                  for (unsigned int i = 0; i < rows; i++)
                    TEFD[make_pair(indx, i)] = FD1[i];
                }
                break;
              case ExternalFunctionNumericalSecondDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                  mxArray *vv = mxCreateString(arg_func_name.c_str());
                  input_arguments[0] = vv;
                  vv = mxCreateDoubleScalar(fc->get_row());
                  input_arguments[1] = vv;
                  vv = mxCreateDoubleScalar(fc->get_col());
                  input_arguments[2] = vv;
                  vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                  for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                    {
                      double rr = Stack.top();
#ifdef DEBUG
                      mexPrintf("i=%d rr = %f\n", i, rr);
#endif
                      mxSetCell(vv, (nb_add_input_arguments - 1) - i, mxCreateDoubleScalar(rr));
                      Stack.pop();
                    }
                  input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
#ifdef DEBUG
                  mexCallMATLAB(0, NULL, 1, &input_arguments[0], "disp");
                  mexCallMATLAB(0, NULL, 1, &input_arguments[1], "disp");
                  mexCallMATLAB(0, NULL, 1, &input_arguments[2], "celldisp");
                  mexPrintf("OK\n");
                  mexEvalString("drawnow;");
#endif
                  nb_input_arguments = 3;
                  mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionSecondDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[i] = vv;
                      Stack.pop();
                    }
                  mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                  unsigned int indx = fc->get_indx();
                  double *FD2 = mxGetPr(output_arguments[2]);
                  unsigned int rows = mxGetM(output_arguments[0]);
                  unsigned int cols = mxGetN(output_arguments[0]);
                  unsigned int k = 0;
                  for (unsigned int j = 0; j < cols; j++)
                    for (unsigned int i = 0; i < rows; i++)
                      TEFDD[make_pair(indx, make_pair(i, j))] = FD2[k++];
                }
                break;
              }
          }
          break;
        case FSTPTEF:
          var = ((FSTPTEF_ *) it_code->second)->get_number();
#ifdef DEBUG
          mexPrintf("FSTPTEF\n");
          mexPrintf("var=%d Stack.size()=%d\n", var, Stack.size());
#endif
          TEF[var-1] = Stack.top();
#ifdef DEBUG
          mexPrintf("FSTP TEF[var-1]=%f done\n", TEF[var-1]);
          mexEvalString("drawnow;");
#endif
          Stack.pop();
          break;
        case FLDTEF:
          var = ((FLDTEF_ *) it_code->second)->get_number();
#ifdef DEBUG
          mexPrintf("FLDTEF\n");
          mexPrintf("var=%d Stack.size()=%d\n", var, Stack.size());
          mexPrintf("FLD TEF[var-1]=%f done\n", TEF[var-1]);
          mexEvalString("drawnow;");
#endif
          Stack.push(TEF[var-1]);
          break;
        case FSTPTEFD:
          {
            unsigned int indx = ((FSTPTEFD_ *) it_code->second)->get_indx();
            unsigned int row = ((FSTPTEFD_ *) it_code->second)->get_row();
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (function_type == ExternalFunctionNumericalFirstDerivative)
              {
                TEFD[make_pair(indx, row-1)] = Stack.top();
#ifdef DEBUG
                mexPrintf("FSTP TEFD[make_pair(indx, row)]=%f done\n", TEFD[make_pair(indx, row-1)]);
                mexEvalString("drawnow;");
#endif
                Stack.pop();
              }
          }

          break;
        case FLDTEFD:
          {
            unsigned int indx = ((FLDTEFD_ *) it_code->second)->get_indx();
            unsigned int row = ((FLDTEFD_ *) it_code->second)->get_row();
#ifdef DEBUG
            mexPrintf("FLDTEFD\n");
            mexPrintf("indx=%d row=%d Stack.size()=%d\n", indx, row, Stack.size());
            mexPrintf("FLD TEFD[make_pair(indx, row)]=%f done\n", TEFD[make_pair(indx, row-1)]);
            mexEvalString("drawnow;");
#endif
            Stack.push(TEFD[make_pair(indx, row-1)]);
          }
          break;
        case FSTPTEFDD:
          {
            unsigned int indx = ((FSTPTEFDD_ *) it_code->second)->get_indx();
            unsigned int row = ((FSTPTEFDD_ *) it_code->second)->get_row();
            unsigned int col = ((FSTPTEFDD_ *) it_code->second)->get_col();
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (function_type == ExternalFunctionNumericalSecondDerivative)
              {
                TEFDD[make_pair(indx, make_pair(row-1, col-1))] = Stack.top();
#ifdef DEBUG
                mexPrintf("FSTP TEFDD[make_pair(indx, make_pair(row, col))]=%f done\n", TEFDD[make_pair(indx, make_pair(row, col))]);
                mexEvalString("drawnow;");
#endif
                Stack.pop();
              }
          }

          break;
        case FLDTEFDD:
          {
            unsigned int indx = ((FLDTEFDD_ *) it_code->second)->get_indx();
            unsigned int row = ((FLDTEFDD_ *) it_code->second)->get_row();
            unsigned int col = ((FSTPTEFDD_ *) it_code->second)->get_col();
#ifdef DEBUG
            mexPrintf("FLDTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
            mexPrintf("FLD TEFD[make_pair(indx, make_pair(row, col))]=%f done\n", TEFDD[make_pair(indx, make_pair(row, col))]);
            mexEvalString("drawnow;");
#endif
            Stack.push(TEFDD[make_pair(indx, make_pair(row-1, col-1))]);
          }
          break;
        case FCUML:
          v1 = Stack.top();
          Stack.pop();
          v2 = Stack.top();
          Stack.pop();
          Stack.push(v1+v2);
          break;
        case FENDBLOCK:
          //it's the block end
#ifdef DEBUG
          mexPrintf("FENDBLOCK\n");
#endif
          go_on = false;
          break;
        case FENDEQU:
          break;
        case FJMPIFEVAL:
          if (evaluate)
            {
#ifdef DEBUG
              mexPrintf("FJMPIFEVAL length=%d\n", ((FJMPIFEVAL_ *) it_code->second)->get_pos());
              mexEvalString("drawnow;");
#endif
              it_code += ((FJMPIFEVAL_ *) it_code->second)->get_pos() /* - 1*/;
            }
          break;
        case FJMP:
#ifdef DEBUG
          mexPrintf("FJMP length=%d\n", ((FJMP_ *) it_code->second)->get_pos());
          mexEvalString("drawnow;");
#endif
          it_code += ((FJMP_ *) it_code->second)->get_pos() /*- 1 */;
          break;
        case FOK:
          op = ((FOK_ *) it_code->second)->get_arg();
          if (Stack.size() > 0)
            {
              ostringstream tmp;
              tmp << " in compute_block_time, stack not empty\n";
              throw FatalExceptionHandling(tmp.str());
            }
          break;
        default:
          ostringstream tmp;
          tmp << " in compute_block_time, unknown opcode " << it_code->first << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
      it_code++;
    }
#ifdef DEBUG
  mexPrintf("==> end of compute_block_time Block = %d\n", block_num);
  mexEvalString("drawnow;");
#endif
}

void
Interpreter::evaluate_a_block(const int size, const int type, string bin_basename, bool steady_state, int block_num,
                              const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag,
                              const int Block_List_Max_Lead, const int u_count_int, int block)
{
  it_code_type begining;

  switch (type)
    {
    case EVALUATE_FORWARD:
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          if (block >= 0)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              {
                //mexPrintf("=>residual[Block_Contain[%d].Equation = %d]=%g (y = %g, ya = %g)\n", j, Block_Contain[j].Equation, y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable], y[Block_Contain[j].Variable], ya[Block_Contain[j].Variable]);
                residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
              }
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              if (block >= 0)
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
            }
        }
      break;
    case SOLVE_FORWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          if (block < 0)
            {
              for (int j = 0; j < size; j++)
                {
                  //mexPrintf("residual[Block_Contain[%d].Equation = %d]=%g\n", j, Block_Contain[j].Equation, r[j]);
                  residual[Block_Contain[j].Equation] = r[j];
                }
            }
          else
            {
              for (int j = 0; j < size; j++)
                residual[j] = r[j];
            }
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              if (block < 0)
                {
                  for (int j = 0; j < size; j++)
                    {
                      //mexPrintf("residual[Per_y + Block_Contain[%d].Equation = %d]=%g\n", j, Per_y_ + Block_Contain[j].Equation, r[j]);
                      residual[Per_y_+Block_Contain[j].Equation] = r[j];
                    }
                }
              else
                {
                  for (int j = 0; j < size; j++)
                    residual[it_*size+j] = r[j];
                }
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_FORWARD_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false, stack_solve_algo, solve_algo);
#ifdef DEBUG
      mexPrintf("in SOLVE_FORWARD_COMPLETE r = mxMalloc(%d*sizeof(double))\n", size);
#endif
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          if (block < 0)
            for (int j = 0; j < size; j++)
              {
                //mexPrintf("residual[Block_Contain[%d].Equation = %d]=%g\n", j, Block_Contain[j].Equation, r[j]);
                residual[Block_Contain[j].Equation] = r[j];
              }
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  {
                    //mexPrintf("residual[it_*y_size+Block_Contain[%d].Equation = %d]=%g\n", j, it_*y_size+Block_Contain[j].Equation, r[j]);
                    residual[it_*y_size+Block_Contain[j].Equation] = r[j];
                  }
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(r);
      break;
    case EVALUATE_BACKWARD:
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          if (block >= 0)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              {
                //mexPrintf("residual[Block_Contain[%d].Equation = %d]=%g\n", j, Block_Contain[j].Equation, y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable]);
                residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
              }
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              if (block >= 0)
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
            }
        }
      break;
    case SOLVE_BACKWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          if (block < 0)
            {
              for (int j = 0; j < size; j++)
                {
                  //mexPrintf("residual[Block_Contain[%d].Equation = %d]=%g\n", j, Block_Contain[j].Equation, r[j]);
                  residual[Block_Contain[j].Equation] = r[j];
                }
            }
          else
            {
              for (int j = 0; j < size; j++)
                residual[j] = r[j];
            }
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              if (block < 0)
                {
                  for (int j = 0; j < size; j++)
                    {
                      //mexPrintf("residual[Per_y_+Block_Contain[%d].Equation = %d]=%g\n", j, Per_y_+Block_Contain[j].Equation, r[j]);
                      residual[Per_y_+Block_Contain[j].Equation] = r[j];
                    }
                }
              else
                {
                  for (int j = 0; j < size; j++)
                    residual[it_*size+j] = r[j];
                }
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_BACKWARD_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false, stack_solve_algo, solve_algo);
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          if (block < 0)
            for (int j = 0; j < size; j++)
              {
                //mexPrintf("residual[Block_Contain[%d].Equation = %d]=%g\n", j, Block_Contain[j].Equation, r[j]);
                residual[Block_Contain[j].Equation] = r[j];
              }
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  {
                    //mexPrintf("residual[Per_y_+Block_Contain[%d].Equation = %d]=%g\n", j, Per_y_+Block_Contain[j].Equation, r[j]);
                    residual[Per_y_+Block_Contain[j].Equation] = r[j];
                  }
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(r);
      break;
    case SOLVE_TWO_BOUNDARIES_SIMPLE:
    case SOLVE_TWO_BOUNDARIES_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, periods, y_kmin, y_kmax, steady_state, true, stack_solve_algo, solve_algo);
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = (double *) mxMalloc(size*sizeof(double));
      begining = it_code;
      for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
        {
          Per_u_ = (it_-y_kmin)*u_count_int;
          Per_y_ = it_*y_size;
          it_code = begining;
          compute_block_time(Per_u_, true, block_num, size, steady_state);
          if (block < 0)
            for (int j = 0; j < size; j++)
              {
                //mexPrintf("residual[it_*y_size+Block_Contain[%d].Equation = %d]=%g\n", j, it_*y_size+Block_Contain[j].Equation, r[j]);
                residual[it_*y_size+Block_Contain[j].Equation] = r[j];
              }
          else
            for (int j = 0; j < size; j++)
              residual[it_*size+j] = r[j];
        }
      mxFree(r);
      break;
    }
}


void
Interpreter::SingularDisplay(int Per_u_, bool evaluate, int Block_Count, int size, bool steady_state, it_code_type begining)
{
  it_code = begining;
  compute_block_time(Per_u_, evaluate, Block_Count, size, steady_state);
  Singular_display(Block_Count, size, steady_state, begining);
}


int
Interpreter::simulate_a_block(const int size, const int type, string file_name, string bin_basename, bool Gaussian_Elimination, bool steady_state, bool print_it, int block_num,
                              const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag, const int Block_List_Max_Lead, const int u_count_int)
{
  it_code_type begining;
  int i;
  bool cvg;
  bool result = true;
  bool singular_system;
  double *y_save;
  res1 = 0;
#ifdef DEBUG
  mexPrintf("simulate_a_block\n");
#endif
  switch (type)
    {
    case EVALUATE_FORWARD:
#ifdef DEBUG
      mexPrintf("EVALUATE_FORWARD\n");
#endif
      if (steady_state)
        compute_block_time(0, false, block_num, size, steady_state);
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, false, block_num, size, steady_state);
            }
        }
      break;
    case EVALUATE_BACKWARD:
#ifdef DEBUG
      mexPrintf("EVALUATE_BACKWARD\n");
#endif
      if (steady_state)
        compute_block_time(0, false, block_num, size, steady_state);
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, false, block_num, size, steady_state);
            }
        }
      break;
    case SOLVE_FORWARD_SIMPLE:
#ifdef DEBUG
      mexPrintf("SOLVE_FORWARD_SIMPLE size=%d\n", size);
#endif
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      begining = it_code;
      if (steady_state)
        {
          cvg = false;
          iter = 0;
          while (!(cvg || (iter > maxit_)))
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, false, block_num, size, steady_state);
              double rr;
              rr = r[0];
              cvg = (fabs(rr) < solve_tolf);
              if (cvg)
                continue;

              try
                {
                  y[Block_Contain[0].Variable] += -divide(r[0], g1[0]);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      \n", fpeh.GetErrorMsg().c_str());
                  mexPrintf("      Singularity in block %d", block_num+1);
                }
              iter++;
            }
          if (!cvg)
            {
              ostringstream tmp;
              tmp << " in Solve Forward simple, convergence not achieved in block " << Block_Count+1 << ", after " << iter << " iterations\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      else
        {
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              cvg = false;
              iter = 0;
              Per_y_ = it_*y_size;
              while (!(cvg || (iter > maxit_)))
                {
                  it_code = begining;
                  Per_y_ = it_*y_size;
                  compute_block_time(0, false, block_num, size, steady_state);
                  double rr;
                  if (fabs(1+y[Per_y_+Block_Contain[0].Variable]) > eps)
                    rr = r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
                  else
                    rr = r[0];
                  cvg = (fabs(rr) < solve_tolf);
                  if (cvg)
                    continue;
                  try
                    {
                      y[Per_y_+Block_Contain[0].Variable] += -divide(r[0], g1[0]);
                    }
                  catch (FloatingPointExceptionHandling &fpeh)
                    {
                      mexPrintf("%s      \n", fpeh.GetErrorMsg().c_str());
                      mexPrintf("      Singularity in block %d", block_num+1);
                    }
                  iter++;
                }
              if (!cvg)
                {
                  ostringstream tmp;
                  tmp << " in Solve Forward simple, convergence not achieved in block " << Block_Count+1 << ", at time " << it_ << ", after " << iter << " iterations\n";
                  throw FatalExceptionHandling(tmp.str());
                }
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_BACKWARD_SIMPLE:
#ifdef DEBUG
      mexPrintf("SOLVE_BACKWARD_SIMPLE\n");
#endif
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      begining = it_code;
      if (steady_state)
        {
          cvg = false;
          iter = 0;
          while (!(cvg || (iter > maxit_)))
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, false, block_num, size, steady_state);
              double rr;
              rr = r[0];
              cvg = (fabs(rr) < solve_tolf);
              if (cvg)
                continue;
              try
                {
                  y[Block_Contain[0].Variable] += -divide(r[0], g1[0]);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      \n", fpeh.GetErrorMsg().c_str());
                  mexPrintf("      Singularity in block %d", block_num+1);
                }
              iter++;
            }
          if (!cvg)
            {
              ostringstream tmp;
              tmp << " in Solve Backward simple, convergence not achieved in block " << Block_Count+1 << ", after " << iter << " iterations\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              cvg = false;
              iter = 0;
              Per_y_ = it_*y_size;
              while (!(cvg || (iter > maxit_)))
                {
                  it_code = begining;
                  Per_y_ = it_*y_size;
                  compute_block_time(0, false, block_num, size, steady_state);
                  double rr;
                  if (fabs(1+y[Per_y_+Block_Contain[0].Variable]) > eps)
                    rr = r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
                  else
                    rr = r[0];
                  cvg = (fabs(rr) < solve_tolf);
                  if (cvg)
                    continue;
                  try
                    {
                      y[Per_y_+Block_Contain[0].Variable] += -divide(r[0], g1[0]);
                    }
                  catch (FloatingPointExceptionHandling &fpeh)
                    {
                      mexPrintf("%s      \n", fpeh.GetErrorMsg().c_str());
                      mexPrintf("      Singularity in block %d", block_num+1);
                    }

                  iter++;
                }
              if (!cvg)
                {
                  ostringstream tmp;
                  tmp << " in Solve Backward simple, convergence not achieved in block " << Block_Count+1 << ", at time " << it_ << ", after " << iter << " iterations\n";
                  throw FatalExceptionHandling(tmp.str());
                }
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_FORWARD_COMPLETE:
#ifdef DEBUG
      mexPrintf("SOLVE_FORWARD_COMPLETE\n");
#endif
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false, stack_solve_algo, solve_algo);
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      begining = it_code;
      Per_u_ = 0;
      if (steady_state)
        {
          if (!is_linear)
            {
              max_res_idx = 0;
              cvg = false;
              iter = 0;
              glambda2 = g0 = very_big;
              try_at_iteration = 0;
              while (!(cvg || (iter > maxit_)))
                {
                  it_code = begining;
                  error_not_printed = true;
                  res2 = 0;
                  res1 = 0;
                  max_res = 0;
                  compute_block_time(0, false, block_num, size, steady_state);
                  if (!(isnan(res1) || isinf(res1)))
                    {
                      for (i = 0; i < size; i++)
                        {
                          double rr;
                          rr = r[i];
                          if (max_res < fabs(rr))
                            {
                              max_res = fabs(rr);
                              max_res_idx = i;
                            }
                          res2 += rr*rr;
                          res1 += fabs(rr);
                        }
                      cvg = (max_res < solve_tolf);
                    }
                  else
                    cvg = false;
                  if (cvg)
                    continue;
                  int prev_iter = iter;
                  singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, print_it, cvg, iter, true, stack_solve_algo, solve_algo);
                  if (singular_system)
                    SingularDisplay(0, false, block_num, size, steady_state, begining);

                  iter++;
                  if (iter > prev_iter)
                    {
                      g0 = res2;
                      gp0 = -res2;
                      try_at_iteration = 0;
                    }
                }
              if (!cvg || !result)
                {
                  ostringstream tmp;
                  tmp << " in Solve Forward complete, convergence not achieved in block " << Block_Count+1 << ", after " << iter << " iterations\n";
                  throw FatalExceptionHandling(tmp.str());
                }
            }
          else
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              iter = 0;
              res1 = 0;
              res2 = 0;
              max_res = 0;
              max_res_idx = 0;
              error_not_printed = true;
              compute_block_time(0, false, block_num, size, steady_state);
              if (!(isnan(res1) || isinf(res1)))
                {
                  for (i = 0; i < size; i++)
                    {
                      double rr;
                      rr = r[i];
                      if (max_res < fabs(rr))
                        {
                          max_res = fabs(rr);
                          max_res_idx = i;
                        }
                      res2 += rr*rr;
                      res1 += fabs(rr);
                    }
                  cvg = (max_res < solve_tolf);
                }
              else
                cvg = false;
              singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, print_it, cvg, iter, true, stack_solve_algo, solve_algo);
              if (singular_system)
                SingularDisplay(0, false, block_num, size, steady_state, begining);
              if (!result)
                {
                  mexPrintf(" in Solve Forward complete, convergence not achieved in block %d\n", Block_Count+1);
                  return ERROR_ON_EXIT;
                }
            }
        }
      else
        {
          if (!is_linear)
            {
              max_res_idx = 0;
              for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
                {
                  cvg = false;
                  iter = 0;
                  glambda2 = g0 = very_big;
                  try_at_iteration = 0;
                  Per_y_ = it_*y_size;
                  while (!(cvg || (iter > maxit_)))
                    {
                      it_code = begining;
                      error_not_printed = true;
                      res2 = 0;
                      res1 = 0;
                      max_res = 0;
                      compute_block_time(0, false, block_num, size, steady_state);
                      if (!(isnan(res1) || isinf(res1)))
                        {
                          for (i = 0; i < size; i++)
                            {
                              double rr;
                              if (fabs(1+y[Per_y_+Block_Contain[i].Variable]) > eps)
                                rr = r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                              else
                                rr = r[i];
                              if (max_res < fabs(rr))
                                {
                                  max_res = fabs(rr);
                                  max_res_idx = i;
                                }
                              res2 += rr*rr;
                              res1 += fabs(rr);
                            }
                          cvg = (max_res < solve_tolf);
                        }
                      else
                        cvg = false;
                      if (cvg)
                        continue;
                      int prev_iter = iter;
                      singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, print_it, cvg, iter, false, stack_solve_algo, solve_algo);
                      if (singular_system)
                        SingularDisplay(0, false, block_num, size, steady_state, begining);
                      iter++;
                      if (iter > prev_iter)
                        {
                          g0 = res2;
                          gp0 = -res2;
                          try_at_iteration = 0;
                        }
                    }
                  if (!cvg)
                    {
                      ostringstream tmp;
                      tmp << " in Solve Forward complete, convergence not achieved in block " << Block_Count+1 << ", at time " << it_ << ", after " << iter << " iterations\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                }
            }
          else
            {
              for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
                {
                  it_code = begining;
                  Per_y_ = it_*y_size;
                  iter = 0;
                  res1 = 0;
                  res2 = 0;
                  max_res = 0; max_res_idx = 0;
                  error_not_printed = true;
                  compute_block_time(0, false, block_num, size, steady_state);
                  if (!(isnan(res1) || isinf(res1)))
                    {
                      for (i = 0; i < size; i++)
                        {
                          double rr;
                          rr = r[i];
                          if (max_res < fabs(rr))
                            {
                              max_res = fabs(rr);
                              max_res_idx = i;
                            }
                          res2 += rr*rr;
                          res1 += fabs(rr);
                        }
                      cvg = (max_res < solve_tolf);
                    }
                  else
                    cvg = false;
                  singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, print_it, cvg, iter, false, stack_solve_algo, solve_algo);
                  if (singular_system)
                    SingularDisplay(0, false, block_num, size, steady_state, begining);
                }
            }
        }
      mxFree(index_equa);
      mxFree(index_vara);
      memset(direction, 0, size_of_direction);
      mxFree(g1);
      mxFree(r);
      mxFree(u);
      break;
    case SOLVE_BACKWARD_COMPLETE:
#ifdef DEBUG
      mexPrintf("SOLVE_BACKWARD_COMPLETE\n");
#endif
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false, stack_solve_algo, solve_algo);
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      begining = it_code;
      if (steady_state)
        {
          if (!is_linear)
            {
              max_res_idx = 0;
              cvg = false;
              iter = 0;
              glambda2 = g0 = very_big;
              try_at_iteration = 0;
              while (!(cvg || (iter > maxit_)))
                {
                  it_code = begining;
                  error_not_printed = true;
                  res2 = 0;
                  res1 = 0;
                  max_res = 0;
                  compute_block_time(0, false, block_num, size, steady_state);
                  if (!(isnan(res1) || isinf(res1)))
                    {
                      for (i = 0; i < size; i++)
                        {
                          double rr;
                          rr = r[i];
                          if (max_res < fabs(rr))
                            {
                              max_res = fabs(rr);
                              max_res_idx = i;
                            }
                          res2 += rr*rr;
                          res1 += fabs(rr);
                        }
                      cvg = (max_res < solve_tolf);
                    }
                  else
                    cvg = false;
                  if (cvg)
                    continue;
                  int prev_iter = iter;
                  singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, print_it, cvg, iter, true, stack_solve_algo, solve_algo);
                  if (singular_system)
                    SingularDisplay(0, false, block_num, size, steady_state, begining);
                  iter++;
                  if (iter > prev_iter)
                    {
                      g0 = res2;
                      gp0 = -res2;
                      try_at_iteration = 0;
                    }
                }
              if (!cvg || !result)
                {
                  ostringstream tmp;
                  tmp << " in Solve Backward complete, convergence not achieved in block " << Block_Count+1 << ", after " << iter << " iterations\n";
                  throw FatalExceptionHandling(tmp.str());
                }
            }
          else
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              iter = 0;
              res1 = 0;
              res2 = 0;
              max_res = 0; max_res_idx = 0;
              error_not_printed = true;
              compute_block_time(0, false, block_num, size, steady_state);
              if (!(isnan(res1) || isinf(res1)))
                {
                  for (i = 0; i < size; i++)
                    {
                      double rr;
                      rr = r[i];
                      if (max_res < fabs(rr))
                        {
                          max_res = fabs(rr);
                          max_res_idx = i;
                        }
                      res2 += rr*rr;
                      res1 += fabs(rr);
                    }
                  cvg = (max_res < solve_tolf);
                }
              else
                cvg = false;
              singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, print_it, cvg, iter, true, stack_solve_algo, solve_algo);
              if (singular_system)
                SingularDisplay(0, false, block_num, size, steady_state, begining);
              if (!result)
                {
                  mexPrintf(" in Solve Backward complete, convergence not achieved in block %d\n", Block_Count+1);
                  return ERROR_ON_EXIT;
                }
            }
        }
      else
        {
          if (!is_linear)
            {
              max_res_idx = 0;
              for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
                {
                  cvg = false;
                  iter = 0;
                  glambda2 = g0 = very_big;
                  try_at_iteration = 0;
                  Per_y_ = it_*y_size;
                  while (!(cvg || (iter > maxit_)))
                    {
                      it_code = begining;
                      error_not_printed = true;
                      res2 = 0;
                      res1 = 0;
                      max_res = 0;
                      compute_block_time(0, false, block_num, size, steady_state);
                      if (!(isnan(res1) || isinf(res1)))
                        {
                          for (i = 0; i < size; i++)
                            {
                              double rr;
                              if (fabs(1+y[Per_y_+Block_Contain[i].Variable]) > eps)
                                rr = r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                              else
                                rr = r[i];
                              if (max_res < fabs(rr))
                                {
                                  max_res = fabs(rr);
                                  max_res_idx = i;
                                }
                              res2 += rr*rr;
                              res1 += fabs(rr);
                            }
                          cvg = (max_res < solve_tolf);
                        }
                      else
                        cvg = false;
                      if (cvg)
                        continue;
                      int prev_iter = iter;
                      singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, print_it, cvg, iter, false, stack_solve_algo, solve_algo);
                      if (singular_system)
                        SingularDisplay(0, false, block_num, size, steady_state, begining);
                      iter++;
                      if (iter > prev_iter)
                        {
                          g0 = res2;
                          gp0 = -res2;
                          try_at_iteration = 0;
                        }
                    }
                  if (!cvg)
                    {
                      ostringstream tmp;
                      tmp << " in Solve Backward complete, convergence not achieved in block " << Block_Count+1 << ", at time " << it_ << ", after " << iter << " iterations\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                }
            }
          else
            {
              for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
                {
                  it_code = begining;
                  Per_y_ = it_*y_size;
                  error_not_printed = true;
                  compute_block_time(0, false, block_num, size, steady_state);
                  if (!(isnan(res1) || isinf(res1)))
                    {
                      for (i = 0; i < size; i++)
                        {
                          double rr;
                          rr = r[i];
                          if (max_res < fabs(rr))
                            {
                              max_res = fabs(rr);
                              max_res_idx = i;
                            }
                          res2 += rr*rr;
                          res1 += fabs(rr);
                        }
                      cvg = (max_res < solve_tolf);
                    }
                  else
                    cvg = false;
                  singular_system = Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, print_it, cvg, iter, false, stack_solve_algo, solve_algo);
                  if (singular_system)
                    SingularDisplay(0, false, block_num, size, steady_state, begining);
                }
            }
        }
      mxFree(index_equa);
      mxFree(index_vara);
      memset(direction, 0, size_of_direction);
      mxFree(g1);
      mxFree(r);
      mxFree(u);
      break;
    case SOLVE_TWO_BOUNDARIES_SIMPLE:
    case SOLVE_TWO_BOUNDARIES_COMPLETE:
#ifdef DEBUG
      mexPrintf("SOLVE_TWO_BOUNDARIES\n");
#endif
      if (steady_state)
        {
          mexPrintf("SOLVE_TWO_BOUNDARIES in a steady state model: impossible case\n");
          return ERROR_ON_EXIT;
        }
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, periods, y_kmin, y_kmax, steady_state, true, stack_solve_algo, solve_algo);
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = (double *) mxMalloc(size*sizeof(double));
      y_save = (double *) mxMalloc(y_size*sizeof(double)*(periods+y_kmax+y_kmin));
      begining = it_code;
      iter = 0;
      if (!is_linear)
        {
          cvg = false;
          glambda2 = g0 = very_big;
          try_at_iteration = 0;
          int u_count_saved = u_count;
          while (!(cvg || (iter > maxit_)))
            {

              res2 = 0;
              res1 = 0;
              max_res = 0;
              max_res_idx = 0;
              memcpy(y_save, y, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
              for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
                {
                  Per_u_ = (it_-y_kmin)*u_count_int;
                  Per_y_ = it_*y_size;
                  it_code = begining;
                  error_not_printed = true;
                  compute_block_time(Per_u_, false, block_num, size, steady_state);
                  if (isnan(res1) || isinf(res1))
                    {
                      memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                      break;
                    }
                  for (i = 0; i < size; i++)
                    {
                      double rr;
                      if (fabs(1+y[Per_y_+Block_Contain[i].Variable]) > eps)
                        rr = r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                      else
                        rr = r[i];
                      if (max_res < fabs(rr))
                        {
                          max_res = fabs(rr);
                          max_res_idx = i;
                        }
                      res2 += rr*rr;
                      res1 += fabs(rr);
                    }
                }
              if (isnan(res1) || isinf(res1))
                cvg = false;
              else
                cvg = (max_res < solve_tolf);
              u_count = u_count_saved;
              int prev_iter = iter;
              Simulate_Newton_Two_Boundaries(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, print_it, cvg, iter, minimal_solving_periods, stack_solve_algo, endo_name_length, P_endo_names);
              iter++;
              if (iter > prev_iter)
                {
                  g0 = res2;
                  gp0 = -res2;
                  try_at_iteration = 0;
                  slowc_save = slowc;
                }
            }
          if (!cvg)
            {
              ostringstream tmp;
              tmp << " in Solve two boundaries, convergence not achieved in block " << Block_Count+1 << ", after " << iter << " iterations\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      else
        {
          res1 = 0;
          res2 = 0;
          max_res = 0; max_res_idx = 0;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              Per_u_ = (it_-y_kmin)*u_count_int;
              Per_y_ = it_*y_size;
              it_code = begining;
              compute_block_time(Per_u_, false, block_num, size, steady_state);
              for (i = 0; i < size; i++)
                {
                  double rr;
                  rr = r[i];
                  if (max_res < fabs(rr))
                    {
                      max_res = fabs(rr);
                      max_res_idx = i;
                    }
                  res2 += rr*rr;
                  res1 += fabs(rr);
                }
            }
          cvg = false;
          Simulate_Newton_Two_Boundaries(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, print_it, cvg, iter, minimal_solving_periods, stack_solve_algo, endo_name_length, P_endo_names);
        }
      mxFree(r);
      mxFree(y_save);
      mxFree(u);
      mxFree(index_vara);
      mxFree(index_equa);
      memset(direction, 0, size_of_direction);
      break;
    default:
      ostringstream tmp;
      tmp << " in simulate_a_block, Unknown type = " << type << "\n";
      throw FatalExceptionHandling(tmp.str());
      return ERROR_ON_EXIT;
    }
  return NO_ERROR_ON_EXIT;
}

void
Interpreter::print_a_block(const int size, const int type, string bin_basename, bool steady_state, int block_num,
                           const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag,
                           const int Block_List_Max_Lead, const int u_count_int, int block)
{
  it_code_type begining;
  mexPrintf("\nBlock %d\n", block_num+1);
  mexPrintf("----------\n");
  if (steady_state)
    residual = vector<double>(size);
  else
    residual = vector<double>(size*(periods+y_kmin));
  bool go_on = true;
  bool space = false;
  while (go_on)
    {
      if (it_code->first == FENDBLOCK)
        {
          go_on = false;
          it_code++;
        }
      else
        {
          string s = print_expression(it_code, false, size, block_num, steady_state, Per_u_, it_, it_code, false);
          if (s == "if (evaluate)" || s == "else")
            space = false;
          if (s.length() > 0)
            {
              if (space)
                mexPrintf("  %s\n", s.c_str());
              else
                mexPrintf("%s\n", s.c_str());
              mexEvalString("drawnow;");
            }
          if (s == "if (evaluate)" || s == "else")
            space = true;
        }
    }
}

bool
Interpreter::compute_blocks(string file_name, string bin_basename, bool steady_state, bool evaluate, int block, int &nb_blocks, bool print_it)
{
  bool result = true;

  int var;
  if (steady_state)
    file_name += "_static";
  else
    file_name += "_dynamic";
  CodeLoad code;
  //First read and store in memory the code
  code_liste = code.get_op_code(file_name);
  EQN_block_number = code.get_block_number();
  if (!code_liste.size())
    {
      ostringstream tmp;
      tmp << " in compute_blocks, " << file_name.c_str() << " cannot be opened\n";
      throw FatalExceptionHandling(tmp.str());
    }
  if (block >= (int) code.get_block_number())
    {
      ostringstream tmp;
      tmp << " in compute_blocks, input argument block = " << block+1 << " is greater than the number of blocks in the model (" << code.get_block_number() << " see M_.blocksMFS)\n";
      throw FatalExceptionHandling(tmp.str());
    }
  //The big loop on intructions
  Block_Count = -1;
  bool go_on = true;

  it_code = code_liste.begin();
  it_code_type Init_Code = it_code;
  if (block < 0)
    {
      if (steady_state)
        residual = vector<double>(y_size);
      else
        residual = vector<double>(y_size*(periods+y_kmin));
    }

  while (go_on)
    {
      switch (it_code->first)
        {
        case FBEGINBLOCK:
          Block_Count++;
#ifdef DEBUG
          mexPrintf("---------------------------------------------------------\n");
          if (block < 0)
            mexPrintf("FBEGINBLOCK Block_Count=%d\n", Block_Count+1);
          else
            mexPrintf("FBEGINBLOCK block=%d\n", block+1);
#endif
          //it's a new block
          {
            FBEGINBLOCK_ *fb = (FBEGINBLOCK_ *) it_code->second;
            Block_Contain = fb->get_Block_Contain();
            it_code++;
            if (print)
              print_a_block(fb->get_size(), fb->get_type(), bin_basename, steady_state, Block_Count,
                            fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int(), block);
            else if (evaluate)
              {
#ifdef DEBUG
                mexPrintf("jacobian_block=mxCreateDoubleMatrix(%d, %d, mxREAL)\n", fb->get_size(), fb->get_nb_col_jacob());
#endif
                jacobian_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_jacob(), mxREAL));
                if (!steady_state)
                  {
#ifdef DEBUG
                    mexPrintf("allocates jacobian_exo_block( %d, %d, mxREAL)\n", fb->get_size(), fb->get_exo_size());
#endif
                    jacobian_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_exo_size(), mxREAL));
                    jacobian_det_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_det_exo_size(), mxREAL));
                    jacobian_other_endo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_other_endo_jacob(), mxREAL));
                  }
                if (block >= 0)
                  {
                    if (steady_state)
                      residual = vector<double>(fb->get_size());
                    else
                      residual = vector<double>(fb->get_size()*(periods+y_kmin));
                  }
                evaluate_a_block(fb->get_size(), fb->get_type(), bin_basename, steady_state, Block_Count,
                                 fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int(), block);
              }
            else
              {
#ifdef DEBUG
                mexPrintf("endo in block=%d, type=%d, steady_state=%d, print_it=%d, Block_Count=%d, fb->get_is_linear()=%d, fb->get_endo_nbr()=%d, fb->get_Max_Lag()=%d, fb->get_Max_Lead()=%d, fb->get_u_count_int()=%d\n",
                          fb->get_size(), fb->get_type(), steady_state, print_it, Block_Count, fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
#endif
                result = simulate_a_block(fb->get_size(), fb->get_type(), file_name, bin_basename, true, steady_state, print_it,Block_Count,
                                          fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
                if (result == ERROR_ON_EXIT)
                  return ERROR_ON_EXIT;
              }
            delete fb;
          }
          if (block >= 0)
            {

              go_on = false;
            }

          break;
        case FEND:
#ifdef DEBUG
          mexPrintf("FEND\n");
#endif
          go_on = false;
          it_code++;
          break;
        case FDIMT:
#ifdef DEBUG
          mexPrintf("FDIMT size=%d\n", ((FDIMT_ *) it_code->second)->get_size());
#endif
          var = ((FDIMT_ *) it_code->second)->get_size();
          if (T)
            mxFree(T);
          T = (double *) mxMalloc(var*(periods+y_kmin+y_kmax)*sizeof(double));
          if (block >= 0)
            {
              it_code = code_liste.begin() + code.get_begin_block(block);
            }
          else
            it_code++;
          break;
        case FDIMST:
#ifdef DEBUG
          mexPrintf("FDIMST size=%d\n", ((FDIMST_ *) it_code->second)->get_size());
#endif
          var = ((FDIMST_ *) it_code->second)->get_size();
          if (T)
            mxFree(T);
          if (global_temporary_terms)
            {
              if (GlobalTemporaryTerms == NULL)
                {
                  mexPrintf("GlobalTemporaryTerms is NULL\n"); mexEvalString("drawnow;");
                }
              if (var != (int) mxGetNumberOfElements(GlobalTemporaryTerms))
                GlobalTemporaryTerms = mxCreateDoubleMatrix(var, 1, mxREAL);
              T = mxGetPr(GlobalTemporaryTerms);
            }
          else
            T = (double *) mxMalloc(var*sizeof(double));

          if (block >= 0)
            it_code = code_liste.begin() + code.get_begin_block(block);
          else
            it_code++;
          break;
        default:
          mexPrintf("Error\n");
          ostringstream tmp;
          tmp << " in compute_blocks, unknown command " << it_code->first << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
    }

  mxFree(Init_Code->second);
  nb_blocks = Block_Count+1;
  if (T and !global_temporary_terms)
    mxFree(T);
  return result;
}
