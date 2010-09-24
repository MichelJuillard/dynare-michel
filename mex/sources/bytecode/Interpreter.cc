/*
 * Copyright (C) 2007-2010 Dynare Team
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
                         string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg)
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
  mxArray *M_ = mexGetVariable("global", "M_");
  nb_endo = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
  endo_name_length = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
  P_endo_names = (char*) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
  nb_exo = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
  exo_name_length = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
  P_exo_names = (char*) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
  nb_param = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
  param_name_length = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
  P_param_names = (char*) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
}

string
Interpreter::add_underscore_to_fpe(const string &str)
{
  string temp;
  int pos1 = -1, pos2 = -1;
  string tmp_n(str.length(), ' ');
  for (unsigned int i = 0; i < str.length(); i++)
    {
      if (str[i] != '$' && str[i] != '£')
        temp += str[i];
      else
        {
          if (str[i] == '$')
            pos1 = temp.length();
          else
            pos2 = temp.length();
          if (pos1>=0 && pos2 >=0)
            {
              tmp_n.erase(pos1, pos2-pos1+1);
              tmp_n.insert(pos1, pos2-pos1, '~');
              pos1 = pos2 = -1;
            }
        }
    }
  temp += "\n" + tmp_n ;
  return temp;
}


string
Interpreter::get_variable(const SymbolType variable_type, const unsigned int variable_num)
{
  ostringstream res;
  switch(variable_type)
    {
    case eEndogenous:
      if (variable_num < nb_endo)
        {
          for (unsigned int i = 0; i < endo_name_length; i++)
            if (P_endo_names[2*(variable_num+i*nb_endo)] != ' ')
              res << P_endo_names[2*(variable_num+i*nb_endo)];
        }
      else
        mexPrintf("=> Unknown endogenous variable n° %d",variable_num);
      break;
    case eExogenous:
    case eExogenousDet:
      if (variable_num < nb_exo)
        {
          for (unsigned int i = 0; i < exo_name_length; i++)
            if (P_exo_names[2*(variable_num+i*nb_exo)] != ' ')
              res << P_exo_names[2*(variable_num+i*nb_exo)];
        }
      else
        mexPrintf("=> Unknown exogenous variable n° %d",variable_num);
      break;
    case eParameter:
      if (variable_num < nb_param)
        {
          for (unsigned int i = 0; i < param_name_length; i++)
            if (P_param_names[2*(variable_num+i*nb_param)] != ' ')
              res << P_param_names[2*(variable_num+i*nb_param)];
        }
      else
        mexPrintf("=> Unknown parameter n° %d",variable_num);
      break;
    default:
      break;
    }
  return(res.str());
}


string
Interpreter::error_location(bool evaluate, bool steady_state)
{
  stringstream Error_loc("");
  if (!steady_state)
    switch(EQN_type)
      {
        case TemporaryTerm:
          if (EQN_block_number > 1)
            Error_loc << "temporary term " << EQN_equation+1 << " in block " << EQN_block+1 << " at time " << it_;
          else
            Error_loc << "temporary term " << EQN_equation+1 << " at time " << it_;
          break;
        case ModelEquation:
          if (EQN_block_number > 1)
            Error_loc << "equation " << EQN_equation+1 << " in block " << EQN_block+1 << " at time " << it_;
          else
            Error_loc << "equation " << EQN_equation+1 << " at time " << it_;
          break;
        case FirstEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to endogenous variable " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to endogenous variable " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          break;
        default:
          return("???");
          break;
      }
  else
    switch(EQN_type)
      {
        case TemporaryTerm:
          if (EQN_block_number > 1)
            Error_loc << "temporary term " << EQN_equation+1 << " in block " << EQN_block+1;
          else
            Error_loc << "temporary term " << EQN_equation+1;
          break;
        case ModelEquation:
          if (EQN_block_number > 1)
            Error_loc << "equation " << EQN_equation+1 << " in block " << EQN_block+1;
          else
            Error_loc << "equation " << EQN_equation+1;
          break;
        case FirstEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to endogenous variable "  << get_variable(eEndogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to endogenous variable " << get_variable(eEndogenous, EQN_dvar1);
          break;
        default:
          return("???");
          break;
      }
  Error_loc << endl << add_underscore_to_fpe("      " + print_expression(it_code_expr, evaluate));
  return(Error_loc.str());
}

double
Interpreter::pow1(double a, double b)
{
  double r = pow_(a, b);
  if (isnan(r))
    {
      res1 = NAN;
      r = 0.0000000000000000000000001;
      throw PowExceptionHandling(a, b);
    }
  return r;
}

double
Interpreter::divide(double a, double b)
{
  double r = a/ b;
  if (isinf(r))
    {
      res1 = NAN;
      r = 1e70;
      throw PowExceptionHandling(a, b);
    }
  return r;
}

double
Interpreter::log1(double a)
{
  double r = log(a);
  if (isnan(r))
    {
      res1 = NAN;
      r = -1e70;
      throw LogExceptionHandling(a);
    }
  return r;
}

double
Interpreter::log10_1(double a)
{
  double r = log(a);
  if (isnan(r))
    {
      res1 = NAN;
      r = -1e70;
      throw Log10ExceptionHandling(a);
    }
  return r;
}


string
Interpreter::print_expression(it_code_type it_code, bool evaluate)
{
  int var, lag = 0, op;
  stack<string> Stack;
  stack<double> Stackf;
  ostringstream tmp_out, tmp_out2;
  string v1, v2, v3;
  double v1f, v2f, v3f;
  bool go_on = true;
  double ll;
  ExpressionType equation_type;
  unsigned int equation_num;
  unsigned int dvar1, dvar2, dvar3;
  int lag1, lag2, lag3;
  size_t found;

  while (go_on)
    {
      switch (it_code->first)
        {
        case FNUMEXPR:
          switch (((FNUMEXPR_ *) it_code->second)->get_expression_type())
            {
            case TemporaryTerm:
              equation_type = TemporaryTerm;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              break;
            case ModelEquation:
              equation_type = ModelEquation;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              break;
            case FirstEndoDerivative:
              equation_type = FirstEndoDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExoDerivative:
              equation_type = FirstExoDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExodetDerivative:
              equation_type = FirstExodetDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstParamDerivative:
              equation_type = FirstParamDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              break;
            case SecondEndoDerivative:
              equation_type = SecondEndoDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExoDerivative:
              equation_type = SecondExoDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExodetDerivative:
              equation_type = SecondExodetDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondParamDerivative:
              equation_type = SecondExodetDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              break;
            case ThirdEndoDerivative:
              equation_type = ThirdEndoDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdExoDerivative:
              equation_type = ThirdExoDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdExodetDerivative:
              equation_type = ThirdExodetDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdParamDerivative:
              equation_type = ThirdExodetDerivative;
              equation_num = ((FNUMEXPR_ *) it_code->second)->get_equation();
              dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              break;
            default:
              ostringstream tmp;
              tmp << " in print_expression, derivatives " << it_code->first << " not implemented yet\n";
              throw FatalExceptionHandling(tmp.str());
            }
          break;
        case FLDV:
          //load a variable in the processor
          switch (((FLDV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDV_ *) it_code->second)->get_pos();
              Stack.push(get_variable(eParameter, var));
              Stackf.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
              tmp_out.str("");
              if(lag != 0)
                tmp_out << get_variable(eEndogenous, var) << "(" << lag << ")";
              else
                tmp_out << get_variable(eEndogenous, var);
              Stack.push(tmp_out.str());
              if (evaluate)
                Stackf.push(ya[(it_+lag)*y_size+var]);
              else
                Stackf.push(y[(it_+lag)*y_size+var]);
              break;
            case eExogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
              tmp_out.str("");
              if(lag != 0)
                tmp_out << get_variable(eExogenous, var) << "(" << lag << ")";
              else
                tmp_out << get_variable(eExogenous, var);
              Stack.push(tmp_out.str());
              Stackf.push(x[it_+lag+var*nb_row_x]);
              break;
            case eExogenousDet:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
              tmp_out.str("");
              if(lag != 0)
                tmp_out << get_variable(eExogenousDet, var) << "(" << lag << ")";
              else
                tmp_out << get_variable(eExogenousDet, var);
              Stack.push(tmp_out.str());
              Stackf.push(x[it_+lag+var*nb_row_xd]);
              break;
            case eModelLocalVariable:
              break;
            default:
              mexPrintf("FLDV: Unknown variable type\n");
            }
          break;
        case FLDSV:
        case FLDVS:
          //load a variable in the processor
          switch (((FLDSV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              Stack.push(get_variable(eParameter, var));
              Stackf.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              Stack.push(get_variable(eEndogenous, var));
              if (it_code->first == FLDSV)
                {
                  if (evaluate)
                    Stackf.push(ya[var]);
                  else
                    Stackf.push(y[var]);
                }
              else
                Stackf.push(steady_y[var]);
              break;
            case eExogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              Stack.push(get_variable(eExogenous, var));
              break;
              Stackf.push(x[var]);
            case eExogenousDet:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              Stack.push(get_variable(eExogenousDet, var));
              Stackf.push(x[var]);
              break;
            case eModelLocalVariable:
              break;
            default:
              mexPrintf("FLDSV: Unknown variable type\n");
            }
          break;
        case FLDT:
          //load a temporary variable in the processor
          var = ((FLDT_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "T" << var+1;
          Stack.push(tmp_out.str());
          Stackf.push(T[var*(periods+y_kmin+y_kmax)+it_]);
          break;
        case FLDST:
          //load a temporary variable in the processor
          var = ((FLDT_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "T" << var+1;
          Stack.push(tmp_out.str());
          Stackf.push(T[var]);
          break;
        case FLDU:
          //load u variable in the processor
          var = ((FLDU_ *) it_code->second)->get_pos();
          var += Per_u_;
          tmp_out.str("");
          tmp_out << "u[" << var+1 << "]";
          Stack.push(tmp_out.str());
          Stackf.push(u[var]);
          break;
        case FLDSU:
          //load u variable in the processor
          var = ((FLDSU_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "u[" << var+1 << "]";
          Stack.push(tmp_out.str());
          Stackf.push(u[var]);
          break;
        case FLDR:
          var = ((FLDR_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "residual[" << var+1 << "]";
          Stack.push(tmp_out.str());
          Stackf.push(r[var]);
          break;
        case FLDZ:
          //load 0 in the processor
          tmp_out.str("");
          tmp_out << 0;
          Stack.push(tmp_out.str());
          Stackf.push(0.0);
          break;
        case FLDC:
          //load a numerical constant in the processor
          ll = ((FLDC_ *) it_code->second)->get_value();
          tmp_out.str("");
          tmp_out << ll;
          Stack.push(tmp_out.str());
          Stackf.push(ll);
          break;
        case FSTPV:
          //load a variable in the processor
          go_on = false;
          switch (((FSTPV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              tmp_out2.str("");
              tmp_out2 << Stack.top();
              tmp_out.str("");
              tmp_out << get_variable(eParameter, var) << " = " << tmp_out2.str();
              Stack.pop();
              params[var] = Stackf.top();
              Stackf.pop();
              break;
            case eEndogenous:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              tmp_out2.str("");
              tmp_out2 << Stack.top();
              tmp_out.str("");
              tmp_out << get_variable(eEndogenous, var);
              if (lag != 0)
                tmp_out << "(" << lag << ")";
              tmp_out << " = " << tmp_out2.str();
              Stack.pop();
              y[(it_+lag)*y_size+var] = Stackf.top();
              Stackf.pop();
              break;
            case eExogenous:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              tmp_out2.str("");
              tmp_out2 << Stack.top();
              tmp_out.str("");
              tmp_out << get_variable(eExogenous, var);
              if (lag != 0)
                tmp_out << "(" << lag << ")";
              tmp_out << " = " << tmp_out2.str();
              Stack.pop();
              x[it_+lag+var*nb_row_x]  = Stackf.top();
              Stackf.pop();
              break;
            case eExogenousDet:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              tmp_out2.str("");
              tmp_out2 << Stack.top();
              tmp_out.str("");
              tmp_out << get_variable(eExogenousDet, var);
              if (lag != 0)
                tmp_out << "(" << lag << ")";
              tmp_out << " = " << tmp_out2.str();
              Stack.pop();
              x[it_+lag+var*nb_row_xd] = Stackf.top();
              Stackf.pop();
              break;
            default:
              mexPrintf("FSTPV: Unknown variable type\n");
            }
          break;
        case FSTPSV:
          go_on = false;
          //load a variable in the processor
          switch (((FSTPSV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              tmp_out2.str("");
              tmp_out2 << Stack.top();
              tmp_out.str("");
              tmp_out << get_variable(eParameter, var);
              tmp_out << " = " << tmp_out2.str();
              Stack.pop();
              params[var] = Stackf.top();
              Stackf.pop();
              break;
            case eEndogenous:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              tmp_out2.str("");
              tmp_out2 << Stack.top();
              tmp_out.str("");
              tmp_out << get_variable(eEndogenous, var);
              tmp_out << " = " << tmp_out2.str();
              Stack.pop();
              y[var] = Stackf.top();
              Stackf.pop();
              break;
            case eExogenous:
            case eExogenousDet:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              tmp_out2.str("");
              tmp_out2 << Stack.top();
              tmp_out.str("");
              tmp_out << get_variable(eExogenous, var);
              tmp_out << " = " << tmp_out2.str();
              Stack.pop();
              x[var]  = Stackf.top();
              Stackf.pop();
              break;
            default:
              mexPrintf("FSTPSV: Unknown variable type\n");
            }
          break;
        case FSTPT:
          go_on = false;
          //store in a temporary variable from the processor
          var = ((FSTPT_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "T" << var+1 << " = " << Stack.top();
          Stack.pop();
          T[var*(periods+y_kmin+y_kmax)+it_] = Stackf.top();
          Stackf.pop();
          break;
        case FSTPST:
          go_on = false;
          //store in a temporary variable from the processor
          var = ((FSTPT_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "T" << var+1 << " = " << Stack.top();
          Stack.pop();
          T[var] = Stackf.top();
          Stackf.pop();
          break;
        case FSTPU:
          go_on = false;
          //store in u variable from the processor
          var = ((FSTPU_ *) it_code->second)->get_pos();
          var += Per_u_;
          tmp_out.str("");
          tmp_out << "u[" << var+1 << "] = " << Stack.top();
          Stack.pop();
          u[var] = Stackf.top();
          Stackf.pop();
          break;
        case FSTPSU:
          go_on = false;
          //store in u variable from the processor
          var = ((FSTPU_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "u[" << var+1 << "] = " << Stack.top();
          Stack.pop();
          u[var] = Stackf.top();
          Stackf.pop();
          break;
        case FSTPR:
          go_on = false;
          //store in residual variable from the processor
          var = ((FSTPR_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "residual[" << var+1 << "] = " << Stack.top();
          Stack.pop();
          r[var] = Stackf.top();
          Stackf.pop();
          break;
        case FSTPG:
          go_on = false;
          //store in derivative (g) variable from the processor
          var = ((FSTPG_ *) it_code->second)->get_pos();
          tmp_out.str("");
          tmp_out << "g1[" << var+1 << "] = " << Stack.top();
          Stack.pop();
          g1[var] = Stackf.top();
          Stackf.pop();
          break;
        case FBINARY:
          op = ((FBINARY_ *) it_code->second)->get_op_type();
          v2 = Stack.top();
          Stack.pop();
          v1 = Stack.top();
          Stack.pop();
          v2f = Stackf.top();
          Stackf.pop();
          v1f = Stackf.top();
          Stackf.pop();
          switch (op)
            {
            case oPlus:
              Stackf.push(v1f + v2f);
              tmp_out.str("");
              tmp_out << v1 << " + " << v2;
              Stack.push(tmp_out.str());
              break;
            case oMinus:
              Stackf.push(v1f - v2f);
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " - ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oTimes:
              Stackf.push(v1f * v2f);
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " * ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oDivide:
              double r;
              r = v1f / v2f;
              Stackf.push(r);
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              if (isinf(r))
                tmp_out << "$";
              tmp_out << " / ";
              if (isinf(r))
                tmp_out << "£";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oLess:
              Stackf.push(double (v1f < v2f));
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " < ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oGreater:
              Stackf.push(double (v1f > v2f));
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " > ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oLessEqual:
              Stackf.push(double (v1f <= v2f));
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " <= ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oGreaterEqual:
              Stackf.push(double (v1f >= v2f));
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " >= ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oEqualEqual:
              Stackf.push(double (v1f == v2f));
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " == ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oDifferent:
              Stackf.push(double (v1f != v2f));
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              tmp_out << " != ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oPower:
              r = pow(v1f, v2f);
              Stackf.push(r);
              tmp_out.str("");
              found = v1.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v1;
              if (found != string::npos)
                tmp_out << ")";
              if(isnan(r))
                tmp_out << "$ ^ £";
              else
                tmp_out << " ^ ";
              found = v2.find(" ");
              if (found != string::npos)
                tmp_out << "(";
              tmp_out << v2;
              if (found != string::npos)
                tmp_out << ")";
              Stack.push(tmp_out.str());
              break;
            case oMax:
              Stackf.push(max(v1f, v2f));
              tmp_out.str("");
              tmp_out << "max(" << v1 << ", " << v2 << ")";
              Stack.push(tmp_out.str());
              break;
            case oMin:
              Stackf.push(min(v1f, v2f));
              tmp_out.str("");
              tmp_out << "min(" << v1 << ", " << v2 << ")";
              Stack.push(tmp_out.str());
              break;
            case oEqual:
            default:
              /*throw EvalException();*/
              ;
            }
          break;
        case FUNARY:
          op = ((FUNARY_ *) it_code->second)->get_op_type();
          v1 = Stack.top();
          Stack.pop();
          v1f = Stackf.top();
          Stackf.pop();
          double r;
          switch (op)
            {
            case oUminus:
              Stackf.push(-v1f);
              tmp_out.str("");
              tmp_out << " -" << v1;
              Stack.push(tmp_out.str());
              break;
            case oExp:
              Stackf.push(exp(v1f));
              tmp_out.str("");
              tmp_out << "exp(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oLog:
              r = log(v1f);
              Stackf.push(r);
              tmp_out.str("");
              if (isnan(r))
                tmp_out << "$log£(" << v1 << ")";
              else
                tmp_out << "log(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oLog10:
              r = log10(v1f);
              Stackf.push(r);
              tmp_out.str("");
              if (isnan(r))
                tmp_out << "$log10£(" << v1 << ")";
              else
                tmp_out << "log10(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oCos:
              Stackf.push(cos(v1f));
              tmp_out.str("");
              tmp_out << "cos(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oSin:
              Stackf.push(sin(v1f));
              tmp_out.str("");
              tmp_out << "sin(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oTan:
              Stackf.push(tan(v1f));
              tmp_out.str("");
              tmp_out << "tan(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oAcos:
              Stackf.push(acos(v1f));
              tmp_out.str("");
              tmp_out << "acos(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oAsin:
              Stackf.push(asin(v1f));
              tmp_out.str("");
              tmp_out << "asin(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oAtan:
              Stackf.push(atan(v1f));
              tmp_out.str("");
              tmp_out << "atan(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oCosh:
              Stackf.push(cosh(v1f));
              tmp_out.str("");
              tmp_out << "cosh(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oSinh:
              Stackf.push(sinh(v1f));
              tmp_out.str("");
              tmp_out << "sinh(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oTanh:
              Stackf.push(tanh(v1f));
              tmp_out.str("");
              tmp_out << "tanh(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oAcosh:
              Stackf.push(acosh(v1f));
              tmp_out.str("");
              tmp_out << "acosh(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oAsinh:
              Stackf.push(asinh(v1f));
              tmp_out.str("");
              tmp_out << "asinh(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oAtanh:
              Stackf.push(atanh(v1f));
              tmp_out.str("");
              tmp_out << "atanh(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oSqrt:
              Stackf.push(sqrt(v1f));
              tmp_out.str("");
              tmp_out << "sqrt(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            case oErf:
              Stackf.push(erf(v1f));
              tmp_out.str("");
              tmp_out << "erf(" << v1 << ")";
              Stack.push(tmp_out.str());
              break;
            default:
              ;
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
          v3f = Stackf.top();
          Stackf.pop();
          v2f = Stackf.top();
          Stackf.pop();
          v1f = Stackf.top();
          Stackf.pop();
          switch (op)
            {
              case oNormcdf:
                Stackf.push(0.5*(1+erf((v1f-v2f)/v3f/M_SQRT2)));
                tmp_out.str("");
                tmp_out << "normcdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
              case oNormpdf:
                Stackf.push(1/(v3f*sqrt(2*M_PI)*exp(pow((v1f-v2f)/v3f,2)/2)));
                tmp_out.str("");
                tmp_out << "normpdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
            }
          break;
        case FCUML:
        case FENDBLOCK:
        case FOK:
        case FENDEQU:
          go_on = false;
          break;
        default:
          ostringstream tmp;
          tmp << " in print_expression, unknown opcode " << it_code->first << "!! FENDEQU=" << FENDEQU << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
      it_code++;
    }
  return(tmp_out.str());
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
#ifdef DEBUG
  mexPrintf("compute_block_time\n");
#endif
  #ifndef DEBUG_EX
  if (evaluate && !steady_state)
    {
      jacob = mxGetPr(jacobian_block[block_num]);
      mexPrintf("jacobian_block[%d]=%x\n",block_num, jacobian_block[block_num]);
      jacob_other_endo = mxGetPr(jacobian_other_endo_block[block_num]);
      jacob_exo = mxGetPr(jacobian_exo_block[block_num]);
      jacob_exo_det = mxGetPr(jacobian_det_exo_block[block_num]);
    }
  #endif

  //feclearexcept (FE_ALL_EXCEPT);
  while (go_on)
    {
      //tmp_it_code = it_code;
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
              mexPrintf("FLDV Param[var=%d]\n",var);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV y[var=%d, lag=%d, it_=%d], y_size=%d evaluate=%d\n",var, lag, it_, y_size, evaluate);
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
              mexPrintf("FLDV x[var=%d, lag=%d, it_=%d], nb_row_x=%d evaluate=%d\n",var, lag, it_, nb_row_x, evaluate);
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
              mexPrintf("FLDSV Param[var=%d]\n",var);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV y[var=%d]\n",var);
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
              mexPrintf("FLDSV x[var=%d]\n",var);
              tmp_out << " x[" << var << "](" << x[var] << ")";
#endif
              Stack.push(x[var]);
              break;
            case eExogenousDet:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV xd[var=%d]\n",var);
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
          mexPrintf("FLDST T[%d]\n",var);
          tmp_out << " T[" << var << "](" << T[var] << ")";
#endif
          Stack.push(T[var]);
          break;
        case FLDU:
          //load u variable in the processor
          var = ((FLDU_ *) it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FLDU u[%d]\n",var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDSU:
          //load u variable in the processor
          var = ((FLDSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDSU u[%d]\n",var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDR:
          //load u variable in the processor
          var = ((FLDR_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDR r[%d]\n",var);
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
          mexPrintf("FLDC = %f\n",ll);
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
              mexPrintf("FSTPV params[%d]\n",var);
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
          var = ((FSTPST_ *) it_code->second)->get_pos();
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
          u[var] = Stack.top();
#ifdef DEBUG
          mexPrintf("FSTPU\n");
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
          if(var>=u_count_alloc || var<0)
            mexPrintf("Erreur var=%d\n",var);
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
          r[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" r[%d](%f)=%s\n", var, r[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPG:
          //store in derivative (g) variable from the processor
          var = ((FSTPG_ *) it_code->second)->get_pos();
          g1[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;

        case FSTPG3:
          //store in derivative (g) variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPG3\n");
          mexEvalString("drawnow;");
#endif
          rr = Stack.top();
          switch(EQN_type)
            {
            case FirstEndoDerivative:
              eq = ((FSTPG3_ *) it_code->second)->get_row();
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
              mexPrintf("jacob[%d(size=%d*pos_col=%d + eq=%d )]=%f\n",eq + size*pos_col, size, pos_col, eq, rr);
              jacob[eq + size*pos_col] = rr;
              break;
            case FirstOtherEndoDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
              mexPrintf("jacob_other_endo[%d(size=%d*pos_col=%d + eq=%d)]=%f\n",size*pos_col + eq, size, pos_col, eq, rr);
              jacob_other_endo[eq + size*pos_col] = rr;
              break;
            case FirstExoDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
              mexPrintf("jacob_exo[%d(size=%d*pos_col=%dr + eq=%d)]=%f\n",size*pos_col+eq, size, pos_col, eq, rr);
              jacob_exo[eq + size*pos_col] = rr;
              break;
            case FirstExodetDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
              mexPrintf("jacob_exo_det[%d(size=%d*pos_col=%dr + eq=%d)]=%f\n",size*pos_col+eq, size, pos_col, eq, rr);
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
              try
                {
                  tmp = divide(v1 , v2);
                }
              catch(FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n",fpeh.GetErrorMsg().c_str(),error_location(evaluate, steady_state).c_str());
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
              try
                {
                  tmp = pow1(v1, v2);
                }
              catch(FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n",fpeh.GetErrorMsg().c_str(),error_location(evaluate, steady_state).c_str());
                  go_on = false;
                }
              Stack.push(tmp);

#ifdef DEBUG
              tmp_out << " |" << v1 << "^" << v2 << "|";
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
            default:
              /*throw EvalException();*/
              ;
            }
          break;
        case FUNARY:
          op = ((FUNARY_ *) it_code->second)->get_op_type();
          v1 = Stack.top();
          Stack.pop();
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
              catch(FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n",fpeh.GetErrorMsg().c_str(),error_location(evaluate, steady_state).c_str());
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
              catch(FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n",fpeh.GetErrorMsg().c_str(),error_location(evaluate, steady_state).c_str());
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

# endif
              break;
            default:
              ;
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
                Stack.push(1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3,2)/2)));
#ifdef DEBUG
                tmp_out << " |normpdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
                break;
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
              mexPrintf("FJMPIFEVAL length=%d\n",((FJMPIFEVAL_ *) it_code->second)->get_pos());
              mexEvalString("drawnow;");
#endif
              it_code += ((FJMPIFEVAL_ *) it_code->second)->get_pos()/* - 1*/;
            }
          break;
        case FJMP:
#ifdef DEBUG
          mexPrintf("FJMP length=%d\n",((FJMP_ *) it_code->second)->get_pos());
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
  mexPrintf("==> end of compute_block_time Block = %d\n",block_num);
  mexEvalString("drawnow;");
#endif
}

void
Interpreter::evaluate_a_block(const int size, const int type, string bin_basename, bool steady_state, int block_num,
                              const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag,
                              const int Block_List_Max_Lead, const int u_count_int)
{
  it_code_type begining;
  switch (type)
    {
    case EVALUATE_FORWARD:
      if (steady_state)
        compute_block_time(0, true, block_num, size, steady_state);
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
            }
        }
      break;
    case SOLVE_FORWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          for (int j = 0; j < size; j++)
            y[Block_Contain[j].Variable] += r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_FORWARD_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false, stack_solve_algo, solve_algo);
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          for (int j = 0; j < size; j++)
            y[Block_Contain[j].Variable] += r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
            }
        }
      mxFree(r);
      break;
    case EVALUATE_BACKWARD:
      if (steady_state)
        compute_block_time(0, true, block_num, size, steady_state);
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
            }
        }
      break;
    case SOLVE_BACKWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num, size, steady_state);
          for (int j = 0; j < size; j++)
            y[Block_Contain[j].Variable] += r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
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
          for (int j = 0; j < size; j++)
            y[Block_Contain[j].Variable] += r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num, size, steady_state);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
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
          for (int j = 0; j < size; j++)
            y[it_*y_size+Block_Contain[j].Variable] += r[j];
        }
      mxFree(r);
      break;
    }
}

int
Interpreter::simulate_a_block(const int size, const int type, string file_name, string bin_basename, bool Gaussian_Elimination, bool steady_state, int block_num,
                              const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag, const int Block_List_Max_Lead, const int u_count_int)
{
  it_code_type begining;
  int i;
  bool cvg;
  int giter;
  bool result = true;
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
      mexPrintf("SOLVE_FORWARD_SIMPLE size=%d\n",size);
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
                  y[Block_Contain[0].Variable] += - divide(r[0],g1[0]);
                }
              catch(FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      \n",fpeh.GetErrorMsg().c_str());
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
                      y[Per_y_+Block_Contain[0].Variable] += - divide(r[0], g1[0]);
                    }
                  catch(FloatingPointExceptionHandling &fpeh)
                    {
                      mexPrintf("%s      \n",fpeh.GetErrorMsg().c_str());
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
                  y[Block_Contain[0].Variable] += - divide(r[0], g1[0]);
                }
              catch(FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      \n",fpeh.GetErrorMsg().c_str());
                  mexPrintf("      Singularity in block %d", block_num+1);
                }
              /*if (isinf(y[Block_Contain[0].Variable]))
                {
                  if (error_not_printed)
                    {
                      mexPrintf("--------------------------------------\n  Error: Divide by zero with %5.15f/%5.15f\nSingularity in block %d\n--------------------------------------\n", r[0], g1[0], block_num);
                      error_not_printed = false;
                    }
                  res1 = NAN;
                }*/
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
          for (it_ = periods+y_kmin; it_ > y_kmin; it_--)
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
                      y[Per_y_+Block_Contain[0].Variable] += - divide(r[0], g1[0]);
                    }
                  catch(FloatingPointExceptionHandling &fpeh)
                    {
                      mexPrintf("%s      \n",fpeh.GetErrorMsg().c_str());
                      mexPrintf("      Singularity in block %d", block_num+1);
                    }

                  /*if (isinf(y[Per_y_+Block_Contain[0].Variable]))
                    {
                      if (error_not_printed)
                        {
                          mexPrintf("--------------------------------------\n  Error: Divide by zero with %5.15f/%5.15f\nSingularity in block %d\n--------------------------------------\n", r[0], g1[0], block_num);
                          error_not_printed = false;
                        }
                      res1 = NAN;
                    }*/
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
                  Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, solve_algo);
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
              cvg = false;
              Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, solve_algo);
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
                      Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, solve_algo);
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
                  cvg = false;
                  Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, solve_algo);
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
                  Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, solve_algo);
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
              cvg = false;
              Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, solve_algo);
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
              for (it_ = periods+y_kmin; it_ > y_kmin; it_--)
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
                      Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, solve_algo);
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
              for (it_ = periods+y_kmin; it_ > y_kmin; it_--)
                {
                  it_code = begining;
                  Per_y_ = it_*y_size;
                  error_not_printed = true;
                  compute_block_time(0, false, block_num, size, steady_state);
                  cvg = false;
                  Simulate_Newton_One_Boundary(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, solve_algo);
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
      giter = 0;
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
              Simulate_Newton_Two_Boundaries(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter, minimal_solving_periods, stack_solve_algo, endo_name_length, P_endo_names);
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
          Simulate_Newton_Two_Boundaries(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter, minimal_solving_periods, stack_solve_algo, endo_name_length, P_endo_names);
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

bool
Interpreter::compute_blocks(string file_name, string bin_basename, bool steady_state, bool evaluate, bool block, int &nb_blocks)
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
  //The big loop on intructions
  Block_Count = -1;
  bool go_on = true;
  it_code = code_liste.begin();
  it_code_type Init_Code = it_code;
  while (go_on)
    {
      switch (it_code->first)
        {
        case FBEGINBLOCK:
          Block_Count++;
#ifdef DEBUG
          mexPrintf("FBEGINBLOCK %d\n", Block_Count+1);
#endif
          //it's a new block
          {
            FBEGINBLOCK_ *fb = (FBEGINBLOCK_ *) it_code->second;
            Block_Contain = fb->get_Block_Contain();
            it_code++;
            if (evaluate)
              {
                if (!steady_state)
                  {
                    jacobian_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_jacob(),mxREAL));
                    mexPrintf("at block = %d, mxCreateDoubleMatrix(%d, %d, mxREAL) jacobian_block.size()=%d\n",Block_Count, fb->get_size(), fb->get_nb_col_jacob(), sizeof(jacobian_block[Block_Count]));
                    jacobian_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_exo_size(),mxREAL));
                    mexPrintf("at block = %d, mxCreateDoubleMatrix(%d, %d, mxREAL) jacobian_exo_block.size()=%d fb->get_exo_size()=%d\n",Block_Count, fb->get_size(), fb->get_exo_size(), sizeof(jacobian_exo_block[Block_Count]), fb->get_exo_size());
                    jacobian_det_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_det_exo_size(),mxREAL));
                    mexPrintf("at block = %d, mxCreateDoubleMatrix(%d, %d, mxREAL) jacobian_det_exo_block.size()=%d\n",Block_Count, fb->get_size(), fb->get_det_exo_size(), sizeof(jacobian_det_exo_block[Block_Count]));
                    jacobian_other_endo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_other_endo_jacob(),mxREAL));
                    mexPrintf("at block = %d, mxCreateDoubleMatrix(%d, %d, mxREAL) jacobian_other_endo_block.size()=%d\n",Block_Count, fb->get_size(), fb->get_nb_col_other_endo_jacob(), sizeof(jacobian_other_endo_block[Block_Count]));
                  }
                evaluate_a_block(fb->get_size(), fb->get_type(), bin_basename, steady_state, Block_Count,
                                 fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
              }
            else
              {
                result = simulate_a_block(fb->get_size(), fb->get_type(), file_name, bin_basename, true, steady_state, Block_Count,
                                          fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
                if (result == ERROR_ON_EXIT)
                  return ERROR_ON_EXIT;
              }
            delete fb;
          }
          if (result == ERROR_ON_EXIT)
            go_on = false;
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
          it_code++;
          break;
        case FDIMST:
#ifdef DEBUG
          mexPrintf("FDIMST\n");
#endif
          var = ((FDIMST_ *) it_code->second)->get_size();
          if (T)
            mxFree(T);
          T = (double *) mxMalloc(var*sizeof(double));
          it_code++;
          break;
        default:
          ostringstream tmp;
          tmp << " in compute_blocks, unknown command " << it_code->first << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  mxFree(Init_Code->second);
  nb_blocks = Block_Count+1;
  if (T)
    mxFree(T);
  return result;
}
