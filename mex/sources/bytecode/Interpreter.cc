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
                         string &filename_arg, int minimal_solving_periods_arg)
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
}

string
Interpreter::remove_white(string str)
{
  string temp;
  for (unsigned int i = 0; i < str.length(); i++)
    if (str[i] != ' ')
      temp += str[i];
  return(temp);
}

string
Interpreter::add_underscore_to_fpe(string str)
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
  return(temp + "\n " + tmp_n );
}


string
Interpreter::get_variable(SymbolType variable_type, int variable_num)
{
  ostringstream res;
#ifndef DEBUG_EX
  mxArray *M_;
  char * P;
  unsigned int R, C;
  M_ = mexGetVariable("global", "M_");
  switch(variable_type)
    {
    case eEndogenous:
      C = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
      R = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
      P = (char*) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
      if (variable_num < R)
        for (unsigned int i = 0; i < C; i++)
          res << P[2*(variable_num+i*R)];
      else
        mexPrintf("=> Unknown endogenous variable n° %d",EQN_dvar1);
      break;
    case eExogenous:
    case eExogenousDet:
      C = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
      R = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
      P = (char*) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
      if (variable_num < R)
        for (unsigned int i = 0; i < C; i++)
          res << P[2*(variable_num+i*R)];
      else
        mexPrintf("=> Unknown exogenous variable n° %d",EQN_dvar1);
      break;
    case eParameter:
      C = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
      R = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
      P = (char*) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
      if (variable_num < R)
        for (unsigned int i = 0; i < C; i++)
          res << P[2*(variable_num+i*R)];
      else
        mexPrintf("=> Unknown parameter n° %d",EQN_dvar1);
      break;
    default:
      break;
    }
#endif
  return(remove_white(res.str()));
}


string
Interpreter::error_location(bool evaluate)
{
  string tmp;
  stringstream Error_loc("");
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
          Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to endogenous variable ";
        else
          Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to endogenous variable " << get_variable(eEndogenous, EQN_dvar1);
        break;
      default:
        return("???");
        break;
    }
  Error_loc << "\n " << add_underscore_to_fpe(print_expression(it_code_expr, evaluate));
  return(Error_loc.str());
}

double
Interpreter::pow1(double a, double b, bool evaluate)
{
  double r = pow_(a, b);
  if (isnan(r))
    {
      if (error_not_printed)
        {
          mexPrintf("--------------------------------------\nError: X^a with X=%5.15f\n in %s \n--------------------------------------\n", a,error_location(evaluate).c_str());
          error_not_printed = false;
          r = 0.0000000000000000000000001;
        }
      res1 = NAN;
      return r;
    }
  return r;
}

double
Interpreter::log1(double a, bool evaluate)
{
  double r = log(a);
  if (isnan(r))
    {
      if (error_not_printed)
        {
          mexPrintf("--------------------------------------\nError: log(X) with X=%5.15f\n in %s \n--------------------------------------\n", a,error_location(evaluate).c_str());
          error_not_printed = false;
          r = 0.0000000000000000000000001;
        }
      res1 = NAN;
      return r;
    }
  return r;
}

double
Interpreter::log10_1(double a, bool evaluate)
{
  double r = log(a);
  if (isnan(r))
    {
      if (error_not_printed)
        {
          mexPrintf("--------------------------------------\nError: log10(X) with X=%5.15f\n in %s \n--------------------------------------\n", a,error_location(evaluate).c_str());
          error_not_printed = false;
          r = 0.0000000000000000000000001;
        }
      res1 = NAN;
      return r;
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
  string v1, v2;
  double v1f, v2f;
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
            }
          break;
        case FLDV:
          //load a variable in the processor
          switch (((FLDV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDV_ *) it_code->second)->get_pos();
              tmp_out.str("");
              tmp_out << get_variable(eParameter, var);
              Stack.push(tmp_out.str());
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
              tmp_out.str("");
              tmp_out << get_variable(eParameter, var);
              Stack.push(tmp_out.str());
              Stackf.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              tmp_out.str("");
              tmp_out << get_variable(eEndogenous, var);
              Stack.push(tmp_out.str());
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
              tmp_out.str("");
              tmp_out << get_variable(eExogenous, var);
              Stack.push(tmp_out.str());
              break;
              Stackf.push(x[var]);
            case eExogenousDet:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              tmp_out.str("");
              tmp_out << get_variable(eExogenousDet, var);
              Stack.push(tmp_out.str());
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
            default:
              ;
            }
          break;
        case FCUML:
        case FENDBLOCK:
        case FOK:
        case FENDEQU:
          go_on = false;
          break;
        default:
          mexPrintf("Unknown opcode %d!! FENDEQU=%d\n", it_code->first, FENDEQU);
          mexErrMsgTxt("End of bytecode");
          break;
        }
      it_code++;
    }
  return(tmp_out.str());
}

void
Interpreter::compute_block_time(int Per_u_, bool evaluate, int block_num)
{
  int var, lag = 0, op;
  ostringstream tmp_out;
  double v1, v2;
  bool go_on = true;
  double ll;

  EQN_block = block_num;
  //feclearexcept (FE_ALL_EXCEPT);
  while (go_on)
    {
      //tmp_it_code = it_code;
      switch (it_code->first)
        {
        case FNUMEXPR:
          it_code_expr = it_code;
          switch (((FNUMEXPR_ *) it_code->second)->get_expression_type())
            {
            case TemporaryTerm:
              EQN_type = TemporaryTerm;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              break;
            case ModelEquation:
              EQN_type = ModelEquation;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              break;
            case FirstEndoDerivative:
              EQN_type = FirstEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExoDerivative:
              EQN_type = FirstExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExodetDerivative:
              EQN_type = FirstExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstParamDerivative:
              EQN_type = FirstParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              break;
            case SecondEndoDerivative:
              EQN_type = SecondEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExoDerivative:
              EQN_type = SecondExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExodetDerivative:
              EQN_type = SecondExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondParamDerivative:
              EQN_type = SecondParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              break;
            case ThirdEndoDerivative:
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
              Stack.push(params[var]);
#ifdef DEBUG
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              break;
            case eEndogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
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
              Stack.push(x[it_+lag+var*nb_row_x]);
#ifdef DEBUG
              tmp_out << " x[" << it_+lag << ", " << var << "](" << x[it_+lag+var*nb_row_x] << ")";
#endif
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
              Stack.push(params[var]);
#ifdef DEBUG
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              break;
            case eEndogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              if (evaluate)
                Stack.push(ya[var]);
              else
                Stack.push(y[var]);
#ifdef DEBUG
              tmp_out << " y[" << var << "](" << y[var] << ")";
              if(var<0 || var>= y_size)
                {
                  mexPrintf("y[%d]=",var);
                  mexErrMsgTxt("End of bytecode");
                }
#endif
              break;
            case eExogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
              Stack.push(x[var]);
#ifdef DEBUG
              tmp_out << " x[" << var << "](" << x[var] << ")";
#endif
              break;
            case eExogenousDet:
              var = ((FLDSV_ *) it_code->second)->get_pos();
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
              mexPrintf("params[%d]=%f\n", var, params[var]);
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf(" steady_y[%d]=%f\n", var, steady_y[var]);
#endif
              Stack.push(steady_y[var]);
              break;
            case eExogenous:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf(" x[%d] = %f\n", var, x[var]);
#endif
              Stack.push(x[var]);
              break;
            case eExogenousDet:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf(" xd[%d] = %f\n", var, x[var]);
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
          tmp_out << " T[" << var << "](" << T[var] << ")";
#endif
          Stack.push(T[var]);
          break;
        case FLDU:
          //load u variable in the processor
          var = ((FLDU_ *) it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDSU:
          //load u variable in the processor
          var = ((FLDSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDR:
          //load u variable in the processor
          var = ((FLDR_ *) it_code->second)->get_pos();
          Stack.push(r[var]);
          break;
        case FLDZ:
          //load 0 in the processor
          Stack.push(0.0);
#ifdef DEBUG
          tmp_out << " 0";
#endif
          break;
        case FLDC:
          //load a numerical constant in the processor
          ll = ((FLDC_ *) it_code->second)->get_value();
#ifdef DEBUG
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
              double r;
              r = v1 / v2;
              if (isinf(r))
                {
                  if (error_not_printed)
                    {
                      mexPrintf("--------------------------------------\n  Error: Divide by zero with %5.15f/%5.15f\n in %s \n--------------------------------------\n", v1, v2,error_location(evaluate).c_str());
                      error_not_printed = false;
                      r = 1e70;
                    }
                  res1 = NAN;
                  go_on = false;
                }
              Stack.push(r);
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
              r = pow1(v1, v2, evaluate);
              if (isnan(res1))
                go_on = false;
              Stack.push(r);

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
              Stack.push(log1(v1, evaluate));
              if (isnan(res1))
                go_on = false;

#ifdef DEBUG
              tmp_out << " |log(" << v1 << ")|";
#endif
              break;
            case oLog10:
              Stack.push(log10_1(v1, evaluate));
              if (isnan(res1))
                go_on = false;
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
            default:
              ;
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
          go_on = false;
          break;
        case FENDEQU:
          break;
        case FOK:
          op = ((FOK_ *) it_code->second)->get_arg();
          if (Stack.size() > 0)
            {
              mexPrintf("error: Stack not empty!\n");
              mexErrMsgTxt("End of bytecode");
            }
          break;
        default:
          mexPrintf("Unknown opcode %d!! FENDEQU=%d\n", it_code->first, FENDEQU);
          mexErrMsgTxt("End of bytecode");
          break;
        }
      it_code++;
    }
}

void
Interpreter::evaluate_a_block(const int size, const int type, string bin_basename, bool steady_state, int block_num,
                              const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag, const int Block_List_Max_Lead, const int u_count_int)
{
  it_code_type begining;
  switch (type)
    {
    case EVALUATE_FORWARD:
      if (steady_state)
        compute_block_time(0, true, block_num);
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num);
            }
        }
      break;
    case SOLVE_FORWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num);
          for (int j = 0; j < size; j++)
            y[Block_Contain[j].Variable] += r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_FORWARD_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false);
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num);
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
              compute_block_time(0, true, block_num);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
            }
        }
      mxFree(r);
      break;
    case EVALUATE_BACKWARD:
      if (steady_state)
        compute_block_time(0, true, block_num);
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, block_num);
            }
        }
      break;
    case SOLVE_BACKWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num);
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
              compute_block_time(0, true, block_num);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_BACKWARD_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false);
      r = (double *) mxMalloc(size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, block_num);
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
              compute_block_time(0, true, block_num);
              for (int j = 0; j < size; j++)
                y[it_*y_size+Block_Contain[j].Variable] += r[j];
            }
        }
      mxFree(r);
      break;
    case SOLVE_TWO_BOUNDARIES_SIMPLE:
    case SOLVE_TWO_BOUNDARIES_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, periods, y_kmin, y_kmax, steady_state, true);
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = (double *) mxMalloc(size*sizeof(double));
      begining = it_code;
      for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
        {
          Per_u_ = (it_-y_kmin)*u_count_int;
          Per_y_ = it_*y_size;
          it_code = begining;
          compute_block_time(Per_u_, true, block_num);
          for (int j = 0; j < size; j++)
            y[it_*y_size+Block_Contain[j].Variable] += r[j];
        }
      mxFree(r);
      break;
    }
}

bool
Interpreter::simulate_a_block(const int size, const int type, string file_name, string bin_basename, bool Gaussian_Elimination, bool steady_state, int block_num,
                              const bool is_linear, const int symbol_table_endo_nbr, const int Block_List_Max_Lag, const int Block_List_Max_Lead, const int u_count_int)
{
  it_code_type begining;
  int i;
  bool cvg;
  int giter;
  bool result = true;
  double *y_save;
  switch (type)
    {
    case EVALUATE_FORWARD:
      if (steady_state)
        compute_block_time(0, false, block_num);
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, false, block_num);
            }
        }
      break;
    case EVALUATE_BACKWARD:
      if (steady_state)
        compute_block_time(0, false, block_num);
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, false, block_num);
            }
        }
      break;
    case SOLVE_FORWARD_SIMPLE:
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
              compute_block_time(0, false, block_num);
              double rr;
              rr = r[0];
              cvg = (fabs(rr) < solve_tolf);
              if (cvg)
                continue;
              y[Block_Contain[0].Variable] += -r[0]/g1[0];
              if (isinf(y[Block_Contain[0].Variable]))
                {
                  if (error_not_printed)
                    {
                      mexPrintf("--------------------------------------\n  Error: Divide by zero with %5.15f/%5.15f\nSingularity in block %d\n--------------------------------------\n", r[0], g1[0], block_num);
                      error_not_printed = false;
                    }
                  res1 = NAN;
                }
              iter++;
            }
          if (!cvg)
            {
              mexPrintf("Convergence not achieved in block %d, after %d iterations\n", Block_Count, iter);
              mexPrintf("r[0]= %f\n", r[0]);
              return false;
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
                  compute_block_time(0, false, block_num);
                  double rr;
                  if (fabs(1+y[Per_y_+Block_Contain[0].Variable]) > eps)
                    rr = r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
                  else
                    rr = r[0];
                  cvg = (fabs(rr) < solve_tolf);
                  if (cvg)
                    continue;
                  y[Per_y_+Block_Contain[0].Variable] += -r[0]/g1[0];
                  if (isinf(y[Per_y_+Block_Contain[0].Variable]))
                    {
                      if (error_not_printed)
                        {
                          mexPrintf("--------------------------------------\n  Error: Divide by zero with %5.15f/%5.15f\nSingularity in block %d\n--------------------------------------\n", r[0], g1[0], block_num);
                          error_not_printed = false;
                        }
                      res1 = NAN;
                    }
                  iter++;
                }
              if (!cvg)
                {
                  mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                  mexErrMsgTxt("End of bytecode");
                }
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_BACKWARD_SIMPLE:
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
              compute_block_time(0, false, block_num);
              double rr;
              rr = r[0];
              cvg = (fabs(rr) < solve_tolf);
              if (cvg)
                continue;
              y[Block_Contain[0].Variable] += -r[0]/g1[0];
              if (isinf(y[Block_Contain[0].Variable]))
                {
                  if (error_not_printed)
                    {
                      mexPrintf("--------------------------------------\n  Error: Divide by zero with %5.15f/%5.15f\nSingularity in block %d\n--------------------------------------\n", r[0], g1[0], block_num);
                      error_not_printed = false;
                    }
                  res1 = NAN;
                }
              iter++;
            }
          if (!cvg)
            {
              mexPrintf("Convergence not achieved in block %d, after %d iterations\n", Block_Count, iter);
              return false;
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
                  compute_block_time(0, false, block_num);
                  double rr;
                  if (fabs(1+y[Per_y_+Block_Contain[0].Variable]) > eps)
                    rr = r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
                  else
                    rr = r[0];
                  cvg = (fabs(rr) < solve_tolf);
                  if (cvg)
                    continue;
                  y[Per_y_+Block_Contain[0].Variable] += -r[0]/g1[0];
                  if (isinf(y[Per_y_+Block_Contain[0].Variable]))
                    {
                      if (error_not_printed)
                        {
                          mexPrintf("--------------------------------------\n  Error: Divide by zero with %5.15f/%5.15f\nSingularity in block %d\n--------------------------------------\n", r[0], g1[0], block_num);
                          error_not_printed = false;
                        }
                      res1 = NAN;
                    }
                  iter++;
                }
              if (!cvg)
                {
                  mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                  mexErrMsgTxt("End of bytecode");
                }
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_FORWARD_COMPLETE:
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false);
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
                  compute_block_time(0, false, block_num);
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
                  result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, EQN_block_number);
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
                  mexPrintf("Convergence not achieved in block %d, after %d iterations\n", Block_Count, iter);
                  return false;
                }
            }
          else
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              iter = 0;
              res1 = res2 = max_res = 0; max_res_idx = 0;
              error_not_printed = true;
              compute_block_time(0, false, block_num);
              cvg = false;
              result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, EQN_block_number);
              if (!result)
                {
                  mexPrintf("Convergence not achieved in block %d\n", Block_Count);
                  return false;
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
                      compute_block_time(0, false, block_num);
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
                      result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, EQN_block_number);
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
                      mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                      mexErrMsgTxt("End of bytecode");
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
                  res1 = res2 = max_res = 0; max_res_idx = 0;
                  error_not_printed = true;
                  compute_block_time(0, false, block_num);
                  cvg = false;
                  result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, EQN_block_number);
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
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state, false);
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
                  compute_block_time(0, false, block_num);
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
                  result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, EQN_block_number);
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
                  mexPrintf("Convergence not achieved in block %d, after %d iterations\n", Block_Count, iter);
                  return false;
                }
            }
          else
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              iter = 0;
              res1 = res2 = max_res = 0; max_res_idx = 0;
              error_not_printed = true;
              compute_block_time(0, false, block_num);
              cvg = false;
              result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true, EQN_block_number);
              if (!result)
                {
                  mexPrintf("Convergence not achieved in block %d\n", Block_Count);
                  return false;
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
                      compute_block_time(0, false, block_num);
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
                      result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, EQN_block_number);
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
                      mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                      mexErrMsgTxt("End of bytecode");
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
                  compute_block_time(0, false, block_num);
                  cvg = false;
                  result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false, EQN_block_number);
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
      if (steady_state)
        {
          mexPrintf("SOLVE_TWO_BOUNDARIES in a steady state model: impossible case\n");
          return false;
        }
      fixe_u(&u, u_count_int, u_count_int);
      Read_SparseMatrix(bin_basename, size, periods, y_kmin, y_kmax, steady_state, true);
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = (double *) mxMalloc(size*sizeof(double));
      y_save = (double *) mxMalloc(y_size*sizeof(double)*(periods+y_kmax+y_kmin));
      begining = it_code;
      if (!Gaussian_Elimination)
        {
        }
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
                  compute_block_time(Per_u_, false, block_num);
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
              simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter, minimal_solving_periods, EQN_block_number);
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
              mexPrintf("Convergence not achieved in block %d, after %d iterations\n", Block_Count, iter);
              mexErrMsgTxt("End of bytecode");
            }
        }
      else
        {
          res1 = res2 = max_res = 0; max_res_idx = 0;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              Per_u_ = (it_-y_kmin)*u_count_int;
              Per_y_ = it_*y_size;
              it_code = begining;
              compute_block_time(Per_u_, false, block_num);
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
          simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter, minimal_solving_periods, EQN_block_number);
        }
      mxFree(r);
      mxFree(y_save);
      mxFree(u);
      mxFree(index_vara);
      mxFree(index_equa);
      memset(direction, 0, size_of_direction);
      break;
    default:
      mexPrintf("Unknown type =%d\n", type);
      mexEvalString("drawnow;");
      mexErrMsgTxt("End of bytecode");
    }
  return true;
}

bool
Interpreter::compute_blocks(string file_name, string bin_basename, bool steady_state, bool evaluate)
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
      mexPrintf("%s.cod Cannot be opened\n", file_name.c_str());
      mexEvalString("drawnow;");
      filename += " stopped";
      mexEvalString("drawnow;");
      mexErrMsgTxt(filename.c_str());
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
              evaluate_a_block(fb->get_size(), fb->get_type(), bin_basename, steady_state, Block_Count,
                               fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
            else
              result = simulate_a_block(fb->get_size(), fb->get_type(), file_name, bin_basename, true, steady_state, Block_Count,
                                        fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
            delete fb;
          }
          if (!result)
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
          mexPrintf("Unknown command \n");
          mexEvalString("drawnow;");
          mexErrMsgTxt("End of bytecode");
          break;
        }
    }
  mxFree(Init_Code->second);
  if (T)
    mxFree(T);
  return result;
}
