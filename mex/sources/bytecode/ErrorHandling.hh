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

#ifndef ERROR_HANDLING
#define ERROR_HANDLING

#include <cstring>
#include <iostream>
#include <sstream>
#include "CodeInterpreter.hh"

using namespace std;

const int NO_ERROR_ON_EXIT = 0;
const int ERROR_ON_EXIT = 1;

typedef vector<pair<Tags, void * > > code_liste_type;
typedef code_liste_type::const_iterator it_code_type;

class GeneralExceptionHandling
{
  string ErrorMsg;
public:
  GeneralExceptionHandling(string ErrorMsg_arg) : ErrorMsg(ErrorMsg_arg){};
  inline string
  GetErrorMsg()
    {
      return ErrorMsg;
    }
  inline void
  completeErrorMsg(string ErrorMsg_arg)
    {
      ErrorMsg += ErrorMsg_arg;
    }
};

class FloatingPointExceptionHandling : public GeneralExceptionHandling
{
public:
  FloatingPointExceptionHandling(string value) : GeneralExceptionHandling(string("Floating point error in bytecode: " + value))
   {};
};


class LogExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  LogExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log(X)"),
                                           value(value_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value << "\n";
     completeErrorMsg(tmp.str());
   };
};

class Log10ExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  Log10ExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log10(X)"),
                                           value(value_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value << "\n";
     completeErrorMsg(tmp.str());
   };
};

class DivideExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  DivideExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("a/X"),
                                   value1(value1_arg),
                                   value2(value2_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value2 << "\n";
     completeErrorMsg(tmp.str());
   };
};


class PowExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  PowExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("X^a"),
                                   value1(value1_arg),
                                   value2(value2_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value1 << "\n";
     completeErrorMsg(tmp.str());
   };
};


class FatalExceptionHandling : public GeneralExceptionHandling
{
public:
  FatalExceptionHandling(string ErrorMsg_arg) : GeneralExceptionHandling("Fatal error in bytecode:")
    {
      completeErrorMsg(ErrorMsg_arg);
    };
  FatalExceptionHandling() : GeneralExceptionHandling("")
    {
    };
};



class ErrorMsg
{
public:
double *T;
int nb_row_xd, nb_row_x, y_size;
int y_kmin, y_kmax, periods;
double *x, *params;
double *u, *y, *ya;
double *steady_y, *steady_x;
double *g2, *g1, *r;
vector<mxArray*> jacobian_block, jacobian_other_endo_block, jacobian_exo_block, jacobian_det_exo_block;
map<unsigned int,double> TEF;
map<pair<unsigned int, unsigned int>, double > TEFD;
map<pair<unsigned int, pair<unsigned int, unsigned int> >, double > TEFDD;

ExpressionType EQN_type;
it_code_type it_code_expr;
unsigned int nb_endo, nb_exo, nb_param;
char *P_endo_names, *P_exo_names, *P_param_names;
unsigned int endo_name_length, exo_name_length, param_name_length;
unsigned int EQN_equation, EQN_block, EQN_block_number;
unsigned int EQN_dvar1, EQN_dvar2, EQN_dvar3;


inline ErrorMsg()
{
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

inline string
add_underscore_to_fpe(const string &str)
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


inline string
get_variable(const SymbolType variable_type, const unsigned int variable_num)
{
  ostringstream res;
  switch(variable_type)
    {
    case eEndogenous:
      if (variable_num < nb_endo)
        {
          for (unsigned int i = 0; i < endo_name_length; i++)
            if (P_endo_names[CHAR_LENGTH*(variable_num+i*nb_endo)] != ' ')
              res << P_endo_names[CHAR_LENGTH*(variable_num+i*nb_endo)];
        }
      else
        mexPrintf("=> Unknown endogenous variable n° %d",variable_num);
      break;
    case eExogenous:
    case eExogenousDet:
      if (variable_num < nb_exo)
        {
          for (unsigned int i = 0; i < exo_name_length; i++)
            if (P_exo_names[CHAR_LENGTH*(variable_num+i*nb_exo)] != ' ')
              res << P_exo_names[CHAR_LENGTH*(variable_num+i*nb_exo)];
        }
      else
        mexPrintf("=> Unknown exogenous variable n° %d",variable_num);
      break;
    case eParameter:
      if (variable_num < nb_param)
        {
          for (unsigned int i = 0; i < param_name_length; i++)
            if (P_param_names[CHAR_LENGTH*(variable_num+i*nb_param)] != ' ')
              res << P_param_names[CHAR_LENGTH*(variable_num+i*nb_param)];
        }
      else
        mexPrintf("=> Unknown parameter n° %d",variable_num);
      break;
    default:
      break;
    }
  return(res.str());
}


inline string
error_location(bool evaluate, bool steady_state, int size, int block_num, int it_, int Per_u_)
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
  Error_loc << endl << add_underscore_to_fpe("      " + print_expression(it_code_expr, evaluate, size, block_num, steady_state, Per_u_, it_));
  return(Error_loc.str());
}

inline string
print_expression(it_code_type it_code, bool evaluate, int size, int block_num, bool steady_state, int Per_u_, int it_)
{
  int var, lag = 0, op, eq;
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
  double *jacob = NULL, *jacob_other_endo = NULL, *jacob_exo = NULL, *jacob_exo_det = NULL;
  if (evaluate)
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
#ifdef DEBUG
              mexPrintf("FLDV_ Param var=%d",var);
#endif
              Stack.push(get_variable(eParameter, var));
              Stackf.push(params[var]);
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case eEndogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV_ endo var=%d, lag=%d",var, lag);
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case eExogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV_ exo var=%d, lag=%d",var, lag);
#endif
              tmp_out.str("");
              if(lag != 0)
                tmp_out << get_variable(eExogenous, var) << "(" << lag << ")";
              else
                tmp_out << get_variable(eExogenous, var);
              Stack.push(tmp_out.str());
              Stackf.push(x[it_+lag+var*nb_row_x]);
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case eExogenousDet:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV_ exo_det var=%d, lag=%d",var, lag);
#endif
              tmp_out.str("");
              if(lag != 0)
                tmp_out << get_variable(eExogenousDet, var) << "(" << lag << ")";
              else
                tmp_out << get_variable(eExogenousDet, var);
              Stack.push(tmp_out.str());
              Stackf.push(x[it_+lag+var*nb_row_xd]);
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
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
#ifdef DEBUG
              mexPrintf("FLDSV_ param var=%d",var);
#endif
              Stack.push(get_variable(eParameter, var));
              Stackf.push(params[var]);
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case eEndogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV_ endo var=%d",var);
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case eExogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV_ exo var=%d",var);
#endif
              Stack.push(get_variable(eExogenous, var));
              Stackf.push(x[var]);
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case eExogenousDet:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV_ exo_det var=%d",var);
#endif
              Stack.push(get_variable(eExogenousDet, var));
              Stackf.push(x[var]);
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
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
        case FSTPG2:
          go_on = false;
          //store in derivative (g) variable from the processor
          eq = ((FSTPG2_ *) it_code->second)->get_row();
          var = ((FSTPG2_ *) it_code->second)->get_col();
          tmp_out.str("");
          tmp_out << "jacob[" << eq+size*var+1 << "] = " << Stack.top();
          Stack.pop();
          jacob[eq + size*var] = Stackf.top();
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
#ifdef DEBUG
              mexPrintf("+");
#endif
              Stackf.push(v1f + v2f);
              tmp_out.str("");
              tmp_out << v1 << " + " << v2;
              Stack.push(tmp_out.str());
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oMinus:
#ifdef DEBUG
              mexPrintf("-");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oTimes:
#ifdef DEBUG
              mexPrintf("*");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oDivide:
#ifdef DEBUG
              mexPrintf("/");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oLess:
#ifdef DEBUG
              mexPrintf("<");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oGreater:
#ifdef DEBUG
              mexPrintf(">");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oLessEqual:
#ifdef DEBUG
              mexPrintf("<=");
#endif
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
 #ifdef DEBUG
              mexPrintf("ok\n");
#endif
             break;
            case oGreaterEqual:
#ifdef DEBUG
              mexPrintf(">=");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oEqualEqual:
#ifdef DEBUG
              mexPrintf("==");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
             break;
            case oDifferent:
#ifdef DEBUG
              mexPrintf("!=");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oPower:
#ifdef DEBUG
              mexPrintf("^");
#endif
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
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oMax:
#ifdef DEBUG
              mexPrintf("max");
#endif
              Stackf.push(max(v1f, v2f));
              tmp_out.str("");
              tmp_out << "max(" << v1 << ", " << v2 << ")";
              Stack.push(tmp_out.str());
#ifdef DEBUG
              mexPrintf("ok\n");
#endif
              break;
            case oMin:
#ifdef DEBUG
              mexPrintf("min");
#endif
              Stackf.push(min(v1f, v2f));
              tmp_out.str("");
              tmp_out << "min(" << v1 << ", " << v2 << ")";
              Stack.push(tmp_out.str());
#ifdef DEBUG
              mexPrintf("ok\n");
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

};

#endif
