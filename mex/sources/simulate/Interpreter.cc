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


Interpreter::Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *direction_arg, int y_size_arg,
                         int nb_row_x_arg, int nb_row_xd_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
                         int maxit_arg_, double solve_tolf_arg, int size_of_direction_arg, double slowc_arg, int y_decal_arg, double markowitz_c_arg,
                         string &filename_arg)
{
  params=params_arg;
  y=y_arg;
  ya=ya_arg;
  x=x_arg;
  direction=direction_arg;
  y_size=y_size_arg;
  nb_row_x=nb_row_x_arg;
  nb_row_xd=nb_row_xd_arg;
  periods=periods_arg;
  y_kmax=y_kmax_arg;
  y_kmin=y_kmin_arg;
  maxit_=maxit_arg_;
  solve_tolf=solve_tolf_arg;
  size_of_direction=size_of_direction_arg;
  slowc=slowc_arg;
  slowc_save = slowc;
  y_decal=y_decal_arg;
  markowitz_c=markowitz_c_arg;
  filename=filename_arg;
  T=NULL;
  error_not_printed = true;
}

double
Interpreter::pow1(double a, double b)
{
	/*double r;
	if(a>=0)
    r=pow_(a,b);
	else
	  {
	     //r=0;
	     //max_res=res1=res2=BIG;
	     if(error_not_printed)
	       {
	         mexPrintf("Error: X^a with X<0\n");
	         error_not_printed = false;
	       }
	     //r = BIG;
	     //r = -pow_(-a, b);
	     //r = 0;
	     //r = SMALL;
	     //r = pow_(-a, b);
	  }*/
	double r = pow_(a, b);
  if (isnan(r) || isinf(r))
    {
    	if(a<0 && error_not_printed)
	       {
	         mexPrintf("Error: X^a with X=%5.25f\n",a);
	         error_not_printed = false;
	         r = 0.0000000000000000000000001;
	       }
      //res1=NAN;
      return(r);
    }
  else
    return r;
}

double
Interpreter::log1(double a)
{
	/*double r;
	if(a>=0)
    r=pow_(a,b);
	else
	  {
	     //r=0;
	     //max_res=res1=res2=BIG;
	     if(error_not_printed)
	       {
	         mexPrintf("Error: X^a with X<0\n");
	         error_not_printed = false;
	       }
	     //r = BIG;
	     //r = -pow_(-a, b);
	     //r = 0;
	     //r = SMALL;
	     //r = pow_(-a, b);
	  }*/
	double r = log(a);
  if (isnan(r) || isinf(r))
    {
    	if(a<=0 && error_not_printed)
	       {
	         mexPrintf("Error: log(X) with X<=0\n");
	         error_not_printed = false;
	       }
      res1=NAN;
      return(r);
    }
  else
    return r;
}



void
Interpreter::compute_block_time(int Per_u_) /*throw(EvalException)*/
{
  int var, lag, op;
  ostringstream tmp_out;
  double v1, v2;
  char/*uint8_t*/ cc;
  bool go_on=true;
  double *ll;
  while (go_on)
    {
      switch (cc=get_code_char)
        {
          case FLDV :
            //load a variable in the processor
            switch (get_code_char)
              {
                case eParameter :
                  var=get_code_int;
                  Stack.push(params[var]);
#ifdef DEBUG
                  tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
                  break;
                case eEndogenous :
                  var=get_code_int;
                  lag=get_code_int;
                  Stack.push(y[(it_+lag)*y_size+var]);
#ifdef DEBUG
                  tmp_out << " y[" << it_+lag << ", " << var << "](" << y[(it_+lag)*y_size+var] << ")";
#endif
                  break;
                case eExogenous :
                  var=get_code_int;
                  lag=get_code_int;
                  Stack.push(x[it_+lag+var*nb_row_x]);
#ifdef DEBUG
                  tmp_out << " x[" << it_+lag << ", " << var << "](" << x[it_+lag+var*nb_row_x] << ")";
#endif
                  break;
                case eExogenousDet :
                  var=get_code_int;
                  lag=get_code_int;
                  Stack.push(x[it_+lag+var*nb_row_xd]);
                  break;
                default:
                  mexPrintf("Unknown variable type\n");
              }
            break;
          case FLDSV :
            //load a variable in the processor
            switch (get_code_char)
              {
                case eParameter :
                  var=get_code_int;
                  Stack.push(params[var]);
#ifdef DEBUG
                  tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
                  break;
                case eEndogenous :
                  var=get_code_int;
                  Stack.push(y[var]);
#ifdef DEBUG
                  tmp_out << " y[" << var << "](" << y[var] << ")";
#endif
                  break;
                case eExogenous :
                  var=get_code_int;
                  Stack.push(x[var]);
#ifdef DEBUG
                  tmp_out << " x[" << var << "](" << x[var] << ")";
#endif
                  break;
                case eExogenousDet :
                  var=get_code_int;
                  Stack.push(x[var]);
                  break;
                default:
                  mexPrintf("Unknown variable type\n");
              }
            break;
          case FLDT :
            //load a temporary variable in the processor
            var=get_code_int;
#ifdef DEBUG
            tmp_out << " T[" << it_ << ", " << var << "](" << T[var*(periods+y_kmin+y_kmax)+it_] << ")";
#endif
            Stack.push(T[var*(periods+y_kmin+y_kmax)+it_]);
            break;
          case FLDST :
            //load a temporary variable in the processor
            var=get_code_int;
#ifdef DEBUG
            tmp_out << " T[" << var << "](" << T[var] << ")";
#endif
            Stack.push(T[var]);
            break;
          case FLDU :
            //load u variable in the processor
            var=get_code_int;
            var+=Per_u_;
#ifdef DEBUG
            tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
            Stack.push(u[var]);
            break;
          case FLDSU :
            //load u variable in the processor
            var=get_code_int;
#ifdef DEBUG
            tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
            Stack.push(u[var]);
            break;
          case FLDR :
            //load u variable in the processor
            var=get_code_int;
            Stack.push(r[var]);
            break;
          case FLDZ :
            //load 0 in the processor
            Stack.push(0);
#ifdef DEBUG
            tmp_out << " 0";
#endif
            break;
          case FLDC :
            //load a numerical constant in the processor
            /*asm("fld\n\t"
                "fstp %%st" : "=t" (ll) : "0" ((double)(*Code)));*/
            ll=get_code_pdouble;
#ifdef DEBUG
            tmp_out << " " << *ll;
#endif

            Stack.push(*ll);
            break;
          case FSTPV :
            //load a variable in the processor
            switch (get_code_char)
              {
                case eParameter :
                  var=get_code_int;
                  params[var] = Stack.top();
                  Stack.pop();
                  break;
                case eEndogenous :
                  var=get_code_int;
                  lag=get_code_int;
                  y[(it_+lag)*y_size+var] = Stack.top();
#ifdef DEBUG
                  tmp_out << "=>";
                  //mexPrintf(" y[%d, %d](%f)=%s\n", it_+lag, var, y[(it_+lag)*y_size+var], tmp_out.str().c_str());
                  tmp_out.str("");
#endif
                  Stack.pop();
                  break;
                case eExogenous :
                  var=get_code_int;
                  lag=get_code_int;
                  x[it_+lag+var*nb_row_x]  = Stack.top();
                  Stack.pop();
                  break;
                case eExogenousDet :
                  var=get_code_int;
                  lag=get_code_int;
                  x[it_+lag+var*nb_row_xd] = Stack.top();
                  Stack.pop();
                  break;
                default:
                  mexPrintf("Unknown vraibale type\n");
              }
            break;
					case FSTPSV :
            //load a variable in the processor
            switch (get_code_char)
              {
                case eParameter :
                  var=get_code_int;
                  params[var] = Stack.top();
                  Stack.pop();
                  break;
                case eEndogenous :
                  var=get_code_int;
                  y[var] = Stack.top();
#ifdef DEBUG
                  tmp_out << "=>";
                  //mexPrintf(" y%d](%f)=%s\n", var, y[var], tmp_out.str().c_str());
                  tmp_out.str("");
#endif
                  Stack.pop();
                  break;
                case eExogenous :
                  var=get_code_int;
                  x[var]  = Stack.top();
                  Stack.pop();
                  break;
                case eExogenousDet :
                  var=get_code_int;
                  x[var] = Stack.top();
                  Stack.pop();
                  break;
                default:
                  mexPrintf("Unknown vraibale type\n");
              }
            break;
          case FSTPT :
            //store in a temporary variable from the processor
            var=get_code_int;
            T[var*(periods+y_kmin+y_kmax)+it_] = Stack.top();
#ifdef DEBUG
						tmp_out << "=>";
            //mexPrintf(" T[%d, %d](%f)=%s\n", it_, var, T[var*(periods+y_kmin+y_kmax)+it_], tmp_out.str().c_str());
            tmp_out.str("");
#endif
            Stack.pop();
            break;
          case FSTPST :
            //store in a temporary variable from the processor
            var=get_code_int;
            T[var] = Stack.top();
#ifdef DEBUG
						tmp_out << "=>";
            //mexPrintf(" T%d](%f)=%s\n", var, T[var], tmp_out.str().c_str());
            tmp_out.str("");
#endif
            Stack.pop();
            break;
          case FSTPU :
            //store in u variable from the processor
            var=get_code_int;
            var+=Per_u_;
            u[var] = Stack.top();
#ifdef DEBUG
						tmp_out << "=>";
						if(var==308)
              mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
            tmp_out.str("");
#endif
            Stack.pop();
            break;
          case FSTPSU :
            //store in u variable from the processor
            var=get_code_int;
            u[var] = Stack.top();
#ifdef DEBUG
						tmp_out << "=>";
						if(var==308)
              mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
            tmp_out.str("");
#endif
            Stack.pop();
            break;
          case FSTPR :
            //store in residual variable from the processor
            var=get_code_int;
            r[var] = Stack.top();
#ifdef DEBUG
            tmp_out << "=>";
            //mexPrintf(" r[%d](%f)=%s\n", var, r[var], tmp_out.str().c_str());
            tmp_out.str("");
#endif
            Stack.pop();
            break;
          case FSTPG :
            //store in derivative (g) variable from the processor
            var=get_code_int;
            g1[var] = Stack.top();
#ifdef DEBUG
            tmp_out << "=>";
            //mexPrintf(" r[%d](%f)=%s\n", var, r[var], tmp_out.str().c_str());
            tmp_out.str("");
#endif
            Stack.pop();
            break;
          case FBINARY :
            op=get_code_int;
            v2=Stack.top();
            Stack.pop();
            v1=Stack.top();
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
                  Stack.push(v1 / v2);
#ifdef DEBUG
                  tmp_out << " |" << v1 << "/" << v2 << "|";
#endif
                  break;
                case oLess:
                  Stack.push(double(v1<v2));
#ifdef DEBUG
                  tmp_out << " |" << v1 << "<" << v2 << "|";
#endif
                  break;
                case oGreater:
                  Stack.push(double(v1>v2));
#ifdef DEBUG
                  tmp_out << " |" << v1 << ">" << v2 << "|";
#endif
                  break;
                case oLessEqual:
                  Stack.push(double(v1<=v2));
#ifdef DEBUG
                  tmp_out << " |" << v1 << "<=" << v2 << "|";
#endif
                  break;
                case oGreaterEqual:
                  Stack.push(double(v1>=v2));
#ifdef DEBUG
                  tmp_out << " |" << v1 << ">=" << v2 << "|";
#endif
                  break;
                case oEqualEqual:
                  Stack.push(double(v1==v2));
#ifdef DEBUG
                  tmp_out << " |" << v1 << "==" << v2 << "|";
#endif
                  break;
                case oDifferent:
                  Stack.push(double(v1!=v2));
#ifdef DEBUG
                  tmp_out << " |" << v1 << "!=" << v2 << "|";
#endif
                  break;
                case oPower:
                  Stack.push(pow1(v1, v2));
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
          case FUNARY :
            op=get_code_int;
            v1=Stack.top();
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
                  Stack.push(log1(v1));
#ifdef DEBUG
                  tmp_out << " |log(" << v1 << ")|";
#endif
                  break;
                case oLog10:
                  Stack.push(log10(v1));
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
                  /*throw EvalException();*/
                  ;
              }
            break;
          case FCUML :
            v1=Stack.top();
            Stack.pop();
            v2=Stack.top();
            Stack.pop();
            Stack.push(v1+v2);
            break;
          case FENDBLOCK :
            //it's the block end
            go_on=false;
            break;
          case FENDEQU :
            break;
          case FOK :
            op=get_code_int;
            if (Stack.size()>0)
              {
                mexPrintf("error: Stack not empty!\n");
                mexEvalString("st=fclose('all');clear all;");
                mexErrMsgTxt("End of simulate");
              }
            break;
          default :
            mexPrintf("Unknow opcode %d!! FENDEQU=%d\n",cc,FENDEQU);
            mexEvalString("st=fclose('all');clear all;");
            mexErrMsgTxt("End of simulate");
            break;
        }
    }
}

bool
Interpreter::simulate_a_block(int size,int type, string file_name, string bin_basename, bool Gaussian_Elimination, bool steady_state, int block_num)
{
  char *begining;
  int i;
  bool is_linear, cvg;
  int max_lag_plus_max_lead_plus_1;
  int symbol_table_endo_nbr;
  int Block_List_Max_Lag;
  int Block_List_Max_Lead;
  int giter;
  int u_count_int;
  bool result = true;
  double *y_save;

  switch (type)
    {
      case EVALUATE_FORWARD :
        if(steady_state)
					compute_block_time(0);
				else
				  {
            begining=get_code_pointer;
            for (it_=y_kmin;it_<periods+y_kmin;it_++)
              {
                set_code_pointer(begining);
                Per_y_=it_*y_size;
                compute_block_time(0);
              }
				  }
        break;
      case EVALUATE_BACKWARD :
        if(steady_state)
          compute_block_time(0);
				else
				  {
            begining=get_code_pointer;
            for (it_=periods+y_kmin-1;it_>=y_kmin;it_--)
              {
                set_code_pointer(begining);
                Per_y_=it_*y_size;
                compute_block_time(0);
              }
				  }
        break;
      case SOLVE_FORWARD_SIMPLE :
        g1=(double*)mxMalloc(size*size*sizeof(double));
				r=(double*)mxMalloc(size*sizeof(double));
				begining=get_code_pointer;
        if(steady_state)
          {
          	cvg=false;
            iter=0;
            while (!(cvg||(iter>maxit_)))
              {
                set_code_pointer(begining);
                Per_y_=it_*y_size;
                    compute_block_time(0);
                    y[Block_Contain[0].Variable] += -r[0]/g1[0];
                    double rr;
							    	rr=r[0];
                    cvg=((rr*rr)<solve_tolf);
                    iter++;
                  }
                if (!cvg)
                  {
                    mexPrintf("Convergence not achieved in block %d, after %d iterations\n",Block_Count,iter);
                    /*mexEvalString("st=fclose('all');clear all;");
                    mexErrMsgTxt("End of simulate");*/
                    return false;
                  }
              }
				   else
				     {
                for (it_=y_kmin;it_<periods+y_kmin;it_++)
                  {
                    cvg=false;
                    iter=0;
                    Per_y_=it_*y_size;
                    while (!(cvg||(iter>maxit_)))
                      {
                        set_code_pointer(begining);
                        Per_y_=it_*y_size;
                        compute_block_time(0);
                        y[Per_y_+Block_Contain[0].Variable] += -r[0]/g1[0];
                        double rr;
                        if(fabs(1+y[Per_y_+Block_Contain[0].Variable])>eps)
				    					    rr=r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
						    		    else
								    	    rr=r[0];
                        cvg=((rr*rr)<solve_tolf);
                        iter++;
                      }
                    if (!cvg)
                      {
                        mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n",Block_Count,it_,iter);
                        mexEvalString("st=fclose('all');clear all;");
                        mexErrMsgTxt("End of simulate");
                     }
								}
					 }
        mxFree(g1);
        mxFree(r);
        break;
      case SOLVE_BACKWARD_SIMPLE :
        g1=(double*)mxMalloc(size*size*sizeof(double));
				r=(double*)mxMalloc(size*sizeof(double));
				begining=get_code_pointer;
        if(steady_state)
          {
          	cvg=false;
            iter=0;
            while (!(cvg||(iter>maxit_)))
              {
                set_code_pointer(begining);
                Per_y_=it_*y_size;
                    compute_block_time(0);
                    y[Block_Contain[0].Variable] += -r[0]/g1[0];
                    double rr;
    					    	rr=r[0];
                    cvg=((rr*rr)<solve_tolf);
                    iter++;
                  }
                if (!cvg)
                  {
                    mexPrintf("Convergence not achieved in block %d, after %d iterations\n",Block_Count,iter);
                    return false;
                    /*mexEvalString("st=fclose('all');clear all;");
                    mexErrMsgTxt("End of simulate");*/
                  }
              }
				    else
				      {
                for (it_=periods+y_kmin;it_>y_kmin;it_--)
                  {
                    cvg=false;
                    iter=0;
                    Per_y_=it_*y_size;
                    while (!(cvg||(iter>maxit_)))
                      {
                        set_code_pointer(begining);
                        Per_y_=it_*y_size;
                        compute_block_time(0);
                        y[Per_y_+Block_Contain[0].Variable] += -r[0]/g1[0];
                        double rr;
                        if(fabs(1+y[Per_y_+Block_Contain[0].Variable])>eps)
				    					    rr=r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
						    		    else
								    	    rr=r[0];
                        cvg=((rr*rr)<solve_tolf);
                        iter++;
                      }
                    if (!cvg)
                      {
                        mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n",Block_Count,it_,iter);
                        mexEvalString("st=fclose('all');clear all;");
                        mexErrMsgTxt("End of simulate");
                      }
                  }
            }
        mxFree(g1);
        mxFree(r);
        break;
      case SOLVE_FORWARD_COMPLETE :
        is_linear=get_code_bool;
        max_lag_plus_max_lead_plus_1=get_code_int;
        symbol_table_endo_nbr=get_code_int;
        Block_List_Max_Lag=get_code_int;
        Block_List_Max_Lead=get_code_int;
        u_count_int=get_code_int;
        fixe_u(&u, u_count_int, u_count_int);
        Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state);
        g1=(double*)mxMalloc(size*size*sizeof(double));
        r=(double*)mxMalloc(size*sizeof(double));
        begining=get_code_pointer;
        Per_u_ = 0;
				if(steady_state)
				  {
			      if (!is_linear)
              {
                max_res_idx=0;
                cvg=false;
                iter=0;
                //Per_y_=it_*y_size;
                while (!(cvg||(iter>maxit_)))
                  {
                 	  /*for (int j = 0; j < y_size; j++)
              	    mexPrintf("   variable %d at time %d and %d = %f\n", j+1, it_, it_+1, y[j+it_*y_size]);*/
                    set_code_pointer(begining);
                    error_not_printed = true;
                    res2=0;
						    		res1=0;
                    max_res=0;
                    compute_block_time(0);
                    /*if (isnan(res1)||isinf(res1))
                    {
                    memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                    }*/
                    if (!(isnan(res1)||isinf(res1)))
										  {
                        for (i=0; i<size ;i++)
                          {
                            double rr;
                            rr=r[i];
                            if (max_res<fabs(rr))
                              {
                                max_res=fabs(rr);
                                max_res_idx=i;
                              }
                            res2+=rr*rr;
                            res1+=fabs(rr);
													}

                        cvg=(max_res<solve_tolf);
											}
		    						else
				    					cvg=false;
                    result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true);
                    iter++;
                  }
                if (!cvg)
                  {
                    mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                    /*mexEvalString("st=fclose('all');clear all;");
                    mexErrMsgTxt("End of simulate");*/
                    return false;
                  }
              }
            else
              {
                set_code_pointer(begining);
                Per_y_=it_*y_size;
                iter = 0;
                res1=res2=max_res=0;max_res_idx=0;
                error_not_printed = true;
                compute_block_time(0);
                cvg=false;
                result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true);
              }
				  }
				else
				  {
			      if (!is_linear)
              {
                max_res_idx=0;
                for (it_=y_kmin;it_<periods+y_kmin;it_++)
                 {
                    cvg=false;
                    iter=0;
                    Per_y_=it_*y_size;
                    while (!(cvg||(iter>maxit_)))
                      {
                      	/*for (int j = 0; j < y_size; j++)
                  	    	mexPrintf("   variable %d at time %d and %d = %f\n", j+1, it_, it_+1, y[j+it_*y_size]);*/
                        set_code_pointer(begining);
                        error_not_printed = true;
                        res2=0;
						    				res1=0;
                        max_res=0;
                        compute_block_time(0);
                        /*if (isnan(res1)||isinf(res1))
                          {
                            memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                          }*/
                        if (!(isnan(res1)||isinf(res1)))
										      {
                            for (i=0; i<size ;i++)
                              {
                                double rr;
                                if(fabs(1+y[Per_y_+Block_Contain[i].Variable])>eps)
                                  rr=r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                                else
                                  rr=r[i];
                                if (max_res<fabs(rr))
                                  {
                                    max_res=fabs(rr);
                                    max_res_idx=i;
                                  }
                                res2+=rr*rr;
                                res1+=fabs(rr);
                              }
                            cvg=(max_res<solve_tolf);
    										  }
		    								else
				    						  cvg=false;
                        result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false);
                        iter++;
                      }
                    if (!cvg)
                      {
                        mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                        mexEvalString("st=fclose('all');clear all;");
                        mexErrMsgTxt("End of simulate");
                      }
                  }
              }
            else
              {
                for (it_=y_kmin;it_<periods+y_kmin;it_++)
                  {
                    set_code_pointer(begining);
                    Per_y_=it_*y_size;
                    iter = 0;
                    res1=res2=max_res=0;max_res_idx=0;
                    error_not_printed = true;
                    compute_block_time(0);
                    cvg=false;
                    result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false);
                  }
              }
				  }
        mxFree(index_equa);
        mxFree(index_vara);
        memset(direction,0,size_of_direction);
        mxFree(g1);
        mxFree(r);
        mxFree(u);
        break;
      case SOLVE_BACKWARD_COMPLETE :
        is_linear=get_code_bool;
        max_lag_plus_max_lead_plus_1=get_code_int;
        symbol_table_endo_nbr=get_code_int;
        Block_List_Max_Lag=get_code_int;
        Block_List_Max_Lead=get_code_int;
        u_count_int=get_code_int;
        fixe_u(&u, u_count_int, u_count_int);
        Read_SparseMatrix(bin_basename, size, 1, 0, 0, steady_state);
        g1=(double*)mxMalloc(size*size*sizeof(double));
        r=(double*)mxMalloc(size*sizeof(double));
        begining=get_code_pointer;
        if(steady_state)
				  {
			      if (!is_linear)
              {
                max_res_idx=0;
                cvg=false;
                iter=0;
                //Per_y_=it_*y_size;
                while (!(cvg||(iter>maxit_)))
                  {
                 	  /*for (int j = 0; j < y_size; j++)
              	    mexPrintf("   variable %d at time %d and %d = %f\n", j+1, it_, it_+1, y[j+it_*y_size]);*/
                    set_code_pointer(begining);
                    error_not_printed = true;
                    res2=0;
						    		res1=0;
                    max_res=0;
                    compute_block_time(0);
                    /*if (isnan(res1)||isinf(res1))
                    {
                    memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                    }*/
                    if (!(isnan(res1)||isinf(res1)))
										  {
                        for (i=0; i<size ;i++)
                          {
                            double rr;
                            rr=r[i];
                            if (max_res<fabs(rr))
                              {
                                max_res=fabs(rr);
                                max_res_idx=i;
                              }
                            res2+=rr*rr;
                            res1+=fabs(rr);
													}
                        cvg=(max_res<solve_tolf);
											}
		    						else
				    					cvg=false;
                    result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true);
                    iter++;
                  }
                if (!cvg)
                  {
                    mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                    /*mexEvalString("st=fclose('all');clear all;");
                    mexErrMsgTxt("End of simulate");*/
                    return false;
                  }
              }
            else
              {
                set_code_pointer(begining);
                Per_y_=it_*y_size;
                iter = 0;
                res1=res2=max_res=0;max_res_idx=0;
                error_not_printed = true;
                compute_block_time(0);
                cvg=false;
                result = simulate_NG(Block_Count, symbol_table_endo_nbr, 0, 0, 0, size, false, cvg, iter, true);
              }
				  }
				else
				 {
            if (!is_linear)
              {
                max_res_idx=0;
                for (it_=periods+y_kmin;it_>y_kmin;it_--)
                  {
                    cvg=false;
                    iter=0;
                    Per_y_=it_*y_size;
                    while (!(cvg||(iter>maxit_)))
                      {
                        set_code_pointer(begining);
                        error_not_printed = true;
                        res2=0;
                        res1=0;
                        max_res=0;
                        compute_block_time(0);
                        /*if (isnan(res1)||isinf(res1))
                          {
                            memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                          }*/
										    if (!(isnan(res1)||isinf(res1)))
										      {
                            for (i=0; i<size ;i++)
                              {
																double rr;
                                if(fabs(1+y[Per_y_+Block_Contain[i].Variable])>eps)
                                  rr=r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                                else
                                  rr=r[i];
                                if (max_res<fabs(rr))
                                  {
                                    max_res=fabs(rr);
                                    max_res_idx=i;
                                  }
                                res2+=rr*rr;
                                res1+=fabs(rr);
                              }
												    cvg=(max_res<solve_tolf);
										      }
										    else
										      cvg=false;
										    result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false);
                        iter++;
                      }
                    if (!cvg)
                      {
                        mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n", Block_Count, it_, iter);
                        mexEvalString("st=fclose('all');clear all;");
                        mexErrMsgTxt("End of simulate");
                      }
                  }
              }
            else
              {
                for (it_=periods+y_kmin;it_>y_kmin;it_--)
                  {
                    set_code_pointer(begining);
                    Per_y_=it_*y_size;
                    error_not_printed = true;
                    compute_block_time(0);
                    cvg=false;
                    result = simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter, false);
                  }
              }
				  }
				mxFree(index_equa);
        mxFree(index_vara);
        memset(direction,0,size_of_direction);
        mxFree(g1);
        mxFree(r);
        mxFree(u);
        break;
      case SOLVE_TWO_BOUNDARIES_SIMPLE :
      case SOLVE_TWO_BOUNDARIES_COMPLETE:
        if(steady_state)
          {
            mexPrintf("SOLVE_TXO_BOUNDARIES in a steady state model: impossible case\n");
            return false;
          }
        is_linear=get_code_bool;
        max_lag_plus_max_lead_plus_1=get_code_int;
        symbol_table_endo_nbr=get_code_int;
        Block_List_Max_Lag=get_code_int;
        Block_List_Max_Lead=get_code_int;
        u_count_int=get_code_int;
        fixe_u(&u, u_count_int, u_count_int);
        Read_SparseMatrix(bin_basename, size, periods, y_kmin, y_kmax, steady_state);
        u_count=u_count_int*(periods+y_kmax+y_kmin);
        r=(double*)mxMalloc(size*sizeof(double));
        y_save=(double*)mxMalloc(y_size*sizeof(double)*(periods+y_kmax+y_kmin));
        begining=get_code_pointer;
        if(!Gaussian_Elimination)
          {
          }
        giter=0;
        iter=0;
        if (!is_linear)
          {
            cvg=false;
            int u_count_saved=u_count;
            while (!(cvg||(iter>maxit_)))
              {
                res2=0;
                res1=0;
                max_res=0;
                max_res_idx=0;
                memcpy(y_save, y, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                for (it_=y_kmin;it_<periods+y_kmin;it_++)
                  {
                    Per_u_=(it_-y_kmin)*u_count_int;
                    Per_y_=it_*y_size;
                    set_code_pointer(begining);
                    compute_block_time(Per_u_);
                    if (isnan(res1)||isinf(res1))
                      {
                        memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                        break;
                      }
                    for (i=0; i< size; i++)
                      {
                        double rr;
                        if(fabs(1+y[Per_y_+Block_Contain[i].Variable])>eps)
                          rr=r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                        else
                          rr=r[i];
                        if (max_res<fabs(rr))
                          {
                            max_res=fabs(rr);
                            max_res_idx=i;
                          }
                        res2+=rr*rr;
                        res1+=fabs(rr);
                      }
                  }
								if (isnan(res1)||isinf(res1))
								  cvg = false;
								else
                  cvg=(max_res<solve_tolf);
                u_count=u_count_saved;
                simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter);
                iter++;
              }
            if (!cvg)
              {
                mexPrintf("Convergence not achieved in block %d, after %d iterations\n",Block_Count, iter);
                mexEvalString("st=fclose('all');clear all;");
                mexErrMsgTxt("End of simulate");
              }
          }
        else
          {
          	res1=res2=max_res=0;max_res_idx=0;
            for (it_=y_kmin;it_<periods+y_kmin;it_++)
              {
                Per_u_=(it_-y_kmin)*u_count_int;
                Per_y_=it_*y_size;
                set_code_pointer(begining);
                compute_block_time(Per_u_);
                for (i=0; i< size; i++)
                  {
                    double rr;
                    rr=r[i];
                    if (max_res<fabs(rr))
                      {
                        max_res=fabs(rr);
                        max_res_idx=i;
                      }
                    res2+=rr*rr;
                    res1+=fabs(rr);
                  }
              }
            cvg = false;
            simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter);
          }
        mxFree(r);
        mxFree(y_save);
        mxFree(u);
        mxFree(index_vara);
        mxFree(index_equa);
        memset(direction,0,size_of_direction);
        break;
      default:
        mexPrintf("Unknow type =%d\n",type);
        mexEvalString("st=fclose('all');clear all;");
        mexEvalString("drawnow;");
        mexErrMsgTxt("End of simulate");
    }
	return true;
}

bool
Interpreter::compute_blocks(string file_name, string bin_basename, bool steady_state)
{
  ifstream CompiledCode;
  bool result = true;
  int Code_Size, var;
  if(steady_state)
    file_name += "_static";
	else
	  file_name += "_dynamic";
  //First read and store in memory the code
  CompiledCode.open((file_name + ".cod").c_str(),std::ios::in | std::ios::binary| std::ios::ate);
  if (!CompiledCode.is_open())
    {
      mexPrintf("%s.cod Cannot be opened\n",file_name.c_str());
      mexEvalString("drawnow;");
      mexEvalString("st=fclose('all');clear all;");
      filename+=" stopped";
      mexEvalString("drawnow;");
      mexErrMsgTxt(filename.c_str());
    }
  Code_Size=CompiledCode.tellg();

  CompiledCode.seekg(std::ios::beg);
  Code=(char*)mxMalloc(Code_Size);
  CompiledCode.seekg(0);
  CompiledCode.read(reinterpret_cast<char *>(Code), Code_Size);
  CompiledCode.close();
  char *Init_Code=Code;

  //The big loop on intructions
  Block_Count=-1;
  bool go_on=true;
  while (go_on)
    {
      char code=get_code_char;
      switch (code)
        {
          case FBEGINBLOCK :
            //it's a new block

            Block_Count++;
            Block_type lBlock;
            Block.clear();
            Block_Contain.clear();
            Block_contain_type lBlock_Contain;
            lBlock.begin=get_code_pos-(uint64_t)Init_Code;
            lBlock.size=get_code_int;
            lBlock.type=get_code_int;
            Block.push_back(lBlock);
            for (int i=0;i<lBlock.size;i++)
              {
                lBlock_Contain.Variable=get_code_int;
                lBlock_Contain.Equation=get_code_int;
                lBlock_Contain.Own_Derivative=get_code_int;
                Block_Contain.push_back(lBlock_Contain);
              }
            result = simulate_a_block(lBlock.size,lBlock.type, file_name, bin_basename,true, steady_state, Block_Count);
            if(!result)
              go_on = false;
            break;
          case FEND :
            go_on=false;
            break;
          case FDIMT :
            var=get_code_int;
            if(T)
              mxFree(T);
            T=(double*)mxMalloc(var*(periods+y_kmin+y_kmax)*sizeof(double));
            break;
					case FDIMST :
            var=get_code_int;
            if(T)
              mxFree(T);
            T=(double*)mxMalloc(var*sizeof(double));
            break;
          default :
            mexPrintf("Unknow command : %d at pos %d !!\n",(long int)(code),(uint64_t*)(get_code_pos)-(uint64_t*)(Init_Code));
            mexEvalString("st=fclose('all');clear all;");
            mexEvalString("drawnow;");
            mexErrMsgTxt("End of simulate");
            break;
        }
    }
  mxFree(Init_Code);
  if(T)
    mxFree(T);
	return result;
}
