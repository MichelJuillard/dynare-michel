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
//#define DEBUGC
#include <cstring>
#include "Interpreter.hh"

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
  //GaussSeidel=true;
}

double
Interpreter::pow1(double a, double b)
{
  double r=pow_(a,b);
  if (isnan(r) || isinf(r))
    {
      //mexPrintf("pow(%f, %f)=%f\n",a, b, r);
      max_res=res1=res2=r;
      return(r);
    }
  else
    return r;
}


void
Interpreter::compute_block_time(int Per_u_) /*throw(EvalException)*/
{
  int var, lag, op;
  double v1, v2;
  char/*uint8_t*/ cc;
  bool go_on=true;
  double *ll;
  while (go_on)
    {
      //mexPrintf("it_=%d",it_);
      switch (cc=get_code_char)
        {
          case FLDV :
            //load a variable in the processor
#ifdef DEBUGC
            if(Block_Count==2)
              {
                mexPrintf("FLDV\n");
                mexEvalString("drawnow;");
              }
#endif
            switch (get_code_char)
              {
                case eParameter :
                  var=get_code_int;
#ifdef DEBUGC
                  if(Block_Count==2)
                    {
                      mexPrintf(" params[%d]=%f\n",var,params[var]);
                      mexEvalString("drawnow;");
                    }
#endif
                  Stack.push(params[var]);
                  break;
                case eEndogenous :
                  var=get_code_int;
                  lag=get_code_int;
#ifdef DEBUGC
                  if(Block_Count==2)
                    {
                      mexPrintf("y[%d, %d]=%f\n",it_+lag, var+1, y[(it_+lag)*y_size+var]);
                    }
#endif
                  Stack.push(y[(it_+lag)*y_size+var]);
                  break;
                case eExogenous :
#ifdef DEBUGC
                  mexPrintf("Exogenous\n");
#endif
                  var=get_code_int;
                  lag=get_code_int;
#ifdef DEBUGC
                  if(Block_Count==2)
                    {
                      mexPrintf("x[%d, %d]\n",it_+lag, var+1);
                    }
#endif
                  Stack.push(x[it_+lag+var*nb_row_x]);
                  break;
                case eExogenousDet :
#ifdef DEBUGC
                  mexPrintf("ExogenousDet\n");
#endif

                  var=get_code_int;
                  lag=get_code_int;
#ifdef DEBUGC
                  mexPrintf(" x(det)[%d]=%f\n",it_+lag+var*nb_row_xd,x[it_+lag+var*nb_row_xd]);
                  mexEvalString("drawnow;");
#endif
                  Stack.push(x[it_+lag+var*nb_row_xd]);
                  break;
                default:
                  mexPrintf("Unknown variable type\n");
              }
            break;
          case FLDT :
            //load a temporary variable in the processor
            var=get_code_int;
#ifdef DEBUGC
            mexPrintf("FLDT %d [%d]=%f Stack.size()=%d\n",var,var*(periods+y_kmin+y_kmax)+it_,T[var*(periods+y_kmin+y_kmax)+it_], Stack.size());
            mexEvalString("drawnow;");
#endif
            Stack.push(T[var*(periods+y_kmin+y_kmax)+it_]);
            break;
          case FLDU :
            //load u variable in the processor
#ifdef DEBUGC
            mexPrintf("FLDU\n");
            mexEvalString("drawnow;");
#endif
            var=get_code_int;
            var+=Per_u_;
            Stack.push(u[var]);
            break;
          case FLDR :
            //load u variable in the processor
            var=get_code_int;
#ifdef DEBUGC
            mexPrintf("FLDR r[%d]=%f\n",var, r[var] );
            mexEvalString("drawnow;");
#endif
            Stack.push(r[var]);
            break;
          case FLDZ :
            //load 0 in the processor
#ifdef DEBUGC
            mexPrintf("FLDZ\n");
            mexEvalString("drawnow;");
#endif
            Stack.push(0);
            break;
          case FLDC :
            //load a numerical constant in the processor
            /*asm("fld\n\t"
                "fstp %%st" : "=t" (ll) : "0" ((double)(*Code)));*/
            ll=get_code_pdouble;
#ifdef DEBUGC
            mexPrintf("FLDC %f\n",*ll);
            mexEvalString("drawnow;");
#endif
            Stack.push(*ll);
            break;
          case FSTPV :
            //load a variable in the processor
#ifdef DEBUGC
            mexPrintf("FSTPV\n");
            mexEvalString("drawnow;");
#endif
            switch (get_code_char)
              {
                case eParameter :
                  var=get_code_int;
                  params[var] = Stack.top();
#ifdef DEBUGC
                  mexPrintf("FSTP params[%d]=%f\n", var, params[var]);
                  mexEvalString("drawnow;");
#endif
                  Stack.pop();
                  break;
                case eEndogenous :
#ifdef DEBUGC
                  mexPrintf("FSTP Endogenous\n");
#endif
                  var=get_code_int;
                  lag=get_code_int;
#ifdef DEBUGC
                  //mexPrintf("y[%d(it_=%d, lag=%d, y_size=%d, var=%d)](%d)=",(it_+lag)*y_size+var,it_, lag, y_size, var, Stack.size());
                  mexPrintf("FSTP y[it_=%d, lag=%d, y_size=%d, var=%d, block=%d)]=",it_, lag, y_size, var+1, Block_Count+1);
                  mexEvalString("drawnow;");
#endif
                  y[(it_+lag)*y_size+var] = Stack.top();
#ifdef DEBUGC
                  mexPrintf("%f\n",y[(it_+lag)*y_size+var]);
                  mexEvalString("drawnow;");
#endif
                  Stack.pop();
                  break;
                case eExogenous :
#ifdef DEBUGC
                  mexPrintf("Exogenous\n");
#endif
                  //var=get_code_int;
                  var=get_code_int;
                  lag=get_code_int;
                  /*mexPrintf("FSTP x[it_=%d, lag=%d, y_size=%d, var=%d, block=%d)]=",it_, lag, y_size, var+1, Block_Count+1);
                  mexEvalString("drawnow;");*/
                  x[it_+lag+var*nb_row_x]  = Stack.top();
                  Stack.pop();
                  /*mexPrintf("%f\n",x[it_+lag+var*nb_row_x]);
                  mexEvalString("drawnow;");*/
                  break;
                case eExogenousDet :
#ifdef DEBUGC
                  mexPrintf("ExogenousDet\n");
#endif
                  var=get_code_int;
                  var=get_code_int;
                  lag=get_code_int;
                  x[it_+lag+var*nb_row_xd] = Stack.top();
                  Stack.pop();
                  break;
                default:
                  mexPrintf("Unknown vraibale type\n");
              }
            break;
          case FSTPT :
            //load a temporary variable in the processor
            var=get_code_int;
#ifdef DEBUGC
            mexPrintf("FSTPT T[(var=%d, it_=%d, periods=%d, y_kmin=%d, y_kmax=%d)%d]=", var, it_, periods, y_kmin, y_kmax, var*(periods+y_kmin+y_kmax)+it_);
            mexEvalString("drawnow;");
#endif
            T[var*(periods+y_kmin+y_kmax)+it_] = Stack.top();
#ifdef DEBUGC
            mexPrintf("%f\n",T[var*(periods+y_kmin+y_kmax)+it_]);
            mexEvalString("drawnow;");
#endif
            Stack.pop();
            break;
          case FSTPU :
            //load u variable in the processor
            var=get_code_int;
            var+=Per_u_;
#ifdef DEBUGC
            mexPrintf("FSTPU u[%d]",var);
            mexEvalString("drawnow;");
#endif
            u[var] = Stack.top();
#ifdef DEBUGC
            mexPrintf("=%f\n",u[var]);
            mexEvalString("drawnow;");
#endif
            Stack.pop();
            break;
          case FSTPR :
            //load u variable in the processor
            var=get_code_int;
#ifdef DEBUGC
            mexPrintf("FSTPR residual[%d]=",var);
            mexEvalString("drawnow;");
#endif
            r[var] = Stack.top();
#ifdef DEBUGC
            mexPrintf("%f\n",r[var]);
            mexEvalString("drawnow;");
#endif
            Stack.pop();
            break;
          case FSTPG :
            //load u variable in the processor
            var=get_code_int;
#ifdef DEBUGC
            mexPrintf("FSTPG g1[%d]=",var);
            mexEvalString("drawnow;");
#endif
            g1[var] = Stack.top();
#ifdef DEBUGC
            mexPrintf("%f\n",g1[var]);
            mexEvalString("drawnow;");
#endif
            Stack.pop();
            break;
          case FBINARY :
#ifdef DEBUGC
            mexPrintf("FBINARY\n");
            mexEvalString("drawnow;");
#endif
            op=get_code_int;
            v2=Stack.top();
            Stack.pop();
            v1=Stack.top();
            Stack.pop();
            switch (op)
              {
                case oPlus:
                  Stack.push(v1 + v2);
#ifdef DEBUGC
                  mexPrintf("+\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oMinus:
                  Stack.push(v1 - v2);
#ifdef DEBUGC
                  mexPrintf("-\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oTimes:
                  Stack.push(v1 * v2);
#ifdef DEBUGC
                  mexPrintf("*\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oDivide:
                  Stack.push(v1 / v2);
#ifdef DEBUGC
                  mexPrintf("/\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oLess:
                  Stack.push(double(v1<v2));
#ifdef DEBUGC
                  mexPrintf("%f < %f\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oGreater:
                  Stack.push(double(v1>v2));
#ifdef DEBUGC
                  mexPrintf("%f > %f\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oLessEqual:
                  Stack.push(double(v1<=v2));
#ifdef DEBUGC
                  mexPrintf("%f <= %f\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oGreaterEqual:
                  Stack.push(double(v1>=v2));
#ifdef DEBUGC
                  mexPrintf("%f >= %f\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oEqualEqual:
                  Stack.push(double(v1==v2));
#ifdef DEBUGC
                  mexPrintf("%f == %f\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oDifferent:
                  Stack.push(double(v1!=v2));
#ifdef DEBUGC
                  mexPrintf("%f > %f\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oPower:
#ifdef DEBUGC
                  mexPrintf("pow(%f, %f)\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  Stack.push(pow1(v1, v2));
                  break;
                case oMax:
                  Stack.push(max(v1, v2));
#ifdef DEBUGC
                  mexPrintf("max(%f, %f)\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oMin:
                  Stack.push(min(v1, v2));
#ifdef DEBUGC
                  mexPrintf("min(%f, %f)\n",v1,v2);
                  mexEvalString("drawnow;");
#endif
                  break;
                case oEqual:
                default:
                  /*throw EvalException();*/
                  ;
              }
            break;
          case FUNARY :
#ifdef DEBUGC
            mexPrintf("FUNARY\n");
            mexEvalString("drawnow;");
#endif
            op=get_code_int;
               v1=Stack.top();
            Stack.pop();
            switch (op)
              {
                case oUminus:
                  Stack.push(-v1);
#ifdef DEBUGC
                  mexPrintf("-\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oExp:
                  Stack.push(exp(v1));
#ifdef DEBUGC
                  mexPrintf("exp\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oLog:
                  Stack.push(log(v1));
#ifdef DEBUGC
                  mexPrintf("log\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oLog10:
                  Stack.push(log10(v1));
#ifdef DEBUGC
                  mexPrintf("log10\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oCos:
                  Stack.push(cos(v1));
#ifdef DEBUGC
                  mexPrintf("cos\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oSin:
                  Stack.push(sin(v1));
#ifdef DEBUGC
                  mexPrintf("sin\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oTan:
                  Stack.push(tan(v1));
#ifdef DEBUGC
                  mexPrintf("tan\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oAcos:
                  Stack.push(acos(v1));
#ifdef DEBUGC
                  mexPrintf("acos\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oAsin:
                  Stack.push(asin(v1));
#ifdef DEBUGC
                  mexPrintf("asin\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oAtan:
                  Stack.push(atan(v1));
#ifdef DEBUGC
                  mexPrintf("atan\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oCosh:
                  Stack.push(cosh(v1));
#ifdef DEBUGC
                  mexPrintf("cosh\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oSinh:
                  Stack.push(sinh(v1));
#ifdef DEBUGC
                  mexPrintf("sinh\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oTanh:
                  Stack.push(tanh(v1));
#ifdef DEBUGC
                  mexPrintf("tanh\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oAcosh:
                  Stack.push(acosh(v1));
#ifdef DEBUGC
                  mexPrintf("acosh\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oAsinh:
                  Stack.push(asinh(v1));
#ifdef DEBUGC
                  mexPrintf("asinh\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oAtanh:
                  Stack.push(atanh(v1));
#ifdef DEBUGC
                  mexPrintf("atanh\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                case oSqrt:
                  Stack.push(sqrt(v1));
#ifdef DEBUGC
                  mexPrintf("sqrt\n");
                  mexEvalString("drawnow;");
#endif
                  break;
                default:
                  /*throw EvalException();*/
                  ;
              }
            break;
          case FCUML :
#ifdef DEBUGC
            mexPrintf("FCUML\n");
            mexEvalString("drawnow;");
#endif
            v1=Stack.top();
            Stack.pop();
            v2=Stack.top();
            Stack.pop();
            Stack.push(v1+v2);
            break;
          case FENDBLOCK :
            //it's the block end
#ifdef DEBUGC
            //mexPrintf("FENDBLOCK\n");
            //mexEvalString("drawnow;");
#endif
            //Block[Block_Count].end=get_code_pos;
            go_on=false;
            break;
          case FENDEQU :
#ifdef DEBUGC
            mexPrintf("FENDEQU\n");
            mexEvalString("drawnow;");
#endif
            /*if (GaussSeidel)
              return;*/
            break;
          case FOK :
#ifdef DEBUGC
            mexPrintf("FOK\n");
            mexEvalString("drawnow;");
#endif
            op=get_code_int;
#ifdef DEBUGC
               mexPrintf("var=%d\n",op);
#endif
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

void
Interpreter::simulate_a_block(int size,int type, string file_name, string bin_basename, bool Gaussian_Elimination)
{
  /*mexPrintf("simulate_a_block\n");
  mexEvalString("drawnow;");*/

  char *begining;
  int i;
  bool is_linear, cvg;
  int max_lag_plus_max_lead_plus_1;
  int symbol_table_endo_nbr;
  int Block_List_Max_Lag;
  int Block_List_Max_Lead;
  int giter;
  int u_count_int;
  double *y_save;
#ifdef LINBCG
  LinBCG linbcg;
  Mat_DP a;
  Vec_INT indx;
#endif
    //SparseMatrix sparse_matrix;

  //int nb_endo, u_count_init;


  //mexPrintf("simulate_a_block\n");
  //mexEvalString("drawnow;");
  //mexPrintf("%d\n",debile);

  //GaussSeidel=false;
  //slowc_save=slowc/2;
  //mexPrintf("simulate_a_block size=%d type=%d\n",size,type);
  switch (type)
    {
      case EVALUATE_FORWARD :
      //case EVALUATE_FORWARD_R :
#ifdef DEBUGC
        mexPrintf("EVALUATE_FORWARD\n");
        mexEvalString("drawnow;");
#endif
        begining=get_code_pointer;
        //mexPrintf("y_kmin=%d periods=%d",y_kmin, periods);
        for (it_=y_kmin;it_<periods+y_kmin;it_++)
          {
            //mexPrintf("begining=%x\n",begining);
            set_code_pointer(begining);
            Per_y_=it_*y_size;
            //mexPrintf("bef compute_block_time()\n");
            compute_block_time(0);
            //mexPrintf("x(%d, 638)=%1.14f/x(%d, 116)=%1.14f=%1.14f\n",it_,x[it_+638*nb_row_x],it_,x[it_+116*nb_row_x],x[it_+638*nb_row_x]/x[it_+116*nb_row_x]);
#ifdef PRINT_OUT
            for (int j = 0; j<size; j++)
              mexPrintf("y[%d, %d] = %1.14f\n", Block_Contain[j].Variable, it_, y[Per_y_ + Block_Contain[j].Variable]);
#endif
          }
        break;
      case EVALUATE_BACKWARD :
      //case EVALUATE_BACKWARD_R :
#ifdef DEBUGC
        mexPrintf("EVALUATE_BACKWARD\n");
        mexEvalString("drawnow;");
#endif
        begining=get_code_pointer;
        for (it_=periods+y_kmin-1;it_>=y_kmin;it_--)
          {
            set_code_pointer(begining);
            Per_y_=it_*y_size;
            compute_block_time(0);
#ifdef PRINT_OUT
            for (int j = 0; j<size; j++)
              mexPrintf("y[%d, %d] = %1.14f\n", Block_Contain[j].Variable, it_, y[Per_y_ + Block_Contain[j].Variable]);
#endif
          }
        break;
      case SOLVE_FORWARD_SIMPLE :
#ifdef DEBUGC
        mexPrintf("SOLVE_FORWARD_SIMPLE\n");
        mexEvalString("drawnow;");
#endif
        g1=(double*)mxMalloc(size*size*sizeof(double));
        r=(double*)mxMalloc(size*sizeof(double));
        begining=get_code_pointer;
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
                //mexPrintf("y[%d] += -r[0] (%f) / g1[0] (%f) = %f\n",Per_y_+Block_Contain[0].Variable,r[0],g1[0],y[Per_y_+Block_Contain[0].Variable]);
                double rr;
                rr=r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
                cvg=((rr*rr)<solve_tolf);
                iter++;
              }
            if (!cvg)
              {
                mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n",Block_Count,it_,iter);
                mexEvalString("st=fclose('all');clear all;");
                mexErrMsgTxt("End of simulate");
              }
#ifdef PRINT_OUT
            mexPrintf("y[%d, %d]=%f \n",it_, Block_Contain[0].Variable ,y[Per_y_ + Block_Contain[0].Variable]);
#endif
          }
        mxFree(g1);
        mxFree(r);
        break;
      case SOLVE_BACKWARD_SIMPLE :
#ifdef DEBUGC
        mexPrintf("SOLVE_BACKWARD_SIMPLE\n");
        mexEvalString("drawnow;");
#endif
        g1=(double*)mxMalloc(size*size*sizeof(double));
        r=(double*)mxMalloc(size*sizeof(double));
        begining=get_code_pointer;
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
                /*mexPrintf("Compute_block_time=> in SOLVE_BACKWARD_SIMPLE : OK\n");
                mexEvalString("drawnow;");*/
                y[Per_y_+Block_Contain[0].Variable] += -r[0]/g1[0];
                double rr;
                rr=r[0]/(1+y[Per_y_+Block_Contain[0].Variable]);
                cvg=((rr*rr)<solve_tolf);
                iter++;
              }
            if (!cvg)
              {
                mexPrintf("Convergence not achieved in block %d, at time %d after %d iterations\n",Block_Count,it_,iter);
                mexEvalString("st=fclose('all');clear all;");
                mexErrMsgTxt("End of simulate");
              }
#ifdef PRINT_OUT
            mexPrintf("y[%d, %d]=%f \n",it_, Block_Contain[0].Variable ,y[Per_y_ + Block_Contain[0].Variable]);
#endif
          }
        mxFree(g1);
        mxFree(r);
        break;
      case SOLVE_FORWARD_COMPLETE :
#ifdef DEBUGC
        mexPrintf("SOLVE FORWARD_COMPLETE\n");
        mexEvalString("drawnow;");
#endif
        is_linear=get_code_bool;
        /*mexPrintf("is_linear=%d\n",is_linear);
        mexEvalString("drawnow;");*/
        max_lag_plus_max_lead_plus_1=get_code_int;
        /*mexPrintf("max_lag_plus_max_lead_plus_1=%d\n",max_lag_plus_max_lead_plus_1);
        mexEvalString("drawnow;");*/
        symbol_table_endo_nbr=get_code_int;
        /*mexPrintf("symbol_table_endo_nbr=%d\n",symbol_table_endo_nbr);
        mexEvalString("drawnow;");*/
        Block_List_Max_Lag=get_code_int;
        /*mexPrintf("Block_List_Max_Lag=%d\n",Block_List_Max_Lag);
        mexEvalString("drawnow;");*/
        Block_List_Max_Lead=get_code_int;
        /*mexPrintf("Block_List_Max_Lead=%d\n",Block_List_Max_Lead);
        mexEvalString("drawnow;");*/
        u_count_int=get_code_int;
        /*mexPrintf("u_count_int=%d\n",u_count_int);
        mexEvalString("drawnow;");*/
        fixe_u(&u, u_count_int, u_count_int);
        //Read_file(file_name, periods, 0, symbol_table_endo_nbr, Block_List_Max_Lag, Block_List_Max_Lead, nb_endo, u_count, u_count_init, u);
        //sparse_matrix.initialize(periods, nb_endo, y_kmin, y_kmax, y_size, u_count, u_count_init, u, y, ya, slowc, y_decal, markowitz_c, res1, res2, max_res);
        Read_SparseMatrix(bin_basename, size, 1, 0, 0);
        //mexPrintf("size=%d\n",size);
        g1=(double*)mxMalloc(size*size*sizeof(double));
        r=(double*)mxMalloc(size*sizeof(double));
        begining=get_code_pointer;
        Per_u_ = 0;

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
                    set_code_pointer(begining);
                    compute_block_time(0);
                    /*mexPrintf("Compute_block_time=> in SOLVE_FORWARD_COMPLETE : OK it_%d\n",it_);
                    mexEvalString("drawnow;");*/
                    //Direct_Simulate(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, 0, false, iter);
                    simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter);
                    res2=0;
                    res1=0;
                    max_res=0;
                    for (i=0; i<size ;i++)
                      {
                        double rr;
                        rr=r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                        if (max_res<fabs(rr))
                          {
                            max_res=fabs(rr);
                            max_res_idx=i;
                          }
                        res2+=rr*rr;
                        res1+=fabs(rr);
                      }
                    cvg=(max_res<solve_tolf);
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
                /*mexPrintf("Compute_block_time=> in SOLVE_FORWARD_COMPLETE before compute_block_time OK\n");
                mexEvalString("drawnow;");*/
                compute_block_time(0);
                //mexPrintf("Compute_block_time=> in SOLVE_FORWARD_COMPLETE : OK\n");
                /*mexPrintf("Compute_block_time=> in SOLVE_FORWARD_COMPLETE : %d OK\n",it_);
                mexEvalString("drawnow;");*/
                //Direct_Simulate(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, 0, false, iter);
                //Direct_Simulate(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, 0, false, iter);
                //simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, 0, /*true*/false, cvg, iter);
                cvg=false;
                /*mexPrintf("Compute_block_time=> in SOLVE_FORWARD_COMPLETE : OK it_%d\n",it_);
                mexEvalString("drawnow;");*/
                simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter);
                /*mexPrintf("End of simulate_NG\n");
                mexEvalString("drawnow;");*/
                //simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter);
              }
            /*mexPrintf("solve forward complete simulation\n");
            for (it_=y_kmin;it_<periods+y_kmin;it_++)
              {
                mexPrintf("it_=%d ",it_);
                for(i=0; i<y_size; i++)
                  mexPrintf(" %f",y[i+it_*y_size]);
                mexPrintf("\n");
              }*/
          }
        /*mexPrintf("end of forward solve\n");
        mexEvalString("drawnow;");*/
        mxFree(index_equa);
        mxFree(index_vara);
        memset(direction,0,size_of_direction);
        mxFree(g1);
        mxFree(r);
        mxFree(u);
        break;
      case SOLVE_BACKWARD_COMPLETE :
#ifdef DEBUGC
        mexPrintf("SOLVE_BACKWARD_COMPLETE\n");
        mexEvalString("drawnow;");
#endif
        is_linear=get_code_bool;
        //mexPrintf("is_linear=%d\n",is_linear);
        //mexEvalString("drawnow;");
        max_lag_plus_max_lead_plus_1=get_code_int;
        //mexPrintf("max_lag_plus_max_lead_plus_1=%d\n",max_lag_plus_max_lead_plus_1);
        //mexEvalString("drawnow;");
        symbol_table_endo_nbr=get_code_int;
        //mexPrintf("symbol_table_endo_nbr=%d\n",symbol_table_endo_nbr);
        //mexEvalString("drawnow;");
        Block_List_Max_Lag=get_code_int;
        //mexPrintf("Block_List_Max_Lag=%d\n",Block_List_Max_Lag);
        //mexEvalString("drawnow;");
        Block_List_Max_Lead=get_code_int;
        //mexPrintf("Block_List_Max_Lead=%d\n",Block_List_Max_Lead);
        //mexEvalString("drawnow;");
        u_count_int=get_code_int;
        /*mexPrintf("u_count_int=%d\n",u_count_int);
        mexEvalString("drawnow;");*/
        fixe_u(&u, u_count_int, u_count_int);
        Read_SparseMatrix(bin_basename, size, 1, 0, 0);
        //mexPrintf("size=%d\n",size);
        //mexEvalString("drawnow;");
        g1=(double*)mxMalloc(size*size*sizeof(double));
        r=(double*)mxMalloc(size*sizeof(double));
        begining=get_code_pointer;
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
                    compute_block_time(0);
                    /*mexPrintf("Compute_block_time=> in SOLVE_BACKWARD_COMPLETE : OK\n");
                    mexEvalString("drawnow;");*/
                    //Direct_Simulate(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, 0, false, iter);
                    simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter);
                    res2=0;
                    res1=0;
                    max_res=0;
                    for (i=0; i<size ;i++)
                      {
                        double rr;
                        rr=r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                        if (max_res<fabs(rr))
                          {
                            max_res=fabs(rr);
                            max_res_idx=i;
                          }
                        res2+=rr*rr;
                        res1+=fabs(rr);
                      }
                    cvg=(max_res<solve_tolf);
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
                compute_block_time(0);
                cvg=false;
                //Direct_Simulate(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, 1, false, iter);
                simulate_NG(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, false, cvg, iter);
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
#if GNUVER >= 432
        //mexPrintf("omp_get_max_threads=%d\n",omp_get_max_threads());
#endif
#ifdef DEBUGC
        mexPrintf("SOLVE_TWO_BOUNDARIES_COMPLETE\n");
        mexEvalString("drawnow;");
#endif
        is_linear=get_code_bool;
#ifdef DEBUGC
        mexPrintf("is_linear=%d\n",is_linear);
        mexEvalString("drawnow;");
#endif
        max_lag_plus_max_lead_plus_1=get_code_int;
#ifdef DEBUGC
        mexPrintf("max_lag_plus_max_lead_plus_1=%d\n",max_lag_plus_max_lead_plus_1);
        mexEvalString("drawnow;");
#endif
        symbol_table_endo_nbr=get_code_int;
#ifdef DEBUGC
        mexPrintf("symbol_table_endo_nbr=%d\n",symbol_table_endo_nbr);
        mexEvalString("drawnow;");
#endif
        Block_List_Max_Lag=get_code_int;
#ifdef DEBUGC
        mexPrintf("Block_List_Max_Lag=%d\n",Block_List_Max_Lag);
        mexEvalString("drawnow;");
#endif
        Block_List_Max_Lead=get_code_int;
#ifdef DEBUGC
        mexPrintf("Block_List_Max_Lead=%d\n",Block_List_Max_Lead);
        mexEvalString("drawnow;");
#endif
        u_count_int=get_code_int;
#ifdef DEBUGC
        mexPrintf("u_count_int=%d\n",u_count_int);
        mexPrintf("periods=%d\n",periods);
        mexEvalString("drawnow;");
#endif

        //sparse_matrix.initialize(periods, nb_endo, y_kmin, y_kmax, y_size, u_count, u_count_init, u, y, ya, slowc, y_decal, markowitz_c, res1, res2, max_res);

        //fixe_u(&u, u_count_int, max_lag_plus_max_lead_plus_1);
        fixe_u(&u, u_count_int, u_count_int);
#ifdef DEBUGC
        mexPrintf("u=%x\n",u);
        mexPrintf("size=%d\n",size);
#endif
        Read_SparseMatrix(bin_basename, size, periods, y_kmin, y_kmax);
#ifdef DEBUGC
        mexPrintf("size=%d\n",size);
        mexEvalString("drawnow;");
#endif
        //mexPrintf("aft reading_sparse_matrix\n");
        //mexEvalString("drawnow;");
        u_count=u_count_int*(periods+y_kmax+y_kmin);
        //g1=(double*)mxMalloc(size*size*sizeof(double));
        //mexPrintf("r=(double*)mxMalloc(%d)=",size*sizeof(double));
        r=(double*)mxMalloc(size*sizeof(double));
        //mexPrintf("%x\n",r);
        y_save=(double*)mxMalloc(y_size*sizeof(double)*(periods+y_kmax+y_kmin));
#ifdef DEBUGC
        mexPrintf("u_count=%d\n",u_count);
        mexEvalString("drawnow;");
#endif
        begining=get_code_pointer;
        if(!Gaussian_Elimination)
          {
#ifdef LINBCG
            it_=y_kmin;
            Per_u_=0;
            Per_y_=it_*y_size;
            set_code_pointer(begining);
            compute_block_time(0);
            linbcg.Initialize(filename, res1, res2, max_res, slowc, ya, direction, iter);
            linbcg.Preconditioner(periods, y_kmin, y_kmax, size, IM_i, index_vara, index_equa, y_size, y, true, 0, a, indx);
#endif
          }
        //GaussSeidel=false;
        giter=0;
        iter=0;
        //mexPrintf("2 boudaries problem\n");
        //mexEvalString("drawnow;");
        //mexPrintf("GaussSeidel=%d\n",GaussSeidel);
        if (!is_linear)
          {
            //double res1a=0;
            cvg=false;
            int u_count_saved=u_count;
            while (!(cvg||(iter>maxit_)))
              {
                res2=0;
                res1=0;
                max_res=0;
                max_res_idx=0;
                memcpy(y_save, y, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                //double res[size][13];
                for (it_=y_kmin;it_<periods+y_kmin;it_++)
                  {
                    Per_u_=(it_-y_kmin)*/*max_lag_plus_max_lead_plus_1*/u_count_int;
                    //mexPrintf("Per_u_=%d\n",Per_u_);
                    Per_y_=it_*y_size;
                    //mexPrintf("ok\n");
                    /*mexPrintf("compute_block_time it_=%d periods=%d y_kmin=%d\n",it_, periods, y_kmin);
                    mexEvalString("drawnow;");*/
                    set_code_pointer(begining);
                    compute_block_time(Per_u_);
                    /*mexPrintf("periods it_=%d ",it_);
                    for (i=0; i< size; i++)
                      mexPrintf(" %f ",r[i]);
								    mexPrintf("\n");*/
                    /*mexPrintf("end of compute_block_time it_=%d\n",it_);*/
                    /*if(Gaussian_Elimination)
                      initialize(periods, nb_endo, y_kmin, y_kmax, y_size, u_count, u_count_init, u, y, ya, slowc, y_decal, markowitz_c, res1, res2, max_res);*/
                    //mexPrintf("ok1\n");
                    /*mexEvalString("drawnow;");*/
                    if (isnan(res1)||isinf(res1))
                      {
                        memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
                        break;
                      }
                    for (i=0; i< size; i++)
                      {
                      	/*if(it_<13)
													res[i][it_-y_kmin]=r[i];*/
                        double rr;
                        /*if(fabs(y[Per_y_+Block_Contain[i].Variable])>solve_tolf)*/
                        //mexPrintf("res[%d]=%f\n",i,r[i]);
                        if(fabs(1+y[Per_y_+Block_Contain[i].Variable])>eps)
                          rr=r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                        else
                          rr=r[i];
                        /*else
                          rr=r[i];*/
                        if (max_res<fabs(rr))
                          {
                            max_res=fabs(rr);
                            max_res_idx=i;
                          }
                        res2+=rr*rr;
                        res1+=fabs(rr);
                        /*if (GaussSeidel && giter)
                          {
                            //mexPrintf("y[%d]-=r[%d]/u[%d]\n",Block_Contain[i].Variable,i,Block_Contain[i].Own_Derivative,);
                            y[Per_y_+Block_Contain[i].Variable]-=r[i]/u[Per_u_+Block_Contain[i].Own_Derivative];
                            //mexPrintf("y[%d]-=r[%d](%f)/u[%d](%f)=%f\n",Block_Contain[i].Variable,i,r[i],Block_Contain[i].Own_Derivative,u[Per_u_+Block_Contain[i].Own_Derivative], y[Per_y_+Block_Contain[i].Variable]);
                          }*/
                        /*mexPrintf("r[%d] (i=%d)",i+size*(it_-y_kmin),i);
                        mexPrintf("=%f\n",r[i]);*/
											  //mexPrintf("u[b[%d]=%d]=%f\n",i+(it_-y_kmin)*u_count_int,b[i+(it_-y_kmin)*u_count_int],u[b[i+(it_-y_kmin)*u_count_int]]);
                      }
										//mexPrintf("---------------------------------------------------------------------------------\n");
                    //mexPrintf("(log(y(%d, 151))) - (x(%d, 338)=%f\n",it_, it_,y[Per_y_+ 150] - exp(x[it_+337*nb_row_x]));
                  }
								/*for(i=0;i<size;i++)
								  {
    							   for(int j=0;j<5;j++)
									     mexPrintf(" % 1.6f ",res[i][j]);
										 mexPrintf("\n");
								  }*/
                cvg=(max_res<solve_tolf);
                //mexPrintf("it_=%d\n",it_);
                if(Gaussian_Elimination)
                  {
                    /*mexPrintf("bef simulate_NG1\n");
                    mexEvalString("drawnow;");*/
                    u_count=u_count_saved;
                    /*mexPrintf("u_count=%d &u_count=%x\n",u_count,&u_count);
                    mexEvalString("drawnow;");*/
                    /*for(int t=0;t<periods;t++)
											{
				  	            mexPrintf("%5.14f %5.14f %5.14f %5.14f \n",T[437*(periods+y_kmin+y_kmax)+t], T[72*(periods+y_kmin+y_kmax)+t], T[473*(periods+y_kmin+y_kmax)+t],
				  	            T[437*(periods+y_kmin+y_kmax)+t] * T[72*(periods+y_kmin+y_kmax)+t] * T[473*(periods+y_kmin+y_kmax)+t]);
				              }
										filename=" stopped";
                    mexErrMsgTxt(filename.c_str());*/
                    simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter);
                    /*mexPrintf("after simulate_NG1\n");
                    mexEvalString("drawnow;");*/
                    /*mexPrintf("aft simulate_NG1\n");
                    mexEvalString("drawnow;");*/
                  }
                else
                  {
#ifdef LINBCG
                    linbcg.Initialize(filename, res1, res2, max_res, slowc, ya, direction, iter);
                    linbcg.SolveLinear(periods, y_kmin, y_kmax, size, IM_i, index_vara, index_equa,y_size,y, true, cvg, a, indx);
#endif
                  }
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
                Per_u_=(it_-y_kmin)*/*max_lag_plus_max_lead_plus_1*/u_count_int;
                Per_y_=it_*y_size;
                set_code_pointer(begining);
                compute_block_time(Per_u_);
                /*mexPrintf("periods it_=%d ",it_);
                for (i=0; i< size; i++)
                  mexPrintf(" %f ",r[i]);
								mexPrintf("\n");*/
#ifdef PRINT_OUT
                for (j=0; j<max_lag_plus_max_lead_plus_1; j++)
                  {
                    mexPrintf(" %f",u[Per_u_+j]);
                  }
                mexPrintf("\n");
#endif
                /*mexPrintf("it_=%d ",it_);
                for(i=0; i<y_size; i++)
                  mexPrintf(" %f",y[i]);
                mexPrintf("\n");*/
                for (i=0; i< size; i++)
                {
                  double rr;
                  /*if(fabs(y[Per_y_+Block_Contain[i].Variable])>solve_tolf)*/
                  //mexPrintf("res[%d]=%f\n",i,r[i]);
                  /*if(fabs(1+y[Per_y_+Block_Contain[i].Variable])>eps)
                    rr=r[i]/(1+y[Per_y_+Block_Contain[i].Variable]);
                  else*/
                    rr=r[i];
                  /*else
                    rr=r[i];*/
                  if (max_res<fabs(rr))
                    {
                      max_res=fabs(rr);
                      max_res_idx=i;
                    }
                  res2+=rr*rr;
                  res1+=fabs(rr);
                }
              }
            //res1=res2=max_res=0;max_res_idx=0;

            cvg = false;
            //mexPrintf("it_=%d\n",it_);
            if(Gaussian_Elimination)
              {
                simulate_NG1(Block_Count, symbol_table_endo_nbr, it_, y_kmin, y_kmax, size, periods, true, cvg, iter);
                /*mexPrintf("after simulate_NG1\n");
                mexEvalString("drawnow;");*/
              }
            else
              {
#ifdef LINBCG
                linbcg.Initialize(filename, res1, res2, max_res, slowc, ya, direction, iter);
                linbcg.SolveLinear(periods, y_kmin, y_kmax, size, IM_i, index_vara, index_equa, y_size, y, true, cvg, a, indx);
#endif
              }
            /*mexPrintf("Two boundaries simulation\n");
            for (it_=y_kmin;it_<periods+y_kmin;it_++)
              {
                mexPrintf("it_=%d ",it_);
                for(i=0; i<y_size; i++)
                  mexPrintf(" %f",y[i+it_*y_size]);
                mexPrintf("\n");
              }*/
          }
#ifdef  DEBUGC
        //mexErrMsgTxt("End of simulate");
#endif

        //mxFree(g1);
        /*mexPrintf("end of simulate_a_block\n");
        mexEvalString("drawnow;");*/
        mxFree(r);
        mxFree(y_save);
        mxFree(u);
        mxFree(index_vara);
        mxFree(index_equa);
        memset(direction,0,size_of_direction);
        /*mexPrintf("after free\n");
        mexEvalString("drawnow;");*/
        //GaussSeidel=false;
        break;
      default:
        mexPrintf("Unknow type =%d\n",type);
        mexEvalString("st=fclose('all');clear all;");
        mexEvalString("drawnow;");
        mexErrMsgTxt("End of simulate");
    }
  /*mexPrintf("finish simulate_a_block\n");
  mexEvalString("drawnow;");*/
}

void
Interpreter::compute_blocks(string file_name, string bin_basename)
{
  ifstream CompiledCode;
  int Code_Size, var;
  //printf("open(%s)\n",(file_name + ".cod").c_str());
  //First read and store inn memory the code
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
#ifdef DEBUGC
  mexPrintf("Code_Size=%d\n",Code_Size);
  mexEvalString("drawnow;");
#endif
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
#ifdef DEBUGC
      mexPrintf("pos=%d\n",int(get_code_pos)-int(Init_Code));
      mexEvalString("drawnow;");
#endif
      char code=get_code_char;
#ifdef DEBUGC
      int icode=(int)code;
      mexPrintf("code=%d\n",icode);
      mexEvalString("drawnow;");
#endif
      switch (code)
        {
          case FBEGINBLOCK :
            //it's a new block

            Block_Count++;
            Block_type lBlock;
            Block.clear();
            Block_Contain.clear();
            Block_contain_type lBlock_Contain;
#ifdef DEBUGC
            mexPrintf("FBEGINBLOCK\n");
            mexEvalString("drawnow;");
#endif
            /*uint64_t *Init_Code_64;
            memcpy(Init_Code_64, Init_Code, sizeof(*Init_Code_64));*/
            lBlock.begin=get_code_pos-(uint64_t)Init_Code;
#ifdef DEBUGC
            mexPrintf("Block[%d].begin=%d\n",Block_Count, lBlock.begin);
            mexEvalString("drawnow;");
#endif
            lBlock.size=get_code_int;
#ifdef DEBUGC
            mexPrintf("Block[Block_Count].size=%d\n",lBlock.size);
            mexEvalString("drawnow;");
#endif
            lBlock.type=get_code_int;
#ifdef DEBUGC
            mexPrintf("Block[Block_Count].type=%d\n",lBlock.type);
            mexEvalString("drawnow;");
#endif
            Block.push_back(lBlock);
            for (int i=0;i</*Block[Block_Count].size*/lBlock.size;i++)
              {
                lBlock_Contain.Variable=get_code_int;
#ifdef DEBUGC
                mexPrintf("Block_Contain[%d].Variable=%d\n",i,lBlock_Contain.Variable);
                mexEvalString("drawnow;");
#endif
                lBlock_Contain.Equation=get_code_int;
#ifdef DEBUGC
                mexPrintf("Block_Contain[%d].Equation=%d\n",i,lBlock_Contain.Equation);
                mexEvalString("drawnow;");
#endif
                lBlock_Contain.Own_Derivative=get_code_int;
                //mexPrintf("Block_Contain[%d].Own_Derivative=%d\n",i,lBlock_Contain.Own_Derivative);
                Block_Contain.push_back(lBlock_Contain);
              }
#ifdef DEBUGC
            mexPrintf("Block Completed\n");
            mexEvalString("drawnow;");
#endif
            simulate_a_block(lBlock.size,lBlock.type, file_name, bin_basename,true);
            /*mexPrintf("after simulate_a_block\n");
            mexEvalString("drawnow;");*/
            //simulate_a_block(lBlock.size,lBlock.type, file_name, bin_basename,false);
            break;
          case FEND :
#ifdef DEBUGC
            mexPrintf("FEND\n");
            mexEvalString("drawnow;");
#endif
            go_on=false;
            break;
          case FDIMT :
            var=get_code_int;
#ifdef DEBUGC
                mexPrintf("FDIMT var=%d mxMalloc(%d)\n",var,var*(periods+y_kmin+y_kmax)*sizeof(double));
                mexEvalString("drawnow;");
#endif
            if(T)
              mxFree(T);
            T=(double*)mxMalloc(var*(periods+y_kmin+y_kmax)*sizeof(double));
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
  /*mexPrintf("compute_blocks\n");
  mexEvalString("drawnow;");*/
}
