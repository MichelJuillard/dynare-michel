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

////////////////////////////////////////////////////////////////////////
//                           simulate.cc                              //
//              simulate file designed for GNU GCC C++ compiler       //
//                 use NO_COMPILER option in MODEL command            //
////////////////////////////////////////////////////////////////////////

#include <cstring>

#include "simulate.hh"
#include "Interpreter.hh"
#include "Mem_Mngr.hh"
#include "linbcg.hh"



int
max(int a, int b)
{
  if (a>b)
    return(a);
  else
    return(b);
}







/*class EvalException
{
};*/


//#define DEBUG
/* The gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *M_, *oo_, *options_;
  int i, row_y, col_y, row_x, col_x;
  double * pind ;
  double *direction;
  //mexPrintf("mexFunction\n");
  //mexEvalString("drawnow;");
  /* Gets model parameters from global workspace of Matlab */
  //mexPrintf("starting simulation\n");
  M_ = mexGetVariable("global","M_");
  if (M_ == NULL )
    {
      mexPrintf("Global variable not found : ");
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("M_ \n");
    }
  /* Gets variables and parameters from global workspace of Matlab */
  oo_ = mexGetVariable("global","oo_");
  if (oo_ == NULL )
    {
      mexPrintf("Global variable not found : ");
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("oo_ \n");
    }
  options_ = mexGetVariable("global","options_");
  if (options_ == NULL )
    {
      mexPrintf("Global variable not found : ");
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("options_ \n");
    }
  //mexPrintf("ok0\n");
  params = mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"params")));
  double *yd, *xd;
  yd= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"endo_simul")));
  row_y=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"endo_simul")));
  xd= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_simul")));
  row_x=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_simul")));
  col_x=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_simul")));
  y_kmin=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"maximum_lag"))))));
  y_kmax=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"maximum_lead"))))));
  y_decal=max(0,y_kmin-int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"maximum_endo_lag")))))));
  periods=int(floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"periods"))))));
  maxit_=int(floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"maxit_"))))));
  slowc=double(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"slowc")))));
  slowc_save=slowc;
  markowitz_c=double(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"markowitz")))));
  nb_row_xd=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"exo_det_nbr"))))));
  mxArray *mxa=mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"fname"));
  int buflen=mxGetM(mxa) * mxGetN(mxa) + 1;
  char *fname;
  fname=(char*)mxCalloc(buflen, sizeof(char));
  string file_name=fname;
  int status = mxGetString(mxa, fname, buflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space. Filename is truncated.");
  //mexPrintf("fname=%s\n",fname);
  col_y=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"endo_simul")));;
  if (col_y<row_x)
    {
      row_y=row_y/row_x;
      col_y=row_x;
    }
  solve_tolf=*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"dynatol"))));
  //mexPrintf("col_y=%d row_y=%d\n",col_y, row_y);
  size_of_direction=col_y*row_y*sizeof(double);
  y=(double*)mxMalloc(size_of_direction);
  ya=(double*)mxMalloc(size_of_direction);
  direction=(double*)mxMalloc(size_of_direction);
  memset(direction,0,size_of_direction);
  x=(double*)mxMalloc(col_x*row_x*sizeof(double));
  for (i=0;i<row_x*col_x;i++)
    x[i]=double(xd[i]);
  for (i=0;i<row_y*col_y;i++)
    y[i]=double(yd[i]);
  y_size=row_y;
  x_size=row_x;
  nb_row_x=row_x;

  /* Call the C subroutines. */
  //mexPrintf("Call subroutines\n");
  //mexEvalString("drawnow;");

  //t0= pctimer();
  t0= clock();
  Interpreter interprete(params, y, ya, x, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal, markowitz_c, file_name);
  string f(fname);
  interprete.compute_blocks(f+"_dynamic", f);
  //t1= pctimer();
  t1= clock();
  mexPrintf("Simulation Time=%f milliseconds\n",1000.0*(double(t1)-double(t0))/double(CLOCKS_PER_SEC));
  //mexPrintf("SaveCode.is_open()=%d nlhs=%d \n",SaveCode.is_open(),nlhs);
  //interprete.sparse_matrix.close_SaveCode();

  //mexPrintf("End all-1\n");
  if (nlhs>0)
    {
      plhs[0] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);
      pind = mxGetPr(plhs[0]);
      for (i=0;i<row_y*col_y;i++)
        pind[i]=y[i];
    }
  //mexPrintf("End all0\n");
  if(x)
    mxFree(x);
  if(y)
    mxFree(y);
  if(ya)
    mxFree(ya);
  if(direction)
    mxFree(direction);
  //mexPrintf("End all\n");
}
