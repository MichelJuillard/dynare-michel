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



int
max(int a, int b)
{
  if (a>b)
    return(a);
  else
    return(b);
}



#ifdef DEBUG_EX
/*The Matlab c++ interface*/

using namespace std;
#include <sstream>
#include "mex_interface.hh"

int
main( int argc, const char* argv[] )
{
  FILE *fid;
  printf("argc=%d\n",argc);
  float f_tmp;
  ostringstream tmp_out("");
  tmp_out << argv[1] << "_options.txt";
  cout << tmp_out.str().c_str() << "\n";
  int nb_params;
  int i, row_y, col_y, row_x, col_x;
  double *yd, *xd;
  double *direction;

  string file_name(argv[1]);
  //mexPrintf("file_name=%s\n",file_name.c_str());

  fid = fopen(tmp_out.str().c_str(),"r");
  fscanf(fid,"%d",&periods);
  fscanf(fid,"%d",&maxit_);
  fscanf(fid,"%f",&f_tmp);
  slowc = f_tmp;
  //mexPrintf("slowc_save=%f\n",slowc_save);
  fscanf(fid,"%f",&f_tmp);
  markowitz_c = f_tmp;
  fscanf(fid,"%f",&f_tmp);
  solve_tolf = f_tmp;
  fclose(fid);

  tmp_out.str("");
  tmp_out << argv[1] << "_M.txt";
  //printf("%s\n",tmp_out.str().c_str());
  fid = fopen(tmp_out.str().c_str(),"r");
  fscanf(fid,"%d",&y_kmin);
  //printf("y_kmin=%d\n",y_kmin);
  fscanf(fid,"%d",&y_kmax);
  //printf("y_kmax=%d\n",y_kmax);
  fscanf(fid,"%d",&y_decal);
  //printf("y_decal=%d\n",y_decal);
  fscanf(fid,"%d",&nb_params);
  //printf("nb_params=%d\n",nb_params);
  fscanf(fid,"%d",&row_x);
  //printf("row_x=%d\n",row_x);
  fscanf(fid,"%d",&col_x);
  //printf("col_x=%d\n",col_x);
  fscanf(fid,"%d",&row_y);
  //printf("row_y=%d\n",row_y);
  fscanf(fid,"%d",&col_y);
  //printf("col_y=%d\n",col_y);
  fscanf(fid,"%d",&nb_row_xd);
  //printf("nb_row_xd=%d\n",nb_row_xd);
  params = (double*)malloc(nb_params*sizeof(params[0]));
  //printf("OK1\n");
  for(i=0; i < nb_params; i++)
    {
      fscanf(fid,"%f",&f_tmp);
      params[i] = f_tmp;
      //printf("param[%d]=%f\n",i,params[i]);
    }
  //printf("OK2\n");
  fclose(fid);
  yd = (double*)malloc(row_y*col_y*sizeof(yd[0]));
  xd = (double*)malloc(row_x*col_x*sizeof(xd[0]));
  tmp_out.str("");
  tmp_out << argv[1] << "_oo.txt";
  fid = fopen(tmp_out.str().c_str(),"r");
  for(i=0; i < col_y*row_y; i++)
    {
      fscanf(fid,"%f",&f_tmp);
      yd[i] = f_tmp;
    }
  for(i=0; i < col_x*row_x; i++)
    {
      fscanf(fid,"%f",&f_tmp);
      xd[i] = f_tmp;
    }
  fclose(fid);

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
  free(yd);
  free(xd);

  y_size=row_y;
  x_size=col_x/*row_x*/;
  nb_row_x=row_x;
  /*for(int i=0; i<y_size; i++)
    {
      for(int it_=0; it_<8;it_++)
        mexPrintf("y[t=%d, var=%d]=%f  ",it_+1, i+1, y[(it_)*y_size+i]);
      mexPrintf("\n");
    }

  for(int i=0; i<col_x; i++)
    {
      for(int it_=0; it_<8;it_++)
        mexPrintf("x[t=%d, var=%d]=%f  ",it_, i+1, x[it_+i*nb_row_x]);
      mexPrintf("\n");
    }*/

  t0= clock();
  Interpreter interprete(params, y, ya, x, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal, markowitz_c, file_name);
  string f(file_name);
  interprete.compute_blocks(f+"_dynamic", f);
  t1= clock();
  mexPrintf("Simulation Time=%f milliseconds\n",1000.0*(double(t1)-double(t0))/double(CLOCKS_PER_SEC));
  /*if (nlhs>0)
    {
      plhs[0] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);
      pind = mxGetPr(plhs[0]);
      for (i=0;i<row_y*col_y;i++)
        pind[i]=y[i];
    }*/
  if(x)
    mxFree(x);
  if(y)
    mxFree(y);
  if(ya)
    mxFree(ya);
  if(direction)
    mxFree(direction);
  free(params);
}

#else
/* The gateway routine */
void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *M_, *oo_, *options_;
  int i, row_y, col_y, row_x, col_x;
  double * pind ;
  double *direction;
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
  //slowc_save=slowc;
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
  col_y=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"endo_simul")));;
  solve_tolf=*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"dynatol"))));
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
  x_size=col_x/*row_x*/;
  nb_row_x=row_x;

  /*for(int i=0; i<y_size; i++)
    {
      for(int it_=0; it_<8;it_++)
        mexPrintf("y[t=%d, var=%d]=%f  ",it_+1, i+1, y[(it_)*y_size+i]);
      mexPrintf("\n");
    }

  for(int i=0; i<col_x; i++)
    {
      for(int it_=0; it_<8;it_++)
        mexPrintf("x[t=%d, var=%d]=%f  ",it_, i+1, x[it_+i*nb_row_x]);
      mexPrintf("\n");
    }*/

  t0= clock();
  Interpreter interprete(params, y, ya, x, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal, markowitz_c, file_name);
  string f(fname);
  interprete.compute_blocks(f+"_dynamic", f);
  t1= clock();
  mexPrintf("Simulation Time=%f milliseconds\n",1000.0*(double(t1)-double(t0))/double(CLOCKS_PER_SEC));
  if (nlhs>0)
    {
      plhs[0] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);
      pind = mxGetPr(plhs[0]);
      for (i=0;i<row_y*col_y;i++)
        pind[i]=y[i];
    }
  if(x)
    mxFree(x);
  if(y)
    mxFree(y);
  if(ya)
    mxFree(ya);
  if(direction)
    mxFree(direction);
}
#endif
