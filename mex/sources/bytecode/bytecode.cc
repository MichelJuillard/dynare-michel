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
#include "bytecode.hh"
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

using namespace std;
#include <sstream>
#include "mex_interface.hh"

int
main( int argc, const char* argv[] )
{
  FILE *fid;
  bool steady_state = false;
  printf("argc=%d\n",argc);
  if(argc<2)
    {
    	mexPrintf("model filename expected\n");
    	mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("Exit from Dynare");
    }
  float f_tmp;
  ostringstream tmp_out("");
  tmp_out << argv[1] << "_options.txt";
  cout << tmp_out.str().c_str() << "\n";
  int nb_params;
  int i, row_y, col_y, row_x, col_x;
  double *yd, *xd;
  double *direction;

  string file_name(argv[1]);

  if(argc>2)
    {
      string f(argv[1]);
      if(f == "steady_state")
        steady_state = true;
    }
  fid = fopen(tmp_out.str().c_str(),"r");
  int periods;
  fscanf(fid,"%d",&periods);
  int maxit_;
  fscanf(fid,"%d",&maxit_);
  fscanf(fid,"%f",&f_tmp);
  double slowc = f_tmp;
  fscanf(fid,"%f",&f_tmp);
  double markowitz_c = f_tmp;
  fscanf(fid,"%f",&f_tmp);
  double solve_tolf = f_tmp;
  fclose(fid);

  tmp_out.str("");
  tmp_out << argv[1] << "_M.txt";
  fid = fopen(tmp_out.str().c_str(),"r");
  int y_kmin;
  fscanf(fid,"%d",&y_kmin);
  int y_kmax;
  fscanf(fid,"%d",&y_kmax);
  int y_decal;
  fscanf(fid,"%d",&y_decal);
  fscanf(fid,"%d",&nb_params);
  fscanf(fid,"%d",&row_x);
  fscanf(fid,"%d",&col_x);
  fscanf(fid,"%d",&row_y);
  fscanf(fid,"%d",&col_y);
  int nb_row_xd;
  fscanf(fid,"%d",&nb_row_xd);
  double * params = (double*)malloc(nb_params*sizeof(params[0]));
  for(i=0; i < nb_params; i++)
    {
      fscanf(fid,"%f",&f_tmp);
      params[i] = f_tmp;
    }
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

  int size_of_direction=col_y*row_y*sizeof(double);
  double * y=(double*)mxMalloc(size_of_direction);
  double * ya=(double*)mxMalloc(size_of_direction);
  direction=(double*)mxMalloc(size_of_direction);
  memset(direction,0,size_of_direction);
  double * x=(double*)mxMalloc(col_x*row_x*sizeof(double));
  for (i=0;i<row_x*col_x;i++)
     x[i]=double(xd[i]);
  for (i=0;i<row_y*col_y;i++)
    y[i]=double(yd[i]);
  free(yd);
  free(xd);

  int y_size=row_y;
  int nb_row_x=row_x;
  clock_t t0= clock();
  Interpreter interprete(params, y, ya, x, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal, markowitz_c, file_name);
  string f(file_name);
  interprete.compute_blocks(f, f, steady_state);
  clock_t t1= clock();
  mexPrintf("Simulation Time=%f milliseconds\n",1000.0*(double(t1)-double(t0))/double(CLOCKS_PER_SEC));
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
  int i, row_y, col_y, row_x, col_x, nb_row_xd;
  int y_kmin=0, y_kmax=0, y_decal=0, periods=1;
  double * pind ;
  double *direction;
  bool steady_state = false;
  if(nrhs>0)
    {
    	const mxArray *mxa = prhs[0];
      int buflen=mxGetM(mxa) * mxGetN(mxa) + 1;
      char *first_argument;
      first_argument=(char*)mxCalloc(buflen, sizeof(char));
      int status = mxGetString(mxa, first_argument, buflen);
      if (status != 0)
        mexWarnMsgTxt("Not enough space. The first argument is truncated.");
      string f(first_argument);
      if(f == "steady_state")
        steady_state = true;
    }
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
  double * params = mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"params")));
  double *yd, *xd;
  if(!steady_state)
    {
      yd= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"endo_simul")));
      row_y=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"endo_simul")));
      col_y=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"endo_simul")));;
      xd= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_simul")));
      row_x=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_simul")));
      col_x=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_simul")));
      nb_row_xd=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"exo_det_nbr"))))));

      y_kmin=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"maximum_lag"))))));
      y_kmax=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"maximum_lead"))))));
      y_decal=max(0,y_kmin-int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"maximum_endo_lag")))))));
      periods=int(floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"periods"))))));
    }
	else
	  {
	  	yd= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"steady_state")));
      row_y=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"steady_state")));
      col_y=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"steady_state")));;
      xd= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_steady_state")));
      row_x=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_steady_state")));
      col_x=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,"exo_steady_state")));
      nb_row_xd=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"exo_det_nbr"))))));
	  }
	int maxit_=int(floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"maxit_"))))));
  double slowc=double(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"slowc")))));
  double markowitz_c=double(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"markowitz")))));
  double solve_tolf;
  if(steady_state)
    solve_tolf=*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"solve_tolf"))));
	else
    solve_tolf=*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,"dynatol"))));
  mxArray *mxa=mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,"fname"));
  int buflen=mxGetM(mxa) * mxGetN(mxa) + 1;
  char *fname;
  fname=(char*)mxCalloc(buflen, sizeof(char));
  string file_name=fname;
  int status = mxGetString(mxa, fname, buflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space. Filename is truncated.");


  int size_of_direction=col_y*row_y*sizeof(double);
  double * y=(double*)mxMalloc(size_of_direction);
  double * ya=(double*)mxMalloc(size_of_direction);
  direction=(double*)mxMalloc(size_of_direction);
  memset(direction,0,size_of_direction);
  double * x=(double*)mxMalloc(col_x*row_x*sizeof(double));
  for (i=0;i<row_x*col_x;i++)
     x[i]=double(xd[i]);
  for (i=0;i<row_y*col_y;i++)
    y[i]=double(yd[i]);
  int y_size=row_y;
  int nb_row_x=row_x;

  /*int it_ = y_kmin;
  for (int j = 0; j < y_size; j++)
		mexPrintf("   variable %d at time %d and %d = %f\n", j+1, it_, it_+1, y[j+it_*y_size]);*/

  clock_t t0= clock();
  Interpreter interprete(params, y, ya, x, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal, markowitz_c, file_name);
  string f(fname);
  bool result = interprete.compute_blocks(f, f, steady_state);
  clock_t t1= clock();
  if(!steady_state)
    mexPrintf("Simulation Time=%f milliseconds\n",1000.0*(double(t1)-double(t0))/double(CLOCKS_PER_SEC));
  if (nlhs>0)
    {
      plhs[0] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);
      pind = mxGetPr(plhs[0]);
      for (i=0;i<row_y*col_y;i++)
        pind[i]=y[i];
			if(nlhs>1)
			  {
			    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
			    pind = mxGetPr(plhs[1]);
			    if(result)
			      pind[0] = 0;
					else
					  pind[0] = 1;
			  }
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
