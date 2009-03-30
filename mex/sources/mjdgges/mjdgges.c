/*
 * Copyright (C) 2006-2008 Dynare Team
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

#include <string.h>
#include "mex.h"

#ifdef MWSIGNEDINDEX_NOT_DEFINED
#  define mwSignedIndex int
#endif

#ifdef NO_LAPACK_H
# if defined(__linux__) || defined(OCTAVE)
#  define dgges dgges_
# endif
void dgges(char *, char *, char *, int (*)(), int *, double *, int *, double *, int *, int *, double *, double *, double *,
           double *, int *, double *, int *, double *, int *, int *, int *);
#else /* NO_LAPACK_H */
# include "lapack.h"
#endif /* NO_LAPACK_H */

double criterium;

mwSignedIndex my_criteria(double *alphar, double *alphai, double *beta)
{
  return( (*alphar * *alphar + *alphai * *alphai) < criterium * *beta * *beta);
}

void mjdgges(double *a, double *b, double *z, double *n, double *sdim, double *eval_r, double *eval_i, double *info)
{
  mwSignedIndex i_n, i_info, i_sdim, one, lwork;
  double *alphar, *alphai, *beta, *work, *par, *pai, *pb, *per, *pei;
  double *junk;
  mwSignedIndex *bwork;

  one = 1;
  i_n = (mwSignedIndex)*n;
  alphar = mxCalloc(i_n,sizeof(double));
  alphai = mxCalloc(i_n,sizeof(double));
  beta = mxCalloc(i_n,sizeof(double));
  lwork = 16*i_n+16;
  work = mxCalloc(lwork,sizeof(double));
  bwork = mxCalloc(i_n,sizeof(mwSignedIndex));
  /* made necessary by bug in Lapack */
  junk = mxCalloc(i_n*i_n,sizeof(double));

  dgges("N", "V", "S", my_criteria, &i_n, a, &i_n, b, &i_n, &i_sdim, alphar, alphai, beta, junk, &i_n, z, &i_n, work, &lwork, bwork, &i_info);
  
  *sdim = i_sdim;
  *info = i_info;

  par = alphar;
  pai = alphai;
  pb = beta;
  pei = eval_i;
  for(per = eval_r; per <= &eval_r[i_n-1]; ++per)
    {
      *per = *par / *pb;
      *pei = *pai / *pb;
      ++par;
      ++pai;
      ++pb;
      ++pei;
    }
}

/* MATLAB interface */
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
  unsigned int m1,n1,m2,n2; 
  double *s, *t, *z, *sdim, *eval_r, *eval_i, *info, *a, *b; 
  double n;

  /* Check for proper number of arguments */
    
  if (nrhs < 2 || nrhs > 3) { 
    mexErrMsgTxt("MJDGGES: two or three input arguments are required."); 
  } else if (nlhs > 6) {
    mexErrMsgTxt("MJDGGES: too many output arguments."); 
  } 
    
  /* Check that A and B are real matrices of the same dimension.*/ 
    
  m1 = mxGetM(prhs[0]); 
  n1 = mxGetN(prhs[0]);
  m2 = mxGetM(prhs[1]); 
  n2 = mxGetN(prhs[1]);
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || 
      !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || 
      (m1 != n1) || (m2!= n1) || (m2 != n2)) { 
    mexErrMsgTxt("MJDGGES requires two square real matrices of the same dimension."); 
  } 
    
  /* Create a matrix for the return argument */ 
  plhs[0] = mxCreateDoubleMatrix(n1, n1, mxREAL); 
  plhs[1] = mxCreateDoubleMatrix(n1, n1, mxREAL); 
  plhs[2] = mxCreateDoubleMatrix(n1, n1, mxREAL); 
  plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL); 
  plhs[4] = mxCreateDoubleMatrix(n1, 1, mxCOMPLEX); 
  plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL); 
    
  /* Assign pointers to the various parameters */ 
  s = mxGetPr(plhs[0]);
  t = mxGetPr(plhs[1]);
  z = mxGetPr(plhs[2]);
  sdim = mxGetPr(plhs[3]);
  eval_r = mxGetPr(plhs[4]);
  eval_i = mxGetPi(plhs[4]);
  info = mxGetPr(plhs[5]);
    
  a = mxGetPr(prhs[0]);
  b = mxGetPr(prhs[1]);

  /* set criterium for stable eigenvalues */
  if ( nrhs == 3)
    {
      criterium = *mxGetPr(prhs[2]);
    }
  else
    {
      criterium = 1+1e-6;
    }

  /* keep a and b intact */
  memcpy(s,a,sizeof(double)*n1*n1);
  memcpy(t,b,sizeof(double)*n1*n1);

  n = n1;

  /* Do the actual computations in a subroutine */
  mjdgges(s, t, z, &n, sdim, eval_r, eval_i, info);


  return;
    
}

/*
07/30/03 MJ added user set criterium for stable eigenvalues
            corrected error messages in mexfunction()
*/
