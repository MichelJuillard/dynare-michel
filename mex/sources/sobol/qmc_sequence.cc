/*
** Computes Quasi Monte-Carlo sequence.
**
** Copyright (C) 2010-2012 Dynare Team
**
** This file is part of Dynare (can be used outside Dynare).
**
** Dynare is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** Dynare is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
**
** AUTHOR(S): stephane DOT adjemian AT univ DASH lemans DOT fr
**/

#include <string.h>
#include <stdint.h>
#include <dynmex.h>

#include "sobol.hh"
#include "gaussian.hh"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* 
  ** INPUTS:
  ** prhs[0]    [integer]    scalar, dimension.
  ** prhs[1]    [integer]    scalar, seed.
  ** prhs[2]    [integer]    scalar, sequence type:
  **                                   0 ==> uniform,
  **                                   1 ==> gaussian,
  **                                   2 ==> uniform on an hypershere.
  ** prhs[3]    [integer]    scalar, sequence size.
  ** prhs[4]    [double]     dimension*2 array, lower and upper bounds of the hypercube (default is 0-1 in all dimensions) if prhs[2]==0,
  **                         dimension*dimension array, covariance of the multivariate gaussian distribution of prhs[2]==1 (default is the identity matrix),
  **                         scalar, radius of the hypershere if prhs[2]==2 (default is one).
  **
  ** OUTPUTS:
  ** plhs[0]    [double]     sequence_size*dimension array, the Sobol sequence.
  ** plhs[1]    [integer]    scalar, seed.
  ** plhs[2]    [integer]    zero in case of success, one in case of error
  ** 
  */
  /*
  ** Check the number of input and output arguments.
  */
  if ( !( (nrhs==5) | (nrhs==4) | (nrhs==3) )  )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: Five, four or three input arguments are required!");
    }
  if (nlhs == 0)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: At least one output argument is required!");
    }
  /*
  ** Test the first input argument and assign it to dimension.
  */
  if (  !( mxIsNumeric(prhs[0]) ) )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: First input (dimension) has to be a positive integer!");
    }
  int dimension = ( int ) mxGetScalar(prhs[0]);
  /*
  ** Test the second input argument and assign it to seed.
  */
  if ( !( mxIsNumeric(prhs[1]) && mxIsClass(prhs[1],"int64") ) )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: Second input (seed) has to be an integer [int64]!");
    }
  int64_T seed = ( int64_T ) mxGetScalar( prhs[1] );
  /*
  ** Test the third input argument and assign it to type (kind of QMC sequence).
  */
  int error_flag_3 = 0;
  if (  !(mxIsNumeric(prhs[2])) )
    {
      error_flag_3 = 1;
    }
  int type = ( int ) mxGetScalar(prhs[2]);
  if ( !(type==0 || type==1 || type==2) )
    {
      error_flag_3 = 1;
    }
  if (error_flag_3==1)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: Third input (type of QMC sequence) has to be an integer equal to 0, 1 or 2!");
    }
  /*
  ** Test dimension>=2 when type==2
  */
  if ( (type==2) && (dimension<2) )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: First input (dimension) has to be greater than 1 for a uniform QMC on an hypershere!");
    }
  /*
  ** Test the optional fourth input argument and assign it to sequence_size.
  */
  if  ( ( nrhs>3 )  &&  !mxIsNumeric(prhs[3])  ) 
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: Fourth input (qmc sequence size) has to be a positive integer!");
    }
  int sequence_size;
  if ( nrhs>3)
    {
      sequence_size = ( int ) mxGetScalar( prhs[3] );
    }
  else
    {
      sequence_size = 1;
    }
  /*
  ** Test the optional fifth input argument and assign it to lower_and_upper_bounds.
  */
  if  (  ( nrhs>4 )  &&  (type==0) && ( !( mxGetN(prhs[4])==2) )  )// Sequence of uniformly distributed numbers in an hypercube.
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: The fifth input argument must be an array with two columns!");
    }
  if  ( (nrhs>4)   &&  (type==0) &&  ( !( (int)mxGetM(prhs[4])==dimension) ) )
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: The fourth input argument must be an array with a number of lines equal to dimension (first input argument)!");
    }
  if  ( ( nrhs>4 )  &&  (type==1) && ( !( ((int)mxGetN(prhs[4])==dimension) && ((int)mxGetM(prhs[4])==dimension) ) ) )// Sequence of normally distributed numbers.
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: The fifth input argument must be a squared matrix (whose dimension is given by the first input argument)!");
    }
  if  ( ( nrhs>4 )  &&  (type==2) && ( !( (mxGetN(prhs[4])==1) && (mxGetM(prhs[4])==1) ) ) )// Sequence of uniformly distributed numbers on an hypershere.
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("qmc_sequence:: The fifth input argument must be a positive scalar!");
    }
  double *lower_bounds = NULL, *upper_bounds = NULL;
  int unit_hypercube_flag = 1;
  if ( (type==0) && (nrhs>4) )
    {
      lower_bounds = (double *) mxCalloc(dimension,sizeof(double));
      upper_bounds = (double *) mxCalloc(dimension,sizeof(double));
      double *tmp;
      tmp = (double *) mxCalloc(dimension*2,sizeof(double));
      memcpy(tmp,mxGetPr(prhs[4]),dimension*2*sizeof(double));
      lower_bounds = &tmp[0];
      upper_bounds = &tmp[dimension];
      unit_hypercube_flag = 0;
    }
  double *cholcov;
  int identity_covariance_matrix = 1;
  if ( (type==1) && (nrhs>4) )
    {
      cholcov = (double *) mxCalloc(dimension*dimension,sizeof(double));
      double *tmp;
      tmp = (double *) mxCalloc(dimension*dimension,sizeof(double));
      memcpy(tmp,mxGetPr(prhs[4]),dimension*dimension*sizeof(double));
      cholcov = &tmp[0];
      identity_covariance_matrix = 0;
    }
  double radius = 1.0;
  int unit_radius = 1;
  if ( (type==2) && (nrhs>4) )
    {
      double *tmp;
      tmp = (double *) mxCalloc(1,sizeof(double));
      memcpy(tmp,mxGetPr(prhs[4]),sizeof(double));
      radius = tmp[0];
      unit_radius = 0;
    }
  /*
  ** Initialize outputs of the mex file.
  */
  double *qmc_draws;
  plhs[0] = mxCreateDoubleMatrix(dimension,sequence_size,mxREAL);
  qmc_draws = mxGetPr(plhs[0]);
  int64_T seed_out;

  if (sequence_size==1)
    {
      next_sobol ( dimension, &seed, qmc_draws );
      seed_out = seed;
    }
  else
    seed_out = sobol_block( dimension, sequence_size, seed, qmc_draws);

  if (type==0 && unit_hypercube_flag==0)// Uniform QMC sequence in an hypercube.
    expand_unit_hypercube( dimension, sequence_size, qmc_draws, lower_bounds, upper_bounds);
  else if (type==1)// Normal QMC sequance in R^n.
    {
      if (identity_covariance_matrix==1)
        icdfm(dimension*sequence_size, qmc_draws);
      else
        icdfmSigma(dimension,sequence_size, qmc_draws, cholcov);
    }
  else if (type==2)// Uniform QMC sequence on an hypershere.
    {
      if (unit_radius==1)
        usphere(dimension, sequence_size, qmc_draws);
      else
        usphereRadius(dimension, sequence_size, radius, qmc_draws);
    }
  
  if (nlhs >= 2)
    {
      plhs[1] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *((int64_T *) mxGetData(plhs[1])) = seed_out;
    }
  if (nlhs >= 3)
    plhs[2] = mxCreateDoubleScalar(0);
}
