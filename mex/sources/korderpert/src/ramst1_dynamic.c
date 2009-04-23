/*
 * ramst1_dynamic.c : Computes dynamic model for Dynare
 *
 * Warning : this file is generated automatically by Dynare
 *           from model file (.mod)

 */
#include <math.h>
#include "mex.h"
void Dynamic(double *y, double *x, int nb_row_x, double *params, int it_, double *residual, double *g1, double *g2)
{
  double lhs, rhs;

  /* Residual equations */
double
T19 = pow(y[3],params[0]-1),
T23 = 1+y[6]*params[0]*T19-params[1],
T29 = pow(y[0],params[0]),
T55 = (params[0]-1)*pow(y[3],params[0]-1-1),
T56 = y[6]*params[0]*T55,
T62 = params[0]*pow(y[0],params[0]-1);
lhs =1/y[2];
rhs =params[2]*1/y[5]*T23;
residual[0]= lhs-rhs;
lhs =y[2]+y[3];
rhs =y[4]*T29+y[0]*(1-params[1]);
residual[1]= lhs-rhs;
lhs =log(y[4]);
rhs =params[3]*log(y[1])+x[it_+0*nb_row_x];
residual[2]= lhs-rhs;
  /* Jacobian  */
  if (g1 == NULL)
    return;
  else
    {
  g1[6]=  g1[6]+(-1)/(y[2]*y[2]);
  g1[15]=  g1[15]+(-(T23*params[2]*(-1)/(y[5]*y[5])));
  g1[18]=  g1[18]+(-(params[2]*1/y[5]*params[0]*T19));
  g1[9]=  g1[9]+(-(params[2]*1/y[5]*T56));
  g1[7]=  g1[7]+1;
  g1[10]=  g1[10]+1;
  g1[13]=  g1[13]+(-T29);
  g1[1]=  g1[1]+(-(1-params[1]+y[4]*T62));
  g1[14]=  g1[14]+1/y[4];
  g1[5]=  g1[5]+(-(params[3]*1/y[1]));
  g1[23]=  g1[23]+(-1);
    }
  /* Hessian for endogenous and exogenous variables */
  if (g2 == NULL)
    return;
  else
    {
  g2[54] = (y[2]+y[2])/(y[2]*y[2]*y[2]*y[2]);
  g2[135] = (-(T23*params[2]*(y[5]+y[5])/(y[5]*y[5]*y[5]*y[5])));
  g2[159] = (-(params[2]*(-1)/(y[5]*y[5])*params[0]*T19));
  g2[87] = (-(params[2]*(-1)/(y[5]*y[5])*T56));
  g2[90] = (-(params[2]*1/y[5]*params[0]*T55));
  g2[81] = (-(params[2]*1/y[5]*y[6]*params[0]*(params[0]-1)*(params[0]-1-1)*pow(y[3],params[0]-1-1-1)));
  g2[13] = (-T62);
  g2[1] = (-(y[4]*params[0]*(params[0]-1)*pow(y[0],params[0]-1-1)));
  g2[110] = (-1)/(y[4]*y[4]);
  g2[29] = (-(params[3]*(-1)/(y[1]*y[1])));
  g2[138] = g2[159];
  g2[129] = g2[87];
  g2[153] = g2[90];
  g2[97] = g2[13];
    }
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *y, *x, *params;
  double *residual, *g1, *g2;
  int nb_row_x, it_;

  /* Create a pointer to the input matrix y. */
  y = mxGetPr(prhs[0]);

  /* Create a pointer to the input matrix x. */
  x = mxGetPr(prhs[1]);

  /* Create a pointer to the input matrix params. */
  params = mxGetPr(prhs[2]);

  /* Fetch time index */
  it_ = (int) mxGetScalar(prhs[3]) - 1;

  /* Gets number of rows of matrix x. */
  nb_row_x = mxGetM(prhs[1]);

  residual = NULL;
  if (nlhs >= 1)
  {
     /* Set the output pointer to the output matrix residual. */
     plhs[0] = mxCreateDoubleMatrix(3,1, mxREAL);
     /* Create a C pointer to a copy of the output matrix residual. */
     residual = mxGetPr(plhs[0]);
  }

  g1 = NULL;
  if (nlhs >= 2)
  {
     /* Set the output pointer to the output matrix g1. */
     plhs[1] = mxCreateDoubleMatrix(3, 8, mxREAL);
     /* Create a C pointer to a copy of the output matrix g1. */
     g1 = mxGetPr(plhs[1]);
  }

  g2 = NULL;
 if (nlhs >= 3)
  {
     /* Set the output pointer to the output matrix g2. */
     plhs[2] = mxCreateDoubleMatrix(3, 64, mxREAL);
     /* Create a C pointer to a copy of the output matrix g1. */
     g2 = mxGetPr(plhs[2]);
  }

  /* Call the C subroutines. */
  Dynamic(y, x, nb_row_x, params, it_, residual, g1, g2);
}
