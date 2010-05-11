

#include "swzmatrix.h"
#include "dw_ascii.h"
#include "dw_rand.h"
#include "dw_matrix_rand.h"
#include "dw_matrix_sort.h"

#include "dw_parse_cmd.h"

#include "switch.h"
#include "switchio.h"
#include "VARbase.h"
#include "VARio.h"
#include "mhm_VAR.h"

#include "spherical.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "modify_for_mex.h"

/*
   Returns ln(exp(a) + exp(b)) computed to avoid overflow.  If
   a = ln(c) and b = ln(d), as is usually the case, then the
   routine returns ln(c + d).

*/
static PRECISION AddLogs_static(PRECISION a, PRECISION b)
{
  return (a > b) ? a + log(1.0 + exp(b-a)) : b + log(exp(a-b) + 1.0);
}

/*
   Let h(x) and f(x) be probability density functions and let c be an unknown 
   constant.  In applications, the following will usually be true. 
      x      - parameters
      y      - data
      f(x)   - posterior distribution = p(x|y)
      h(x)   - proposal distribution
      c      - marginal distribution = p(y)
      c*f(x) - likelihood*prior = p(y|x)*p(x)

   Assumes:
     posterior : M x 2 matrix with posterior[i][0] = ln(h(x(i))) and 
                 posterior[i][1] = ln(c*f(x(i))) where x(i) is sampled from f(x).

   Returns:
     Estimate of c.

   Notes:
     Uses the fact that c = 1/I(L) where I(L) is the integral over x of

                        h(x)/(c*f(x)) * f(x) 

     I(L) can be computed from the posterior draws.
    
*/
PRECISION ComputeMarginalDensity_Standard(TMatrix posterior)
{
  PRECISION I=MINUS_INFINITY;
  int i;

  for (i=RowM(posterior)-1; i >= 0; i--)
    I=AddLogs_static(I,ElementM(posterior,i,0) - ElementM(posterior,i,1));

  I-=log((PRECISION)RowM(posterior));

  return -I;
}


/*
   Let h(x) and f(x) be probability density functions and let c be an unknown 
   constant.  In applications, the following will usually be true. 
      x      - parameters
      y      - data
      f(x)   - posterior distribution = p(x|y)
      h(x)   - proposal distribution
      c      - marginal distribution = p(y)
      c*f(x) - likelihood*prior = p(y|x)*p(x)

   Assumes:
     posterior : M x 2 matrix with posterior[i][0] = ln(h(x(i))) and 
                 posterior[i][1] = ln(c*f(x(i))) where x(i) is sampled from f(x).
     L         : cutoff value in logs.
     q         : estimate of the probability that x sampled from h(x) satisfies
                 c*f(x) > exp(L).

   Returns:
     Estimate of c or MINUS_INFINITY if no proposal draws satisfied the
     restriction given by the cutoff value L.

   Notes:
     Uses the fact that c = q/I(L) where I(L) is the integral over x of

                      1{c*f(x)>exp(L)}*h(x)/(c*f(x)) * f(x) 

     I(L) can be computed from the posterior draws.
    
*/
PRECISION ComputeMarginalDensity_WZ1_q(TMatrix posterior, PRECISION L, PRECISION q, int *in_posterior)
{
  PRECISION I=MINUS_INFINITY;
  int i;

  for (*in_posterior=0, i=RowM(posterior)-1; i >= 0; i--)
    if (ElementM(posterior,i,1) >= L)
      {
	(*in_posterior)++;
	I=AddLogs_static(I,ElementM(posterior,i,0) - ElementM(posterior,i,1));
      }

  if ((*in_posterior) > 0) I-=log((PRECISION)RowM(posterior));

  return (q == 0.0) ? MINUS_INFINITY : log(q) - I;
}

/*
   Let h(x) and f(x) be probability density functions and let c be an unknown 
   constant.  In applications, the following will usually be true. 
      x      - parameters
      y      - data
      f(x)   - posterior distribution = p(x|y)
      h(x)   - proposal distribution
      c      - marginal distribution = p(y)
      c*f(x) - likelihood*prior = p(y|x)*p(x)

   Assumes:
     proposal  : N x 2 matrix with proposal[i][0] = ln(h(x(i))) and 
                 proposal[i][1] = ln(c*f(x(i))) where x(i) is sampled from h(x).
     posterior : M x 2 matrix with posterior[i][0] = ln(h(x(i))) and 
                 posterior[i][1] = ln(c*f(x(i))) where x(i) is sampled from f(x).
     L         : cutoff value.

   Returns:
     Estimate of c or MINUS_INFINITY if no proposal draws satisfied the
     restriction given by the cutoff value L.

   Notes:
     Uses the fact that c = P(L)/I(L) where p(L) is the probability that x 
     sampled from h(x) satisfies c*f(x) > exp(L) and I(L) is the integral over x 
     of

                      1{c*f(x)>exp(L)}*h(x)/(c*f(x)) * f(x) 

     P(L) can be computed from the proposal draws and I(L) can be computed from 
     the posterior draws.
    
*/
PRECISION ComputeMarginalDensity_WZ1(TMatrix proposal, TMatrix posterior, PRECISION L, int *in_proposal, int *in_posterior)
{
  PRECISION I=MINUS_INFINITY;
  int i;

  for (*in_proposal=0, i=RowM(proposal)-1; i >= 0; i--)
    if (ElementM(proposal,i,1) >= L) (*in_proposal)++;

  for (*in_posterior=0, i=RowM(posterior)-1; i >= 0; i--)
    if (ElementM(posterior,i,1) >= L)
      {
	(*in_posterior)++;
	I=AddLogs_static(I,ElementM(posterior,i,0) - ElementM(posterior,i,1));
      }

  if ((*in_posterior) > 0) I-=log((PRECISION)RowM(posterior));

  return ((*in_proposal) == 0) ? MINUS_INFINITY : log((PRECISION)(*in_proposal)/(PRECISION)RowM(proposal)) - I;
}

/*
   ess - effective sample size
    This is the sum of h(x)/(c*f(x)) divided by the max of h(x)/(c*f(x)).  See 
    the header for ComputeMarginalDensity_WZ1() for the definitions of h(x) and c*f(x).
*/
PRECISION ComputeEffectiveSampleSize_WZ1(TMatrix posterior, PRECISION L)
{
  PRECISION sum=MINUS_INFINITY, max=MINUS_INFINITY, tmp;
  int i;

  for (i=RowM(posterior)-1; i >= 0; i--)
    if (ElementM(posterior,i,1) >= L)
      {
	tmp=ElementM(posterior,i,0) - ElementM(posterior,i,1);
	sum=AddLogs_static(sum,tmp);
	if (tmp > max) max=tmp;
      }

  return exp(sum - max);
}


/*
   Let h(x) and f(x) be probability density functions and let c be an unknown 
   constant.  In applications, the following will usually be true. 
      x      - parameters
      y      - data
      f(x)   - posterior distribution = p(x|y)
      h(x)   - proposal distribution
      c      - marginal distribution = p(y)
      c*f(x) - likelihood*prior = p(y|x)*p(x)

   Assumes:
     proposal  : N x 2 matrix with proposal[i][0] = ln(h(x(i))) and 
                 proposal[i][1] = ln(c*f(x(i))) where x(i) is sampled from h(x).
     posterior : M x 2 matrix with posterior[i][0] = ln(h(x(i))) and 
                 posterior[i][1] = ln(c*f(x(i))) where x(i) is sampled from f(x).
     L1        : cutoff value (likelihood(x)*prior(x) > L1)
     L2        : cutoff value for (proposal(x) < L2)

   Returns:
     Estimate of c or MINUS_INFINITY if no proposal draws satisfied the
     restriction given by the cutoff values of L1 and L2.

   Notes:
     Uses the fact that c = P(L1,L2)/I(L1,L2) where p(L1,L2) is the probability
     that x sampled from h(x) satisfies c*f(x) > exp(L1) and h(x) < exp(L2) and 
     I(L1,L2) is the integral over x of

              1{c*f(x)>exp(L1) and h(x)<exp(L2)}*h(x)/(c*f(x)) * f(x) 

     P(L) can be computed from the proposal draws and I(L) can be computed from 
     the posterior draws.  The function 1{.} denotes the 
    
*/
PRECISION ComputeMarginalDensity_WZ2(TMatrix proposal, TMatrix posterior, PRECISION L1, PRECISION L2, int *in_proposal, int *in_posterior)
{
  PRECISION I=MINUS_INFINITY;
  int i;

  for (*in_proposal=0, i=RowM(proposal)-1; i >= 0; i--)
    if ((ElementM(proposal,i,1) >= L1) && (ElementM(proposal,i,0) <= L2))
      (*in_proposal)++;

  for (*in_posterior=0, i=RowM(posterior)-1; i >= 0; i--)
    if ((ElementM(posterior,i,1) >= L1) && (ElementM(posterior,i,0) <= L2))
      {
	(*in_posterior)++;
	I=AddLogs_static(I,ElementM(posterior,i,0) - ElementM(posterior,i,1));
      }

  if ((*in_posterior) > 0) I-=log((PRECISION)RowM(posterior));

  return ((*in_proposal) == 0) ? MINUS_INFINITY : log((PRECISION)(*in_proposal)/(PRECISION)RowM(proposal)) - I;
}

PRECISION ComputeMarginalDensity_WZ3(TMatrix proposal, TMatrix posterior, PRECISION L1, PRECISION L2, int *in)
{
  PRECISION I=MINUS_INFINITY, tmp;
  int i;

  for (i=RowM(proposal)-1; i >= 0; i--)
    {
      tmp=ElementM(proposal,i,0) - ElementM(proposal,i,1);
      if ((L1 < tmp) && (tmp < L2)) (*in)++;
    }

  for (i=RowM(posterior)-1; i >= 0; i--)
    {
      tmp=ElementM(posterior,i,0) - ElementM(posterior,i,1);
      if ((L1 < tmp) && (tmp < L2))
	I=AddLogs_static(I,tmp);
    }

  if (I > MINUS_INFINITY) I-=log((PRECISION)RowM(posterior));

  return ((*in) == 0) ? MINUS_INFINITY : log((PRECISION)(*in)/(PRECISION)RowM(proposal)) - I;
}

/*******************************************************************************/
/************************ Mueller's Method (alternate ) ************************/
/*******************************************************************************/
#define MAX_C 1E50
#define TOL   1E-5
PRECISION ComputeDifference(TMatrix proposal, TMatrix posterior, PRECISION c, int *in1, int *in2)
{
  int i, n;
  PRECISION sum1=0.0, sum2=0.0, tmp;
  for (*in1=0, i=RowM(proposal)-1; i >= 0; i--)
    if (c < (tmp=ElementM(proposal,i,0)-ElementM(proposal,i,1)))
      {
	sum1+=1-exp(c-tmp);
	(*in1)++;
      }
  for (*in2=0, i=RowM(posterior)-1; i >= 0; i--)
    if (c > (tmp=ElementM(posterior,i,0)-ElementM(posterior,i,1)))
      {
	sum2+=1-exp(tmp-c);
	(*in2)++;
      }
  return sum1/(PRECISION)RowM(proposal) - sum2/(PRECISION)RowM(posterior);
}
PRECISION ComputeMarginalDensity_Mueller_Check(TMatrix proposal, TMatrix posterior,int *in1, int *in2)
{
  PRECISION min_c, max_c, mid_c, diff;
  int i;

  if ((diff=ComputeDifference(proposal,posterior,mid_c=0.0,in1,in2)) < 0.0)
    {
printf("%lf %le - %d %d\n",mid_c,diff,*in1,*in2);                                            //**
      max_c=mid_c;
      for (min_c=-1.0; min_c > -MAX_C; max_c=min_c, min_c*=10)
	if ((diff=ComputeDifference(proposal,posterior,min_c,in1,in2)) > 0) break;
else printf("%lf %lf %le - %d %d\n",min_c,max_c,diff,*in1,*in2);                             //**
      if (min_c <= -MAX_C) return min_c;
    }
  else
    {
printf("%lf %le - %d %d\n",mid_c,diff,*in1,*in2);                                            //**
      min_c=mid_c;
      for (max_c=1.0; max_c < MAX_C; min_c=max_c, max_c*=10)
	if ((diff=ComputeDifference(proposal,posterior,max_c,in1,in2)) < 0) break;
else printf("%lf %lf %le - %d %d\n",min_c,max_c,diff,*in1,*in2);                             //**
      if (max_c >= MAX_C) return max_c;
    }
printf("%lf %lf %le - %d %d\n",min_c,max_c,diff,*in1,*in2);                                  //**
  diff=ComputeDifference(proposal,posterior,mid_c=(min_c + max_c)/2.0,in1,in2);
  for (i=0; i < 200; i++)
    {
      if (diff > 0)
	min_c=mid_c;
      else
	max_c=mid_c;
printf("%lf %lf %le %d - %d %d\n",min_c,max_c,diff,i,*in1,*in2);                             //**
      if ((fabs(diff=ComputeDifference(proposal,posterior,mid_c=(min_c + max_c)/2.0,in1,in2)) < TOL) && (i > 20)) break;
    }
printf("%lf %le - %d %d\n",mid_c,diff,*in1,*in2);                                            //**
  return -mid_c;
}
#undef MAX_C
#undef TOL
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/****************************** Mueller's Method *******************************/
/*******************************************************************************/
#define MAX_C 1E50
//#define VERBOSE
PRECISION ComputeLinear(TMatrix proposal, PRECISION log_c, int *intercept)
{
  int i;
  PRECISION slope=0.0, tmp;
  for (*intercept=0, i=RowM(proposal)-1; i >= 0; i--)
    if (log_c < (tmp=ElementM(proposal,i,0) - ElementM(proposal,i,1)))
      {
	(*intercept)++;
	slope+=exp(log_c-tmp);
      }
  return slope;
}
PRECISION ComputeInverseLinear(TMatrix posterior, PRECISION log_c, int *intercept)
{
  int i;
  PRECISION slope=0.0, tmp;
  for (*intercept=0, i=RowM(posterior)-1; i >= 0; i--)
    if (log_c > (tmp=ElementM(posterior,i,0)-ElementM(posterior,i,1)))
      {
	(*intercept)++;
	slope+=exp(tmp-log_c);
      }
  return slope;
}
PRECISION ComputeMarginalDensity_Mueller(TMatrix proposal, TMatrix posterior,int *in1, int *in2)
{
  PRECISION log_c=0.0, slope1, slope2, N1=1.0/(PRECISION)RowM(proposal), N2=1.0/(PRECISION)RowM(posterior),
    intercept1, intercept2, max_log_c, min_log_c, diff, tmp;
  int min_in1, max_in1, min_in2, max_in2;

  slope1=ComputeLinear(proposal,log_c,in1);
  slope2=ComputeInverseLinear(posterior,log_c,in2);
  diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);

  // Bracket the intersection
  if (diff < 0.0)
    {
      do
	if (log_c < -MAX_C) 
	  return log_c;
	else
	  {
	    max_in1=*in1;
	    max_in2=*in2;
	    max_log_c=log_c;
	    log_c=10*(log_c-1);
	    slope1=ComputeLinear(proposal,log_c,in1);
	    slope2=ComputeInverseLinear(posterior,log_c,in2);
	    diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);
#ifdef VERBOSE
printf("%lf %lf %le - %d %d\n",log_c,max_log_c,diff,*in1,*in2);                                             
#endif
	  }
      while (diff < 0.0);
      min_in1=*in1;
      min_in2=*in2;
      min_log_c=log_c;
    }
  else
    {
      do
	if (log_c > MAX_C) 
	  return log_c;
	else
	  {
	    min_in1=*in1;
	    min_in2=*in2;
	    min_log_c=log_c;
	    log_c=10*(log_c+1);
	    slope1=ComputeLinear(proposal,log_c,in1);
	    slope2=ComputeInverseLinear(posterior,log_c,in2);
	    diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);
#ifdef VERBOSE
printf("%lf %lf %le - %d %d\n",min_log_c,log_c,diff,*in1,*in2);                                            
#endif
	  }
      while (diff >= 0.0);
      max_in1=*in1;
      max_in2=*in2;
      max_log_c=log_c;
    }

  // At this point diff(min_log_c) >= 0 and diff(max_log_c) < 0.
  while ((min_in1 != max_in1) || (min_in2 != max_in2))
    {
      log_c=(min_log_c + max_log_c)/2.0;
      slope1=ComputeLinear(proposal,log_c,in1);
      slope2=ComputeInverseLinear(posterior,log_c,in2);
      diff=N1*((PRECISION)(*in1) - slope1) - N2*((PRECISION)(*in2) - slope2);
      if (diff > 0)
	{
	  min_in1=*in1;
	  min_in2=*in2;
	  min_log_c=log_c;
	}
      else
	{
	  max_in1=*in1;
	  max_in2=*in2;
	  max_log_c=log_c;
	}
#ifdef VERBOSE
printf("%lf %lf %le - %d %d\n",min_log_c,max_log_c,diff,*in1,*in2);                                        
#endif
    }

  slope1=N1*ComputeLinear(proposal,min_log_c,in1);
  intercept1=N1*(PRECISION)(*in1);
  slope2=N2*ComputeInverseLinear(posterior,min_log_c,in2);
  intercept2=N2*(PRECISION)(*in2);

  tmp=intercept1-intercept2;
  if (slope1 > 0)
    {
      tmp+=sqrt(tmp*tmp + 4*slope1*slope2);
      if (tmp >= 2.0*slope1)
	tmp=min_log_c + log(tmp) - log(2.0*slope1);
      else
	return -min_log_c;
    }
  else
    {
#ifdef VERBOSE
printf("Flat linear slope\n");                                       
#endif
    if (tmp > 0)
      if (slope2 > tmp)
	tmp=min_log_c + log(slope2) - log(tmp);
      else
	return -min_log_c;
    else
      return -max_log_c;
    }
  return (tmp > max_log_c) ? -max_log_c : -tmp;
}

void ComputeDiagnostics_Mueller(TMatrix proposal, TMatrix posterior, PRECISION log_c, 
				PRECISION *log_slope1, PRECISION *ess1, PRECISION *log_slope2, PRECISION *ess2)
{
  PRECISION slope, sum, max, tmp;
  int i;

  // proposal - piecewise linear
  sum=max=slope=0.0;
  for (i=RowM(proposal)-1; i >= 0; i--)
    if (log_c < (tmp=ElementM(proposal,i,0) - ElementM(proposal,i,1)))
      {
	tmp=exp(log_c-tmp);
	slope+=tmp;
	sum+=1.0-tmp;
	if (1.0-tmp > max) max=1.0-tmp;
      }

  *log_slope1=log(slope) - log_c - log(RowM(proposal));
  *ess1=(max > 0.0) ? sum/max : 0.0;

  // posterior - piecewise hyperbolic
  sum=max=slope=0.0;
  for (i=RowM(posterior)-1; i >= 0; i--)
    if (log_c > (tmp=ElementM(posterior,i,0)-ElementM(posterior,i,1)))
      {
	tmp=exp(tmp-log_c);
	slope-=tmp;
	sum+=1.0-tmp;
	if (1.0-tmp > max) max=1.0-tmp;
      }
 
  *log_slope2=log(slope) - log_c - log(RowM(posterior));
  *ess2=(max > 0.0) ? sum/max : 0.0;
} 
#undef MAX_C
#undef VERBOSE
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
TVector Sorted_Level_Cuts(TMatrix X, int n_cuts)
{
  int step, start, i;
  PRECISION min, max, inc, tmp;
  TVector cuts=CreateVector(n_cuts), Y=ColumnVector((TVector)NULL,X,1);

  SortVectorAscending(Y,Y);

  step=RowM(X)/(n_cuts-1);
  start=(RowM(X) - (n_cuts-1)*step)/2;

  for (i=0; i < n_cuts; i++)
    ElementV(cuts,i)=ElementV(Y,start+i*step);

  FreeVector(Y);
  return cuts;
}

TVector Sorted_Center_Cuts(TMatrix X, int n_cuts)
{
  int step, start, i;
  PRECISION min, max, inc, tmp;
  TVector cuts=CreateVector(n_cuts), Y=ColumnVector((TVector)NULL,X,1);

  step=RowM(X)/(n_cuts-1);
  start=(RowM(X) - (n_cuts-1)*step)/2;

  SortVectorAscending(Y,Y);
  for (i=0; i < n_cuts; i++)
    ElementV(cuts,i)=ElementV(Y,start+i*step);

  FreeVector(Y);
  return cuts;
}


TVector Linear_Level_Cuts(TMatrix X, int n_cuts)
{
  int k;
  PRECISION min, max, inc, tmp;
  TVector cuts=CreateVector(n_cuts);

  max=min=ElementM(X,0,1);
  for (k=RowM(X)-1; k >= 0; k--)
    {
      tmp=ElementM(X,k,1);
      if (tmp < min)
	min=tmp;
      else
	if (tmp > max)
	  max=tmp;
    }
  min-=0.1;
  max+=0.1;

  ElementV(cuts,0)=min;
  inc=(max - min)/(n_cuts-1);

  for (k=1; k < n_cuts; k++)
    ElementV(cuts,k)=ElementV(cuts,k-1)+inc;

  return cuts;
}

TVector Exponential_Level_Cuts(TMatrix X, int n_cuts)
{
  int k;
  PRECISION min, max, inc, tmp;
  TVector cuts=CreateVector(n_cuts);

  max=min=ElementM(X,0,1);
  for (k=RowM(X)-1; k >= 0; k--)
    {
      tmp=ElementM(X,k,1);
      if (tmp < min)
	min=tmp;
      else
	if (tmp > max)
	  max=tmp;
    }

  ElementV(cuts,0)=min;
  inc=max+log(1-exp(min-max))-log(n_cuts-1);

  for (k=1; k < n_cuts; k++)
    ElementV(cuts,k)=AddLogs_static(ElementV(cuts,k-1),inc);

  return cuts;
}

/*
   Assumes:
     center : center point for spherical distribution
     scale  : linear transformation for spherical distribution
     model  : valid point to TStateModel structure
     alpha  : array of parameters of independent draws from the Dirichlet 
              distribution

   Results:
     Draws theta from affine transformation of spherical distribution and Q from 
     independent Dirichlet distribution.  Computes the scaled log posterior 
     (likelihood times prior), with the states integerated out.  Returns the
     vector q.  q[i] is the proportion of the draws that are greater than or 
     equal to level_cuts[i].
*/
TVector Create_q(int ndraws, TVector center, TMatrix scale, TStateModel* model, TVector* alpha, TVector level_cuts)
{
  TVector free_parameters, q;
  PRECISION density;
  int i, j, k, begin_time, end_time;

  // VAR specific
  PRECISION *zeta;
  int n_zeta;

  InitializeVector(q=CreateVector(DimV(level_cuts)),0.0);
  free_parameters=CreateVector(NumberFreeParametersTheta(model));

  // VAR specific
  zeta=pElementV(free_parameters) + ZetaIndex((T_VAR_Parameters*)(model->theta));
  n_zeta=ZetaLength((T_VAR_Parameters*)(model->theta));

  // timings
  begin_time=time((time_t*)NULL);

  for (i=ndraws-1; i >= 0; i--)
    {
      // draw uniform on sphere
      DrawSpherical(free_parameters);

      // scale and center to get draw for theta
      ProductMV(free_parameters,scale,free_parameters);
      AddVV(free_parameters,free_parameters,center);

      // draw Q from Dirichet
      DrawIndependentDirichletVector(model->sv->ba,alpha);

      // VAR specific - perform checks on theta parameters to ensure they are valid
      for (j=n_zeta-1; j >= 0; j--)
	if (zeta[j] <= 0)
	  {
	    // non-positive values of zeta - log posterior is minus infinity
	    break;
	  }

      if (j < 0)
	{
	  // force free parameters into model
	  ConvertFreeParametersToTheta(model,pElementV(free_parameters));
	  Update_Q_from_B_SV(model->sv);
	  if (!(model->sv->valid_transition_matrix)) ValidateTransitionMatrices_SV(model->sv);
	  TransitionMatricesChanged(model);

	  // Normalize
	  if (!IsNormalized_VAR(model->theta))
	    {
	      // parameters changed - log posterior is minus infinity
	      break;
	    }
	  else
	    {
	      // compute log posterior
	      density=LogPosterior_StatesIntegratedOut(model);
	      for (k=DimV(level_cuts)-1; k >= 0; k--)
		if (density >= ElementV(level_cuts,k))
		  ElementV(q,k)+=1.0;
	    }
	}
    }

  for (k=DimV(level_cuts)-1; k >= 0; k--)
    ElementV(q,k)/=(PRECISION)ndraws;

  // timings
  end_time=time((time_t*)NULL);
  printf("Elapsed Time: %d seconds\n",end_time - begin_time);

  FreeVector(free_parameters);

  return q;
}

/*
   Assumes:
     center : center point for spherical distribution
     scale  : linear transformation for spherical distribution
     model  : valid pointer to TStateModel structure
     alpha  : array of parameters of independent draws from the Dirichlet 
              distribution

   Results:
     Draws theta from affine transformation of spherical distribution and Q from 
     independent Dirichlet distribution.  Computes the log density, properly 
     scaled, of the proposal distribution and stores this in the first column of 
     proposal.  Computes the scaled log posterior (likelihood time prior), with 
     the states integerated out, and stores this in the second column of 
     proposal.
*/
TMatrix CreateProposal(int ndraws, TVector center, TMatrix scale, TStateModel* model, TVector* alpha)
{
  TMatrix proposal;
  TVector free_parameters;
  PRECISION r, Jacobian;
  int i, j, begin_time, end_time;

  // VAR specific
  PRECISION *zeta;
  int n_zeta;

  proposal=CreateMatrix(ndraws,2);
  free_parameters=CreateVector(NumberFreeParametersTheta(model));
  Jacobian=-LogAbsDeterminant_LU(scale);

  // VAR specific
  zeta=pElementV(free_parameters) + ZetaIndex((T_VAR_Parameters*)(model->theta));
  n_zeta=ZetaLength((T_VAR_Parameters*)(model->theta));

  // timings
  begin_time=time((time_t*)NULL);

  for (i=ndraws-1; i >= 0; i--)
    {
      // draw uniform on sphere
      r=DrawSpherical(free_parameters);

      // scale and center to get draw for theta
      ProductMV(free_parameters,scale,free_parameters);
      AddVV(free_parameters,free_parameters,center);

      // draw Q from Dirichet
      DrawIndependentDirichletVector(model->sv->ba,alpha);

      // compute log proposal density
      ElementM(proposal,i,0)=LogSphericalDensity(r) + Jacobian + LogIndependentDirichlet_pdf(model->sv->ba,alpha);

      // VAR specific - perform checks on theta parameters to ensure they are valid
      for (j=n_zeta-1; j >= 0; j--)
	if (zeta[j] <= 0)
	  {
	    ElementM(proposal,i,1)=MINUS_INFINITY;
	    break;
	  }

      if (j < 0)
	{
	  // force free parameters into model
	  ConvertFreeParametersToTheta(model,pElementV(free_parameters));
	  Update_Q_from_B_SV(model->sv);
	  if (!(model->sv->valid_transition_matrix)) ValidateTransitionMatrices_SV(model->sv);
	  TransitionMatricesChanged(model);

	  // compute log posterior
	  ElementM(proposal,i,1)=LogPosterior_StatesIntegratedOut(model);
	}
    }

  // timings
  end_time=time((time_t*)NULL);
  printf("Elapsed Time: %d seconds\n",end_time - begin_time);

  FreeVector(free_parameters);

  return proposal;
}

/*
   Assumes:
     center : center point for spherical distribution
     scale  : linear transformation for spherical distribution
     model  : valid point to TStateModel structure
     alpha  : array of parameters of independent draws from the Dirichlet 
              distribution

   Results:
     Draws theta from affine transformation of spherical distribution and Q from 
     independent Dirichlet distribution.  Computes the log density, properly 
     scaled, of the proposal distribution and stores this in the first column of 
     proposal.  Computes the scaled log posterior (likelihood time prior), with 
     the states integerated out, and stores this in the second column of 
     proposal.
*/
TMatrix CreateProposal_Radius(int ndraws, TVector center, TMatrix scale, TStateModel* model, TVector* alpha)
{
  TMatrix proposal;
  TVector free_parameters;
  PRECISION r, Jacobian;
  int i, j, begin_time, end_time;

  // VAR specific
  PRECISION *zeta;
  int n_zeta;

  proposal=CreateMatrix(ndraws,2);
  free_parameters=CreateVector(NumberFreeParametersTheta(model));
  Jacobian=-LogAbsDeterminant_LU(scale);

  // VAR specific
  zeta=pElementV(free_parameters) + ZetaIndex((T_VAR_Parameters*)(model->theta));
  n_zeta=ZetaLength((T_VAR_Parameters*)(model->theta));

  // timings
  begin_time=time((time_t*)NULL);

  for (i=ndraws-1; i >= 0; i--)
    {
      // draw uniform on sphere
      r=DrawSpherical(free_parameters);

      // scale and center to get draw for theta
      ProductMV(free_parameters,scale,free_parameters);
      AddVV(free_parameters,free_parameters,center);

      // draw Q from Dirichet
      DrawIndependentDirichletVector(model->sv->ba,alpha);

      // save radius
      ElementM(proposal,i,0)=r;

      // VAR specific - perform checks on theta parameters to ensure they are valid
      for (j=n_zeta-1; j >= 0; j--)
	if (zeta[j] <= 0)
	  {
	    ElementM(proposal,i,1)=MINUS_INFINITY;
	    continue;
	  }

      // force free parameters into model
      ConvertFreeParametersToTheta(model,pElementV(free_parameters));
      Update_Q_from_B_SV(model->sv);
      if (!(model->sv->valid_transition_matrix)) ValidateTransitionMatrices_SV(model->sv);
      TransitionMatricesChanged(model);

      // compute log posterior
      ElementM(proposal,i,1)=LogPosterior_StatesIntegratedOut(model);
    }

  // timings
  end_time=time((time_t*)NULL);
  printf("Elapsed Time: %d seconds\n",end_time - begin_time);

  FreeVector(free_parameters);

  return proposal;
}

/*
    Assumes:
      X is a m x k matrix with k > 2+idx_alpha
        Column 1   = log posterior density of (theta,Q) (not properly scaled).
        Column 2   = (theta - center)' * Inverse(scale * scale') * (theta - center).
        Column 2+i = log Dirichlet density of Q.

    Returns:
      m x 2 matrix.  The first column is the log of the proposal density of the
      posterior draw and the second column is the log of the posterior density.

    Notes:
      The parameters are theta and Q.  The proposal density is independent across
      these with a Dirichlet distribution on Q and a spherical density on theta.
      The scale for spherical distribution is (rescale_factor * base_scale). 
*/
TMatrix CreatePosterior(TMatrix X, PRECISION rescale_factor, TMatrix base_scale, int idx_alpha, int dim)
{
  int i;
  PRECISION factor=1.0/rescale_factor, Jacobian=-LogAbsDeterminant_LU(base_scale) - dim*log(rescale_factor);
  TMatrix Y=CreateMatrix(RowM(X),3);

  for (i=RowM(X)-1; i >= 0; i--)
    {
      ElementM(Y,i,0)=LogSphericalDensity(factor*sqrt(ElementM(X,i,1))) + Jacobian + ElementM(X,i,2+idx_alpha);
      ElementM(Y,i,1)=ElementM(X,i,0);
      ElementM(Y,i,2)=ElementM(X,i,5);
    }

  return Y;
}

void ComputeMarginal(int ndraws_proposal, TMatrix X, TStateModel *model, T_MHM *mhm, int idx_alpha, PRECISION rescale_factor, 
		     char *tag, char *proposal_tag)
{
  TMatrix proposal, posterior, base_scale, scale;
  TVector *alpha=(TVector*)NULL, level_cuts;
  PRECISION diff, marginal_density, ess, slope1, ess1, slope2, ess2;
  int i, in1, in2;
  char filename[256];
  FILE *f_out;

  // compute base scale and scale for proposal
  base_scale=CholeskyUT((TMatrix)NULL,mhm->variance);
  Transpose(base_scale,base_scale);
  scale=ProductSM((TMatrix)NULL,rescale_factor,base_scale);

  // allocate alpha
  alpha=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    alpha[i]=CreateVector(DimV(mhm->BaseAlpha[i]));

  // Dirichlet index
  idx_alpha=1;
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    ProductVS(alpha[i],mhm->BaseAlpha[i],ElementV(mhm->alpha_scales,idx_alpha));

  // Create posterior matrix
  printf("Creating posterior (%d draws)\n",RowM(X));
  posterior=CreatePosterior(X,rescale_factor,base_scale,idx_alpha,NumberFreeParametersTheta(model));

  // Create proposal matrix
  printf("Creating proposal (%d draws)\n",ndraws_proposal);
  proposal=CreateProposal(ndraws_proposal,mhm->center,scale,model,alpha);

  // Create output file
  sprintf(filename,"%s_md_%s.dat",proposal_tag,tag);
  f_out=dw_AppendTextFile(filename);

  // Compute marginal density using standard technique
  marginal_density=ComputeMarginalDensity_Standard(posterior);
  fprintf(f_out,"Standard MHM method: %0.14lg\n",marginal_density);

  // Compute marginal density using WZ1
  level_cuts=Sorted_Level_Cuts(posterior,10);
  fprintf(f_out,"WZ method - no center cut\nmarginal density,level,number proposal draws,number posterior draws,effective sample size\n");
  for (i=0; i < DimV(level_cuts); i++)
    {
      marginal_density=ComputeMarginalDensity_WZ1(proposal,posterior,ElementV(level_cuts,i),&in1,&in2);
      ess=ComputeEffectiveSampleSize_WZ1(posterior,ElementV(level_cuts,i));

      fprintf(f_out,"%0.14lg,%lg,%d,%d,%lg\n",marginal_density,ElementV(level_cuts,i),in1,in2,ess);
    }
  FreeVector(level_cuts);

  // Compute marginal density Mueller
  fprintf(f_out,"\nMueller method\nmarginal density,difference,# proposal draws,proposal log slope,proposal ess,# posterior draws,posterior log slope,posterior ess\n");
  marginal_density=ComputeMarginalDensity_Mueller(proposal,posterior,&in1,&in2);
  ComputeDiagnostics_Mueller(proposal,posterior,-marginal_density,&slope1,&ess1,&slope2,&ess2);
  diff=ComputeDifference(proposal,posterior,-marginal_density,&in1,&in2);
  fprintf(f_out,"%0.14lg,%lg,%d,%lg,%lg,%d,%lg,%lg\n\n",marginal_density,diff,in1,slope1,ess1,in2,slope2,ess2);

  fclose(f_out);

  sprintf(filename,"%s_md_proposal_%s.dat",proposal_tag,tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,proposal,"%le,");
  fclose(f_out);

  sprintf(filename,"%s_md_posterior_%s.dat",proposal_tag,tag);
  f_out=dw_CreateTextFile(filename);
  //dw_PrintMatrix(f_out,posterior,"%le,");
  fclose(f_out);

  FreeMatrix(posterior);
  FreeMatrix(proposal);

  dw_FreeArray(alpha);
  FreeMatrix(base_scale);
  FreeMatrix(scale);
}


/*******************************************************************************/
PRECISION SetupSphericalFromPosterior_Table(TMatrix X, int n)
{
  TVector table, Y;
  PRECISION cut_point=0.01, p, inc;
  int n_table=100, i, j, B, N;

  Y=ColumnVector((TVector)NULL,X,1);
  SortVectorAscending(Y,Y);

  table=CreateVector(n_table+1);

  B=(int)floor(cut_point*(PRECISION)DimV(Y));
  N=DimV(Y) - B;

  ElementV(table,0)=0.0;
  inc=(PRECISION)N/(PRECISION)n_table;
  for (i=1; i <= n_table; i++)
    {
      j=(int)floor(inc*(PRECISION)i);
      p=(inc*(PRECISION)i - (PRECISION)j)/inc; 
      j+=B;
      ElementV(table,i)=(j < DimV(Y) - 1) ? (1 - p)*sqrt(ElementV(Y,j)) + p*sqrt(ElementV(Y,j+1)) : sqrt(ElementV(Y,DimV(Y)-1));
    }

  FreeVector(Y);

  SetupSpherical_Table(n,pElementV(table),n_table);

  FreeVector(table);

  return 1.0;
}

void ComputeMarginal_Table(int ndraws_proposal, TMatrix X, TStateModel *model, T_MHM *mhm, char *tag)
{
  int idx_alpha=1;
  char filename[256], *proposal_tag="table";
  FILE *f_out;

  SetupSphericalFromPosterior_Table(X,NumberFreeParametersTheta(model));

  // Create output file amd write headers
  sprintf(filename,"%s_md_%s.dat",proposal_tag,tag);
  f_out=dw_CreateTextFile(filename);

  fprintf(f_out,"tag: %s\n",tag);
  fprintf(f_out,"Table proposal\n");

  fprintf(f_out,"number draws posterior (kept): %d\n",RowM(X));
  fprintf(f_out,"thinning factor: %d\n",mhm->n_thin);
  fprintf(f_out,"total posterior draws: %d\n",RowM(X) * mhm->n_thin);
  fprintf(f_out,"number draws proposal: %d\n",ndraws_proposal);

  fclose(f_out);

  ComputeMarginal(ndraws_proposal,X,model,mhm,idx_alpha,1.0,tag,proposal_tag);
}

void PlotCumulative_Table(TMatrix X, int n, char *tag)
{
  PRECISION cut_point=0.01, inc, max;
  int bins=1000, n_table=100, i, j;
  char filename[256];
  FILE *f_out;
  TVector cumulative;

  SetupSphericalFromPosterior_Table(X,n);
  max=sqrt(ElementM(X,RowM(X)-1,1));
  cumulative=SphericalCumulativeDensity(max,bins,20);
  sprintf(filename,"cumulative_table_%s.csv",tag);
  f_out=dw_CreateTextFile(filename);
  fprintf(f_out,"tag: %s\n",tag);
  fprintf(f_out,"number table entries: %d\n",n_table);
  fprintf(f_out,"number draws posterior: %d\n",RowM(X));
  fprintf(f_out,"cumulative densities\n");
  inc=max/(PRECISION)bins;
  for (i=j=0; i < bins; i++)
    {
      while ((j < RowM(X)) && (sqrt(ElementM(X,j,1)) <= (i+1)*inc)) j++;
      fprintf(f_out,"%lf,%lf,%lf\n",(PRECISION)(i+1)*inc,(PRECISION)j/(PRECISION)RowM(X),ElementV(cumulative,i));
    }
  fclose(f_out);
  FreeVector(cumulative);
}
/*******************************************************************************/

/*******************************************************************************/
PRECISION SetupSphericalFromPosterior_TruncatedPower(TMatrix X, int n)
{
  TVector Y;
  PRECISION cut_point=0.01, min_point=0.1, max_point=0.9, a, b, k, truncate, rescale_factor;
  int i;

  Y=ColumnVector((TVector)NULL,X,1);
  SortVectorAscending(Y,Y);

  i=floor(min_point*(PRECISION)DimV(Y));
  a=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));
  i=floor(max_point*(PRECISION)DimV(Y));
  b=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));

  k=log(min_point/max_point)/log(a/b);
  rescale_factor=b/pow(max_point,1.0/k);
  i=floor(cut_point*(PRECISION)DimV(Y));
  truncate=((i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i)))/rescale_factor;

  FreeVector(Y);

  SetupSpherical_TruncatedPower(n,k,truncate);

  return rescale_factor;
}
void ComputeMarginal_TruncatedPowerProposal(int ndraws_proposal, TMatrix X, TStateModel *model, T_MHM *mhm, char *tag)
{
  TMatrix proposal, posterior, base_scale, scale;
  TVector Y, *alpha=(TVector*)NULL, level_cuts;
  PRECISION rescale_factor, p, s=1.0/(PRECISION)RowM(X), cut_point=0.01, min_point=0.1, max_point=0.9, a, b, x, y, inc, k, truncate, 
    diff, marginal_density, ess, slope1, ess1, slope2, ess2;
  int idx_alpha=1, i, j, in1, in2, bins=1000, ia, ib;
  char filename[256];
  FILE *f_out;

  Y=ColumnVector((TVector)NULL,X,1);
  SortVectorAscending(Y,Y);

  i=floor(min_point*(PRECISION)DimV(Y));
  a=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));
  i=floor(max_point*(PRECISION)DimV(Y));
  b=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));

  k=log(min_point/max_point)/log(a/b);
  rescale_factor=b/pow(max_point,1.0/k);
  i=floor(cut_point*(PRECISION)DimV(Y));
  truncate=((i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i)))/rescale_factor;

  FreeVector(Y);
  SetupSpherical_TruncatedPower(NumberFreeParametersTheta(model),k,truncate);

/*   // Plot cummulative densities */
/*   f_out=dw_CreateTextFile("tmp.csv"); */
/*   fprintf(f_out,"tag: %s\n",tag); */
/*   fprintf(f_out,"scale factor: %lf\n",rescale_factor); */
/*   fprintf(f_out,"power: %lf\n",k); */
/*   fprintf(f_out,"truncation: %lf\n",truncate); */
/*   fprintf(f_out,"number draws posterior: %d\n",RowM(Y)); */
/*   inc=1.0/(PRECISION)bins; */
/*   p=pow(truncate,k); */
/*   for (i=j=0; i < bins; i++) */
/*     { */
/*       while ((j < RowM(Y)) && (sqrt(ElementM(Y,j,1)) <= (i+1)*inc*rescale_factor)) j++; */
/*       x=((PRECISION)i+0.5)*inc; */
/*       y=(x < truncate) ? 0.0 : (pow(x,k) - p)/(1.0 - p); */
/*       fprintf(f_out,"%lf,%lf,%lf\n",x,(PRECISION)j/(PRECISION)RowM(Y),y); */
/*     } */
/*   fclose(f_out); */
/*   return; */

  // compute base scale and scale
  base_scale=CholeskyUT((TMatrix)NULL,mhm->variance);
  Transpose(base_scale,base_scale);

  // allocate alpha
  alpha=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    alpha[i]=CreateVector(DimV(mhm->BaseAlpha[i]));

  // Dirichlet index
  idx_alpha=1;
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    ProductVS(alpha[i],mhm->BaseAlpha[i],ElementV(mhm->alpha_scales,idx_alpha));

  // Create posterior matrix
  printf("Creating posterior (%d draws)\n",RowM(X));
  posterior=CreatePosterior(X,rescale_factor,base_scale,idx_alpha,NumberFreeParametersTheta(model));

  // Create proposal matrix
  printf("Creating proposal (%d draws)\n",ndraws_proposal);
  scale=ProductSM((TMatrix)NULL,rescale_factor,base_scale);
  proposal=CreateProposal(ndraws_proposal,mhm->center,scale,model,alpha);

  // Create output file
  sprintf(filename,"truncatedpower_md_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);

  fprintf(f_out,"tag: %s\n",tag);
  fprintf(f_out,"Truncated power proposal\n");
  fprintf(f_out,"scale factor: %lf\n",rescale_factor);
  fprintf(f_out,"power: %lf\n",k);
  fprintf(f_out,"truncation: %lf\n",truncate);

  fprintf(f_out,"number draws posterior (kept): %d\n",RowM(X));
  fprintf(f_out,"thinning factor: %d\n",mhm->n_thin);
  fprintf(f_out,"total posterior draws: %d\n",RowM(X) * mhm->n_thin);
  fprintf(f_out,"number draws proposal: %d\n",ndraws_proposal);

  // Compute marginal density using standard technique
  marginal_density=ComputeMarginalDensity_Standard(posterior);
  fprintf(f_out,"Standard MHM method: %0.14lg\n",marginal_density);

  // Compute marginal density using WZ1
  level_cuts=Sorted_Level_Cuts(posterior,10);
  fprintf(f_out,"WZ method - no center cut\nmarginal density,level,number proposal draws,number posterior draws,effective sample size\n");
  for (i=0; i < DimV(level_cuts); i++)
    {
      marginal_density=ComputeMarginalDensity_WZ1(proposal,posterior,ElementV(level_cuts,i),&in1,&in2);
      ess=ComputeEffectiveSampleSize_WZ1(posterior,ElementV(level_cuts,i));

      fprintf(f_out,"%0.14lg,%lg,%d,%d,%lg\n",marginal_density,ElementV(level_cuts,i),in1,in2,ess);
    }
  FreeVector(level_cuts);

  // Compute marginal density Mueller
  fprintf(f_out,"\nMueller method\nmarginal density,difference,# proposal draws,proposal log slope,proposal ess,# posterior draws,posterior log slope,posterior ess\n");
  marginal_density=ComputeMarginalDensity_Mueller(proposal,posterior,&in1,&in2);
  ComputeDiagnostics_Mueller(proposal,posterior,-marginal_density,&slope1,&ess1,&slope2,&ess2);
  diff=ComputeDifference(proposal,posterior,-marginal_density,&in1,&in2);
  fprintf(f_out,"%0.14lg,%lg,%d,%lg,%lg,%d,%lg,%lg\n\n",marginal_density,diff,in1,slope1,ess1,in2,slope2,ess2);

  fclose(f_out);

  sprintf(filename,"truncatedpower_md_proposal_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,proposal,"%le,");
  fclose(f_out);

  sprintf(filename,"truncatedpower_md_posterior_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,posterior,"%le,");
  fclose(f_out);

  FreeMatrix(posterior);
  FreeMatrix(proposal);

  dw_FreeArray(alpha);
  FreeMatrix(base_scale);
  FreeMatrix(scale);
}

PRECISION SetupSphericalFromPosterior_Power(TMatrix X, int n)
{
  int i;
  PRECISION a, b, k, rescale_factor, min_point=0.1, max_point=0.9;
  TVector Y=ColumnVector((TVector)NULL,X,1);

  SortVectorAscending(Y,Y);

  i=floor(min_point*(PRECISION)DimV(Y));
  a=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));
  i=floor(max_point*(PRECISION)DimV(Y));
  b=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));

  k=log(min_point/max_point)/log(a/b);
  rescale_factor=b/pow(max_point,1.0/k);

  FreeVector(Y);

  SetupSpherical_Power(n,k);

  return rescale_factor;
}
void ComputeMarginal_PowerProposal(int ndraws_proposal, TMatrix X, TStateModel *model, T_MHM *mhm, char *tag)
{
  TMatrix proposal, posterior, base_scale, scale;
  TVector Y, *alpha=(TVector*)NULL, level_cuts;
  PRECISION rescale_factor, p, s=1.0/RowM(X), min_point=0.1, max_point=0.9, a, b, k, 
    diff, marginal_density, ess, slope1, ess1, slope2, ess2, inc;
  int idx_alpha=1, i, j, in1, in2, bins=1000;
  char filename[256];
  FILE *f_out;

  Y=ColumnVector((TVector)NULL,X,1);
  SortVectorAscending(Y,Y);

  i=floor(min_point*(PRECISION)DimV(Y));
  a=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));
  i=floor(max_point*(PRECISION)DimV(Y));
  b=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i));

  k=log(min_point/max_point)/log(a/b);
  rescale_factor=b/pow(max_point,1.0/k);

  FreeVector(Y);

  SetupSpherical_Power(NumberFreeParametersTheta(model),k);

/*   // Plot cummulative densities */
/*   f_out=dw_CreateTextFile("tmp.csv"); */
/*   fprintf(f_out,"tag: %s\n",tag); */
/*   fprintf(f_out,"scale factor: %lf\n",rescale_factor); */
/*   fprintf(f_out,"power: %lf\n",k); */
/*   fprintf(f_out,"number draws posterior: %d\n",RowM(Y)); */
/*   inc=1.0/(PRECISION)bins; */
/*   for (i=j=0; i < bins; i++) */
/*     { */
/*       while ((j < RowM(Y)) && (sqrt(ElementM(Y,j,1)) <= (i+1)*inc*rescale_factor)) j++; */
/*       fprintf(f_out,"%lf,%lf,%lf\n",((PRECISION)i+0.5)*inc,(PRECISION)j/(PRECISION)RowM(Y),pow(((PRECISION)i+0.5)*inc,k)); */
/*     } */
/*   fclose(f_out); */
/*   return; */

  // compute base scale and scale
  base_scale=CholeskyUT((TMatrix)NULL,mhm->variance);
  Transpose(base_scale,base_scale);

  // allocate alpha
  alpha=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    alpha[i]=CreateVector(DimV(mhm->BaseAlpha[i]));

  // Dirichlet index
  idx_alpha=1;
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    ProductVS(alpha[i],mhm->BaseAlpha[i],ElementV(mhm->alpha_scales,idx_alpha));

  // Create posterior matrix
  printf("Creating posterior (%d draws)\n",RowM(X));
  posterior=CreatePosterior(X,rescale_factor,base_scale,idx_alpha,NumberFreeParametersTheta(model));

  // Create proposal matrix
  printf("Creating proposal (%d draws)\n",ndraws_proposal);
  scale=ProductSM((TMatrix)NULL,rescale_factor,base_scale);
  proposal=CreateProposal(ndraws_proposal,mhm->center,scale,model,alpha);

  // Create output file
  sprintf(filename,"power_md_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);

  fprintf(f_out,"tag: %s\n",tag);
  fprintf(f_out,"Power proposal\n");
  fprintf(f_out,"scale factor: %lf\n",rescale_factor);
  fprintf(f_out,"power: %lf\n",k);

  fprintf(f_out,"number draws posterior (kept): %d\n",RowM(X));
  fprintf(f_out,"thinning factor: %d\n",mhm->n_thin);
  fprintf(f_out,"total posterior draws: %d\n",RowM(X) * mhm->n_thin);
  fprintf(f_out,"number draws proposal: %d\n",ndraws_proposal);

  // Compute marginal density using standard technique
  marginal_density=ComputeMarginalDensity_Standard(posterior);
  fprintf(f_out,"Standard MHM method: %0.14lg\n",marginal_density);

  // Compute marginal density using WZ1
  level_cuts=Sorted_Level_Cuts(posterior,10);
  fprintf(f_out,"WZ method - no center cut\nmarginal density,level,number proposal draws,number posterior draws,effective sample size\n");
  for (i=0; i < DimV(level_cuts); i++)
    {
      marginal_density=ComputeMarginalDensity_WZ1(proposal,posterior,ElementV(level_cuts,i),&in1,&in2);
      ess=ComputeEffectiveSampleSize_WZ1(posterior,ElementV(level_cuts,i));

      fprintf(f_out,"%0.14lg,%lg,%d,%d,%lg\n",marginal_density,ElementV(level_cuts,i),in1,in2,ess);
    }
  FreeVector(level_cuts);

  // Compute marginal density Mueller
  fprintf(f_out,"\nMueller method\nmarginal density,difference,# proposal draws,proposal log slope,proposal ess,# posterior draws,posterior log slope,posterior ess\n");
  marginal_density=ComputeMarginalDensity_Mueller(proposal,posterior,&in1,&in2);
  ComputeDiagnostics_Mueller(proposal,posterior,-marginal_density,&slope1,&ess1,&slope2,&ess2);
  diff=ComputeDifference(proposal,posterior,-marginal_density,&in1,&in2);
  fprintf(f_out,"%0.14lg,%lg,%d,%lg,%lg,%d,%lg,%lg\n\n",marginal_density,diff,in1,slope1,ess1,in2,slope2,ess2);

  fclose(f_out);

  sprintf(filename,"power_md_proposal_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,proposal,"%le,");
  fclose(f_out);

  sprintf(filename,"power_md_posterior_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,posterior,"%le,");
  fclose(f_out);

  FreeMatrix(posterior);
  FreeMatrix(proposal);

  dw_FreeArray(alpha);
  FreeMatrix(base_scale);
  FreeMatrix(scale);
}

void ComputeMarginal_GaussianProposal(int ndraws_proposal, TMatrix X, TStateModel *model, T_MHM *mhm, char *tag, PRECISION variance)
{
  TMatrix proposal, posterior, base_scale, scale;
  TVector *alpha=(TVector*)NULL, level_cuts;
  PRECISION rescale_factor, diff, marginal_density, ess, slope1, ess1, slope2, ess2;
  int idx_alpha=1, i, in1, in2;
  char filename[256];
  FILE *f_out;

  // set spherical type
  SetupSpherical_Gaussian(NumberFreeParametersTheta(model));

  // compute base scale
  base_scale=CholeskyUT((TMatrix)NULL,mhm->variance);
  Transpose(base_scale,base_scale);

  // allocate alpha
  alpha=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    alpha[i]=CreateVector(DimV(mhm->BaseAlpha[i]));

  // Dirichlet index
  idx_alpha=1;
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    ProductVS(alpha[i],mhm->BaseAlpha[i],ElementV(mhm->alpha_scales,idx_alpha));

  // rescale factor
  rescale_factor=sqrt(variance);

  // Create posterior matrix
  printf("Creating posterior (%d draws)\n",RowM(X));
  posterior=CreatePosterior(X,rescale_factor,base_scale,idx_alpha,NumberFreeParametersTheta(model));

  // Create proposal matrix
  printf("Creating proposal (%d draws)\n",ndraws_proposal);
  scale=ProductSM((TMatrix)NULL,rescale_factor,base_scale);
  proposal=CreateProposal(ndraws_proposal,mhm->center,scale,model,alpha);

  // Create output file
  sprintf(filename,"gaussian_md_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);

  fprintf(f_out,"tag: %s\n",tag);
  fprintf(f_out,"Gaussian proposal\n");
  fprintf(f_out,"variance: %lf\n",variance);

  fprintf(f_out,"number draws posterior (kept): %d\n",RowM(X));
  fprintf(f_out,"thinning factor: %d\n",mhm->n_thin);
  fprintf(f_out,"total posterior draws: %d\n",RowM(X) * mhm->n_thin);
  fprintf(f_out,"number draws proposal: %d\n",ndraws_proposal);

  // Compute marginal density using standard technique
  marginal_density=ComputeMarginalDensity_Standard(posterior);
  fprintf(f_out,"Standard MHM method: %0.14lg\n",marginal_density);

  // Compute marginal density using WZ1
  level_cuts=Sorted_Level_Cuts(posterior,10);
  fprintf(f_out,"WZ method - no center cut\nmarginal density,level,number proposal draws,number posterior draws,effective sample size\n");
  for (i=0; i < DimV(level_cuts); i++)
    {
      marginal_density=ComputeMarginalDensity_WZ1(proposal,posterior,ElementV(level_cuts,i),&in1,&in2);
      ess=ComputeEffectiveSampleSize_WZ1(posterior,ElementV(level_cuts,i));

      fprintf(f_out,"%0.14lg,%lg,%d,%d,%lg\n",marginal_density,ElementV(level_cuts,i),in1,in2,ess);
    }
  FreeVector(level_cuts);

  // Compute marginal density Mueller
  fprintf(f_out,"\nMueller method\nmarginal density,difference,# proposal draws,proposal log slope,proposal ess,# posterior draws,posterior log slope,posterior ess\n");
  marginal_density=ComputeMarginalDensity_Mueller(proposal,posterior,&in1,&in2);
  ComputeDiagnostics_Mueller(proposal,posterior,-marginal_density,&slope1,&ess1,&slope2,&ess2);
  diff=ComputeDifference(proposal,posterior,-marginal_density,&in1,&in2);
  fprintf(f_out,"%0.14lg,%lg,%d,%lg,%lg,%d,%lg,%lg\n\n",marginal_density,diff,in1,slope1,ess1,in2,slope2,ess2);

  fclose(f_out);

  sprintf(filename,"gaussian_md_proposal_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,proposal,"%le,");
  fclose(f_out);

  sprintf(filename,"gaussian_md_posterior_%s.dat",tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,posterior,"%le,");
  fclose(f_out);

  FreeMatrix(posterior);
  FreeMatrix(proposal);

  dw_FreeArray(alpha);
  FreeMatrix(base_scale);
  FreeMatrix(scale);
}

void ComputeMarginal_TruncatedGaussianProposal(int ndraws_proposal, TMatrix X, TStateModel *model, T_MHM *mhm, char *tag, PRECISION p1, PRECISION p2)
{
  TMatrix proposal, posterior, base_scale, scale;
  TVector Y, *alpha=(TVector*)NULL, level_cuts;
  PRECISION rescale_factor, diff, marginal_density, ess, slope1, ess1, slope2, ess2, min_point=0.05, max_point=0.95, r1, r2;
  int idx_alpha=1, i, in1, in2;
  char filename[256];
  FILE *f_out;

  // set spherical type

/*   Y=ColumnVector((TVector)NULL,X,1); */
/*   SortVectorAscending(Y,Y); */

/*   i=floor(min_point*(PRECISION)DimV(Y)); */
/*   r1=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i)); */
/*   i=floor(max_point*(PRECISION)DimV(Y)); */
/*   r2=(i > 0) ? 0.5*(sqrt(ElementV(Y,i-1)) + sqrt(ElementV(Y,i))) : sqrt(ElementV(Y,i)); */

/*   FreeVector(Y); */

  r1=sqrt(dw_chi_square_invcdf(p1,NumberFreeParametersTheta(model)));
  r2=sqrt(dw_chi_square_invcdf(1.0 - p2,NumberFreeParametersTheta(model)));

  // rescale factor
  rescale_factor=1;

  SetupSpherical_TruncatedGaussian(NumberFreeParametersTheta(model),r1,r2);

/*   // Plot cummulative densities */
/*   f_out=dw_CreateTextFile("tmp.csv"); */
/*   fprintf(f_out,"tag: %s\n",tag); */
/*   fprintf(f_out,"scale factor: %lf\n",rescale_factor); */
/*   fprintf(f_out,"power: %lf\n",k); */
/*   fprintf(f_out,"truncation: %lf\n",truncate); */
/*   fprintf(f_out,"number draws posterior: %d\n",RowM(Y)); */
/*   inc=1.0/(PRECISION)bins; */
/*   p=pow(truncate,k); */
/*   for (i=j=0; i < bins; i++) */
/*     { */
/*       while ((j < RowM(Y)) && (sqrt(ElementM(Y,j,1)) <= (i+1)*inc*rescale_factor)) j++; */
/*       x=((PRECISION)i+0.5)*inc; */
/*       y=(x < truncate) ? 0.0 : (pow(x,k) - p)/(1.0 - p); */
/*       fprintf(f_out,"%lf,%lf,%lf\n",x,(PRECISION)j/(PRECISION)RowM(Y),y); */
/*     } */
/*   fclose(f_out); */
/*   return; */

  // compute base scale
  base_scale=CholeskyUT((TMatrix)NULL,mhm->variance);
  Transpose(base_scale,base_scale);

  // allocate alpha
  alpha=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    alpha[i]=CreateVector(DimV(mhm->BaseAlpha[i]));

  // Dirichlet index
  idx_alpha=1;
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    ProductVS(alpha[i],mhm->BaseAlpha[i],ElementV(mhm->alpha_scales,idx_alpha));

  // Create posterior matrix
  printf("Creating posterior (%d draws)\n",RowM(X));
  posterior=CreatePosterior(X,rescale_factor,base_scale,idx_alpha,NumberFreeParametersTheta(model));

  // Create proposal matrix
  printf("Creating proposal (%d draws)\n",ndraws_proposal);
  scale=ProductSM((TMatrix)NULL,rescale_factor,base_scale);
  proposal=CreateProposal(ndraws_proposal,mhm->center,scale,model,alpha);

  // Create output file
  sprintf(filename,"truncatedgaussian_%.2lf_%.2lf_md_%s.dat",p1,p2,tag);
  f_out=dw_CreateTextFile(filename);

  fprintf(f_out,"tag: %s\n",tag);
  fprintf(f_out,"Truncated Gaussian proposal\n");
  fprintf(f_out,"center cut: %lf  (radius: %lf)\n",p1,r1);
  fprintf(f_out,"tail cut: %lf   (radius: %lf)\n",p2,r2);

  fprintf(f_out,"number draws posterior (kept): %d\n",RowM(X));
  fprintf(f_out,"thinning factor: %d\n",mhm->n_thin);
  fprintf(f_out,"total posterior draws: %d\n",RowM(X) * mhm->n_thin);
  fprintf(f_out,"number draws proposal: %d\n",ndraws_proposal);

  // Compute marginal density using standard technique
  marginal_density=ComputeMarginalDensity_Standard(posterior);
  fprintf(f_out,"Standard MHM method: %0.14lg\n",marginal_density);

  // Compute marginal density using WZ1
  level_cuts=Sorted_Level_Cuts(posterior,10);
  fprintf(f_out,"WZ method - no center cut\nmarginal density,level,number proposal draws,number posterior draws,effective sample size\n");
  for (i=0; i < DimV(level_cuts); i++)
    {
      marginal_density=ComputeMarginalDensity_WZ1(proposal,posterior,ElementV(level_cuts,i),&in1,&in2);
      ess=ComputeEffectiveSampleSize_WZ1(posterior,ElementV(level_cuts,i));

      fprintf(f_out,"%0.14lg,%lg,%d,%d,%lg\n",marginal_density,ElementV(level_cuts,i),in1,in2,ess);
    }
  FreeVector(level_cuts);

  // Compute marginal density Mueller
  fprintf(f_out,"\nMueller method\nmarginal density,difference,# proposal draws,proposal log slope,proposal ess,# posterior draws,posterior log slope,posterior ess\n");
  marginal_density=ComputeMarginalDensity_Mueller(proposal,posterior,&in1,&in2);
  ComputeDiagnostics_Mueller(proposal,posterior,-marginal_density,&slope1,&ess1,&slope2,&ess2);
  diff=ComputeDifference(proposal,posterior,-marginal_density,&in1,&in2);
  fprintf(f_out,"%0.14lg,%lg,%d,%lg,%lg,%d,%lg,%lg\n\n",marginal_density,diff,in1,slope1,ess1,in2,slope2,ess2);

  fclose(f_out);

  sprintf(filename,"truncatedgaussian_%.2lf_%.2lf_md_proposal_%s.dat",p1,p2,tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,proposal,"%le,");
  fclose(f_out);

  sprintf(filename,"truncatedgaussian_%.2lf_%.2lf_md_posterior_%s.dat",p1,p2,tag);
  f_out=dw_CreateTextFile(filename);
  dw_PrintMatrix(f_out,posterior,"%le,");
  fclose(f_out);

  FreeMatrix(posterior);
  FreeMatrix(proposal);

  dw_FreeArray(alpha);
  FreeMatrix(base_scale);
  FreeMatrix(scale);
}


/* PRECISION SetupSphericalFromPosterior(TMatrix X, int n, int type) */
/* { */
/*   switch (type) */
/*     { */
/*       //case TYPE_GAUSSIAN: */
/*     case TYPE_POWER: */
/*       return SetupSphericalFromPosterior_Power(X,n); */
/*     case TYPE_TRUNCATED_POWER: */
/*       return SetupSphericalFromPosterior_TruncatedPower(X,n); */
/*     case TYPE_TABLE: */
/*       return SetupSphericalFromPosterior_Table(X,n); */
/*     default: */
/*       swz_fprintf_err("Unknown proposal type\n"); */
/*       exit(0); */
/*     } */
/* } */

/* FILE* CreateOutputFile(int type, char *tag, PRECISION rescale_factor) */
/* { */
/*   FILE *f_out; */
/*   char filename[256]; */

/*   // Create output file */
/*   switch (type) */
/*     { */
/*     case TYPE_GAUSSIAN: */
/*       sprintf(filename,"gaussian_md_%s.dat",tag); */
/*       break; */
/*     case TYPE_POWER: */
/*       sprintf(filename,"power_md_%s.dat",tag); */
/*       break; */
/*     case TYPE_TRUNCATED_POWER: */
/*       sprintf(filename,"truncatedpower_md_%s.dat",tag); */
/*       break; */
/*     case TYPE_TABLE: */
/*       sprintf(filename,"table_md_%s.dat",tag); */
/*       break; */
/*     default: */
/*       swz_fprintf_err("Unknown proposal type\n"); */
/*       exit(0); */
/*     } */
/*   f_out=dw_CreateTextFile(filename); */

/*   // Write header information */
/*   fprintf(f_out,"tag: %s\n",tag); */
/*   fprintf(f_out,"%s\n",SphericalType()); */

/*   return f_out; */
/* } */

/*
   Computes marginal data density using the WZ technique.  A small footprint technique is used to minimize 
   memory and file usage.
*/
/* void ComputeMarginal_Small(int type, int ndraws_proposal, TMatrix X, TStateModel *model, T_MHM *mhm, char *tag, PRECISION variance) */
/* { */
/*   TMatrix  posterior, base_scale, scale; */
/*   TVector *alpha=(TVector*)NULL, level_cuts, q; */
/*   PRECISION rescale_factor, marginal_density, ess; */
/*   int idx_alpha=1, i, in2; */
/*   FILE *f_out; */

/*   // Setup spherical distribution */
/*   rescale_factor=SetupSphericalFromPosterior(X,NumberFreeParametersTheta(model),type); */
 
/*   // compute base scale and scale */
/*   base_scale=CholeskyUT((TMatrix)NULL,mhm->variance); */
/*   Transpose(base_scale,base_scale); */
/*   scale=ProductSM((TMatrix)NULL,rescale_factor,base_scale); */

/*   // allocate alpha */
/*   alpha=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha)); */
/*   for (i=dw_DimA(alpha)-1; i >= 0; i--) */
/*     alpha[i]=CreateVector(DimV(mhm->BaseAlpha[i])); */

/*   // Dirichlet index */
/*   idx_alpha=1; */
/*   for (i=dw_DimA(alpha)-1; i >= 0; i--) */
/*     ProductVS(alpha[i],mhm->BaseAlpha[i],ElementV(mhm->alpha_scales,idx_alpha)); */

/*   // Create posterior matrix */
/*   printf("Creating posterior (%d draws)\n",RowM(X)); */
/*   posterior=CreatePosterior(X,rescale_factor,base_scale,idx_alpha,NumberFreeParametersTheta(model)); */

/*   // Create output file and write header */
/*   f_out=CreateOutputFile(type,tag,rescale_factor); */
/*   fprintf(f_out,"number draws posterior: %d (trimming factor: %d)\n",RowM(X),mhm->n_thin); */
/*   fprintf(f_out,"number draws proposal: %d\n",ndraws_proposal); */

/*   // Compute marginal density using WZ1 */
/*   level_cuts=Sorted_Level_Cuts(posterior,10); */
/*   printf("Creating proposal (%d draws)\n",ndraws_proposal); */
/*   q=Create_q(ndraws_proposal,mhm->center,scale,model,alpha,level_cuts); */
/*   fprintf(f_out,"WZ method\nmarginal density,level,number proposal draws,number posterior draws,effective sample size\n"); */
/*   for (i=0; i < DimV(level_cuts); i++) */
/*     { */
/*       marginal_density=ComputeMarginalDensity_WZ1_q(posterior,ElementV(level_cuts,i),ElementV(q,i),&in2); */
/*       ess=ComputeEffectiveSampleSize_WZ1(posterior,ElementV(level_cuts,i)); */

/*       fprintf(f_out,"%0.14lg,%lg,%d,%d,%lg\n",marginal_density,ElementV(level_cuts,i),(int)floor(ElementV(q,i)*ndraws_proposal),in2,ess); */
/*     } */

/*   // Clean up */
/*   FreeVector(level_cuts); */
/*   FreeVector(q); */

/*   fclose(f_out); */

/*   FreeMatrix(posterior); */

/*   dw_FreeArray(alpha); */
/*   FreeMatrix(base_scale); */
/*   FreeMatrix(scale); */
/* } */

/*
   Creates a cumulative distribution of the radius from the posterior 
   distribution.  The square of the radius is the second column of X.
   The larger the bin size, the finer the plot.  The first column is the
   radius and the second column is the cumulative distribution.
*/
void CreateCumulativeDistributionPosteriorRadius(TMatrix X, char *filename, int bins)
{
  FILE *fout=fopen(filename,"wt");
  TVector r=CreateVector(RowM(X));
  TMatrix cumulative;
  PRECISION inc, max;
  int i, j;

  for (i=DimV(r)-1; i >= 0; i--) ElementV(r,i)=sqrt(ElementM(X,i,1));
  SortVectorAscending(r,r);
  max=ElementV(r,DimV(r)-1);
  inc=max/bins;

  cumulative=CreateMatrix(bins,2);
  for (i=j=0; i < bins; i++)
    {
      while ((j < DimV(r)) && (ElementV(r,j) <= (i+1)*inc)) j++;
      ElementM(cumulative,i,0)=(i+0.5)*inc;
      ElementM(cumulative,i,1)=(PRECISION)j/(PRECISION)DimV(r);
    }

  dw_PrintMatrix(fout,cumulative,"%lf,");
  fclose(fout);
  FreeMatrix(cumulative);
  FreeVector(r);
}

/*
   Creates a cumulative distribution of the radius from the posterior 
   distribution.  The square of the radius is the second column of X.
   The larger the bin size, the finer the plot.  The first column is the
   radius and the second column is the cumulative distribution.
*/
void CreateCumulativeDistributionSphericalRadius(int ndraws, int dim, char *filename, int bins)
{
  FILE *fout=fopen(filename,"wt");
  TMatrix cumulative;
  TVector x=CreateVector(dim);
  PRECISION inc, max=1, r, s=1.0/(PRECISION)ndraws;
  int i, j;

  inc=max/bins;
  cumulative=CreateMatrix(bins,2);
  for (i=0; i < bins; i++)
    {
      ElementM(cumulative,i,0)=(i+0.5)*inc;
      ElementM(cumulative,i,1)=0.0;
    }

  for (i=ndraws; i > 0; i--)
    {
      r=DrawSpherical(x);
      j=(int)floor(bins*r/max);
      if (j >= bins) j=bins-1;
      ElementM(cumulative,j,1)+=1.0;
    }
  ElementM(cumulative,0,1)*=s;
  for (i=1; i < bins; i++) 
    ElementM(cumulative,i,1)=s*ElementM(cumulative,i,1) + ElementM(cumulative,i-1,1);

  dw_PrintMatrix(fout,cumulative,"%lf,");
  fclose(fout);
  FreeMatrix(cumulative);
  FreeVector(x);
}

/*
   Outputs  
*/
void Plot_Posterior_vs_Posterior_Radius(TMatrix X, char* filename)
{
  FILE *f_out=dw_CreateTextFile(filename);
  int bins=1000, inc=(RowM(X)-1)/(bins-1), start=(RowM(X) - inc*(bins-1))/2, stop, i, j;
  TMatrix Y;

  SortMatrixRowsAscending(X,X,1);
  Y=CreateMatrix(bins,4);
  for (i=0; i < bins; i++)
    {
      stop=start+inc;
      if (stop > RowM(X)) stop=RowM(X);
      ElementM(Y,i,0)=(ElementM(X,start,1) + ElementM(X,stop-1,1))/2;
      ElementM(Y,i,1)=ElementM(Y,i,2)=ElementM(Y,i,3)=ElementM(X,start,0);
      for (j=start+1; j < stop; j++)
	{
	  ElementM(Y,i,2)+=ElementM(X,j,0);
	  if (ElementM(X,j,0) < ElementM(Y,i,1))
	    ElementM(Y,i,1)=ElementM(X,j,0);
	  else
	    if (ElementM(X,j,0) > ElementM(Y,i,3))
	      ElementM(Y,i,3)=ElementM(X,j,0);
	}
      ElementM(Y,i,2)/=(PRECISION)(stop - start);
      start+=inc;
    }

  dw_PrintMatrix(f_out,Y,"%lf,");
  FreeMatrix(Y);
  fclose(f_out);
}

void Plot_Posterior_vs_Proposal_Radius(int ndraws, TStateModel *model, T_MHM *mhm, char* filename, PRECISION rescale_factor)
{
  FILE *f_out=dw_CreateTextFile(filename);
  int bins=1000, inc=(ndraws-1)/(bins-1), start=(ndraws - inc*(bins-1))/2, idx_alpha=1, stop, i, j;
  TMatrix Y, X, base_scale, scale;
  TVector *alpha;
  
  // compute base scale and scale
  base_scale=CholeskyUT((TMatrix)NULL,mhm->variance);
  Transpose(base_scale,base_scale);

  // allocate scale and alpha
  alpha=dw_CreateArray_vector(dw_DimA(mhm->BaseAlpha));
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    alpha[i]=CreateVector(DimV(mhm->BaseAlpha[i]));

  // Dirichlet index
  idx_alpha=1;
  for (i=dw_DimA(alpha)-1; i >= 0; i--)
    ProductVS(alpha[i],mhm->BaseAlpha[i],ElementV(mhm->alpha_scales,idx_alpha));

  // Create proposal matrix
  printf("Creating proposal (%d draws)\n",ndraws);
  scale=ProductSM((TMatrix)NULL,rescale_factor,base_scale);
  X=CreateProposal_Radius(ndraws,mhm->center,scale,model,alpha);
  FreeMatrix(scale);
  FreeMatrix(base_scale);
  dw_FreeArray(alpha);

  SortMatrixRowsAscending(X,X,0);
  Y=CreateMatrix(bins,4);
  for (i=0; i < bins; i++)
    {
      stop=start+inc;
      if (stop > RowM(X)) stop=RowM(X);
      ElementM(Y,i,0)=rescale_factor*(ElementM(X,start,0) + ElementM(X,stop-1,0))/2;
      ElementM(Y,i,1)=ElementM(Y,i,2)=ElementM(Y,i,3)=ElementM(X,start,1);
      for (j=start+1; j < stop; j++)
	{
	  ElementM(Y,i,2)+=ElementM(X,j,1);
	  if (ElementM(X,j,1) < ElementM(Y,i,1))
	    ElementM(Y,i,1)=ElementM(X,j,1);
	  else
	    if (ElementM(X,j,1) > ElementM(Y,i,3))
	      ElementM(Y,i,3)=ElementM(X,j,1);
	}
      ElementM(Y,i,2)/=(PRECISION)(stop - start);
      start+=inc;
    }
  FreeMatrix(X);

  dw_PrintMatrix(f_out,Y,"%lf,");
  FreeMatrix(Y);
  fclose(f_out);
}

/*
   Overview:
     filename contains the following:

       (A) Information required to create and setup T_MHM structure.
       (B) Data - the columns are

           (1) - log posterior value for each draw
           (2) - fp_theta' * Inverse(variance) * fp_theta for each draw
           (3) - log Dirichlet value for each draw under several different values
                 of the parameters for the Dirichlet parameter

           The parameters of the model are Q and theta.  Q is the matrix of
           transition probabilities and theta contains the model specific
           parameters.  fp_theta denotes the free parametes in theta.

     spec_filename contains information to create and setup the TStateModel
     structure and parameter values for the posterior mode.

   Output:
     for each draw in filename computes the following

                       h(theta,Q)/posterior(theta,Q)

     where h(theta,Q) = f(theta)*g(Q) where f is a spherical distribution on
     x = sqrt(Inverse(variance))*fp_theta and g is Dirichlet distribution on Q.
     Recall that a spherical distribution is one in which the density of x depends
     only on the norm of x and sqrt(X) = Y only if Y'*Y = X.
*/
int main(int nargs, char **args)
{
  char *id, *fmt, *tag, *mhm_filename, *spec_filename;
  FILE *f_in, *f_out=stdout;
  int n_fields, ndraws_proposal, ndraws_posterior, proposal_type;
  TMatrix X;
  TStateModel *model;
  T_MHM *mhm;
  PRECISION p1, p2;

  /*** Test Spherical distribution ***
  int n=120;

  PRECISION r1, r2;
  p1=0.0;
  p2=0.7;
  r1=sqrt(dw_chi_square_invcdf(p1,n));
  r2=sqrt(dw_chi_square_invcdf(1.0 - p2,n));
  printf("p1: %lf   r1: %lf  p1 (computed): %lf\n",p1,r1,dw_chi_square_cdf(r1*r1,n));
  printf("p2: %lf   r2: %lf  p2 (computed): %lf\n",p2,r2,1.0-dw_chi_square_cdf(r2*r2,n));
  SetupSpherical_TruncatedGaussian(n,r1,r2);

  TestSpherical((FILE*)NULL,"tmp.csv",10.0);
  return 0;
  /************************************/

  /*****************************************************************************/
  /* Read command line and input files                                         */
  /*****************************************************************************/
  // help
  if (dw_FindArgument(nargs,args,'h') >= 0)
    {
      printf("\nSyntax: marginal_VAR -d <number posterior draws>"
             "\n                     -ft <file tag name>"
             "\n                     -t <proposal type>"
             "\n                         1: gaussian"
             "\n                         2: power"
             "\n                         3: truncated power"
             "\n                         4: table"
             "\n                         5: truncated gaussian"
             "\n                     -u <upper tail truncation>"
             "\n                     -l <center region truncation>"
             "\n\n");
      exit(0);
    }

  // setup filenames
  if (!(tag=dw_ParseString_String(nargs,args,"ft",(char*)NULL)))
    {
      printf("Tag not specified.  Usage: -ft <file tag name>\n\n");
      exit(0);
    }

  fmt="mhm_draws_%s.dat";
  sprintf(mhm_filename=(char*)malloc(strlen(tag)+strlen(fmt)-1),fmt,tag);
  fmt="est_final_%s.dat";
  sprintf(spec_filename=(char*)malloc(strlen(tag)+strlen(fmt)-1),fmt,tag);

  // get number of proposal draws - default 100000
  ndraws_proposal=dw_ParseInteger_String(nargs,args,"d",100000);

  // Read and create TStateModel strucure and setup normalization
  printf("Creating TStateModel\n");
  model=Read_VAR_Specification((FILE*)NULL,spec_filename);
  ReadTransitionMatrices((FILE*)NULL,spec_filename,"Posterior mode: ",model);
  Read_VAR_Parameters((FILE*)NULL,spec_filename,"Posterior mode: ",model);
  Setup_WZ_Normalization(model->theta,((T_VAR_Parameters*)(model->theta))->A0);

  // Read and create T_MHM structure
  printf("Creating T_MHM structure\n");
  f_in=dw_OpenTextFile(mhm_filename);
  mhm=ReadMHM_Input(f_in,(char*)NULL,(T_MHM*)NULL);
  AddStateModel(model,mhm);
  ReadMeanVariance(f_in,mhm);

  // Read draws
  printf("Reading draws\n");
  ndraws_posterior=mhm->n_mhm;
  n_fields=DimV(mhm->alpha_scales)+3;
  id="//== Draws ==//";
  if (!dw_SetFilePosition(f_in,id) || !dw_ReadMatrix(f_in,X=CreateMatrix(ndraws_posterior,n_fields)))
    {
      if (!dw_SetFilePosition(f_in,id))
	printf("Error in parsing header\n  %s\n",id);
      else
	{
	  char **line;
	  int i;
	  printf("Error in reading data matrix - checking data...\n");
	  rewind(f_in);
	  dw_SetFilePosition(f_in,id);
	  for (i=0; i < ndraws_posterior; i++)
	    {
	      line=dw_ReadDelimitedLine(f_in,' ',REMOVE_EMPTY_FIELDS | STRIP_WHITESPACE);
	      if (!line)
		{
		  printf("Not enough lines - %d\n",i+1);
		  exit(0);
		}
	      if (dw_DimA(line) != n_fields)
		{
		  printf("Error on line %d - incorrect number of fields\n  ",i+1);
		  dw_PrintDelimitedArray(stdout,line,' ');
		  fgetc(stdin);
		}
	      dw_FreeArray(line);
	    }
	}
      exit(0);
    }
  fclose(f_in);

  // set proposal type - default power
  proposal_type=dw_ParseInteger_String(nargs,args,"t",2);

  /*****************************************************************************/
  /*****************************************************************************/

  // Initial generator from clock
  dw_initialize_generator(0);
                    
  // Compute marginal using power functions
  switch (proposal_type)
    {
    case 1:
      ComputeMarginal_GaussianProposal(ndraws_proposal,X,model,mhm,tag,1.0);
      break;
    case 2:
      ComputeMarginal_PowerProposal(ndraws_proposal,X,model,mhm,tag);
      break;
    case 3:
      ComputeMarginal_TruncatedPowerProposal(ndraws_proposal,X,model,mhm,tag);
      break;
    case 4:
      ComputeMarginal_Table(ndraws_proposal,X,model,mhm,tag);
      break;
    case 5:
      p1=dw_ParseFloating_String(nargs,args,"l",0.0);
      p2=dw_ParseFloating_String(nargs,args,"u",0.0);
      ComputeMarginal_TruncatedGaussianProposal(ndraws_proposal,X,model,mhm,tag,p1,p2);
    default:
      printf("\nUnknown proposal density."
             "\n  Usage: -t <proposal type>"
             "\n             1: gaussian"
             "\n             2: power"
             "\n             3: truncated power"
             "\n             4: table"
             "\n             5: truncated gaussian"
             "\n\n");
      break;
    }
  
  FreeMatrix(X);
  FreeStateModel(model);
  FreeMHM(mhm);
}

