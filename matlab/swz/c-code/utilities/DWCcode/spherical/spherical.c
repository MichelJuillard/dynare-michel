
#include "spherical.h"
#include "dw_rand.h"
#include "dw_matrix_rand.h"
#include "dw_error.h"
#include "dw_ascii.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "modify_for_mex.h"

#define SPHERICAL_GAUSSIAN           1
#define SPHERICAL_UNIFORM            2
#define SPHERICAL_POWER              3
#define SPHERICAL_TRUNCATED_POWER    4
#define SPHERICAL_TABLE              5
#define SPHERICAL_TRUNCATED_GAUSSIAN 6

#define PI 3.141592653589793

static int SPHERICAL_TYPE=0;
static int SPHERICAL_DIM=0;
static PRECISION SPHERICAL_CONSTANT=0.0;

static PRECISION SPHERICAL_POWER_EXP=0.0;
static PRECISION SPHERICAL_LOWER_TRUNCATE=0.0;
static PRECISION SPHERICAL_UPPER_TRUNCATE=0.0;

static PRECISION *SPHERICAL_TABLE_VALUES=(PRECISION*)NULL;
static int SPHERICAL_TABLE_LENGTH=0;

/*
   Returns ln(exp(a) + exp(b)) computed to avoid overflow.  If
   a = ln(c) and b = ln(d), as is usually the case, then the
   routine returns ln(c + d).

*/
static PRECISION AddLogs_static(PRECISION a, PRECISION b)
{
  return (a > b) ? a + log(1.0 + exp(b-a)) : b + log(exp(a-b) + 1.0);
}

char* SphericalType(void)
{
  static char buffer[128];
  switch (SPHERICAL_TYPE)
    {
    case SPHERICAL_GAUSSIAN:
      return "Gaussian";
    case SPHERICAL_UNIFORM:
      return "Uniform";
    case SPHERICAL_POWER:
      sprintf(buffer,"Power(%lg)",SPHERICAL_POWER_EXP);
      return buffer;
    case SPHERICAL_TRUNCATED_POWER:
      sprintf(buffer,"TruncatedPower(%lg,%lg)",SPHERICAL_POWER_EXP,SPHERICAL_LOWER_TRUNCATE);
      return buffer;
    case SPHERICAL_TABLE:
      sprintf(buffer,"Table(%d)",SPHERICAL_TABLE_LENGTH);
      return buffer;
    case SPHERICAL_TRUNCATED_GAUSSIAN:
      sprintf(buffer,"TruncatedGaussian(%lg,%lg)",SPHERICAL_LOWER_TRUNCATE,SPHERICAL_UPPER_TRUNCATE);
      return buffer;
    default:
      return "Spherical type not set";
    }
}

void SetupSpherical_Gaussian(int n)
{
  SPHERICAL_TYPE=SPHERICAL_GAUSSIAN;
  SPHERICAL_DIM=n;
  SPHERICAL_CONSTANT=-0.5*n*log(2.0*PI);
}

void SetupSpherical_TruncatedGaussian(int n, PRECISION r1, PRECISION r2)
{
  SPHERICAL_TYPE=SPHERICAL_TRUNCATED_GAUSSIAN;
  SPHERICAL_DIM=n;
  SPHERICAL_CONSTANT=-0.5*n*log(2.0*PI) - log(dw_chi_square_cdf(r2*r2,n) - dw_chi_square_cdf(r1*r1,n));
  SPHERICAL_LOWER_TRUNCATE=r1;
  SPHERICAL_UPPER_TRUNCATE=r2;
}

void SetupSpherical_Uniform(int n)
{
  SPHERICAL_TYPE=SPHERICAL_UNIFORM;
  SPHERICAL_DIM=n;
  SPHERICAL_CONSTANT=log(0.5*n) + dw_log_gamma(0.5*n) - 0.5*n*log(PI);
}

/*
   See the function PowerUnitBall() below for the description of the
   distribution.
*/
void SetupSpherical_Power(int n, PRECISION k)
{
  SPHERICAL_TYPE=SPHERICAL_POWER;
  SPHERICAL_DIM=n;
  SPHERICAL_CONSTANT=log(0.5*k) + dw_log_gamma(0.5*n) - 0.5*n*log(PI);
  SPHERICAL_POWER_EXP=k;
}

void SetupSpherical_TruncatedPower(int n, PRECISION k, PRECISION a)
{
  SPHERICAL_TYPE=SPHERICAL_TRUNCATED_POWER;
  SPHERICAL_DIM=n;
  SPHERICAL_CONSTANT=log(0.5*k/(1.0 - pow(a,k))) + dw_log_gamma(0.5*n) - 0.5*n*log(PI);
  SPHERICAL_POWER_EXP=k;
  SPHERICAL_LOWER_TRUNCATE=a;
}

void SetupSpherical_Table(int n, PRECISION *table, int m)
{
  int i;
  SPHERICAL_TYPE=SPHERICAL_TABLE;
  SPHERICAL_DIM=n;
  SPHERICAL_CONSTANT=log(0.5) + dw_log_gamma(0.5*n) - 0.5*n*log(PI);
  if (SPHERICAL_TABLE_VALUES) swzFree(SPHERICAL_TABLE_VALUES);
  SPHERICAL_TABLE_VALUES=(PRECISION*)swzMalloc((m+1)*sizeof(PRECISION));
  SPHERICAL_TABLE_LENGTH=m;
  memcpy(SPHERICAL_TABLE_VALUES,table,(m+1)*sizeof(PRECISION));

/*    // Check   ansi-c*/
  if (SPHERICAL_TABLE_VALUES[0] != 0.0)
    {
      printf("First entry of inverse cumulative spherical table must be zero\n");
      swzExit(0);
    }
  for (i=1; i < SPHERICAL_TABLE_LENGTH; i++)
    if (SPHERICAL_TABLE_VALUES[i-1] >= SPHERICAL_TABLE_VALUES[i])
      {
    printf("Inverse cumulative spherical table must be strictly increasing\n");
    for (i=0; i <= m; i++) printf("%lf\n",table[i]);
    swzExit(0);
      }
}

PRECISION DrawSpherical(TVector x)
{
  PRECISION r;
  switch (SPHERICAL_TYPE)
    {
    case SPHERICAL_GAUSSIAN:
      dw_NormalVector(x);
      return Norm(x);
    case SPHERICAL_UNIFORM:
      return UniformUnitBall(x);
    case SPHERICAL_POWER:
      return PowerUnitBall(x,SPHERICAL_POWER_EXP);
    case SPHERICAL_TRUNCATED_POWER:
      return TruncatedPowerUnitBall(x,SPHERICAL_POWER_EXP,SPHERICAL_LOWER_TRUNCATE);
    case SPHERICAL_TABLE:
      return SphericalTable(x,SPHERICAL_TABLE_VALUES,SPHERICAL_TABLE_LENGTH);
    case SPHERICAL_TRUNCATED_GAUSSIAN:
      do
    {
      dw_NormalVector(x);
      r=Norm(x);
    }
      while ((r < SPHERICAL_LOWER_TRUNCATE) || (SPHERICAL_UPPER_TRUNCATE < r));
      return r;
    default:
      swz_fprintf_err("Unknown spherical type\n");
      swzExit(0);
    }
}

PRECISION LogSphericalDensity(PRECISION r)
{
  switch (SPHERICAL_TYPE)
    {
    case SPHERICAL_GAUSSIAN:
      return -0.5*r*r + SPHERICAL_CONSTANT;
    case SPHERICAL_UNIFORM:
      return (r > 1.0) ? MINUS_INFINITY : SPHERICAL_CONSTANT;
    case SPHERICAL_POWER:
      return (r > 1.0) ? MINUS_INFINITY : SPHERICAL_CONSTANT + (SPHERICAL_POWER_EXP - SPHERICAL_DIM)*log(r);
    case SPHERICAL_TRUNCATED_POWER:
      return ((r < SPHERICAL_LOWER_TRUNCATE) || (r > 1.0)) ? MINUS_INFINITY
                                                               : SPHERICAL_CONSTANT + (SPHERICAL_POWER_EXP - SPHERICAL_DIM)*log(r);
    case SPHERICAL_TABLE:
      return SPHERICAL_CONSTANT - (SPHERICAL_DIM - 1)*log(r) + LogSphericalTableDensity(r,SPHERICAL_TABLE_VALUES,SPHERICAL_TABLE_LENGTH);
    case SPHERICAL_TRUNCATED_GAUSSIAN:
      return ((r < SPHERICAL_LOWER_TRUNCATE) || (r > SPHERICAL_UPPER_TRUNCATE)) ? MINUS_INFINITY : -0.5*r*r + SPHERICAL_CONSTANT;
    default:
      swz_fprintf_err("Unknown spherical type\n");
      swzExit(0);
    }
}

/*
   The ith entry of the returned vector is the cumulative density evaluated at
   (i + 1) * max / bins.  The integer cum_bins controls the accuracy of the
   estimation.  The larger the value, the more accuate the estimate.
*/
TVector SphericalCumulativeDensity(PRECISION max, int bins, int cum_bins)
{
  TVector cumulative=CreateVector(bins);
  int i, j;
  PRECISION r, z, inc=max/(PRECISION)bins, cum_inc=inc/(PRECISION)cum_bins;
  for (i=0; i < bins; i++)
    {
      for (r=(PRECISION)i*inc + 0.5*cum_inc, z=MINUS_INFINITY, j=0; j < cum_bins; r+=cum_inc, j++)
    z=AddLogs_static(z,LogSphericalDensity(r) + (SPHERICAL_DIM - 1)*log(r));
      ElementV(cumulative,i)=exp(z - log(0.5) + 0.5*SPHERICAL_DIM*log(PI) - dw_log_gamma(0.5*SPHERICAL_DIM) + log(cum_inc));
    }
  for (i=1; i < bins; i++)
    ElementV(cumulative,i)+=ElementV(cumulative,i-1);
  return cumulative;
}

void TestSpherical(FILE *f, char *filename, PRECISION max)
{
  TMatrix cumulative;
  TVector x;
  int i, j, bins=1000, cum_bins=20, ndraws=1000000;
  PRECISION r, z, inc=max/(PRECISION)bins, cum_inc=inc/(PRECISION)cum_bins, s=1.0/(PRECISION)ndraws;
  FILE *f_out;

  cumulative=CreateMatrix(bins,3);
  for (i=0; i < bins; i++)
    {
      ElementM(cumulative,i,0)=(PRECISION)(i+1)*inc;
      for (r=(PRECISION)i*inc + 0.5*cum_inc, z=MINUS_INFINITY, j=0; j < cum_bins; r+=cum_inc, j++)
    z=AddLogs_static(z,LogSphericalDensity(r) + (SPHERICAL_DIM - 1)*log(r));
      ElementM(cumulative,i,1)=exp(z - log(0.5) + 0.5*SPHERICAL_DIM*log(PI) - dw_log_gamma(0.5*SPHERICAL_DIM) + log(cum_inc));
      ElementM(cumulative,i,2)=0.0;
    }

  x=CreateVector(SPHERICAL_DIM);
  for (i=ndraws; i > 0; i--)
    {
      r=DrawSpherical(x);
      if ((j=(int)floor(r/inc)) < bins)
    ElementM(cumulative,j,2)+=s;
    }
  FreeVector(x);

  for (i=1; i < bins; i++)
    {
      ElementM(cumulative,i,1)+=ElementM(cumulative,i-1,1);
      ElementM(cumulative,i,2)+=ElementM(cumulative,i-1,2);
    }

  f_out=!f ? dw_CreateTextFile(filename) : f;
  dw_PrintMatrix(f_out,cumulative,"%lf,");
  if (!f) fclose(f_out);
  FreeMatrix(cumulative);
}
#undef PI

#undef SPHERICAL_GAUSSIAN
#undef SPHERICAL_UNIFORM
#undef SPHERICAL_POWER
#undef SPHERICAL_TURNCATED_POWER
#undef SPHERICAL_TRUNCATED_GAUSSIAN
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*
   Assumes:
     x : m-vector

   Results:
     The vector x is filled with a vector drawn from the uniform distribution on
     the m dimensional solid unit sphere.

   Returns:
     Upon success, returns the norm of x, upon failure returns negative value.

   Notes:
     The vector is drawn by drawing a m-vector from the standard normal
     distribution and a real number u from the uniform distribution on [0,1], and
     normalizing the vector so its length equal to u^(1/m).
*/
PRECISION UniformUnitBall(TVector x)
{
  PRECISION r, s;
  if (!x)
    {
      dw_Error(NULL_ERR);
      return -1.0;
    }

  do
    dw_NormalVector(x);
  while ((s=Norm(x)) == 0.0);

  ProductSV(x,(r=pow(dw_uniform_rnd(),1.0/DimV(x)))/s,x);

  return r;
}

/*
   Assumes:
     x : n-vector

   Results:
     The vector x is filled with a vector drawn from the distribution

            0.5 * k * Gamma(n/2) * pi^(-n/2) * norm(x)^(k-n)

   Returns:
     norm(x) upon success and a negative value upon failure.

   Notes:
     If x is obtained by drawing y from the standard n-dimensional Gaussian
     distribtuion and r from the distribution on [0,1] with density

                                  k * r^(k-1)

     Since the cumulative density of r is

                                      r^k

     a draw of r can be obtained by drawing u from the uniform on [0,1] and
     defining r = u^(1/k).  This assumes that k > 0.
*/
PRECISION PowerUnitBall(TVector x, PRECISION k)
{
  PRECISION r, s;
  if (!x)
    {
      dw_Error(NULL_ERR);
      return -1.0;
    }

  do
    dw_NormalVector(x);
  while ((s=Norm(x)) == 0.0);

  ProductSV(x,(r=pow(dw_uniform_rnd(),1.0/k))/s,x);

  return r;
}

/*
   Assumes:
     x : n-vector

   Results:
     The vector x is filled with a vector drawn from the distribution

         0.5 * k * Gamma(n/2) * pi^(-n/2) * norm(x)^(k-n) / (1 - a^k)

   Returns:
     norm(x) upon success and a negative value upon failure.

   Notes:
     If x is obtained by drawing y from the standard n-dimensional Gaussian
     distribtuion and r from the distribution on [a,1] with density

                            k * r^(k-1) / (1 - a^k)

     Since the cumulative density of r is

                            (r^k - a^k) / (1 - a^k)

     a draw of r can be obtained by drawing u from the uniform on [0,1] and
     defining r = (u(1-a^k) + a^k)^(1/k).  This assumes that k > 0.
*/
PRECISION TruncatedPowerUnitBall(TVector x, PRECISION k, PRECISION a)
{
  PRECISION r, s, t;
  if (!x)
    {
      dw_Error(NULL_ERR);
      return -1.0;
    }

  do
    dw_NormalVector(x);
  while ((s=Norm(x)) == 0.0);

  t=pow(a,k);
  ProductSV(x,(r=pow(dw_uniform_rnd()*(1.0 - t) + t,1.0/k))/s,x);

  return r;
}


PRECISION SphericalTable(TVector x, PRECISION *table, int m)
{
  PRECISION r, s;
  int i, j;
  if (!x)
    {
      dw_Error(NULL_ERR);
      return -1.0;
    }

  do
    dw_NormalVector(x);
  while ((s=Norm(x)) == 0.0);

  j=(int)floor(dw_uniform_rnd()*(PRECISION)m);
  r=(j < m) ? table[j] + dw_uniform_rnd()*(table[j+1] - table[j]) : table[m];
  ProductSV(x,r/s,x);

  return r;
}

PRECISION LogSphericalTableDensity(PRECISION r, PRECISION *table, int m)
{
  int min=0, max=m, mid;
  if (r > table[m]) return MINUS_INFINITY;
  while (max - min > 1)
    if (r > table[mid=(min + max)/2])
      min=mid;
    else
      max=mid;
  return -log((PRECISION)m*(table[max] - table[min]));
}

