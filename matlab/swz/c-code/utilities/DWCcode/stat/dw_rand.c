
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include "prcsn.h"
#include "dw_rand.h"
#include "dw_error.h"

#include "modify_for_mex.h"

/*  //=== Static routines ===   ansi-c*/
static void gser(PRECISION *gamser, PRECISION a, PRECISION x, PRECISION *gln);
static void gcf(PRECISION *gammcf, PRECISION a, PRECISION x, PRECISION *gln);
static PRECISION gammp(PRECISION a, PRECISION x);

/*******************************************************************************/
/*************************** Uniform Random Numbers ****************************/
/*******************************************************************************/
/*
   Flag controling which uniform random number to choose
*/
/*  //#define USE_NR1_RNG   ansi-c*/
#define USE_NR2_RNG
/*  //#define USE_IMSL_RNG   ansi-c*/

#if defined (USE_IMSL_RNG)
#include <imsls.h>
#elif defined(USE_NR1_RNG)
#define NTAB 32
static int idum=-1;
static int iy=0;
static int iv[NTAB];
#elif defined(USE_NR2_RNG)
#define NTAB 32
static int idum=-1;
static int idum2=123456789;
static int iy=0;
static int iv[NTAB];
#endif

/*
   Initializes seed value for uniform random number generator.  The seed value
   can be any integer.  A value of 0 will initialize the seed from the system
   clock for the Numerical Recipies algorithms.
*/
void dw_initialize_generator(int init)
{
#ifdef USE_IMSL_RNG
 imsls_random_option(7);
 imsls_random_seed_set((init < 0) ? -init : init);
#else
 if (init)
   idum=(init > 0) ? -init : init;
 else
   {
     idum=0;
     idum=(int)(-INT_MAX*dw_uniform_rnd());
   }
#endif
}

/*
   Allocates memory and saves the state of the random number generator.  The
   calling routine is responsible for freeing the returned memory.
*/
void* dw_get_generator_state(void)
{
#if defined(USE_IMSL_RNG)
  int *state=(int*)NULL;
  if (state=(int*)swzMalloc(1566*sizeof(int)))
    {
      imsls_random_GFSR_table_get(&state,IMSLS_RETURN_USER,state,0);
      state[1565]=imsls_random_seed_get();
    }
  return state;
#elif defined (USE_NR1_RNG)
  int *state=(int*)NULL;
  if (state=(int*)swzMalloc((NTAB+2)*sizeof(int)))
    {
      memcpy(state,iv,NTAB*sizeof(int));
      state[NTAB]=iy;
      state[NTAB+1]=idum;
    }
  return state;
#elif defined (USE_NR2_RNG)
  int *state=(int*)NULL;
  if (state=(int*)swzMalloc((NTAB+3)*sizeof(int)))
    {
      memcpy(state,iv,NTAB*sizeof(int));
      state[NTAB]=iy;
      state[NTAB+1]=idum;
      state[NTAB+2]=idum2;
    }
  return state;
#endif
}

/*
   Returns the size in bytes of the void pointer returned by
   dw_get_generator_state().
*/
int dw_get_generator_state_size(void)
{
#if defined(USE_IMSL_RNG)
  return 1566*sizeof(int);
#elif defined (USE_NR1_RNG)
  return (NTAB+2)*sizeof(int);
#elif defined (USE_NR2_RNG)
  return (NTAB+3)*sizeof(int);
#endif
}

/*
   Sets the state of the random number generator.  The void pointer must have
   been obtained via a call to dw_get_generator_state().
*/
void dw_set_generator_state(void *state)
{
#if defined(USE_IMSL_RNG)
  imsls_random_GFSR_table_set((int*)state);
  imsls_random_seed_set(((int*)state)[1565]);
#elif defined (USE_NR1_RNG)
  memcpy(iv,state,NTAB*sizeof(int));
  iy=((int*)state)[NTAB];
  idum=((int*)state)[NTAB+1];
#elif defined (USE_NR2_RNG)
  memcpy(iv,state,NTAB*sizeof(int));
  iy=((int*)state)[NTAB];
  idum=((int*)state)[NTAB+1];
  idum2=((int*)state)[NTAB+2];
#endif
}

void dw_print_generator_state(FILE *f)
{
  if (f)
    {
#if defined(USE_IMSL_RNG)
      int i, *state;
      if (state=dw_get_generator_state())
    {
      for (i=0; i < 1566; i++) fprintf(f,"%d ",state[i]);
      fprintf(f,"\n");
      free(state);
    }
#elif defined (USE_NR1_RNG)
      int i, *state;
      if (state=dw_get_generator_state())
    {
      for (i=0; i < NTAB+2; i++) fprintf(f,"%d ",state[i]);
      fprintf(f,"\n");
      free(state);
    }
#elif defined (USE_NR2_RNG)
      int i, *state;
      if (state=dw_get_generator_state())
    {
      for (i=0; i < NTAB+3; i++) fprintf(f,"%d ",state[i]);
      fprintf(f,"\n");
      free(state);
    }
#endif
  }
}
void dw_read_generator_state(FILE *f)
{
  if (f)
    {
#if defined(USE_IMSL_RNG)
  int i, *state;
  if (state=(int*)swzMalloc(1566*sizeof(int)))
    {
      for (i=0; i < 1566; i++) fscanf(f," %d ",state+i);
      dw_set_generator_state(state);
      free(state);
    }
#elif defined (USE_NR1_RNG)
  int i, *state;
  if (state=(int*)swzMalloc((NTAB+2)*sizeof(int)))
    {
      for (i=0; i < NTAB+2; i++) fscanf(f," %d ",state+i);
      dw_set_generator_state(state);
      free(state);
    }
#elif defined (USE_NR2_RNG)
  int i, *state;
  if (state=(int*)swzMalloc((NTAB+3)*sizeof(int)))
    {
      for (i=0; i < NTAB+3; i++) fscanf(f," %d ",state+i);
      dw_set_generator_state(state);
      free(state);
    }
#endif
  }
}

/*
   The following code is adapted from Numerical Recipes in C.  This rnd1() from
   that text.
*/
#ifdef USE_NR1_RNG
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-MACHINE_EPSILON)
PRECISION dw_uniform_rnd(void)
{
  int j, k;
  PRECISION temp;

  if (idum <= 0)
    {
      if (idum == 0)
        {
          if(constant_seed==0)
            idum=abs((int)time((time_t *)NULL));
          else
            {
              srand(constant_seed);
              idum=rand();
            }
          if (idum == 0) idum=1;
        }
      else
        idum=-idum;

      for (j=NTAB+7; j >= 0; j--)
    {
      k=idum/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if (idum < 0) idum+=IM;
      if (j < NTAB) iv[j]=idum;
    }
      iy=iv[0];
    }
  k=idum/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=idum;
  return ((temp=(PRECISION)(AM*iy)) > RNMX) ? (PRECISION)RNMX : temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NDIV
#undef RNMX
#endif

/*
   The following code is adapted from Numerical Recipes in C.  This rnd2() from
   that text.
*/
#ifdef USE_NR2_RNG
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-MACHINE_EPSILON)
PRECISION dw_uniform_rnd(void)
{
  int j, k;
  PRECISION temp;

  if (idum <= 0)
    {
      if (idum == 0)
        {
          if(constant_seed==0)
            idum=abs((int)time((time_t *)NULL));
          else
            {
              srand(constant_seed);
              idum=rand();
            }
          if (idum == 0) idum=1;
        }
      else
        idum=-idum;

      idum2=idum;
      for (j=NTAB+7; j>=0; j--)
    {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum += IM1;
      if (j < NTAB) iv[j] = idum;
    }
      iy=iv[0];
    }
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  return ((temp=AM*iy) > RNMX) ? RNMX : temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NDIV
#undef RNMX
#endif

#ifdef USE_IMSL_RNG
PRECISION dw_uniform_rnd(void)
{
  PRECISION x;
#if PRECISION_SIZE == 8
  imsls_d_random_uniform(1,IMSLS_RETURN_USER,&x,0);
#else
  imsls_f_random_uniform(1,IMSLS_RETURN_USER,&x,0);
#endif
  return x;
}
#endif

#if defined (USE_IMSL_RNG)
#undef USE_IMSL_RNG
#elif defined (USE_NR1_RNG)
#undef NTAB
#undef USE_NR1_RNG
#elif defined (USE_NR2_RNG)
#undef NTAB
#undef USE_NR2_RNG
#endif

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*
   Returns a standard gaussian deviate.  The density function for the
   standard gaussian is

                          1
                     ----------- exp(-0.5*x^2)
                      sqrt(2*Pi)

*/
PRECISION dw_gaussian_rnd(void)
{
  static int iset=0;
  static PRECISION gset;
  PRECISION fac,r,v1,v2;

  if  (iset == 0)
    {
      do
    {
      v1=2.0*dw_uniform_rnd()-1.0;
      v2=2.0*dw_uniform_rnd()-1.0;
      r=v1*v1+v2*v2;
    }
      while (r >= 1.0);
      fac=sqrt(-2.0*log(r)/r);
      gset=v1*fac;
      iset=1;
      return v2*fac;
    }
  else
    {
      iset=0;
      return gset;
    }
}

#undef PI
/*
   Returns a standard gamma deviate.  The density function for a standard gamma
   distribution is

                                           x^(a-1)*exp(-x)
                   gamma_density(x;a) =   ----------------
                                              gamma(a)

   for a > 0.  The function gamma(a) is the integral with from 0 to infinity of
   exp(-t)*t^(a-1).

   When a = 1.0, then gamma is exponential. (Devroye, page 405).
   When a < 1.0, Johnk's generator (Devroye, page 418).
   When a > 1.0, a rejection method or Best's algorithm (Devroye, page 410).

   A general gamma variate can be obtained as follows.  Let z=b*x.  Then,
   z is drawn from a general gamma distribution whose density is

                                        z^(a-1)*exp(-z/b)
                gamma_density(z;a,b) = ------------------
                                          gamma(a)*b^a

   Uses algorithm translated by Iskander Karibzhanov from the Matlab function
   gamrnd.m, which follows Johnk's generator in Devroye ("Non-Uniform Random
   Variate Generation", Springer-Verlag, 1986, page 418).

   Notes:
    Does not check if a > 0.
*/
PRECISION dw_gamma_rnd(PRECISION a)
{
  PRECISION b, u, v, w, x, y, z;

   if (a == 1.0) return -log(dw_uniform_rnd());

   if (a < 1.0)
     {
       u=1.0/a;
       v=1.0/(1.0-a);
       do
     {
       x=pow(dw_uniform_rnd(),u);
       y=pow(dw_uniform_rnd(),v);
     }
       while (x+y > 1.0);
       return -log(dw_uniform_rnd())*x/(x+y);
     }

   b=a - 1.0;
   while(1)
     {
       u=dw_uniform_rnd();
       w=u*(1.0 - u);
       y=sqrt((3.0*a - 0.75)/w)*(u - 0.5);
       x=b + y;
       if (x > 0.0)
     {
       v=dw_uniform_rnd();
       z=64.0*w*w*w*v*v;
       if ((z <= 1.0 - 2.0*y*y/x) || (log(z) <= 2.0*(b*log(x/b) - y)))
         return x;
     }
     }
}

/*
   Returns a lognormal deviate.  The mean and standard deviations of the
   underlying normal distributions are passed.
*/
PRECISION dw_lognormal_rnd(PRECISION mean, PRECISION standard_deviation)
{
  return exp(standard_deviation * dw_gaussian_rnd() + mean);
}


/*
   Returns the integral from -infinity to x of 1/sqrt(2*PI)*exp(-y^2/2).
   Routine adapted from Numerical Recipes in C.
*/
double dw_normal_cdf(double x)
{
 double z=fabs(0.7071067811865*x), t=2.0/(2.0+z);

 return (x > 0) ?
           1.0-0.5*t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
             t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
               t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
                :
           0.5*t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
             t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
               t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));

}

PRECISION dw_chi_square_cdf(PRECISION x, int df)
{
  return gammp(0.5*df,0.5*x);
}

#define MAXITER 1000
PRECISION dw_chi_square_invcdf(PRECISION p, int df)
{
  int i;
  PRECISION p_lo=p-SQRT_MACHINE_EPSILON, p_hi=p+SQRT_MACHINE_EPSILON, hi, lo=0.0, mid, cdf;
  if (p <= 0)
    {
      if (p < 0) dw_Error(ARG_ERR);
      return 0.0;
    }
  else
    if (p >= 1)
      {
        if (p > 1) dw_Error(ARG_ERR);
        return PLUS_INFINITY;
      }
  if ((cdf=dw_chi_square_cdf(hi=2*df,df)) < p_lo)
    {
      for (lo=hi, i=MAXITER; (i > 0) && ((cdf=dw_chi_square_cdf(hi*=2,df)) < p_lo); lo=hi, i--);
      if (i == 0)
        {
          dw_Error(ITERATION_ERR);
          return PLUS_INFINITY;
        }
    }
  if (cdf < p_hi) return hi;
  for (i=MAXITER; i > 0; i--)
    if ((cdf=dw_chi_square_cdf(mid=0.5*(lo+hi),df)) < p_lo)
      lo=mid;
    else
      if (cdf > p_hi)
        hi=mid;
      else
        return mid;
  return 0.5*(lo+hi);
}
#undef MAXITER

/*
   Returns the natural logrithm of the gamma function applied to x.  The gamma
   function of x is the integral from 0 to infinity of t^(x-1)*exp(-t)dt.

   Routine adapted from the gammln routine from Numerical Recipes in C.
*/
PRECISION dw_log_gamma(PRECISION x)
{
  static PRECISION cof[6]={ 76.18009172947146,  -86.50532032941677,
                24.01409824083091,  -1.231739572450155,
                0.1208650973866179e-2, -0.5395239384953e-5};
  PRECISION y, z, ser;
  int j;
  z=x+5.5;
  z-=(x+0.5)*log(z);
  ser=1.000000000190015;
  for (y=x, j=0; j <= 5; j++) ser+=cof[j]/++y;
  return -z+log(2.5066282746310005*ser/x);
}

/******************************************************************************/
/************************** Numerical Recipies in C ***************************/
/******************************************************************************/
#define ITMAX 1000
#define EPS 3.0e-7
static void gser(PRECISION *gamser, PRECISION a, PRECISION x, PRECISION *gln)
{
  int n;
  PRECISION sum,del,ap;

  dw_ClearError();
  *gln=dw_log_gamma(a);
  if (x <= 0.0)
    {
      if (x < 0.0)
        dw_Error(ARG_ERR);
      else
        *gamser=0.0;
    }
  else
    {
      ap=a;
      del=sum=1.0/a;
      for (n=1; n <= ITMAX; n++)
        {
          ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS)
        {
          *gamser=sum*exp(-x+a*log(x)-(*gln));
          return;
        }
    }
      dw_Error(ITERATION_ERR);
    }
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software */

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
static void gcf(PRECISION *gammcf, PRECISION a, PRECISION x, PRECISION *gln)
{
  int i;
  PRECISION an,b,c,d,del,h;

  *gln=dw_log_gamma(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1; i <= ITMAX; i++)
    {
      an = -i*(i-a);
      b += 2.0;
      d=an*d+b;
      if (fabs(d) < FPMIN) d=FPMIN;
      c=b+an/c;
      if (fabs(c) < FPMIN) c=FPMIN;
      d=1.0/d;
      del=d*c;
      h *= del;
      if (fabs(del-1.0) < EPS) break;
    }
  if (i > ITMAX)
    dw_Error(ITERATION_ERR);
  else
    {
      *gammcf=exp(-x+a*log(x)-(*gln))*h;
      dw_ClearError();
    }
}
#undef ITMAX
#undef EPS
#undef FPMIN
/* (C) Copr. 1986-92 Numerical Recipes Software */

static PRECISION gammp(PRECISION a, PRECISION x)
{
  PRECISION gamser,gammcf,gln;

  if (x < 0.0 || a <= 0.0)
    {
      dw_Error(ARG_ERR);
      return 0.0;
    }
  dw_ClearError();
  if (x < (a+1.0))
    {
      gser(&gamser,a,x,&gln);
      return gamser;
    }
   else
    {
      gcf(&gammcf,a,x,&gln);
      return 1.0-gammcf;
    }
}
/* (C) Copr. 1986-92 Numerical Recipes Software */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
