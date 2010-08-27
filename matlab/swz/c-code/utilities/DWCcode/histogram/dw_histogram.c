
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "dw_histogram.h"
#include "dw_error.h"

static void Resize(PRECISION x, int *h, PRECISION *min, PRECISION *max, int intervals);
static void AddObservationVariable(PRECISION x, int *h, PRECISION *min, PRECISION *max, int intervals);
static void AddObservationFixed(PRECISION x, int *low, int *h, int *high, PRECISION min, PRECISION max, int intervals);

static PRECISION Cumulative(PRECISION level, int low, int *h, PRECISION min, PRECISION max, int intervals, int sample_size);
static PRECISION Percentile(PRECISION percentile, int low, int *h, PRECISION min, PRECISION max, int intervals, int sample_size);
static TMatrix MakeHistogram(int low, int *h, PRECISION min, PRECISION max,int intervals, int sample_size,
			     PRECISION min_out, PRECISION max_out, int bins);
static TMatrix MakeHistogramAuto(int low, int *h, int high, PRECISION min, PRECISION max, int intervals, int sample_size, int bins);


/*******************************************************************************
  The following set of routines create a matrix of histograms on the fly.
*******************************************************************************/
/*
   Assumes
     rows > 0
     cols > 0
     intervals > 0
     type = HISTOGRAM_FIXED or  HISTOGRAM_VARIABLE

   Results
     Creates and returns a matrix histogram data structure.  The size of the
     matrix is m x n and the number of intervals is intrvls.

*/
TMatrixHistogram *CreateMatrixHistogram(int rows, int cols, int intervals, int type)
{
 int i, j;
 TMatrixHistogram *h;

 if (!(h=(TMatrixHistogram *)malloc(sizeof(TMatrixHistogram)))) dw_Error(MEM_ERR);

 if (!(h->freq=(int***)malloc(rows*sizeof(int**)))) dw_Error(MEM_ERR);
 for (i=rows-1; i >= 0; i--)
  {
   if (!(h->freq[i]=(int**)malloc(cols*sizeof(int*)))) dw_Error(MEM_ERR);
   for (j=cols-1; j >= 0; j--)
    if (!(h->freq[i][j]=(int*)malloc(intervals*sizeof(int)))) dw_Error(MEM_ERR);
  }

 if (!(h->low=(int**)malloc(rows*sizeof(int*)))) dw_Error(MEM_ERR);
  for (i=rows-1; i >= 0; i--)
   if (!(h->low[i]=(int*)malloc(cols*sizeof(int)))) dw_Error(MEM_ERR);

 if (!(h->high=(int**)malloc(rows*sizeof(int*)))) dw_Error(MEM_ERR);
  for (i=rows-1; i >= 0; i--)
   if (!(h->high[i]=(int*)malloc(cols*sizeof(int)))) dw_Error(MEM_ERR);

 h->Min=CreateMatrix(rows,cols);
 h->Max=CreateMatrix(rows,cols);

 h->rows=rows;
 h->cols=cols;
 h->intervals=intervals;
 h->sample_size=0;
 h->type=type;

 return h;
}

void SetMaxMinMatrixHistogram(TMatrix Min, TMatrix Max, TMatrixHistogram *h)
{
 EquateMatrix(h->Min,Min);
 EquateMatrix(h->Max,Max);
 h->sample_size=0;
}

void FreeMatrixHistogram(TMatrixHistogram *h)
{
 int i, j;
 for (i=h->rows-1; i >= 0; i--)
  {
   for (j=h->cols-1; j >= 0; j--) free(h->freq[i][j]);
   free(h->freq[i]);
  }
 free(h->freq);
 for (i=h->rows-1; i >= 0; i--) free(h->low[i]);
 free(h->low);
 for (i=h->rows-1; i >= 0; i--) free(h->high[i]);
 free(h->high);
 FreeMatrix(h->Min);
 FreeMatrix(h->Max);
 free(h);
}

void AddMatrixObservation(TMatrix X, TMatrixHistogram *h)
{
 int i, j, k;

 if ((h->rows != RowM(X)) || (h->cols != ColM(X))) dw_Error(SIZE_ERR);

 if (h->sample_size <= 0)
   {
     for (i=h->rows-1; i >= 0; i--)
       for (j=h->cols-1; j >= 0; j--)
	 {
	   h->low[i][j]=h->high[i][j]=0;
	   for (k=h->intervals-1; k >= 0; k--) h->freq[i][j][k]=0;
	 }
     if (h->type == HISTOGRAM_VARIABLE)
       for (i=h->rows-1; i >= 0; i--)
	 for (j=h->cols-1; j >= 0; j--)
	   ElementM(h->Min,i,j)=ElementM(h->Max,i,j)=ElementM(X,i,j);
   }

 if (h->type == HISTOGRAM_FIXED)
   for (i=h->rows-1; i >= 0; i--)
     for (j=h->cols-1; j >= 0; j--)
       AddObservationFixed(ElementM(X,i,j),h->low[i]+j,h->freq[i][j],h->high[i]+j,ElementM(h->Min,i,j),ElementM(h->Max,i,j),h->intervals);
 else
   for (i=h->rows-1; i >= 0; i--)
     for (j=h->cols-1; j >= 0; j--)
       AddObservationVariable(ElementM(X,i,j),h->freq[i][j],&ElementM(h->Min,i,j),&ElementM(h->Max,i,j),h->intervals);

 h->sample_size++;
}

void MatrixPercentile(TMatrix X, PRECISION percentile, TMatrixHistogram *h)
{
 int i, j;

 if ((h->rows != RowM(X)) || (h->cols != ColM(X))) dw_Error(SIZE_ERR);

 for (i=h->rows-1; i >= 0; i--)
  for (j=h->cols-1; j >= 0; j--)
   ElementM(X,i,j)=Percentile(percentile,h->low[i][j],h->freq[i][j],ElementM(h->Min,i,j),ElementM(h->Max,i,j),h->intervals,h->sample_size);
}

/*
   Returns the probability that an observation is less than or equal to
   level.

   Assumes
     For 0 <= i < h->rows and 0 <= j < h->cols, let

       I[i][j][k]=(h->min[i][j] + k*inc[i][j], h->min[i][j] + (k+1)*inc[i][j]),

     where inc[i][j]=(h->max[i][j] - h->min[i][j])/h->samples_size.  The
     distribution is uniform on I[i][k][j] and

      P(h->min[i][j] + k*inc[i][j] < x[i][j] < h->min[i][j] + (k+1)*inc[i][j])
                          = h->freq[i][j][k]/h->sample_size.

     Furthermore,

          P(x[i][j] < h->min[i][j]) = 0 and P(x[i][j] > h->min[i][j]) = 0.

     In addition, if h->type == FIXED, then

              P(x[i][j] = h->min[i][j]) = h->low[i][j]/h->sample_size

     and

             P(x[i][j] = h->min[i][j]) = h->high[i][j]/h->sample_size.
*/
void MatrixCumulative(TMatrix P, TMatrix Level, TMatrixHistogram *h)
{
 int i, j;

 if ((h->rows != RowM(P)) || (h->cols != ColM(P)) ||
        (h->rows != RowM(Level)) || (h->cols != ColM(Level)))
  dw_Error(SIZE_ERR);

 for (i=h->rows-1; i >= 0; i--)
  for (j=h->cols-1; j >= 0; j--)
   ElementM(P,i,j)=Cumulative(ElementM(Level,i,j),h->low[i][j],h->freq[i][j],ElementM(h->Min,i,j),ElementM(h->Max,i,j),h->intervals,h->sample_size);
}

TMatrix PlotMatrixHistogramAuto(int i, int j, int bins, TMatrixHistogram *h)
{
  return MakeHistogramAuto(h->low[i][j],h->freq[i][j],h->high[i][j],ElementM(h->Min,i,j),ElementM(h->Max,i,j),h->intervals,h->sample_size,bins);
}

TMatrix PlotMatrixHistogram(int i, int j, PRECISION min, PRECISION max, int bins, TMatrixHistogram *h)
{
  return MakeHistogram(h->low[i][j],h->freq[i][j],ElementM(h->Min,i,j),ElementM(h->Max,i,j),h->intervals,h->sample_size,min,max,bins);
}

/*******************************************************************************
  The following set of routines create a vector of histograms on the fly.
*******************************************************************************/
TVectorHistogram *CreateVectorHistogram(int dim, int intervals, int type)
{
 int i;
 TVectorHistogram *h;

 if (!(h=(TVectorHistogram *)malloc(sizeof(TVectorHistogram))))
  dw_Error(MEM_ERR);

 if (!(h->freq=(int**)malloc(dim*sizeof(int*)))) dw_Error(MEM_ERR);
 for (i=dim-1; i >= 0; i--)
  if (!(h->freq[i]=(int*)malloc(intervals*sizeof(int)))) dw_Error(MEM_ERR);

 if (!(h->low=(int*)malloc(dim*sizeof(int)))) dw_Error(MEM_ERR);
 if (!(h->high=(int*)malloc(dim*sizeof(int)))) dw_Error(MEM_ERR);

 h->Min=CreateVector(dim);
 h->Max=CreateVector(dim);

 h->dim=dim;
 h->intervals=intervals;
 h->sample_size=0;
 h->type=type;

 return h;
}

void SetMaxMinVectorHistogram(TVector Min, TVector Max, TVectorHistogram *h)
{
 EquateVector(h->Min,Min);
 EquateVector(h->Max,Max);
 h->sample_size=0;
}

void FreeVectorHistogram(TVectorHistogram *h)
{
 int i;
 for (i=h->dim-1; i >= 0; i--) free(h->freq[i]);
 free(h->freq);
 free(h->low);
 free(h->high);
 FreeVector(h->Min);
 FreeVector(h->Max);
 free(h);
}

void AddVectorObservation(TVector x, TVectorHistogram *h)
{
 int i, k;

 if (h->dim != DimV(x)) dw_Error(SIZE_ERR);

 if (h->sample_size <= 0)
   {
     for (i=h->dim-1; i >= 0; i--)
       {
	 h->low[i]=h->high[i]=0;
	 for (k=h->intervals-1; k >= 0; k--) h->freq[i][k]=0;
       }
     if (h->type == HISTOGRAM_VARIABLE)
       for (i=h->dim-1; i >= 0; i--)
	 ElementV(h->Min,i)=ElementV(h->Max,i)=ElementV(x,i);
   }

 if (h->type == HISTOGRAM_FIXED)
   for (i=h->dim-1; i >= 0; i--)
     AddObservationFixed(ElementV(x,i),h->low+i,h->freq[i],h->high+i,ElementV(h->Min,i),ElementV(h->Max,i),h->intervals);
 else
   for (i=h->dim-1; i >= 0; i--)
     AddObservationVariable(ElementV(x,i),h->freq[i],&ElementV(h->Min,i),&ElementV(h->Max,i),h->intervals);

 h->sample_size++;
}

void VectorPercentile(TVector x, PRECISION percentile, TVectorHistogram *h)
{
 int i;

 if (h->dim != DimV(x)) dw_Error(SIZE_ERR);

 for (i=h->dim-1; i >= 0; i--)
  ElementV(x,i)=Percentile(percentile,h->low[i],h->freq[i],ElementV(h->Min,i),ElementV(h->Max,i),h->intervals,h->sample_size);
}



/*
   Returns the probability that an observation is less than or equal to
   level.

   Assumes
     For 0 <= i < h->dim, let

             I[i][k]=(h->min[i] + k*inc[i], h->min[i] + (k+1)*inc[i]),

     where inc[i]=(h->max[i] - h->min[i])/h->samples_size.  The distribution
     is uniform on I[i][k] and

       P(h->min[i] + k*inc[i] < x[i] < h->min[i] + (k+1)*inc[i])
                                              = h->freq[i][k]/h->sample_size.

     Furthermore,

            P(x[i] < h->min[i]) = 0 and P(x[i] > h->min[i]) = 0.

     In addition, if h->type == FIXED, then

                P(x[i] = h->min[i]) = h->low[i]/h->sample_size

     and

                P(x[i] = h->min[i]) = h->high[i]/h->sample_size.
*/
void VectorCumulative(TVector p, TVector level, TVectorHistogram *h)
{
 int i;

 if (h->dim != DimV(p) || (h->dim != DimV(level)))
  dw_Error(SIZE_ERR);

 for (i=h->dim-1; i >= 0; i--)
  ElementV(p,i)=Cumulative(ElementV(level,i),h->low[i],h->freq[i],ElementV(h->Min,i),ElementV(h->Max,i),h->intervals,h->sample_size);
}

TMatrix PlotVectorHistogramAuto(int i, int bins, TVectorHistogram *h)
{
  return MakeHistogramAuto(h->low[i],h->freq[i],h->high[i],ElementV(h->Min,i),ElementV(h->Max,i),h->intervals,h->sample_size,bins);
}

TMatrix PlotVectorHistogram(int i, PRECISION min, PRECISION max, int bins, TVectorHistogram *h)
{
  return MakeHistogram(h->low[i],h->freq[i],ElementV(h->Min,i),ElementV(h->Max,i),h->intervals,h->sample_size,min,max,bins);
}
/*******************************************************************************
  The following set of routines create a scalar histogram on the fly.
*******************************************************************************/
/*
   Assumes

   Results
     Creates and returns a scalar histogram data structure.
*/
TScalarHistogram *CreateScalarHistogram(int intervals, int type)
{
 TScalarHistogram *h;

 if (!(h=(TScalarHistogram *)malloc(sizeof(TScalarHistogram)))) dw_Error(MEM_ERR);

 if (!(h->freq=(int*)malloc(intervals*sizeof(int)))) dw_Error(MEM_ERR);

 h->intervals=intervals;
 h->sample_size=0;
 h->type=type;

 return h;
}

void SetMaxMinScalarHistogram(PRECISION Min, PRECISION Max, TScalarHistogram *h)
{
 h->Min=Min;
 h->Max=Max;
 h->sample_size=0;
}

void FreeScalarHistogram(TScalarHistogram *h)
{
 free(h->freq);
 free(h);
}

void AddScalarObservation(PRECISION x, TScalarHistogram *h)
{
 int k;

 if (h->sample_size <= 0)
  {
   h->low=h->high=0;
   for (k=h->intervals-1; k >= 0; k--) h->freq[k]=0;
   if (h->type == HISTOGRAM_VARIABLE) h->Min=h->Max=x;
  }

 if (h->type == HISTOGRAM_FIXED)
   AddObservationFixed(x,&(h->low),h->freq,&(h->high),h->Min,h->Max,h->intervals);
  else
   AddObservationVariable(x,h->freq,&(h->Min),&(h->Max),h->intervals);

 h->sample_size++;
}

PRECISION ScalarPercentile(PRECISION percentile, TScalarHistogram *h)
{
 return Percentile(percentile,h->low,h->freq,h->Min,h->Max,h->intervals,h->sample_size);
}

/*
   Returns the probability that an observation is less than or equal to
   level.

   Assumes
     Let

                I[k]=(h->min + k*inc, h->min + (k+1)*inc),

     where inc=(h->max - h->min)/h->samples_size.  The distribution
     is uniform on I[k] and

       P(h->min + k*inc < x < h->min + (k+1)*inc) = h->freq[k]/h->sample_size.

     Furthermore,

                 P(x < h->min) = 0 and P(x > h->min) = 0.

     In addition, if h->type == FIXED, then

                   P(x = h->min) = h->low/h->sample_size

     and

                  P(x = h->min) = h->high/h->sample_size.
*/
PRECISION ScalarCumulative(PRECISION level, TScalarHistogram *h)
{
 return Cumulative(level,h->low,h->freq,h->Min,h->Max,h->intervals,h->sample_size);
}

TMatrix PlotScalarHistogramAuto(int bins, TScalarHistogram *h)
{
  return MakeHistogramAuto(h->low,h->freq,h->Min,h->high,h->Max,h->intervals,h->sample_size,bins);
}

TMatrix PlotScalarHistogram(PRECISION min, PRECISION max, int bins, TScalarHistogram *h)
{
  return MakeHistogram(h->low,h->freq,h->Min,h->Max,h->intervals,h->sample_size,min,max,bins);
}

/*******************************************************************************/
/***************************** Low Level Routines ******************************/
/*******************************************************************************/
/*
   Resizes the histogram.  After resizing, it is guaranteed that *min <= x <= *max.
   The type of the histogram must be HISTOGRAM_VARIABLE.
*/
static void Resize(PRECISION x, int *h, PRECISION *min, PRECISION *max, int intervals)
{
 int i, j, k, m;
 if (x > *max)
   if (x - *min >= (PRECISION)intervals*(*max - *min))
     {
      for (i=1; i < intervals; i++)
       {
        h[0]+=h[i];
        h[i]=0;
       }
      *max=x;
     }
    else
     {
      m=(int)ceil((x - *min)/(*max - *min));
      for (i=j=0; i < intervals; j++)
       for(h[j]=h[i++], k=1; (k < m) && (i < intervals); k++)
        h[j]+=h[i++];
      for ( ; j < intervals; j++) h[j]=0;
      *max=*min + m*(*max - *min);
      if (x > *max) *max=x;
     }
  else
   if (x < *min)
    if (*max - x >= (PRECISION)intervals*(*max - *min))
      {
       for (j=intervals-1, i=intervals-2; i >= 0; i--)
        {
         h[j]+=h[i];
         h[i]=0;
        }
       *min=x;
      }
    else
     {
      m=(int)ceil((*max - x)/(*max - *min));
      for (i=j=intervals-1; i >= 0; j--)
       for(h[j]=h[i--], k=1; (k < m) && (i >= 0); k++)
        h[j]+=h[i--];
      for ( ; j >= 0; j--) h[j]=0;
      *min=*max - m*(*max - *min);
      if (x < *min) *min=x;
     }
}

/*
   Adds a observation to the histogram.  The type of the histogram must
   be HISTOGRAM_VARIABLE.
*/
static void AddObservationVariable(PRECISION x, int *h, PRECISION *min, PRECISION *max, int intervals)
{
 int i;

 if ((x < *min) || (x > *max)) Resize(x,h,min,max,intervals);

 if (*max > *min)
   {
    i=(int)(intervals*(x - *min)/(*max - *min));
    h[(i < intervals) ? i : intervals-1]++;
   }
  else
   h[0]++;
}

/*
   Adds a observation to the histogram.  The type of the histogram must
   be HISTOGRAM_FIXED.
*/
static void AddObservationFixed(PRECISION x, int *low, int *h, int *high, PRECISION min, PRECISION max, int intervals)
{
 PRECISION y=floor(intervals*(x - min)/(max - min));
 if (y < 0)
   (*low)++;
  else
   if (y < intervals)
     h[(int)y]++;
    else
     (*high)++;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
   Returns the level such that the probability of observing an observation
   less than or equal to level is percentile.  If there is a point mass at
   x, and P(y < x) <= percentile <= P(y <= x), then x is returned.

   Assumes
     Both intervals and sample_size are poitive and low and h[i] are
     non-negative.  Also if

           high = sample_size - (low + h[0] + ... + h[intervals - 1]),

     then high is non-negative.

     If min < max, let inc=(max - min)/intervals and define

                     I[k]=(min + k*inc, min + (k+1)*inc),

     The distribution is uniform on I[k] and

           P(min + k*inc < x < min + (k+1)*inc) = h[k]/sample_size.

     Furthermore, there are point masses at min and max with probability

                   P(x = min) = low/sample_size
     and
                   P(x = max) = high/sample_size.

     If min = max, then there is a single point mass at this point.
*/
static PRECISION Percentile(PRECISION percentile, int low, int *h, PRECISION min, PRECISION max, int intervals, int sample_size)
{
 int i;
 percentile=percentile*sample_size - low;
 if (percentile <= 0) return min;
 for (i=0; i < intervals; i++)
  if (h[i] && (percentile-=h[i]) <= 0)
   return min + ((PRECISION)(i+1) + percentile/(PRECISION)h[i])*(max - min)/(PRECISION)intervals;
 return max;
}

/*
   Returns the probability that an observation is less than or equal to
   level.

    Assumes
      Both intervals and sample_size are poitive and low and h[i] are
      non-negative.  Also, if

            high = sample_size - (low + h[0] + ... + h[intervals - 1]),

      then high is non-negative.

      If min < max, let inc=(max - min)/intervals and define

                      I[k]=(min + k*inc, min + (k+1)*inc),

      The distribution is uniform on I[k] and

            P(min + k*inc < x < min + (k+1)*inc) = h[k]/sample_size.

      Furthermore, there are point masses at min and max with probability

                    P(x = min) = low/sample_size
      and
                    P(x = max) = high/sample_size.

      If min = max, then there is a single point mass at this point
*/
static PRECISION Cumulative(PRECISION level, int low, int *h, PRECISION min, PRECISION max, int intervals, int sample_size)
{
 PRECISION inc=(max-min)/(PRECISION)intervals;
 int i, count;

 if (level < min) return 0.0;
 if (level >= max) return 1.0;

 for (count=low, i=0; i < intervals; count+=h[i++])
  if ((min+=inc) >= level)
   return ((PRECISION)count + (PRECISION)h[i]*(level - min + inc)/inc)/(PRECISION)sample_size;
 return 1.0;
}

/*
    Returns a histogram over the interval I=[min_out,max_out].  The matrix returned
    has bins rows and 2 columns.  If inc=(max_out - min_out)/bins, then the first
    element of the ith row is

                                 min + (i + 0.5)*inc,

    which is the mid-point of the ith interval.  The second element is

                        P(min + i*inc < x <= min + (i + 1)*inc)/inc,

    which is the average density over the ith interval.

    Assumes
     Both intervals and sample_size are poitive and low and h[i] are
     non-negative.  Also if

           high = sample_size - (low + h[0] + ... + h[intervals - 1]),

     then high is non-negative.

     If min < max, let inc=(max - min)/intervals and define

                     I[k]=(min + k*inc, min + (k+1)*inc),

     The distribution is uniform on I[k] and

           P(min + k*inc < x < min + (k+1)*inc) = h[k]/sample_size.

     Furthermore, there are point masses at min and max with probability

                   P(x = min) = low/sample_size
     and
                   P(x = max) = high/sample_size.

     If min = max, then there is a single point mass at this point.
*/
static TMatrix MakeHistogram(int low, int *h, PRECISION min, PRECISION max,
                      int intervals, int sample_size, PRECISION min_out, PRECISION max_out, int bins)
{
 int i;
 PRECISION inc, x, cdf_lower, cdf_upper;
 TMatrix X;

 inc=(max_out-min_out)/(PRECISION)bins;

 if (inc > 0)
   {
    X=CreateMatrix(bins,2);
    x=min_out+inc;

    cdf_lower=Cumulative(min_out,low,h,min,max,intervals,sample_size);

    for (i=0; i < bins; i++)
     {
      cdf_upper=Cumulative(x,low,h,min,max,intervals,sample_size);

      ElementM(X,i,0)=x - 0.5*inc;
      ElementM(X,i,1)=(cdf_upper-cdf_lower)/inc;

      cdf_lower=cdf_upper;
      x+=inc;
     }
   }
 else
   return (TMatrix)NULL;

 return X;
}

/*
    Automatically chooses lenth of interval over which to produce histogram and
    then calls MakeHistogram().
*/
static TMatrix MakeHistogramAuto(int low, int *h, int high, PRECISION min, PRECISION max, int intervals, int sample_size, int bins)
{
 PRECISION inc=(max-min)/intervals, max_out, min_out;
 int lo, hi;

 if ((low == sample_size) || (inc <= 0))
   {
     min_out=min-1.0;
     max_out=min+1.0;
   }
 else
   {
     if (low > 0)
       lo=-1;
     else
       for (lo=0; (lo < intervals) && !h[lo]; lo++);

     if (lo == intervals)
       {
	 min_out=max-1.0;
	 max_out=max+1.0;
       }
     else
       {
	 if (high > 0)
	   hi=intervals;
	 else
	   for (hi=intervals-1; !h[hi]; hi--);

	 if (lo >= 0)
	   if (hi < intervals)
	     {
	       min_out=min+lo*inc;
	       max_out=min+(hi+1)*inc;
	     }
	   else
	     {
	       min_out=min+lo*inc;
	       if (bins == 1)
		 max_out=(1+SQRT_MACHINE_EPSILON)*max;
	       else
		 {
		   inc=(1-SQRT_MACHINE_EPSILON)*(max - min_out)/(PRECISION)(bins-1);
		   max_out=max + inc;
		 }
	     }
	 else
	   if (hi < intervals)
	     {
	       max_out=min+(hi+1)*inc;
	       if (bins == 1)
		 min_out=(1-SQRT_MACHINE_EPSILON)*min;
	       else
		 {
		   inc=(1-SQRT_MACHINE_EPSILON)*(max_out - min)/(PRECISION)(bins-1);
		   min_out=min - inc;
		 }
	     }
	   else
	     if (bins <= 2)
	       {
		 min_out=(1-SQRT_MACHINE_EPSILON)*min;
		 max_out=(1+SQRT_MACHINE_EPSILON)*max;
	       }
	     else
	       {
		 inc=(1-SQRT_MACHINE_EPSILON)*(max_out - min)/(PRECISION)(bins-2);
		 min_out=min - inc;
		 max_out=max +inc;
	       }
       }
   }

 return MakeHistogram(low,h,min,max,intervals,sample_size,min_out,max_out,bins);
}


