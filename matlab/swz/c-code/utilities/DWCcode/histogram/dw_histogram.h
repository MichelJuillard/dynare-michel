#ifndef __HISTOGRAMS__
#define __HISTOGRAMS__

#include "matrix.h"

#define HISTOGRAM_FIXED      1
#define HISTOGRAM_VARIABLE   2

/* Matrix histograms */
typedef struct
{
  TMatrix Min;
  TMatrix Max;
  int   **low;
  int   **high;
  int  ***freq;
  int     rows;
  int     cols;
  int     intervals;
  int     sample_size;
  int     type;
} TMatrixHistogram;

/* Vector histograms */
typedef struct
{
  TVector Min;
  TVector Max;
  int    *low;
  int    *high;
  int   **freq;
  int     dim;
  int     intervals;
  int     sample_size;
  int     type;
} TVectorHistogram;

/* Scalar histograms */
typedef struct
{
  PRECISION  Min;
  PRECISION  Max;
  int        low;
  int        high;
  int       *freq;
  int        intervals;
  int        sample_size;
  int        type;
} TScalarHistogram;

TMatrixHistogram *CreateMatrixHistogram(int rows, int cols, int intervals, int type);
void SetMaxMinMatrixHistogram(TMatrix Min, TMatrix Max, TMatrixHistogram *h);
void FreeMatrixHistogram(TMatrixHistogram *h);
void AddMatrixObservation(TMatrix X, TMatrixHistogram *h);
void MatrixPercentile(TMatrix X, PRECISION percentile, TMatrixHistogram *h);
void MatrixCumulative(TMatrix P, TMatrix Level, TMatrixHistogram *h);
TMatrix PlotMatrixHistogramAuto(int i, int j, int bins, TMatrixHistogram *h);
TMatrix PlotMatrixHistogram(int i, int j, PRECISION min, PRECISION max, int bins, TMatrixHistogram *h);

TVectorHistogram *CreateVectorHistogram(int dim, int intervals, int type);
void SetMaxMinVectorHistogram(TVector Min, TVector Max, TVectorHistogram *h);
void FreeVectorHistogram(TVectorHistogram *h);
void AddVectorObservation(TVector X, TVectorHistogram *h);
void VectorPercentile(TVector X, PRECISION percentile, TVectorHistogram *h);
void VectorCumulative(TVector p, TVector level, TVectorHistogram *h);
TMatrix PlotVectorHistogramAuto(int i, int bins, TVectorHistogram *h);
TMatrix PlotVectorHistogram(int i, PRECISION min, PRECISION max, int bins, TVectorHistogram *h);

TScalarHistogram *CreateScalarHistogram(int intervals, int type);
void SetMaxMinScalarHistogram(PRECISION Min, PRECISION Max, TScalarHistogram *h);
void FreeScalarHistogram(TScalarHistogram *h);
void AddScalarObservation(PRECISION x, TScalarHistogram *h);
PRECISION ScalarPercentile(PRECISION percentile, TScalarHistogram *h);
PRECISION ScalarCumulative(PRECISION level, TScalarHistogram *h);
TMatrix PlotScalarHistogramAuto(int bins, TScalarHistogram *h);
TMatrix PlotScalarHistogram(PRECISION min, PRECISION max, int bins, TScalarHistogram *h);
#endif
