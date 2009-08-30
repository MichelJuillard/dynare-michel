
#include "dw_matrix_sort.h"
#include "dw_error.h"
#include <stdlib.h>
#include <string.h>

static void b_qsort_array_ascending_real(PRECISION *x, int m);
static void b_qsort_array_descending_real(PRECISION *x, int m);
static void b_qsort_matrix_columns_ascending_real(PRECISION *x, int m, int n, int idx);
static void b_qsort_matrix_columns_descending_real(PRECISION *x, int m, int n, int idx);
static void b_qsort_matrix_rows_ascending_real(PRECISION *x, int m, int n, int br, int er, int idx);
static void b_qsort_matrix_rows_descending_real(PRECISION *x, int m, int n, int br, int er, int idx);

/*
   Assumes
     X : m x n matrix or null 
     Y : m x n matrix
     j : column to sort

   Results
     The rows of X are sorted in ascending order on the ith column.  The matrix 
     X is created if null.

   Returns
     Returns X upon success and null on failure.  Call GetError() to
     determine the cause of failure.

   Notes
     X and Y do not have to be distinct matrices.  Uses the quick sort algorithm,
*/
TMatrix SortMatrixRowsAscending(TMatrix X, TMatrix Y, int j)
{
  if (!Y)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  if ((X != Y) && !(X=EquateMatrix(X,Y)))
    return (TMatrix)NULL;
  if (MajorForm(X) == ROW_MAJOR)
    b_qsort_matrix_columns_ascending_real(pElementM(X),ColM(X),RowM(X),j);
  else
    b_qsort_matrix_rows_ascending_real(pElementM(X),RowM(X),ColM(X),0,RowM(X)-1,j*RowM(X));
  return X;
}

/*
   Assumes
     X : m x n matrix or null 
     Y : m x n matrix
     j : column to sort

   Results
     The rows of X are sorted in descending order on the ith column.  The matrix 
     X is created if null.

   Returns
     Returns X upon success and null on failure.  Call GetError() to
     determine the cause of failure.

   Notes
     X and Y do not have to be distinct matrices.  Uses the quick sort algorithm,
*/
TMatrix SortMatrixRowsDescending(TMatrix X, TMatrix Y, int j)
{
  if (!Y)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  if ((X != Y) && !(X=EquateMatrix(X,Y)))
    return (TMatrix)NULL;
  if (MajorForm(X) == ROW_MAJOR)
    b_qsort_matrix_columns_descending_real(pElementM(X),ColM(X),RowM(X),j);
  else
    b_qsort_matrix_rows_descending_real(pElementM(X),RowM(X),ColM(X),0,RowM(X)-1,j*RowM(X));
  return X;
}

/*
   Assumes
     X : m x n matrix or null 
     Y : m x n matrix
     i : row to sort

   Results
     The columns of X are sorted in ascending order on the ith row.  The matrix X
     is created if null.

   Returns
     Returns X upon success and null on failure.  Call GetError() to
     determine the cause of failure.

   Notes
     X and Y do not have to be distinct matrices.  Uses the quick sort algorithm,
*/
TMatrix SortMatrixColumnsAscending(TMatrix X, TMatrix Y, int i)
{
  if (!Y)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  if ((X != Y) && !(X=EquateMatrix(X,Y)))
    return (TMatrix)NULL;
  if (MajorForm(X) == ROW_MAJOR)
    b_qsort_matrix_rows_ascending_real(pElementM(X),ColM(X),RowM(X),0,ColM(X)-1,i*RowM(X));
  else 
    b_qsort_matrix_columns_ascending_real(pElementM(X),RowM(X),ColM(X),i);
  return X;
}

/*
   Assumes
     X : m x n matrix or null 
     Y : m x n matrix
     i : row to sort

   Results
     The columns of X are sorted in descending order on the ith row.  The matrix 
     X is created if null.

   Returns
     Returns X upon success and null on failure.  Call GetError() to
     determine the cause of failure.

   Notes
     X and Y do not have to be distinct matrices.  Uses the quick sort algorithm,
*/
TMatrix SortMatrixColumnsDescending(TMatrix X, TMatrix Y, int i)
{
  if (!Y)
    {
      dw_Error(NULL_ERR);
      return (TMatrix)NULL;
    }
  if ((X != Y) && !(X=EquateMatrix(X,Y)))
    return (TMatrix)NULL;
  if (MajorForm(X) == ROW_MAJOR)
    b_qsort_matrix_rows_descending_real(pElementM(X),ColM(X),RowM(X),0,ColM(X)-1,i*RowM(X));
  else 
    b_qsort_matrix_columns_descending_real(pElementM(X),RowM(X),ColM(X),i);
  return X;
}

/*
   Assumes
     x : m vector or null 
     y : m vector

   Results
     The vector x is sorted in ascending order.  The vector x is created if 
     null.

   Returns
     Returns x upon success and null on failure.  Call GetError() to
     determine the cause of failure.

   Notes
     x and x do not have to be distinct vectors.  Uses the quick sort algorithm,
*/
TVector SortVectorAscending(TVector x, TVector y)
{
  if (!y)
    {
      dw_Error(NULL_ERR);
      return (TVector)NULL;
    }
  if ((x != x) && !(x=EquateVector(x,y)))
    return (TVector)NULL;
  b_qsort_array_ascending_real(pElementV(x),DimV(x));
  return x;
}

/*
   Assumes
     x : m vector or null 
     y : m vector

   Results
     The vector x is sorted in descending order.  The vector x is created if 
     null.

   Returns
     Returns x upon success and null on failure.  Call GetError() to
     determine the cause of failure.

   Notes
     x and x do not have to be distinct vectors.  Uses the quick sort algorithm,
*/
TVector SortVectorDescending(TVector x, TVector y)
{
  if (!y)
    {
      dw_Error(NULL_ERR);
      return (TVector)NULL;
    }
  if ((x != x) && !(x=EquateVector(x,y)))
    return (TVector)NULL;
  b_qsort_array_descending_real(pElementV(x),DimV(x));
  return x;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/*
   Assumes:
     x - array of length m

   Results:
     x is sorted in ascending order

   Notes:
     Uses the quick sort mean algorithm.  Switches to insertion sort when the
     size of the list is 10 or less.
*/
static void b_qsort_array_ascending_real(PRECISION *x, int m)
{
  PRECISION y, c;
  int j, k;
  if (m > 10)
    {
      // quick sort
      m--;

      if (x[0] == x[m])
	c=x[0];
      else
	{
	  if (x[0] > x[m])
	    { y=x[m]; x[m]=x[0]; x[0]=y; }
	  c=0.5*(x[0] + x[m]);
	}

      for (j=1; (j < m) && (x[j] <= c); j++);
      for (k=m-1; (k > 0) && (x[k] >= c); k--);
      while (j < k)
	{
	  y=x[j]; x[j]=x[k]; x[k]=y;
	  while (x[j] <= c) j++;
	  while (x[k] >= c) k--;
	}
      if (k > 0)
	b_qsort_array_ascending_real(x,k+1);
      if (j < m)
	b_qsort_array_ascending_real(x+j,m-j+1);
    }
  else
    {
      // insertion sort
      for (j=1; j < m; j++)
	{
	  y=x[j];
	  for (k=j-1; k >= 0; k--)
	    if (x[k] <= y)
	      break;
	    else
	      x[k+1]=x[k];
	  x[k+1]=y;
	}
    }
}

/*
   Assumes:
     x - array of length m

   Results:
     x is sorted in descending order

   Notes:
     Uses the quick sort mean algorithm.  Switches to insertion sort when the
     size of the list is 10 or less.
*/
static void b_qsort_array_descending_real(PRECISION *x, int m)
{
  PRECISION y, c;
  int j, k;
  if (m > 10)
    {
      // quick sort
      m--;

      if (x[0] == x[m])
	c=x[0];
      else
	{
	  if (x[0] < x[m])
	    { y=x[m]; x[m]=x[0]; x[0]=y; }
	  c=0.5*(x[0] + x[m]);
	}

      for (j=1; (j < m) && (x[j] >= c); j++);
      for (k=m-1; (k > 0) && (x[k] <= c); k--);
      while (j < k)
	{
	  y=x[j]; x[j]=x[k]; x[k]=y;
	  while (x[j] >= c) j++;
	  while (x[k] <= c) k--;
	}
      if (k > 0)
	b_qsort_array_descending_real(x,k+1);
      if (j < m)
	b_qsort_array_descending_real(x+j,m-j+1);
    }
  else
    {
      // insertion sort
      for (j=1; j < m; j++)
	{
	  y=x[j];
	  for (k=j-1; k >= 0; k--)
	    if (x[k] >= y)
	      break;
	    else
	      x[k+1]=x[k];
	  x[k+1]=y;
	}
    }

}

/*
   Assumes:
     x - array of length m*n in colum major format.
     m - number of rows
     n - number of columns

   Results:
     The columns of x are sorted in ascending order on row idx.

   Notes:
     Uses the quick sort mean algorithm.  Switches to insertion sort when the
     size of the list is 10 or less.  If the matrix is in row major format, then
     m is the number of columns, n is the number of rows, and the rows of x are
     sorted in ascending order on column idx. 
*/
static void b_qsort_matrix_columns_ascending_real(PRECISION *x, int m, int n, int idx)
{
  PRECISION *y, c;
  int j, k, p, s;
  y=(PRECISION*)malloc(s=m*sizeof(PRECISION));
  if (n > 10)
    {
      // quick sort
      p=(n-1)*m;
      k=p+idx;

      if (x[idx] == x[k])
	c=x[idx];
      else
	{
	  if (x[idx] > x[k])
	    { memcpy(y,x+p,s); memcpy(x+p,x,s); memcpy(x,y,s); }
	  c=0.5*(x[idx] + x[k]);
	}

      for (j=m+idx; (j < p) && (x[j] <= c); j+=m);
      for (k-=m; (k > idx) && (x[k] >= c); k-=m);
      while (j < k)
	{
	  memcpy(y,x+j-idx,s); memcpy(x+j-idx,x+k-idx,s); memcpy(x+k-idx,y,s);
	  while (x[j] <= c) j+=m;
	  while (x[k] >= c) k-=m;
	}
      if (k > idx)
	b_qsort_matrix_columns_ascending_real(x,m,(k-idx)/m+1,idx);
      if (j < p)
	b_qsort_matrix_columns_ascending_real(x+j-idx,m,n-(j-idx)/m,idx);
    }
  else
    {
      // insertion sort
      p=n*m;
      for (j=m+idx; j < p; j+=m)
	if (x[j-m] > x[j])
	  {
	    memcpy(y,x+j-idx,s);
	    memcpy(x+j-idx,x+j-m-idx,s);
	    for (k=j-m-m; k >= 0; k-=m)
	      if (x[k] <= y[idx])
		break;
	      else
		memcpy(x+k+m-idx,x+k-idx,s);
	    memcpy(x+k+m-idx,y,s);
	  }
    }
  free(y);
}

/*
   Assumes:
     x - array of length m*n in colum major format.
     m - number of rows
     n - number of columns

   Results:
     The columns of x are sorted in ascending order on row idx.

   Notes:
     Uses the quick sort mean algorithm.  Switches to insertion sort when the
     size of the list is 10 or less.  If the matrix is in row major format, then
     m is the number of columns, n is the number of rows, and the rows of x are
     sorted in ascending order on column idx. 
*/
static void b_qsort_matrix_columns_descending_real(PRECISION *x, int m, int n, int idx)
{
  PRECISION *y, c;
  int j, k, p, s;
  y=(PRECISION*)malloc(s=m*sizeof(PRECISION));
  if (n > 10)
    {
      // quick sort
      p=(n-1)*m;
      k=p+idx;

      if (x[idx] == x[k])
	c=x[idx];
      else
	{
	  if (x[idx] < x[k])
	    { memcpy(y,x+p,s); memcpy(x+p,x,s); memcpy(x,y,s); }
	  c=0.5*(x[idx] + x[k]);
	}

      for (j=m+idx; (j < p) && (x[j] >= c); j+=m);
      for (k-=m; (k > idx) && (x[k] <= c); k-=m);
      while (j < k)
	{
	  memcpy(y,x+j-idx,s); memcpy(x+j-idx,x+k-idx,s); memcpy(x+k-idx,y,s);
	  while (x[j] >= c) j+=m;
	  while (x[k] <= c) k-=m;
	}
      if (k > idx)
	b_qsort_matrix_columns_descending_real(x,m,(k-idx)/m+1,idx);
      if (j < p)
	b_qsort_matrix_columns_descending_real(x+j-idx,m,n-(j-idx)/m,idx);
    }
  else
    {
      // insertion sort
      p=n*m;
      for (j=m+idx; j < p; j+=m)
	if (x[j-m] < x[j])
	  {
	    memcpy(y,x+j-idx,s);
	    memcpy(x+j-idx,x+j-m-idx,s);
	    for (k=j-m-m; k >= 0; k-=m)
	      if (x[k] >= y[idx])
		break;
	      else
		memcpy(x+k+m-idx,x+k-idx,s);
	    memcpy(x+k+m-idx,y,s);
	  }
    }
  free(y);
}

/*
   Assumes:
     x   - array of length m*n in colum major format.
     m   - number of rows
     n   - number of columns
     br  - first row in block to sort
     er  - last row in block to sort to sort
     idx - idx/m is column to sort 

   Results:
     The rows of x are sorted in ascending order on column idx/m.

   Notes:
     Uses the quick sort mean algorithm.  Switches to insertion sort when the
     size of the list is 10 or less.  If the matrix is in row major format, then
     m is the number of columns, n is the number of rows, and the columns of x 
     are sorted in ascending order on row idx. 
*/
static void b_qsort_matrix_rows_ascending_real(PRECISION *x, int m, int n, int br, int er, int idx)
{
  PRECISION y, c;
  int i, j, k;
  if (er-br+1 > 10)
    {
      // quick sort
      if (x[idx+br] == x[idx+er])
	c=x[idx+br];
      else
	{
	  if (x[idx+br] > x[idx+er])
	    for (i=(n-1)*m; i >= 0; i-=m)
	      { y=x[i+br]; x[i+br]=x[i+er]; x[i+er]=y; }
	  c=0.5*(x[idx+br] + x[idx+er]);
	}

      for (j=br+1; (j < er) && (x[idx+j] <= c); j++);
      for (k=er-1; (k > br) && (x[idx+k] >= c); k--);
      while (j < k)
	{
	  for (i=(n-1)*m; i >= 0; i-=m)
	    { y=x[i+j]; x[i+j]=x[i+k]; x[i+k]=y; }
	  while (x[idx+j] <= c) j++;
	  while (x[idx+k] >= c) k--;
	}
      if (k > br)
	b_qsort_matrix_rows_ascending_real(x,m,n,br,k,idx);
      if (j < er)
	b_qsort_matrix_rows_ascending_real(x,m,n,j,er,idx);
    }
  else
    {
      // insertion sort
      int r;
      for (j=br+1; j <= er; j++)
	{
	  for (k=j-1; k >= br; k--)
	    if (x[idx+k] <= x[idx+j]) break;
          if (++k < j)
	    for (i=(n-1)*m; i >= 0; i-=m)
	      {
		y=x[i+j];
		for (r=j; r > k; r--) x[i+r]=x[i+r-1];
		x[i+k]=y;
	      }
	}
    }
}

/*
   Assumes:
     x   - array of length m*n in colum major format.
     m   - number of rows
     n   - number of columns
     br  - first row in block to sort
     er  - last row in block to sort to sort
     idx - idx/m is column to sort 

   Results:
     The rows of x are sorted in ascending order on column idx/m.

   Notes:
     Uses the quick sort mean algorithm.  Switches to insertion sort when the
     size of the list is 10 or less.  If the matrix is in row major format, then
     m is the number of columns, n is the number of rows, and the columns of x 
     are sorted in ascending order on row idx. 
*/
static void b_qsort_matrix_rows_descending_real(PRECISION *x, int m, int n, int br, int er, int idx)
{
  PRECISION y, c;
  int i, j, k;
  if (er-br+1 > 10)
    {
      // quick sort
      if (x[idx+br] == x[idx+er])
	c=x[idx+br];
      else
	{
	  if (x[idx+br] < x[idx+er])
	    for (i=(n-1)*m; i >= 0; i-=m)
	      { y=x[i+br]; x[i+br]=x[i+er]; x[i+er]=y; }
	  c=0.5*(x[idx+br] + x[idx+er]);
	}

      for (j=br+1; (j < er) && (x[idx+j] >= c); j++);
      for (k=er-1; (k > br) && (x[idx+k] <= c); k--);
      while (j < k)
	{
	  for (i=(n-1)*m; i >= 0; i-=m)
	    { y=x[i+j]; x[i+j]=x[i+k]; x[i+k]=y; }
	  while (x[idx+j] >= c) j++;
	  while (x[idx+k] <= c) k--;
	}
      if (k > br)
	b_qsort_matrix_rows_descending_real(x,m,n,br,k,idx);
      if (j < er)
	b_qsort_matrix_rows_descending_real(x,m,n,j,er,idx);
    }
  else
    {
      // insertion sort
      int r;
      for (j=br+1; j <= er; j++)
	{
	  for (k=j-1; k >= br; k--)
	    if (x[idx+k] >= x[idx+j]) break;
          if (++k < j)
	    for (i=(n-1)*m; i >= 0; i-=m)
	      {
		y=x[i+j];
		for (r=j; r > k; r--) x[i+r]=x[i+r-1];
		x[i+k]=y;
	      }
	}
    }
}

/*
   Assumes:
     x - array of length m

   Results:
     x is sorted in ascending order

   Notes:
     Uses the quick sort median of three algorithm
*/
static void b_median_qsort_array_ascending(PRECISION *x, int m)
{
  PRECISION y;
  int j, k;
  if (m > 10)
    {
      // Quick sort
      j=(m--)/2;

      y=x[j]; x[j]=x[1]; x[1]=y;

      if (x[1] > x[m])
	if (x[0] > x[m])
	  if (x[0] > x[1])
	    { y=x[0]; x[0]=x[m]; x[m]=y; }
	  else
	    { y=x[0]; x[0]=x[m]; x[m]=x[1]; x[1]=y; }
	else
	  { y=x[1]; x[1]=x[m]; x[m]=y; }
      else
	if (x[0] > x[1])
	  if (x[0] > x[m])
	    { y=x[0]; x[0]=x[1]; x[1]=x[m]; x[m]=y; }
	  else
	    { y=x[0]; x[0]=x[1]; x[1]=y; };


      for (j=2; (j < m) && (x[j] <= x[1]); j++);
      for (k=m-1; (k > 1) && (x[k] >= x[1]); k--);
      while (j < k)
	{
	  y=x[j]; x[j]=x[k]; x[k]=y;
	  while (x[j] <= x[1]) j++;
	  while (x[k] >= x[1]) k--;
	}
      if (k > 1)
	{
	  y=x[k]; x[k]=x[1]; x[1]=y;
	  b_median_qsort_array_ascending(x,k);
	}
      if (j < m)
	b_median_qsort_array_ascending(x+j,m-j+1);
    }
  else
    {
      // Insertion sort
      for (j=1; j < m; j++)
	{
	  for (k=j-1; k >= 0; k--)
	    if (x[j] >= x[k]) break;
	  if (++k < j)
	    {
	      y=x[j];
	      memmove(x+k+1,x+k,(j-k)*sizeof(PRECISION));
	      x[k]=y;
	    }
	}
    }
}
