/** Example:
   #if defined (USE_DEBUG_FILE)
   fprintf(FPTR_DEBUG, "\nWARNING: .../mathlib.c/TransposeSquare():  the matrix is already both SU and SL, so there is no need to transpose.\n");
   fflush(FPTR_DEBUG);
   #else
   printf("\nWARNING: .../mathlib.c/TransposeSquare():  the matrix is already both SU and SL, so there is no need to transpose.\n");
   fflush(stdout);
   #endif
/**/


#include "tzmatlab.h"

#include "modify_for_mex.h"

FILE *FPTR_DEBUG = (FILE *)NULL;  //Debug output file, to be opened by main.c.
FILE *FPTR_OPT = (FILE *)NULL;  //Optimization output file, to be opened by main.c.

//-----------------
// Some high-level functions.
//-----------------
int fn_locofyearqm(int q_m, int yrstart, int qmstart, int yrend, int qmend)
{
   //Returns the (base 0) location of a specified year and month (quarter) for the time series.
   //All the other inputs take the usual (base-1) numbers, I guess 01/17/05.  For example, yrstart = 1960 means the year 1960.
   int tmpi, loc;

   if ( q_m != 12 )
      if ( q_m != 4 )  fn_DisplayError(".../tzmatlab.c/fn_locofyearqm(): This function only works for monthly or quarterly data");

   if ( (tmpi=yrend - yrstart) < 0 )
      fn_DisplayError(".../cstz.c/fn_locofyearqm(): the end year must be greater than or equal to the start year");
   else if ( (loc = (tmpi==0) ? qmend-qmstart : tmpi*q_m+qmend-qmstart) < 0 )
      fn_DisplayError(".../tzmatlab.c/fn_locofyearqm(): the end month or quarter must be greater than or equal to the start month or quarter for the given year");

   return (loc);
}



//-----------------
// Function to display errors.
//-----------------
void fn_DisplayError(char *msg_s)
{
   #if defined (USE_DEBUG_FILE)
      fprintf(FPTR_DEBUG, "\nFatal Error:\n"
             "  %s!\n", msg_s);
      fflush(FPTR_DEBUG);
   #else
      printf("\nFatal Error:\n"
             "\t %s!\n", msg_s);
      fflush(stdout);
   #endif

   #ifdef WIN_MATLABAPI
      mexErrMsgTxt(".");
   #else
      //getchar();
      exit( EXIT_FAILURE );                          // This exits the entire C program.
   #endif
}


//-----------------
// Error-checking memory allocators
//-----------------
void *m_alloc(size_t size) {
   void *new_mem;
   if ( (new_mem = malloc(size)) == NULL )  fn_DisplayError("Out of Memory!");
   return(new_mem);
}
//+
void *c_alloc(size_t elt_count, size_t elt_size) {
   void *new_mem;
   if ( (new_mem = calloc(elt_count, elt_size)) == NULL )  fn_DisplayError("Out of Memory!");
   return(new_mem);
}


//-----------------
// Creat and destroy vectors, matrices, and cells.
//-----------------
/**
TSvoidvector *CreateVector_void(int _n)
{
   TSvoidvector *x_voidv = tzMalloc(1, TSvoidvector);
   x_voidv->flag = V_UNDEF;
   x_voidv->n = _n;
   x_voidv->v = tzMalloc(_n, void);
   return(x_voidv);
}
TSvoidvector *DestroyVector_void(TSvoidvector *x_voidv)
{
   if (x_voidv) {
      free(x_voidv->v);
      free(x_voidv);
      return ((TSvoidvector *)NULL);
   }
   else  return (x_voidv);
}
/**/


TScvector *CreateVector_c(int _n)
{
   TScvector *x_cv = tzMalloc(1, TScvector);
   x_cv->flag = V_UNDEF;
   x_cv->n = _n;
   if (_n<1)  fn_DisplayError(".../tzmatlab.c/CreateVector_c(): dimension input _n must be a positive integer");
   x_cv->v = tzMalloc(_n, char);
   return( x_cv );
}
TScvector *DestroyVector_c(TScvector *x_cv)
{
   if (x_cv) {
      free(x_cv->v);
      free(x_cv);
      return ((TScvector *)NULL);
   }
   else  return (x_cv);
}

TSivector *CreateVector_int(int _n)
{
   TSivector *x_iv=tzMalloc(1, TSivector);
   x_iv->flag = V_UNDEF;
   x_iv->n = _n;
   if (_n<1)  fn_DisplayError(".../tzmatlab.c/CreateVector_int(): dimension input _n must be a positive integer");
   x_iv->v = tzMalloc(_n, int);
   return(x_iv);
}
TSivector *DestroyVector_int(TSivector *x_iv)
{
   if (x_iv) {
      free(x_iv->v);
      free(x_iv);
      return ((TSivector *)NULL);
   }
   else  return (x_iv);
}

TSimatrix *CreateMatrix_int(int nrows, int ncols)
{
   TSimatrix *x_im=tzMalloc(1, TSimatrix);
   x_im->nrows = nrows;
   x_im->ncols = ncols;
   if (nrows<1 || ncols<1)  fn_DisplayError(".../tzmatlab.c/CreateMatrix_int(): dimension inputs nrows and ncols must both be positive integers");
   x_im->M = tzMalloc(nrows*ncols, int);
   return (x_im);
}
TSimatrix *DestroyMatrix_int(TSimatrix *x_im)
{
   if (x_im) {
      free(x_im->M);
      free(x_im);
		return ((TSimatrix *)NULL);
   }
	else  return  (x_im);
}

TSicellvec *CreateCellvec_int(TSivector *n_iv)
{
   int _i,
       ncells;
	TSicellvec *x_icv = tzMalloc(1, TSicellvec);

	if (!n_iv || !n_iv->flag)  fn_DisplayError(".../CreateCellvec_int( ): Dimension vector n_iv must (1) created and (2) assigned legal values");
   x_icv->ncells = ncells = n_iv->n;
   x_icv->C = tzMalloc(ncells, TSivector *);
   for (_i=ncells-1; _i>-0; _i--)   *(x_icv->C + _i) = CreateVector_int(n_iv->v[_i]);
   return(x_icv);
}
TSicellvec *DestroyCellvec_int(TSicellvec *x_icv)
{
   int _i;
   if (x_icv) {
      for (_i=0; _i<x_icv->ncells; _i++)  DestroyVector_int(x_icv->C[_i]);
      free(x_icv->C);
      free(x_icv);
      return ((TSicellvec *)NULL);
   }
   else  return (x_icv);
}

TSicell *CreateCell_int(TSivector *row_iv, TSivector *col_iv)
{
   int _i,
       ncells;
   TSicell *x_ic=NULL;
   if (!row_iv || !col_iv || !row_iv->flag || !col_iv->flag)  fn_DisplayError(".../CreateCell_int( ): Dimension vectors row_iv and col_iv must (1) created and (2) assigned legal values");
   if ((ncells = row_iv->n) != col_iv->n)  fn_DisplayError(".../CreateCell_int( ): the lengths of row_iv and col_iv (i.e., numbers of cells) must be the same");
   x_ic = tzMalloc(1, TSicell);
   x_ic->ncells = ncells;
   x_ic->C = tzMalloc(ncells, TSimatrix *);
   for (_i=ncells-1; _i>=0; _i--) {
      *(x_ic->C + _i) = CreateMatrix_int(row_iv->v[_i], col_iv->v[_i]);
   }
   return(x_ic);
}
TSicell *DestroyCell_int(TSicell *x_ic)
{
   int _i;
   if (x_ic) {
      for (_i=x_ic->ncells-1; _i>=0; _i--)   x_ic->C[_i] = DestroyMatrix_int(x_ic->C[_i]);
      tzDestroy(x_ic->C);
      free(x_ic);
      return ((TSicell *)NULL);
   }
   else  return (x_ic);
}




TSdvector *CreateVector_lf(int _n)
{
   TSdvector *x_dv=tzMalloc(1, TSdvector);
   x_dv->flag = V_UNDEF;
   x_dv->n = _n;
   if (_n<1)  fn_DisplayError(".../tzmatlab.c/CreateVector_lf(): dimension input _n must be a positive integers");
   x_dv->v = tzMalloc(_n, double);
   return(x_dv);
}
TSdvector *DestroyVector_lf(TSdvector *x_dv)
{
   if (x_dv) {
      free(x_dv->v);
      free(x_dv);
      return ((TSdvector *)NULL);
   }
   else  return (x_dv);
}

TSdmatrix *CreateMatrix_lf(int nrows, int ncols)
{
   TSdmatrix *x_dm=tzMalloc(1, TSdmatrix);
   x_dm->flag = M_UNDEF;
   x_dm->nrows = nrows;
   x_dm->ncols = ncols;
      if (nrows<1 || ncols<1)  fn_DisplayError(".../tzmatlab.c/CreateMatrix_lf(): dimension inputs nrows and ncols must both be positive integers");
   x_dm->M = tzMalloc(nrows*ncols, double);
   return(x_dm);
}
TSdmatrix *DestroyMatrix_lf(TSdmatrix *x_dm)
{
   if (x_dm) {
      free(x_dm->M);
      free(x_dm);
      return ((TSdmatrix *)NULL);
   }
   else  return (x_dm);
}

TSdcell *CreateCell_lf(TSivector *row_iv, TSivector *col_iv)
{
   int _i,
       ncells;
   TSdcell *x_dc=NULL;
   //-------------- The following line must be enacted when we produce new code in the future.   ---------------------
   //-------------- In old code I forgot to set the flags for row_iv and col_iv but change them in all places are too time-consuming at this point.  ---------------------
   //if (!row_iv || !col_iv || !row_iv->flag || !col_iv->flag)  fn_DisplayError(".../CreateCell_lf( ): Dimension vectors row_iv and col_iv must (1) created and (2) assigned legal values");
   if ((ncells = row_iv->n) != col_iv->n)  fn_DisplayError(".../CreateCell_lf( ): the lengths of row_iv and col_iv (i.e., numbers of cells) must be the same");
   x_dc = tzMalloc(1, TSdcell);
   x_dc->ncells = ncells;
   x_dc->C = tzMalloc(ncells, TSdmatrix *);
   for (_i=ncells-1; _i>=0; _i--) {
      *(x_dc->C + _i) = CreateMatrix_lf(row_iv->v[_i], col_iv->v[_i]);
   }
   return(x_dc);
}
TSdcell *DestroyCell_lf(TSdcell *x_dc)
{
   int _i;
   if (x_dc) {
      for (_i=x_dc->ncells-1; _i>=0; _i--)  x_dc->C[_i] = DestroyMatrix_lf(x_dc->C[_i]);
      tzDestroy(x_dc->C);
      free(x_dc);
      return ((TSdcell *)NULL);
   }
   else  return (x_dc);
}

TSdcellvec *CreateCellvec_lf(TSivector *n_iv) {
   TSdcellvec *x_dcv = tzMalloc(1, TSdcellvec);
   int _i,
       ncells;
   //-------------- The following line must be enacted when we produce new code in the future.   ---------------------
   //-------------- In old code I forgot to set the flag for n_iv but change it in all places are too time-consuming at this point.  ---------------------
   //if (!n_iv || !n_iv->flag)  fn_DisplayError(".../CreateCellvec_lf( ): Dimension vector n_iv must (1) created and (2) assigned legal values");
   x_dcv->ncells = ncells = n_iv->n;
   x_dcv->C = tzMalloc(ncells, TSdvector *);
   for (_i=0; _i<ncells; _i++)   *(x_dcv->C + _i) = CreateVector_lf(n_iv->v[_i]);
   return(x_dcv);
}
TSdcellvec *DestroyCellvec_lf(TSdcellvec *x_dcv) {
   int _i;
   if (x_dcv) {
      for (_i=x_dcv->ncells-1; _i>=0; _i--)  DestroyVector_lf(x_dcv->C[_i]);
      free(x_dcv->C);
      free(x_dcv);
      return ((TSdcellvec *)NULL);
   }
   else  return (x_dcv);
}

TSdfourth *CreateFourth_lf(int ndims, TSivector *row_iv, TSivector *col_iv) {
   int _i;
   TSdfourth *x_d4 = NULL;
   //if (row_iv->n != col_iv->n) fn_DisplayError(".../CreateFourth_lf( ): the lengths of row_iv and col_iv (i.e., sizes of dimensions) must be the same");

   x_d4 = tzMalloc(1, TSdfourth);
   x_d4->ndims = ndims;
   x_d4->F = tzMalloc(ndims, TSdcell *);
   for (_i=ndims-1; _i>=0; _i--) {
      *(x_d4->F + _i) = CreateCell_lf(row_iv, col_iv);
   }
   return(x_d4);
}
TSdfourth *DestroyFourth_lf(TSdfourth *x_d4) {
   int _i;
   if (x_d4) {
      for (_i=x_d4->ndims-1; _i>=0; _i--)  DestroyCell_lf(x_d4->F[_i]);
      free(x_d4->F);
      free(x_d4);
      return ((TSdfourth *)NULL);
   }
	else  return (x_d4);
}

TSdfourthvec *CreateFourthvec_lf(int ndims, TSivector *n_iv)
{
   int _i;
   TSdfourthvec *x_d4v = NULL;
   //if (n_iv->n != col_iv->n) fn_DisplayError(".../CreateFourth_lf( ): the lengths of n_iv and col_iv (i.e., sizes of dimensions) must be the same");

   x_d4v = tzMalloc(1, TSdfourthvec);
   x_d4v->ndims = ndims;
   x_d4v->F = tzMalloc(ndims, TSdcellvec *);
   for (_i=ndims-1; _i>=0; _i--) {
      *(x_d4v->F + _i) = CreateCellvec_lf(n_iv);
   }
   return(x_d4v);
}
TSdfourthvec *DestroyFourthvec_lf(TSdfourthvec *x_d4v)
{
   int _i;
   if (x_d4v) {
      for (_i=x_d4v->ndims-1; _i>=0; _i--)  DestroyCellvec_lf(x_d4v->F[_i]);
      free(x_d4v->F);
      free(x_d4v);
      return ((TSdfourthvec *)NULL);
   }
	else  return (x_d4v);
}

TSdzvector *CreateVector_dz(int _n)
{
   TSdzvector *x_dzv=tzMalloc(1, TSdzvector);
   x_dzv->real = CreateVector_lf(_n);
   x_dzv->imag = CreateVector_lf(_n);
   return( x_dzv );
}
TSdzvector *DestroyVector_dz(TSdzvector *x_dzv)
{
   if (x_dzv) {
      DestroyVector_lf(x_dzv->real);
      DestroyVector_lf(x_dzv->imag);
      free(x_dzv);
      return ((TSdzvector *)NULL);
   }
   else  return (x_dzv);
}

TSdzmatrix *CreateMatrix_dz(int nrows, int ncols) {
   TSdzmatrix *x_dzm=tzMalloc(1, TSdzmatrix);
   x_dzm->real = CreateMatrix_lf(nrows, ncols);
   x_dzm->imag = CreateMatrix_lf(nrows, ncols);
   return( x_dzm );
}
TSdzmatrix *DestroyMatrix_dz(TSdzmatrix *x_dzm)
{
   if (x_dzm) {
      DestroyMatrix_lf(x_dzm->real);
      DestroyMatrix_lf(x_dzm->imag);
      free(x_dzm);
      return ((TSdzmatrix *)NULL);
   }
   else  return (x_dzm);
}



//-----------------
// Creates special vectors, matrices, and cells but uses the same destroy utilities as above.
//-----------------
//=== Creates two special matrices: zeros and identity.  Use DestroyMatrix_lf to free the memory allocated to these functions.
TSdmatrix *CreateZeroMatrix_lf(const int nrows, const int ncols) {
   int _i;
   TSdmatrix *x_dm=CreateMatrix_lf(nrows, ncols);
   //x_dm->flag = M_GE | M_SU | M_SL | M_UT | M_LT;
   x_dm->flag = M_GE;
   for (_i=nrows*ncols-1; _i>=0; _i--)
      x_dm->M[_i] = 0.0;
   return(x_dm);
}
TSdmatrix *CreateIdentityMatrix_lf(const int nrows, const int ncols) {
   int _i;
   TSdmatrix *x_dm=CreateZeroMatrix_lf(nrows, ncols);
   if (nrows==ncols) {
      //x_dm->flag = M_GE | M_SU | M_SL | M_UT | M_LT;
      //x_dm->flag = M_GE;
      for (_i=square(nrows)-1; _i>=0; _i -= nrows+1)  x_dm->M[_i] = 1.0;
      x_dm->flag = M_GE | M_SU | M_SL | M_UT | M_LT;
   }
   else if (nrows<ncols) {
      //x_dm->flag = M_GE | M_SU | M_UT;
      //x_dm->flag = M_GE;
      for (_i=square(nrows)-1; _i>=0; _i -= nrows+1)  x_dm->M[_i] = 1.0;
      x_dm->flag = M_GE | M_UT | M_LT;
   }
   else {
      //x_dm->flag = M_GE | M_SL | M_LT;
      //x_dm->flag = M_GE;
      for (_i=(ncols-1)*(nrows+1); _i>=0; _i -= nrows+1)  x_dm->M[_i] = 1.0;
      x_dm->flag = M_GE | M_UT | M_LT;
   }
   return(x_dm);
}

//=== Other speicial matrices.
TSivector *CreateConstantVector_int(const int _n, const int _k) {
   //Inputs:
   //  _k:  Integer constant;
   //  _n: Dimension of the vector.
   int _i;
   TSivector *x_iv=CreateVector_int(_n);
   for (_i=_n-1; _i>=0; _i--)
      x_iv->v[_i] = _k;
   x_iv->flag = V_DEF;
   return(x_iv);
}

TSimatrix *CreateConstantMatrix_int(const int nrows, const int ncols, const int _n)
{
   int _i;
   TSimatrix *x_im=CreateMatrix_int(nrows, ncols);

   for (_i=nrows*ncols-1; _i>=0; _i--)  x_im->M[_i] = _n;
   if ( nrows==ncols )   x_im->flag = M_GE | M_SU | M_SL | M_CN;
   else  x_im->flag = M_GE | M_CN;
   return(x_im);
}

TSicellvec *CreateConstantCellvec_int(TSivector *n_iv, const int _n)
{
   int _i,
       ncells;
	TSicellvec *x_icv = tzMalloc(1, TSicellvec);

	if (!n_iv || !n_iv->flag)  fn_DisplayError(".../CreateCellvec_int( ): Dimension vector n_iv must (1) created and (2) assigned legal values");
   x_icv->ncells = ncells = n_iv->n;
   x_icv->C = tzMalloc(ncells, TSivector *);
   for (_i=ncells-1; _i>=0; _i--)   *(x_icv->C + _i) = CreateConstantVector_int(n_iv->v[_i], _n);
   return(x_icv);
}

TSicell *CreateConstantCell_int(TSivector *row_iv, TSivector *col_iv, const int _n)
{
   int _i,
       ncells;
   TSicell *x_ic=NULL;
   if (!row_iv || !col_iv || !row_iv->flag || !col_iv->flag)  fn_DisplayError(".../CreateConstantCell_int( ): Dimension vectors row_iv and col_iv must (1) created and (2) assigned legal values");
   if ((ncells = row_iv->n) != col_iv->n)  fn_DisplayError(".../CreateCell_int( ): the lengths of row_iv and col_iv (i.e., numbers of cells) must be the same");

   x_ic = tzMalloc(1, TSicell);
   x_ic->ncells = ncells;
   x_ic->C = tzMalloc(ncells, TSimatrix *);
   for (_i=ncells-1; _i>=0; _i--)   *(x_ic->C + _i) = CreateConstantMatrix_int(row_iv->v[_i], col_iv->v[_i], _n);
   return(x_ic);
}


TSdvector *CreateConstantVector_lf(const int _n, const double _alpha) {
   int _i;
   TSdvector *x_dv=CreateVector_lf(_n);
   for (_i=_n-1; _i>=0; _i--)   x_dv->v[_i] = _alpha;
   x_dv->flag = V_DEF;
   return(x_dv);
}

TSdmatrix *CreateConstantMatrix_lf(const int nrows, const int ncols, const double _alpha) {
   //Inputs:
   //  _alpha:  Double constant;
   //  nrows and ncols: Dimensions of the matrix.
   int _i;
   TSdmatrix *x_dm=CreateMatrix_lf(nrows, ncols);

   for (_i=nrows*ncols-1; _i>=0; _i--)  x_dm->M[_i] = _alpha;
   if ( nrows==ncols )   x_dm->flag = M_GE | M_SU | M_SL | M_CN;
   else  x_dm->flag = M_GE | M_CN;
   return(x_dm);
}

TSdcellvec *CreateConstantCellvec_lf(TSivector *n_iv, const double _alpha) {
   //Inputs:
   //  _alpha:  Double constant;
   //  _n: Length (dimension) of the vector.
   int _i,
       ncells;
   TSdcellvec *x_dcv = tzMalloc(1, TSdcellvec);
   //-------------- The following line must be enacted when we produce new code in the future.   ---------------------
   //-------------- In old code I forgot to set the flag for n_iv but change it in all places are too time-consuming at this point.  ---------------------
   //if (!n_iv || !n_iv->flag)  fn_DisplayError(".../CreateConstantCellvec_lf( ): Dimension vector n_iv must (1) created and (2) assigned legal values");
   x_dcv->ncells = ncells = n_iv->n;
   x_dcv->C = tzMalloc(ncells, TSdvector *);
   for (_i=ncells-1; _i>=0; _i--)   *(x_dcv->C + _i) = CreateConstantVector_lf(n_iv->v[_i], _alpha);
   return(x_dcv);
}

TSdcell *CreateConstantCell_lf(TSivector *row_iv, TSivector *col_iv, const double _alpha) {
   //Inputs:
   //  _alpha:  Double constant;
   //  nrows: Number of rows;
   //  ncols: Number of columns.
   int _i,
       ncells;
   TSdcell *x_dc=NULL;
   //-------------- The following line must be enacted when we produce new code in the future.   ---------------------
   //-------------- In old code I forgot to set the flags for row_iv and col_iv but change them in all places are too time-consuming at this point.  ---------------------
   //if (!row_iv || !col_iv || !row_iv->flag || !col_iv->flag)  fn_DisplayError(".../CreateConstantCell_lf( ): Dimension vectors row_iv and col_iv must (1) created and (2) assigned legal values");
   if ((ncells = row_iv->n) != col_iv->n)  fn_DisplayError(".../CreateCell_lf( ): the lengths of row_iv and col_iv (i.e., numbers of cells) must be the same");

   x_dc = tzMalloc(1, TSdcell);
   x_dc->ncells = ncells;
   x_dc->C = tzMalloc(ncells, TSdmatrix *);
   for (_i=ncells-1; _i>=0; _i--)   *(x_dc->C + _i) = CreateConstantMatrix_lf(row_iv->v[_i], col_iv->v[_i], _alpha);
   return(x_dc);
}


TSdvector *CreateDatesVector_lf(int nq_m, int yrstart, int qmstart, int yrend, int qmend)
{
   //If nq_m==4, quarterly data; nq_m==12, monthly data.
   //All the other inputs take the usual (base-1) numbers, I guess 01/17/05.  For example, yrstart = 1960 means the year 1960.
   int _t;
   int samplesize = 1+fn_locofyearqm(nq_m, yrstart, qmstart, yrend, qmend);  //1+ because fn_locofyearqm() returns a 0-based integer.
   //
   TSdvector *dates_dv = tzMalloc(1, TSdvector);
   dates_dv->n = samplesize;
   dates_dv->v = tzMalloc(samplesize, double);

   if (nq_m==4 || nq_m==12) {
      for (_t=samplesize-1; _t>=0; _t--)  dates_dv->v[_t] = (double)yrstart + (double)(qmstart+_t-1)/(double)nq_m;
      dates_dv->flag = V_DEF;
   }
   else  fn_DisplayError(".../tzmatlab.c/CreateDatesVector_lf(): Dates have to be either monthly or quarterly");


   return (dates_dv);
}



//-----------------
// Initializes already-created special vectors, matrices, and cells.
//-----------------
void InitializeConstantVector_lf(TSdvector *x_dv, const double _alpha)
{
   //Ouputs:
   //  x_dv: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_dv:  Memory allocated already.
   //  _alpha:  Double constant;
   int _i, _n;

   if (!x_dv)  fn_DisplayError(".../tzmatlab.c/InitializeConstantVector_lf():  Input vector must be created (memory-allocated)");
   else {
      _n=x_dv->n;
   }
   for (_i=_n-1; _i>=0; _i--)   x_dv->v[_i] = _alpha;
   x_dv->flag = V_DEF;
}

void InitializeConstantVector_int(TSivector *x_iv, const int _k)
{
   //Ouputs:
   //  x_iv: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_iv:  Memory allocated already.
   //  _alpha:  Integer constant;
   int _i, _n;

   if (!x_iv)  fn_DisplayError(".../tzmatlab.c/InitializeConstantVector_int():  Input vector must be created (memory-allocated)");
   else {
      _n=x_iv->n;
   }
   for (_i=_n-1; _i>=0; _i--)  x_iv->v[_i] = _k;
   x_iv->flag = V_DEF;
}

void InitializeConstantMatrix_lf(TSdmatrix *x_dm, const double _alpha)
{
   //Ouputs:
   //  x_dm: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_dm:  Memory allocated already.
   //  _alpha:  Double constant;
   //See Kenneth Reek, pp.202-212.

//   int _i;
//   for (_i=x_dm->nrows*x_dm->ncols-1; _i>=0; _i--)
//      x_dm->M[_i] = _alpha;
//   int nrows, ncols;
   double *ptrcnt, *lastptr;

   if ( !x_dm) fn_DisplayError(".../tzmathlab.c/InitializeConstantMatrix_int(): Input matrix must be created (memory-allocated)");
   else {
//      nrows = x_dm->nrows;
//      ncols = x_dm->ncols;
      lastptr = (ptrcnt = x_dm->M) + x_dm->nrows * x_dm->ncols;
   }

//   if (nrows==ncols)  x_dm->flag = M_GE | M_SU | M_SL;
//   else if (nrows<ncols)  x_dm->flag = M_GE | M_SU;
//   else  x_dm->flag = M_GE | M_SL;
   x_dm->flag = M_GE | M_CN;
   for ( ; ptrcnt<lastptr; ptrcnt++ )  *ptrcnt = _alpha;
}

void InitializeDiagonalMatrix_lf(TSdmatrix *x_dm, const double _alpha) {
   int _i, n2, nrows, ncols;
   double *M;

   if ( !x_dm )  fn_DisplayError(".../tzmathlab.c/InitializeIdentiyMatrix_lf(): (1) Input matrix must be created (memory-allocated)");
   else {
      nrows = x_dm->nrows;
      ncols = x_dm->ncols;
      M = x_dm->M;
   }

   if (nrows==ncols) {
      for (_i=(n2=square(nrows))-1; _i>=0; _i--)  M[_i] = 0.0;
      for (_i=n2-1; _i>=0; _i -= nrows+1)  M[_i] = _alpha;
      x_dm->flag = M_GE | M_SU | M_SL | M_UT | M_LT;
   }
   else if (nrows<ncols) {
      for (_i=nrows*ncols-1; _i>=0; _i--)  M[_i] = 0.0;
      for (_i=square(nrows)-1; _i>=0; _i -= nrows+1)  M[_i] = _alpha;
      x_dm->flag = M_GE | M_UT | M_LT;
   }
   else {
      for (_i=nrows*ncols-1; _i>=0; _i--)  M[_i] = 0.0;
      for (_i=(ncols-1)*(nrows+1); _i>=0; _i -= nrows+1)  M[_i] = _alpha;
      x_dm->flag = M_GE | M_UT | M_LT;
   }
}

void InitializeConstantMatrix_int(TSimatrix *x_im, const int _alpha) {
   //Ouputs:
   //  x_im: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_im:  Memory allocated already.
   //  _alpha:  Integer constant;
   //
   //See Kenneth Reek, pp.202-212.


//   int _i;
//   for (_i=x_im->nrows*x_im->ncols-1; _i>=0; _i--)
//      x_im->M[_i] = _alpha;

   int *ptrcnt, *lastptr;

   if ( !x_im) fn_DisplayError(".../tzmathlab.c/InitializeConstantMatrix_int(): Input matrix must be created (memory-allocated)");
   else lastptr = (ptrcnt = x_im->M) + x_im->nrows * x_im->ncols;

   for ( ; ptrcnt<lastptr; ptrcnt++ ) *ptrcnt = _alpha;
}

void InitializeConstantCellvec_lf(TSdcellvec *x_dcv, const double _alpha) {
   //Ouputs:
   //  x_dcv: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_dcv:  Memory allocated already.
   //  _alpha:  Double constant;
   int _i, _k, _n;
   double *v;

   if ( !x_dcv ) fn_DisplayError(".../tzmatlab.c/InitializeConstantCellvec_lf(): Input cell vector must be created (memory-allocated)");


   for (_i=x_dcv->ncells-1; _i>=0; _i--) {
      v = x_dcv->C[_i]->v;
      _n = x_dcv->C[_i]->n;
      for (_k=_n-1; _k>=0; _k--) v[_k] = _alpha;
      x_dcv->C[_i]->flag = V_DEF;
   }
}

void InitializeConstantCell_lf(TSdcell *x_dc, const double _alpha)
{
   //Ouputs:
   //  x_dc: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_dc:  Memory allocated already.
   //  _alpha:  Double constant;
   int _i, _k, nrows, ncols;
   double *M;

   if ( !x_dc ) fn_DisplayError(".../tzmatlab.c/InitializeConstantCell_lf(): Input cell must be created (memory-allocated)");


   for (_i=x_dc->ncells-1; _i>=0; _i--) {
      M = x_dc->C[_i]->M;
      nrows = x_dc->C[_i]->nrows;
      ncols = x_dc->C[_i]->ncols;
//      if (nrows==ncols)  x_dc->C[_i]->flag = M_GE | M_SU | M_SL;
//      else if (nrows<ncols)  x_dc->C[_i]->flag = M_GE | M_SU;
//      else  x_dc->C[_i]->flag = M_GE | M_SL;
      for (_k=nrows*ncols-1; _k>=0; _k--) M[_k] = _alpha;
      x_dc->C[_i]->flag = M_GE | M_CN;
   }
}



void InitializeConstantFourthvec_lf(TSdfourthvec *x_d4v, const double _alpha) {
   //Ouputs:
   //  x_d4v: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_d4v:  Memory allocated already.
   //  _alpha:  Double constant;
   int _j, _i, _k;
   double *v;

   if ( !x_d4v ) fn_DisplayError(".../tzmatlab.c/InitializeConstantFourthvec_lf(): Input fourth must be created (memory-allocated)");

   for (_j=x_d4v->ndims-1; _j>=0; _j--) {
      for (_i=x_d4v->F[_j]->ncells-1; _i>=0; _i--) {
         v = x_d4v->F[_j]->C[_i]->v;
         for (_k=x_d4v->F[_j]->C[_i]->n-1; _k>=0; _k--)  v[_k] = _alpha;
         x_d4v->F[_j]->C[_i]->flag = V_DEF;
      }
   }
}
void InitializeConstantFourth_lf(TSdfourth *x_d4, const double _alpha) {
   //Ouputs:
   //  x_d4: Initialized to a constant value _alpha for all elements.
   //Inputs:
   //  x_d4:  Memory allocated already.
   //  _alpha:  Double constant;
   int _j, _i, _k, nrows, ncols;
   double *M;

   if ( !x_d4 ) fn_DisplayError(".../tzmatlab.c/InitializeConstantFourth_lf(): Input fourth must be created (memory-allocated)");

   for (_j=x_d4->ndims-1; _j>=0; _j--) {
      for (_i=x_d4->F[_j]->ncells-1; _i>=0; _i--) {
         M = x_d4->F[_j]->C[_i]->M;
         nrows = x_d4->F[_j]->C[_i]->nrows;
         ncols = x_d4->F[_j]->C[_i]->ncols;
         for (_k=nrows*ncols-1; _k>=0; _k--)  M[_k] = _alpha;
         x_d4->F[_j]->C[_i]->flag = M_GE | M_CN;
      }
   }
}


void NegateColofMatrix_lf(TSdvector *y_dv, TSdmatrix *X_dm, int jx) {
   //Ouputs:
   //  If y_dv!=NULL, y_dv is the negative of the jx_th column of X_dm (i.e., multiplied by -1.0).
   //  If !y_dv, the jx_th column of X_dm is replaced by its negated value (i.e., multiplied by -1.0).
   //Inputs:
   //  X_dm:  Memory allocated and legal values given already.
   //  jx: The jx_th column of X_dm.

   int _i, nrows_x;
   double *M, *v;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../tzmathlab.c/NegateColumnofMatrix_lf(): (1) input matrix must be created (memory-allocated); (2) legal values must be given");
   if (jx >= X_dm->ncols)  fn_DisplayError(".../tzmathlab.c/NegateColumnofMatrix_lf(): The jx_th column specified exceeds the column dimension of the input matrix");

   M = X_dm->M + (jx+1)*(nrows_x=X_dm->nrows) - 1;  //Points to the end of the jx_th column.
   if ( !y_dv )
      for (_i=nrows_x-1; _i>=0; _i--, M--)  *M = -(*M);
   else {
      for (_i=nrows_x-1, v=y_dv->v+_i; _i>=0; _i--, M--, v--)  *v = -(*M);
      y_dv->flag = V_DEF;
   }
}


void InitializeConstantColofMatrix_lf(TSdmatrix *X_dm, int jx, double _alpha) {
   //Ouputs:
   //  The jx_th column of X_dm is replaced by its original value multiplied by _alpha.
   //Inputs:
   //  X_dm:  Memory allocated and legal values given already.
   //  jx: The jx_th column of X_dm.
   //  _alpha: A double constant.

   int _i, nrows_x;
   double *M;

   if ( !X_dm || !X_dm->flag )  fn_DisplayError(".../tzmathlab.c/NegateColumnofMatrix_lf(): (1) input matrix must be created (memory-allocated); (2) legal values must be given");
   if (jx >= X_dm->ncols)  fn_DisplayError(".../tzmathlab.c/NegateColumnofMatrix_lf(): The jx_th column specified exceeds the column dimension of the input matrix");

   M = X_dm->M + (jx+1)*(nrows_x=X_dm->nrows) - 1;  //Points to the end of the jx_th column.
   for (_i=nrows_x-1; _i>=0; _i--, M--)  *M = _alpha;
}




//-----------------
// Open files.
//-----------------
FILE *tzFopen(char *filename, char *mode) {
   FILE *fptr_dummy;

   if (filename)
   {
      if ( !(fptr_dummy = fopen(filename,mode)) ) {
         printf("\n\n...tzmatlab.c/tzFopen(): Fatal Error -- unable to write, read, or append the file %s!\n", filename);
         //getchar();
         exit(EXIT_FAILURE);
      }
   }
   else  fn_DisplayError(".../tzmatlab.c/tzFopen(): the input filename must exit");

   return (fptr_dummy);
}
