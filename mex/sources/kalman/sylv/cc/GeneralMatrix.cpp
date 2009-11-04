/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/GeneralMatrix.cpp,v 1.4 2004/11/24 20:41:59 kamenik Exp $ */

/* Tag $Name:  $ */


#include "SylvException.h"
#include "GeneralMatrix.h"

#include "cppblas.h"
#include "cpplapack.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <limits>

//vector<int>nullVec(0);
int GeneralMatrix::md_length = 32;

GeneralMatrix::GeneralMatrix(const GeneralMatrix& m)
: data(m.rows*m.cols), rows(m.rows), cols(m.cols), ld(m.rows)
  {
  copy(m);
  }

GeneralMatrix::GeneralMatrix(const ConstGeneralMatrix& m)
: data(m.rows*m.cols), rows(m.rows), cols(m.cols), ld(m.rows)
  {
  copy(m);
  }

GeneralMatrix::GeneralMatrix(const GeneralMatrix& m, const char* dummy)
: data(m.rows*m.cols), rows(m.cols), cols(m.rows), ld(m.cols)
  {
  for (int i = 0; i < m.rows; i++)
    for (int j = 0; j < m.cols; j++)
      get(j,i) = m.get(i,j);
  }

GeneralMatrix::GeneralMatrix(const ConstGeneralMatrix& m, const char* dummy)
: data(m.rows*m.cols), rows(m.cols), cols(m.rows), ld(m.cols)
  {
  for (int i = 0; i < m.rows; i++)
    for (int j = 0; j < m.cols; j++)
      get(j,i) = m.get(i,j);
  }


GeneralMatrix::GeneralMatrix(const GeneralMatrix& m, int i, int j, int nrows, int ncols)
: data(nrows*ncols), rows(nrows), cols(ncols), ld(nrows)
  {
  copy(m, i, j);
  }

GeneralMatrix::GeneralMatrix(GeneralMatrix& m, int i, int j, int nrows, int ncols)
: data(m.base()+m.ld*j+i, m.ld*(ncols-1)+nrows), rows(nrows), cols(ncols), ld(m.ld)
  {}

GeneralMatrix::GeneralMatrix(const GeneralMatrix& a, const GeneralMatrix& b)
: data(a.rows*b.cols), rows(a.rows), cols(b.cols), ld(a.rows)
  {
  gemm("N", a, "N", b, 1.0, 0.0);
  }

GeneralMatrix::GeneralMatrix(const GeneralMatrix& a, const GeneralMatrix& b, const char* dum)
: data(a.rows*b.rows), rows(a.rows), cols(b.rows), ld(a.rows)
  {
  gemm("N", a, "T", b, 1.0, 0.0);
  }

GeneralMatrix::GeneralMatrix(const GeneralMatrix& a, const char* dum, const GeneralMatrix& b)
: data(a.cols*b.cols), rows(a.cols), cols(b.cols), ld(a.cols)
  {
  gemm("T", a, "N", b, 1.0, 0.0);
  }

GeneralMatrix::GeneralMatrix(const GeneralMatrix& a, const char* dum1,
                             const GeneralMatrix& b, const char* dum2)
                             : data(a.cols*b.rows), rows(a.cols), cols(b.rows), ld(a.cols)
  {
  gemm("T", a, "T", b, 1.0, 0.0);
  }

/* generate new matrix b as subset or whole of matrix a but reordered by 
vrows and vcols as Matlab b=a(vrows,vcols) where vectors vrows and vcols start from 1.
It ignores non-positive elements passing zero length vector is equivalent 
to Matlab operator ":" = all elements of that dimension in its order */
GeneralMatrix::GeneralMatrix(const GeneralMatrix& a, const vector<int>&vrows, const vector<int>&vcols)
  {
  int nrows=0, ncols=0;
  if (vrows.size()==0 && vcols.size()==0)
    {
    *this=a;
    return;
    }
  else
    {
    if (vrows.size()==0)
      nrows=a.numRows();
    else
      {
      for (int i=0;i<vrows.size();++i)
        {
        if (vrows[i]>0) 
          nrows++;
        else 
          throw SYLV_MES_EXCEPTION("Non-positive indices in construction by vector.");
        }
      }
    if (nrows>a.numRows())
      throw SYLV_MES_EXCEPTION("Wrong dimensions for construction by vector.");

    if (vcols.size()==0)
      ncols=a.numCols();
    else
      {
      for (int i=0;i<vcols.size();++i)
        {
        if (vcols[i]>0) 
          ncols++;
        else 
          throw SYLV_MES_EXCEPTION("Non-positive indices in construction by vector.");
        }
      }
    if (ncols>a.numCols())
      throw SYLV_MES_EXCEPTION("Wrong dimensions for construction by vector.");

    data= *(new Vector(nrows*ncols));
    rows=nrows; 
    cols=ncols;
    if(nrows*ncols==0) return;

    for (int i=0;i<nrows;++i)
      {
      for (int j=0;j<nrows;++j)
        if (vrows.size()!=0 && vcols.size()!=0)
          {
          if (vrows[i]>0 && vcols[j]  >0)
            get(i,j)=a.get(vrows[i]-1, vcols[j]-1);
          }
        else if (vrows.size()!=0 && vcols.size()==0)
          {
          if (vrows[i]>0 )
            get(i,j)=a.get(vrows[i]-1, j);
          }
        else if (vrows.size()==0 && vcols.size()!=0)
          {
          if (vcols[j] >0)
            get(i,j)=a.get(i, vcols[j]-1);
          }
      }
    }
  }



GeneralMatrix::~GeneralMatrix()
  {
  }


/* Matlab element product: this = this .*m */
void 
GeneralMatrix::multElements(const GeneralMatrix& m)
  {
  if(cols!=m.numCols() || rows!=m.numRows())
     throw SYLV_MES_EXCEPTION("multiply Element porduct: matrices must be same dimension.");
  for (int i=0;i<cols;++i)
    for (int j=0;j<rows;++j)
      get(i,j)*=m.get(i,j);
  };


/* emulates Matlab repmat: new matrix = multv*multh*this */
GeneralMatrix& 
GeneralMatrix::repmat(int multv, int multh)
  {
  GeneralMatrix* repMat=(new GeneralMatrix ( multv*rows,  multh*cols));
  for (int i=0;i<multv;++i)
    for (int j=0;j<multh;++j)
      (*repMat).place(*this, multv*i,  multh*j);
    return *repMat;
  };


void GeneralMatrix::place(const ConstGeneralMatrix& m, int i, int j)
  {
  if (i + m.numRows() > numRows() ||
    j + m.numCols() > numCols())
    throw SYLV_MES_EXCEPTION("Bad submatrix placement, matrix dimensions exceeded.");
  
  GeneralMatrix tmpsub(*this, i, j, m.numRows(), m.numCols());
  tmpsub.copy(m);
  }

/* this = a*b */
void GeneralMatrix::mult(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b)
  {
  gemm("N", a, "N", b, 1.0, 0.0);
  }

/* this = this + scalar*a*b */
void GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b,
                               double mult)
  {
  gemm("N", a, "N", b, mult, 1.0);
  }


/* this = this + scalar*a*b' */
void GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, const ConstGeneralMatrix& b,
                               const char* dum, double mult)
  {
  gemm("N", a, "T", b, mult, 1.0);
  }

/* this = this + scalar*a'*b */
void GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, const char* dum,
                               const ConstGeneralMatrix& b, double mult)
  {
  gemm("T", a, "N", b, mult, 1.0);
  }

/* this = this + scalar*a'*b' */
void GeneralMatrix::multAndAdd(const ConstGeneralMatrix& a, const char* dum1,
                               const ConstGeneralMatrix& b, const char* dum2, double mult)
  {
  gemm("T", a, "T", b, mult, 1.0);
  }

void GeneralMatrix::addOuter(const ConstVector& a, double mult)
  {
  if (numRows() != numCols())
    throw SYLV_MES_EXCEPTION("Matrix is not square in GeneralMatrix::addOuter.");
  if (numRows() != a.length())
    throw SYLV_MES_EXCEPTION("Wrong length of a vector in GeneralMatrix::addOuter.");
  
  // since BLAS dsyr (symmetric rank 1 update) assumes symmetricity, we do this manually
  for (int i = 0; i < numRows(); i++)
    for (int j = i; j < numRows(); j++) {
      double s = mult*a[i]*a[j];
      get(i,j) = get(i,j) + s;
      if (i != j)
        get(j,i) = get(j,i) + s;
      }
  }


void GeneralMatrix::multRight(const ConstGeneralMatrix& m)
  {
  gemm_partial_right("N", m, 1.0, 0.0);
  }

void GeneralMatrix::multLeft(const ConstGeneralMatrix& m)
  {
  gemm_partial_left("N", m, 1.0, 0.0);
  }

void GeneralMatrix::multRightTrans(const ConstGeneralMatrix& m)
  {
  gemm_partial_right("T", m, 1.0, 0.0);
  }

void GeneralMatrix::multLeftTrans(const ConstGeneralMatrix& m)
  {
  gemm_partial_left("T", m, 1.0, 0.0);
  }


/* this = this * A^(-1)  */
void 
GeneralMatrix::multInvRight( GeneralMatrix&A)
  {
  // check or allocate tmp space for Transpose *this
  /**
  if (tmpGMp)
    {
    if (tmpGMp->numCols()!=rows || tmpGMp->numRows()!=cols)
      delete (tmpGMp);
    }
  if (!tmpGMp) 
  ********/
  tmpGMp= new GeneralMatrix(cols,rows); // allocate space only once if and when needed!

  // tmpGMp=(*this)' i.e. Transpose (*this)
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      tmpGMp->get(j,i) = get(i,j);
  // check A and this suiability
  const int mcols=A.numCols();
  if (A.numRows() != mcols) {
    throw SYLV_MES_EXCEPTION("The matrix is not square for inversion.");
    }
  if (cols != A.numRows()) {
    throw SYLV_MES_EXCEPTION("Wrong dimensions for matrix inverse mutliply.");
    }
  
  if (rows > 0) 
    {
    /* out =tmpGMp out = inv(A')*(*this)' = inv(A')*(*tmpGMp) */
    int* ipiv = new int[cols];
    int info;
    const int mld=A.getLD();
    LAPACK_dgetrf(&cols, &cols, A.base(), &cols, ipiv, &info);
    LAPACK_dgetrs("T", &cols, &mcols, A.base(), &cols, ipiv, tmpGMp->base(),
      &mld, &info);
    delete [] ipiv;
    // *this= Transpose(tmpGMp out)
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        get(i,j) = tmpGMp->get(j,i);
    }
  delete tmpGMp;  
  }


// here we must be careful for ld
void GeneralMatrix::zeros()
  {
  if (ld == rows)
    data.zeros();
  else {
    for (int i = 0; i < rows; i++) 
      for (int j = 0; j < cols; j++)
        get(i,j) = 0;
    }
  }

void GeneralMatrix::unit()
  {
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      if (i == j)
        get(i,j) = 1.0;
      else
        get(i,j) = 0.0;
  }

void GeneralMatrix::nans()
  {
  for (int i = 0; i < rows; i++) 
    for (int j = 0; j < cols; j++)
      get(i,j) = std::numeric_limits<double>::quiet_NaN();
  }

void GeneralMatrix::infs()
  {
  for (int i = 0; i < rows; i++) 
    for (int j = 0; j < cols; j++)
      get(i,j) = std::numeric_limits<double>::infinity();
  }


// here we must be careful for ld
void GeneralMatrix::mult(double a)
  {
  if (ld == rows)
    data.mult(a);
  else {
    for (int i = 0; i < rows; i++) 
      for (int j = 0; j < cols; j++)
        get(i,j) *= a;
    }
  }

// here we must be careful for ld
void GeneralMatrix::add(double a, const ConstGeneralMatrix& m)
  {
  if (m.numRows() != rows || m.numCols() != cols)
    throw SYLV_MES_EXCEPTION("Matrix has different size in GeneralMatrix::add.");
  
  if (ld == rows && m.ld == m.rows)
    data.add(a, m.data);
  else {
    for (int i = 0; i < rows; i++) 
      for (int j = 0; j < cols; j++)
        get(i,j) += a*m.get(i,j);
    }
  }

void GeneralMatrix::add(double a, const ConstGeneralMatrix& m, const char* dum)
  {
  if (m.numRows() != cols || m.numCols() != rows)
    throw SYLV_MES_EXCEPTION("Matrix has different size in GeneralMatrix::add.");
  
  for (int i = 0; i < rows; i++) 
    for (int j = 0; j < cols; j++)
      get(i,j) += a*m.get(j,i);
  }


bool GeneralMatrix::isDiff(const GeneralMatrix& m, const double tol=0.0)const
  {
  if (m.numRows() != rows || m.numCols() != cols)
    throw SYLV_MES_EXCEPTION("Matrix has different size in GeneralMatrix::isDiff.");
  for (int i = 0; i < rows; i++) 
    for (int j = 0; j < cols; j++)
      if (fabs(get(i,j) - m.get(i,j))>tol)
        return true;
    return false;
  }

bool GeneralMatrix::isDiffSym(const GeneralMatrix& m, const double tol=0.0)const
  {
  if (m.numRows() != rows || m.numCols() != cols || m.numRows() != cols || m.numCols() != rows)
    throw SYLV_MES_EXCEPTION("Matrix has different size or not square in GeneralMatrix::isDiffSym.");
  for (int i = 0; i < cols; i++) 
    for (int j = 0; i+j < cols ; j++) // traverse the upper triangle only
      if (fabs(get(j,j+i) - m.get(j,j+i))>tol) // along diagonals where higher changes occur
        return true;
    return false;
  }


/* x = scalar(a)*x + scalar(b)*this*d */
void GeneralMatrix::multVec(double a, Vector& x, double b, const ConstVector& d) const
  {
  if (x.length() != rows || cols != d.length()) {
    throw SYLV_MES_EXCEPTION("Wrong dimensions for vector multiply.");
    }
  if (rows > 0) {
    int mm = rows;
    int nn = cols;
    double alpha = b;
    int lda = ld;
    int incx = d.skip();
    double beta = a;
    int incy = x.skip();
    BLAS_dgemv("N", &mm, &nn, &alpha, data.base(), &lda, d.base(), &incx,
      &beta, x.base(), &incy);
    }
  
  }

void GeneralMatrix::copy(const ConstGeneralMatrix& m, int ioff, int joff)
  {
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      get(i,j) = m.get(i+ioff,j+joff);
  }


void 
GeneralMatrix::copy(const ConstGeneralMatrix& m)
  {
  memcpy(data.base() ,m.getData().base() ,m.numCols()*m.numRows()*sizeof(double));
  };

void 
GeneralMatrix::copy(const GeneralMatrix& m)
  {
  memcpy(data.base(),m.getData().base(),m.numCols()*m.numRows()*sizeof(double));
  };

void 
GeneralMatrix::copyColumns(const GeneralMatrix& m, int istart, int iend, int ito)
  {
  if ((rows!=m.numRows())|| istart<iend|| istart> m.numCols()-1 || iend> m.numCols()-1
    || ito> cols-1 )
    throw SYLV_MES_EXCEPTION("Wrong dimensions for copying matrix columns.");

  memcpy(data.base()+ito*rows*sizeof(double)
    ,m.getData().base()+istart*rows*sizeof(double)
    ,(iend-istart+1)*rows*sizeof(double));
  };



/* emulates Matlab command A(a,b)=B(c,d) where a,b,c,d are vectors or ":")*/
void 
GeneralMatrix::AssignByVectors(GeneralMatrix& a, const vector<int>& vToRows, const vector<int>& vToCols
      , const GeneralMatrix& b, const vector<int>& vrows, const vector<int>& vcols)
  {
  int nrows=0, ncols=0,  tonrows=0, toncols=0;
  const vector<int> *vpToCols=0, *vpToRows=0, *vpRows=0, *vpCols=0;
  vector<int> *tmpvpToCols=0, *tmpvpToRows=0, *tmpvpRows=0, *tmpvpCols=0;

  if (vToRows.size()==0 && vToCols.size()==0 &&vrows.size()==0 && vcols.size()==0)
    a=b;
  else
    {
    if (vToRows.size()==0)
      {
      tonrows=a.numRows();
//      vpToRows=new vector<int>(tonrows);
      tmpvpToRows=new vector<int>(tonrows);
      for (int i=0;i<tonrows;++i)
        (*tmpvpToRows)[i]=i+1;
      vpToRows=(const vector<int>*)tmpvpToRows;
      }
    else
      {
      for (int i=0;i<vToRows.size();++i)
        {
        if (vToRows[i]>0) 
          tonrows++;
        else 
          throw SYLV_MES_EXCEPTION("Non-positive indices in assignment by vector.");
        }
      vpToRows=&vToRows;
      }
    if (tonrows>a.numRows())
      throw SYLV_MES_EXCEPTION("Wrong dimensions for assignment by vector.");

    if (vToCols.size()==0)
      {
      toncols=a.numCols();
      tmpvpToCols=new vector<int>(toncols);
      for (int i=0;i<toncols;++i)
        (*tmpvpToCols)[i]=i+1;
      vpToCols=(const vector<int>*)tmpvpToCols;
      }
    else
      {
      for (int i=0;i<vToCols.size();++i)
        {
        if (vToCols[i]>0) 
          toncols++;
        else 
          throw SYLV_MES_EXCEPTION("Non-positive indices in assignment by vector.");
        }
      vpToCols=&vToCols;
      }
    if (toncols>a.numCols())
      throw SYLV_MES_EXCEPTION("Wrong dimensions for assignment by vector.");

    if (vrows.size()==0)
      {
      nrows=b.numRows();
      tmpvpRows=new vector<int>(nrows);
      for (int i=0;i<nrows;++i)
        (*tmpvpToRows)[i]=i+1;
      vpRows=(const vector<int>*)tmpvpRows;
      }
    else
      {
      for (int i=0;i<vrows.size();++i)
        {
        if (vrows[i]>0) 
          nrows++;
        else 
          throw SYLV_MES_EXCEPTION("Non-positive indices in assignment by vector.");
        }
      vpRows=&vrows;
      }
    if (nrows>b.numRows())
      throw SYLV_MES_EXCEPTION("Wrong dimensions for assignment by vector.");

    if (vcols.size()==0)
      {
      ncols=b.numCols();
      tmpvpCols=new vector<int>(ncols);
      for (int i=0;i<ncols;++i)
        (*tmpvpCols)[i]=i+1;
      vpCols=(const vector<int>*)tmpvpCols;
      }
    else
      {
      for (int i=0;i<vcols.size();++i)
        {
        if (vcols[i]>0) 
          ncols++;
        else 
          throw SYLV_MES_EXCEPTION("Non-positive indices in assignment by vector.");
        }
      vpCols=&vcols;
      }
    if (ncols>b.numCols())
      throw SYLV_MES_EXCEPTION("Wrong dimensions for assignment by vector.");

    if (tonrows!=nrows || toncols!=ncols)
      throw SYLV_MES_EXCEPTION("Wrong indices dimensions for assignment by vector.");

    if(!(nrows*ncols==0 || tonrows*toncols==0))
      {
      for (int i=0;i<nrows;++i)
        {
        for (int j=0;j<nrows;++j)
              a.get((*vpToRows)[i]-1,(*vpToCols)[j]-1)=b.get((*vpRows)[i]-1, (*vpCols)[j]-1);
        }
      }
    if (tmpvpToCols) delete(tmpvpToCols); 
    if (tmpvpToRows) delete(tmpvpToRows);
    if (tmpvpRows) delete(tmpvpRows);
    if (tmpvpCols) delete(tmpvpCols);
    }
  }


void GeneralMatrix::gemm(const char* transa, const ConstGeneralMatrix& a,
                         const char* transb, const ConstGeneralMatrix& b,
                         double alpha, double beta)
  {
  int opa_rows = a.numRows();
  int opa_cols = a.numCols();
  if (!strcmp(transa, "T")) {
    opa_rows = a.numCols();
    opa_cols = a.numRows();
    }
  int opb_rows = b.numRows();
  int opb_cols = b.numCols();
  if (!strcmp(transb, "T")) {
    opb_rows = b.numCols();
    opb_cols = b.numRows();
    }
  
  if (opa_rows != numRows() ||
    opb_cols != numCols() ||
    opa_cols != opb_rows) {
    throw SYLV_MES_EXCEPTION("Wrong dimensions for matrix multiplication.");
    }
  
  int m = opa_rows;
  int n = opb_cols;
  int k = opa_cols;
  int lda = a.ld;
  int ldb = b.ld;
  int ldc = ld;
  if (lda > 0 && ldb > 0 && ldc > 0) {
    BLAS_dgemm(transa, transb, &m, &n, &k, &alpha, a.data.base(), &lda,
      b.data.base(), &ldb, &beta, data.base(), &ldc); 
    } else if (numRows()*numCols() > 0) {
    if (beta == 0.0)
      zeros();
    else
      mult(beta);
      }
  }

void GeneralMatrix::gemm_partial_left(const char* trans, const ConstGeneralMatrix& m,
                                      double alpha, double beta)
  {
  int icol;
  for (icol = 0; icol + md_length < cols; icol += md_length) {
    GeneralMatrix incopy((const GeneralMatrix&)*this, 0, icol, rows, md_length);
    GeneralMatrix inplace((GeneralMatrix&)*this, 0, icol, rows, md_length);
    inplace.gemm(trans, m, "N", ConstGeneralMatrix(incopy), alpha, beta);
    }
  if (cols > icol) {
    GeneralMatrix incopy((const GeneralMatrix&)*this, 0, icol, rows, cols - icol);
    GeneralMatrix inplace((GeneralMatrix&)*this, 0, icol, rows, cols - icol);
    inplace.gemm(trans, m, "N", ConstGeneralMatrix(incopy), alpha, beta);
    }
  }

void GeneralMatrix::gemm_partial_right(const char* trans, const ConstGeneralMatrix& m,
                                       double alpha, double beta)
  {
  int irow;
  for (irow = 0; irow + md_length < rows; irow += md_length) {
    GeneralMatrix incopy((const GeneralMatrix&)*this, irow, 0, md_length, cols);
    GeneralMatrix inplace((GeneralMatrix&)*this, irow, 0, md_length, cols);
    inplace.gemm("N", ConstGeneralMatrix(incopy), trans, m, alpha, beta);
    }
  if (rows > irow) {
    GeneralMatrix incopy((const GeneralMatrix&)*this, irow, 0, rows - irow, cols);
    GeneralMatrix inplace((GeneralMatrix&)*this, irow, 0, rows - irow, cols);
    inplace.gemm("N", ConstGeneralMatrix(incopy), trans, m, alpha, beta);
    }
  }

ConstGeneralMatrix::ConstGeneralMatrix(const GeneralMatrix& m, int i, int j, int nrows, int ncols)
: data(m.getData(), j*m.getLD()+i, (ncols-1)*m.getLD()+nrows), rows(nrows), cols(ncols), ld(m.getLD())
  {
  // can check that the submatirx is fully in the matrix
  }

ConstGeneralMatrix::ConstGeneralMatrix(const ConstGeneralMatrix& m, int i, int j, int nrows, int ncols)
: data(m.getData(), j*m.getLD()+i, (ncols-1)*m.getLD()+nrows), rows(nrows), cols(ncols), ld(m.getLD())
  {
  // can check that the submatirx is fully in the matrix
  }


ConstGeneralMatrix::ConstGeneralMatrix(const GeneralMatrix& m)
: data(m.data), rows(m.rows), cols(m.cols), ld(m.ld) {}


double ConstGeneralMatrix::getNormInf() const
  {
  double norm = 0.0;
  for (int i = 0; i < numRows(); i++) {
    ConstVector rowi(data.base()+i, ld, cols);
    double normi = rowi.getNorm1();
    if (norm < normi)
      norm = normi;
    }
  return norm;
  }

double ConstGeneralMatrix::getNorm1() const
  {
  double norm = 0.0;
  for (int j = 0; j < numCols(); j++) {
    ConstVector colj(data.base()+ld*j, 1, rows);
    double normj = colj.getNorm1();
    if (norm < normj)
      norm = normj;
    }
  return norm;
  }

/* x = scalar(a)*x + scalar(b)*this*d */
void ConstGeneralMatrix::multVec(double a, Vector& x, double b, const ConstVector& d) const
  {
  if (x.length() != rows || cols != d.length()) {
    throw SYLV_MES_EXCEPTION("Wrong dimensions for vector multiply.");
    }
  if (rows > 0) {
    int mm = rows;
    int nn = cols;
    double alpha = b;
    int lda = ld;
    int incx = d.skip();
    double beta = a;
    int incy = x.skip();
    BLAS_dgemv("N", &mm, &nn, &alpha, data.base(), &lda, d.base(), &incx,
      &beta, x.base(), &incy);
    }
  
  }

void ConstGeneralMatrix::multVecTrans(double a, Vector& x, double b,
                                      const ConstVector& d) const
  {
  if (x.length() != cols || rows != d.length()) {
    throw SYLV_MES_EXCEPTION("Wrong dimensions for vector multiply.");
    }
  if (rows > 0) {
    int mm = rows;
    int nn = cols;
    double alpha = b;
    int lda = rows;
    int incx = d.skip();
    double beta = a;
    int incy = x.skip();
    BLAS_dgemv("T", &mm, &nn, &alpha, data.base(), &lda, d.base(), &incx,
      &beta, x.base(), &incy);
    }
  }

/* m = inv(this)*m */
void ConstGeneralMatrix::multInvLeft(const char* trans, int mrows, int mcols, int mld, double* d) const
  {
  if (rows != cols) {
    throw SYLV_MES_EXCEPTION("The matrix is not square for inversion.");
    }
  if (cols != mrows) {
    throw SYLV_MES_EXCEPTION("Wrong dimensions for matrix inverse mutliply.");
    }
  
  if (rows > 0) {
    GeneralMatrix inv(*this);
    int* ipiv = new int[rows];
    int info;
    LAPACK_dgetrf(&rows, &rows, inv.getData().base(), &rows, ipiv, &info);
    LAPACK_dgetrs(trans, &rows, &mcols, inv.base(), &rows, ipiv, d,
      &mld, &info);
    delete [] ipiv;
    }
  }

/* m = inv(this)*m */
void ConstGeneralMatrix::multInvLeft(GeneralMatrix& m) const
  {
  multInvLeft("N", m.numRows(), m.numCols(), m.getLD(), m.getData().base());
  }

/* m = inv(this')*m */
void ConstGeneralMatrix::multInvLeftTrans(GeneralMatrix& m) const
  {
  multInvLeft("T", m.numRows(), m.numCols(), m.getLD(), m.getData().base());
  }

/* d = inv(this)*d */
void ConstGeneralMatrix::multInvLeft(Vector& d) const
  {
  if (d.skip() != 1) {
    throw SYLV_MES_EXCEPTION("Skip!=1 not implemented in ConstGeneralMatrix::multInvLeft(Vector&)");
    }
  
  multInvLeft("N", d.length(), 1, d.length(), d.base());
  }

/* d = inv(this')*d */
void ConstGeneralMatrix::multInvLeftTrans(Vector& d) const
  {
  if (d.skip() != 1) {
    throw SYLV_MES_EXCEPTION("Skip!=1 not implemented in ConstGeneralMatrix::multInvLeft(Vector&)");
    }
  
  multInvLeft("T", d.length(), 1, d.length(), d.base());
  }


bool ConstGeneralMatrix::isFinite() const
  {
  for (int i = 0; i < numRows(); i++)
    for (int j = 0; j < numCols(); j++)
      if (! std::isfinite(get(i,j)))
        return false;
      return true;
  }

bool ConstGeneralMatrix::isZero() const
  {
  for (int i = 0; i < numRows(); i++)
    for (int j = 0; j < numCols(); j++)
      if (get(i,j) != 0.0)
        return false;
      return true;
  }

void ConstGeneralMatrix::print() const
  {
  printf("rows=%d, cols=%d\n",rows, cols);
  for (int i = 0; i < rows; i++) {
    printf("row %d:\n",i);
    for (int j = 0; j < cols; j++) {
      printf("%6.3g ",get(i,j));
      }
    printf("\n");
    }
  }
