/*
* Copyright (C) 2008-2009 Dynare Team
*
* This file is part of Dynare.
*
* Dynare is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Dynare is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
*/

/* derived from c++kalman_filter library by O. Kamenik */

#include "utils.h"
#include "ts_exception.h"
#include "cppblas.h"
#include "cpplapack.h"

#include <math.h> 
#include <cmath> 
#include <float.h> 


LowerTriangle::LowerTriangle(const GeneralMatrix&m)
:GeneralMatrix(m)
  {
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "The matrix is not square in LowerTriangle constructor");
  }


void LowerTriangle::multInvLeft(GeneralMatrix&m)const
  {
  TS_RAISE_IF(numCols()!=m.numRows(),
    "Wrong dimensions of the matrix for LowerTriangle::multInvLeft");
  int mrows= m.numRows();
  int mcols= m.numCols();
  double alpha= 1.0;
  int ld= getLD();
  int ldm= m.getLD();
  BLAS_dtrsm("L","L","N","N",&mrows,&mcols,&alpha,getData().base(),
    &ld,m.getData().base(),&ldm);
  }

;

void LowerTriangle::multInvLeftUnit(GeneralMatrix&m)const
  {
  TS_RAISE_IF(numCols()!=m.numRows(),
    "Wrong dimensions of the matrix for LowerTriangle::multInvLeftUnit");
  int mrows= m.numRows();
  int mcols= m.numCols();
  double alpha= 1.0;
  int ld= getLD();
  int ldm= m.getLD();
  BLAS_dtrsm("L","L","N","U",&mrows,&mcols,&alpha,getData().base(),
    &ld,m.getData().base(),&ldm);
  }


;

NormCholesky::NormCholesky(const GeneralMatrix&a)
:L(a),D(a.numRows())
  {
  TS_RAISE_IF(a.numRows()!=a.numCols(),
    "The matrix is not square in NormCholesky constructor");
  
  int lrows= L.numRows();
  int ldl= L.getLD();
  int info;
  LAPACK_dpotrf("L",&lrows,L.getData().base(),&ldl,&info);
  TS_RAISE_IF(info> 0,
    "The matrix is not positive definite in NormCholesky constructor");
  TS_RAISE_IF(info<0,
    "Internal error in NormCholesky constructor");
  
  for(int i= 0;i<L.numRows();i++)
    for(int j= i+1;j<L.numCols();j++)
      L.get(i,j)= 0.0;
    
  for(int j= 0;j<L.numCols();j++){
    double d= L.get(j,j);
    Vector Lj(L,j);
    Lj.mult(1.0/d);
    D[j]= d*d;
    }
  
  }

;

PLUFact::PLUFact(const GeneralMatrix&m)
:inv(m.numRows()*m.numCols()),ipiv(new int[m.numRows()]),
rows(m.numRows())
  {
  TS_RAISE_IF(!m.isFinite(),
    "Matrix is not finite in PLUFact constructor");
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "Matrix not square in PLUFact constructor");
  inv= m.getData();
  LAPACK_dgetrf(&rows,&rows,inv.base(),&rows,ipiv,&info);
  TS_RAISE_IF(info<0,
    "Internal error in PLUFact constructor");
  
  double mnorm= m.getNormInf();
  double*work= new double[4*rows];
  int*iwork= new int[rows];
  int infotmp;
  LAPACK_dgecon("I",&rows,inv.base(),&rows,&mnorm,&rcond,work,
    iwork,&infotmp);
  delete[]iwork;
  delete[]work;
  TS_RAISE_IF(infotmp<0,
    "Internal error in PLUFact constructor");
  
  ;
  calcDetSign();
  }


PLUFact::PLUFact(const PLUFact&fact)
:inv(fact.inv),ipiv(new int[fact.rows]),
rows(fact.rows),rcond(fact.rcond),detsign(fact.detsign),info(fact.info)
  {
  memcpy(ipiv,fact.ipiv,rows*sizeof(int));
  }

PLUFact::PLUFact(const int nc,const int nr )
    :inv(nr*nc),ipiv(new int[nr]),rows(nr)
  {
  TS_RAISE_IF(nr!=nc,
    "Matrix not square in PLUFact constructor");
  }

const PLUFact&
PLUFact::operator = (const GeneralMatrix&m)
  {
  TS_RAISE_IF(!m.isFinite(),
    "Matrix is not finite in PLUFact assignement");
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "Matrix not square in PLUFact assignement");
  TS_RAISE_IF(m.numRows()!=rows,
    "Matrix not matching PLUFact size for assignement");
  inv= m.getData();
  LAPACK_dgetrf(&rows,&rows,inv.base(),&rows,ipiv,&info);
  TS_RAISE_IF(info<0,
    "Internal error in PLUFact assignement");
  
  double mnorm= m.getNormInf();
  double*work= new double[4*rows];
  int*iwork= new int[rows];
  int infotmp;
  LAPACK_dgecon("I",&rows,inv.base(),&rows,&mnorm,&rcond,work,
    iwork,&infotmp);
  delete[]iwork;
  delete[]work;
  TS_RAISE_IF(infotmp<0,
    "Internal error in PLUFact assignement");
  calcDetSign();
  return *this;
  }

;

void PLUFact::PL_dgetrs(const char*trans,double*b,int ldb,int bcols)const
  {
  if(rows> 0){
    int info;
    LAPACK_dgetrs(trans,&rows,&bcols,inv.base(),&rows,ipiv,b,&ldb,&info);
    TS_RAISE_IF(info<0,
      "Internal error in PLUFact::dgetrs");
    }
  }

;

void PLUFact::multInvLeft(GeneralMatrix&a)const
  {
  TS_RAISE_IF(rows!=a.numRows(),
    "Wrong dimension of the matrix in PLUFact::multInvLeft");
  PL_dgetrs("N",a.getData().base(),a.getLD(),a.numCols());
  }

;

void PLUFact::multInvRight(GeneralMatrix&a)const
  {
  GeneralMatrix atrans(a,"trans");
  TS_RAISE_IF(rows!=atrans.numRows(),
    "Wrong dimension of the matrix in PLUFact::multInvLeft");
  PL_dgetrs("T",atrans.getData().base(),atrans.getLD(),atrans.numCols());
  for(int i= 0;i<a.numRows();i++)
    for(int j= 0;j<a.numCols();j++)
      a.get(i,j)= atrans.get(j,i);
  }

;

void PLUFact::multInvLeft(Vector&a)const
  {
  TS_RAISE_IF(rows!=a.length(),
    "Wrong dimension of the vector in PLUFact::multInvLeft");
  TS_RAISE_IF(a.skip()!=1,
    "Not implemented error in PLUFact::multInvLeft");
  PL_dgetrs("N",a.base(),a.length(),1);
  }

;

void PLUFact::multInvRight(Vector&a)const
  {
  TS_RAISE_IF(rows!=a.length(),
    "Wrong dimension of the vector in PLUFact::multInvLeft");
  TS_RAISE_IF(a.skip()!=1,
    "Not implemented error in PLUFact::multInvLeft");
  PL_dgetrs("T",a.base(),a.length(),1);
  }


;

double PLUFact::getDeterminant()const
  {
  double res= 1;
  for(int i= 0;i<rows;i++)
    res*= std::abs(inv[(rows+1)*i]);
  return detsign*res;
  }

;

double PLUFact::getLogDeterminant()const
  {
  double res= 0;
  for(int i= 0;i<rows;i++)
    res+= log(std::abs(inv[(rows+1)*i]));
  TS_RAISE_IF(detsign==-1,
    "Negative determinant in PLUFact::getLogDeterminant");
  return res;
  }

;

void PLUFact::calcDetSign()
  {
  detsign= 1;
  
  for(int i= 0;i<rows;i++)
    if(ipiv[i]!=i+1)
      detsign*= -1;
    
  for(int i= 0;i<rows;i++)
    if(inv[i*(rows+1)]<0)
      detsign*= -1;
      
  }

;

void PLUFact::print()const
  {
  for(int i= 0;i<rows;i++)
    printf(" %d",ipiv[i]);
  printf("\n");
  for(int i= 0;i<rows;i++){
    for(int j= 0;j<rows;j++)
      printf(" %15.12g",inv[j*rows+i]);
    printf("\n");
    }
  }

;

VDVFact::VDVFact(const GeneralMatrix&m)
:V(m),D(m.numRows())
  {
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "Matrix is not square in VDVFact constructor");
  
  int n= m.numRows();
  int lda= V.getLD();
  double tmpwork;
  int lwork= -1;
  int info;
  LAPACK_dsyev("V","U",&n,V.base(),&lda,D.base(),&tmpwork,&lwork,&info);
  lwork= (int)tmpwork;
  double*work= new double[lwork];
  LAPACK_dsyev("V","U",&n,V.base(),&lda,D.base(),work,&lwork,&info);
  delete[]work;
  
  TS_RAISE_IF(info<0,
    "Internal error in VDVFact constructor");
  converged= true;
  if(info)
    converged= false;
  }


;

bool TSUtils::isDiagonal(const ConstGeneralMatrix&m)
  {
  bool res= (m.numCols()==m.numRows());
  for(int i= 0;i<m.numRows()&&res;i++)
    for(int j= i+1;j<m.numCols()&&res;j++)
      if(m.get(i,j)!=0.0||m.get(j,i)!=0.0)
        res= false;
      return res;
  }

;

bool TSUtils::isZero(const ConstGeneralMatrix&m)
  {
  bool res= true;
  for(int i= 0;i<m.numRows()&&res;i++)
    for(int j= 0;j<m.numCols()&&res;j++)
      if(m.get(i,j)!=0.0)
        res= false;
      return res;
  }

;

bool TSUtils::hasNegativeDiagonal(const ConstGeneralMatrix&m)
  {
  int r= m.numRows()<m.numCols()?m.numRows():m.numCols();
  bool res= false;
  for(int i= 0;i<r&&!res;i++)
    res= m.get(i,i)<0.0;
  return res;
  }

;

bool TSUtils::isSymDiagDominant(const ConstGeneralMatrix&m)
  {
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "The matrix is not square in TSUtils::isSymDiagDominant");
  
  bool res= true;
  for(int i= 0;i<m.numRows()&&res;i++)
    for(int j= i+1;j<m.numCols()&&res;j++)
      res= 2*std::abs(m.get(i,j))<=
      std::abs(m.get(i,i))+std::abs(m.get(j,j));
    return res;
  }

;

double TSUtils::correctDefinitness(GeneralMatrix&m)
  {
  VDVFact f(m);
  if(!f.hasConverged())
    return-1;
  
  Vector d(f.getD());
  double correct= 0;
  int i= 0;
  while(i<d.length()&&d[i]<2*DBL_EPSILON){
    correct+= d[i]*d[i];
    d[i]= 0.0;
    i++;
    }
  
  m= f.getV();
  for(int i= 0;i<d.length();i++){
    Vector mi(m,i);
    mi.mult(d[i]);
    }
  m.multRightTrans(f.getV());
  
  return sqrt(correct);
  }

;

void TSUtils::correctSymmetricity(GeneralMatrix&m)
  {
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "Matrix is not square in TSUtils::correctSymmetricity");
  GeneralMatrix tmp((const GeneralMatrix&)m,"trans");
  m.add(1.0,tmp);
  m.mult(0.5);
  }


;

