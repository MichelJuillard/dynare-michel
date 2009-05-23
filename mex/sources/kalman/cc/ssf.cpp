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

#include "ssf.h"
#include "ts_exception.h"
#include "utils.h"

#include <vector> 


TMatrixCycle::TMatrixCycle(int n,int nr,int nc)
:matrices(new GeneralMatrix*[n]),num(n),nrows(nr),ncols(nc)
  {
  for(int i= 0;i<num;i++)
    matrices[i]= NULL;
  }

;

TMatrixCycle::TMatrixCycle(const TMatrixCycle&m)
:matrices(new GeneralMatrix*[m.num]),num(m.num),
nrows(m.nrows),ncols(m.ncols)
  {
  for(int i= 0;i<num;i++)
    if(m.matrices[i])
      matrices[i]= new GeneralMatrix(*(m.matrices[i]));
    else
      matrices[i]= NULL;
  }

;

TMatrixCycle::TMatrixCycle(const GeneralMatrix&m)
:matrices(new GeneralMatrix*[m.numRows()]),num(m.numRows()),
nrows(1),ncols(m.numCols())
  {
  for(int i= 0;i<num;i++)
    matrices[i]= new GeneralMatrix(m,i,0,1,ncols);
  }

;

TMatrixCycle::TMatrixCycle(const TMatrix&m,const char*dummy)
:matrices(new GeneralMatrix*[m.numRows()*m.period()]),
num(m.numRows()*m.period()),nrows(1),ncols(m.numCols())
  {
  TS_RAISE_IF(m.period()==0,
    "Infinite period in TMatrixCycle constructor");
  for(int i= 0;i<m.period();i++)
    for(int j= 0;j<m.numRows();j++)
      matrices[i*m.numRows()+j]
      = new GeneralMatrix(m[i],j,0,1,ncols);
  }

;

TMatrixCycle::~TMatrixCycle()
  {
  for(int i= 0;i<num;i++)
    delete matrices[i];
  delete[]matrices;
  }

;

const GeneralMatrix&TMatrixCycle::operator[](int t)const
  {
  int i= (t-1)%num;
  TS_RAISE_IF(matrices[i]==NULL,
    "The matrix has not ben set in TMatrixCycle::operator[]");
  return*(matrices[i]);
  }

GeneralMatrix&TMatrixCycle::operator[](int t)
  {
  int i= (t-1)%num;
  TS_RAISE_IF(matrices[i]==NULL,
    "The matrix has not ben set in TMatrixCycle::operator[]");
  return*(matrices[i]);
  }

;

void TMatrixCycle::set(int t,const GeneralMatrix&m)
  {
  TS_RAISE_IF(m.numRows()!=numRows()||m.numCols()!=numCols(),
    "Wrong matrix dimensions for TMatrixCycle::set");
  int i= (t-1)%num;
  if(matrices[i])
    delete matrices[i];
  matrices[i]= new GeneralMatrix(m);
  }

;

TMatrixPadUnit::TMatrixPadUnit(const TMatrix&m,int s)
:tmat(m.clone()),skip(s),unit(m.numRows(),m.numRows())
  {
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "TMatrix not square in TMatrixPadUnit constructor");
  unit.zeros();
  for(int i= 0;i<numRows();i++)
    unit.get(i,i)= 1.0;
  }

;

const GeneralMatrix&TMatrixPadUnit::operator[](int t)const
  {
  if(isUnit(t))
    return unit;
  else
    return(*tmat)[t/skip];
  }


GeneralMatrix&TMatrixPadUnit::operator[](int t)
  {
  TS_RAISE_IF(isUnit(t),
    "Attempt to return non-const unit in TMatrixPadUnit::operator[]");
  return(*tmat)[t/skip];
  }


TMatrixPadZero::TMatrixPadZero(const TMatrix&m,int s)
:tmat(m.clone()),skip(s),zero(m.numRows(),m.numCols())
  {
  zero.zeros();
  }


const GeneralMatrix&TMatrixPadZero::operator[](int t)const
  {
  if(isZero(t))
    return zero;
  else
    return(*tmat)[t/skip];
  }


GeneralMatrix&TMatrixPadZero::operator[](int t)
  {
  TS_RAISE_IF(isZero(t),
    "Attempt to return non-const zero in TMatrixPadZero::operator[]");
  return(*tmat)[t/skip];
  }



SSForm::SSForm(const TMatrix&zz,const TMatrix&hh,const TMatrix&tt,
               const TMatrix&rr,const TMatrix&qq)
               :Z(zz.clone()),
               H(hh.clone()),
               T(tt.clone()),
               R(rr.clone()),
               Q(qq.clone()),
               p(zz.numRows()),m(zz.numCols()),r(qq.numRows())
  {
  TS_RAISE_IF(T->numRows()!=m||T->numCols()!=m||
    H->numRows()!=p||H->numCols()!=p||
    R->numRows()!=m||R->numCols()!=r||
    Q->numCols()!=r,
    "Wrong TMatrix dimension in SSForm constructor");
  }


SSForm::SSForm(const GeneralMatrix&zz,const GeneralMatrix&hh,
               const GeneralMatrix&tt,const GeneralMatrix&rr,
               const GeneralMatrix&qq)
               :Z(new TMatrixInvariant(zz)),
               H(new TMatrixInvariant(hh)),
               T(new TMatrixInvariant(tt)),
               R(new TMatrixInvariant(rr)),
               Q(new TMatrixInvariant(qq)),
               p(zz.numRows()),m(zz.numCols()),r(qq.numRows())
  {
  TS_RAISE_IF(T->numRows()!=m||T->numCols()!=m||
    H->numRows()!=p||H->numCols()!=p||
    R->numRows()!=m||R->numCols()!=r||
    Q->numCols()!=r,
    "Wrong TMatrix dimension in SSForm constructor");
  }


SSForm::SSForm(const SSForm&f)
:Z(f.Z->clone()),
H(f.H->clone()),
T(f.T->clone()),
R(f.R->clone()),
Q(f.Q->clone()),
p(f.p),m(f.m),r(f.r)
  {}



SSForm::~SSForm()
  {
  delete Z;
  delete H;
  delete T;
  delete R;
  delete Q;
  }


MesEquation::MesEquation(const GeneralMatrix&data,const TMatrix&zz,
                         const TMatrix&hh)
                         :y(data),
                         Z((zz.period()*hh.period()==1)?(TMatrix*)new TMatrixInvariant(zz[1]):
(zz.period()*hh.period()==0)?(TMatrix*)new TMatrixCycle(y.numCols(),
                                                        zz.numRows(),zz.numCols())
                                                        :(TMatrix*)new TMatrixCycle(zz.period()*hh.period(),zz.numRows(),zz.numCols())),
                                                        H((zz.period()*hh.period()==1)?(TMatrix*)new TMatrixInvariant(hh[1]):
(zz.period()*hh.period()==0)?(TMatrix*)new TMatrixCycle(y.numCols(),
                                                        hh.numRows(),hh.numCols())
                                                        :(TMatrix*)new TMatrixCycle(zz.period()*hh.period(),hh.numRows(),hh.numCols()))
  {
  TS_RAISE_IF(y.numRows()!=Z->numRows()||y.numRows()!=H->numRows()||
    y.numRows()!=H->numCols(),
    "Incompatible dimension in MesEquation constructor");
  
  int mper= zz.period()*hh.period();
  if(mper==1){
    construct_invariant();
    }else{
    std::vector<NormCholesky*> chols;
    int per= (mper==0)?y.numCols():mper;
    for(int t= 1;t<=per;t++){
      
      GeneralMatrix ycol(y,0,t-1,y.numRows(),1);
      int hi= t;
      if(hh.period()> 0)
        hi= (t-1)%hh.period()+1;
      
      ;
      
      NormCholesky*ch;
      if(hh.period()==0){
        ch= new NormCholesky(hh[t]);
        }else if(hi-1>=(int)chols.size()){
        ch= new NormCholesky(hh[t]);
        chols.push_back(ch);
          }else{
          ch= chols[hi-1];
            }
          
          ;
          
          ch->getL().multInvLeftUnit(ycol);
          if(t-1<mper){
            GeneralMatrix Zt(zz[t]);
            ch->getL().multInvLeftUnit(Zt);
            ((TMatrixCycle*)Z)->set(t,Zt);
            GeneralMatrix Ht(hh.numRows(),hh.numRows());
            Ht.zeros();
            for(int i= 0;i<Ht.numRows();i++)
              Ht.get(i,i)= ch->getD()[i];
            ((TMatrixCycle*)H)->set(t,Ht);
            }
          
          
          
          ;
          if(hh.period()==0)
            delete ch;
      }
    for(unsigned int i= 0;i<chols.size();i++)
      delete chols[i];
      }
  }

MesEquation::MesEquation(const GeneralMatrix&data,const GeneralMatrix&zz,
                         const GeneralMatrix&hh)
                         :y(data),
                         Z(new TMatrixInvariant(zz)),
                         H(new TMatrixInvariant(hh))
  {
  TS_RAISE_IF(y.numRows()!=Z->numRows()||y.numRows()!=H->numRows()||
    y.numRows()!=H->numCols(),
    "Incompatible dimension in MesEquation constructor");
  
  construct_invariant();
  }


MesEquation::~MesEquation()
  {
  delete Z;
  delete H;
  }

void MesEquation::construct_invariant()
  {
  if(!TSUtils::isDiagonal((*H)[1])){
    NormCholesky chol((*H)[1]);
    chol.getL().multInvLeftUnit(y);
    chol.getL().multInvLeftUnit((*Z)[1]);
    (*H)[1].zeros();
    for(int i= 0;i<H->numRows();i++)
      (*H)[1].get(i,i)= chol.getD()[i];
    }
  }

;

