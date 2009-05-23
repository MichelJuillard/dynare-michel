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

#include "ssf_uni.h"
#include "ts_exception.h"


TScalarCycle::TScalarCycle(int n)
:ss(new double[n]),flags(new bool[n]),num(n)
  {
  for(int i= 0;i<num;i++)
    flags[i]= false;
  }

;

TScalarCycle::TScalarCycle(const TScalarCycle&c)
:ss(new double[c.num]),flags(new bool[c.num]),num(c.num)
  {
  for(int i= 0;i<num;i++){
    flags[i]= c.flags[i];
    ss[i]= c.ss[i];
    }
  }

;

TScalarCycle::TScalarCycle(const GeneralMatrix&m)
:ss(new double[m.numRows()]),flags(new bool[m.numRows()]),
num(m.numRows())
  {
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "Matrix is not diagonal in TScalarCycle diagonal constructor");
  for(int i= 0;i<m.numRows();i++){
    ss[i]= m.get(i,i);
    flags[i]= true;
    }
  }

;

TScalarCycle::TScalarCycle(const TMatrix&m)
:ss(new double[m.numRows()*m.period()]),
flags(new bool[m.numRows()*m.period()]),
num(m.numRows()*m.period())
  {
  TS_RAISE_IF(m.period()==0,
    "Infinite period in TScalarCycle diagonal constructor");
  TS_RAISE_IF(m.numRows()!=m.numCols(),
    "TMatrix is not diagonal in TScalarCycle diagonal constructor");
  for(int i= 0;i<m.period();i++)
    for(int j= 0;j<m.numRows();j++){
      ss[i*m.numRows()+j]= m[i].get(j,j);
      flags[i*m.numRows()+j]= true;
      }
  }


;

TScalarCycle::~TScalarCycle()
  {
  delete[]flags;
  delete[]ss;
  }

;

const double&TScalarCycle::operator[](int t)const
  {
  int i= (t-1)%num;
  TS_RAISE_IF(!flags[i],
    "The scalar has not been set in TScalarCycle::operator[]");
  return ss[i];
  }


;

void TScalarCycle::set(int t,double s)
  {
  int i= (t-1)%num;
  flags[i]= true;
  ss[i]= s;
  }

;

SSFormUni::SSFormUni(const TMatrix&zz,const TScalar&hh,const TMatrix&tt,
                     const TMatrix&rr,const TMatrix&qq)
                     :Z(zz.clone()),
                     H(hh.clone()),
                     T(tt.clone()),
                     R(rr.clone()),
                     Q(qq.clone()),
                     m(zz.numCols()),r(qq.numRows())
  {
  TS_RAISE_IF(T->numRows()!=m||T->numCols()!=m||
    R->numRows()!=m||R->numCols()!=r||
    Q->numCols()!=r,
    "Wrong TMatrix dimension in SSFormUni constructor");
  TS_RAISE_IF(Z->numRows()!=1,
    "Z is not univariate in SSFormUni constructor");
  }

;

SSFormUni::SSFormUni(const GeneralMatrix&zz,double hh,
                     const GeneralMatrix&tt,const GeneralMatrix&rr,
                     const GeneralMatrix&qq)
                     :Z(new TMatrixInvariant(zz)),
                     H(new TScalarInvariant(hh)),
                     T(new TMatrixInvariant(tt)),
                     R(new TMatrixInvariant(rr)),
                     Q(new TMatrixInvariant(qq)),
                     m(zz.numCols()),r(qq.numRows())
  {
  TS_RAISE_IF(T->numRows()!=m||T->numCols()!=m||
    R->numRows()!=m||R->numCols()!=r||
    Q->numCols()!=r,
    "Wrong TMatrix dimension in SSFormUni constructor");
  TS_RAISE_IF(Z->numRows()!=1,
    "Z is not univariate in SSFormUni constructor");
  }

SSFormUni::SSFormUni(const SSFormUni&ssfu)
:Z(ssfu.Z->clone()),
H(ssfu.H->clone()),
T(ssfu.T->clone()),
R(ssfu.R->clone()),
Q(ssfu.Q->clone()),
m(ssfu.m),r(ssfu.r)
  {}



SSFormUni::~SSFormUni()
  {
  delete Z;
  delete H;
  delete T;
  delete R;
  delete Q;
  }


;

