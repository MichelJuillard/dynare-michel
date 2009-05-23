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

#ifndef SSF_H
#define SSF_H

#include "GeneralMatrix.h"


class TMatrix{
  public:
    virtual const GeneralMatrix&operator[](int t)const= 0;
    virtual GeneralMatrix&operator[](int t)= 0;
    virtual int numRows()const= 0;
    virtual int numCols()const= 0;
    virtual int period()const= 0;
    virtual bool isZero(int t)const= 0;
    virtual bool isUnit(int t)const= 0;
    virtual~TMatrix(){}
    virtual TMatrix*clone()const= 0;
  };


class TMatrixInvariant:public TMatrix,public GeneralMatrix{
  public:
    TMatrixInvariant(const GeneralMatrix&m)
      :GeneralMatrix(m){}
    TMatrixInvariant(const TMatrixInvariant&m)
      :GeneralMatrix(m){}
    const GeneralMatrix&operator[](int t)const
      {return*this;}
    GeneralMatrix&operator[](int t)
      {return*this;}
    int numRows()const
      {return GeneralMatrix::numRows();}
    int numCols()const
      {return GeneralMatrix::numCols();}
    int period()const
      {return 1;}
    bool isZero(int t)const
      {return false;}
    bool isUnit(int t)const
      {return false;}
    TMatrix*clone()const
      {return new TMatrixInvariant(*this);}
  };



class TMatrixCycle:public TMatrix{
  protected:
    GeneralMatrix**const matrices;
    int num;
    int nrows;
    int ncols;
  public:
    TMatrixCycle(int n,int nr,int nc);
    TMatrixCycle(const TMatrixCycle&m);
    
    TMatrixCycle(const GeneralMatrix&m);
    
    TMatrixCycle(const TMatrix&m,const char*dummy);
    ~TMatrixCycle();
    const GeneralMatrix&operator[](int t)const;
    GeneralMatrix&operator[](int t);
    int numRows()const
      {return nrows;}
    int numCols()const
      {return ncols;}
    int period()const
      {return num;}
    bool isZero(int t)const
      {return false;}
    bool isUnit(int t)const
      {return false;}
    TMatrix*clone()const
      {return new TMatrixCycle(*this);}
    void set(int t,const GeneralMatrix&m);
  };


class TMatrixPadUnit:public TMatrix{
  TMatrix*const tmat;
  int skip;
  GeneralMatrix unit;
  public:
    TMatrixPadUnit(const TMatrix&m,int s);
    TMatrixPadUnit(const TMatrixPadUnit&m)
      :tmat(m.tmat->clone()),skip(m.skip),unit(m.unit){}
    ~TMatrixPadUnit()
      {delete tmat;}
    const GeneralMatrix&operator[](int t)const;
    GeneralMatrix&operator[](int t);
    int numRows()const
      {return tmat->numRows();}
    int numCols()const
      {return tmat->numCols();}
    int period()const
      {return skip*tmat->period();}
    bool isZero(int t)const
      {return false;}
    bool isUnit(int t)const
      {return(t/skip)*skip!=t;}
    TMatrix*clone()const
      {return new TMatrixPadUnit(*this);}
  };


class TMatrixPadZero:public TMatrix{
  TMatrix*const tmat;
  int skip;
  GeneralMatrix zero;
  public:
    TMatrixPadZero(const TMatrix&m,int s);
    TMatrixPadZero(const TMatrixPadZero&m)
      :tmat(m.tmat->clone()),skip(m.skip),zero(m.zero){}
    ~TMatrixPadZero()
      {delete tmat;}
    const GeneralMatrix&operator[](int t)const;
    GeneralMatrix&operator[](int t);
    int numRows()const
      {return tmat->numRows();}
    int numCols()const
      {return tmat->numCols();}
    int period()const
      {return skip*tmat->period();}
    bool isUnit(int t)const
      {return false;}
    bool isZero(int t)const
      {return(t/skip)*skip!=t;}
    TMatrix*clone()const
      {return new TMatrixPadZero(*this);}
  };



struct SSForm{
  TMatrix*const Z;
  TMatrix*const H;
  TMatrix*const T;
  TMatrix*const R;
  TMatrix*const Q;
  const int p;
  const int m;
  const int r;
  
  SSForm(const TMatrix&zz,const TMatrix&hh,const TMatrix&tt,
    const TMatrix&rr,const TMatrix&qq);
  SSForm(const GeneralMatrix&zz,const GeneralMatrix&hh,
    const GeneralMatrix&tt,const GeneralMatrix&rr,
    const GeneralMatrix&qq);
  SSForm::SSForm(const SSForm&f);
  
  ~SSForm();
  };


struct MesEquation{
  GeneralMatrix y;
  TMatrix*const Z;
  TMatrix*const H;
  
  MesEquation(const GeneralMatrix&data,const GeneralMatrix&zz,
    const GeneralMatrix&hh);
  MesEquation(const GeneralMatrix&data,const TMatrix&zz,
    const TMatrix&hh);
  ~MesEquation();
  protected:
    void construct_invariant();
  };


;
#endif

