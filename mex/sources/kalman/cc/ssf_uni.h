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

#ifndef SSF_UNI_H
#define SSF_UNI_H

#include "ssf.h"


class TScalar{
  public:
    virtual const double&operator[](int t)const= 0;
    virtual~TScalar(){}
    virtual int period()const= 0;
    virtual TScalar*clone()const= 0;
  };

;

class TScalarInvariant:public TScalar{
  protected:
    double s;
  public:
    TScalarInvariant(double ss)
      :s(ss){}
    TScalarInvariant(const TScalarInvariant&c)
      :s(c.s){}
    const double&operator[](int t)const
      {return s;}
    
    int period()const
      {return 1;}
    TScalar*clone()const
      {return new TScalarInvariant(*this);}
  };

;

class TScalarCycle:public TScalar{
  protected:
    double*const ss;
    bool*const flags;
    int num;
  public:
    TScalarCycle(int n);
    TScalarCycle(const TScalarCycle&c);
    
    TScalarCycle(const GeneralMatrix&m);
    
    TScalarCycle(const TMatrix&m);
    ~TScalarCycle();
    const double&operator[](int t)const;
    int period()const
      {return num;}
    TScalar*clone()const
      {return new TScalarCycle(*this);}
    void set(int t,double s);
  };

;

struct SSFormUni{
  TMatrix*const Z;
  TScalar*const H;
  TMatrix*const T;
  TMatrix*const R;
  TMatrix*const Q;
  const int m;
  const int r;
  
  SSFormUni(const TMatrix&zz,const TScalar&hh,const TMatrix&tt,
    const TMatrix&rr,const TMatrix&qq);
  SSFormUni(const GeneralMatrix&zz,double hh,
    const GeneralMatrix&tt,const GeneralMatrix&rr,
    const GeneralMatrix&qq);
  SSFormUni(const SSFormUni&ssfu);
  ~SSFormUni();
  };

;

#endif

