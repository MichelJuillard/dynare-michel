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

#ifndef STATE_INIT_H
#define STATE_INIT_H

#include "GeneralMatrix.h"


class StateInit{
  const int m;
  int ndiffuse;
  GeneralMatrix Pstar;
  GeneralMatrix Pinf;
  Vector a;
  public:
    StateInit(const GeneralMatrix&PPstar,const Vector&aa);
    StateInit(const GeneralMatrix&PPstar,const GeneralMatrix&PPinf,
      const Vector&aa);
    StateInit(const StateInit&init)
      :m(init.m),ndiffuse(init.ndiffuse),Pstar(init.Pstar),
      Pinf(init.Pinf),a(init.a){}
    virtual~StateInit(){}
    int getM()const
      {return m;}
    bool isDiffuse()const
      {return ndiffuse> 0;}
    const Vector&getA()const
      {return a;}
    const GeneralMatrix&getPstar()const
      {return Pstar;}
    const GeneralMatrix&getPinf()const
      {return Pinf;}
    int getNDiff()const
      {return ndiffuse;}
  };


;

#endif

