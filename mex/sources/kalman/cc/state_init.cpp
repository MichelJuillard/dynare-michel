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

#include "state_init.h"
#include "ts_exception.h"
#include "utils.h"


StateInit::StateInit(const GeneralMatrix&PPstar,const Vector&aa)
:m(PPstar.numRows()),ndiffuse(0),Pstar(PPstar),
Pinf(m,m),a(aa)
  {
  TS_RAISE_IF(Pstar.numRows()!=Pstar.numCols(),
    "Pstar not square in StateInit non-diffuse constructor");
  TS_RAISE_IF(m!=a.length(),
    "Bad length of initial state vector in StateInit non-diffuse constructor");
  Pinf.zeros();
  }

;

StateInit::StateInit(const GeneralMatrix&PPstar,const GeneralMatrix&PPinf,
                     const Vector&aa)
                     :m(PPstar.numRows()),ndiffuse(0),Pstar(PPstar),
                     Pinf(PPinf),a(aa)
  {
  TS_RAISE_IF(m!=Pstar.numCols()||m!=Pinf.numRows()||
    m!=Pinf.numCols()||m!=a.length(),
    "Wrong dimensions for StateInit diffuse constructor");
  TS_RAISE_IF(!TSUtils::isDiagonal(Pinf),
    "Pinf is not diagonal in StateInit diffuse constructor");
  
  for(int i= 0;i<m;i++)
    if(Pinf.get(i,i)!=0.0)
      ndiffuse++;
  }



;

