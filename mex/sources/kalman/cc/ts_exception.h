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

#ifndef TS_EXCEPTION_H
#define TS_EXCEPTION_H

#include <stdio.h> 
#include <string.h> 

#ifdef MATLAB		
#include "mex.h"
#endif

#define TS_RAISE(mes) \
throw TSException(__FILE__, __LINE__, mes);

#define TS_RAISE_IF(expr, mes) \
if (expr) throw TSException(__FILE__, __LINE__, mes);


class TSException{
  char fname[50];
  int lnum;
  char message[500];
  public:
    TSException(const char*f,int l,const char*mes)
      {
      strncpy(fname,f,50);fname[49]= '\0';
      strncpy(message,mes,500);message[499]= '\0';
      lnum= l;
      }
    virtual void print()const
      {
      printf("At %s:%d:%s\n",fname,lnum,message);
#ifdef MATLAB		
      mexPrintf("At %s:%d:%s\n",fname,lnum,message);
#endif
      }
    
    virtual const char*getMessage()const
      {return message;}
  };

;
#endif

