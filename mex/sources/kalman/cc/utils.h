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

#ifndef UTILS_H
#define UTILS_H

#include "GeneralMatrix.h"


class LowerTriangle:public GeneralMatrix{
  public:
    LowerTriangle(const GeneralMatrix&m);
    LowerTriangle(const LowerTriangle&t)
      :GeneralMatrix(t){}
    void multInvLeft(GeneralMatrix&m)const;
    void multInvLeftUnit(GeneralMatrix&m)const;
  };

;

class NormCholesky{
  LowerTriangle L;
  Vector D;
  public:
    NormCholesky(const GeneralMatrix&m);
    NormCholesky(const NormCholesky&chol)
      :L(chol.L),D(chol.D){}
    const LowerTriangle&getL()const
      {return L;}
    const Vector&getD()const
      {return D;}
  };

;

class PLUFact{
  Vector inv;
  int*ipiv;
  int rows;
  double rcond;
  int detsign;
  int info;
  public:
    PLUFact(const GeneralMatrix&m);
    PLUFact(const PLUFact&plu);
    virtual~PLUFact()
      {delete[]ipiv;}
    void multInvLeft(GeneralMatrix&a)const;
    void multInvRight(GeneralMatrix&a)const;
    void multInvLeft(Vector&a)const;
    void multInvRight(Vector&a)const;
    bool isRegular()const
      {return info==0;}
    double getDeterminant()const;
    double getLogDeterminant()const;
    int getDetSign()const
      {return detsign;}
    int numRows()const
      {return rows;}
    double getRcond()const
      {return rcond;}
    void print()const;
  private:
    void PL_dgetrs(const char*trans,double*b,int ldb,int bcols)const;
    void calcDetSign();
  };

;

class VDVFact{
  GeneralMatrix V;
  Vector D;
  bool converged;
  public:
    VDVFact(const GeneralMatrix&m);
    const GeneralMatrix&getV()const
      {return V;}
    const Vector&getD()const
      {return D;}
    bool hasConverged()const
      {return converged;}
  };


;

struct TSUtils{
  static bool isDiagonal(const ConstGeneralMatrix&m);
  static bool isZero(const ConstGeneralMatrix&m);
  static bool hasNegativeDiagonal(const ConstGeneralMatrix&m);
  static bool isSymDiagDominant(const ConstGeneralMatrix&m);
  static double correctDefinitness(GeneralMatrix&m);
  static void correctSymmetricity(GeneralMatrix&m);
  };

;

#endif

