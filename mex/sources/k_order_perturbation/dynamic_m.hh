/*
 * Copyright (C) 2010-2011 Dynare Team
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

#ifndef _DYNAMIC_M_HH
#define _DYNAMIC_M_HH

#include "dynamic_abstract_class.hh"

#include "mex.h"
#include <dynmex.h>

/**
 * handles calls to <model>_dynamic.m
 *
 **/
class DynamicModelMFile : public DynamicModelAC
{
private:
  const string DynamicMFilename;
  const static int nlhs_dynamic = 4;
  const static int nrhs_dynamic = 5;
public:
  DynamicModelMFile(const string &modName) throw (DynareException);
  virtual ~DynamicModelMFile();
  void eval(const Vector &y, const Vector &x, const Vector &params, const Vector &ySteady,
            Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) throw (DynareException);
};
#endif
