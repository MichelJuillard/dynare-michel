/*
 * Copyright (C) 2007-2009 Dynare Team
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

#ifndef ERROR_HANDLING
#define ERROR_HANDLING

#include <cstring>
#include <iostream>
#include <sstream>

using namespace std;

const int NO_ERROR_ON_EXIT = 0;
const int ERROR_ON_EXIT = 1;
class GeneralExceptionHandling
{
  string ErrorMsg;
public:
  GeneralExceptionHandling(string ErrorMsg_arg) : ErrorMsg(ErrorMsg_arg){};
  inline string
  GetErrorMsg()
    {
      return ErrorMsg;
    }
  inline void
  completeErrorMsg(string ErrorMsg_arg)
    {
      ErrorMsg += ErrorMsg_arg;
    }
};

class FloatingPointExceptionHandling : public GeneralExceptionHandling
{
public:
  FloatingPointExceptionHandling(string value) : GeneralExceptionHandling(string("Floating point error in bytecode: " + value))
   {};
};


class LogExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  LogExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log(X)"),
                                           value(value_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value << "\n";
     completeErrorMsg(tmp.str());
   };
};

class Log10ExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  Log10ExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log10(X)"),
                                           value(value_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value << "\n";
     completeErrorMsg(tmp.str());
   };
};

class DivideExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  DivideExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("a/X"),
                                   value1(value1_arg),
                                   value2(value2_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value2 << "\n";
     completeErrorMsg(tmp.str());
   };
};


class PowExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  PowExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("X^a"),
                                   value1(value1_arg),
                                   value2(value2_arg)
   {
     ostringstream tmp;
     tmp << " with X=" << value1 << "\n";
     completeErrorMsg(tmp.str());
   };
};


class FatalExceptionHandling : public GeneralExceptionHandling
{
public:
  FatalExceptionHandling(string ErrorMsg_arg) : GeneralExceptionHandling("Fatal error in bytecode:")
    {
      completeErrorMsg(ErrorMsg_arg);
    };
  FatalExceptionHandling() : GeneralExceptionHandling("")
    {
    };
};
#endif
