/*
 * Copyright (C) 2010-2013 Dynare Team
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

#include "static_dll.hh"

#include <sstream>

using namespace std;

StaticModelDLL::StaticModelDLL(const std::string &basename) throw (TSException)
{
  std::string fName;
#if !defined(__CYGWIN32__) && !defined(_WIN32)
  if (basename[0] != '/')
    fName = "./";
#endif
  fName += basename + "_static" + MEXEXT;

  try
    {
#if defined(__CYGWIN32__) || defined(_WIN32)
      staticHinstance = LoadLibrary(fName.c_str());
      if (staticHinstance == NULL)
        throw 1;
      Static = (StaticFn) GetProcAddress(staticHinstance, "Static");
      if (Static == NULL)
        {
          FreeLibrary(staticHinstance); // Free the library
          throw 2;
        }
#else // Linux or Mac
      staticHinstance = dlopen(fName.c_str(), RTLD_NOW);
      if ((staticHinstance == NULL) || dlerror())
        {
          cerr << dlerror() << endl;
          throw 1;
        }
      Static = (StaticFn) dlsym(staticHinstance, "Static");
      if ((Static  == NULL) || dlerror())
        {
          dlclose(staticHinstance); // Free the library
          cerr << dlerror() << endl;
          throw 2;
        }
#endif

    }
  catch (int i)
    {
      std::ostringstream msg;
      msg << "Error when loading " << fName << " (";
      if (i == 1)
        msg << "can't dynamically load the file";
      if (i == 2)
        msg << "can't locate the 'Static' symbol";
      msg << ")";
      throw TSException(__FILE__, __LINE__, msg.str());
    }
  catch (...)
    {
      throw TSException(__FILE__, __LINE__, std::string("Can't find Static function in ") + fName);
    }
}

StaticModelDLL::~StaticModelDLL()
{
#if defined(__CYGWIN32__) || defined(_WIN32)
  bool result = FreeLibrary(staticHinstance);
  if (result == 0)
    throw TSException(__FILE__, __LINE__, std::string("Can't free the *_static DLL"));
#else
  dlclose(staticHinstance);
#endif
}
