/*
 * Copyright (C) 2003-2008 Dynare Team
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

using namespace std;

#include "ParsingDriver.hh"
#include "ModFile.hh"

/*!
  \brief Main function of Dynare.
  \param argc Number of command line argumetns from runtime system
  \param argv Command line arguments from runtime system
*/
int
main(int argc, char** argv)
{
  if (argc < 2)
    {
      cerr << "Missing model file" << endl;
      cerr << "Dynare usage: dynare model_file [debug]" << endl;
      exit(-1);
    }

  ParsingDriver p;

  bool clear_all = true;

  // Parse options
  for (int arg = 2; arg < argc; arg++)
    {
      if (string(argv[arg]) == string("debug"))
        {
          p.trace_scanning = true;
          p.trace_parsing = true;
        }
      else
        if (string(argv[arg]) == string("noclearall"))
          clear_all = false;
    }

  cout << "Starting Dynare ..." << endl;
  cout << "Parsing your model file ..." << endl;

  // Launch parsing
  ModFile *mod_file = p.parse(argv[1]);

  // Run checking pass
  mod_file->checkPass();

  // Do computations
  mod_file->computingPass();

  // FIXME
  string basename = argv[1];
  basename.erase(basename.size() - 4, 4);

  mod_file->writeOutputFiles(basename, clear_all);

  delete mod_file;

  cout << "Parsing done" << endl;
  cout << "Starting Matlab computing ..." << endl;
  return 0;
}
