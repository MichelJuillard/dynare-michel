/*!
  \file DynareMain.cc
  \par Main file of Dynare Parser.
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

  // FIXME
  string basename = argv[1];
  basename.erase(basename.size() - 4, 4);

  mod_file->writeOutputFiles(basename, clear_all);

  delete mod_file;

  cout << "Parsing done" << endl;
  cout << "Starting Matlab computing ..." << endl;
  return 0;
}
