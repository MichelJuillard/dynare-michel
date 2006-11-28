/*!
  \file DynareMain.cc
  \par Main file of Dynare Parser.
*/

using namespace std;

#include "ParsingDriver.hh"
#include "OutputFile.hh"
#include "ModFile.hh"

/*!
  \brief Main function of Dynare.
  \param argc Number of command line argumetns from runtime system
  \param argv Command line arguments from runtime system
*/
int
main(int argc, char** argv)
{
  OutputFile output_file;
  ostringstream output;

  if (argc < 2)
    {
      cerr << "Missing model file" << endl;
      cerr << "Dynare usage: dynare model_file [debug]" << endl;
      exit(-1);
    }

  ParsingDriver p;

  // Sets string output of parser
  p.setoutput(&output);

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
          output_file.clear_all = false;
    }

  cout << "Starting Dynare ..." << endl;
  cout << "Parsing your model file ..." << endl;

  // Launch parsing
  ModFile *mod_file = p.parse(argv[1]);

  // Execute final instructions
  p.finish();

  string name = argv[1];
  name.erase(name.size() - 4,4);
  // Opening and init main Output file (.m or .sci file)
  output_file.Open(name, mod_file);
  // Writing remaining string output to output file
  output_file.Save(output, mod_file);

  delete mod_file;

  cout << "Parsing done" << endl;
  cout << "Starting Matlab computing ..." << endl;
  return 0;
}
