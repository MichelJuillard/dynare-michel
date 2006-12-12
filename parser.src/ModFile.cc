#include "ModFile.hh"
#include "Interface.hh"

ModFile::ModFile() : symbol_table(model_parameters),
                     variable_table(symbol_table, model_parameters),
                     model_tree(symbol_table, variable_table, model_parameters, num_constants),
                     order(-1), linear(-1)
{
}

ModFile::~ModFile()
{
  for(vector<Statement *>::iterator it = statements.begin();
      it != statements.end(); it++)
    delete (*it);
}

void
ModFile::addStatement(Statement *st)
{
  statements.push_back(st);
}

void
ModFile::writeOutputFiles(const string &basename, bool clear_all)
{
  ofstream mOutputFile;

  if (basename.size())
    {
      string fname(basename);
      fname += interfaces::function_file_extension();
      mOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!mOutputFile.is_open())
        {
          cerr << "Error: Can't open file " << fname
               << " for writing" << endl;
          exit(-1);
        }
    }
  else
    {
      cerr << "Error: Missing file name" << endl;
      exit(-1);
    }

  mOutputFile << interfaces::comment() << endl;
  mOutputFile << interfaces::comment() << "Status : main Dynare file " << endl;
  mOutputFile << interfaces::comment() << endl;
  mOutputFile << interfaces::comment() << "Warning : this file is generated automatically by Dynare" << endl;
  mOutputFile << interfaces::comment() << "          from model file (.mod)" << endl << endl;

  if (clear_all)
    mOutputFile << "clear all" << endl;
  mOutputFile << "tic;" << endl;
  mOutputFile << "global M_ oo_ exedet_ exdet_ recur_ recurs_ " << endl;
  mOutputFile << "global options_ endval_" << endl;
  mOutputFile << "global ys0_ recurs0_ ex0_ ct_" << endl;
  mOutputFile << "options_ = [];" << endl;
  mOutputFile << "M_.fname = '" << basename << "';" << endl;
  mOutputFile << interfaces::comment() << endl;
  mOutputFile << interfaces::comment() << "Some global variables initialisation" << endl;
  mOutputFile << interfaces::comment() << endl;
  mOutputFile << "global_initialization;" << endl;
  mOutputFile << "diary off;" << endl << "warning off;" << endl << endl;
  mOutputFile << interfaces::delete_file(basename + ".log") << ";" << endl;
  mOutputFile << "warning on;" << endl << "warning backtrace;" << endl;
  mOutputFile << "logname_ = '" << basename << ".log';" << endl;
  mOutputFile << "diary '" << basename << ".log';" << endl;

  if (model_tree.offset == 0)
    {
      mOutputFile << "if ";
      mOutputFile << interfaces::file_exist(basename + "_static.c)") << endl;
      mOutputFile << "   clear " << basename << "_static" << endl;
      mOutputFile << "   " << interfaces::compile(basename +"_static.c") << endl;
      mOutputFile << "end" << endl;
      mOutputFile << "if ";
      mOutputFile << interfaces::file_exist(basename + "_dynamic.c)") << endl;
      mOutputFile << "   clear " << basename << "_dynamic" << endl;
      mOutputFile << "   " + interfaces::compile(basename+"_dynamic.c") << endl;
      mOutputFile << "end" << endl;
    }
  else
    {
      mOutputFile << "erase_compiled_function('" + basename +"_static');" << endl;
      mOutputFile << "erase_compiled_function('" + basename +"_dynamic');" << endl;
      mOutputFile << interfaces::load_model_function_files(basename);
    }

  symbol_table.writeOutput(mOutputFile);

  if (linear == 1)
    mOutputFile << "options_.linear = 1;" << endl;

  model_tree.writeOutput(mOutputFile, basename, order, linear);

  // Print statements
  for(vector<Statement *>::iterator it = statements.begin();
      it != statements.end(); it++)
    (*it)->writeOutput(mOutputFile);

  mOutputFile << "save('" << basename << "_results', 'oo_');" << endl;
  mOutputFile << "diary off" << endl;

  mOutputFile << endl << "disp(['Total computing time : ' sec2hms(round(toc)) ]);" << endl;
  mOutputFile.close();
}
