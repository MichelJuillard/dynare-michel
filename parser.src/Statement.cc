#include "Statement.hh"

Statement::~Statement()
{
}

NativeStatement::NativeStatement(const string &native_statement_arg) :
  native_statement(native_statement_arg)
{
}

void
NativeStatement::writeOutput(ostream &output) const
{
  output << native_statement << endl;
}

void
OptionsList::writeOutput(ostream &output) const
{
  for(num_options_type::const_iterator it = num_options.begin();
      it != num_options.end(); it++)
    output << "options_." << it->first << " = " << it->second << ";" << endl;

  for(paired_num_options_type::const_iterator it = paired_num_options.begin();
      it != paired_num_options.end(); it++)
    output << "options_." << it->first << " = [" << it->second.first << "; "
           << it->second.second << "];" << endl;

  for(string_options_type::const_iterator it = string_options.begin();
      it != string_options.end(); it++)
    output << "options_." << it->first << " = '" << it->second << "';" << endl;
}

void
OptionsList::clear()
{
  num_options.clear();
  paired_num_options.clear();
  string_options.clear();
}
