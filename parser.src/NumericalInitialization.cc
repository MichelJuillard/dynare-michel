#include "NumericalInitialization.hh"
#include "Interface.hh"

InitParamStatement::InitParamStatement(const string &param_name_arg,
                                       const NodeID param_value_arg,
                                       const SymbolTable &symbol_table_arg) :
  param_name(param_name_arg),
  param_value(param_value_arg),
  symbol_table(symbol_table_arg)
{
}

void
InitParamStatement::writeOutput(ostream &output, const string &basename) const
{
  int id = symbol_table.getID(param_name) + 1;
  output << "M_.params( " << id << " ) = ";
  param_value->writeOutput(output);
  output << ";" << endl;
  output << param_name << " = M_.params( " << id << " );\n";
}

InitOrEndValStatement::InitOrEndValStatement(const init_values_type &init_values_arg,
                                             const SymbolTable &symbol_table_arg) :
  init_values(init_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
InitOrEndValStatement::writeInitValues(ostream &output) const
{
  for(init_values_type::const_iterator it = init_values.begin();
      it != init_values.end(); it++)
    {
      const string &name = it->first;
      const NodeID expression = it->second;

      Type type = symbol_table.getType(name);
      int id = symbol_table.getID(name) + 1;

      if (type == eEndogenous)
        output << "oo_.steady_state";
      else if (type == eExogenous)
        output << "oo_.exo_steady_state";
      else if (type == eExogenousDet)
        output << "oo_.exo_det_steady_state";

      output << "( " << id << " ) = ";
      expression->writeOutput(output);
      output << ";" << endl;
    }
}

InitValStatement::InitValStatement(const init_values_type &init_values_arg,
                                   const SymbolTable &symbol_table_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg)
{
}

void
InitValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << interfaces::comment() << "\n" << interfaces::comment() << "INITVAL instructions \n"
         << interfaces::comment() << "\n";
  // Writing initval block to set initial values for variables
  output << "options_.initval_file = 0;\nendval_=0;\n";

  if (symbol_table.recur_nbr > 0)
    output << "recurs_ = zeros(" << symbol_table.recur_nbr << ", 1);\n";

  writeInitValues(output);

  output << "oo_.y_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];\n";
  output << "if M_.exo_nbr > 0;\n";
  output << "\too_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];\n";
  output <<"end;\n";
  output << "if M_.exo_det_nbr > 0;\n";
  output << "\too_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];\n";
  output <<"end;\n";
}


EndValStatement::EndValStatement(const init_values_type &init_values_arg,
                                 const SymbolTable &symbol_table_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg)
{
}


void
EndValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << interfaces::comment() << "\n" << interfaces::comment() << "ENDVAL instructions\n"
         << interfaces::comment() << "\n";
  // Writing endval block to set terminal values for variables
  output << "ys0_= oo_.steady_state;\nex0_ = oo_.exo_steady_state;\nrecurs0_ = recurs_;\nendval_ = 1;\n";

  writeInitValues(output);
}

HistValStatement::HistValStatement(const hist_values_type &hist_values_arg,
                                   const SymbolTable &symbol_table_arg) :
  hist_values(hist_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
HistValStatement::writeOutput(ostream &output, const string &basename) const
{
  output << interfaces::comment() << "\n" << interfaces::comment() << "HISTVAL instructions\n"
         << interfaces::comment() << "\n";

  for(hist_values_type::const_iterator it = hist_values.begin();
      it != hist_values.end(); it++)
    {
      const string &name = it->first.first;
      const int &lag = it->first.second;
      const NodeID expression = it->second;

      Type type = symbol_table.getType(name);
      int id = symbol_table.getID(name) + 1;

      if (type == eEndogenous)
        output << "oo_.endo_simul( " << id << ", M_.maximum_lag + " << lag + 1 << ") = ";
      else if (type == eExogenous)
        output << "oo_.exo_simul( M_.maximum_lag + " << lag + 1 << ", " << id << " ) = ";
      else if (type != eExogenousDet)
        output << "oo_.exo_det_simul( M_.maximum_lag + " << lag + 1 << ", " << id << " ) = ";

      expression->writeOutput(output);
      output << ";" << endl;
    }
}
