/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the SymbolTable class methodes.
*/

#include <iostream>
#include <algorithm>
#include <sstream>

#include "SymbolTable.hh"
#include "Interface.hh"
using namespace std;

SymbolTable::SymbolTable() : endo_nbr(0), exo_nbr(0), exo_det_nbr(0), parameter_nbr(0),
                             model_local_variable_nbr(0), modfile_local_variable_nbr(0),
                             recur_nbr(0), unknown_function_nbr(0)
{
  name_table.resize(20);
  tex_name_table.resize(20);
}

int SymbolTable::AddSymbol(string name,Type type, string tex_name)
{
  symboltable[name].type = type;
  name_table[(int) type].push_back(name);
  tex_name_table[(int) type].push_back(tex_name);

  switch (type)
    {
    case eExogenous:
      symboltable[name].id = exo_nbr;
      return exo_nbr++;
    case eExogenousDet:
      symboltable[name].id = exo_det_nbr;
      return exo_det_nbr++;
    case eEndogenous:
      symboltable[name].id = endo_nbr;
      return endo_nbr++;
    case eParameter:
      symboltable[name].id = parameter_nbr;
      return parameter_nbr++;
    case eRecursiveVariable:
      symboltable[name].id = recur_nbr;
      return recur_nbr++;
    case eModelLocalVariable:
      symboltable[name].id = model_local_variable_nbr;
      return model_local_variable_nbr++;
    case eModFileLocalVariable:
      symboltable[name].id = modfile_local_variable_nbr;
      return modfile_local_variable_nbr++;
    case eUnknownFunction:
      symboltable[name].id = unknown_function_nbr;
      return unknown_function_nbr++;
    }
  // should never happen
  return -1;
}

int SymbolTable::AddSymbolDeclar(string name,Type type, string tex_name)
{
  //Testing if the symbol exist in the map
  if ( !Exist(name) )
    {
      //The symbol dosn't exist, adding it
      return AddSymbol(name,type, tex_name);
    }
  else
    {
      //The symbol exists, testing its type
      if (symboltable[name].type == type)
        {
          cout << "Warning : symbol " << name << " declared more than once.\n";
          return getID(name);
        }
      else
        {
          string msg = "symbol " + name + " declared more than once with different types.";
          (* error) (msg.c_str());
          return -1;
        }
    }

}

void SymbolTable::AddSymbolRange(string name,int nbr,Type type, string tex_name)
{
}

void  SymbolTable::ResetType(string name,Type new_type)
{
  symboltable[name].type = new_type;
}

void
SymbolTable::writeOutput(ostream &output) const
{
  if (exo_nbr > 0)
    {
      output << "M_.exo_names = '" << getNameByID(eExogenous, 0) << "';\n";
      output << "M_.exo_names_tex = '" << getTexNameByID(eExogenous, 0) << "';\n";
      for (int id = 1; id < exo_nbr; id++)
        {
          output << "M_.exo_names = " + interfaces::strvcat("M_.exo_names","'"+getNameByID(eExogenous, id)+"'") + ";\n";
          output << "M_.exo_names_tex = " + interfaces::strvcat("M_.exo_names_tex","'"+getTexNameByID(eExogenous, id)+"'") + ";\n";
        }
    }
  if (exo_det_nbr > 0)
    {
      output << "lgxdet_ = '" << getNameByID(eExogenousDet, 0) << "';\n";
      output << "lgxdet_tex_ = '" << getTexNameByID(eExogenousDet, 0) << "';\n";
      for (int id = 1; id < exo_det_nbr; id++)
        {
          output << "lgxdet_ = " + interfaces::strvcat("lgxdet_","'"+getNameByID(eExogenousDet, id)+"'") + ";\n";
          output << "lgxdet_tex_ = " + interfaces::strvcat("lgxdet_tex_","'"+getTexNameByID(eExogenousDet, id)+"'") + ";\n";
        }
    }
  if (endo_nbr > 0)
    {
      output << "M_.endo_names = '" << getNameByID(eEndogenous, 0) << "';\n";
      output << "M_.endo_names_tex = '" << getTexNameByID(eEndogenous, 0) << "';\n";
      for (int id = 1; id < endo_nbr; id++)
        {
          output << "M_.endo_names = " + interfaces::strvcat("M_.endo_names","'"+getNameByID(eEndogenous, id)+"'") + ";\n";
          output << "M_.endo_names_tex = " + interfaces::strvcat("M_.endo_names_tex","'"+getTexNameByID(eEndogenous, id)+"'") + ";\n";
        }
    }
  if (recur_nbr > 0)
    {
      output << "M_.recur_names = '" << getNameByID(eRecursiveVariable, 0) << "';\n";
      output << "M_.recur_names_tex = '" << getTexNameByID(eRecursiveVariable, 0) << "';\n";
      for (int id = 1; id < recur_nbr; id++)
        {
          output << "M_.recur_names = " + interfaces::strvcat("M_.recur_names","'"+getNameByID(eRecursiveVariable, id)+"'") + ";\n";
          output << "M_.recur_names_tex = " + interfaces::strvcat("M_.recur_names_tex","'"+getTexNameByID(eRecursiveVariable, id)+"'") + ";\n";
        }
    }
  if (parameter_nbr > 0)
    {
      output << "M_.param_names = '" << getNameByID(eParameter, 0) << "';\n";
      output << "M_.param_names_tex = '" << getTexNameByID(eParameter, 0) << "';\n";
      for (int id = 1; id < parameter_nbr; id++)
        {
          output << "M_.param_names = " + interfaces::strvcat("M_.param_names","'"+getNameByID(eParameter, id)+"'") + ";\n";
          output << "M_.param_names_tex = " + interfaces::strvcat("M_.param_names_tex","'"+getTexNameByID(eParameter, id)+"'") + ";\n";
        }
    }

  output << "M_.exo_det_nbr = " << exo_det_nbr << ";\n";
  output << "M_.exo_nbr = " << exo_nbr << ";\n";
  output << "M_.Sigma_e = zeros(" << exo_nbr
         << ", " << exo_nbr << ");\n";
  output << "M_.endo_nbr = " << endo_nbr << ";\n";
  output << "M_.recur_nbr = " << recur_nbr << ";\n";
  output << "M_.param_nbr = " << parameter_nbr << ";\n";
}
