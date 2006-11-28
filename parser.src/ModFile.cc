#include "ModFile.hh"

ModFile::ModFile() : symbol_table(model_parameters),
                     variable_table(symbol_table, model_parameters),
                     numerical_initialization(symbol_table, model_parameters),
                     computing_tasks(symbol_table),
                     model_tree(symbol_table, variable_table, model_parameters, num_constants)
{
}
