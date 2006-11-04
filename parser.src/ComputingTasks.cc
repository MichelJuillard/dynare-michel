/** \file
 * \version 1.0
 * \date 12/16/2003
 * \par This file implements the ComputingTasks class methodes.
 */
#include <iostream>
#include <sstream>

using namespace std;
//------------------------------------------------------------------------------
#include "ComputingTasks.hh"
#include "Interface.hh"
//------------------------------------------------------------------------------
//ostringstream	ComputingTasks::output;
//------------------------------------------------------------------------------
ComputingTasks::ComputingTasks()
{
  // Empty
}

//------------------------------------------------------------------------------
ComputingTasks::~ComputingTasks()
{
  // Empty
}

//------------------------------------------------------------------------------
void ComputingTasks::setOutput(ostringstream* iOutput)
{
  output = iOutput;
}

//------------------------------------------------------------------------------
void ComputingTasks::set(void)
{
  // Empty
}

//------------------------------------------------------------------------------
void ComputingTasks::setSteady(void)
{
  *output << "steady;\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setCheck(void)
{
  *output << "check;\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setSimul(void)
{
  *output << "simul(oo_.dr);\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setStochSimul(string tmp1)

{
  *output << tmp1;
  *output << "stoch_simul(var_list_);\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setOption(string name, string value)
{
  *output << "options_." << name << " = " << value << ";\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setOption(string name, string value1, string value2)
{
  *output << "options_." << name << " = [" << value1 << "; " << value2 << "];\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::runEstimation(string tmp1)
{
  *output << tmp1;
  *output << "dynare_estimation(var_list_);\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::runRplot(string tmp1)
{
  *output << tmp1;
  *output << "rplot(var_list_);\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setEstimationInit(void)
{
  *output << "global estim_params_\n";
  *output << "var_list_ = [];\n";
  *output << "estim_params_.var_exo = [];\n";
  *output << "estim_params_.var_endo = [];\n";
  *output << "estim_params_.corrx = [];\n";
  *output << "estim_params_.corrn = [];\n";
  *output << "estim_params_.param_names = [];\n";
  *output << "estim_params_.user_param_names = [];\n";
  *output << "estim_params_.param_vals = [];\n";
  *output << "M_.H = 0;\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setOptimOptions(string str1, string str2, int task)
{
  static string optim_string;
  static int start;
  switch(task)
    {
    case 1:
      optim_string = "options_.optim_opt = '";
      start = 0;
      return;
    case 2:
      if (start > 0)
        {
          optim_string += ",";
        }
      else
        {
          start = 1;
        }
      optim_string += "''";
      optim_string += str1;
      optim_string += "'',";
      if (str2[0] >= 'A' && str2[0] <= 'z')
        {
          optim_string += "''";
          optim_string += str2;
          optim_string += "''";
        }
      else
        {
          optim_string += str2;
        }
      return;
    case 3:
      optim_string += "';\n";
      *output << optim_string;
    }
}

//------------------------------------------------------------------------------
void ComputingTasks::setEstimatedElements(void)
{
  if (!SymbolTable::Exist(EstimParams->name))
    {
      string msg = "Unknown symbol: "+EstimParams->name;
      error(msg.c_str());
    }
  if (SymbolTable::isReferenced(EstimParams->name) == eNotReferenced & EstimParams->name != "dsge_prior_weight")
    {
      return;
    }
  if ((EstimParams->init_val).size() == 0)
    {
      EstimParams->init_val = EstimParams->mean;
    }
  switch(EstimParams->type)
    {
    case 1:
      if( SymbolTable::getType(EstimParams->name) == eExogenous)
        {
          *output << "estim_params_.var_exo = [estim_params_.var_exo; ";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eEndogenous)
        {
          *output << "estim_params_.var_endo = [estim_params_.var_endo; ";
        }
      *output << SymbolTable::getID(EstimParams->name)+1;
      break;
    case 2:
      *output << "estim_params_.param_vals = [estim_params_.param_vals; ";
      *output << SymbolTable::getID(EstimParams->name)+1;
      break;
    case 3:
      if( SymbolTable::getType(EstimParams->name) == eExogenous)
        {
          *output << "estim_params_.corrx = [estim_params_.corrx; ";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eEndogenous)
        {
          *output << "estim_params_.corrn = [estim_params_.corrn; ";
        }
      *output << SymbolTable::getID(EstimParams->name)+1;
      *output << " " << SymbolTable::getID(EstimParams->name2)+1;
      break;
    }
  *output << " " << EstimParams->init_val << " " <<  EstimParams->low_bound << " " <<
    EstimParams->up_bound << " " <<  EstimParams->prior << " ";
  *output <<  EstimParams->mean << " " <<  EstimParams->std << " " <<
    EstimParams->p3 << " " <<  EstimParams->p4  << " " <<  EstimParams->jscale << "];\n";
  EstimParams->clear();
}

//------------------------------------------------------------------------------
void ComputingTasks::setEstimatedInitElements(void)
{
  if (!SymbolTable::Exist(EstimParams->name))
    {
      string msg = "Unknown symbol: "+EstimParams->name;
      error(msg.c_str());
    }
  if (SymbolTable::isReferenced(EstimParams->name) == eNotReferenced)
    {
      return;
    }
  if ((EstimParams->init_val).size() == 0)
    {
      EstimParams->init_val = EstimParams->mean;
    }
  if (EstimParams->type < 3)
    {
      if( SymbolTable::getType(EstimParams->name) == eExogenous)
        {
          *output << "tmp1 = find(estim_params_.var_exo(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ");\n";
          *output << "estim_params_.var_exo(tmp1,2) = " << EstimParams->init_val << ";\n";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eEndogenous)
        {
          *output << "tmp1 = find(estim_params_.var_endo(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ");\n";
          *output << "estim_params_.var_endo(tmp1,2) = " << EstimParams->init_val << ";\n";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eParameter)
        {
          *output << "tmp1 = find(estim_params_.param_vals(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ");\n";
          *output << "estim_params_.param_vals(tmp1,2) = " << EstimParams->init_val << ";\n";
        }
    }
  else
    {
      if( SymbolTable::getType(EstimParams->name) == eExogenous)
        {
          *output << "tmp1 = find((estim_params_.corrx(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ")) & (estim_params_.corrx(:,2)==SymbolTable::getID(EstimParams->name2)+1);\n";
          *output << "estim_params_.corrx(tmp1,3) = " << EstimParams->init_val << ";\n";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eEndogenous)
        {
          *output << "tmp1 = find((estim_params_.corrn(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ")) & (estim_params_.corrn(:,2)==" << SymbolTable::getID(EstimParams->name2)+1 << ";\n";
          *output << "estim_params_.corrx(tmp1,3) = " << EstimParams->init_val << ";\n";
        }
    }
  EstimParams->clear();
}

//------------------------------------------------------------------------------
void ComputingTasks::setEstimatedBoundsElements(void)
{
  if (!SymbolTable::Exist(EstimParams->name))
    {
      string msg = "Unknown symbol: "+EstimParams->name;
      error(msg.c_str());
    }
  if (SymbolTable::isReferenced(EstimParams->name) == eNotReferenced)
    {
      return;
    }
  if ((EstimParams->init_val).size() == 0)
    {
      EstimParams->init_val = EstimParams->mean;
    }
  if (EstimParams->type < 3)
    {
      if( SymbolTable::getType(EstimParams->name) == eExogenous)
        {
          *output << "tmp1 = find(estim_params_.var_exo(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ");\n";
          *output << "estim_params_.var_exo(tmp1,3) = " << EstimParams->low_bound << ";\n";
          *output << "estim_params_.var_exo(tmp1,4) = " << EstimParams->up_bound << ";\n";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eEndogenous)
        {
          *output << "tmp1 = find(estim_params_.var_endo(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ");\n";
          *output << "estim_params_.var_endo(tmp1,3) = " << EstimParams->low_bound << ";\n";
          *output << "estim_params_.var_endo(tmp1,4) = " << EstimParams->up_bound << ";\n";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eParameter)
        {
          *output << "tmp1 = find(estim_params_.param_vals(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ");\n";
          *output << "estim_params_.param_vals(tmp1,3) = " << EstimParams->low_bound << ";\n";
          *output << "estim_params_.param_vals(tmp1,4) = " << EstimParams->up_bound << ";\n";
        }
    }
  else
    {
      if( SymbolTable::getType(EstimParams->name) == eExogenous)
        {
          *output << "tmp1 = find((estim_params_.corrx(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ")) & (estim_params_.corrx(:,2)==SymbolTable::getID(EstimParams->name2)+1);\n";
          *output << "estim_params_.corrx(tmp1,4) = " << EstimParams->low_bound << ";\n";
          *output << "estim_params_.corrx(tmp1,5) = " << EstimParams->up_bound << ";\n";
        }
      else if ( SymbolTable::getType(EstimParams->name) == eEndogenous)
        {
          *output << "tmp1 = find((estim_params_.corrn(:,1)==" << SymbolTable::getID(EstimParams->name)+1 << ")) & (estim_params_.corrn(:,2)==" << SymbolTable::getID(EstimParams->name2)+1 << ";\n";
          *output << "estim_params_.corrx(tmp1,4) = " << EstimParams->low_bound << ";\n";
          *output << "estim_params_.corrx(tmp1,5) = " << EstimParams->up_bound << ";\n";
        }
    }
  EstimParams->clear();
}

//-----------------------------------------------------------------------
void ComputingTasks::set_trend_element (string name, string expression)
{
  //Testing if symbol exists
  if (!SymbolTable::Exist(name))
    {
      string msg = "Unknown variable: " + name;
      (* error) (msg.c_str());
    }
  Type  type = SymbolTable::getType(name);
  //    int 	id = SymbolTable::getID(name);
  if (type == eEndogenous)
    {
      *output << "tmp1 = strmatch('" << name << "',options_.varobs,'exact');\n";
      *output << "options_.trend_coeffs{tmp1} = '" << expression << "';\n";
    }
  else
    {
      cout << "Error : Non-variable symbol used in TREND_COEFF: " << name << endl;
    }
}

//------------------------------------------------------------------------------
void ComputingTasks::BeginCalibVar(void)
{

  *output << interfaces::comment() << "\n" << interfaces::comment() << "CALIB_VAR \n"
          << interfaces::comment() << "\n";
  for(int i=1;i<4;++i)
    {
      *output << "calib_var_index{" << i << "} = [];\n";
      *output << "calib_targets{" << i << "} = [];\n";
      *output << "calib_weights{" << i << "}=[];\n";
    }
}

//------------------------------------------------------------------------------
void ComputingTasks::setCalibVar(string name, string weight, string expression)
{
  if (!SymbolTable::Exist(name))
    {
      string msg = "calib_var: " + name + " doesn't exist";
      error(msg.c_str());
    }
  int id = SymbolTable::getID(name)+1;
  if (SymbolTable::getType(name) == eEndogenous)
    {
      *output << "calib_var_index{1} = [calib_var_index{1};" <<  id << "," << id << "];\n";
      *output << "calib_weights{1} = [calib_weights{1}; " << weight << "];\n";
      *output << "calib_targets{1} =[calib_targets{1}; " << expression << "];\n";
    }
  else if (SymbolTable::getType(name) == eExogenous)
    {
      *output << "calib_var_index{3} = [calib_var_index{3};" <<  id << "," << id << "];\n";
      *output << "calib_weights{3} = [calib_weights{3}; " << weight << "];\n";
      *output << "calib_targets{3} =[calib_targets{3}; " << expression << "];\n";
    }
  else
    {
      string msg = "calib_var: " + name + "isn't a endogenous or an exogenous variable";
      error(msg.c_str());
    }
}

//------------------------------------------------------------------------------
void ComputingTasks::setCalibVar(string name1, string name2, string weight, string expression)
{
  if (!SymbolTable::Exist(name1))
    {
      string msg = "calib_var: " + name1 + " doesn't exist";
      error(msg.c_str());
    }
  if (!SymbolTable::Exist(name2))
    {
      string msg = "calib_var: " + name2 + " doesn't exist";
      error(msg.c_str());
    }
  if (SymbolTable::getType(name1) != SymbolTable::getType(name2))
    {
      string msg = "calib_var: " + name1 + " and " + name2 + " don't have the same type";
      error(msg.c_str());
    }
  int id1 = SymbolTable::getID(name1)+1;
  int id2 = SymbolTable::getID(name2)+1;
  if (SymbolTable::getType(name1) == eEndogenous)
    {
      *output << "calib_var_index{1} = [calib_var_index{1};" <<  id1 << "," << id2 << "];\n";
      *output << "calib_weights{1} = [calib_weights{1}; " << weight << "];\n";
      *output << "calib_targets{1} =[calib_targets{1}; " << expression << "];\n";
    }
  else if (SymbolTable::getType(name1) == eExogenous)
    {
      *output << "calib_var_index{3} = [calib_var_index{3};" <<  id1 << "," << id2 << "];\n";
      *output << "calib_weights{3} = [calib_weights{3}; " << weight << "];\n";
      *output << "calib_targets{3} =[calib_targets{3}; " << expression << "];\n";
    }
  else
    {
      string msg = "calib_var: " + name1 + " and " + name2 + "aren't endogenous or exogenous variables";
      error(msg.c_str());
    }
}

void ComputingTasks::setCalibAc(string name, string ar, string weight, string expression)
{
  static int max_iar = 3;
  if (!SymbolTable::Exist(name))
    {
      string msg = "calib_var: " + name + " doesn't exist";
      error(msg.c_str());
    }
  int id = SymbolTable::getID(name)+1;
  int iar = atoi(ar.c_str())+3;
  if (iar > max_iar)
    {
      // creates new variables
      for(int i=max_iar+1; i <= iar; ++i)
        {
          *output << "calib_var_index{" << i << "} = [];\n";
          *output << "calib_targets{" << i << "} = [];\n";
          *output << "calib_weights{" << i << "}=[];\n";
        }
      max_iar = iar;
    }
  if (SymbolTable::getType(name) == eEndogenous)
    {
      *output << "calib_var_index{" << iar << "} = [calib_var_index{" << iar << "};" <<  id << "];\n";
      *output << "calib_weights{" << iar << "} = [calib_weights{" << iar << "}; " << weight << "];\n";
      *output << "calib_targets{" << iar << "} =[calib_targets{" << iar << "}; " << expression << "];\n";
    }
  else
    {
      string msg = "calib_var: " + name + "isn't a endogenous variable";
      error(msg.c_str());
    }
}

//------------------------------------------------------------------------------
void ComputingTasks::runCalib(int cova)
{
  *output << "M_.Sigma_e=calib(calib_var_index,calib_targets,calib_weights," << cova << ",Sigma_e_);\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setOsrParams(string tmp)
{
  *output << tmp;
}

//------------------------------------------------------------------------------
void ComputingTasks::runOsr(string tmp1)
{
  *output << tmp1;
  *output << "osr(var_list_,osr_params_,obj_var_,optim_weights_);\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setOlrInst(string tmp)
{
  *output << tmp;
}

//------------------------------------------------------------------------------
void ComputingTasks::runOlr(string tmp1)
{
  *output << tmp1;
  *output << "options_.olr = 1;\n";
  *output << "options_.olr_w = optim_weights_;\n";
  *output << "options_.olr_inst = olr_inst_;\n";
  *output << "info = stoch_simul(var_list_);\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::BeginOptimWeights(void)
{
  *output << interfaces::comment() << "OPTIM_WEIGHTS\n\n";
  *output << "optim_weights_ = sparse(M_.endo_nbr,M_.endo_nbr);\n";
  *output << "obj_var_ = [];\n\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setOptimWeights(string name, string exp)
{
  if (!SymbolTable::Exist(name) || SymbolTable::getType(name) != eEndogenous)
    {
      string msg = "optim_weights: " + name + " isn't an endogenous variable";
      error(msg.c_str());
    }
  int id = SymbolTable::getID(name)+1;
  *output <<  "optim_weights_(" << id << "," << id << ") = " << exp << ";\n";
  *output << "obj_var_ = [obj_var_; " << id << "];\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::setOptimWeights(string name1, string name2, string exp)
{
  if (!SymbolTable::Exist(name1) || SymbolTable::getType(name1) != eEndogenous)
    {
      string msg = "optim_weights: " + name1 + " isn't an endogenous variable";
      error(msg.c_str());
    }
  if (!SymbolTable::Exist(name2) || SymbolTable::getType(name2) != eEndogenous)
    {
      string msg = "optim_weights: " + name2 + " isn't an endogenous variable";
      error(msg.c_str());
    }
  int id1 = SymbolTable::getID(name1)+1;
  int id2 = SymbolTable::getID(name2)+1;
  *output <<  "optim_weights_(" << id1 << "," << id2 << ") = " << exp << ";\n";
  *output << "obj_var_ = [obj_var_; " << id1 << " " << id2 << "];\n";
}

//------------------------------------------------------------------------------
void ComputingTasks::runDynasave(string filename, string ext, string varlist)
{
  *output << varlist;
  *output << "dynasave(" << filename;
  if (ext.size() > 0)
    {
      *output << "," << ext;
    }
  *output << ",varlist_);\n";
}

void ComputingTasks::runDynatype(string filename, string ext, string varlist)
{
  *output << varlist;
  *output << "dynatype(" << filename;
  if (ext.size() > 0)
    {
      *output << "," << ext;
    }
  *output << ",varlist_);\n";
}

void ComputingTasks::beginModelComparison(void)
{
  *output << "ModelNames_ = {};\n";
  *output << "ModelPriors_ = {};\n";
}

void ComputingTasks::addMcFilename(string filename, string prior)
{
  *output << "ModelNames_ = { ModelNames_{:} '" << filename << "};\n";
  *output << "ModelPriors_ = { ModelPriors_{:} '" << prior << "};\n";
}

void ComputingTasks::runModelComparison(void)
{
  *output << "model_comparison(ModelNames_,ModelPriors_);\n";
}

/*
  string ComputingTasks::get(void)
  {
  return output.str();
  }
*/
