/** \file 
 * \version 1.0
 * \date 12/16/2003
 * \par This file implements the ComputingTasks class methodes.
 */
 #include <iostream>
#include <sstream>

using namespace std;
//------------------------------------------------------------------------------
#include "ComputingTasks.h"
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
void ComputingTasks::runEstimation(string tmp1)			
{
  *output << tmp1;
  *output << "dynare_estimation(var_list_);\n";
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
}
//------------------------------------------------------------------------------
void ComputingTasks::setOptimOptions(string str1, string str2, int task)
{
  static string optim_string;
  static int start;
  switch(task){
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
		string msg = "Unknown symbol : "+EstimParams->name;
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
	if (EstimParams->type == 2)
	{
			// Setting user parameter name
			*output << "estim_params_.user_param_names = strvcat(estim_params_.user_param_names, ";
			*output << "'" << EstimParams->name << "');\n";
			// Setting Dynare parameter name
			*output << "estim_params_.param_names = strvcat(estim_params_.param_names, ";
			*output << "'M_.params(" << SymbolTable::getID(EstimParams->name)+1 << ")');\n";
	}
	if( SymbolTable::getType(EstimParams->name) == eExogenous)
	{
		*output << "estim_params_.var_exo = [estim_params_.var_exo; ";
	}
	else if ( SymbolTable::getType(EstimParams->name) == eEndogenous)
	{
		*output << "estim_params_.var_endo = [estim_params_.var_endo; ";
	}
	else if ( SymbolTable::getType(EstimParams->name) == eParameter)
	{
		*output << "estim_params_.param_vals = [estim_params_.param_vals; ";
	}
	if (EstimParams->type == 1)
	{
		*output << SymbolTable::getID(EstimParams->name)+1 << " " << 
			EstimParams->init_val << " " <<  EstimParams->low_bound << " " <<
			EstimParams->up_bound << " " <<  EstimParams->prior << " ";
		*output <<  EstimParams->mean << " " <<  EstimParams->std << " " <<
			EstimParams->p3 << " " <<  EstimParams->p4  << " " <<  EstimParams->jscale << "];\n";	
	}
	else if (EstimParams->type == 2)
	{
		*output << EstimParams->init_val << " " <<  EstimParams->low_bound << " " <<
			EstimParams->up_bound << " " <<  EstimParams->prior << " ";
		*output <<  EstimParams->mean << " " <<  EstimParams->std << " " <<  
			EstimParams->p3 << " " <<  EstimParams->p4  << " " <<  EstimParams->jscale << "];\n";
	}
	EstimParams->clear();
}
void ComputingTasks::set_trend_element (string name, string expression)
{
  //Testing if symbol exists 
  if (!SymbolTable::Exist(name))
    {
      string msg = "Unknown variable: " + name;
      (* error) (msg.c_str());
    }
    Type 	type = SymbolTable::getType(name);
    int 	id = SymbolTable::getID(name);
    if (type == eEndogenous)
      {
	*output << "tmp1 = strmatch(" << name << ",options_.varobs,'exact');\n";
	*output << "options_.trend_coeffs{tmp1} = " << expression << ";\n";
      }
    else
      {
    	cout << "Error : Non-variable symbol used in TREND_COEFF: " << name << endl;
      }
}

//------------------------------------------------------------------------------
void ComputingTasks::setCalibInit(void)
{
	  
	*output << "%\n% CALIB_VAR \n%\n";
	for(int i=1;i<4;++i)
	{
		*output << "calib_var_index{" << i << "} = [];\n";
    	*output << "calib_targets{" << i << "} = [];\n";
		*output << "calib_weights{" << i << "}=[];\n";
	}
}
//------------------------------------------------------------------------------
void ComputingTasks::setCalibVariance(void)
{
/*	
  char buffer[200];
  if (p_t->endo_exo == 1 || p_t->endo_exo == 3)
    {
      sprintf(buffer,"calib_var_index{1} = [calib_var_index{1};%d];\n",p_t->nbr+1);
      str_output(buffer);
      str_output("calib_targets{1} =[calib_targets{1}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{1} = [calib_weights{1}; %s];\n",weight);
      str_output(buffer);
    }
  else if (p_t->endo_exo == 0)
    {
      sprintf(buffer,"calib_var_index{3} = [calib_var_index{3};%d %d];\n",p_t->nbr+1,p_t->nbr+1);
      str_output(buffer);
      str_output("calib_targets{3} =[calib_targets{3}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{3} = [calib_weights{3}; %s];\n",weight);
      str_output(buffer);
    }
  else
    {
      printf("ERROR in CALIB: one of the targets isn't legitimate\n");
    }	
*/
}
//------------------------------------------------------------------------------
void ComputingTasks::setCalibCovariance(void)
{
/*
0 : exo
1 endo
  char buffer[200];
  if ((p_t1->endo_exo == 0 && p_t2->endo_exo == 1)|| (p_t1->endo_exo == 1 && p_t2->endo_exo == 0))
    {
      printf("ERROR in CALIB: can't target correlation between an Endogenousous and an Exogenousous variable\n");
      exit(1);
    }
  else if (p_t1->endo_exo == 1 || p_t1->endo_exo == 3)
    {
      sprintf(buffer,"calib_var_index{2} = [calib_var_index{2};%d %d];\n",p_t1->nbr+1,p_t2->nbr+1);
      str_output(buffer);
      str_output("calib_targets{2} =[calib_targets{2}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{2} = [calib_weights{2}; %s];\n",weight);
      str_output(buffer);
    }
  else if (p_t1->endo_exo == 0)
    {
      sprintf(buffer,"calib_var_index{3} = [calib_var_index{2};%d %d];\n",p_t1->nbr+1,p_t2->nbr+1);
      str_output(buffer);
      str_output("calib_targets{3} =[calib_targets{3}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{3} = [calib_weights{3}; %s];\n",weight);
      str_output(buffer);
    }
  else
    {
      printf("ERROR in CALIB: one of the targets isn''t legitimate\n");
      exit(1);
    }
*/
}
//------------------------------------------------------------------------------
void ComputingTasks::setCalibAutoCorrelation(void)
{
/*
  char buffer[200];
  int i, iar;

  iar = atoi(ar)+3;
  if (iar > max_iar)
    {
      for(i=max_iar+1; i <= iar; ++i)
	{
	  sprintf(buffer,"calib_var_index{%d} = [];\ncalib_targets{%d} = [];\ncalib_weights{%d}=[];\n",i,i,i);
	  str_output(buffer);
	}
      max_iar = iar;
    }
  sprintf(buffer,"calib_var_index{%d} = [calib_var_index{%d};%d];\n",iar,iar,p_t->nbr+1);
  str_output(buffer);
  sprintf(buffer,"calib_targets{%d} =[calib_targets{%d}; ",iar,iar);
  str_output(buffer);
  p_expression(p_q);
  str_output("];\n");
  sprintf(buffer,"calib_weights{%d} = [calib_weights{%d}; %s];\n",iar,iar,weight);
  str_output(buffer);
*/
}
//------------------------------------------------------------------------------
void ComputingTasks::setCalib(void)
{
/*	
	sprintf(buffer,"M_.Sigma_e=calib(calib_var_index,calib_targets,calib_weights,%d,%d,Sigma_e_);\n",max_iar-3,cova);
  	str_output(buffer);
*/
}
//------------------------------------------------------------------------------
void ComputingTasks::setOsr(string tmp1)
{
	*output << tmp1;
	*output << "osr(var1_list_,osr_params_,optim_weights_);\n";
}
//------------------------------------------------------------------------------
void ComputingTasks::setOlr(string tmp1, string tmp2)
{
	*output << tmp1 << tmp2;
	*output << "olr(var_list_,olr_inst_,obj_var_,optim_weights_);\n";
	
}
//------------------------------------------------------------------------------
void ComputingTasks::setOptimWeightsInit(void)
{
	*output << "% OPTIM_WEIGHTS\n\n";
	*output << "optim_weights_ = sparse(endo_nbr,endo_nbr);\n";
	*output << "obj_var_ = [];\n\n";
}
//------------------------------------------------------------------------------
void ComputingTasks::setOptimWeights1(void)
{
/*	
  char buffer[200];

  if (p_t->endo_exo != 1 && p_t->endo_exo != 3)
    {
      fprintf(stdout,"OPTIM_WEIGHTS ERROR: only Endogenousous variables can have weights" );
      exit(1);
    }

  sprintf(buffer,"optim_weights_(%d,%d) = ",p_t->nbr+1,p_t->nbr+1);
  str_output(buffer);
  p_expression(p_q);
  str_output(";\n");
  sprintf(buffer,"obj_var_ = [obj_var_; %d];\n",p_t->nbr+1);
  str_output(buffer);
*/  
}
//------------------------------------------------------------------------------
void ComputingTasks::setOptimWeights2(void)
{
/*	
   char buffer[200];

  if ((p_t1->endo_exo != 1 && p_t1->endo_exo != 3) || (p_t2->endo_exo != 1 && p_t2->endo_exo != 3))
    {
      fprintf(stdout,"OPTIM_WEIGHTS ERROR: only Endogenousous variables can have weights" );
      exit(1);
    }

  sprintf(buffer,"optim_weights_(%d,%d) = ",p_t1->nbr+1,p_t2->nbr+1);
  str_output(buffer);
  p_expression(p_q);
  str_output(";\n");
  sprintf(buffer,"optim_weights_(%d,%d) = optim_weights_(%d,%d);\n",p_t2->nbr+1,p_t1->nbr+1,p_t1->nbr+1,p_t2->nbr+1);
  str_output(buffer);
  sprintf(buffer,"obj_var_ = [obj_var_; %d];\n",p_t1->nbr+1);
  sprintf(buffer,"obj_var_ = [obj_var_; %d];\n",p_t2->nbr+1);
*/  
}
//------------------------------------------------------------------------------
/*
string ComputingTasks::get(void)
{
	return output.str();
}
*/
