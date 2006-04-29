/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the parser class methodes.
*/
//------------------------------------------------------------------------------
#include "ModelParameters.h"
#include "SymbolTable.h"
#include "Expression.h"
#include "NumericalInitialization.h"
#include "ModelTree.h"
#include "VariableTable.h"
#include "Shocks.h"
#include "SigmaeInitialization.h"
#include "ComputingTasks.h"
#include "TmpSymbolTable.h"
#include "DynareParser.h"

string dynare::parser::file_name = "";
void dynare::parser::set_file_name(string fname)
{
	file_name = fname;
}

void dynare::parser::setoutput(ostringstream* ostr) 
{
	output = ostr;	
	numerical_initialization.setOutput(ostr);
	shocks.setOutput(ostr);
	sigmae.setOutput(ostr);
	computing_tasks.setOutput(ostr);
}

dynare::Objects* dynare::parser::add_endogenous(Objects* obj, Objects* tex_name)
{
  //cout << "add_endogenous \n";
  
  obj->ID = (NodeID) symbol_table.AddSymbolDeclar(obj->symbol,eEndogenous, tex_name->symbol);
  obj->type = eEndogenous;
  return (obj);
}
dynare::Objects* dynare::parser::add_exogenous(Objects* obj, Objects* tex_name)
{
  obj->ID = (NodeID) symbol_table.AddSymbolDeclar(obj->symbol,eExogenous, tex_name->symbol);
  obj->type = eExogenous;
  return (obj);
}
dynare::Objects* dynare::parser::add_exogenous_det(Objects* obj, Objects* tex_name)                                    
{
  obj->ID = (NodeID) symbol_table.AddSymbolDeclar(obj->symbol,eExogenousDet, tex_name->symbol);
  obj->type = eExogenousDet;
  return (obj);
}
dynare::Objects* dynare::parser::add_parameter(Objects* obj, Objects* tex_name)
{
  obj->ID = (NodeID) symbol_table.AddSymbolDeclar(obj->symbol,eParameter, tex_name->symbol);
  obj->type = eParameter;
  return (obj);
}
dynare::Objects* dynare::parser::add_local_parameter(Objects* obj)
{
  obj->ID = (NodeID) symbol_table.AddSymbolDeclar(obj->symbol,eLocalParameter, obj->symbol);
  obj->type = eLocalParameter;
  NodeID id = model_tree.AddTerminal(obj->symbol);
  return new Objects("", id, eTempResult);
}
dynare::Objects* dynare::parser::add_constant(Objects* obj)
{
	obj->ID = (NodeID) num_constants.AddConstant(obj->symbol);
	obj->type = eNumericalConstant;
	return obj;
}
dynare::Objects* dynare::parser::add_model_constant(Objects* constant)
{
	constant = add_constant(constant);
	NodeID id = model_tree.AddTerminal(constant->ID, eNumericalConstant);
	return new Objects("", id, eTempResult);
}
dynare::Objects* dynare::parser::add_variable(Objects* var)
{
	//cout << "add_variable1 : " << var->symbol << endl;
	var = get_symbol(var);
	if((var->type == eEndogenous) 
	   || (var->type == eExogenous)
	   || (var->type == eExogenousDet))
		variable_table.AddVariable(var->symbol,0);
	//cout   << "add_model_token : " << var->ID << endl;
	NodeID id = model_tree.AddTerminal(var->symbol);
	return new Objects("", id, eTempResult);
}
dynare::Objects* dynare::parser::add_variable(Objects* var,Objects* olag)
{
	//cout << "add_variable2\n";
	
	var = get_symbol(var);
	int lag = atoi((olag->symbol).c_str());
	//cout << "symbol = " << olag->symbol << endl;
	//cout << "lag = " << lag << endl;
	if ((var->type == eEndogenous) || (var->type == eExogenous))
		variable_table.AddVariable(var->symbol,lag);
	//cout   << "add_model_token : " << var->ID << endl;
	NodeID id = model_tree.AddTerminal(var->symbol,lag);
	return new Objects("", id, eTempResult);
}
dynare::Objects* dynare::parser::get_symbol(Objects* obj)
{
	if (!symbol_table.Exist(obj->symbol))
	{
		string msg = "Unknown symbol : "+obj->symbol;
		error(msg.c_str());
	}
	obj->ID = (NodeID) symbol_table.getID(obj->symbol);
	obj->type = symbol_table.getType(obj->symbol);
	return obj;
}
dynare::Objects* dynare::parser::translate_symbol(Objects* obj)
{
	if (!symbol_table.Exist(obj->symbol))
	{
		string msg = "Unknown symbol : "+obj->symbol;
		error(msg.c_str());
	}
	obj->ID = (NodeID) symbol_table.getID(obj->symbol);
	obj->type = symbol_table.getType(obj->symbol);
	ostringstream symbol;
	if (obj->type == eEndogenous)
	  {
	    symbol << "oo_.steady_state(" << (int)obj->ID+1 << ")";
	    obj->symbol = symbol.str();
	  }
	else if (obj->type == eExogenous)
	  {
	    symbol << "oo_.exo_steady_state(" << (int)obj->ID+1 << ")";
	    obj->symbol = symbol.str();
	  }
	else if (obj->type == eExogenousDet)
	  {
	    symbol << "oo_.exo_det_steady_state(" << (int)obj->ID+1 << ")";
	    obj->symbol = symbol.str();
	  }
	else if (obj->type == eParameter)
	  {
	    symbol << "M_.params(" << (int)obj->ID+1 << ")";
	    obj->symbol = symbol.str();
	  }
	else if (obj->type == eLocalParameter)
	  {
	    symbol << obj->symbol;
	  }
	return obj;
}

dynare::Objects* dynare::parser::add_expression_token( Objects* arg1,  Objects* arg2,  Objects* op)
{
	int id = expression.AddToken((int) arg1->ID,arg1->type,
						(int) arg2->ID,arg2->type,
						op->opcode);
	//cout << "after add_expression_token\n";
	return new Objects("", (NodeID) id, eTempResult);
}
dynare::Objects* dynare::parser::add_expression_token( Objects* arg1, Objects* op)
{
  int id;
  if (op->opcode != NAME)
    {
      id = expression.AddToken((int) arg1->ID,arg1->type,
						op->opcode);
    }
  else
    {
      id = expression.AddToken((int) arg1->ID,arg1->type,
						op->symbol);
    }

  //cout << "after add_expression_token\n";
  return new Objects("", (NodeID) id, eTempResult);
}
dynare::Objects* dynare::parser::get_expression(Objects* exp) 
{
	if (exp->type == eTempResult)
	{
		expression.set();
		string exp = expression.get();
		expression.clear();
		return new Objects(exp);
	}
	else
		return exp;
}
dynare::Objects* dynare::parser::cat(Objects* string1, Objects* string2)
{
  dynare::Objects* result = new dynare::Objects;
  result->symbol = string1->symbol+string2->symbol;
  return result;
}
dynare::Objects* dynare::parser::cat_with_space(Objects* string1, Objects* string2)
{
  dynare::Objects* result = new dynare::Objects;
  result->symbol = string1->symbol+" "+string2->symbol;
  return result;
}
void dynare::parser::init_param(Objects* lhs,  Objects* rhs)
{
	numerical_initialization.SetConstant(lhs->symbol, rhs->symbol);
}
void dynare::parser::init_param(Objects* lhs)
{
	//cout << "Befor set\n";
	expression.set();   	
	numerical_initialization.SetConstant(lhs->symbol, expression.get());
	expression.clear();
}
void dynare::parser::init_val(Objects* lhs,  Objects* rhs)
{
	numerical_initialization.SetInit(lhs->symbol, rhs->symbol);
}
void dynare::parser::hist_val(Objects* lhs, Objects* slag, Objects* rhs)
{
	int lag = atoi((slag->symbol).c_str());
	numerical_initialization.SetHist(lhs->symbol, lag, rhs->symbol);
}
void dynare::parser::hist_val(Objects* lhs, Objects* slag)
{
	int lag = atoi((slag->symbol).c_str());
	expression.set();  	
	numerical_initialization.SetHist(lhs->symbol, lag, expression.get());
	expression.clear();
}
void dynare::parser::initialize_model(void)
{
}
void dynare::parser::use_dll(void)
{
	// Seetting variable momber offset to use C outputs
	model_tree.offset = 0;
}
void dynare::parser::check_model(void)
{
	symbol_table.clean();
}
void dynare::parser::finish(void)
{
	
  string model_file_name(file_name);   

  // Setting flags to compute what is necessary	
  if (order == 1 || linear == 1)
    {
      model_tree.computeJacobianExo = true;
      model_tree.computeJacobian = false;
    }
  else if (order != -1 && linear != -1)
    {
      model_tree.computeHessian = true;
      model_tree.computeJacobianExo = true;
    }
  // Removing extension chars
  model_file_name.erase(model_file_name.size()-4,4);
  model_tree.ModelInitialization();

  if ( model_tree.computeHessian )
    {
      model_tree.derive(2);
    }
  else
    {
      model_tree.derive(1);
    }
	
  cout << "Processing outputs ...\n";
  model_tree.setStaticModel();
  model_tree.setDynamicModel();

  if (model_tree.offset == 0)
    {
      model_tree.OpenCFiles(model_file_name+"_static", model_file_name+"_dynamic");
      model_tree.SaveCFiles();
    }
  else
    {
      model_tree.OpenMFiles(model_file_name+"_static", model_file_name+"_dynamic");
      model_tree.SaveMFiles();
    }
  
  *output << "save('" << model_file_name << "_results', 'oo_');\n";
  *output << "diary off\n";

  symbol_table.erase_local_parameters();
}
void dynare::parser::begin_initval(void)
{
	numerical_initialization.BeginInitval();
}
void dynare::parser::end_initval(void)
{
	numerical_initialization.EndInitval();
}
void dynare::parser::begin_endval(void)
{
	numerical_initialization.BeginEndval();
}
void dynare::parser::end_endval(void)
{
	numerical_initialization.EndEndval();
}
void dynare::parser::begin_histval(void)
{
	numerical_initialization.BeginHistval();
}
void dynare::parser::begin_shocks(void)
{
	shocks.BeginShocks();
}
void dynare::parser::begin_mshocks(void)
{
	shocks.BeginMShocks();
}
void dynare::parser::end_shocks(void)
{
  shocks.EndShocks();
}
void dynare::parser::add_det_shock(Objects* var)
{
  if (!symbol_table.Exist(var->symbol))
    {
      string msg = "Unknown symbol : "+var->symbol;
      error(msg.c_str());
    }
  int id = symbol_table.getID(var->symbol);
  switch (symbol_table.getType(var->symbol))
    {
    case eExogenous:
      shocks.AddDetShockExo(id);
      return;
    case eExogenousDet:
      shocks.AddDetShockExoDet(id);
      return;
    default:	    
      error("Shocks can only be applied to exogenous variables");
    }
}
void dynare::parser::add_stderr_shock(Objects* var, Objects* value)
{
	if (!symbol_table.Exist(var->symbol))
	{
		string msg = "Unknown symbol : "+var->symbol;
		error(msg.c_str());
	}
	int id = symbol_table.getID(var->symbol);
	shocks.AddSTDShock(id, value->symbol);
}
void dynare::parser::add_var_shock(Objects* var, Objects* value)
{
	if (!symbol_table.Exist(var->symbol))
	{
		string msg = "Unknown symbol : "+var->symbol;
		error(msg.c_str());
	}
	int id = symbol_table.getID(var->symbol);
	shocks.AddVARShock(id, value->symbol);
}
void dynare::parser::add_covar_shock(Objects* var1, Objects* var2, Objects* value)
{
	if (!symbol_table.Exist(var1->symbol))
	{
		string msg = "Unknown symbol : "+var1->symbol;
		error(msg.c_str());
	}
	if (!symbol_table.Exist(var2->symbol))
	{
		string msg = "Unknown symbol : "+var2->symbol;
		error(msg.c_str());
	}
	int id1 = symbol_table.getID(var1->symbol);
	int id2 = symbol_table.getID(var2->symbol);
	shocks.AddCOVAShock(id1, id2, value->symbol);
}
void dynare::parser::add_correl_shock(Objects* var1, Objects* var2, Objects* value)
{
	if (!symbol_table.Exist(var1->symbol))
	{
		string msg = "Unknown symbol : "+var1->symbol;
		error(msg.c_str());
	}
	if (!symbol_table.Exist(var2->symbol))
	{
		string msg = "Unknown symbol : "+var2->symbol;
		error(msg.c_str());
	}
	int id1 = symbol_table.getID(var1->symbol);
	int id2 = symbol_table.getID(var2->symbol);
	shocks.AddCORRShock(id1, id2, value->symbol);
}
void dynare::parser::add_period(Objects* p1, Objects* p2)
{
	shocks.AddPeriod(p1->symbol, p2->symbol);
}
void dynare::parser::add_period(Objects* p1)
{
	shocks.AddPeriod(p1->symbol, p1->symbol);
}
void dynare::parser::add_value(Objects* value)
{
	shocks.AddValue(value->symbol);
}
void dynare::parser::do_sigma_e(void)
{
	sigmae.set();
}
void dynare::parser::end_of_row(void)
{
	sigmae.EndOfRow();
}
void dynare::parser::add_to_row(Objects* s)
{
	sigmae.AddExpression(s->symbol);
}
void dynare::parser::steady(void)
{
	computing_tasks.setSteady();
	model_tree.computeJacobian = true;
}
void dynare::parser::option_num(string name_option, Objects* opt)
{
	computing_tasks.setOption(name_option, opt->symbol);
	if (name_option == "order")
		order = atoi((opt->symbol).c_str());
	else if (name_option == "linear")
		linear = atoi((opt->symbol).c_str());
}
void dynare::parser::option_num(string name_option, Objects* opt1, Objects* opt2)
{
  computing_tasks.setOption(name_option, opt1->symbol, opt2->symbol);
}
void dynare::parser::option_num(string name_option, string opt)
{
	computing_tasks.setOption(name_option, opt);
	if (name_option == "order")
		order = atoi(opt.c_str());
	else if (name_option == "linear")
		linear = atoi(opt.c_str());
}
void dynare::parser::option_str(string name_option, Objects* opt)
{
	opt->symbol = "'"+opt->symbol;
	opt->symbol += "'";
	computing_tasks.setOption(name_option, opt->symbol);
}
void dynare::parser::option_str(string name_option, string opt)
{
	opt = "'"+opt;
	opt += "'";
	computing_tasks.setOption(name_option, opt);
}
void dynare::parser::add_tmp_var(Objects* tmp_var1, Objects* tmp_var2)
{
	tmp_symbol_table.AddTempSymbol(tmp_var1->symbol, tmp_var2->symbol);
}
void dynare::parser::add_tmp_var(Objects* tmp_var)
{
	tmp_symbol_table.AddTempSymbol(tmp_var->symbol);
}
// dynare::Objects* get_tmp_var(string)
// {
// 	//string str = tmp_symbol_table.get
// }
void dynare::parser::rplot()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runRplot(tmp);
}
void dynare::parser::stoch_simul()
{
	// If order and linear not set, then set them to default values
	if (order == -1)
	{
		order = 2;
	}
	if (linear == -1)
	{
		linear = 0;
	}

	tmp_symbol_table.set("var_list_");
	string tmp = tmp_symbol_table.get();
	computing_tasks.setStochSimul(tmp);
}
void dynare::parser::simul()
{
	computing_tasks.setSimul();
	model_tree.computeJacobian = true;
}
void dynare::parser::check()
{
	computing_tasks.setCheck();
	model_tree.computeJacobian = true;
}
void dynare::parser::estimation_init()
{
	computing_tasks.EstimParams = &estim_params;
	computing_tasks.setEstimationInit();
	model_tree.computeJacobianExo = true;
}
void dynare::parser::set_estimated_elements(void)
{
	computing_tasks.setEstimatedElements();
}
void dynare::parser::set_estimated_init_elements(void)
{
  computing_tasks.setEstimatedInitElements();
}
void dynare::parser::set_estimated_bounds_elements(void)
{
  computing_tasks.setEstimatedBoundsElements();
}
void dynare::parser::set_unit_root_vars()
{
  tmp_symbol_table.set("options_.unit_root_vars");
  *output << tmp_symbol_table.get();
}
void dynare::parser::run_estimation()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runEstimation(tmp);
}
void dynare::parser::optim_options(Objects* str1, Objects* str2, int task)
{
	computing_tasks.setOptimOptions(str1->symbol, str2->symbol, task);
}
void dynare::parser::optim_options(int task)
{
	computing_tasks.setOptimOptions("", "", task);
}
void dynare::parser::set_varobs()
{
  tmp_symbol_table.set("options_.varobs");
  *output << tmp_symbol_table.get();
}
void dynare::parser::set_trend_init()
{
  *output << "options_.trend_coeff_ = {};\n";
}
void dynare::parser::set_trend_element(Objects* arg1, Objects* arg2)
{
  computing_tasks.set_trend_element(arg1->symbol, arg2->symbol);
}
void dynare::parser::begin_optim_weights(void)
{
  computing_tasks.BeginOptimWeights();
}
void dynare::parser::set_optim_weights(Objects* arg1, Objects* arg2)
{
  computing_tasks.setOptimWeights(arg1->symbol, arg2->symbol);
}
void dynare::parser::set_optim_weights(Objects* arg1, Objects* arg2, Objects* arg3)
{
  computing_tasks.setOptimWeights(arg1->symbol, arg2->symbol, arg3->symbol);
}
void dynare::parser::set_osr_params(void)
{
  tmp_symbol_table.set("osr_params_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.setOsrParams(tmp);
}
void dynare::parser::run_osr(void)
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runOsr(tmp);
}
void dynare::parser::set_olr_inst(void)
{
  tmp_symbol_table.set("options_.olr_inst");
  string tmp = tmp_symbol_table.get();
  computing_tasks.setOlrInst(tmp);
}
void dynare::parser::run_olr(void)
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runOlr(tmp);
}
void dynare::parser::begin_calib_var(void)
{
  computing_tasks.BeginCalibVar();
}
void dynare::parser::set_calib_var(Objects* name, Objects* weight, Objects* expression)
{
  Objects* exp = get_expression(expression);
  computing_tasks.setCalibVar(name->symbol,weight->symbol,exp->symbol);
}
void dynare::parser::set_calib_var(Objects* name1, Objects* name2, Objects* weight, Objects* expression)
{
  Objects* exp = get_expression(expression);
  computing_tasks.setCalibVar(name1->symbol,name2->symbol,weight->symbol,exp->symbol);
}
void dynare::parser::set_calib_ac(Objects* name, Objects* ar, Objects* weight, Objects* expression)
{
  Objects* exp = get_expression(expression);
  computing_tasks.setCalibAc(name->symbol,ar->symbol,weight->symbol,exp->symbol);
}
void dynare::parser::run_calib(int flag)
{
  computing_tasks.runCalib(flag);
}
void dynare::parser::run_dynatype(Objects* filename, Objects* ext)
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runDynatype(filename->symbol,ext->symbol,tmp);
}
void dynare::parser::run_dynasave(Objects* filename, Objects* ext)
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runDynasave(filename->symbol,ext->symbol,tmp);
}
void dynare::parser::begin_model_comparison(void)
{
  computing_tasks.beginModelComparison();
}
void dynare::parser::add_mc_filename(Objects* filename, Objects* prior)
{
  computing_tasks.addMcFilename(filename->symbol, prior->symbol);
}
void dynare::parser::run_model_comparison(void)
{
  computing_tasks.runModelComparison();
}
dynare::Objects*	dynare::parser::add_equal(Objects* arg1,  Objects* arg2)
{
	NodeID id = model_tree.AddEqual(arg1->ID, arg2->ID);
	model_parameters.eq_nbr++;
	return new Objects("", id, eTempResult);
}
		
dynare::Objects*	dynare::parser::init_local_parameter(Objects* arg1,  Objects* arg2)
{
	NodeID id = model_tree.AddEqual(arg1->ID, arg2->ID);
	return new Objects("", id, eTempResult);
}
		
dynare::Objects*	dynare::parser::add_plus(Objects* arg1,  Objects* arg2)
{
	NodeID id = model_tree.AddPlus(arg1->ID, arg2->ID);
	return new Objects("", id, eTempResult);
}

dynare::Objects*	dynare::parser::add_minus(Objects* arg1,  Objects* arg2)
{
	NodeID id = model_tree.AddMinus(arg1->ID, arg2->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_uminus(Objects* arg1)
{
	NodeID id = model_tree.AddUMinus(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_times(Objects* arg1,  Objects* arg2)
{
	NodeID id = model_tree.AddTimes(arg1->ID, arg2->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_divide(Objects* arg1,  Objects* arg2)
{
	NodeID id = model_tree.AddDivide(arg1->ID, arg2->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_power(Objects* arg1,  Objects* arg2)
{
	NodeID id = model_tree.AddPower(arg1->ID, arg2->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_exp(Objects* arg1)
{
	NodeID id = model_tree.AddExp(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_log(Objects* arg1)
{
	NodeID id = model_tree.AddLog(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_log10(Objects* arg1)
{
	NodeID id = model_tree.AddLog10(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_cos(Objects* arg1)
{
	NodeID id = model_tree.AddCos(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_sin(Objects* arg1)
{
	NodeID id = model_tree.AddSin(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_tan(Objects* arg1)
{
	NodeID id = model_tree.AddTan(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_acos(Objects* arg1)
{
	NodeID id = model_tree.AddACos(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_asin(Objects* arg1)
{
	NodeID id = model_tree.AddASin(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_atan(Objects* arg1)
{
	NodeID id = model_tree.AddATan(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_cosh(Objects* arg1)
{
	NodeID id = model_tree.AddCosH(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_sinh(Objects* arg1)
{
	NodeID id = model_tree.AddSinH(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_tanh(Objects* arg1)
{
	NodeID id = model_tree.AddTanH(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_acosh(Objects* arg1)
{
	NodeID id = model_tree.AddACosH(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_asinh(Objects* arg1)
{
	NodeID id = model_tree.AddASinH(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_atanh(Objects* arg1)
{
	NodeID id = model_tree.AddATanH(arg1->ID);
	return new Objects("", id, eTempResult);
}
dynare::Objects*	dynare::parser::add_sqrt(Objects* arg1)
{
	NodeID id = model_tree.AddSqRt(arg1->ID);
	return new Objects("", id, eTempResult);
}
