#include "ParsingDriver.hh"

ParsingDriver::ParsingDriver() : trace_scanning(false), trace_parsing(false)
{
  order = -1;
  linear = -1;
  model_tree.error = error;
  symbol_table.error = error;
  variable_table.error = error;
  shocks.error = error;
  numerical_initialization.error = error;
  computing_tasks.error = error;
  tmp_symbol_table.error = error;
}

ParsingDriver::~ParsingDriver()
{
}

void
ParsingDriver::check_symbol_existence(const string &name)
{
  if (!symbol_table.Exist(name))
    error("Unknown symbol: " + name);
}

string
ParsingDriver::get_expression(ExpObj *exp)
{
  // Here we don't call "delete exp", since this will be done by the calling function
  if (exp->second == eTempResult)
    {
      expression.set();
      string sexp = expression.get();
      expression.clear();
      return sexp;
    }
  else
    return Expression::getArgument(exp->second, exp->first);
}

void
ParsingDriver::parse(const string &f)
{
  file = f;
  scan_begin();
  yy::parser parser(*this);
  parser.set_debug_level(trace_parsing);
  parser.parse();
  scan_end();
}

void
ParsingDriver::error(const yy::parser::location_type &l, const string &m)
{
  cerr << l << ": " << m << endl;
  exit(-1);
}

void
ParsingDriver::error(const string &m)
{
  extern int yylineno;
  cerr << file << ":" << yylineno << ": " << m << endl;
  exit(-1);
}

void
ParsingDriver::setoutput(ostringstream *ostr)
{
  output = ostr;
  numerical_initialization.setOutput(ostr);
  shocks.setOutput(ostr);
  sigmae.setOutput(ostr);
  computing_tasks.setOutput(ostr);
}

void
ParsingDriver::declare_endogenous(string *name, string *tex_name)
{
  symbol_table.AddSymbolDeclar(*name, eEndogenous, *tex_name);
  delete name;
  delete tex_name;
}

void
ParsingDriver::declare_exogenous(string *name, string *tex_name)
{
  symbol_table.AddSymbolDeclar(*name, eExogenous, *tex_name);
  delete name;
  delete tex_name;
}

void
ParsingDriver::declare_exogenous_det(string *name, string *tex_name)
{
  symbol_table.AddSymbolDeclar(*name, eExogenousDet, *tex_name);
  delete name;
  delete tex_name;
}

void
ParsingDriver::declare_parameter(string *name, string *tex_name)
{
  symbol_table.AddSymbolDeclar(*name, eParameter, *tex_name);
  delete name;
  delete tex_name;
}

ExpObj *
ParsingDriver::add_expression_constant(string *constant)
{
  int id = num_constants.AddConstant(*constant);
  delete constant;
  return new ExpObj(id, eNumericalConstant);
}

NodeID
ParsingDriver::add_model_constant(string *constant)
{
  int id = num_constants.AddConstant(*constant);
  delete constant;
  return model_tree.AddTerminal((NodeID) id, eNumericalConstant);
}

NodeID
ParsingDriver::add_model_variable(string *name)
{
  check_symbol_existence(*name);
  Type type = symbol_table.getType(*name);

  if ((type == eEndogenous)
      || (type == eExogenous)
      || (type == eExogenousDet))
    variable_table.AddVariable(*name, 0);

  NodeID id = model_tree.AddTerminal(*name);
  delete name;
  return id;
}

NodeID
ParsingDriver::add_model_variable(string *name, string *olag)
{
  check_symbol_existence(*name);
  Type type = symbol_table.getType(*name);
  int lag = atoi(olag->c_str());

  if ((type == eExogenous) && lag != 0)
    {
      cout << "Warning: exogenous variable "
           << *name
           << " has lag " << lag << endl;
    }

  if ((type == eEndogenous) || (type == eExogenous))
    variable_table.AddVariable(*name, lag);

  NodeID id = model_tree.AddTerminal(*name, lag);
  delete name;
  delete olag;
  return id;
}

ExpObj *
ParsingDriver::add_expression_variable(string *name)
{
  check_symbol_existence(*name);
  int id = symbol_table.getID(*name);
  Type type = symbol_table.getType(*name);
  delete name;
  return new ExpObj(id, type);
}

ExpObj *
ParsingDriver::add_expression_token(ExpObj *arg1, ExpObj *arg2, int op)
{
  int id = expression.AddToken(arg1->first, arg1->second,
                               arg2->first, arg2->second,
                               op);
  delete arg1;
  delete arg2;
  return new ExpObj(id, eTempResult);
}

ExpObj *
ParsingDriver::add_expression_token(ExpObj *arg1, int op)
{
  int id = expression.AddToken(arg1->first, arg1->second, op);
  delete arg1;
  return new ExpObj(id, eTempResult);
}

ExpObj *
ParsingDriver::add_expression_token(ExpObj *arg1, string *op_name)
{
  int id = expression.AddToken(arg1->first, arg1->second, *op_name);
  delete arg1;
  delete op_name;
  return new ExpObj(id, eTempResult);
}

void
ParsingDriver::init_param(string *name, ExpObj *rhs)
{
  numerical_initialization.SetConstant(*name, get_expression(rhs));
  delete name;
  delete rhs;
}

void
ParsingDriver::init_val(string *name, ExpObj *rhs)
{
  numerical_initialization.SetInit(*name, get_expression(rhs));
  delete name;
  delete rhs;
}

void
ParsingDriver::hist_val(string *name, string *lag, ExpObj *rhs)
{
  int ilag = atoi(lag->c_str());
  numerical_initialization.SetHist(*name, ilag, get_expression(rhs));
  delete name;
  delete lag;
  delete rhs;
}

void
ParsingDriver::use_dll()
{
  // Seetting variable momber offset to use C outputs
  model_tree.offset = 0;
}

void
ParsingDriver::check_model()
{
  // creates too many problems MJ 11/12/06
  //  symbol_table.clean();
}

void
ParsingDriver::finish()
{
  string model_file_name(file);

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

  cout << "Processing outputs ..." << endl;
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

  //  symbol_table.erase_local_parameters();
}

void
ParsingDriver::begin_initval()
{
  numerical_initialization.BeginInitval();
}

void
ParsingDriver::end_initval()
{
  numerical_initialization.EndInitval();
}

void
ParsingDriver::begin_endval()
{
  numerical_initialization.BeginEndval();
}

void
ParsingDriver::end_endval()
{
  numerical_initialization.EndEndval();
}

void
ParsingDriver::begin_histval()
{
  numerical_initialization.BeginHistval();
}

void
ParsingDriver::begin_shocks()
{
  shocks.BeginShocks();
}

void
ParsingDriver::begin_mshocks()
{
  shocks.BeginMShocks();
}

void
ParsingDriver::end_shocks()
{
  shocks.EndShocks();
}

void
ParsingDriver::add_det_shock(string *var)
{
  check_symbol_existence(*var);
  int id = symbol_table.getID(*var);
  switch (symbol_table.getType(*var))
    {
    case eExogenous:
      shocks.AddDetShockExo(id);
      break;
    case eExogenousDet:
      shocks.AddDetShockExoDet(id);
      break;
    default:
      error("Shocks can only be applied to exogenous variables");
    }
  delete var;
}

void
ParsingDriver::add_stderr_shock(string *var, ExpObj *value)
{
  check_symbol_existence(*var);
  int id = symbol_table.getID(*var);
  shocks.AddSTDShock(id, get_expression(value));
  delete var;
  delete value;
}

void
ParsingDriver::add_var_shock(string *var, ExpObj *value)
{
  check_symbol_existence(*var);
  int id = symbol_table.getID(*var);
  shocks.AddVARShock(id, get_expression(value));
  delete var;
  delete value;
}

void
ParsingDriver::add_covar_shock(string *var1, string *var2, ExpObj *value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);
  int id1 = symbol_table.getID(*var1);
  int id2 = symbol_table.getID(*var2);
  shocks.AddCOVAShock(id1, id2, get_expression(value));
  delete var1;
  delete var2;
  delete value;
}

void
ParsingDriver::add_correl_shock(string *var1, string *var2, ExpObj *value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);
  int id1 = symbol_table.getID(*var1);
  int id2 = symbol_table.getID(*var2);
  shocks.AddCORRShock(id1, id2, get_expression(value));
  delete var1;
  delete var2;
  delete value;
}

void
ParsingDriver::add_period(string *p1, string *p2)
{
  shocks.AddPeriod(*p1, *p2);
  delete p1;
  delete p2;
}

void
ParsingDriver::add_period(string *p1)
{
  shocks.AddPeriod(*p1, *p1);
  delete p1;
}

void
ParsingDriver::add_value(string *value)
{
  shocks.AddValue(*value);
  delete value;
}

void
ParsingDriver::add_value(ExpObj *value)
{
  shocks.AddValue(get_expression(value));
  delete value;
}

void
ParsingDriver::do_sigma_e()
{
  sigmae.set();
}

void
ParsingDriver::end_of_row()
{
  sigmae.EndOfRow();
}

void
ParsingDriver::add_to_row(string *s)
{
  sigmae.AddExpression(*s);
  delete s;
}

void
ParsingDriver::add_to_row(ExpObj *v)
{
  sigmae.AddExpression(get_expression(v));
  delete v;
}

void
ParsingDriver::steady()
{
  computing_tasks.setSteady();
  model_tree.computeJacobian = true;
}

void
ParsingDriver::option_num(const string &name_option, string *opt1, string *opt2)
{
  computing_tasks.setOption(name_option, *opt1, *opt2);
  delete opt1;
  delete opt2;
}

void
ParsingDriver::option_num(const string &name_option, string *opt)
{
  option_num(name_option, *opt);
  delete opt;
}

void
ParsingDriver::option_num(const string &name_option, const string &opt)
{
  computing_tasks.setOption(name_option, opt);
  if (name_option == "order")
    order = atoi(opt.c_str());
  else if (name_option == "linear")
    linear = atoi(opt.c_str());
}

void
ParsingDriver::option_str(const string &name_option, string *opt)
{
  option_str(name_option, *opt);
  delete opt;
}

void
ParsingDriver::option_str(const string &name_option, const string &opt)
{
  computing_tasks.setOption(name_option, "'" + opt + "'");
}

void
ParsingDriver::add_tmp_var(string *tmp_var1, string *tmp_var2)
{
  tmp_symbol_table.AddTempSymbol(*tmp_var1, *tmp_var2);
  delete tmp_var1;
  delete tmp_var2;
}

void
ParsingDriver::add_tmp_var(string *tmp_var)
{
  tmp_symbol_table.AddTempSymbol(*tmp_var);
  delete tmp_var;
}

void ParsingDriver::rplot()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runRplot(tmp);
}

void ParsingDriver::stoch_simul()
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

void
ParsingDriver::simul()
{
  computing_tasks.setSimul();
  model_tree.computeJacobian = true;
}

void
ParsingDriver::check()
{
  computing_tasks.setCheck();
  model_tree.computeJacobian = true;
}

void
ParsingDriver::estimation_init()
{
  computing_tasks.EstimParams = &estim_params;
  computing_tasks.setEstimationInit();
  model_tree.computeJacobianExo = true;
}

void
ParsingDriver::set_estimated_elements()
{
  computing_tasks.setEstimatedElements();
}

void
ParsingDriver::set_estimated_init_elements()
{
  computing_tasks.setEstimatedInitElements();
}

void
ParsingDriver::set_estimated_bounds_elements()
{
  computing_tasks.setEstimatedBoundsElements();
}

void
ParsingDriver::set_unit_root_vars()
{
  tmp_symbol_table.set("options_.unit_root_vars");
  *output << tmp_symbol_table.get();
}

void
ParsingDriver::run_estimation()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runEstimation(tmp);
}

void
ParsingDriver::run_prior_analysis()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runPriorAnalysis(tmp);
}

void
ParsingDriver::run_posterior_analysis()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runPosteriorAnalysis(tmp);
}

void
ParsingDriver::optim_options(string *str1, string *str2, int task)
{
  computing_tasks.setOptimOptions(*str1, *str2, task);
  delete str1;
  delete str2;
}

void
ParsingDriver::set_varobs()
{
  tmp_symbol_table.set("options_.varobs");
  *output << tmp_symbol_table.get();
}

void
ParsingDriver::set_trend_init()
{
  *output << "options_.trend_coeff_ = {};" << endl;
}

void
ParsingDriver::set_trend_element(string *arg1, ExpObj *arg2)
{
  computing_tasks.set_trend_element(*arg1, get_expression(arg2));
  delete arg1;
  delete arg2;
}

void
ParsingDriver::begin_optim_weights()
{
  computing_tasks.BeginOptimWeights();
}

void
ParsingDriver::set_optim_weights(string *arg1, ExpObj *arg2)
{
  computing_tasks.setOptimWeights(*arg1, get_expression(arg2));
  delete arg1;
  delete arg2;
}

void
ParsingDriver::set_optim_weights(string *arg1, string *arg2, ExpObj *arg3)
{
  computing_tasks.setOptimWeights(*arg1, *arg2, get_expression(arg3));
  delete arg1;
  delete arg2;
  delete arg3;
}

void
ParsingDriver::set_osr_params()
{
  tmp_symbol_table.set("osr_params_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.setOsrParams(tmp);
}

void
ParsingDriver::run_osr()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runOsr(tmp);
  model_tree.computeJacobianExo = true;
}

void
ParsingDriver::set_olr_inst()
{
  tmp_symbol_table.set("options_.olr_inst");
  string tmp = tmp_symbol_table.get();
  computing_tasks.setOlrInst(tmp);
}

void
ParsingDriver::run_olr()
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runOlr(tmp);
}

void
ParsingDriver::begin_calib_var()
{
  computing_tasks.BeginCalibVar();
}

void
ParsingDriver::set_calib_var(string *name, string *weight, ExpObj *expression)
{
  computing_tasks.setCalibVar(*name, *weight, get_expression(expression));
  delete name;
  delete weight;
  delete expression;
}

void
ParsingDriver::set_calib_var(string *name1, string *name2,
                             string *weight, ExpObj *expression)
{
  computing_tasks.setCalibVar(*name1, *name2, *weight, get_expression(expression));
  delete name1;
  delete name2;
  delete weight;
  delete expression;
}

void
ParsingDriver::set_calib_ac(string *name, string *ar,
                            string *weight, ExpObj *expression)
{
  computing_tasks.setCalibAc(*name, *ar, *weight, get_expression(expression));
  delete name;
  delete ar;
  delete weight;
  delete expression;
}

void
ParsingDriver::run_calib(int flag)
{
  computing_tasks.runCalib(flag);
}

void
ParsingDriver::run_dynatype(string *filename, string *ext)
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runDynatype(*filename, *ext, tmp);
  delete filename;
  delete ext;
}

void
ParsingDriver::run_dynasave(string *filename, string *ext)
{
  tmp_symbol_table.set("var_list_");
  string tmp = tmp_symbol_table.get();
  computing_tasks.runDynasave(*filename, *ext, tmp);
  delete filename;
  delete ext;
}

void
ParsingDriver::begin_model_comparison()
{
  computing_tasks.beginModelComparison();
}

void
ParsingDriver::add_mc_filename(string *filename, string *prior)
{
  computing_tasks.addMcFilename(*filename, *prior);
  delete filename;
  delete prior;
}

void
ParsingDriver::run_model_comparison()
{
  computing_tasks.runModelComparison();
}

NodeID
ParsingDriver::add_model_equal(NodeID arg1, NodeID arg2)
{
  NodeID id = model_tree.AddEqual(arg1, arg2);
  model_parameters.eq_nbr++;
  return id;
}

void
ParsingDriver::declare_and_init_local_parameter(string *name, NodeID rhs)
{
  symbol_table.AddSymbolDeclar(*name, eLocalParameter, *name);
  NodeID id = model_tree.AddTerminal(*name);
  model_tree.AddAssign(id, rhs);
  delete name;
}

NodeID
ParsingDriver::add_model_plus(NodeID arg1, NodeID arg2)
{
  return model_tree.AddPlus(arg1, arg2);
}

NodeID
ParsingDriver::add_model_minus(NodeID arg1, NodeID arg2)
{
  return model_tree.AddMinus(arg1, arg2);
}

NodeID
ParsingDriver::add_model_uminus(NodeID arg1)
{
  return model_tree.AddUMinus(arg1);
}

NodeID
ParsingDriver::add_model_times(NodeID arg1, NodeID arg2)
{
  return model_tree.AddTimes(arg1, arg2);
}

NodeID
ParsingDriver::add_model_divide(NodeID arg1, NodeID arg2)
{
  return model_tree.AddDivide(arg1, arg2);
}

NodeID
ParsingDriver::add_model_power(NodeID arg1, NodeID arg2)
{
  return model_tree.AddPower(arg1, arg2);
}

NodeID
ParsingDriver::add_model_exp(NodeID arg1)
{
  return model_tree.AddExp(arg1);
}

NodeID
ParsingDriver::add_model_log(NodeID arg1)
{
  return model_tree.AddLog(arg1);
}

NodeID
ParsingDriver::add_model_log10(NodeID arg1)
{
  return model_tree.AddLog10(arg1);
}

NodeID
ParsingDriver::add_model_cos(NodeID arg1)
{
  return model_tree.AddCos(arg1);
}

NodeID
ParsingDriver::add_model_sin(NodeID arg1)
{
  return model_tree.AddSin(arg1);
}

NodeID
ParsingDriver::add_model_tan(NodeID arg1)
{
  return model_tree.AddTan(arg1);
}

NodeID
ParsingDriver::add_model_acos(NodeID arg1)
{
  return model_tree.AddACos(arg1);
}

NodeID
ParsingDriver::add_model_asin(NodeID arg1)
{
  return model_tree.AddASin(arg1);
}

NodeID
ParsingDriver::add_model_atan(NodeID arg1)
{
  return model_tree.AddATan(arg1);
}

NodeID
ParsingDriver::add_model_cosh(NodeID arg1)
{
  return model_tree.AddCosH(arg1);
}

NodeID
ParsingDriver::add_model_sinh(NodeID arg1)
{
  return model_tree.AddSinH(arg1);
}

NodeID
ParsingDriver::add_model_tanh(NodeID arg1)
{
  return model_tree.AddTanH(arg1);
}

NodeID
ParsingDriver::add_model_acosh(NodeID arg1)
{
  return model_tree.AddACosH(arg1);
}

NodeID
ParsingDriver::add_model_asinh(NodeID arg1)
{
  return model_tree.AddASinH(arg1);
}

NodeID
ParsingDriver::add_model_atanh(NodeID arg1)
{
  return model_tree.AddATanH(arg1);
}

NodeID
ParsingDriver::add_model_sqrt(NodeID arg1)
{
  return model_tree.AddSqRt(arg1);
}
