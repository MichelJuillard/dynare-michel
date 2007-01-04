#include "ParsingDriver.hh"
#include "Statement.hh"

ParsingDriver::ParsingDriver() : trace_scanning(false), trace_parsing(false)
{
}

ParsingDriver::~ParsingDriver()
{
}

bool
ParsingDriver::exists_symbol(const char *s)
{
  return mod_file->symbol_table.Exist(s);
}

void
ParsingDriver::check_symbol_existence(const string &name)
{
  if (!mod_file->symbol_table.Exist(name))
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
    return expression.getArgument(exp->second, exp->first);
}

ModFile *
ParsingDriver::parse(const string &f)
{
  mod_file = new ModFile();

  mod_file->model_tree.error = error;
  mod_file->symbol_table.error = error;

  expression.setNumericalConstants(&mod_file->num_constants);
  tmp_symbol_table = new TmpSymbolTable(mod_file->symbol_table);

  file = f;
  scan_begin();
  yy::parser parser(*this);
  parser.set_debug_level(trace_parsing);
  parser.parse();
  scan_end();

  delete tmp_symbol_table;

  return mod_file;
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
ParsingDriver::declare_endogenous(string *name, string *tex_name)
{
  mod_file->symbol_table.AddSymbolDeclar(*name, eEndogenous, *tex_name);
  delete name;
  delete tex_name;
}

void
ParsingDriver::declare_exogenous(string *name, string *tex_name)
{
  mod_file->symbol_table.AddSymbolDeclar(*name, eExogenous, *tex_name);
  delete name;
  delete tex_name;
}

void
ParsingDriver::declare_exogenous_det(string *name, string *tex_name)
{
  mod_file->symbol_table.AddSymbolDeclar(*name, eExogenousDet, *tex_name);
  delete name;
  delete tex_name;
}

void
ParsingDriver::declare_parameter(string *name, string *tex_name)
{
  mod_file->symbol_table.AddSymbolDeclar(*name, eParameter, *tex_name);
  delete name;
  delete tex_name;
}

ExpObj *
ParsingDriver::add_expression_constant(string *constant)
{
  int id = mod_file->num_constants.AddConstant(*constant);
  delete constant;
  return new ExpObj(id, eNumericalConstant);
}

NodeID
ParsingDriver::add_model_constant(string *constant)
{
  int id = mod_file->num_constants.AddConstant(*constant);
  delete constant;
  return model_tree->AddTerminal((NodeID) id, eNumericalConstant);
}

NodeID
ParsingDriver::add_model_variable(string *name)
{
  check_symbol_existence(*name);
  NodeID id = model_tree->AddTerminal(*name);
  delete name;
  return id;
}

NodeID
ParsingDriver::add_model_variable(string *name, string *olag)
{
  check_symbol_existence(*name);
  Type type = mod_file->symbol_table.getType(*name);
  int lag = atoi(olag->c_str());

  if ((type == eExogenous) && lag != 0)
    {
      cout << "Warning: exogenous variable "
           << *name
           << " has lag " << lag << endl;
    }
  NodeID id = model_tree->AddTerminal(*name, lag);
  delete name;
  delete olag;
  return id;
}

ExpObj *
ParsingDriver::add_expression_variable(string *name)
{
  check_symbol_existence(*name);
  int id = mod_file->symbol_table.getID(*name);
  Type type = mod_file->symbol_table.getType(*name);
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
ParsingDriver::periods(string *periods)
{
  int periods_val = atoi(periods->c_str());
  mod_file->addStatement(new PeriodsStatement(periods_val));
  delete periods;
}

void
ParsingDriver::dsample(string *arg1)
{
  int arg1_val = atoi(arg1->c_str());
  mod_file->addStatement(new DsampleStatement(arg1_val));
  delete arg1;
}

void
ParsingDriver::dsample(string *arg1, string *arg2)
{
  int arg1_val = atoi(arg1->c_str());
  int arg2_val = atoi(arg2->c_str());
  mod_file->addStatement(new DsampleStatement(arg1_val, arg2_val));
  delete arg1;
  delete arg2;
}

void
ParsingDriver::init_param(string *name, ExpObj *rhs)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eParameter)
    error(*name + " is not a parameter");

  mod_file->addStatement(new InitParamStatement(*name, get_expression(rhs), mod_file->symbol_table));
  delete name;
  delete rhs;
}

void
ParsingDriver::init_val(string *name, ExpObj *rhs)
{
  check_symbol_existence(*name);
  Type type = mod_file->symbol_table.getType(*name);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("initval/endval: " + *name + " should be an endogenous or exogenous variable");

  init_values.push_back(pair<string, string>(*name, get_expression(rhs)));

  delete name;
  delete rhs;
}

void
ParsingDriver::hist_val(string *name, string *lag, ExpObj *rhs)
{
  check_symbol_existence(*name);
  Type type = mod_file->symbol_table.getType(*name);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("hist_val: " + *name + " should be an endogenous or exogenous variable");

  int ilag = atoi(lag->c_str());
  pair<string, int> key(*name, ilag);

  if (hist_values.find(key) != hist_values.end())
    error("hist_val: (" + *name + ", " + *lag + ") declared twice");

  hist_values[key] = get_expression(rhs);

  delete name;
  delete lag;
  delete rhs;
}

void
ParsingDriver::use_dll()
{
  // Seetting variable momber offset to use C outputs
  mod_file->model_tree.offset = 0;
}

void
ParsingDriver::end_initval()
{
  mod_file->addStatement(new InitValStatement(init_values, mod_file->symbol_table));
  init_values.clear();
}

void
ParsingDriver::end_endval()
{
  mod_file->addStatement(new EndValStatement(init_values, mod_file->symbol_table));
  init_values.clear();
}

void
ParsingDriver::end_histval()
{
  mod_file->addStatement(new HistValStatement(hist_values, mod_file->symbol_table));
  hist_values.clear();
}

void
ParsingDriver::begin_model()
{
  model_tree = &mod_file->model_tree;
}

void
ParsingDriver::end_shocks()
{
  mod_file->addStatement(new ShocksStatement(det_shocks, var_shocks, std_shocks,
                                             covar_shocks, corr_shocks, mod_file->symbol_table));
  det_shocks.clear();
  var_shocks.clear();
  std_shocks.clear();
  covar_shocks.clear();
  corr_shocks.clear();
}

void
ParsingDriver::end_mshocks()
{
  mod_file->addStatement(new MShocksStatement(det_shocks, var_shocks, std_shocks,
                                              covar_shocks, corr_shocks, mod_file->symbol_table));
  det_shocks.clear();
  var_shocks.clear();
  std_shocks.clear();
  covar_shocks.clear();
  corr_shocks.clear();
}

void
ParsingDriver::add_det_shock(string *var)
{
  check_symbol_existence(*var);
  Type type = mod_file->symbol_table.getType(*var);
  if (type != eExogenous && type != eExogenousDet)
    error("shocks: shocks can only be applied to exogenous variables");

  if (det_shocks.find(*var) != det_shocks.end())
    error("shocks: variable " + *var + " declared twice");

  if (det_shocks_periods.size() != det_shocks_values.size())
    error("shocks: variable " + *var + ": number of periods is different from number of shock values");

  vector<ShocksStatement::DetShockElement> v;

  for(unsigned int i = 0; i < det_shocks_periods.size(); i++)
    {
      ShocksStatement::DetShockElement dse;
      dse.period1 = det_shocks_periods[i].first;
      dse.period2 = det_shocks_periods[i].second;
      dse.value = det_shocks_values[i];
      v.push_back(dse);
    }

  det_shocks[*var] = v;

  det_shocks_periods.clear();
  det_shocks_values.clear();
  delete var;
}

void
ParsingDriver::add_stderr_shock(string *var, ExpObj *value)
{
  check_symbol_existence(*var);
  if (var_shocks.find(*var) != var_shocks.end()
      || std_shocks.find(*var) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + *var + " declared twice");

  std_shocks[*var] = get_expression(value);

  delete var;
  delete value;
}

void
ParsingDriver::add_var_shock(string *var, ExpObj *value)
{
  check_symbol_existence(*var);
  if (var_shocks.find(*var) != var_shocks.end()
      || std_shocks.find(*var) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + *var + " declared twice");

  var_shocks[*var] = get_expression(value);

  delete var;
  delete value;
}

void
ParsingDriver::add_covar_shock(string *var1, string *var2, ExpObj *value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);

  pair<string, string> key(*var1, *var2), key_inv(*var2, *var1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + *var1 + ", "
          + *var2 + ") declared twice");
  
  covar_shocks[key] = get_expression(value);

  delete var1;
  delete var2;
  delete value;
}

void
ParsingDriver::add_correl_shock(string *var1, string *var2, ExpObj *value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);

  pair<string, string> key(*var1, *var2), key_inv(*var2, *var1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + *var1 + ", "
          + *var2 + ") declared twice");
  
  corr_shocks[key] = get_expression(value);

  delete var1;
  delete var2;
  delete value;
}

void
ParsingDriver::add_period(string *p1, string *p2)
{
  int p1_val = atoi(p1->c_str());
  int p2_val = atoi(p2->c_str());
  det_shocks_periods.push_back(pair<int, int>(p1_val, p2_val));
  delete p1;
  delete p2;
}

void
ParsingDriver::add_period(string *p1)
{
  int p1_val = atoi(p1->c_str());
  det_shocks_periods.push_back(pair<int, int>(p1_val, p1_val));
  delete p1;
}

void
ParsingDriver::add_value(string *value)
{
  det_shocks_values.push_back(*value);
  delete value;
}

void
ParsingDriver::add_value(ExpObj *value)
{
  det_shocks_values.push_back(get_expression(value));
  delete value;
}

void
ParsingDriver::do_sigma_e()
{
  try
    {
      mod_file->addStatement(new SigmaeStatement(sigmae_matrix));
    }
  catch(SigmaeStatement::MatrixFormException &e)
    {
      error("Sigma_e: matrix is neither upper triangular nor lower triangular");
    }
  sigmae_matrix.clear();
}

void
ParsingDriver::end_of_row()
{
  sigmae_matrix.push_back(sigmae_row);
  sigmae_row.clear();
}

void
ParsingDriver::add_to_row(string *s)
{
  sigmae_row.push_back(*s);
  delete s;
}

void
ParsingDriver::add_to_row(ExpObj *v)
{
  sigmae_row.push_back(get_expression(v));
  delete v;
}

void
ParsingDriver::steady()
{
  mod_file->addStatement(new SteadyStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::option_num(const string &name_option, string *opt1, string *opt2)
{
  if (options_list.paired_num_options.find(name_option)
      != options_list.paired_num_options.end())
    error("option " + name_option + " declared twice");

  options_list.paired_num_options[name_option] = pair<string, string>(*opt1, *opt2);
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
  if (options_list.num_options.find(name_option)
      != options_list.num_options.end())
    error("option " + name_option + " declared twice");

  options_list.num_options[name_option] = opt;
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
  if (options_list.string_options.find(name_option)
      != options_list.string_options.end())
    error("option " + name_option + " declared twice");

  options_list.string_options[name_option] = opt;
}

void
ParsingDriver::linear()
{
  mod_file->linear = true;
}

void
ParsingDriver::add_tmp_var(string *tmp_var1, string *tmp_var2)
{
  check_symbol_existence(*tmp_var1);
  check_symbol_existence(*tmp_var2);
  tmp_symbol_table->AddTempSymbol(*tmp_var1, *tmp_var2);
  delete tmp_var1;
  delete tmp_var2;
}

void
ParsingDriver::add_tmp_var(string *tmp_var)
{
  check_symbol_existence(*tmp_var);
  tmp_symbol_table->AddTempSymbol(*tmp_var);
  delete tmp_var;
}

void ParsingDriver::rplot()
{
  mod_file->addStatement(new RplotStatement(*tmp_symbol_table, options_list));
  options_list.clear();
  tmp_symbol_table->clear();
}

void ParsingDriver::stoch_simul()
{
  mod_file->addStatement(new StochSimulStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::simul()
{
  mod_file->addStatement(new SimulStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::check()
{
  mod_file->addStatement(new CheckStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::add_estimated_params_element()
{
  check_symbol_existence(estim_params.name);
  if (estim_params.name2.size() > 0)
    check_symbol_existence(estim_params.name2);

  estim_params_list.push_back(estim_params);
  estim_params.clear();
}

void
ParsingDriver::estimated_params()
{
  mod_file->addStatement(new EstimatedParamsStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_init()
{
  mod_file->addStatement(new EstimatedParamsInitStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_bounds()
{
  mod_file->addStatement(new EstimatedParamsBoundsStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::set_unit_root_vars()
{
  mod_file->addStatement(new UnitRootVarsStatement(*tmp_symbol_table));
  tmp_symbol_table->clear();
}

void
ParsingDriver::run_estimation()
{
  mod_file->addStatement(new EstimationStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::run_prior_analysis()
{
  mod_file->addStatement(new PriorAnalysisStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::run_posterior_analysis()
{
  mod_file->addStatement(new PosteriorAnalysisStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::optim_options_helper(const string &name)
{
  if (options_list.string_options.find("optim_opt") == options_list.string_options.end())
    options_list.string_options["optim_opt"] = "";
  else
    options_list.string_options["optim_opt"] += ",";
  options_list.string_options["optim_opt"] += "''" + name + "'',";
}

void
ParsingDriver::optim_options_string(string *name, string *value)
{
  optim_options_helper(*name);
  options_list.string_options["optim_opt"] += "''" + *value + "''";
  delete name;
  delete value;
}

void
ParsingDriver::optim_options_num(string *name, string *value)
{
  optim_options_helper(*name);
  options_list.string_options["optim_opt"] += *value;
  delete name;
  delete value;
}

void
ParsingDriver::set_varobs()
{
  mod_file->addStatement(new VarobsStatement(*tmp_symbol_table));
  tmp_symbol_table->clear();
}

void
ParsingDriver::set_trends()
{
  mod_file->addStatement(new ObservationTrendsStatement(trend_elements, mod_file->symbol_table));
  trend_elements.clear();
}

void
ParsingDriver::set_trend_element(string *arg1, ExpObj *arg2)
{
  check_symbol_existence(*arg1);
  if (trend_elements.find(*arg1) != trend_elements.end())
    error("observation_trends: " + *arg1 + " declared twice");
  trend_elements[*arg1] = get_expression(arg2);
  delete arg1;
  delete arg2;
}

void
ParsingDriver::set_optim_weights(string *name, ExpObj *value)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eEndogenous)
    error("optim_weights: " + *name + " isn't an endogenous variable");
  if (var_weights.find(*name) != var_weights.end())
    error("optim_weights: " + *name + " declared twice");
  var_weights[*name] = get_expression(value);
  delete name;
  delete value;
}

void
ParsingDriver::set_optim_weights(string *name1, string *name2, ExpObj *value)
{
  check_symbol_existence(*name1);
  if (mod_file->symbol_table.getType(*name1) != eEndogenous)
    error("optim_weights: " + *name1 + " isn't an endogenous variable");

  check_symbol_existence(*name2);
  if (mod_file->symbol_table.getType(*name2) != eEndogenous)
    error("optim_weights: " + *name2 + " isn't an endogenous variable");

  pair<string, string> covar_key(*name1, *name2);

  if (covar_weights.find(covar_key) != covar_weights.end())
    error("optim_weights: pair of variables (" + *name1 + ", " + *name2
          + ") declared twice");

  covar_weights[covar_key] = get_expression(value);
  delete name1;
  delete name2;
  delete value;
}

void
ParsingDriver::optim_weights()
{
  mod_file->addStatement(new OptimWeightsStatement(var_weights, covar_weights, mod_file->symbol_table));
  var_weights.clear();
  covar_weights.clear();
}

void
ParsingDriver::set_osr_params()
{
  mod_file->addStatement(new OsrParamsStatement(*tmp_symbol_table));
  tmp_symbol_table->clear();
}

void
ParsingDriver::run_osr()
{
  mod_file->addStatement(new OsrStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::set_olr_inst()
{
  mod_file->addStatement(new OlrInstStatement(*tmp_symbol_table));
  tmp_symbol_table->clear();
}

void
ParsingDriver::run_olr()
{
  mod_file->addStatement(new OlrStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

void
ParsingDriver::set_calib_var(string *name, string *weight, ExpObj *expression)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eEndogenous
      && mod_file->symbol_table.getType(*name) != eExogenous)
    error("calib_var: " + *name + " isn't an endogenous or exogenous variable");

  if (calib_var.find(*name) != calib_var.end())
    error("calib_var: " + *name + " declared twice");

  calib_var[*name] = pair<string, string>(*weight, get_expression(expression));

  delete name;
  delete weight;
  delete expression;
}

void
ParsingDriver::set_calib_covar(string *name1, string *name2,
                               string *weight, ExpObj *expression)
{
  check_symbol_existence(*name1);
  check_symbol_existence(*name2);
  if (mod_file->symbol_table.getType(*name1) != mod_file->symbol_table.getType(*name2))
    error("calib_var: " + *name1 + " and " + *name2 + "dont't have the same type");
  if (mod_file->symbol_table.getType(*name1) != eEndogenous
      && mod_file->symbol_table.getType(*name1) != eExogenous)
    error("calib_var: " + *name1 + " and " + *name2 + "aren't endogenous or exogenous variables");

  pair<string, string> covar_key(*name1, *name2);

  if (calib_covar.find(covar_key) != calib_covar.end())
    error("calib_var: pair of variables (" + *name1 + ", " + *name2
          + ") declared twice");

  calib_covar[covar_key] = pair<string, string>(*weight, get_expression(expression));

  delete name1;
  delete name2;
  delete weight;
  delete expression;
}

void
ParsingDriver::set_calib_ac(string *name, string *ar,
                            string *weight, ExpObj *expression)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eEndogenous)
    error("calib_var: " + *name + "isn't an endogenous variable");

  int iar = atoi(ar->c_str());
  pair<string, int> ac_key(*name, iar);

  if (calib_ac.find(ac_key) != calib_ac.end())
    error("calib_var: autocorr " + *name + "(" + *ar + ") declared twice");

  calib_ac[ac_key] = pair<string, string>(*weight, get_expression(expression));

  delete name;
  delete ar;
  delete weight;
  delete expression;
}

void
ParsingDriver::run_calib_var()
{
  mod_file->addStatement(new CalibVarStatement(calib_var, calib_covar, calib_ac,
                                               mod_file->symbol_table));
  calib_var.clear();
  calib_covar.clear();
  calib_ac.clear();
}

void
ParsingDriver::run_calib(int covar)
{
  mod_file->addStatement(new CalibStatement(covar));
}

void
ParsingDriver::run_dynatype(string *filename, string *ext)
{
  mod_file->addStatement(new DynaTypeStatement(*tmp_symbol_table, *filename, *ext));
  tmp_symbol_table->clear();
  delete filename;
  delete ext;
}

void
ParsingDriver::run_dynasave(string *filename, string *ext)
{
  mod_file->addStatement(new DynaSaveStatement(*tmp_symbol_table, *filename, *ext));
  tmp_symbol_table->clear();
  delete filename;
  delete ext;
}

void
ParsingDriver::add_mc_filename(string *filename, string *prior)
{
  if (filename_list.find(*filename) != filename_list.end())
    error("model_comparison: filename " + *filename + " declared twice");
  filename_list[*filename] = *prior;
  delete filename;
  delete prior;
}

void
ParsingDriver::run_model_comparison()
{
  mod_file->addStatement(new ModelComparisonStatement(filename_list, options_list));
  filename_list.clear();
  options_list.clear();
}

void
ParsingDriver::begin_planner_objective()
{
  model_tree = new ModelTree(mod_file->symbol_table, mod_file->num_constants);
}

void
ParsingDriver::end_planner_objective(NodeID expr)
{
  // Add equation corresponding to expression
  model_tree->AddEqual(expr, model_tree->Zero);
  model_tree->eq_nbr++;

  mod_file->addStatement(new PlannerObjectiveStatement(model_tree));
}

void
ParsingDriver::ramsey_policy()
{
  mod_file->addStatement(new RamseyPolicyStatement(*tmp_symbol_table, options_list));
  tmp_symbol_table->clear();
  options_list.clear();
}

NodeID
ParsingDriver::add_model_equal(NodeID arg1, NodeID arg2)
{
  NodeID id = model_tree->AddEqual(arg1, arg2);
  model_tree->eq_nbr++;
  return id;
}

NodeID
ParsingDriver::add_model_equal_with_zero_rhs(NodeID arg)
{
  return add_model_equal(arg, model_tree->Zero);
}

void
ParsingDriver::declare_and_init_local_parameter(string *name, NodeID rhs)
{
  mod_file->symbol_table.AddSymbolDeclar(*name, eLocalParameter, *name);
  NodeID id = model_tree->AddTerminal(*name);
  model_tree->AddAssign(id, rhs);
  delete name;
}

NodeID
ParsingDriver::add_model_plus(NodeID arg1, NodeID arg2)
{
  return model_tree->AddPlus(arg1, arg2);
}

NodeID
ParsingDriver::add_model_minus(NodeID arg1, NodeID arg2)
{
  return model_tree->AddMinus(arg1, arg2);
}

NodeID
ParsingDriver::add_model_uminus(NodeID arg1)
{
  return model_tree->AddUMinus(arg1);
}

NodeID
ParsingDriver::add_model_times(NodeID arg1, NodeID arg2)
{
  return model_tree->AddTimes(arg1, arg2);
}

NodeID
ParsingDriver::add_model_divide(NodeID arg1, NodeID arg2)
{
  return model_tree->AddDivide(arg1, arg2);
}

NodeID
ParsingDriver::add_model_power(NodeID arg1, NodeID arg2)
{
  return model_tree->AddPower(arg1, arg2);
}

NodeID
ParsingDriver::add_model_exp(NodeID arg1)
{
  return model_tree->AddExp(arg1);
}

NodeID
ParsingDriver::add_model_log(NodeID arg1)
{
  return model_tree->AddLog(arg1);
}

NodeID
ParsingDriver::add_model_log10(NodeID arg1)
{
  return model_tree->AddLog10(arg1);
}

NodeID
ParsingDriver::add_model_cos(NodeID arg1)
{
  return model_tree->AddCos(arg1);
}

NodeID
ParsingDriver::add_model_sin(NodeID arg1)
{
  return model_tree->AddSin(arg1);
}

NodeID
ParsingDriver::add_model_tan(NodeID arg1)
{
  return model_tree->AddTan(arg1);
}

NodeID
ParsingDriver::add_model_acos(NodeID arg1)
{
  return model_tree->AddACos(arg1);
}

NodeID
ParsingDriver::add_model_asin(NodeID arg1)
{
  return model_tree->AddASin(arg1);
}

NodeID
ParsingDriver::add_model_atan(NodeID arg1)
{
  return model_tree->AddATan(arg1);
}

NodeID
ParsingDriver::add_model_cosh(NodeID arg1)
{
  return model_tree->AddCosH(arg1);
}

NodeID
ParsingDriver::add_model_sinh(NodeID arg1)
{
  return model_tree->AddSinH(arg1);
}

NodeID
ParsingDriver::add_model_tanh(NodeID arg1)
{
  return model_tree->AddTanH(arg1);
}

NodeID
ParsingDriver::add_model_acosh(NodeID arg1)
{
  return model_tree->AddACosH(arg1);
}

NodeID
ParsingDriver::add_model_asinh(NodeID arg1)
{
  return model_tree->AddASinH(arg1);
}

NodeID
ParsingDriver::add_model_atanh(NodeID arg1)
{
  return model_tree->AddATanH(arg1);
}

NodeID
ParsingDriver::add_model_sqrt(NodeID arg1)
{
  return model_tree->AddSqRt(arg1);
}

void
ParsingDriver::add_native(const char *s)
{
  mod_file->addStatement(new NativeStatement(s));
}
