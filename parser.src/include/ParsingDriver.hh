#ifndef _PARSING_DRIVER_HH
#define _PARSING_DRIVER_HH

#include <iostream>

#include "ModFile.hh"
#include "Expression.hh"
#include "TmpSymbolTable.hh"
#include "DynareBison.hh"
#include "ComputingTasks.hh"
#include "Shocks.hh"
#include "SigmaeInitialization.hh"
#include "NumericalInitialization.hh"
#include "ModelTree.hh"

using namespace std;

// Announce to Flex the prototype we want for lexing function, ...
#define YY_DECL                                              \
  yy::parser::token_type                                     \
  yylex(yy::parser::semantic_type *yylval,                   \
        yy::parser::location_type *yylloc,                   \
        ParsingDriver &driver)
// ... and declare it for the parser's sake.
YY_DECL;

//! Drives the scanning and parsing of the .mod file, and constructs its abstract representation
/*! It is built along the guidelines given in Bison 2.3 manual. */
class ParsingDriver
{
private:
  //! Start scanning
  /*! Body defined at the end of DynareFlex.ll, for convenience reasons. */
  void scan_begin();

  //! Stop scanning
  /*! Body defined at the end of DynareFlex.ll, for convenience reasons. */
  void scan_end();

  //! Returns string output of last expression parsed
  string get_expression(ExpObj *exp);

  //! Checks that a given symbol exists, and stops with an error message if it doesn't
  void check_symbol_existence(const string &name);

  //! Creates option "optim_opt" in OptionsList if it doesn't exist, else add a comma, and adds the option name
  void optim_options_helper(const string &name);

  //! Stores expressions
  Expression expression;
  //! Stores temporary symbol table
  TmpSymbolTable *tmp_symbol_table;

  //! The model tree in which to add expressions currently parsed
  ModelTree *model_tree;

  //! Stores options lists
  OptionsList options_list;
  //! Temporary storage for trend elements
  ObservationTrendsStatement::trend_elements_type trend_elements;
  //! Temporary storage for filename list of ModelComparison
  ModelComparisonStatement::filename_list_type filename_list;
  //! Temporary storage for list of EstimationParams (from estimated_params* statements)
  vector<EstimationParams> estim_params_list;
  //! Temporary storage of variances from optim_weights
  OptimWeightsStatement::var_weights_type var_weights;
  //! Temporary storage of covariances from optim_weights
  OptimWeightsStatement::covar_weights_type covar_weights;
  //! Temporary storage of variances from calib_var
  CalibVarStatement::calib_var_type calib_var;
  //! Temporary storage of covariances from calib_var
  CalibVarStatement::calib_covar_type calib_covar;
  //! Temporary storage of autocorrelations from calib_var
  CalibVarStatement::calib_ac_type calib_ac;
  //! Temporary storage for deterministic shocks
  ShocksStatement::det_shocks_type det_shocks;
  //! Temporary storage for periods of deterministic shocks
  vector<pair<int, int> > det_shocks_periods;
  //! Temporary storage for values of deterministic shocks
  vector<string> det_shocks_values;
  //! Temporary storage for variances of shocks
  ShocksStatement::var_and_std_shocks_type var_shocks;
  //! Temporary storage for standard errors of shocks
  ShocksStatement::var_and_std_shocks_type std_shocks;
  //! Temporary storage for covariances of shocks
  ShocksStatement::covar_and_corr_shocks_type covar_shocks;
  //! Temporary storage for correlations of shocks
  ShocksStatement::covar_and_corr_shocks_type corr_shocks;
  //! Temporary storage for Sigma_e rows
  SigmaeStatement::row_type sigmae_row;
  //! Temporary storage for Sigma_e matrix
  SigmaeStatement::matrix_type sigmae_matrix;
  //! Temporary storage for initval/endval blocks
  InitOrEndValStatement::init_values_type init_values;
  //! Temporary storage for histval blocks
  HistValStatement::hist_values_type hist_values;

  //! The mod file representation constructed by this ParsingDriver
  ModFile *mod_file;

public:
  //! Constructor
  ParsingDriver();
  //! Destructor
  virtual ~ParsingDriver();

  //! Starts parsing, and constructs the MOD file representation
  /*! \param f Name of file to parse

      The returned pointer should be deleted after use.
   */
  ModFile *parse(const string &f);

  //! Name of file being parsed
  string file;

  //! Trace scanning ?
  /*! If set to true before calling parse(), the flex scanner will dump a lot of debugging information. Defaults to false.
  */
  bool trace_scanning;

  //! Trace parsing ?
  /*! If set to true before calling parse(), the bison parser will dump debugging information. Defaults to false. */
  bool trace_parsing;

  //! Estimation parameters
  EstimationParams estim_params;

  //! Error handler with location
  void error(const yy::parser::location_type &l, const string &m);
  //! Error handler without location
  void error(const string &m);

  //! Static error handler
  /*! To be removed in the future. */
  static void error(const char *m)
  {
    extern int yylineno;
    cerr << "Error at line " << yylineno << ": " << m << endl;
    exit(-1);
  }

  //! Check if a given symbol exists in the parsing context
  bool exists_symbol(const char *s);
  //! Sets variable offset of ModelTree class to use C output
  void use_dll();
  //! Declares an endogenous variable by adding it to SymbolTable
  void declare_endogenous(string *name, string *tex_name = new string);
  //! Declares an exogenous variable by adding it to SymbolTable
  void declare_exogenous(string *name, string *tex_name = new string);
  //! Declares an exogenous deterministic variable by adding it to SymbolTable
  void declare_exogenous_det(string *name, string *tex_name = new string);
  //! Declares a parameter by adding it to SymbolTable
  void declare_parameter(string *name, string *tex_name = new string);
  //! Declares and initializes a local parameter
  void declare_and_init_local_parameter(string *name, NodeID rhs);
  //! Adds an Expression's numerical constant
  ExpObj *add_expression_constant(string *constant);
  //! Adds a model constant to ModelTree
  NodeID add_model_constant(string *constant);
  //! Adds a model variable to ModelTree and VariableTable
  NodeID add_model_variable(string *name);
  //! Adds a model lagged variable to ModelTree and VariableTable
  NodeID add_model_variable(string *name, string *olag);
  //! Adds an Expression's variable
  ExpObj *add_expression_variable(string *name);
  //! Adds a binary token to an expression
  ExpObj *add_expression_token(ExpObj *arg1, ExpObj *arg2, int op);
  //! Adds an unary token to an expression
  ExpObj *add_expression_token(ExpObj *arg1, int op);
  //! Adds a unary token to an expression, with function name unknown
  ExpObj *add_expression_token(ExpObj *arg1, string *op_name);
  //! Adds a "periods" statement
  void periods(string *periods);
  //! Adds a "dsample" statement
  void dsample(string *arg1);
  //! Adds a "dsample" statement
  void dsample(string *arg1, string *arg2);
  //! Writes parameter intitialisation expression
  void init_param(string *name, ExpObj *rhs);
  //! Writes an initval block
  void init_val(string *name, ExpObj *rhs);
  //! Writes an histval block
  void hist_val(string *name, string *lag, ExpObj *rhs);
  //! Writes end of an initval block
  void end_initval();
  //! Writes end of an endval block
  void end_endval();
  //! Writes end of an histval block
  void end_histval();
  //! Begin a model block
  void begin_model();
  //! Writes a shocks statement
  void end_shocks();
  //! Writes a mshocks statement
  void end_mshocks();
  //! Adds a deterministic chock
  void add_det_shock(string *var);
  //! Adds a std error chock
  void add_stderr_shock(string *var, ExpObj *value);
  //! Adds a variance chock
  void add_var_shock(string *var, ExpObj *value);
  //! Adds a covariance chock
  void add_covar_shock(string *var1, string *var2, ExpObj *value);
  //! Adds a correlated chock
  void add_correl_shock(string *var1, string *var2, ExpObj *value);
  //! Adds a shock period range 
  void add_period(string *p1, string *p2);
  //! Adds a shock period 
  void add_period(string *p1);
  //! Adds a shock value
  void add_value(string *value);
  //! Adds a shock value
  void add_value(ExpObj *value);
  //! Writes a Sigma_e block
  void do_sigma_e();
  //! Ends row of Sigma_e block
  void end_of_row();
  //! Adds an element to current row of Sigma_e
  void add_to_row(string *s);
  //! Adds an element to current row of Sigma_e
  void add_to_row(ExpObj *v);
  //! Write a steady command
  void steady();
  //! Sets an option to a numerical value
  void option_num(const string &name_option, string *opt);
  //! Sets an option to a numerical value
  void option_num(const string &name_option, const string &opt);
  //! Sets an option to a numerical value
  void option_num(const string &name_option, string *opt1, string *opt2);
  //! Sets an option to a string value
  void option_str(const string &name_option, string *opt);
  //! Sets an option to a string value
  void option_str(const string &name_option, const string &opt);
  //! Indicates that the model is linear
  void linear();
  //! Adds a variable to temp symbol table and sets its value
  void add_tmp_var(string *tmp_var1, string *tmp_var2);
  //! Adds a variable to temp symbol table
  void add_tmp_var(string *tmp_var);
  //! Writes a rplot() command
  void rplot();
  //! Writes a stock_simul command
  void stoch_simul();
  //! Writes a simul command
  void simul();
  //! Writes check command
  void check();
  //! Writes estimated params command
  void estimated_params();
  //! Writes estimated params init command
  void estimated_params_init();
  //! Writes estimated params bound command
  void estimated_params_bounds();
  //! Add a line in an estimated params block
  void add_estimated_params_element();
  //! Runs estimation process
  void run_estimation();
  //! Runs prior_analysis();
  void run_prior_analysis();
  //! Runs posterior_analysis();
  void run_posterior_analysis();
  //! Adds an optimization option (string value)
  void optim_options_string(string *name, string *value);
  //! Adds an optimization option (numeric value)
  void optim_options_num(string *name, string *value);
  //! Prints varops instructions
  void set_varobs();
  void set_trends();
  void set_trend_element(string *arg1, ExpObj *arg2);
  void set_unit_root_vars();
  void optim_weights();
  void set_optim_weights(string *name, ExpObj *value);
  void set_optim_weights(string *name1, string *name2, ExpObj *value);
  void set_osr_params();
  void run_osr();
  void set_olr_inst();
  void run_olr();
  void run_calib_var();
  void set_calib_var(string *name, string *weight, ExpObj *expression);
  void set_calib_covar(string *name1, string *name2, string *weight, ExpObj *expression);
  void set_calib_ac(string *name, string *ar, string *weight, ExpObj *expression);
  void run_calib(int covar);
  void run_dynasave(string *arg1, string *arg2 = new string);
  void run_dynatype(string *arg1, string *arg2 = new string);
  void add_mc_filename(string *filename, string *prior = new string("1"));
  void run_model_comparison();
  //! Begin a planner_objective statement
  void begin_planner_objective();
  //! End a planner objective statement
  void end_planner_objective(NodeID expr);
  //! ramsey policy statement
  void ramsey_policy();
  //! Writes token "arg1=arg2" to model tree
  NodeID add_model_equal(NodeID arg1, NodeID arg2);
  //! Writes token "arg=0" to model tree
  NodeID add_model_equal_with_zero_rhs(NodeID arg);
  //! Writes token "arg1+arg2" to model tree
  NodeID add_model_plus(NodeID arg1, NodeID arg2);
  //! Writes token "arg1-arg2" to model tree
  NodeID add_model_minus(NodeID arg1,  NodeID arg2);
  //! Writes token "-arg1" to model tree
  NodeID add_model_uminus(NodeID arg1);
  //! Writes token "arg1*arg2" to model tree
  NodeID add_model_times(NodeID arg1,  NodeID arg2);
  //! Writes token "arg1/arg2" to model tree
  NodeID add_model_divide(NodeID arg1,  NodeID arg2);
  //! Writes token "arg1^arg2" to model tree
  NodeID add_model_power(NodeID arg1,  NodeID arg2);
  //! Writes token "exp(arg1)" to model tree
  NodeID add_model_exp(NodeID arg1);
  //! Writes token "log(arg1)" to model tree
  NodeID add_model_log(NodeID arg1);
  //! Writes token "log10(arg1)" to model tree
  NodeID add_model_log10(NodeID arg1);
  //! Writes token "cos(arg1)" to model tree
  NodeID add_model_cos(NodeID arg1);
  //! Writes token "sin(arg1)" to model tree
  NodeID add_model_sin(NodeID arg1);
  //! Writes token "tan(arg1)" to model tree
  NodeID add_model_tan(NodeID arg1);
  //! Writes token "acos(arg1)" to model tree
  NodeID add_model_acos(NodeID arg1);
  //! Writes token "asin(arg1)" to model tree
  NodeID add_model_asin(NodeID arg1);
  //! Writes token "atan(arg1)" to model tree
  NodeID add_model_atan(NodeID arg1);
  //! Writes token "cosh(arg1)" to model tree
  NodeID add_model_cosh(NodeID arg1);
  //! Writes token "sinh(arg1)" to model tree
  NodeID add_model_sinh(NodeID arg1);
  //! Writes token "tanh(arg1)" to model tree
  NodeID add_model_tanh(NodeID arg1);
  //! Writes token "acosh(arg1)" to model tree
  NodeID add_model_acosh(NodeID arg1);
  //! Writes token "asin(arg1)" to model tree
  NodeID add_model_asinh(NodeID arg1);
  //! Writes token "atanh(arg1)" to model tree
  NodeID add_model_atanh(NodeID arg1);
  //! Writes token "sqrt(arg1)" to model tree
  NodeID add_model_sqrt(NodeID arg1);
  //! Adds a native statement
  void add_native(const char *s);
};

#endif // ! PARSING_DRIVER_HH
