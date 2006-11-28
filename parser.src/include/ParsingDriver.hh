#ifndef _PARSING_DRIVER_HH
#define _PARSING_DRIVER_HH

#include <iostream>

#include "ModFile.hh"
#include "Expression.hh"
#include "TmpSymbolTable.hh"
#include "DynareBison.hh"
#include "ComputingTasks.hh"

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

  //! Stores expressions
  Expression expression;
  //! Stores temporary symbol table
  TmpSymbolTable tmp_symbol_table;
  //! Stores operator table
  OperatorTable op_table;

  //! The mod file representation constructed by this ParsingDriver
  ModFile *mod_file;

  //! Reference to output string
  ostringstream *output;

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

  //! Sets reference to output string
  void setoutput(ostringstream *ostr);
  //! Executes final instructions
  void finish();
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
  //! Writes parameter intitialisation expression
  void init_param(string *name, ExpObj *rhs);
  //! Writes an initval block
  void init_val(string *name, ExpObj *rhs);
  //! Writes an histval block
  void hist_val(string *name, string *lag, ExpObj *rhs);
  //! Writes begining of an initval block
  void begin_initval();
  //! Writes end of an initval block
  void end_initval();
  //! Writes begining of an endval block
  void begin_endval();
  //! Writes end of an endval block
  void end_endval();
  //! Writes begining of an histval block
  void begin_histval();
  //! Write begining of a shock block
  void begin_shocks();
  void begin_mshocks();
  void end_shocks();
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
  //! Writes instructions for estimation initialization
  void estimation_init();
  //! Writes instructions for estimated elements
  void set_estimated_elements();
  void set_estimated_init_elements();
  void set_estimated_bounds_elements();
  //! Runs estimation process
  void run_estimation();
  //! Runs prior_analysis();
  void run_prior_analysis();
  //! Runs posterior_analysis();
  void run_posterior_analysis();
  //! Prints optimization options
  void optim_options(string *str1, string *str2, int task);
  //! Prints varops instructions
  void set_varobs();
  void set_trend_init();
  void set_trend_element(string *arg1, ExpObj *arg2);
  void set_unit_root_vars();
  void begin_optim_weights();
  void set_optim_weights(string *arg1, ExpObj *arg2);
  void set_optim_weights(string *arg1, string *arg2, ExpObj *arg3);
  void set_osr_params();
  void run_osr();
  void set_olr_inst();
  void run_olr();
  void begin_calib_var();
  void set_calib_var(string *name, string *weight, ExpObj *expression);
  void set_calib_var(string *name1, string *name2, string *weight, ExpObj *expression);
  void set_calib_ac(string *name, string *ar, string *weight, ExpObj *expression);
  void run_calib(int);
  void run_dynasave(string *arg1, string *arg2 = new string);
  void run_dynatype(string *arg1, string *arg2 = new string);
  void begin_model_comparison();
  void add_mc_filename(string *filename, string *prior = new string("1"));
  void run_model_comparison();
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
