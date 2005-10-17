#ifndef DYNAREPARSER_H
#define DYNAREPARSER_H
//------------------------------------------------------------------------------
/*! \file 
 \version 1.0
 \date 04/27/2004
 \par This file defines the parser class.
*/
//------------------------------------------------------------------------------
#include <sstream>
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
#include "Objects.h"
#ifndef YLMM_basic_parser 
#include "ylmm/basic_parser.hh"
#endif
#ifndef YLMM_basic_scanner
#include "ylmm/basic_scanner.hh"
#endif
#ifndef __IOSTREAM__
#include <iostream>
#endif

//------------------------------------------------------------------------------
extern int yylineno;
/*! \namespace dynare
 */
namespace dynare 
{
  /*!
    \class  parser
    \brief  Member functions of this class are called from DyanreBison.y  
  */
  class parser : public ylmm::basic_parser<Objects*> 
  {
  private:
    ylmm::basic_scanner<Objects*>& _scanner;
    /*! Output file name */
	static string				file_name;
  public :
    /*! Refrence to output string */
    ostringstream				*output;
    /*! Stores model parameters */
    ModelParameters			model_parameters;
    /*! Stores symbol table */
  	SymbolTable				symbol_table;
    /*! Stores expressions */
  	Expression  			expression;
    /*! Stores numerical constants*/
  	NumericalConstants 		num_constants;
    /*! Handles numerical initalisations*/
  	NumericalInitialization numerical_initialization;
    /*! Handles shock command */
  	Shocks					shocks;
    /*! Handles sigma_e command */
  	SigmaeInitialization    sigmae;
    /*! Handles computing tasks commands*/
  	ComputingTasks			computing_tasks;
    /*! Stores temporary symbol table */
  	TmpSymbolTable			tmp_symbol_table;
    /*! Stores model tree */
  	ModelTree				model_tree;
    /*! Stores variable table */
  	VariableTable			variable_table;
    /*! Stores operator table */
  	OperatorTable			op_table;
	/*! Value of option order */
	int 					order;
	/*! Value of option linear */
	int 					linear;
	EstimationParams		estim_params;
    /*! Prints an arror to stdout */
	static void error(const char* m)
	{ 
		std::cout << file_name << " : Error in line " << yylineno << " : " << m << endl;
		exit(-1);
	}
    /*!
     Constuctor
	 \param s reference to scanner  
	 */
    parser(ylmm::basic_scanner<Objects*>& s) : _scanner(s) 
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
    /*! Destructor */
    virtual ~parser() {}                    
    /*!
     Scan input                          
	 \param arg Optional argument            
	 \return The next token ID 
	 */
    int scan(void* arg=0) {                 
      if (need_where())                     
		return _scanner.next(*_token,*_location); 
      return _scanner.next(*_token);
      
    }
    /*! Sets output file name */
    void set_file_name(string fname);
    /*! Sets reference to output string */
    void setoutput(ostringstream* ostr);
    /*! Remove unused symbol from symbol table */
    void check_model(void);
    /*! Executes final instructions */
    void finish(void);
    /*! Sets variable offset of ModelTree class to use C output */
    void use_dll(void);
    /*! Adds an endogenous variable to SymbolTable*/
    Objects* add_endogenous(Objects* name, Objects* tex_name = new Objects("",NULL, eUNDEF));
    /*! Adds an exogenous variable to SymbolTable*/
    Objects* add_exogenous(Objects* name, Objects* tex_name = new Objects("",NULL, eUNDEF));
    /*! Adds an exogenous determinist variable to SymbolTable*/
    Objects* add_exogenous_det(Objects* name, Objects* tex_name = new Objects("",NULL, eUNDEF));
    /*! Adds a parameter to SymbolTable*/
    Objects* add_parameter(Objects* name, Objects* tex_name = new Objects("",NULL, eUNDEF));
    /*! Adds a constant to NumericalConstants */
    Objects* add_constant(Objects* obj);
    /*! Adds a constant to ModelTree */
    Objects* add_model_constant(Objects* constant);
    /*! Adds a variable to VariableTable */
    Objects* add_variable(Objects* var);
    /*! Adds a lag variable to VariableTable */
    Objects* add_variable(Objects* var,Objects* olag);
    /*! Sets type and ID flag of an object */
    Objects* get_symbol(Objects* obj);
    /*! translates symbols for expressions */
    Objects* translate_symbol(Objects* obj);
    /*! Adds a binary token to an expression */
    Objects* add_expression_token( Objects* arg1,  Objects* arg2,  Objects* op);
    /*! Adds an unary token to an expression */
    Objects* add_expression_token( Objects* arg1, Objects* op);
    /*! Gets literal expression string */
    Objects* get_expression(Objects* exp);
    /* Concatenates two string objects */
    Objects* cat(Objects* string1, Objects* string2);
    /*! Writes parameter intitialisation expression */
	void init_param(Objects* lhs,  Objects* rhs);
    /*! Writes parameter intitialisation expression */
	void init_param(Objects* lhs);
    /*! Writes an initval block */
	void init_val(Objects* lhs,  Objects* rhs);
    /*! Writes an histval block */
	void hist_val(Objects* lhs, Objects* lag, Objects* rhs);
    /*! Writes an histval block */
	void hist_val(Objects* lhs, Objects* lag);
    /*! Writes begining of an initval block */
	void begin_initval(void);
    /*! Writes end of an initval block */
	void end_initval(void);
    /*! Writes begining of an endval block */
	void begin_endval(void);
    /*! Writes end of an endval block */
	void end_endval(void);
    /*! Writes begining of an histval block */
	void begin_histval(void);
	/*! Write begining of a shock block */
	void begin_shocks(void);
    /*! Adds a deterministic chock */
	void add_det_shock(Objects* var);
    /*! Adds a std error chock */
	void add_stderr_shock(Objects* var, Objects* value);
    /*! Adds a varriance chock */
	void add_var_shock(Objects* var, Objects* value);
    /*! Adds a covariance chock */
	void add_covar_shock(Objects* var1, Objects* var2, Objects* value);
    /*! Adds a correlated chock */
	void add_correl_shock(Objects* var1, Objects* var2, Objects* value);
    /*! Adds a shock period range  */
	void add_period(Objects* p1, Objects* p2);
    /*! Adds a shock period  */
	void add_period(Objects* p1);
    /*! Adds a shock value */
	void add_value(Objects* value);
    /*! Writes a Sigma_e block */
	void do_sigma_e(void);
    /*! Ends row of Sigma_e block */
	void end_of_row(void);
    /*! Adds an element to current row of Sigma_e */
	void add_to_row(Objects* s);
    /*! Write a steady command */
	void steady(void);
    /*! Sets an option to a numerical value */
	void option_num(string name_option, Objects* opt);
	void option_num(string name_option, Objects* opt1, Objects* opt2);
    /*! Sets an option to a string value */
	void option_str(string name_option, Objects* opt);
    /*! Sets an option to a numerical value */
	void option_num(string name_option, string opt);
    /*! Sets an option (string value) */
	void option_str(string name_option, string opt);
    /*! Adds a variable to temp symbol table and sets its value */
	void add_tmp_var(Objects* tmp_var1, Objects* tmp_var2);
    /*! Adds a variable to temp symbol table */
	void add_tmp_var(Objects* tmp_var);
    /*! Gets temp symbol table output */
	Objects* get_tmp_var(void);
	/*! Writes a rplot() command */
	void rplot(void);
	/*! Writes a stock_simul command */
	void stoch_simul(void);
    /*! Writes a simul command */
	void simul(void);
    /*! Writes check command */
	void check(void);
	/*! Writes instructions for estimation initialization */
	void estimation_init(void);
	/*! Writes instructions for estimated elements */
	void set_estimated_elements(void);
	void set_estimated_init_elements(void);
	void set_estimated_bounds_elements(void);
	/*! Runs estimation process */
	void run_estimation(void);
	/*! Prints optimization options */
	void optim_options(Objects* str1, Objects* str2, int task);
	void optim_options(int task);
	/*! Prints varops instructions */
	void set_varobs(void);
	void set_trend_init(void);
	void set_trend_element(Objects*, Objects*);
	void set_unit_root_vars(void);
	void begin_optim_weights(void);
	void set_optim_weights(Objects*,Objects*);
	void set_optim_weights(Objects*,Objects*,Objects*);
	void set_osr_params(void);
	void run_osr(void);
	void set_olr_inst(void);
	void run_olr(void);
	void begin_calib_var(void);
	void set_calib_var(Objects*,Objects*,Objects*);
	void set_calib_var(Objects*,Objects*,Objects*,Objects*);
	void set_calib_ac(Objects*,Objects*,Objects*,Objects*);
	void run_calib(int);
	void run_dynasave(Objects* arg1,Objects* arg2 = new Objects(""));
	void run_dynatype(Objects* arg1,Objects* arg2 = new Objects(""));
	void begin_model_comparison(void);
	void add_mc_filename(Objects* filename, Objects* prior = new Objects("1"));
	void run_model_comparison();
    /*! Writes token "arg1=arg2" to model tree */
	Objects*	add_equal(Objects* arg1,  Objects* arg2 = new Objects("0.0",ModelTree::Zero, eTempResult));
    /*! Writes token "arg1+arg2" to model tree */
	Objects*	add_plus(Objects* arg1,  Objects* arg2);
    /*! Writes token "arg1-arg2" to model tree */
	Objects*	add_minus(Objects* arg1,  Objects* arg2);
    /*! Writes token "-arg12" to model tree */
	Objects*	add_uminus(Objects* arg1);
    /*! Writes token "arg1*arg2" to model tree */
	Objects*	add_times(Objects* arg1,  Objects* arg2);
    /*! Writes token "arg1/arg2" to model tree */
	Objects*	add_divide(Objects* arg1,  Objects* arg2);
    /*! Writes token "arg1^arg2" to model tree */
	Objects*	add_power(Objects* arg1,  Objects* arg2);
    /*! Writes token "exp(arg1)" to model tree */
	Objects*	add_exp(Objects* arg1);
    /*! Writes token "log(arg1)" to model tree */
	Objects*	add_log(Objects* arg1);
    /*! Writes token "log10(arg1)" to model tree */
	Objects*	add_log10(Objects* arg1);
    /*! Writes token "cos(arg1)" to model tree */
	Objects*	add_cos(Objects* arg1);
    /*! Writes token "sin(arg1)" to model tree */
	Objects*	add_sin(Objects* arg1);
    /*! Writes token "tan(arg1)" to model tree */
	Objects*	add_tan(Objects* arg1);
    /*! Writes token "acos(arg1)" to model tree */
	Objects*	add_acos(Objects* arg1);
    /*! Writes token "asin(arg1)" to model tree */
	Objects*	add_asin(Objects* arg1);
    /*! Writes token "atan(arg1)" to model tree */
	Objects*	add_atan(Objects* arg1);
    /*! Writes token "cosh(arg1)" to model tree */
	Objects*	add_cosh(Objects* arg1);
    /*! Writes token "sinh(arg1)" to model tree */
	Objects*	add_sinh(Objects* arg1);
    /*! Writes token "tanh(arg1)" to model tree */
	Objects*	add_tanh(Objects* arg1);
    /*! Writes token "acosh(arg1)" to model tree */
	Objects*	add_acosh(Objects* arg1);
    /*! Writes token "asin(arg1)" to model tree */
	Objects*	add_asinh(Objects* arg1);
    /*! Writes token "atanh(arg1)" to model tree */
	Objects*	add_atanh(Objects* arg1);
    /*! Writes token "sqrt(arg1)" to model tree */
	Objects*	add_sqrt(Objects* arg1);
  };
}
//------------------------------------------------------------------------------
#endif
