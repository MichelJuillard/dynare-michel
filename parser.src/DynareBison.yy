%skeleton "lalr1.cc"
%require "2.3"
%defines

%{
namespace dynare
{
  class Objects;
}

class ParsingDriver;
%}

%parse-param { ParsingDriver& driver }
%lex-param { ParsingDriver& driver }

%locations
%initial-action
{
  // Initialize the location filenames
  @$.begin.filename = @$.end.filename = &driver.file;
};

%debug
%error-verbose

%union
{
  dynare::Objects *obj;
};

%{
#include "ParsingDriver.hh"
%}

%token AR AUTOCORR
%token BAYESIAN_IRF BETA_PDF
%token CALIB CALIB_VAR CHECK CONF_SIG CORR COVAR
%token DATAFILE DIAGNOSTIC DIFFUSE_D DOLLAR DR_ALGO DROP DSAMPLE DYN2VEC DYNASAVE DYNATYPE 
%token END ENDVAL EQUAL ESTIMATION ESTIMATED_PARAMS ESTIMATED_PARAMS_BOUNDS ESTIMATED_PARAMS_INIT
%token FILTERED_VARS FIRST_OBS
%token <obj> FLOAT_NUMBER
%token FORECAST FUNCTIONS
%token GAMMA_PDF GRAPH
%token HISTVAL HP_FILTER HP_NGRID    
%token INITVAL INITVALF
%token <obj> INT_NUMBER
%token INV_GAMMA_PDF INV_GAMMA1_PDF INV_GAMMA2_PDF IRF
%token KALMAN_ALGO KALMAN_TOL CONSTANT NOCONSTANT
%token LAPLACE LIK_ALGO LIK_INIT LINEAR LOAD_MH_FILE LOGLINEAR
%token MH_DROP MH_INIT_SCALE MH_JSCALE MH_MODE MH_NBLOCKS MH_REPLIC MODE_CHECK MODE_COMPUTE MODE_FILE MODEL MODEL_COMPARISON MODEL_COMPARISON_APPROXIMATION MODIFIEDHARMONICMEAN MOMENTS MOMENTS_VARENDO MSHOCKS
%token <obj> NAME
%token NOBS NOCORR NODIAGNOSTIC NOFUNCTIONS NOGRAPH XLS_SHEET XLS_RANGE
%token NOMOMENTS NOPRINT NORMAL_PDF
%token OBSERVATION_TRENDS OLR OLR_INST OLR_BETA OPTIM OPTIM_WEIGHTS ORDER OSR OSR_PARAMS 
%token PARAMETERS PERIODS PREFILTER PRESAMPLE PRINT PRIOR_TRUNC FILTER_STEP_AHEAD
%token QZ_CRITERIUM
%token RELATIVE_IRF REPLIC RESOL RPLOT
%token SHOCKS SIGMA_E SIMUL SIMUL_ALGO SIMUL_SEED SMOOTHER SOLVE_ALGO STDERR STEADY STOCH_SIMUL  
%token TEX
%token <obj> TEX_NAME
%token UNIFORM_PDF UNIT_ROOT_VARS USE_DLL
%token VALUES VAR VAREXO VAREXO_DET VAROBS
%token XTICK XTICKLABEL
%left <obj> COMMA
%left <obj> PLUS MINUS
%left <obj> TIMES DIVIDE
%left UMINUS
%right <obj> POWER 
%token <obj> EXP LOG LOG10 SIN COS TAN ASIN ACOS ATAN SINH COSH TANH ASINH ACOSH ATANH SQRT
/* isn't parsed from the *.mod file, but used to distinguish EQUAL in equation and EQUAL in assignment in    operation codes
*/
%token ASSIGN

%type <obj> var_list varexo_list varexo_det_list parameter_list equality_expression
%type <obj> expression comma_expression initval_elem histval_elem equation hand_side
%type <obj> pound_expression model_var shock_elem value_list triangular_row signed_integer
%type <obj> signed_float prior value trend_element optim_weights_list calib_arg2 filename
%type <obj> filename_elem range vec_int_elem vec_int_1 vec_int

%%

%start statement_list;

 statement_list
 	: statement
        | statement_list statement
        ;

 statement
  	: declaration
 	| periods
 	| model
 	| initval
 	| endval
 	| histval
 	| equality_expression
 	| shocks
 	| mshocks
 	| sigma_e
 	| steady
 	| check
 	| simul
 	| stoch_simul
 	| estimation
	| estimated_params
	| estimated_params_bounds
	| estimated_params_init
	| varobs
	| observation_trends
	| unit_root_vars
	| dsample
	| rplot
	| optim_weights
	| osr_params
	| osr
	| calib_var
	| calib
	| dynatype
	| dynasave
	| olr
	| olr_inst
        | model_comparison
	;

    
 declaration
 	: parameters
 	| var
 	| varexo
 	| varexo_det
 	;

 	
 dsample : DSAMPLE INT_NUMBER ';' {driver.option_num("dsample", $2);}
         | DSAMPLE INT_NUMBER INT_NUMBER ';' {driver.option_num("dsample", $2, $3);}
         ; 

 rplot : RPLOT tmp_var_list ';' {driver.rplot();}
      ; 
 var 
 	: VAR var_list ';' 
 	;

 varexo 
 	: VAREXO varexo_list ';'
	;

 varexo_det
 	: VAREXO_DET varexo_det_list ';'
	;
      	
 parameters
 	: PARAMETERS parameter_list ';'
 	;
 
 var_list
 	: var_list NAME  
 		{$$  = driver.add_endogenous($2);}
 	| var_list COMMA NAME  
 		{$$  = driver.add_endogenous($3);}
 	| NAME
 		{$$  = driver.add_endogenous($1);}
 	| var_list NAME  TEX_NAME 
 		{$$  = driver.add_endogenous($2,$3);}
 	| var_list COMMA NAME  TEX_NAME
 		{$$  = driver.add_endogenous($3, $4);}
 	| NAME TEX_NAME
 		{$$  = driver.add_endogenous($1, $2);}
	; 
	
 varexo_list
 	: varexo_list NAME
 		{$$  = driver.add_exogenous($2);}              
 	| varexo_list COMMA NAME
 		{$$  = driver.add_exogenous($3);}              
 	| NAME
 		{$$  = driver.add_exogenous($1);}
 	| varexo_list NAME TEX_NAME
 		{$$  = driver.add_exogenous($2, $3);}              
 	| varexo_list COMMA NAME TEX_NAME
 		{$$  = driver.add_exogenous($3, $4);}              
 	| NAME TEX_NAME
 		{$$  = driver.add_exogenous($1, $2);}
	; 

 varexo_det_list
 	: varexo_det_list NAME
 		{$$  = driver.add_exogenous_det($2);}              
 	| varexo_det_list COMMA NAME
 		{$$  = driver.add_exogenous_det($3);}              
 	| NAME
 		{$$  = driver.add_exogenous_det($1);}
 	| varexo_det_list NAME TEX_NAME
 		{$$  = driver.add_exogenous_det($2, $3);}              
 	| varexo_det_list COMMA NAME TEX_NAME
 		{$$  = driver.add_exogenous_det($3, $4);}              
 	| NAME TEX_NAME
 		{$$  = driver.add_exogenous_det($1, $2);}
	; 

 parameter_list
 	: parameter_list NAME
 		{$$  = driver.add_parameter($2);}
 	| parameter_list COMMA NAME
 		{$$  = driver.add_parameter($3);}
 	| NAME
 		{$$  = driver.add_parameter($1);}
 	| parameter_list NAME TEX_NAME
 		{$$  = driver.add_parameter($2, $3);}
 	| parameter_list COMMA NAME TEX_NAME
 		{$$  = driver.add_parameter($3, $4);}
 	| NAME TEX_NAME
 		{$$  = driver.add_parameter($1, $2);}
	; 	

 periods 
 	: PERIODS INT_NUMBER ';'
 		{
		 driver.option_num("periods", $2);
		 driver.option_num("simul", "1");
		}
    | PERIODS EQUAL INT_NUMBER ';'
 		{
		 driver.option_num("periods", $3);
		 driver.option_num("simul", "1");
		}
    ;

 		   
 equality_expression
 	: NAME EQUAL expression {$$ = driver.get_expression($3);} ';' 
    {driver.init_param($1, $<obj>4);} 
	;
 	
 expression
	: '(' expression ')'
		{ $$ = $2;}
	| NAME
		{$$ = driver.translate_symbol($1);}
	| FLOAT_NUMBER
		{$$ = driver.add_constant($1);}
	| INT_NUMBER
		{$$ = driver.add_constant($1);}
	| expression PLUS expression 
    	{$$ = driver.add_expression_token($1, $3, $2);}
	| expression MINUS expression
    	{$$ = driver.add_expression_token($1, $3, $2);}
	| expression DIVIDE expression	 
    	{$$ = driver.add_expression_token($1, $3, $2);}
	| expression TIMES expression 
    	{$$ = driver.add_expression_token($1, $3, $2);}
	| expression POWER expression 
    	{$$ = driver.add_expression_token($1, $3, $2);}	
	| MINUS expression %prec UMINUS
      {$1->opcode = token::UMINUS; $$ = driver.add_expression_token($2, $1);}
	| PLUS expression
	{$$ = $2;}
	| EXP '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| LOG '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| LOG10 '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| SIN '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| COS '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| TAN '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| ASIN '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| ACOS '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| ATAN '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
	| SQRT '(' expression ')'
    	{$$ = driver.add_expression_token($3, $1);}
        | NAME '(' expression ')' 
    	{$$ = driver.add_expression_token($3, $1);}
        | NAME '(' comma_expression ')' 
    	{$$ = driver.add_expression_token($3, $1);}
	; 

 comma_expression : expression COMMA expression {$$ = driver.add_expression_token($1, $3, $2);}
        | comma_expression COMMA expression {$$ = driver.add_expression_token($1, $3, $2);}

 initval
 	: INITVAL ';' {driver.begin_initval();} initval_list END
 		{driver.end_initval();}
 	;
 	
 endval
 	: ENDVAL ';' {driver.begin_endval();} initval_list END
 		{driver.end_endval();}
 	; 	

 initval_list
 	: initval_list initval_elem
	| initval_elem
	;

 initval_elem 
 	: NAME EQUAL expression {$$ = driver.get_expression($3);} ';'
    {driver.init_val($1, $<obj>4);}
 	; 	
 	
 histval
 	: HISTVAL ';' {driver.begin_histval();} histval_list END 

 histval_list
 	: histval_list histval_elem
	| histval_elem
	;
	
 histval_elem
 	: NAME '(' signed_integer ')' EQUAL expression 
                {$$ = driver.get_expression($6);} ';'
    {driver.hist_val($1, $3, $<obj>7);}
 	;
	
 model
 	: MODEL ';' equation_list END 
 		{driver.check_model();}
 	| MODEL '(' LINEAR ')' ';' {driver.option_num("linear","1");} 
		equation_list END {driver.check_model();}
 	| MODEL '(' USE_DLL ')' ';' {driver.use_dll();} 
		equation_list END {driver.check_model();}
 	;

 equation_list
    : equation_list equation     	
    | equation_list pound_expression
    | equation
    | pound_expression
    ;
    
 equation
 	: hand_side EQUAL hand_side ';'
 		{$$ = driver.add_equal($1, $3);}
 	| hand_side ';'
 		{$$ = driver.add_equal($1);}
	;
 
 hand_side
	: '(' hand_side ')' {$$ = $2;}
	| model_var
	| FLOAT_NUMBER
		{$$ = driver.add_model_constant($1);}
	| INT_NUMBER
		{$1->symbol += ".0"; $$ = driver.add_model_constant($1);}
	| hand_side PLUS hand_side 
    	{$$ = driver.add_plus($1, $3);}
	| hand_side MINUS hand_side
    	{$$ = driver.add_minus($1, $3);}
	| hand_side DIVIDE hand_side	 
    	{$$ = driver.add_divide($1, $3);}
	| hand_side TIMES hand_side 
    	{$$ = driver.add_times($1, $3);}
	| hand_side POWER hand_side 
    	{$$ = driver.add_power($1, $3);}	
        | MINUS hand_side %prec UMINUS
      {$1->opcode = token::UMINUS; $$ = driver.add_uminus($2);}
	| PLUS hand_side
	{$$ = $2;}
	| EXP '(' hand_side ')'
    	{$$ = driver.add_exp($3);}
	| LOG '(' hand_side ')'
    	{$$ = driver.add_log($3);}
	| LOG10 '(' hand_side ')'
    	{$$ = driver.add_log10($3);}
	| SIN '(' hand_side ')'
    	{$$ = driver.add_sin($3);}
	| COS '(' hand_side ')'
    	{$$ = driver.add_cos($3);}
	| TAN '(' hand_side ')'
    	{$$ = driver.add_tan($3);}
	| ASIN '(' hand_side ')'
    	{$$ = driver.add_asin($3);}
	| ACOS '(' hand_side ')'
    	{$$ = driver.add_acos($3);}
	| ATAN '(' hand_side ')'
    	{$$ = driver.add_atan($3);}
	| SQRT '(' hand_side ')'
    	{$$ = driver.add_sqrt($3);}
	;
	
 pound_expression: '#' NAME 
                   {$$ = driver.add_local_parameter($2);}
		   EQUAL hand_side ';'
                   {$$ = driver.init_local_parameter($<obj>3, $5);}

 model_var
 	: NAME 
 		{$$ = driver.add_variable($1);}
	| NAME '(' signed_integer ')'
		{$$ = driver.add_variable($1, $3);}
	;
	
 shocks
 	: SHOCKS  ';' {driver.begin_shocks();} shock_list END {driver.end_shocks();} 
 	;

 mshocks
 	: MSHOCKS  ';' {driver.begin_mshocks();} shock_list END {driver.end_shocks();} 
 	;

 shock_list 
 	: shock_list shock_elem
	| shock_elem
    ;

 shock_elem 
	: VAR NAME ';' PERIODS period_list ';' VALUES value_list ';'
		{driver.add_det_shock($2);}
	| VAR NAME ';' STDERR expression {$$ = driver.get_expression($5);} ';'
    {driver.add_stderr_shock($2, $<obj>6);}
	| VAR NAME EQUAL expression {$$ = driver.get_expression($4);} ';'
    {driver.add_var_shock($2, $<obj>5);}	
	| VAR NAME COMMA NAME EQUAL expression {$$ = driver.get_expression($6);} ';'
    {driver.add_covar_shock($2, $4, $<obj>7);}
	| CORR NAME COMMA NAME EQUAL expression {$$ = driver.get_expression($6);} ';'
    {driver.add_correl_shock($2, $4, $<obj>7);}
	;

 period_list
 	: period_list INT_NUMBER
 		{driver.add_period($2);}
	| period_list INT_NUMBER ':' INT_NUMBER 
 		{driver.add_period($2,$4);}
        | period_list COMMA INT_NUMBER
 		{driver.add_period($3);}
	| period_list COMMA INT_NUMBER ':' INT_NUMBER 
 		{driver.add_period($3, $5);}
	| INT_NUMBER ':' INT_NUMBER
 		{driver.add_period($1, $3);}
	| INT_NUMBER
 		{driver.add_period($1);}
	;


 value_list
 	: value_list signed_float  
		{driver.add_value($2);}
	| value_list signed_integer
		{driver.add_value($2);}
	| value_list NAME
		{driver.add_value($2);}
	| signed_float
		{driver.add_value($1);}
	| signed_integer
		{driver.add_value($1);}
	| NAME
		{driver.add_value($1);}
	| value_list '(' expression {$$ = driver.get_expression($3);} ')'
    {driver.add_value($<obj>4);}
	| '(' expression {$$ = driver.get_expression($2);} ')'
    {driver.add_value($<obj>3);}
	;
	
 sigma_e 
 	: SIGMA_E EQUAL '[' triangular_matrix ']' ';' 
 		{driver.do_sigma_e();}
    ;

 triangular_matrix 
 	: triangular_matrix ';' triangular_row 
 		{driver.end_of_row();}
    | triangular_row 
    	{driver.end_of_row();}
	;
	
 triangular_row 
 	: triangular_row COMMA '(' expression {$$ = driver.get_expression($4);} ')' 
    {driver.add_to_row($<obj>5);}
	| triangular_row COMMA FLOAT_NUMBER 
		{driver.add_to_row($3);}
	| triangular_row COMMA INT_NUMBER 
		{driver.add_to_row($3);}
	| triangular_row '(' expression {$$ = driver.get_expression($3);} ')' 
    {driver.add_to_row($<obj>4);}
	| triangular_row FLOAT_NUMBER 
		{driver.add_to_row($2);}
	| triangular_row INT_NUMBER 
		{driver.add_to_row($2);}
	| '(' expression {$$ = driver.get_expression($2);} ')' 
    {driver.add_to_row($<obj>3);}
	| FLOAT_NUMBER 
		{driver.add_to_row($1);}
	| INT_NUMBER 
		{driver.add_to_row($1);}
	;
	
 steady 
 	: STEADY ';' 
 		{
 			driver.steady();
 		}
    | STEADY '(' steady_options_list ')' ';'
    	 {driver.steady();}
    ;
    
 steady_options_list : steady_options_list COMMA steady_options
                     | steady_options
                     ;

 steady_options: o_solve_algo
               ;

 check 
 	: CHECK ';' 
 		{driver.check();}
	| CHECK '(' check_options_list ')' ';' 
		{driver.check();} 
	;

 check_options_list : check_options_list COMMA check_options
        | check_options
        ;

 check_options : o_solve_algo
               ;

 simul 
	: SIMUL ';' 
		{driver.simul();}
	| SIMUL '(' simul_options_list ')' ';'
		{driver.simul();}
        ;

 simul_options_list: simul_options_list COMMA simul_options
                   | simul_options
                   ;

 simul_options: o_periods
              ;

 stoch_simul
 	: STOCH_SIMUL ';'
 		{driver.stoch_simul();}
	| STOCH_SIMUL '(' stoch_simul_options_list ')' ';' 
		{driver.stoch_simul();}
	| STOCH_SIMUL tmp_var_list ';'
		{driver.stoch_simul();}
	| STOCH_SIMUL '(' stoch_simul_options_list ')' tmp_var_list ';' 
		{driver.stoch_simul();}
	;
	
 stoch_simul_options_list: stoch_simul_options_list COMMA stoch_simul_options
                         | stoch_simul_options
                         ;

 stoch_simul_options : o_dr_algo
                     | o_solve_algo
         	     | o_simul_algo
		     | o_linear
		     | o_order			   
                     | o_replic
                     | o_drop
                     | o_ar
                     | o_nocorr
                     | o_nofunction
                     | o_nomoments
                     | o_irf
                     | o_relative_irf
                     | o_hp_filter
                     | o_hp_ngrid
                     | o_periods
                     | o_simul
                     | o_simul_seed
                     | o_qz_criterium
                     ;

 tmp_var_list
	: tmp_var_list NAME
		{driver.add_tmp_var($2);}
	| tmp_var_list NAME EQUAL NAME
		{driver.add_tmp_var($2, $4);}
	| tmp_var_list COMMA NAME
		{driver.add_tmp_var($3);}
	| tmp_var_list COMMA NAME EQUAL NAME
		{driver.add_tmp_var($3, $5);}
	| NAME
		{driver.add_tmp_var($1);}
	| NAME EQUAL NAME
		{driver.add_tmp_var($1, $3);}
	;

 signed_integer
 	: PLUS INT_NUMBER
 		{$$ = $2;}
 	| MINUS INT_NUMBER
 		{$2->symbol.insert(0, "-"); $$ = $2;}
 	| INT_NUMBER
 		{$$ = $1;}
 	;
 	
 signed_float
 	: PLUS FLOAT_NUMBER
 		{$$ = $2;}
 	| MINUS FLOAT_NUMBER
 		{$2->symbol.insert(0, "-"); $$ = $2;}
 	| FLOAT_NUMBER
 		{$$ = $1;}
 	;

 estimated_params 
	: ESTIMATED_PARAMS ';' {driver.estimation_init();} estimated_list END
	;
	
 estimated_list 
	: estimated_list estimated_elem 
		{driver.set_estimated_elements();}
	| estimated_elem 
		{driver.set_estimated_elements();}
	;

 estimated_elem 
	: estimated_elem1 COMMA estimated_elem2 ';'
	;

 estimated_elem1 
	: STDERR NAME 
		{driver.estim_params.type = 1;
		 driver.estim_params.name = *$2;
		}
	| NAME
		{driver.estim_params.type = 2;
		 driver.estim_params.name = *$1;
		}
	| CORR NAME COMMA NAME
		{driver.estim_params.type = 3;
		 driver.estim_params.name = *$2;
		 driver.estim_params.name2 = *$4;
		}
	;

 estimated_elem2 
	: prior COMMA estimated_elem3 
		{driver.estim_params.prior=*$1;}
	| value COMMA prior COMMA estimated_elem3 
		{driver.estim_params.init_val=*$1;
		 driver.estim_params.prior=*$3;
		}
	| value COMMA value COMMA value COMMA prior COMMA estimated_elem3 
		{driver.estim_params.init_val=*$1;
		 driver.estim_params.low_bound=*$3;
		 driver.estim_params.up_bound=*$5;
		 driver.estim_params.prior=*$7;
		}
	| value 
		{driver.estim_params.init_val=*$1;}
	| value COMMA value COMMA value 
		{driver.estim_params.init_val=*$1;
		 driver.estim_params.low_bound=*$3;
		 driver.estim_params.up_bound=*$5;
		}
	;

 estimated_elem3 
 	: value COMMA value 
 		{driver.estim_params.mean=*$1;
 		 driver.estim_params.std=*$3;
 		}
	| value COMMA value COMMA value
		{driver.estim_params.mean=*$1;
		 driver.estim_params.std=*$3;
		 driver.estim_params.p3=*$5;
		}	
	| value COMMA value COMMA value COMMA value 
		{driver.estim_params.mean=*$1;
		 driver.estim_params.std=*$3;
		 driver.estim_params.p3=*$5;
		 driver.estim_params.p4=*$7;
		}
	| value COMMA value COMMA value COMMA value COMMA value 
		{driver.estim_params.mean=*$1;
		 driver.estim_params.std=*$3;
		 driver.estim_params.p3=*$5;
		 driver.estim_params.p4=*$7;
		 driver.estim_params.jscale=*$9;
		}
	;

 estimated_params_init: ESTIMATED_PARAMS_INIT ';' estimated_init_list END
                      ;

 estimated_init_list : estimated_init_list estimated_init_elem
                       {driver.set_estimated_init_elements();}
                     | estimated_init_elem
                       {driver.set_estimated_init_elements();}
                     ;

 estimated_init_elem : STDERR NAME COMMA value ';'
 		     		{driver.estim_params.type = 1;
				 driver.estim_params.name = *$2;
				 driver.estim_params.init_val=*$4;
				}
                     | CORR NAME COMMA NAME COMMA value ';'
 		     		{driver.estim_params.type = 3;
				 driver.estim_params.name = *$2;
				 driver.estim_params.name2 = *$4;
				 driver.estim_params.init_val=*$6;
				}
                     | NAME COMMA value ';'
 		     		{driver.estim_params.type = 2;
				 driver.estim_params.name = *$1;
				 driver.estim_params.init_val=*$3;
				}
                     ;

 estimated_params_bounds: ESTIMATED_PARAMS_BOUNDS ';' estimated_bounds_list END
                      ;

 estimated_bounds_list : estimated_bounds_list estimated_bounds_elem
                       {driver.set_estimated_bounds_elements();}
                     | estimated_bounds_elem
                       {driver.set_estimated_bounds_elements();}
                     ;

 estimated_bounds_elem : STDERR NAME COMMA value COMMA value ';'
 		     		{driver.estim_params.type = 1;
				 driver.estim_params.name = *$2;
				 driver.estim_params.low_bound=*$4;
				 driver.estim_params.up_bound=*$6;
				}
                     | CORR NAME COMMA NAME COMMA value COMMA value ';'
 		     		{driver.estim_params.type = 3;
				 driver.estim_params.name = *$2;
				 driver.estim_params.name2 = *$4;
				 driver.estim_params.low_bound=*$6;
				 driver.estim_params.up_bound=*$8;
				}
                     | NAME COMMA value COMMA value ';'
 		     		{driver.estim_params.type = 2;
				 driver.estim_params.name = *$1;
				 driver.estim_params.low_bound=*$3;
				 driver.estim_params.up_bound=*$5;
				}
                     ;

 prior
	: BETA_PDF 
		{$$ = new dynare::Objects("1");}
	| GAMMA_PDF 
		{$$ = new dynare::Objects("2");}
	| NORMAL_PDF 
		{$$ = new dynare::Objects("3");}
	| INV_GAMMA_PDF 
		{$$ = new dynare::Objects("4");}
	| UNIFORM_PDF 
		{$$ = new dynare::Objects("5");}
	;

 value 
	: {$$ = new dynare::Objects("NaN");} 
        | INT_NUMBER
	| FLOAT_NUMBER
	| NAME
        | MINUS INT_NUMBER
 		{$2->symbol.insert(0, "-"); $$ = $2;}
        | MINUS FLOAT_NUMBER
 		{$2->symbol.insert(0, "-"); $$ = $2;}
	;

 
       
 estimation 
	: ESTIMATION ';'
		{driver.run_estimation();}
	| ESTIMATION '(' estimation_options_list ')' ';' 
		{driver.run_estimation();}
	| ESTIMATION tmp_var_list ';'
		{driver.run_estimation();}
	| ESTIMATION '(' estimation_options_list ')' tmp_var_list ';' 
		{driver.run_estimation();}
	;

 estimation_options_list 
	: estimation_options_list COMMA estimation_options
	| estimation_options
	;

 estimation_options : o_datafile
                   | o_nobs
                   | o_first_obs
                   | o_prefilter
                   | o_presample
                   | o_lik_algo 
                   | o_lik_init 
                   | o_nograph
                   | o_conf_sig 
                   | o_mh_replic
                   | o_mh_drop
                   | o_mh_jscale
                   | o_optim
                   | o_mh_init_scale 
                   | o_mode_file 
                   | o_mode_compute 
                   | o_mode_check
                   | o_prior_trunc 
                   | o_mh_mode 
                   | o_mh_nblcks 
                   | o_load_mh_file 
                   | o_loglinear
                   | o_nodiagnostic
                   | o_bayesian_irf
                   | o_irf
                   | o_tex
                   | o_forecast
                   | o_smoother
                   | o_moments_varendo
                   | o_filtered_vars
                   | o_kalman_algo
                   | o_kalman_tol
                   | o_xls_sheet
                   | o_xls_range
                   | o_filter_step_ahead
                   | o_solve_algo
                   | o_constant
                   | o_noconstant
                   ;
	
 list_optim_option
 	: '\'' NAME '\'' COMMA '\'' NAME '\'' {driver.optim_options($2, $6, 2);}
	| '\'' NAME '\'' COMMA value {driver.optim_options($2, $5, 2);}
	;

 optim_options
	: list_optim_option
	| optim_options COMMA list_optim_option;
	;
	
 varobs 
 	: VAROBS tmp_var_list ';' 
 		{driver.set_varobs();}
	;

 observation_trends
        : OBSERVATION_TRENDS ';' {driver.set_trend_init();} trend_list END
	;

 trend_list 
        : trend_list trend_element
	| trend_element
	;

 trend_element :  NAME '(' expression {$$ = driver.get_expression($3);} ')' ';' {driver.set_trend_element($1, $<obj>4);}
               ;

 unit_root_vars : UNIT_ROOT_VARS tmp_var_list ';' {driver.set_unit_root_vars();}
                ;

 optim_weights : OPTIM_WEIGHTS ';' {driver.begin_optim_weights();} optim_weights_list END
               ;

 optim_weights_list : optim_weights_list NAME expression {$$ = driver.get_expression($3);} ';'  {driver.set_optim_weights($2, $<obj>4);}
                    | optim_weights_list NAME COMMA NAME expression {$$ = driver.get_expression($3);} ';' {driver.set_optim_weights($2, $4, $<obj>6);}
                    | NAME expression {$$ = driver.get_expression($2);} ';' {driver.set_optim_weights($1, $<obj>3);}
                    | NAME COMMA NAME expression {$$ = driver.get_expression($3);} ';' {driver.set_optim_weights($1, $3, $<obj>5);}
                    ;

 osr_params : OSR_PARAMS tmp_var_list ';' {driver.set_osr_params();}
            ;

 osr : OSR ';' {driver.run_osr();}
     | OSR '(' olr_options ')' ';' {driver.run_osr();}
     | OSR tmp_var_list ';' {driver.run_osr();}
     | OSR '(' olr_options ')' tmp_var_list ';' {driver.run_osr();}
     ;
 
 olr : OLR ';' {driver.run_olr();}
     | OLR '(' olr_options ')' ';' {driver.run_olr();}
     | OLR tmp_var_list ';' {driver.run_olr();}
     | OLR '(' olr_options ')' tmp_var_list ';' {driver.run_olr();}
     ;
 
 olr_option : o_olr_beta
     | stoch_simul_options
     ;
 
 olr_options : olr_option
     | olr_options COMMA olr_option
     ;

 olr_inst : OLR_INST tmp_var_list ';' {driver.set_olr_inst();}
          ;

 calib_var : CALIB_VAR ';' {driver.begin_calib_var();} calib_var_list END
           ;

 calib_var_list : calib_var_list calib_arg1
                | calib_arg1
                ;

 calib_arg1 : NAME calib_arg2 EQUAL expression ';' {driver.set_calib_var($1, $2, $4);}
            | NAME COMMA NAME calib_arg2 EQUAL expression ';' {driver.set_calib_var($1, $3, $4, $6);}
            | AUTOCORR NAME '(' INT_NUMBER ')' calib_arg2 EQUAL expression ';' {driver.set_calib_ac($2, $4, $6, $8);}
            ;

 calib_arg2 : {$$ = new dynare::Objects("1");}
            | '(' INT_NUMBER ')' {$$ = $2;}
            | '(' FLOAT_NUMBER ')' {$$ = $2;}
            ;

 calib : CALIB ';' {driver.run_calib(0);}
       | CALIB '(' COVAR ')' ';' {driver.run_calib(1);}
       ;

 dynatype : DYNATYPE '(' NAME ')'';' {driver.run_dynatype($3);}
          | DYNATYPE '(' NAME ')' tmp_var_list ';' {driver.run_dynatype($3);}
          | DYNATYPE NAME ';' {driver.run_dynatype($2);}
          | DYNATYPE '(' NAME '.' NAME ')'';' {driver.run_dynatype($3, $5);}
          | DYNATYPE '(' NAME '.' NAME ')' tmp_var_list ';' {driver.run_dynatype($3, $5);}
          | DYNATYPE NAME '.' NAME ';' {driver.run_dynatype($2,$4);};

 dynasave : DYNASAVE '(' NAME ')'';' {driver.run_dynasave($3);}
          | DYNASAVE '(' NAME ')' tmp_var_list ';' {driver.run_dynasave($3);}
          | DYNASAVE NAME ';' {driver.run_dynasave($2);}
          | DYNASAVE '(' NAME '.' NAME ')'';' {driver.run_dynasave($3, $5);}
          | DYNASAVE '(' NAME '.' NAME ')' tmp_var_list ';' {driver.run_dynasave($3, $5);}
          | DYNASAVE NAME '.' NAME ';' {driver.run_dynasave($2, $4);};

 model_comparison : MODEL_COMPARISON '(' model_comparison_options ')' {driver.begin_model_comparison();} 
                       filename_list ';' {driver.run_model_comparison();}
                  ;

 model_comparison_options: model_comparison_options ',' model_comparison_option
              | model_comparison_option
              ;

 model_comparison_option : o_model_comparison_approximation
              | o_print
              | o_noprint
              ;

 filename_list : filename {driver.add_mc_filename($1);}
        | filename_list ',' filename {driver.add_mc_filename($3);}
	| filename '(' value ')' {driver.add_mc_filename($1, $3);}
        | filename_list ',' filename '(' value ')' {driver.add_mc_filename($3, $5);}
        ;

 filename : filename_elem {$$ = $1;}
        | filename filename_elem {$$ = driver.cat($1, $2);}
        ;

 filename_elem : NAME
        | '\\' {$$ = new dynare::Objects("\\");}
        | '/' {$$ = new dynare::Objects("/");}
        | ':' {$$ = new dynare::Objects(":");}
        | '.' {$$ = new dynare::Objects(".");}
        ;

 o_dr_algo: DR_ALGO EQUAL INT_NUMBER {driver.option_num("dr_algo", $3);};
 o_solve_algo: SOLVE_ALGO EQUAL INT_NUMBER {driver.option_num("solve_algo", $3);};
 o_simul_algo: SIMUL_ALGO EQUAL INT_NUMBER {driver.option_num("simul_algo", $3);};
 o_linear: LINEAR {driver.option_num("linear", "1");};
 o_order: ORDER EQUAL INT_NUMBER {driver.option_num("order", $3);};
 o_replic: REPLIC EQUAL INT_NUMBER {driver.option_num("replic", $3);};
 o_drop: DROP EQUAL INT_NUMBER {driver.option_num("drop", $3);};
 o_ar: AR EQUAL INT_NUMBER {driver.option_num("ar", $3);};
 o_nocorr: NOCORR {driver.option_num("nocorr", "1");};
 o_nofunction: NOFUNCTIONS {driver.option_num("nofunctions", "1");};
 o_nomoments: NOMOMENTS {driver.option_num("nomoments", "1");};
 o_irf: IRF EQUAL INT_NUMBER {driver.option_num("irf", $3);};
 o_hp_filter: HP_FILTER EQUAL INT_NUMBER {driver.option_num("hp_filter", $3);};
 o_hp_ngrid: HP_NGRID EQUAL INT_NUMBER {driver.option_num("hp_ngrid", $3);};
 o_periods: PERIODS EQUAL INT_NUMBER {driver.option_num("periods", $3); driver.option_num("simul", "1");};
 o_simul: SIMUL {driver.option_num("simul", "1");};
 o_simul_seed: SIMUL_SEED EQUAL INT_NUMBER { driver.option_num("simul_seed", $3)};
 o_qz_criterium: QZ_CRITERIUM EQUAL INT_NUMBER { driver.option_num("qz_criterium", $3)}
               | QZ_CRITERIUM EQUAL FLOAT_NUMBER { driver.option_num("qz_criterium", $3)}
               ;
 o_datafile: DATAFILE EQUAL NAME {driver.option_str("datafile", $3);};
 o_nobs: NOBS EQUAL vec_int {driver.option_num("nobs", $3);}
       | NOBS EQUAL INT_NUMBER {driver.option_num("nobs", $3);}
       ;	      
 o_first_obs: FIRST_OBS EQUAL INT_NUMBER {driver.option_num("first_obs", $3);};
 o_prefilter: PREFILTER EQUAL INT_NUMBER {driver.option_num("prefilter", $3);};
 o_presample: PRESAMPLE EQUAL INT_NUMBER {driver.option_num("presample", $3);};
 o_lik_algo: LIK_ALGO EQUAL INT_NUMBER {driver.option_num("lik_algo", $3);}; 
 o_lik_init: LIK_INIT EQUAL INT_NUMBER {driver.option_num("lik_init", $3);}; 
 o_nograph: NOGRAPH {driver.option_num("nograph","1");}; 
 o_conf_sig: CONF_SIG EQUAL FLOAT_NUMBER {driver.option_num("conf_sig", $3);}; 
 o_mh_replic: MH_REPLIC EQUAL INT_NUMBER {driver.option_num("mh_replic", $3);}; 
 o_mh_drop: MH_DROP EQUAL FLOAT_NUMBER {driver.option_num("mh_drop", $3);}; 
 o_mh_jscale: MH_JSCALE EQUAL FLOAT_NUMBER {driver.option_num("mh_jscale", $3);}; 
 o_optim: OPTIM  EQUAL '(' optim_options ')';
 o_mh_init_scale :MH_INIT_SCALE EQUAL FLOAT_NUMBER {driver.option_num("mh_init_scale", $3);};
 o_mh_init_scale :MH_INIT_SCALE EQUAL INT_NUMBER {driver.option_num("mh_init_scale", $3);};
 o_mode_file : MODE_FILE EQUAL NAME {driver.option_str("mode_file", $3);};
 o_mode_compute : MODE_COMPUTE EQUAL INT_NUMBER {driver.option_num("mode_compute", $3);};
 o_mode_check : MODE_CHECK {driver.option_num("mode_check", "1");};
 o_prior_trunc : PRIOR_TRUNC EQUAL FLOAT_NUMBER {driver.option_num("prior_trunc", $3);};
 o_mh_mode : MH_MODE EQUAL INT_NUMBER {driver.option_num("mh_mode", $3);};
 o_mh_nblcks : MH_NBLOCKS EQUAL INT_NUMBER {driver.option_num("mh_nblck", $3);};
 o_load_mh_file : LOAD_MH_FILE {driver.option_num("load_mh_file", "1");};
 o_loglinear : LOGLINEAR {driver.option_num("loglinear", "1");};
 o_nodiagnostic : NODIAGNOSTIC {driver.option_num("nodiagnostic", "1");};
 o_bayesian_irf : BAYESIAN_IRF {driver.option_num("bayesian_irf", "1");};
 o_tex : TEX {driver.option_num("TeX", "1");};
 o_forecast : FORECAST EQUAL INT_NUMBER {driver.option_num("forecast", $3);};
 o_smoother : SMOOTHER {driver.option_num("smoother", "1");};
 o_moments_varendo : MOMENTS_VARENDO {driver.option_num("moments_varendo", "1");};
 o_filtered_vars : FILTERED_VARS {driver.option_num("filtered_vars", "1");};
 o_relative_irf : RELATIVE_IRF {driver.option_num("relative_irf", "1");};
 o_kalman_algo : KALMAN_ALGO EQUAL INT_NUMBER {driver.option_num("kalman_algo", $3);};
 o_kalman_tol : KALMAN_TOL EQUAL INT_NUMBER {driver.option_num("kalman_tol", $3);};
 o_olr_beta: OLR_BETA EQUAL value {driver.option_num("olr_beta", $3);};
 o_model_comparison_approximation: MODEL_COMPARISON_APPROXIMATION EQUAL LAPLACE {dynare::Objects* tmp = new dynare::Objects("Laplace"); driver.option_str("model_comparison_approximation", tmp);}
   | MODEL_COMPARISON_APPROXIMATION EQUAL MODIFIEDHARMONICMEAN {dynare::Objects* tmp = new dynare::Objects("MODIFIEDHARMONICMEAN"); driver.option_str("model_comparison_approximation", tmp);}
   ;
 o_print : PRINT {driver.option_num("noprint", "0");};
 o_noprint : NOPRINT {driver.option_num("noprint", "1");};
 o_nograph : GRAPH {driver.option_num("nograph", "0");};	
 o_xls_sheet : XLS_SHEET EQUAL NAME {driver.option_str("xls_sheet", $3);} 
 o_xls_range : XLS_RANGE EQUAL range {driver.option_str("xls_range", $3);} 
 o_filter_step_ahead : FILTER_STEP_AHEAD EQUAL vec_int {driver.option_num("filter_step_ahead", $3);}
 o_constant : CONSTANT {driver.option_num("noconstant", "0");}
 o_noconstant : NOCONSTANT {driver.option_num("noconstant", "1");}

 range : NAME ':' NAME {$$ = new dynare::Objects(":");$$ = driver.cat($1, $$);$$ = driver.cat($$, $3);}
 vec_int_elem : INT_NUMBER
 	      | INT_NUMBER ':' {$$ = new dynare::Objects(":"); $$ = driver.cat($1, $$); }
          INT_NUMBER {$$ = driver.cat($<obj>3,$4);}   		    	        
	      ;

 vec_int_1 : '[' vec_int_elem {$$ = new dynare::Objects("["); $$ = driver.cat($$, $2);}
           | vec_int_1 vec_int_elem {$$ = driver.cat_with_space($1, $2);}
           ;

 vec_int : vec_int_1 ']' {$$ = new dynare::Objects("]"); $$ = driver.cat($1, $$);};

%%

void
yy::parser::error(const yy::parser::location_type& l,
                  const std::string& m)
{
  driver.error(l, m);
}

/*
  Local variables:
  mode: C++
  End:
*/
