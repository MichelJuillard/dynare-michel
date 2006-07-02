%{
  /* Declarations */
#include "DynareParser.h"
#define yyparse tcparse 
#define yylex tclex
//#define yyerror tcerror
#define yylval tclval
#define yychar tcchar 
//#define yydebug 
#define YYDEBUG 1
#define yynerrs tcnerrs
//#define YYLSP_NEEDED 1
#define YYERROR_VERBOSE
#define YLMM_PARSER_CLASS dynare::parser 
#define YLMM_LEX_STATIC
#include "ylmm/yaccmm.hh"
#ifdef yyerror
#undef yyerror
#define yyerror _parser->error 
#endif
%}
/* %pure_parser */
%token AR AUTOCORR
%token BAYESIAN_IRF BETA_PDF
%token CALIB CALIB_VAR CHECK CONF_SIG CORR COVAR
%token DATAFILE DIAGNOSTIC DIFFUSE_D DOLLAR DR_ALGO DROP DSAMPLE DYN2VEC DYNASAVE DYNATYPE 
%token END ENDVAL EQUAL ESTIMATION ESTIMATED_PARAMS ESTIMATED_PARAMS_BOUNDS ESTIMATED_PARAMS_INIT
%token FILTERED_VARS FIRST_OBS FLOAT_NUMBER FORECAST FUNCTIONS
%token GAMMA_PDF GRAPH
%token HISTVAL HP_FILTER HP_NGRID    
%token INITVAL INITVALF INT_NUMBER INV_GAMMA_PDF INV_GAMMA1_PDF INV_GAMMA2_PDF IRF
%token KALMAN_ALGO KALMAN_TOL
%token LAPLACE LIK_ALGO LIK_INIT LINEAR LOAD_MH_FILE LOGLINEAR
%token MH_DROP MH_INIT_SCALE MH_JSCALE MH_MODE MH_NBLOCKS MH_REPLIC MODE_CHECK MODE_COMPUTE MODE_FILE MODEL MODEL_COMPARISON MODEL_COMPARISON_APPROXIMATION MODIFIEDHARMONICMEAN MOMENTS MOMENTS_VARENDO MSHOCKS
%token NAME NOBS NOCORR NODIAGNOSTIC NOFUNCTIONS NOGRAPH XLS_SHEET XLS_RANGE
%token NOMOMENTS NOPRINT NORMAL_PDF
%token OBSERVATION_TRENDS OLR OLR_INST OLR_BETA OPTIM OPTIM_WEIGHTS ORDER OSR OSR_PARAMS 
%token PARAMETERS PERIODS PREFILTER PRESAMPLE PRINT PRIOR_TRUNC FILTER_STEP_AHEAD
%token QZ_CRITERIUM
%token RELATIVE_IRF REPLIC RESOL RPLOT
%token SHOCKS SIGMA_E SIMUL SIMUL_ALGO SIMUL_SEED SMOOTHER SOLVE_ALGO STDERR STEADY STOCH_SIMUL  
%token TEX TEX_NAME
%token UNIFORM_PDF UNIT_ROOT_VARS USE_DLL
%token VALUES VAR VAREXO VAREXO_DET VAROBS
%token XTICK XTICKLABEL  
%left COMMA
%left PLUS MINUS 
%left TIMES DIVIDE
%left UMINUS
%right POWER 
%nonassoc FACTORIAL 
%token EXP LOG LOG10 LN SIN COS TAN ASIN ACOS ATAN SINH COSH TANH ASINH ACOSH ATANH SQRT
%%

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

 	
 dsample : DSAMPLE INT_NUMBER ';' {_parser->option_num("dsample",$2);}
         | DSAMPLE INT_NUMBER INT_NUMBER ';' {_parser->option_num("dsample",$2,$3);}
         ; 

 rplot : RPLOT tmp_var_list ';' {_parser->rplot();}
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
 		{$$  = _parser->add_endogenous($2);}
 	| var_list COMMA NAME  
 		{$$  = _parser->add_endogenous($3);}
 	| NAME
 		{$$  = _parser->add_endogenous($1);}
 	| var_list NAME  TEX_NAME 
 		{$$  = _parser->add_endogenous($2,$3);}
 	| var_list COMMA NAME  TEX_NAME
 		{$$  = _parser->add_endogenous($3,$4);}
 	| NAME TEX_NAME
 		{$$  = _parser->add_endogenous($1,$2);}
	; 
	
 varexo_list
 	: varexo_list NAME
 		{$$  = _parser->add_exogenous($2);}              
 	| varexo_list COMMA NAME
 		{$$  = _parser->add_exogenous($3);}              
 	| NAME
 		{$$  = _parser->add_exogenous($1);}
 	| varexo_list NAME TEX_NAME
 		{$$  = _parser->add_exogenous($2, $3);}              
 	| varexo_list COMMA NAME TEX_NAME
 		{$$  = _parser->add_exogenous($3, $4);}              
 	| NAME TEX_NAME
 		{$$  = _parser->add_exogenous($1, $2);}
	; 

 varexo_det_list
 	: varexo_det_list NAME
 		{$$  = _parser->add_exogenous_det($2);}              
 	| varexo_det_list COMMA NAME
 		{$$  = _parser->add_exogenous_det($3);}              
 	| NAME
 		{$$  = _parser->add_exogenous_det($1);}
 	| varexo_det_list NAME TEX_NAME
 		{$$  = _parser->add_exogenous_det($2, $3);}              
 	| varexo_det_list COMMA NAME TEX_NAME
 		{$$  = _parser->add_exogenous_det($3, $4);}              
 	| NAME TEX_NAME
 		{$$  = _parser->add_exogenous_det($1, $2);}
	; 

 parameter_list
 	: parameter_list NAME
 		{$$  = _parser->add_parameter($2);}
 	| parameter_list COMMA NAME
 		{$$  = _parser->add_parameter($3);}
 	| NAME
 		{$$  = _parser->add_parameter($1);}
 	| parameter_list NAME TEX_NAME
 		{$$  = _parser->add_parameter($2, $3);}
 	| parameter_list COMMA NAME TEX_NAME
 		{$$  = _parser->add_parameter($3,$4);}
 	| NAME TEX_NAME
 		{$$  = _parser->add_parameter($1, $2);}
	; 	

 periods 
 	: PERIODS INT_NUMBER ';'
 		{
		 _parser->option_num("periods",$2);
		 _parser->option_num("simul","1");
		}
    | PERIODS EQUAL INT_NUMBER ';'
 		{
		 _parser->option_num("periods",$3);
		 _parser->option_num("simul","1");
		}
    ;

 		   
 equality_expression
 	: NAME EQUAL expression {$$ = _parser->get_expression($3);} ';' 
    	{_parser->init_param($1, $4);} 
	;
 	
 expression
	: '(' expression ')'
		{ $$ = $2;}
	| NAME
		{$$ = _parser->translate_symbol($1);}
	| FLOAT_NUMBER
		{$$ = _parser->add_constant($1);}
	| INT_NUMBER
		{$$ = _parser->add_constant($1);}
	| expression PLUS expression 
    	{$$ = _parser->add_expression_token($1,$3, $2);}
	| expression MINUS expression
    	{$$ = _parser->add_expression_token($1,$3, $2);}
	| expression DIVIDE expression	 
    	{$$ = _parser->add_expression_token($1,$3, $2);}
	| expression TIMES expression 
    	{$$ = _parser->add_expression_token($1,$3, $2);}
	| expression POWER expression 
    	{$$ = _parser->add_expression_token($1,$3, $2);}	
	| MINUS expression %prec UMINUS
    	{$1->opcode = UMINUS;$$ = _parser->add_expression_token($2,$1);}
	| PLUS expression
	{$$ = $2;}
	| EXP '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| LOG '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| LOG10 '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| LN '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| SIN '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| COS '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| TAN '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| ASIN '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| ACOS '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| ATAN '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
	| SQRT '(' expression ')'
    	{$$ = _parser->add_expression_token($3, $1);}
        | NAME '(' expression ')' 
    	{$$ = _parser->add_expression_token($3, $1);}
        | NAME '(' comma_expression ')' 
    	{$$ = _parser->add_expression_token($3, $1);}
	; 

 comma_expression : expression COMMA expression {$$ = _parser->add_expression_token($1,$3,$2);}
        | comma_expression COMMA expression {$$ = _parser->add_expression_token($1,$3,$2);}

 initval
 	: INITVAL ';' {_parser->begin_initval();} initval_list END
 		{_parser->end_initval();}
 	;
 	
 endval
 	: ENDVAL ';' {_parser->begin_endval();} initval_list END
 		{_parser->end_endval();}
 	; 	

 initval_list
 	: initval_list initval_elem
	| initval_elem
	;

 initval_elem 
 	: NAME EQUAL expression {$$ = _parser->get_expression($3);} ';'
 		{_parser->init_val($1, $4);}
 	; 	
 	
 histval
 	: HISTVAL ';' {_parser->begin_histval();} histval_list END 

 histval_list
 	: histval_list histval_elem
	| histval_elem
	;
	
 histval_elem
 	: NAME '(' signed_integer ')' EQUAL expression ';'
 		{_parser->hist_val($1, $3, $6);}
 	;
	
 model
 	: MODEL {_parser->initialize_model();} ';' equation_list END 
 		{_parser->check_model();}
 	| MODEL '(' LINEAR ')' ';' {_parser->option_num("linear","1");
	  	    	       	   _parser->initialize_model();} 
		equation_list END {_parser->check_model();}
 	| MODEL '(' USE_DLL ')' ';' {_parser->use_dll();
	  	    	    	    _parser->initialize_model();} 
		equation_list END {_parser->check_model();}
 	;

 equation_list
    : equation_list equation     	
    | equation_list pound_expression
    | equation
    | pound_expression
    ;
    
 equation
 	: hand_side EQUAL hand_side ';'
 		{$$ = _parser->add_equal($1,$3);}
 	| hand_side ';'
 		{$$ = _parser->add_equal($1);}
	;
 
 hand_side
	: '(' hand_side ')' {$$ = $2;}
	| model_var
	| FLOAT_NUMBER
		{$$ = _parser->add_model_constant($1);}
	| INT_NUMBER
		{$1->symbol += ".0";$$ = _parser->add_model_constant($1);}
	| hand_side PLUS hand_side 
    	{$$ = _parser->add_plus($1,$3);}
	| hand_side MINUS hand_side
    	{$$ = _parser->add_minus($1,$3);}
	| hand_side DIVIDE hand_side	 
    	{$$ = _parser->add_divide($1,$3);}
	| hand_side TIMES hand_side 
    	{$$ = _parser->add_times($1,$3);}
	| hand_side POWER hand_side 
    	{$$ = _parser->add_power($1,$3);}	
        | MINUS hand_side %prec UMINUS
    	{$1->opcode = UMINUS;$$ = _parser->add_uminus($2);}
	| PLUS hand_side
	{$$ = $2;}
	| EXP '(' hand_side ')'
    	{$$ = _parser->add_exp($3);}
	| LOG '(' hand_side ')'
    	{$$ = _parser->add_log($3);}
	| LOG10 '(' hand_side ')'
    	{$$ = _parser->add_log10($3);}
	| LN '(' hand_side ')'
    	{$$ = _parser->add_log($3);}
	| SIN '(' hand_side ')'
    	{$$ = _parser->add_sin($3);}
	| COS '(' hand_side ')'
    	{$$ = _parser->add_cos($3);}
	| TAN '(' hand_side ')'
    	{$$ = _parser->add_tan($3);}
	| ASIN '(' hand_side ')'
    	{$$ = _parser->add_asin($3);}
	| ACOS '(' hand_side ')'
    	{$$ = _parser->add_acos($3);}
	| ATAN '(' hand_side ')'
    	{$$ = _parser->add_atan($3);}
	| SQRT '(' hand_side ')'
    	{$$ = _parser->add_sqrt($3);}
	;
	
 pound_expression: '#' NAME 
                   {$$ = _parser->add_local_parameter($2);}
		   EQUAL hand_side ';'
                   {$$ = _parser->init_local_parameter($3,$5);}

 model_var
 	: NAME 
 		{$$ = _parser->add_variable($1);}
	| NAME '(' signed_integer ')'
		{$$ = _parser->add_variable($1,$3);}
	;
	
 shocks
 	: SHOCKS  ';' {_parser->begin_shocks();} shock_list END {_parser->end_shocks();} 
 	;

 mshocks
 	: MSHOCKS  ';' {_parser->begin_mshocks();} shock_list END {_parser->end_shocks();} 
 	;

 shock_list 
 	: shock_list shock_elem
	| shock_elem
    ;

 shock_elem 
	: VAR NAME ';' PERIODS period_list ';' VALUES value_list ';'
		{_parser->add_det_shock($2);}
	| VAR NAME ';' STDERR expression {$$ = _parser->get_expression($5);} ';'
		{_parser->add_stderr_shock($2, $6);}
	| VAR NAME EQUAL expression {$$ = _parser->get_expression($4);} ';'
		{_parser->add_var_shock($2, $5);}	
	| VAR NAME COMMA NAME EQUAL expression {$$ = _parser->get_expression($6);} ';'
		{_parser->add_covar_shock($2, $4, $7);}
	| CORR NAME COMMA NAME EQUAL expression {$$ = _parser->get_expression($6);} ';'
		{_parser->add_correl_shock($2, $4, $7);}
	;

 period_list
 	: period_list INT_NUMBER
 		{_parser->add_period($2);}
	| period_list INT_NUMBER ':' INT_NUMBER 
 		{_parser->add_period($2,$4);}
        | period_list COMMA INT_NUMBER
 		{_parser->add_period($3);}
	| period_list COMMA INT_NUMBER ':' INT_NUMBER 
 		{_parser->add_period($3,$5);}
	| INT_NUMBER ':' INT_NUMBER
 		{_parser->add_period($1,$3);}
	| INT_NUMBER
 		{_parser->add_period($1);}
	;


 value_list
 	: value_list signed_float  
		{_parser->add_value($2);}
	| value_list signed_integer
		{_parser->add_value($2);}
	| value_list NAME
		{_parser->add_value($2);}
	| signed_float
		{_parser->add_value($1);}
	| signed_integer
		{_parser->add_value($1);}
	| NAME
		{_parser->add_value($1);}
	| value_list '(' expression {$$ = _parser->get_expression($3);} ')'
		{_parser->add_value($4);}
	| '(' expression {$$ = _parser->get_expression($2);} ')'
		{_parser->add_value($3);}
	;
	
 sigma_e 
 	: SIGMA_E EQUAL '[' triangular_matrix ']' ';' 
 		{_parser->do_sigma_e();}
    ;

 triangular_matrix 
 	: triangular_matrix ';' triangular_row 
 		{_parser->end_of_row();}
    | triangular_row 
    	{_parser->end_of_row();}
	;
	
 triangular_row 
 	: triangular_row COMMA '(' expression {$$ = _parser->get_expression($4);} ')' 
 		{_parser->add_to_row($5);}
	| triangular_row COMMA FLOAT_NUMBER 
		{_parser->add_to_row($3);}
	| triangular_row COMMA INT_NUMBER 
		{_parser->add_to_row($3);}
	| triangular_row '(' expression {$$ = _parser->get_expression($3);} ')' 
		{_parser->add_to_row($4);}
	| triangular_row FLOAT_NUMBER 
		{_parser->add_to_row($2);}
	| triangular_row INT_NUMBER 
		{_parser->add_to_row($2);}
	| '(' expression {$$ = _parser->get_expression($2);} ')' 
		{_parser->add_to_row($3);}
	| FLOAT_NUMBER 
		{_parser->add_to_row($1);}
	| INT_NUMBER 
		{_parser->add_to_row($1);}
	;
	
 steady 
 	: STEADY ';' 
 		{
 			_parser->steady();
 		}
    | STEADY '(' steady_options_list ')' ';'
    	 {_parser->steady();}
    ;
    
 steady_options_list : steady_options_list COMMA steady_options
                     | steady_options
                     ;

 steady_options: o_solve_algo
               ;

 check 
 	: CHECK ';' 
 		{_parser->check();}
	| CHECK '(' check_options_list ')' ';' 
		{_parser->check();} 
	;

 check_options_list : check_options_list COMMA check_options
        | check_options
        ;

 check_options : o_solve_algo
               ;

 simul 
	: SIMUL ';' 
		{_parser->simul();}
	| SIMUL '(' simul_options_list ')' ';'
		{_parser->simul();}
        ;

 simul_options_list: simul_options_list COMMA simul_options
                   | simul_options
                   ;

 simul_options: o_periods
              ;

 stoch_simul
 	: STOCH_SIMUL ';'
 		{_parser->stoch_simul();}
	| STOCH_SIMUL '(' stoch_simul_options_list ')' ';' 
		{_parser->stoch_simul();}
	| STOCH_SIMUL tmp_var_list ';'
		{_parser->stoch_simul();}
	| STOCH_SIMUL '(' stoch_simul_options_list ')' tmp_var_list ';' 
		{_parser->stoch_simul();}
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
		{_parser->add_tmp_var($2);}
	| tmp_var_list NAME EQUAL NAME
		{_parser->add_tmp_var($2,$4);}
	| tmp_var_list COMMA NAME
		{_parser->add_tmp_var($3);}
	| tmp_var_list COMMA NAME EQUAL NAME
		{_parser->add_tmp_var($3,$5);}
	| NAME
		{_parser->add_tmp_var($1);}
	| NAME EQUAL NAME
		{_parser->add_tmp_var($1,$3);}
	;

 signed_integer
 	: PLUS INT_NUMBER
 		{$$ = $2;}
 	| MINUS INT_NUMBER
 		{$2->symbol.insert(0,"-"); $$ = $2;}
 	| INT_NUMBER
 		{$$ = $1;}
 	;
 	
 signed_float
 	: PLUS FLOAT_NUMBER
 		{$$ = $2;}
 	| MINUS FLOAT_NUMBER
 		{$2->symbol.insert(0,"-"); $$ = $2;}
 	| FLOAT_NUMBER
 		{$$ = $1;}
 	;

 estimated_params 
	: ESTIMATED_PARAMS ';' {_parser->estimation_init();} estimated_list END
	;
	
 estimated_list 
	: estimated_list estimated_elem 
		{_parser->set_estimated_elements();}
	| estimated_elem 
		{_parser->set_estimated_elements();}
	;

 estimated_elem 
	: estimated_elem1 COMMA estimated_elem2 ';'
	;

 estimated_elem1 
	: STDERR NAME 
		{_parser->estim_params.type = 1;
		 _parser->estim_params.name = *$2;
		}
	| NAME
		{_parser->estim_params.type = 2;
		 _parser->estim_params.name = *$1;
		}
	| CORR NAME COMMA NAME
		{_parser->estim_params.type = 3;
		 _parser->estim_params.name = *$2;
		 _parser->estim_params.name2 = *$4;
		}
	;

 estimated_elem2 
	: prior COMMA estimated_elem3 
		{_parser->estim_params.prior=*$1;}
	| value COMMA prior COMMA estimated_elem3 
		{_parser->estim_params.init_val=*$1;
		 _parser->estim_params.prior=*$3;
		}
	| value COMMA value COMMA value COMMA prior COMMA estimated_elem3 
		{_parser->estim_params.init_val=*$1;
		 _parser->estim_params.low_bound=*$3;
		 _parser->estim_params.up_bound=*$5;
		 _parser->estim_params.prior=*$7;
		}
	| value 
		{_parser->estim_params.init_val=*$1;}
	| value COMMA value COMMA value 
		{_parser->estim_params.init_val=*$1;
		 _parser->estim_params.low_bound=*$3;
		 _parser->estim_params.up_bound=*$5;
		}
	;

 estimated_elem3 
 	: value COMMA value 
 		{_parser->estim_params.mean=*$1;
 		 _parser->estim_params.std=*$3;
 		}
	| value COMMA value COMMA value COMMA value 
		{_parser->estim_params.mean=*$1;
		 _parser->estim_params.std=*$3;
		 _parser->estim_params.p3=*$5;
		 _parser->estim_params.p4=*$7;
		}
	| value COMMA value COMMA value COMMA value COMMA value 
		{_parser->estim_params.mean=*$1;
		 _parser->estim_params.std=*$3;
		 _parser->estim_params.p3=*$5;
		 _parser->estim_params.p4=*$7;
		 _parser->estim_params.jscale=*$9;
		}
	;

 estimated_params_init: ESTIMATED_PARAMS_INIT ';' estimated_init_list END
                      ;

 estimated_init_list : estimated_init_list estimated_init_elem
                       {_parser->set_estimated_init_elements();}
                     | estimated_init_elem
                       {_parser->set_estimated_init_elements();}
                     ;

 estimated_init_elem : STDERR NAME COMMA value ';'
 		     		{_parser->estim_params.type = 1;
				 _parser->estim_params.name = *$2;
				 _parser->estim_params.init_val=*$4;
				}
                     | CORR NAME COMMA NAME COMMA value ';'
 		     		{_parser->estim_params.type = 3;
				 _parser->estim_params.name = *$2;
				 _parser->estim_params.name2 = *$4;
				 _parser->estim_params.init_val=*$6;
				}
                     | NAME COMMA value ';'
 		     		{_parser->estim_params.type = 2;
				 _parser->estim_params.name = *$1;
				 _parser->estim_params.init_val=*$3;
				}
                     ;

 estimated_params_bounds: ESTIMATED_PARAMS_BOUNDS ';' estimated_bounds_list END
                      ;

 estimated_bounds_list : estimated_bounds_list estimated_bounds_elem
                       {_parser->set_estimated_bounds_elements();}
                     | estimated_bounds_elem
                       {_parser->set_estimated_bounds_elements();}
                     ;

 estimated_bounds_elem : STDERR NAME COMMA value COMMA value ';'
 		     		{_parser->estim_params.type = 1;
				 _parser->estim_params.name = *$2;
				 _parser->estim_params.low_bound=*$4;
				 _parser->estim_params.up_bound=*$6;
				}
                     | CORR NAME COMMA NAME COMMA value COMMA value ';'
 		     		{_parser->estim_params.type = 3;
				 _parser->estim_params.name = *$2;
				 _parser->estim_params.name2 = *$4;
				 _parser->estim_params.low_bound=*$6;
				 _parser->estim_params.up_bound=*$8;
				}
                     | NAME COMMA value COMMA value ';'
 		     		{_parser->estim_params.type = 2;
				 _parser->estim_params.name = *$1;
				 _parser->estim_params.low_bound=*$3;
				 _parser->estim_params.up_bound=*$5;
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
 		{$2->symbol.insert(0,"-"); $$ = $2;}
        | MINUS FLOAT_NUMBER
 		{$2->symbol.insert(0,"-"); $$ = $2;}
	;

 
       
 estimation 
	: ESTIMATION ';'
		{_parser->run_estimation();}
	| ESTIMATION '(' estimation_options_list ')' ';' 
		{_parser->run_estimation();}
	| ESTIMATION tmp_var_list ';'
		{_parser->run_estimation();}
	| ESTIMATION '(' estimation_options_list ')' tmp_var_list ';' 
		{_parser->run_estimation();}
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
                   ;
	
 list_optim_option
 	: '\'' NAME '\'' COMMA '\'' NAME '\'' {_parser->optim_options($2,$6,2);}
	| '\'' NAME '\'' COMMA value {_parser->optim_options($2,$5,2);}
	;

 optim_options
	: list_optim_option
	| optim_options COMMA list_optim_option;
	;
	
 varobs 
 	: VAROBS tmp_var_list ';' 
 		{_parser->set_varobs();}
	;

 observation_trends
        : OBSERVATION_TRENDS ';' {_parser->set_trend_init();} trend_list END
	;

 trend_list 
        : trend_list trend_element
	| trend_element
	;

 trend_element :  NAME '(' expression {$$ = _parser->get_expression($3);} ')' ';' {_parser->set_trend_element($1,$4);}
               ;

 unit_root_vars : UNIT_ROOT_VARS tmp_var_list ';' {_parser->set_unit_root_vars();}
                ;

 optim_weights : OPTIM_WEIGHTS ';' {_parser->begin_optim_weights();} optim_weights_list END
               ;

 optim_weights_list : optim_weights_list NAME expression {$$=_parser->get_expression($3);} ';'  {_parser->set_optim_weights($2, $4);}
                    | optim_weights_list NAME COMMA NAME expression {$$=_parser->get_expression($3);} ';' {_parser->set_optim_weights($2, $4, $6);}
                    | NAME expression {$$=_parser->get_expression($2);} ';' {_parser->set_optim_weights($1, $3);}
                    | NAME COMMA NAME expression {$$=_parser->get_expression($3);} ';' {_parser->set_optim_weights($1, $3, $5);}
                    ;

 osr_params : OSR_PARAMS tmp_var_list ';' {_parser->set_osr_params();}
            ;

 osr : OSR ';' {_parser->run_osr();}
     | OSR '(' olr_options ')' ';' {_parser->run_osr();}
     | OSR tmp_var_list ';' {_parser->run_osr();}
     | OSR '(' olr_options ')' tmp_var_list ';' {_parser->run_osr();}
     ;
 
 olr : OLR ';' {_parser->run_olr();}
     | OLR '(' olr_options ')' ';' {_parser->run_olr();}
     | OLR tmp_var_list ';' {_parser->run_olr();}
     | OLR '(' olr_options ')' tmp_var_list ';' {_parser->run_olr();}
     ;
 
 olr_option : o_olr_beta
     | stoch_simul_options_list
     ;
 
 olr_options : olr_option
     | olr_options COMMA olr_option
     ;

 olr_inst : OLR_INST tmp_var_list ';' {_parser->set_olr_inst();}
          ;

 calib_var : CALIB_VAR ';' {_parser->begin_calib_var();} calib_var_list END
           ;

 calib_var_list : calib_var_list calib_arg1
                | calib_arg1
                ;

 calib_arg1 : NAME calib_arg2 EQUAL expression ';' {_parser->set_calib_var($1,$2,$4);}
            | NAME COMMA NAME calib_arg2 EQUAL expression ';' {_parser->set_calib_var($1,$3,$4,$6);}
            | AUTOCORR NAME '(' INT_NUMBER ')' calib_arg2 EQUAL expression ';' {_parser->set_calib_ac($2,$4,$6,$8);}
            ;

 calib_arg2 : {$$ = new dynare::Objects("1");}
            | '(' INT_NUMBER ')' {$$ = $2;}
            | '(' FLOAT_NUMBER ')' {$$ = $2;}
            ;

 calib : CALIB ';' {_parser->run_calib(0);}
       | CALIB '(' COVAR ')' ';' {_parser->run_calib(1);}
       ;

 dynatype : DYNATYPE '(' NAME ')'';' {_parser->run_dynatype($3);}
          | DYNATYPE '(' NAME ')' tmp_var_list ';' {_parser->run_dynatype($3);}
          | DYNATYPE NAME ';' {_parser->run_dynatype($2);}
          | DYNATYPE '(' NAME '.' NAME ')'';' {_parser->run_dynatype($3,$5);}
          | DYNATYPE '(' NAME '.' NAME ')' tmp_var_list ';' {_parser->run_dynatype($3,$5);}
          | DYNATYPE NAME '.' NAME ';' {_parser->run_dynatype($2,$4);};

 dynasave : DYNASAVE '(' NAME ')'';' {_parser->run_dynasave($3);}
          | DYNASAVE '(' NAME ')' tmp_var_list ';' {_parser->run_dynasave($3);}
          | DYNASAVE NAME ';' {_parser->run_dynasave($2);}
          | DYNASAVE '(' NAME '.' NAME ')'';' {_parser->run_dynasave($3,$5);}
          | DYNASAVE '(' NAME '.' NAME ')' tmp_var_list ';' {_parser->run_dynasave($3,$5);}
          | DYNASAVE NAME '.' NAME ';' {_parser->run_dynasave($2,$4);};

 model_comparison : MODEL_COMPARISON '(' model_comparison_options ')' {_parser->begin_model_comparison();} 
                       filename_list ';' {_parser->run_model_comparison();}
                  ;

 model_comparison_options: model_comparison_options ',' model_comparison_option
              | model_comparison_option
              ;

 model_comparison_option : o_model_comparison_approximation
              | o_print
              | o_noprint
              ;

 filename_list : filename {_parser->add_mc_filename($1);}
        | filename_list ',' filename {_parser->add_mc_filename($3);}
	| filename '(' value ')' {_parser->add_mc_filename($1,$3);}
        | filename_list ',' filename '(' value ')' {_parser->add_mc_filename($3,$5);}
        ;

 filename : filename_elem {$$=$1;}
        | filename filename_elem {$$ = _parser->cat($1,$2);}
        ;

 filename_elem : NAME
        | '\\' {$$ = new dynare::Objects("\\");}
        | '/' {$$ = new dynare::Objects("/");}
        | ':' {$$ = new dynare::Objects(":");}
        | '.' {$$ = new dynare::Objects(".");}
        ;

 o_dr_algo: DR_ALGO EQUAL INT_NUMBER {_parser->option_num("dr_algo",$3);};
 o_solve_algo: SOLVE_ALGO EQUAL INT_NUMBER {_parser->option_num("solve_algo",$3);};
 o_simul_algo: SIMUL_ALGO EQUAL INT_NUMBER {_parser->option_num("simul_algo",$3);};
 o_linear: LINEAR {_parser->option_num("linear","1");};
 o_order: ORDER EQUAL INT_NUMBER {_parser->option_num("order",$3);};
 o_replic: REPLIC EQUAL INT_NUMBER {_parser->option_num("replic",$3);};
 o_drop: DROP EQUAL INT_NUMBER {_parser->option_num("drop",$3);};
 o_ar: AR EQUAL INT_NUMBER {_parser->option_num("ar",$3);};
 o_nocorr: NOCORR {_parser->option_num("nocorr","1");};
 o_nofunction: NOFUNCTIONS {_parser->option_num("nofunctions","1");};
 o_nomoments: NOMOMENTS {_parser->option_num("nomoments","1");};
 o_irf: IRF EQUAL INT_NUMBER {_parser->option_num("irf",$3);};
 o_hp_filter: HP_FILTER EQUAL INT_NUMBER {_parser->option_num("hp_filter",$3);};
 o_hp_ngrid: HP_NGRID EQUAL INT_NUMBER {_parser->option_num("hp_ngrid",$3);};
 o_periods: PERIODS EQUAL INT_NUMBER {_parser->option_num("periods",$3);_parser->option_num("simul","1");};
 o_simul: SIMUL {_parser->option_num("simul","1");};
 o_simul_seed: SIMUL_SEED EQUAL INT_NUMBER { _parser->option_num("simul_seed",$3)};
 o_qz_criterium: QZ_CRITERIUM EQUAL INT_NUMBER { _parser->option_num("qz_criterium",$3)}
               | QZ_CRITERIUM EQUAL FLOAT_NUMBER { _parser->option_num("qz_criterium",$3)}
               ;
 o_datafile: DATAFILE EQUAL NAME {_parser->option_str("datafile",$3);};   
 o_nobs: NOBS EQUAL vec_int {_parser->option_num("nobs",$3);}
       | NOBS EQUAL INT_NUMBER {_parser->option_num("nobs",$3);}
       ;	      
 o_first_obs: FIRST_OBS EQUAL INT_NUMBER {_parser->option_num("first_obs",$3);};
 o_prefilter: PREFILTER EQUAL INT_NUMBER {_parser->option_num("prefilter",$3);};
 o_presample: PRESAMPLE EQUAL INT_NUMBER {_parser->option_num("presample",$3);};
 o_lik_algo: LIK_ALGO EQUAL INT_NUMBER {_parser->option_num("lik_algo",$3);}; 
 o_lik_init: LIK_INIT EQUAL INT_NUMBER {_parser->option_num("lik_init",$3);}; 
 o_nograph: NOGRAPH {_parser->option_num("nograph","1");}; 
 o_conf_sig: CONF_SIG EQUAL FLOAT_NUMBER {_parser->option_num("conf_sig",$3);}; 
 o_mh_replic: MH_REPLIC EQUAL INT_NUMBER {_parser->option_num("mh_replic",$3);}; 
 o_mh_drop: MH_DROP EQUAL FLOAT_NUMBER {_parser->option_num("mh_drop",$3);}; 
 o_mh_jscale: MH_JSCALE EQUAL FLOAT_NUMBER {_parser->option_num("mh_jscale",$3);}; 
 o_optim: OPTIM  EQUAL '(' optim_options ')';
 o_mh_init_scale :MH_INIT_SCALE EQUAL FLOAT_NUMBER {_parser->option_num("mh_init_scale",$3);};
 o_mh_init_scale :MH_INIT_SCALE EQUAL INT_NUMBER {_parser->option_num("mh_init_scale",$3);};
 o_mode_file : MODE_FILE EQUAL NAME {_parser->option_str("mode_file",$3);};
 o_mode_compute : MODE_COMPUTE EQUAL INT_NUMBER {_parser->option_num("mode_compute",$3);};
 o_mode_check : MODE_CHECK {_parser->option_num("mode_check","1");};
 o_prior_trunc : PRIOR_TRUNC EQUAL FLOAT_NUMBER {_parser->option_num("prior_trunc",$3);};
 o_mh_mode : MH_MODE EQUAL INT_NUMBER {_parser->option_num("mh_mode",$3);};
 o_mh_nblcks : MH_NBLOCKS EQUAL INT_NUMBER {_parser->option_num("mh_nblck",$3);};
 o_load_mh_file : LOAD_MH_FILE {_parser->option_num("load_mh_file","1");};
 o_loglinear : LOGLINEAR {_parser->option_num("loglinear","1");};
 o_nodiagnostic : NODIAGNOSTIC {_parser->option_num("nodiagnostic","1");};
 o_bayesian_irf : BAYESIAN_IRF {_parser->option_num("bayesian_irf","1");};
 o_tex : TEX {_parser->option_num("TeX","1");};
 o_forecast : FORECAST EQUAL INT_NUMBER {_parser->option_num("forecast",$3);};
 o_smoother : SMOOTHER {_parser->option_num("smoother","1");};
 o_moments_varendo : MOMENTS_VARENDO {_parser->option_num("moments_varendo","1");};
 o_filtered_vars : FILTERED_VARS {_parser->option_num("filtered_vars","1");};
 o_relative_irf : RELATIVE_IRF {_parser->option_num("relative_irf","1");};
 o_kalman_algo : KALMAN_ALGO EQUAL INT_NUMBER {_parser->option_num("kalman_algo",$3);};
 o_kalman_tol : KALMAN_TOL EQUAL INT_NUMBER {_parser->option_num("kalman_tol",$3);};
 o_olr_beta: OLR_BETA EQUAL value {_parser->option_num("olr_beta",$3);};
 o_model_comparison_approximation: MODEL_COMPARISON_APPROXIMATION EQUAL LAPLACE {dynare::Objects* tmp = new dynare::Objects("Laplace");_parser->option_str("model_comparison_approximation",tmp);}
   | MODEL_COMPARISON_APPROXIMATION EQUAL MODIFIEDHARMONICMEAN {dynare::Objects* tmp = new dynare::Objects("MODIFIEDHARMONICMEAN");_parser->option_str("model_comparison_approximation",tmp);}
   ;
 o_print : PRINT {_parser->option_num("noprint","0");};
 o_noprint : NOPRINT {_parser->option_num("noprint","1");};
 o_nograph : GRAPH {_parser->option_num("nograph","0");};	
 o_xls_sheet : XLS_SHEET EQUAL NAME {_parser->option_str("xls_sheet",$3);} 
 o_xls_range : XLS_RANGE EQUAL range {_parser->option_str("xls_range",$3);} 
 o_filter_step_ahead : FILTER_STEP_AHEAD EQUAL vec_int {_parser->option_num("filter_step_ahead",$3);}

 range : NAME ':' NAME {$$ = new dynare::Objects(":");$$ = _parser->cat($1,$$);$$ = _parser->cat($$,$3);}
 vec_int_elem : INT_NUMBER
 	      | INT_NUMBER ':' {$$ = new dynare::Objects(":"); $$ = _parser->cat($1,$$); }
	        INT_NUMBER {$$ = _parser->cat($3,$4);}   		    	        
	      ;

 vec_int_1 : '[' vec_int_elem {$$ = new dynare::Objects("["); $$ =_parser->cat($$,$2);}
           | vec_int_1 vec_int_elem {$$ = _parser->cat_with_space($1,$2);}
           ;

 vec_int : vec_int_1 ']' {$$ = new dynare::Objects("]"); $$ = _parser->cat($1,$$);};

 %%
