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
#include "yaccmm.h"
#ifdef yyerror
#undef yyerror
#define yyerror _parser->error 
#endif
%}
/* %pure_parser */

%token SIGMA_E INITVAL ENDVAL HISTVAL SHOCKS  VALUES STDERR CORR MODEL  END
%token VAR VAREXO VAREXO_DET PARAMETER PERIODS
%token  NAME TEX_NAME INT_NUMBER FLOAT_NUMBER EQUAL
%token DR_ALGO SOLVE_ALGO STOCH_SIMUL STEADY SIMUL CHECK
%token ORDER NOFUNCTIONS REPLIC DROP IRF AR LINEAR HP_FILTER NOCORR SIMUL_SEED USE_DLL
%token HP_NGRID HP_NGRIDDR_ALGO SIMUL_ALGO QZ_CRITERIUM NOMOMENTS
%token VAROBS ESTIMATION OBSERVATION_TRENDS UNIT_ROOT_VARS
%token PERC DATAFILE NOBS FIRST_OBS
%token MH_REPLIC MH_DROP
%token MH_JSCALE OPTIM MH_INIT_SCALE MODE_FILE MODE_COMPUTE MODE_CHECK
%token PRIOR_TRUNC MH_MODE MH_NBLOCKS LOAD_MH_FILE LOGLINEAR
%token NODIAGNOSTIC NOGRAPH PRESAMPLE CONF_SIG LIK_INIT PREFILTER LIK_ALGO
%token ESTIMATED_PARAMS GAMMA_PDF BETA_PDF NORMAL_PDF INV_GAMMA_PDF UNIFORM_PDF
%token INV_GAMMA1_PDF INV_GAMMA2_PDF
%token KALMAN_ALGO KALMAN_TOL FORECAST SMOOTHER BAYESIAN_IRF MOMENTS_VARENDO
%token FILTERED_VARS RELATIVE_IRF TEX
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
 	| sigma_e
 	| steady
 	| check
 	| simul
 	| stoch_simul
 	| estimation
	| estimated_params
	| varobs
	| observation_trends
	| unit_root_vars
	  /*
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
	  */
	;

    
 declaration
 	: parameters
 	| var
 	| varexo
 	| varexo_det
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
 	: PARAMETER parameter_list ';'
 	;
 
 var_list
 	: var_list NAME  
 		{$$  = _parser->add_endogenous($2);}
 	| var_list ',' NAME  
 		{$$  = _parser->add_endogenous($3);}
 	| NAME
 		{$$  = _parser->add_endogenous($1);}
 	| var_list NAME  TEX_NAME 
 		{$$  = _parser->add_endogenous($2,$3);}
 	| var_list ',' NAME  TEX_NAME
 		{$$  = _parser->add_endogenous($3,$4);}
 	| NAME TEX_NAME
 		{$$  = _parser->add_endogenous($1,$2);}
	; 
	
 varexo_list
 	: varexo_list NAME
 		{$$  = _parser->add_exogenous($2);}              
 	| varexo_list ',' NAME
 		{$$  = _parser->add_exogenous($3);}              
 	| NAME
 		{$$  = _parser->add_exogenous($1);}
 	| varexo_list NAME TEX_NAME
 		{$$  = _parser->add_exogenous($2, $3);}              
 	| varexo_list ',' NAME TEX_NAME
 		{$$  = _parser->add_exogenous($3, $4);}              
 	| NAME TEX_NAME
 		{$$  = _parser->add_exogenous($1, $2);}
	; 

 varexo_det_list
 	: varexo_det_list NAME
 		{$$  = _parser->add_exogenous_det($2);}              
 	| varexo_det_list ',' NAME
 		{$$  = _parser->add_exogenous_det($3);}              
 	| NAME
 		{$$  = _parser->add_exogenous_det($1);}
 	| varexo_det_list NAME TEX_NAME
 		{$$  = _parser->add_exogenous_det($2, $3);}              
 	| varexo_det_list ',' NAME TEX_NAME
 		{$$  = _parser->add_exogenous_det($3, $4);}              
 	| NAME TEX_NAME
 		{$$  = _parser->add_exogenous_det($1, $2);}
	; 

 parameter_list
 	: parameter_list NAME
 		{$$  = _parser->add_parameter($2);}
 	| parameter_list ',' NAME
 		{$$  = _parser->add_parameter($3);}
 	| NAME
 		{$$  = _parser->add_parameter($1);}
 	| parameter_list NAME TEX_NAME
 		{$$  = _parser->add_parameter($2, $3);}
 	| parameter_list ',' NAME TEX_NAME
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
	; 
	
 initval
 	: INITVAL ';' {_parser->begin_initval();} initval_list END ';'
 		{_parser->end_initval();}
 	;
 	
 endval
 	: ENDVAL ';' {_parser->begin_endval();} initval_list END ';'
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
 	: HISTVAL ';' {_parser->begin_histval();} histval_list END ';' 

 histval_list
 	: histval_list histval_elem
	| histval_elem
	;
	
 histval_elem
 	: NAME '(' signed_integer ')' EQUAL expression ';'
 		{_parser->hist_val($1, $3, $6);}
 	;
	
 model
 	: MODEL ';' equation_list END ';' 
 		{_parser->check_model();}
 	| MODEL '(' LINEAR ')' ';' {_parser->option_num("linear","1");} equation_list END ';'
 		{_parser->check_model();}
 	| MODEL '(' USE_DLL ')' ';' {_parser->use_dll();} equation_list END ';'
 		{_parser->check_model();}
 	;

 equation_list
    : equation_list equation     	
    | equation
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
	
 model_var
 	: NAME 
 		{$$ = _parser->add_variable($1);}
	| NAME '(' signed_integer ')'
		{$$ = _parser->add_variable($1,$3);}
	;
	
 shocks
 	: SHOCKS  ';' {_parser->begin_shocks();} shock_list END ';'
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
	| VAR NAME ',' NAME EQUAL expression {$$ = _parser->get_expression($6);} ';'
		{_parser->add_covar_shock($2, $4, $7);}
	| CORR NAME ',' NAME EQUAL expression {$$ = _parser->get_expression($6);} ';'
		{_parser->add_correl_shock($2, $4, $7);}
	;

 period_list
 	: period_list INT_NUMBER
 		{_parser->add_period($2);}
	| period_list INT_NUMBER ':' INT_NUMBER 
 		{_parser->add_period($2,$4);}
        | period_list ',' INT_NUMBER
 		{_parser->add_period($3);}
	| period_list ',' INT_NUMBER ':' INT_NUMBER 
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
 	: triangular_row ',' '(' expression {$$ = _parser->get_expression($4);} ')' 
 		{_parser->add_to_row($5);}
	| triangular_row ',' FLOAT_NUMBER 
		{_parser->add_to_row($3);}
	| triangular_row ',' INT_NUMBER 
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
 			_parser->option_num("solve_algo","0");
 			_parser->steady();
 		}
    | STEADY '(' steady_options_list ')' ';'
    	 {_parser->steady();}
    ;
    
 steady_options_list : steady_options_list ',' steady_options
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

 check_options_list : check_options_list ',' check_options
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

 simul_options_list: simul_options_list ',' simul_options
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
	
 stoch_simul_options_list: stoch_simul_options_list ',' stoch_simul_options
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
	| tmp_var_list ',' NAME
		{_parser->add_tmp_var($3);}
	| tmp_var_list ',' NAME EQUAL NAME
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
	: ESTIMATED_PARAMS ';' {_parser->estimation_init();} estimated_list END ';'
	;
	
 estimated_list 
	: estimated_list estimated_elem 
		{_parser->set_estimated_elements();}
	| estimated_elem 
		{_parser->set_estimated_elements();}
	;

 estimated_elem 
	: estimated_elem1 ',' estimated_elem2 ';'
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
	;

 estimated_elem2 
	: prior ',' estimated_elem3 
		{_parser->estim_params.prior=*$1;}
	| value ',' prior ',' estimated_elem3 
		{_parser->estim_params.init_val=*$1;
		 _parser->estim_params.prior=*$3;
		}
	| value ',' value ',' value ',' prior ',' estimated_elem3 
		{_parser->estim_params.init_val=*$1;
		 _parser->estim_params.low_bound=*$3;
		 _parser->estim_params.up_bound=*$5;
		 _parser->estim_params.prior=*$7;
		}
	| value 
		{_parser->estim_params.init_val=*$1;}
	| value ',' value ',' value 
		{_parser->estim_params.init_val=*$1;
		 _parser->estim_params.low_bound=*$3;
		 _parser->estim_params.up_bound=*$5;
		}
	;

 estimated_elem3 
 	: value ',' value 
 		{_parser->estim_params.mean=*$1;
 		 _parser->estim_params.std=*$3;
 		}
	| value ',' value ',' value ',' value 
		{_parser->estim_params.mean=*$1;
		 _parser->estim_params.std=*$3;
		 _parser->estim_params.p3=*$5;
		 _parser->estim_params.p4=*$7;
		}
	| value ',' value ',' value ',' value ',' value 
		{_parser->estim_params.mean=*$1;
		 _parser->estim_params.std=*$3;
		 _parser->estim_params.p3=*$5;
		 _parser->estim_params.p4=*$7;
		 _parser->estim_params.jscale=*$9;
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
	: INT_NUMBER
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
	: estimation_options_list ',' estimation_options
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
                   | o_tex
                   | o_forecast
                   | o_smoother
                   | o_moments_varendo
                   | o_filtered_vars
                   | o_kalman_algo
                   | o_kalman_tol
                   ;
	
 list_optim_option
 	: '\'' NAME '\'' ',' '\'' NAME '\'' {_parser->optim_options($2,$6,2);}
	| '\'' NAME '\'' ',' value {_parser->optim_options($2,$5,2);}
	;

 optim_options
	: list_optim_option
	| optim_options ',' list_optim_option;
	;
	
 varobs 
 	: VAROBS tmp_var_list ';' 
 		{_parser->set_varobs();}
	;

 observation_trends
        : OBSERVATION_TRENDS ';' {_parser->set_trend_init();} trend_list END ';'
	;

 trend_list 
        : trend_list trend_element
	| trend_element
	;

 trend_element :  NAME '(' expression {$$ = _parser->get_expression($3);} ')' ';' {_parser->set_trend_element($1,$4);}
               ;

 unit_root_vars : UNIT_ROOT_VARS tmp_var_list ';' {_parser->set_unit_root_vars();}
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
 o_nobs: NOBS EQUAL INT_NUMBER {_parser->option_num("nobs",$3);};	      
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
 %%
 
