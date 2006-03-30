%{
#define YY_BUF_SIZE 1000000
#include <unistd.h>
#include <string.h>
#include "DynareScanner.h" 
#ifdef HAVE_CONFIG_H
# include "config.hh"
#endif
#define YLMM_SCANNER_CLASS dynare::scanner
#define LEXDEBUG 1
#define FINISH 0
#define YY_READ_BUF_SIZE 1000000
#include "ylmm/lexmm.hh"
int lineno = 1; 
int comment_caller;
/* Particular value : when sigma_e command is found 
 this flag is set to 1, when command finished it is set to 0
 */
int sigma_e = 0;
%} 


%option stack
%option yylineno
%x COMMENT
%x SET_STATEMENT
%x DYNARE_STATEMENT
%x DYNARE_BLOCK
%x NATIVE

%%
 /* spaces, tabs and EOL are ignored */
<*>[ \t\n\r\f] ;

 /* Comments */
<INITIAL,SET_STATEMENT,DYNARE_STATEMENT,DYNARE_BLOCK>["%"].*
<INITIAL,SET_STATEMENT,DYNARE_STATEMENT,DYNARE_BLOCK>["/"]["/"].* 
<INITIAL,SET_STATEMENT,DYNARE_STATEMENT,DYNARE_BLOCK>"/*"   {comment_caller = YYSTATE; BEGIN COMMENT;}

<COMMENT>[^*\n]* 		
<COMMENT>"*"+[^*/\n]*	
<COMMENT>"*"+"/"        {BEGIN comment_caller;}

 /* Begin of a Dynare statement */
<INITIAL>var {BEGIN DYNARE_STATEMENT; return VAR;}
<INITIAL>varexo {BEGIN DYNARE_STATEMENT; return VAREXO;}
<INITIAL>varexo_det {BEGIN DYNARE_STATEMENT; return VAREXO_DET;}
<INITIAL>parameters {BEGIN DYNARE_STATEMENT; return PARAMETERS;}
<INITIAL>periods 	{BEGIN DYNARE_STATEMENT; return PERIODS;}
<INITIAL>initvalf 	{BEGIN DYNARE_STATEMENT; return INITVALF;}
<INITIAL>estimation {BEGIN DYNARE_STATEMENT; return ESTIMATION;}
<INITIAL>varobs 	{BEGIN DYNARE_STATEMENT; return VAROBS;}
<INITIAL>unit_root_vars	{BEGIN DYNARE_STATEMENT; return UNIT_ROOT_VARS;}
<INITIAL>dyn2vec 	{BEGIN DYNARE_STATEMENT; return DYN2VEC;}
<INITIAL>rplot	 	{BEGIN DYNARE_STATEMENT; return RPLOT;}
<INITIAL>osr_params 	{BEGIN DYNARE_STATEMENT; return OSR_PARAMS;}
<INITIAL>osr	 	{BEGIN DYNARE_STATEMENT; return OSR;}
<INITIAL>calib_var 	{BEGIN DYNARE_STATEMENT; return CALIB_VAR;}
<INITIAL>dynatype	{BEGIN DYNARE_STATEMENT; return DYNATYPE;}
<INITIAL>dynasave 	{BEGIN DYNARE_STATEMENT; return DYNASAVE;}
<INITIAL>olr	 	{BEGIN DYNARE_STATEMENT; return OLR;}
<INITIAL>olr_inst	 	{BEGIN DYNARE_STATEMENT; return OLR_INST;}
<INITIAL>model_comparison 	{BEGIN DYNARE_STATEMENT; return MODEL_COMPARISON;}

<INITIAL>steady {BEGIN DYNARE_STATEMENT; return STEADY;}
<INITIAL>check {BEGIN DYNARE_STATEMENT; return CHECK;}
<INITIAL>simul {BEGIN DYNARE_STATEMENT; return SIMUL;}
<INITIAL>stoch_simul {BEGIN DYNARE_STATEMENT; return STOCH_SIMUL;}
<INITIAL>dsample {BEGIN DYNARE_STATEMENT; return DSAMPLE;}
<INITIAL>Sigma_e {BEGIN DYNARE_STATEMENT;sigma_e = 1; return SIGMA_E;}

 /* End of a Dynare statement */
<DYNARE_STATEMENT>; {if (!sigma_e) BEGIN INITIAL; return yytext[0];}


 /* Begin of a Dynare block */
<INITIAL>model {BEGIN DYNARE_BLOCK;return MODEL;}
<INITIAL>initval {BEGIN DYNARE_BLOCK;return INITVAL;}
<INITIAL>endval {BEGIN DYNARE_BLOCK;return ENDVAL;}
<INITIAL>histval {BEGIN DYNARE_BLOCK;return HISTVAL;}
<INITIAL>shocks {BEGIN DYNARE_BLOCK;return SHOCKS;}
<INITIAL>estimated_params {BEGIN DYNARE_BLOCK;return ESTIMATED_PARAMS;}
<INITIAL>estimated_params_init 		{BEGIN DYNARE_BLOCK;return ESTIMATED_PARAMS_INIT;}
<INITIAL>estimated_params_bounds 	{BEGIN DYNARE_BLOCK;return ESTIMATED_PARAMS_BOUNDS;}
<INITIAL>observation_trends {BEGIN DYNARE_BLOCK;return OBSERVATION_TRENDS;}
<INITIAL>optim_weights {BEGIN DYNARE_BLOCK;return OPTIM_WEIGHTS;}

 /* End of a Dynare block */
<DYNARE_BLOCK>end[ \t\n]*; 	{BEGIN INITIAL;return END;}   

 /* Inside  of a Dynare statement */
<DYNARE_STATEMENT>datafile 		{return DATAFILE;}
<DYNARE_STATEMENT>nobs 			{return NOBS;}
<DYNARE_STATEMENT>first_obs 		{return FIRST_OBS;}
<DYNARE_STATEMENT>prefilter 		{return PREFILTER;} 
<DYNARE_STATEMENT>presample 		{return PRESAMPLE;} 
<DYNARE_STATEMENT>lik_algo  		{return LIK_ALGO;}  
<DYNARE_STATEMENT>lik_init  		{return LIK_INIT;}  
<DYNARE_STATEMENT>graph   		{return GRAPH;}  	  
<DYNARE_STATEMENT>nograph   		{return NOGRAPH;}  	  
<DYNARE_STATEMENT>print   		{return PRINT;}  	  
<DYNARE_STATEMENT>noprint   		{return NOPRINT;}  	  
<DYNARE_STATEMENT>conf_sig  		{return CONF_SIG;}  
<DYNARE_STATEMENT>mh_replic 		{return MH_REPLIC;} 
<DYNARE_STATEMENT>mh_drop   		{return MH_DROP;}   
<DYNARE_STATEMENT>mh_jscale   		{return MH_JSCALE;}   
<DYNARE_STATEMENT>mh_init_scale 	{return MH_INIT_SCALE;}
<DYNARE_STATEMENT>mode_file 		{return MODE_FILE;}
<DYNARE_STATEMENT>mode_compute 	{return MODE_COMPUTE;}
<DYNARE_STATEMENT>mode_check 		{return MODE_CHECK;}
<DYNARE_STATEMENT>prior_trunc 	{return PRIOR_TRUNC;}
<DYNARE_STATEMENT>mh_mode 		{return MH_MODE;}
<DYNARE_STATEMENT>mh_nblocks 		{return MH_NBLOCKS;}
<DYNARE_STATEMENT>load_mh_file 	{return LOAD_MH_FILE;}
<DYNARE_STATEMENT>loglinear 		{return LOGLINEAR;}
<DYNARE_STATEMENT>nodiagnostic 	{return NODIAGNOSTIC;}
<DYNARE_STATEMENT>kalman_algo 	{return KALMAN_ALGO;}
<DYNARE_STATEMENT>kalman_tol 	{return KALMAN_TOL;}
<DYNARE_STATEMENT>forecast 	{return FORECAST;}
<DYNARE_STATEMENT>smoother 	{return SMOOTHER;}
<DYNARE_STATEMENT>bayesian_irf 	{return BAYESIAN_IRF;}
<DYNARE_STATEMENT>moments_varendo {return MOMENTS_VARENDO;}
<DYNARE_STATEMENT>filtered_vars	{return FILTERED_VARS;}
<DYNARE_STATEMENT>filter_step_ahead	{return FILTER_STEP_AHEAD;}
<DYNARE_STATEMENT>relative_irf 	{return RELATIVE_IRF;}
<DYNARE_STATEMENT>tex		{return TEX;}
<DYNARE_STATEMENT>shock_size	{return SHOCK_SIZE;}
<DYNARE_STATEMENT>moments	{return MOMENTS;}
<DYNARE_STATEMENT>nomoments	{return NOMOMENTS;}
<DYNARE_STATEMENT>corr		{return CORR;}
<DYNARE_STATEMENT>nocorr	{return NOCORR;}
<DYNARE_STATEMENT>optim		{return OPTIM;}
<DYNARE_STATEMENT>periods	{return PERIODS;}
<DYNARE_STATEMENT>diffuse_d {return DIFFUSE_D;}
<DYNARE_STATEMENT>model_comparison_approximation {return MODEL_COMPARISON;}
<DYNARE_STATEMENT>laplace {return LAPLACE;}
<DYNARE_STATEMENT>modifiedharmonicmean {return MODIFIEDHARMONICMEAN;}

<DYNARE_STATEMENT>[\$][^$]*[\$] {
	strtok(yytext+1,"$");
	_scanner->do_name(yytext+1); 
	return TEX_NAME;}

 /* Inside a Dynare block */
<DYNARE_BLOCK>var {return VAR;}
<DYNARE_BLOCK>stderr {return STDERR;}
<DYNARE_BLOCK>values {return VALUES;}
<DYNARE_BLOCK>corr {return CORR;}
<DYNARE_BLOCK>periods {return PERIODS;}

<DYNARE_BLOCK>gamma_pdf {return GAMMA_PDF;}
<DYNARE_BLOCK>beta_pdf {return BETA_PDF;}
<DYNARE_BLOCK>normal_pdf {return NORMAL_PDF;}
<DYNARE_BLOCK>inv_gamma_pdf {return INV_GAMMA_PDF;}
<DYNARE_BLOCK>inv_gamma1_pdf {return INV_GAMMA_PDF;}
<DYNARE_BLOCK>inv_gamma2_pdf {return INV_GAMMA_PDF;}
<DYNARE_BLOCK>uniform_pdf {return UNIFORM_PDF;}

<DYNARE_BLOCK>; {return yytext[0];}


 /* Inside Dynare statement */
<DYNARE_STATEMENT>solve_algo {return SOLVE_ALGO;}
<DYNARE_STATEMENT>dr_algo {return DR_ALGO;}
<DYNARE_STATEMENT>simul_algo {return SIMUL_ALGO;}
<DYNARE_STATEMENT>drop {return DROP;}
<DYNARE_STATEMENT>order {return ORDER;}
<DYNARE_STATEMENT>replic {return REPLIC;}
<DYNARE_STATEMENT>ar {return AR;}
<DYNARE_STATEMENT>nofunctions {return NOFUNCTIONS;}
<DYNARE_STATEMENT>irf {return IRF;}
<DYNARE_STATEMENT>hp_filter {return HP_FILTER;}
<DYNARE_STATEMENT>hp_ngrid {return HP_NGRID;}
<DYNARE_STATEMENT>simul_seed {return SIMUL_SEED;}
<DYNARE_STATEMENT>qz_criterium {return QZ_CRITERIUM;}
<DYNARE_STATEMENT>simul {return SIMUL;}
<DYNARE_STATEMENT>autocorr {return AUTOCORR;}
<DYNARE_STATEMENT>olr_beta {return OLR_BETA;}
<DYNARE_STATEMENT>xtick   		{return XTICK;}  	  
<DYNARE_STATEMENT>xticklabel   		{return XTICKLABEL;}  	  
<DYNARE_STATEMENT>xls_sheet {return XLS_SHEET;}
<DYNARE_STATEMENT>xls_range {return XLS_RANGE;}

<DYNARE_STATEMENT,DYNARE_BLOCK>use_dll {return USE_DLL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>linear {return LINEAR;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[,] {_scanner->do_operator(COMMA); return COMMA;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\(\)] {return yytext[0];} 
<DYNARE_STATEMENT,DYNARE_BLOCK>[\[] {return yytext[0];}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\]] {if (sigma_e) sigma_e=0; return yytext[0];}
<DYNARE_STATEMENT,DYNARE_BLOCK>[+] {_scanner->do_operator(PLUS); return PLUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[-] {_scanner->do_operator(MINUS);return MINUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[*] {_scanner->do_operator(TIMES);return TIMES;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[/] {_scanner->do_operator(DIVIDE);return DIVIDE;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[=] {_scanner->do_operator(EQUAL); return EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\^] {_scanner->do_operator(POWER);return POWER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[!] {_scanner->do_operator(FACTORIAL);return FACTORIAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>exp {_scanner->do_operator(EXP);return EXP;}
<DYNARE_STATEMENT,DYNARE_BLOCK>log {_scanner->do_operator(LOG);return LOG;}
<DYNARE_STATEMENT,DYNARE_BLOCK>log10 {_scanner->do_operator(LOG10);return LOG10;}
<DYNARE_STATEMENT,DYNARE_BLOCK>ln {_scanner->do_operator(LN);return LN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sin {_scanner->do_operator(SIN);return SIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>cos {_scanner->do_operator(COS);return COS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>tan {_scanner->do_operator(TAN);return TAN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>asin {_scanner->do_operator(ASIN);return ASIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>acos {_scanner->do_operator(ACOS);return ACOS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>atan {_scanner->do_operator(ATAN);return ATAN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sinh {_scanner->do_operator(SINH);return SINH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>cosh {_scanner->do_operator(COSH);return COSH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>tanh {_scanner->do_operator(TANH);return TANH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>asinh {_scanner->do_operator(ASINH);return ASINH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>acosh {_scanner->do_operator(ACOSH);return ACOSH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>atanh {_scanner->do_operator(ATANH);return ATANH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sqrt  {_scanner->do_operator(SQRT);return SQRT;}

<DYNARE_STATEMENT,DYNARE_BLOCK>[A-Za-z_][A-Za-z0-9_]* {
	_scanner->do_name(yytext); 
	return NAME;}

<DYNARE_STATEMENT,DYNARE_BLOCK>((([0-9]*\.[0-9]+)|([0-9]+\.))([edED][-+]?[0-9]+)?)|([0-9]+[edED][-+]?[0-9]+) {
	_scanner->do_num_constant(yytext);
	 return FLOAT_NUMBER;}

<DYNARE_STATEMENT,DYNARE_BLOCK>[0-9]+ {
	_scanner->do_num_constant(yytext);
   	return INT_NUMBER;}
   	
 /* an instruction starting with a recognized symbol is passed as NAME,
    otherwise it is a native statement until the end of the line
 */
<INITIAL>[A-Za-z_][A-Za-z0-9_]* {	      
		if (SymbolTable::getID(yytext) != -1)
		{
		 	BEGIN DYNARE_STATEMENT;
	 		_scanner->do_name(yytext);
			return NAME;
	 	}
	 	else
	 	{
	 		BEGIN NATIVE;
	 		_scanner->do_as_is(yytext);
	 	}
 }
<INITIAL>. {BEGIN NATIVE; _scanner->do_as_is(yytext);}

 /* NATIVE Block */
<NATIVE>.* {BEGIN INITIAL;_scanner->do_as_is(yytext);_scanner->do_as_is("\n");}

