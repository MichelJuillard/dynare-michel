%{
using namespace std;

/* The prototype for yylex() is defined in ParsingDriver.hh */

#include "ParsingDriver.hh"
#include "DynareBison.hh"

/* Shortcut to access tokens defined by Bison */
typedef yy::parser::token token;

/* By default yylex returns int, we use token_type.
   Unfortunately yyterminate by default returns 0, which is
   not of token_type.  */
#define yyterminate() return yy::parser::token_type (0);

int comment_caller;
/* Particular value : when sigma_e command is found 
 this flag is set to 1, when command finished it is set to 0
 */
int sigma_e = 0;
%} 

%option case-insensitive noyywrap nounput batch debug never-interactive

/* To be removed in the future (see static ParsingDriver::error) */
%option yylineno

%x COMMENT
%x SET_STATEMENT
%x DYNARE_STATEMENT
%x DYNARE_BLOCK
%x NATIVE

%{
// Increments location counter for every token read
#define YY_USER_ACTION yylloc->columns(yyleng);
%}
%%
%{
  // Reset location before reading token
  yylloc->step();
%}

 /* spaces, tabs and EOL are ignored */
<*>[ \t\r\f]+  { yylloc->step(); }
<*>[\n]+       { yylloc->lines(yyleng); yylloc->step(); }

 /* Comments */
<INITIAL,SET_STATEMENT,DYNARE_STATEMENT,DYNARE_BLOCK>["%"].*
<INITIAL,SET_STATEMENT,DYNARE_STATEMENT,DYNARE_BLOCK>["/"]["/"].* 
<INITIAL,SET_STATEMENT,DYNARE_STATEMENT,DYNARE_BLOCK>"/*"   {comment_caller = YY_START; BEGIN COMMENT;}

<COMMENT>[^*\n]* 		
<COMMENT>"*"+[^/\n]
<COMMENT>"*"+"/"        {BEGIN comment_caller;}

 /* Begin of a Dynare statement */
<INITIAL>var {BEGIN DYNARE_STATEMENT; return token::VAR;}
<INITIAL>varexo {BEGIN DYNARE_STATEMENT; return token::VAREXO;}
<INITIAL>varexo_det {BEGIN DYNARE_STATEMENT; return token::VAREXO_DET;}
<INITIAL>parameters {BEGIN DYNARE_STATEMENT; return token::PARAMETERS;}
<INITIAL>periods 	{BEGIN DYNARE_STATEMENT; return token::PERIODS;}
<INITIAL>initvalf 	{BEGIN DYNARE_STATEMENT; return token::INITVALF;}
<INITIAL>estimation {BEGIN DYNARE_STATEMENT; return token::ESTIMATION;}
<INITIAL>varobs 	{BEGIN DYNARE_STATEMENT; return token::VAROBS;}
<INITIAL>unit_root_vars	{BEGIN DYNARE_STATEMENT; return token::UNIT_ROOT_VARS;}
<INITIAL>dyn2vec 	{BEGIN DYNARE_STATEMENT; return token::DYN2VEC;}
<INITIAL>rplot	 	{BEGIN DYNARE_STATEMENT; return token::RPLOT;}
<INITIAL>osr_params 	{BEGIN DYNARE_STATEMENT; return token::OSR_PARAMS;}
<INITIAL>osr	 	{BEGIN DYNARE_STATEMENT; return token::OSR;}
<INITIAL>calib_var 	{BEGIN DYNARE_STATEMENT; return token::CALIB_VAR;}
<INITIAL>dynatype	{BEGIN DYNARE_STATEMENT; return token::DYNATYPE;}
<INITIAL>dynasave 	{BEGIN DYNARE_STATEMENT; return token::DYNASAVE;}
<INITIAL>olr	 	{BEGIN DYNARE_STATEMENT; return token::OLR;}
<INITIAL>olr_inst	 	{BEGIN DYNARE_STATEMENT; return token::OLR_INST;}
<INITIAL>model_comparison 	{BEGIN DYNARE_STATEMENT; return token::MODEL_COMPARISON;}

<INITIAL>steady {BEGIN DYNARE_STATEMENT; return token::STEADY;}
<INITIAL>check {BEGIN DYNARE_STATEMENT; return token::CHECK;}
<INITIAL>simul {BEGIN DYNARE_STATEMENT; return token::SIMUL;}
<INITIAL>stoch_simul {BEGIN DYNARE_STATEMENT; return token::STOCH_SIMUL;}
<INITIAL>dsample {BEGIN DYNARE_STATEMENT; return token::DSAMPLE;}
<INITIAL>Sigma_e {BEGIN DYNARE_STATEMENT; sigma_e = 1; return token::SIGMA_E;}

 /* End of a Dynare statement */
<DYNARE_STATEMENT>; {
  if (!sigma_e)
    BEGIN INITIAL;
  return yy::parser::token_type (yytext[0]);
}


 /* Begin of a Dynare block */
<INITIAL>model {BEGIN DYNARE_BLOCK; return token::MODEL;}
<INITIAL>initval {BEGIN DYNARE_BLOCK; return token::INITVAL;}
<INITIAL>endval {BEGIN DYNARE_BLOCK; return token::ENDVAL;}
<INITIAL>histval {BEGIN DYNARE_BLOCK; return token::HISTVAL;}
<INITIAL>shocks {BEGIN DYNARE_BLOCK; return token::SHOCKS;}
<INITIAL>estimated_params {BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS;}
<INITIAL>estimated_params_init 		{BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS_INIT;}
<INITIAL>estimated_params_bounds 	{BEGIN DYNARE_BLOCK; return token::ESTIMATED_PARAMS_BOUNDS;}
<INITIAL>observation_trends {BEGIN DYNARE_BLOCK; return token::OBSERVATION_TRENDS;}
<INITIAL>optim_weights {BEGIN DYNARE_BLOCK; return token::OPTIM_WEIGHTS;}

 /* End of a Dynare block */
<DYNARE_BLOCK>end[ \t\n]*; 	{BEGIN INITIAL; return token::END;}   

 /* Inside  of a Dynare statement */
<DYNARE_STATEMENT>datafile 		{return token::DATAFILE;}
<DYNARE_STATEMENT>nobs 			{return token::NOBS;}
<DYNARE_STATEMENT>first_obs 		{return token::FIRST_OBS;}
<DYNARE_STATEMENT>prefilter 		{return token::PREFILTER;} 
<DYNARE_STATEMENT>presample 		{return token::PRESAMPLE;} 
<DYNARE_STATEMENT>lik_algo  		{return token::LIK_ALGO;}  
<DYNARE_STATEMENT>lik_init  		{return token::LIK_INIT;}  
<DYNARE_STATEMENT>graph   		{return token::GRAPH;}  	  
<DYNARE_STATEMENT>nograph   		{return token::NOGRAPH;}  	  
<DYNARE_STATEMENT>print   		{return token::PRINT;}  	  
<DYNARE_STATEMENT>noprint   		{return token::NOPRINT;}  	  
<DYNARE_STATEMENT>conf_sig  		{return token::CONF_SIG;}  
<DYNARE_STATEMENT>mh_replic 		{return token::MH_REPLIC;} 
<DYNARE_STATEMENT>mh_drop   		{return token::MH_DROP;}   
<DYNARE_STATEMENT>mh_jscale   		{return token::MH_JSCALE;}   
<DYNARE_STATEMENT>mh_init_scale 	{return token::MH_INIT_SCALE;}
<DYNARE_STATEMENT>mode_file 		{return token::MODE_FILE;}
<DYNARE_STATEMENT>mode_compute 	{return token::MODE_COMPUTE;}
<DYNARE_STATEMENT>mode_check 		{return token::MODE_CHECK;}
<DYNARE_STATEMENT>prior_trunc 	{return token::PRIOR_TRUNC;}
<DYNARE_STATEMENT>mh_mode 		{return token::MH_MODE;}
<DYNARE_STATEMENT>mh_nblocks 		{return token::MH_NBLOCKS;}
<DYNARE_STATEMENT>load_mh_file 	{return token::LOAD_MH_FILE;}
<DYNARE_STATEMENT>loglinear 		{return token::LOGLINEAR;}
<DYNARE_STATEMENT>nodiagnostic 	{return token::NODIAGNOSTIC;}
<DYNARE_STATEMENT>kalman_algo 	{return token::KALMAN_ALGO;}
<DYNARE_STATEMENT>kalman_tol 	{return token::KALMAN_TOL;}
<DYNARE_STATEMENT>forecast 	{return token::FORECAST;}
<DYNARE_STATEMENT>smoother 	{return token::SMOOTHER;}
<DYNARE_STATEMENT>bayesian_irf 	{return token::BAYESIAN_IRF;}
<DYNARE_STATEMENT>moments_varendo {return token::MOMENTS_VARENDO;}
<DYNARE_STATEMENT>filtered_vars	{return token::FILTERED_VARS;}
<DYNARE_STATEMENT>filter_step_ahead	{return token::FILTER_STEP_AHEAD;}
<DYNARE_STATEMENT>relative_irf 	{return token::RELATIVE_IRF;}
<DYNARE_STATEMENT>tex		{return token::TEX;}
<DYNARE_STATEMENT>moments	{return token::MOMENTS;}
<DYNARE_STATEMENT>nomoments	{return token::NOMOMENTS;}
<DYNARE_STATEMENT>corr		{return token::CORR;}
<DYNARE_STATEMENT>nocorr	{return token::NOCORR;}
<DYNARE_STATEMENT>optim		{return token::OPTIM;}
<DYNARE_STATEMENT>periods	{return token::PERIODS;}
<DYNARE_STATEMENT>diffuse_d {return token::DIFFUSE_D;}
<DYNARE_STATEMENT>model_comparison_approximation {return token::MODEL_COMPARISON;}
<DYNARE_STATEMENT>laplace {return token::LAPLACE;}
<DYNARE_STATEMENT>modifiedharmonicmean {return token::MODIFIEDHARMONICMEAN;}
<DYNARE_STATEMENT>constant	{return token::CONSTANT;}
<DYNARE_STATEMENT>noconstant	{return token::NOCONSTANT;}

<DYNARE_STATEMENT>[\$][^$]*[\$] {
  strtok(yytext+1, "$");
  yylval->string_val = new string(yytext + 1); 
  return token::TEX_NAME;
}

 /* Inside a Dynare block */
<DYNARE_BLOCK>var {return token::VAR;}
<DYNARE_BLOCK>stderr {return token::STDERR;}
<DYNARE_BLOCK>values {return token::VALUES;}
<DYNARE_BLOCK>corr {return token::CORR;}
<DYNARE_BLOCK>periods {return token::PERIODS;}

<DYNARE_BLOCK>gamma_pdf {return token::GAMMA_PDF;}
<DYNARE_BLOCK>beta_pdf {return token::BETA_PDF;}
<DYNARE_BLOCK>normal_pdf {return token::NORMAL_PDF;}
<DYNARE_BLOCK>inv_gamma_pdf {return token::INV_GAMMA_PDF;}
<DYNARE_BLOCK>inv_gamma1_pdf {return token::INV_GAMMA_PDF;}
<DYNARE_BLOCK>inv_gamma2_pdf {return token::INV_GAMMA_PDF;}
<DYNARE_BLOCK>uniform_pdf {return token::UNIFORM_PDF;}

<DYNARE_BLOCK>; {return yy::parser::token_type (yytext[0]);}


 /* Inside Dynare statement */
<DYNARE_STATEMENT>solve_algo {return token::SOLVE_ALGO;}
<DYNARE_STATEMENT>dr_algo {return token::DR_ALGO;}
<DYNARE_STATEMENT>simul_algo {return token::SIMUL_ALGO;}
<DYNARE_STATEMENT>drop {return token::DROP;}
<DYNARE_STATEMENT>order {return token::ORDER;}
<DYNARE_STATEMENT>replic {return token::REPLIC;}
<DYNARE_STATEMENT>ar {return token::AR;}
<DYNARE_STATEMENT>nofunctions {return token::NOFUNCTIONS;}
<DYNARE_STATEMENT>irf {return token::IRF;}
<DYNARE_STATEMENT>hp_filter {return token::HP_FILTER;}
<DYNARE_STATEMENT>hp_ngrid {return token::HP_NGRID;}
<DYNARE_STATEMENT>simul_seed {return token::SIMUL_SEED;}
<DYNARE_STATEMENT>qz_criterium {return token::QZ_CRITERIUM;}
<DYNARE_STATEMENT>simul {return token::SIMUL;}
<DYNARE_STATEMENT>autocorr {return token::AUTOCORR;}
<DYNARE_STATEMENT>olr_beta {return token::OLR_BETA;}
<DYNARE_STATEMENT>xtick   		{return token::XTICK;}  	  
<DYNARE_STATEMENT>xticklabel   		{return token::XTICKLABEL;}  	  
<DYNARE_STATEMENT>xls_sheet {return token::XLS_SHEET;}
<DYNARE_STATEMENT>xls_range {return token::XLS_RANGE;}

<DYNARE_STATEMENT,DYNARE_BLOCK>use_dll {return token::USE_DLL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>linear {return token::LINEAR;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[,] {return token::COMMA;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\(\)] {return yy::parser::token_type (yytext[0]);} 
<DYNARE_STATEMENT,DYNARE_BLOCK>[\[] {return yy::parser::token_type (yytext[0]);}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\]] {
  if (sigma_e)
    sigma_e = 0;
  return yy::parser::token_type (yytext[0]);
}
<DYNARE_STATEMENT,DYNARE_BLOCK>[+] {return token::PLUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[-] {return token::MINUS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[*] {return token::TIMES;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[/] {return token::DIVIDE;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[=] {return token::EQUAL;}
<DYNARE_STATEMENT,DYNARE_BLOCK>[\^] {return token::POWER;}
<DYNARE_STATEMENT,DYNARE_BLOCK>exp {return token::EXP;}
<DYNARE_STATEMENT,DYNARE_BLOCK>log {return token::LOG;}
<DYNARE_STATEMENT,DYNARE_BLOCK>log10 {return token::LOG10;}
<DYNARE_STATEMENT,DYNARE_BLOCK>ln {return token::LOG;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sin {return token::SIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>cos {return token::COS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>tan {return token::TAN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>asin {return token::ASIN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>acos {return token::ACOS;}
<DYNARE_STATEMENT,DYNARE_BLOCK>atan {return token::ATAN;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sinh {return token::SINH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>cosh {return token::COSH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>tanh {return token::TANH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>asinh {return token::ASINH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>acosh {return token::ACOSH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>atanh {return token::ATANH;}
<DYNARE_STATEMENT,DYNARE_BLOCK>sqrt  {return token::SQRT;}

<DYNARE_STATEMENT,DYNARE_BLOCK>[A-Za-z_][A-Za-z0-9_]* {
  yylval->string_val = new string(yytext);
  return token::NAME;
}

<DYNARE_STATEMENT,DYNARE_BLOCK>((([0-9]*\.[0-9]+)|([0-9]+\.))([edED][-+]?[0-9]+)?)|([0-9]+[edED][-+]?[0-9]+) {
  yylval->string_val = new string(yytext);
  return token::FLOAT_NUMBER;
}

<DYNARE_STATEMENT,DYNARE_BLOCK>[0-9]+ {
  yylval->string_val = new string(yytext);
  return token::INT_NUMBER;
}
   	
 /* an instruction starting with a recognized symbol is passed as NAME,
    otherwise it is a native statement until the end of the line
 */
<INITIAL>[A-Za-z_][A-Za-z0-9_]* {	      
  if (SymbolTable::getID(yytext) != -1)
    {
      BEGIN DYNARE_STATEMENT;
      yylval->string_val = new string(yytext);
      return token::NAME;
    }
  else
    {
      BEGIN NATIVE;
      *(driver.output) << yytext;
    }
}

<INITIAL>. {BEGIN NATIVE; *(driver.output) << yytext;}

 /* NATIVE Block */
<NATIVE>.* {BEGIN INITIAL; *(driver.output) << yytext << endl;}

<*>. {return yy::parser::token_type (yytext[0]);}

%%

void
ParsingDriver::scan_begin()
{
  yy_flex_debug = trace_scanning;
  if (!(yyin = fopen(file.c_str(), "r")))
    error(string("cannot open file"));
}

void
ParsingDriver::scan_end()
{
  fclose(yyin);
}

/*
  Local variables:
  mode: C++
  End:
*/
