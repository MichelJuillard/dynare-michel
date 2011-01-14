%{
// Copyright (C) 2006-2011, Ondra Kamenik

#include "parser/cc/location.h"
#include "dynare_model.h"
#include "dynglob_tab.hh"

#include <stdio.h>

	void dynglob_error(const char*);
	int dynglob_lex(void);
	extern int dynglob_lineno;
	extern ogdyn::DynareParser* dynare_parser;
	int symblist_flag;

	static void print_token_value1 (FILE *, int, YYSTYPE);
#define YYPRINT(file, type, value) print_token_value1 (file, type, value)

%}

%union {
	int integer;
	char *string;
	char character;
}

%token  END INITVAL MODEL PARAMETERS VAR VAREXO SEMICOLON COMMA EQUAL_SIGN CHARACTER
%token  VCOV LEFT_BRACKET RIGHT_BRACKET ORDER PLANNEROBJECTIVE PLANNERDISCOUNT
%token <string> NAME;

%name-prefix="dynglob_"

%locations
%error-verbose

%%

dynare_file : preamble paramset model rest {
	dynare_parser->set_paramset_pos(@2.off, @3.off);}
  | preamble model rest {
	dynare_parser->set_paramset_pos(0, 0);}
  | preamble paramset planner model rest {
	dynare_parser->set_paramset_pos(@2.off, @3.off);}
  ;

preamble : preamble preamble_statement | preamble_statement;

preamble_statement : var | varexo | parameters;

var : VAR {symblist_flag=1;} symblist SEMICOLON;

varexo : VAREXO {symblist_flag=2;} symblist SEMICOLON;

parameters : PARAMETERS {symblist_flag=3;} symblist SEMICOLON;


symblist : symblist NAME          {dynare_parser->add_name($2,symblist_flag);}
     | symblist COMMA NAME        {dynare_parser->add_name($3,symblist_flag);}
     | NAME                       {dynare_parser->add_name($1,symblist_flag);}
     ;

paramset : recnameset;

recnameset : recnameset onenameset | onenameset;

onenameset : NAME EQUAL_SIGN material SEMICOLON;

material : material CHARACTER | material NAME | NAME | CHARACTER;

model : MODEL SEMICOLON equations END SEMICOLON {
	dynare_parser->set_model_pos(@3.off, @4.off);
};

equations : equations equation | equation;

equation : material EQUAL_SIGN material SEMICOLON | material SEMICOLON;

rest : rest_statement | rest rest_statement;

rest_statement : initval | vcov | order | planner;

initval : INITVAL SEMICOLON recnameset END SEMICOLON {
	dynare_parser->set_initval_pos(@3.off, @4.off);
};

vcov : VCOV EQUAL_SIGN LEFT_BRACKET m_material RIGHT_BRACKET SEMICOLON {
	dynare_parser->set_vcov_pos(@4.off, @5.off);
};

m_material : m_material CHARACTER | m_material NAME | m_material SEMICOLON | m_material COMMA | CHARACTER | NAME | SEMICOLON | COMMA; 

order : ORDER EQUAL_SIGN material SEMICOLON {
    dynare_parser->set_order_pos(@3.off, @4.off);
};

planner : planner_objective planner_discount
  | planner_discount planner_objective
;

planner_objective : PLANNEROBJECTIVE material SEMICOLON {
	dynare_parser->set_pl_objective_pos(@2.off, @3.off);
};

planner_discount : PLANNERDISCOUNT NAME SEMICOLON {
	dynare_parser->set_pl_discount_pos(@2.off, @3.off);
};

%%

void dynglob_error(const char* mes)
{
	dynare_parser->error(mes);
}

static void print_token_value1(FILE* file, int type, YYSTYPE value)
{
	if (type == NAME)
		fprintf(file, "%s", value.string);
	if (type == CHARACTER)
		fprintf(file, "%c", value.character);
}
