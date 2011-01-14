%{
/* Copyright (C) 2006-2011, Ondra Kamenik */

#include <cstdio>

#include "location.h"
#include "formula_parser.h" 
#include "formula_tab.hh"

	void fmla_error(const char*);
	int fmla_lex(void);
	extern int fmla_lineno;
	extern ogp::FormulaParser* fparser;
	extern YYLTYPE fmla_lloc;

  //  static void print_token_value (FILE *, int, YYSTYPE);
  //  #define YYPRINT(file, type, value) print_token_value (file, type, value)

%}

%union {
	char* string;
	double dvalue;
	int integer;
}

%token EQUAL_SIGN
%left YPLUS YMINUS
%left YTIMES YDIVIDE
%left YUMINUS YUPLUS
%right YPOWER
%token YEXP YLOG YSIN YCOS YTAN YSQRT YERF YERFC YDIFF
%token <string> DNUMBER NAME
%type <integer> expression

%name-prefix="fmla_"

%locations
%error-verbose

%%
 root : equation_list
      | expression
                                  {fparser->add_formula($1);}
      ; 

 equation_list : equation_list equation | equation ;

 equation : expression EQUAL_SIGN expression ';' 
                                  {fparser->add_formula(fparser->add_binary(ogp::MINUS,$1,$3));}
      | expression ';'
                                  {fparser->add_formula($1);}
      ;

  expression : '(' expression ')' { $$ = $2;}
      | expression YPLUS expression {$$=fparser->add_binary(ogp::PLUS,$1,$3);}
      | expression YMINUS expression {$$=fparser->add_binary(ogp::MINUS,$1,$3);}
      | expression YTIMES expression {$$=fparser->add_binary(ogp::TIMES,$1,$3);}
      | expression YDIVIDE expression {$$=fparser->add_binary(ogp::DIVIDE,$1,$3);}
      | expression YPOWER expression {$$=fparser->add_binary(ogp::POWER,$1,$3);}
      | YMINUS expression %prec YUMINUS {$$=fparser->add_unary(ogp::UMINUS,$2);}
      | YPLUS expression %prec YUPLUS {$$ = $2;}
      | YSIN '(' expression ')' {$$=fparser->add_unary(ogp::SIN,$3);}
      | YCOS '(' expression ')' {$$=fparser->add_unary(ogp::COS,$3);}
      | YTAN '(' expression ')' {$$=fparser->add_unary(ogp::TAN,$3);}
      | YEXP '(' expression ')' {$$=fparser->add_unary(ogp::EXP,$3);}
      | YLOG '(' expression ')' {$$=fparser->add_unary(ogp::LOG,$3);}
      | YSQRT '(' expression ')' {$$=fparser->add_unary(ogp::SQRT,$3);}
      | YERF '(' expression ')' {$$=fparser->add_unary(ogp::ERF,$3);}
      | YERFC '(' expression ')' {$$=fparser->add_unary(ogp::ERFC,$3);}
      | YDIFF '(' expression ',' NAME ')' {$$=fparser->add_derivative($3, fparser->add_nulary($5));}
      | NAME {$$=fparser->add_nulary($1);}
      | DNUMBER {$$=fparser->add_nulary($1);}
      ;

%%

void fmla_error(const char* s)
{
	fparser->error(s);
}

/*
static void print_token_value(FILE* file, int type, YYSTYPE value)
{
	if (type == NAME)
		fprintf(file, "%s", value.string);
}
*/
