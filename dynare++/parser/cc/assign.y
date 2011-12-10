%{
/* Copyright (C) 2006-2011, Ondra Kamenik */

#include "location.h"
#include "atom_assignings.h"
#include "assign_tab.hh"

#include <stdio.h>

	void asgn_error(const char*);
	int asgn_lex(void);
	extern int asgn_lineno;
	extern ogp::AtomAssignings* aparser;

%}

%union {
	int integer;
	char *string;
	char character;
}

%token EQUAL_SIGN SEMICOLON CHARACTER BLANK
%token <string> NAME;

%name-prefix="asgn_"

%locations
%error-verbose

%%

root : assignments | ;

assignments : assignments BLANK | assignments assignment | assignment | BLANK;

assignment : NAME EQUAL_SIGN material SEMICOLON {
	aparser->add_assignment(@1.off, $1, @1.ll, @3.off-@1.off, @3.ll + @4.ll);}
  | NAME space EQUAL_SIGN material SEMICOLON {
	aparser->add_assignment(@1.off, $1, @1.ll, @4.off-@1.off, @4.ll + @5.ll);}
  ;

material : material CHARACTER | material NAME | material BLANK | NAME | CHARACTER | BLANK;

space : space BLANK | BLANK;

%%

void asgn_error(const char* mes)
{
	aparser->error(mes);
}
