// Copyright (C) 2007-2011, Ondra Kamenik

%{
#include "location.h"
#include "namelist.h"
#include "namelist_tab.hh"

	void namelist_error(const char*);
	int namelist_lex(void);
	extern ogp::NameListParser* name_list_parser;

%}

%union {
	int integer;
	char *string;
	char character;
}

%token COMMA CHARACTER
%token <string> NAME;

%name-prefix="namelist_"

%locations
%error-verbose

%%

namelist : namelist NAME       {name_list_parser->add_name($2);}
         | namelist COMMA NAME {name_list_parser->add_name($3);}
         | NAME                {name_list_parser->add_name($1);}
         ;

%%

void namelist_error(const char* mes)
{
	name_list_parser->namelist_error(mes);
}
