%{
// Copyright (C) 2007-2011, Ondra Kamenik

#include "location.h"
#include "csv_tab.hh"

	extern YYLTYPE csv_lloc;

#define YY_USER_ACTION SET_LLOC(csv_);
%}

%option nounput
%option noyy_top_state
%option yylineno
%option prefix="csv_"
%option never-interactive

%%

,                     {return COMMA;}
\n                    {return NEWLINE;}
\r\n                  {return NEWLINE;}
[^,\n\r]+             {return ITEM;}

%%

int csv_wrap()
{
	return 1;
}

void csv__destroy_buffer(void* p)
{
	csv__delete_buffer((YY_BUFFER_STATE)p);
}
