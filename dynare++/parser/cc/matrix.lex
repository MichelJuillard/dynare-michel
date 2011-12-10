%{
// Copyright (C) 2006-2011, Ondra Kamenik

#include "location.h"
#include "matrix_tab.hh"

	extern YYLTYPE matrix_lloc;
	extern void matrix_error(const char*);

#define YY_USER_ACTION SET_LLOC(matrix_);
%}

%option nounput
%option noyy_top_state
%option stack
%option yylineno
%option prefix="matrix_"
%option never-interactive
%x CMT

%%

 /* comments */
<*>"/*"              {yy_push_state(CMT);}
<CMT>[^*\n]*
<CMT>"*"+[^*/\n]*
<CMT>"*"+"/"         {yy_pop_state();}
<CMT>[\n]
"//".*\n

 /* ignore spaces and commas */
[ \t,]
 /* new row */
\r\n                 {return NEW_ROW;}
\n                   {return NEW_ROW;}
;[ \t]*\n            {return NEW_ROW;}
;[ \t]*\r\n          {return NEW_ROW;}
;                    {return NEW_ROW;}

[+-]?(([0-9]*\.?[0-9]+)|([0-9]+\.))([edED][-+]?[0-9]+)? {
	matrix_lval.val = strtod(matrix_text, NULL);
	return DNUMBER;
}

. {
	char mes[300];
	sprintf(mes, "Unrecognized character %s", matrix_text);
	matrix_error(mes); 
}

%%

int matrix_wrap()
{
	return 1;
}

void matrix__destroy_buffer(void* p)
{
	matrix__delete_buffer((YY_BUFFER_STATE)p);
}
