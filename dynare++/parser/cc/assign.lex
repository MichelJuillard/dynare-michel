%{
#include "location.h"
#include "assign_tab.hh"

	extern YYLTYPE asgn_lloc;

#define YY_USER_ACTION SET_LLOC(asgn_);
%}

%option nounput
%option noyy_top_state
%option stack
%option yylineno
%option prefix="asgn_"
%option never-interactive
%x CMT

%%

 /* comments */
<*>"/*"            {yy_push_state(CMT);}
<CMT>[^*\n]*
<CMT>"*"+[^*/\n]*
<CMT>"*"+"/"       {yy_pop_state();}
<CMT>[\n]
"//".*\n

 /* spaces */
[ \t\r\n]          {return BLANK;}

 /* names */
[A-Za-z_][A-Za-z0-9_]* {
	asgn_lval.string = asgn_text;
	return NAME;
}

;                  {return SEMICOLON;}
=                  {return EQUAL_SIGN;}
. {
	asgn_lval.character = asgn_text[0];
	return CHARACTER;
}

%%

int asgn_wrap()
{
	return 1;
}

void asgn__destroy_buffer(void* p)
{
	asgn__delete_buffer((YY_BUFFER_STATE)p);
}
