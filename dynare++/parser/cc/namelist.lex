%{
#include "location.h"
#include "namelist_tab.hh"

	extern YYLTYPE namelist_lloc;

#define YY_USER_ACTION SET_LLOC(namelist_);
%}

%option nounput
%option noyy_top_state
%option stack
%option prefix="namelist_"
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

 /* initial spaces or tabs are ignored */
[ \t\r\n\0]

 /* names */
[A-Za-z_][A-Za-z0-9_]* {
	namelist_lval.string = namelist_text;
	return NAME;
}

,                  {return COMMA;}
. {
	namelist_lval.character = namelist_text[0];
	return CHARACTER;
}

%%

int namelist_wrap()
{
	return 1;
}

void namelist__destroy_buffer(void* p)
{
	namelist__delete_buffer((YY_BUFFER_STATE)p);
}
