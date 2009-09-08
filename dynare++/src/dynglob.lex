%{
#include "parser/cc/location.h"
#include "dynglob_tab.hh"

	extern YYLTYPE dynglob_lloc;

#define YY_USER_ACTION SET_LLOC(dynglob_);
%}

%option nounput
%option noyy_top_state
%option stack
%option yylineno
%option prefix="dynglob_"
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
var                {return VAR;}
varexo             {return VAREXO;}
parameters         {return PARAMETERS;}
model              {return MODEL;}
end                {return END;}
initval            {return INITVAL;}
order              {return ORDER;}
vcov               {return VCOV;}
planner_objective  {return PLANNEROBJECTIVE;}
planner_discount   {return PLANNERDISCOUNT;}

 /* names */
[A-Za-z_][A-Za-z0-9_]* {
	dynglob_lval.string = dynglob_text;
	return NAME;
}

;                  {return SEMICOLON;}
,                  {return COMMA;}
=                  {return EQUAL_SIGN;}
\[                 {return LEFT_BRACKET;}
\]                 {return RIGHT_BRACKET;}
. {
	dynglob_lval.character = dynglob_text[0];
	return CHARACTER;
}

%%

int dynglob_wrap()
{
	return 1;
}

void dynglob__destroy_buffer(void* p)
{
	dynglob__delete_buffer((YY_BUFFER_STATE)p);
}
