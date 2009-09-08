%{
#include "location.h"
#include "formula_tab.hh"

	extern YYLTYPE fmla_lloc;

#define YY_USER_ACTION SET_LLOC(fmla_);
%}

%option nounput
%option noyy_top_state
%option stack
%option yylineno
%option prefix="fmla_"
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

 /* initial spaces or tabs are ignored */

[ \t\r\n]
[+]                  {return YPLUS;}
[-]                  {return YMINUS;}
[*]                  {return YTIMES;}
[/]                  {return YDIVIDE;}
[\^]                 {return YPOWER;}
exp                  {return YEXP;}
log                  {return YLOG;}
sin                  {return YSIN;}
cos                  {return YCOS;}
tan                  {return YTAN;}
sqrt                 {return YSQRT;}
erf                  {return YERF;}
erfc                 {return YERFC;}
diff                 {return YDIFF;}

 /* names: parameters, variables (lagged/leaded) */
[A-Za-z_][A-Za-z0-9_]*([\(\{][+-]?[0-9]+[\)\}])? {
	fmla_lval.string=fmla_text;
	return NAME;
}

 /* floating point numbers */
(([0-9]*\.?[0-9]+)|([0-9]+\.))([edED][-+]?[0-9]+)? {
	fmla_lval.string=fmla_text;
	return DNUMBER;
}

=                    {return EQUAL_SIGN;}

.                    {return fmla_text[0];}

%%

int fmla_wrap()
{
	return 1;
}

void fmla__destroy_buffer(void* p)
{
	fmla__delete_buffer((YY_BUFFER_STATE)p);
}
