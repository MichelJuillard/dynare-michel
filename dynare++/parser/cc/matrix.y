%{
#include "location.h"
#include "matrix_parser.h" 
#include "matrix_tab.hh"

	void matrix_error(char*);
	int matrix_lex(void);
	extern int matrix_lineno;
	extern ogp::MatrixParser* mparser;
	extern YYLTYPE matrix_lloc;

//	static void print_token_value (FILE *, int, YYSTYPE);
//#define YYPRINT(file, type, value) print_token_value (file, type, value)

%}

%union {
	double val;
	int integer;
}

%token NEW_ROW
%token <val> DNUMBER

%name-prefix="matrix_";

%locations
%error-verbose

%%

matrix : first_row other_rows
    | first_row other_rows empty_rows
    | first_row empty_rows other_rows empty_rows
    | first_row empty_rows other_rows
    | empty_rows first_row other_rows
    | empty_rows first_row other_rows empty_rows
    | empty_rows first_row empty_rows other_rows empty_rows
    | empty_rows first_row empty_rows
    | first_row empty_rows
    | empty_rows first_row
    | first_row
    | empty_rows
    ;

empty_rows : empty_rows NEW_ROW | NEW_ROW;

lod : DNUMBER {mparser->add_item($1);}
    | lod DNUMBER {mparser->add_item($2);}
    ;

first_row : lod;

other_rows : other_rows one_row | other_rows empty_rows one_row |one_row ;

one_row : NEW_ROW {mparser->start_row();} lod;


%%

void matrix_error(char* s)
{
	mparser->error(s);
}


