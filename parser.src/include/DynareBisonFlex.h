#ifndef DYNARE_BISON_FLEX_H
#define DYNARE_BISON_FLEX_H

#include <iostream>
#include <string>
#include <FlexLexer.h>

using namespace std;

#include "SymbolTable.h"

// Parse function of 'bison' is defined externally
extern "C" int yyparse(void *);

// The error function that 'bison' calls
inline void yyerror(char const *what_error) { cerr << what_error << endl; }


// Control class for a flex-scanner
// It has a pointer to the compiler in which it is contained
class FlexerClass//sem_flexer
  : public yyFlexLexer
{
public:
  FlexerClass(void *_compiler) : Compiler(_compiler) {}

  void * getCompiler() const { return Compiler; }

private:
  void *Compiler;
};


// Base class for compilers which use 'flex' and 'bison' as scanner/parser
class FlexBisonClass//flex_bison_compiler
{
public:
  friend int yylex(string *new_token, void *compiler);

  FlexBisonClass() : flexer(this) {}

  // Compiles program from a stream
  void Compile(istream &program_strm = cin)
  {
    // Tells 'flex' the right stream
    flexer.switch_streams(&program_strm, 0);

    // If there is an error parsing the source
    if (! parse())
      {
        cerr << "How bad! A parse error has occurred!" << endl;
        exit(1);
      }
  }

private:
  // The scanner used in this compiler (it is a flex-scanner)
  FlexerClass flexer;

  int scan() { return flexer.yylex(); }
  bool parse() { return ! yyparse((void *) this); }
};


// Directs the call from 'bison' to the scanner in the right compiler
// "new_token" is not needed in this case...
inline int yylex(string *new_token, void *compiler)
{
  FlexBisonClass &what_compiler
    = *static_cast<FlexBisonClass *>(compiler); 
  cout << *new_token; 
  return what_compiler.scan();
}


// Definitions for 'flex' and 'bison'
// Read the manuals of the two beasties to see what these macros mean!

#define yywrap() 1
#define YY_SKIP_YYWRAP
#define YYPARSE_PARAM parm
#define YYLEX_PARAM parm
#define YYSTYPE string

/*
// A test compiler
class some_compiler
  : public FlexBisonClass
{
public:
  void just_a_flex_test(string const &some_string) const
  {
    cout << some_string << endl; 
  }

  void just_a_bison_test() const
  {
    cout << "Cool! Even bison has found something!" << endl;
  }
};


// Important! These are the "shortcuts" which you can use in your
// ".l"- and ".y"-files to access the corresponding compiler-object!
#define my_flex_compiler (*static_cast<some_compiler *> (static_cast<FlexerClass *>(this)->getCompiler()))
#define my_bison_compiler (*static_cast<some_compiler *> (parm))

*/

#endif
