/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison LALR(1) parsers in C++

   Copyright (C) 2002, 2003, 2004, 2005, 2006 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C++ LALR(1) parser skeleton written by Akim Demaille.  */

#ifndef PARSER_HEADER_H
# define PARSER_HEADER_H

#include <string>
#include <iostream>
#include "stack.hh"

namespace yy
{
  class position;
  class location;
}

/* First part of user declarations.  */
#line 5 "DynareBison.yy"

using namespace std;

class ParsingDriver;

#include "SymbolTableTypes.hh"
#include "ExprNode.hh"

//! Type for semantic value of non-derivable expressions
typedef pair<int, Type> ExpObj;


/* Line 303 of lalr1.cc.  */
#line 66 "DynareBison.hh"

#include "location.hh"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)		\
do {							\
  if (N)						\
    {							\
      (Current).begin = (Rhs)[1].begin;			\
      (Current).end   = (Rhs)[N].end;			\
    }							\
  else							\
    {							\
      (Current).begin = (Current).end = (Rhs)[0].end;	\
    }							\
} while (false)
#endif

namespace yy
{

  /// A Bison parser.
  class parser
  {
  public:
    /// Symbol semantic values.
#ifndef YYSTYPE
    union semantic_type
#line 31 "DynareBison.yy"
{
  string *string_val;
  ExpObj *exp_val;
  NodeID model_val;
}
/* Line 303 of lalr1.cc.  */
#line 124 "DynareBison.hh"
	;
#else
    typedef YYSTYPE semantic_type;
#endif
    /// Symbol locations.
    typedef location location_type;
    /// Tokens.
    struct token
    {
      /* Tokens.  */
   enum yytokentype {
     AR = 258,
     AUTOCORR = 259,
     BAYESIAN_IRF = 260,
     BETA_PDF = 261,
     CALIB = 262,
     CALIB_VAR = 263,
     CHECK = 264,
     CONF_SIG = 265,
     CONSTANT = 266,
     CORR = 267,
     COVAR = 268,
     DATAFILE = 269,
     DR_ALGO = 270,
     DROP = 271,
     DSAMPLE = 272,
     DYNASAVE = 273,
     DYNATYPE = 274,
     END = 275,
     ENDVAL = 276,
     EQUAL = 277,
     ESTIMATION = 278,
     ESTIMATED_PARAMS = 279,
     ESTIMATED_PARAMS_BOUNDS = 280,
     ESTIMATED_PARAMS_INIT = 281,
     FILTER_STEP_AHEAD = 282,
     FILTERED_VARS = 283,
     FIRST_OBS = 284,
     FLOAT_NUMBER = 285,
     FORECAST = 286,
     GAMMA_PDF = 287,
     GRAPH = 288,
     HISTVAL = 289,
     HP_FILTER = 290,
     HP_NGRID = 291,
     INITVAL = 292,
     INT_NUMBER = 293,
     INV_GAMMA_PDF = 294,
     IRF = 295,
     KALMAN_ALGO = 296,
     KALMAN_TOL = 297,
     LAPLACE = 298,
     LIK_ALGO = 299,
     LIK_INIT = 300,
     LINEAR = 301,
     LOAD_MH_FILE = 302,
     LOGLINEAR = 303,
     MH_DROP = 304,
     MH_INIT_SCALE = 305,
     MH_JSCALE = 306,
     MH_MODE = 307,
     MH_NBLOCKS = 308,
     MH_REPLIC = 309,
     MH_RECOVER = 310,
     MODE_CHECK = 311,
     MODE_COMPUTE = 312,
     MODE_FILE = 313,
     MODEL = 314,
     MODEL_COMPARISON = 315,
     MSHOCKS = 316,
     MODEL_COMPARISON_APPROXIMATION = 317,
     MODIFIEDHARMONICMEAN = 318,
     MOMENTS_VARENDO = 319,
     NAME = 320,
     NOBS = 321,
     NOCONSTANT = 322,
     NOCORR = 323,
     NODIAGNOSTIC = 324,
     NOFUNCTIONS = 325,
     NOGRAPH = 326,
     NOMOMENTS = 327,
     NOPRINT = 328,
     NORMAL_PDF = 329,
     OBSERVATION_TRENDS = 330,
     OLR = 331,
     OLR_INST = 332,
     OLR_BETA = 333,
     OPTIM = 334,
     OPTIM_WEIGHTS = 335,
     ORDER = 336,
     OSR = 337,
     OSR_PARAMS = 338,
     PARAMETERS = 339,
     PERIODS = 340,
     PLANNER_OBJECTIVE = 341,
     PREFILTER = 342,
     PRESAMPLE = 343,
     PRINT = 344,
     PRIOR_TRUNC = 345,
     PRIOR_ANALYSIS = 346,
     POSTERIOR_ANALYSIS = 347,
     QZ_CRITERIUM = 348,
     RELATIVE_IRF = 349,
     REPLIC = 350,
     RPLOT = 351,
     SHOCKS = 352,
     SIGMA_E = 353,
     SIMUL = 354,
     SIMUL_ALGO = 355,
     SIMUL_SEED = 356,
     SMOOTHER = 357,
     SOLVE_ALGO = 358,
     STDERR = 359,
     STEADY = 360,
     STOCH_SIMUL = 361,
     TEX = 362,
     RAMSEY_POLICY = 363,
     PLANNER_DISCOUNT = 364,
     TEX_NAME = 365,
     UNIFORM_PDF = 366,
     UNIT_ROOT_VARS = 367,
     USE_DLL = 368,
     VALUES = 369,
     VAR = 370,
     VAREXO = 371,
     VAREXO_DET = 372,
     VAROBS = 373,
     XLS_SHEET = 374,
     XLS_RANGE = 375,
     COMMA = 376,
     MINUS = 377,
     PLUS = 378,
     DIVIDE = 379,
     TIMES = 380,
     UMINUS = 381,
     POWER = 382,
     EXP = 383,
     LOG = 384,
     LOG10 = 385,
     SIN = 386,
     COS = 387,
     TAN = 388,
     ASIN = 389,
     ACOS = 390,
     ATAN = 391,
     SINH = 392,
     COSH = 393,
     TANH = 394,
     ASINH = 395,
     ACOSH = 396,
     ATANH = 397,
     SQRT = 398
   };

    };
    /// Token type.
    typedef token::yytokentype token_type;

    /// Build a parser object.
    parser (ParsingDriver &driver_yyarg);
    virtual ~parser ();

    /// Parse.
    /// \returns  0 iff parsing succeeded.
    virtual int parse ();

    /// The current debugging stream.
    std::ostream& debug_stream () const;
    /// Set the current debugging stream.
    void set_debug_stream (std::ostream &);

    /// Type for debugging levels.
    typedef int debug_level_type;
    /// The current debugging level.
    debug_level_type debug_level () const;
    /// Set the current debugging level.
    void set_debug_level (debug_level_type l);

  private:
    /// Report a syntax error.
    /// \param loc    where the syntax error is found.
    /// \param msg    a description of the syntax error.
    virtual void error (const location_type& loc, const std::string& msg);

    /// Generate an error message.
    /// \param state   the state where the error occurred.
    /// \param tok     the look-ahead token.
    virtual std::string yysyntax_error_ (int yystate, int tok);

#if YYDEBUG
    /// \brief Report a symbol value on the debug stream.
    /// \param yytype       The token type.
    /// \param yyvaluep     Its semantic value.
    /// \param yylocationp  Its location.
    virtual void yy_symbol_value_print_ (int yytype,
					 const semantic_type* yyvaluep,
					 const location_type* yylocationp);
    /// \brief Report a symbol on the debug stream.
    /// \param yytype       The token type.
    /// \param yyvaluep     Its semantic value.
    /// \param yylocationp  Its location.
    virtual void yy_symbol_print_ (int yytype,
				   const semantic_type* yyvaluep,
				   const location_type* yylocationp);
#endif /* ! YYDEBUG */


    /// State numbers.
    typedef int state_type;
    /// State stack type.
    typedef stack<state_type>    state_stack_type;
    /// Semantic value stack type.
    typedef stack<semantic_type> semantic_stack_type;
    /// location stack type.
    typedef stack<location_type> location_stack_type;

    /// The state stack.
    state_stack_type yystate_stack_;
    /// The semantic value stack.
    semantic_stack_type yysemantic_stack_;
    /// The location stack.
    location_stack_type yylocation_stack_;

    /// Internal symbol numbers.
    typedef unsigned char token_number_type;
    /* Tables.  */
    /// For a state, the index in \a yytable_ of its portion.
    static const short int yypact_[];
    static const short int yypact_ninf_;

    /// For a state, default rule to reduce.
    /// Unless\a  yytable_ specifies something else to do.
    /// Zero means the default is an error.
    static const unsigned short int yydefact_[];

    static const short int yypgoto_[];
    static const short int yydefgoto_[];

    /// What to do in a state.
    /// \a yytable_[yypact_[s]]: what to do in state \a s.
    /// - if positive, shift that token.
    /// - if negative, reduce the rule which number is the opposite.
    /// - if zero, do what YYDEFACT says.
    static const short int yytable_[];
    static const signed char yytable_ninf_;

    static const unsigned short int yycheck_[];

    /// For a state, its accessing symbol.
    static const unsigned short int yystos_[];

    /// For a rule, its LHS.
    static const unsigned short int yyr1_[];
    /// For a rule, its RHS length.
    static const unsigned char yyr2_[];

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
    /// For a symbol, its name in clear.
    static const char* const yytname_[];
#endif

#if YYERROR_VERBOSE
    /// Convert the symbol name \a n to a form suitable for a diagnostic.
    virtual std::string yytnamerr_ (const char *n);
#endif

#if YYDEBUG
    /// A type to store symbol numbers and -1.
    typedef short int rhs_number_type;
    /// A `-1'-separated list of the rules' RHS.
    static const rhs_number_type yyrhs_[];
    /// For each rule, the index of the first RHS symbol in \a yyrhs_.
    static const unsigned short int yyprhs_[];
    /// For each rule, its source line number.
    static const unsigned short int yyrline_[];
    /// For each scanner token number, its symbol number.
    static const unsigned short int yytoken_number_[];
    /// Report on the debug stream that the rule \a r is going to be reduced.
    virtual void yy_reduce_print_ (int r);
    /// Print the state stack on the debug stream.
    virtual void yystack_print_ ();
#endif

    /// Convert a scanner token number \a t to a symbol number.
    token_number_type yytranslate_ (int t);

    /// \brief Reclaim the memory associated to a symbol.
    /// \param yymsg        Why this token is reclaimed.
    /// \param yytype       The symbol type.
    /// \param yyvaluep     Its semantic value.
    /// \param yylocationp  Its location.
    inline void yydestruct_ (const char* yymsg,
			     int yytype,
			     semantic_type* yyvaluep,
			     location_type* yylocationp);

    /// Pop \a n symbols the three stacks.
    inline void yypop_ (unsigned int n = 1);

    /* Constants.  */
    static const int yyeof_;
    /* LAST_ -- Last index in TABLE_.  */
    static const int yylast_;
    static const int yynnts_;
    static const int yyempty_;
    static const int yyfinal_;
    static const int yyterror_;
    static const int yyerrcode_;
    static const int yyntokens_;
    static const unsigned int yyuser_token_number_max_;
    static const token_number_type yyundef_token_;

    /* Debugging.  */
    int yydebug_;
    std::ostream* yycdebug_;


    /* User arguments.  */
    ParsingDriver &driver;
  };
}


#endif /* ! defined PARSER_HEADER_H */
