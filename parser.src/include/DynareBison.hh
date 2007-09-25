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

#include "ExprNode.hh"


/* Line 35 of lalr1.cc.  */
#line 62 "DynareBison.hh"

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
#line 27 "DynareBison.yy"
{
  string *string_val;
  NodeID node_val;
}
/* Line 35 of lalr1.cc.  */
#line 119 "DynareBison.hh"
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
     BVAR_DENSITY = 262,
     BVAR_FORECAST = 263,
     BVAR_PRIOR_DECAY = 264,
     BVAR_PRIOR_FLAT = 265,
     BVAR_PRIOR_LAMBDA = 266,
     BVAR_PRIOR_MU = 267,
     BVAR_PRIOR_OMEGA = 268,
     BVAR_PRIOR_TAU = 269,
     BVAR_PRIOR_TRAIN = 270,
     BVAR_REPLIC = 271,
     CALIB = 272,
     CALIB_VAR = 273,
     CHECK = 274,
     CONF_SIG = 275,
     CONSTANT = 276,
     CORR = 277,
     COVAR = 278,
     CUTOFF = 279,
     DATAFILE = 280,
     DR_ALGO = 281,
     DROP = 282,
     DSAMPLE = 283,
     DYNASAVE = 284,
     DYNATYPE = 285,
     END = 286,
     ENDVAL = 287,
     EQUAL = 288,
     ESTIMATION = 289,
     ESTIMATED_PARAMS = 290,
     ESTIMATED_PARAMS_BOUNDS = 291,
     ESTIMATED_PARAMS_INIT = 292,
     FILENAME = 293,
     FILTER_STEP_AHEAD = 294,
     FILTERED_VARS = 295,
     FIRST_OBS = 296,
     FLOAT_NUMBER = 297,
     FORECAST = 298,
     GAMMA_PDF = 299,
     GCC_COMPILER = 300,
     GRAPH = 301,
     HISTVAL = 302,
     HP_FILTER = 303,
     HP_NGRID = 304,
     INITVAL = 305,
     INT_NUMBER = 306,
     INV_GAMMA_PDF = 307,
     IRF = 308,
     KALMAN_ALGO = 309,
     KALMAN_TOL = 310,
     LAPLACE = 311,
     LCC_COMPILER = 312,
     LIK_ALGO = 313,
     LIK_INIT = 314,
     LINEAR = 315,
     LOAD_MH_FILE = 316,
     LOGLINEAR = 317,
     MARKOWITZ = 318,
     MH_DROP = 319,
     MH_INIT_SCALE = 320,
     MH_JSCALE = 321,
     MH_MODE = 322,
     MH_NBLOCKS = 323,
     MH_REPLIC = 324,
     MH_RECOVER = 325,
     MODE_CHECK = 326,
     MODE_COMPUTE = 327,
     MODE_FILE = 328,
     MODEL = 329,
     MODEL_COMPARISON = 330,
     MSHOCKS = 331,
     MODEL_COMPARISON_APPROXIMATION = 332,
     MODIFIEDHARMONICMEAN = 333,
     MOMENTS_VARENDO = 334,
     NAME = 335,
     NOBS = 336,
     NOCONSTANT = 337,
     NOCORR = 338,
     NODIAGNOSTIC = 339,
     NOFUNCTIONS = 340,
     NOGRAPH = 341,
     NOMOMENTS = 342,
     NOPRINT = 343,
     NORMAL_PDF = 344,
     OBSERVATION_TRENDS = 345,
     OPTIM = 346,
     OPTIM_WEIGHTS = 347,
     ORDER = 348,
     OSR = 349,
     OSR_PARAMS = 350,
     PARAMETERS = 351,
     PERIODS = 352,
     PLANNER_OBJECTIVE = 353,
     PREFILTER = 354,
     PRESAMPLE = 355,
     PRINT = 356,
     PRIOR_TRUNC = 357,
     PRIOR_ANALYSIS = 358,
     POSTERIOR_ANALYSIS = 359,
     QZ_CRITERIUM = 360,
     RELATIVE_IRF = 361,
     REPLIC = 362,
     RPLOT = 363,
     SHOCKS = 364,
     SIGMA_E = 365,
     SIMUL = 366,
     SIMUL_ALGO = 367,
     SIMUL_SEED = 368,
     SMOOTHER = 369,
     SOLVE_ALGO = 370,
     SPARSE_DLL = 371,
     STDERR = 372,
     STEADY = 373,
     STOCH_SIMUL = 374,
     TEX = 375,
     RAMSEY_POLICY = 376,
     PLANNER_DISCOUNT = 377,
     TEX_NAME = 378,
     UNIFORM_PDF = 379,
     UNIT_ROOT_VARS = 380,
     USE_DLL = 381,
     VALUES = 382,
     VAR = 383,
     VAREXO = 384,
     VAREXO_DET = 385,
     VAROBS = 386,
     XLS_SHEET = 387,
     XLS_RANGE = 388,
     COMMA = 389,
     MINUS = 390,
     PLUS = 391,
     DIVIDE = 392,
     TIMES = 393,
     UMINUS = 394,
     POWER = 395,
     EXP = 396,
     LOG = 397,
     LOG10 = 398,
     SIN = 399,
     COS = 400,
     TAN = 401,
     ASIN = 402,
     ACOS = 403,
     ATAN = 404,
     SINH = 405,
     COSH = 406,
     TANH = 407,
     ASINH = 408,
     ACOSH = 409,
     ATANH = 410,
     SQRT = 411
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

    static const short int yycheck_[];

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
