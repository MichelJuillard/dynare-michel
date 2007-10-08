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
     DUMMY = 284,
     DYNASAVE = 285,
     DYNATYPE = 286,
     END = 287,
     ENDVAL = 288,
     EQUAL = 289,
     ESTIMATION = 290,
     ESTIMATED_PARAMS = 291,
     ESTIMATED_PARAMS_BOUNDS = 292,
     ESTIMATED_PARAMS_INIT = 293,
     FILENAME = 294,
     FILTER_STEP_AHEAD = 295,
     FILTERED_VARS = 296,
     FIRST_OBS = 297,
     FLOAT_NUMBER = 298,
     FORECAST = 299,
     GAMMA_PDF = 300,
     GCC_COMPILER = 301,
     GRAPH = 302,
     HISTVAL = 303,
     HP_FILTER = 304,
     HP_NGRID = 305,
     INITVAL = 306,
     INT_NUMBER = 307,
     INV_GAMMA_PDF = 308,
     IRF = 309,
     KALMAN_ALGO = 310,
     KALMAN_TOL = 311,
     LAPLACE = 312,
     LCC_COMPILER = 313,
     LIK_ALGO = 314,
     LIK_INIT = 315,
     LINEAR = 316,
     LOAD_MH_FILE = 317,
     LOGLINEAR = 318,
     MARKOWITZ = 319,
     MAX = 320,
     MH_DROP = 321,
     MH_INIT_SCALE = 322,
     MH_JSCALE = 323,
     MH_MODE = 324,
     MH_NBLOCKS = 325,
     MH_REPLIC = 326,
     MH_RECOVER = 327,
     MIN = 328,
     MODE_CHECK = 329,
     MODE_COMPUTE = 330,
     MODE_FILE = 331,
     MODEL = 332,
     MODEL_COMPARISON = 333,
     MSHOCKS = 334,
     MODEL_COMPARISON_APPROXIMATION = 335,
     MODIFIEDHARMONICMEAN = 336,
     MOMENTS_VARENDO = 337,
     NAME = 338,
     NO_COMPILER = 339,
     NOBS = 340,
     NOCONSTANT = 341,
     NOCORR = 342,
     NODIAGNOSTIC = 343,
     NOFUNCTIONS = 344,
     NOGRAPH = 345,
     NOMOMENTS = 346,
     NOPRINT = 347,
     NORMAL_PDF = 348,
     OBSERVATION_TRENDS = 349,
     OPTIM = 350,
     OPTIM_WEIGHTS = 351,
     ORDER = 352,
     OSR = 353,
     OSR_PARAMS = 354,
     PARAMETERS = 355,
     PERIODS = 356,
     PLANNER_OBJECTIVE = 357,
     PREFILTER = 358,
     PRESAMPLE = 359,
     PRINT = 360,
     PRIOR_TRUNC = 361,
     PRIOR_ANALYSIS = 362,
     POSTERIOR_ANALYSIS = 363,
     QZ_CRITERIUM = 364,
     RELATIVE_IRF = 365,
     REPLIC = 366,
     RPLOT = 367,
     SHOCKS = 368,
     SIGMA_E = 369,
     SIMUL = 370,
     SIMUL_ALGO = 371,
     SIMUL_SEED = 372,
     SMOOTHER = 373,
     SOLVE_ALGO = 374,
     SPARSE_DLL = 375,
     STDERR = 376,
     STEADY = 377,
     STOCH_SIMUL = 378,
     TEX = 379,
     RAMSEY_POLICY = 380,
     PLANNER_DISCOUNT = 381,
     TEX_NAME = 382,
     UNIFORM_PDF = 383,
     UNIT_ROOT_VARS = 384,
     USE_DLL = 385,
     VALUES = 386,
     VAR = 387,
     VAREXO = 388,
     VAREXO_DET = 389,
     VAROBS = 390,
     XLS_SHEET = 391,
     XLS_RANGE = 392,
     COMMA = 393,
     GREATER = 394,
     LESS = 395,
     EXCLAMATION_EQUAL = 396,
     EXCLAMATION = 397,
     EQUAL_EQUAL = 398,
     GREATER_EQUAL = 399,
     LESS_EQUAL = 400,
     MINUS = 401,
     PLUS = 402,
     DIVIDE = 403,
     TIMES = 404,
     UMINUS = 405,
     POWER = 406,
     EXP = 407,
     LOG = 408,
     LOG10 = 409,
     SIN = 410,
     COS = 411,
     TAN = 412,
     ASIN = 413,
     ACOS = 414,
     ATAN = 415,
     SINH = 416,
     COSH = 417,
     TANH = 418,
     ASINH = 419,
     ACOSH = 420,
     ATANH = 421,
     SQRT = 422,
     DYNARE_SENSITIVITY = 423,
     IDENTIFICATION = 424,
     MORRIS = 425,
     STAB = 426,
     REDFORM = 427,
     PPRIOR = 428,
     PRIOR_RANGE = 429,
     PPOST = 430,
     ILPTAU = 431,
     GLUE = 432,
     MORRIS_NLIV = 433,
     MORRIS_NTRA = 434,
     NSAM = 435,
     LOAD_REDFORM = 436,
     LOAD_RMSE = 437,
     LOAD_STAB = 438,
     ALPHA2_STAB = 439,
     KSSTAT = 440,
     LOGTRANS_REDFORM = 441,
     THRESHOLD_REDFORM = 442,
     KSSTAT_REDFORM = 443,
     ALPHA2_REDFORM = 444,
     NAMENDO = 445,
     NAMLAGENDO = 446,
     NAMEXO = 447,
     RMSE = 448,
     LIK_ONLY = 449,
     VAR_RMSE = 450,
     PFILT_RMSE = 451,
     ISTART_RMSE = 452,
     ALPHA_RMSE = 453,
     ALPHA2_RMSE = 454
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
