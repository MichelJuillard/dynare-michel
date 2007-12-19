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
     BICGSTAB = 262,
     BVAR_DENSITY = 263,
     BVAR_FORECAST = 264,
     BVAR_PRIOR_DECAY = 265,
     BVAR_PRIOR_FLAT = 266,
     BVAR_PRIOR_LAMBDA = 267,
     BVAR_PRIOR_MU = 268,
     BVAR_PRIOR_OMEGA = 269,
     BVAR_PRIOR_TAU = 270,
     BVAR_PRIOR_TRAIN = 271,
     BVAR_REPLIC = 272,
     CALIB = 273,
     CALIB_VAR = 274,
     CHECK = 275,
     CONF_SIG = 276,
     CONSTANT = 277,
     CORR = 278,
     COVAR = 279,
     CUTOFF = 280,
     DATAFILE = 281,
     DR_ALGO = 282,
     DROP = 283,
     DSAMPLE = 284,
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
     GAUSSIAN_ELIMINATION = 301,
     GCC_COMPILER = 302,
     GMRES = 303,
     GRAPH = 304,
     HISTVAL = 305,
     HP_FILTER = 306,
     HP_NGRID = 307,
     INITVAL = 308,
     INT_NUMBER = 309,
     INV_GAMMA_PDF = 310,
     IRF = 311,
     KALMAN_ALGO = 312,
     KALMAN_TOL = 313,
     LAPLACE = 314,
     LCC_COMPILER = 315,
     LIK_ALGO = 316,
     LIK_INIT = 317,
     LINEAR = 318,
     LOAD_MH_FILE = 319,
     LOGLINEAR = 320,
     LU = 321,
     MARKOWITZ = 322,
     MAX = 323,
     METHOD = 324,
     MH_DROP = 325,
     MH_INIT_SCALE = 326,
     MH_JSCALE = 327,
     MH_MODE = 328,
     MH_NBLOCKS = 329,
     MH_REPLIC = 330,
     MH_RECOVER = 331,
     MIN = 332,
     MODE_CHECK = 333,
     MODE_COMPUTE = 334,
     MODE_FILE = 335,
     MODEL = 336,
     MODEL_COMPARISON = 337,
     MSHOCKS = 338,
     MODEL_COMPARISON_APPROXIMATION = 339,
     MODIFIEDHARMONICMEAN = 340,
     MOMENTS_VARENDO = 341,
     NAME = 342,
     NO_COMPILER = 343,
     NOBS = 344,
     NOCONSTANT = 345,
     NOCORR = 346,
     NODIAGNOSTIC = 347,
     NOFUNCTIONS = 348,
     NOGRAPH = 349,
     NOMOMENTS = 350,
     NOPRINT = 351,
     NORMAL_PDF = 352,
     OBSERVATION_TRENDS = 353,
     OPTIM = 354,
     OPTIM_WEIGHTS = 355,
     ORDER = 356,
     OSR = 357,
     OSR_PARAMS = 358,
     PARAMETERS = 359,
     PERIODS = 360,
     PLANNER_OBJECTIVE = 361,
     PREFILTER = 362,
     PRESAMPLE = 363,
     PRINT = 364,
     PRIOR_TRUNC = 365,
     PRIOR_ANALYSIS = 366,
     POSTERIOR_ANALYSIS = 367,
     QZ_CRITERIUM = 368,
     RELATIVE_IRF = 369,
     REPLIC = 370,
     RPLOT = 371,
     SHOCKS = 372,
     SIGMA_E = 373,
     SIMUL = 374,
     SIMUL_ALGO = 375,
     SIMUL_SEED = 376,
     SMOOTHER = 377,
     SOLVE_ALGO = 378,
     SPARSE = 379,
     SPARSE_DLL = 380,
     STDERR = 381,
     STEADY = 382,
     STOCH_SIMUL = 383,
     TEX = 384,
     RAMSEY_POLICY = 385,
     PLANNER_DISCOUNT = 386,
     TEX_NAME = 387,
     UNIFORM_PDF = 388,
     UNIT_ROOT_VARS = 389,
     USE_DLL = 390,
     VALUES = 391,
     VAR = 392,
     VAREXO = 393,
     VAREXO_DET = 394,
     VAROBS = 395,
     XLS_SHEET = 396,
     XLS_RANGE = 397,
     NORMCDF = 398,
     HOMOTOPY_SETUP = 399,
     HOMOTOPY_MODE = 400,
     HOMOTOPY_STEPS = 401,
     EXCLAMATION_EQUAL = 402,
     EXCLAMATION = 403,
     EQUAL_EQUAL = 404,
     GREATER_EQUAL = 405,
     LESS_EQUAL = 406,
     GREATER = 407,
     LESS = 408,
     COMMA = 409,
     MINUS = 410,
     PLUS = 411,
     DIVIDE = 412,
     TIMES = 413,
     UMINUS = 414,
     POWER = 415,
     EXP = 416,
     LOG = 417,
     LN = 418,
     LOG10 = 419,
     SIN = 420,
     COS = 421,
     TAN = 422,
     ASIN = 423,
     ACOS = 424,
     ATAN = 425,
     SINH = 426,
     COSH = 427,
     TANH = 428,
     ASINH = 429,
     ACOSH = 430,
     ATANH = 431,
     SQRT = 432,
     DYNARE_SENSITIVITY = 433,
     IDENTIFICATION = 434,
     MORRIS = 435,
     STAB = 436,
     REDFORM = 437,
     PPRIOR = 438,
     PRIOR_RANGE = 439,
     PPOST = 440,
     ILPTAU = 441,
     GLUE = 442,
     MORRIS_NLIV = 443,
     MORRIS_NTRA = 444,
     NSAM = 445,
     LOAD_REDFORM = 446,
     LOAD_RMSE = 447,
     LOAD_STAB = 448,
     ALPHA2_STAB = 449,
     KSSTAT = 450,
     LOGTRANS_REDFORM = 451,
     THRESHOLD_REDFORM = 452,
     KSSTAT_REDFORM = 453,
     ALPHA2_REDFORM = 454,
     NAMENDO = 455,
     NAMLAGENDO = 456,
     NAMEXO = 457,
     RMSE = 458,
     LIK_ONLY = 459,
     VAR_RMSE = 460,
     PFILT_RMSE = 461,
     ISTART_RMSE = 462,
     ALPHA_RMSE = 463,
     ALPHA2_RMSE = 464
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
