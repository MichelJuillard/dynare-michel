/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison LALR(1) parsers in C++

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


#include "DynareBison.hh"

/* User implementation prologue.  */
#line 37 "DynareBison.yy"

#include "ParsingDriver.hh"


/* Line 317 of lalr1.cc.  */
#line 46 "DynareBison.cc"

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* FIXME: INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#define YYUSE(e) ((void) (e))

/* A pseudo ostream that takes yydebug_ into account.  */
# define YYCDEBUG							\
  for (bool yydebugcond_ = yydebug_; yydebugcond_; yydebugcond_ = false)	\
    (*yycdebug_)

/* Enable debugging if requested.  */
#if YYDEBUG

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)	\
do {							\
  if (yydebug_)						\
    {							\
      *yycdebug_ << Title << ' ';			\
      yy_symbol_print_ ((Type), (Value), (Location));	\
      *yycdebug_ << std::endl;				\
    }							\
} while (false)

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug_)				\
    yy_reduce_print_ (Rule);		\
} while (false)

# define YY_STACK_PRINT()		\
do {					\
  if (yydebug_)				\
    yystack_print_ ();			\
} while (false)

#else /* !YYDEBUG */

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_REDUCE_PRINT(Rule)
# define YY_STACK_PRINT()

#endif /* !YYDEBUG */

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab

namespace yy
{
#if YYERROR_VERBOSE

  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  parser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr = "";
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              /* Fall through.  */
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }

#endif

  /// Build a parser object.
  parser::parser (ParsingDriver &driver_yyarg)
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
      driver (driver_yyarg)
  {
  }

  parser::~parser ()
  {
  }

#if YYDEBUG
  /*--------------------------------.
  | Print this symbol on YYOUTPUT.  |
  `--------------------------------*/

  inline void
  parser::yy_symbol_value_print_ (int yytype,
			   const semantic_type* yyvaluep, const location_type* yylocationp)
  {
    YYUSE (yylocationp);
    YYUSE (yyvaluep);
    switch (yytype)
      {
         default:
	  break;
      }
  }


  void
  parser::yy_symbol_print_ (int yytype,
			   const semantic_type* yyvaluep, const location_type* yylocationp)
  {
    *yycdebug_ << (yytype < yyntokens_ ? "token" : "nterm")
	       << ' ' << yytname_[yytype] << " ("
	       << *yylocationp << ": ";
    yy_symbol_value_print_ (yytype, yyvaluep, yylocationp);
    *yycdebug_ << ')';
  }
#endif /* ! YYDEBUG */

  void
  parser::yydestruct_ (const char* yymsg,
			   int yytype, semantic_type* yyvaluep, location_type* yylocationp)
  {
    YYUSE (yylocationp);
    YYUSE (yymsg);
    YYUSE (yyvaluep);

    YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

    switch (yytype)
      {
  
	default:
	  break;
      }
  }

  void
  parser::yypop_ (unsigned int n)
  {
    yystate_stack_.pop (n);
    yysemantic_stack_.pop (n);
    yylocation_stack_.pop (n);
  }

  std::ostream&
  parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  parser::debug_level_type
  parser::debug_level () const
  {
    return yydebug_;
  }

  void
  parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }


  int
  parser::parse ()
  {
    /// Look-ahead and look-ahead in internal form.
    int yychar = yyempty_;
    int yytoken = 0;

    /* State.  */
    int yyn;
    int yylen = 0;
    int yystate = 0;

    /* Error handling.  */
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// Semantic value of the look-ahead.
    semantic_type yylval;
    /// Location of the look-ahead.
    location_type yylloc;
    /// The locations where the error started and ended.
    location yyerror_range[2];

    /// $$.
    semantic_type yyval;
    /// @$.
    location_type yyloc;

    int yyresult;

    YYCDEBUG << "Starting parse" << std::endl;


    /* User initialization code.  */
    #line 22 "DynareBison.yy"
{
  // Initialize the location filenames
  yylloc.begin.filename = yylloc.end.filename = &driver.file;
}
  /* Line 555 of yacc.c.  */
#line 283 "DynareBison.cc"
    /* Initialize the stacks.  The initial state will be pushed in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystate_stack_ = state_stack_type (0);
    yysemantic_stack_ = semantic_stack_type (0);
    yylocation_stack_ = location_stack_type (0);
    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yylloc);

    /* New state.  */
  yynewstate:
    yystate_stack_.push (yystate);
    YYCDEBUG << "Entering state " << yystate << std::endl;
    goto yybackup;

    /* Backup.  */
  yybackup:

    /* Try to take a decision without look-ahead.  */
    yyn = yypact_[yystate];
    if (yyn == yypact_ninf_)
      goto yydefault;

    /* Read a look-ahead token.  */
    if (yychar == yyempty_)
      {
	YYCDEBUG << "Reading a token: ";
	yychar = yylex (&yylval, &yylloc, driver);
      }


    /* Convert token to internal form.  */
    if (yychar <= yyeof_)
      {
	yychar = yytoken = yyeof_;
	YYCDEBUG << "Now at end of input." << std::endl;
      }
    else
      {
	yytoken = yytranslate_ (yychar);
	YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
      }

    /* If the proper action on seeing token YYTOKEN is to reduce or to
       detect an error, take that action.  */
    yyn += yytoken;
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yytoken)
      goto yydefault;

    /* Reduce or error.  */
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
	if (yyn == 0 || yyn == yytable_ninf_)
	goto yyerrlab;
	yyn = -yyn;
	goto yyreduce;
      }

    /* Accept?  */
    if (yyn == yyfinal_)
      goto yyacceptlab;

    /* Shift the look-ahead token.  */
    YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

    /* Discard the token being shifted unless it is eof.  */
    if (yychar != yyeof_)
      yychar = yyempty_;

    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yylloc);

    /* Count tokens shifted since error; after three, turn off error
       status.  */
    if (yyerrstatus_)
      --yyerrstatus_;

    yystate = yyn;
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystate];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    /* If YYLEN is nonzero, implement the default value of the action:
       `$$ = $1'.  Otherwise, use the top of the stack.

       Otherwise, the following line sets YYVAL to garbage.
       This behavior is undocumented and Bison
       users should not rely upon it.  */
    if (yylen)
      yyval = yysemantic_stack_[yylen - 1];
    else
      yyval = yysemantic_stack_[0];

    {
      slice<location_type, location_stack_type> slice (yylocation_stack_, yylen);
      YYLLOC_DEFAULT (yyloc, slice, yylen);
    }
    YY_REDUCE_PRINT (yyn);
    switch (yyn)
      {
	  case 45:
#line 142 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 46:
#line 143 "DynareBison.yy"
    {driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 47:
#line 146 "DynareBison.yy"
    {driver.rplot();;}
    break;

  case 52:
#line 166 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 53:
#line 168 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 54:
#line 170 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 55:
#line 172 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 56:
#line 174 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 57:
#line 176 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 58:
#line 181 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 59:
#line 183 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 60:
#line 185 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 61:
#line 187 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 62:
#line 189 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 63:
#line 191 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 64:
#line 196 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 65:
#line 198 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 66:
#line 200 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 67:
#line 202 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 68:
#line 204 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 69:
#line 206 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 70:
#line 211 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 71:
#line 213 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 72:
#line 215 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 73:
#line 217 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 74:
#line 219 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 75:
#line 221 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 76:
#line 226 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 77:
#line 230 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 78:
#line 238 "DynareBison.yy"
    {driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 79:
#line 243 "DynareBison.yy"
    { (yyval.exp_val) = (yysemantic_stack_[(3) - (2)].exp_val);;}
    break;

  case 80:
#line 245 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 81:
#line 247 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 82:
#line 249 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 83:
#line 251 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::PLUS);;}
    break;

  case 84:
#line 253 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::MINUS);;}
    break;

  case 85:
#line 255 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::DIVIDE);;}
    break;

  case 86:
#line 257 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::TIMES);;}
    break;

  case 87:
#line 259 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::POWER);;}
    break;

  case 88:
#line 261 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(2) - (2)].exp_val), token::UMINUS);;}
    break;

  case 89:
#line 263 "DynareBison.yy"
    {(yyval.exp_val) = (yysemantic_stack_[(2) - (2)].exp_val);;}
    break;

  case 90:
#line 265 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::EXP);;}
    break;

  case 91:
#line 267 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::LOG);;}
    break;

  case 92:
#line 269 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::LOG10);;}
    break;

  case 93:
#line 271 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::SIN);;}
    break;

  case 94:
#line 273 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::COS);;}
    break;

  case 95:
#line 275 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::TAN);;}
    break;

  case 96:
#line 277 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::ASIN);;}
    break;

  case 97:
#line 279 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::ACOS);;}
    break;

  case 98:
#line 281 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::ATAN);;}
    break;

  case 99:
#line 283 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::SQRT);;}
    break;

  case 100:
#line 285 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), (yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 101:
#line 287 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), (yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 102:
#line 292 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::COMMA);;}
    break;

  case 103:
#line 294 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::COMMA);;}
    break;

  case 104:
#line 298 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 105:
#line 303 "DynareBison.yy"
    {driver.end_endval();;}
    break;

  case 108:
#line 313 "DynareBison.yy"
    {driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 109:
#line 318 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 112:
#line 328 "DynareBison.yy"
    {driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 113:
#line 332 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 115:
#line 333 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 117:
#line 335 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 123:
#line 348 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].model_val), (yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 124:
#line 350 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].model_val));;}
    break;

  case 125:
#line 354 "DynareBison.yy"
    {(yyval.model_val) = (yysemantic_stack_[(3) - (2)].model_val);;}
    break;

  case 127:
#line 357 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 128:
#line 359 "DynareBison.yy"
    {(yysemantic_stack_[(1) - (1)].string_val)->append(".0"); (yyval.model_val) = driver.add_model_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 129:
#line 361 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_plus((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 130:
#line 363 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_minus((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 131:
#line 365 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_divide((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 132:
#line 367 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_times((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 133:
#line 369 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_power((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 134:
#line 371 "DynareBison.yy"
    { (yyval.model_val) = driver.add_model_uminus((yysemantic_stack_[(2) - (2)].model_val));;}
    break;

  case 135:
#line 373 "DynareBison.yy"
    {(yyval.model_val) = (yysemantic_stack_[(2) - (2)].model_val);;}
    break;

  case 136:
#line 375 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_exp((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 137:
#line 377 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_log((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 138:
#line 379 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_log10((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 139:
#line 381 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_sin((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 140:
#line 383 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_cos((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 141:
#line 385 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_tan((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 142:
#line 387 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_asin((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 143:
#line 389 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_acos((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 144:
#line 391 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_atan((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 145:
#line 393 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_sqrt((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 146:
#line 397 "DynareBison.yy"
    {driver.declare_and_init_local_parameter((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].model_val));;}
    break;

  case 147:
#line 401 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 148:
#line 403 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 149:
#line 407 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 150:
#line 411 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 153:
#line 421 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 154:
#line 423 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].exp_val));;}
    break;

  case 155:
#line 425 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 156:
#line 427 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 157:
#line 429 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 158:
#line 434 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 159:
#line 436 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 160:
#line 438 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 161:
#line 440 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 162:
#line 442 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 163:
#line 444 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 164:
#line 450 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 165:
#line 452 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 166:
#line 454 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 167:
#line 456 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 168:
#line 458 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 169:
#line 460 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 170:
#line 462 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 171:
#line 464 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (2)].exp_val));;}
    break;

  case 172:
#line 469 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 173:
#line 474 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 174:
#line 476 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 175:
#line 481 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 176:
#line 483 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 177:
#line 485 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 178:
#line 487 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 179:
#line 489 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 180:
#line 491 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 181:
#line 493 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].exp_val));;}
    break;

  case 182:
#line 495 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 183:
#line 497 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 184:
#line 502 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 185:
#line 506 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 189:
#line 518 "DynareBison.yy"
    {driver.check();;}
    break;

  case 190:
#line 520 "DynareBison.yy"
    {driver.check();;}
    break;

  case 194:
#line 532 "DynareBison.yy"
    {driver.simul();;}
    break;

  case 195:
#line 534 "DynareBison.yy"
    {driver.simul();;}
    break;

  case 199:
#line 546 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 200:
#line 548 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 201:
#line 550 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 202:
#line 552 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 224:
#line 582 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 225:
#line 584 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 226:
#line 586 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 227:
#line 588 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 228:
#line 590 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 229:
#line 592 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 230:
#line 597 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 231:
#line 599 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 232:
#line 601 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 233:
#line 606 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 234:
#line 608 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 235:
#line 610 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 236:
#line 615 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 237:
#line 620 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 238:
#line 622 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 240:
#line 631 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 241:
#line 636 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 242:
#line 641 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 243:
#line 651 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 244:
#line 656 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 245:
#line 662 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(9) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(9) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(9) - (5)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(9) - (7)].string_val);
     delete (yysemantic_stack_[(9) - (1)].string_val);
     delete (yysemantic_stack_[(9) - (3)].string_val);
     delete (yysemantic_stack_[(9) - (5)].string_val);
     delete (yysemantic_stack_[(9) - (7)].string_val);
		;}
    break;

  case 246:
#line 672 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 247:
#line 677 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 248:
#line 688 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 249:
#line 694 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 250:
#line 702 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(7) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(7) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(7) - (5)].string_val);
		 driver.estim_params.p4=*(yysemantic_stack_[(7) - (7)].string_val);
     delete (yysemantic_stack_[(7) - (1)].string_val);
     delete (yysemantic_stack_[(7) - (3)].string_val);
     delete (yysemantic_stack_[(7) - (5)].string_val);
     delete (yysemantic_stack_[(7) - (7)].string_val);
		;}
    break;

  case 251:
#line 712 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(9) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(9) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(9) - (5)].string_val);
		 driver.estim_params.p4=*(yysemantic_stack_[(9) - (7)].string_val);
		 driver.estim_params.jscale=*(yysemantic_stack_[(9) - (9)].string_val);
     delete (yysemantic_stack_[(9) - (1)].string_val);
     delete (yysemantic_stack_[(9) - (3)].string_val);
     delete (yysemantic_stack_[(9) - (5)].string_val);
     delete (yysemantic_stack_[(9) - (7)].string_val);
     delete (yysemantic_stack_[(9) - (9)].string_val);
		;}
    break;

  case 252:
#line 726 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 253:
#line 730 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 254:
#line 732 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 255:
#line 736 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 256:
#line 743 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 257:
#line 752 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 258:
#line 761 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 259:
#line 765 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 260:
#line 767 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 261:
#line 771 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 262:
#line 780 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(9) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(9) - (4)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(9) - (6)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(9) - (8)].string_val);
         delete (yysemantic_stack_[(9) - (2)].string_val);
         delete (yysemantic_stack_[(9) - (4)].string_val);
         delete (yysemantic_stack_[(9) - (6)].string_val);
         delete (yysemantic_stack_[(9) - (8)].string_val);
				;}
    break;

  case 263:
#line 791 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 264:
#line 803 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 265:
#line 805 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 266:
#line 807 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 267:
#line 809 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 268:
#line 811 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 269:
#line 815 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 273:
#line 820 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 274:
#line 822 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 275:
#line 829 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 276:
#line 831 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 277:
#line 833 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 278:
#line 835 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 320:
#line 886 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 321:
#line 888 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 337:
#line 914 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 338:
#line 916 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 339:
#line 920 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 340:
#line 921 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 343:
#line 931 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 344:
#line 936 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 347:
#line 945 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].exp_val));;}
    break;

  case 348:
#line 948 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 349:
#line 952 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 350:
#line 956 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 351:
#line 958 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].exp_val));;}
    break;

  case 352:
#line 960 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].exp_val));;}
    break;

  case 353:
#line 962 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 354:
#line 965 "DynareBison.yy"
    {driver.set_osr_params();;}
    break;

  case 355:
#line 968 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 356:
#line 969 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 357:
#line 970 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 358:
#line 971 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 359:
#line 974 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 360:
#line 975 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 361:
#line 976 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 362:
#line 977 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 367:
#line 988 "DynareBison.yy"
    {driver.set_olr_inst();;}
    break;

  case 368:
#line 992 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 371:
#line 999 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 372:
#line 1000 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 373:
#line 1001 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].exp_val));;}
    break;

  case 374:
#line 1004 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 375:
#line 1005 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 376:
#line 1006 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 377:
#line 1009 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 378:
#line 1010 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 379:
#line 1013 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 380:
#line 1014 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 381:
#line 1015 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 382:
#line 1016 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 383:
#line 1017 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 384:
#line 1018 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 385:
#line 1020 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 386:
#line 1021 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 387:
#line 1022 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 388:
#line 1023 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 389:
#line 1024 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 390:
#line 1025 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 391:
#line 1028 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 397:
#line 1040 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 398:
#line 1041 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 399:
#line 1042 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 400:
#line 1043 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 401:
#line 1046 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 402:
#line 1047 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 404:
#line 1051 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 405:
#line 1052 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 406:
#line 1053 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 407:
#line 1054 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 408:
#line 1057 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 409:
#line 1057 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].model_val)); ;}
    break;

  case 411:
#line 1061 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 412:
#line 1063 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 413:
#line 1065 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 414:
#line 1067 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 419:
#line 1079 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 420:
#line 1080 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 421:
#line 1081 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 422:
#line 1082 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 423:
#line 1083 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 424:
#line 1084 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 425:
#line 1085 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 426:
#line 1086 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 427:
#line 1087 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 428:
#line 1088 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 429:
#line 1089 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 430:
#line 1090 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 431:
#line 1091 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 432:
#line 1092 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 433:
#line 1093 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 434:
#line 1094 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 435:
#line 1095 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 436:
#line 1096 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 437:
#line 1097 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 438:
#line 1099 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 439:
#line 1100 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 440:
#line 1101 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 441:
#line 1103 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 442:
#line 1104 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 443:
#line 1105 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 444:
#line 1106 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 445:
#line 1107 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 446:
#line 1108 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 447:
#line 1109 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 448:
#line 1110 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 449:
#line 1111 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 450:
#line 1112 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 451:
#line 1113 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 453:
#line 1115 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 454:
#line 1116 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 455:
#line 1117 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 456:
#line 1118 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 457:
#line 1119 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 458:
#line 1120 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 459:
#line 1121 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 460:
#line 1122 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 461:
#line 1123 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 462:
#line 1124 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 463:
#line 1125 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 464:
#line 1126 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 465:
#line 1127 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 466:
#line 1128 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 467:
#line 1129 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 468:
#line 1130 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 469:
#line 1131 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 470:
#line 1132 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 471:
#line 1133 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 472:
#line 1134 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 473:
#line 1135 "DynareBison.yy"
    {driver.option_num("olr_beta", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 474:
#line 1138 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 475:
#line 1140 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 476:
#line 1142 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 477:
#line 1143 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 478:
#line 1144 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 479:
#line 1145 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 480:
#line 1146 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 481:
#line 1147 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 482:
#line 1148 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 483:
#line 1149 "DynareBison.yy"
    {driver.option_num("load_mh_file", "-1");;}
    break;

  case 484:
#line 1150 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 485:
#line 1153 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 487:
#line 1162 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 488:
#line 1165 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 489:
#line 1167 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 490:
#line 1175 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2094 "DynareBison.cc"
	default: break;
      }
    YY_SYMBOL_PRINT ("-> $$ =", yyr1_[yyn], &yyval, &yyloc);

    yypop_ (yylen);
    yylen = 0;
    YY_STACK_PRINT ();

    yysemantic_stack_.push (yyval);
    yylocation_stack_.push (yyloc);

    /* Shift the result of the reduction.  */
    yyn = yyr1_[yyn];
    yystate = yypgoto_[yyn - yyntokens_] + yystate_stack_[0];
    if (0 <= yystate && yystate <= yylast_
	&& yycheck_[yystate] == yystate_stack_[0])
      yystate = yytable_[yystate];
    else
      yystate = yydefgoto_[yyn - yyntokens_];
    goto yynewstate;

  /*------------------------------------.
  | yyerrlab -- here on detecting error |
  `------------------------------------*/
  yyerrlab:
    /* If not already recovering from an error, report this error.  */
    if (!yyerrstatus_)
      {
	++yynerrs_;
	error (yylloc, yysyntax_error_ (yystate, yytoken));
      }

    yyerror_range[0] = yylloc;
    if (yyerrstatus_ == 3)
      {
	/* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

	if (yychar <= yyeof_)
	  {
	  /* Return failure if at end of input.  */
	  if (yychar == yyeof_)
	    YYABORT;
	  }
	else
	  {
	    yydestruct_ ("Error: discarding", yytoken, &yylval, &yylloc);
	    yychar = yyempty_;
	  }
      }

    /* Else will try to reuse look-ahead token after shifting the error
       token.  */
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (false)
      goto yyerrorlab;

    yyerror_range[0] = yylocation_stack_[yylen - 1];
    /* Do not reclaim the symbols of the rule which action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    yystate = yystate_stack_[0];
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;	/* Each real token shifted decrements this.  */

    for (;;)
      {
	yyn = yypact_[yystate];
	if (yyn != yypact_ninf_)
	{
	  yyn += yyterror_;
	  if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_)
	    {
	      yyn = yytable_[yyn];
	      if (0 < yyn)
		break;
	    }
	}

	/* Pop the current state because it cannot handle the error token.  */
	if (yystate_stack_.height () == 1)
	YYABORT;

	yyerror_range[0] = yylocation_stack_[0];
	yydestruct_ ("Error: popping",
		     yystos_[yystate],
		     &yysemantic_stack_[0], &yylocation_stack_[0]);
	yypop_ ();
	yystate = yystate_stack_[0];
	YY_STACK_PRINT ();
      }

    if (yyn == yyfinal_)
      goto yyacceptlab;

    yyerror_range[1] = yylloc;
    // Using YYLLOC is tempting, but would change the location of
    // the look-ahead.  YYLOC is available though.
    YYLLOC_DEFAULT (yyloc, (yyerror_range - 1), 2);
    yysemantic_stack_.push (yylval);
    yylocation_stack_.push (yyloc);

    /* Shift the error token.  */
    YY_SYMBOL_PRINT ("Shifting", yystos_[yyn],
		   &yysemantic_stack_[0], &yylocation_stack_[0]);

    yystate = yyn;
    goto yynewstate;

    /* Accept.  */
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    /* Abort.  */
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (yychar != yyeof_ && yychar != yyempty_)
      yydestruct_ ("Cleanup: discarding lookahead", yytoken, &yylval, &yylloc);

    /* Do not reclaim the symbols of the rule which action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (yystate_stack_.height () != 1)
      {
	yydestruct_ ("Cleanup: popping",
		   yystos_[yystate_stack_[0]],
		   &yysemantic_stack_[0],
		   &yylocation_stack_[0]);
	yypop_ ();
      }

    return yyresult;
  }

  // Generate an error message.
  std::string
  parser::yysyntax_error_ (int yystate, int tok)
  {
    std::string res;
    YYUSE (yystate);
#if YYERROR_VERBOSE
    int yyn = yypact_[yystate];
    if (yypact_ninf_ < yyn && yyn <= yylast_)
      {
	/* Start YYX at -YYN if negative to avoid negative indexes in
	   YYCHECK.  */
	int yyxbegin = yyn < 0 ? -yyn : 0;

	/* Stay within bounds of both yycheck and yytname.  */
	int yychecklim = yylast_ - yyn + 1;
	int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
	int count = 0;
	for (int x = yyxbegin; x < yyxend; ++x)
	  if (yycheck_[x + yyn] == x && x != yyterror_)
	    ++count;

	// FIXME: This method of building the message is not compatible
	// with internationalization.  It should work like yacc.c does it.
	// That is, first build a string that looks like this:
	// "syntax error, unexpected %s or %s or %s"
	// Then, invoke YY_ on this string.
	// Finally, use the string as a format to output
	// yytname_[tok], etc.
	// Until this gets fixed, this message appears in English only.
	res = "syntax error, unexpected ";
	res += yytnamerr_ (yytname_[tok]);
	if (count < 5)
	  {
	    count = 0;
	    for (int x = yyxbegin; x < yyxend; ++x)
	      if (yycheck_[x + yyn] == x && x != yyterror_)
		{
		  res += (!count++) ? ", expecting " : " or ";
		  res += yytnamerr_ (yytname_[x]);
		}
	  }
      }
    else
#endif
      res = YY_("syntax error");
    return res;
  }


  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
  const short int parser::yypact_ninf_ = -864;
  const short int
  parser::yypact_[] =
  {
       803,   146,   -48,   190,    71,   182,   188,   -31,   -14,   -26,
     -24,    -9,    82,   133,   309,   185,   150,   103,   237,    77,
     286,   239,   204,   286,   329,    84,  -864,   282,   297,   286,
     311,   493,   462,   486,   209,   234,   286,   453,   474,   488,
     286,   517,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,   543,    64,  -864,   460,    24,
     131,   514,   137,   526,   556,   617,  -864,   643,    -6,    32,
     196,   197,   586,   556,  -864,   174,    14,    81,   363,   590,
    -864,   901,   217,   291,   591,  -864,   901,   292,   302,   548,
     367,   625,   520,   673,   412,   412,   368,    81,   518,  -864,
     583,  -864,   460,  -864,   952,   379,  -864,   857,   439,   440,
     559,   444,   563,   445,   567,   451,   483,  -864,  -864,   535,
     640,   -51,   260,  -864,   686,   -30,  -864,  -864,   571,  -864,
    -864,   654,   193,  -864,   655,   246,   701,    46,  -864,   659,
    -864,   704,  -864,   705,   706,  -864,   710,   712,  -864,   713,
     714,   715,   720,   721,  -864,  -864,   722,   726,   727,   733,
     735,   736,  -864,  -864,   737,   738,  -864,   739,  -864,  -864,
    -864,   742,   743,   744,   745,  -864,  -864,   746,   751,   -20,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
     767,   682,  -864,   729,  -864,   732,    52,  -864,   671,   748,
     677,   749,   203,  -864,   750,   678,   752,   240,  -864,   674,
      49,  -864,    54,   181,  -864,   679,   689,   778,  -864,  -864,
     127,  -864,  -864,  -864,  -864,   758,   773,    68,  -864,  -864,
    -864,   694,   363,   363,   700,   708,   724,   730,   731,   753,
     760,   761,   762,   765,   363,   821,   768,    58,  -864,   829,
     832,   844,   845,   848,  -864,  -864,  -864,   852,   855,   859,
     860,  -864,   862,  -864,   868,   869,  -864,  -864,   142,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,   241,    63,   149,  -864,  -864,  -864,   786,   847,  -864,
     770,  -864,  -864,  -864,   777,   673,   673,   779,   781,   783,
     789,   794,   795,   804,   814,   816,   817,   673,   663,  -864,
     172,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,   201,  -864,    78,    18,   219,
    -864,  -864,   336,  -864,  -864,   406,  -864,  -864,   931,  -864,
     423,  -864,  -864,  -864,  -864,  -864,   813,   898,  -864,  -864,
     854,   905,  -864,  -864,   864,   907,  -864,  -864,   837,   838,
     919,   254,   963,  -864,  -864,   951,   460,   846,  -864,   849,
     -11,   926,   865,   255,   932,   363,  -864,  -864,  -864,   973,
     950,   870,   979,   980,   983,   985,   993,   996,   998,  1002,
     357,  1012,  1006,  1010,  1016,  1022,  1008,    12,   922,  1037,
    1044,  1020,  1030,  1050,   643,   256,  1057,  1080,  1081,  -864,
    -864,  -864,    47,  1082,   220,  1083,  -864,  -864,  1103,   220,
    1105,  -864,  -864,   175,  -864,  -864,  -864,  1070,    17,  -864,
     132,  -864,  1011,  1018,   395,    14,   278,  1106,    36,  -864,
    -864,   363,  1015,   513,   363,   363,   363,   363,   363,   363,
     363,   363,   363,   363,   629,   363,   363,   363,   363,   363,
    -864,   363,  -864,  -864,  1137,  1144,  1157,  1168,  1175,   220,
    1190,  1196,   499,  1199,  1210,  1212,   901,   276,  1186,  1061,
    -864,   648,   281,  -864,  1147,  -864,   175,  1132,   555,   673,
     673,   673,   673,   673,   673,   673,   673,   673,   673,   719,
     673,   673,   673,   673,   673,  1117,   412,   293,   315,  -864,
    -864,  -864,   363,   461,    22,   583,  1128,   460,  1130,   952,
     326,  1250,   857,   369,  -864,  1172,  -864,  1173,  -864,  1189,
    -864,  1247,  1151,  1155,  1156,   363,  -864,  -864,  -864,  -864,
    -864,   497,  1158,  -864,  -864,   503,  1159,  1067,  -864,  -864,
    1262,     5,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  1152,  -864,  -864,  -864,  -864,  1161,  -864,  -864,  -864,
     505,  -864,  1241,  1242,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,   537,  1166,  1191,  1192,  1246,  1194,   220,  1251,
    1174,   220,  -864,  1279,  1281,  1176,  1298,  -864,  -864,  -864,
     673,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,   401,   184,  -864,  1256,   363,  1258,    20,   709,
     434,   725,   808,   853,   883,   889,   903,   916,   934,   941,
     947,  -864,   513,   513,  1015,  1015,  1197,   954,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,   521,   363,  -864,  1260,  1074,  -864,   522,
    -864,  1180,   961,   967,   974,   981,   987,   994,  1001,  1007,
    1014,  1021,  -864,   555,   555,  1132,  1132,  1197,  -864,  -864,
    -864,   531,  -864,   581,  1027,    18,  1183,  -864,  -864,    37,
     363,  -864,  -864,  -864,  -864,  -864,  -864,   585,  -864,  -864,
    -864,   595,  -864,  -864,  -864,  1182,  1307,  -864,  -864,  1085,
    -864,   376,  -864,   386,  -864,  1184,  -864,  -864,  -864,  1265,
    -864,   452,  1266,  -864,  -864,  -864,  -864,  -864,  -864,   220,
      47,  1213,   220,  1214,  1215,  -864,  1193,  -864,  -864,  1311,
     673,  1092,   181,   181,   278,  -864,   220,  -864,  1316,  1098,
    1317,  1302,   363,   363,  -864,   363,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  1198,  -864,  1108,
     363,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,    22,  -864,  -864,
    -864,   363,  1034,  -864,  -864,  1151,   363,  -864,  -864,   596,
    -864,   597,  1303,  1195,  1152,  -864,  -864,  -864,  1222,  1223,
    1224,   220,  1203,   220,   220,  -864,   363,  1116,  -864,    43,
      69,   207,  1202,   363,  -864,   363,  1201,    66,  1122,   734,
     734,  -864,  -864,  1131,  1041,  -864,  1328,  1140,  -864,  -864,
    -864,  1230,  -864,   220,   220,   220,  1231,  -864,  1209,  1211,
    1146,  -864,  -864,  -864,   220,  -864,  1154,  1164,  1318,  1206,
    1319,  1244,  -864,  -864,  -864,   363,  -864,    19,  1238,  -864,
    1239,   220,  -864,  -864,  -864,  1216,  -864,  -864,  -864,  1323,
    1217,    72,  1170,  1299,  -864,   220,   212,  1219,  -864,  -864,
    1329,  -864,  -864,   582,   598,   363,    62,  -864,  -864,  -864,
    1218,  1245,  1249,  -864,  -864,  -864,  -864,  1047,  -864,  -864,
     363,  -864,  -864,  -864,   220,   220,  -864,  1054,  1252,  -864,
    -864,   220,  -864
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   408,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     2,     4,    27,    28,    42,    43,    44,    41,
       5,    10,     7,     8,     9,     6,    11,    12,    13,    14,
      15,    16,    17,    21,    23,    22,    18,    19,    20,    24,
      25,    26,    29,    30,    31,    36,    37,    32,    33,    34,
      35,    38,    39,    40,   377,     0,     0,   189,     0,     0,
       0,     0,     0,     0,     0,   228,   275,     0,     0,     0,
       0,     0,     0,     0,   113,     0,     0,     0,     0,     0,
     359,     0,     0,     0,     0,   355,     0,     0,     0,    72,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   194,
       0,   184,     0,   199,     0,     0,   411,     0,     0,     0,
      54,     0,    60,     0,    66,     0,     0,     1,     3,     0,
       0,   374,     0,   370,     0,     0,   192,   193,     0,    45,
     387,     0,     0,   381,     0,     0,     0,     0,   107,     0,
     464,     0,   481,     0,     0,   469,     0,     0,   447,     0,
       0,     0,     0,     0,   461,   462,     0,     0,     0,     0,
       0,     0,   483,   457,     0,     0,   468,     0,   482,   463,
     446,     0,     0,     0,     0,   467,   465,     0,     0,     0,
     280,   316,   305,   281,   282,   283,   284,   285,   286,   287,
     288,   289,   290,   291,   292,   293,   294,   295,   296,   297,
     298,   299,   300,   301,   302,   303,   304,   306,   307,   308,
     309,   310,   311,   312,   313,   314,   315,   317,   318,   319,
     224,     0,   277,     0,   241,     0,     0,   238,     0,     0,
       0,     0,     0,   260,     0,     0,     0,     0,   254,     0,
       0,   111,     0,     0,   422,     0,     0,     0,   477,   476,
       0,   393,   394,   395,   396,     0,     0,     0,   152,    81,
      82,    80,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   346,     0,
       0,     0,     0,     0,   427,   428,   429,     0,     0,     0,
       0,   470,     0,   434,     0,     0,   364,   365,     0,   205,
     206,   207,   208,   209,   210,   211,   212,   213,   214,   215,
     216,   218,   219,   220,   221,   222,   223,   217,   363,   361,
     367,     0,     0,     0,   357,   354,    75,    70,     0,    51,
       0,    76,   127,   128,   147,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   409,   126,
       0,   323,   328,   324,   325,   326,   327,   329,   330,   331,
     332,   333,   334,   335,   336,     0,    47,     0,     0,     0,
     197,   198,     0,   187,   188,     0,   204,   201,     0,   417,
       0,   416,   418,   413,   348,    57,    52,     0,    48,    63,
      58,     0,    49,    69,    64,     0,    50,   343,     0,     0,
       0,     0,     0,   368,   369,     0,     0,     0,    46,     0,
       0,     0,     0,     0,     0,     0,   105,   106,   229,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   226,     0,   240,
     236,   237,   269,     0,   269,     0,   258,   259,     0,   269,
       0,   252,   253,     0,   109,   110,   104,     0,     0,   121,
       0,   122,     0,     0,     0,     0,     0,     0,     0,   150,
     151,     0,    88,    89,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      78,     0,   344,   345,     0,     0,     0,     0,     0,   269,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     349,     0,     0,    73,    71,    77,     0,   134,   135,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   149,
     182,   183,     0,     0,   174,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    55,    53,    61,    59,    67,    65,
     378,     0,   374,     0,     0,     0,   420,   191,   190,   390,
     385,     0,     0,   384,   379,     0,     0,     0,   448,   438,
       0,     0,   480,   441,   466,   430,   471,   472,   444,   445,
     450,   453,   454,   451,   459,   460,   449,   456,   455,   440,
     439,     0,   442,   443,   458,   478,     0,   479,   279,   276,
       0,   225,     0,     0,   264,   271,   265,   270,   267,   272,
     266,   268,     0,     0,     0,   246,     0,     0,   269,     0,
       0,   269,   232,     0,     0,     0,     0,   114,   119,   120,
       0,   124,   117,   115,   474,   475,   392,   403,   405,   406,
     407,   404,     0,   397,   401,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    79,    84,    83,    85,    86,    87,     0,   426,   419,
     425,   431,   432,   473,   423,   433,   437,   436,   424,   421,
     435,   366,   360,     0,     0,   352,     0,     0,   356,     0,
      74,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   125,   130,   129,   131,   132,   133,   410,   322,
     320,     0,   337,     0,     0,     0,     0,   179,   180,     0,
       0,   196,   195,   186,   185,   203,   200,     0,   484,   415,
     412,     0,    56,    62,    68,     0,     0,   376,   375,     0,
     386,     0,   380,     0,   108,   486,   488,   490,   489,     0,
     341,     0,     0,   278,   227,   242,   274,   273,   239,   269,
     269,     0,   269,     0,     0,   257,     0,   231,   230,     0,
       0,     0,     0,     0,     0,   391,   269,   402,     0,     0,
       0,     0,     0,     0,   100,     0,   101,    90,    91,    92,
      93,    94,    95,    96,    97,    98,    99,     0,   362,     0,
       0,   350,   358,   148,   136,   137,   138,   139,   140,   141,
     142,   143,   144,   145,   321,   338,   181,   173,   172,   176,
     177,     0,     0,   202,   414,   374,     0,   371,   388,     0,
     382,     0,     0,     0,     0,   452,   485,   243,     0,     0,
       0,   269,     0,   269,   269,   255,     0,     0,   123,     0,
       0,   398,     0,     0,   155,     0,   163,     0,     0,   102,
     103,   347,   353,     0,     0,   178,     0,     0,   389,   383,
     487,     0,   342,   269,   269,   269,     0,   263,     0,     0,
       0,   146,   118,   116,   269,   399,     0,     0,     0,   158,
       0,     0,   154,   351,   175,     0,   372,   269,   248,   244,
     247,   269,   261,   256,   112,     0,   157,   156,   162,     0,
     160,     0,     0,     0,   340,   269,     0,     0,   400,   159,
       0,   235,   169,     0,     0,     0,     0,   168,   167,   373,
       0,   249,     0,   262,   161,   234,   233,     0,   166,   153,
       0,   165,   164,   339,   269,   269,   171,     0,   250,   245,
     170,   269,   251
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -864,  -864,  1327,  -864,  -864,  -864,  -864,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -286,  -864,  -864,  -864,
    1268,  -128,  -864,  -864,  1102,  -864,  -864,  -864,  -864,  -178,
    -493,   -99,  -484,  -864,  -864,  -864,  1248,  -216,  -864,  -864,
    -864,  -864,   609,  -864,  -864,   790,  -864,  -864,   940,  -864,
    -864,   793,  -864,  -864,  -112,   -19,  -547,   403,  -864,  -864,
    1124,  -864,  -864,  -863,  -864,  -864,  1114,  -864,  -864,  1120,
    -777,  -444,  -864,  -864,   909,  -864,  1259,   809,  -864,   502,
    -864,  -864,  -864,  -864,  1084,  -864,  -864,  -864,  -864,  -864,
    -864,   841,  1272,  -864,  -864,  -864,  1237,  -573,  -864,  -864,
    -864,  -864,  -864,   885,  -864,   568,  -663,  -864,  -864,  -864,
    -864,  -864,   801,  -864,   -86,  -864,  1289,  -864,  -864,  -864,
    -864,  -864,  -864,  -864,   -89,  -864,  -864,  -107,  -864,  -864,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,   -85,   -84,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,  -864,   -78,  -864,
    -864,  -864,  -864,  -864,   -77,   -71,   -70,   -69,   -66,   -65,
    -864,  -864,  -864,  -864,  -864,  -864,  -864,   -63,   -56,   -55,
    -864,  -864,  -864,  -864,  -864,   774,  -864,   929
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    41,    42,    43,    44,    45,    46,    47,    48,    49,
     141,   143,   145,   120,    50,    51,   305,   700,    52,    53,
     167,   168,    54,   270,   271,    55,   273,   823,   822,   498,
     499,   500,   501,   379,    56,    57,   287,   288,   907,   976,
      58,   583,   584,    59,   402,   403,    60,   155,   156,    61,
     399,   400,    62,   405,   326,    98,   675,   978,    63,   256,
     257,   258,   663,   887,    64,   267,   268,    65,   262,   263,
     664,   888,    66,   209,   210,    67,   380,   381,    68,   800,
     801,    69,    70,   307,   308,    71,    72,   352,    73,    74,
      75,   327,   328,    76,    77,   152,   153,   432,    78,    79,
      80,    81,   280,   281,   692,   693,   694,    82,   123,   575,
      83,   410,   411,   329,   330,   331,   332,   333,   334,   335,
     336,   337,   338,   339,   340,   341,   342,   343,   344,   345,
     346,   213,   214,   215,   216,   217,   218,   219,   383,   384,
     222,   223,   224,   225,   226,   227,   228,   229,   385,   231,
     232,   233,   234,   235,   386,   387,   388,   389,   390,   391,
     347,   242,   243,   348,   282,   283,   284,   392,   393,   394,
     247,   248,   249,   412,   647,   796,   621,   622
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       112,   113,   157,   117,   118,   678,   512,   513,   212,   741,
     126,   211,   220,   221,   679,   135,   138,   139,   524,   230,
     236,   146,   406,   401,   378,   409,   237,   238,   239,   786,
     827,   240,   241,   889,   244,   382,   382,   677,   665,   447,
     667,   245,   246,   795,   253,   670,   404,   362,   580,   655,
     639,    95,   767,   654,    95,   363,   581,   657,   696,   250,
     768,   949,   158,   932,   253,   549,   446,   869,   150,   494,
     430,   510,   480,   362,   496,   870,   277,   655,   532,   656,
     285,   363,   364,   550,   659,   657,   658,   278,   509,   933,
     285,   436,   971,   285,   431,   723,    86,   254,   579,   362,
     672,   474,   971,   279,   939,   831,   121,   363,   364,    89,
     672,   166,   659,    94,   269,   251,   437,   254,    99,   166,
     100,   660,   122,   306,   832,   108,   475,   988,   551,   151,
      96,    97,   999,   610,   364,   101,   255,   972,   252,   365,
     366,   662,    95,   769,   447,   367,   368,   369,   370,   371,
     372,   373,   374,   375,   680,   797,   255,   697,   661,   617,
     376,   620,   377,   582,   497,   365,   366,   770,   159,   662,
     963,   367,   368,   369,   370,   371,   372,   373,   374,   375,
     698,   510,   871,   286,   973,   974,   376,   940,   377,   982,
     497,   365,   366,   286,   973,   974,   286,   367,   368,   369,
     370,   371,   372,   373,   374,   375,   989,   990,   259,   264,
     941,   362,   376,   672,   377,   259,   497,   975,   654,   363,
     274,   110,   111,   486,   813,   699,   102,   816,   701,   702,
     703,   704,   705,   706,   707,   708,   709,   710,   827,   712,
     713,   714,   715,   716,   656,   717,   364,    90,   505,   687,
     655,   658,   264,    92,   570,   571,   572,   573,   657,   574,
     491,   260,   265,   546,   150,   737,   557,   558,   260,    95,
     546,   289,   687,   506,    95,   160,   681,   103,   569,   290,
     433,   163,   250,   161,   603,   659,   660,   275,   547,   164,
      84,    85,   604,   576,   107,   552,   764,   673,   674,    95,
     261,   266,   916,   365,   366,   265,   291,   261,   688,   367,
     368,   369,   370,   371,   372,   373,   374,   375,   577,   789,
      95,    95,   576,   661,   376,   151,   377,    91,   497,   826,
     106,   688,   689,    93,    87,    88,   690,   691,   251,   440,
     585,    95,   662,   687,   266,   441,    95,   578,   115,   116,
     157,    95,   934,   133,   134,   689,   250,   250,    95,   690,
     691,   349,   548,   292,   293,   586,   890,   250,   892,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   136,   137,
      95,   109,   902,   114,   303,   212,   304,   631,   211,   220,
     221,    95,   443,   289,   119,   632,   230,   236,   444,   614,
     649,   290,   688,   237,   238,   239,   678,   678,   240,   241,
     829,   244,   251,   251,   977,   679,   679,   170,   245,   246,
     732,   611,   171,   251,   615,   738,   689,   124,   291,   991,
     690,   691,   357,   250,    95,   350,   354,   760,   684,   174,
     175,    95,   125,   177,   250,   178,   355,   926,   849,   928,
     929,    95,   179,   104,   105,   127,   650,   587,   685,   762,
     742,   743,   744,   745,   746,   747,   748,   749,   750,   751,
     776,   753,   754,   755,   756,   757,   196,   775,   401,   948,
     409,   950,   588,   200,   872,   292,   293,   382,   358,   251,
     955,   294,   295,   296,   297,   298,   299,   300,   301,   302,
     251,   404,   204,   964,   250,   250,   303,   967,   304,   416,
     420,   359,   396,   780,   205,   128,   424,   147,   140,   206,
     878,   981,   824,   407,     1,     2,     3,   589,   733,   726,
     880,   207,   208,   739,     4,     5,     6,   727,     7,   142,
       8,     9,    10,    11,   592,   825,   908,   909,   250,   910,
     998,    12,   590,   144,    13,   835,   149,  1002,   761,   763,
     251,   251,   250,   154,   913,   417,   421,   806,   250,   593,
     250,   777,   425,   884,   781,   807,    14,    15,    16,   162,
     836,   821,    17,   413,   414,   914,   250,   250,   418,   422,
     917,   165,    18,    19,    20,   426,   250,    21,   885,    22,
      23,    24,    25,    26,   251,   765,   129,   130,    27,    28,
     930,   766,   985,    29,    30,    31,    32,   936,   251,   937,
     817,   166,    33,    34,   251,    35,   251,   427,   986,    36,
     131,   132,    37,    38,    39,    40,   818,   527,   528,   169,
     529,   790,   251,   251,   899,   900,   250,   792,   170,   803,
     250,   269,   251,   171,   172,   306,   351,   173,   356,   962,
     250,   250,   250,   360,   361,   848,   852,   398,   319,   415,
     174,   175,   176,   419,   177,   864,   178,   423,   289,   572,
     573,   428,   574,   179,   180,   181,   290,   182,   183,   987,
     184,   185,   186,   187,   188,   189,   190,   191,   192,   193,
     194,   195,   251,   362,   997,   429,   251,   196,   435,   197,
     198,   363,   199,   291,   200,   438,   251,   251,   251,   439,
     442,   897,   201,   445,   448,   865,   449,   450,   451,   873,
     202,   203,   452,   204,   453,   454,   455,   456,   364,   874,
     918,   919,   457,   458,   459,   205,   154,   477,   460,   461,
     206,   525,   526,   527,   528,   462,   529,   463,   464,   465,
     466,   467,   207,   208,   468,   469,   470,   471,   472,   736,
     292,   293,   879,   473,   881,   711,   294,   295,   296,   297,
     298,   299,   300,   301,   302,   570,   571,   572,   573,   476,
     574,   303,   482,   304,   478,   365,   366,   479,   484,   489,
     504,   367,   368,   369,   370,   371,   372,   373,   374,   375,
       1,     2,     3,   483,   485,   488,   376,   490,   377,   493,
       4,     5,     6,   507,     7,   502,     8,     9,    10,    11,
     833,   525,   526,   527,   528,   503,   529,    12,   508,   511,
      13,   570,   571,   572,   573,   514,   574,   525,   526,   527,
     528,   534,   529,   515,   535,   834,   525,   526,   527,   528,
     309,   529,    14,    15,    16,   752,   536,   537,    17,   516,
     538,   837,   310,   311,   539,   517,   518,   540,    18,    19,
      20,   541,   542,    21,   543,    22,    23,    24,    25,    26,
     544,   545,   312,   313,    27,    28,   553,   179,   519,    29,
      30,    31,    32,   274,   309,   520,   521,   522,    33,    34,
     523,    35,   554,   531,   555,    36,   310,   311,    37,    38,
      39,    40,   556,   594,   559,   314,   560,   315,   561,   316,
     525,   526,   527,   528,   562,   529,   312,   313,   318,   563,
     564,   179,   319,   525,   526,   527,   528,   274,   529,   565,
     320,   321,   322,   591,   838,   309,   323,   324,   325,   566,
     154,   567,   568,   595,   596,   530,   408,   310,   311,   314,
     597,   315,   599,   316,   598,   525,   526,   527,   528,   317,
     529,   600,   318,   601,   602,   605,   319,   312,   313,   606,
     608,   612,   179,   609,   320,   321,   322,   616,   274,   839,
     323,   324,   325,   618,   154,   525,   526,   527,   528,   613,
     529,   525,   526,   527,   528,   619,   529,   623,   624,   620,
     314,   625,   315,   626,   316,   525,   526,   527,   528,   840,
     529,   627,   630,   318,   628,   841,   629,   319,   525,   526,
     527,   528,   633,   529,   634,   320,   321,   322,   635,   842,
     644,   323,   324,   325,   636,   154,   525,   526,   527,   528,
     637,   529,   843,   525,   526,   527,   528,   641,   529,   525,
     526,   527,   528,   638,   529,   642,   525,   526,   527,   528,
     844,   529,   643,   570,   571,   572,   573,   845,   574,   570,
     571,   572,   573,   846,   574,   645,   570,   571,   572,   573,
     847,   574,   652,   570,   571,   572,   573,   854,   574,   570,
     571,   572,   573,   855,   574,   646,   570,   571,   572,   573,
     856,   574,   651,   570,   571,   572,   573,   857,   574,   570,
     571,   572,   573,   858,   574,   676,   570,   571,   572,   573,
     859,   574,   529,   570,   571,   572,   573,   860,   574,   525,
     526,   527,   528,   861,   529,   682,   525,   526,   527,   528,
     862,   529,   683,   525,   526,   527,   528,   863,   529,   525,
     526,   527,   528,   866,   529,   718,   525,   526,   527,   528,
     915,   529,   719,   525,   526,   527,   528,   944,   529,   525,
     526,   527,   528,   996,   529,   720,   525,   526,   527,   528,
    1000,   529,   653,   666,   668,   735,   721,   525,   526,   527,
     528,   794,   529,   722,   570,   571,   572,   573,   851,   574,
     525,   526,   527,   528,   669,   529,   671,   695,   724,   877,
     525,   526,   527,   528,   725,   529,   898,   728,   570,   571,
     572,   573,   904,   574,   525,   526,   527,   528,   729,   529,
     730,   734,   912,   525,   526,   527,   528,   740,   529,   574,
     931,   758,   525,   526,   527,   528,   942,   529,   525,   526,
     527,   528,   772,   529,   774,   943,   525,   526,   527,   528,
     778,   529,   782,   783,   946,   785,   525,   526,   527,   528,
     954,   529,   525,   526,   527,   528,   431,   529,   956,   784,
     795,   787,   788,   799,   791,   793,   804,   805,   957,   802,
     808,   811,   809,   810,   979,   812,   814,   817,   815,   818,
     820,   828,   819,   830,    -1,   850,   853,   868,   875,   876,
     883,   886,   882,   896,   891,   893,   894,   895,   903,   905,
     906,   920,   911,   923,   924,   925,   921,   927,   935,   938,
     945,   947,   951,   952,   959,   953,   958,   960,   961,   965,
     966,   969,   968,   983,   980,   970,   994,   984,   148,   993,
     995,   272,   495,  1001,   867,   397,   607,   773,   771,   992,
     481,   492,   487,   648,   395,   759,   922,   731,   353,   434,
     686,   533,   901,   779,   276,   798,   640
  };

  /* YYCHECK.  */
  const unsigned short int
  parser::yycheck_[] =
  {
        19,    20,    88,    22,    23,   498,   292,   293,    97,   556,
      29,    97,    97,    97,   498,    34,    35,    36,   304,    97,
      97,    40,   134,   130,   123,   137,    97,    97,    97,   602,
     693,    97,    97,   810,    97,   124,   125,    20,   482,   167,
     484,    97,    97,    38,    12,   489,   132,    30,    30,    30,
      38,    65,    30,     6,    65,    38,    38,    38,    22,    65,
      38,   924,    38,    20,    12,   351,    20,    30,     4,    20,
     121,   287,    20,    30,    20,    38,    62,    30,    20,    32,
      12,    38,    65,    20,    65,    38,    39,    73,    20,    20,
      12,   121,    30,    12,   145,   539,   144,    65,    20,    30,
      38,   121,    30,    89,    38,    85,    22,    38,    65,    38,
      38,    65,    65,   144,    65,   121,   146,    65,   144,    65,
     144,    74,    38,    65,   104,    22,   146,    65,    65,    65,
     144,   145,   995,   144,    65,   144,   104,    65,   144,   122,
     123,   122,    65,   121,   272,   128,   129,   130,   131,   132,
     133,   134,   135,   136,    22,   150,   104,   121,   111,   445,
     143,   149,   145,   145,   147,   122,   123,   145,   144,   122,
     151,   128,   129,   130,   131,   132,   133,   134,   135,   136,
     144,   397,   145,   115,   122,   123,   143,   121,   145,   966,
     147,   122,   123,   115,   122,   123,   115,   128,   129,   130,
     131,   132,   133,   134,   135,   136,   144,   145,    12,    12,
     144,    30,   143,    38,   145,    12,   147,   145,     6,    38,
      46,   144,   145,    20,   668,   511,   144,   671,   514,   515,
     516,   517,   518,   519,   520,   521,   522,   523,   901,   525,
     526,   527,   528,   529,    32,   531,    65,    65,   121,    65,
      30,    39,    12,    65,   122,   123,   124,   125,    38,   127,
      20,    65,    65,   121,     4,   551,   365,   366,    65,    65,
     121,    30,    65,   146,    65,   144,   144,   144,   377,    38,
      20,   144,    65,   152,    30,    65,    74,   113,   146,   152,
     144,   145,    38,   121,   144,   146,   582,   122,   123,    65,
     104,   104,   875,   122,   123,    65,    65,   104,   124,   128,
     129,   130,   131,   132,   133,   134,   135,   136,   146,   605,
      65,    65,   121,   111,   143,    65,   145,   145,   147,   145,
     145,   124,   148,   145,   144,   145,   152,   153,   121,   146,
     121,    65,   122,    65,   104,   152,    65,   146,   144,   145,
     436,    65,   145,   144,   145,   148,    65,    65,    65,   152,
     153,   144,   121,   122,   123,   146,   810,    65,   812,   128,
     129,   130,   131,   132,   133,   134,   135,   136,   144,   145,
      65,   144,   826,   144,   143,   474,   145,    30,   474,   474,
     474,    65,   146,    30,    65,    38,   474,   474,   152,   144,
     144,    38,   124,   474,   474,   474,   899,   900,   474,   474,
     696,   474,   121,   121,   961,   899,   900,     5,   474,   474,
     144,   440,    10,   121,   443,   144,   148,   145,    65,   976,
     152,   153,    65,    65,    65,   144,   144,   144,    43,    27,
      28,    65,   145,    31,    65,    33,   144,   891,   734,   893,
     894,    65,    40,   144,   145,   144,   475,   121,    63,   144,
     559,   560,   561,   562,   563,   564,   565,   566,   567,   568,
     144,   570,   571,   572,   573,   574,    64,   589,   585,   923,
     592,   925,   146,    71,   770,   122,   123,   576,   121,   121,
     934,   128,   129,   130,   131,   132,   133,   134,   135,   136,
     121,   587,    90,   947,    65,    65,   143,   951,   145,    65,
      65,   144,   144,   144,   102,    22,    65,     0,    65,   107,
     144,   965,   121,   144,     7,     8,     9,   121,   547,    30,
     144,   119,   120,   552,    17,    18,    19,    38,    21,    65,
      23,    24,    25,    26,   121,   144,   832,   833,    65,   835,
     994,    34,   146,    65,    37,   121,    13,  1001,   577,   578,
     121,   121,    65,   103,   850,   121,   121,    30,    65,   146,
      65,   590,   121,   121,   593,    38,    59,    60,    61,    65,
     146,   680,    65,   144,   144,   871,    65,    65,   144,   144,
     876,    65,    75,    76,    77,   144,    65,    80,   146,    82,
      83,    84,    85,    86,   121,   144,   144,   145,    91,    92,
     896,   150,    30,    96,    97,    98,    99,   903,   121,   905,
      38,    65,   105,   106,   121,   108,   121,   144,    30,   112,
     144,   145,   115,   116,   117,   118,    38,   124,   125,    22,
     127,   144,   121,   121,   822,   823,    65,   144,     5,   144,
      65,    65,   121,    10,    11,    65,    65,    14,   110,   945,
      65,    65,    65,    38,   144,   144,   144,   149,    85,   110,
      27,    28,    29,   110,    31,   144,    33,   110,    30,   124,
     125,   146,   127,    40,    41,    42,    38,    44,    45,   975,
      47,    48,    49,    50,    51,    52,    53,    54,    55,    56,
      57,    58,   121,    30,   990,    65,   121,    64,    22,    66,
      67,    38,    69,    65,    71,   144,   121,   121,   121,    65,
      65,   820,    79,    22,    65,   144,    22,    22,    22,   144,
      87,    88,    22,    90,    22,    22,    22,    22,    65,   144,
     144,   144,    22,    22,    22,   102,   103,    65,    22,    22,
     107,   122,   123,   124,   125,    22,   127,    22,    22,    22,
      22,    22,   119,   120,    22,    22,    22,    22,    22,   121,
     122,   123,   791,    22,   793,   146,   128,   129,   130,   131,
     132,   133,   134,   135,   136,   122,   123,   124,   125,    22,
     127,   143,   121,   145,    65,   122,   123,    65,   121,   121,
      22,   128,   129,   130,   131,   132,   133,   134,   135,   136,
       7,     8,     9,    65,    65,    65,   143,    65,   145,   145,
      17,    18,    19,    65,    21,   146,    23,    24,    25,    26,
     121,   122,   123,   124,   125,   146,   127,    34,    65,   145,
      37,   122,   123,   124,   125,   145,   127,   122,   123,   124,
     125,    22,   127,   145,    22,   146,   122,   123,   124,   125,
       3,   127,    59,    60,    61,   146,    22,    22,    65,   145,
      22,   146,    15,    16,    22,   145,   145,    22,    75,    76,
      77,    22,    22,    80,    22,    82,    83,    84,    85,    86,
      22,    22,    35,    36,    91,    92,   110,    40,   145,    96,
      97,    98,    99,    46,     3,   145,   145,   145,   105,   106,
     145,   108,    65,   145,   144,   112,    15,    16,   115,   116,
     117,   118,   145,   110,   145,    68,   145,    70,   145,    72,
     122,   123,   124,   125,   145,   127,    35,    36,    81,   145,
     145,    40,    85,   122,   123,   124,   125,    46,   127,   145,
      93,    94,    95,    22,   146,     3,    99,   100,   101,   145,
     103,   145,   145,    65,   110,   144,   109,    15,    16,    68,
      65,    70,    65,    72,   110,   122,   123,   124,   125,    78,
     127,   144,    81,   145,    65,    22,    85,    35,    36,    38,
     144,    65,    40,   144,    93,    94,    95,    65,    46,   146,
      99,   100,   101,    30,   103,   122,   123,   124,   125,   144,
     127,   122,   123,   124,   125,    65,   127,    38,    38,   149,
      68,    38,    70,    38,    72,   122,   123,   124,   125,   146,
     127,    38,    30,    81,    38,   146,    38,    85,   122,   123,
     124,   125,    30,   127,    38,    93,    94,    95,    38,   146,
      30,    99,   100,   101,    38,   103,   122,   123,   124,   125,
      38,   127,   146,   122,   123,   124,   125,   145,   127,   122,
     123,   124,   125,    65,   127,    38,   122,   123,   124,   125,
     146,   127,    38,   122,   123,   124,   125,   146,   127,   122,
     123,   124,   125,   146,   127,    65,   122,   123,   124,   125,
     146,   127,    22,   122,   123,   124,   125,   146,   127,   122,
     123,   124,   125,   146,   127,    65,   122,   123,   124,   125,
     146,   127,    65,   122,   123,   124,   125,   146,   127,   122,
     123,   124,   125,   146,   127,    65,   122,   123,   124,   125,
     146,   127,   127,   122,   123,   124,   125,   146,   127,   122,
     123,   124,   125,   146,   127,   144,   122,   123,   124,   125,
     146,   127,   144,   122,   123,   124,   125,   146,   127,   122,
     123,   124,   125,   146,   127,    38,   122,   123,   124,   125,
     146,   127,    38,   122,   123,   124,   125,   146,   127,   122,
     123,   124,   125,   146,   127,    38,   122,   123,   124,   125,
     146,   127,   121,   121,   121,   144,    38,   122,   123,   124,
     125,   144,   127,    38,   122,   123,   124,   125,   144,   127,
     122,   123,   124,   125,   121,   127,   121,   121,    38,   144,
     122,   123,   124,   125,    38,   127,   144,    38,   122,   123,
     124,   125,   144,   127,   122,   123,   124,   125,    38,   127,
      38,    65,   144,   122,   123,   124,   125,   110,   127,   127,
     144,   144,   122,   123,   124,   125,   144,   127,   122,   123,
     124,   125,   144,   127,   144,   144,   122,   123,   124,   125,
      30,   127,   110,   110,   144,    38,   122,   123,   124,   125,
     144,   127,   122,   123,   124,   125,   145,   127,   144,   110,
      38,   146,   146,   151,   146,   146,    65,    65,   144,   148,
     144,    65,   121,   121,   144,   121,    65,    38,   144,    38,
      22,    65,   146,    65,   127,    65,   146,   144,   146,    22,
      65,    65,   148,    22,   121,   121,   121,   144,    22,    22,
      38,    38,   144,   121,   121,   121,   151,   144,   146,   148,
      22,   121,   121,   144,   148,   144,    38,    38,   114,   121,
     121,    38,   146,   144,    65,   148,   121,    38,    41,   151,
     121,   103,   270,   121,   765,   127,   436,   587,   585,   976,
     256,   267,   262,   474,   125,   576,   884,   546,   116,   152,
     505,   307,   824,   592,   105,   621,   467
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,     9,    17,    18,    19,    21,    23,    24,
      25,    26,    34,    37,    59,    60,    61,    65,    75,    76,
      77,    80,    82,    83,    84,    85,    86,    91,    92,    96,
      97,    98,    99,   105,   106,   108,   112,   115,   116,   117,
     118,   155,   156,   157,   158,   159,   160,   161,   162,   163,
     168,   169,   172,   173,   176,   179,   188,   189,   194,   197,
     200,   203,   206,   212,   218,   221,   226,   229,   232,   235,
     236,   239,   240,   242,   243,   244,   247,   248,   252,   253,
     254,   255,   261,   264,   144,   145,   144,   144,   145,    38,
      65,   145,    65,   145,   144,    65,   144,   145,   209,   144,
     144,   144,   144,   144,   144,   145,   145,   144,    22,   144,
     144,   145,   209,   209,   144,   144,   145,   209,   209,    65,
     167,    22,    38,   262,   145,   145,   209,   144,    22,   144,
     145,   144,   145,   144,   145,   209,   144,   145,   209,   209,
      65,   164,    65,   165,    65,   166,   209,     0,   156,    13,
       4,    65,   249,   250,   103,   201,   202,   268,    38,   144,
     144,   152,    65,   144,   152,    65,    65,   174,   175,    22,
       5,    10,    11,    14,    27,    28,    29,    31,    33,    40,
      41,    42,    44,    45,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    56,    57,    58,    64,    66,    67,    69,
      71,    79,    87,    88,    90,   102,   107,   119,   120,   227,
     228,   268,   278,   285,   286,   287,   288,   289,   290,   291,
     292,   293,   294,   295,   296,   297,   298,   299,   300,   301,
     302,   303,   304,   305,   306,   307,   308,   309,   310,   311,
     312,   313,   315,   316,   321,   322,   323,   324,   325,   326,
      65,   121,   144,    12,    65,   104,   213,   214,   215,    12,
      65,   104,   222,   223,    12,    65,   104,   219,   220,    65,
     177,   178,   174,   180,    46,   113,   270,    62,    73,    89,
     256,   257,   318,   319,   320,    12,   115,   190,   191,    30,
      38,    65,   122,   123,   128,   129,   130,   131,   132,   133,
     134,   135,   136,   143,   145,   170,    65,   237,   238,     3,
      15,    16,    35,    36,    68,    70,    72,    78,    81,    85,
      93,    94,    95,    99,   100,   101,   208,   245,   246,   267,
     268,   269,   270,   271,   272,   273,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   284,   314,   317,   144,
     144,    65,   241,   246,   144,   144,   110,    65,   121,   144,
      38,   144,    30,    38,    65,   122,   123,   128,   129,   130,
     131,   132,   133,   134,   135,   136,   143,   145,   185,   187,
     230,   231,   278,   292,   293,   302,   308,   309,   310,   311,
     312,   313,   321,   322,   323,   230,   144,   190,   149,   204,
     205,   281,   198,   199,   268,   207,   208,   144,   109,   208,
     265,   266,   327,   144,   144,   110,    65,   121,   144,   110,
      65,   121,   144,   110,    65,   121,   144,   144,   146,    65,
     121,   145,   251,    20,   250,    22,   121,   146,   144,    65,
     146,   152,    65,   146,   152,    22,    20,   175,    65,    22,
      22,    22,    22,    22,    22,    22,    22,    22,    22,    22,
      22,    22,    22,    22,    22,    22,    22,    22,    22,    22,
      22,    22,    22,    22,   121,   146,    22,    65,    65,    65,
      20,   214,   121,    65,   121,    65,    20,   223,    65,   121,
      65,    20,   220,   145,    20,   178,    20,   147,   183,   184,
     185,   186,   146,   146,    22,   121,   146,    65,    65,    20,
     191,   145,   170,   170,   145,   145,   145,   145,   145,   145,
     145,   145,   145,   145,   170,   122,   123,   124,   125,   127,
     144,   145,    20,   238,    22,    22,    22,    22,    22,    22,
      22,    22,    22,    22,    22,    22,   121,   146,   121,   170,
      20,    65,   146,   110,    65,   144,   145,   185,   185,   145,
     145,   145,   145,   145,   145,   145,   145,   145,   145,   185,
     122,   123,   124,   125,   127,   263,   121,   146,   146,    20,
      30,    38,   145,   195,   196,   121,   146,   121,   146,   121,
     146,    22,   121,   146,   110,    65,   110,    65,   110,    65,
     144,   145,    65,    30,    38,    22,    38,   202,   144,   144,
     144,   209,    65,   144,   144,   209,    65,   170,    30,    65,
     149,   330,   331,    38,    38,    38,    38,    38,    38,    38,
      30,    30,    38,    30,    38,    38,    38,    38,    65,    38,
     331,   145,    38,    38,    30,    65,    65,   328,   228,   144,
     209,    65,    22,   121,     6,    30,    32,    38,    39,    65,
      74,   111,   122,   216,   224,   225,   121,   225,   121,   121,
     225,   121,    38,   122,   123,   210,    65,    20,   184,   186,
      22,   144,   144,   144,    43,    63,   257,    65,   124,   148,
     152,   153,   258,   259,   260,   121,    22,   121,   144,   170,
     171,   170,   170,   170,   170,   170,   170,   170,   170,   170,
     170,   146,   170,   170,   170,   170,   170,   170,    38,    38,
      38,    38,    38,   225,    38,    38,    30,    38,    38,    38,
      38,   245,   144,   209,    65,   144,   121,   170,   144,   209,
     110,   210,   185,   185,   185,   185,   185,   185,   185,   185,
     185,   185,   146,   185,   185,   185,   185,   185,   144,   231,
     144,   209,   144,   209,   170,   144,   150,    30,    38,   121,
     145,   205,   144,   199,   144,   208,   144,   209,    30,   266,
     144,   209,   110,   110,   110,    38,   251,   146,   146,   170,
     144,   146,   144,   146,   144,    38,   329,   150,   329,   151,
     233,   234,   148,   144,    65,    65,    30,    38,   144,   121,
     121,    65,   121,   225,    65,   144,   225,    38,    38,   146,
      22,   185,   182,   181,   121,   144,   145,   260,    65,   170,
      65,    85,   104,   121,   146,   121,   146,   146,   146,   146,
     146,   146,   146,   146,   146,   146,   146,   146,   144,   170,
      65,   144,   144,   146,   146,   146,   146,   146,   146,   146,
     146,   146,   146,   146,   144,   144,   146,   196,   144,    30,
      38,   145,   170,   144,   144,   146,    22,   144,   144,   209,
     144,   209,   148,    65,   121,   146,    65,   217,   225,   224,
     225,   121,   225,   121,   121,   144,    22,   185,   144,   183,
     183,   259,   225,    22,   144,    22,    38,   192,   170,   170,
     170,   144,   144,   170,   170,   146,   251,   170,   144,   144,
      38,   151,   233,   121,   121,   121,   225,   144,   225,   225,
     170,   144,    20,    20,   145,   146,   170,   170,   148,    38,
     121,   144,   144,   144,   146,    22,   144,   121,   225,   217,
     225,   121,   144,   144,   144,   225,   144,   144,    38,   148,
      38,   114,   170,   151,   225,   121,   121,   225,   146,    38,
     148,    30,    65,   122,   123,   145,   193,   210,   211,   144,
      65,   225,   224,   144,    38,    30,    30,   170,    65,   144,
     145,   210,   211,   151,   121,   121,   146,   170,   225,   217,
     146,   121,   225
  };

#if YYDEBUG
  /* TOKEN_NUMBER_[YYLEX-NUM] -- Internal symbol number corresponding
     to YYLEX-NUM.  */
  const unsigned short int
  parser::yytoken_number_[] =
  {
         0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,    59,    40,    41,    35,    58,    91,
      93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   154,   155,   155,   156,   156,   156,   156,   156,   156,
     156,   156,   156,   156,   156,   156,   156,   156,   156,   156,
     156,   156,   156,   156,   156,   156,   156,   156,   156,   156,
     156,   156,   156,   156,   156,   156,   156,   156,   156,   156,
     156,   157,   157,   157,   157,   158,   158,   159,   160,   161,
     162,   163,   164,   164,   164,   164,   164,   164,   165,   165,
     165,   165,   165,   165,   166,   166,   166,   166,   166,   166,
     167,   167,   167,   167,   167,   167,   168,   168,   169,   170,
     170,   170,   170,   170,   170,   170,   170,   170,   170,   170,
     170,   170,   170,   170,   170,   170,   170,   170,   170,   170,
     170,   170,   171,   171,   172,   173,   174,   174,   175,   176,
     177,   177,   178,   180,   179,   181,   179,   182,   179,   183,
     183,   183,   183,   184,   184,   185,   185,   185,   185,   185,
     185,   185,   185,   185,   185,   185,   185,   185,   185,   185,
     185,   185,   185,   185,   185,   185,   186,   187,   187,   188,
     189,   190,   190,   191,   191,   191,   191,   191,   192,   192,
     192,   192,   192,   192,   193,   193,   193,   193,   193,   193,
     193,   193,   194,   195,   195,   196,   196,   196,   196,   196,
     196,   196,   196,   196,   197,   197,   198,   198,   199,   200,
     200,   201,   201,   202,   203,   203,   204,   204,   205,   206,
     206,   206,   206,   207,   207,   208,   208,   208,   208,   208,
     208,   208,   208,   208,   208,   208,   208,   208,   208,   208,
     208,   208,   208,   208,   209,   209,   209,   209,   209,   209,
     210,   210,   210,   211,   211,   211,   212,   213,   213,   214,
     215,   215,   215,   216,   216,   216,   216,   216,   217,   217,
     217,   217,   218,   219,   219,   220,   220,   220,   221,   222,
     222,   223,   223,   223,   224,   224,   224,   224,   224,   225,
     225,   225,   225,   225,   225,   226,   226,   226,   226,   227,
     227,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     229,   229,   230,   230,   231,   231,   231,   231,   231,   231,
     231,   231,   231,   231,   231,   231,   231,   232,   232,   233,
     233,   234,   234,   235,   236,   237,   237,   238,   239,   240,
     241,   241,   241,   241,   242,   243,   243,   243,   243,   244,
     244,   244,   244,   245,   245,   246,   246,   247,   248,   249,
     249,   250,   250,   250,   251,   251,   251,   252,   252,   253,
     253,   253,   253,   253,   253,   254,   254,   254,   254,   254,
     254,   255,   256,   256,   257,   257,   257,   258,   258,   258,
     258,   259,   259,   260,   260,   260,   260,   260,   262,   263,
     261,   264,   264,   264,   264,   265,   265,   266,   266,   267,
     268,   269,   270,   271,   272,   273,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   284,   284,   285,   286,
     286,   287,   288,   289,   290,   291,   292,   292,   293,   294,
     295,   296,   297,   298,   298,   299,   300,   301,   302,   303,
     304,   305,   306,   307,   308,   309,   310,   311,   312,   313,
     314,   315,   316,   317,   318,   318,   319,   320,   321,   322,
     323,   324,   325,   326,   327,   328,   329,   329,   330,   330,
     331
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  parser::yyr2_[] =
  {
         0,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     4,     3,     3,     3,
       3,     3,     2,     3,     1,     3,     4,     2,     2,     3,
       1,     3,     4,     2,     2,     3,     1,     3,     4,     2,
       2,     3,     1,     3,     4,     2,     3,     4,     4,     3,
       1,     1,     1,     3,     3,     3,     3,     3,     2,     2,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     3,     3,     4,     4,     2,     1,     4,     4,
       2,     1,     7,     0,     5,     0,     8,     0,     8,     2,
       2,     1,     1,     4,     2,     3,     1,     1,     1,     3,
       3,     3,     3,     3,     2,     2,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     5,     1,     4,     4,
       4,     2,     1,     9,     6,     5,     7,     7,     2,     4,
       3,     5,     3,     1,     2,     2,     2,     1,     1,     1,
       4,     3,     6,     3,     1,     5,     3,     3,     4,     2,
       2,     3,     1,     1,     2,     5,     3,     1,     1,     2,
       5,     3,     1,     1,     2,     5,     3,     1,     1,     2,
       5,     3,     6,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     2,     4,     3,     5,     1,     3,
       2,     2,     1,     2,     2,     1,     4,     2,     1,     4,
       2,     1,     4,     3,     5,     9,     1,     5,     3,     5,
       7,     9,     4,     2,     1,     5,     7,     4,     4,     2,
       1,     7,     9,     6,     1,     1,     1,     1,     1,     0,
       1,     1,     1,     2,     2,     2,     5,     3,     6,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       5,     6,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     5,     6,     7,
       5,     1,     3,     3,     4,     2,     1,     5,     3,     4,
       4,     6,     3,     5,     3,     2,     5,     3,     6,     2,
       5,     3,     6,     1,     1,     1,     3,     3,     4,     2,
       1,     5,     7,     9,     0,     3,     3,     2,     5,     5,
       6,     3,     7,     8,     5,     5,     6,     3,     7,     8,
       5,     6,     3,     1,     1,     1,     1,     1,     3,     4,
       6,     1,     2,     1,     1,     1,     1,     1,     0,     0,
       5,     2,     5,     3,     6,     3,     1,     1,     1,     3,
       3,     3,     1,     3,     3,     3,     3,     1,     1,     1,
       3,     3,     3,     3,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     1,     3,     3,
       3,     3,     5,     3,     3,     3,     3,     1,     3,     3,
       3,     1,     1,     1,     1,     1,     3,     1,     1,     1,
       1,     3,     3,     3,     3,     3,     1,     1,     3,     3,
       3,     1,     1,     1,     3,     3,     1,     3,     2,     2,
       2
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const parser::yytname_[] =
  {
    "$end", "error", "$undefined", "AR", "AUTOCORR", "BAYESIAN_IRF",
  "BETA_PDF", "CALIB", "CALIB_VAR", "CHECK", "CONF_SIG", "CONSTANT",
  "CORR", "COVAR", "DATAFILE", "DR_ALGO", "DROP", "DSAMPLE", "DYNASAVE",
  "DYNATYPE", "END", "ENDVAL", "EQUAL", "ESTIMATION", "ESTIMATED_PARAMS",
  "ESTIMATED_PARAMS_BOUNDS", "ESTIMATED_PARAMS_INIT", "FILTER_STEP_AHEAD",
  "FILTERED_VARS", "FIRST_OBS", "FLOAT_NUMBER", "FORECAST", "GAMMA_PDF",
  "GRAPH", "HISTVAL", "HP_FILTER", "HP_NGRID", "INITVAL", "INT_NUMBER",
  "INV_GAMMA_PDF", "IRF", "KALMAN_ALGO", "KALMAN_TOL", "LAPLACE",
  "LIK_ALGO", "LIK_INIT", "LINEAR", "LOAD_MH_FILE", "LOGLINEAR", "MH_DROP",
  "MH_INIT_SCALE", "MH_JSCALE", "MH_MODE", "MH_NBLOCKS", "MH_REPLIC",
  "MH_RECOVER", "MODE_CHECK", "MODE_COMPUTE", "MODE_FILE", "MODEL",
  "MODEL_COMPARISON", "MSHOCKS", "MODEL_COMPARISON_APPROXIMATION",
  "MODIFIEDHARMONICMEAN", "MOMENTS_VARENDO", "NAME", "NOBS", "NOCONSTANT",
  "NOCORR", "NODIAGNOSTIC", "NOFUNCTIONS", "NOGRAPH", "NOMOMENTS",
  "NOPRINT", "NORMAL_PDF", "OBSERVATION_TRENDS", "OLR", "OLR_INST",
  "OLR_BETA", "OPTIM", "OPTIM_WEIGHTS", "ORDER", "OSR", "OSR_PARAMS",
  "PARAMETERS", "PERIODS", "PLANNER_OBJECTIVE", "PREFILTER", "PRESAMPLE",
  "PRINT", "PRIOR_TRUNC", "PRIOR_ANALYSIS", "POSTERIOR_ANALYSIS",
  "QZ_CRITERIUM", "RELATIVE_IRF", "REPLIC", "RPLOT", "SHOCKS", "SIGMA_E",
  "SIMUL", "SIMUL_ALGO", "SIMUL_SEED", "SMOOTHER", "SOLVE_ALGO", "STDERR",
  "STEADY", "STOCH_SIMUL", "TEX", "RAMSEY_POLICY", "PLANNER_DISCOUNT",
  "TEX_NAME", "UNIFORM_PDF", "UNIT_ROOT_VARS", "USE_DLL", "VALUES", "VAR",
  "VAREXO", "VAREXO_DET", "VAROBS", "XLS_SHEET", "XLS_RANGE", "COMMA",
  "MINUS", "PLUS", "DIVIDE", "TIMES", "UMINUS", "POWER", "EXP", "LOG",
  "LOG10", "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN", "SINH", "COSH",
  "TANH", "ASINH", "ACOSH", "ATANH", "SQRT", "';'", "'('", "')'", "'#'",
  "':'", "'['", "']'", "'''", "'.'", "'\\\\'", "$accept", "statement_list",
  "statement", "declaration", "dsample", "rplot", "var", "varexo",
  "varexo_det", "parameters", "var_list", "varexo_list", "varexo_det_list",
  "parameter_list", "periods", "init_param", "expression",
  "comma_expression", "initval", "endval", "initval_list", "initval_elem",
  "histval", "histval_list", "histval_elem", "model", "@1", "@2", "@3",
  "equation_list", "equation", "hand_side", "pound_expression",
  "model_var", "shocks", "mshocks", "shock_list", "shock_elem",
  "period_list", "value_list", "sigma_e", "triangular_matrix",
  "triangular_row", "steady", "steady_options_list", "steady_options",
  "check", "check_options_list", "check_options", "simul",
  "simul_options_list", "simul_options", "stoch_simul",
  "stoch_simul_options_list", "stoch_simul_options", "tmp_var_list",
  "signed_integer", "signed_float", "estimated_params", "estimated_list",
  "estimated_elem", "estimated_elem1", "estimated_elem2",
  "estimated_elem3", "estimated_params_init", "estimated_init_list",
  "estimated_init_elem", "estimated_params_bounds",
  "estimated_bounds_list", "estimated_bounds_elem", "prior", "value",
  "estimation", "estimation_options_list", "estimation_options",
  "prior_analysis", "prior_posterior_options_list",
  "prior_posterior_options", "posterior_analysis", "list_optim_option",
  "optim_options", "varobs", "observation_trends", "trend_list",
  "trend_element", "unit_root_vars", "optim_weights", "optim_weights_list",
  "osr_params", "osr", "olr", "olr_option", "olr_options", "olr_inst",
  "calib_var", "calib_var_list", "calib_arg1", "calib_arg2", "calib",
  "dynatype", "dynasave", "model_comparison", "model_comparison_options",
  "model_comparison_option", "filename_list", "filename", "filename_elem",
  "planner_objective", "@4", "@5", "ramsey_policy",
  "ramsey_policy_options_list", "ramsey_policy_options", "o_dr_algo",
  "o_solve_algo", "o_simul_algo", "o_linear", "o_order", "o_replic",
  "o_drop", "o_ar", "o_nocorr", "o_nofunctions", "o_nomoments", "o_irf",
  "o_hp_filter", "o_hp_ngrid", "o_periods", "o_simul", "o_simul_seed",
  "o_qz_criterium", "o_datafile", "o_nobs", "o_first_obs", "o_prefilter",
  "o_presample", "o_lik_algo", "o_lik_init", "o_nograph", "o_conf_sig",
  "o_mh_replic", "o_mh_drop", "o_mh_jscale", "o_optim", "o_mh_init_scale",
  "o_mode_file", "o_mode_compute", "o_mode_check", "o_prior_trunc",
  "o_mh_mode", "o_mh_nblcks", "o_load_mh_file", "o_loglinear",
  "o_nodiagnostic", "o_bayesian_irf", "o_tex", "o_forecast", "o_smoother",
  "o_moments_varendo", "o_filtered_vars", "o_relative_irf",
  "o_kalman_algo", "o_kalman_tol", "o_olr_beta",
  "o_model_comparison_approximation", "o_print", "o_noprint",
  "o_xls_sheet", "o_xls_range", "o_filter_step_ahead", "o_constant",
  "o_noconstant", "o_mh_recover", "o_planner_discount", "range",
  "vec_int_elem", "vec_int_1", "vec_int", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const parser::rhs_number_type
  parser::yyrhs_[] =
  {
       155,     0,    -1,   156,    -1,   155,   156,    -1,   157,    -1,
     168,    -1,   179,    -1,   172,    -1,   173,    -1,   176,    -1,
     169,    -1,   188,    -1,   189,    -1,   194,    -1,   197,    -1,
     200,    -1,   203,    -1,   206,    -1,   226,    -1,   229,    -1,
     232,    -1,   212,    -1,   221,    -1,   218,    -1,   235,    -1,
     236,    -1,   239,    -1,   158,    -1,   159,    -1,   240,    -1,
     242,    -1,   243,    -1,   248,    -1,   252,    -1,   253,    -1,
     254,    -1,   244,    -1,   247,    -1,   255,    -1,   261,    -1,
     264,    -1,   163,    -1,   160,    -1,   161,    -1,   162,    -1,
      17,    38,   144,    -1,    17,    38,    38,   144,    -1,    96,
     209,   144,    -1,   115,   164,   144,    -1,   116,   165,   144,
      -1,   117,   166,   144,    -1,    84,   167,   144,    -1,   164,
      65,    -1,   164,   121,    65,    -1,    65,    -1,   164,    65,
     110,    -1,   164,   121,    65,   110,    -1,    65,   110,    -1,
     165,    65,    -1,   165,   121,    65,    -1,    65,    -1,   165,
      65,   110,    -1,   165,   121,    65,   110,    -1,    65,   110,
      -1,   166,    65,    -1,   166,   121,    65,    -1,    65,    -1,
     166,    65,   110,    -1,   166,   121,    65,   110,    -1,    65,
     110,    -1,   167,    65,    -1,   167,   121,    65,    -1,    65,
      -1,   167,    65,   110,    -1,   167,   121,    65,   110,    -1,
      65,   110,    -1,    85,    38,   144,    -1,    85,    22,    38,
     144,    -1,    65,    22,   170,   144,    -1,   145,   170,   146,
      -1,    65,    -1,    30,    -1,    38,    -1,   170,   123,   170,
      -1,   170,   122,   170,    -1,   170,   124,   170,    -1,   170,
     125,   170,    -1,   170,   127,   170,    -1,   122,   170,    -1,
     123,   170,    -1,   128,   145,   170,   146,    -1,   129,   145,
     170,   146,    -1,   130,   145,   170,   146,    -1,   131,   145,
     170,   146,    -1,   132,   145,   170,   146,    -1,   133,   145,
     170,   146,    -1,   134,   145,   170,   146,    -1,   135,   145,
     170,   146,    -1,   136,   145,   170,   146,    -1,   143,   145,
     170,   146,    -1,    65,   145,   170,   146,    -1,    65,   145,
     171,   146,    -1,   170,   121,   170,    -1,   171,   121,   170,
      -1,    37,   144,   174,    20,    -1,    21,   144,   174,    20,
      -1,   174,   175,    -1,   175,    -1,    65,    22,   170,   144,
      -1,    34,   144,   177,    20,    -1,   177,   178,    -1,   178,
      -1,    65,   145,   210,   146,    22,   170,   144,    -1,    -1,
      59,   144,   180,   183,    20,    -1,    -1,    59,   145,   270,
     146,   144,   181,   183,    20,    -1,    -1,    59,   145,   113,
     146,   144,   182,   183,    20,    -1,   183,   184,    -1,   183,
     186,    -1,   184,    -1,   186,    -1,   185,    22,   185,   144,
      -1,   185,   144,    -1,   145,   185,   146,    -1,   187,    -1,
      30,    -1,    38,    -1,   185,   123,   185,    -1,   185,   122,
     185,    -1,   185,   124,   185,    -1,   185,   125,   185,    -1,
     185,   127,   185,    -1,   122,   185,    -1,   123,   185,    -1,
     128,   145,   185,   146,    -1,   129,   145,   185,   146,    -1,
     130,   145,   185,   146,    -1,   131,   145,   185,   146,    -1,
     132,   145,   185,   146,    -1,   133,   145,   185,   146,    -1,
     134,   145,   185,   146,    -1,   135,   145,   185,   146,    -1,
     136,   145,   185,   146,    -1,   143,   145,   185,   146,    -1,
     147,    65,    22,   185,   144,    -1,    65,    -1,    65,   145,
     210,   146,    -1,    97,   144,   190,    20,    -1,    61,   144,
     190,    20,    -1,   190,   191,    -1,   191,    -1,   115,    65,
     144,    85,   192,   144,   114,   193,   144,    -1,   115,    65,
     144,   104,   170,   144,    -1,   115,    65,    22,   170,   144,
      -1,   115,    65,   121,    65,    22,   170,   144,    -1,    12,
      65,   121,    65,    22,   170,   144,    -1,   192,    38,    -1,
     192,    38,   148,    38,    -1,   192,   121,    38,    -1,   192,
     121,    38,   148,    38,    -1,    38,   148,    38,    -1,    38,
      -1,   193,   211,    -1,   193,   210,    -1,   193,    65,    -1,
     211,    -1,   210,    -1,    65,    -1,   193,   145,   170,   146,
      -1,   145,   170,   146,    -1,    98,    22,   149,   195,   150,
     144,    -1,   195,   144,   196,    -1,   196,    -1,   196,   121,
     145,   170,   146,    -1,   196,   121,    30,    -1,   196,   121,
      38,    -1,   196,   145,   170,   146,    -1,   196,    30,    -1,
     196,    38,    -1,   145,   170,   146,    -1,    30,    -1,    38,
      -1,   105,   144,    -1,   105,   145,   198,   146,   144,    -1,
     198,   121,   199,    -1,   199,    -1,   268,    -1,     9,   144,
      -1,     9,   145,   201,   146,   144,    -1,   201,   121,   202,
      -1,   202,    -1,   268,    -1,    99,   144,    -1,    99,   145,
     204,   146,   144,    -1,   204,   121,   205,    -1,   205,    -1,
     281,    -1,   106,   144,    -1,   106,   145,   207,   146,   144,
      -1,   106,   209,   144,    -1,   106,   145,   207,   146,   209,
     144,    -1,   207,   121,   208,    -1,   208,    -1,   267,    -1,
     268,    -1,   269,    -1,   270,    -1,   271,    -1,   272,    -1,
     273,    -1,   274,    -1,   275,    -1,   276,    -1,   277,    -1,
     278,    -1,   314,    -1,   279,    -1,   280,    -1,   281,    -1,
     282,    -1,   283,    -1,   284,    -1,   209,    65,    -1,   209,
      65,    22,    65,    -1,   209,   121,    65,    -1,   209,   121,
      65,    22,    65,    -1,    65,    -1,    65,    22,    65,    -1,
     123,    38,    -1,   122,    38,    -1,    38,    -1,   123,    30,
      -1,   122,    30,    -1,    30,    -1,    24,   144,   213,    20,
      -1,   213,   214,    -1,   214,    -1,   215,   121,   216,   144,
      -1,   104,    65,    -1,    65,    -1,    12,    65,   121,    65,
      -1,   224,   121,   217,    -1,   225,   121,   224,   121,   217,
      -1,   225,   121,   225,   121,   225,   121,   224,   121,   217,
      -1,   225,    -1,   225,   121,   225,   121,   225,    -1,   225,
     121,   225,    -1,   225,   121,   225,   121,   225,    -1,   225,
     121,   225,   121,   225,   121,   225,    -1,   225,   121,   225,
     121,   225,   121,   225,   121,   225,    -1,    26,   144,   219,
      20,    -1,   219,   220,    -1,   220,    -1,   104,    65,   121,
     225,   144,    -1,    12,    65,   121,    65,   121,   225,   144,
      -1,    65,   121,   225,   144,    -1,    25,   144,   222,    20,
      -1,   222,   223,    -1,   223,    -1,   104,    65,   121,   225,
     121,   225,   144,    -1,    12,    65,   121,    65,   121,   225,
     121,   225,   144,    -1,    65,   121,   225,   121,   225,   144,
      -1,     6,    -1,    32,    -1,    74,    -1,    39,    -1,   111,
      -1,    -1,    38,    -1,    30,    -1,    65,    -1,   122,    38,
      -1,   122,    30,    -1,    23,   144,    -1,    23,   145,   227,
     146,   144,    -1,    23,   209,   144,    -1,    23,   145,   227,
     146,   209,   144,    -1,   227,   121,   228,    -1,   228,    -1,
     285,    -1,   286,    -1,   287,    -1,   288,    -1,   289,    -1,
     290,    -1,   291,    -1,   292,    -1,   293,    -1,   294,    -1,
     295,    -1,   296,    -1,   297,    -1,   298,    -1,   299,    -1,
     300,    -1,   301,    -1,   302,    -1,   303,    -1,   304,    -1,
     305,    -1,   306,    -1,   307,    -1,   308,    -1,   278,    -1,
     309,    -1,   310,    -1,   311,    -1,   312,    -1,   313,    -1,
     315,    -1,   316,    -1,   321,    -1,   322,    -1,   323,    -1,
     268,    -1,   324,    -1,   325,    -1,   326,    -1,    91,   145,
     230,   146,   144,    -1,    91,   145,   230,   146,   209,   144,
      -1,   230,   121,   231,    -1,   231,    -1,   292,    -1,   293,
      -1,   302,    -1,   308,    -1,   278,    -1,   309,    -1,   310,
      -1,   311,    -1,   312,    -1,   313,    -1,   321,    -1,   322,
      -1,   323,    -1,    92,   145,   230,   146,   144,    -1,    92,
     145,   230,   146,   209,   144,    -1,   151,    65,   151,   121,
     151,    65,   151,    -1,   151,    65,   151,   121,   225,    -1,
     233,    -1,   234,   121,   233,    -1,   118,   209,   144,    -1,
      75,   144,   237,    20,    -1,   237,   238,    -1,   238,    -1,
      65,   145,   170,   146,   144,    -1,   112,   209,   144,    -1,
      80,   144,   241,    20,    -1,   241,    65,   170,   144,    -1,
     241,    65,   121,    65,   170,   144,    -1,    65,   170,   144,
      -1,    65,   121,    65,   170,   144,    -1,    83,   209,   144,
      -1,    82,   144,    -1,    82,   145,   246,   146,   144,    -1,
      82,   209,   144,    -1,    82,   145,   246,   146,   209,   144,
      -1,    76,   144,    -1,    76,   145,   246,   146,   144,    -1,
      76,   209,   144,    -1,    76,   145,   246,   146,   209,   144,
      -1,   317,    -1,   208,    -1,   245,    -1,   246,   121,   245,
      -1,    77,   209,   144,    -1,     8,   144,   249,    20,    -1,
     249,   250,    -1,   250,    -1,    65,   251,    22,   170,   144,
      -1,    65,   121,    65,   251,    22,   170,   144,    -1,     4,
      65,   145,    38,   146,   251,    22,   170,   144,    -1,    -1,
     145,    38,   146,    -1,   145,    30,   146,    -1,     7,   144,
      -1,     7,   145,    13,   146,   144,    -1,    19,   145,    65,
     146,   144,    -1,    19,   145,    65,   146,   209,   144,    -1,
      19,    65,   144,    -1,    19,   145,    65,   152,    65,   146,
     144,    -1,    19,   145,    65,   152,    65,   146,   209,   144,
      -1,    19,    65,   152,    65,   144,    -1,    18,   145,    65,
     146,   144,    -1,    18,   145,    65,   146,   209,   144,    -1,
      18,    65,   144,    -1,    18,   145,    65,   152,    65,   146,
     144,    -1,    18,   145,    65,   152,    65,   146,   209,   144,
      -1,    18,    65,   152,    65,   144,    -1,    60,   145,   256,
     146,   258,   144,    -1,   256,   121,   257,    -1,   257,    -1,
     318,    -1,   319,    -1,   320,    -1,   259,    -1,   258,   121,
     259,    -1,   259,   145,   225,   146,    -1,   258,   121,   259,
     145,   225,   146,    -1,   260,    -1,   259,   260,    -1,    65,
      -1,   153,    -1,   124,    -1,   148,    -1,   152,    -1,    -1,
      -1,    86,   262,   185,   263,   144,    -1,   108,   144,    -1,
     108,   145,   265,   146,   144,    -1,   108,   209,   144,    -1,
     108,   145,   265,   146,   209,   144,    -1,   265,   121,   266,
      -1,   266,    -1,   208,    -1,   327,    -1,    15,    22,    38,
      -1,   103,    22,    38,    -1,   100,    22,    38,    -1,    46,
      -1,    81,    22,    38,    -1,    95,    22,    38,    -1,    16,
      22,    38,    -1,     3,    22,    38,    -1,    68,    -1,    70,
      -1,    72,    -1,    40,    22,    38,    -1,    35,    22,    38,
      -1,    36,    22,    38,    -1,    85,    22,    38,    -1,    99,
      -1,   101,    22,    38,    -1,    93,    22,    38,    -1,    93,
      22,    30,    -1,    14,    22,    65,    -1,    66,    22,   331,
      -1,    66,    22,    38,    -1,    29,    22,    38,    -1,    87,
      22,    38,    -1,    88,    22,    38,    -1,    44,    22,    38,
      -1,    45,    22,    38,    -1,    71,    -1,    33,    -1,    10,
      22,    30,    -1,    54,    22,    38,    -1,    49,    22,    30,
      -1,    51,    22,    30,    -1,    79,    22,   145,   234,   146,
      -1,    50,    22,    30,    -1,    50,    22,    38,    -1,    58,
      22,    65,    -1,    57,    22,    38,    -1,    56,    -1,    90,
      22,    30,    -1,    52,    22,    38,    -1,    53,    22,    38,
      -1,    47,    -1,    48,    -1,    69,    -1,     5,    -1,   107,
      -1,    31,    22,    38,    -1,   102,    -1,    64,    -1,    28,
      -1,    94,    -1,    41,    22,    38,    -1,    42,    22,    38,
      -1,    78,    22,   225,    -1,    62,    22,    43,    -1,    62,
      22,    63,    -1,    89,    -1,    73,    -1,   119,    22,    65,
      -1,   120,    22,   328,    -1,    27,    22,   331,    -1,    11,
      -1,    67,    -1,    55,    -1,   109,    22,    30,    -1,    65,
     148,    65,    -1,    38,    -1,    38,   148,    38,    -1,   149,
     329,    -1,   330,   329,    -1,   330,   150,    -1
  };

  /* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
     YYRHS.  */
  const unsigned short int
  parser::yyprhs_[] =
  {
         0,     0,     3,     5,     8,    10,    12,    14,    16,    18,
      20,    22,    24,    26,    28,    30,    32,    34,    36,    38,
      40,    42,    44,    46,    48,    50,    52,    54,    56,    58,
      60,    62,    64,    66,    68,    70,    72,    74,    76,    78,
      80,    82,    84,    86,    88,    90,    94,    99,   103,   107,
     111,   115,   119,   122,   126,   128,   132,   137,   140,   143,
     147,   149,   153,   158,   161,   164,   168,   170,   174,   179,
     182,   185,   189,   191,   195,   200,   203,   207,   212,   217,
     221,   223,   225,   227,   231,   235,   239,   243,   247,   250,
     253,   258,   263,   268,   273,   278,   283,   288,   293,   298,
     303,   308,   313,   317,   321,   326,   331,   334,   336,   341,
     346,   349,   351,   359,   360,   366,   367,   376,   377,   386,
     389,   392,   394,   396,   401,   404,   408,   410,   412,   414,
     418,   422,   426,   430,   434,   437,   440,   445,   450,   455,
     460,   465,   470,   475,   480,   485,   490,   496,   498,   503,
     508,   513,   516,   518,   528,   535,   541,   549,   557,   560,
     565,   569,   575,   579,   581,   584,   587,   590,   592,   594,
     596,   601,   605,   612,   616,   618,   624,   628,   632,   637,
     640,   643,   647,   649,   651,   654,   660,   664,   666,   668,
     671,   677,   681,   683,   685,   688,   694,   698,   700,   702,
     705,   711,   715,   722,   726,   728,   730,   732,   734,   736,
     738,   740,   742,   744,   746,   748,   750,   752,   754,   756,
     758,   760,   762,   764,   766,   769,   774,   778,   784,   786,
     790,   793,   796,   798,   801,   804,   806,   811,   814,   816,
     821,   824,   826,   831,   835,   841,   851,   853,   859,   863,
     869,   877,   887,   892,   895,   897,   903,   911,   916,   921,
     924,   926,   934,   944,   951,   953,   955,   957,   959,   961,
     962,   964,   966,   968,   971,   974,   977,   983,   987,   994,
     998,  1000,  1002,  1004,  1006,  1008,  1010,  1012,  1014,  1016,
    1018,  1020,  1022,  1024,  1026,  1028,  1030,  1032,  1034,  1036,
    1038,  1040,  1042,  1044,  1046,  1048,  1050,  1052,  1054,  1056,
    1058,  1060,  1062,  1064,  1066,  1068,  1070,  1072,  1074,  1076,
    1078,  1084,  1091,  1095,  1097,  1099,  1101,  1103,  1105,  1107,
    1109,  1111,  1113,  1115,  1117,  1119,  1121,  1123,  1129,  1136,
    1144,  1150,  1152,  1156,  1160,  1165,  1168,  1170,  1176,  1180,
    1185,  1190,  1197,  1201,  1207,  1211,  1214,  1220,  1224,  1231,
    1234,  1240,  1244,  1251,  1253,  1255,  1257,  1261,  1265,  1270,
    1273,  1275,  1281,  1289,  1299,  1300,  1304,  1308,  1311,  1317,
    1323,  1330,  1334,  1342,  1351,  1357,  1363,  1370,  1374,  1382,
    1391,  1397,  1404,  1408,  1410,  1412,  1414,  1416,  1418,  1422,
    1427,  1434,  1436,  1439,  1441,  1443,  1445,  1447,  1449,  1450,
    1451,  1457,  1460,  1466,  1470,  1477,  1481,  1483,  1485,  1487,
    1491,  1495,  1499,  1501,  1505,  1509,  1513,  1517,  1519,  1521,
    1523,  1527,  1531,  1535,  1539,  1541,  1545,  1549,  1553,  1557,
    1561,  1565,  1569,  1573,  1577,  1581,  1585,  1587,  1589,  1593,
    1597,  1601,  1605,  1611,  1615,  1619,  1623,  1627,  1629,  1633,
    1637,  1641,  1643,  1645,  1647,  1649,  1651,  1655,  1657,  1659,
    1661,  1663,  1667,  1671,  1675,  1679,  1683,  1685,  1687,  1691,
    1695,  1699,  1701,  1703,  1705,  1709,  1713,  1715,  1719,  1722,
    1725
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    89,    89,    90,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   122,   123,   124,   125,   126,   127,   128,   129,
     130,   135,   136,   137,   138,   142,   143,   146,   149,   153,
     157,   161,   165,   167,   169,   171,   173,   175,   180,   182,
     184,   186,   188,   190,   195,   197,   199,   201,   203,   205,
     210,   212,   214,   216,   218,   220,   225,   229,   237,   242,
     244,   246,   248,   250,   252,   254,   256,   258,   260,   262,
     264,   266,   268,   270,   272,   274,   276,   278,   280,   282,
     284,   286,   291,   293,   297,   302,   307,   308,   312,   317,
     322,   323,   327,   332,   332,   333,   333,   335,   335,   340,
     341,   342,   343,   347,   349,   354,   355,   356,   358,   360,
     362,   364,   366,   368,   370,   372,   374,   376,   378,   380,
     382,   384,   386,   388,   390,   392,   396,   400,   402,   407,
     411,   415,   416,   420,   422,   424,   426,   428,   433,   435,
     437,   439,   441,   443,   449,   451,   453,   455,   457,   459,
     461,   463,   468,   473,   475,   480,   482,   484,   486,   488,
     490,   492,   494,   496,   501,   505,   509,   510,   513,   517,
     519,   523,   524,   527,   531,   533,   537,   538,   541,   545,
     547,   549,   551,   555,   556,   559,   560,   561,   562,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,   581,   583,   585,   587,   589,   591,
     596,   598,   600,   605,   607,   609,   614,   619,   621,   626,
     630,   635,   640,   650,   655,   661,   671,   676,   687,   693,
     701,   711,   725,   729,   731,   735,   742,   751,   760,   764,
     766,   770,   779,   790,   802,   804,   806,   808,   810,   815,
     816,   817,   818,   819,   821,   828,   830,   832,   834,   839,
     840,   843,   844,   845,   846,   847,   848,   849,   850,   851,
     852,   853,   854,   855,   856,   857,   858,   859,   860,   861,
     862,   863,   864,   865,   866,   867,   868,   869,   870,   871,
     872,   873,   874,   875,   876,   877,   878,   879,   880,   881,
     885,   887,   892,   893,   897,   898,   899,   900,   901,   902,
     903,   904,   905,   906,   907,   908,   909,   913,   915,   920,
     921,   925,   926,   930,   935,   940,   941,   944,   948,   951,
     955,   957,   959,   961,   965,   968,   969,   970,   971,   974,
     975,   976,   977,   980,   981,   984,   985,   988,   991,   995,
     996,   999,  1000,  1001,  1004,  1005,  1006,  1009,  1010,  1013,
    1014,  1015,  1016,  1017,  1018,  1020,  1021,  1022,  1023,  1024,
    1025,  1027,  1031,  1032,  1035,  1036,  1037,  1040,  1041,  1042,
    1043,  1046,  1047,  1050,  1051,  1052,  1053,  1054,  1057,  1057,
    1057,  1060,  1062,  1064,  1066,  1071,  1072,  1075,  1076,  1079,
    1080,  1081,  1082,  1083,  1084,  1085,  1086,  1087,  1088,  1089,
    1090,  1091,  1092,  1093,  1094,  1095,  1096,  1097,  1099,  1100,
    1101,  1103,  1104,  1105,  1106,  1107,  1108,  1109,  1110,  1111,
    1112,  1113,  1114,  1115,  1116,  1117,  1118,  1119,  1120,  1121,
    1122,  1123,  1124,  1125,  1126,  1127,  1128,  1129,  1130,  1131,
    1132,  1133,  1134,  1135,  1137,  1139,  1142,  1143,  1144,  1145,
    1146,  1147,  1148,  1149,  1150,  1152,  1160,  1161,  1165,  1166,
    1175
  };

  // Print the state stack on the debug stream.
  void
  parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (state_stack_type::const_iterator i = yystate_stack_.begin ();
	 i != yystate_stack_.end (); ++i)
      *yycdebug_ << ' ' << *i;
    *yycdebug_ << std::endl;
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  parser::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    /* Print the symbols being reduced, and their result.  */
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
	       << " (line " << yylno << "), ";
    /* The symbols being reduced.  */
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
		       yyrhs_[yyprhs_[yyrule] + yyi],
		       &(yysemantic_stack_[(yynrhs) - (yyi + 1)]),
		       &(yylocation_stack_[(yynrhs) - (yyi + 1)]));
  }
#endif // YYDEBUG

  /* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
  parser::token_number_type
  parser::yytranslate_ (int t)
  {
    static
    const token_number_type
    translate_table[] =
    {
           0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,   147,     2,     2,     2,   151,
     145,   146,     2,     2,     2,     2,   152,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   148,   144,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   149,   153,   150,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1396;
  const int parser::yynnts_ = 178;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 147;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 154;

  const unsigned int parser::yyuser_token_number_max_ = 398;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1177 "DynareBison.yy"


void
yy::parser::error(const yy::parser::location_type &l,
                  const string &m)
{
  driver.error(l, m);
}

/*
  Local variables:
  mode: C++
  End:
*/

