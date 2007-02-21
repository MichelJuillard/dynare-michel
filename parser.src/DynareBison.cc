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
  /* Line 547 of yacc.c.  */
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
	  case 46:
#line 143 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 47:
#line 144 "DynareBison.yy"
    {driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 48:
#line 147 "DynareBison.yy"
    {driver.rplot();;}
    break;

  case 53:
#line 167 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 54:
#line 169 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 55:
#line 171 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 56:
#line 173 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 57:
#line 175 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 58:
#line 177 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 59:
#line 182 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 60:
#line 184 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 61:
#line 186 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 62:
#line 188 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 63:
#line 190 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 64:
#line 192 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 65:
#line 197 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 66:
#line 199 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 67:
#line 201 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 68:
#line 203 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 69:
#line 205 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 70:
#line 207 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 71:
#line 212 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 72:
#line 214 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 73:
#line 216 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 74:
#line 218 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 75:
#line 220 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 76:
#line 222 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 77:
#line 227 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 78:
#line 231 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 79:
#line 238 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 80:
#line 242 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 81:
#line 249 "DynareBison.yy"
    {driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 82:
#line 254 "DynareBison.yy"
    { (yyval.exp_val) = (yysemantic_stack_[(3) - (2)].exp_val);;}
    break;

  case 83:
#line 256 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 84:
#line 258 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 85:
#line 260 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 86:
#line 262 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::PLUS);;}
    break;

  case 87:
#line 264 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::MINUS);;}
    break;

  case 88:
#line 266 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::DIVIDE);;}
    break;

  case 89:
#line 268 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::TIMES);;}
    break;

  case 90:
#line 270 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::POWER);;}
    break;

  case 91:
#line 272 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(2) - (2)].exp_val), token::UMINUS);;}
    break;

  case 92:
#line 274 "DynareBison.yy"
    {(yyval.exp_val) = (yysemantic_stack_[(2) - (2)].exp_val);;}
    break;

  case 93:
#line 276 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::EXP);;}
    break;

  case 94:
#line 278 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::LOG);;}
    break;

  case 95:
#line 280 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::LOG10);;}
    break;

  case 96:
#line 282 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::SIN);;}
    break;

  case 97:
#line 284 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::COS);;}
    break;

  case 98:
#line 286 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::TAN);;}
    break;

  case 99:
#line 288 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::ASIN);;}
    break;

  case 100:
#line 290 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::ACOS);;}
    break;

  case 101:
#line 292 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::ATAN);;}
    break;

  case 102:
#line 294 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), token::SQRT);;}
    break;

  case 103:
#line 296 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), (yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 104:
#line 298 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(4) - (3)].exp_val), (yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 105:
#line 303 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::COMMA);;}
    break;

  case 106:
#line 305 "DynareBison.yy"
    {(yyval.exp_val) = driver.add_expression_token((yysemantic_stack_[(3) - (1)].exp_val), (yysemantic_stack_[(3) - (3)].exp_val), token::COMMA);;}
    break;

  case 107:
#line 309 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 108:
#line 311 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 109:
#line 315 "DynareBison.yy"
    {driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 110:
#line 320 "DynareBison.yy"
    {driver.end_endval();;}
    break;

  case 113:
#line 330 "DynareBison.yy"
    {driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 114:
#line 335 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 117:
#line 345 "DynareBison.yy"
    {driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 120:
#line 353 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 121:
#line 354 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 123:
#line 359 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 125:
#line 360 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 127:
#line 362 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 129:
#line 364 "DynareBison.yy"
    { driver.sparse_dll(); driver.begin_model(); ;}
    break;

  case 131:
#line 366 "DynareBison.yy"
    { driver.sparse_dll(); driver.begin_model(); ;}
    break;

  case 137:
#line 379 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].model_val), (yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 138:
#line 381 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].model_val));;}
    break;

  case 139:
#line 385 "DynareBison.yy"
    {(yyval.model_val) = (yysemantic_stack_[(3) - (2)].model_val);;}
    break;

  case 141:
#line 388 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 142:
#line 390 "DynareBison.yy"
    {(yysemantic_stack_[(1) - (1)].string_val)->append(".0"); (yyval.model_val) = driver.add_model_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 143:
#line 392 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_plus((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 144:
#line 394 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_minus((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 145:
#line 396 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_divide((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 146:
#line 398 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_times((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 147:
#line 400 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_power((yysemantic_stack_[(3) - (1)].model_val), (yysemantic_stack_[(3) - (3)].model_val));;}
    break;

  case 148:
#line 402 "DynareBison.yy"
    { (yyval.model_val) = driver.add_model_uminus((yysemantic_stack_[(2) - (2)].model_val));;}
    break;

  case 149:
#line 404 "DynareBison.yy"
    {(yyval.model_val) = (yysemantic_stack_[(2) - (2)].model_val);;}
    break;

  case 150:
#line 406 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_exp((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 151:
#line 408 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_log((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 152:
#line 410 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_log10((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 153:
#line 412 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_sin((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 154:
#line 414 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_cos((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 155:
#line 416 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_tan((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 156:
#line 418 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_asin((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 157:
#line 420 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_acos((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 158:
#line 422 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_atan((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 159:
#line 424 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_sqrt((yysemantic_stack_[(4) - (3)].model_val));;}
    break;

  case 160:
#line 428 "DynareBison.yy"
    {driver.declare_and_init_local_parameter((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].model_val));;}
    break;

  case 161:
#line 432 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 162:
#line 434 "DynareBison.yy"
    {(yyval.model_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 163:
#line 438 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 164:
#line 442 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 167:
#line 452 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 168:
#line 454 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].exp_val));;}
    break;

  case 169:
#line 456 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 170:
#line 458 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 171:
#line 460 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 172:
#line 465 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 173:
#line 467 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 174:
#line 469 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 175:
#line 471 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 176:
#line 473 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 177:
#line 475 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 178:
#line 481 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 179:
#line 483 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 180:
#line 485 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 181:
#line 487 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 182:
#line 489 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 183:
#line 491 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 184:
#line 493 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 185:
#line 495 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (2)].exp_val));;}
    break;

  case 186:
#line 500 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 187:
#line 505 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 188:
#line 507 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 189:
#line 512 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 190:
#line 514 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 191:
#line 516 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 192:
#line 518 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 193:
#line 520 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 194:
#line 522 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 195:
#line 524 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].exp_val));;}
    break;

  case 196:
#line 526 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 197:
#line 528 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 198:
#line 533 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 199:
#line 537 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 203:
#line 549 "DynareBison.yy"
    {driver.check();;}
    break;

  case 204:
#line 551 "DynareBison.yy"
    {driver.check();;}
    break;

  case 208:
#line 563 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 209:
#line 565 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 213:
#line 577 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 214:
#line 579 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 215:
#line 581 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 216:
#line 583 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 239:
#line 614 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 240:
#line 616 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 241:
#line 618 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 242:
#line 620 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 243:
#line 622 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 244:
#line 624 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 245:
#line 629 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 246:
#line 631 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 247:
#line 633 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 248:
#line 638 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 249:
#line 640 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 250:
#line 642 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 251:
#line 647 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 252:
#line 652 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 253:
#line 654 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 255:
#line 663 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 256:
#line 668 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 257:
#line 673 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 258:
#line 683 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 259:
#line 688 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 260:
#line 694 "DynareBison.yy"
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

  case 261:
#line 704 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 262:
#line 709 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 263:
#line 720 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 264:
#line 726 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 265:
#line 734 "DynareBison.yy"
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

  case 266:
#line 744 "DynareBison.yy"
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

  case 267:
#line 758 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 268:
#line 762 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 269:
#line 764 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 270:
#line 768 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 271:
#line 775 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 272:
#line 784 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 273:
#line 793 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 274:
#line 797 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 275:
#line 799 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 276:
#line 803 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 277:
#line 812 "DynareBison.yy"
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

  case 278:
#line 823 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 279:
#line 835 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 280:
#line 837 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 281:
#line 839 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 282:
#line 841 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 283:
#line 843 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 284:
#line 847 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 288:
#line 852 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 289:
#line 854 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 290:
#line 861 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 291:
#line 863 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 292:
#line 865 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 293:
#line 867 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 335:
#line 918 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 336:
#line 920 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 352:
#line 946 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 353:
#line 948 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 354:
#line 952 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 355:
#line 953 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 358:
#line 963 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 359:
#line 968 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 362:
#line 977 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].exp_val));;}
    break;

  case 363:
#line 980 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 364:
#line 984 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 365:
#line 988 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].exp_val));;}
    break;

  case 366:
#line 990 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].exp_val));;}
    break;

  case 367:
#line 992 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].exp_val));;}
    break;

  case 368:
#line 994 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 369:
#line 997 "DynareBison.yy"
    {driver.set_osr_params();;}
    break;

  case 370:
#line 1000 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 371:
#line 1001 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 372:
#line 1002 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 373:
#line 1003 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 374:
#line 1006 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 375:
#line 1007 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 376:
#line 1008 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 377:
#line 1009 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 382:
#line 1020 "DynareBison.yy"
    {driver.set_olr_inst();;}
    break;

  case 383:
#line 1024 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 386:
#line 1031 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].exp_val));;}
    break;

  case 387:
#line 1032 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].exp_val));;}
    break;

  case 388:
#line 1033 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].exp_val));;}
    break;

  case 389:
#line 1036 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 390:
#line 1037 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 391:
#line 1038 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 392:
#line 1041 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 393:
#line 1042 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 394:
#line 1045 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 395:
#line 1046 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 396:
#line 1047 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 397:
#line 1048 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 398:
#line 1049 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 399:
#line 1050 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 400:
#line 1052 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 401:
#line 1053 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 402:
#line 1054 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 403:
#line 1055 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 404:
#line 1056 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 405:
#line 1057 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 406:
#line 1060 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 412:
#line 1072 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 413:
#line 1073 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 414:
#line 1074 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 415:
#line 1075 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 416:
#line 1078 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 417:
#line 1079 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 419:
#line 1083 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 420:
#line 1084 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 421:
#line 1085 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 422:
#line 1086 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 423:
#line 1089 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 424:
#line 1089 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].model_val)); ;}
    break;

  case 426:
#line 1093 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 427:
#line 1095 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 428:
#line 1097 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 429:
#line 1099 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 434:
#line 1111 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 435:
#line 1112 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 436:
#line 1113 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 437:
#line 1114 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 438:
#line 1115 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 439:
#line 1116 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 440:
#line 1117 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 441:
#line 1118 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 442:
#line 1119 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 443:
#line 1120 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 444:
#line 1121 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 445:
#line 1122 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 446:
#line 1123 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 447:
#line 1124 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 448:
#line 1125 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 449:
#line 1126 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 450:
#line 1127 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 451:
#line 1128 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 452:
#line 1129 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 453:
#line 1130 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 454:
#line 1132 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 455:
#line 1133 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 456:
#line 1134 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 457:
#line 1136 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 458:
#line 1137 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 459:
#line 1138 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 460:
#line 1139 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 461:
#line 1140 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 462:
#line 1141 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 463:
#line 1142 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 464:
#line 1143 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 465:
#line 1144 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 466:
#line 1145 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 467:
#line 1146 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 469:
#line 1148 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 470:
#line 1149 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 471:
#line 1150 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 472:
#line 1151 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 473:
#line 1152 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 474:
#line 1153 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 475:
#line 1154 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 476:
#line 1155 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 477:
#line 1156 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 478:
#line 1157 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 479:
#line 1158 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 480:
#line 1159 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 481:
#line 1160 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 482:
#line 1161 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 483:
#line 1162 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 484:
#line 1163 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 485:
#line 1164 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 486:
#line 1165 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 487:
#line 1166 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 488:
#line 1167 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 489:
#line 1168 "DynareBison.yy"
    {driver.option_num("olr_beta", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 490:
#line 1171 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 491:
#line 1173 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 492:
#line 1175 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 493:
#line 1176 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 494:
#line 1177 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1178 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 496:
#line 1179 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 497:
#line 1180 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 498:
#line 1181 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 499:
#line 1182 "DynareBison.yy"
    {driver.option_num("load_mh_file", "-1");;}
    break;

  case 500:
#line 1183 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 501:
#line 1186 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 503:
#line 1195 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 504:
#line 1198 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 505:
#line 1200 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 506:
#line 1208 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2143 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -911;
  const short int
  parser::yypact_[] =
  {
       964,   185,   -38,   411,   268,    35,    -6,    62,   -23,   117,
     -10,    81,   114,   119,   475,   478,   -34,   148,   105,   166,
     158,   225,   189,   160,   225,   291,    99,  -911,   214,   344,
     225,   371,   400,   505,   538,   187,   209,   225,   462,   468,
     474,   225,   658,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,   535,    52,  -911,
     479,   552,   441,    12,   108,   527,   430,   533,   543,   577,
    -911,   772,   164,    40,   208,   219,   560,   543,   576,  -911,
     289,   226,    53,   243,   569,  -911,  1088,   261,   262,   572,
    -911,  1088,   263,   276,   493,   303,   616,   494,   829,   250,
     250,   337,    53,   507,  -911,   570,  -911,   479,  -911,  1133,
     373,  -911,  1018,   379,   380,   548,   386,   559,   389,   574,
     390,   392,  -911,  -911,   540,   612,   -16,   439,  -911,   673,
     -33,  -911,  -911,   561,  -911,   562,  -911,  -911,   628,   277,
    -911,   639,   299,   686,    58,  -911,   643,  -911,   692,  -911,
     694,   696,  -911,   697,   701,  -911,   702,   703,   705,   712,
     713,  -911,  -911,   718,   720,   726,   727,   729,   732,  -911,
    -911,   733,   734,  -911,   743,  -911,  -911,  -911,   747,   749,
     750,   751,  -911,  -911,   761,   765,   -31,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,   773,   707,  -911,
     728,  -911,   730,   193,  -911,   678,   737,   681,   740,   205,
    -911,   741,   685,   744,   207,  -911,   662,    61,  -911,    64,
     795,   668,   753,  -911,    -9,   695,   706,   825,  -911,  -911,
     138,  -911,  -911,  -911,  -911,   784,   789,    91,  -911,  -911,
    -911,   715,   243,   243,   717,   719,   721,   724,   725,   735,
     745,   752,   754,   756,   243,   602,   757,   251,  -911,   833,
     836,   839,   845,   853,   854,  -911,  -911,  -911,   859,   860,
     876,   877,  -911,   889,  -911,   895,   897,  -911,  -911,   145,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,   566,   369,   172,  -911,  -911,  -911,   810,
     858,  -911,   779,  -911,  -911,  -911,   783,   829,   829,   785,
     796,   797,   798,   799,   800,   802,   804,   805,   808,   829,
     541,  -911,   173,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,   175,  -911,    93,
      29,   266,  -911,  -911,   279,  -911,  -911,   385,  -911,  -911,
     911,  -911,   401,  -911,  -911,  -911,  -911,  -911,   844,   891,
    -911,  -911,   846,   905,  -911,  -911,   861,   906,  -911,  -911,
     831,   835,   912,   425,   973,  -911,  -911,   946,   479,   849,
    -911,  -911,   850,   -14,   931,   856,    -1,   933,   243,  -911,
    -911,  -911,   971,   941,   864,   978,   979,   981,   983,   984,
     985,   990,   991,   517,  1005,   997,   998,   999,  1000,   977,
       6,   892,  1006,  1008,  1023,   989,   993,   772,    54,   994,
    1046,   944,  -911,  -911,  -911,    60,   945,   184,   947,  -911,
    -911,   950,   184,   952,  -911,  -911,    19,  -911,  -911,  -911,
    1003,   930,  1011,    30,  -911,    23,  -911,    43,  -911,   934,
     939,   239,   226,   -21,   956,    49,  -911,  -911,   243,   957,
     -49,   243,   243,   243,   243,   243,   243,   243,   243,   243,
     243,   635,   243,   243,   243,   243,   243,  -911,   243,  -911,
    -911,  1058,  1061,  1060,  1065,  1067,  1072,   184,  1077,  1079,
     531,  1087,  1089,  1093,  1088,   112,  1055,  1194,  -911,   803,
     218,  -911,  1014,  -911,    19,  1013,   465,   829,   829,   829,
     829,   829,   829,   829,   829,   829,   829,   709,   829,   829,
     829,   829,   829,   986,   250,   236,   252,  -911,  -911,  -911,
     243,   341,    59,   570,   988,   479,  1009,  1133,   260,  1111,
    1018,   273,  -911,  1031,  -911,  1033,  -911,  1040,  -911,  1116,
    1015,  1010,  1016,   243,  -911,  -911,  -911,  -911,  -911,   393,
    1017,  -911,  -911,   398,  1027,  1201,  -911,  -911,  1118,     4,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  1032,
    -911,  -911,  -911,  -911,  1021,  -911,  -911,  -911,   428,  -911,
    1097,  1100,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
     532,  1041,  1063,  1068,  1131,  1076,   184,  1135,  1057,   184,
    -911,  1173,  1175,  1066,  -911,   543,  1185,  -911,  -911,  -911,
     829,  -911,  -911,  -911,   404,  -911,  -911,  1070,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,   -30,   274,
    -911,  1152,   243,  1158,   -25,   663,   408,   722,   781,   794,
     865,   879,   885,   968,   982,  1012,  1024,  -911,   -49,   -49,
     957,   957,  1101,  1052,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
     440,   243,  -911,  1160,  1207,  -911,   452,  -911,  1084,  1069,
    1083,  1096,  1114,  1121,  1127,  1134,  1141,  1147,  1154,  -911,
     465,   465,  1013,  1013,  1101,  -911,  -911,  -911,   456,  -911,
     477,  1161,    29,  1090,  -911,  -911,    57,   243,  -911,  -911,
    -911,  -911,  -911,  -911,   487,  -911,  -911,  -911,   488,  -911,
    -911,  -911,  1094,  1229,  -911,  -911,  1217,  -911,   286,  -911,
     288,  -911,  1105,  -911,  -911,  -911,  1191,  -911,   414,  1198,
    -911,  -911,  -911,  -911,  -911,  -911,   184,    60,  1214,   184,
    1215,  1216,  -911,  1138,  -911,  -911,  1257,   375,   829,  1225,
      43,  -911,   753,   753,   753,   -21,  -911,   184,  -911,  1277,
    1231,  1284,  1307,   243,   243,  -911,   243,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  1171,  -911,
    1240,   243,  -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -911,    59,  -911,
    -911,  -911,   243,  1167,  -911,  -911,  1015,   243,  -911,  -911,
     530,  -911,   564,  1310,  1206,  1032,  -911,  -911,  -911,  1238,
    1239,  1245,   184,  1178,   184,   184,  -911,   243,  -911,  1249,
    -911,  -911,  1224,    56,   213,   482,   335,  1235,   243,  -911,
     243,  1222,    28,  1255,   787,   787,  -911,  -911,  1263,  1174,
    -911,  1365,  1273,  -911,  -911,  -911,  1268,  -911,   184,   184,
     184,  1270,  -911,  1248,  1250,  1279,  -911,   753,  -911,  -911,
    -911,   184,  -911,  1286,  1297,  1369,  1264,  1378,  1301,  -911,
    -911,  -911,   243,  -911,    18,  1295,  -911,  1311,   184,  -911,
    -911,  -911,   512,  1272,  -911,  -911,  -911,  1389,  1285,    74,
    1304,  1370,  -911,   184,   217,  1291,  -911,  -911,  -911,  1400,
    -911,  -911,   539,   542,   243,    72,  -911,  -911,  -911,  1287,
    1316,  1318,  -911,  -911,  -911,  -911,  1181,  -911,  -911,   243,
    -911,  -911,  -911,   184,   184,  -911,  1187,  1319,  -911,  -911,
     184,  -911
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   423,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     2,     4,    28,    29,    43,    44,    45,
      42,     5,     6,    11,     8,     9,    10,     7,    12,    13,
      14,    15,    16,    17,    18,    22,    24,    23,    19,    20,
      21,    25,    26,    27,    30,    31,    32,    37,    38,    33,
      34,    35,    36,    39,    40,    41,   392,     0,     0,   203,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   243,
     290,     0,     0,     0,     0,     0,     0,     0,     0,   123,
       0,     0,     0,     0,     0,   374,     0,     0,     0,     0,
     370,     0,     0,     0,    73,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   208,     0,   198,     0,   213,     0,
       0,   426,     0,     0,     0,    55,     0,    61,     0,    67,
       0,     0,     1,     3,     0,     0,   389,     0,   385,     0,
       0,   206,   207,     0,    79,     0,    46,   402,     0,     0,
     396,     0,     0,     0,     0,   112,     0,   480,     0,   497,
       0,     0,   485,     0,     0,   463,     0,     0,     0,     0,
       0,   477,   478,     0,     0,     0,     0,     0,     0,   499,
     473,     0,     0,   484,     0,   498,   479,   462,     0,     0,
       0,     0,   483,   481,     0,     0,     0,   295,   331,   320,
     296,   297,   298,   299,   300,   301,   302,   303,   304,   305,
     306,   307,   308,   309,   310,   311,   312,   313,   314,   315,
     316,   317,   318,   319,   321,   322,   323,   324,   325,   326,
     327,   328,   329,   330,   332,   333,   334,   239,     0,   292,
       0,   256,     0,     0,   253,     0,     0,     0,     0,     0,
     275,     0,     0,     0,     0,   269,     0,     0,   116,     0,
       0,     0,     0,   437,     0,     0,     0,     0,   493,   492,
       0,   408,   409,   410,   411,     0,     0,     0,   166,    84,
      85,    83,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   361,     0,
       0,     0,     0,     0,     0,   442,   443,   444,     0,     0,
       0,     0,   486,     0,   450,     0,     0,   379,   380,     0,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   232,   233,   234,   235,   236,   237,   238,   231,
     378,   376,   382,     0,     0,     0,   372,   369,    76,    71,
       0,    52,     0,    77,   141,   142,   161,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     424,   140,     0,   338,   343,   339,   340,   341,   342,   344,
     345,   346,   347,   348,   349,   350,   351,     0,    48,     0,
       0,     0,   211,   212,     0,   201,   202,     0,   218,   215,
       0,   432,     0,   431,   433,   428,   363,    58,    53,     0,
      49,    64,    59,     0,    50,    70,    65,     0,    51,   358,
       0,     0,     0,     0,     0,   383,   384,     0,     0,     0,
      80,    47,     0,     0,     0,     0,     0,     0,     0,   110,
     111,   244,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     241,     0,   255,   251,   252,   284,     0,   284,     0,   273,
     274,     0,   284,     0,   267,   268,     0,   114,   115,   107,
       0,     0,     0,     0,   135,     0,   136,     0,   131,     0,
       0,     0,     0,     0,     0,     0,   164,   165,     0,    91,
      92,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    81,     0,   359,
     360,     0,     0,     0,     0,     0,     0,   284,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   364,     0,
       0,    74,    72,    78,     0,   148,   149,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   163,   196,   197,
       0,     0,   188,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    56,    54,    62,    60,    68,    66,   393,     0,
     389,     0,     0,     0,   435,   205,   204,   405,   400,     0,
       0,   399,   394,     0,     0,     0,   464,   454,     0,     0,
     496,   457,   482,   445,   487,   488,   460,   461,   466,   469,
     470,   467,   475,   476,   465,   472,   471,   456,   455,     0,
     458,   459,   474,   494,     0,   495,   294,   291,     0,   240,
       0,     0,   279,   286,   280,   285,   282,   287,   281,   283,
       0,     0,     0,   261,     0,     0,   284,     0,     0,   284,
     247,     0,     0,     0,   109,     0,     0,   124,   133,   134,
       0,   138,   121,   120,     0,   119,   122,     0,   127,   125,
     490,   491,   407,   418,   420,   421,   422,   419,     0,   412,
     416,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    82,    87,    86,
      88,    89,    90,     0,   441,   449,   434,   440,   446,   447,
     489,   438,   448,   453,   452,   439,   436,   451,   381,   375,
       0,     0,   367,     0,     0,   371,     0,    75,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   139,
     144,   143,   145,   146,   147,   425,   337,   335,     0,   352,
       0,     0,     0,     0,   193,   194,     0,     0,   210,   209,
     200,   199,   217,   214,     0,   500,   430,   427,     0,    57,
      63,    69,     0,     0,   391,   390,     0,   401,     0,   395,
       0,   113,   502,   504,   506,   505,     0,   356,     0,     0,
     293,   242,   257,   289,   288,   254,   284,   284,     0,   284,
       0,     0,   272,     0,   246,   245,     0,     0,     0,     0,
       0,   129,     0,     0,     0,     0,   406,   284,   417,     0,
       0,     0,     0,     0,     0,   103,     0,   104,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,     0,   377,
       0,     0,   365,   373,   162,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   336,   353,   195,   187,   186,
     190,   191,     0,     0,   216,   429,   389,     0,   386,   403,
       0,   397,     0,     0,     0,     0,   468,   501,   258,     0,
       0,     0,   284,     0,   284,   284,   270,     0,   108,     0,
     137,   118,     0,     0,     0,     0,   413,     0,     0,   169,
       0,   177,     0,     0,   105,   106,   362,   368,     0,     0,
     192,     0,     0,   404,   398,   503,     0,   357,   284,   284,
     284,     0,   278,     0,     0,     0,   160,     0,   132,   128,
     126,   284,   414,     0,     0,     0,   172,     0,     0,   168,
     366,   189,     0,   387,   284,   263,   259,   262,   284,   276,
     271,   117,     0,     0,   171,   170,   176,     0,   174,     0,
       0,     0,   355,   284,     0,     0,   130,   415,   173,     0,
     250,   183,     0,     0,     0,     0,   182,   181,   388,     0,
     264,     0,   277,   175,   249,   248,     0,   180,   167,     0,
     179,   178,   354,   284,   284,   185,     0,   265,   260,   184,
     284,   266
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -911,  -911,  1405,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,  -911,  -911,  -296,  -911,  -911,
    -911,  -911,  -102,  -172,  -911,  -911,  1172,  -911,   598,  -911,
    -911,  -911,  -911,  -911,  -911,  -779,  -494,  -108,  -488,  -911,
    -911,  -911,  1320,  -253,  -911,  -911,  -911,  -911,   659,  -911,
    -911,   851,  -911,  -911,  1002,  -911,  -911,   852,  -911,  -911,
    -100,   -20,  -520,   442,  -911,  -911,  1195,  -911,  -911,  -910,
    -911,  -911,  1180,  -911,  -911,  1190,  -794,  -471,  -911,  -911,
     974,  -911,  1330,   868,  -911,   549,  -911,  -911,  -911,  -911,
    1146,  -911,  -911,  -911,  -911,  -911,  -911,   901,  1345,  -911,
    -911,  -911,  1312,  -584,  -911,  -911,  -911,  -911,  -911,   948,
    -911,   613,  -678,  -911,  -911,  -911,  -911,  -911,   857,  -911,
     -79,  -911,  1361,  -911,  -911,  -911,  -911,  -911,  -911,  -911,
     -92,  -911,  -911,  -112,  -477,  -911,  -911,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,   -93,   -89,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,  -911,   -88,  -911,  -911,  -911,  -911,
    -911,   -87,   -74,   -73,   -72,   -71,   -69,  -911,  -911,  -911,
    -911,  -911,  -911,  -911,   -68,   -67,   -66,  -911,  -911,  -911,
    -911,  -911,   834,  -911,   992
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    42,    43,    44,    45,    46,    47,    48,    49,    50,
     146,   148,   150,   125,    51,    52,    53,   315,   726,    54,
     281,    55,   174,   175,    56,   277,   278,   704,   705,    57,
     282,   854,   853,   932,   707,   513,   514,   515,   516,   391,
      58,    59,   297,   298,   942,  1015,    60,   601,   602,    61,
     414,   415,    62,   160,   161,    63,   411,   412,    64,   417,
     337,   102,   693,  1017,    65,   263,   264,   265,   681,   918,
      66,   274,   275,    67,   269,   270,   682,   919,    68,   216,
     217,    69,   392,   393,    70,   827,   828,    71,    72,   317,
     318,    73,    74,   364,    75,    76,    77,   338,   339,    78,
      79,   157,   158,   444,    80,    81,    82,    83,   290,   291,
     718,   719,   720,    84,   128,   593,    85,   422,   423,   340,
     341,   342,   343,   344,   345,   346,   347,   348,   349,   350,
     351,   352,   353,   354,   355,   356,   357,   358,   220,   221,
     222,   223,   224,   225,   226,   395,   396,   229,   230,   231,
     232,   233,   234,   235,   236,   397,   238,   239,   240,   241,
     242,   398,   399,   400,   401,   402,   403,   359,   249,   250,
     360,   292,   293,   294,   404,   405,   406,   254,   255,   256,
     424,   665,   823,   639,   640
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       117,   118,   460,   122,   123,   279,   529,   530,   227,   219,
     131,   162,   228,   237,   243,   140,   143,   144,   541,   698,
     390,   151,   218,   413,   683,   699,   685,   244,   245,   246,
     247,   688,   248,   251,   252,   253,   813,   394,   394,   418,
     706,   858,   421,   920,   527,   822,   700,   657,   713,   986,
     673,   697,   260,   165,   768,    99,   155,   320,   416,   675,
     690,   598,   374,    94,   862,   295,   672,   567,    99,   976,
     599,   375,   722,   933,   934,   935,    93,   968,   702,   459,
     544,   545,   507,   546,   863,   509,   750,   677,   374,   900,
     703,   794,   673,   448,   674,   487,   855,   375,   901,   376,
     795,   675,   676,   295,  1010,   295,  1010,   460,   714,   261,
     442,    88,   526,   690,   597,   690,   111,   517,   449,   856,
     488,   156,   126,    99,  1038,   376,    98,   173,   113,   677,
     276,    96,   715,   173,   443,   628,   716,   717,   678,   103,
     127,  1027,   518,  1011,    95,   680,   691,   692,   632,   262,
     588,   589,   590,   591,   977,   592,   527,   377,   378,   824,
     638,   166,   635,   379,   380,   381,   382,   383,   384,   385,
     386,   387,   701,   296,  1001,   723,   679,   978,   388,   600,
     389,    99,   512,   377,   378,   796,    99,   680,   992,   379,
     380,   381,   382,   383,   384,   385,   386,   387,   724,  1012,
    1013,  1012,  1013,   667,   388,   260,   389,   902,   512,   797,
    1021,   296,    97,   296,   493,   840,   673,   266,   843,   271,
     266,  1028,  1029,   672,  1014,   675,   499,    99,   504,    99,
     104,   271,   725,   257,   969,   727,   728,   729,   730,   731,
     732,   733,   734,   735,   736,   374,   738,   739,   740,   741,
     742,   674,   743,   677,   375,   177,    99,   167,   858,   676,
     178,   759,   261,   105,   522,   168,   100,   101,   106,   575,
     576,   564,   549,   764,   267,   299,   272,   267,    99,   181,
     182,   587,   376,   184,   300,   710,   185,    99,   272,   523,
     258,    91,   287,   186,    99,   678,   565,   112,   564,   594,
      92,   594,   262,   288,   791,    99,   711,   115,   116,   120,
     121,   680,   301,   259,   268,   114,   273,   268,   203,   289,
     316,    99,   951,   570,   595,   207,   596,   816,   273,    99,
     257,   257,   257,   679,    86,    87,   138,   139,   119,   283,
     377,   378,    99,   713,   211,   257,   379,   380,   381,   382,
     383,   384,   385,   386,   387,    99,   212,    99,   141,   142,
     124,   388,   213,   389,   129,   512,   921,   765,   923,   162,
     302,   303,   369,   706,   214,   215,   304,   305,   306,   307,
     308,   309,   310,   311,   312,   787,   937,   258,   258,   258,
     568,   313,   603,   314,   227,   219,   928,   284,   228,   237,
     243,   789,   258,   714,   713,   605,   257,   285,   218,   803,
     361,   362,   366,   244,   245,   246,   247,   604,   248,   251,
     252,   253,   807,   133,   857,   367,   860,   715,   453,   370,
     606,   716,   717,   629,   454,   909,   633,   911,   569,   698,
     698,   698,   257,   155,   173,   699,   699,   699,   257,   257,
     456,   961,   371,   963,   964,   428,   457,   621,   432,   436,
     445,   257,   257,   258,   714,   880,   622,   257,   668,   769,
     770,   771,   772,   773,   774,   775,   776,   777,   778,  1016,
     780,   781,   782,   783,   784,   971,   408,   985,   715,   987,
     792,   413,   716,   717,   130,  1030,   793,   257,   698,   258,
     993,   903,   394,   970,   699,   258,   258,   802,   156,   257,
     421,   607,   429,  1002,   374,   433,   437,  1005,   258,   258,
     132,   257,   419,   375,   258,   257,   416,   610,   425,   426,
     850,   145,  1020,  1006,   866,   430,   608,   147,   434,   438,
     915,   439,   817,   149,   374,   760,   257,   819,   154,   649,
     766,   376,   611,   375,   258,   851,   257,   257,   650,   867,
      89,    90,  1037,   753,   833,   916,   258,   943,   944,  1041,
     945,  1024,   754,   834,  1025,   788,   790,   830,   258,   170,
     844,   376,   258,   845,   163,   948,   159,   171,   804,   879,
     164,   808,   849,   847,   590,   591,   169,   592,   299,   257,
     176,   883,   172,   258,   280,   895,   949,   300,   368,   377,
     378,   952,   173,   258,   258,   379,   380,   381,   382,   383,
     384,   385,   386,   387,   107,   108,   896,   109,   110,   276,
     388,   965,   389,   257,   512,   301,   904,   905,   316,   377,
     378,   363,   973,   373,   974,   379,   380,   381,   382,   383,
     384,   385,   386,   387,   134,   135,   258,   372,   152,   330,
     388,   410,   389,   427,   512,     1,     2,     3,   588,   589,
     590,   591,     4,   592,   431,   460,     5,     6,     7,   953,
       8,   441,     9,    10,    11,    12,  1000,   136,   137,   435,
     258,   440,   566,   302,   303,    13,   447,   452,    14,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   455,   458,
     450,   451,   461,   954,   313,   462,   314,   463,  1026,   464,
     465,    15,    16,    17,   466,   467,   468,    18,   469,   542,
     543,   544,   545,  1036,   546,   470,   471,    19,    20,    21,
     929,   472,    22,   473,    23,    24,    25,    26,    27,   474,
     475,   547,   476,    28,    29,   477,   478,   479,    30,    31,
      32,    33,   542,   543,   544,   545,   480,   546,    34,    35,
     481,    36,   482,   483,   484,    37,   490,   177,    38,    39,
      40,    41,   178,   179,   485,   374,   737,   180,   486,   864,
     542,   543,   544,   545,   375,   546,   489,   491,   910,   492,
     912,   181,   182,   183,   495,   184,   496,   497,   185,   498,
     501,   502,   506,   503,   865,   186,   187,   188,   510,   511,
     189,   190,   376,   191,   192,   193,   194,   195,   196,   197,
     198,   199,   200,   201,   202,   299,   588,   589,   590,   591,
     203,   592,   204,   205,   300,   206,   519,   207,   521,   542,
     543,   544,   545,   524,   546,   208,   551,   520,   525,   552,
     779,   374,   553,   209,   210,   528,   211,   531,   554,   532,
     375,   533,   301,   868,   534,   535,   555,   556,   212,   159,
     377,   378,   557,   558,   213,   536,   379,   380,   381,   382,
     383,   384,   385,   386,   387,   537,   214,   215,   376,   559,
     560,   388,   538,   389,   539,   512,   540,   548,   542,   543,
     544,   545,   561,   546,   542,   543,   544,   545,   562,   546,
     563,   542,   543,   544,   545,   571,   546,   572,   573,   763,
     302,   303,   869,   574,   609,   577,   304,   305,   306,   307,
     308,   309,   310,   311,   312,   870,   578,   579,   580,   581,
     582,   313,   583,   314,   584,   585,   377,   378,   586,   612,
     613,   614,   379,   380,   381,   382,   383,   384,   385,   386,
     387,     1,     2,     3,   615,   617,   616,   388,     4,   389,
     618,   620,     5,     6,     7,   619,     8,   624,     9,    10,
      11,    12,   542,   543,   544,   545,   623,   546,   626,   627,
     630,    13,   634,   636,    14,   631,   542,   543,   544,   545,
     637,   546,   542,   543,   544,   545,   871,   546,   638,   641,
     642,   319,   643,   648,   644,   645,   646,    15,    16,    17,
     872,   647,   320,    18,   321,   322,   873,   651,   652,   653,
     654,   655,   659,    19,    20,    21,   656,   660,    22,   661,
      23,    24,    25,    26,    27,   662,   323,   324,   663,    28,
      29,   186,   664,   669,    30,    31,    32,    33,   283,   670,
     671,   684,   694,   686,    34,    35,   687,    36,   689,   695,
     696,    37,   721,   708,    38,    39,    40,    41,   709,   546,
     325,   319,   326,   745,   327,   542,   543,   544,   545,   744,
     546,   746,   320,   329,   321,   322,   747,   330,   748,   542,
     543,   544,   545,   749,   546,   331,   332,   333,   751,   874,
     752,   334,   335,   336,   761,   159,   323,   324,   755,   767,
     756,   186,   420,   875,   757,   785,   319,   799,   283,   542,
     543,   544,   545,   805,   546,   592,   809,   320,   810,   321,
     322,   542,   543,   544,   545,   811,   546,   812,   801,   822,
     325,   814,   326,   876,   327,   443,   831,   815,   818,   832,
     328,   323,   324,   329,   829,   877,   186,   330,   820,   542,
     543,   544,   545,   283,   546,   331,   332,   333,   826,   836,
     835,   334,   335,   336,   837,   159,   588,   589,   590,   591,
     838,   592,   839,   878,   841,   325,   842,   326,   848,   327,
     588,   589,   590,   591,   844,   592,   845,   846,   329,   852,
     885,   859,   330,   588,   589,   590,   591,   861,   592,   881,
     331,   332,   333,    -1,   886,   884,   334,   335,   336,   899,
     159,   588,   589,   590,   591,   906,   592,   887,   588,   589,
     590,   591,   907,   592,   588,   589,   590,   591,   913,   592,
     914,   588,   589,   590,   591,   888,   592,   917,   588,   589,
     590,   591,   889,   592,   588,   589,   590,   591,   890,   592,
     927,   588,   589,   590,   591,   891,   592,   926,   542,   543,
     544,   545,   892,   546,   542,   543,   544,   545,   893,   546,
     938,   542,   543,   544,   545,   894,   546,   940,   542,   543,
     544,   545,   897,   546,   542,   543,   544,   545,   950,   546,
     946,   542,   543,   544,   545,   981,   546,   962,   542,   543,
     544,   545,  1035,   546,   542,   543,   544,   545,  1039,   546,
     922,   924,   925,   762,   542,   543,   544,   545,   941,   546,
     821,   955,   588,   589,   590,   591,   882,   592,   542,   543,
     544,   545,   956,   546,   958,   959,   908,   542,   543,   544,
     545,   960,   546,   967,   930,   975,   588,   589,   590,   591,
     939,   592,   542,   543,   544,   545,   972,   546,   982,   947,
     542,   543,   544,   545,   984,   546,   988,   989,   966,   990,
     542,   543,   544,   545,   979,   546,   542,   543,   544,   545,
     996,   546,   980,   542,   543,   544,   545,   997,   546,   998,
     999,  1003,   983,  1007,   542,   543,   544,   545,   991,   546,
    1008,   542,   543,   544,   545,   994,   546,  1004,  1009,  1019,
    1022,  1023,  1033,  1032,  1034,  1040,   995,   153,   931,   508,
     625,   898,   409,  1018,   505,   798,   800,  1031,   494,   500,
     407,   666,   786,   550,   957,   758,   365,   806,   936,   446,
     712,   286,   658,   825
  };

  /* YYCHECK.  */
  const unsigned short int
  parser::yycheck_[] =
  {
        20,    21,   174,    23,    24,   107,   302,   303,   101,   101,
      30,    90,   101,   101,   101,    35,    36,    37,   314,   513,
     128,    41,   101,   135,   495,   513,   497,   101,   101,   101,
     101,   502,   101,   101,   101,   101,   620,   129,   130,   139,
     517,   719,   142,   837,   297,    41,    23,    41,    69,   959,
      32,    21,    12,    41,   574,    69,     4,    14,   137,    41,
      41,    32,    32,    69,    89,    12,     6,   363,    69,    41,
      41,    41,    23,   852,   853,   854,    41,    21,    35,    21,
     129,   130,    21,   132,   109,    21,   557,    69,    32,    32,
      47,    32,    32,   126,    34,   126,   126,    41,    41,    69,
      41,    41,    42,    12,    32,    12,    32,   279,   129,    69,
     126,   149,    21,    41,    21,    41,   150,   126,   151,   149,
     151,    69,    23,    69,  1034,    69,   149,    69,    23,    69,
      69,    69,   153,    69,   150,   149,   157,   158,    78,   149,
      41,    69,   151,    69,   150,   127,   127,   128,   149,   109,
     127,   128,   129,   130,   126,   132,   409,   127,   128,   155,
     154,   149,   458,   133,   134,   135,   136,   137,   138,   139,
     140,   141,   149,   120,   156,   126,   116,   149,   148,   150,
     150,    69,   152,   127,   128,   126,    69,   127,   967,   133,
     134,   135,   136,   137,   138,   139,   140,   141,   149,   127,
     128,   127,   128,   149,   148,    12,   150,   150,   152,   150,
    1004,   120,   150,   120,    21,   686,    32,    12,   689,    12,
      12,   149,   150,     6,   150,    41,    21,    69,    21,    69,
     149,    12,   528,    69,    21,   531,   532,   533,   534,   535,
     536,   537,   538,   539,   540,    32,   542,   543,   544,   545,
     546,    34,   548,    69,    41,     5,    69,   149,   936,    42,
      10,   149,    69,   149,   126,   157,   149,   150,   149,   377,
     378,   126,    21,   569,    69,    32,    69,    69,    69,    29,
      30,   389,    69,    33,    41,    46,    36,    69,    69,   151,
     126,    23,    66,    43,    69,    78,   151,   149,   126,   126,
      32,   126,   109,    77,   600,    69,    67,   149,   150,   149,
     150,   127,    69,   149,   109,   149,   109,   109,    68,    93,
      69,    69,   906,   151,   151,    75,   151,   623,   109,    69,
      69,    69,    69,   116,   149,   150,   149,   150,   149,    50,
     127,   128,    69,    69,    94,    69,   133,   134,   135,   136,
     137,   138,   139,   140,   141,    69,   106,    69,   149,   150,
      69,   148,   112,   150,   150,   152,   837,   149,   839,   448,
     127,   128,    69,   850,   124,   125,   133,   134,   135,   136,
     137,   138,   139,   140,   141,   149,   857,   126,   126,   126,
      21,   148,   126,   150,   487,   487,    21,   108,   487,   487,
     487,   149,   126,   129,    69,   126,    69,   118,   487,   149,
     149,   149,   149,   487,   487,   487,   487,   151,   487,   487,
     487,   487,   149,    23,   150,   149,   722,   153,   151,   126,
     151,   157,   158,   453,   157,   149,   456,   149,    69,   933,
     934,   935,    69,     4,    69,   933,   934,   935,    69,    69,
     151,   922,   149,   924,   925,    69,   157,    32,    69,    69,
      21,    69,    69,   126,   129,   761,    41,    69,   488,   577,
     578,   579,   580,   581,   582,   583,   584,   585,   586,   999,
     588,   589,   590,   591,   592,   150,   149,   958,   153,   960,
     149,   603,   157,   158,   150,  1015,   155,    69,   992,   126,
     971,   797,   594,    21,   992,   126,   126,   607,    69,    69,
     610,   126,   126,   984,    32,   126,   126,   988,   126,   126,
     149,    69,   149,    41,   126,    69,   605,   126,   149,   149,
     126,    69,  1003,    21,   126,   149,   151,    69,   149,   149,
     126,   149,   149,    69,    32,   565,    69,   149,    13,    32,
     570,    69,   151,    41,   126,   151,    69,    69,    41,   151,
     149,   150,  1033,    32,    32,   151,   126,   863,   864,  1040,
     866,    32,    41,    41,    32,   595,   596,   149,   126,   149,
      41,    69,   126,    41,    32,   881,   107,   157,   608,   149,
     149,   611,   700,   695,   129,   130,    69,   132,    32,    69,
      23,   149,    69,   126,    28,   149,   902,    41,   115,   127,
     128,   907,    69,   126,   126,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   149,   150,   149,   149,   150,    69,
     148,   927,   150,    69,   152,    69,   149,   149,    69,   127,
     128,    69,   938,   149,   940,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   149,   150,   126,    41,     0,    89,
     148,   154,   150,   115,   152,     7,     8,     9,   127,   128,
     129,   130,    14,   132,   115,   847,    18,    19,    20,   149,
      22,    69,    24,    25,    26,    27,   982,   149,   150,   115,
     126,   151,   126,   127,   128,    37,    23,    69,    40,   133,
     134,   135,   136,   137,   138,   139,   140,   141,    69,    23,
     149,   149,    69,   149,   148,    23,   150,    23,  1014,    23,
      23,    63,    64,    65,    23,    23,    23,    69,    23,   127,
     128,   129,   130,  1029,   132,    23,    23,    79,    80,    81,
     848,    23,    84,    23,    86,    87,    88,    89,    90,    23,
      23,   149,    23,    95,    96,    23,    23,    23,   100,   101,
     102,   103,   127,   128,   129,   130,    23,   132,   110,   111,
      23,   113,    23,    23,    23,   117,    69,     5,   120,   121,
     122,   123,    10,    11,    23,    32,   151,    15,    23,   126,
     127,   128,   129,   130,    41,   132,    23,    69,   818,    69,
     820,    29,    30,    31,   126,    33,    69,   126,    36,    69,
      69,   126,   150,    69,   151,    43,    44,    45,    23,   151,
      48,    49,    69,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    32,   127,   128,   129,   130,
      68,   132,    70,    71,    41,    73,   151,    75,    23,   127,
     128,   129,   130,    69,   132,    83,    23,   151,    69,    23,
     151,    32,    23,    91,    92,   150,    94,   150,    23,   150,
      41,   150,    69,   151,   150,   150,    23,    23,   106,   107,
     127,   128,    23,    23,   112,   150,   133,   134,   135,   136,
     137,   138,   139,   140,   141,   150,   124,   125,    69,    23,
      23,   148,   150,   150,   150,   152,   150,   150,   127,   128,
     129,   130,    23,   132,   127,   128,   129,   130,    23,   132,
      23,   127,   128,   129,   130,   115,   132,    69,   149,   126,
     127,   128,   151,   150,    23,   150,   133,   134,   135,   136,
     137,   138,   139,   140,   141,   151,   150,   150,   150,   150,
     150,   148,   150,   150,   150,   150,   127,   128,   150,   115,
      69,   115,   133,   134,   135,   136,   137,   138,   139,   140,
     141,     7,     8,     9,    69,    69,   115,   148,    14,   150,
     149,    69,    18,    19,    20,   150,    22,    41,    24,    25,
      26,    27,   127,   128,   129,   130,    23,   132,   149,   149,
      69,    37,    69,    32,    40,   149,   127,   128,   129,   130,
      69,   132,   127,   128,   129,   130,   151,   132,   154,    41,
      41,     3,    41,    32,    41,    41,    41,    63,    64,    65,
     151,    41,    14,    69,    16,    17,   151,    32,    41,    41,
      41,    41,   150,    79,    80,    81,    69,    41,    84,    41,
      86,    87,    88,    89,    90,    32,    38,    39,    69,    95,
      96,    43,    69,    69,   100,   101,   102,   103,    50,    23,
     126,   126,    69,   126,   110,   111,   126,   113,   126,   149,
      69,   117,   126,   149,   120,   121,   122,   123,   149,   132,
      72,     3,    74,    32,    76,   127,   128,   129,   130,    41,
     132,    41,    14,    85,    16,    17,    41,    89,    41,   127,
     128,   129,   130,    41,   132,    97,    98,    99,    41,   151,
      41,   103,   104,   105,    69,   107,    38,    39,    41,   115,
      41,    43,   114,   151,    41,   149,     3,   149,    50,   127,
     128,   129,   130,    32,   132,   132,   115,    14,   115,    16,
      17,   127,   128,   129,   130,   115,   132,    41,   149,    41,
      72,   151,    74,   151,    76,   150,    69,   151,   151,    69,
      82,    38,    39,    85,   153,   151,    43,    89,   151,   127,
     128,   129,   130,    50,   132,    97,    98,    99,   156,   126,
     149,   103,   104,   105,   126,   107,   127,   128,   129,   130,
      69,   132,   126,   151,    69,    72,   149,    74,    23,    76,
     127,   128,   129,   130,    41,   132,    41,   151,    85,   149,
     151,    69,    89,   127,   128,   129,   130,    69,   132,    69,
      97,    98,    99,   132,   151,   151,   103,   104,   105,   149,
     107,   127,   128,   129,   130,   151,   132,   151,   127,   128,
     129,   130,    23,   132,   127,   128,   129,   130,   153,   132,
      69,   127,   128,   129,   130,   151,   132,    69,   127,   128,
     129,   130,   151,   132,   127,   128,   129,   130,   151,   132,
      23,   127,   128,   129,   130,   151,   132,   149,   127,   128,
     129,   130,   151,   132,   127,   128,   129,   130,   151,   132,
      23,   127,   128,   129,   130,   151,   132,    23,   127,   128,
     129,   130,   151,   132,   127,   128,   129,   130,   151,   132,
     149,   127,   128,   129,   130,   151,   132,   149,   127,   128,
     129,   130,   151,   132,   127,   128,   129,   130,   151,   132,
     126,   126,   126,   149,   127,   128,   129,   130,    41,   132,
     149,    41,   127,   128,   129,   130,   149,   132,   127,   128,
     129,   130,   156,   132,   126,   126,   149,   127,   128,   129,
     130,   126,   132,   149,   149,   153,   127,   128,   129,   130,
     149,   132,   127,   128,   129,   130,   151,   132,    23,   149,
     127,   128,   129,   130,   126,   132,   126,   149,   149,   149,
     127,   128,   129,   130,   149,   132,   127,   128,   129,   130,
      41,   132,   149,   127,   128,   129,   130,   153,   132,    41,
     119,   126,   149,   151,   127,   128,   129,   130,   149,   132,
      41,   127,   128,   129,   130,   149,   132,   126,   153,    69,
     149,    41,   126,   156,   126,   126,   149,    42,   850,   277,
     448,   792,   132,   149,   274,   603,   605,  1015,   263,   269,
     130,   487,   594,   317,   915,   564,   121,   610,   855,   157,
     522,   110,   480,   639
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,     9,    14,    18,    19,    20,    22,    24,
      25,    26,    27,    37,    40,    63,    64,    65,    69,    79,
      80,    81,    84,    86,    87,    88,    89,    90,    95,    96,
     100,   101,   102,   103,   110,   111,   113,   117,   120,   121,
     122,   123,   160,   161,   162,   163,   164,   165,   166,   167,
     168,   173,   174,   175,   178,   180,   183,   188,   199,   200,
     205,   208,   211,   214,   217,   223,   229,   232,   237,   240,
     243,   246,   247,   250,   251,   253,   254,   255,   258,   259,
     263,   264,   265,   266,   272,   275,   149,   150,   149,   149,
     150,    23,    32,    41,    69,   150,    69,   150,   149,    69,
     149,   150,   220,   149,   149,   149,   149,   149,   150,   149,
     150,   150,   149,    23,   149,   149,   150,   220,   220,   149,
     149,   150,   220,   220,    69,   172,    23,    41,   273,   150,
     150,   220,   149,    23,   149,   150,   149,   150,   149,   150,
     220,   149,   150,   220,   220,    69,   169,    69,   170,    69,
     171,   220,     0,   161,    13,     4,    69,   260,   261,   107,
     212,   213,   279,    32,   149,    41,   149,   149,   157,    69,
     149,   157,    69,    69,   181,   182,    23,     5,    10,    11,
      15,    29,    30,    31,    33,    36,    43,    44,    45,    48,
      49,    51,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    68,    70,    71,    73,    75,    83,    91,
      92,    94,   106,   112,   124,   125,   238,   239,   279,   289,
     297,   298,   299,   300,   301,   302,   303,   304,   305,   306,
     307,   308,   309,   310,   311,   312,   313,   314,   315,   316,
     317,   318,   319,   320,   321,   322,   323,   324,   325,   327,
     328,   333,   334,   335,   336,   337,   338,    69,   126,   149,
      12,    69,   109,   224,   225,   226,    12,    69,   109,   233,
     234,    12,    69,   109,   230,   231,    69,   184,   185,   181,
      28,   179,   189,    50,   108,   118,   281,    66,    77,    93,
     267,   268,   330,   331,   332,    12,   120,   201,   202,    32,
      41,    69,   127,   128,   133,   134,   135,   136,   137,   138,
     139,   140,   141,   148,   150,   176,    69,   248,   249,     3,
      14,    16,    17,    38,    39,    72,    74,    76,    82,    85,
      89,    97,    98,    99,   103,   104,   105,   219,   256,   257,
     278,   279,   280,   281,   282,   283,   284,   285,   286,   287,
     288,   289,   290,   291,   292,   293,   294,   295,   296,   326,
     329,   149,   149,    69,   252,   257,   149,   149,   115,    69,
     126,   149,    41,   149,    32,    41,    69,   127,   128,   133,
     134,   135,   136,   137,   138,   139,   140,   141,   148,   150,
     196,   198,   241,   242,   289,   304,   305,   314,   320,   321,
     322,   323,   324,   325,   333,   334,   335,   241,   149,   201,
     154,   215,   216,   292,   209,   210,   279,   218,   219,   149,
     114,   219,   276,   277,   339,   149,   149,   115,    69,   126,
     149,   115,    69,   126,   149,   115,    69,   126,   149,   149,
     151,    69,   126,   150,   262,    21,   261,    23,   126,   151,
     149,   149,    69,   151,   157,    69,   151,   157,    23,    21,
     182,    69,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,    23,    23,    23,   126,   151,    23,
      69,    69,    69,    21,   225,   126,    69,   126,    69,    21,
     234,    69,   126,    69,    21,   231,   150,    21,   185,    21,
      23,   151,   152,   194,   195,   196,   197,   126,   151,   151,
     151,    23,   126,   151,    69,    69,    21,   202,   150,   176,
     176,   150,   150,   150,   150,   150,   150,   150,   150,   150,
     150,   176,   127,   128,   129,   130,   132,   149,   150,    21,
     249,    23,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,   126,   151,   126,   176,    21,    69,
     151,   115,    69,   149,   150,   196,   196,   150,   150,   150,
     150,   150,   150,   150,   150,   150,   150,   196,   127,   128,
     129,   130,   132,   274,   126,   151,   151,    21,    32,    41,
     150,   206,   207,   126,   151,   126,   151,   126,   151,    23,
     126,   151,   115,    69,   115,    69,   115,    69,   149,   150,
      69,    32,    41,    23,    41,   213,   149,   149,   149,   220,
      69,   149,   149,   220,    69,   176,    32,    69,   154,   342,
     343,    41,    41,    41,    41,    41,    41,    41,    32,    32,
      41,    32,    41,    41,    41,    41,    69,    41,   343,   150,
      41,    41,    32,    69,    69,   340,   239,   149,   220,    69,
      23,   126,     6,    32,    34,    41,    42,    69,    78,   116,
     127,   227,   235,   236,   126,   236,   126,   126,   236,   126,
      41,   127,   128,   221,    69,   149,    69,    21,   195,   197,
      23,   149,    35,    47,   186,   187,   293,   193,   149,   149,
      46,    67,   268,    69,   129,   153,   157,   158,   269,   270,
     271,   126,    23,   126,   149,   176,   177,   176,   176,   176,
     176,   176,   176,   176,   176,   176,   176,   151,   176,   176,
     176,   176,   176,   176,    41,    32,    41,    41,    41,    41,
     236,    41,    41,    32,    41,    41,    41,    41,   256,   149,
     220,    69,   149,   126,   176,   149,   220,   115,   221,   196,
     196,   196,   196,   196,   196,   196,   196,   196,   196,   151,
     196,   196,   196,   196,   196,   149,   242,   149,   220,   149,
     220,   176,   149,   155,    32,    41,   126,   150,   216,   149,
     210,   149,   219,   149,   220,    32,   277,   149,   220,   115,
     115,   115,    41,   262,   151,   151,   176,   149,   151,   149,
     151,   149,    41,   341,   155,   341,   156,   244,   245,   153,
     149,    69,    69,    32,    41,   149,   126,   126,    69,   126,
     236,    69,   149,   236,    41,    41,   151,   181,    23,   196,
     126,   151,   149,   191,   190,   126,   149,   150,   271,    69,
     176,    69,    89,   109,   126,   151,   126,   151,   151,   151,
     151,   151,   151,   151,   151,   151,   151,   151,   151,   149,
     176,    69,   149,   149,   151,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   149,   149,   151,   207,   149,
      32,    41,   150,   176,   149,   149,   151,    23,   149,   149,
     220,   149,   220,   153,    69,   126,   151,    69,   228,   236,
     235,   236,   126,   236,   126,   126,   149,    23,    21,   196,
     149,   187,   192,   194,   194,   194,   270,   236,    23,   149,
      23,    41,   203,   176,   176,   176,   149,   149,   176,   176,
     151,   262,   176,   149,   149,    41,   156,   244,   126,   126,
     126,   236,   149,   236,   236,   176,   149,   149,    21,    21,
      21,   150,   151,   176,   176,   153,    41,   126,   149,   149,
     149,   151,    23,   149,   126,   236,   228,   236,   126,   149,
     149,   149,   194,   236,   149,   149,    41,   153,    41,   119,
     176,   156,   236,   126,   126,   236,    21,   151,    41,   153,
      32,    69,   127,   128,   150,   204,   221,   222,   149,    69,
     236,   235,   149,    41,    32,    32,   176,    69,   149,   150,
     221,   222,   156,   126,   126,   151,   176,   236,   228,   151,
     126,   236
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
     395,   396,   397,   398,   399,   400,   401,   402,   403,    59,
      40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   159,   160,   160,   161,   161,   161,   161,   161,   161,
     161,   161,   161,   161,   161,   161,   161,   161,   161,   161,
     161,   161,   161,   161,   161,   161,   161,   161,   161,   161,
     161,   161,   161,   161,   161,   161,   161,   161,   161,   161,
     161,   161,   162,   162,   162,   162,   163,   163,   164,   165,
     166,   167,   168,   169,   169,   169,   169,   169,   169,   170,
     170,   170,   170,   170,   170,   171,   171,   171,   171,   171,
     171,   172,   172,   172,   172,   172,   172,   173,   173,   174,
     174,   175,   176,   176,   176,   176,   176,   176,   176,   176,
     176,   176,   176,   176,   176,   176,   176,   176,   176,   176,
     176,   176,   176,   176,   176,   177,   177,   178,   178,   179,
     180,   181,   181,   182,   183,   184,   184,   185,   186,   186,
     187,   187,   187,   189,   188,   190,   188,   191,   188,   192,
     188,   193,   188,   194,   194,   194,   194,   195,   195,   196,
     196,   196,   196,   196,   196,   196,   196,   196,   196,   196,
     196,   196,   196,   196,   196,   196,   196,   196,   196,   196,
     197,   198,   198,   199,   200,   201,   201,   202,   202,   202,
     202,   202,   203,   203,   203,   203,   203,   203,   204,   204,
     204,   204,   204,   204,   204,   204,   205,   206,   206,   207,
     207,   207,   207,   207,   207,   207,   207,   207,   208,   208,
     209,   209,   210,   211,   211,   212,   212,   213,   214,   214,
     215,   215,   216,   217,   217,   217,   217,   218,   218,   219,
     219,   219,   219,   219,   219,   219,   219,   219,   219,   219,
     219,   219,   219,   219,   219,   219,   219,   219,   219,   220,
     220,   220,   220,   220,   220,   221,   221,   221,   222,   222,
     222,   223,   224,   224,   225,   226,   226,   226,   227,   227,
     227,   227,   227,   228,   228,   228,   228,   229,   230,   230,
     231,   231,   231,   232,   233,   233,   234,   234,   234,   235,
     235,   235,   235,   235,   236,   236,   236,   236,   236,   236,
     237,   237,   237,   237,   238,   238,   239,   239,   239,   239,
     239,   239,   239,   239,   239,   239,   239,   239,   239,   239,
     239,   239,   239,   239,   239,   239,   239,   239,   239,   239,
     239,   239,   239,   239,   239,   239,   239,   239,   239,   239,
     239,   239,   239,   239,   239,   240,   240,   241,   241,   242,
     242,   242,   242,   242,   242,   242,   242,   242,   242,   242,
     242,   242,   243,   243,   244,   244,   245,   245,   246,   247,
     248,   248,   249,   250,   251,   252,   252,   252,   252,   253,
     254,   254,   254,   254,   255,   255,   255,   255,   256,   256,
     257,   257,   258,   259,   260,   260,   261,   261,   261,   262,
     262,   262,   263,   263,   264,   264,   264,   264,   264,   264,
     265,   265,   265,   265,   265,   265,   266,   267,   267,   268,
     268,   268,   269,   269,   269,   269,   270,   270,   271,   271,
     271,   271,   271,   273,   274,   272,   275,   275,   275,   275,
     276,   276,   277,   277,   278,   279,   280,   281,   282,   283,
     284,   285,   286,   287,   288,   289,   290,   291,   292,   293,
     294,   295,   296,   296,   297,   298,   298,   299,   300,   301,
     302,   303,   304,   304,   305,   306,   307,   308,   309,   310,
     310,   311,   312,   313,   314,   315,   316,   317,   318,   319,
     320,   321,   322,   323,   324,   325,   326,   327,   328,   329,
     330,   330,   331,   332,   333,   334,   335,   336,   337,   338,
     339,   340,   341,   341,   342,   342,   343
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  parser::yyr2_[] =
  {
         0,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     4,     3,     3,
       3,     3,     3,     2,     3,     1,     3,     4,     2,     2,
       3,     1,     3,     4,     2,     2,     3,     1,     3,     4,
       2,     2,     3,     1,     3,     4,     2,     3,     4,     3,
       4,     4,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     2,     2,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     3,     3,     4,     7,     3,
       4,     2,     1,     4,     4,     2,     1,     7,     3,     1,
       1,     1,     1,     0,     5,     0,     8,     0,     8,     0,
      10,     0,     8,     2,     2,     1,     1,     4,     2,     3,
       1,     1,     1,     3,     3,     3,     3,     3,     2,     2,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       5,     1,     4,     4,     4,     2,     1,     9,     6,     5,
       7,     7,     2,     4,     3,     5,     3,     1,     2,     2,
       2,     1,     1,     1,     4,     3,     6,     3,     1,     5,
       3,     3,     4,     2,     2,     3,     1,     1,     2,     5,
       3,     1,     1,     2,     5,     3,     1,     1,     2,     5,
       3,     1,     1,     2,     5,     3,     6,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     2,
       4,     3,     5,     1,     3,     2,     2,     1,     2,     2,
       1,     4,     2,     1,     4,     2,     1,     4,     3,     5,
       9,     1,     5,     3,     5,     7,     9,     4,     2,     1,
       5,     7,     4,     4,     2,     1,     7,     9,     6,     1,
       1,     1,     1,     1,     0,     1,     1,     1,     2,     2,
       2,     5,     3,     6,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     5,     6,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     5,     6,     7,     5,     1,     3,     3,     4,
       2,     1,     5,     3,     4,     4,     6,     3,     5,     3,
       2,     5,     3,     6,     2,     5,     3,     6,     1,     1,
       1,     3,     3,     4,     2,     1,     5,     7,     9,     0,
       3,     3,     2,     5,     5,     6,     3,     7,     8,     5,
       5,     6,     3,     7,     8,     5,     6,     3,     1,     1,
       1,     1,     1,     3,     4,     6,     1,     2,     1,     1,
       1,     1,     1,     0,     0,     5,     2,     5,     3,     6,
       3,     1,     1,     1,     3,     3,     3,     1,     3,     3,
       3,     3,     1,     1,     1,     3,     3,     3,     3,     3,
       1,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     3,     3,     3,     3,     5,     3,
       3,     3,     3,     1,     3,     3,     3,     1,     1,     1,
       1,     1,     3,     1,     1,     1,     1,     3,     3,     3,
       3,     3,     1,     1,     3,     3,     3,     1,     1,     1,
       3,     3,     1,     3,     2,     2,     2
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const parser::yytname_[] =
  {
    "$end", "error", "$undefined", "AR", "AUTOCORR", "BAYESIAN_IRF",
  "BETA_PDF", "CALIB", "CALIB_VAR", "CHECK", "CONF_SIG", "CONSTANT",
  "CORR", "COVAR", "CUTOFF", "DATAFILE", "DR_ALGO", "DROP", "DSAMPLE",
  "DYNASAVE", "DYNATYPE", "END", "ENDVAL", "EQUAL", "ESTIMATION",
  "ESTIMATED_PARAMS", "ESTIMATED_PARAMS_BOUNDS", "ESTIMATED_PARAMS_INIT",
  "FILENAME", "FILTER_STEP_AHEAD", "FILTERED_VARS", "FIRST_OBS",
  "FLOAT_NUMBER", "FORECAST", "GAMMA_PDF", "GCC_COMPILER", "GRAPH",
  "HISTVAL", "HP_FILTER", "HP_NGRID", "INITVAL", "INT_NUMBER",
  "INV_GAMMA_PDF", "IRF", "KALMAN_ALGO", "KALMAN_TOL", "LAPLACE",
  "LCC_COMPILER", "LIK_ALGO", "LIK_INIT", "LINEAR", "LOAD_MH_FILE",
  "LOGLINEAR", "MH_DROP", "MH_INIT_SCALE", "MH_JSCALE", "MH_MODE",
  "MH_NBLOCKS", "MH_REPLIC", "MH_RECOVER", "MODE_CHECK", "MODE_COMPUTE",
  "MODE_FILE", "MODEL", "MODEL_COMPARISON", "MSHOCKS",
  "MODEL_COMPARISON_APPROXIMATION", "MODIFIEDHARMONICMEAN",
  "MOMENTS_VARENDO", "NAME", "NOBS", "NOCONSTANT", "NOCORR",
  "NODIAGNOSTIC", "NOFUNCTIONS", "NOGRAPH", "NOMOMENTS", "NOPRINT",
  "NORMAL_PDF", "OBSERVATION_TRENDS", "OLR", "OLR_INST", "OLR_BETA",
  "OPTIM", "OPTIM_WEIGHTS", "ORDER", "OSR", "OSR_PARAMS", "PARAMETERS",
  "PERIODS", "PLANNER_OBJECTIVE", "PREFILTER", "PRESAMPLE", "PRINT",
  "PRIOR_TRUNC", "PRIOR_ANALYSIS", "POSTERIOR_ANALYSIS", "QZ_CRITERIUM",
  "RELATIVE_IRF", "REPLIC", "RPLOT", "SHOCKS", "SIGMA_E", "SIMUL",
  "SIMUL_ALGO", "SIMUL_SEED", "SMOOTHER", "SOLVE_ALGO", "SPARSE_DLL",
  "STDERR", "STEADY", "STOCH_SIMUL", "TEX", "RAMSEY_POLICY",
  "PLANNER_DISCOUNT", "TEX_NAME", "UNIFORM_PDF", "UNIT_ROOT_VARS",
  "USE_DLL", "VALUES", "VAR", "VAREXO", "VAREXO_DET", "VAROBS",
  "XLS_SHEET", "XLS_RANGE", "COMMA", "MINUS", "PLUS", "DIVIDE", "TIMES",
  "UMINUS", "POWER", "EXP", "LOG", "LOG10", "SIN", "COS", "TAN", "ASIN",
  "ACOS", "ATAN", "SINH", "COSH", "TANH", "ASINH", "ACOSH", "ATANH",
  "SQRT", "';'", "'('", "')'", "'#'", "':'", "'['", "']'", "'''", "'.'",
  "'\\\\'", "$accept", "statement_list", "statement", "declaration",
  "dsample", "rplot", "var", "varexo", "varexo_det", "parameters",
  "var_list", "varexo_list", "varexo_det_list", "parameter_list",
  "periods", "cutoff", "init_param", "expression", "comma_expression",
  "initval", "initval_option", "endval", "initval_list", "initval_elem",
  "histval", "histval_list", "histval_elem", "model_sparse_options_list",
  "model_sparse_options", "model", "@1", "@2", "@3", "@4", "@5",
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
  "planner_objective", "@6", "@7", "ramsey_policy",
  "ramsey_policy_options_list", "ramsey_policy_options", "o_dr_algo",
  "o_solve_algo", "o_simul_algo", "o_linear", "o_order", "o_replic",
  "o_drop", "o_ar", "o_nocorr", "o_nofunctions", "o_nomoments", "o_irf",
  "o_hp_filter", "o_hp_ngrid", "o_periods", "o_cutoff", "o_simul",
  "o_simul_seed", "o_qz_criterium", "o_datafile", "o_nobs", "o_first_obs",
  "o_prefilter", "o_presample", "o_lik_algo", "o_lik_init", "o_nograph",
  "o_conf_sig", "o_mh_replic", "o_mh_drop", "o_mh_jscale", "o_optim",
  "o_mh_init_scale", "o_mode_file", "o_mode_compute", "o_mode_check",
  "o_prior_trunc", "o_mh_mode", "o_mh_nblcks", "o_load_mh_file",
  "o_loglinear", "o_nodiagnostic", "o_bayesian_irf", "o_tex", "o_forecast",
  "o_smoother", "o_moments_varendo", "o_filtered_vars", "o_relative_irf",
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
       160,     0,    -1,   161,    -1,   160,   161,    -1,   162,    -1,
     173,    -1,   174,    -1,   188,    -1,   178,    -1,   180,    -1,
     183,    -1,   175,    -1,   199,    -1,   200,    -1,   205,    -1,
     208,    -1,   211,    -1,   214,    -1,   217,    -1,   237,    -1,
     240,    -1,   243,    -1,   223,    -1,   232,    -1,   229,    -1,
     246,    -1,   247,    -1,   250,    -1,   163,    -1,   164,    -1,
     251,    -1,   253,    -1,   254,    -1,   259,    -1,   263,    -1,
     264,    -1,   265,    -1,   255,    -1,   258,    -1,   266,    -1,
     272,    -1,   275,    -1,   168,    -1,   165,    -1,   166,    -1,
     167,    -1,    18,    41,   149,    -1,    18,    41,    41,   149,
      -1,   100,   220,   149,    -1,   120,   169,   149,    -1,   121,
     170,   149,    -1,   122,   171,   149,    -1,    88,   172,   149,
      -1,   169,    69,    -1,   169,   126,    69,    -1,    69,    -1,
     169,    69,   115,    -1,   169,   126,    69,   115,    -1,    69,
     115,    -1,   170,    69,    -1,   170,   126,    69,    -1,    69,
      -1,   170,    69,   115,    -1,   170,   126,    69,   115,    -1,
      69,   115,    -1,   171,    69,    -1,   171,   126,    69,    -1,
      69,    -1,   171,    69,   115,    -1,   171,   126,    69,   115,
      -1,    69,   115,    -1,   172,    69,    -1,   172,   126,    69,
      -1,    69,    -1,   172,    69,   115,    -1,   172,   126,    69,
     115,    -1,    69,   115,    -1,    89,    41,   149,    -1,    89,
      23,    41,   149,    -1,    14,    32,   149,    -1,    14,    23,
      32,   149,    -1,    69,    23,   176,   149,    -1,   150,   176,
     151,    -1,    69,    -1,    32,    -1,    41,    -1,   176,   128,
     176,    -1,   176,   127,   176,    -1,   176,   129,   176,    -1,
     176,   130,   176,    -1,   176,   132,   176,    -1,   127,   176,
      -1,   128,   176,    -1,   133,   150,   176,   151,    -1,   134,
     150,   176,   151,    -1,   135,   150,   176,   151,    -1,   136,
     150,   176,   151,    -1,   137,   150,   176,   151,    -1,   138,
     150,   176,   151,    -1,   139,   150,   176,   151,    -1,   140,
     150,   176,   151,    -1,   141,   150,   176,   151,    -1,   148,
     150,   176,   151,    -1,    69,   150,   176,   151,    -1,    69,
     150,   177,   151,    -1,   176,   126,   176,    -1,   177,   126,
     176,    -1,    40,   149,   181,    21,    -1,    40,   150,   179,
     151,   149,   181,    21,    -1,    28,    23,    69,    -1,    22,
     149,   181,    21,    -1,   181,   182,    -1,   182,    -1,    69,
      23,   176,   149,    -1,    37,   149,   184,    21,    -1,   184,
     185,    -1,   185,    -1,    69,   150,   221,   151,    23,   176,
     149,    -1,   186,   126,   187,    -1,   187,    -1,    47,    -1,
      35,    -1,   293,    -1,    -1,    63,   149,   189,   194,    21,
      -1,    -1,    63,   150,   281,   151,   149,   190,   194,    21,
      -1,    -1,    63,   150,   118,   151,   149,   191,   194,    21,
      -1,    -1,    63,   150,   108,   126,   186,   151,   192,   149,
     194,    21,    -1,    -1,    63,   150,   108,   151,   193,   149,
     194,    21,    -1,   194,   195,    -1,   194,   197,    -1,   195,
      -1,   197,    -1,   196,    23,   196,   149,    -1,   196,   149,
      -1,   150,   196,   151,    -1,   198,    -1,    32,    -1,    41,
      -1,   196,   128,   196,    -1,   196,   127,   196,    -1,   196,
     129,   196,    -1,   196,   130,   196,    -1,   196,   132,   196,
      -1,   127,   196,    -1,   128,   196,    -1,   133,   150,   196,
     151,    -1,   134,   150,   196,   151,    -1,   135,   150,   196,
     151,    -1,   136,   150,   196,   151,    -1,   137,   150,   196,
     151,    -1,   138,   150,   196,   151,    -1,   139,   150,   196,
     151,    -1,   140,   150,   196,   151,    -1,   141,   150,   196,
     151,    -1,   148,   150,   196,   151,    -1,   152,    69,    23,
     196,   149,    -1,    69,    -1,    69,   150,   221,   151,    -1,
     101,   149,   201,    21,    -1,    65,   149,   201,    21,    -1,
     201,   202,    -1,   202,    -1,   120,    69,   149,    89,   203,
     149,   119,   204,   149,    -1,   120,    69,   149,   109,   176,
     149,    -1,   120,    69,    23,   176,   149,    -1,   120,    69,
     126,    69,    23,   176,   149,    -1,    12,    69,   126,    69,
      23,   176,   149,    -1,   203,    41,    -1,   203,    41,   153,
      41,    -1,   203,   126,    41,    -1,   203,   126,    41,   153,
      41,    -1,    41,   153,    41,    -1,    41,    -1,   204,   222,
      -1,   204,   221,    -1,   204,    69,    -1,   222,    -1,   221,
      -1,    69,    -1,   204,   150,   176,   151,    -1,   150,   176,
     151,    -1,   102,    23,   154,   206,   155,   149,    -1,   206,
     149,   207,    -1,   207,    -1,   207,   126,   150,   176,   151,
      -1,   207,   126,    32,    -1,   207,   126,    41,    -1,   207,
     150,   176,   151,    -1,   207,    32,    -1,   207,    41,    -1,
     150,   176,   151,    -1,    32,    -1,    41,    -1,   110,   149,
      -1,   110,   150,   209,   151,   149,    -1,   209,   126,   210,
      -1,   210,    -1,   279,    -1,     9,   149,    -1,     9,   150,
     212,   151,   149,    -1,   212,   126,   213,    -1,   213,    -1,
     279,    -1,   103,   149,    -1,   103,   150,   215,   151,   149,
      -1,   215,   126,   216,    -1,   216,    -1,   292,    -1,   111,
     149,    -1,   111,   150,   218,   151,   149,    -1,   111,   220,
     149,    -1,   111,   150,   218,   151,   220,   149,    -1,   218,
     126,   219,    -1,   219,    -1,   278,    -1,   279,    -1,   280,
      -1,   281,    -1,   282,    -1,   283,    -1,   284,    -1,   285,
      -1,   286,    -1,   287,    -1,   288,    -1,   289,    -1,   326,
      -1,   290,    -1,   291,    -1,   292,    -1,   293,    -1,   294,
      -1,   295,    -1,   296,    -1,   220,    69,    -1,   220,    69,
      23,    69,    -1,   220,   126,    69,    -1,   220,   126,    69,
      23,    69,    -1,    69,    -1,    69,    23,    69,    -1,   128,
      41,    -1,   127,    41,    -1,    41,    -1,   128,    32,    -1,
     127,    32,    -1,    32,    -1,    25,   149,   224,    21,    -1,
     224,   225,    -1,   225,    -1,   226,   126,   227,   149,    -1,
     109,    69,    -1,    69,    -1,    12,    69,   126,    69,    -1,
     235,   126,   228,    -1,   236,   126,   235,   126,   228,    -1,
     236,   126,   236,   126,   236,   126,   235,   126,   228,    -1,
     236,    -1,   236,   126,   236,   126,   236,    -1,   236,   126,
     236,    -1,   236,   126,   236,   126,   236,    -1,   236,   126,
     236,   126,   236,   126,   236,    -1,   236,   126,   236,   126,
     236,   126,   236,   126,   236,    -1,    27,   149,   230,    21,
      -1,   230,   231,    -1,   231,    -1,   109,    69,   126,   236,
     149,    -1,    12,    69,   126,    69,   126,   236,   149,    -1,
      69,   126,   236,   149,    -1,    26,   149,   233,    21,    -1,
     233,   234,    -1,   234,    -1,   109,    69,   126,   236,   126,
     236,   149,    -1,    12,    69,   126,    69,   126,   236,   126,
     236,   149,    -1,    69,   126,   236,   126,   236,   149,    -1,
       6,    -1,    34,    -1,    78,    -1,    42,    -1,   116,    -1,
      -1,    41,    -1,    32,    -1,    69,    -1,   127,    41,    -1,
     127,    32,    -1,    24,   149,    -1,    24,   150,   238,   151,
     149,    -1,    24,   220,   149,    -1,    24,   150,   238,   151,
     220,   149,    -1,   238,   126,   239,    -1,   239,    -1,   297,
      -1,   298,    -1,   299,    -1,   300,    -1,   301,    -1,   302,
      -1,   303,    -1,   304,    -1,   305,    -1,   306,    -1,   307,
      -1,   308,    -1,   309,    -1,   310,    -1,   311,    -1,   312,
      -1,   313,    -1,   314,    -1,   315,    -1,   316,    -1,   317,
      -1,   318,    -1,   319,    -1,   320,    -1,   289,    -1,   321,
      -1,   322,    -1,   323,    -1,   324,    -1,   325,    -1,   327,
      -1,   328,    -1,   333,    -1,   334,    -1,   335,    -1,   279,
      -1,   336,    -1,   337,    -1,   338,    -1,    95,   150,   241,
     151,   149,    -1,    95,   150,   241,   151,   220,   149,    -1,
     241,   126,   242,    -1,   242,    -1,   304,    -1,   305,    -1,
     314,    -1,   320,    -1,   289,    -1,   321,    -1,   322,    -1,
     323,    -1,   324,    -1,   325,    -1,   333,    -1,   334,    -1,
     335,    -1,    96,   150,   241,   151,   149,    -1,    96,   150,
     241,   151,   220,   149,    -1,   156,    69,   156,   126,   156,
      69,   156,    -1,   156,    69,   156,   126,   236,    -1,   244,
      -1,   245,   126,   244,    -1,   123,   220,   149,    -1,    79,
     149,   248,    21,    -1,   248,   249,    -1,   249,    -1,    69,
     150,   176,   151,   149,    -1,   117,   220,   149,    -1,    84,
     149,   252,    21,    -1,   252,    69,   176,   149,    -1,   252,
      69,   126,    69,   176,   149,    -1,    69,   176,   149,    -1,
      69,   126,    69,   176,   149,    -1,    87,   220,   149,    -1,
      86,   149,    -1,    86,   150,   257,   151,   149,    -1,    86,
     220,   149,    -1,    86,   150,   257,   151,   220,   149,    -1,
      80,   149,    -1,    80,   150,   257,   151,   149,    -1,    80,
     220,   149,    -1,    80,   150,   257,   151,   220,   149,    -1,
     329,    -1,   219,    -1,   256,    -1,   257,   126,   256,    -1,
      81,   220,   149,    -1,     8,   149,   260,    21,    -1,   260,
     261,    -1,   261,    -1,    69,   262,    23,   176,   149,    -1,
      69,   126,    69,   262,    23,   176,   149,    -1,     4,    69,
     150,    41,   151,   262,    23,   176,   149,    -1,    -1,   150,
      41,   151,    -1,   150,    32,   151,    -1,     7,   149,    -1,
       7,   150,    13,   151,   149,    -1,    20,   150,    69,   151,
     149,    -1,    20,   150,    69,   151,   220,   149,    -1,    20,
      69,   149,    -1,    20,   150,    69,   157,    69,   151,   149,
      -1,    20,   150,    69,   157,    69,   151,   220,   149,    -1,
      20,    69,   157,    69,   149,    -1,    19,   150,    69,   151,
     149,    -1,    19,   150,    69,   151,   220,   149,    -1,    19,
      69,   149,    -1,    19,   150,    69,   157,    69,   151,   149,
      -1,    19,   150,    69,   157,    69,   151,   220,   149,    -1,
      19,    69,   157,    69,   149,    -1,    64,   150,   267,   151,
     269,   149,    -1,   267,   126,   268,    -1,   268,    -1,   330,
      -1,   331,    -1,   332,    -1,   270,    -1,   269,   126,   270,
      -1,   270,   150,   236,   151,    -1,   269,   126,   270,   150,
     236,   151,    -1,   271,    -1,   270,   271,    -1,    69,    -1,
     158,    -1,   129,    -1,   153,    -1,   157,    -1,    -1,    -1,
      90,   273,   196,   274,   149,    -1,   113,   149,    -1,   113,
     150,   276,   151,   149,    -1,   113,   220,   149,    -1,   113,
     150,   276,   151,   220,   149,    -1,   276,   126,   277,    -1,
     277,    -1,   219,    -1,   339,    -1,    16,    23,    41,    -1,
     107,    23,    41,    -1,   104,    23,    41,    -1,    50,    -1,
      85,    23,    41,    -1,    99,    23,    41,    -1,    17,    23,
      41,    -1,     3,    23,    41,    -1,    72,    -1,    74,    -1,
      76,    -1,    43,    23,    41,    -1,    38,    23,    41,    -1,
      39,    23,    41,    -1,    89,    23,    41,    -1,    14,    23,
      32,    -1,   103,    -1,   105,    23,    41,    -1,    97,    23,
      41,    -1,    97,    23,    32,    -1,    15,    23,    69,    -1,
      70,    23,   343,    -1,    70,    23,    41,    -1,    31,    23,
      41,    -1,    91,    23,    41,    -1,    92,    23,    41,    -1,
      48,    23,    41,    -1,    49,    23,    41,    -1,    75,    -1,
      36,    -1,    10,    23,    32,    -1,    58,    23,    41,    -1,
      53,    23,    32,    -1,    55,    23,    32,    -1,    83,    23,
     150,   245,   151,    -1,    54,    23,    32,    -1,    54,    23,
      41,    -1,    62,    23,    69,    -1,    61,    23,    41,    -1,
      60,    -1,    94,    23,    32,    -1,    56,    23,    41,    -1,
      57,    23,    41,    -1,    51,    -1,    52,    -1,    73,    -1,
       5,    -1,   112,    -1,    33,    23,    41,    -1,   106,    -1,
      68,    -1,    30,    -1,    98,    -1,    44,    23,    41,    -1,
      45,    23,    41,    -1,    82,    23,   236,    -1,    66,    23,
      46,    -1,    66,    23,    67,    -1,    93,    -1,    77,    -1,
     124,    23,    69,    -1,   125,    23,   340,    -1,    29,    23,
     343,    -1,    11,    -1,    71,    -1,    59,    -1,   114,    23,
      32,    -1,    69,   153,    69,    -1,    41,    -1,    41,   153,
      41,    -1,   154,   341,    -1,   342,   341,    -1,   342,   155,
      -1
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
      80,    82,    84,    86,    88,    90,    92,    96,   101,   105,
     109,   113,   117,   121,   124,   128,   130,   134,   139,   142,
     145,   149,   151,   155,   160,   163,   166,   170,   172,   176,
     181,   184,   187,   191,   193,   197,   202,   205,   209,   214,
     218,   223,   228,   232,   234,   236,   238,   242,   246,   250,
     254,   258,   261,   264,   269,   274,   279,   284,   289,   294,
     299,   304,   309,   314,   319,   324,   328,   332,   337,   345,
     349,   354,   357,   359,   364,   369,   372,   374,   382,   386,
     388,   390,   392,   394,   395,   401,   402,   411,   412,   421,
     422,   433,   434,   443,   446,   449,   451,   453,   458,   461,
     465,   467,   469,   471,   475,   479,   483,   487,   491,   494,
     497,   502,   507,   512,   517,   522,   527,   532,   537,   542,
     547,   553,   555,   560,   565,   570,   573,   575,   585,   592,
     598,   606,   614,   617,   622,   626,   632,   636,   638,   641,
     644,   647,   649,   651,   653,   658,   662,   669,   673,   675,
     681,   685,   689,   694,   697,   700,   704,   706,   708,   711,
     717,   721,   723,   725,   728,   734,   738,   740,   742,   745,
     751,   755,   757,   759,   762,   768,   772,   779,   783,   785,
     787,   789,   791,   793,   795,   797,   799,   801,   803,   805,
     807,   809,   811,   813,   815,   817,   819,   821,   823,   825,
     828,   833,   837,   843,   845,   849,   852,   855,   857,   860,
     863,   865,   870,   873,   875,   880,   883,   885,   890,   894,
     900,   910,   912,   918,   922,   928,   936,   946,   951,   954,
     956,   962,   970,   975,   980,   983,   985,   993,  1003,  1010,
    1012,  1014,  1016,  1018,  1020,  1021,  1023,  1025,  1027,  1030,
    1033,  1036,  1042,  1046,  1053,  1057,  1059,  1061,  1063,  1065,
    1067,  1069,  1071,  1073,  1075,  1077,  1079,  1081,  1083,  1085,
    1087,  1089,  1091,  1093,  1095,  1097,  1099,  1101,  1103,  1105,
    1107,  1109,  1111,  1113,  1115,  1117,  1119,  1121,  1123,  1125,
    1127,  1129,  1131,  1133,  1135,  1137,  1143,  1150,  1154,  1156,
    1158,  1160,  1162,  1164,  1166,  1168,  1170,  1172,  1174,  1176,
    1178,  1180,  1182,  1188,  1195,  1203,  1209,  1211,  1215,  1219,
    1224,  1227,  1229,  1235,  1239,  1244,  1249,  1256,  1260,  1266,
    1270,  1273,  1279,  1283,  1290,  1293,  1299,  1303,  1310,  1312,
    1314,  1316,  1320,  1324,  1329,  1332,  1334,  1340,  1348,  1358,
    1359,  1363,  1367,  1370,  1376,  1382,  1389,  1393,  1401,  1410,
    1416,  1422,  1429,  1433,  1441,  1450,  1456,  1463,  1467,  1469,
    1471,  1473,  1475,  1477,  1481,  1486,  1493,  1495,  1498,  1500,
    1502,  1504,  1506,  1508,  1509,  1510,  1516,  1519,  1525,  1529,
    1536,  1540,  1542,  1544,  1546,  1550,  1554,  1558,  1560,  1564,
    1568,  1572,  1576,  1578,  1580,  1582,  1586,  1590,  1594,  1598,
    1602,  1604,  1608,  1612,  1616,  1620,  1624,  1628,  1632,  1636,
    1640,  1644,  1648,  1650,  1652,  1656,  1660,  1664,  1668,  1674,
    1678,  1682,  1686,  1690,  1692,  1696,  1700,  1704,  1706,  1708,
    1710,  1712,  1714,  1718,  1720,  1722,  1724,  1726,  1730,  1734,
    1738,  1742,  1746,  1748,  1750,  1754,  1758,  1762,  1764,  1766,
    1768,  1772,  1776,  1778,  1782,  1785,  1788
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    89,    89,    90,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,   108,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   122,   123,   124,   125,   126,   127,   128,   129,
     130,   131,   136,   137,   138,   139,   143,   144,   147,   150,
     154,   158,   162,   166,   168,   170,   172,   174,   176,   181,
     183,   185,   187,   189,   191,   196,   198,   200,   202,   204,
     206,   211,   213,   215,   217,   219,   221,   226,   230,   237,
     241,   248,   253,   255,   257,   259,   261,   263,   265,   267,
     269,   271,   273,   275,   277,   279,   281,   283,   285,   287,
     289,   291,   293,   295,   297,   302,   304,   308,   310,   315,
     319,   324,   325,   329,   334,   339,   340,   344,   348,   349,
     353,   354,   355,   359,   359,   360,   360,   362,   362,   364,
     364,   366,   366,   371,   372,   373,   374,   378,   380,   385,
     386,   387,   389,   391,   393,   395,   397,   399,   401,   403,
     405,   407,   409,   411,   413,   415,   417,   419,   421,   423,
     427,   431,   433,   438,   442,   446,   447,   451,   453,   455,
     457,   459,   464,   466,   468,   470,   472,   474,   480,   482,
     484,   486,   488,   490,   492,   494,   499,   504,   506,   511,
     513,   515,   517,   519,   521,   523,   525,   527,   532,   536,
     540,   541,   544,   548,   550,   554,   555,   558,   562,   564,
     568,   569,   572,   576,   578,   580,   582,   586,   587,   590,
     591,   592,   593,   594,   595,   596,   597,   598,   599,   600,
     601,   602,   603,   604,   605,   606,   607,   608,   609,   613,
     615,   617,   619,   621,   623,   628,   630,   632,   637,   639,
     641,   646,   651,   653,   658,   662,   667,   672,   682,   687,
     693,   703,   708,   719,   725,   733,   743,   757,   761,   763,
     767,   774,   783,   792,   796,   798,   802,   811,   822,   834,
     836,   838,   840,   842,   847,   848,   849,   850,   851,   853,
     860,   862,   864,   866,   871,   872,   875,   876,   877,   878,
     879,   880,   881,   882,   883,   884,   885,   886,   887,   888,
     889,   890,   891,   892,   893,   894,   895,   896,   897,   898,
     899,   900,   901,   902,   903,   904,   905,   906,   907,   908,
     909,   910,   911,   912,   913,   917,   919,   924,   925,   929,
     930,   931,   932,   933,   934,   935,   936,   937,   938,   939,
     940,   941,   945,   947,   952,   953,   957,   958,   962,   967,
     972,   973,   976,   980,   983,   987,   989,   991,   993,   997,
    1000,  1001,  1002,  1003,  1006,  1007,  1008,  1009,  1012,  1013,
    1016,  1017,  1020,  1023,  1027,  1028,  1031,  1032,  1033,  1036,
    1037,  1038,  1041,  1042,  1045,  1046,  1047,  1048,  1049,  1050,
    1052,  1053,  1054,  1055,  1056,  1057,  1059,  1063,  1064,  1067,
    1068,  1069,  1072,  1073,  1074,  1075,  1078,  1079,  1082,  1083,
    1084,  1085,  1086,  1089,  1089,  1089,  1092,  1094,  1096,  1098,
    1103,  1104,  1107,  1108,  1111,  1112,  1113,  1114,  1115,  1116,
    1117,  1118,  1119,  1120,  1121,  1122,  1123,  1124,  1125,  1126,
    1127,  1128,  1129,  1130,  1132,  1133,  1134,  1136,  1137,  1138,
    1139,  1140,  1141,  1142,  1143,  1144,  1145,  1146,  1147,  1148,
    1149,  1150,  1151,  1152,  1153,  1154,  1155,  1156,  1157,  1158,
    1159,  1160,  1161,  1162,  1163,  1164,  1165,  1166,  1167,  1168,
    1170,  1172,  1175,  1176,  1177,  1178,  1179,  1180,  1181,  1182,
    1183,  1185,  1193,  1194,  1198,  1199,  1208
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
       2,     2,     2,     2,     2,   152,     2,     2,     2,   156,
     150,   151,     2,     2,     2,     2,   157,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   153,   149,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   154,   158,   155,     2,     2,     2,     2,     2,     2,
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
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1473;
  const int parser::yynnts_ = 185;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 152;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 159;

  const unsigned int parser::yyuser_token_number_max_ = 403;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1210 "DynareBison.yy"


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

