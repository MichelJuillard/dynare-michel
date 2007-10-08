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
#line 32 "DynareBison.yy"

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
    #line 18 "DynareBison.yy"
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
	  case 48:
#line 152 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 49:
#line 154 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 50:
#line 157 "DynareBison.yy"
    { driver.rplot(); ;}
    break;

  case 55:
#line 168 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 56:
#line 170 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 57:
#line 172 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 58:
#line 174 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 59:
#line 176 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 60:
#line 178 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 61:
#line 182 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 62:
#line 184 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 63:
#line 186 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 64:
#line 188 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 65:
#line 190 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 66:
#line 192 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 67:
#line 196 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 68:
#line 198 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 69:
#line 200 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 70:
#line 202 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 71:
#line 204 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 72:
#line 206 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 73:
#line 210 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 74:
#line 212 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 75:
#line 214 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 76:
#line 216 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 77:
#line 218 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 78:
#line 220 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 79:
#line 224 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 80:
#line 226 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 81:
#line 230 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 82:
#line 232 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 83:
#line 236 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 84:
#line 238 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 85:
#line 242 "DynareBison.yy"
    { driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 86:
#line 245 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 87:
#line 247 "DynareBison.yy"
    { (yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 88:
#line 249 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 89:
#line 251 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 90:
#line 253 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 91:
#line 255 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 92:
#line 257 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 93:
#line 259 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 94:
#line 261 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 95:
#line 263 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 96:
#line 265 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 97:
#line 267 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 98:
#line 269 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 99:
#line 271 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 100:
#line 273 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 101:
#line 275 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 102:
#line 277 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 103:
#line 279 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 104:
#line 281 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 105:
#line 283 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 106:
#line 285 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 107:
#line 287 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 108:
#line 289 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 109:
#line 291 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 110:
#line 293 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 111:
#line 295 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 112:
#line 297 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 113:
#line 299 "DynareBison.yy"
    { (yyval.node_val) = driver.add_dummy((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 114:
#line 301 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 115:
#line 303 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 116:
#line 305 "DynareBison.yy"
    { (yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val)); ;}
    break;

  case 117:
#line 309 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 118:
#line 311 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 119:
#line 315 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 120:
#line 317 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 121:
#line 320 "DynareBison.yy"
    { driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 122:
#line 322 "DynareBison.yy"
    { driver.end_endval(); ;}
    break;

  case 125:
#line 328 "DynareBison.yy"
    { driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 126:
#line 330 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 129:
#line 336 "DynareBison.yy"
    { driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 132:
#line 343 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 133:
#line 345 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 134:
#line 347 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 137:
#line 352 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 138:
#line 353 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 139:
#line 354 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 140:
#line 355 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 141:
#line 356 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 142:
#line 357 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 143:
#line 359 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 144:
#line 360 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 145:
#line 361 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 146:
#line 362 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 151:
#line 372 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 152:
#line 374 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val)); ;}
    break;

  case 153:
#line 378 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 155:
#line 381 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 156:
#line 383 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 157:
#line 385 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 158:
#line 387 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 159:
#line 389 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 160:
#line 391 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 161:
#line 393 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 162:
#line 395 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 163:
#line 397 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 164:
#line 399 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 165:
#line 401 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 166:
#line 403 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 167:
#line 405 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 168:
#line 407 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 169:
#line 409 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 170:
#line 411 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 171:
#line 413 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 172:
#line 415 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 173:
#line 417 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 174:
#line 419 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 175:
#line 421 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 176:
#line 423 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 177:
#line 425 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 178:
#line 427 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 179:
#line 429 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 180:
#line 431 "DynareBison.yy"
    { (yyval.node_val) = driver.add_dummy((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 181:
#line 433 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 182:
#line 435 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 183:
#line 439 "DynareBison.yy"
    { driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 184:
#line 442 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 185:
#line 444 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 186:
#line 447 "DynareBison.yy"
    { driver.end_shocks(); ;}
    break;

  case 187:
#line 449 "DynareBison.yy"
    { driver.end_mshocks(); ;}
    break;

  case 190:
#line 456 "DynareBison.yy"
    { driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val)); ;}
    break;

  case 191:
#line 458 "DynareBison.yy"
    { driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 192:
#line 460 "DynareBison.yy"
    { driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 193:
#line 462 "DynareBison.yy"
    { driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 194:
#line 464 "DynareBison.yy"
    { driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 195:
#line 468 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 196:
#line 470 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 197:
#line 472 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 198:
#line 474 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 199:
#line 476 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 200:
#line 478 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 201:
#line 482 "DynareBison.yy"
    { driver.add_value((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 202:
#line 484 "DynareBison.yy"
    { driver.add_value((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 203:
#line 487 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 204:
#line 490 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 205:
#line 492 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 206:
#line 496 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 207:
#line 498 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 208:
#line 500 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 209:
#line 502 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 210:
#line 504 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 211:
#line 506 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 212:
#line 508 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 213:
#line 510 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 214:
#line 512 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 215:
#line 516 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 216:
#line 518 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 220:
#line 528 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 221:
#line 530 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 225:
#line 540 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 226:
#line 542 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 231:
#line 554 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 232:
#line 556 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 233:
#line 558 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 234:
#line 560 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 260:
#line 593 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 261:
#line 595 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 262:
#line 597 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 263:
#line 599 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 264:
#line 601 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 265:
#line 603 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 266:
#line 607 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 267:
#line 609 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 268:
#line 611 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 269:
#line 615 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 270:
#line 617 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 271:
#line 619 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 272:
#line 622 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 273:
#line 625 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 274:
#line 627 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 276:
#line 633 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 277:
#line 639 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 278:
#line 645 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 279:
#line 655 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 280:
#line 660 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 281:
#line 667 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(9) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(9) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(9) - (5)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(9) - (7)].string_val);
                    delete (yysemantic_stack_[(9) - (1)].string_val);
                    delete (yysemantic_stack_[(9) - (3)].string_val);
                    delete (yysemantic_stack_[(9) - (5)].string_val);
                    delete (yysemantic_stack_[(9) - (7)].string_val);
                  ;}
    break;

  case 282:
#line 678 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 283:
#line 683 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 284:
#line 694 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 285:
#line 701 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 286:
#line 710 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(7) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(7) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(7) - (5)].string_val);
                    driver.estim_params.p4 = *(yysemantic_stack_[(7) - (7)].string_val);
                    delete (yysemantic_stack_[(7) - (1)].string_val);
                    delete (yysemantic_stack_[(7) - (3)].string_val);
                    delete (yysemantic_stack_[(7) - (5)].string_val);
                    delete (yysemantic_stack_[(7) - (7)].string_val);
                  ;}
    break;

  case 287:
#line 721 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(9) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(9) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(9) - (5)].string_val);
                    driver.estim_params.p4 = *(yysemantic_stack_[(9) - (7)].string_val);
                    driver.estim_params.jscale = *(yysemantic_stack_[(9) - (9)].string_val);
                    delete (yysemantic_stack_[(9) - (1)].string_val);
                    delete (yysemantic_stack_[(9) - (3)].string_val);
                    delete (yysemantic_stack_[(9) - (5)].string_val);
                    delete (yysemantic_stack_[(9) - (7)].string_val);
                    delete (yysemantic_stack_[(9) - (9)].string_val);
                  ;}
    break;

  case 288:
#line 736 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 289:
#line 739 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 290:
#line 741 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 291:
#line 745 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 292:
#line 753 "DynareBison.yy"
    {
                        driver.estim_params.type = 3;
                        driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
                        driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(7) - (6)].string_val);
                        delete (yysemantic_stack_[(7) - (2)].string_val);
                        delete (yysemantic_stack_[(7) - (4)].string_val);
                        delete (yysemantic_stack_[(7) - (6)].string_val);
                      ;}
    break;

  case 293:
#line 763 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 294:
#line 773 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 295:
#line 776 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 296:
#line 778 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 297:
#line 782 "DynareBison.yy"
    {
                          driver.estim_params.type = 1;
                          driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
                          driver.estim_params.low_bound = *(yysemantic_stack_[(7) - (4)].string_val);
                          driver.estim_params.up_bound = *(yysemantic_stack_[(7) - (6)].string_val);
                          delete (yysemantic_stack_[(7) - (2)].string_val);
                          delete (yysemantic_stack_[(7) - (4)].string_val);
                          delete (yysemantic_stack_[(7) - (6)].string_val);
                        ;}
    break;

  case 298:
#line 792 "DynareBison.yy"
    {
                          driver.estim_params.type = 3;
                          driver.estim_params.name = *(yysemantic_stack_[(9) - (2)].string_val);
                          driver.estim_params.name2 = *(yysemantic_stack_[(9) - (4)].string_val);
                          driver.estim_params.low_bound = *(yysemantic_stack_[(9) - (6)].string_val);
                          driver.estim_params.up_bound = *(yysemantic_stack_[(9) - (8)].string_val);
                          delete (yysemantic_stack_[(9) - (2)].string_val);
                          delete (yysemantic_stack_[(9) - (4)].string_val);
                          delete (yysemantic_stack_[(9) - (6)].string_val);
                          delete (yysemantic_stack_[(9) - (8)].string_val);
                        ;}
    break;

  case 299:
#line 804 "DynareBison.yy"
    {
                          driver.estim_params.type = 2;
                          driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
                          driver.estim_params.low_bound = *(yysemantic_stack_[(6) - (3)].string_val);
                          driver.estim_params.up_bound = *(yysemantic_stack_[(6) - (5)].string_val);
                          delete (yysemantic_stack_[(6) - (1)].string_val);
                          delete (yysemantic_stack_[(6) - (3)].string_val);
                          delete (yysemantic_stack_[(6) - (5)].string_val);
                        ;}
    break;

  case 300:
#line 816 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 301:
#line 818 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 302:
#line 820 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 303:
#line 822 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 304:
#line 824 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 305:
#line 827 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 309:
#line 832 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 310:
#line 834 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 311:
#line 838 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 312:
#line 840 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 313:
#line 842 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 314:
#line 844 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 356:
#line 893 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 357:
#line 895 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 373:
#line 918 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 374:
#line 920 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 375:
#line 924 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 376:
#line 926 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 379:
#line 933 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 380:
#line 935 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 383:
#line 941 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 384:
#line 943 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 385:
#line 945 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 386:
#line 948 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 387:
#line 950 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 388:
#line 952 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 389:
#line 954 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 390:
#line 957 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 391:
#line 960 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 392:
#line 962 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 393:
#line 964 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 394:
#line 966 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 395:
#line 969 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 398:
#line 976 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 399:
#line 978 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 400:
#line 980 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 401:
#line 983 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 402:
#line 985 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 403:
#line 987 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 404:
#line 991 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 405:
#line 993 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 406:
#line 997 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 407:
#line 999 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 408:
#line 1001 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 409:
#line 1003 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 410:
#line 1005 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 411:
#line 1007 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 412:
#line 1011 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 413:
#line 1013 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 414:
#line 1015 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 415:
#line 1017 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 416:
#line 1019 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 417:
#line 1021 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 418:
#line 1025 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 424:
#line 1037 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 425:
#line 1039 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 426:
#line 1041 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 427:
#line 1043 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 428:
#line 1047 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 429:
#line 1049 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 431:
#line 1054 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 432:
#line 1056 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 433:
#line 1058 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 434:
#line 1060 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 435:
#line 1063 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 436:
#line 1064 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 438:
#line 1067 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 439:
#line 1069 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 440:
#line 1071 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 441:
#line 1073 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 465:
#line 1110 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 466:
#line 1112 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 473:
#line 1126 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 474:
#line 1128 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 475:
#line 1132 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 476:
#line 1134 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 508:
#line 1174 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 509:
#line 1175 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 510:
#line 1176 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 511:
#line 1177 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 512:
#line 1178 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 513:
#line 1179 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 514:
#line 1180 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 515:
#line 1181 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 516:
#line 1182 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 517:
#line 1183 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 518:
#line 1184 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 519:
#line 1185 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 520:
#line 1186 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 521:
#line 1187 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 522:
#line 1189 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 523:
#line 1190 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 524:
#line 1191 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 525:
#line 1192 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 526:
#line 1193 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 527:
#line 1194 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 528:
#line 1195 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1197 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1199 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1201 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1202 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1203 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 534:
#line 1204 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 535:
#line 1205 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 536:
#line 1207 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 537:
#line 1209 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 538:
#line 1211 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1212 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 540:
#line 1213 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1214 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 543:
#line 1216 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 544:
#line 1217 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 545:
#line 1218 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 546:
#line 1219 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 547:
#line 1220 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 548:
#line 1221 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 549:
#line 1222 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 550:
#line 1223 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 551:
#line 1224 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 552:
#line 1225 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 553:
#line 1226 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 554:
#line 1227 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 555:
#line 1228 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 556:
#line 1229 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 557:
#line 1230 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 558:
#line 1231 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 559:
#line 1232 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 560:
#line 1233 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 561:
#line 1234 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 562:
#line 1236 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 563:
#line 1238 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 564:
#line 1240 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 565:
#line 1241 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 566:
#line 1242 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1243 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 568:
#line 1244 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 569:
#line 1245 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 570:
#line 1246 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 571:
#line 1247 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 572:
#line 1248 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1250 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1251 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 575:
#line 1252 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 576:
#line 1253 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 577:
#line 1254 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1255 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 579:
#line 1256 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1257 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 581:
#line 1259 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 582:
#line 1260 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1261 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1262 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1263 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1264 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1265 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1266 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1267 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1268 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1269 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 592:
#line 1270 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 593:
#line 1271 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 594:
#line 1272 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 595:
#line 1273 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 596:
#line 1274 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 597:
#line 1275 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 598:
#line 1276 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 599:
#line 1277 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 600:
#line 1279 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 601:
#line 1280 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 602:
#line 1286 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 603:
#line 1287 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 604:
#line 1291 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 605:
#line 1292 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 606:
#line 1293 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 607:
#line 1294 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 608:
#line 1298 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 610:
#line 1307 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); 
               ;}
    break;

  case 611:
#line 1316 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 612:
#line 1318 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 613:
#line 1326 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2419 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1112;
  const short int
  parser::yypact_[] =
  {
      1257,    14,    30,    84,   -59,   250,   109,    70,   -12,    40,
     -52,     5,   -39,    86,   102,   104,   263,   110,   310,   -94,
     152,   124,   154,   167,    46,   192,   332,   103, -1112,   285,
     337,   192,   291,   471,   365,   369,    51,    53,   192,   434,
     458,   460,   192,   410,  1138, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
     380,   716,   385,  1473, -1112,   592,   127, -1112,   498,   569,
     419,    48,   367,   538,   368,   539,   544,   595, -1112,  1373,
      79,    55,   356,   413,   547,   544,   598,   588,   440, -1112,
      72,   427,    50,   951,   562,   563, -1112,  1557,    92,   173,
     527,   174,   604,   457,   977,   363,   363,   176,    50,   454,
   -1112,    69, -1112,   498, -1112,  1557,   200, -1112,  1477,   222,
     223,   531,   224,   533,   232,   536,   233,   234, -1112,  2046,
   -1112, -1112, -1112,   630, -1112,   636,   637,   638,   640,   641,
   -1112,   643,   644,   646, -1112,   647,   659,   661,   662, -1112,
     529,   495, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,   665,
     666,   667, -1112,   566,   503, -1112, -1112, -1112,   509,   631,
     -21,   382, -1112,   679,   -30, -1112, -1112,   517, -1112,   523,
   -1112, -1112,   649,   -81, -1112,   650,   218,   700,    71, -1112,
     652, -1112,   702, -1112, -1112,   704,   705,   706,   710,   712,
   -1112, -1112,   713,   727,   728,   731,   733,   737, -1112, -1112,
     738,   740, -1112, -1112, -1112,   741,   742, -1112, -1112,   175,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
     743,   681, -1112,   695, -1112,   698,   399, -1112,   645,   711,
     658,   715,   453, -1112,   721,   668,   722,   456, -1112,   606,
      94, -1112,   397,   775,   608,   611, -1112,   730, -1112,   217,
     612,   614,   784, -1112, -1112,   241, -1112, -1112, -1112, -1112,
     734,   739,   113, -1112,   620, -1112, -1112,   622,   623,   624,
     951,   951,   626,   627,   628,   629,   639,   642,   648,   653,
     654,   656,   951,   892,   657,   414, -1112,   783,   415,   797,
     798,   804,   805,   807,   810, -1112, -1112, -1112,   816,   825,
     827, -1112,   828, -1112,   829,   830,   663, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112,   744,   759, -1112,   669, -1112,   671, -1112,
   -1112,   673,   674,   677,   977,   977,   678,   691,   692,   693,
     694,   708,   723,   725,   726,   746,   977,  2135, -1112,   246,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112,   247, -1112,   117,    35,   260, -1112,
   -1112, -1112,   287, -1112, -1112,   290, -1112, -1112,   833, -1112,
     292, -1112, -1112, -1112, -1112, -1112,   753,   785, -1112, -1112,
     772,   787, -1112, -1112,   780,   840, -1112, -1112,   894,   898,
     900,   910,   911,   914,   915,   917,   919,   924,   925,   935,
     937,   938,   940,   942,   943,   944,   963,   964,   965,   967,
     968,   970,   971,   973,   974,   750,   873, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112,   431,   326,   431,   960,   326,   961,
     931,   966,    23,   969,   975,   932,   934,   716,   976,   978,
     431,   993,  1473,   994,   819,   821,   972,   520,   992, -1112,
   -1112,   995,   498,   848, -1112, -1112,   849,    63,   985,   851,
      73,   987,   951, -1112, -1112, -1112,   847,  1001,  1002,  1004,
    1019,  1020,   431,   431,   431,  1023,  1024,  1025,  1026,  1005,
     888,   431,  1373,    74,  1007,  1057,   955, -1112, -1112, -1112,
     313,   957,   375,   958, -1112, -1112,   962,   375,   979, -1112,
   -1112,    27, -1112, -1112, -1112,  1016,   901, -1112,  1029,    41,
   -1112,   150, -1112,   589, -1112,   902,   913,    61,   427,   155,
     981,    29, -1112, -1112,   951,   951,   951,   951,   988,   493,
     951,   951,   951,   951,   951,   951,   951,   951,   951,   951,
     543,   951,   951,   951,   951,   951,   951,   951,   951,   951,
     951,   951, -1112,   951, -1112, -1112,  1031,  1079, -1112,   808,
    1063,   431,  1064,  1069,  1070,  1073,  1074,  1088,   431,  1089,
    1090,  1091,    76, -1112,  1021, -1112,   977,   977,   977,    27,
     996,   501,   977,   977,   977,   977,   977,   977,   977,   977,
     977,   977,   771,   977,   977,   977,   977,   977,   977,   977,
     977,   977,   977,   977,   949,   363,    77,    80, -1112, -1112,
   -1112,   951,   -91,    31,    69,   950,   498,   953,  1557,    81,
     431,  1477,    82, -1112,  1027, -1112,  1032, -1112,  1033,  1099,
    1106,  1109,  1111,  1112,  1113,  1115,  1118,  1120,  1136,  1142,
    1152,  1153,  1154,  1155,   431,   431,  1156,   847,   431,   431,
    1157,  1158,   431,  1159,   431,   431,   990,  2046, -1112, -1112,
   -1112, -1112,  1169,  1170, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112,  1162,    12, -1112, -1112, -1112, -1112,  1037, -1112,
   -1112,  1035, -1112, -1112, -1112, -1112,  1042, -1112,  1177,  1030,
    1041,  1045,   951, -1112, -1112, -1112, -1112, -1112,   237,  1046,
   -1112, -1112,   270,  1047,  1827, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,  1048,
   -1112, -1112, -1112,   314, -1112,  1150,  1161, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112,   521,  1056,  1119,  1121,  1175,
    1124,   375,  1183,  1068,   375, -1112,  1217,  1225,  1076, -1112,
     544,  1246, -1112, -1112, -1112,   977, -1112, -1112, -1112,  1248,
   -1112,   302, -1112, -1112, -1112,  1083, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112,   -34,   -18, -1112,  1201,
     951,  1203,     4,   842,  2032,  2047,  2148,   304,   918,  1052,
    1171,  1198,  1396,  1476,  1489,  1541,  1554,  1567, -1112,   757,
     757,   561,   561,   561,   561,   493,   493,   988,   988,  1140,
    1580,   951, -1112,  1206,  1840, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,   316, -1112,
    1593,  2108,  2122,  1094,  1606,  1619,  1632,  1645,  1658,  1671,
    1684,  1697,  1710,  1723, -1112,  1036,  1036,   573,   573,   573,
     573,   501,   501,   996,   996,  1140, -1112, -1112, -1112,   318,
   -1112,   322,  1736,    35,  1097, -1112, -1112,    49,   951, -1112,
   -1112, -1112, -1112, -1112, -1112,   328, -1112, -1112, -1112,   329,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112,  1095, -1112, -1112, -1112,  1215, -1112,
   -1112,  1098,  1267, -1112, -1112,  1853, -1112,    88, -1112,   126,
   -1112,  1219, -1112,   306, -1112, -1112, -1112, -1112, -1112, -1112,
     375,   313,  1165,   375,  1166,  1185, -1112,  1107, -1112, -1112,
    1275,   429,   977,  1866,   431,   589, -1112,   730,   730,   730,
     155, -1112,   375, -1112,  1279,  1885,  1290,  1273,   951, -1112,
     951,   951,   951, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112,  1126,  1898,   951, -1112, -1112, -1112,
     977,   977, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112,    31, -1112, -1112, -1112,
     951,  1749, -1112, -1112,  1276, -1112,  1030,   951, -1112, -1112,
     330, -1112,   333,  1122,  1048, -1112, -1112,  1189,  1192,  1193,
     375,  1132,   375,   375, -1112,   951, -1112,  1911, -1112, -1112,
   -1112,  1133,    67,   235,   596,     3,  1146,   951, -1112,   951,
    1148,    43,  1924,  1762,  1775,  2148, -1112, -1112,  1943,  1788,
    1801,  1814, -1112, -1112,  1316,  1956, -1112, -1112,  1216, -1112,
     375,   375,   375,  1222, -1112,  1163,  1167,  1969, -1112,   730,
   -1112, -1112, -1112,   375, -1112,  1982,  2001,  1309,  1164,  1310,
    1235, -1112, -1112, -1112, -1112, -1112, -1112, -1112,   951, -1112,
      33,  1236, -1112,  1237,   375, -1112, -1112, -1112,   633,  1174,
   -1112, -1112, -1112,  1325,  1179,   951,  2014,  1298, -1112,   375,
      75,  1184, -1112, -1112, -1112,  1333,  2148,   927, -1112,  1180,
    1250,  1258, -1112, -1112, -1112,  2148, -1112,   375,   375,  1259,
   -1112,   375, -1112
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   435,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     4,    29,    30,    45,
      46,    47,    44,     5,     6,     7,    12,     9,    10,    11,
       8,    13,    14,    15,    16,    17,    18,    19,    23,    25,
      24,    20,    21,    22,    26,    27,    28,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,   404,     0,     0,   220,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   264,   311,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   137,
       0,     0,     0,     0,     0,     0,   391,     0,     0,     0,
      75,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     225,     0,   215,     0,   231,     0,     0,   438,     0,     0,
       0,    57,     0,    63,     0,    69,     0,     0,   475,     0,
       1,     3,   465,     0,   578,     0,     0,     0,     0,     0,
     569,     0,     0,     0,   570,     0,     0,     0,     0,   453,
     464,     0,   454,   459,   457,   460,   458,   455,   456,   461,
     462,   446,   447,   448,   449,   450,   451,   452,   473,     0,
       0,     0,   467,   472,     0,   469,   468,   470,     0,     0,
     401,     0,   397,     0,     0,   223,   224,     0,    81,     0,
      48,   414,     0,     0,   408,     0,     0,     0,     0,   124,
       0,   553,     0,   558,   537,     0,     0,     0,     0,     0,
     550,   551,     0,     0,     0,     0,     0,     0,   571,   546,
       0,     0,   557,   552,   536,     0,     0,   556,   554,     0,
     316,   352,   341,   317,   318,   319,   320,   321,   322,   323,
     324,   325,   326,   327,   328,   329,   330,   331,   332,   333,
     334,   335,   336,   337,   338,   339,   340,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   353,   354,   355,
     260,     0,   313,     0,   277,     0,     0,   274,     0,     0,
       0,     0,     0,   296,     0,     0,     0,     0,   290,     0,
       0,   128,     0,     0,     0,     0,    83,     0,   511,     0,
       0,     0,     0,   565,   564,     0,   420,   421,   422,   423,
       0,     0,     0,   189,     0,    88,    89,     0,     0,    87,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   382,     0,     0,     0,
       0,     0,     0,     0,     0,   516,   517,   518,     0,     0,
       0,   559,     0,   525,     0,     0,     0,   237,   238,   239,
     240,   241,   242,   243,   244,   245,   246,   247,   249,   251,
     252,   253,   254,   255,   256,   257,   248,   250,   258,   259,
     393,   390,    78,    73,     0,    54,     0,    79,     0,   155,
     156,     0,     0,   184,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   436,   154,     0,
     359,   364,   360,   361,   362,   363,   365,   366,   367,   368,
     369,   370,   371,   372,     0,    50,     0,     0,     0,   228,
     229,   230,     0,   218,   219,     0,   236,   233,     0,   444,
       0,   443,   445,   440,   384,    60,    55,     0,    51,    66,
      61,     0,    52,    72,    67,     0,    53,   379,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   478,   479,   480,   481,
     482,   483,   484,   485,   486,   487,   488,   489,   490,   491,
     492,   493,   494,   495,   496,   505,   497,   498,   499,   500,
     501,   502,   503,   504,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   395,
     396,     0,     0,     0,    82,    49,     0,     0,     0,     0,
       0,     0,     0,   122,   123,   265,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   262,     0,   276,   272,   273,
     305,     0,   305,     0,   294,   295,     0,   305,     0,   288,
     289,     0,   126,   127,   119,     0,     0,    84,     0,     0,
     149,     0,   150,     0,   145,     0,     0,     0,     0,     0,
       0,     0,   187,   188,     0,     0,     0,     0,   101,   102,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    85,     0,   380,   381,     0,     0,   385,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    76,    74,    80,     0,     0,     0,     0,
     168,   169,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   186,   213,
     214,     0,     0,   205,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    58,    56,    64,    62,    70,    68,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   507,   506,
     574,   271,     0,     0,   575,   576,   577,   573,   579,   528,
     531,   530,     0,     0,   529,   532,   533,   566,     0,   567,
     463,     0,   580,   538,   555,   471,     0,   405,     0,   401,
       0,     0,     0,   509,   222,   221,   417,   412,     0,     0,
     411,   406,     0,     0,     0,   568,   519,   560,   561,   534,
     535,   540,   543,   541,   548,   549,   539,   545,   544,     0,
     547,   315,   312,     0,   261,     0,     0,   300,   307,   301,
     306,   303,   308,   302,   304,     0,     0,     0,   282,     0,
       0,   305,     0,     0,   305,   268,     0,     0,     0,   121,
       0,     0,   138,   147,   148,     0,   152,   133,   132,     0,
     134,     0,   131,   135,   136,     0,   141,   139,   562,   563,
     419,   430,   432,   433,   434,   431,     0,   424,   428,     0,
       0,     0,     0,     0,     0,     0,   117,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    86,    96,
      95,   100,    99,    98,    97,    91,    90,    92,    93,    94,
       0,     0,   388,     0,     0,   515,   523,   508,   514,   520,
     521,   512,   522,   527,   513,   510,   526,   392,     0,    77,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   153,   162,   161,   166,   165,   164,
     163,   158,   157,   159,   160,   167,   437,   358,   356,     0,
     373,     0,     0,     0,     0,   210,   211,     0,     0,   227,
     226,   217,   216,   235,   232,     0,   572,   442,   439,     0,
      59,    65,    71,   581,   582,   583,   584,   585,   586,   587,
     588,   589,   590,   591,   592,   593,   594,   595,   596,   597,
     598,   599,   600,   601,   602,   603,   604,   605,   606,   607,
     476,   477,   270,   269,   609,   611,   613,   612,     0,   466,
     474,     0,     0,   403,   402,     0,   413,     0,   407,     0,
     125,     0,   377,     0,   314,   263,   278,   310,   309,   275,
     305,   305,     0,   305,     0,     0,   293,     0,   267,   266,
       0,     0,     0,     0,     0,     0,   143,     0,     0,     0,
       0,   418,   305,   429,     0,     0,     0,     0,     0,   113,
       0,     0,     0,   116,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,     0,     0,     0,   386,   394,   180,
       0,     0,   185,   170,   171,   172,   173,   174,   175,   176,
     177,   178,   179,   357,   374,   212,   204,   203,   207,   208,
       0,     0,   234,   441,     0,   608,   401,     0,   398,   415,
       0,   409,     0,     0,     0,   542,   279,     0,     0,     0,
     305,     0,   305,   305,   291,     0,   120,     0,   151,   524,
     130,     0,     0,     0,     0,   425,     0,     0,   192,     0,
     200,     0,     0,     0,     0,   118,   383,   389,     0,     0,
       0,     0,   209,   610,     0,     0,   416,   410,     0,   378,
     305,   305,   305,     0,   299,     0,     0,     0,   183,     0,
     146,   142,   140,   305,   426,     0,     0,     0,   195,     0,
       0,   191,   114,   115,   387,   181,   182,   206,     0,   399,
     305,   284,   280,   283,   305,   297,   292,   129,     0,     0,
     194,   193,   199,     0,   197,     0,     0,     0,   376,   305,
       0,     0,   144,   427,   196,     0,   202,     0,   400,     0,
     285,     0,   298,   198,   190,   201,   375,   305,   305,   286,
     281,   305,   287
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1112, -1112,  1351, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,  -322, -1112,
   -1112, -1112, -1112,  -113,  -225, -1112, -1112,  1081, -1112,   324,
   -1112, -1112, -1112, -1112, -1112, -1112,  -966,  -612,  -115,  -605,
   -1112, -1112, -1112,  1264,  -240, -1112, -1112, -1112, -1112,   420,
   -1112, -1112,   670, -1112, -1112,   832, -1112, -1112,   675, -1112,
   -1112,  -106,   -24,   709,   857, -1112, -1112,  1101, -1112, -1112,
   -1111, -1112, -1112,  1093, -1112, -1112,  1100, -1001,  -606, -1112,
   -1112,   809, -1112,  1280,   696, -1112,   274, -1112, -1112, -1112,
   -1112,  1054, -1112, -1112, -1112, -1112, -1112, -1112, -1112,  1211,
    -756, -1112, -1112, -1112, -1112, -1112,   786, -1112,   343,  -838,
   -1112, -1112, -1112, -1112, -1112,   685, -1112,   -35,   874, -1112,
   -1112,   868, -1112, -1112,   660, -1112,  -503, -1112,   -93, -1112,
    1314, -1112, -1112, -1112, -1112, -1112, -1112, -1112,  -101, -1112,
   -1112,  -114,  -582, -1112, -1112, -1112, -1112,  -100,   -79,   -77,
     -76,   -71, -1112, -1112,   -99,   -78, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112,   -65, -1112, -1112, -1112, -1112, -1112,
     -57,   -54,   -67,   -53,   -48,   -47, -1112, -1112, -1112, -1112,
     -98,   -96,   -89,   -87,   -42,   -41,   -40, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112, -1112,
   -1112, -1112, -1112, -1112, -1112,   664, -1112,  -532
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    44,    45,    46,    47,    48,    49,    50,    51,    52,
     152,   154,   156,   131,    53,    54,    55,    56,   363,   907,
      57,   324,    58,   228,   229,    59,   320,   321,   881,   882,
      60,   327,  1079,  1078,  1161,   885,   629,   630,   631,   632,
     438,    61,    62,   342,   343,  1171,  1247,    63,   732,   733,
      64,   462,   463,    65,   214,   215,    66,   458,   459,    67,
     465,   469,   110,   868,   784,    68,   306,   307,   308,   856,
    1146,    69,   317,   318,    70,   312,   313,   857,  1147,    71,
     259,   260,    72,   439,   440,    73,  1052,  1053,    74,    75,
     365,   366,    76,    77,   368,    78,    79,    80,   211,   212,
     568,    81,    82,    83,    84,   335,   336,   896,   897,   898,
      85,   134,   724,    86,   470,   471,   179,   180,   181,    87,
     203,   204,    88,    89,   515,   516,   780,   387,   388,   389,
     390,   391,   392,   393,   394,   395,   396,   397,   398,   399,
     400,   401,   402,   884,   403,   404,   405,   182,   183,   184,
     185,   186,   268,   269,   406,   443,   272,   273,   274,   275,
     276,   277,   278,   279,   444,   281,   282,   283,   284,   285,
     445,   446,   447,   448,   449,   450,   407,   292,   293,   337,
     408,   409,   187,   188,   453,   189,   190,   299,   472,   191,
     192,   193,   194,   195,   196,   197,   207,   517,   518,   519,
     520,   521,   522,   523,   524,   525,   526,   527,   528,   529,
     530,   531,   532,   533,   534,   535,   536,   537,   538,   539,
     540,   541,   542,   543,   799,  1035,   793,   794
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       128,   129,   322,   584,   858,   216,   860,   137,   262,   263,
     270,   863,   146,   149,   150,   205,   261,   873,   157,   437,
     294,   386,   295,   338,   874,   339,   206,   460,   648,   649,
     264,   271,   265,   266,   441,   441,   442,   442,   267,   466,
     660,   461,   288,   785,   280,   677,   451,   451,   452,   452,
     464,   883,   286,  1042,   825,   287,   289,   803,   202,  1083,
    1148,   290,   291,   900,  1034,   891,    90,   296,   297,   298,
     418,   102,   340,   872,   985,   791,   848,   303,   729,   865,
    1222,   847,    92,   986,   419,   850,   891,   730,   107,   831,
     832,   833,  1128,   420,   171,  1208,   418,   584,   840,  1200,
     219,  1129,   643,   583,  1080,  1087,   421,   121,   572,   983,
     419,  1162,  1163,  1164,   422,   984,   852,   566,   888,   420,
     849,   577,   101,   104,   423,  1088,   622,   578,   851,   107,
     892,   209,   421,   328,   107,   340,   107,   132,   304,   340,
     422,    96,   889,    99,   117,   642,   107,  1260,   106,   728,
     423,   892,   100,   118,   227,   133,   107,   107,   123,   107,
     107,   111,   300,   107,   107,   107,  1081,   901,   853,   987,
     379,   107,   573,   866,   867,   300,   305,   319,   936,   855,
     567,  1209,   341,  1082,   875,   943,   893,   424,   425,   103,
     894,   895,   329,   426,   427,   428,   429,   430,   431,   432,
     433,   434,   330,   854,  1203,   108,   109,   893,   435,   107,
     210,   894,   895,   424,   425,    91,   643,   301,  1036,   426,
     427,   428,   429,   430,   431,   432,   433,   434,   792,   902,
     301,    93,   988,  1228,   435,  1021,   731,   996,   891,  1251,
    1237,   105,   436,  1210,   628,   341,   126,   127,   220,   341,
    1130,   144,   145,   147,   148,  1064,   300,   413,  1067,   300,
     824,  1018,  1019,   817,   418,  1022,  1023,  1201,   436,  1026,
     628,  1028,  1029,   821,   842,   107,   947,   978,   419,   302,
     980,   994,   998,   300,    94,    95,   112,   420,  1139,   713,
     714,   715,   410,   716,   717,   718,   719,   720,   721,   722,
     421,   723,   113,   892,   114,   300,   300,   476,   422,   700,
     701,   301,   414,   602,   301,   480,   484,   300,   423,   847,
     300,   712,   903,   904,   905,   906,  1141,  1083,   908,   909,
     910,   911,   912,   913,   914,   915,   916,   917,   301,   919,
     920,   921,   922,   923,   924,   925,   926,   927,   928,   929,
     876,   930,   122,   300,   124,   633,   848,   934,   849,   893,
     301,   301,   477,   894,   895,   850,   851,   125,   231,   781,
     481,   485,   301,   411,   415,   301,   455,   603,   309,   638,
    1184,   424,   425,   200,   725,   725,   209,   426,   427,   428,
     429,   430,   431,   432,   433,   434,   852,   300,   734,   300,
     467,   300,   435,   232,   233,   300,   853,   201,   301,   982,
     234,   300,   300,   300,   569,   130,   300,   235,   848,   634,
     580,   303,   473,   474,   478,   736,   581,   850,   738,   624,
     741,   608,   482,   486,   487,   314,   436,  1046,   628,   310,
    1075,   854,  1092,   639,  1144,   252,   674,   678,   726,   727,
      97,    98,   301,   254,   301,  1149,   301,  1151,   852,   855,
     301,  1156,   735,   115,   116,   210,   301,   301,   301,   256,
    1048,   301,   782,   783,   778,   309,  1166,   311,   314,   216,
     227,   257,   304,   779,   205,   614,   135,   258,   619,   737,
    1045,   138,   739,   883,   742,   206,   315,   364,   679,   177,
     178,   262,   263,   270,  1076,   139,  1093,   332,  1145,   261,
     119,   120,   227,   294,  1054,   295,  1108,   151,  1123,   333,
     305,   855,  1124,   264,   271,   265,   266,   202,  1132,  1133,
    1186,   267,   334,  1187,   316,   288,   310,   280,   136,   315,
     338,   153,   339,   155,  1193,   286,  1195,  1196,   287,   289,
     873,   873,   873,   818,   290,   291,   822,   874,   874,   874,
     296,   297,   298,   810,  1057,   140,   141,   221,   224,   142,
     143,  1159,   811,  1058,   311,   222,   225,   316,  1085,   843,
     162,   950,   951,   952,  1221,   198,  1223,   954,   955,   956,
     957,   958,   959,   960,   961,   962,   963,  1229,   965,   966,
     967,   968,   969,   970,   971,   972,   973,   974,   975,  1105,
     158,   159,   217,   370,  1238,   208,   873,   213,  1241,   218,
     460,   223,   226,   874,   441,   418,   442,   227,  1202,   230,
     319,   325,   993,  1250,   461,   877,   451,   323,   452,   419,
     326,   669,   670,   464,   671,   364,   367,   878,   420,   721,
     722,  1259,   723,   879,   412,  1262,   416,   417,   475,   457,
     479,   421,   418,   483,   544,  1242,  1131,   557,   948,   422,
     545,   546,   547,   880,   548,   549,   419,   550,   551,   423,
     552,   553,   661,   662,   663,   420,   664,   665,   666,   667,
     668,   669,   670,   554,   671,   555,   556,   558,   421,   559,
     560,   561,   979,   981,   562,   563,   422,   667,   668,   669,
     670,   564,   671,   571,   565,   995,   423,   574,   999,   719,
     720,   721,   722,   575,   723,   163,   164,   165,   166,   167,
     168,   169,   576,   579,   582,   585,   586,   170,   587,   588,
     589,   171,   424,   425,   590,   918,   591,   592,   426,   427,
     428,   429,   430,   431,   432,   433,   434,  1071,   172,   418,
    1073,   593,   594,   435,   605,   595,  1172,   596,  1173,  1174,
    1175,   597,   598,   419,   599,   600,   601,   604,   606,   424,
     425,   607,   420,   610,  1178,   426,   427,   428,   429,   430,
     431,   432,   433,   434,   611,   421,   612,   436,   613,   628,
     435,   173,   174,   422,   616,   618,   617,   621,  1181,   625,
     626,   627,   344,   423,   635,  1185,   636,   640,   637,   175,
     176,   644,   641,   645,   646,   647,   345,   650,   651,   652,
     653,   680,   681,  1197,   436,   346,   628,   344,   682,   683,
     654,   684,   694,   655,   685,  1205,   584,  1206,   347,   656,
     686,   345,   177,   178,   657,   658,   348,   659,   673,   687,
     346,   688,   689,   690,   691,   692,   349,   740,   744,   695,
     746,   693,   696,   347,   697,   698,   424,   425,   699,   702,
     743,   348,   426,   427,   428,   429,   430,   431,   432,   433,
     434,   349,   703,   704,   705,   706,  1236,   435,   663,   745,
     664,   665,   666,   667,   668,   669,   670,   747,   671,   707,
     713,   714,   715,  1246,   716,   717,   718,   719,   720,   721,
     722,   676,   723,   748,   708,  1255,   709,   710,   749,   350,
     351,   436,   750,   628,   751,   352,   353,   354,   355,   356,
     357,   358,   359,   360,   752,   753,   933,   711,   754,   755,
     361,   756,   776,   757,   350,   351,   344,  1157,   758,   759,
     352,   353,   354,   355,   356,   357,   358,   359,   360,   760,
     345,   761,   762,   964,   763,   361,   764,   765,   766,   346,
     344,   661,   662,   663,   362,   664,   665,   666,   667,   668,
     669,   670,   347,   671,   345,  1179,  1180,   767,   768,   769,
     348,   770,   771,   346,   772,   773,   418,   774,   775,   362,
     349,   777,   786,   788,   789,   797,   347,   798,   790,   807,
     419,   795,   808,  1140,   348,  1142,   812,   796,   801,   420,
     802,   661,   662,   663,   349,   664,   665,   666,   667,   668,
     669,   670,   421,   671,  1089,   804,   806,   813,   815,   816,
     422,   820,   792,   826,   827,   809,   828,   661,   662,   663,
     423,   664,   665,   666,   667,   668,   669,   670,   819,   671,
     823,   829,   830,   350,   351,   834,   835,   836,   837,   352,
     353,   354,   355,   356,   357,   358,   359,   360,   838,   839,
     844,   845,   672,   846,   361,   859,   861,   350,   351,   869,
     862,   870,   886,   352,   353,   354,   355,   356,   357,   358,
     359,   360,   871,   887,   931,   935,   937,   864,   361,   899,
    1094,   938,   939,   424,   425,   940,   941,  1254,   362,   426,
     427,   428,   429,   430,   431,   432,   433,   434,   160,   671,
     942,   944,   945,   946,   435,     1,     2,   723,   949,   976,
     990,  1003,   362,   992,  1000,     3,     4,     5,  1004,  1001,
    1002,  1005,     6,  1006,  1007,  1008,     7,  1009,     8,     9,
    1010,    10,  1011,    11,    12,    13,    14,   715,   436,   716,
     717,   718,   719,   720,   721,   722,    15,   723,  1012,    16,
    1030,   661,   662,   663,  1013,   664,   665,   666,   667,   668,
     669,   670,    17,   671,  1014,  1015,  1016,  1017,  1020,  1024,
    1025,  1027,  1032,  1033,  1034,    18,    19,    20,   661,   662,
     663,    21,   664,   665,   666,   667,   668,   669,   670,  1041,
     671,   567,    22,  1055,    23,  1039,    24,    25,    26,    27,
      28,  1038,  1040,  1043,  1056,    29,    30,  1044,  1047,  1049,
      31,    32,    33,    34,  1095,  1051,  1059,  1060,  1062,  1061,
      35,    36,  1063,    37,     1,     2,  1065,    38,  1066,  1068,
      39,    40,    41,    42,     3,     4,     5,  1069,  1070,   932,
    1072,     6,  1074,  1077,  1084,     7,  1086,     8,     9,  1106,
      10,    -1,    11,    12,    13,    14,  1112,  1127,  1135,  1134,
    1136,  1137,  1143,  1150,  1152,    15,    43,  1154,    16,  1155,
     661,   662,   663,  1167,   664,   665,   666,   667,   668,   669,
     670,    17,   671,  1153,  1169,  1170,  1176,  1190,  1183,  1188,
    1191,  1192,  1194,  1199,    18,    19,    20,   661,   662,   663,
      21,   664,   665,   666,   667,   668,   669,   670,  1204,   671,
    1218,    22,  1207,    23,  1220,    24,    25,    26,    27,    28,
    1224,  1232,  1234,  1225,    29,    30,  1235,  1226,  1233,    31,
      32,    33,    34,  1096,  1239,  1240,  1243,  1244,   231,    35,
      36,  1249,    37,  1245,  1252,  1253,    38,  1256,  1257,    39,
      40,    41,    42,   200,   170,   161,  1258,  1261,   171,  1160,
    1097,   623,   456,  1126,   814,   787,   991,   609,   953,   989,
     620,   841,   615,   232,   233,   172,   454,   201,  1189,   675,
     234,   977,   570,  1165,   890,    43,   997,   235,   236,   237,
     805,   800,   238,   239,   331,   240,   241,  1031,     0,   242,
     243,   244,   245,   246,   247,   248,     0,   249,   250,   251,
       0,     0,     0,     0,     0,   252,     0,  1037,   173,   174,
       0,   253,     0,   254,     0,     0,     0,     0,   255,     0,
       0,     0,     0,     0,     0,     0,   175,   176,     0,   256,
     369,     0,   163,   164,   165,   166,   167,   168,   169,   199,
       0,   257,   213,   200,   170,     0,     0,   258,   171,     0,
       0,   370,     0,   371,   372,     0,     0,     0,     0,   177,
     178,     0,     0,     0,     0,   172,     0,   201,     0,     0,
       0,     0,     0,     0,   234,     0,   373,   374,     0,     0,
       0,   235,     0,     0,     0,   661,   662,   663,   328,   664,
     665,   666,   667,   668,   669,   670,     0,   671,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   173,   174,
     369,     0,     0,     0,   375,     0,   376,   254,   377,   333,
       0,     0,     0,     0,   378,     0,   175,   176,   379,     0,
       0,   370,   334,   371,   372,     0,   380,   381,   382,     0,
       0,     0,   383,   384,   385,     0,   213,     0,  1098,     0,
       0,     0,     0,   468,   234,     0,   373,   374,     0,   177,
     178,   235,     0,     0,     0,   661,   662,   663,   328,   664,
     665,   666,   667,   668,   669,   670,     0,   671,   661,   662,
     663,     0,   664,   665,   666,   667,   668,   669,   670,     0,
     671,     0,     0,     0,   375,     0,   376,   254,   377,   333,
       0,     0,     0,     0,   378,     0,     0,     0,   379,     0,
       0,     0,   334,     0,     0,     0,   380,   381,   382,     0,
       0,     0,   383,   384,   385,     0,   213,     0,  1099,     0,
     661,   662,   663,     0,   664,   665,   666,   667,   668,   669,
     670,  1100,   671,   661,   662,   663,     0,   664,   665,   666,
     667,   668,   669,   670,     0,   671,   661,   662,   663,     0,
     664,   665,   666,   667,   668,   669,   670,     0,   671,   661,
     662,   663,     0,   664,   665,   666,   667,   668,   669,   670,
       0,   671,   713,   714,   715,     0,   716,   717,   718,   719,
     720,   721,   722,  1101,   723,   713,   714,   715,     0,   716,
     717,   718,   719,   720,   721,   722,  1102,   723,   713,   714,
     715,     0,   716,   717,   718,   719,   720,   721,   722,  1103,
     723,   713,   714,   715,     0,   716,   717,   718,   719,   720,
     721,   722,  1104,   723,   713,   714,   715,     0,   716,   717,
     718,   719,   720,   721,   722,  1109,   723,   713,   714,   715,
       0,   716,   717,   718,   719,   720,   721,   722,  1113,   723,
     713,   714,   715,     0,   716,   717,   718,   719,   720,   721,
     722,  1114,   723,   713,   714,   715,     0,   716,   717,   718,
     719,   720,   721,   722,  1115,   723,   713,   714,   715,     0,
     716,   717,   718,   719,   720,   721,   722,  1116,   723,   713,
     714,   715,     0,   716,   717,   718,   719,   720,   721,   722,
    1117,   723,   713,   714,   715,     0,   716,   717,   718,   719,
     720,   721,   722,  1118,   723,   661,   662,   663,     0,   664,
     665,   666,   667,   668,   669,   670,  1119,   671,   661,   662,
     663,     0,   664,   665,   666,   667,   668,   669,   670,  1120,
     671,   661,   662,   663,     0,   664,   665,   666,   667,   668,
     669,   670,  1121,   671,   661,   662,   663,     0,   664,   665,
     666,   667,   668,   669,   670,  1122,   671,   713,   714,   715,
       0,   716,   717,   718,   719,   720,   721,   722,  1125,   723,
     713,   714,   715,     0,   716,   717,   718,   719,   720,   721,
     722,  1182,   723,   661,   662,   663,     0,   664,   665,   666,
     667,   668,   669,   670,  1212,   671,   661,   662,   663,     0,
     664,   665,   666,   667,   668,   669,   670,  1213,   671,   661,
     662,   663,     0,   664,   665,   666,   667,   668,   669,   670,
    1215,   671,   661,   662,   663,     0,   664,   665,   666,   667,
     668,   669,   670,  1216,   671,   713,   714,   715,     0,   716,
     717,   718,   719,   720,   721,   722,  1217,   723,     0,     0,
       0,     0,     0,     0,   661,   662,   663,  1050,   664,   665,
     666,   667,   668,   669,   670,     0,   671,   661,   662,   663,
    1107,   664,   665,   666,   667,   668,   669,   670,     0,   671,
     713,   714,   715,  1138,   716,   717,   718,   719,   720,   721,
     722,     0,   723,   661,   662,   663,  1158,   664,   665,   666,
     667,   668,   669,   670,     0,   671,     0,     0,     0,     0,
       0,     0,   661,   662,   663,  1168,   664,   665,   666,   667,
     668,   669,   670,     0,   671,   661,   662,   663,  1177,   664,
     665,   666,   667,   668,   669,   670,     0,   671,   661,   662,
     663,  1198,   664,   665,   666,   667,   668,   669,   670,     0,
     671,   661,   662,   663,  1211,   664,   665,   666,   667,   668,
     669,   670,     0,   671,     0,     0,     0,     0,     0,     0,
     661,   662,   663,  1214,   664,   665,   666,   667,   668,   669,
     670,     0,   671,   661,   662,   663,  1219,   664,   665,   666,
     667,   668,   669,   670,     0,   671,     0,     0,     0,  1227,
    1090,   661,   662,   663,     0,   664,   665,   666,   667,   668,
     669,   670,  1230,   671,     0,  1091,   661,   662,   663,     0,
     664,   665,   666,   667,   668,   669,   670,     0,   671,     0,
       0,  1231,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1248,   488,   489,   490,   491,   492,
     493,   494,   495,   496,   497,   498,   499,   500,   501,   502,
     503,   504,   505,   506,   507,   508,     0,     0,     0,   509,
     510,     0,   511,   512,   513,   514,  1110,   713,   714,   715,
       0,   716,   717,   718,   719,   720,   721,   722,     0,   723,
    1111,   713,   714,   715,     0,   716,   717,   718,   719,   720,
     721,   722,     0,   723,   713,   714,   715,     0,   716,   717,
     718,   719,   720,   721,   722,     0,   723,   661,   662,   663,
       0,   664,   665,   666,   667,   668,   669,   670,     0,   671
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   115,   228,   610,    98,   612,    31,   109,   109,
     109,   617,    36,    37,    38,    93,   109,   629,    42,   134,
     109,   127,   109,   121,   629,   121,    93,   141,   350,   351,
     109,   109,   109,   109,   135,   136,   135,   136,   109,   145,
     362,   141,   109,   546,   109,   367,   135,   136,   135,   136,
     143,   633,   109,   809,   586,   109,   109,   560,    93,   897,
    1061,   109,   109,    34,    52,    83,    52,   109,   109,   109,
      29,    83,    22,    32,    43,    52,    43,    22,    43,    52,
    1191,     6,    52,    52,    43,    52,    83,    52,    83,   592,
     593,   594,    43,    52,    25,    52,    29,   322,   601,    32,
      52,    52,   342,    32,   138,   101,    65,   201,   138,   200,
      43,  1077,  1078,  1079,    73,   206,    83,   138,    57,    52,
      45,   202,    52,    83,    83,   121,    32,   208,    53,    83,
     148,     4,    65,    61,    83,    22,    83,    34,    83,    22,
      73,   200,    81,    34,    34,    32,    83,  1258,   200,    32,
      83,   148,    43,    43,    83,    52,    83,    83,    34,    83,
      83,   200,    83,    83,    83,    83,   200,   138,    93,   138,
     101,    83,   202,   146,   147,    83,   121,    83,   681,   146,
     201,   138,   132,   201,    34,   688,   204,   146,   147,   201,
     208,   209,   120,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   130,   128,   201,   200,   201,   204,   167,    83,
      83,   208,   209,   146,   147,   201,   456,   138,   206,   152,
     153,   154,   155,   156,   157,   158,   159,   160,   205,   200,
     138,   201,   201,  1199,   167,   767,   201,   740,    83,  1240,
     207,   201,   201,   200,   203,   132,   200,   201,   200,   132,
     201,   200,   201,   200,   201,   861,    83,    83,   864,    83,
     582,   764,   765,   200,    29,   768,   769,    32,   201,   772,
     203,   774,   775,   200,   200,    83,   200,   200,    43,   200,
     200,   200,   200,    83,   200,   201,   200,    52,   200,   139,
     140,   141,   200,   143,   144,   145,   146,   147,   148,   149,
      65,   151,   200,   148,   200,    83,    83,    83,    73,   424,
     425,   138,   138,   138,   138,    83,    83,    83,    83,     6,
      83,   436,   644,   645,   646,   647,   200,  1165,   650,   651,
     652,   653,   654,   655,   656,   657,   658,   659,   138,   661,
     662,   663,   664,   665,   666,   667,   668,   669,   670,   671,
     200,   673,   200,    83,   200,   138,    43,   679,    45,   204,
     138,   138,   138,   208,   209,    52,    53,   200,     5,    43,
     138,   138,   138,   200,   200,   138,   200,   202,    22,   138,
    1136,   146,   147,    20,   138,   138,     4,   152,   153,   154,
     155,   156,   157,   158,   159,   160,    83,    83,   138,    83,
     200,    83,   167,    40,    41,    83,    93,    44,   138,   731,
      47,    83,    83,    83,    32,    83,    83,    54,    43,   202,
     202,    22,   200,   200,   200,   138,   208,    52,   138,    32,
     138,    32,   200,   200,   200,    22,   201,   200,   203,    83,
     138,   128,   138,   202,   138,    82,    32,    32,   202,   202,
     200,   201,   138,    90,   138,  1061,   138,  1063,    83,   146,
     138,    32,   202,   200,   201,    83,   138,   138,   138,   106,
     200,   138,   146,   147,    43,    22,  1082,   121,    22,   572,
      83,   118,    83,    52,   562,    32,   201,   124,    32,   202,
     812,   200,   202,  1075,   202,   562,    83,    83,    83,   136,
     137,   602,   602,   602,   202,    34,   202,    80,   202,   602,
     200,   201,    83,   602,   200,   602,   200,    83,   200,    92,
     121,   146,   200,   602,   602,   602,   602,   562,   200,   200,
     200,   602,   105,   200,   121,   602,    83,   602,   201,    83,
     638,    83,   638,    83,  1150,   602,  1152,  1153,   602,   602,
    1162,  1163,  1164,   577,   602,   602,   580,  1162,  1163,  1164,
     602,   602,   602,    43,    43,   200,   201,   200,   200,   200,
     201,  1074,    52,    52,   121,   208,   208,   121,   900,   603,
     200,   696,   697,   698,  1190,   200,  1192,   702,   703,   704,
     705,   706,   707,   708,   709,   710,   711,  1203,   713,   714,
     715,   716,   717,   718,   719,   720,   721,   722,   723,   931,
     200,   201,    43,    24,  1220,    23,  1228,   119,  1224,   200,
     734,    83,    83,  1228,   725,    29,   725,    83,    32,    34,
      83,    43,   738,  1239,   734,    46,   725,    39,   725,    43,
     200,   148,   149,   736,   151,    83,    83,    58,    52,   148,
     149,  1257,   151,    64,   127,  1261,    52,   200,   127,   205,
     127,    65,    29,   127,    34,    32,   988,   138,   692,    73,
      34,    34,    34,    84,    34,    34,    43,    34,    34,    83,
      34,    34,   139,   140,   141,    52,   143,   144,   145,   146,
     147,   148,   149,    34,   151,    34,    34,   202,    65,    34,
      34,    34,   726,   727,   138,   202,    73,   146,   147,   148,
     149,   202,   151,    34,    83,   739,    83,   200,   742,   146,
     147,   148,   149,   200,   151,     9,    10,    11,    12,    13,
      14,    15,    83,    83,    34,    83,    34,    21,    34,    34,
      34,    25,   146,   147,    34,   202,    34,    34,   152,   153,
     154,   155,   156,   157,   158,   159,   160,   870,    42,    29,
     875,    34,    34,   167,    83,    34,  1088,    34,  1090,  1091,
    1092,    34,    34,    43,    34,    34,    34,    34,    83,   146,
     147,    83,    52,   138,  1106,   152,   153,   154,   155,   156,
     157,   158,   159,   160,    83,    65,   138,   201,    83,   203,
     167,    85,    86,    73,    83,    83,   138,   201,  1130,    34,
     202,   200,    29,    83,   202,  1137,   202,    83,    34,   103,
     104,   201,    83,   201,   201,   201,    43,   201,   201,   201,
     201,    34,    34,  1155,   201,    52,   203,    29,    34,    34,
     201,    34,    83,   201,    34,  1167,  1071,  1169,    65,   201,
      34,    43,   136,   137,   201,   201,    73,   201,   201,    34,
      52,    34,    34,    34,    34,   202,    83,    34,    83,   200,
      83,   127,   201,    65,   201,   201,   146,   147,   201,   201,
     127,    73,   152,   153,   154,   155,   156,   157,   158,   159,
     160,    83,   201,   201,   201,   201,  1218,   167,   141,   127,
     143,   144,   145,   146,   147,   148,   149,   127,   151,   201,
     139,   140,   141,  1235,   143,   144,   145,   146,   147,   148,
     149,   138,   151,    83,   201,  1247,   201,   201,    34,   146,
     147,   201,    34,   203,    34,   152,   153,   154,   155,   156,
     157,   158,   159,   160,    34,    34,   138,   201,    34,    34,
     167,    34,   202,    34,   146,   147,    29,  1072,    34,    34,
     152,   153,   154,   155,   156,   157,   158,   159,   160,    34,
      43,    34,    34,   202,    34,   167,    34,    34,    34,    52,
      29,   139,   140,   141,   201,   143,   144,   145,   146,   147,
     148,   149,    65,   151,    43,  1110,  1111,    34,    34,    34,
      73,    34,    34,    52,    34,    34,    29,    34,    34,   201,
      83,   138,    52,    52,    83,    83,    65,    83,    52,   200,
      43,    52,   201,  1047,    73,  1049,    34,    52,    52,    52,
      52,   139,   140,   141,    83,   143,   144,   145,   146,   147,
     148,   149,    65,   151,   202,    52,    52,    52,   200,   200,
      73,   200,   205,    52,    52,    83,    52,   139,   140,   141,
      83,   143,   144,   145,   146,   147,   148,   149,    83,   151,
      83,    52,    52,   146,   147,    52,    52,    52,    52,   152,
     153,   154,   155,   156,   157,   158,   159,   160,    83,   201,
      83,    34,   200,   138,   167,   138,   138,   146,   147,    83,
     138,   200,   200,   152,   153,   154,   155,   156,   157,   158,
     159,   160,    83,   200,    83,    52,    52,   138,   167,   138,
     202,    52,    52,   146,   147,    52,    52,   200,   201,   152,
     153,   154,   155,   156,   157,   158,   159,   160,     0,   151,
      52,    52,    52,    52,   167,     7,     8,   151,   127,   200,
     200,    52,   201,   200,   127,    17,    18,    19,    52,   127,
     127,    52,    24,    52,    52,    52,    28,    52,    30,    31,
      52,    33,    52,    35,    36,    37,    38,   141,   201,   143,
     144,   145,   146,   147,   148,   149,    48,   151,    52,    51,
     200,   139,   140,   141,    52,   143,   144,   145,   146,   147,
     148,   149,    64,   151,    52,    52,    52,    52,    52,    52,
      52,    52,    43,    43,    52,    77,    78,    79,   139,   140,
     141,    83,   143,   144,   145,   146,   147,   148,   149,    52,
     151,   201,    94,    83,    96,   200,    98,    99,   100,   101,
     102,   204,   200,   202,    83,   107,   108,   202,   202,   202,
     112,   113,   114,   115,   202,   207,   200,   138,    83,   138,
     122,   123,   138,   125,     7,     8,    83,   129,   200,    52,
     132,   133,   134,   135,    17,    18,    19,    52,   202,   200,
      34,    24,    34,   200,    83,    28,    83,    30,    31,    83,
      33,   151,    35,    36,    37,    38,   202,   200,    83,   204,
     202,    34,    83,   138,   138,    48,   168,   200,    51,    34,
     139,   140,   141,    34,   143,   144,   145,   146,   147,   148,
     149,    64,   151,   138,    34,    52,   200,   138,    52,   207,
     138,   138,   200,   200,    77,    78,    79,   139,   140,   141,
      83,   143,   144,   145,   146,   147,   148,   149,   202,   151,
      34,    94,   204,    96,   138,    98,    99,   100,   101,   102,
     138,    52,    52,   200,   107,   108,   131,   200,   204,   112,
     113,   114,   115,   202,   138,   138,   202,    52,     5,   122,
     123,    83,   125,   204,   200,    52,   129,   207,   138,   132,
     133,   134,   135,    20,    21,    44,   138,   138,    25,  1075,
     202,   320,   138,   983,   572,   548,   736,   306,   699,   734,
     317,   602,   312,    40,    41,    42,   136,    44,  1144,   365,
      47,   725,   211,  1080,   638,   168,   741,    54,    55,    56,
     562,   557,    59,    60,   120,    62,    63,   777,    -1,    66,
      67,    68,    69,    70,    71,    72,    -1,    74,    75,    76,
      -1,    -1,    -1,    -1,    -1,    82,    -1,   793,    85,    86,
      -1,    88,    -1,    90,    -1,    -1,    -1,    -1,    95,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   103,   104,    -1,   106,
       3,    -1,     9,    10,    11,    12,    13,    14,    15,    16,
      -1,   118,   119,    20,    21,    -1,    -1,   124,    25,    -1,
      -1,    24,    -1,    26,    27,    -1,    -1,    -1,    -1,   136,
     137,    -1,    -1,    -1,    -1,    42,    -1,    44,    -1,    -1,
      -1,    -1,    -1,    -1,    47,    -1,    49,    50,    -1,    -1,
      -1,    54,    -1,    -1,    -1,   139,   140,   141,    61,   143,
     144,   145,   146,   147,   148,   149,    -1,   151,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    85,    86,
       3,    -1,    -1,    -1,    87,    -1,    89,    90,    91,    92,
      -1,    -1,    -1,    -1,    97,    -1,   103,   104,   101,    -1,
      -1,    24,   105,    26,    27,    -1,   109,   110,   111,    -1,
      -1,    -1,   115,   116,   117,    -1,   119,    -1,   202,    -1,
      -1,    -1,    -1,   126,    47,    -1,    49,    50,    -1,   136,
     137,    54,    -1,    -1,    -1,   139,   140,   141,    61,   143,
     144,   145,   146,   147,   148,   149,    -1,   151,   139,   140,
     141,    -1,   143,   144,   145,   146,   147,   148,   149,    -1,
     151,    -1,    -1,    -1,    87,    -1,    89,    90,    91,    92,
      -1,    -1,    -1,    -1,    97,    -1,    -1,    -1,   101,    -1,
      -1,    -1,   105,    -1,    -1,    -1,   109,   110,   111,    -1,
      -1,    -1,   115,   116,   117,    -1,   119,    -1,   202,    -1,
     139,   140,   141,    -1,   143,   144,   145,   146,   147,   148,
     149,   202,   151,   139,   140,   141,    -1,   143,   144,   145,
     146,   147,   148,   149,    -1,   151,   139,   140,   141,    -1,
     143,   144,   145,   146,   147,   148,   149,    -1,   151,   139,
     140,   141,    -1,   143,   144,   145,   146,   147,   148,   149,
      -1,   151,   139,   140,   141,    -1,   143,   144,   145,   146,
     147,   148,   149,   202,   151,   139,   140,   141,    -1,   143,
     144,   145,   146,   147,   148,   149,   202,   151,   139,   140,
     141,    -1,   143,   144,   145,   146,   147,   148,   149,   202,
     151,   139,   140,   141,    -1,   143,   144,   145,   146,   147,
     148,   149,   202,   151,   139,   140,   141,    -1,   143,   144,
     145,   146,   147,   148,   149,   202,   151,   139,   140,   141,
      -1,   143,   144,   145,   146,   147,   148,   149,   202,   151,
     139,   140,   141,    -1,   143,   144,   145,   146,   147,   148,
     149,   202,   151,   139,   140,   141,    -1,   143,   144,   145,
     146,   147,   148,   149,   202,   151,   139,   140,   141,    -1,
     143,   144,   145,   146,   147,   148,   149,   202,   151,   139,
     140,   141,    -1,   143,   144,   145,   146,   147,   148,   149,
     202,   151,   139,   140,   141,    -1,   143,   144,   145,   146,
     147,   148,   149,   202,   151,   139,   140,   141,    -1,   143,
     144,   145,   146,   147,   148,   149,   202,   151,   139,   140,
     141,    -1,   143,   144,   145,   146,   147,   148,   149,   202,
     151,   139,   140,   141,    -1,   143,   144,   145,   146,   147,
     148,   149,   202,   151,   139,   140,   141,    -1,   143,   144,
     145,   146,   147,   148,   149,   202,   151,   139,   140,   141,
      -1,   143,   144,   145,   146,   147,   148,   149,   202,   151,
     139,   140,   141,    -1,   143,   144,   145,   146,   147,   148,
     149,   202,   151,   139,   140,   141,    -1,   143,   144,   145,
     146,   147,   148,   149,   202,   151,   139,   140,   141,    -1,
     143,   144,   145,   146,   147,   148,   149,   202,   151,   139,
     140,   141,    -1,   143,   144,   145,   146,   147,   148,   149,
     202,   151,   139,   140,   141,    -1,   143,   144,   145,   146,
     147,   148,   149,   202,   151,   139,   140,   141,    -1,   143,
     144,   145,   146,   147,   148,   149,   202,   151,    -1,    -1,
      -1,    -1,    -1,    -1,   139,   140,   141,   200,   143,   144,
     145,   146,   147,   148,   149,    -1,   151,   139,   140,   141,
     200,   143,   144,   145,   146,   147,   148,   149,    -1,   151,
     139,   140,   141,   200,   143,   144,   145,   146,   147,   148,
     149,    -1,   151,   139,   140,   141,   200,   143,   144,   145,
     146,   147,   148,   149,    -1,   151,    -1,    -1,    -1,    -1,
      -1,    -1,   139,   140,   141,   200,   143,   144,   145,   146,
     147,   148,   149,    -1,   151,   139,   140,   141,   200,   143,
     144,   145,   146,   147,   148,   149,    -1,   151,   139,   140,
     141,   200,   143,   144,   145,   146,   147,   148,   149,    -1,
     151,   139,   140,   141,   200,   143,   144,   145,   146,   147,
     148,   149,    -1,   151,    -1,    -1,    -1,    -1,    -1,    -1,
     139,   140,   141,   200,   143,   144,   145,   146,   147,   148,
     149,    -1,   151,   139,   140,   141,   200,   143,   144,   145,
     146,   147,   148,   149,    -1,   151,    -1,    -1,    -1,   200,
     138,   139,   140,   141,    -1,   143,   144,   145,   146,   147,
     148,   149,   200,   151,    -1,   138,   139,   140,   141,    -1,
     143,   144,   145,   146,   147,   148,   149,    -1,   151,    -1,
      -1,   200,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   200,   169,   170,   171,   172,   173,
     174,   175,   176,   177,   178,   179,   180,   181,   182,   183,
     184,   185,   186,   187,   188,   189,    -1,    -1,    -1,   193,
     194,    -1,   196,   197,   198,   199,   138,   139,   140,   141,
      -1,   143,   144,   145,   146,   147,   148,   149,    -1,   151,
     138,   139,   140,   141,    -1,   143,   144,   145,   146,   147,
     148,   149,    -1,   151,   139,   140,   141,    -1,   143,   144,
     145,   146,   147,   148,   149,    -1,   151,   139,   140,   141,
      -1,   143,   144,   145,   146,   147,   148,   149,    -1,   151
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,    17,    18,    19,    24,    28,    30,    31,
      33,    35,    36,    37,    38,    48,    51,    64,    77,    78,
      79,    83,    94,    96,    98,    99,   100,   101,   102,   107,
     108,   112,   113,   114,   115,   122,   123,   125,   129,   132,
     133,   134,   135,   168,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   224,   225,   226,   227,   230,   232,   235,
     240,   251,   252,   257,   260,   263,   266,   269,   275,   281,
     284,   289,   292,   295,   298,   299,   302,   303,   305,   306,
     307,   311,   312,   313,   314,   320,   323,   329,   332,   333,
      52,   201,    52,   201,   200,   201,   200,   200,   201,    34,
      43,    52,    83,   201,    83,   201,   200,    83,   200,   201,
     272,   200,   200,   200,   200,   200,   201,    34,    43,   200,
     201,   201,   200,    34,   200,   200,   200,   201,   272,   272,
      83,   223,    34,    52,   321,   201,   201,   272,   200,    34,
     200,   201,   200,   201,   200,   201,   272,   200,   201,   272,
     272,    83,   220,    83,   221,    83,   222,   272,   200,   201,
       0,   212,   200,     9,    10,    11,    12,    13,    14,    15,
      21,    25,    42,    85,    86,   103,   104,   136,   137,   326,
     327,   328,   357,   358,   359,   360,   361,   392,   393,   395,
     396,   399,   400,   401,   402,   403,   404,   405,   200,    16,
      20,    44,   327,   330,   331,   365,   382,   406,    23,     4,
      83,   308,   309,   119,   264,   265,   338,    43,   200,    52,
     200,   200,   208,    83,   200,   208,    83,    83,   233,   234,
      34,     5,    40,    41,    47,    54,    55,    56,    59,    60,
      62,    63,    66,    67,    68,    69,    70,    71,    72,    74,
      75,    76,    82,    88,    90,    95,   106,   118,   124,   290,
     291,   338,   348,   357,   358,   359,   360,   361,   362,   363,
     364,   365,   366,   367,   368,   369,   370,   371,   372,   373,
     374,   375,   376,   377,   378,   379,   380,   381,   382,   383,
     384,   385,   387,   388,   392,   393,   394,   395,   396,   397,
      83,   138,   200,    22,    83,   121,   276,   277,   278,    22,
      83,   121,   285,   286,    22,    83,   121,   282,   283,    83,
     236,   237,   233,    39,   231,    43,   200,   241,    61,   120,
     130,   340,    80,    92,   105,   315,   316,   389,   390,   391,
      22,   132,   253,   254,    29,    43,    52,    65,    73,    83,
     146,   147,   152,   153,   154,   155,   156,   157,   158,   159,
     160,   167,   201,   228,    83,   300,   301,    83,   304,     3,
      24,    26,    27,    49,    50,    87,    89,    91,    97,   101,
     109,   110,   111,   115,   116,   117,   271,   337,   338,   339,
     340,   341,   342,   343,   344,   345,   346,   347,   348,   349,
     350,   351,   352,   354,   355,   356,   364,   386,   390,   391,
     200,   200,   127,    83,   138,   200,    52,   200,    29,    43,
      52,    65,    73,    83,   146,   147,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   167,   201,   248,   250,   293,
     294,   348,   364,   365,   374,   380,   381,   382,   383,   384,
     385,   392,   393,   394,   293,   200,   253,   205,   267,   268,
     351,   357,   261,   262,   338,   270,   271,   200,   126,   271,
     324,   325,   398,   200,   200,   127,    83,   138,   200,   127,
      83,   138,   200,   127,    83,   138,   200,   200,   169,   170,
     171,   172,   173,   174,   175,   176,   177,   178,   179,   180,
     181,   182,   183,   184,   185,   186,   187,   188,   189,   193,
     194,   196,   197,   198,   199,   334,   335,   407,   408,   409,
     410,   411,   412,   413,   414,   415,   416,   417,   418,   419,
     420,   421,   422,   423,   424,   425,   426,   427,   428,   429,
     430,   431,   432,   433,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,   138,   202,    34,
      34,    34,   138,   202,   202,    83,   138,   201,   310,    32,
     309,    34,   138,   202,   200,   200,    83,   202,   208,    83,
     202,   208,    34,    32,   234,    83,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,   138,   202,    34,    83,    83,    83,    32,   277,
     138,    83,   138,    83,    32,   286,    83,   138,    83,    32,
     283,   201,    32,   237,    32,    34,   202,   200,   203,   246,
     247,   248,   249,   138,   202,   202,   202,    34,   138,   202,
      83,    83,    32,   254,   201,   201,   201,   201,   228,   228,
     201,   201,   201,   201,   201,   201,   201,   201,   201,   201,
     228,   139,   140,   141,   143,   144,   145,   146,   147,   148,
     149,   151,   200,   201,    32,   301,   138,   228,    32,    83,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,   202,   127,    83,   200,   201,   201,   201,   201,
     248,   248,   201,   201,   201,   201,   201,   201,   201,   201,
     201,   201,   248,   139,   140,   141,   143,   144,   145,   146,
     147,   148,   149,   151,   322,   138,   202,   202,    32,    43,
      52,   201,   258,   259,   138,   202,   138,   202,   138,   202,
      34,   138,   202,   127,    83,   127,    83,   127,    83,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,   202,   138,    43,    52,
     336,    43,   146,   147,   274,   336,    52,   274,    52,    83,
      52,    52,   205,   436,   437,    52,    52,    83,    83,   434,
     328,    52,    52,   336,    52,   331,    52,   200,   201,    83,
      43,    52,    34,    52,   265,   200,   200,   200,   272,    83,
     200,   200,   272,    83,   228,   437,    52,    52,    52,    52,
      52,   336,   336,   336,    52,    52,    52,    52,    83,   201,
     336,   291,   200,   272,    83,    34,   138,     6,    43,    45,
      52,    53,    83,    93,   128,   146,   279,   287,   288,   138,
     288,   138,   138,   288,   138,    52,   146,   147,   273,    83,
     200,    83,    32,   247,   249,    34,   200,    46,    58,    64,
      84,   238,   239,   352,   353,   245,   200,   200,    57,    81,
     316,    83,   148,   204,   208,   209,   317,   318,   319,   138,
      34,   138,   200,   228,   228,   228,   228,   229,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   202,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,    83,   200,   138,   228,    52,   336,    52,    52,    52,
      52,    52,    52,   336,    52,    52,    52,   200,   272,   127,
     248,   248,   248,   273,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   202,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   200,   294,   200,   272,
     200,   272,   228,   200,   206,    43,    52,   138,   201,   268,
     200,   262,   200,   271,   200,   272,   336,   325,   200,   272,
     127,   127,   127,    52,    52,    52,    52,    52,    52,    52,
      52,    52,    52,    52,    52,    52,    52,    52,   336,   336,
      52,   437,   336,   336,    52,    52,   336,    52,   336,   336,
     200,   334,    43,    43,    52,   435,   206,   435,   204,   200,
     200,    52,   310,   202,   202,   228,   200,   202,   200,   202,
     200,   207,   296,   297,   200,    83,    83,    43,    52,   200,
     138,   138,    83,   138,   288,    83,   200,   288,    52,    52,
     202,   233,    34,   248,    34,   138,   202,   200,   243,   242,
     138,   200,   201,   319,    83,   228,    83,   101,   121,   202,
     138,   138,   138,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   228,    83,   200,   200,   202,
     138,   138,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   200,   200,   202,   259,   200,    43,    52,
     201,   228,   200,   200,   204,    83,   202,    34,   200,   200,
     272,   200,   272,    83,   138,   202,   280,   288,   287,   288,
     138,   288,   138,   138,   200,    34,    32,   248,   200,   336,
     239,   244,   246,   246,   246,   318,   288,    34,   200,    34,
      52,   255,   228,   228,   228,   228,   200,   200,   228,   248,
     248,   228,   202,    52,   310,   228,   200,   200,   207,   296,
     138,   138,   138,   288,   200,   288,   288,   228,   200,   200,
      32,    32,    32,   201,   202,   228,   228,   204,    52,   138,
     200,   200,   202,   202,   200,   202,   202,   202,    34,   200,
     138,   288,   280,   288,   138,   200,   200,   200,   246,   288,
     200,   200,    52,   204,    52,   131,   228,   207,   288,   138,
     138,   288,    32,   202,    52,   204,   228,   256,   200,    83,
     288,   287,   200,    52,   200,   228,   207,   138,   138,   288,
     280,   138,   288
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
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
      59,    40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   210,   211,   211,   212,   212,   212,   212,   212,   212,
     212,   212,   212,   212,   212,   212,   212,   212,   212,   212,
     212,   212,   212,   212,   212,   212,   212,   212,   212,   212,
     212,   212,   212,   212,   212,   212,   212,   212,   212,   212,
     212,   212,   212,   212,   213,   213,   213,   213,   214,   214,
     215,   216,   217,   218,   219,   220,   220,   220,   220,   220,
     220,   221,   221,   221,   221,   221,   221,   222,   222,   222,
     222,   222,   222,   223,   223,   223,   223,   223,   223,   224,
     224,   225,   225,   226,   226,   227,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   229,   229,   230,
     230,   231,   232,   233,   233,   234,   235,   236,   236,   237,
     238,   238,   239,   239,   239,   239,   239,   241,   240,   242,
     240,   243,   240,   244,   240,   245,   240,   246,   246,   246,
     246,   247,   247,   248,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   249,   250,   250,   251,   252,   253,   253,
     254,   254,   254,   254,   254,   255,   255,   255,   255,   255,
     255,   256,   256,   257,   258,   258,   259,   259,   259,   259,
     259,   259,   259,   259,   259,   260,   260,   261,   261,   262,
     263,   263,   264,   264,   265,   266,   266,   267,   267,   268,
     268,   269,   269,   269,   269,   270,   270,   271,   271,   271,
     271,   271,   271,   271,   271,   271,   271,   271,   271,   271,
     271,   271,   271,   271,   271,   271,   271,   271,   271,   271,
     272,   272,   272,   272,   272,   272,   273,   273,   273,   274,
     274,   274,   275,   276,   276,   277,   278,   278,   278,   279,
     279,   279,   279,   279,   280,   280,   280,   280,   281,   282,
     282,   283,   283,   283,   284,   285,   285,   286,   286,   286,
     287,   287,   287,   287,   287,   288,   288,   288,   288,   288,
     288,   289,   289,   289,   289,   290,   290,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   292,   292,   293,   293,
     294,   294,   294,   294,   294,   294,   294,   294,   294,   294,
     294,   294,   294,   295,   295,   296,   296,   297,   297,   298,
     299,   300,   300,   301,   302,   303,   304,   304,   304,   304,
     305,   306,   306,   306,   306,   307,   308,   308,   309,   309,
     309,   310,   310,   310,   311,   311,   312,   312,   312,   312,
     312,   312,   313,   313,   313,   313,   313,   313,   314,   315,
     315,   316,   316,   316,   317,   317,   317,   317,   318,   318,
     319,   319,   319,   319,   319,   321,   322,   320,   323,   323,
     323,   323,   324,   324,   325,   325,   326,   326,   326,   326,
     326,   326,   326,   327,   327,   327,   327,   327,   327,   327,
     327,   327,   327,   328,   328,   329,   329,   330,   330,   330,
     330,   331,   331,   332,   332,   333,   333,   334,   334,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   336,   336,   337,   338,
     339,   340,   341,   342,   343,   344,   345,   346,   347,   348,
     349,   350,   351,   352,   353,   354,   355,   356,   357,   358,
     358,   359,   360,   361,   362,   363,   364,   364,   365,   366,
     367,   368,   369,   370,   371,   372,   373,   374,   375,   376,
     377,   378,   379,   380,   381,   382,   383,   384,   385,   386,
     387,   388,   389,   389,   390,   391,   392,   393,   394,   395,
     396,   397,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,   412,   413,   414,   415,
     416,   417,   418,   419,   420,   421,   422,   423,   424,   425,
     426,   427,   428,   429,   430,   431,   432,   433,   434,   435,
     435,   436,   436,   437
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  parser::yyr2_[] =
  {
         0,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     4,
       3,     3,     3,     3,     3,     2,     3,     1,     3,     4,
       2,     2,     3,     1,     3,     4,     2,     2,     3,     1,
       3,     4,     2,     2,     3,     1,     3,     4,     2,     3,
       4,     3,     4,     3,     4,     4,     3,     1,     1,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     2,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     6,     6,     4,     1,     3,     4,
       7,     3,     4,     2,     1,     4,     4,     2,     1,     7,
       3,     1,     1,     1,     1,     1,     1,     0,     5,     0,
       8,     0,     8,     0,    10,     0,     8,     2,     2,     1,
       1,     4,     2,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     2,     2,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     6,     6,     5,     1,     4,     4,     4,     2,     1,
       9,     6,     5,     7,     7,     2,     4,     3,     5,     3,
       1,     2,     1,     6,     3,     1,     5,     3,     3,     4,
       2,     2,     3,     1,     1,     2,     5,     3,     1,     1,
       2,     5,     3,     1,     1,     2,     5,     3,     1,     1,
       1,     2,     5,     3,     6,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       2,     4,     3,     5,     1,     3,     2,     2,     1,     2,
       2,     1,     4,     2,     1,     4,     2,     1,     4,     3,
       5,     9,     1,     5,     3,     5,     7,     9,     4,     2,
       1,     5,     7,     4,     4,     2,     1,     7,     9,     6,
       1,     1,     1,     1,     1,     0,     1,     1,     1,     2,
       2,     2,     5,     3,     6,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     5,     6,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     5,     6,     7,     5,     1,     3,     3,
       4,     2,     1,     5,     3,     4,     4,     6,     3,     5,
       3,     2,     5,     3,     6,     4,     2,     1,     5,     7,
       9,     0,     3,     3,     2,     5,     5,     6,     3,     7,
       8,     5,     5,     6,     3,     7,     8,     5,     6,     3,
       1,     1,     1,     1,     1,     3,     4,     6,     1,     2,
       1,     1,     1,     1,     1,     0,     0,     5,     2,     5,
       3,     6,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     3,     1,     3,     6,     1,     1,     1,
       1,     3,     1,     3,     6,     2,     5,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     3,
       3,     1,     3,     3,     3,     3,     1,     1,     1,     3,
       3,     3,     3,     3,     3,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     1,     3,     3,
       3,     3,     5,     3,     3,     3,     1,     3,     3,     3,
       1,     1,     1,     1,     1,     3,     1,     1,     1,     1,
       3,     3,     3,     3,     1,     1,     3,     3,     3,     1,
       1,     1,     3,     3,     3,     3,     3,     3,     1,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       3,     2,     2,     2
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const parser::yytname_[] =
  {
    "$end", "error", "$undefined", "AR", "AUTOCORR", "BAYESIAN_IRF",
  "BETA_PDF", "BVAR_DENSITY", "BVAR_FORECAST", "BVAR_PRIOR_DECAY",
  "BVAR_PRIOR_FLAT", "BVAR_PRIOR_LAMBDA", "BVAR_PRIOR_MU",
  "BVAR_PRIOR_OMEGA", "BVAR_PRIOR_TAU", "BVAR_PRIOR_TRAIN", "BVAR_REPLIC",
  "CALIB", "CALIB_VAR", "CHECK", "CONF_SIG", "CONSTANT", "CORR", "COVAR",
  "CUTOFF", "DATAFILE", "DR_ALGO", "DROP", "DSAMPLE", "DUMMY", "DYNASAVE",
  "DYNATYPE", "END", "ENDVAL", "EQUAL", "ESTIMATION", "ESTIMATED_PARAMS",
  "ESTIMATED_PARAMS_BOUNDS", "ESTIMATED_PARAMS_INIT", "FILENAME",
  "FILTER_STEP_AHEAD", "FILTERED_VARS", "FIRST_OBS", "FLOAT_NUMBER",
  "FORECAST", "GAMMA_PDF", "GCC_COMPILER", "GRAPH", "HISTVAL", "HP_FILTER",
  "HP_NGRID", "INITVAL", "INT_NUMBER", "INV_GAMMA_PDF", "IRF",
  "KALMAN_ALGO", "KALMAN_TOL", "LAPLACE", "LCC_COMPILER", "LIK_ALGO",
  "LIK_INIT", "LINEAR", "LOAD_MH_FILE", "LOGLINEAR", "MARKOWITZ", "MAX",
  "MH_DROP", "MH_INIT_SCALE", "MH_JSCALE", "MH_MODE", "MH_NBLOCKS",
  "MH_REPLIC", "MH_RECOVER", "MIN", "MODE_CHECK", "MODE_COMPUTE",
  "MODE_FILE", "MODEL", "MODEL_COMPARISON", "MSHOCKS",
  "MODEL_COMPARISON_APPROXIMATION", "MODIFIEDHARMONICMEAN",
  "MOMENTS_VARENDO", "NAME", "NO_COMPILER", "NOBS", "NOCONSTANT", "NOCORR",
  "NODIAGNOSTIC", "NOFUNCTIONS", "NOGRAPH", "NOMOMENTS", "NOPRINT",
  "NORMAL_PDF", "OBSERVATION_TRENDS", "OPTIM", "OPTIM_WEIGHTS", "ORDER",
  "OSR", "OSR_PARAMS", "PARAMETERS", "PERIODS", "PLANNER_OBJECTIVE",
  "PREFILTER", "PRESAMPLE", "PRINT", "PRIOR_TRUNC", "PRIOR_ANALYSIS",
  "POSTERIOR_ANALYSIS", "QZ_CRITERIUM", "RELATIVE_IRF", "REPLIC", "RPLOT",
  "SHOCKS", "SIGMA_E", "SIMUL", "SIMUL_ALGO", "SIMUL_SEED", "SMOOTHER",
  "SOLVE_ALGO", "SPARSE_DLL", "STDERR", "STEADY", "STOCH_SIMUL", "TEX",
  "RAMSEY_POLICY", "PLANNER_DISCOUNT", "TEX_NAME", "UNIFORM_PDF",
  "UNIT_ROOT_VARS", "USE_DLL", "VALUES", "VAR", "VAREXO", "VAREXO_DET",
  "VAROBS", "XLS_SHEET", "XLS_RANGE", "COMMA", "GREATER", "LESS",
  "EXCLAMATION_EQUAL", "EXCLAMATION", "EQUAL_EQUAL", "GREATER_EQUAL",
  "LESS_EQUAL", "MINUS", "PLUS", "DIVIDE", "TIMES", "UMINUS", "POWER",
  "EXP", "LOG", "LOG10", "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN",
  "SINH", "COSH", "TANH", "ASINH", "ACOSH", "ATANH", "SQRT",
  "DYNARE_SENSITIVITY", "IDENTIFICATION", "MORRIS", "STAB", "REDFORM",
  "PPRIOR", "PRIOR_RANGE", "PPOST", "ILPTAU", "GLUE", "MORRIS_NLIV",
  "MORRIS_NTRA", "NSAM", "LOAD_REDFORM", "LOAD_RMSE", "LOAD_STAB",
  "ALPHA2_STAB", "KSSTAT", "LOGTRANS_REDFORM", "THRESHOLD_REDFORM",
  "KSSTAT_REDFORM", "ALPHA2_REDFORM", "NAMENDO", "NAMLAGENDO", "NAMEXO",
  "RMSE", "LIK_ONLY", "VAR_RMSE", "PFILT_RMSE", "ISTART_RMSE",
  "ALPHA_RMSE", "ALPHA2_RMSE", "';'", "'('", "')'", "'#'", "':'", "'['",
  "']'", "'''", "'.'", "'\\\\'", "$accept", "statement_list", "statement",
  "declaration", "dsample", "rplot", "var", "varexo", "varexo_det",
  "parameters", "var_list", "varexo_list", "varexo_det_list",
  "parameter_list", "periods", "cutoff", "markowitz", "init_param",
  "expression", "comma_expression", "initval", "initval_option", "endval",
  "initval_list", "initval_elem", "histval", "histval_list",
  "histval_elem", "model_sparse_options_list", "model_sparse_options",
  "model", "@1", "@2", "@3", "@4", "@5", "equation_list", "equation",
  "hand_side", "pound_expression", "model_var", "shocks", "mshocks",
  "shock_list", "shock_elem", "period_list", "value_list", "sigma_e",
  "triangular_matrix", "triangular_row", "steady", "steady_options_list",
  "steady_options", "check", "check_options_list", "check_options",
  "simul", "simul_options_list", "simul_options", "stoch_simul",
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
  "osr_params", "osr", "calib_var", "calib_var_list", "calib_arg1",
  "calib_arg2", "calib", "dynatype", "dynasave", "model_comparison",
  "model_comparison_options", "model_comparison_option", "filename_list",
  "filename", "filename_elem", "planner_objective", "@6", "@7",
  "ramsey_policy", "ramsey_policy_options_list", "ramsey_policy_options",
  "bvar_prior_option", "bvar_common_option", "bvar_density_options_list",
  "bvar_density", "bvar_forecast_option", "bvar_forecast_options_list",
  "bvar_forecast", "dynare_sensitivity", "dynare_sensitivity_options_list",
  "dynare_sensitivity_option", "number", "o_dr_algo", "o_solve_algo",
  "o_simul_algo", "o_linear", "o_order", "o_replic", "o_drop", "o_ar",
  "o_nocorr", "o_nofunctions", "o_nomoments", "o_irf", "o_hp_filter",
  "o_hp_ngrid", "o_periods", "o_cutoff", "o_markowitz", "o_simul",
  "o_simul_seed", "o_qz_criterium", "o_datafile", "o_nobs", "o_first_obs",
  "o_prefilter", "o_presample", "o_lik_algo", "o_lik_init", "o_nograph",
  "o_conf_sig", "o_mh_replic", "o_mh_drop", "o_mh_jscale", "o_optim",
  "o_mh_init_scale", "o_mode_file", "o_mode_compute", "o_mode_check",
  "o_prior_trunc", "o_mh_mode", "o_mh_nblcks", "o_load_mh_file",
  "o_loglinear", "o_nodiagnostic", "o_bayesian_irf", "o_tex", "o_forecast",
  "o_smoother", "o_moments_varendo", "o_filtered_vars", "o_relative_irf",
  "o_kalman_algo", "o_kalman_tol", "o_model_comparison_approximation",
  "o_print", "o_noprint", "o_xls_sheet", "o_xls_range",
  "o_filter_step_ahead", "o_constant", "o_noconstant", "o_mh_recover",
  "o_planner_discount", "o_bvar_prior_tau", "o_bvar_prior_decay",
  "o_bvar_prior_lambda", "o_bvar_prior_mu", "o_bvar_prior_omega",
  "o_bvar_prior_flat", "o_bvar_prior_train", "o_bvar_replic",
  "o_gsa_identification", "o_gsa_morris", "o_gsa_stab", "o_gsa_redform",
  "o_gsa_pprior", "o_gsa_prior_range", "o_gsa_ppost", "o_gsa_ilptau",
  "o_gsa_glue", "o_gsa_morris_nliv", "o_gsa_morris_ntra", "o_gsa_nsam",
  "o_gsa_load_redform", "o_gsa_load_rmse", "o_gsa_load_stab",
  "o_gsa_alpha2_stab", "o_gsa_ksstat", "o_gsa_logtrans_redform",
  "o_gsa_threshold_redform", "o_gsa_ksstat_redform",
  "o_gsa_alpha2_redform", "o_gsa_rmse", "o_gsa_lik_only",
  "o_gsa_pfilt_rmse", "o_gsa_istart_rmse", "o_gsa_alpha_rmse",
  "o_gsa_alpha2_rmse", "range", "vec_int_elem", "vec_int_1", "vec_int", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const parser::rhs_number_type
  parser::yyrhs_[] =
  {
       211,     0,    -1,   212,    -1,   211,   212,    -1,   213,    -1,
     224,    -1,   225,    -1,   226,    -1,   240,    -1,   230,    -1,
     232,    -1,   235,    -1,   227,    -1,   251,    -1,   252,    -1,
     257,    -1,   260,    -1,   263,    -1,   266,    -1,   269,    -1,
     289,    -1,   292,    -1,   295,    -1,   275,    -1,   284,    -1,
     281,    -1,   298,    -1,   299,    -1,   302,    -1,   214,    -1,
     215,    -1,   303,    -1,   305,    -1,   306,    -1,   307,    -1,
     311,    -1,   312,    -1,   313,    -1,   314,    -1,   320,    -1,
     323,    -1,   329,    -1,   332,    -1,   333,    -1,   219,    -1,
     216,    -1,   217,    -1,   218,    -1,    28,    52,   200,    -1,
      28,    52,    52,   200,    -1,   112,   272,   200,    -1,   132,
     220,   200,    -1,   133,   221,   200,    -1,   134,   222,   200,
      -1,   100,   223,   200,    -1,   220,    83,    -1,   220,   138,
      83,    -1,    83,    -1,   220,    83,   127,    -1,   220,   138,
      83,   127,    -1,    83,   127,    -1,   221,    83,    -1,   221,
     138,    83,    -1,    83,    -1,   221,    83,   127,    -1,   221,
     138,    83,   127,    -1,    83,   127,    -1,   222,    83,    -1,
     222,   138,    83,    -1,    83,    -1,   222,    83,   127,    -1,
     222,   138,    83,   127,    -1,    83,   127,    -1,   223,    83,
      -1,   223,   138,    83,    -1,    83,    -1,   223,    83,   127,
      -1,   223,   138,    83,   127,    -1,    83,   127,    -1,   101,
      52,   200,    -1,   101,    34,    52,   200,    -1,    24,    43,
     200,    -1,    24,    34,    43,   200,    -1,    64,    43,   200,
      -1,    64,    34,    43,   200,    -1,    83,    34,   228,   200,
      -1,   201,   228,   202,    -1,    83,    -1,    43,    -1,    52,
      -1,   228,   147,   228,    -1,   228,   146,   228,    -1,   228,
     148,   228,    -1,   228,   149,   228,    -1,   228,   151,   228,
      -1,   228,   140,   228,    -1,   228,   139,   228,    -1,   228,
     145,   228,    -1,   228,   144,   228,    -1,   228,   143,   228,
      -1,   228,   141,   228,    -1,   146,   228,    -1,   147,   228,
      -1,   152,   201,   228,   202,    -1,   153,   201,   228,   202,
      -1,   154,   201,   228,   202,    -1,   155,   201,   228,   202,
      -1,   156,   201,   228,   202,    -1,   157,   201,   228,   202,
      -1,   158,   201,   228,   202,    -1,   159,   201,   228,   202,
      -1,   160,   201,   228,   202,    -1,   167,   201,   228,   202,
      -1,    29,   201,   228,   202,    -1,    65,   201,   228,   138,
     228,   202,    -1,    73,   201,   228,   138,   228,   202,    -1,
      83,   201,   229,   202,    -1,   228,    -1,   229,   138,   228,
      -1,    51,   200,   233,    32,    -1,    51,   201,   231,   202,
     200,   233,    32,    -1,    39,    34,    83,    -1,    33,   200,
     233,    32,    -1,   233,   234,    -1,   234,    -1,    83,    34,
     228,   200,    -1,    48,   200,   236,    32,    -1,   236,   237,
      -1,   237,    -1,    83,   201,   273,   202,    34,   228,   200,
      -1,   238,   138,   239,    -1,   239,    -1,    58,    -1,    46,
      -1,    84,    -1,   352,    -1,   353,    -1,    -1,    77,   200,
     241,   246,    32,    -1,    -1,    77,   201,   340,   202,   200,
     242,   246,    32,    -1,    -1,    77,   201,   130,   202,   200,
     243,   246,    32,    -1,    -1,    77,   201,   120,   138,   238,
     202,   244,   200,   246,    32,    -1,    -1,    77,   201,   120,
     202,   245,   200,   246,    32,    -1,   246,   247,    -1,   246,
     249,    -1,   247,    -1,   249,    -1,   248,    34,   248,   200,
      -1,   248,   200,    -1,   201,   248,   202,    -1,   250,    -1,
      43,    -1,    52,    -1,   248,   147,   248,    -1,   248,   146,
     248,    -1,   248,   148,   248,    -1,   248,   149,   248,    -1,
     248,   140,   248,    -1,   248,   139,   248,    -1,   248,   145,
     248,    -1,   248,   144,   248,    -1,   248,   143,   248,    -1,
     248,   141,   248,    -1,   248,   151,   248,    -1,   146,   248,
      -1,   147,   248,    -1,   152,   201,   248,   202,    -1,   153,
     201,   248,   202,    -1,   154,   201,   248,   202,    -1,   155,
     201,   248,   202,    -1,   156,   201,   248,   202,    -1,   157,
     201,   248,   202,    -1,   158,   201,   248,   202,    -1,   159,
     201,   248,   202,    -1,   160,   201,   248,   202,    -1,   167,
     201,   248,   202,    -1,    29,   201,   248,   202,    -1,    65,
     201,   248,   138,   248,   202,    -1,    73,   201,   248,   138,
     248,   202,    -1,   203,    83,    34,   248,   200,    -1,    83,
      -1,    83,   201,   273,   202,    -1,   113,   200,   253,    32,
      -1,    79,   200,   253,    32,    -1,   253,   254,    -1,   254,
      -1,   132,    83,   200,   101,   255,   200,   131,   256,   200,
      -1,   132,    83,   200,   121,   228,   200,    -1,   132,    83,
      34,   228,   200,    -1,   132,    83,   138,    83,    34,   228,
     200,    -1,    22,    83,   138,    83,    34,   228,   200,    -1,
     255,    52,    -1,   255,    52,   204,    52,    -1,   255,   138,
      52,    -1,   255,   138,    52,   204,    52,    -1,    52,   204,
      52,    -1,    52,    -1,   256,   228,    -1,   228,    -1,   114,
      34,   205,   258,   206,   200,    -1,   258,   200,   259,    -1,
     259,    -1,   259,   138,   201,   228,   202,    -1,   259,   138,
      43,    -1,   259,   138,    52,    -1,   259,   201,   228,   202,
      -1,   259,    43,    -1,   259,    52,    -1,   201,   228,   202,
      -1,    43,    -1,    52,    -1,   122,   200,    -1,   122,   201,
     261,   202,   200,    -1,   261,   138,   262,    -1,   262,    -1,
     338,    -1,    19,   200,    -1,    19,   201,   264,   202,   200,
      -1,   264,   138,   265,    -1,   265,    -1,   338,    -1,   115,
     200,    -1,   115,   201,   267,   202,   200,    -1,   267,   138,
     268,    -1,   268,    -1,   351,    -1,   357,    -1,   123,   200,
      -1,   123,   201,   270,   202,   200,    -1,   123,   272,   200,
      -1,   123,   201,   270,   202,   272,   200,    -1,   270,   138,
     271,    -1,   271,    -1,   337,    -1,   338,    -1,   339,    -1,
     340,    -1,   341,    -1,   342,    -1,   343,    -1,   344,    -1,
     345,    -1,   346,    -1,   347,    -1,   364,    -1,   348,    -1,
     386,    -1,   349,    -1,   350,    -1,   351,    -1,   352,    -1,
     354,    -1,   355,    -1,   356,    -1,   390,    -1,   391,    -1,
     272,    83,    -1,   272,    83,    34,    83,    -1,   272,   138,
      83,    -1,   272,   138,    83,    34,    83,    -1,    83,    -1,
      83,    34,    83,    -1,   147,    52,    -1,   146,    52,    -1,
      52,    -1,   147,    43,    -1,   146,    43,    -1,    43,    -1,
      36,   200,   276,    32,    -1,   276,   277,    -1,   277,    -1,
     278,   138,   279,   200,    -1,   121,    83,    -1,    83,    -1,
      22,    83,   138,    83,    -1,   287,   138,   280,    -1,   288,
     138,   287,   138,   280,    -1,   288,   138,   288,   138,   288,
     138,   287,   138,   280,    -1,   288,    -1,   288,   138,   288,
     138,   288,    -1,   288,   138,   288,    -1,   288,   138,   288,
     138,   288,    -1,   288,   138,   288,   138,   288,   138,   288,
      -1,   288,   138,   288,   138,   288,   138,   288,   138,   288,
      -1,    38,   200,   282,    32,    -1,   282,   283,    -1,   283,
      -1,   121,    83,   138,   288,   200,    -1,    22,    83,   138,
      83,   138,   288,   200,    -1,    83,   138,   288,   200,    -1,
      37,   200,   285,    32,    -1,   285,   286,    -1,   286,    -1,
     121,    83,   138,   288,   138,   288,   200,    -1,    22,    83,
     138,    83,   138,   288,   138,   288,   200,    -1,    83,   138,
     288,   138,   288,   200,    -1,     6,    -1,    45,    -1,    93,
      -1,    53,    -1,   128,    -1,    -1,    52,    -1,    43,    -1,
      83,    -1,   146,    52,    -1,   146,    43,    -1,    35,   200,
      -1,    35,   201,   290,   202,   200,    -1,    35,   272,   200,
      -1,    35,   201,   290,   202,   272,   200,    -1,   290,   138,
     291,    -1,   291,    -1,   357,    -1,   358,    -1,   359,    -1,
     360,    -1,   361,    -1,   362,    -1,   363,    -1,   364,    -1,
     365,    -1,   366,    -1,   367,    -1,   368,    -1,   369,    -1,
     370,    -1,   371,    -1,   372,    -1,   373,    -1,   374,    -1,
     375,    -1,   376,    -1,   377,    -1,   378,    -1,   379,    -1,
     380,    -1,   348,    -1,   381,    -1,   382,    -1,   383,    -1,
     384,    -1,   385,    -1,   387,    -1,   388,    -1,   392,    -1,
     393,    -1,   394,    -1,   338,    -1,   395,    -1,   396,    -1,
     397,    -1,   107,   201,   293,   202,   200,    -1,   107,   201,
     293,   202,   272,   200,    -1,   293,   138,   294,    -1,   294,
      -1,   364,    -1,   365,    -1,   374,    -1,   380,    -1,   348,
      -1,   381,    -1,   382,    -1,   383,    -1,   384,    -1,   385,
      -1,   392,    -1,   393,    -1,   394,    -1,   108,   201,   293,
     202,   200,    -1,   108,   201,   293,   202,   272,   200,    -1,
     207,    83,   207,   138,   207,    83,   207,    -1,   207,    83,
     207,   138,   288,    -1,   296,    -1,   297,   138,   296,    -1,
     135,   272,   200,    -1,    94,   200,   300,    32,    -1,   300,
     301,    -1,   301,    -1,    83,   201,   228,   202,   200,    -1,
     129,   272,   200,    -1,    96,   200,   304,    32,    -1,   304,
      83,   228,   200,    -1,   304,    83,   138,    83,   228,   200,
      -1,    83,   228,   200,    -1,    83,   138,    83,   228,   200,
      -1,    99,   272,   200,    -1,    98,   200,    -1,    98,   201,
     271,   202,   200,    -1,    98,   272,   200,    -1,    98,   201,
     271,   202,   272,   200,    -1,    18,   200,   308,    32,    -1,
     308,   309,    -1,   309,    -1,    83,   310,    34,   228,   200,
      -1,    83,   138,    83,   310,    34,   228,   200,    -1,     4,
      83,   201,    52,   202,   310,    34,   228,   200,    -1,    -1,
     201,    52,   202,    -1,   201,    43,   202,    -1,    17,   200,
      -1,    17,   201,    23,   202,   200,    -1,    31,   201,    83,
     202,   200,    -1,    31,   201,    83,   202,   272,   200,    -1,
      31,    83,   200,    -1,    31,   201,    83,   208,    83,   202,
     200,    -1,    31,   201,    83,   208,    83,   202,   272,   200,
      -1,    31,    83,   208,    83,   200,    -1,    30,   201,    83,
     202,   200,    -1,    30,   201,    83,   202,   272,   200,    -1,
      30,    83,   200,    -1,    30,   201,    83,   208,    83,   202,
     200,    -1,    30,   201,    83,   208,    83,   202,   272,   200,
      -1,    30,    83,   208,    83,   200,    -1,    78,   201,   315,
     202,   317,   200,    -1,   315,   138,   316,    -1,   316,    -1,
     389,    -1,   390,    -1,   391,    -1,   318,    -1,   317,   138,
     318,    -1,   318,   201,   288,   202,    -1,   317,   138,   318,
     201,   288,   202,    -1,   319,    -1,   318,   319,    -1,    83,
      -1,   209,    -1,   148,    -1,   204,    -1,   208,    -1,    -1,
      -1,   102,   321,   248,   322,   200,    -1,   125,   200,    -1,
     125,   201,   324,   202,   200,    -1,   125,   272,   200,    -1,
     125,   201,   324,   202,   272,   200,    -1,   324,   138,   325,
      -1,   325,    -1,   271,    -1,   398,    -1,   399,    -1,   400,
      -1,   401,    -1,   402,    -1,   403,    -1,   404,    -1,   405,
      -1,   326,    -1,   357,    -1,   392,    -1,   393,    -1,   359,
      -1,   361,    -1,   358,    -1,   360,    -1,   395,    -1,   396,
      -1,   327,   138,   328,    -1,   327,    -1,     7,    52,   200,
      -1,     7,   201,   328,   202,    52,   200,    -1,   327,    -1,
     382,    -1,   365,    -1,   406,    -1,   330,   138,   331,    -1,
     330,    -1,     8,    52,   200,    -1,     8,   201,   331,   202,
      52,   200,    -1,   168,   200,    -1,   168,   201,   334,   202,
     200,    -1,   335,   138,   334,    -1,   335,    -1,   407,    -1,
     408,    -1,   409,    -1,   410,    -1,   411,    -1,   412,    -1,
     413,    -1,   414,    -1,   415,    -1,   416,    -1,   417,    -1,
     418,    -1,   419,    -1,   420,    -1,   421,    -1,   422,    -1,
     423,    -1,   424,    -1,   426,    -1,   427,    -1,   428,    -1,
     429,    -1,   430,    -1,   431,    -1,   432,    -1,   433,    -1,
     425,    -1,    52,    -1,    43,    -1,    26,    34,    52,    -1,
     119,    34,    52,    -1,   116,    34,    52,    -1,    61,    -1,
      97,    34,    52,    -1,   111,    34,    52,    -1,    27,    34,
      52,    -1,     3,    34,    52,    -1,    87,    -1,    89,    -1,
      91,    -1,    54,    34,    52,    -1,    49,    34,    52,    -1,
      50,    34,    52,    -1,   101,    34,    52,    -1,    24,    34,
     336,    -1,    64,    34,   336,    -1,   115,    -1,   117,    34,
      52,    -1,   109,    34,   336,    -1,    25,    34,    83,    -1,
      85,    34,   437,    -1,    85,    34,    52,    -1,    42,    34,
      52,    -1,   103,    34,    52,    -1,   104,    34,    52,    -1,
      59,    34,    52,    -1,    60,    34,    52,    -1,    90,    -1,
      47,    -1,    20,    34,   336,    -1,    71,    34,    52,    -1,
      66,    34,   336,    -1,    68,    34,   336,    -1,    95,    34,
     201,   297,   202,    -1,    67,    34,   336,    -1,    76,    34,
      83,    -1,    75,    34,    52,    -1,    74,    -1,   106,    34,
     336,    -1,    69,    34,    52,    -1,    70,    34,    52,    -1,
      62,    -1,    63,    -1,    88,    -1,     5,    -1,   124,    -1,
      44,    34,    52,    -1,   118,    -1,    82,    -1,    41,    -1,
     110,    -1,    55,    34,    52,    -1,    56,    34,    52,    -1,
      80,    34,    57,    -1,    80,    34,    81,    -1,   105,    -1,
      92,    -1,   136,    34,    83,    -1,   137,    34,   434,    -1,
      40,    34,   437,    -1,    21,    -1,    86,    -1,    72,    -1,
     126,    34,   336,    -1,    14,    34,   274,    -1,     9,    34,
     336,    -1,    11,    34,   274,    -1,    12,    34,   336,    -1,
      13,    34,    52,    -1,    10,    -1,    15,    34,    52,    -1,
      16,    34,    52,    -1,   169,    34,    52,    -1,   170,    34,
      52,    -1,   171,    34,    52,    -1,   172,    34,    52,    -1,
     173,    34,    52,    -1,   174,    34,    52,    -1,   175,    34,
      52,    -1,   176,    34,    52,    -1,   177,    34,    52,    -1,
     178,    34,    52,    -1,   179,    34,    52,    -1,   180,    34,
      52,    -1,   181,    34,    52,    -1,   182,    34,    52,    -1,
     183,    34,    52,    -1,   184,    34,   336,    -1,   185,    34,
     336,    -1,   186,    34,    52,    -1,   187,    34,   437,    -1,
     188,    34,   336,    -1,   189,    34,   336,    -1,   193,    34,
      52,    -1,   194,    34,    52,    -1,   196,    34,   336,    -1,
     197,    34,    52,    -1,   198,    34,   336,    -1,   199,    34,
     336,    -1,    83,   204,    83,    -1,    52,    -1,    52,   204,
      52,    -1,   205,   435,    -1,   436,   435,    -1,   436,   206,
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
      80,    82,    84,    86,    88,    90,    92,    94,    96,   100,
     105,   109,   113,   117,   121,   125,   128,   132,   134,   138,
     143,   146,   149,   153,   155,   159,   164,   167,   170,   174,
     176,   180,   185,   188,   191,   195,   197,   201,   206,   209,
     213,   218,   222,   227,   231,   236,   241,   245,   247,   249,
     251,   255,   259,   263,   267,   271,   275,   279,   283,   287,
     291,   295,   298,   301,   306,   311,   316,   321,   326,   331,
     336,   341,   346,   351,   356,   363,   370,   375,   377,   381,
     386,   394,   398,   403,   406,   408,   413,   418,   421,   423,
     431,   435,   437,   439,   441,   443,   445,   447,   448,   454,
     455,   464,   465,   474,   475,   486,   487,   496,   499,   502,
     504,   506,   511,   514,   518,   520,   522,   524,   528,   532,
     536,   540,   544,   548,   552,   556,   560,   564,   568,   571,
     574,   579,   584,   589,   594,   599,   604,   609,   614,   619,
     624,   629,   636,   643,   649,   651,   656,   661,   666,   669,
     671,   681,   688,   694,   702,   710,   713,   718,   722,   728,
     732,   734,   737,   739,   746,   750,   752,   758,   762,   766,
     771,   774,   777,   781,   783,   785,   788,   794,   798,   800,
     802,   805,   811,   815,   817,   819,   822,   828,   832,   834,
     836,   838,   841,   847,   851,   858,   862,   864,   866,   868,
     870,   872,   874,   876,   878,   880,   882,   884,   886,   888,
     890,   892,   894,   896,   898,   900,   902,   904,   906,   908,
     910,   913,   918,   922,   928,   930,   934,   937,   940,   942,
     945,   948,   950,   955,   958,   960,   965,   968,   970,   975,
     979,   985,   995,   997,  1003,  1007,  1013,  1021,  1031,  1036,
    1039,  1041,  1047,  1055,  1060,  1065,  1068,  1070,  1078,  1088,
    1095,  1097,  1099,  1101,  1103,  1105,  1106,  1108,  1110,  1112,
    1115,  1118,  1121,  1127,  1131,  1138,  1142,  1144,  1146,  1148,
    1150,  1152,  1154,  1156,  1158,  1160,  1162,  1164,  1166,  1168,
    1170,  1172,  1174,  1176,  1178,  1180,  1182,  1184,  1186,  1188,
    1190,  1192,  1194,  1196,  1198,  1200,  1202,  1204,  1206,  1208,
    1210,  1212,  1214,  1216,  1218,  1220,  1222,  1228,  1235,  1239,
    1241,  1243,  1245,  1247,  1249,  1251,  1253,  1255,  1257,  1259,
    1261,  1263,  1265,  1267,  1273,  1280,  1288,  1294,  1296,  1300,
    1304,  1309,  1312,  1314,  1320,  1324,  1329,  1334,  1341,  1345,
    1351,  1355,  1358,  1364,  1368,  1375,  1380,  1383,  1385,  1391,
    1399,  1409,  1410,  1414,  1418,  1421,  1427,  1433,  1440,  1444,
    1452,  1461,  1467,  1473,  1480,  1484,  1492,  1501,  1507,  1514,
    1518,  1520,  1522,  1524,  1526,  1528,  1532,  1537,  1544,  1546,
    1549,  1551,  1553,  1555,  1557,  1559,  1560,  1561,  1567,  1570,
    1576,  1580,  1587,  1591,  1593,  1595,  1597,  1599,  1601,  1603,
    1605,  1607,  1609,  1611,  1613,  1615,  1617,  1619,  1621,  1623,
    1625,  1627,  1629,  1631,  1635,  1637,  1641,  1648,  1650,  1652,
    1654,  1656,  1660,  1662,  1666,  1673,  1676,  1682,  1686,  1688,
    1690,  1692,  1694,  1696,  1698,  1700,  1702,  1704,  1706,  1708,
    1710,  1712,  1714,  1716,  1718,  1720,  1722,  1724,  1726,  1728,
    1730,  1732,  1734,  1736,  1738,  1740,  1742,  1744,  1746,  1750,
    1754,  1758,  1760,  1764,  1768,  1772,  1776,  1778,  1780,  1782,
    1786,  1790,  1794,  1798,  1802,  1806,  1808,  1812,  1816,  1820,
    1824,  1828,  1832,  1836,  1840,  1844,  1848,  1850,  1852,  1856,
    1860,  1864,  1868,  1874,  1878,  1882,  1886,  1888,  1892,  1896,
    1900,  1902,  1904,  1906,  1908,  1910,  1914,  1916,  1918,  1920,
    1922,  1926,  1930,  1934,  1938,  1940,  1942,  1946,  1950,  1954,
    1956,  1958,  1960,  1964,  1968,  1972,  1976,  1980,  1984,  1986,
    1990,  1994,  1998,  2002,  2006,  2010,  2014,  2018,  2022,  2026,
    2030,  2034,  2038,  2042,  2046,  2050,  2054,  2058,  2062,  2066,
    2070,  2074,  2078,  2082,  2086,  2090,  2094,  2098,  2102,  2106,
    2108,  2112,  2115,  2118
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    98,    98,    99,   102,   103,   104,   105,   106,   107,
     108,   109,   110,   111,   112,   113,   114,   115,   116,   117,
     118,   119,   120,   121,   122,   123,   124,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   144,   145,   146,   147,   151,   153,
     157,   159,   161,   163,   165,   167,   169,   171,   173,   175,
     177,   181,   183,   185,   187,   189,   191,   195,   197,   199,
     201,   203,   205,   209,   211,   213,   215,   217,   219,   223,
     225,   229,   231,   235,   237,   242,   244,   246,   248,   250,
     252,   254,   256,   258,   260,   262,   264,   266,   268,   270,
     272,   274,   276,   278,   280,   282,   284,   286,   288,   290,
     292,   294,   296,   298,   300,   302,   304,   308,   310,   314,
     316,   320,   322,   324,   325,   328,   330,   332,   333,   336,
     338,   339,   342,   344,   346,   348,   349,   352,   352,   354,
     354,   356,   356,   359,   358,   361,   361,   365,   366,   367,
     368,   371,   373,   377,   379,   380,   382,   384,   386,   388,
     390,   392,   394,   396,   398,   400,   402,   404,   406,   408,
     410,   412,   414,   416,   418,   420,   422,   424,   426,   428,
     430,   432,   434,   438,   441,   443,   447,   449,   451,   452,
     455,   457,   459,   461,   463,   467,   469,   471,   473,   475,
     477,   481,   483,   487,   489,   491,   495,   497,   499,   501,
     503,   505,   507,   509,   511,   515,   517,   521,   522,   525,
     527,   529,   533,   534,   537,   539,   541,   545,   546,   549,
     550,   553,   555,   557,   559,   563,   564,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,   582,   583,   584,   585,   586,   587,   588,   589,
     592,   594,   596,   598,   600,   602,   606,   608,   610,   614,
     616,   618,   622,   624,   626,   630,   632,   638,   644,   654,
     659,   666,   677,   682,   693,   700,   709,   720,   735,   738,
     740,   744,   752,   762,   772,   775,   777,   781,   791,   803,
     815,   817,   819,   821,   823,   827,   828,   829,   830,   831,
     833,   837,   839,   841,   843,   847,   848,   851,   852,   853,
     854,   855,   856,   857,   858,   859,   860,   861,   862,   863,
     864,   865,   866,   867,   868,   869,   870,   871,   872,   873,
     874,   875,   876,   877,   878,   879,   880,   881,   882,   883,
     884,   885,   886,   887,   888,   889,   892,   894,   898,   899,
     902,   903,   904,   905,   906,   907,   908,   909,   910,   911,
     912,   913,   914,   917,   919,   923,   925,   929,   930,   933,
     935,   937,   938,   941,   943,   945,   947,   949,   951,   953,
     957,   959,   961,   963,   965,   969,   971,   972,   975,   977,
     979,   983,   984,   986,   990,   992,   996,   998,  1000,  1002,
    1004,  1006,  1010,  1012,  1014,  1016,  1018,  1020,  1024,  1027,
    1028,  1031,  1032,  1033,  1036,  1038,  1040,  1042,  1046,  1048,
    1052,  1053,  1055,  1057,  1059,  1063,  1064,  1063,  1066,  1068,
    1070,  1072,  1076,  1077,  1080,  1081,  1084,  1085,  1086,  1087,
    1088,  1089,  1090,  1093,  1094,  1095,  1096,  1097,  1098,  1099,
    1100,  1101,  1102,  1105,  1106,  1109,  1111,  1115,  1116,  1117,
    1118,  1121,  1122,  1125,  1127,  1131,  1133,  1137,  1138,  1141,
    1142,  1143,  1144,  1145,  1146,  1147,  1148,  1149,  1150,  1151,
    1152,  1153,  1154,  1155,  1156,  1157,  1158,  1159,  1160,  1161,
    1162,  1163,  1164,  1165,  1166,  1167,  1170,  1171,  1174,  1175,
    1176,  1177,  1178,  1179,  1180,  1181,  1182,  1183,  1184,  1185,
    1186,  1187,  1188,  1190,  1191,  1192,  1193,  1194,  1195,  1196,
    1198,  1201,  1202,  1203,  1204,  1205,  1206,  1208,  1211,  1212,
    1213,  1214,  1215,  1216,  1217,  1218,  1219,  1220,  1221,  1222,
    1223,  1224,  1225,  1226,  1227,  1228,  1229,  1230,  1231,  1232,
    1233,  1234,  1235,  1237,  1240,  1241,  1242,  1243,  1244,  1245,
    1246,  1247,  1248,  1250,  1251,  1252,  1253,  1254,  1255,  1256,
    1257,  1259,  1260,  1261,  1262,  1263,  1264,  1265,  1266,  1267,
    1268,  1269,  1270,  1271,  1272,  1273,  1274,  1275,  1276,  1277,
    1279,  1280,  1286,  1287,  1291,  1292,  1293,  1294,  1297,  1305,
    1306,  1315,  1317,  1326
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
       2,     2,     2,     2,     2,   203,     2,     2,     2,   207,
     201,   202,     2,     2,     2,     2,   208,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   204,   200,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   205,   209,   206,     2,     2,     2,     2,     2,     2,
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
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 2299;
  const int parser::yynnts_ = 228;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 160;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 210;

  const unsigned int parser::yyuser_token_number_max_ = 454;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1328 "DynareBison.yy"


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

