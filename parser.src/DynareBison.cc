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
#line 51 "DynareBison.yy"

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
    #line 37 "DynareBison.yy"
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
	  case 49:
#line 173 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 50:
#line 175 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 51:
#line 178 "DynareBison.yy"
    { driver.rplot(); ;}
    break;

  case 56:
#line 189 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 57:
#line 191 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 58:
#line 193 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 59:
#line 195 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 60:
#line 197 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 61:
#line 199 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 62:
#line 203 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 63:
#line 205 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 64:
#line 207 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 65:
#line 209 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 66:
#line 211 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 67:
#line 213 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 68:
#line 217 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 69:
#line 219 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 70:
#line 221 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 71:
#line 223 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 72:
#line 225 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 73:
#line 227 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 74:
#line 231 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 75:
#line 233 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 76:
#line 235 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 77:
#line 237 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 78:
#line 239 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 79:
#line 241 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 80:
#line 245 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 81:
#line 247 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 82:
#line 251 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 83:
#line 253 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 84:
#line 257 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 85:
#line 259 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 86:
#line 263 "DynareBison.yy"
    { driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 87:
#line 266 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 88:
#line 268 "DynareBison.yy"
    { (yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 89:
#line 270 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 90:
#line 272 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 91:
#line 274 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 92:
#line 276 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 93:
#line 278 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 94:
#line 280 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 95:
#line 282 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 96:
#line 284 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 97:
#line 286 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 98:
#line 288 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 99:
#line 290 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 100:
#line 292 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 101:
#line 294 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 102:
#line 296 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 103:
#line 298 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 104:
#line 300 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 105:
#line 302 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 106:
#line 304 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 107:
#line 306 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 108:
#line 308 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 109:
#line 310 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 110:
#line 312 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 111:
#line 314 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 112:
#line 316 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 113:
#line 318 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 114:
#line 320 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 115:
#line 322 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 116:
#line 324 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 117:
#line 326 "DynareBison.yy"
    { (yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val)); ;}
    break;

  case 118:
#line 328 "DynareBison.yy"
    { (yyval.node_val) = driver.add_normcdf((yysemantic_stack_[(8) - (3)].node_val),(yysemantic_stack_[(8) - (5)].node_val),(yysemantic_stack_[(8) - (7)].node_val));;}
    break;

  case 119:
#line 332 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 120:
#line 334 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 121:
#line 338 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 122:
#line 340 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 123:
#line 343 "DynareBison.yy"
    { driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 124:
#line 345 "DynareBison.yy"
    { driver.end_endval(); ;}
    break;

  case 127:
#line 351 "DynareBison.yy"
    { driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 128:
#line 353 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 131:
#line 359 "DynareBison.yy"
    { driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 138:
#line 374 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 139:
#line 376 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 140:
#line 378 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 143:
#line 385 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 144:
#line 386 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 145:
#line 387 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 146:
#line 388 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 147:
#line 389 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 148:
#line 390 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 149:
#line 392 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 150:
#line 393 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 151:
#line 394 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 152:
#line 395 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 153:
#line 397 "DynareBison.yy"
    { driver.begin_model(); driver.sparse(); ;}
    break;

  case 154:
#line 398 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 155:
#line 399 "DynareBison.yy"
    { driver.begin_model(); driver.sparse(); ;}
    break;

  case 156:
#line 400 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 161:
#line 410 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 162:
#line 412 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val)); ;}
    break;

  case 163:
#line 416 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 165:
#line 419 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 166:
#line 421 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 167:
#line 423 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 168:
#line 425 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 169:
#line 427 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 170:
#line 429 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 171:
#line 431 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 172:
#line 433 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 173:
#line 435 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 174:
#line 437 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 175:
#line 439 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 176:
#line 441 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 177:
#line 443 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 178:
#line 445 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 179:
#line 447 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 180:
#line 449 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 181:
#line 451 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 182:
#line 453 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 183:
#line 455 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 184:
#line 457 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 185:
#line 459 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 186:
#line 461 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 187:
#line 463 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 188:
#line 465 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 189:
#line 467 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 190:
#line 469 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 191:
#line 471 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 192:
#line 473 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 193:
#line 475 "DynareBison.yy"
    { (yyval.node_val) = driver.add_normcdf((yysemantic_stack_[(8) - (3)].node_val),(yysemantic_stack_[(8) - (5)].node_val),(yysemantic_stack_[(8) - (7)].node_val));;}
    break;

  case 194:
#line 479 "DynareBison.yy"
    { driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 195:
#line 482 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 196:
#line 484 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 197:
#line 487 "DynareBison.yy"
    { driver.end_shocks(); ;}
    break;

  case 198:
#line 489 "DynareBison.yy"
    { driver.end_mshocks(); ;}
    break;

  case 201:
#line 496 "DynareBison.yy"
    { driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val)); ;}
    break;

  case 202:
#line 498 "DynareBison.yy"
    { driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 203:
#line 500 "DynareBison.yy"
    { driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 204:
#line 502 "DynareBison.yy"
    { driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 205:
#line 504 "DynareBison.yy"
    { driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 206:
#line 508 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 207:
#line 510 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 208:
#line 512 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 209:
#line 514 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 210:
#line 516 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 211:
#line 518 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 212:
#line 522 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 213:
#line 526 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 214:
#line 528 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 215:
#line 530 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].node_val));;}
    break;

  case 216:
#line 534 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 217:
#line 536 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 218:
#line 540 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 219:
#line 542 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 220:
#line 544 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 221:
#line 546 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 222:
#line 548 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 223:
#line 550 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 224:
#line 552 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 225:
#line 554 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 226:
#line 556 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 227:
#line 560 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 228:
#line 562 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 234:
#line 575 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 235:
#line 577 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 239:
#line 587 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 240:
#line 589 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 246:
#line 602 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 247:
#line 604 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 248:
#line 606 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 249:
#line 608 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 275:
#line 641 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 276:
#line 643 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 277:
#line 645 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 278:
#line 647 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 279:
#line 649 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 280:
#line 651 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 281:
#line 655 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 282:
#line 657 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 283:
#line 659 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 284:
#line 663 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 285:
#line 665 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 286:
#line 667 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 287:
#line 670 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 288:
#line 673 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 289:
#line 675 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 291:
#line 681 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 292:
#line 687 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 293:
#line 693 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 294:
#line 703 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 295:
#line 708 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 296:
#line 715 "DynareBison.yy"
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

  case 297:
#line 726 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 298:
#line 731 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 299:
#line 742 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 300:
#line 749 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 301:
#line 758 "DynareBison.yy"
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

  case 302:
#line 769 "DynareBison.yy"
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

  case 303:
#line 784 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 304:
#line 787 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 305:
#line 789 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 306:
#line 793 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 307:
#line 801 "DynareBison.yy"
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

  case 308:
#line 811 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 309:
#line 821 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 310:
#line 824 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 311:
#line 826 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 312:
#line 830 "DynareBison.yy"
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

  case 313:
#line 840 "DynareBison.yy"
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

  case 314:
#line 852 "DynareBison.yy"
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

  case 315:
#line 864 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 316:
#line 866 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 317:
#line 868 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 318:
#line 870 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 319:
#line 872 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 320:
#line 875 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 324:
#line 880 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 325:
#line 882 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 326:
#line 886 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 327:
#line 888 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 328:
#line 890 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 329:
#line 892 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 371:
#line 941 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 372:
#line 943 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 388:
#line 966 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 389:
#line 968 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 390:
#line 972 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 391:
#line 974 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 394:
#line 981 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 395:
#line 983 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 398:
#line 989 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 399:
#line 991 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 400:
#line 993 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 401:
#line 996 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 402:
#line 998 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 403:
#line 1000 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 404:
#line 1002 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 405:
#line 1005 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 406:
#line 1008 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 407:
#line 1010 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 408:
#line 1012 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 409:
#line 1014 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 410:
#line 1017 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 413:
#line 1024 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 414:
#line 1026 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 415:
#line 1028 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 416:
#line 1031 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 417:
#line 1033 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 418:
#line 1035 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 419:
#line 1039 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 420:
#line 1041 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 421:
#line 1045 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 422:
#line 1047 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 423:
#line 1049 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 424:
#line 1051 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 425:
#line 1053 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 426:
#line 1055 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 427:
#line 1059 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 428:
#line 1061 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 429:
#line 1063 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 430:
#line 1065 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 431:
#line 1067 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 432:
#line 1069 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 433:
#line 1073 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 439:
#line 1085 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 440:
#line 1087 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 441:
#line 1089 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 442:
#line 1091 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 443:
#line 1095 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 444:
#line 1097 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 446:
#line 1102 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 447:
#line 1104 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 448:
#line 1106 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 449:
#line 1108 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 450:
#line 1111 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 451:
#line 1112 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 453:
#line 1115 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 454:
#line 1117 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 455:
#line 1119 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 456:
#line 1121 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 480:
#line 1158 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 481:
#line 1160 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 488:
#line 1174 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 489:
#line 1176 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 490:
#line 1180 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 491:
#line 1182 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 521:
#line 1219 "DynareBison.yy"
    { driver.end_homotopy();;}
    break;

  case 524:
#line 1226 "DynareBison.yy"
    { driver.homotopy_val((yysemantic_stack_[(6) - (1)].string_val),(yysemantic_stack_[(6) - (3)].node_val),(yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 527:
#line 1232 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 528:
#line 1233 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1234 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1235 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 531:
#line 1236 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1237 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1238 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 534:
#line 1239 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 535:
#line 1240 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 536:
#line 1241 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 537:
#line 1242 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 538:
#line 1243 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1244 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 540:
#line 1245 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1247 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 542:
#line 1248 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 543:
#line 1249 "DynareBison.yy"
    { driver.option_num("simulation_method",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 544:
#line 1250 "DynareBison.yy"
    { driver.option_num("simulation_method", "0"); ;}
    break;

  case 545:
#line 1251 "DynareBison.yy"
    { driver.option_num("simulation_method", "1"); ;}
    break;

  case 546:
#line 1252 "DynareBison.yy"
    { driver.option_num("simulation_method", "2"); ;}
    break;

  case 547:
#line 1253 "DynareBison.yy"
    { driver.option_num("simulation_method", "3"); ;}
    break;

  case 548:
#line 1254 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 549:
#line 1255 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 550:
#line 1256 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 551:
#line 1257 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 552:
#line 1258 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 553:
#line 1260 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 554:
#line 1262 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 555:
#line 1264 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 556:
#line 1265 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 557:
#line 1266 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 558:
#line 1267 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 559:
#line 1268 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 560:
#line 1270 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 561:
#line 1272 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 562:
#line 1274 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 563:
#line 1275 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 564:
#line 1276 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1277 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1279 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 568:
#line 1280 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 569:
#line 1281 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 570:
#line 1282 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 571:
#line 1283 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1284 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1285 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1286 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 575:
#line 1287 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 576:
#line 1288 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 577:
#line 1289 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 578:
#line 1290 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 579:
#line 1291 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1292 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 581:
#line 1293 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 582:
#line 1294 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 583:
#line 1295 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 584:
#line 1296 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1297 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1299 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 587:
#line 1301 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 588:
#line 1303 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 589:
#line 1304 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 590:
#line 1305 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1306 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 592:
#line 1307 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 593:
#line 1308 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 594:
#line 1309 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 595:
#line 1310 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 596:
#line 1311 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 597:
#line 1313 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 598:
#line 1314 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 599:
#line 1315 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 600:
#line 1316 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 601:
#line 1317 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 602:
#line 1318 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 603:
#line 1319 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 604:
#line 1320 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 605:
#line 1322 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 606:
#line 1323 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 607:
#line 1324 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 608:
#line 1325 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 609:
#line 1326 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 610:
#line 1327 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 611:
#line 1328 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 612:
#line 1329 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 613:
#line 1330 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 614:
#line 1331 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 615:
#line 1332 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 616:
#line 1333 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 617:
#line 1334 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 618:
#line 1335 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 619:
#line 1336 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 620:
#line 1337 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 621:
#line 1338 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 622:
#line 1339 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 623:
#line 1340 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 624:
#line 1342 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 625:
#line 1343 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 626:
#line 1349 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 627:
#line 1350 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 628:
#line 1354 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 629:
#line 1355 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 630:
#line 1356 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 631:
#line 1357 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 632:
#line 1359 "DynareBison.yy"
    {driver.option_num("homotopy_mode",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 633:
#line 1360 "DynareBison.yy"
    {driver.option_num("homotopy_steps",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 634:
#line 1363 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 636:
#line 1372 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
               ;}
    break;

  case 637:
#line 1381 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 638:
#line 1383 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 639:
#line 1391 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2499 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1155;
  const short int
  parser::yypact_[] =
  {
      1186,    20,    25,   -99,  -103,   105,   491,    60,    21,    47,
     -81,    59,   -72,   -12,    10,    82,   461,   531,   463,   -54,
      94,   126,   229,   233,    61,   168,   373,    87, -1155,   251,
     275,   168,   263,   489,   474,   485,    66,    85,   168,   439,
     466,   473,   168,   357,   495,  1063, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155,   365,   667,   405,   935, -1155,   559,    86, -1155,
     503,   593,   437,    44,   -55,   562,   256,   564,   569,   625,
   -1155,  1600,    -4,    95,    96,   394,   574,   569,   628,   626,
     465, -1155,   438,   380,    71,  1345,   599,   616, -1155,  1761,
      13,    98,   576,   218,   657,   502,  1448,  1662,  1662,   255,
      71,   498, -1155,   411, -1155,   -22, -1155,  1761,   257, -1155,
    1706,   258,   261,   582,   262,   583,   268,   584,   269,   281,
     630, -1155,  2294, -1155, -1155, -1155,   685, -1155,   686,   687,
     691,   693,   694, -1155,   696,   697,   699, -1155,   700,   702,
     705,   706, -1155,   572,   530, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155,   710,   711,   713, -1155,   594,   537, -1155, -1155,
   -1155,   538,   668,   -62,   401, -1155,   717,   -50, -1155, -1155,
     548, -1155,   549, -1155, -1155,   674,   340, -1155,   676,   396,
     720,    90, -1155,   678, -1155,   732, -1155, -1155,   733,   734,
     736,   737,   738, -1155, -1155,   739,   742,   743,   744,   746,
     753, -1155, -1155,   755,   768, -1155, -1155, -1155,   769,   770,
   -1155, -1155,   -48, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155,   771,   719, -1155,   723, -1155,   725,   524,
   -1155,   653,   740,   669,   741,   532, -1155,   745,   672,   747,
     534, -1155,   635,   113, -1155,   241,   815,   638,   654, -1155,
    1184, -1155,   -17,   235,   651,   655,   831, -1155, -1155,   238,
   -1155, -1155, -1155, -1155,   783,   785,   138, -1155, -1155, -1155,
     662,   663,   666,   670,  1345,  1345,   673,   679,   680,   681,
     683,   684,   689,   690,   695,   698,   701,  1345,  1259,   703,
     397, -1155,  1288,   398,   841,   844,   845,   846,   851,   853,
   -1155, -1155, -1155,   855,   862,   863, -1155,   864, -1155,   868,
     870,   707, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,   776,   824,
   -1155,   712, -1155, -1155, -1155,   704,   709,   716,   718,  1448,
    1448,   729,   731,   748,   761,   763,   764,   767,   772,   773,
     774,   775,  1448,   813, -1155,   270, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
     278, -1155,   142,    42,   879,   280, -1155, -1155, -1155, -1155,
     883,   889,   287, -1155, -1155, -1155, -1155,   288, -1155, -1155,
     890, -1155,   290, -1155, -1155, -1155, -1155, -1155,   821,   871,
   -1155, -1155,   823,   893, -1155, -1155,   849,   900, -1155, -1155,
     834,   431, -1155,   955,   956,   957,   958,   962,   963,   964,
     965,   977,   978,   979,   980,   981,   985,   987,   992,   993,
    1004,  1005,  1006,  1007,  1011,  1012,  1013,  1014,  1015,  1017,
     842,   899, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,   147,
     264,   147,  1001,   264,  1003,   969,  1008,    28,  1016,  1019,
     972,   973,   667,  1020,  1021,   147,  1025,   935,  1030,   848,
     850,   982,   415,  1031, -1155, -1155,  1032,   503,   857, -1155,
   -1155,   858,    18,   998,   877,    53,  1010,  1345, -1155, -1155,
   -1155,   874,  1036,  1037,  1041,  1048,  1049,   147,   147,   147,
    1050,  1051,  1052,  1053,  1022,   897,   147,  1600,    56,  1023,
    1077,   960, -1155, -1155, -1155,    72,   961,   364,   983, -1155,
   -1155,   984,   364,   986, -1155, -1155,    41, -1155, -1155, -1155,
    1038,   902, -1155,  1045,    48, -1155,   246, -1155,   122, -1155,
     111, -1155,   907,   909,    54,   380,    69,   988,    38, -1155,
   -1155,  1345,  1345,  1345,  1345,   974,   470,  1345,  1345,  1345,
    1345,  1345,  1345,  1345,  1345,  1345,  1345,  1345,   482,  1345,
    1345,  1345,  1345,  1345,  1345,  1345,  1345,  1345,  1345,  1345,
   -1155,  1345, -1155, -1155,  1046,  2056, -1155,  1316,  1081,   147,
    1085,  1087,  1089,  1093,  1094,  1102,   147,  1110,  1116,  1117,
      88, -1155,  1040, -1155,  1448,  1448,    41,  1448,  1018,   505,
    1448,  1448,  1448,  1448,  1448,  1448,  1448,  1448,  1448,  1448,
    1448,   664,  1448,  1448,  1448,  1448,  1448,  1448,  1448,  1448,
    1448,  1448,  1448,   966,  1662,    89,    92, -1155, -1155, -1155,
    1345,   366,    30,   429,   411,   967,  1119,  1130,   -22,   975,
    1761,    93,   147,  1706,    99, -1155,  1054, -1155,  1055, -1155,
    1056,  1345, -1155, -1155,  1135,  1138,  1142,  1144,  1145,  1154,
    1155,  1156,  1158,  1159,  1164,  1166,  1171,  1172,  1174,   147,
     147,  1175,   874,   147,   147,  1176,  1177,   147,  1178,   147,
     147,  1024,  2294, -1155, -1155, -1155, -1155,  1190,  1192, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155,  1183,    17, -1155,
   -1155, -1155, -1155,  1026, -1155, -1155,  1033, -1155, -1155, -1155,
   -1155,  1034, -1155,  1188,  1035,  1039,  1042,  1345, -1155, -1155,
   -1155, -1155, -1155,   284,  1043, -1155, -1155,   286,  1044,  2068,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155,  1028, -1155, -1155, -1155,   298, -1155,
    1160,  1161, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
     416,  1047,  1095,  1096,  1173,  1104,   364,  1179,  1060,   364,
   -1155,  1205,  1208,  1067, -1155,   569,  1229, -1155, -1155, -1155,
    1448, -1155,  1230,   291, -1155, -1155, -1155,  1083, -1155, -1155,
   -1155,   292, -1155, -1155, -1155,  1084, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155,   -44,    67, -1155,  1200,
    1345,  1209,     4,   641,  2298,  2280,   294,  2357,   781,   971,
    1002,  1125,  1381,  1395,  1411,  1423,  1683,  1736,  1750, -1155,
     488,   488,   488,   488,   488,   488,   470,   470,   974,   974,
    1105,  1762,  1345, -1155,  1212,  2082, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,   299,
   -1155,  2369,  2381,  1088,  2393,  1774,  1786,  1800,  1820,  1832,
    1844,  1858,  1870,  1890,  1902,  1916, -1155,   544,   544,   544,
     544,   544,   544,   505,   505,  1018,  1018,  1105, -1155, -1155,
   -1155,   300, -1155,   301,  1928,    42,  1091, -1155, -1155,    45,
    1345, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155,   303, -1155, -1155, -1155,   307, -1155,
   -1155, -1155,  2405, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155,  1092, -1155, -1155, -1155,  1220, -1155,
   -1155,  1097,  1261, -1155, -1155,  2094, -1155,   100, -1155,   203,
   -1155,  1221, -1155,   295, -1155, -1155, -1155, -1155, -1155, -1155,
     364,    72,  1157,   364,  1163,  1165, -1155,  1100, -1155, -1155,
    1278,   490,  1448,  2112,   147,   122, -1155,  1184,   111, -1155,
    1184,  1184,  1184,    69, -1155,   364, -1155,  1281,  2124,  1284,
    1267,  1345,  1345,  1345,  1345, -1155,  1345, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,  1112,  2138,
    1345, -1155, -1155,  1448,  1448, -1155,  1448, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155,    30, -1155, -1155, -1155,  1345,  1940, -1155, -1155,  1345,
    1274, -1155,  1035,  1345, -1155, -1155,   339, -1155,   344,  1115,
    1028, -1155, -1155,  1180,  1181,  1182,   364,  1123,   364,   364,
   -1155,  1345, -1155,  2150, -1155, -1155, -1155,  1128,   197, -1155,
    1131,   214,   675,   692,   132,  1132,  1345, -1155,  1345,  1129,
      39,  2168,  1960,  1974,  2280,  2417, -1155, -1155,  2180,  1986,
    1998,  2429,  2010, -1155,  2194, -1155,  1295,  2206, -1155, -1155,
    1201, -1155,   364,   364,   364,  1203, -1155,  1148,  1150,  2224,
   -1155,  1184, -1155,  1184, -1155, -1155, -1155,   364, -1155,  2236,
    2250,  1308,  1149,  1312,  1231, -1155, -1155, -1155,  1345, -1155,
   -1155, -1155,  1448, -1155, -1155,  1345, -1155,    33,  1214, -1155,
    1215,   364, -1155, -1155, -1155,   839,   867,  1162, -1155, -1155,
   -1155,  1317,  1167,  1345,  2024,  2044,  2262,  1285, -1155,   364,
      97,  1168, -1155, -1155, -1155, -1155,  1319,  2280,    27, -1155,
   -1155, -1155,  1169,  1222,  1223, -1155, -1155,  1345, -1155, -1155,
   -1155,   364,   364,  2280,  1225, -1155,   364, -1155
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   450,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     2,     4,    29,    30,
      46,    47,    48,    45,     5,     6,     7,    12,     9,    10,
      11,     8,    13,    14,    15,    16,    17,    18,    19,    23,
      25,    24,    20,    21,    22,    26,    27,    28,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,     0,     0,     0,     0,   419,     0,     0,   234,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   279,
     326,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   143,     0,     0,     0,     0,     0,     0,   406,     0,
       0,     0,    76,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   239,     0,   227,     0,   246,     0,     0,   453,
       0,     0,     0,    58,     0,    64,     0,    70,     0,     0,
       0,   490,     0,     1,     3,   480,     0,   602,     0,     0,
       0,     0,     0,   593,     0,     0,     0,   594,     0,     0,
       0,     0,   468,   479,     0,   469,   474,   472,   475,   473,
     470,   471,   476,   477,   461,   462,   463,   464,   465,   466,
     467,   488,     0,     0,     0,   482,   487,     0,   484,   483,
     485,     0,     0,   416,     0,   412,     0,     0,   237,   238,
       0,    82,     0,    49,   429,     0,     0,   423,     0,     0,
       0,     0,   126,     0,   577,     0,   582,   561,     0,     0,
       0,     0,     0,   574,   575,     0,     0,     0,     0,     0,
       0,   595,   570,     0,     0,   581,   576,   560,     0,     0,
     580,   578,     0,   331,   367,   356,   332,   333,   334,   335,
     336,   337,   338,   339,   340,   341,   342,   343,   344,   345,
     346,   347,   348,   349,   350,   351,   352,   353,   354,   355,
     357,   358,   359,   360,   361,   362,   363,   364,   365,   366,
     368,   369,   370,   275,     0,   328,     0,   292,     0,     0,
     289,     0,     0,     0,     0,     0,   311,     0,     0,     0,
       0,   305,     0,     0,   130,     0,     0,     0,     0,    84,
       0,   530,     0,     0,     0,     0,     0,   589,   588,     0,
     435,   436,   437,   438,     0,     0,     0,   200,    89,    90,
       0,     0,    88,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   397,     0,     0,     0,     0,     0,     0,     0,     0,
     535,   536,   537,     0,     0,     0,   583,     0,   549,     0,
       0,     0,   252,   253,   254,   255,   256,   257,   258,   259,
     260,   261,   262,   264,   266,   267,   268,   269,   270,   271,
     272,   263,   265,   273,   274,   408,   405,    79,    74,     0,
      55,     0,    80,   165,   166,     0,     0,   195,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   451,   164,     0,   374,   379,   375,   376,
     377,   378,   380,   381,   382,   383,   384,   385,   386,   387,
       0,    51,     0,     0,     0,     0,   242,   243,   245,   244,
       0,     0,     0,   230,   231,   232,   233,     0,   251,   248,
       0,   459,     0,   458,   460,   455,   399,    61,    56,     0,
      52,    67,    62,     0,    53,    73,    68,     0,    54,   394,
       0,     0,   522,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   493,   494,   495,   496,   497,   498,   499,   500,   501,
     502,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     520,   512,   513,   514,   515,   516,   517,   518,   519,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   410,   411,     0,     0,     0,    83,
      50,     0,     0,     0,     0,     0,     0,     0,   124,   125,
     280,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     277,     0,   291,   287,   288,   320,     0,   320,     0,   309,
     310,     0,   320,     0,   303,   304,     0,   128,   129,   121,
       0,     0,    85,     0,     0,   159,     0,   160,     0,   155,
       0,   151,     0,     0,     0,     0,     0,     0,     0,   198,
     199,     0,     0,     0,     0,   102,   103,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      86,     0,   395,   396,     0,     0,   400,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    77,    75,    81,     0,     0,     0,     0,   178,   179,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   197,   225,   226,
       0,     0,   217,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    59,    57,    65,    63,    71,
      69,     0,   521,   523,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   526,   525,   598,   286,     0,     0,   599,
     600,   601,   597,   603,   552,   555,   554,     0,     0,   553,
     556,   557,   590,     0,   591,   478,     0,   604,   562,   579,
     486,     0,   420,     0,   416,     0,     0,     0,   528,   236,
     235,   432,   427,     0,     0,   426,   421,     0,     0,     0,
     592,   538,   584,   585,   558,   559,   564,   567,   565,   572,
     573,   563,   569,   568,     0,   571,   330,   327,     0,   276,
       0,     0,   315,   322,   316,   321,   318,   323,   317,   319,
       0,     0,     0,   297,     0,     0,   320,     0,     0,   320,
     283,     0,     0,     0,   123,     0,     0,   144,   157,   158,
       0,   162,     0,     0,   135,   141,   142,     0,   139,   138,
     140,     0,   133,   136,   137,     0,   147,   145,   586,   587,
     434,   445,   447,   448,   449,   446,     0,   439,   443,     0,
       0,     0,     0,     0,     0,   119,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    87,
     101,   100,    99,    98,    97,    96,    92,    91,    93,    94,
      95,     0,     0,   403,     0,     0,   534,   542,   527,   533,
     539,   540,   531,   541,   551,   532,   529,   550,   407,     0,
      78,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   163,   176,   175,   174,
     173,   172,   171,   168,   167,   169,   170,   177,   452,   373,
     371,     0,   388,     0,     0,     0,     0,   222,   223,     0,
       0,   547,   545,   546,   543,   544,   241,   240,   632,   633,
     229,   228,   250,   247,     0,   596,   457,   454,     0,    60,
      66,    72,     0,   605,   606,   607,   608,   609,   610,   611,
     612,   613,   614,   615,   616,   617,   618,   619,   620,   621,
     622,   623,   624,   625,   626,   627,   628,   629,   630,   631,
     491,   492,   285,   284,   635,   637,   639,   638,     0,   481,
     489,     0,     0,   418,   417,     0,   428,     0,   422,     0,
     127,     0,   392,     0,   329,   278,   293,   325,   324,   290,
     320,   320,     0,   320,     0,     0,   308,     0,   282,   281,
       0,     0,     0,     0,     0,     0,   153,     0,     0,   149,
       0,     0,     0,     0,   433,   320,   444,     0,     0,     0,
       0,     0,     0,     0,     0,   117,     0,   104,   105,   106,
     107,   108,   109,   110,   111,   112,   113,   114,     0,     0,
       0,   401,   409,     0,     0,   196,     0,   180,   181,   182,
     183,   184,   185,   186,   187,   188,   189,   190,   372,   389,
     224,   216,   212,   219,   220,     0,     0,   249,   456,     0,
       0,   634,   416,     0,   413,   430,     0,   424,     0,     0,
       0,   566,   294,     0,     0,     0,   320,     0,   320,   320,
     306,     0,   122,     0,   161,   548,   134,     0,     0,   132,
       0,     0,     0,     0,   440,     0,     0,   203,     0,   211,
       0,     0,     0,     0,   120,     0,   398,   404,     0,     0,
       0,     0,     0,   221,     0,   636,     0,     0,   431,   425,
       0,   393,   320,   320,   320,     0,   314,     0,     0,     0,
     194,     0,   156,     0,   152,   148,   146,   320,   441,     0,
       0,     0,   207,     0,     0,   202,   115,   116,     0,   402,
     191,   192,     0,   218,   524,     0,   414,   320,   299,   295,
     298,   320,   312,   307,   131,     0,     0,     0,   205,   204,
     210,     0,   206,     0,     0,     0,     0,     0,   391,   320,
       0,     0,   154,   150,   442,   209,     0,   215,     0,   118,
     193,   415,     0,   300,     0,   313,   208,     0,   201,   214,
     390,   320,   320,   213,   301,   296,   320,   302
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1155, -1155,  1335, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,  -350, -1155,
   -1155, -1155, -1155,  -110,  -228, -1155, -1155,  1059, -1155, -1155,
     265, -1155,  -618, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155,  -989,  -595,  -136,  -588, -1155, -1155, -1155,  1245,  -260,
   -1155, -1155, -1155, -1155,   372, -1155, -1155,   631, -1155, -1155,
     803, -1155, -1155,   637, -1155, -1155,  -119,   -23,   682,   829,
   -1155, -1155,  1098, -1155, -1155, -1154, -1155, -1155,  1074, -1155,
   -1155,  1086, -1033,  -589, -1155, -1155,   779, -1155,  1262,   658,
   -1155,   215, -1155, -1155, -1155, -1155,  1057, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155,  1204,  -768, -1155, -1155, -1155, -1155,
   -1155,   749, -1155,   297,  -860, -1155, -1155, -1155, -1155, -1155,
     660, -1155,   -53,   852, -1155, -1155,   856, -1155, -1155,   619,
   -1155, -1155, -1155,   924,  -555, -1155,   -82, -1155,  1304, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155,   -90, -1155, -1155,  -132,
    -573, -1155, -1155, -1155, -1155, -1155,  -102,   -88,   -86,   -85,
     -84, -1155, -1155,   -87,   -56, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155,   -79, -1155, -1155, -1155, -1155, -1155,   -77,
     -76,   -49,   -74,   -71,   -67, -1155, -1155, -1155, -1155,  -111,
    -107,   -80,   -78,   -66,   -47,   -46, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155, -1155,
   -1155, -1155, -1155, -1155, -1155, -1155,   610, -1155,  -532
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    45,    46,    47,    48,    49,    50,    51,    52,    53,
     154,   156,   158,   133,    54,    55,    56,    57,   368,   936,
      58,   327,    59,   231,   232,    60,   323,   324,   911,   903,
     912,   913,   914,    61,   330,  1122,  1121,  1210,   915,  1207,
     907,   644,   645,   646,   647,   444,    62,    63,   346,   347,
    1220,    64,  1308,   751,   752,    65,   472,   473,    66,   217,
     218,    67,   465,   466,    68,   477,   481,   112,   893,   809,
      69,   309,   310,   311,   881,  1192,    70,   320,   321,    71,
     315,   316,   882,  1193,    72,   262,   263,    73,   445,   446,
      74,  1092,  1093,    75,    76,   370,   371,    77,    78,   373,
      79,    80,    81,   214,   215,   583,    82,    83,    84,    85,
     339,   340,   926,   927,   928,    86,   136,   743,    87,   482,
     483,   182,   183,   184,    88,   206,   207,    89,    90,   530,
     531,    91,   501,   502,   805,   392,   393,   394,   395,   396,
     397,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   468,   906,   408,   409,   410,   185,   186,   187,   188,
     189,   271,   272,   411,   449,   275,   276,   277,   278,   279,
     280,   281,   282,   450,   284,   285,   286,   287,   288,   451,
     452,   453,   454,   455,   456,   412,   295,   296,   341,   413,
     414,   190,   191,   459,   192,   193,   302,   484,   194,   195,
     196,   197,   198,   199,   200,   210,   532,   533,   534,   535,
     536,   537,   538,   539,   540,   541,   542,   543,   544,   545,
     546,   547,   548,   549,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   475,   476,   824,  1075,   818,   819
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       443,   130,   131,   599,   665,   666,   810,   325,   139,   266,
     391,   467,   342,   148,   151,   152,   343,   678,   219,   159,
     828,   265,   695,   267,   273,   268,   269,   270,   478,   264,
     904,   297,   283,   298,   289,   290,   883,   292,   885,   208,
     293,   469,   205,   888,   294,   299,   209,   447,   447,   898,
     448,   448,   856,   857,   858,   274,   899,   457,   457,   458,
     458,   865,   291,   474,   300,   301,  1082,  1126,  1194,   850,
     803,  1074,   930,  1017,    92,   905,   873,   905,   872,    94,
     897,   804,   816,   303,  1018,   748,   660,   875,  1173,  1279,
     212,   423,   581,  1262,   344,   890,   749,   599,   222,  1174,
     303,   216,   424,   872,   587,   109,   617,    98,   104,  1130,
    1123,    96,    97,   918,   103,   873,   425,   874,   306,   312,
     877,   134,   598,   470,   471,   426,   875,   876,  1208,   108,
    1131,  1211,  1212,  1213,   106,   427,   375,   648,   113,   919,
     109,   135,   874,   109,   967,   637,   109,   375,   109,   582,
     304,   974,   876,   109,   921,   224,   921,   123,   908,   877,
     125,   344,   588,   225,   618,   344,  1124,   304,  1325,   878,
     659,   909,   109,   213,   747,   109,   109,   230,   902,   109,
     109,  1317,   307,   313,  1019,   303,   109,   109,   880,   902,
     803,   428,   931,  1263,   878,   649,   891,   892,   114,   910,
     322,   804,   660,   429,   430,   879,   305,  1035,   345,   431,
     432,   433,   434,   435,   436,   437,   438,   439,   440,   921,
     115,   308,   314,   415,   922,   441,   922,   880,   842,  1252,
     879,    93,   105,  1076,  1058,  1059,    95,  1318,  1062,  1063,
     423,  1020,  1066,   817,  1068,  1069,  1254,   849,   932,  1264,
    1297,   424,   304,   750,   223,   109,  1175,   423,   107,   442,
    1061,   643,  1285,   846,  1286,   425,   867,  1314,   424,   110,
     111,   128,   129,   639,   426,   345,   146,   147,  1125,   345,
     900,   923,   425,   923,   427,   924,   925,   924,   925,   922,
     109,   426,   116,   718,   719,   149,   150,  1104,   978,  1010,
    1107,   427,  1012,  1033,   124,   418,   731,   806,   416,  1037,
    1185,   933,   934,   935,   937,    99,   100,   938,   939,   940,
     941,   942,   943,   944,   945,   946,   947,   948,   230,   950,
     951,   952,   953,   954,   955,   956,   957,   958,   959,   960,
     428,   961,   303,  1257,   303,   303,   923,   965,   303,   488,
     924,   925,   429,   430,  1126,   492,   496,   428,   431,   432,
     433,   434,   435,   436,   437,   438,   439,   440,   303,   429,
     430,   303,   419,   303,   441,   431,   432,   433,   434,   435,
     436,   437,   438,   439,   440,   303,   303,   303,   303,   650,
     303,   441,   655,   732,   303,   733,   734,   735,   736,   737,
    1014,   738,   739,   740,   741,   212,   742,   873,   442,   304,
     643,   304,   304,  1187,  1236,   304,   489,   317,   875,   807,
     808,  1042,   493,   497,   744,   442,   303,   643,   420,   692,
     696,   303,   744,   584,   754,   304,  1021,   174,   304,   126,
     304,   758,   760,   127,   763,  1115,  1118,   651,  1134,  1190,
     656,   877,   304,   304,   304,   304,   901,   304,   835,  1097,
     132,   304,   137,   772,   336,   461,   227,   479,   485,   836,
    1098,   486,   490,   140,   228,  1022,   337,  1023,   494,   498,
     464,   318,   745,  1024,   369,   697,   138,  1085,   213,   338,
     746,   499,   755,   304,  1086,  1025,  1088,  1206,   304,   759,
     761,   331,   764,  1116,  1119,   219,  1135,  1191,  1094,  1152,
    1168,  1169,  1195,  1177,  1197,   266,   384,  1178,   500,   880,
     319,   208,  1202,   141,   205,   101,   153,   265,   209,   267,
     273,   268,   269,   270,   102,   264,  1215,   297,   283,   298,
     289,   290,   905,   292,   342,   905,   293,   306,   343,  1238,
     294,   299,   592,   155,  1239,   312,   623,   317,   593,  1205,
     157,   274,   332,   333,   629,   119,   634,   160,   291,   843,
     300,   301,   847,   334,   120,   165,  1015,   230,   981,   982,
    1128,   984,  1016,   211,   985,   986,   987,   988,   989,   990,
     991,   992,   993,   994,   995,   868,   997,   998,   999,  1000,
    1001,  1002,  1003,  1004,  1005,  1006,  1007,  1245,   595,  1247,
    1248,   307,  1149,   898,   596,   201,   898,   898,   898,   313,
     899,   318,   467,   899,   899,   899,   216,   687,   688,   679,
     689,   680,   681,   682,   683,   684,   220,   685,   686,   687,
     688,  1032,   689,   685,   686,   687,   688,   221,   689,   226,
     308,   229,   469,  1278,   447,  1280,   230,   448,   314,   233,
     319,   322,   740,   741,   457,   742,   458,   326,  1287,   328,
    1176,   117,   118,   121,   122,   329,   474,   166,   167,   168,
     169,   170,   171,   172,   142,   143,   369,   979,  1298,   173,
     898,   898,  1301,   174,   949,   144,   145,   899,   899,   738,
     739,   740,   741,   372,   742,   161,   162,  1255,   417,   175,
    1313,   421,   422,   463,   487,   491,   495,   500,   423,   559,
     560,   561,  1011,  1013,  1256,   562,   572,   563,   564,   424,
     565,   566,  1324,   567,   568,   423,   569,  1327,  1034,   570,
     571,  1038,   573,   425,   574,   575,   424,   576,   577,   578,
     579,   586,   426,  1319,   597,   580,   176,   177,   589,   590,
     425,   591,   427,   594,  1113,   600,   601,   602,   603,   426,
     604,   605,   606,   607,   178,   179,   608,   609,   610,   427,
     611,  1221,  1222,  1223,  1224,  1111,  1225,   612,   679,   613,
     680,   681,   682,   683,   684,  1132,   685,   686,   687,   688,
    1228,   689,   614,   615,   616,   619,   620,   625,   180,   181,
     621,   732,   622,   733,   734,   735,   736,   737,   428,   738,
     739,   740,   741,   627,   742,  1232,   632,   626,   628,  1234,
     429,   430,   631,  1237,   633,   428,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   636,   429,   430,   640,
     641,  1249,   441,   431,   432,   433,   434,   435,   436,   437,
     438,   439,   440,   652,   642,   654,  1259,   653,  1260,   441,
     657,  1302,   658,   661,   662,   698,   996,   663,   699,   700,
     701,   664,   423,   599,   667,   702,   442,   703,   643,   704,
     668,   669,   670,   424,   671,   672,   705,   706,   707,  1303,
     673,   674,   708,   442,   709,   643,   675,   425,   711,   676,
     423,   712,   677,   753,   691,   714,   426,   756,  1294,   710,
     715,   424,   713,   757,   762,  1296,   427,   716,   679,   717,
     680,   681,   682,   683,   684,   425,   685,   686,   687,   688,
     720,   689,   721,  1307,   426,   166,   167,   168,   169,   170,
     171,   172,   202,   765,   427,   767,   203,   173,   766,   722,
     732,   174,   733,   734,   735,   736,   737,  1323,   738,   739,
     740,   741,   723,   742,   724,   725,  1203,   175,   726,   204,
     768,   769,   428,   727,   728,   729,   730,   770,   771,   774,
     775,   776,   777,  1137,   429,   430,   778,   779,   780,   781,
     431,   432,   433,   434,   435,   436,   437,   438,   439,   440,
     428,   782,   783,   784,   785,   786,   441,  1229,  1230,   787,
    1231,   788,   429,   430,   176,   177,   789,   790,   431,   432,
     433,   434,   435,   436,   437,   438,   439,   440,   791,   792,
     793,   794,   178,   179,   441,   795,   796,   797,   798,   799,
     442,   800,   643,   802,   801,   811,   814,   813,   832,   822,
     823,   833,   815,   163,  1186,   837,  1188,   840,   841,   834,
     820,     1,     2,   821,   826,   827,   180,   181,   442,   829,
     643,     3,     4,     5,   831,   844,   838,   845,     6,   817,
     851,   852,     7,     8,     9,   853,    10,   848,    11,    12,
      13,    14,   854,   855,   859,   860,   861,   862,   864,   863,
     869,   870,   895,    15,   871,   884,    16,   916,   679,   917,
     680,   681,   682,   683,   684,   894,   685,   686,   687,   688,
      17,   689,   896,   962,   689,   966,  1295,   886,   887,   968,
     889,   969,   929,   970,    18,    19,    20,   971,   972,   679,
      21,   680,   681,   682,   683,   684,   973,   685,   686,   687,
     688,    22,   689,    23,   975,    24,    25,    26,    27,    28,
     976,   977,   980,  1028,    29,    30,  1008,  1027,   742,    31,
      32,    33,    34,  1138,  1029,  1031,  1039,  1040,  1041,  1043,
      35,    36,  1044,    37,     1,     2,  1045,    38,  1046,  1047,
      39,    40,    41,    42,     3,     4,     5,    43,  1048,  1049,
    1050,     6,  1051,  1052,  1139,     7,     8,     9,  1053,    10,
    1054,    11,    12,    13,    14,  1055,  1056,   423,  1057,  1060,
    1064,  1065,  1067,  1072,  1070,  1073,    15,  1074,   424,    16,
    1078,    44,  1081,  1079,  1080,  1091,   582,  1095,  1096,  1100,
    1101,  1083,   425,    17,  1084,  1087,  1089,  1099,  1103,  1108,
    1102,   426,  1109,  1112,  1114,    -1,  1105,    18,    19,    20,
    1106,   427,   679,    21,   680,   681,   682,   683,   684,  1110,
     685,   686,   687,   688,    22,   689,    23,  1127,    24,    25,
      26,    27,    28,  1117,  1120,  1183,  1129,    29,    30,  1150,
    1155,  1172,    31,    32,    33,    34,  1180,  1181,  1189,  1182,
    1200,  1196,  1201,    35,    36,  1216,    37,  1198,  1218,  1199,
      38,  1219,  1226,    39,    40,    41,    42,   428,  1235,  1275,
      43,   348,  1240,  1246,  1242,  1243,  1244,  1140,  1251,   429,
     430,  1253,   349,  1261,  1258,   431,   432,   433,   434,   435,
     436,   437,   438,   439,   440,  1277,   350,  1281,  1282,   348,
    1283,   441,  1290,  1291,    44,   351,  1292,  1293,  1299,  1300,
     349,  1305,  1312,  1316,  1304,   352,  1321,  1322,  1315,  1326,
     164,  1306,   638,  1209,   350,   462,  1320,  1171,   348,  1030,
     839,  1026,   812,   351,   635,   442,   866,   643,   983,   349,
     460,   630,  1009,   352,   920,  1241,   679,   624,   680,   681,
     682,   683,   684,   350,   685,   686,   687,   688,   585,   689,
    1214,  1071,   351,  1036,   825,   773,   335,   693,  1077,     0,
       0,   353,   352,   830,     0,     0,     0,     0,     0,     0,
       0,     0,   694,   354,   355,     0,     0,     0,     0,   356,
     357,   358,   359,   360,   361,   362,   363,   364,   365,   353,
       0,     0,     0,     0,     0,   366,     0,     0,     0,   690,
     964,   354,   355,     0,     0,     0,     0,   356,   357,   358,
     359,   360,   361,   362,   363,   364,   365,     0,   353,     0,
       0,   423,     0,   366,     0,     0,     0,     0,     0,   367,
     354,   355,   424,     0,     0,     0,   356,   357,   358,   359,
     360,   361,   362,   363,   364,   365,   425,     0,     0,     0,
       0,     0,   366,     0,     0,   426,     0,   367,   679,     0,
     680,   681,   682,   683,   684,   427,   685,   686,   687,   688,
       0,   689,   679,     0,   680,   681,   682,   683,   684,     0,
     685,   686,   687,   688,     0,   689,   367,     0,   679,     0,
     680,   681,   682,   683,   684,     0,   685,   686,   687,   688,
     679,   689,   680,   681,   682,   683,   684,     0,   685,   686,
     687,   688,     0,   689,     0,     0,     0,     0,     0,     0,
       0,   428,     0,  1141,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   429,   430,   234,     0,  1142,     0,   431,
     432,   433,   434,   435,   436,   437,   438,   439,   440,     0,
       0,   203,   173,  1143,     0,   441,   174,     0,     0,     0,
       0,     0,     0,     0,     0,  1144,     0,     0,     0,     0,
     235,   236,   175,     0,   204,     0,     0,     0,     0,   237,
       0,     0,     0,     0,     0,     0,   238,   239,   240,   442,
       0,   241,   242,     0,   243,   244,     0,   234,     0,     0,
     245,   246,   247,   248,   249,   250,   251,     0,   252,   253,
     254,     0,     0,   203,     0,     0,   255,     0,     0,   176,
     177,     0,   256,     0,   257,     0,     0,     0,     0,   258,
       0,     0,   235,   236,     0,     0,   204,   178,   179,   374,
     259,   237,     0,     0,     0,     0,     0,     0,   238,     0,
       0,     0,   260,   216,     0,     0,     0,     0,     0,   261,
       0,   375,     0,   376,   377,     0,     0,     0,     0,     0,
       0,   180,   181,     0,     0,     0,     0,     0,   255,     0,
       0,     0,     0,     0,     0,   237,   257,   378,   379,     0,
       0,     0,   238,     0,   374,     0,     0,     0,     0,   331,
       0,     0,   259,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   260,     0,   375,     0,   376,   377,
       0,   261,     0,     0,     0,     0,     0,   380,     0,   381,
     257,   382,   337,   180,   181,     0,     0,   383,     0,     0,
     237,   384,   378,   379,     0,   338,     0,   238,     0,   385,
     386,   387,     0,     0,   331,   388,   389,   390,     0,   216,
     679,     0,   680,   681,   682,   683,   684,   480,   685,   686,
     687,   688,     0,   689,     0,     0,     0,     0,     0,     0,
       0,     0,   380,     0,   381,   257,   382,   337,     0,     0,
       0,     0,   383,     0,     0,     0,   384,     0,     0,     0,
     338,     0,     0,     0,   385,   386,   387,     0,     0,     0,
     388,   389,   390,   679,   216,   680,   681,   682,   683,   684,
       0,   685,   686,   687,   688,  1145,   689,   679,     0,   680,
     681,   682,   683,   684,     0,   685,   686,   687,   688,   679,
     689,   680,   681,   682,   683,   684,     0,   685,   686,   687,
     688,   732,   689,   733,   734,   735,   736,   737,     0,   738,
     739,   740,   741,   732,   742,   733,   734,   735,   736,   737,
       0,   738,   739,   740,   741,     0,   742,   732,  1146,   733,
     734,   735,   736,   737,     0,   738,   739,   740,   741,     0,
     742,     0,  1147,     0,     0,     0,     0,   732,     0,   733,
     734,   735,   736,   737,  1148,   738,   739,   740,   741,   732,
     742,   733,   734,   735,   736,   737,  1157,   738,   739,   740,
     741,   732,   742,   733,   734,   735,   736,   737,  1158,   738,
     739,   740,   741,     0,   742,   732,     0,   733,   734,   735,
     736,   737,  1159,   738,   739,   740,   741,   732,   742,   733,
     734,   735,   736,   737,     0,   738,   739,   740,   741,     0,
     742,     0,  1160,     0,     0,     0,     0,   732,     0,   733,
     734,   735,   736,   737,  1161,   738,   739,   740,   741,   732,
     742,   733,   734,   735,   736,   737,  1162,   738,   739,   740,
     741,     0,   742,   732,     0,   733,   734,   735,   736,   737,
    1163,   738,   739,   740,   741,   679,   742,   680,   681,   682,
     683,   684,  1164,   685,   686,   687,   688,   679,   689,   680,
     681,   682,   683,   684,     0,   685,   686,   687,   688,     0,
     689,     0,  1165,     0,     0,     0,     0,   679,     0,   680,
     681,   682,   683,   684,  1166,   685,   686,   687,   688,     0,
     689,   679,     0,   680,   681,   682,   683,   684,  1167,   685,
     686,   687,   688,   732,   689,   733,   734,   735,   736,   737,
    1170,   738,   739,   740,   741,   732,   742,   733,   734,   735,
     736,   737,  1233,   738,   739,   740,   741,   679,   742,   680,
     681,   682,   683,   684,     0,   685,   686,   687,   688,     0,
     689,   679,  1266,   680,   681,   682,   683,   684,     0,   685,
     686,   687,   688,     0,   689,     0,  1267,     0,     0,     0,
       0,   732,     0,   733,   734,   735,   736,   737,  1270,   738,
     739,   740,   741,   679,   742,   680,   681,   682,   683,   684,
    1271,   685,   686,   687,   688,   679,   689,   680,   681,   682,
     683,   684,  1273,   685,   686,   687,   688,     0,   689,   679,
       0,   680,   681,   682,   683,   684,  1309,   685,   686,   687,
     688,   679,   689,   680,   681,   682,   683,   684,     0,   685,
     686,   687,   688,     0,   689,     0,  1310,     0,     0,   732,
       0,   733,   734,   735,   736,   737,   963,   738,   739,   740,
     741,   679,   742,   680,   681,   682,   683,   684,  1090,   685,
     686,   687,   688,     0,   689,   679,     0,   680,   681,   682,
     683,   684,  1151,   685,   686,   687,   688,   732,   689,   733,
     734,   735,   736,   737,  1184,   738,   739,   740,   741,     0,
     742,     0,     0,     0,     0,   679,     0,   680,   681,   682,
     683,   684,  1204,   685,   686,   687,   688,   679,   689,   680,
     681,   682,   683,   684,  1217,   685,   686,   687,   688,     0,
     689,   679,     0,   680,   681,   682,   683,   684,  1227,   685,
     686,   687,   688,   679,   689,   680,   681,   682,   683,   684,
    1250,   685,   686,   687,   688,     0,   689,     0,     0,     0,
       0,   679,     0,   680,   681,   682,   683,   684,  1265,   685,
     686,   687,   688,   679,   689,   680,   681,   682,   683,   684,
    1269,   685,   686,   687,   688,     0,   689,   679,     0,   680,
     681,   682,   683,   684,  1274,   685,   686,   687,   688,   679,
     689,   680,   681,   682,   683,   684,  1276,   685,   686,   687,
     688,     0,   689,     0,     0,     0,     0,   679,     0,   680,
     681,   682,   683,   684,  1284,   685,   686,   687,   688,     0,
     689,     0,     0,     0,     0,   679,  1288,   680,   681,   682,
     683,   684,  1133,   685,   686,   687,   688,     0,   689,     0,
    1289,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1311,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,   517,   518,   519,
     520,   521,   522,   523,     0,     0,     0,   524,   525,     0,
     526,   527,   528,   529,   679,     0,   680,   681,   682,   683,
     684,  1136,   685,   686,   687,   688,   732,   689,   733,   734,
     735,   736,   737,  1153,   738,   739,   740,   741,   732,   742,
     733,   734,   735,   736,   737,  1154,   738,   739,   740,   741,
     732,   742,   733,   734,   735,   736,   737,  1156,   738,   739,
     740,   741,   679,   742,   680,   681,   682,   683,   684,  1179,
     685,   686,   687,   688,   679,   689,   680,   681,   682,   683,
     684,  1268,   685,   686,   687,   688,   732,   689,   733,   734,
     735,   736,   737,  1272,   738,   739,   740,   741,     0,   742
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
       136,    24,    25,   231,   354,   355,   561,   117,    31,   111,
     129,   143,   123,    36,    37,    38,   123,   367,   100,    42,
     575,   111,   372,   111,   111,   111,   111,   111,   147,   111,
     648,   111,   111,   111,   111,   111,   625,   111,   627,    95,
     111,   143,    95,   632,   111,   111,    95,   137,   138,   644,
     137,   138,   607,   608,   609,   111,   644,   137,   138,   137,
     138,   616,   111,   145,   111,   111,   834,   927,  1101,   601,
      43,    54,    34,    43,    54,   648,    43,   650,     6,    54,
      32,    54,    54,    87,    54,    43,   346,    54,    43,  1243,
       4,    43,   154,    54,    23,    54,    54,   325,    54,    54,
      87,   123,    54,     6,   154,    87,   154,   210,    87,   105,
     154,   210,   211,    59,    54,    43,    68,    45,    23,    23,
      87,    34,    32,   145,   146,    77,    54,    55,  1117,   210,
     126,  1120,  1121,  1122,    87,    87,    25,   154,   210,    85,
      87,    54,    45,    87,   699,    32,    87,    25,    87,   211,
     154,   706,    55,    87,    87,   210,    87,   211,    47,    87,
      34,    23,   212,   218,   212,    23,   210,   154,  1322,    97,
      32,    60,    87,    87,    32,    87,    87,    87,    67,    87,
      87,   154,    87,    87,   154,    87,    87,    87,   155,    67,
      43,   143,   154,   154,    97,   212,   155,   156,   210,    88,
      87,    54,   462,   155,   156,   133,   210,   762,   137,   161,
     162,   163,   164,   165,   166,   167,   168,   169,   170,    87,
     210,   126,   126,   210,   157,   177,   157,   155,   210,    32,
     133,   211,   211,   216,   789,   790,   211,   210,   793,   794,
      43,   211,   797,   215,   799,   800,    32,   597,   210,   210,
     217,    54,   154,   211,   210,    87,   211,    43,   211,   211,
     792,   213,  1251,   210,  1253,    68,   210,  1300,    54,   210,
     211,   210,   211,    32,    77,   137,   210,   211,   211,   137,
      34,   214,    68,   214,    87,   218,   219,   218,   219,   157,
      87,    77,   210,   429,   430,   210,   211,   886,   210,   210,
     889,    87,   210,   210,   210,    87,   442,    43,   210,   210,
     210,   661,   662,   663,   664,   210,   211,   667,   668,   669,
     670,   671,   672,   673,   674,   675,   676,   677,    87,   679,
     680,   681,   682,   683,   684,   685,   686,   687,   688,   689,
     143,   691,    87,   211,    87,    87,   214,   697,    87,    87,
     218,   219,   155,   156,  1214,    87,    87,   143,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,    87,   155,
     156,    87,   154,    87,   177,   161,   162,   163,   164,   165,
     166,   167,   168,   169,   170,    87,    87,    87,    87,   154,
      87,   177,   154,   147,    87,   149,   150,   151,   152,   153,
     750,   155,   156,   157,   158,     4,   160,    43,   211,   154,
     213,   154,   154,   210,  1182,   154,   154,    23,    54,   155,
     156,   771,   154,   154,   154,   211,    87,   213,   210,    32,
      32,    87,   154,    32,   154,   154,     7,    26,   154,   210,
     154,   154,   154,   210,   154,   154,   154,   212,   154,   154,
     212,    87,   154,   154,   154,   154,   210,   154,    43,    43,
      87,   154,   211,    32,    84,   210,   210,   210,   210,    54,
      54,   210,   210,   210,   218,    46,    96,    48,   210,   210,
      69,    87,   212,    54,    87,    87,   211,   837,    87,   109,
     212,   210,   212,   154,   210,    66,   210,  1115,   154,   212,
     212,    63,   212,   212,   212,   587,   212,   212,   210,   210,
     210,   210,  1101,   210,  1103,   617,   105,   210,    87,   155,
     126,   577,    32,    34,   577,    34,    87,   617,   577,   617,
     617,   617,   617,   617,    43,   617,  1125,   617,   617,   617,
     617,   617,  1115,   617,   655,  1118,   617,    23,   655,   210,
     617,   617,   212,    87,   210,    23,    32,    23,   218,  1114,
      87,   617,   124,   125,    32,    34,    32,   210,   617,   592,
     617,   617,   595,   135,    43,   210,   210,    87,   714,   715,
     930,   717,   216,    24,   720,   721,   722,   723,   724,   725,
     726,   727,   728,   729,   730,   618,   732,   733,   734,   735,
     736,   737,   738,   739,   740,   741,   742,  1196,   212,  1198,
    1199,    87,   962,  1208,   218,   210,  1211,  1212,  1213,    87,
    1208,    87,   754,  1211,  1212,  1213,   123,   157,   158,   147,
     160,   149,   150,   151,   152,   153,    43,   155,   156,   157,
     158,   760,   160,   155,   156,   157,   158,   210,   160,    87,
     126,    87,   754,  1242,   744,  1244,    87,   744,   126,    34,
     126,    87,   157,   158,   744,   160,   744,    39,  1257,    43,
    1020,   210,   211,   210,   211,   210,   758,    10,    11,    12,
      13,    14,    15,    16,   210,   211,    87,   710,  1277,    22,
    1285,  1286,  1281,    26,   212,   210,   211,  1285,  1286,   155,
     156,   157,   158,    87,   160,   210,   211,    32,   132,    42,
    1299,    54,   210,   215,   132,   132,   132,    87,    43,    34,
      34,    34,   745,   746,    32,    34,   154,    34,    34,    54,
      34,    34,  1321,    34,    34,    43,    34,  1326,   761,    34,
      34,   764,   212,    68,    34,    34,    54,    34,   154,   212,
     212,    34,    77,  1308,    34,    87,    89,    90,   210,   210,
      68,    87,    87,    87,   900,    87,    34,    34,    34,    77,
      34,    34,    34,    34,   107,   108,    34,    34,    34,    87,
      34,  1131,  1132,  1133,  1134,   895,  1136,    34,   147,    34,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
    1150,   160,    34,    34,    34,    34,    87,   154,   141,   142,
      87,   147,    87,   149,   150,   151,   152,   153,   143,   155,
     156,   157,   158,   154,   160,  1175,   154,    87,    87,  1179,
     155,   156,    87,  1183,    87,   143,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   211,   155,   156,    34,
     212,  1201,   177,   161,   162,   163,   164,   165,   166,   167,
     168,   169,   170,   212,   210,    34,  1216,   212,  1218,   177,
      87,    32,    87,   211,   211,    34,   212,   211,    34,    34,
      34,   211,    43,  1111,   211,    34,   211,    34,   213,    34,
     211,   211,   211,    54,   211,   211,    34,    34,    34,    32,
     211,   211,    34,   211,    34,   213,   211,    68,   132,   211,
      43,    87,   211,    34,   211,   211,    77,    34,  1268,   212,
     211,    54,   210,    34,    34,  1275,    87,   211,   147,   211,
     149,   150,   151,   152,   153,    68,   155,   156,   157,   158,
     211,   160,   211,  1293,    77,    10,    11,    12,    13,    14,
      15,    16,    17,   132,    87,   132,    21,    22,    87,   211,
     147,    26,   149,   150,   151,   152,   153,  1317,   155,   156,
     157,   158,   211,   160,   211,   211,  1112,    42,   211,    44,
      87,   132,   143,   211,   211,   211,   211,    87,   154,    34,
      34,    34,    34,   212,   155,   156,    34,    34,    34,    34,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   170,
     143,    34,    34,    34,    34,    34,   177,  1153,  1154,    34,
    1156,    34,   155,   156,    89,    90,    34,    34,   161,   162,
     163,   164,   165,   166,   167,   168,   169,   170,    34,    34,
      34,    34,   107,   108,   177,    34,    34,    34,    34,    34,
     211,    34,   213,   154,   212,    54,    87,    54,   210,    87,
      87,   211,    54,     0,  1087,    34,  1089,   210,   210,    87,
      54,     8,     9,    54,    54,    54,   141,   142,   211,    54,
     213,    18,    19,    20,    54,    87,    54,   210,    25,   215,
      54,    54,    29,    30,    31,    54,    33,    87,    35,    36,
      37,    38,    54,    54,    54,    54,    54,    54,   211,    87,
      87,    34,   210,    50,   154,   154,    53,   210,   147,   210,
     149,   150,   151,   152,   153,    87,   155,   156,   157,   158,
      67,   160,    87,    87,   160,    54,  1272,   154,   154,    54,
     154,    54,   154,    54,    81,    82,    83,    54,    54,   147,
      87,   149,   150,   151,   152,   153,    54,   155,   156,   157,
     158,    98,   160,   100,    54,   102,   103,   104,   105,   106,
      54,    54,   132,    54,   111,   112,   210,   210,   160,   116,
     117,   118,   119,   212,    54,   210,   132,   132,   132,    54,
     127,   128,    54,   130,     8,     9,    54,   134,    54,    54,
     137,   138,   139,   140,    18,    19,    20,   144,    54,    54,
      54,    25,    54,    54,   212,    29,    30,    31,    54,    33,
      54,    35,    36,    37,    38,    54,    54,    43,    54,    54,
      54,    54,    54,    43,   210,    43,    50,    54,    54,    53,
     214,   178,    54,   210,   210,   217,   211,    87,    87,   154,
     154,   212,    68,    67,   212,   212,   212,   210,   154,    54,
      87,    77,    54,    34,    34,   160,    87,    81,    82,    83,
     210,    87,   147,    87,   149,   150,   151,   152,   153,   212,
     155,   156,   157,   158,    98,   160,   100,    87,   102,   103,
     104,   105,   106,   210,   210,    34,    87,   111,   112,    87,
     212,   210,   116,   117,   118,   119,   214,    87,    87,   212,
     210,   154,    34,   127,   128,    34,   130,   154,    34,   154,
     134,    54,   210,   137,   138,   139,   140,   143,    54,    34,
     144,    43,   217,   210,   154,   154,   154,   212,   210,   155,
     156,   210,    54,   214,   212,   161,   162,   163,   164,   165,
     166,   167,   168,   169,   170,   154,    68,   154,   210,    43,
     210,   177,    54,   214,   178,    77,    54,   136,   154,   154,
      54,    54,    87,    54,   212,    87,   154,   154,   210,   154,
      45,   214,   323,  1118,    68,   140,   217,  1015,    43,   758,
     587,   754,   563,    77,   320,   211,   617,   213,   716,    54,
     138,   315,   744,    87,   655,  1190,   147,   309,   149,   150,
     151,   152,   153,    68,   155,   156,   157,   158,   214,   160,
    1123,   802,    77,   763,   572,   501,   122,   370,   818,    -1,
      -1,   143,    87,   577,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   154,   155,   156,    -1,    -1,    -1,    -1,   161,
     162,   163,   164,   165,   166,   167,   168,   169,   170,   143,
      -1,    -1,    -1,    -1,    -1,   177,    -1,    -1,    -1,   210,
     154,   155,   156,    -1,    -1,    -1,    -1,   161,   162,   163,
     164,   165,   166,   167,   168,   169,   170,    -1,   143,    -1,
      -1,    43,    -1,   177,    -1,    -1,    -1,    -1,    -1,   211,
     155,   156,    54,    -1,    -1,    -1,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,    68,    -1,    -1,    -1,
      -1,    -1,   177,    -1,    -1,    77,    -1,   211,   147,    -1,
     149,   150,   151,   152,   153,    87,   155,   156,   157,   158,
      -1,   160,   147,    -1,   149,   150,   151,   152,   153,    -1,
     155,   156,   157,   158,    -1,   160,   211,    -1,   147,    -1,
     149,   150,   151,   152,   153,    -1,   155,   156,   157,   158,
     147,   160,   149,   150,   151,   152,   153,    -1,   155,   156,
     157,   158,    -1,   160,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   143,    -1,   212,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   155,   156,     5,    -1,   212,    -1,   161,
     162,   163,   164,   165,   166,   167,   168,   169,   170,    -1,
      -1,    21,    22,   212,    -1,   177,    26,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   212,    -1,    -1,    -1,    -1,
      40,    41,    42,    -1,    44,    -1,    -1,    -1,    -1,    49,
      -1,    -1,    -1,    -1,    -1,    -1,    56,    57,    58,   211,
      -1,    61,    62,    -1,    64,    65,    -1,     5,    -1,    -1,
      70,    71,    72,    73,    74,    75,    76,    -1,    78,    79,
      80,    -1,    -1,    21,    -1,    -1,    86,    -1,    -1,    89,
      90,    -1,    92,    -1,    94,    -1,    -1,    -1,    -1,    99,
      -1,    -1,    40,    41,    -1,    -1,    44,   107,   108,     3,
     110,    49,    -1,    -1,    -1,    -1,    -1,    -1,    56,    -1,
      -1,    -1,   122,   123,    -1,    -1,    -1,    -1,    -1,   129,
      -1,    25,    -1,    27,    28,    -1,    -1,    -1,    -1,    -1,
      -1,   141,   142,    -1,    -1,    -1,    -1,    -1,    86,    -1,
      -1,    -1,    -1,    -1,    -1,    49,    94,    51,    52,    -1,
      -1,    -1,    56,    -1,     3,    -1,    -1,    -1,    -1,    63,
      -1,    -1,   110,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   122,    -1,    25,    -1,    27,    28,
      -1,   129,    -1,    -1,    -1,    -1,    -1,    91,    -1,    93,
      94,    95,    96,   141,   142,    -1,    -1,   101,    -1,    -1,
      49,   105,    51,    52,    -1,   109,    -1,    56,    -1,   113,
     114,   115,    -1,    -1,    63,   119,   120,   121,    -1,   123,
     147,    -1,   149,   150,   151,   152,   153,   131,   155,   156,
     157,   158,    -1,   160,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    91,    -1,    93,    94,    95,    96,    -1,    -1,
      -1,    -1,   101,    -1,    -1,    -1,   105,    -1,    -1,    -1,
     109,    -1,    -1,    -1,   113,   114,   115,    -1,    -1,    -1,
     119,   120,   121,   147,   123,   149,   150,   151,   152,   153,
      -1,   155,   156,   157,   158,   212,   160,   147,    -1,   149,
     150,   151,   152,   153,    -1,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,    -1,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,    -1,   155,
     156,   157,   158,   147,   160,   149,   150,   151,   152,   153,
      -1,   155,   156,   157,   158,    -1,   160,   147,   212,   149,
     150,   151,   152,   153,    -1,   155,   156,   157,   158,    -1,
     160,    -1,   212,    -1,    -1,    -1,    -1,   147,    -1,   149,
     150,   151,   152,   153,   212,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,   212,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,   212,   155,
     156,   157,   158,    -1,   160,   147,    -1,   149,   150,   151,
     152,   153,   212,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,    -1,   155,   156,   157,   158,    -1,
     160,    -1,   212,    -1,    -1,    -1,    -1,   147,    -1,   149,
     150,   151,   152,   153,   212,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,   212,   155,   156,   157,
     158,    -1,   160,   147,    -1,   149,   150,   151,   152,   153,
     212,   155,   156,   157,   158,   147,   160,   149,   150,   151,
     152,   153,   212,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,    -1,   155,   156,   157,   158,    -1,
     160,    -1,   212,    -1,    -1,    -1,    -1,   147,    -1,   149,
     150,   151,   152,   153,   212,   155,   156,   157,   158,    -1,
     160,   147,    -1,   149,   150,   151,   152,   153,   212,   155,
     156,   157,   158,   147,   160,   149,   150,   151,   152,   153,
     212,   155,   156,   157,   158,   147,   160,   149,   150,   151,
     152,   153,   212,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,    -1,   155,   156,   157,   158,    -1,
     160,   147,   212,   149,   150,   151,   152,   153,    -1,   155,
     156,   157,   158,    -1,   160,    -1,   212,    -1,    -1,    -1,
      -1,   147,    -1,   149,   150,   151,   152,   153,   212,   155,
     156,   157,   158,   147,   160,   149,   150,   151,   152,   153,
     212,   155,   156,   157,   158,   147,   160,   149,   150,   151,
     152,   153,   212,   155,   156,   157,   158,    -1,   160,   147,
      -1,   149,   150,   151,   152,   153,   212,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,    -1,   155,
     156,   157,   158,    -1,   160,    -1,   212,    -1,    -1,   147,
      -1,   149,   150,   151,   152,   153,   210,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,   210,   155,
     156,   157,   158,    -1,   160,   147,    -1,   149,   150,   151,
     152,   153,   210,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,   210,   155,   156,   157,   158,    -1,
     160,    -1,    -1,    -1,    -1,   147,    -1,   149,   150,   151,
     152,   153,   210,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,   210,   155,   156,   157,   158,    -1,
     160,   147,    -1,   149,   150,   151,   152,   153,   210,   155,
     156,   157,   158,   147,   160,   149,   150,   151,   152,   153,
     210,   155,   156,   157,   158,    -1,   160,    -1,    -1,    -1,
      -1,   147,    -1,   149,   150,   151,   152,   153,   210,   155,
     156,   157,   158,   147,   160,   149,   150,   151,   152,   153,
     210,   155,   156,   157,   158,    -1,   160,   147,    -1,   149,
     150,   151,   152,   153,   210,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,   210,   155,   156,   157,
     158,    -1,   160,    -1,    -1,    -1,    -1,   147,    -1,   149,
     150,   151,   152,   153,   210,   155,   156,   157,   158,    -1,
     160,    -1,    -1,    -1,    -1,   147,   210,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,    -1,   160,    -1,
     210,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   210,   179,   180,   181,   182,   183,   184,   185,
     186,   187,   188,   189,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   199,    -1,    -1,    -1,   203,   204,    -1,
     206,   207,   208,   209,   147,    -1,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   147,   160,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,   147,   160,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
     147,   160,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,   147,   160,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   147,   160,   149,   150,   151,   152,
     153,   154,   155,   156,   157,   158,   147,   160,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,    -1,   160
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     8,     9,    18,    19,    20,    25,    29,    30,    31,
      33,    35,    36,    37,    38,    50,    53,    67,    81,    82,
      83,    87,    98,   100,   102,   103,   104,   105,   106,   111,
     112,   116,   117,   118,   119,   127,   128,   130,   134,   137,
     138,   139,   140,   144,   178,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   234,   235,   236,   237,   240,   242,
     245,   253,   266,   267,   271,   275,   278,   281,   284,   290,
     296,   299,   304,   307,   310,   313,   314,   317,   318,   320,
     321,   322,   326,   327,   328,   329,   335,   338,   344,   347,
     348,   351,    54,   211,    54,   211,   210,   211,   210,   210,
     211,    34,    43,    54,    87,   211,    87,   211,   210,    87,
     210,   211,   287,   210,   210,   210,   210,   210,   211,    34,
      43,   210,   211,   211,   210,    34,   210,   210,   210,   211,
     287,   287,    87,   233,    34,    54,   336,   211,   211,   287,
     210,    34,   210,   211,   210,   211,   210,   211,   287,   210,
     211,   287,   287,    87,   230,    87,   231,    87,   232,   287,
     210,   210,   211,     0,   222,   210,    10,    11,    12,    13,
      14,    15,    16,    22,    26,    42,    89,    90,   107,   108,
     141,   142,   341,   342,   343,   376,   377,   378,   379,   380,
     411,   412,   414,   415,   418,   419,   420,   421,   422,   423,
     424,   210,    17,    21,    44,   342,   345,   346,   384,   401,
     425,    24,     4,    87,   323,   324,   123,   279,   280,   356,
      43,   210,    54,   210,   210,   218,    87,   210,   218,    87,
      87,   243,   244,    34,     5,    40,    41,    49,    56,    57,
      58,    61,    62,    64,    65,    70,    71,    72,    73,    74,
      75,    76,    78,    79,    80,    86,    92,    94,    99,   110,
     122,   129,   305,   306,   356,   366,   376,   377,   378,   379,
     380,   381,   382,   383,   384,   385,   386,   387,   388,   389,
     390,   391,   392,   393,   394,   395,   396,   397,   398,   399,
     400,   401,   402,   403,   404,   406,   407,   411,   412,   413,
     414,   415,   416,    87,   154,   210,    23,    87,   126,   291,
     292,   293,    23,    87,   126,   300,   301,    23,    87,   126,
     297,   298,    87,   246,   247,   243,    39,   241,    43,   210,
     254,    63,   124,   125,   135,   358,    84,    96,   109,   330,
     331,   408,   409,   410,    23,   137,   268,   269,    43,    54,
      68,    77,    87,   143,   155,   156,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   177,   211,   238,    87,
     315,   316,    87,   319,     3,    25,    27,    28,    51,    52,
      91,    93,    95,   101,   105,   113,   114,   115,   119,   120,
     121,   286,   355,   356,   357,   358,   359,   360,   361,   362,
     363,   364,   365,   366,   367,   368,   369,   370,   373,   374,
     375,   383,   405,   409,   410,   210,   210,   132,    87,   154,
     210,    54,   210,    43,    54,    68,    77,    87,   143,   155,
     156,   161,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   177,   211,   263,   265,   308,   309,   366,   383,   384,
     393,   399,   400,   401,   402,   403,   404,   411,   412,   413,
     308,   210,   268,   215,    69,   282,   283,   369,   371,   376,
     145,   146,   276,   277,   356,   453,   454,   285,   286,   210,
     131,   286,   339,   340,   417,   210,   210,   132,    87,   154,
     210,   132,    87,   154,   210,   132,    87,   154,   210,   210,
      87,   352,   353,   179,   180,   181,   182,   183,   184,   185,
     186,   187,   188,   189,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   199,   203,   204,   206,   207,   208,   209,
     349,   350,   426,   427,   428,   429,   430,   431,   432,   433,
     434,   435,   436,   437,   438,   439,   440,   441,   442,   443,
     444,   445,   446,   447,   448,   449,   450,   451,   452,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,   154,   212,    34,    34,    34,   154,   212,   212,
      87,   154,   211,   325,    32,   324,    34,   154,   212,   210,
     210,    87,   212,   218,    87,   212,   218,    34,    32,   244,
      87,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,   154,   212,    34,
      87,    87,    87,    32,   292,   154,    87,   154,    87,    32,
     301,    87,   154,    87,    32,   298,   211,    32,   247,    32,
      34,   212,   210,   213,   261,   262,   263,   264,   154,   212,
     154,   212,   212,   212,    34,   154,   212,    87,    87,    32,
     269,   211,   211,   211,   211,   238,   238,   211,   211,   211,
     211,   211,   211,   211,   211,   211,   211,   211,   238,   147,
     149,   150,   151,   152,   153,   155,   156,   157,   158,   160,
     210,   211,    32,   316,   154,   238,    32,    87,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
     212,   132,    87,   210,   211,   211,   211,   211,   263,   263,
     211,   211,   211,   211,   211,   211,   211,   211,   211,   211,
     211,   263,   147,   149,   150,   151,   152,   153,   155,   156,
     157,   158,   160,   337,   154,   212,   212,    32,    43,    54,
     211,   273,   274,    34,   154,   212,    34,    34,   154,   212,
     154,   212,    34,   154,   212,   132,    87,   132,    87,   132,
      87,   154,    32,   353,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,   212,   154,    43,    54,   354,    43,   155,   156,   289,
     354,    54,   289,    54,    87,    54,    54,   215,   457,   458,
      54,    54,    87,    87,   455,   343,    54,    54,   354,    54,
     346,    54,   210,   211,    87,    43,    54,    34,    54,   280,
     210,   210,   210,   287,    87,   210,   210,   287,    87,   238,
     458,    54,    54,    54,    54,    54,   354,   354,   354,    54,
      54,    54,    54,    87,   211,   354,   306,   210,   287,    87,
      34,   154,     6,    43,    45,    54,    55,    87,    97,   133,
     155,   294,   302,   303,   154,   303,   154,   154,   303,   154,
      54,   155,   156,   288,    87,   210,    87,    32,   262,   264,
      34,   210,    67,   249,   252,   370,   372,   260,    47,    60,
      88,   248,   250,   251,   252,   258,   210,   210,    59,    85,
     331,    87,   157,   214,   218,   219,   332,   333,   334,   154,
      34,   154,   210,   238,   238,   238,   239,   238,   238,   238,
     238,   238,   238,   238,   238,   238,   238,   238,   238,   212,
     238,   238,   238,   238,   238,   238,   238,   238,   238,   238,
     238,   238,    87,   210,   154,   238,    54,   354,    54,    54,
      54,    54,    54,    54,   354,    54,    54,    54,   210,   287,
     132,   263,   263,   288,   263,   263,   263,   263,   263,   263,
     263,   263,   263,   263,   263,   263,   212,   263,   263,   263,
     263,   263,   263,   263,   263,   263,   263,   263,   210,   309,
     210,   287,   210,   287,   238,   210,   216,    43,    54,   154,
     211,     7,    46,    48,    54,    66,   283,   210,    54,    54,
     277,   210,   286,   210,   287,   354,   340,   210,   287,   132,
     132,   132,   238,    54,    54,    54,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,   354,   354,
      54,   458,   354,   354,    54,    54,   354,    54,   354,   354,
     210,   349,    43,    43,    54,   456,   216,   456,   214,   210,
     210,    54,   325,   212,   212,   238,   210,   212,   210,   212,
     210,   217,   311,   312,   210,    87,    87,    43,    54,   210,
     154,   154,    87,   154,   303,    87,   210,   303,    54,    54,
     212,   243,    34,   263,    34,   154,   212,   210,   154,   212,
     210,   256,   255,   154,   210,   211,   334,    87,   238,    87,
     105,   126,   154,   154,   154,   212,   154,   212,   212,   212,
     212,   212,   212,   212,   212,   212,   212,   212,   212,   238,
      87,   210,   210,   154,   154,   212,   154,   212,   212,   212,
     212,   212,   212,   212,   212,   212,   212,   212,   210,   210,
     212,   274,   210,    43,    54,   211,   238,   210,   210,   154,
     214,    87,   212,    34,   210,   210,   287,   210,   287,    87,
     154,   212,   295,   303,   302,   303,   154,   303,   154,   154,
     210,    34,    32,   263,   210,   354,   252,   259,   261,   250,
     257,   261,   261,   261,   333,   303,    34,   210,    34,    54,
     270,   238,   238,   238,   238,   238,   210,   210,   238,   263,
     263,   263,   238,   212,   238,    54,   325,   238,   210,   210,
     217,   311,   154,   154,   154,   303,   210,   303,   303,   238,
     210,   210,    32,   210,    32,    32,    32,   211,   212,   238,
     238,   214,    54,   154,   210,   210,   212,   212,   154,   210,
     212,   212,   154,   212,   210,    34,   210,   154,   303,   295,
     303,   154,   210,   210,   210,   261,   261,   303,   210,   210,
      54,   214,    54,   136,   238,   263,   238,   217,   303,   154,
     154,   303,    32,    32,   212,    54,   214,   238,   272,   212,
     212,   210,    87,   303,   302,   210,    54,   154,   210,   354,
     217,   154,   154,   238,   303,   295,   154,   303
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
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
      59,    40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   220,   221,   221,   222,   222,   222,   222,   222,   222,
     222,   222,   222,   222,   222,   222,   222,   222,   222,   222,
     222,   222,   222,   222,   222,   222,   222,   222,   222,   222,
     222,   222,   222,   222,   222,   222,   222,   222,   222,   222,
     222,   222,   222,   222,   222,   223,   223,   223,   223,   224,
     224,   225,   226,   227,   228,   229,   230,   230,   230,   230,
     230,   230,   231,   231,   231,   231,   231,   231,   232,   232,
     232,   232,   232,   232,   233,   233,   233,   233,   233,   233,
     234,   234,   235,   235,   236,   236,   237,   238,   238,   238,
     238,   238,   238,   238,   238,   238,   238,   238,   238,   238,
     238,   238,   238,   238,   238,   238,   238,   238,   238,   238,
     238,   238,   238,   238,   238,   238,   238,   238,   238,   239,
     239,   240,   240,   241,   242,   243,   243,   244,   245,   246,
     246,   247,   248,   248,   249,   249,   250,   250,   251,   251,
     251,   252,   252,   254,   253,   255,   253,   256,   253,   257,
     253,   258,   253,   259,   253,   260,   253,   261,   261,   261,
     261,   262,   262,   263,   263,   263,   263,   263,   263,   263,
     263,   263,   263,   263,   263,   263,   263,   263,   263,   263,
     263,   263,   263,   263,   263,   263,   263,   263,   263,   263,
     263,   263,   263,   263,   264,   265,   265,   266,   267,   268,
     268,   269,   269,   269,   269,   269,   270,   270,   270,   270,
     270,   270,   271,   272,   272,   272,   273,   273,   274,   274,
     274,   274,   274,   274,   274,   274,   274,   275,   275,   276,
     276,   277,   277,   277,   278,   278,   279,   279,   280,   281,
     281,   282,   282,   283,   283,   283,   284,   284,   284,   284,
     285,   285,   286,   286,   286,   286,   286,   286,   286,   286,
     286,   286,   286,   286,   286,   286,   286,   286,   286,   286,
     286,   286,   286,   286,   286,   287,   287,   287,   287,   287,
     287,   288,   288,   288,   289,   289,   289,   290,   291,   291,
     292,   293,   293,   293,   294,   294,   294,   294,   294,   295,
     295,   295,   295,   296,   297,   297,   298,   298,   298,   299,
     300,   300,   301,   301,   301,   302,   302,   302,   302,   302,
     303,   303,   303,   303,   303,   303,   304,   304,   304,   304,
     305,   305,   306,   306,   306,   306,   306,   306,   306,   306,
     306,   306,   306,   306,   306,   306,   306,   306,   306,   306,
     306,   306,   306,   306,   306,   306,   306,   306,   306,   306,
     306,   306,   306,   306,   306,   306,   306,   306,   306,   306,
     306,   307,   307,   308,   308,   309,   309,   309,   309,   309,
     309,   309,   309,   309,   309,   309,   309,   309,   310,   310,
     311,   311,   312,   312,   313,   314,   315,   315,   316,   317,
     318,   319,   319,   319,   319,   320,   321,   321,   321,   321,
     322,   323,   323,   324,   324,   324,   325,   325,   325,   326,
     326,   327,   327,   327,   327,   327,   327,   328,   328,   328,
     328,   328,   328,   329,   330,   330,   331,   331,   331,   332,
     332,   332,   332,   333,   333,   334,   334,   334,   334,   334,
     336,   337,   335,   338,   338,   338,   338,   339,   339,   340,
     340,   341,   341,   341,   341,   341,   341,   341,   342,   342,
     342,   342,   342,   342,   342,   342,   342,   342,   343,   343,
     344,   344,   345,   345,   345,   345,   346,   346,   347,   347,
     348,   348,   349,   349,   350,   350,   350,   350,   350,   350,
     350,   350,   350,   350,   350,   350,   350,   350,   350,   350,
     350,   350,   350,   350,   350,   350,   350,   350,   350,   350,
     350,   351,   352,   352,   353,   354,   354,   355,   356,   357,
     358,   359,   360,   361,   362,   363,   364,   365,   366,   367,
     368,   369,   370,   371,   371,   371,   371,   371,   372,   373,
     374,   375,   376,   377,   377,   378,   379,   380,   381,   382,
     383,   383,   384,   385,   386,   387,   388,   389,   390,   391,
     392,   393,   394,   395,   396,   397,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   408,   409,   410,
     411,   412,   413,   414,   415,   416,   417,   418,   419,   420,
     421,   422,   423,   424,   425,   426,   427,   428,   429,   430,
     431,   432,   433,   434,   435,   436,   437,   438,   439,   440,
     441,   442,   443,   444,   445,   446,   447,   448,   449,   450,
     451,   452,   453,   454,   455,   456,   456,   457,   457,   458
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  parser::yyr2_[] =
  {
         0,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       4,     3,     3,     3,     3,     3,     2,     3,     1,     3,
       4,     2,     2,     3,     1,     3,     4,     2,     2,     3,
       1,     3,     4,     2,     2,     3,     1,     3,     4,     2,
       3,     4,     3,     4,     3,     4,     4,     3,     1,     1,
       1,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     2,     2,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     6,     6,     4,     8,     1,
       3,     4,     7,     3,     4,     2,     1,     4,     4,     2,
       1,     7,     3,     1,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     0,     5,     0,     8,     0,     8,     0,
      10,     0,     8,     0,    10,     0,     8,     2,     2,     1,
       1,     4,     2,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     2,     2,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     6,     6,     8,     5,     1,     4,     4,     4,     2,
       1,     9,     6,     5,     7,     7,     3,     2,     5,     4,
       3,     1,     6,     3,     2,     1,     3,     1,     5,     3,
       3,     4,     2,     2,     3,     1,     1,     2,     5,     3,
       1,     1,     1,     1,     2,     5,     3,     1,     1,     2,
       5,     3,     1,     1,     1,     1,     2,     5,     3,     6,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     2,     4,     3,     5,     1,
       3,     2,     2,     1,     2,     2,     1,     4,     2,     1,
       4,     2,     1,     4,     3,     5,     9,     1,     5,     3,
       5,     7,     9,     4,     2,     1,     5,     7,     4,     4,
       2,     1,     7,     9,     6,     1,     1,     1,     1,     1,
       0,     1,     1,     1,     2,     2,     2,     5,     3,     6,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     5,     6,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     5,     6,
       7,     5,     1,     3,     3,     4,     2,     1,     5,     3,
       4,     4,     6,     3,     5,     3,     2,     5,     3,     6,
       4,     2,     1,     5,     7,     9,     0,     3,     3,     2,
       5,     5,     6,     3,     7,     8,     5,     5,     6,     3,
       7,     8,     5,     6,     3,     1,     1,     1,     1,     1,
       3,     4,     6,     1,     2,     1,     1,     1,     1,     1,
       0,     0,     5,     2,     5,     3,     6,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     1,
       3,     6,     1,     1,     1,     1,     3,     1,     3,     6,
       2,     5,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     4,     1,     2,     6,     1,     1,     3,     3,     3,
       1,     3,     3,     3,     3,     1,     1,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       1,     1,     3,     3,     3,     3,     5,     3,     3,     3,
       1,     3,     3,     3,     1,     1,     1,     1,     1,     3,
       1,     1,     1,     1,     3,     3,     3,     3,     1,     1,
       3,     3,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     3,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     1,     3,     2,     2,     2
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const parser::yytname_[] =
  {
    "$end", "error", "$undefined", "AR", "AUTOCORR", "BAYESIAN_IRF",
  "BETA_PDF", "BICGSTAB", "BVAR_DENSITY", "BVAR_FORECAST",
  "BVAR_PRIOR_DECAY", "BVAR_PRIOR_FLAT", "BVAR_PRIOR_LAMBDA",
  "BVAR_PRIOR_MU", "BVAR_PRIOR_OMEGA", "BVAR_PRIOR_TAU",
  "BVAR_PRIOR_TRAIN", "BVAR_REPLIC", "CALIB", "CALIB_VAR", "CHECK",
  "CONF_SIG", "CONSTANT", "CORR", "COVAR", "CUTOFF", "DATAFILE", "DR_ALGO",
  "DROP", "DSAMPLE", "DYNASAVE", "DYNATYPE", "END", "ENDVAL", "EQUAL",
  "ESTIMATION", "ESTIMATED_PARAMS", "ESTIMATED_PARAMS_BOUNDS",
  "ESTIMATED_PARAMS_INIT", "FILENAME", "FILTER_STEP_AHEAD",
  "FILTERED_VARS", "FIRST_OBS", "FLOAT_NUMBER", "FORECAST", "GAMMA_PDF",
  "GAUSSIAN_ELIMINATION", "GCC_COMPILER", "GMRES", "GRAPH", "HISTVAL",
  "HP_FILTER", "HP_NGRID", "INITVAL", "INT_NUMBER", "INV_GAMMA_PDF", "IRF",
  "KALMAN_ALGO", "KALMAN_TOL", "LAPLACE", "LCC_COMPILER", "LIK_ALGO",
  "LIK_INIT", "LINEAR", "LOAD_MH_FILE", "LOGLINEAR", "LU", "MARKOWITZ",
  "MAX", "METHOD", "MH_DROP", "MH_INIT_SCALE", "MH_JSCALE", "MH_MODE",
  "MH_NBLOCKS", "MH_REPLIC", "MH_RECOVER", "MIN", "MODE_CHECK",
  "MODE_COMPUTE", "MODE_FILE", "MODEL", "MODEL_COMPARISON", "MSHOCKS",
  "MODEL_COMPARISON_APPROXIMATION", "MODIFIEDHARMONICMEAN",
  "MOMENTS_VARENDO", "NAME", "NO_COMPILER", "NOBS", "NOCONSTANT", "NOCORR",
  "NODIAGNOSTIC", "NOFUNCTIONS", "NOGRAPH", "NOMOMENTS", "NOPRINT",
  "NORMAL_PDF", "OBSERVATION_TRENDS", "OPTIM", "OPTIM_WEIGHTS", "ORDER",
  "OSR", "OSR_PARAMS", "PARAMETERS", "PERIODS", "PLANNER_OBJECTIVE",
  "PREFILTER", "PRESAMPLE", "PRINT", "PRIOR_TRUNC", "PRIOR_ANALYSIS",
  "POSTERIOR_ANALYSIS", "QZ_CRITERIUM", "RELATIVE_IRF", "REPLIC", "RPLOT",
  "SHOCKS", "SIGMA_E", "SIMUL", "SIMUL_ALGO", "SIMUL_SEED", "SMOOTHER",
  "SOLVE_ALGO", "SPARSE", "SPARSE_DLL", "STDERR", "STEADY", "STOCH_SIMUL",
  "TEX", "RAMSEY_POLICY", "PLANNER_DISCOUNT", "TEX_NAME", "UNIFORM_PDF",
  "UNIT_ROOT_VARS", "USE_DLL", "VALUES", "VAR", "VAREXO", "VAREXO_DET",
  "VAROBS", "XLS_SHEET", "XLS_RANGE", "NORMCDF", "HOMOTOPY_SETUP",
  "HOMOTOPY_MODE", "HOMOTOPY_STEPS", "EXCLAMATION_EQUAL", "EXCLAMATION",
  "EQUAL_EQUAL", "GREATER_EQUAL", "LESS_EQUAL", "GREATER", "LESS", "COMMA",
  "MINUS", "PLUS", "DIVIDE", "TIMES", "UMINUS", "POWER", "EXP", "LOG",
  "LN", "LOG10", "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN", "SINH",
  "COSH", "TANH", "ASINH", "ACOSH", "ATANH", "SQRT", "DYNARE_SENSITIVITY",
  "IDENTIFICATION", "MORRIS", "STAB", "REDFORM", "PPRIOR", "PRIOR_RANGE",
  "PPOST", "ILPTAU", "GLUE", "MORRIS_NLIV", "MORRIS_NTRA", "NSAM",
  "LOAD_REDFORM", "LOAD_RMSE", "LOAD_STAB", "ALPHA2_STAB", "KSSTAT",
  "LOGTRANS_REDFORM", "THRESHOLD_REDFORM", "KSSTAT_REDFORM",
  "ALPHA2_REDFORM", "NAMENDO", "NAMLAGENDO", "NAMEXO", "RMSE", "LIK_ONLY",
  "VAR_RMSE", "PFILT_RMSE", "ISTART_RMSE", "ALPHA_RMSE", "ALPHA2_RMSE",
  "';'", "'('", "')'", "'#'", "':'", "'['", "']'", "'''", "'.'", "'\\\\'",
  "$accept", "statement_list", "statement", "declaration", "dsample",
  "rplot", "var", "varexo", "varexo_det", "parameters", "var_list",
  "varexo_list", "varexo_det_list", "parameter_list", "periods", "cutoff",
  "markowitz", "init_param", "expression", "comma_expression", "initval",
  "initval_option", "endval", "initval_list", "initval_elem", "histval",
  "histval_list", "histval_elem", "model_sparse_dll_options_list",
  "model_sparse_options_list", "model_sparse_dll_options",
  "model_compiler_options", "model_sparse_common_options", "model", "@1",
  "@2", "@3", "@4", "@5", "@6", "@7", "equation_list", "equation",
  "hand_side", "pound_expression", "model_var", "shocks", "mshocks",
  "shock_list", "shock_elem", "period_list", "sigma_e", "value_list",
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
  "filename", "filename_elem", "planner_objective", "@8", "@9",
  "ramsey_policy", "ramsey_policy_options_list", "ramsey_policy_options",
  "bvar_prior_option", "bvar_common_option", "bvar_density_options_list",
  "bvar_density", "bvar_forecast_option", "bvar_forecast_options_list",
  "bvar_forecast", "dynare_sensitivity", "dynare_sensitivity_options_list",
  "dynare_sensitivity_option", "homotopy_setup", "homotopy_list",
  "homotopy_item", "number", "o_dr_algo", "o_solve_algo", "o_simul_algo",
  "o_linear", "o_order", "o_replic", "o_drop", "o_ar", "o_nocorr",
  "o_nofunctions", "o_nomoments", "o_irf", "o_hp_filter", "o_hp_ngrid",
  "o_periods", "o_cutoff", "o_method", "o_markowitz", "o_simul",
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
  "o_gsa_alpha2_rmse", "o_homotopy_mode", "o_homotopy_steps", "range",
  "vec_int_elem", "vec_int_1", "vec_int", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const parser::rhs_number_type
  parser::yyrhs_[] =
  {
       221,     0,    -1,   222,    -1,   221,   222,    -1,   223,    -1,
     234,    -1,   235,    -1,   236,    -1,   253,    -1,   240,    -1,
     242,    -1,   245,    -1,   237,    -1,   266,    -1,   267,    -1,
     271,    -1,   275,    -1,   278,    -1,   281,    -1,   284,    -1,
     304,    -1,   307,    -1,   310,    -1,   290,    -1,   299,    -1,
     296,    -1,   313,    -1,   314,    -1,   317,    -1,   224,    -1,
     225,    -1,   318,    -1,   320,    -1,   321,    -1,   322,    -1,
     326,    -1,   327,    -1,   328,    -1,   329,    -1,   335,    -1,
     338,    -1,   344,    -1,   347,    -1,   348,    -1,   351,    -1,
     229,    -1,   226,    -1,   227,    -1,   228,    -1,    29,    54,
     210,    -1,    29,    54,    54,   210,    -1,   116,   287,   210,
      -1,   137,   230,   210,    -1,   138,   231,   210,    -1,   139,
     232,   210,    -1,   104,   233,   210,    -1,   230,    87,    -1,
     230,   154,    87,    -1,    87,    -1,   230,    87,   132,    -1,
     230,   154,    87,   132,    -1,    87,   132,    -1,   231,    87,
      -1,   231,   154,    87,    -1,    87,    -1,   231,    87,   132,
      -1,   231,   154,    87,   132,    -1,    87,   132,    -1,   232,
      87,    -1,   232,   154,    87,    -1,    87,    -1,   232,    87,
     132,    -1,   232,   154,    87,   132,    -1,    87,   132,    -1,
     233,    87,    -1,   233,   154,    87,    -1,    87,    -1,   233,
      87,   132,    -1,   233,   154,    87,   132,    -1,    87,   132,
      -1,   105,    54,   210,    -1,   105,    34,    54,   210,    -1,
      25,    43,   210,    -1,    25,    34,    43,   210,    -1,    67,
      43,   210,    -1,    67,    34,    43,   210,    -1,    87,    34,
     238,   210,    -1,   211,   238,   212,    -1,    87,    -1,    43,
      -1,    54,    -1,   238,   156,   238,    -1,   238,   155,   238,
      -1,   238,   157,   238,    -1,   238,   158,   238,    -1,   238,
     160,   238,    -1,   238,   153,   238,    -1,   238,   152,   238,
      -1,   238,   151,   238,    -1,   238,   150,   238,    -1,   238,
     149,   238,    -1,   238,   147,   238,    -1,   155,   238,    -1,
     156,   238,    -1,   161,   211,   238,   212,    -1,   162,   211,
     238,   212,    -1,   163,   211,   238,   212,    -1,   164,   211,
     238,   212,    -1,   165,   211,   238,   212,    -1,   166,   211,
     238,   212,    -1,   167,   211,   238,   212,    -1,   168,   211,
     238,   212,    -1,   169,   211,   238,   212,    -1,   170,   211,
     238,   212,    -1,   177,   211,   238,   212,    -1,    68,   211,
     238,   154,   238,   212,    -1,    77,   211,   238,   154,   238,
     212,    -1,    87,   211,   239,   212,    -1,   143,   211,   238,
     154,   238,   154,   238,   212,    -1,   238,    -1,   239,   154,
     238,    -1,    53,   210,   243,    32,    -1,    53,   211,   241,
     212,   210,   243,    32,    -1,    39,    34,    87,    -1,    33,
     210,   243,    32,    -1,   243,   244,    -1,   244,    -1,    87,
      34,   238,   210,    -1,    50,   210,   246,    32,    -1,   246,
     247,    -1,   247,    -1,    87,   211,   288,   212,    34,   238,
     210,    -1,   248,   154,   250,    -1,   250,    -1,   249,   154,
     252,    -1,   252,    -1,   251,    -1,   252,    -1,    60,    -1,
      47,    -1,    88,    -1,   370,    -1,   372,    -1,    -1,    81,
     210,   254,   261,    32,    -1,    -1,    81,   211,   358,   212,
     210,   255,   261,    32,    -1,    -1,    81,   211,   135,   212,
     210,   256,   261,    32,    -1,    -1,    81,   211,   125,   154,
     248,   212,   257,   210,   261,    32,    -1,    -1,    81,   211,
     125,   212,   258,   210,   261,    32,    -1,    -1,    81,   211,
     124,   154,   249,   212,   259,   210,   261,    32,    -1,    -1,
      81,   211,   124,   212,   260,   210,   261,    32,    -1,   261,
     262,    -1,   261,   264,    -1,   262,    -1,   264,    -1,   263,
      34,   263,   210,    -1,   263,   210,    -1,   211,   263,   212,
      -1,   265,    -1,    43,    -1,    54,    -1,   263,   156,   263,
      -1,   263,   155,   263,    -1,   263,   157,   263,    -1,   263,
     158,   263,    -1,   263,   153,   263,    -1,   263,   152,   263,
      -1,   263,   151,   263,    -1,   263,   150,   263,    -1,   263,
     149,   263,    -1,   263,   147,   263,    -1,   263,   160,   263,
      -1,   155,   263,    -1,   156,   263,    -1,   161,   211,   263,
     212,    -1,   162,   211,   263,   212,    -1,   163,   211,   263,
     212,    -1,   164,   211,   263,   212,    -1,   165,   211,   263,
     212,    -1,   166,   211,   263,   212,    -1,   167,   211,   263,
     212,    -1,   168,   211,   263,   212,    -1,   169,   211,   263,
     212,    -1,   170,   211,   263,   212,    -1,   177,   211,   263,
     212,    -1,    68,   211,   263,   154,   263,   212,    -1,    77,
     211,   263,   154,   263,   212,    -1,   143,   211,   263,   154,
     263,   154,   263,   212,    -1,   213,    87,    34,   263,   210,
      -1,    87,    -1,    87,   211,   288,   212,    -1,   117,   210,
     268,    32,    -1,    83,   210,   268,    32,    -1,   268,   269,
      -1,   269,    -1,   137,    87,   210,   105,   270,   210,   136,
     272,   210,    -1,   137,    87,   210,   126,   238,   210,    -1,
     137,    87,    34,   238,   210,    -1,   137,    87,   154,    87,
      34,   238,   210,    -1,    23,    87,   154,    87,    34,   238,
     210,    -1,   270,   154,    54,    -1,   270,    54,    -1,   270,
     154,    54,   214,    54,    -1,   270,    54,   214,    54,    -1,
      54,   214,    54,    -1,    54,    -1,   118,    34,   215,   273,
     216,   210,    -1,   272,   154,   238,    -1,   272,   354,    -1,
     238,    -1,   273,   210,   274,    -1,   274,    -1,   274,   154,
     211,   238,   212,    -1,   274,   154,    43,    -1,   274,   154,
      54,    -1,   274,   211,   238,   212,    -1,   274,    43,    -1,
     274,    54,    -1,   211,   238,   212,    -1,    43,    -1,    54,
      -1,   127,   210,    -1,   127,   211,   276,   212,   210,    -1,
     276,   154,   277,    -1,   277,    -1,   356,    -1,   453,    -1,
     454,    -1,    20,   210,    -1,    20,   211,   279,   212,   210,
      -1,   279,   154,   280,    -1,   280,    -1,   356,    -1,   119,
     210,    -1,   119,   211,   282,   212,   210,    -1,   282,   154,
     283,    -1,   283,    -1,   369,    -1,   376,    -1,   371,    -1,
     128,   210,    -1,   128,   211,   285,   212,   210,    -1,   128,
     287,   210,    -1,   128,   211,   285,   212,   287,   210,    -1,
     285,   154,   286,    -1,   286,    -1,   355,    -1,   356,    -1,
     357,    -1,   358,    -1,   359,    -1,   360,    -1,   361,    -1,
     362,    -1,   363,    -1,   364,    -1,   365,    -1,   383,    -1,
     366,    -1,   405,    -1,   367,    -1,   368,    -1,   369,    -1,
     370,    -1,   373,    -1,   374,    -1,   375,    -1,   409,    -1,
     410,    -1,   287,    87,    -1,   287,    87,    34,    87,    -1,
     287,   154,    87,    -1,   287,   154,    87,    34,    87,    -1,
      87,    -1,    87,    34,    87,    -1,   156,    54,    -1,   155,
      54,    -1,    54,    -1,   156,    43,    -1,   155,    43,    -1,
      43,    -1,    36,   210,   291,    32,    -1,   291,   292,    -1,
     292,    -1,   293,   154,   294,   210,    -1,   126,    87,    -1,
      87,    -1,    23,    87,   154,    87,    -1,   302,   154,   295,
      -1,   303,   154,   302,   154,   295,    -1,   303,   154,   303,
     154,   303,   154,   302,   154,   295,    -1,   303,    -1,   303,
     154,   303,   154,   303,    -1,   303,   154,   303,    -1,   303,
     154,   303,   154,   303,    -1,   303,   154,   303,   154,   303,
     154,   303,    -1,   303,   154,   303,   154,   303,   154,   303,
     154,   303,    -1,    38,   210,   297,    32,    -1,   297,   298,
      -1,   298,    -1,   126,    87,   154,   303,   210,    -1,    23,
      87,   154,    87,   154,   303,   210,    -1,    87,   154,   303,
     210,    -1,    37,   210,   300,    32,    -1,   300,   301,    -1,
     301,    -1,   126,    87,   154,   303,   154,   303,   210,    -1,
      23,    87,   154,    87,   154,   303,   154,   303,   210,    -1,
      87,   154,   303,   154,   303,   210,    -1,     6,    -1,    45,
      -1,    97,    -1,    55,    -1,   133,    -1,    -1,    54,    -1,
      43,    -1,    87,    -1,   155,    54,    -1,   155,    43,    -1,
      35,   210,    -1,    35,   211,   305,   212,   210,    -1,    35,
     287,   210,    -1,    35,   211,   305,   212,   287,   210,    -1,
     305,   154,   306,    -1,   306,    -1,   376,    -1,   377,    -1,
     378,    -1,   379,    -1,   380,    -1,   381,    -1,   382,    -1,
     383,    -1,   384,    -1,   385,    -1,   386,    -1,   387,    -1,
     388,    -1,   389,    -1,   390,    -1,   391,    -1,   392,    -1,
     393,    -1,   394,    -1,   395,    -1,   396,    -1,   397,    -1,
     398,    -1,   399,    -1,   366,    -1,   400,    -1,   401,    -1,
     402,    -1,   403,    -1,   404,    -1,   406,    -1,   407,    -1,
     411,    -1,   412,    -1,   413,    -1,   356,    -1,   414,    -1,
     415,    -1,   416,    -1,   111,   211,   308,   212,   210,    -1,
     111,   211,   308,   212,   287,   210,    -1,   308,   154,   309,
      -1,   309,    -1,   383,    -1,   384,    -1,   393,    -1,   399,
      -1,   366,    -1,   400,    -1,   401,    -1,   402,    -1,   403,
      -1,   404,    -1,   411,    -1,   412,    -1,   413,    -1,   112,
     211,   308,   212,   210,    -1,   112,   211,   308,   212,   287,
     210,    -1,   217,    87,   217,   154,   217,    87,   217,    -1,
     217,    87,   217,   154,   303,    -1,   311,    -1,   312,   154,
     311,    -1,   140,   287,   210,    -1,    98,   210,   315,    32,
      -1,   315,   316,    -1,   316,    -1,    87,   211,   238,   212,
     210,    -1,   134,   287,   210,    -1,   100,   210,   319,    32,
      -1,   319,    87,   238,   210,    -1,   319,    87,   154,    87,
     238,   210,    -1,    87,   238,   210,    -1,    87,   154,    87,
     238,   210,    -1,   103,   287,   210,    -1,   102,   210,    -1,
     102,   211,   286,   212,   210,    -1,   102,   287,   210,    -1,
     102,   211,   286,   212,   287,   210,    -1,    19,   210,   323,
      32,    -1,   323,   324,    -1,   324,    -1,    87,   325,    34,
     238,   210,    -1,    87,   154,    87,   325,    34,   238,   210,
      -1,     4,    87,   211,    54,   212,   325,    34,   238,   210,
      -1,    -1,   211,    54,   212,    -1,   211,    43,   212,    -1,
      18,   210,    -1,    18,   211,    24,   212,   210,    -1,    31,
     211,    87,   212,   210,    -1,    31,   211,    87,   212,   287,
     210,    -1,    31,    87,   210,    -1,    31,   211,    87,   218,
      87,   212,   210,    -1,    31,   211,    87,   218,    87,   212,
     287,   210,    -1,    31,    87,   218,    87,   210,    -1,    30,
     211,    87,   212,   210,    -1,    30,   211,    87,   212,   287,
     210,    -1,    30,    87,   210,    -1,    30,   211,    87,   218,
      87,   212,   210,    -1,    30,   211,    87,   218,    87,   212,
     287,   210,    -1,    30,    87,   218,    87,   210,    -1,    82,
     211,   330,   212,   332,   210,    -1,   330,   154,   331,    -1,
     331,    -1,   408,    -1,   409,    -1,   410,    -1,   333,    -1,
     332,   154,   333,    -1,   333,   211,   303,   212,    -1,   332,
     154,   333,   211,   303,   212,    -1,   334,    -1,   333,   334,
      -1,    87,    -1,   219,    -1,   157,    -1,   214,    -1,   218,
      -1,    -1,    -1,   106,   336,   263,   337,   210,    -1,   130,
     210,    -1,   130,   211,   339,   212,   210,    -1,   130,   287,
     210,    -1,   130,   211,   339,   212,   287,   210,    -1,   339,
     154,   340,    -1,   340,    -1,   286,    -1,   417,    -1,   418,
      -1,   419,    -1,   420,    -1,   421,    -1,   422,    -1,   423,
      -1,   424,    -1,   341,    -1,   376,    -1,   411,    -1,   412,
      -1,   378,    -1,   380,    -1,   377,    -1,   379,    -1,   414,
      -1,   415,    -1,   342,   154,   343,    -1,   342,    -1,     8,
      54,   210,    -1,     8,   211,   343,   212,    54,   210,    -1,
     342,    -1,   401,    -1,   384,    -1,   425,    -1,   345,   154,
     346,    -1,   345,    -1,     9,    54,   210,    -1,     9,   211,
     346,   212,    54,   210,    -1,   178,   210,    -1,   178,   211,
     349,   212,   210,    -1,   350,   154,   349,    -1,   350,    -1,
     426,    -1,   427,    -1,   428,    -1,   429,    -1,   430,    -1,
     431,    -1,   432,    -1,   433,    -1,   434,    -1,   435,    -1,
     436,    -1,   437,    -1,   438,    -1,   439,    -1,   440,    -1,
     441,    -1,   442,    -1,   443,    -1,   445,    -1,   446,    -1,
     447,    -1,   448,    -1,   449,    -1,   450,    -1,   451,    -1,
     452,    -1,   444,    -1,   144,   210,   352,    32,    -1,   353,
      -1,   352,   353,    -1,    87,   154,   238,   154,   238,   210,
      -1,    54,    -1,    43,    -1,    27,    34,    54,    -1,   123,
      34,    54,    -1,   120,    34,    54,    -1,    63,    -1,   101,
      34,    54,    -1,   115,    34,    54,    -1,    28,    34,    54,
      -1,     3,    34,    54,    -1,    91,    -1,    93,    -1,    95,
      -1,    56,    34,    54,    -1,    51,    34,    54,    -1,    52,
      34,    54,    -1,   105,    34,    54,    -1,    25,    34,   354,
      -1,    69,    34,    54,    -1,    69,    34,    66,    -1,    69,
      34,    46,    -1,    69,    34,    48,    -1,    69,    34,     7,
      -1,    67,    34,   354,    -1,   119,    -1,   121,    34,    54,
      -1,   113,    34,   354,    -1,    26,    34,    87,    -1,    89,
      34,   458,    -1,    89,    34,    54,    -1,    42,    34,    54,
      -1,   107,    34,    54,    -1,   108,    34,    54,    -1,    61,
      34,    54,    -1,    62,    34,    54,    -1,    94,    -1,    49,
      -1,    21,    34,   354,    -1,    75,    34,    54,    -1,    70,
      34,   354,    -1,    72,    34,   354,    -1,    99,    34,   211,
     312,   212,    -1,    71,    34,   354,    -1,    80,    34,    87,
      -1,    79,    34,    54,    -1,    78,    -1,   110,    34,   354,
      -1,    73,    34,    54,    -1,    74,    34,    54,    -1,    64,
      -1,    65,    -1,    92,    -1,     5,    -1,   129,    -1,    44,
      34,    54,    -1,   122,    -1,    86,    -1,    41,    -1,   114,
      -1,    57,    34,    54,    -1,    58,    34,    54,    -1,    84,
      34,    59,    -1,    84,    34,    85,    -1,   109,    -1,    96,
      -1,   141,    34,    87,    -1,   142,    34,   455,    -1,    40,
      34,   458,    -1,    22,    -1,    90,    -1,    76,    -1,   131,
      34,   354,    -1,    15,    34,   289,    -1,    10,    34,   354,
      -1,    12,    34,   289,    -1,    13,    34,   354,    -1,    14,
      34,    54,    -1,    11,    -1,    16,    34,    54,    -1,    17,
      34,    54,    -1,   179,    34,    54,    -1,   180,    34,    54,
      -1,   181,    34,    54,    -1,   182,    34,    54,    -1,   183,
      34,    54,    -1,   184,    34,    54,    -1,   185,    34,    54,
      -1,   186,    34,    54,    -1,   187,    34,    54,    -1,   188,
      34,    54,    -1,   189,    34,    54,    -1,   190,    34,    54,
      -1,   191,    34,    54,    -1,   192,    34,    54,    -1,   193,
      34,    54,    -1,   194,    34,   354,    -1,   195,    34,   354,
      -1,   196,    34,    54,    -1,   197,    34,   458,    -1,   198,
      34,   354,    -1,   199,    34,   354,    -1,   203,    34,    54,
      -1,   204,    34,    54,    -1,   206,    34,   354,    -1,   207,
      34,    54,    -1,   208,    34,   354,    -1,   209,    34,   354,
      -1,   145,    34,    54,    -1,   146,    34,    54,    -1,    87,
     214,    87,    -1,    54,    -1,    54,   214,    54,    -1,   215,
     456,    -1,   457,   456,    -1,   457,   216,    -1
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
      80,    82,    84,    86,    88,    90,    92,    94,    96,    98,
     102,   107,   111,   115,   119,   123,   127,   130,   134,   136,
     140,   145,   148,   151,   155,   157,   161,   166,   169,   172,
     176,   178,   182,   187,   190,   193,   197,   199,   203,   208,
     211,   215,   220,   224,   229,   233,   238,   243,   247,   249,
     251,   253,   257,   261,   265,   269,   273,   277,   281,   285,
     289,   293,   297,   300,   303,   308,   313,   318,   323,   328,
     333,   338,   343,   348,   353,   358,   365,   372,   377,   386,
     388,   392,   397,   405,   409,   414,   417,   419,   424,   429,
     432,   434,   442,   446,   448,   452,   454,   456,   458,   460,
     462,   464,   466,   468,   469,   475,   476,   485,   486,   495,
     496,   507,   508,   517,   518,   529,   530,   539,   542,   545,
     547,   549,   554,   557,   561,   563,   565,   567,   571,   575,
     579,   583,   587,   591,   595,   599,   603,   607,   611,   614,
     617,   622,   627,   632,   637,   642,   647,   652,   657,   662,
     667,   672,   679,   686,   695,   701,   703,   708,   713,   718,
     721,   723,   733,   740,   746,   754,   762,   766,   769,   775,
     780,   784,   786,   793,   797,   800,   802,   806,   808,   814,
     818,   822,   827,   830,   833,   837,   839,   841,   844,   850,
     854,   856,   858,   860,   862,   865,   871,   875,   877,   879,
     882,   888,   892,   894,   896,   898,   900,   903,   909,   913,
     920,   924,   926,   928,   930,   932,   934,   936,   938,   940,
     942,   944,   946,   948,   950,   952,   954,   956,   958,   960,
     962,   964,   966,   968,   970,   972,   975,   980,   984,   990,
     992,   996,   999,  1002,  1004,  1007,  1010,  1012,  1017,  1020,
    1022,  1027,  1030,  1032,  1037,  1041,  1047,  1057,  1059,  1065,
    1069,  1075,  1083,  1093,  1098,  1101,  1103,  1109,  1117,  1122,
    1127,  1130,  1132,  1140,  1150,  1157,  1159,  1161,  1163,  1165,
    1167,  1168,  1170,  1172,  1174,  1177,  1180,  1183,  1189,  1193,
    1200,  1204,  1206,  1208,  1210,  1212,  1214,  1216,  1218,  1220,
    1222,  1224,  1226,  1228,  1230,  1232,  1234,  1236,  1238,  1240,
    1242,  1244,  1246,  1248,  1250,  1252,  1254,  1256,  1258,  1260,
    1262,  1264,  1266,  1268,  1270,  1272,  1274,  1276,  1278,  1280,
    1282,  1284,  1290,  1297,  1301,  1303,  1305,  1307,  1309,  1311,
    1313,  1315,  1317,  1319,  1321,  1323,  1325,  1327,  1329,  1335,
    1342,  1350,  1356,  1358,  1362,  1366,  1371,  1374,  1376,  1382,
    1386,  1391,  1396,  1403,  1407,  1413,  1417,  1420,  1426,  1430,
    1437,  1442,  1445,  1447,  1453,  1461,  1471,  1472,  1476,  1480,
    1483,  1489,  1495,  1502,  1506,  1514,  1523,  1529,  1535,  1542,
    1546,  1554,  1563,  1569,  1576,  1580,  1582,  1584,  1586,  1588,
    1590,  1594,  1599,  1606,  1608,  1611,  1613,  1615,  1617,  1619,
    1621,  1622,  1623,  1629,  1632,  1638,  1642,  1649,  1653,  1655,
    1657,  1659,  1661,  1663,  1665,  1667,  1669,  1671,  1673,  1675,
    1677,  1679,  1681,  1683,  1685,  1687,  1689,  1691,  1693,  1697,
    1699,  1703,  1710,  1712,  1714,  1716,  1718,  1722,  1724,  1728,
    1735,  1738,  1744,  1748,  1750,  1752,  1754,  1756,  1758,  1760,
    1762,  1764,  1766,  1768,  1770,  1772,  1774,  1776,  1778,  1780,
    1782,  1784,  1786,  1788,  1790,  1792,  1794,  1796,  1798,  1800,
    1802,  1804,  1809,  1811,  1814,  1821,  1823,  1825,  1829,  1833,
    1837,  1839,  1843,  1847,  1851,  1855,  1857,  1859,  1861,  1865,
    1869,  1873,  1877,  1881,  1885,  1889,  1893,  1897,  1901,  1905,
    1907,  1911,  1915,  1919,  1923,  1927,  1931,  1935,  1939,  1943,
    1947,  1949,  1951,  1955,  1959,  1963,  1967,  1973,  1977,  1981,
    1985,  1987,  1991,  1995,  1999,  2001,  2003,  2005,  2007,  2009,
    2013,  2015,  2017,  2019,  2021,  2025,  2029,  2033,  2037,  2039,
    2041,  2045,  2049,  2053,  2055,  2057,  2059,  2063,  2067,  2071,
    2075,  2079,  2083,  2085,  2089,  2093,  2097,  2101,  2105,  2109,
    2113,  2117,  2121,  2125,  2129,  2133,  2137,  2141,  2145,  2149,
    2153,  2157,  2161,  2165,  2169,  2173,  2177,  2181,  2185,  2189,
    2193,  2197,  2201,  2205,  2209,  2213,  2215,  2219,  2222,  2225
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,   118,   118,   119,   122,   123,   124,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   160,   161,   162,   165,   166,   167,   168,   172,
     174,   178,   180,   182,   184,   186,   188,   190,   192,   194,
     196,   198,   202,   204,   206,   208,   210,   212,   216,   218,
     220,   222,   224,   226,   230,   232,   234,   236,   238,   240,
     244,   246,   250,   252,   256,   258,   263,   265,   267,   269,
     271,   273,   275,   277,   279,   281,   283,   285,   287,   289,
     291,   293,   295,   297,   299,   301,   303,   305,   307,   309,
     311,   313,   315,   317,   319,   321,   323,   325,   327,   331,
     333,   337,   339,   343,   345,   347,   348,   351,   353,   355,
     356,   359,   361,   362,   365,   366,   369,   370,   373,   375,
     377,   381,   382,   385,   385,   387,   387,   389,   389,   392,
     391,   394,   394,   397,   396,   399,   399,   403,   404,   405,
     406,   409,   411,   415,   417,   418,   420,   422,   424,   426,
     428,   430,   432,   434,   436,   438,   440,   442,   444,   446,
     448,   450,   452,   454,   456,   458,   460,   462,   464,   466,
     468,   470,   472,   474,   478,   481,   483,   487,   489,   491,
     492,   495,   497,   499,   501,   503,   507,   509,   511,   513,
     515,   517,   522,   525,   527,   529,   533,   535,   539,   541,
     543,   545,   547,   549,   551,   553,   555,   559,   561,   565,
     566,   569,   570,   571,   574,   576,   580,   581,   584,   586,
     588,   592,   593,   596,   597,   598,   601,   603,   605,   607,
     611,   612,   615,   616,   617,   618,   619,   620,   621,   622,
     623,   624,   625,   626,   627,   628,   629,   630,   631,   632,
     633,   634,   635,   636,   637,   640,   642,   644,   646,   648,
     650,   654,   656,   658,   662,   664,   666,   670,   672,   674,
     678,   680,   686,   692,   702,   707,   714,   725,   730,   741,
     748,   757,   768,   783,   786,   788,   792,   800,   810,   820,
     823,   825,   829,   839,   851,   863,   865,   867,   869,   871,
     875,   876,   877,   878,   879,   881,   885,   887,   889,   891,
     895,   896,   899,   900,   901,   902,   903,   904,   905,   906,
     907,   908,   909,   910,   911,   912,   913,   914,   915,   916,
     917,   918,   919,   920,   921,   922,   923,   924,   925,   926,
     927,   928,   929,   930,   931,   932,   933,   934,   935,   936,
     937,   940,   942,   946,   947,   950,   951,   952,   953,   954,
     955,   956,   957,   958,   959,   960,   961,   962,   965,   967,
     971,   973,   977,   978,   981,   983,   985,   986,   989,   991,
     993,   995,   997,   999,  1001,  1005,  1007,  1009,  1011,  1013,
    1017,  1019,  1020,  1023,  1025,  1027,  1031,  1032,  1034,  1038,
    1040,  1044,  1046,  1048,  1050,  1052,  1054,  1058,  1060,  1062,
    1064,  1066,  1068,  1072,  1075,  1076,  1079,  1080,  1081,  1084,
    1086,  1088,  1090,  1094,  1096,  1100,  1101,  1103,  1105,  1107,
    1111,  1112,  1111,  1114,  1116,  1118,  1120,  1124,  1125,  1128,
    1129,  1132,  1133,  1134,  1135,  1136,  1137,  1138,  1141,  1142,
    1143,  1144,  1145,  1146,  1147,  1148,  1149,  1150,  1153,  1154,
    1157,  1159,  1163,  1164,  1165,  1166,  1169,  1170,  1173,  1175,
    1179,  1181,  1185,  1186,  1189,  1190,  1191,  1192,  1193,  1194,
    1195,  1196,  1197,  1198,  1199,  1200,  1201,  1202,  1203,  1204,
    1205,  1206,  1207,  1208,  1209,  1210,  1211,  1212,  1213,  1214,
    1215,  1218,  1221,  1222,  1225,  1228,  1229,  1232,  1233,  1234,
    1235,  1236,  1237,  1238,  1239,  1240,  1241,  1242,  1243,  1244,
    1245,  1246,  1248,  1249,  1250,  1251,  1252,  1253,  1254,  1255,
    1256,  1257,  1258,  1259,  1261,  1264,  1265,  1266,  1267,  1268,
    1269,  1271,  1274,  1275,  1276,  1277,  1278,  1279,  1280,  1281,
    1282,  1283,  1284,  1285,  1286,  1287,  1288,  1289,  1290,  1291,
    1292,  1293,  1294,  1295,  1296,  1297,  1298,  1300,  1303,  1304,
    1305,  1306,  1307,  1308,  1309,  1310,  1311,  1313,  1314,  1315,
    1316,  1317,  1318,  1319,  1320,  1322,  1323,  1324,  1325,  1326,
    1327,  1328,  1329,  1330,  1331,  1332,  1333,  1334,  1335,  1336,
    1337,  1338,  1339,  1340,  1342,  1343,  1349,  1350,  1354,  1355,
    1356,  1357,  1359,  1360,  1362,  1370,  1371,  1380,  1382,  1391
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
       2,     2,     2,     2,     2,   213,     2,     2,     2,   217,
     211,   212,     2,     2,     2,     2,   218,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   214,   210,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   215,   219,   216,     2,     2,     2,     2,     2,     2,
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
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 2589;
  const int parser::yynnts_ = 239;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 163;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 220;

  const unsigned int parser::yyuser_token_number_max_ = 464;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1393 "DynareBison.yy"


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

