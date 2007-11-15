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
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 114:
#line 301 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 115:
#line 303 "DynareBison.yy"
    { (yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val)); ;}
    break;

  case 116:
#line 305 "DynareBison.yy"
    { (yyval.node_val) = driver.add_normcdf((yysemantic_stack_[(8) - (3)].node_val),(yysemantic_stack_[(8) - (5)].node_val),(yysemantic_stack_[(8) - (7)].node_val));;}
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
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 181:
#line 433 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 182:
#line 435 "DynareBison.yy"
    { (yyval.node_val) = driver.add_normcdf((yysemantic_stack_[(8) - (3)].node_val),(yysemantic_stack_[(8) - (5)].node_val),(yysemantic_stack_[(8) - (7)].node_val));;}
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
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 196:
#line 470 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 197:
#line 472 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 198:
#line 474 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 199:
#line 478 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 200:
#line 482 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 201:
#line 484 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].node_val));;}
    break;

  case 202:
#line 488 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 203:
#line 490 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 204:
#line 494 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 205:
#line 496 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 206:
#line 498 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 207:
#line 500 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 208:
#line 502 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 209:
#line 504 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 210:
#line 506 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 211:
#line 508 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 212:
#line 510 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 213:
#line 514 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 214:
#line 516 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 218:
#line 526 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 219:
#line 528 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 223:
#line 538 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 224:
#line 540 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 229:
#line 552 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 230:
#line 554 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 231:
#line 556 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 232:
#line 558 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 258:
#line 591 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 259:
#line 593 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 260:
#line 595 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 261:
#line 597 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 262:
#line 599 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 263:
#line 601 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 264:
#line 605 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 265:
#line 607 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 266:
#line 609 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 267:
#line 613 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 268:
#line 615 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 269:
#line 617 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 270:
#line 620 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 271:
#line 623 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 272:
#line 625 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 274:
#line 631 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 275:
#line 637 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 276:
#line 643 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 277:
#line 653 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 278:
#line 658 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 279:
#line 665 "DynareBison.yy"
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

  case 280:
#line 676 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 281:
#line 681 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 282:
#line 692 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 283:
#line 699 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 284:
#line 708 "DynareBison.yy"
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

  case 285:
#line 719 "DynareBison.yy"
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

  case 286:
#line 734 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 287:
#line 737 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 288:
#line 739 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 289:
#line 743 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 290:
#line 751 "DynareBison.yy"
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

  case 291:
#line 761 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 292:
#line 771 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 293:
#line 774 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 294:
#line 776 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 295:
#line 780 "DynareBison.yy"
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

  case 296:
#line 790 "DynareBison.yy"
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

  case 297:
#line 802 "DynareBison.yy"
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

  case 298:
#line 814 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 299:
#line 816 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 300:
#line 818 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 301:
#line 820 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 302:
#line 822 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 303:
#line 825 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 307:
#line 830 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 308:
#line 832 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 309:
#line 836 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 310:
#line 838 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 311:
#line 840 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 312:
#line 842 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 354:
#line 891 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 355:
#line 893 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 371:
#line 916 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 372:
#line 918 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 373:
#line 922 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 374:
#line 924 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 377:
#line 931 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 378:
#line 933 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 381:
#line 939 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 382:
#line 941 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 383:
#line 943 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 384:
#line 946 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 385:
#line 948 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 386:
#line 950 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 387:
#line 952 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 388:
#line 955 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 389:
#line 958 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 390:
#line 960 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 391:
#line 962 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 392:
#line 964 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 393:
#line 967 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 396:
#line 974 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 397:
#line 976 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 398:
#line 978 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 399:
#line 981 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 400:
#line 983 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 401:
#line 985 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 402:
#line 989 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 403:
#line 991 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 404:
#line 995 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 405:
#line 997 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 406:
#line 999 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 407:
#line 1001 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 408:
#line 1003 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 409:
#line 1005 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 410:
#line 1009 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 411:
#line 1011 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 412:
#line 1013 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 413:
#line 1015 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 414:
#line 1017 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 415:
#line 1019 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 416:
#line 1023 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 422:
#line 1035 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 423:
#line 1037 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 424:
#line 1039 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 425:
#line 1041 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 426:
#line 1045 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 427:
#line 1047 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 429:
#line 1052 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 430:
#line 1054 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 431:
#line 1056 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 432:
#line 1058 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 433:
#line 1061 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 434:
#line 1062 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 436:
#line 1065 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 437:
#line 1067 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 438:
#line 1069 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 439:
#line 1071 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 463:
#line 1108 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 464:
#line 1110 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 471:
#line 1124 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 472:
#line 1126 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 473:
#line 1130 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 474:
#line 1132 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 506:
#line 1172 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 507:
#line 1173 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 508:
#line 1174 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 509:
#line 1175 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 510:
#line 1176 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 511:
#line 1177 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 512:
#line 1178 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 513:
#line 1179 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 514:
#line 1180 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 515:
#line 1181 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 516:
#line 1182 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 517:
#line 1183 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 518:
#line 1184 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 519:
#line 1185 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 520:
#line 1187 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 521:
#line 1188 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 522:
#line 1189 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 523:
#line 1190 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 524:
#line 1191 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 525:
#line 1192 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 526:
#line 1193 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 527:
#line 1195 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 528:
#line 1197 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1199 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1200 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1201 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1202 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1203 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 534:
#line 1205 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 535:
#line 1207 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 536:
#line 1209 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 537:
#line 1210 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 538:
#line 1211 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1212 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1214 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 542:
#line 1215 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 543:
#line 1216 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 544:
#line 1217 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 545:
#line 1218 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 546:
#line 1219 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 547:
#line 1220 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 548:
#line 1221 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 549:
#line 1222 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 550:
#line 1223 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 551:
#line 1224 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 552:
#line 1225 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 553:
#line 1226 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 554:
#line 1227 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 555:
#line 1228 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 556:
#line 1229 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 557:
#line 1230 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 558:
#line 1231 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 559:
#line 1232 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 560:
#line 1234 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 561:
#line 1236 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 562:
#line 1238 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 563:
#line 1239 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 564:
#line 1240 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1241 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 566:
#line 1242 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1243 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 568:
#line 1244 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 569:
#line 1245 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 570:
#line 1246 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 571:
#line 1248 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1249 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1250 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1251 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 575:
#line 1252 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 576:
#line 1253 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 577:
#line 1254 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1255 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 579:
#line 1257 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1258 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 581:
#line 1259 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 582:
#line 1260 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1261 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1262 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1263 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1264 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1265 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1266 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1267 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1268 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1269 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 592:
#line 1270 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 593:
#line 1271 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 594:
#line 1272 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 595:
#line 1273 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 596:
#line 1274 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 597:
#line 1275 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 598:
#line 1277 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 599:
#line 1278 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 600:
#line 1284 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 601:
#line 1285 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 602:
#line 1289 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 603:
#line 1290 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 604:
#line 1291 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 605:
#line 1292 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 606:
#line 1296 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 608:
#line 1305 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); 
               ;}
    break;

  case 609:
#line 1314 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 610:
#line 1316 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 611:
#line 1324 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2409 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1102;
  const short int
  parser::yypact_[] =
  {
       979,    37,    47,    84,  -118,   334,    87,    71,    25,    27,
     -53,    72,    57,    74,   115,   160,   375,    91,   406,    80,
     183,   269,   228,   235,    78,   267,   360,   102, -1102,   245,
     248,   267,   263,   435,   430,   441,    86,    88,   267,   391,
     405,   409,   267,   446,   844, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
     295,   640,   305,  1341, -1102,   513,   225, -1102,   462,   551,
     419,    51,   154,   544,   333,   556,   566,   582, -1102,  1242,
      12,   358,   368,   403,   574,   566,   619,   621,   464, -1102,
     341,   533,   129,  1074,   585,   590, -1102,  1445,    52,   210,
     548,   227,   626,   479,  1098,   655,   655,   229,   129,   475,
   -1102,    93, -1102,   462, -1102,  1445,   234, -1102,  1345,   264,
     293,   562,   299,   564,   300,   565,   302,   306, -1102,  2094,
   -1102, -1102, -1102,   649, -1102,   659,   660,   664,   666,   667,
   -1102,   671,   673,   676, -1102,   677,   678,   679,   680, -1102,
     571,   494, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,   684,
     686,   687, -1102,   576,   520, -1102, -1102, -1102,   524,   641,
     -46,   109, -1102,   695,   -37, -1102, -1102,   529, -1102,   530,
   -1102, -1102,   650,   415, -1102,   651,   420,   701,   279, -1102,
     656, -1102,   704, -1102, -1102,   706,   722,   724,   726,   733,
   -1102, -1102,   741,   744,   746,   748,   749,   750, -1102, -1102,
     751,   752, -1102, -1102, -1102,   754,   755, -1102, -1102,   -29,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
     756,   710, -1102,   711, -1102,   712,   408, -1102,   653,   714,
     654,   719,   450, -1102,   720,   658,   723,   587, -1102,   605,
     354, -1102,   400,   774,   606,   612, -1102,   663, -1102,   149,
     611,   623,   781, -1102, -1102,   150, -1102, -1102, -1102, -1102,
     742,   745,   172, -1102, -1102, -1102,   625,   627,   630,   631,
    1074,  1074,   632,   633,   634,   635,   636,   637,   638,   639,
     645,   646,  1074,  1816,   647,   402, -1102,   977,   499,   809,
     816,   817,   820,   821,   822, -1102, -1102, -1102,   823,   824,
     825, -1102,   826, -1102,   827,   832,   665, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102,   743,   788, -1102,   675, -1102, -1102, -1102,
     670,   681,   682,   683,  1098,  1098,   685,   689,   691,   696,
     697,   698,   699,   700,   702,   703,  1098,  2105, -1102,   246,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102,   259, -1102,   328,    33,   287, -1102,
   -1102, -1102,   292, -1102, -1102,   309, -1102, -1102,   852, -1102,
     311, -1102, -1102, -1102, -1102, -1102,   761,   795, -1102, -1102,
     762,   814, -1102, -1102,   776,   828, -1102, -1102,   872,   873,
     875,   876,   879,   880,   881,   882,   883,   884,   885,   886,
     890,   891,   892,   894,   895,   896,   897,   898,   899,   901,
     902,   903,   905,   907,   913,   747,   802, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102,    99,   167,    99,   908,   167,   909,
     866,   910,    36,   911,   912,   870,   887,   640,   916,   919,
      99,   920,  1341,   923,   764,   753,   900,   425,   940, -1102,
   -1102,   928,   462,   780, -1102, -1102,   783,    55,   906,   784,
      63,   917,  1074, -1102, -1102, -1102,   787,   930,   934,   938,
     942,   943,    99,    99,    99,   944,   949,   950,   951,   924,
     803,    99,  1242,    66,   935,   972,   865, -1102, -1102, -1102,
      79,   877,    61,   888, -1102, -1102,   889,    61,   893, -1102,
   -1102,   153, -1102, -1102, -1102,   936,   830, -1102,   939,    29,
   -1102,   275, -1102,    81, -1102,   831,   835,    69,   533,    50,
     937,    45, -1102, -1102,  1074,  1074,  1074,  1074,   869,   416,
    1074,  1074,  1074,  1074,  1074,  1074,  1074,  1074,  1074,  1074,
    1017,  1074,  1074,  1074,  1074,  1074,  1074,  1074,  1074,  1074,
    1074,  1074, -1102,  1074, -1102, -1102,   945,  1828, -1102,  1045,
     973,    99,   981,   985,   986,   988,   989,   992,    99,   993,
     994,   995,    70, -1102,   921, -1102,  1098,  1098,   153,  1074,
     932,   522,  1098,  1098,  1098,  1098,  1098,  1098,  1098,  1098,
    1098,  1098,  1266,  1098,  1098,  1098,  1098,  1098,  1098,  1098,
    1098,  1098,  1098,  1098,   848,   655,    76,    90, -1102, -1102,
   -1102,  1074,   434,    32,    93,   850,   462,   851,  1445,    96,
      99,  1345,    97, -1102,   926, -1102,   927, -1102,   955,  1003,
    1007,  1024,  1033,  1037,  1038,  1043,  1044,  1046,  1047,  1048,
    1051,  1053,  1054,  1055,    99,    99,  1057,   787,    99,    99,
    1064,  1067,    99,  1068,    99,    99,   941,  2094, -1102, -1102,
   -1102, -1102,  1078,  1079, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102,  1075,    21, -1102, -1102, -1102, -1102,   946, -1102,
   -1102,   948, -1102, -1102, -1102, -1102,   952, -1102,  1077,   953,
     965,   967,  1074, -1102, -1102, -1102, -1102, -1102,   307,   969,
   -1102, -1102,   320,   970,  1849, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,   966,
   -1102, -1102, -1102,   321, -1102,  1060,  1061, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102,   521,   974,  1000,  1006,  1071,
    1030,    61,  1094,   983,    61, -1102,  1088,  1126,   982, -1102,
     566,  1146, -1102, -1102, -1102,  1098, -1102, -1102, -1102,  1148,
   -1102,   312, -1102, -1102, -1102,   987, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102,   -41,    -9, -1102,  1103,
    1074,  1104,    -3,   922,  2090,  2228,   313,  2156,  1366,  1378,
    1435,  1447,  1459,  1471,  1483,  1504,  1516,  1528, -1102,   538,
     538,   538,   538,   538,   538,   416,   416,   869,   869,  1042,
    1540,  1074, -1102,  1106,  1861, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,   324, -1102,
    2168,  2180,  1004,  2192,  1552,  1573,  1585,  1597,  1609,  1621,
    1642,  1654,  1666,  1678, -1102,   622,   622,   622,   622,   622,
     622,   522,   522,   932,   932,  1042, -1102, -1102, -1102,   326,
   -1102,   329,  1690,    33,   996, -1102, -1102,    35,  1074, -1102,
   -1102, -1102, -1102, -1102, -1102,   338, -1102, -1102, -1102,   371,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102,   990, -1102, -1102, -1102,  1107, -1102,
   -1102,  1005,  1162, -1102, -1102,  1883, -1102,    98, -1102,   123,
   -1102,  1127, -1102,   317, -1102, -1102, -1102, -1102, -1102, -1102,
      61,    79,  1063,    61,  1065,  1069, -1102,  1013, -1102, -1102,
    1182,   500,  1098,  1895,    99,    81, -1102,   663,   663,   663,
      50, -1102,    61, -1102,  1183,  1916,  1184,  1167,  1074,  1074,
    1074,  1074, -1102,  1074, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102,  1022,  1928,  1074, -1102, -1102,  1098,
    1098, -1102,  1074, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102,    32, -1102, -1102, -1102,
    1074,  1711, -1102, -1102,  1172, -1102,   953,  1074, -1102, -1102,
     373, -1102,   414,  1018,   966, -1102, -1102,  1091,  1092,  1093,
      61,  1039,    61,    61, -1102,  1074, -1102,  1950, -1102, -1102,
   -1102,  1040,    64,   211,   240,    60,  1041,  1074, -1102,  1074,
    1020,   -31,  1962,  1723,  1735,  2228,  2204, -1102, -1102,  1983,
    1747,  1759,  2216,  1780, -1102, -1102,  1209,  1995, -1102, -1102,
    1114, -1102,    61,    61,    61,  1115, -1102,  1049,  1066,  2017,
   -1102,   663, -1102, -1102, -1102,    61, -1102,  2029,  2050,  1197,
    1210,  1134, -1102, -1102, -1102,  1074, -1102, -1102, -1102,  1074,
   -1102,  1074, -1102,    28,  1123, -1102,  1124,    61, -1102, -1102,
   -1102,   594,  1070, -1102, -1102, -1102,  1072,  1074,  1792,  1804,
    2062,  1188, -1102,    61,    75,  1073, -1102, -1102,  1220,  2228,
     156, -1102, -1102, -1102,  1080,  1129,  1132, -1102, -1102,  1074,
   -1102, -1102,    61,    61,  2228,  1133, -1102,    61, -1102
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   433,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     4,    29,    30,    45,
      46,    47,    44,     5,     6,     7,    12,     9,    10,    11,
       8,    13,    14,    15,    16,    17,    18,    19,    23,    25,
      24,    20,    21,    22,    26,    27,    28,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,   402,     0,     0,   218,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   262,   309,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   137,
       0,     0,     0,     0,     0,     0,   389,     0,     0,     0,
      75,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     223,     0,   213,     0,   229,     0,     0,   436,     0,     0,
       0,    57,     0,    63,     0,    69,     0,     0,   473,     0,
       1,     3,   463,     0,   576,     0,     0,     0,     0,     0,
     567,     0,     0,     0,   568,     0,     0,     0,     0,   451,
     462,     0,   452,   457,   455,   458,   456,   453,   454,   459,
     460,   444,   445,   446,   447,   448,   449,   450,   471,     0,
       0,     0,   465,   470,     0,   467,   466,   468,     0,     0,
     399,     0,   395,     0,     0,   221,   222,     0,    81,     0,
      48,   412,     0,     0,   406,     0,     0,     0,     0,   124,
       0,   551,     0,   556,   535,     0,     0,     0,     0,     0,
     548,   549,     0,     0,     0,     0,     0,     0,   569,   544,
       0,     0,   555,   550,   534,     0,     0,   554,   552,     0,
     314,   350,   339,   315,   316,   317,   318,   319,   320,   321,
     322,   323,   324,   325,   326,   327,   328,   329,   330,   331,
     332,   333,   334,   335,   336,   337,   338,   340,   341,   342,
     343,   344,   345,   346,   347,   348,   349,   351,   352,   353,
     258,     0,   311,     0,   275,     0,     0,   272,     0,     0,
       0,     0,     0,   294,     0,     0,     0,     0,   288,     0,
       0,   128,     0,     0,     0,     0,    83,     0,   509,     0,
       0,     0,     0,   563,   562,     0,   418,   419,   420,   421,
       0,     0,     0,   189,    88,    89,     0,     0,    87,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   380,     0,     0,     0,
       0,     0,     0,     0,     0,   514,   515,   516,     0,     0,
       0,   557,     0,   523,     0,     0,     0,   235,   236,   237,
     238,   239,   240,   241,   242,   243,   244,   245,   247,   249,
     250,   251,   252,   253,   254,   255,   246,   248,   256,   257,
     391,   388,    78,    73,     0,    54,     0,    79,   155,   156,
       0,     0,   184,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   434,   154,     0,
     357,   362,   358,   359,   360,   361,   363,   364,   365,   366,
     367,   368,   369,   370,     0,    50,     0,     0,     0,   226,
     227,   228,     0,   216,   217,     0,   234,   231,     0,   442,
       0,   441,   443,   438,   382,    60,    55,     0,    51,    66,
      61,     0,    52,    72,    67,     0,    53,   377,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   476,   477,   478,   479,
     480,   481,   482,   483,   484,   485,   486,   487,   488,   489,
     490,   491,   492,   493,   494,   503,   495,   496,   497,   498,
     499,   500,   501,   502,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   393,
     394,     0,     0,     0,    82,    49,     0,     0,     0,     0,
       0,     0,     0,   122,   123,   263,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   260,     0,   274,   270,   271,
     303,     0,   303,     0,   292,   293,     0,   303,     0,   286,
     287,     0,   126,   127,   119,     0,     0,    84,     0,     0,
     149,     0,   150,     0,   145,     0,     0,     0,     0,     0,
       0,     0,   187,   188,     0,     0,     0,     0,   101,   102,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    85,     0,   378,   379,     0,     0,   383,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    76,    74,    80,     0,     0,     0,     0,
     168,   169,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   186,   211,
     212,     0,     0,   203,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    58,    56,    64,    62,    70,    68,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   505,   504,
     572,   269,     0,     0,   573,   574,   575,   571,   577,   526,
     529,   528,     0,     0,   527,   530,   531,   564,     0,   565,
     461,     0,   578,   536,   553,   469,     0,   403,     0,   399,
       0,     0,     0,   507,   220,   219,   415,   410,     0,     0,
     409,   404,     0,     0,     0,   566,   517,   558,   559,   532,
     533,   538,   541,   539,   546,   547,   537,   543,   542,     0,
     545,   313,   310,     0,   259,     0,     0,   298,   305,   299,
     304,   301,   306,   300,   302,     0,     0,     0,   280,     0,
       0,   303,     0,     0,   303,   266,     0,     0,     0,   121,
       0,     0,   138,   147,   148,     0,   152,   133,   132,     0,
     134,     0,   131,   135,   136,     0,   141,   139,   560,   561,
     417,   428,   430,   431,   432,   429,     0,   422,   426,     0,
       0,     0,     0,     0,     0,   117,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    86,   100,
      99,    98,    97,    96,    95,    91,    90,    92,    93,    94,
       0,     0,   386,     0,     0,   513,   521,   506,   512,   518,
     519,   510,   520,   525,   511,   508,   524,   390,     0,    77,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   153,   166,   165,   164,   163,   162,
     161,   158,   157,   159,   160,   167,   435,   356,   354,     0,
     371,     0,     0,     0,     0,   208,   209,     0,     0,   225,
     224,   215,   214,   233,   230,     0,   570,   440,   437,     0,
      59,    65,    71,   579,   580,   581,   582,   583,   584,   585,
     586,   587,   588,   589,   590,   591,   592,   593,   594,   595,
     596,   597,   598,   599,   600,   601,   602,   603,   604,   605,
     474,   475,   268,   267,   607,   609,   611,   610,     0,   464,
     472,     0,     0,   401,   400,     0,   411,     0,   405,     0,
     125,     0,   375,     0,   312,   261,   276,   308,   307,   273,
     303,   303,     0,   303,     0,     0,   291,     0,   265,   264,
       0,     0,     0,     0,     0,     0,   143,     0,     0,     0,
       0,   416,   303,   427,     0,     0,     0,     0,     0,     0,
       0,     0,   115,     0,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,     0,     0,     0,   384,   392,     0,
       0,   185,     0,   170,   171,   172,   173,   174,   175,   176,
     177,   178,   179,   355,   372,   210,   202,   199,   205,   206,
       0,     0,   232,   439,     0,   606,   399,     0,   396,   413,
       0,   407,     0,     0,     0,   540,   277,     0,     0,     0,
     303,     0,   303,   303,   289,     0,   120,     0,   151,   522,
     130,     0,     0,     0,     0,   423,     0,     0,   192,     0,
     198,     0,     0,     0,     0,   118,     0,   381,   387,     0,
       0,     0,     0,     0,   207,   608,     0,     0,   414,   408,
       0,   376,   303,   303,   303,     0,   297,     0,     0,     0,
     183,     0,   146,   142,   140,   303,   424,     0,     0,     0,
       0,     0,   191,   113,   114,     0,   385,   180,   181,     0,
     204,     0,   397,   303,   282,   278,   281,   303,   295,   290,
     129,     0,     0,   194,   193,   197,   195,     0,     0,     0,
       0,     0,   374,   303,     0,     0,   144,   425,     0,   201,
       0,   116,   182,   398,     0,   283,     0,   296,   196,     0,
     190,   373,   303,   303,   200,   284,   279,   303,   285
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1102, -1102,  1235, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,  -326, -1102,
   -1102, -1102, -1102,  -112,  -226, -1102, -1102,   960, -1102,   209,
   -1102, -1102, -1102, -1102, -1102, -1102,  -988,  -620,  -119,  -602,
   -1102, -1102, -1102,  1151,  -242, -1102, -1102, -1102, -1102,   303,
   -1102, -1102,   554, -1102, -1102,   721, -1102, -1102,   557, -1102,
   -1102,  -122,   -24,   596,   757, -1102, -1102,  1008, -1102, -1102,
   -1101, -1102, -1102,   975, -1102, -1102,  1009,  -998,  -584, -1102,
   -1102,   716, -1102,  1166,   573, -1102,   162, -1102, -1102, -1102,
   -1102,   954, -1102, -1102, -1102, -1102, -1102, -1102, -1102,  1109,
    -760, -1102, -1102, -1102, -1102, -1102,   690, -1102,   242,  -831,
   -1102, -1102, -1102, -1102, -1102,   583, -1102,   -59,   768, -1102,
   -1102,   770, -1102, -1102,   553, -1102,  -525, -1102,   -92, -1102,
    1213, -1102, -1102, -1102, -1102, -1102, -1102, -1102,  -105, -1102,
   -1102,  -121,  -583, -1102, -1102, -1102, -1102,  -101,   -80,   -77,
     -65,   -62, -1102, -1102,   -98,   -54, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102,   -57, -1102, -1102, -1102, -1102, -1102,
     -56,   -55,   -45,   -52,   -51,   -50, -1102, -1102, -1102, -1102,
    -111,   -99,   -93,   -90,   -48,   -47,   -44, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102, -1102,
   -1102, -1102, -1102, -1102, -1102,   541, -1102,  -530
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    44,    45,    46,    47,    48,    49,    50,    51,    52,
     152,   154,   156,   131,    53,    54,    55,    56,   363,   906,
      57,   324,    58,   228,   229,    59,   320,   321,   881,   882,
      60,   327,  1079,  1078,  1161,   885,   629,   630,   631,   632,
     438,    61,    62,   342,   343,  1171,    63,  1250,   732,   733,
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
       128,   129,   584,   322,   262,   386,   216,   137,   263,   873,
     338,   270,   146,   149,   150,   437,   294,   261,   157,   295,
     460,   785,   339,   466,   648,   649,   858,   874,   860,   264,
     441,   441,   265,   863,   202,   803,   660,   442,   442,   205,
     461,   677,   451,   451,   266,   452,   452,   267,   206,  1042,
     883,   464,   280,   286,   287,   271,   825,   289,   290,   291,
     872,   296,   297,  1148,   288,   298,  1083,   831,   832,   833,
     848,   418,  1034,   891,   985,   729,   840,  1128,   900,   850,
     419,   847,    96,   986,   730,   847,  1129,   791,    90,  1162,
    1163,  1164,  1225,   420,   300,  1202,   584,  1087,    92,   566,
     643,   421,   219,   848,  1080,   370,   418,   102,   572,   104,
     852,   422,   850,   209,  1210,   419,   602,  1088,   171,   849,
      99,   848,   101,   849,   117,   888,   877,   851,   420,   100,
     850,   851,   891,   118,   300,   132,   421,   107,   878,   892,
     569,   778,   891,   852,   879,   107,   422,   106,   107,   889,
     779,   340,   107,   133,   107,   567,   936,   301,   107,  1081,
     107,   852,  1266,   943,   880,   573,   423,   853,   107,  1211,
     107,   853,   107,   603,   855,   424,   425,   987,   107,   107,
     107,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     901,   210,  1082,   379,   340,   893,   435,   301,   892,   894,
     895,   423,   854,   642,   865,   107,   854,   855,   892,   781,
     424,   425,   302,  1231,   643,   996,   426,   427,   428,   429,
     430,   431,   432,   433,   434,   855,   103,  1036,   105,   209,
     436,   435,   628,   988,   731,  1241,  1130,  1021,    91,  1018,
    1019,   792,  1203,  1022,  1023,   902,  1256,  1026,    93,  1028,
    1029,   220,   410,   418,   893,   817,   824,   111,   894,   895,
     341,  1205,   419,   821,   893,   436,   842,   628,   894,   895,
     947,  1204,   108,   109,   112,   420,   978,  1064,   126,   127,
    1067,   121,   418,   421,    94,    95,   144,   145,   147,   148,
     980,   419,   300,   422,   633,   638,   994,   998,  1139,   866,
     867,  1259,   123,   341,   420,   700,   701,   210,   875,   413,
     583,   300,   421,   782,   783,   113,   300,   712,   903,   904,
     905,   907,   422,  1141,   908,   909,   910,   911,   912,   913,
     914,   915,   916,   917,  1083,   919,   920,   921,   922,   923,
     924,   925,   926,   927,   928,   929,   300,   930,   423,   107,
     340,   634,   639,   934,   221,   301,  1260,   424,   425,   728,
     114,   227,   222,   426,   427,   428,   429,   430,   431,   432,
     433,   434,   414,   953,   301,   300,  1186,   423,   435,   301,
     303,   476,   480,   122,   484,   622,   424,   425,   300,   300,
     309,   725,   426,   427,   428,   429,   430,   431,   432,   433,
     434,   328,   300,   300,   725,   982,   300,   435,   300,   301,
     411,   300,   436,   713,   628,   714,   715,   716,   717,   718,
     300,   719,   720,   721,   722,   314,   723,   415,   124,   455,
     303,   624,   734,   674,   467,   125,   319,   736,   301,   608,
     304,   436,   130,   628,   477,   481,   135,   485,   726,   136,
     310,   301,   301,   300,   738,   300,   741,  1075,  1091,   341,
     329,   727,  1144,   138,   473,   301,   301,   810,   139,   301,
     330,   301,   309,   151,   301,   876,   811,  1149,   305,  1151,
     216,   614,   227,   301,   364,   315,  1045,   153,   311,   735,
     304,   155,   883,   474,   737,   162,   300,   262,  1166,   478,
     482,   263,   486,   202,   270,   198,   487,  1046,   205,   294,
     261,   739,   295,   742,  1076,  1092,   301,   206,   301,  1145,
    1048,  1054,   264,   316,  1108,   265,  1123,   338,   305,  1124,
     678,  1156,   310,   224,    97,    98,   208,   266,  1132,   339,
     267,   225,   873,   873,   873,   280,   286,   287,   271,  1159,
     289,   290,   291,   818,   296,   297,   822,   288,   298,   301,
     874,   874,   874,  1057,   669,   670,  1195,   671,  1197,  1198,
     311,  1133,  1058,  1188,  1085,   115,   116,   950,   951,   843,
     213,   679,   227,   954,   955,   956,   957,   958,   959,   960,
     961,   962,   963,   217,   965,   966,   967,   968,   969,   970,
     971,   972,   973,   974,   975,  1105,   119,   120,  1224,   314,
    1226,   873,   332,   460,  1189,   230,   993,   577,   619,   218,
     441,  1232,   580,   578,   333,  1246,   223,   442,   581,   874,
     140,   141,   451,   461,   983,   452,   418,   334,   226,  1242,
     984,   142,   143,  1245,   464,   419,   158,   159,   227,   163,
     164,   165,   166,   167,   168,   169,   319,   323,   420,  1255,
     231,   170,  1131,   325,   326,   171,   421,   364,   948,   315,
     721,   722,   367,   723,   412,   200,   422,   416,  1265,   417,
     457,   172,   544,  1268,   667,   668,   669,   670,   475,   671,
     479,   483,   545,   546,   232,   233,   558,   547,   201,   548,
     549,   234,   979,   981,   550,   418,   551,   316,   235,   552,
     553,   554,   555,   556,   419,   995,   557,   559,   999,   560,
     561,   562,   563,   565,   173,   174,   564,   420,   571,   574,
     575,   423,   576,   579,   582,   421,   252,   586,   585,   587,
     424,   425,   175,   176,   254,   422,   426,   427,   428,   429,
     430,   431,   432,   433,   434,   588,  1073,   589,  1071,   590,
     256,   435,  1172,  1173,  1174,  1175,   591,  1176,   719,   720,
     721,   722,   257,   723,   592,   177,   178,   593,   258,   594,
    1179,   595,   596,   597,   598,   599,  1182,   600,   601,   604,
     177,   178,   605,   606,   607,   436,   611,   628,   610,   612,
     423,   613,   616,   617,  1183,   618,   621,   625,   626,   424,
     425,  1187,   627,   635,   637,   426,   427,   428,   429,   430,
     431,   432,   433,   434,   640,   636,   644,   641,   645,  1199,
     435,   646,   647,   650,   651,   652,   653,   654,   655,   656,
     657,  1207,   680,  1208,   160,   584,   658,   659,   673,   681,
     682,     1,     2,   683,   684,   685,   686,   687,   688,   689,
     690,     3,     4,     5,   436,   691,   628,   692,     6,   693,
     694,   696,     7,     8,     9,   695,    10,   744,    11,    12,
      13,    14,   697,   698,   699,   740,   702,   743,   745,  1238,
     703,    15,   704,  1239,    16,  1240,   746,   705,   706,   707,
     708,   709,   747,   710,   711,   749,   750,    17,   751,   752,
     748,  1249,   753,   754,   755,   756,   757,   758,   759,   760,
      18,    19,    20,   761,   762,   763,    21,   764,   765,   766,
     767,   768,   769,  1264,   770,   771,   772,    22,   773,    23,
     774,    24,    25,    26,    27,    28,   775,   777,   789,   776,
      29,    30,   797,  1157,   808,    31,    32,    33,    34,   786,
     788,   790,   795,   796,   807,    35,    36,   801,    37,   798,
     802,   804,    38,   812,   806,    39,    40,    41,    42,   813,
     815,   826,   809,   816,   820,   827,     1,     2,   819,   828,
    1180,  1181,   792,   829,   830,   834,     3,     4,     5,   823,
     835,   836,   837,     6,   839,   845,   838,     7,     8,     9,
     846,    10,    43,    11,    12,    13,    14,   844,   869,   344,
     671,   871,   859,  1140,   935,  1142,    15,   931,   345,    16,
     870,   886,   937,   861,   862,   887,   938,   939,   864,   940,
     941,   346,    17,   942,   944,   945,   946,   949,   976,   347,
     990,   992,  1000,  1001,  1003,    18,    19,    20,  1004,   348,
     661,    21,   662,   663,   664,   665,   666,  1089,   667,   668,
     669,   670,    22,   671,    23,  1005,    24,    25,    26,    27,
      28,  1002,   899,   723,  1006,    29,    30,   344,  1007,  1008,
      31,    32,    33,    34,  1009,  1010,   345,  1011,  1012,  1013,
      35,    36,  1014,    37,  1015,  1016,  1017,    38,  1020,   346,
      39,    40,    41,    42,   349,  1024,   344,   347,  1025,  1027,
    1032,  1033,   676,   350,   351,   345,  1034,   348,  1041,   352,
     353,   354,   355,   356,   357,   358,   359,   360,   346,  1068,
     418,  1030,  1055,  1056,   361,  1060,   347,    43,  1039,   419,
    1038,  1061,  1040,  1062,   567,   661,   348,   662,   663,   664,
     665,   666,   420,   667,   668,   669,   670,  1043,   671,  1044,
     421,  1047,  1049,  1051,  1059,  1063,  1065,  1069,   362,  1072,
     422,  1074,   349,  1066,  1070,  1084,  1086,  1077,  1106,  1135,
     933,   350,   351,    -1,  1134,  1137,  1127,   352,   353,   354,
     355,   356,   357,   358,   359,   360,  1111,  1136,  1150,  1143,
    1152,   349,   361,  1154,  1153,  1155,  1167,  1169,  1170,   918,
     350,   351,  1177,  1185,  1209,  1190,   352,   353,   354,   355,
     356,   357,   358,   359,   360,   423,  1192,  1193,  1194,  1196,
    1201,   361,  1221,  1206,   424,   425,   362,   231,  1235,  1228,
     426,   427,   428,   429,   430,   431,   432,   433,   434,  1223,
    1227,  1236,   200,   170,  1237,   435,  1229,   171,  1243,  1244,
    1254,  1258,  1247,  1257,  1262,   362,  1248,  1263,  1267,   161,
     623,   232,   233,   172,  1160,   201,  1126,  1261,   234,   456,
     991,   989,   620,   814,   952,   235,   236,   237,   977,   436,
     238,   239,   454,   240,   241,   787,  1191,   242,   243,   244,
     245,   246,   247,   248,   609,   249,   250,   251,   841,   675,
     570,   615,  1165,   252,   997,   800,   173,   174,   890,   253,
    1031,   254,   805,   331,  1037,     0,   255,     0,     0,     0,
       0,     0,     0,     0,   175,   176,     0,   256,   369,     0,
     163,   164,   165,   166,   167,   168,   169,   199,     0,   257,
     213,   200,   170,     0,     0,   258,   171,     0,     0,   370,
       0,   371,   372,     0,     0,     0,     0,   177,   178,     0,
       0,     0,   172,     0,   201,     0,     0,     0,     0,     0,
       0,   234,     0,   373,   374,     0,     0,     0,   235,     0,
       0,     0,     0,     0,   713,   328,   714,   715,   716,   717,
     718,     0,   719,   720,   721,   722,     0,   723,     0,     0,
       0,     0,     0,     0,     0,   173,   174,     0,     0,     0,
       0,   375,     0,   376,   254,   377,   333,     0,     0,     0,
       0,   378,     0,   175,   176,   379,     0,     0,   369,   334,
       0,     0,     0,   380,   381,   382,     0,     0,     0,   383,
     384,   385,     0,   213,     0,     0,     0,     0,   964,   370,
     468,   371,   372,     0,     0,     0,   177,   178,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   234,     0,   373,   374,     0,     0,     0,   235,     0,
       0,     0,     0,     0,   661,   328,   662,   663,   664,   665,
     666,     0,   667,   668,   669,   670,   661,   671,   662,   663,
     664,   665,   666,     0,   667,   668,   669,   670,     0,   671,
       0,   375,     0,   376,   254,   377,   333,     0,     0,     0,
       0,   378,     0,     0,     0,   379,     0,     0,     0,   334,
       0,     0,     0,   380,   381,   382,     0,     0,     0,   383,
     384,   385,     0,   213,     0,     0,     0,     0,  1094,     0,
       0,     0,     0,   661,     0,   662,   663,   664,   665,   666,
    1095,   667,   668,   669,   670,   661,   671,   662,   663,   664,
     665,   666,     0,   667,   668,   669,   670,   661,   671,   662,
     663,   664,   665,   666,     0,   667,   668,   669,   670,   661,
     671,   662,   663,   664,   665,   666,     0,   667,   668,   669,
     670,   661,   671,   662,   663,   664,   665,   666,     0,   667,
     668,   669,   670,     0,   671,     0,     0,  1096,     0,     0,
       0,     0,   661,     0,   662,   663,   664,   665,   666,  1097,
     667,   668,   669,   670,   661,   671,   662,   663,   664,   665,
     666,  1098,   667,   668,   669,   670,   661,   671,   662,   663,
     664,   665,   666,  1099,   667,   668,   669,   670,   661,   671,
     662,   663,   664,   665,   666,  1100,   667,   668,   669,   670,
     713,   671,   714,   715,   716,   717,   718,     0,   719,   720,
     721,   722,     0,   723,     0,     0,  1101,     0,     0,     0,
       0,   713,     0,   714,   715,   716,   717,   718,  1102,   719,
     720,   721,   722,   713,   723,   714,   715,   716,   717,   718,
    1103,   719,   720,   721,   722,   713,   723,   714,   715,   716,
     717,   718,  1104,   719,   720,   721,   722,   713,   723,   714,
     715,   716,   717,   718,  1113,   719,   720,   721,   722,   713,
     723,   714,   715,   716,   717,   718,     0,   719,   720,   721,
     722,     0,   723,     0,     0,  1114,     0,     0,     0,     0,
     713,     0,   714,   715,   716,   717,   718,  1115,   719,   720,
     721,   722,   713,   723,   714,   715,   716,   717,   718,  1116,
     719,   720,   721,   722,   713,   723,   714,   715,   716,   717,
     718,  1117,   719,   720,   721,   722,   713,   723,   714,   715,
     716,   717,   718,  1118,   719,   720,   721,   722,   661,   723,
     662,   663,   664,   665,   666,     0,   667,   668,   669,   670,
       0,   671,     0,     0,  1119,     0,     0,     0,     0,   661,
       0,   662,   663,   664,   665,   666,  1120,   667,   668,   669,
     670,   661,   671,   662,   663,   664,   665,   666,  1121,   667,
     668,   669,   670,   661,   671,   662,   663,   664,   665,   666,
    1122,   667,   668,   669,   670,   713,   671,   714,   715,   716,
     717,   718,  1125,   719,   720,   721,   722,   713,   723,   714,
     715,   716,   717,   718,     0,   719,   720,   721,   722,     0,
     723,     0,     0,  1184,     0,     0,     0,     0,   661,     0,
     662,   663,   664,   665,   666,  1213,   667,   668,   669,   670,
     661,   671,   662,   663,   664,   665,   666,  1214,   667,   668,
     669,   670,   661,   671,   662,   663,   664,   665,   666,  1217,
     667,   668,   669,   670,   661,   671,   662,   663,   664,   665,
     666,  1218,   667,   668,   669,   670,   661,   671,   662,   663,
     664,   665,   666,     0,   667,   668,   669,   670,     0,   671,
       0,     0,  1220,     0,     0,     0,     0,   661,     0,   662,
     663,   664,   665,   666,  1251,   667,   668,   669,   670,   661,
     671,   662,   663,   664,   665,   666,  1252,   667,   668,   669,
     670,     0,   671,     0,     0,     0,   672,     0,     0,     0,
       0,   661,     0,   662,   663,   664,   665,   666,   932,   667,
     668,   669,   670,   713,   671,   714,   715,   716,   717,   718,
       0,   719,   720,   721,   722,     0,   723,     0,     0,  1050,
       0,     0,     0,     0,   661,     0,   662,   663,   664,   665,
     666,  1107,   667,   668,   669,   670,   661,   671,   662,   663,
     664,   665,   666,     0,   667,   668,   669,   670,     0,   671,
       0,     0,     0,  1138,     0,     0,     0,     0,   713,     0,
     714,   715,   716,   717,   718,  1158,   719,   720,   721,   722,
     661,   723,   662,   663,   664,   665,   666,     0,   667,   668,
     669,   670,     0,   671,     0,     0,  1168,     0,     0,     0,
       0,   661,     0,   662,   663,   664,   665,   666,  1178,   667,
     668,   669,   670,   661,   671,   662,   663,   664,   665,   666,
       0,   667,   668,   669,   670,     0,   671,     0,     0,     0,
    1200,     0,     0,     0,     0,   661,     0,   662,   663,   664,
     665,   666,  1212,   667,   668,   669,   670,   661,   671,   662,
     663,   664,   665,   666,     0,   667,   668,   669,   670,     0,
     671,     0,     0,  1216,     0,     0,     0,     0,   661,     0,
     662,   663,   664,   665,   666,  1222,   667,   668,   669,   670,
     661,   671,   662,   663,   664,   665,   666,     0,   667,   668,
     669,   670,     0,   671,     0,     0,     0,  1230,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   661,  1233,
     662,   663,   664,   665,   666,  1090,   667,   668,   669,   670,
       0,   671,     0,   713,     0,   714,   715,   716,   717,   718,
    1234,   719,   720,   721,   722,     0,   723,     0,     0,     0,
       0,     0,  1253,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,     0,     0,     0,   509,   510,     0,
     511,   512,   513,   514,   661,     0,   662,   663,   664,   665,
     666,  1093,   667,   668,   669,   670,   713,   671,   714,   715,
     716,   717,   718,  1109,   719,   720,   721,   722,   713,   723,
     714,   715,   716,   717,   718,  1110,   719,   720,   721,   722,
     661,   723,   662,   663,   664,   665,   666,  1112,   667,   668,
     669,   670,   661,   671,   662,   663,   664,   665,   666,  1215,
     667,   668,   669,   670,   661,   671,   662,   663,   664,   665,
     666,  1219,   667,   668,   669,   670,   661,   671,   662,   663,
     664,   665,   666,     0,   667,   668,   669,   670,     0,   671
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   228,   115,   109,   127,    98,    31,   109,   629,
     121,   109,    36,    37,    38,   134,   109,   109,    42,   109,
     141,   546,   121,   145,   350,   351,   610,   629,   612,   109,
     135,   136,   109,   617,    93,   560,   362,   135,   136,    93,
     141,   367,   135,   136,   109,   135,   136,   109,    93,   809,
     633,   143,   109,   109,   109,   109,   586,   109,   109,   109,
      31,   109,   109,  1061,   109,   109,   897,   592,   593,   594,
      42,    42,    51,    82,    42,    42,   601,    42,    33,    51,
      51,     6,   200,    51,    51,     6,    51,    51,    51,  1077,
    1078,  1079,  1193,    64,    82,    31,   322,   100,    51,   145,
     342,    72,    51,    42,   145,    24,    42,    82,   145,    82,
      82,    82,    51,     4,   145,    51,   145,   120,    25,    44,
      33,    42,    51,    44,    33,    56,    45,    52,    64,    42,
      51,    52,    82,    42,    82,    33,    72,    82,    57,   148,
      31,    42,    82,    82,    63,    82,    82,   200,    82,    80,
      51,    22,    82,    51,    82,   201,   681,   145,    82,   200,
      82,    82,  1263,   688,    83,   202,   137,    92,    82,   200,
      82,    92,    82,   202,   146,   146,   147,   145,    82,    82,
      82,   152,   153,   154,   155,   156,   157,   158,   159,   160,
     145,    82,   201,   100,    22,   204,   167,   145,   148,   208,
     209,   137,   127,    31,    51,    82,   127,   146,   148,    42,
     146,   147,   200,  1201,   456,   740,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   146,   201,   206,   201,     4,
     201,   167,   203,   201,   201,   207,   201,   767,   201,   764,
     765,   205,    31,   768,   769,   200,  1244,   772,   201,   774,
     775,   200,   200,    42,   204,   200,   582,   200,   208,   209,
     131,   201,    51,   200,   204,   201,   200,   203,   208,   209,
     200,    31,   200,   201,   200,    64,   200,   861,   200,   201,
     864,   201,    42,    72,   200,   201,   200,   201,   200,   201,
     200,    51,    82,    82,   145,   145,   200,   200,   200,   146,
     147,   145,    33,   131,    64,   424,   425,    82,    33,    82,
      31,    82,    72,   146,   147,   200,    82,   436,   644,   645,
     646,   647,    82,   200,   650,   651,   652,   653,   654,   655,
     656,   657,   658,   659,  1165,   661,   662,   663,   664,   665,
     666,   667,   668,   669,   670,   671,    82,   673,   137,    82,
      22,   202,   202,   679,   200,   145,   200,   146,   147,    31,
     200,    82,   208,   152,   153,   154,   155,   156,   157,   158,
     159,   160,   145,   699,   145,    82,  1136,   137,   167,   145,
      22,    82,    82,   200,    82,    31,   146,   147,    82,    82,
      22,   145,   152,   153,   154,   155,   156,   157,   158,   159,
     160,    60,    82,    82,   145,   731,    82,   167,    82,   145,
     200,    82,   201,   138,   203,   140,   141,   142,   143,   144,
      82,   146,   147,   148,   149,    22,   151,   200,   200,   200,
      22,    31,   145,    31,   200,   200,    82,   145,   145,    31,
      82,   201,    82,   203,   145,   145,   201,   145,   202,   201,
      82,   145,   145,    82,   145,    82,   145,   145,   145,   131,
     119,   202,   145,   200,   200,   145,   145,    42,    33,   145,
     129,   145,    22,    82,   145,   200,    51,  1061,   120,  1063,
     572,    31,    82,   145,    82,    82,   812,    82,   120,   202,
      82,    82,  1075,   200,   202,   200,    82,   602,  1082,   200,
     200,   602,   200,   562,   602,   200,   200,   200,   562,   602,
     602,   202,   602,   202,   202,   202,   145,   562,   145,   202,
     200,   200,   602,   120,   200,   602,   200,   638,   120,   200,
      31,    31,    82,   200,   200,   201,    23,   602,   200,   638,
     602,   208,  1162,  1163,  1164,   602,   602,   602,   602,  1074,
     602,   602,   602,   577,   602,   602,   580,   602,   602,   145,
    1162,  1163,  1164,    42,   148,   149,  1150,   151,  1152,  1153,
     120,   200,    51,   200,   900,   200,   201,   696,   697,   603,
     118,    82,    82,   702,   703,   704,   705,   706,   707,   708,
     709,   710,   711,    42,   713,   714,   715,   716,   717,   718,
     719,   720,   721,   722,   723,   931,   200,   201,  1192,    22,
    1194,  1231,    79,   734,   200,    33,   738,   202,    31,   200,
     725,  1205,   202,   208,    91,    31,    82,   725,   208,  1231,
     200,   201,   725,   734,   200,   725,    42,   104,    82,  1223,
     206,   200,   201,  1227,   736,    51,   200,   201,    82,     9,
      10,    11,    12,    13,    14,    15,    82,    38,    64,  1243,
       5,    21,   988,    42,   200,    25,    72,    82,   692,    82,
     148,   149,    82,   151,   126,    20,    82,    51,  1262,   200,
     205,    41,    33,  1267,   146,   147,   148,   149,   126,   151,
     126,   126,    33,    33,    39,    40,   202,    33,    43,    33,
      33,    46,   726,   727,    33,    42,    33,   120,    53,    33,
      33,    33,    33,    33,    51,   739,   145,    33,   742,    33,
      33,   145,   202,    82,    84,    85,   202,    64,    33,   200,
     200,   137,    82,    82,    33,    72,    81,    33,    82,    33,
     146,   147,   102,   103,    89,    82,   152,   153,   154,   155,
     156,   157,   158,   159,   160,    33,   875,    33,   870,    33,
     105,   167,  1088,  1089,  1090,  1091,    33,  1093,   146,   147,
     148,   149,   117,   151,    33,   135,   136,    33,   123,    33,
    1106,    33,    33,    33,    33,    33,  1112,    33,    33,    33,
     135,   136,    82,    82,    82,   201,    82,   203,   145,   145,
     137,    82,    82,   145,  1130,    82,   201,    33,   202,   146,
     147,  1137,   200,   202,    33,   152,   153,   154,   155,   156,
     157,   158,   159,   160,    82,   202,   201,    82,   201,  1155,
     167,   201,   201,   201,   201,   201,   201,   201,   201,   201,
     201,  1167,    33,  1169,     0,  1071,   201,   201,   201,    33,
      33,     7,     8,    33,    33,    33,    33,    33,    33,    33,
      33,    17,    18,    19,   201,    33,   203,   202,    24,   126,
      82,   201,    28,    29,    30,   200,    32,    82,    34,    35,
      36,    37,   201,   201,   201,    33,   201,   126,   126,  1215,
     201,    47,   201,  1219,    50,  1221,    82,   201,   201,   201,
     201,   201,   126,   201,   201,    33,    33,    63,    33,    33,
      82,  1237,    33,    33,    33,    33,    33,    33,    33,    33,
      76,    77,    78,    33,    33,    33,    82,    33,    33,    33,
      33,    33,    33,  1259,    33,    33,    33,    93,    33,    95,
      33,    97,    98,    99,   100,   101,    33,   145,    82,   202,
     106,   107,    82,  1072,   201,   111,   112,   113,   114,    51,
      51,    51,    51,    51,   200,   121,   122,    51,   124,    82,
      51,    51,   128,    33,    51,   131,   132,   133,   134,    51,
     200,    51,    82,   200,   200,    51,     7,     8,    82,    51,
    1109,  1110,   205,    51,    51,    51,    17,    18,    19,    82,
      51,    51,    51,    24,   201,    33,    82,    28,    29,    30,
     145,    32,   168,    34,    35,    36,    37,    82,    82,    42,
     151,    82,   145,  1047,    51,  1049,    47,    82,    51,    50,
     200,   200,    51,   145,   145,   200,    51,    51,   145,    51,
      51,    64,    63,    51,    51,    51,    51,   126,   200,    72,
     200,   200,   126,   126,    51,    76,    77,    78,    51,    82,
     138,    82,   140,   141,   142,   143,   144,   145,   146,   147,
     148,   149,    93,   151,    95,    51,    97,    98,    99,   100,
     101,   126,   145,   151,    51,   106,   107,    42,    51,    51,
     111,   112,   113,   114,    51,    51,    51,    51,    51,    51,
     121,   122,    51,   124,    51,    51,    51,   128,    51,    64,
     131,   132,   133,   134,   137,    51,    42,    72,    51,    51,
      42,    42,   145,   146,   147,    51,    51,    82,    51,   152,
     153,   154,   155,   156,   157,   158,   159,   160,    64,    51,
      42,   200,    82,    82,   167,   145,    72,   168,   200,    51,
     204,   145,   200,    82,   201,   138,    82,   140,   141,   142,
     143,   144,    64,   146,   147,   148,   149,   202,   151,   202,
      72,   202,   202,   207,   200,   145,    82,    51,   201,    33,
      82,    33,   137,   200,   202,    82,    82,   200,    82,    82,
     145,   146,   147,   151,   204,    33,   200,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   202,   202,   145,    82,
     145,   137,   167,   200,   145,    33,    33,    33,    51,   202,
     146,   147,   200,    51,   204,   207,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   137,   145,   145,   145,   200,
     200,   167,    33,   202,   146,   147,   201,     5,    51,   200,
     152,   153,   154,   155,   156,   157,   158,   159,   160,   145,
     145,    51,    20,    21,   130,   167,   200,    25,   145,   145,
      82,    51,   202,   200,   145,   201,   204,   145,   145,    44,
     320,    39,    40,    41,  1075,    43,   983,   207,    46,   138,
     736,   734,   317,   572,   698,    53,    54,    55,   725,   201,
      58,    59,   136,    61,    62,   548,  1144,    65,    66,    67,
      68,    69,    70,    71,   306,    73,    74,    75,   602,   365,
     211,   312,  1080,    81,   741,   557,    84,    85,   638,    87,
     777,    89,   562,   120,   793,    -1,    94,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   102,   103,    -1,   105,     3,    -1,
       9,    10,    11,    12,    13,    14,    15,    16,    -1,   117,
     118,    20,    21,    -1,    -1,   123,    25,    -1,    -1,    24,
      -1,    26,    27,    -1,    -1,    -1,    -1,   135,   136,    -1,
      -1,    -1,    41,    -1,    43,    -1,    -1,    -1,    -1,    -1,
      -1,    46,    -1,    48,    49,    -1,    -1,    -1,    53,    -1,
      -1,    -1,    -1,    -1,   138,    60,   140,   141,   142,   143,
     144,    -1,   146,   147,   148,   149,    -1,   151,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    84,    85,    -1,    -1,    -1,
      -1,    86,    -1,    88,    89,    90,    91,    -1,    -1,    -1,
      -1,    96,    -1,   102,   103,   100,    -1,    -1,     3,   104,
      -1,    -1,    -1,   108,   109,   110,    -1,    -1,    -1,   114,
     115,   116,    -1,   118,    -1,    -1,    -1,    -1,   202,    24,
     125,    26,    27,    -1,    -1,    -1,   135,   136,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    46,    -1,    48,    49,    -1,    -1,    -1,    53,    -1,
      -1,    -1,    -1,    -1,   138,    60,   140,   141,   142,   143,
     144,    -1,   146,   147,   148,   149,   138,   151,   140,   141,
     142,   143,   144,    -1,   146,   147,   148,   149,    -1,   151,
      -1,    86,    -1,    88,    89,    90,    91,    -1,    -1,    -1,
      -1,    96,    -1,    -1,    -1,   100,    -1,    -1,    -1,   104,
      -1,    -1,    -1,   108,   109,   110,    -1,    -1,    -1,   114,
     115,   116,    -1,   118,    -1,    -1,    -1,    -1,   202,    -1,
      -1,    -1,    -1,   138,    -1,   140,   141,   142,   143,   144,
     202,   146,   147,   148,   149,   138,   151,   140,   141,   142,
     143,   144,    -1,   146,   147,   148,   149,   138,   151,   140,
     141,   142,   143,   144,    -1,   146,   147,   148,   149,   138,
     151,   140,   141,   142,   143,   144,    -1,   146,   147,   148,
     149,   138,   151,   140,   141,   142,   143,   144,    -1,   146,
     147,   148,   149,    -1,   151,    -1,    -1,   202,    -1,    -1,
      -1,    -1,   138,    -1,   140,   141,   142,   143,   144,   202,
     146,   147,   148,   149,   138,   151,   140,   141,   142,   143,
     144,   202,   146,   147,   148,   149,   138,   151,   140,   141,
     142,   143,   144,   202,   146,   147,   148,   149,   138,   151,
     140,   141,   142,   143,   144,   202,   146,   147,   148,   149,
     138,   151,   140,   141,   142,   143,   144,    -1,   146,   147,
     148,   149,    -1,   151,    -1,    -1,   202,    -1,    -1,    -1,
      -1,   138,    -1,   140,   141,   142,   143,   144,   202,   146,
     147,   148,   149,   138,   151,   140,   141,   142,   143,   144,
     202,   146,   147,   148,   149,   138,   151,   140,   141,   142,
     143,   144,   202,   146,   147,   148,   149,   138,   151,   140,
     141,   142,   143,   144,   202,   146,   147,   148,   149,   138,
     151,   140,   141,   142,   143,   144,    -1,   146,   147,   148,
     149,    -1,   151,    -1,    -1,   202,    -1,    -1,    -1,    -1,
     138,    -1,   140,   141,   142,   143,   144,   202,   146,   147,
     148,   149,   138,   151,   140,   141,   142,   143,   144,   202,
     146,   147,   148,   149,   138,   151,   140,   141,   142,   143,
     144,   202,   146,   147,   148,   149,   138,   151,   140,   141,
     142,   143,   144,   202,   146,   147,   148,   149,   138,   151,
     140,   141,   142,   143,   144,    -1,   146,   147,   148,   149,
      -1,   151,    -1,    -1,   202,    -1,    -1,    -1,    -1,   138,
      -1,   140,   141,   142,   143,   144,   202,   146,   147,   148,
     149,   138,   151,   140,   141,   142,   143,   144,   202,   146,
     147,   148,   149,   138,   151,   140,   141,   142,   143,   144,
     202,   146,   147,   148,   149,   138,   151,   140,   141,   142,
     143,   144,   202,   146,   147,   148,   149,   138,   151,   140,
     141,   142,   143,   144,    -1,   146,   147,   148,   149,    -1,
     151,    -1,    -1,   202,    -1,    -1,    -1,    -1,   138,    -1,
     140,   141,   142,   143,   144,   202,   146,   147,   148,   149,
     138,   151,   140,   141,   142,   143,   144,   202,   146,   147,
     148,   149,   138,   151,   140,   141,   142,   143,   144,   202,
     146,   147,   148,   149,   138,   151,   140,   141,   142,   143,
     144,   202,   146,   147,   148,   149,   138,   151,   140,   141,
     142,   143,   144,    -1,   146,   147,   148,   149,    -1,   151,
      -1,    -1,   202,    -1,    -1,    -1,    -1,   138,    -1,   140,
     141,   142,   143,   144,   202,   146,   147,   148,   149,   138,
     151,   140,   141,   142,   143,   144,   202,   146,   147,   148,
     149,    -1,   151,    -1,    -1,    -1,   200,    -1,    -1,    -1,
      -1,   138,    -1,   140,   141,   142,   143,   144,   200,   146,
     147,   148,   149,   138,   151,   140,   141,   142,   143,   144,
      -1,   146,   147,   148,   149,    -1,   151,    -1,    -1,   200,
      -1,    -1,    -1,    -1,   138,    -1,   140,   141,   142,   143,
     144,   200,   146,   147,   148,   149,   138,   151,   140,   141,
     142,   143,   144,    -1,   146,   147,   148,   149,    -1,   151,
      -1,    -1,    -1,   200,    -1,    -1,    -1,    -1,   138,    -1,
     140,   141,   142,   143,   144,   200,   146,   147,   148,   149,
     138,   151,   140,   141,   142,   143,   144,    -1,   146,   147,
     148,   149,    -1,   151,    -1,    -1,   200,    -1,    -1,    -1,
      -1,   138,    -1,   140,   141,   142,   143,   144,   200,   146,
     147,   148,   149,   138,   151,   140,   141,   142,   143,   144,
      -1,   146,   147,   148,   149,    -1,   151,    -1,    -1,    -1,
     200,    -1,    -1,    -1,    -1,   138,    -1,   140,   141,   142,
     143,   144,   200,   146,   147,   148,   149,   138,   151,   140,
     141,   142,   143,   144,    -1,   146,   147,   148,   149,    -1,
     151,    -1,    -1,   200,    -1,    -1,    -1,    -1,   138,    -1,
     140,   141,   142,   143,   144,   200,   146,   147,   148,   149,
     138,   151,   140,   141,   142,   143,   144,    -1,   146,   147,
     148,   149,    -1,   151,    -1,    -1,    -1,   200,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   138,   200,
     140,   141,   142,   143,   144,   145,   146,   147,   148,   149,
      -1,   151,    -1,   138,    -1,   140,   141,   142,   143,   144,
     200,   146,   147,   148,   149,    -1,   151,    -1,    -1,    -1,
      -1,    -1,   200,   169,   170,   171,   172,   173,   174,   175,
     176,   177,   178,   179,   180,   181,   182,   183,   184,   185,
     186,   187,   188,   189,    -1,    -1,    -1,   193,   194,    -1,
     196,   197,   198,   199,   138,    -1,   140,   141,   142,   143,
     144,   145,   146,   147,   148,   149,   138,   151,   140,   141,
     142,   143,   144,   145,   146,   147,   148,   149,   138,   151,
     140,   141,   142,   143,   144,   145,   146,   147,   148,   149,
     138,   151,   140,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   138,   151,   140,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   138,   151,   140,   141,   142,   143,
     144,   145,   146,   147,   148,   149,   138,   151,   140,   141,
     142,   143,   144,    -1,   146,   147,   148,   149,    -1,   151
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,    17,    18,    19,    24,    28,    29,    30,
      32,    34,    35,    36,    37,    47,    50,    63,    76,    77,
      78,    82,    93,    95,    97,    98,    99,   100,   101,   106,
     107,   111,   112,   113,   114,   121,   122,   124,   128,   131,
     132,   133,   134,   168,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   224,   225,   226,   227,   230,   232,   235,
     240,   251,   252,   256,   260,   263,   266,   269,   275,   281,
     284,   289,   292,   295,   298,   299,   302,   303,   305,   306,
     307,   311,   312,   313,   314,   320,   323,   329,   332,   333,
      51,   201,    51,   201,   200,   201,   200,   200,   201,    33,
      42,    51,    82,   201,    82,   201,   200,    82,   200,   201,
     272,   200,   200,   200,   200,   200,   201,    33,    42,   200,
     201,   201,   200,    33,   200,   200,   200,   201,   272,   272,
      82,   223,    33,    51,   321,   201,   201,   272,   200,    33,
     200,   201,   200,   201,   200,   201,   272,   200,   201,   272,
     272,    82,   220,    82,   221,    82,   222,   272,   200,   201,
       0,   212,   200,     9,    10,    11,    12,    13,    14,    15,
      21,    25,    41,    84,    85,   102,   103,   135,   136,   326,
     327,   328,   357,   358,   359,   360,   361,   392,   393,   395,
     396,   399,   400,   401,   402,   403,   404,   405,   200,    16,
      20,    43,   327,   330,   331,   365,   382,   406,    23,     4,
      82,   308,   309,   118,   264,   265,   338,    42,   200,    51,
     200,   200,   208,    82,   200,   208,    82,    82,   233,   234,
      33,     5,    39,    40,    46,    53,    54,    55,    58,    59,
      61,    62,    65,    66,    67,    68,    69,    70,    71,    73,
      74,    75,    81,    87,    89,    94,   105,   117,   123,   290,
     291,   338,   348,   357,   358,   359,   360,   361,   362,   363,
     364,   365,   366,   367,   368,   369,   370,   371,   372,   373,
     374,   375,   376,   377,   378,   379,   380,   381,   382,   383,
     384,   385,   387,   388,   392,   393,   394,   395,   396,   397,
      82,   145,   200,    22,    82,   120,   276,   277,   278,    22,
      82,   120,   285,   286,    22,    82,   120,   282,   283,    82,
     236,   237,   233,    38,   231,    42,   200,   241,    60,   119,
     129,   340,    79,    91,   104,   315,   316,   389,   390,   391,
      22,   131,   253,   254,    42,    51,    64,    72,    82,   137,
     146,   147,   152,   153,   154,   155,   156,   157,   158,   159,
     160,   167,   201,   228,    82,   300,   301,    82,   304,     3,
      24,    26,    27,    48,    49,    86,    88,    90,    96,   100,
     108,   109,   110,   114,   115,   116,   271,   337,   338,   339,
     340,   341,   342,   343,   344,   345,   346,   347,   348,   349,
     350,   351,   352,   354,   355,   356,   364,   386,   390,   391,
     200,   200,   126,    82,   145,   200,    51,   200,    42,    51,
      64,    72,    82,   137,   146,   147,   152,   153,   154,   155,
     156,   157,   158,   159,   160,   167,   201,   248,   250,   293,
     294,   348,   364,   365,   374,   380,   381,   382,   383,   384,
     385,   392,   393,   394,   293,   200,   253,   205,   267,   268,
     351,   357,   261,   262,   338,   270,   271,   200,   125,   271,
     324,   325,   398,   200,   200,   126,    82,   145,   200,   126,
      82,   145,   200,   126,    82,   145,   200,   200,   169,   170,
     171,   172,   173,   174,   175,   176,   177,   178,   179,   180,
     181,   182,   183,   184,   185,   186,   187,   188,   189,   193,
     194,   196,   197,   198,   199,   334,   335,   407,   408,   409,
     410,   411,   412,   413,   414,   415,   416,   417,   418,   419,
     420,   421,   422,   423,   424,   425,   426,   427,   428,   429,
     430,   431,   432,   433,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,   145,   202,    33,
      33,    33,   145,   202,   202,    82,   145,   201,   310,    31,
     309,    33,   145,   202,   200,   200,    82,   202,   208,    82,
     202,   208,    33,    31,   234,    82,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   145,   202,    33,    82,    82,    82,    31,   277,
     145,    82,   145,    82,    31,   286,    82,   145,    82,    31,
     283,   201,    31,   237,    31,    33,   202,   200,   203,   246,
     247,   248,   249,   145,   202,   202,   202,    33,   145,   202,
      82,    82,    31,   254,   201,   201,   201,   201,   228,   228,
     201,   201,   201,   201,   201,   201,   201,   201,   201,   201,
     228,   138,   140,   141,   142,   143,   144,   146,   147,   148,
     149,   151,   200,   201,    31,   301,   145,   228,    31,    82,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   202,   126,    82,   200,   201,   201,   201,   201,
     248,   248,   201,   201,   201,   201,   201,   201,   201,   201,
     201,   201,   248,   138,   140,   141,   142,   143,   144,   146,
     147,   148,   149,   151,   322,   145,   202,   202,    31,    42,
      51,   201,   258,   259,   145,   202,   145,   202,   145,   202,
      33,   145,   202,   126,    82,   126,    82,   126,    82,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,   202,   145,    42,    51,
     336,    42,   146,   147,   274,   336,    51,   274,    51,    82,
      51,    51,   205,   436,   437,    51,    51,    82,    82,   434,
     328,    51,    51,   336,    51,   331,    51,   200,   201,    82,
      42,    51,    33,    51,   265,   200,   200,   200,   272,    82,
     200,   200,   272,    82,   228,   437,    51,    51,    51,    51,
      51,   336,   336,   336,    51,    51,    51,    51,    82,   201,
     336,   291,   200,   272,    82,    33,   145,     6,    42,    44,
      51,    52,    82,    92,   127,   146,   279,   287,   288,   145,
     288,   145,   145,   288,   145,    51,   146,   147,   273,    82,
     200,    82,    31,   247,   249,    33,   200,    45,    57,    63,
      83,   238,   239,   352,   353,   245,   200,   200,    56,    80,
     316,    82,   148,   204,   208,   209,   317,   318,   319,   145,
      33,   145,   200,   228,   228,   228,   229,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   202,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,    82,   200,   145,   228,    51,   336,    51,    51,    51,
      51,    51,    51,   336,    51,    51,    51,   200,   272,   126,
     248,   248,   273,   228,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   202,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   200,   294,   200,   272,
     200,   272,   228,   200,   206,    42,    51,   145,   201,   268,
     200,   262,   200,   271,   200,   272,   336,   325,   200,   272,
     126,   126,   126,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    51,    51,   336,   336,
      51,   437,   336,   336,    51,    51,   336,    51,   336,   336,
     200,   334,    42,    42,    51,   435,   206,   435,   204,   200,
     200,    51,   310,   202,   202,   228,   200,   202,   200,   202,
     200,   207,   296,   297,   200,    82,    82,    42,    51,   200,
     145,   145,    82,   145,   288,    82,   200,   288,    51,    51,
     202,   233,    33,   248,    33,   145,   202,   200,   243,   242,
     145,   200,   201,   319,    82,   228,    82,   100,   120,   145,
     145,   145,   202,   145,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   228,    82,   200,   200,   145,
     145,   202,   145,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   200,   200,   202,   259,   200,    42,    51,
     201,   228,   200,   200,   204,    82,   202,    33,   200,   200,
     272,   200,   272,    82,   145,   202,   280,   288,   287,   288,
     145,   288,   145,   145,   200,    33,    31,   248,   200,   336,
     239,   244,   246,   246,   246,   318,   288,    33,   200,    33,
      51,   255,   228,   228,   228,   228,   228,   200,   200,   228,
     248,   248,   228,   228,   202,    51,   310,   228,   200,   200,
     207,   296,   145,   145,   145,   288,   200,   288,   288,   228,
     200,   200,    31,    31,    31,   201,   202,   228,   228,   204,
     145,   200,   200,   202,   202,   145,   200,   202,   202,   145,
     202,    33,   200,   145,   288,   280,   288,   145,   200,   200,
     200,   246,   288,   200,   200,    51,    51,   130,   228,   228,
     228,   207,   288,   145,   145,   288,    31,   202,   204,   228,
     257,   202,   202,   200,    82,   288,   287,   200,    51,   145,
     200,   207,   145,   145,   228,   288,   280,   145,   288
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
     254,   254,   254,   254,   254,   255,   255,   255,   255,   256,
     257,   257,   258,   258,   259,   259,   259,   259,   259,   259,
     259,   259,   259,   260,   260,   261,   261,   262,   263,   263,
     264,   264,   265,   266,   266,   267,   267,   268,   268,   269,
     269,   269,   269,   270,   270,   271,   271,   271,   271,   271,
     271,   271,   271,   271,   271,   271,   271,   271,   271,   271,
     271,   271,   271,   271,   271,   271,   271,   271,   272,   272,
     272,   272,   272,   272,   273,   273,   273,   274,   274,   274,
     275,   276,   276,   277,   278,   278,   278,   279,   279,   279,
     279,   279,   280,   280,   280,   280,   281,   282,   282,   283,
     283,   283,   284,   285,   285,   286,   286,   286,   287,   287,
     287,   287,   287,   288,   288,   288,   288,   288,   288,   289,
     289,   289,   289,   290,   290,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   292,   292,   293,   293,   294,   294,
     294,   294,   294,   294,   294,   294,   294,   294,   294,   294,
     294,   295,   295,   296,   296,   297,   297,   298,   299,   300,
     300,   301,   302,   303,   304,   304,   304,   304,   305,   306,
     306,   306,   306,   307,   308,   308,   309,   309,   309,   310,
     310,   310,   311,   311,   312,   312,   312,   312,   312,   312,
     313,   313,   313,   313,   313,   313,   314,   315,   315,   316,
     316,   316,   317,   317,   317,   317,   318,   318,   319,   319,
     319,   319,   319,   321,   322,   320,   323,   323,   323,   323,
     324,   324,   325,   325,   326,   326,   326,   326,   326,   326,
     326,   327,   327,   327,   327,   327,   327,   327,   327,   327,
     327,   328,   328,   329,   329,   330,   330,   330,   330,   331,
     331,   332,   332,   333,   333,   334,   334,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   336,   336,   337,   338,   339,   340,
     341,   342,   343,   344,   345,   346,   347,   348,   349,   350,
     351,   352,   353,   354,   355,   356,   357,   358,   358,   359,
     360,   361,   362,   363,   364,   364,   365,   366,   367,   368,
     369,   370,   371,   372,   373,   374,   375,   376,   377,   378,
     379,   380,   381,   382,   383,   384,   385,   386,   387,   388,
     389,   389,   390,   391,   392,   393,   394,   395,   396,   397,
     398,   399,   400,   401,   402,   403,   404,   405,   406,   407,
     408,   409,   410,   411,   412,   413,   414,   415,   416,   417,
     418,   419,   420,   421,   422,   423,   424,   425,   426,   427,
     428,   429,   430,   431,   432,   433,   434,   435,   435,   436,
     436,   437
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
       4,     4,     4,     6,     6,     4,     8,     1,     3,     4,
       7,     3,     4,     2,     1,     4,     4,     2,     1,     7,
       3,     1,     1,     1,     1,     1,     1,     0,     5,     0,
       8,     0,     8,     0,    10,     0,     8,     2,     2,     1,
       1,     4,     2,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     2,     2,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       6,     6,     8,     5,     1,     4,     4,     4,     2,     1,
       9,     6,     5,     7,     7,     3,     5,     3,     1,     6,
       3,     1,     3,     1,     5,     3,     3,     4,     2,     2,
       3,     1,     1,     2,     5,     3,     1,     1,     2,     5,
       3,     1,     1,     2,     5,     3,     1,     1,     1,     2,
       5,     3,     6,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     2,     4,
       3,     5,     1,     3,     2,     2,     1,     2,     2,     1,
       4,     2,     1,     4,     2,     1,     4,     3,     5,     9,
       1,     5,     3,     5,     7,     9,     4,     2,     1,     5,
       7,     4,     4,     2,     1,     7,     9,     6,     1,     1,
       1,     1,     1,     0,     1,     1,     1,     2,     2,     2,
       5,     3,     6,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     5,     6,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     5,     6,     7,     5,     1,     3,     3,     4,     2,
       1,     5,     3,     4,     4,     6,     3,     5,     3,     2,
       5,     3,     6,     4,     2,     1,     5,     7,     9,     0,
       3,     3,     2,     5,     5,     6,     3,     7,     8,     5,
       5,     6,     3,     7,     8,     5,     6,     3,     1,     1,
       1,     1,     1,     3,     4,     6,     1,     2,     1,     1,
       1,     1,     1,     0,     0,     5,     2,     5,     3,     6,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     1,     3,     6,     1,     1,     1,     1,     3,
       1,     3,     6,     2,     5,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     3,     3,     1,
       3,     3,     3,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     3,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     1,     1,     3,     3,     3,     3,
       5,     3,     3,     3,     1,     3,     3,     3,     1,     1,
       1,     1,     1,     3,     1,     1,     1,     1,     3,     3,
       3,     3,     1,     1,     3,     3,     3,     1,     1,     1,
       3,     3,     3,     3,     3,     3,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     3,     2,
       2,     2
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
  "CUTOFF", "DATAFILE", "DR_ALGO", "DROP", "DSAMPLE", "DYNASAVE",
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
  "VAROBS", "XLS_SHEET", "XLS_RANGE", "NORMCDF", "EXCLAMATION_EQUAL",
  "EXCLAMATION", "EQUAL_EQUAL", "GREATER_EQUAL", "LESS_EQUAL", "GREATER",
  "LESS", "COMMA", "MINUS", "PLUS", "DIVIDE", "TIMES", "UMINUS", "POWER",
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
     256,    -1,   260,    -1,   263,    -1,   266,    -1,   269,    -1,
     289,    -1,   292,    -1,   295,    -1,   275,    -1,   284,    -1,
     281,    -1,   298,    -1,   299,    -1,   302,    -1,   214,    -1,
     215,    -1,   303,    -1,   305,    -1,   306,    -1,   307,    -1,
     311,    -1,   312,    -1,   313,    -1,   314,    -1,   320,    -1,
     323,    -1,   329,    -1,   332,    -1,   333,    -1,   219,    -1,
     216,    -1,   217,    -1,   218,    -1,    28,    51,   200,    -1,
      28,    51,    51,   200,    -1,   111,   272,   200,    -1,   131,
     220,   200,    -1,   132,   221,   200,    -1,   133,   222,   200,
      -1,    99,   223,   200,    -1,   220,    82,    -1,   220,   145,
      82,    -1,    82,    -1,   220,    82,   126,    -1,   220,   145,
      82,   126,    -1,    82,   126,    -1,   221,    82,    -1,   221,
     145,    82,    -1,    82,    -1,   221,    82,   126,    -1,   221,
     145,    82,   126,    -1,    82,   126,    -1,   222,    82,    -1,
     222,   145,    82,    -1,    82,    -1,   222,    82,   126,    -1,
     222,   145,    82,   126,    -1,    82,   126,    -1,   223,    82,
      -1,   223,   145,    82,    -1,    82,    -1,   223,    82,   126,
      -1,   223,   145,    82,   126,    -1,    82,   126,    -1,   100,
      51,   200,    -1,   100,    33,    51,   200,    -1,    24,    42,
     200,    -1,    24,    33,    42,   200,    -1,    63,    42,   200,
      -1,    63,    33,    42,   200,    -1,    82,    33,   228,   200,
      -1,   201,   228,   202,    -1,    82,    -1,    42,    -1,    51,
      -1,   228,   147,   228,    -1,   228,   146,   228,    -1,   228,
     148,   228,    -1,   228,   149,   228,    -1,   228,   151,   228,
      -1,   228,   144,   228,    -1,   228,   143,   228,    -1,   228,
     142,   228,    -1,   228,   141,   228,    -1,   228,   140,   228,
      -1,   228,   138,   228,    -1,   146,   228,    -1,   147,   228,
      -1,   152,   201,   228,   202,    -1,   153,   201,   228,   202,
      -1,   154,   201,   228,   202,    -1,   155,   201,   228,   202,
      -1,   156,   201,   228,   202,    -1,   157,   201,   228,   202,
      -1,   158,   201,   228,   202,    -1,   159,   201,   228,   202,
      -1,   160,   201,   228,   202,    -1,   167,   201,   228,   202,
      -1,    64,   201,   228,   145,   228,   202,    -1,    72,   201,
     228,   145,   228,   202,    -1,    82,   201,   229,   202,    -1,
     137,   201,   228,   145,   228,   145,   228,   202,    -1,   228,
      -1,   229,   145,   228,    -1,    50,   200,   233,    31,    -1,
      50,   201,   231,   202,   200,   233,    31,    -1,    38,    33,
      82,    -1,    32,   200,   233,    31,    -1,   233,   234,    -1,
     234,    -1,    82,    33,   228,   200,    -1,    47,   200,   236,
      31,    -1,   236,   237,    -1,   237,    -1,    82,   201,   273,
     202,    33,   228,   200,    -1,   238,   145,   239,    -1,   239,
      -1,    57,    -1,    45,    -1,    83,    -1,   352,    -1,   353,
      -1,    -1,    76,   200,   241,   246,    31,    -1,    -1,    76,
     201,   340,   202,   200,   242,   246,    31,    -1,    -1,    76,
     201,   129,   202,   200,   243,   246,    31,    -1,    -1,    76,
     201,   119,   145,   238,   202,   244,   200,   246,    31,    -1,
      -1,    76,   201,   119,   202,   245,   200,   246,    31,    -1,
     246,   247,    -1,   246,   249,    -1,   247,    -1,   249,    -1,
     248,    33,   248,   200,    -1,   248,   200,    -1,   201,   248,
     202,    -1,   250,    -1,    42,    -1,    51,    -1,   248,   147,
     248,    -1,   248,   146,   248,    -1,   248,   148,   248,    -1,
     248,   149,   248,    -1,   248,   144,   248,    -1,   248,   143,
     248,    -1,   248,   142,   248,    -1,   248,   141,   248,    -1,
     248,   140,   248,    -1,   248,   138,   248,    -1,   248,   151,
     248,    -1,   146,   248,    -1,   147,   248,    -1,   152,   201,
     248,   202,    -1,   153,   201,   248,   202,    -1,   154,   201,
     248,   202,    -1,   155,   201,   248,   202,    -1,   156,   201,
     248,   202,    -1,   157,   201,   248,   202,    -1,   158,   201,
     248,   202,    -1,   159,   201,   248,   202,    -1,   160,   201,
     248,   202,    -1,   167,   201,   248,   202,    -1,    64,   201,
     248,   145,   248,   202,    -1,    72,   201,   248,   145,   248,
     202,    -1,   137,   201,   228,   145,   228,   145,   228,   202,
      -1,   203,    82,    33,   248,   200,    -1,    82,    -1,    82,
     201,   273,   202,    -1,   112,   200,   253,    31,    -1,    78,
     200,   253,    31,    -1,   253,   254,    -1,   254,    -1,   131,
      82,   200,   100,   255,   200,   130,   257,   200,    -1,   131,
      82,   200,   120,   228,   200,    -1,   131,    82,    33,   228,
     200,    -1,   131,    82,   145,    82,    33,   228,   200,    -1,
      22,    82,   145,    82,    33,   228,   200,    -1,   255,   145,
      51,    -1,   255,   145,    51,   204,    51,    -1,    51,   204,
      51,    -1,    51,    -1,   113,    33,   205,   258,   206,   200,
      -1,   257,   145,   228,    -1,   228,    -1,   258,   200,   259,
      -1,   259,    -1,   259,   145,   201,   228,   202,    -1,   259,
     145,    42,    -1,   259,   145,    51,    -1,   259,   201,   228,
     202,    -1,   259,    42,    -1,   259,    51,    -1,   201,   228,
     202,    -1,    42,    -1,    51,    -1,   121,   200,    -1,   121,
     201,   261,   202,   200,    -1,   261,   145,   262,    -1,   262,
      -1,   338,    -1,    19,   200,    -1,    19,   201,   264,   202,
     200,    -1,   264,   145,   265,    -1,   265,    -1,   338,    -1,
     114,   200,    -1,   114,   201,   267,   202,   200,    -1,   267,
     145,   268,    -1,   268,    -1,   351,    -1,   357,    -1,   122,
     200,    -1,   122,   201,   270,   202,   200,    -1,   122,   272,
     200,    -1,   122,   201,   270,   202,   272,   200,    -1,   270,
     145,   271,    -1,   271,    -1,   337,    -1,   338,    -1,   339,
      -1,   340,    -1,   341,    -1,   342,    -1,   343,    -1,   344,
      -1,   345,    -1,   346,    -1,   347,    -1,   364,    -1,   348,
      -1,   386,    -1,   349,    -1,   350,    -1,   351,    -1,   352,
      -1,   354,    -1,   355,    -1,   356,    -1,   390,    -1,   391,
      -1,   272,    82,    -1,   272,    82,    33,    82,    -1,   272,
     145,    82,    -1,   272,   145,    82,    33,    82,    -1,    82,
      -1,    82,    33,    82,    -1,   147,    51,    -1,   146,    51,
      -1,    51,    -1,   147,    42,    -1,   146,    42,    -1,    42,
      -1,    35,   200,   276,    31,    -1,   276,   277,    -1,   277,
      -1,   278,   145,   279,   200,    -1,   120,    82,    -1,    82,
      -1,    22,    82,   145,    82,    -1,   287,   145,   280,    -1,
     288,   145,   287,   145,   280,    -1,   288,   145,   288,   145,
     288,   145,   287,   145,   280,    -1,   288,    -1,   288,   145,
     288,   145,   288,    -1,   288,   145,   288,    -1,   288,   145,
     288,   145,   288,    -1,   288,   145,   288,   145,   288,   145,
     288,    -1,   288,   145,   288,   145,   288,   145,   288,   145,
     288,    -1,    37,   200,   282,    31,    -1,   282,   283,    -1,
     283,    -1,   120,    82,   145,   288,   200,    -1,    22,    82,
     145,    82,   145,   288,   200,    -1,    82,   145,   288,   200,
      -1,    36,   200,   285,    31,    -1,   285,   286,    -1,   286,
      -1,   120,    82,   145,   288,   145,   288,   200,    -1,    22,
      82,   145,    82,   145,   288,   145,   288,   200,    -1,    82,
     145,   288,   145,   288,   200,    -1,     6,    -1,    44,    -1,
      92,    -1,    52,    -1,   127,    -1,    -1,    51,    -1,    42,
      -1,    82,    -1,   146,    51,    -1,   146,    42,    -1,    34,
     200,    -1,    34,   201,   290,   202,   200,    -1,    34,   272,
     200,    -1,    34,   201,   290,   202,   272,   200,    -1,   290,
     145,   291,    -1,   291,    -1,   357,    -1,   358,    -1,   359,
      -1,   360,    -1,   361,    -1,   362,    -1,   363,    -1,   364,
      -1,   365,    -1,   366,    -1,   367,    -1,   368,    -1,   369,
      -1,   370,    -1,   371,    -1,   372,    -1,   373,    -1,   374,
      -1,   375,    -1,   376,    -1,   377,    -1,   378,    -1,   379,
      -1,   380,    -1,   348,    -1,   381,    -1,   382,    -1,   383,
      -1,   384,    -1,   385,    -1,   387,    -1,   388,    -1,   392,
      -1,   393,    -1,   394,    -1,   338,    -1,   395,    -1,   396,
      -1,   397,    -1,   106,   201,   293,   202,   200,    -1,   106,
     201,   293,   202,   272,   200,    -1,   293,   145,   294,    -1,
     294,    -1,   364,    -1,   365,    -1,   374,    -1,   380,    -1,
     348,    -1,   381,    -1,   382,    -1,   383,    -1,   384,    -1,
     385,    -1,   392,    -1,   393,    -1,   394,    -1,   107,   201,
     293,   202,   200,    -1,   107,   201,   293,   202,   272,   200,
      -1,   207,    82,   207,   145,   207,    82,   207,    -1,   207,
      82,   207,   145,   288,    -1,   296,    -1,   297,   145,   296,
      -1,   134,   272,   200,    -1,    93,   200,   300,    31,    -1,
     300,   301,    -1,   301,    -1,    82,   201,   228,   202,   200,
      -1,   128,   272,   200,    -1,    95,   200,   304,    31,    -1,
     304,    82,   228,   200,    -1,   304,    82,   145,    82,   228,
     200,    -1,    82,   228,   200,    -1,    82,   145,    82,   228,
     200,    -1,    98,   272,   200,    -1,    97,   200,    -1,    97,
     201,   271,   202,   200,    -1,    97,   272,   200,    -1,    97,
     201,   271,   202,   272,   200,    -1,    18,   200,   308,    31,
      -1,   308,   309,    -1,   309,    -1,    82,   310,    33,   228,
     200,    -1,    82,   145,    82,   310,    33,   228,   200,    -1,
       4,    82,   201,    51,   202,   310,    33,   228,   200,    -1,
      -1,   201,    51,   202,    -1,   201,    42,   202,    -1,    17,
     200,    -1,    17,   201,    23,   202,   200,    -1,    30,   201,
      82,   202,   200,    -1,    30,   201,    82,   202,   272,   200,
      -1,    30,    82,   200,    -1,    30,   201,    82,   208,    82,
     202,   200,    -1,    30,   201,    82,   208,    82,   202,   272,
     200,    -1,    30,    82,   208,    82,   200,    -1,    29,   201,
      82,   202,   200,    -1,    29,   201,    82,   202,   272,   200,
      -1,    29,    82,   200,    -1,    29,   201,    82,   208,    82,
     202,   200,    -1,    29,   201,    82,   208,    82,   202,   272,
     200,    -1,    29,    82,   208,    82,   200,    -1,    77,   201,
     315,   202,   317,   200,    -1,   315,   145,   316,    -1,   316,
      -1,   389,    -1,   390,    -1,   391,    -1,   318,    -1,   317,
     145,   318,    -1,   318,   201,   288,   202,    -1,   317,   145,
     318,   201,   288,   202,    -1,   319,    -1,   318,   319,    -1,
      82,    -1,   209,    -1,   148,    -1,   204,    -1,   208,    -1,
      -1,    -1,   101,   321,   248,   322,   200,    -1,   124,   200,
      -1,   124,   201,   324,   202,   200,    -1,   124,   272,   200,
      -1,   124,   201,   324,   202,   272,   200,    -1,   324,   145,
     325,    -1,   325,    -1,   271,    -1,   398,    -1,   399,    -1,
     400,    -1,   401,    -1,   402,    -1,   403,    -1,   404,    -1,
     405,    -1,   326,    -1,   357,    -1,   392,    -1,   393,    -1,
     359,    -1,   361,    -1,   358,    -1,   360,    -1,   395,    -1,
     396,    -1,   327,   145,   328,    -1,   327,    -1,     7,    51,
     200,    -1,     7,   201,   328,   202,    51,   200,    -1,   327,
      -1,   382,    -1,   365,    -1,   406,    -1,   330,   145,   331,
      -1,   330,    -1,     8,    51,   200,    -1,     8,   201,   331,
     202,    51,   200,    -1,   168,   200,    -1,   168,   201,   334,
     202,   200,    -1,   335,   145,   334,    -1,   335,    -1,   407,
      -1,   408,    -1,   409,    -1,   410,    -1,   411,    -1,   412,
      -1,   413,    -1,   414,    -1,   415,    -1,   416,    -1,   417,
      -1,   418,    -1,   419,    -1,   420,    -1,   421,    -1,   422,
      -1,   423,    -1,   424,    -1,   426,    -1,   427,    -1,   428,
      -1,   429,    -1,   430,    -1,   431,    -1,   432,    -1,   433,
      -1,   425,    -1,    51,    -1,    42,    -1,    26,    33,    51,
      -1,   118,    33,    51,    -1,   115,    33,    51,    -1,    60,
      -1,    96,    33,    51,    -1,   110,    33,    51,    -1,    27,
      33,    51,    -1,     3,    33,    51,    -1,    86,    -1,    88,
      -1,    90,    -1,    53,    33,    51,    -1,    48,    33,    51,
      -1,    49,    33,    51,    -1,   100,    33,    51,    -1,    24,
      33,   336,    -1,    63,    33,   336,    -1,   114,    -1,   116,
      33,    51,    -1,   108,    33,   336,    -1,    25,    33,    82,
      -1,    84,    33,   437,    -1,    84,    33,    51,    -1,    41,
      33,    51,    -1,   102,    33,    51,    -1,   103,    33,    51,
      -1,    58,    33,    51,    -1,    59,    33,    51,    -1,    89,
      -1,    46,    -1,    20,    33,   336,    -1,    70,    33,    51,
      -1,    65,    33,   336,    -1,    67,    33,   336,    -1,    94,
      33,   201,   297,   202,    -1,    66,    33,   336,    -1,    75,
      33,    82,    -1,    74,    33,    51,    -1,    73,    -1,   105,
      33,   336,    -1,    68,    33,    51,    -1,    69,    33,    51,
      -1,    61,    -1,    62,    -1,    87,    -1,     5,    -1,   123,
      -1,    43,    33,    51,    -1,   117,    -1,    81,    -1,    40,
      -1,   109,    -1,    54,    33,    51,    -1,    55,    33,    51,
      -1,    79,    33,    56,    -1,    79,    33,    80,    -1,   104,
      -1,    91,    -1,   135,    33,    82,    -1,   136,    33,   434,
      -1,    39,    33,   437,    -1,    21,    -1,    85,    -1,    71,
      -1,   125,    33,   336,    -1,    14,    33,   274,    -1,     9,
      33,   336,    -1,    11,    33,   274,    -1,    12,    33,   336,
      -1,    13,    33,    51,    -1,    10,    -1,    15,    33,    51,
      -1,    16,    33,    51,    -1,   169,    33,    51,    -1,   170,
      33,    51,    -1,   171,    33,    51,    -1,   172,    33,    51,
      -1,   173,    33,    51,    -1,   174,    33,    51,    -1,   175,
      33,    51,    -1,   176,    33,    51,    -1,   177,    33,    51,
      -1,   178,    33,    51,    -1,   179,    33,    51,    -1,   180,
      33,    51,    -1,   181,    33,    51,    -1,   182,    33,    51,
      -1,   183,    33,    51,    -1,   184,    33,   336,    -1,   185,
      33,   336,    -1,   186,    33,    51,    -1,   187,    33,   437,
      -1,   188,    33,   336,    -1,   189,    33,   336,    -1,   193,
      33,    51,    -1,   194,    33,    51,    -1,   196,    33,   336,
      -1,   197,    33,    51,    -1,   198,    33,   336,    -1,   199,
      33,   336,    -1,    82,   204,    82,    -1,    51,    -1,    51,
     204,    51,    -1,   205,   435,    -1,   436,   435,    -1,   436,
     206,    -1
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
     336,   341,   346,   351,   358,   365,   370,   379,   381,   385,
     390,   398,   402,   407,   410,   412,   417,   422,   425,   427,
     435,   439,   441,   443,   445,   447,   449,   451,   452,   458,
     459,   468,   469,   478,   479,   490,   491,   500,   503,   506,
     508,   510,   515,   518,   522,   524,   526,   528,   532,   536,
     540,   544,   548,   552,   556,   560,   564,   568,   572,   575,
     578,   583,   588,   593,   598,   603,   608,   613,   618,   623,
     628,   635,   642,   651,   657,   659,   664,   669,   674,   677,
     679,   689,   696,   702,   710,   718,   722,   728,   732,   734,
     741,   745,   747,   751,   753,   759,   763,   767,   772,   775,
     778,   782,   784,   786,   789,   795,   799,   801,   803,   806,
     812,   816,   818,   820,   823,   829,   833,   835,   837,   839,
     842,   848,   852,   859,   863,   865,   867,   869,   871,   873,
     875,   877,   879,   881,   883,   885,   887,   889,   891,   893,
     895,   897,   899,   901,   903,   905,   907,   909,   911,   914,
     919,   923,   929,   931,   935,   938,   941,   943,   946,   949,
     951,   956,   959,   961,   966,   969,   971,   976,   980,   986,
     996,   998,  1004,  1008,  1014,  1022,  1032,  1037,  1040,  1042,
    1048,  1056,  1061,  1066,  1069,  1071,  1079,  1089,  1096,  1098,
    1100,  1102,  1104,  1106,  1107,  1109,  1111,  1113,  1116,  1119,
    1122,  1128,  1132,  1139,  1143,  1145,  1147,  1149,  1151,  1153,
    1155,  1157,  1159,  1161,  1163,  1165,  1167,  1169,  1171,  1173,
    1175,  1177,  1179,  1181,  1183,  1185,  1187,  1189,  1191,  1193,
    1195,  1197,  1199,  1201,  1203,  1205,  1207,  1209,  1211,  1213,
    1215,  1217,  1219,  1221,  1223,  1229,  1236,  1240,  1242,  1244,
    1246,  1248,  1250,  1252,  1254,  1256,  1258,  1260,  1262,  1264,
    1266,  1268,  1274,  1281,  1289,  1295,  1297,  1301,  1305,  1310,
    1313,  1315,  1321,  1325,  1330,  1335,  1342,  1346,  1352,  1356,
    1359,  1365,  1369,  1376,  1381,  1384,  1386,  1392,  1400,  1410,
    1411,  1415,  1419,  1422,  1428,  1434,  1441,  1445,  1453,  1462,
    1468,  1474,  1481,  1485,  1493,  1502,  1508,  1515,  1519,  1521,
    1523,  1525,  1527,  1529,  1533,  1538,  1545,  1547,  1550,  1552,
    1554,  1556,  1558,  1560,  1561,  1562,  1568,  1571,  1577,  1581,
    1588,  1592,  1594,  1596,  1598,  1600,  1602,  1604,  1606,  1608,
    1610,  1612,  1614,  1616,  1618,  1620,  1622,  1624,  1626,  1628,
    1630,  1632,  1636,  1638,  1642,  1649,  1651,  1653,  1655,  1657,
    1661,  1663,  1667,  1674,  1677,  1683,  1687,  1689,  1691,  1693,
    1695,  1697,  1699,  1701,  1703,  1705,  1707,  1709,  1711,  1713,
    1715,  1717,  1719,  1721,  1723,  1725,  1727,  1729,  1731,  1733,
    1735,  1737,  1739,  1741,  1743,  1745,  1747,  1751,  1755,  1759,
    1761,  1765,  1769,  1773,  1777,  1779,  1781,  1783,  1787,  1791,
    1795,  1799,  1803,  1807,  1809,  1813,  1817,  1821,  1825,  1829,
    1833,  1837,  1841,  1845,  1849,  1851,  1853,  1857,  1861,  1865,
    1869,  1875,  1879,  1883,  1887,  1889,  1893,  1897,  1901,  1903,
    1905,  1907,  1909,  1911,  1915,  1917,  1919,  1921,  1923,  1927,
    1931,  1935,  1939,  1941,  1943,  1947,  1951,  1955,  1957,  1959,
    1961,  1965,  1969,  1973,  1977,  1981,  1985,  1987,  1991,  1995,
    1999,  2003,  2007,  2011,  2015,  2019,  2023,  2027,  2031,  2035,
    2039,  2043,  2047,  2051,  2055,  2059,  2063,  2067,  2071,  2075,
    2079,  2083,  2087,  2091,  2095,  2099,  2103,  2107,  2109,  2113,
    2116,  2119
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
     455,   457,   459,   461,   463,   467,   469,   471,   473,   478,
     481,   483,   487,   489,   493,   495,   497,   499,   501,   503,
     505,   507,   509,   513,   515,   519,   520,   523,   525,   527,
     531,   532,   535,   537,   539,   543,   544,   547,   548,   551,
     553,   555,   557,   561,   562,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,   582,   583,   584,   585,   586,   587,   590,   592,
     594,   596,   598,   600,   604,   606,   608,   612,   614,   616,
     620,   622,   624,   628,   630,   636,   642,   652,   657,   664,
     675,   680,   691,   698,   707,   718,   733,   736,   738,   742,
     750,   760,   770,   773,   775,   779,   789,   801,   813,   815,
     817,   819,   821,   825,   826,   827,   828,   829,   831,   835,
     837,   839,   841,   845,   846,   849,   850,   851,   852,   853,
     854,   855,   856,   857,   858,   859,   860,   861,   862,   863,
     864,   865,   866,   867,   868,   869,   870,   871,   872,   873,
     874,   875,   876,   877,   878,   879,   880,   881,   882,   883,
     884,   885,   886,   887,   890,   892,   896,   897,   900,   901,
     902,   903,   904,   905,   906,   907,   908,   909,   910,   911,
     912,   915,   917,   921,   923,   927,   928,   931,   933,   935,
     936,   939,   941,   943,   945,   947,   949,   951,   955,   957,
     959,   961,   963,   967,   969,   970,   973,   975,   977,   981,
     982,   984,   988,   990,   994,   996,   998,  1000,  1002,  1004,
    1008,  1010,  1012,  1014,  1016,  1018,  1022,  1025,  1026,  1029,
    1030,  1031,  1034,  1036,  1038,  1040,  1044,  1046,  1050,  1051,
    1053,  1055,  1057,  1061,  1062,  1061,  1064,  1066,  1068,  1070,
    1074,  1075,  1078,  1079,  1082,  1083,  1084,  1085,  1086,  1087,
    1088,  1091,  1092,  1093,  1094,  1095,  1096,  1097,  1098,  1099,
    1100,  1103,  1104,  1107,  1109,  1113,  1114,  1115,  1116,  1119,
    1120,  1123,  1125,  1129,  1131,  1135,  1136,  1139,  1140,  1141,
    1142,  1143,  1144,  1145,  1146,  1147,  1148,  1149,  1150,  1151,
    1152,  1153,  1154,  1155,  1156,  1157,  1158,  1159,  1160,  1161,
    1162,  1163,  1164,  1165,  1168,  1169,  1172,  1173,  1174,  1175,
    1176,  1177,  1178,  1179,  1180,  1181,  1182,  1183,  1184,  1185,
    1186,  1188,  1189,  1190,  1191,  1192,  1193,  1194,  1196,  1199,
    1200,  1201,  1202,  1203,  1204,  1206,  1209,  1210,  1211,  1212,
    1213,  1214,  1215,  1216,  1217,  1218,  1219,  1220,  1221,  1222,
    1223,  1224,  1225,  1226,  1227,  1228,  1229,  1230,  1231,  1232,
    1233,  1235,  1238,  1239,  1240,  1241,  1242,  1243,  1244,  1245,
    1246,  1248,  1249,  1250,  1251,  1252,  1253,  1254,  1255,  1257,
    1258,  1259,  1260,  1261,  1262,  1263,  1264,  1265,  1266,  1267,
    1268,  1269,  1270,  1271,  1272,  1273,  1274,  1275,  1277,  1278,
    1284,  1285,  1289,  1290,  1291,  1292,  1295,  1303,  1304,  1313,
    1315,  1324
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
  const int parser::yylast_ = 2379;
  const int parser::yynnts_ = 228;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 160;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 210;

  const unsigned int parser::yyuser_token_number_max_ = 454;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1326 "DynareBison.yy"


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

