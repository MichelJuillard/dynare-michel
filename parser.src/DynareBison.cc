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
	  case 49:
#line 154 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 50:
#line 156 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 51:
#line 159 "DynareBison.yy"
    { driver.rplot(); ;}
    break;

  case 56:
#line 170 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 57:
#line 172 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 58:
#line 174 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 59:
#line 176 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 60:
#line 178 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 61:
#line 180 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 62:
#line 184 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 63:
#line 186 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 64:
#line 188 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 65:
#line 190 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 66:
#line 192 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 67:
#line 194 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 68:
#line 198 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 69:
#line 200 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 70:
#line 202 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 71:
#line 204 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 72:
#line 206 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 73:
#line 208 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 74:
#line 212 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 75:
#line 214 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 76:
#line 216 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 77:
#line 218 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 78:
#line 220 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 79:
#line 222 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 80:
#line 226 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 81:
#line 228 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 82:
#line 232 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 83:
#line 234 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 84:
#line 238 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 85:
#line 240 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 86:
#line 244 "DynareBison.yy"
    { driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 87:
#line 247 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 88:
#line 249 "DynareBison.yy"
    { (yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 89:
#line 251 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 90:
#line 253 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 91:
#line 255 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 92:
#line 257 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 93:
#line 259 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 94:
#line 261 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 95:
#line 263 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 96:
#line 265 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 97:
#line 267 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 98:
#line 269 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 99:
#line 271 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 100:
#line 273 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 101:
#line 275 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 102:
#line 277 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 103:
#line 279 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 104:
#line 281 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 105:
#line 283 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 106:
#line 285 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 107:
#line 287 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 108:
#line 289 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 109:
#line 291 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 110:
#line 293 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 111:
#line 295 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 112:
#line 297 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 113:
#line 299 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
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
#line 307 "DynareBison.yy"
    { (yyval.node_val) = driver.add_normcdf((yysemantic_stack_[(8) - (3)].node_val),(yysemantic_stack_[(8) - (5)].node_val),(yysemantic_stack_[(8) - (7)].node_val));;}
    break;

  case 118:
#line 311 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 119:
#line 313 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 120:
#line 317 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 121:
#line 319 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 122:
#line 322 "DynareBison.yy"
    { driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 123:
#line 324 "DynareBison.yy"
    { driver.end_endval(); ;}
    break;

  case 126:
#line 330 "DynareBison.yy"
    { driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 127:
#line 332 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 130:
#line 338 "DynareBison.yy"
    { driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 133:
#line 345 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 134:
#line 347 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 135:
#line 349 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 138:
#line 354 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 139:
#line 355 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 140:
#line 356 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 141:
#line 357 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 142:
#line 358 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 143:
#line 359 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 144:
#line 361 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 145:
#line 362 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 146:
#line 363 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 147:
#line 364 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 152:
#line 374 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 153:
#line 376 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val)); ;}
    break;

  case 154:
#line 380 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 156:
#line 383 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 157:
#line 385 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 158:
#line 387 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 159:
#line 389 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 160:
#line 391 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 161:
#line 393 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 162:
#line 395 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 163:
#line 397 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 164:
#line 399 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 165:
#line 401 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 166:
#line 403 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 167:
#line 405 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 168:
#line 407 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 169:
#line 409 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 170:
#line 411 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 171:
#line 413 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 172:
#line 415 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 173:
#line 417 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 174:
#line 419 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 175:
#line 421 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 176:
#line 423 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 177:
#line 425 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 178:
#line 427 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 179:
#line 429 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 180:
#line 431 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
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
#line 437 "DynareBison.yy"
    { (yyval.node_val) = driver.add_normcdf((yysemantic_stack_[(8) - (3)].node_val),(yysemantic_stack_[(8) - (5)].node_val),(yysemantic_stack_[(8) - (7)].node_val));;}
    break;

  case 184:
#line 441 "DynareBison.yy"
    { driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 185:
#line 444 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 186:
#line 446 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 187:
#line 449 "DynareBison.yy"
    { driver.end_shocks(); ;}
    break;

  case 188:
#line 451 "DynareBison.yy"
    { driver.end_mshocks(); ;}
    break;

  case 191:
#line 458 "DynareBison.yy"
    { driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val)); ;}
    break;

  case 192:
#line 460 "DynareBison.yy"
    { driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 193:
#line 462 "DynareBison.yy"
    { driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 194:
#line 464 "DynareBison.yy"
    { driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 195:
#line 466 "DynareBison.yy"
    { driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 196:
#line 470 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 197:
#line 472 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 198:
#line 474 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 199:
#line 476 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 200:
#line 480 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 201:
#line 484 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 202:
#line 486 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].node_val));;}
    break;

  case 203:
#line 490 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 204:
#line 492 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 205:
#line 496 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 206:
#line 498 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 207:
#line 500 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 208:
#line 502 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 209:
#line 504 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 210:
#line 506 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 211:
#line 508 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 212:
#line 510 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 213:
#line 512 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 214:
#line 516 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 215:
#line 518 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 221:
#line 531 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 222:
#line 533 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 226:
#line 543 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 227:
#line 545 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 232:
#line 557 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 233:
#line 559 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 234:
#line 561 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 235:
#line 563 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 261:
#line 596 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 262:
#line 598 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 263:
#line 600 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 264:
#line 602 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 265:
#line 604 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 266:
#line 606 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 267:
#line 610 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 268:
#line 612 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 269:
#line 614 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 270:
#line 618 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 271:
#line 620 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 272:
#line 622 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 273:
#line 625 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 274:
#line 628 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 275:
#line 630 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 277:
#line 636 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 278:
#line 642 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 279:
#line 648 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 280:
#line 658 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 281:
#line 663 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 282:
#line 670 "DynareBison.yy"
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

  case 283:
#line 681 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 284:
#line 686 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 285:
#line 697 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 286:
#line 704 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 287:
#line 713 "DynareBison.yy"
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

  case 288:
#line 724 "DynareBison.yy"
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

  case 289:
#line 739 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 290:
#line 742 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 291:
#line 744 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 292:
#line 748 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 293:
#line 756 "DynareBison.yy"
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

  case 294:
#line 766 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 295:
#line 776 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 296:
#line 779 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 297:
#line 781 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 298:
#line 785 "DynareBison.yy"
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

  case 299:
#line 795 "DynareBison.yy"
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

  case 300:
#line 807 "DynareBison.yy"
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

  case 301:
#line 819 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 302:
#line 821 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 303:
#line 823 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 304:
#line 825 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 305:
#line 827 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 306:
#line 830 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 310:
#line 835 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 311:
#line 837 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 312:
#line 841 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 313:
#line 843 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 314:
#line 845 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 315:
#line 847 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 357:
#line 896 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 358:
#line 898 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 374:
#line 921 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 375:
#line 923 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 376:
#line 927 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 377:
#line 929 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 380:
#line 936 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 381:
#line 938 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 384:
#line 944 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 385:
#line 946 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 386:
#line 948 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 387:
#line 951 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 388:
#line 953 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 389:
#line 955 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 390:
#line 957 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 391:
#line 960 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 392:
#line 963 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 393:
#line 965 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 394:
#line 967 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 395:
#line 969 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 396:
#line 972 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 399:
#line 979 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 400:
#line 981 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 401:
#line 983 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 402:
#line 986 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 403:
#line 988 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 404:
#line 990 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 405:
#line 994 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 406:
#line 996 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 407:
#line 1000 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 408:
#line 1002 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 409:
#line 1004 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 410:
#line 1006 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 411:
#line 1008 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 412:
#line 1010 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 413:
#line 1014 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 414:
#line 1016 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 415:
#line 1018 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 416:
#line 1020 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 417:
#line 1022 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 418:
#line 1024 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 419:
#line 1028 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 425:
#line 1040 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 426:
#line 1042 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 427:
#line 1044 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 428:
#line 1046 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 429:
#line 1050 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 430:
#line 1052 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 432:
#line 1057 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 433:
#line 1059 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 434:
#line 1061 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 435:
#line 1063 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 436:
#line 1066 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 437:
#line 1067 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 439:
#line 1070 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 440:
#line 1072 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 441:
#line 1074 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 442:
#line 1076 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 466:
#line 1113 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 467:
#line 1115 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 474:
#line 1129 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 475:
#line 1131 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 476:
#line 1135 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 477:
#line 1137 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 507:
#line 1174 "DynareBison.yy"
    { driver.end_homotopy();;}
    break;

  case 510:
#line 1181 "DynareBison.yy"
    { driver.homotopy_val((yysemantic_stack_[(6) - (1)].string_val),(yysemantic_stack_[(6) - (3)].node_val),(yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 513:
#line 1187 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 514:
#line 1188 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 515:
#line 1189 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 516:
#line 1190 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 517:
#line 1191 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 518:
#line 1192 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 519:
#line 1193 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 520:
#line 1194 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 521:
#line 1195 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 522:
#line 1196 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 523:
#line 1197 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 524:
#line 1198 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 525:
#line 1199 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 526:
#line 1200 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 527:
#line 1202 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 528:
#line 1203 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1204 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1205 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 531:
#line 1206 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 532:
#line 1207 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 533:
#line 1208 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 534:
#line 1210 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 535:
#line 1212 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 536:
#line 1214 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 537:
#line 1215 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 538:
#line 1216 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1217 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 540:
#line 1218 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1220 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 542:
#line 1222 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 543:
#line 1224 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 544:
#line 1225 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 545:
#line 1226 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 546:
#line 1227 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 548:
#line 1229 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 549:
#line 1230 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 550:
#line 1231 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 551:
#line 1232 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 552:
#line 1233 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 553:
#line 1234 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 554:
#line 1235 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 555:
#line 1236 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 556:
#line 1237 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 557:
#line 1238 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 558:
#line 1239 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 559:
#line 1240 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 560:
#line 1241 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 561:
#line 1242 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 562:
#line 1243 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 563:
#line 1244 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 564:
#line 1245 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 565:
#line 1246 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 566:
#line 1247 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1249 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 568:
#line 1251 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 569:
#line 1253 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 570:
#line 1254 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 571:
#line 1255 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1256 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1257 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1258 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 575:
#line 1259 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 576:
#line 1260 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 577:
#line 1261 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1263 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 579:
#line 1264 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1265 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 581:
#line 1266 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 582:
#line 1267 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1268 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 584:
#line 1269 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1270 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1272 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1273 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1274 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1275 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1276 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1277 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 592:
#line 1278 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 593:
#line 1279 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 594:
#line 1280 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 595:
#line 1281 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 596:
#line 1282 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 597:
#line 1283 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 598:
#line 1284 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 599:
#line 1285 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 600:
#line 1286 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 601:
#line 1287 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 602:
#line 1288 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 603:
#line 1289 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 604:
#line 1290 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 605:
#line 1292 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 606:
#line 1293 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 607:
#line 1299 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 608:
#line 1300 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 609:
#line 1304 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 610:
#line 1305 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 611:
#line 1306 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 612:
#line 1307 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 613:
#line 1309 "DynareBison.yy"
    {driver.option_num("homotopy_mode",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 614:
#line 1310 "DynareBison.yy"
    {driver.option_num("homotopy_steps",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 615:
#line 1313 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 617:
#line 1322 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); 
               ;}
    break;

  case 618:
#line 1331 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 619:
#line 1333 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 620:
#line 1341 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2429 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1086;
  const short int
  parser::yypact_[] =
  {
       978,    10,    36,   -73,  -102,   411,   251,    62,     7,    38,
      78,    -1,    96,   146,   150,   245,   429,   425,   492,   -92,
     250,   112,   261,   276,     4,    65,   131,    89, -1086,   -86,
      59,    65,   284,   312,   495,   497,    72,    74,    65,   418,
     458,   470,    65,   344,   500,   838, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086,   403,  1302,   435,   806, -1086,   588,   100, -1086,
     522,   601,   453,    33,   -62,   567,   -52,   585,   595,   659,
   -1086,  1279,    45,    76,   402,   408,   623,   595,   668,   666,
     506, -1086,   349,    42,    40,  1103,   628,   629, -1086,  1457,
      57,    79,   587,    80,   663,   512,  1132,  1319,  1319,    97,
      40,   508, -1086,    94, -1086,   491, -1086,  1457,   198, -1086,
    1386,   242,   258,   591,   260,   592,   266,   597,   268,   270,
     652, -1086,  2003, -1086, -1086, -1086,   702, -1086,   704,   705,
     708,   709,   710, -1086,   712,   713,   714, -1086,   715,   716,
     717,   718, -1086,   604,   548, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086,   722,   723,   724, -1086,   610,   554, -1086, -1086,
   -1086,   559,   684,   -49,    92, -1086,   735,   -55, -1086, -1086,
     566, -1086,   568, -1086, -1086,   690,   227, -1086,   691,   328,
     742,    84, -1086,   695, -1086,   748, -1086, -1086,   749,   750,
     751,   761,   762, -1086, -1086,   763,   764,   765,   766,   768,
     770, -1086, -1086,   771,   772, -1086, -1086, -1086,   774,   775,
   -1086, -1086,   -32, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086,   776,   729, -1086,   730, -1086,   731,    86,
   -1086,   675,   732,   676,   743,   421, -1086,   747,   682,   753,
     577, -1086,   633,   400, -1086,   443,   799,   634,   638, -1086,
     630, -1086,   -10,   639,   643,   800, -1086, -1086,   142, -1086,
   -1086, -1086, -1086,   769,   779,   122, -1086, -1086, -1086,   646,
     648,   649,   650,  1103,  1103,   654,   656,   660,   661,   665,
     667,   672,   673,   674,   677,  1103,  1019,   678,   463, -1086,
     976,   571,   810,   830,   846,   847,   850,   851, -1086, -1086,
   -1086,   853,   854,   859, -1086,   860, -1086,   861,   862,   692,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086,   773,   816, -1086,   697,
   -1086, -1086, -1086,   698,   699,   700,   701,  1132,  1132,   703,
     706,   707,   736,   739,   744,   752,   754,   757,   759,  1132,
    1990, -1086,   244, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086,   288, -1086,   404,
      29,   296, -1086, -1086, -1086,   873,   880,   303, -1086, -1086,
   -1086, -1086,   314, -1086, -1086,   884, -1086,   318, -1086, -1086,
   -1086, -1086, -1086,   792,   844, -1086, -1086,   820,   865, -1086,
   -1086,   827,   875, -1086, -1086,   807,   576, -1086,   931,   932,
     934,   935,   940,   941,   942,   944,   945,   946,   947,   948,
     949,   951,   954,   955,   956,   957,   960,   965,   966,   967,
     968,   970,   971,   972,   983,   812,   863, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086,   465,   297,   465,   969,   297,   973,
     937,   975,    31,   979,   980,   939,   950,  1302,   982,   984,
     465,   985,   806,   986,   819,   825,   952,   478,   990, -1086,
   -1086,   987,   522,   836, -1086, -1086,   839,    54,   961,   841,
      66,   963,  1103, -1086, -1086, -1086,   843,   995,   996,   998,
    1002,  1006,   465,   465,   465,  1015,  1023,  1029,  1030,  1000,
     879,   465,  1279,    82,  1004,  1054,   953, -1086, -1086, -1086,
     368,   959,    52,   974, -1086, -1086,   981,    52,   993, -1086,
   -1086,   290, -1086, -1086, -1086,  1011,   885, -1086,  1012,    28,
   -1086,   253, -1086,   441, -1086,   892,   893,   446,    42,    61,
     994,    44, -1086, -1086,  1103,  1103,  1103,  1103,   943,   466,
    1103,  1103,  1103,  1103,  1103,  1103,  1103,  1103,  1103,  1103,
     519,  1103,  1103,  1103,  1103,  1103,  1103,  1103,  1103,  1103,
    1103,  1103, -1086,  1103, -1086, -1086,  1016,  1276, -1086,  1076,
    1053,   465,  1057,  1063,  1064,  1066,  1068,  1069,   465,  1070,
    1077,  1079,    85, -1086,  1017, -1086,  1132,  1132,   290,  1132,
     997,   524,  1132,  1132,  1132,  1132,  1132,  1132,  1132,  1132,
    1132,  1132,   539,  1132,  1132,  1132,  1132,  1132,  1132,  1132,
    1132,  1132,  1132,  1132,   902,  1319,    93,   222, -1086, -1086,
   -1086,  1103,   -74,    21,    94,   958,  1093,  1096,   491,   988,
    1457,   225,   465,  1386,   230, -1086,  1024, -1086,  1026, -1086,
    1027,  1103, -1086, -1086,  1104,  1105,  1106,  1108,  1121,  1125,
    1126,  1127,  1128,  1130,  1131,  1147,  1148,  1149,  1150,   465,
     465,  1151,   843,   465,   465,  1152,  1154,   465,  1155,   465,
     465,  1005,  2003, -1086, -1086, -1086, -1086,  1165,  1167, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086,  1159,    17, -1086,
   -1086, -1086, -1086,  1008, -1086, -1086,  1009, -1086, -1086, -1086,
   -1086,  1013, -1086,  1160,  1014,  1018,  1022,  1103, -1086, -1086,
   -1086, -1086, -1086,   273,  1025, -1086, -1086,   289,  1036,  1789,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086,  1007, -1086, -1086, -1086,   293, -1086,
    1137,  1138, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
     562,  1039,  1073,  1080,  1161,  1081,    52,  1162,  1042,    52,
   -1086,  1196,  1198,  1045, -1086,   595,  1218, -1086, -1086, -1086,
    1132, -1086, -1086, -1086,  1221, -1086,   327, -1086, -1086, -1086,
    1052, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086,   -23,    90, -1086,  1174,  1103,  1175,    37,  2065,  2077,
    2002,   329,  2089,   778,   918,  1043,  1375,  1387,  1444,  1456,
    1468,  1480,  1492, -1086,   473,   473,   473,   473,   473,   473,
     466,   466,   943,   943,  1113,  1513,  1103, -1086,  1186,  1801,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086,   294, -1086,  2101,  2113,  1065,  2125,  1525,
    1537,  1549,  1561,  1582,  1594,  1606,  1618,  1630,  1651, -1086,
     496,   496,   496,   496,   496,   496,   524,   524,   997,   997,
    1113, -1086, -1086, -1086,   306, -1086,   307,  1663,    29,  1071,
   -1086, -1086,    46,  1103, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086,   309, -1086, -1086, -1086,   335, -1086, -1086, -1086,
    2137, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086,  1072, -1086, -1086, -1086,  1189, -1086, -1086,  1067,
    1242, -1086, -1086,  1813, -1086,   231, -1086,   256, -1086,  1194,
   -1086,   332, -1086, -1086, -1086, -1086, -1086, -1086,    52,   368,
    1129,    52,  1135,  1153, -1086,  1075, -1086, -1086,  1252,   579,
    1132,  1825,   465,   441, -1086,   630,   630,   630,    61, -1086,
      52, -1086,  1253,  1837,  1263,  1246,  1103,  1103,  1103,  1103,
   -1086,  1103, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086,  1095,  1856,  1103, -1086, -1086,  1132,  1132, -1086,
    1132, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086,    21, -1086, -1086, -1086,  1103,  1675,
   -1086, -1086,  1103,  1254, -1086,  1014,  1103, -1086, -1086,   340,
   -1086,   341,  1098,  1007, -1086, -1086,  1158,  1173,  1178,    52,
    1100,    52,    52, -1086,  1103, -1086,  1868, -1086, -1086, -1086,
    1107,    60,   207,   223,   147,  1123,  1103, -1086,  1103,  1102,
     -22,  1880,  1687,  1699,  2002,  2149, -1086, -1086,  1892,  1720,
    1732,  2161,  1744, -1086,  1904, -1086,  1296,  1923, -1086, -1086,
    1182, -1086,    52,    52,    52,  1183, -1086,  1139,  1164,  1935,
   -1086,   630, -1086, -1086, -1086,    52, -1086,  1947,  1959,  1284,
    1300,  1225, -1086, -1086, -1086,  1103, -1086, -1086, -1086,  1132,
   -1086, -1086,  1103, -1086,    27,  1208, -1086,  1209,    52, -1086,
   -1086, -1086,   570,  1156, -1086, -1086, -1086,  1163,  1103,  1756,
    1768,  1971,  1287, -1086,    52,   108,  1168, -1086, -1086,  1323,
    2002,   158, -1086, -1086, -1086,  1166,  1227,  1229, -1086, -1086,
    1103, -1086, -1086,    52,    52,  2002,  1230, -1086,    52, -1086
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   436,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     2,     4,    29,    30,
      46,    47,    48,    45,     5,     6,     7,    12,     9,    10,
      11,     8,    13,    14,    15,    16,    17,    18,    19,    23,
      25,    24,    20,    21,    22,    26,    27,    28,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,     0,     0,     0,     0,   405,     0,     0,   221,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   265,
     312,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   138,     0,     0,     0,     0,     0,     0,   392,     0,
       0,     0,    76,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   226,     0,   214,     0,   232,     0,     0,   439,
       0,     0,     0,    58,     0,    64,     0,    70,     0,     0,
       0,   476,     0,     1,     3,   466,     0,   583,     0,     0,
       0,     0,     0,   574,     0,     0,     0,   575,     0,     0,
       0,     0,   454,   465,     0,   455,   460,   458,   461,   459,
     456,   457,   462,   463,   447,   448,   449,   450,   451,   452,
     453,   474,     0,     0,     0,   468,   473,     0,   470,   469,
     471,     0,     0,   402,     0,   398,     0,     0,   224,   225,
       0,    82,     0,    49,   415,     0,     0,   409,     0,     0,
       0,     0,   125,     0,   558,     0,   563,   542,     0,     0,
       0,     0,     0,   555,   556,     0,     0,     0,     0,     0,
       0,   576,   551,     0,     0,   562,   557,   541,     0,     0,
     561,   559,     0,   317,   353,   342,   318,   319,   320,   321,
     322,   323,   324,   325,   326,   327,   328,   329,   330,   331,
     332,   333,   334,   335,   336,   337,   338,   339,   340,   341,
     343,   344,   345,   346,   347,   348,   349,   350,   351,   352,
     354,   355,   356,   261,     0,   314,     0,   278,     0,     0,
     275,     0,     0,     0,     0,     0,   297,     0,     0,     0,
       0,   291,     0,     0,   129,     0,     0,     0,     0,    84,
       0,   516,     0,     0,     0,     0,   570,   569,     0,   421,
     422,   423,   424,     0,     0,     0,   190,    89,    90,     0,
       0,    88,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   383,
       0,     0,     0,     0,     0,     0,     0,     0,   521,   522,
     523,     0,     0,     0,   564,     0,   530,     0,     0,     0,
     238,   239,   240,   241,   242,   243,   244,   245,   246,   247,
     248,   250,   252,   253,   254,   255,   256,   257,   258,   249,
     251,   259,   260,   394,   391,    79,    74,     0,    55,     0,
      80,   156,   157,     0,     0,   185,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     437,   155,     0,   360,   365,   361,   362,   363,   364,   366,
     367,   368,   369,   370,   371,   372,   373,     0,    51,     0,
       0,     0,   229,   230,   231,     0,     0,     0,   217,   218,
     219,   220,     0,   237,   234,     0,   445,     0,   444,   446,
     441,   385,    61,    56,     0,    52,    67,    62,     0,    53,
      73,    68,     0,    54,   380,     0,     0,   508,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   479,   480,   481,   482,
     483,   484,   485,   486,   487,   488,   489,   490,   491,   492,
     493,   494,   495,   496,   497,   506,   498,   499,   500,   501,
     502,   503,   504,   505,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   396,
     397,     0,     0,     0,    83,    50,     0,     0,     0,     0,
       0,     0,     0,   123,   124,   266,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   263,     0,   277,   273,   274,
     306,     0,   306,     0,   295,   296,     0,   306,     0,   289,
     290,     0,   127,   128,   120,     0,     0,    85,     0,     0,
     150,     0,   151,     0,   146,     0,     0,     0,     0,     0,
       0,     0,   188,   189,     0,     0,     0,     0,   102,   103,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    86,     0,   381,   382,     0,     0,   386,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    77,    75,    81,     0,     0,     0,     0,
     169,   170,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   187,   212,
     213,     0,     0,   204,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    59,    57,    65,    63,    71,
      69,     0,   507,   509,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   512,   511,   579,   272,     0,     0,   580,
     581,   582,   578,   584,   533,   536,   535,     0,     0,   534,
     537,   538,   571,     0,   572,   464,     0,   585,   543,   560,
     472,     0,   406,     0,   402,     0,     0,     0,   514,   223,
     222,   418,   413,     0,     0,   412,   407,     0,     0,     0,
     573,   524,   565,   566,   539,   540,   545,   548,   546,   553,
     554,   544,   550,   549,     0,   552,   316,   313,     0,   262,
       0,     0,   301,   308,   302,   307,   304,   309,   303,   305,
       0,     0,     0,   283,     0,     0,   306,     0,     0,   306,
     269,     0,     0,     0,   122,     0,     0,   139,   148,   149,
       0,   153,   134,   133,     0,   135,     0,   132,   136,   137,
       0,   142,   140,   567,   568,   420,   431,   433,   434,   435,
     432,     0,   425,   429,     0,     0,     0,     0,     0,     0,
     118,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    87,   101,   100,    99,    98,    97,    96,
      92,    91,    93,    94,    95,     0,     0,   389,     0,     0,
     520,   528,   513,   519,   525,   526,   517,   527,   532,   518,
     515,   531,   393,     0,    78,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   154,
     167,   166,   165,   164,   163,   162,   159,   158,   160,   161,
     168,   438,   359,   357,     0,   374,     0,     0,     0,     0,
     209,   210,     0,     0,   228,   227,   613,   614,   216,   215,
     236,   233,     0,   577,   443,   440,     0,    60,    66,    72,
       0,   586,   587,   588,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   599,   600,   601,   602,   603,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   477,   478,
     271,   270,   616,   618,   620,   619,     0,   467,   475,     0,
       0,   404,   403,     0,   414,     0,   408,     0,   126,     0,
     378,     0,   315,   264,   279,   311,   310,   276,   306,   306,
       0,   306,     0,     0,   294,     0,   268,   267,     0,     0,
       0,     0,     0,     0,   144,     0,     0,     0,     0,   419,
     306,   430,     0,     0,     0,     0,     0,     0,     0,     0,
     116,     0,   104,   105,   106,   107,   108,   109,   110,   111,
     112,   113,     0,     0,     0,   387,   395,     0,     0,   186,
       0,   171,   172,   173,   174,   175,   176,   177,   178,   179,
     180,   358,   375,   211,   203,   200,   206,   207,     0,     0,
     235,   442,     0,     0,   615,   402,     0,   399,   416,     0,
     410,     0,     0,     0,   547,   280,     0,     0,     0,   306,
       0,   306,   306,   292,     0,   121,     0,   152,   529,   131,
       0,     0,     0,     0,   426,     0,     0,   193,     0,   199,
       0,     0,     0,     0,   119,     0,   384,   390,     0,     0,
       0,     0,     0,   208,     0,   617,     0,     0,   417,   411,
       0,   379,   306,   306,   306,     0,   300,     0,     0,     0,
     184,     0,   147,   143,   141,   306,   427,     0,     0,     0,
       0,     0,   192,   114,   115,     0,   388,   181,   182,     0,
     205,   510,     0,   400,   306,   285,   281,   284,   306,   298,
     293,   130,     0,     0,   195,   194,   198,   196,     0,     0,
       0,     0,     0,   377,   306,     0,     0,   145,   428,     0,
     202,     0,   117,   183,   401,     0,   286,     0,   299,   197,
       0,   191,   376,   306,   306,   201,   287,   282,   306,   288
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1086, -1086,  1334, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,  -346, -1086,
   -1086, -1086, -1086,  -115,  -230, -1086, -1086,  1060, -1086,   287,
   -1086, -1086, -1086, -1086, -1086, -1086,  -990,  -633,  -136,  -624,
   -1086, -1086, -1086,  1245,  -255, -1086, -1086, -1086, -1086,   390,
   -1086, -1086,   642, -1086, -1086,   809, -1086, -1086,   651, -1086,
   -1086,   -96,   -15,   685,   834, -1086, -1086,  1085, -1086, -1086,
   -1085, -1086, -1086,  1078, -1086, -1086,  1084, -1022,  -609, -1086,
   -1086,   789, -1086,  1265,   671, -1086,   246, -1086, -1086, -1086,
   -1086,  1048, -1086, -1086, -1086, -1086, -1086, -1086, -1086,  1193,
    -768, -1086, -1086, -1086, -1086, -1086,   781, -1086,   313,  -848,
   -1086, -1086, -1086, -1086, -1086,   680, -1086,   -41,   864, -1086,
   -1086,   868, -1086, -1086,   626, -1086, -1086, -1086,   962,  -528,
   -1086,   -97, -1086,  1321, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086,  -106, -1086, -1086,  -118,  -594, -1086, -1086, -1086, -1086,
    -107,   -85,   -82,   -76,   -70, -1086, -1086,   -99,   -61, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086,   -66, -1086, -1086,
   -1086, -1086, -1086,   -59,   -56,   -58,   -51,   -46,   -45, -1086,
   -1086, -1086, -1086,   -93,   -83,   -94,   -91,   -44,   -38,   -26,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
   -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086, -1086,
     636, -1086,  -538
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    45,    46,    47,    48,    49,    50,    51,    52,    53,
     154,   156,   158,   133,    54,    55,    56,    57,   366,   921,
      58,   327,    59,   231,   232,    60,   323,   324,   896,   897,
      61,   330,  1097,  1096,  1180,   900,   639,   640,   641,   642,
     441,    62,    63,   345,   346,  1190,    64,  1271,   742,   743,
      65,   467,   468,    66,   217,   218,    67,   461,   462,    68,
     472,   476,   112,   883,   799,    69,   309,   310,   311,   871,
    1165,    70,   320,   321,    71,   315,   316,   872,  1166,    72,
     262,   263,    73,   442,   443,    74,  1070,  1071,    75,    76,
     368,   369,    77,    78,   371,    79,    80,    81,   214,   215,
     578,    82,    83,    84,    85,   338,   339,   911,   912,   913,
      86,   136,   734,    87,   477,   478,   182,   183,   184,    88,
     206,   207,    89,    90,   525,   526,    91,   496,   497,   795,
     390,   391,   392,   393,   394,   395,   396,   397,   398,   399,
     400,   401,   402,   403,   404,   405,   899,   406,   407,   408,
     185,   186,   187,   188,   189,   271,   272,   409,   446,   275,
     276,   277,   278,   279,   280,   281,   282,   447,   284,   285,
     286,   287,   288,   448,   449,   450,   451,   452,   453,   410,
     295,   296,   340,   411,   412,   190,   191,   456,   192,   193,
     302,   479,   194,   195,   196,   197,   198,   199,   200,   210,
     527,   528,   529,   530,   531,   532,   533,   534,   535,   536,
     537,   538,   539,   540,   541,   542,   543,   544,   545,   546,
     547,   548,   549,   550,   551,   552,   553,   470,   471,   814,
    1053,   808,   809
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       440,   594,   325,   219,   266,   265,   888,   658,   659,   130,
     131,   873,   273,   875,   264,   889,   139,   297,   878,   670,
     298,   148,   151,   152,   687,   463,   267,   159,   800,   268,
     341,   444,   444,   389,   208,   269,   464,   209,   445,   445,
     342,   270,   818,   454,   454,   283,   455,   455,   469,   898,
     274,   473,   289,   291,   205,   290,  1060,  1167,   840,   887,
     292,    92,   343,  1000,  1101,   293,   294,   299,  1052,   863,
     421,   739,  1001,   300,   846,   847,   848,   915,   865,   422,
     740,   109,   806,   855,   222,   301,   109,    94,  1146,   104,
     653,  1222,   423,   582,   863,   594,   212,  1147,   306,   576,
     424,    98,   421,   865,   212,  1181,  1182,  1183,   306,   867,
     425,   422,   123,   103,   862,   593,   612,   618,   137,   174,
     106,   335,   134,   579,   423,  1098,  1230,   303,  1246,   998,
      96,    97,   424,   336,   867,   999,   109,  1105,   643,   303,
     135,   224,   425,   906,   343,   125,   337,   109,   109,   225,
     583,   227,   864,   652,   109,   577,   109,  1106,   307,   228,
     866,   303,   416,   951,   109,   426,   230,   109,   307,  1002,
     958,   344,   906,   613,   213,   109,   870,   427,   428,   303,
    1099,  1231,   213,   429,   430,   431,   432,   433,   434,   435,
     436,   437,   916,   304,   382,   644,   308,   426,   438,  1287,
     868,   870,   110,   111,   653,   304,   308,   128,   129,   427,
     428,   105,   907,   132,    93,   429,   430,   431,   432,   433,
     434,   435,   436,   437,  1013,  1003,  1054,   304,   417,   906,
     438,  1252,   439,   741,   638,   869,   223,  1262,  1223,   807,
      95,   907,   107,  1277,  1039,   304,   839,   917,   305,   421,
    1148,  1036,  1037,   344,  1224,  1040,  1041,   832,   422,  1044,
     413,  1046,  1047,   138,   439,   421,   638,  1082,   908,   836,
    1085,   423,   909,   910,   422,   146,   147,   149,   150,   424,
     303,   108,   414,   418,   101,   857,   890,   423,   962,   425,
     648,   710,   711,   102,  1100,   424,   993,   908,   907,   113,
     458,   909,   910,   722,   109,   425,  1280,   109,   918,   919,
     920,   922,   109,   109,   923,   924,   925,   926,   927,   928,
     929,   930,   931,   932,   303,   934,   935,   936,   937,   938,
     939,   940,   941,   942,   943,   944,  1101,   945,   109,   796,
     303,   880,   483,   949,   426,   141,   304,   649,   487,   114,
     491,  1225,   303,   115,   908,   303,   427,   428,   909,   910,
     426,  1281,   429,   430,   431,   432,   433,   434,   435,   436,
     437,   303,   427,   428,   862,   303,   303,   438,   429,   430,
     431,   432,   433,   434,   435,   436,   437,  1206,   303,   303,
     304,   303,   735,   438,   723,   997,   724,   725,   726,   727,
     728,   474,   729,   730,   731,   732,   304,   733,   484,   331,
     863,   439,   864,   638,   488,  1020,   492,   303,   304,   865,
     866,   304,   303,   303,   312,   995,   343,   439,  1011,   638,
     317,   632,   587,  1015,  1158,   738,   735,   304,   588,   881,
     882,   304,   304,   312,   744,   480,   797,   798,   116,   736,
     867,   748,   624,   124,   304,   304,   891,   304,   119,  1160,
     868,   481,   750,   485,   126,   373,   753,   120,   332,   489,
    1168,   493,  1170,   494,   634,  1093,  1064,  1109,   333,   127,
    1163,  1063,   322,   304,   313,   219,   892,   140,   304,   304,
     318,  1185,  1066,   737,   684,   869,  1072,  1126,   893,   898,
     153,   745,   903,   313,   894,   266,   265,   793,   749,  1141,
    1142,   208,  1150,   273,   209,   264,   794,   870,   297,   751,
     825,   298,   314,   754,   895,   230,   904,   267,   319,   826,
     268,   205,  1094,   590,  1110,   344,   269,  1164,  1151,   591,
     155,   314,   270,  1208,  1209,   367,   283,   160,   888,   888,
     888,   274,   157,   289,   291,   341,   290,   889,   889,   889,
    1215,   292,  1217,  1218,  1178,   342,   293,   294,   299,  1103,
     965,   966,   833,   968,   300,   837,   969,   970,   971,   972,
     973,   974,   975,   976,   977,   978,   301,   980,   981,   982,
     983,   984,   985,   986,   987,   988,   989,   990,   858,   317,
    1123,  1267,   688,  1245,  1075,  1247,   165,   762,   629,   216,
    1175,   211,   421,  1076,    99,   100,  1253,   679,   680,   888,
     681,   422,   677,   678,   679,   680,   463,   681,   889,   444,
     465,   466,   117,   118,   423,  1263,   445,   464,   201,  1266,
     216,   454,   424,   220,   455,   729,   730,   731,   732,   226,
     733,   469,   425,   689,  1010,  1276,   221,  1149,   495,   318,
     671,   230,   672,   673,   674,   675,   676,   229,   677,   678,
     679,   680,   421,   681,  1286,   731,   732,   230,   733,  1289,
     723,   422,   724,   725,   726,   727,   728,   963,   729,   730,
     731,   732,   233,   733,   423,   121,   122,   319,   142,   143,
     144,   145,   424,   161,   162,   322,   326,   426,   328,   329,
     367,   370,   425,   415,   419,   420,   460,   482,   486,   427,
     428,   994,   996,   490,   933,   429,   430,   431,   432,   433,
     434,   435,   436,   437,   495,   554,  1012,   555,   556,  1016,
     438,   557,   558,   559,   979,   560,   561,   562,   563,   564,
     565,   566,   567,   568,  1091,   569,   570,   571,   572,   573,
    1191,  1192,  1193,  1194,   574,  1195,   575,   426,   581,   584,
    1089,   585,   586,   589,   439,   592,   638,   595,  1198,   427,
     428,   596,   597,   598,   599,   429,   430,   431,   432,   433,
     434,   435,   436,   437,   600,   601,   602,   603,   604,   605,
     438,   606,  1202,   607,   608,   609,  1204,   610,   611,   614,
    1207,   615,   616,   617,   621,   166,   167,   168,   169,   170,
     171,   172,   202,   620,   622,   623,   203,   173,  1219,   626,
     627,   174,   635,   647,   439,   628,   638,   631,   163,   636,
    1227,   637,  1228,   690,   645,     1,     2,   175,   646,   204,
     654,   650,   655,   656,   657,     3,     4,     5,   660,   594,
     661,   651,     6,   691,   662,   663,     7,     8,     9,   664,
      10,   665,    11,    12,    13,    14,   666,   667,   668,   692,
     693,   669,   683,   694,   695,    15,   696,   697,    16,  1259,
     176,   177,   698,   699,   700,   701,  1261,   702,   704,   703,
     705,    17,   706,   707,   708,   709,   746,   712,   178,   179,
     713,   714,  1270,   747,    18,    19,    20,   752,   755,   671,
      21,   672,   673,   674,   675,   676,   756,   677,   678,   679,
     680,    22,   681,    23,  1285,    24,    25,    26,    27,    28,
     715,   180,   181,   716,    29,    30,   757,   758,   717,    31,
      32,    33,    34,   759,  1176,   761,   718,   760,   719,    35,
      36,   720,    37,   721,   764,   765,    38,   766,   767,    39,
      40,    41,    42,   768,   769,   770,    43,   771,   772,   773,
     774,   775,   776,  1112,   777,     1,     2,   778,   779,   780,
     781,  1199,  1200,   782,  1201,     3,     4,     5,   783,   784,
     785,   786,     6,   787,   788,   789,     7,     8,     9,    44,
      10,   792,    11,    12,    13,    14,   790,   791,   347,   804,
     801,   812,   822,   827,   803,    15,   805,   348,    16,   823,
     810,   811,   813,   816,   824,   817,   819,   821,   828,   830,
     349,    17,   831,   834,   835,   838,   841,   842,   350,   843,
    1159,   807,  1161,   844,    18,    19,    20,   845,   351,   671,
      21,   672,   673,   674,   675,   676,   849,   677,   678,   679,
     680,    22,   681,    23,   850,    24,    25,    26,    27,    28,
     851,   852,   853,   854,    29,    30,   859,   860,   885,    31,
      32,    33,    34,   884,   886,   901,   902,   681,   946,    35,
      36,   861,    37,  1260,   950,   991,    38,   874,   952,    39,
      40,    41,    42,   352,   953,   954,    43,   955,   347,   956,
     957,   959,   876,  1113,   686,   353,   354,   348,   960,   877,
     961,   355,   356,   357,   358,   359,   360,   361,   362,   363,
     349,   879,   914,   964,  1006,   347,   364,  1007,   350,    44,
    1017,   733,  1018,  1019,   348,  1021,  1022,  1023,   351,  1024,
     671,  1005,   672,   673,   674,   675,   676,   349,   677,   678,
     679,   680,  1025,   681,   421,   350,  1026,  1027,  1028,  1029,
     365,  1030,  1031,   422,   671,   351,   672,   673,   674,   675,
     676,  1009,   677,   678,   679,   680,   423,   681,  1032,  1033,
    1034,  1035,  1038,  1042,   424,  1043,  1045,  1050,  1048,  1051,
    1052,  1059,  1057,   352,   425,  1056,  1058,  1069,   577,  1073,
    1074,  1078,   682,  1061,   948,   353,   354,  1062,  1079,  1081,
    1065,   355,   356,   357,   358,   359,   360,   361,   362,   363,
     352,  1067,  1077,  1080,  1083,  1084,   364,  1086,  1114,  1087,
    1088,  1090,   353,   354,  1092,  1095,  1102,  1104,   355,   356,
     357,   358,   359,   360,   361,   362,   363,    -1,  1124,   426,
    1129,  1154,  1155,   364,  1145,  1156,  1162,  1169,  1173,  1153,
     365,   427,   428,  1171,   234,  1174,  1186,   429,   430,   431,
     432,   433,   434,   435,   436,   437,  1188,  1189,  1196,   203,
     173,  1172,   438,  1216,   174,  1205,  1212,   365,  1210,  1229,
    1221,   166,   167,   168,   169,   170,   171,   172,   235,   236,
     175,  1213,   204,   173,   234,   237,  1214,   174,  1226,  1242,
    1244,  1248,   238,   239,   240,  1256,   439,   241,   242,   203,
     243,   244,  1249,   175,   245,   246,   247,   248,   249,   250,
     251,  1257,   252,   253,   254,  1258,  1264,  1265,   235,   236,
     255,  1268,   204,   176,   177,   237,   256,  1250,   257,  1275,
    1269,  1278,   238,   258,  1279,  1283,  1282,  1284,  1288,   164,
    1179,   178,   179,   633,   259,   459,   176,   177,  1144,   372,
    1008,   829,   802,   967,   619,  1004,   260,   216,   630,   625,
     255,   856,   261,   457,   178,   179,   992,   580,   257,  1211,
     373,  1184,   374,   375,   180,   181,   685,   671,  1049,   672,
     673,   674,   675,   676,   259,   677,   678,   679,   680,   905,
     681,   815,   237,  1014,   376,   377,   260,   180,   181,   238,
     820,     0,   261,   334,  1055,     0,   331,     0,     0,     0,
       0,     0,     0,     0,   180,   181,     0,     0,   763,     0,
     372,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   378,     0,   379,   257,   380,   336,     0,   947,
       0,   373,   381,   374,   375,     0,   382,     0,     0,     0,
     337,     0,     0,     0,   383,   384,   385,     0,     0,     0,
     386,   387,   388,   237,   216,   376,   377,     0,     0,     0,
     238,   475,     0,     0,     0,     0,   671,   331,   672,   673,
     674,   675,   676,     0,   677,   678,   679,   680,   671,   681,
     672,   673,   674,   675,   676,     0,   677,   678,   679,   680,
       0,   681,     0,   378,     0,   379,   257,   380,   336,     0,
       0,     0,     0,   381,     0,     0,     0,   382,     0,     0,
       0,   337,     0,     0,     0,   383,   384,   385,     0,     0,
       0,   386,   387,   388,     0,   216,     0,     0,     0,     0,
    1115,     0,     0,     0,     0,   671,     0,   672,   673,   674,
     675,   676,  1116,   677,   678,   679,   680,   671,   681,   672,
     673,   674,   675,   676,     0,   677,   678,   679,   680,   671,
     681,   672,   673,   674,   675,   676,     0,   677,   678,   679,
     680,   671,   681,   672,   673,   674,   675,   676,     0,   677,
     678,   679,   680,   671,   681,   672,   673,   674,   675,   676,
       0,   677,   678,   679,   680,     0,   681,     0,     0,  1117,
       0,     0,     0,     0,   671,     0,   672,   673,   674,   675,
     676,  1118,   677,   678,   679,   680,   723,   681,   724,   725,
     726,   727,   728,  1119,   729,   730,   731,   732,   723,   733,
     724,   725,   726,   727,   728,  1120,   729,   730,   731,   732,
     723,   733,   724,   725,   726,   727,   728,  1121,   729,   730,
     731,   732,   723,   733,   724,   725,   726,   727,   728,     0,
     729,   730,   731,   732,     0,   733,     0,     0,  1122,     0,
       0,     0,     0,   723,     0,   724,   725,   726,   727,   728,
    1131,   729,   730,   731,   732,   723,   733,   724,   725,   726,
     727,   728,  1132,   729,   730,   731,   732,   723,   733,   724,
     725,   726,   727,   728,  1133,   729,   730,   731,   732,   723,
     733,   724,   725,   726,   727,   728,  1134,   729,   730,   731,
     732,   723,   733,   724,   725,   726,   727,   728,     0,   729,
     730,   731,   732,     0,   733,     0,     0,  1135,     0,     0,
       0,     0,   723,     0,   724,   725,   726,   727,   728,  1136,
     729,   730,   731,   732,   671,   733,   672,   673,   674,   675,
     676,  1137,   677,   678,   679,   680,   671,   681,   672,   673,
     674,   675,   676,  1138,   677,   678,   679,   680,   671,   681,
     672,   673,   674,   675,   676,  1139,   677,   678,   679,   680,
     671,   681,   672,   673,   674,   675,   676,     0,   677,   678,
     679,   680,     0,   681,     0,     0,  1140,     0,     0,     0,
       0,   723,     0,   724,   725,   726,   727,   728,  1143,   729,
     730,   731,   732,   723,   733,   724,   725,   726,   727,   728,
    1203,   729,   730,   731,   732,   671,   733,   672,   673,   674,
     675,   676,  1233,   677,   678,   679,   680,   671,   681,   672,
     673,   674,   675,   676,  1234,   677,   678,   679,   680,   723,
     681,   724,   725,   726,   727,   728,     0,   729,   730,   731,
     732,     0,   733,     0,     0,  1237,     0,     0,     0,     0,
     671,     0,   672,   673,   674,   675,   676,  1238,   677,   678,
     679,   680,   671,   681,   672,   673,   674,   675,   676,  1240,
     677,   678,   679,   680,   671,   681,   672,   673,   674,   675,
     676,  1272,   677,   678,   679,   680,   723,   681,   724,   725,
     726,   727,   728,  1273,   729,   730,   731,   732,   671,   733,
     672,   673,   674,   675,   676,     0,   677,   678,   679,   680,
       0,   681,  1068,     0,     0,     0,     0,   671,     0,   672,
     673,   674,   675,   676,  1125,   677,   678,   679,   680,   723,
     681,   724,   725,   726,   727,   728,  1157,   729,   730,   731,
     732,   671,   733,   672,   673,   674,   675,   676,  1177,   677,
     678,   679,   680,   671,   681,   672,   673,   674,   675,   676,
    1187,   677,   678,   679,   680,   671,   681,   672,   673,   674,
     675,   676,     0,   677,   678,   679,   680,     0,   681,  1197,
       0,     0,     0,     0,   671,     0,   672,   673,   674,   675,
     676,  1220,   677,   678,   679,   680,   671,   681,   672,   673,
     674,   675,   676,  1232,   677,   678,   679,   680,   671,   681,
     672,   673,   674,   675,   676,  1236,   677,   678,   679,   680,
     671,   681,   672,   673,   674,   675,   676,  1241,   677,   678,
     679,   680,   671,   681,   672,   673,   674,   675,   676,     0,
     677,   678,   679,   680,     0,   681,  1243,     0,     0,     0,
       0,   723,     0,   724,   725,   726,   727,   728,  1251,   729,
     730,   731,   732,   671,   733,   672,   673,   674,   675,   676,
    1254,   677,   678,   679,   680,     0,   681,     0,     0,     0,
       0,     0,  1255,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1274,   498,   499,   500,   501,   502,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,   517,   518,     0,     0,     0,   519,
     520,     0,   521,   522,   523,   524,   671,     0,   672,   673,
     674,   675,   676,  1107,   677,   678,   679,   680,   671,   681,
     672,   673,   674,   675,   676,  1108,   677,   678,   679,   680,
     671,   681,   672,   673,   674,   675,   676,  1111,   677,   678,
     679,   680,   723,   681,   724,   725,   726,   727,   728,  1127,
     729,   730,   731,   732,   723,   733,   724,   725,   726,   727,
     728,  1128,   729,   730,   731,   732,   723,   733,   724,   725,
     726,   727,   728,  1130,   729,   730,   731,   732,   671,   733,
     672,   673,   674,   675,   676,  1152,   677,   678,   679,   680,
     671,   681,   672,   673,   674,   675,   676,  1235,   677,   678,
     679,   680,   723,   681,   724,   725,   726,   727,   728,  1239,
     729,   730,   731,   732,     0,   733
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
       136,   231,   117,   100,   111,   111,   639,   353,   354,    24,
      25,   620,   111,   622,   111,   639,    31,   111,   627,   365,
     111,    36,    37,    38,   370,   143,   111,    42,   556,   111,
     123,   137,   138,   129,    95,   111,   143,    95,   137,   138,
     123,   111,   570,   137,   138,   111,   137,   138,   145,   643,
     111,   147,   111,   111,    95,   111,   824,  1079,   596,    31,
     111,    51,    22,    42,   912,   111,   111,   111,    51,    42,
      42,    42,    51,   111,   602,   603,   604,    33,    51,    51,
      51,    82,    51,   611,    51,   111,    82,    51,    42,    82,
     345,    31,    64,   148,    42,   325,     4,    51,    22,   148,
      72,   203,    42,    51,     4,  1095,  1096,  1097,    22,    82,
      82,    51,   204,    51,     6,    31,   148,    31,   204,    25,
      82,    79,    33,    31,    64,   148,   148,    82,  1213,   203,
     203,   204,    72,    91,    82,   209,    82,   100,   148,    82,
      51,   203,    82,    82,    22,    33,   104,    82,    82,   211,
     205,   203,    44,    31,    82,   204,    82,   120,    82,   211,
      52,    82,    82,   691,    82,   137,    82,    82,    82,   148,
     698,   131,    82,   205,    82,    82,   149,   149,   150,    82,
     203,   203,    82,   155,   156,   157,   158,   159,   160,   161,
     162,   163,   148,   148,   100,   205,   120,   137,   170,  1284,
      92,   149,   203,   204,   459,   148,   120,   203,   204,   149,
     150,   204,   151,    82,   204,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   752,   204,   209,   148,   148,    82,
     170,  1221,   204,   204,   206,   127,   203,   210,    31,   208,
     204,   151,   204,  1265,   782,   148,   592,   203,   203,    42,
     204,   779,   780,   131,    31,   783,   784,   203,    51,   787,
     203,   789,   790,   204,   204,    42,   206,   876,   207,   203,
     879,    64,   211,   212,    51,   203,   204,   203,   204,    72,
      82,   203,   203,   203,    33,   203,    33,    64,   203,    82,
     148,   427,   428,    42,   204,    72,   203,   207,   151,   203,
     203,   211,   212,   439,    82,    82,   148,    82,   654,   655,
     656,   657,    82,    82,   660,   661,   662,   663,   664,   665,
     666,   667,   668,   669,    82,   671,   672,   673,   674,   675,
     676,   677,   678,   679,   680,   681,  1184,   683,    82,    42,
      82,    51,    82,   689,   137,    33,   148,   205,    82,   203,
      82,   204,    82,   203,   207,    82,   149,   150,   211,   212,
     137,   203,   155,   156,   157,   158,   159,   160,   161,   162,
     163,    82,   149,   150,     6,    82,    82,   170,   155,   156,
     157,   158,   159,   160,   161,   162,   163,  1155,    82,    82,
     148,    82,   148,   170,   141,   741,   143,   144,   145,   146,
     147,   203,   149,   150,   151,   152,   148,   154,   148,    60,
      42,   204,    44,   206,   148,   761,   148,    82,   148,    51,
      52,   148,    82,    82,    22,   203,    22,   204,   203,   206,
      22,    31,   205,   203,   203,    31,   148,   148,   211,   149,
     150,   148,   148,    22,   148,   203,   149,   150,   203,   205,
      82,   148,    31,   203,   148,   148,   203,   148,    33,   203,
      92,   203,   148,   203,   203,    24,   148,    42,   119,   203,
    1079,   203,  1081,   203,    31,   148,   203,   148,   129,   203,
     148,   827,    82,   148,    82,   582,    45,   203,   148,   148,
      82,  1100,   203,   205,    31,   127,   203,   203,    57,  1093,
      82,   205,    56,    82,    63,   612,   612,    42,   205,   203,
     203,   572,   203,   612,   572,   612,    51,   149,   612,   205,
      42,   612,   120,   205,    83,    82,    80,   612,   120,    51,
     612,   572,   205,   205,   205,   131,   612,   205,   203,   211,
      82,   120,   612,   203,   203,    82,   612,   203,  1181,  1182,
    1183,   612,    82,   612,   612,   648,   612,  1181,  1182,  1183,
    1169,   612,  1171,  1172,  1092,   648,   612,   612,   612,   915,
     706,   707,   587,   709,   612,   590,   712,   713,   714,   715,
     716,   717,   718,   719,   720,   721,   612,   723,   724,   725,
     726,   727,   728,   729,   730,   731,   732,   733,   613,    22,
     946,    31,    31,  1212,    42,  1214,   203,    31,    31,   118,
      31,    23,    42,    51,   203,   204,  1225,   151,   152,  1252,
     154,    51,   149,   150,   151,   152,   744,   154,  1252,   735,
     139,   140,   203,   204,    64,  1244,   735,   744,   203,  1248,
     118,   735,    72,    42,   735,   149,   150,   151,   152,    82,
     154,   748,    82,    82,   750,  1264,   203,  1003,    82,    82,
     141,    82,   143,   144,   145,   146,   147,    82,   149,   150,
     151,   152,    42,   154,  1283,   151,   152,    82,   154,  1288,
     141,    51,   143,   144,   145,   146,   147,   702,   149,   150,
     151,   152,    33,   154,    64,   203,   204,   120,   203,   204,
     203,   204,    72,   203,   204,    82,    38,   137,    42,   203,
      82,    82,    82,   126,    51,   203,   208,   126,   126,   149,
     150,   736,   737,   126,   205,   155,   156,   157,   158,   159,
     160,   161,   162,   163,    82,    33,   751,    33,    33,   754,
     170,    33,    33,    33,   205,    33,    33,    33,    33,    33,
      33,    33,   148,   205,   890,    33,    33,    33,   148,   205,
    1106,  1107,  1108,  1109,   205,  1111,    82,   137,    33,   203,
     885,   203,    82,    82,   204,    33,   206,    82,  1124,   149,
     150,    33,    33,    33,    33,   155,   156,   157,   158,   159,
     160,   161,   162,   163,    33,    33,    33,    33,    33,    33,
     170,    33,  1148,    33,    33,    33,  1152,    33,    33,    33,
    1156,    82,    82,    82,    82,     9,    10,    11,    12,    13,
      14,    15,    16,   148,   148,    82,    20,    21,  1174,    82,
     148,    25,    33,    33,   204,    82,   206,   204,     0,   205,
    1186,   203,  1188,    33,   205,     7,     8,    41,   205,    43,
     204,    82,   204,   204,   204,    17,    18,    19,   204,  1089,
     204,    82,    24,    33,   204,   204,    28,    29,    30,   204,
      32,   204,    34,    35,    36,    37,   204,   204,   204,    33,
      33,   204,   204,    33,    33,    47,    33,    33,    50,  1235,
      84,    85,    33,    33,    33,    33,  1242,   205,    82,   126,
     203,    63,   204,   204,   204,   204,    33,   204,   102,   103,
     204,   204,  1258,    33,    76,    77,    78,    33,   126,   141,
      82,   143,   144,   145,   146,   147,    82,   149,   150,   151,
     152,    93,   154,    95,  1280,    97,    98,    99,   100,   101,
     204,   135,   136,   204,   106,   107,   126,    82,   204,   111,
     112,   113,   114,   126,  1090,   148,   204,    82,   204,   121,
     122,   204,   124,   204,    33,    33,   128,    33,    33,   131,
     132,   133,   134,    33,    33,    33,   138,    33,    33,    33,
      33,    33,    33,   205,    33,     7,     8,    33,    33,    33,
      33,  1127,  1128,    33,  1130,    17,    18,    19,    33,    33,
      33,    33,    24,    33,    33,    33,    28,    29,    30,   171,
      32,   148,    34,    35,    36,    37,    33,   205,    42,    82,
      51,    82,   203,    33,    51,    47,    51,    51,    50,   204,
      51,    51,    82,    51,    82,    51,    51,    51,    51,   203,
      64,    63,   203,    82,   203,    82,    51,    51,    72,    51,
    1065,   208,  1067,    51,    76,    77,    78,    51,    82,   141,
      82,   143,   144,   145,   146,   147,    51,   149,   150,   151,
     152,    93,   154,    95,    51,    97,    98,    99,   100,   101,
      51,    51,    82,   204,   106,   107,    82,    33,   203,   111,
     112,   113,   114,    82,    82,   203,   203,   154,    82,   121,
     122,   148,   124,  1239,    51,   203,   128,   148,    51,   131,
     132,   133,   134,   137,    51,    51,   138,    51,    42,    51,
      51,    51,   148,   205,   148,   149,   150,    51,    51,   148,
      51,   155,   156,   157,   158,   159,   160,   161,   162,   163,
      64,   148,   148,   126,    51,    42,   170,    51,    72,   171,
     126,   154,   126,   126,    51,    51,    51,    51,    82,    51,
     141,   203,   143,   144,   145,   146,   147,    64,   149,   150,
     151,   152,    51,   154,    42,    72,    51,    51,    51,    51,
     204,    51,    51,    51,   141,    82,   143,   144,   145,   146,
     147,   203,   149,   150,   151,   152,    64,   154,    51,    51,
      51,    51,    51,    51,    72,    51,    51,    42,   203,    42,
      51,    51,   203,   137,    82,   207,   203,   210,   204,    82,
      82,   148,   203,   205,   148,   149,   150,   205,   148,   148,
     205,   155,   156,   157,   158,   159,   160,   161,   162,   163,
     137,   205,   203,    82,    82,   203,   170,    51,   205,    51,
     205,    33,   149,   150,    33,   203,    82,    82,   155,   156,
     157,   158,   159,   160,   161,   162,   163,   154,    82,   137,
     205,    82,   205,   170,   203,    33,    82,   148,   203,   207,
     204,   149,   150,   148,     5,    33,    33,   155,   156,   157,
     158,   159,   160,   161,   162,   163,    33,    51,   203,    20,
      21,   148,   170,   203,    25,    51,   148,   204,   210,   207,
     203,     9,    10,    11,    12,    13,    14,    15,    39,    40,
      41,   148,    43,    21,     5,    46,   148,    25,   205,    33,
     148,   148,    53,    54,    55,    51,   204,    58,    59,    20,
      61,    62,   203,    41,    65,    66,    67,    68,    69,    70,
      71,    51,    73,    74,    75,   130,   148,   148,    39,    40,
      81,   205,    43,    84,    85,    46,    87,   203,    89,    82,
     207,   203,    53,    94,    51,   148,   210,   148,   148,    45,
    1093,   102,   103,   323,   105,   140,    84,    85,   998,     3,
     748,   582,   558,   708,   309,   744,   117,   118,   320,   315,
      81,   612,   123,   138,   102,   103,   735,   214,    89,  1163,
      24,  1098,    26,    27,   135,   136,   368,   141,   792,   143,
     144,   145,   146,   147,   105,   149,   150,   151,   152,   648,
     154,   567,    46,   753,    48,    49,   117,   135,   136,    53,
     572,    -1,   123,   122,   808,    -1,    60,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   135,   136,    -1,    -1,   496,    -1,
       3,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    86,    -1,    88,    89,    90,    91,    -1,   203,
      -1,    24,    96,    26,    27,    -1,   100,    -1,    -1,    -1,
     104,    -1,    -1,    -1,   108,   109,   110,    -1,    -1,    -1,
     114,   115,   116,    46,   118,    48,    49,    -1,    -1,    -1,
      53,   125,    -1,    -1,    -1,    -1,   141,    60,   143,   144,
     145,   146,   147,    -1,   149,   150,   151,   152,   141,   154,
     143,   144,   145,   146,   147,    -1,   149,   150,   151,   152,
      -1,   154,    -1,    86,    -1,    88,    89,    90,    91,    -1,
      -1,    -1,    -1,    96,    -1,    -1,    -1,   100,    -1,    -1,
      -1,   104,    -1,    -1,    -1,   108,   109,   110,    -1,    -1,
      -1,   114,   115,   116,    -1,   118,    -1,    -1,    -1,    -1,
     205,    -1,    -1,    -1,    -1,   141,    -1,   143,   144,   145,
     146,   147,   205,   149,   150,   151,   152,   141,   154,   143,
     144,   145,   146,   147,    -1,   149,   150,   151,   152,   141,
     154,   143,   144,   145,   146,   147,    -1,   149,   150,   151,
     152,   141,   154,   143,   144,   145,   146,   147,    -1,   149,
     150,   151,   152,   141,   154,   143,   144,   145,   146,   147,
      -1,   149,   150,   151,   152,    -1,   154,    -1,    -1,   205,
      -1,    -1,    -1,    -1,   141,    -1,   143,   144,   145,   146,
     147,   205,   149,   150,   151,   152,   141,   154,   143,   144,
     145,   146,   147,   205,   149,   150,   151,   152,   141,   154,
     143,   144,   145,   146,   147,   205,   149,   150,   151,   152,
     141,   154,   143,   144,   145,   146,   147,   205,   149,   150,
     151,   152,   141,   154,   143,   144,   145,   146,   147,    -1,
     149,   150,   151,   152,    -1,   154,    -1,    -1,   205,    -1,
      -1,    -1,    -1,   141,    -1,   143,   144,   145,   146,   147,
     205,   149,   150,   151,   152,   141,   154,   143,   144,   145,
     146,   147,   205,   149,   150,   151,   152,   141,   154,   143,
     144,   145,   146,   147,   205,   149,   150,   151,   152,   141,
     154,   143,   144,   145,   146,   147,   205,   149,   150,   151,
     152,   141,   154,   143,   144,   145,   146,   147,    -1,   149,
     150,   151,   152,    -1,   154,    -1,    -1,   205,    -1,    -1,
      -1,    -1,   141,    -1,   143,   144,   145,   146,   147,   205,
     149,   150,   151,   152,   141,   154,   143,   144,   145,   146,
     147,   205,   149,   150,   151,   152,   141,   154,   143,   144,
     145,   146,   147,   205,   149,   150,   151,   152,   141,   154,
     143,   144,   145,   146,   147,   205,   149,   150,   151,   152,
     141,   154,   143,   144,   145,   146,   147,    -1,   149,   150,
     151,   152,    -1,   154,    -1,    -1,   205,    -1,    -1,    -1,
      -1,   141,    -1,   143,   144,   145,   146,   147,   205,   149,
     150,   151,   152,   141,   154,   143,   144,   145,   146,   147,
     205,   149,   150,   151,   152,   141,   154,   143,   144,   145,
     146,   147,   205,   149,   150,   151,   152,   141,   154,   143,
     144,   145,   146,   147,   205,   149,   150,   151,   152,   141,
     154,   143,   144,   145,   146,   147,    -1,   149,   150,   151,
     152,    -1,   154,    -1,    -1,   205,    -1,    -1,    -1,    -1,
     141,    -1,   143,   144,   145,   146,   147,   205,   149,   150,
     151,   152,   141,   154,   143,   144,   145,   146,   147,   205,
     149,   150,   151,   152,   141,   154,   143,   144,   145,   146,
     147,   205,   149,   150,   151,   152,   141,   154,   143,   144,
     145,   146,   147,   205,   149,   150,   151,   152,   141,   154,
     143,   144,   145,   146,   147,    -1,   149,   150,   151,   152,
      -1,   154,   203,    -1,    -1,    -1,    -1,   141,    -1,   143,
     144,   145,   146,   147,   203,   149,   150,   151,   152,   141,
     154,   143,   144,   145,   146,   147,   203,   149,   150,   151,
     152,   141,   154,   143,   144,   145,   146,   147,   203,   149,
     150,   151,   152,   141,   154,   143,   144,   145,   146,   147,
     203,   149,   150,   151,   152,   141,   154,   143,   144,   145,
     146,   147,    -1,   149,   150,   151,   152,    -1,   154,   203,
      -1,    -1,    -1,    -1,   141,    -1,   143,   144,   145,   146,
     147,   203,   149,   150,   151,   152,   141,   154,   143,   144,
     145,   146,   147,   203,   149,   150,   151,   152,   141,   154,
     143,   144,   145,   146,   147,   203,   149,   150,   151,   152,
     141,   154,   143,   144,   145,   146,   147,   203,   149,   150,
     151,   152,   141,   154,   143,   144,   145,   146,   147,    -1,
     149,   150,   151,   152,    -1,   154,   203,    -1,    -1,    -1,
      -1,   141,    -1,   143,   144,   145,   146,   147,   203,   149,
     150,   151,   152,   141,   154,   143,   144,   145,   146,   147,
     203,   149,   150,   151,   152,    -1,   154,    -1,    -1,    -1,
      -1,    -1,   203,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   203,   172,   173,   174,   175,   176,
     177,   178,   179,   180,   181,   182,   183,   184,   185,   186,
     187,   188,   189,   190,   191,   192,    -1,    -1,    -1,   196,
     197,    -1,   199,   200,   201,   202,   141,    -1,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   141,   154,
     143,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     141,   154,   143,   144,   145,   146,   147,   148,   149,   150,
     151,   152,   141,   154,   143,   144,   145,   146,   147,   148,
     149,   150,   151,   152,   141,   154,   143,   144,   145,   146,
     147,   148,   149,   150,   151,   152,   141,   154,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   141,   154,
     143,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     141,   154,   143,   144,   145,   146,   147,   148,   149,   150,
     151,   152,   141,   154,   143,   144,   145,   146,   147,   148,
     149,   150,   151,   152,    -1,   154
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
     132,   133,   134,   138,   171,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   227,   228,   229,   230,   233,   235,
     238,   243,   254,   255,   259,   263,   266,   269,   272,   278,
     284,   287,   292,   295,   298,   301,   302,   305,   306,   308,
     309,   310,   314,   315,   316,   317,   323,   326,   332,   335,
     336,   339,    51,   204,    51,   204,   203,   204,   203,   203,
     204,    33,    42,    51,    82,   204,    82,   204,   203,    82,
     203,   204,   275,   203,   203,   203,   203,   203,   204,    33,
      42,   203,   204,   204,   203,    33,   203,   203,   203,   204,
     275,   275,    82,   226,    33,    51,   324,   204,   204,   275,
     203,    33,   203,   204,   203,   204,   203,   204,   275,   203,
     204,   275,   275,    82,   223,    82,   224,    82,   225,   275,
     203,   203,   204,     0,   215,   203,     9,    10,    11,    12,
      13,    14,    15,    21,    25,    41,    84,    85,   102,   103,
     135,   136,   329,   330,   331,   363,   364,   365,   366,   367,
     398,   399,   401,   402,   405,   406,   407,   408,   409,   410,
     411,   203,    16,    20,    43,   330,   333,   334,   371,   388,
     412,    23,     4,    82,   311,   312,   118,   267,   268,   344,
      42,   203,    51,   203,   203,   211,    82,   203,   211,    82,
      82,   236,   237,    33,     5,    39,    40,    46,    53,    54,
      55,    58,    59,    61,    62,    65,    66,    67,    68,    69,
      70,    71,    73,    74,    75,    81,    87,    89,    94,   105,
     117,   123,   293,   294,   344,   354,   363,   364,   365,   366,
     367,   368,   369,   370,   371,   372,   373,   374,   375,   376,
     377,   378,   379,   380,   381,   382,   383,   384,   385,   386,
     387,   388,   389,   390,   391,   393,   394,   398,   399,   400,
     401,   402,   403,    82,   148,   203,    22,    82,   120,   279,
     280,   281,    22,    82,   120,   288,   289,    22,    82,   120,
     285,   286,    82,   239,   240,   236,    38,   234,    42,   203,
     244,    60,   119,   129,   346,    79,    91,   104,   318,   319,
     395,   396,   397,    22,   131,   256,   257,    42,    51,    64,
      72,    82,   137,   149,   150,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   170,   204,   231,    82,   303,   304,
      82,   307,     3,    24,    26,    27,    48,    49,    86,    88,
      90,    96,   100,   108,   109,   110,   114,   115,   116,   274,
     343,   344,   345,   346,   347,   348,   349,   350,   351,   352,
     353,   354,   355,   356,   357,   358,   360,   361,   362,   370,
     392,   396,   397,   203,   203,   126,    82,   148,   203,    51,
     203,    42,    51,    64,    72,    82,   137,   149,   150,   155,
     156,   157,   158,   159,   160,   161,   162,   163,   170,   204,
     251,   253,   296,   297,   354,   370,   371,   380,   386,   387,
     388,   389,   390,   391,   398,   399,   400,   296,   203,   256,
     208,   270,   271,   357,   363,   139,   140,   264,   265,   344,
     440,   441,   273,   274,   203,   125,   274,   327,   328,   404,
     203,   203,   126,    82,   148,   203,   126,    82,   148,   203,
     126,    82,   148,   203,   203,    82,   340,   341,   172,   173,
     174,   175,   176,   177,   178,   179,   180,   181,   182,   183,
     184,   185,   186,   187,   188,   189,   190,   191,   192,   196,
     197,   199,   200,   201,   202,   337,   338,   413,   414,   415,
     416,   417,   418,   419,   420,   421,   422,   423,   424,   425,
     426,   427,   428,   429,   430,   431,   432,   433,   434,   435,
     436,   437,   438,   439,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,   148,   205,    33,
      33,    33,   148,   205,   205,    82,   148,   204,   313,    31,
     312,    33,   148,   205,   203,   203,    82,   205,   211,    82,
     205,   211,    33,    31,   237,    82,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   148,   205,    33,    82,    82,    82,    31,   280,
     148,    82,   148,    82,    31,   289,    82,   148,    82,    31,
     286,   204,    31,   240,    31,    33,   205,   203,   206,   249,
     250,   251,   252,   148,   205,   205,   205,    33,   148,   205,
      82,    82,    31,   257,   204,   204,   204,   204,   231,   231,
     204,   204,   204,   204,   204,   204,   204,   204,   204,   204,
     231,   141,   143,   144,   145,   146,   147,   149,   150,   151,
     152,   154,   203,   204,    31,   304,   148,   231,    31,    82,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   205,   126,    82,   203,   204,   204,   204,   204,
     251,   251,   204,   204,   204,   204,   204,   204,   204,   204,
     204,   204,   251,   141,   143,   144,   145,   146,   147,   149,
     150,   151,   152,   154,   325,   148,   205,   205,    31,    42,
      51,   204,   261,   262,   148,   205,    33,    33,   148,   205,
     148,   205,    33,   148,   205,   126,    82,   126,    82,   126,
      82,   148,    31,   341,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,   205,   148,    42,    51,   342,    42,   149,   150,   277,
     342,    51,   277,    51,    82,    51,    51,   208,   444,   445,
      51,    51,    82,    82,   442,   331,    51,    51,   342,    51,
     334,    51,   203,   204,    82,    42,    51,    33,    51,   268,
     203,   203,   203,   275,    82,   203,   203,   275,    82,   231,
     445,    51,    51,    51,    51,    51,   342,   342,   342,    51,
      51,    51,    51,    82,   204,   342,   294,   203,   275,    82,
      33,   148,     6,    42,    44,    51,    52,    82,    92,   127,
     149,   282,   290,   291,   148,   291,   148,   148,   291,   148,
      51,   149,   150,   276,    82,   203,    82,    31,   250,   252,
      33,   203,    45,    57,    63,    83,   241,   242,   358,   359,
     248,   203,   203,    56,    80,   319,    82,   151,   207,   211,
     212,   320,   321,   322,   148,    33,   148,   203,   231,   231,
     231,   232,   231,   231,   231,   231,   231,   231,   231,   231,
     231,   231,   231,   205,   231,   231,   231,   231,   231,   231,
     231,   231,   231,   231,   231,   231,    82,   203,   148,   231,
      51,   342,    51,    51,    51,    51,    51,    51,   342,    51,
      51,    51,   203,   275,   126,   251,   251,   276,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   205,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   203,   297,   203,   275,   203,   275,   231,   203,   209,
      42,    51,   148,   204,   271,   203,    51,    51,   265,   203,
     274,   203,   275,   342,   328,   203,   275,   126,   126,   126,
     231,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,   342,   342,    51,   445,
     342,   342,    51,    51,   342,    51,   342,   342,   203,   337,
      42,    42,    51,   443,   209,   443,   207,   203,   203,    51,
     313,   205,   205,   231,   203,   205,   203,   205,   203,   210,
     299,   300,   203,    82,    82,    42,    51,   203,   148,   148,
      82,   148,   291,    82,   203,   291,    51,    51,   205,   236,
      33,   251,    33,   148,   205,   203,   246,   245,   148,   203,
     204,   322,    82,   231,    82,   100,   120,   148,   148,   148,
     205,   148,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   205,   205,   231,    82,   203,   203,   148,   148,   205,
     148,   205,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   203,   203,   205,   262,   203,    42,    51,   204,   231,
     203,   203,   148,   207,    82,   205,    33,   203,   203,   275,
     203,   275,    82,   148,   205,   283,   291,   290,   291,   148,
     291,   148,   148,   203,    33,    31,   251,   203,   342,   242,
     247,   249,   249,   249,   321,   291,    33,   203,    33,    51,
     258,   231,   231,   231,   231,   231,   203,   203,   231,   251,
     251,   251,   231,   205,   231,    51,   313,   231,   203,   203,
     210,   299,   148,   148,   148,   291,   203,   291,   291,   231,
     203,   203,    31,    31,    31,   204,   205,   231,   231,   207,
     148,   203,   203,   205,   205,   148,   203,   205,   205,   148,
     205,   203,    33,   203,   148,   291,   283,   291,   148,   203,
     203,   203,   249,   291,   203,   203,    51,    51,   130,   231,
     251,   231,   210,   291,   148,   148,   291,    31,   205,   207,
     231,   260,   205,   205,   203,    82,   291,   290,   203,    51,
     148,   203,   210,   148,   148,   231,   291,   283,   148,   291
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
     455,   456,   457,    59,    40,    41,    35,    58,    91,    93,
      39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   213,   214,   214,   215,   215,   215,   215,   215,   215,
     215,   215,   215,   215,   215,   215,   215,   215,   215,   215,
     215,   215,   215,   215,   215,   215,   215,   215,   215,   215,
     215,   215,   215,   215,   215,   215,   215,   215,   215,   215,
     215,   215,   215,   215,   215,   216,   216,   216,   216,   217,
     217,   218,   219,   220,   221,   222,   223,   223,   223,   223,
     223,   223,   224,   224,   224,   224,   224,   224,   225,   225,
     225,   225,   225,   225,   226,   226,   226,   226,   226,   226,
     227,   227,   228,   228,   229,   229,   230,   231,   231,   231,
     231,   231,   231,   231,   231,   231,   231,   231,   231,   231,
     231,   231,   231,   231,   231,   231,   231,   231,   231,   231,
     231,   231,   231,   231,   231,   231,   231,   231,   232,   232,
     233,   233,   234,   235,   236,   236,   237,   238,   239,   239,
     240,   241,   241,   242,   242,   242,   242,   242,   244,   243,
     245,   243,   246,   243,   247,   243,   248,   243,   249,   249,
     249,   249,   250,   250,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   252,   253,   253,   254,   255,   256,
     256,   257,   257,   257,   257,   257,   258,   258,   258,   258,
     259,   260,   260,   261,   261,   262,   262,   262,   262,   262,
     262,   262,   262,   262,   263,   263,   264,   264,   265,   265,
     265,   266,   266,   267,   267,   268,   269,   269,   270,   270,
     271,   271,   272,   272,   272,   272,   273,   273,   274,   274,
     274,   274,   274,   274,   274,   274,   274,   274,   274,   274,
     274,   274,   274,   274,   274,   274,   274,   274,   274,   274,
     274,   275,   275,   275,   275,   275,   275,   276,   276,   276,
     277,   277,   277,   278,   279,   279,   280,   281,   281,   281,
     282,   282,   282,   282,   282,   283,   283,   283,   283,   284,
     285,   285,   286,   286,   286,   287,   288,   288,   289,   289,
     289,   290,   290,   290,   290,   290,   291,   291,   291,   291,
     291,   291,   292,   292,   292,   292,   293,   293,   294,   294,
     294,   294,   294,   294,   294,   294,   294,   294,   294,   294,
     294,   294,   294,   294,   294,   294,   294,   294,   294,   294,
     294,   294,   294,   294,   294,   294,   294,   294,   294,   294,
     294,   294,   294,   294,   294,   294,   294,   295,   295,   296,
     296,   297,   297,   297,   297,   297,   297,   297,   297,   297,
     297,   297,   297,   297,   298,   298,   299,   299,   300,   300,
     301,   302,   303,   303,   304,   305,   306,   307,   307,   307,
     307,   308,   309,   309,   309,   309,   310,   311,   311,   312,
     312,   312,   313,   313,   313,   314,   314,   315,   315,   315,
     315,   315,   315,   316,   316,   316,   316,   316,   316,   317,
     318,   318,   319,   319,   319,   320,   320,   320,   320,   321,
     321,   322,   322,   322,   322,   322,   324,   325,   323,   326,
     326,   326,   326,   327,   327,   328,   328,   329,   329,   329,
     329,   329,   329,   329,   330,   330,   330,   330,   330,   330,
     330,   330,   330,   330,   331,   331,   332,   332,   333,   333,
     333,   333,   334,   334,   335,   335,   336,   336,   337,   337,
     338,   338,   338,   338,   338,   338,   338,   338,   338,   338,
     338,   338,   338,   338,   338,   338,   338,   338,   338,   338,
     338,   338,   338,   338,   338,   338,   338,   339,   340,   340,
     341,   342,   342,   343,   344,   345,   346,   347,   348,   349,
     350,   351,   352,   353,   354,   355,   356,   357,   358,   359,
     360,   361,   362,   363,   364,   364,   365,   366,   367,   368,
     369,   370,   370,   371,   372,   373,   374,   375,   376,   377,
     378,   379,   380,   381,   382,   383,   384,   385,   386,   387,
     388,   389,   390,   391,   392,   393,   394,   395,   395,   396,
     397,   398,   399,   400,   401,   402,   403,   404,   405,   406,
     407,   408,   409,   410,   411,   412,   413,   414,   415,   416,
     417,   418,   419,   420,   421,   422,   423,   424,   425,   426,
     427,   428,   429,   430,   431,   432,   433,   434,   435,   436,
     437,   438,   439,   440,   441,   442,   443,   443,   444,   444,
     445
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
       4,     4,     4,     4,     6,     6,     4,     8,     1,     3,
       4,     7,     3,     4,     2,     1,     4,     4,     2,     1,
       7,     3,     1,     1,     1,     1,     1,     1,     0,     5,
       0,     8,     0,     8,     0,    10,     0,     8,     2,     2,
       1,     1,     4,     2,     3,     1,     1,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     2,
       2,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     6,     6,     8,     5,     1,     4,     4,     4,     2,
       1,     9,     6,     5,     7,     7,     3,     5,     3,     1,
       6,     3,     1,     3,     1,     5,     3,     3,     4,     2,
       2,     3,     1,     1,     2,     5,     3,     1,     1,     1,
       1,     2,     5,     3,     1,     1,     2,     5,     3,     1,
       1,     1,     2,     5,     3,     6,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     4,     3,     5,     1,     3,     2,     2,     1,
       2,     2,     1,     4,     2,     1,     4,     2,     1,     4,
       3,     5,     9,     1,     5,     3,     5,     7,     9,     4,
       2,     1,     5,     7,     4,     4,     2,     1,     7,     9,
       6,     1,     1,     1,     1,     1,     0,     1,     1,     1,
       2,     2,     2,     5,     3,     6,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     5,     6,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     5,     6,     7,     5,     1,     3,
       3,     4,     2,     1,     5,     3,     4,     4,     6,     3,
       5,     3,     2,     5,     3,     6,     4,     2,     1,     5,
       7,     9,     0,     3,     3,     2,     5,     5,     6,     3,
       7,     8,     5,     5,     6,     3,     7,     8,     5,     6,
       3,     1,     1,     1,     1,     1,     3,     4,     6,     1,
       2,     1,     1,     1,     1,     1,     0,     0,     5,     2,
       5,     3,     6,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     1,     3,     6,     1,     1,
       1,     1,     3,     1,     3,     6,     2,     5,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     4,     1,     2,
       6,     1,     1,     3,     3,     3,     1,     3,     3,     3,
       3,     1,     1,     1,     3,     3,     3,     3,     3,     3,
       1,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     1,     1,     3,     3,     3,     3,     5,     3,     3,
       3,     1,     3,     3,     3,     1,     1,     1,     1,     1,
       3,     1,     1,     1,     1,     3,     3,     3,     3,     1,
       1,     3,     3,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     3,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     3,     2,     2,
       2
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
  "VAROBS", "XLS_SHEET", "XLS_RANGE", "NORMCDF", "HOMOTOPY_SETUP",
  "HOMOTOPY_MODE", "HOMOTOPY_STEPS", "EXCLAMATION_EQUAL", "EXCLAMATION",
  "EQUAL_EQUAL", "GREATER_EQUAL", "LESS_EQUAL", "GREATER", "LESS", "COMMA",
  "MINUS", "PLUS", "DIVIDE", "TIMES", "UMINUS", "POWER", "EXP", "LOG",
  "LOG10", "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN", "SINH", "COSH",
  "TANH", "ASINH", "ACOSH", "ATANH", "SQRT", "DYNARE_SENSITIVITY",
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
  "histval_list", "histval_elem", "model_sparse_options_list",
  "model_sparse_options", "model", "@1", "@2", "@3", "@4", "@5",
  "equation_list", "equation", "hand_side", "pound_expression",
  "model_var", "shocks", "mshocks", "shock_list", "shock_elem",
  "period_list", "sigma_e", "value_list", "triangular_matrix",
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
  "osr_params", "osr", "calib_var", "calib_var_list", "calib_arg1",
  "calib_arg2", "calib", "dynatype", "dynasave", "model_comparison",
  "model_comparison_options", "model_comparison_option", "filename_list",
  "filename", "filename_elem", "planner_objective", "@6", "@7",
  "ramsey_policy", "ramsey_policy_options_list", "ramsey_policy_options",
  "bvar_prior_option", "bvar_common_option", "bvar_density_options_list",
  "bvar_density", "bvar_forecast_option", "bvar_forecast_options_list",
  "bvar_forecast", "dynare_sensitivity", "dynare_sensitivity_options_list",
  "dynare_sensitivity_option", "homotopy_setup", "homotopy_list",
  "homotopy_item", "number", "o_dr_algo", "o_solve_algo", "o_simul_algo",
  "o_linear", "o_order", "o_replic", "o_drop", "o_ar", "o_nocorr",
  "o_nofunctions", "o_nomoments", "o_irf", "o_hp_filter", "o_hp_ngrid",
  "o_periods", "o_cutoff", "o_markowitz", "o_simul", "o_simul_seed",
  "o_qz_criterium", "o_datafile", "o_nobs", "o_first_obs", "o_prefilter",
  "o_presample", "o_lik_algo", "o_lik_init", "o_nograph", "o_conf_sig",
  "o_mh_replic", "o_mh_drop", "o_mh_jscale", "o_optim", "o_mh_init_scale",
  "o_mode_file", "o_mode_compute", "o_mode_check", "o_prior_trunc",
  "o_mh_mode", "o_mh_nblcks", "o_load_mh_file", "o_loglinear",
  "o_nodiagnostic", "o_bayesian_irf", "o_tex", "o_forecast", "o_smoother",
  "o_moments_varendo", "o_filtered_vars", "o_relative_irf",
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
       214,     0,    -1,   215,    -1,   214,   215,    -1,   216,    -1,
     227,    -1,   228,    -1,   229,    -1,   243,    -1,   233,    -1,
     235,    -1,   238,    -1,   230,    -1,   254,    -1,   255,    -1,
     259,    -1,   263,    -1,   266,    -1,   269,    -1,   272,    -1,
     292,    -1,   295,    -1,   298,    -1,   278,    -1,   287,    -1,
     284,    -1,   301,    -1,   302,    -1,   305,    -1,   217,    -1,
     218,    -1,   306,    -1,   308,    -1,   309,    -1,   310,    -1,
     314,    -1,   315,    -1,   316,    -1,   317,    -1,   323,    -1,
     326,    -1,   332,    -1,   335,    -1,   336,    -1,   339,    -1,
     222,    -1,   219,    -1,   220,    -1,   221,    -1,    28,    51,
     203,    -1,    28,    51,    51,   203,    -1,   111,   275,   203,
      -1,   131,   223,   203,    -1,   132,   224,   203,    -1,   133,
     225,   203,    -1,    99,   226,   203,    -1,   223,    82,    -1,
     223,   148,    82,    -1,    82,    -1,   223,    82,   126,    -1,
     223,   148,    82,   126,    -1,    82,   126,    -1,   224,    82,
      -1,   224,   148,    82,    -1,    82,    -1,   224,    82,   126,
      -1,   224,   148,    82,   126,    -1,    82,   126,    -1,   225,
      82,    -1,   225,   148,    82,    -1,    82,    -1,   225,    82,
     126,    -1,   225,   148,    82,   126,    -1,    82,   126,    -1,
     226,    82,    -1,   226,   148,    82,    -1,    82,    -1,   226,
      82,   126,    -1,   226,   148,    82,   126,    -1,    82,   126,
      -1,   100,    51,   203,    -1,   100,    33,    51,   203,    -1,
      24,    42,   203,    -1,    24,    33,    42,   203,    -1,    63,
      42,   203,    -1,    63,    33,    42,   203,    -1,    82,    33,
     231,   203,    -1,   204,   231,   205,    -1,    82,    -1,    42,
      -1,    51,    -1,   231,   150,   231,    -1,   231,   149,   231,
      -1,   231,   151,   231,    -1,   231,   152,   231,    -1,   231,
     154,   231,    -1,   231,   147,   231,    -1,   231,   146,   231,
      -1,   231,   145,   231,    -1,   231,   144,   231,    -1,   231,
     143,   231,    -1,   231,   141,   231,    -1,   149,   231,    -1,
     150,   231,    -1,   155,   204,   231,   205,    -1,   156,   204,
     231,   205,    -1,   157,   204,   231,   205,    -1,   158,   204,
     231,   205,    -1,   159,   204,   231,   205,    -1,   160,   204,
     231,   205,    -1,   161,   204,   231,   205,    -1,   162,   204,
     231,   205,    -1,   163,   204,   231,   205,    -1,   170,   204,
     231,   205,    -1,    64,   204,   231,   148,   231,   205,    -1,
      72,   204,   231,   148,   231,   205,    -1,    82,   204,   232,
     205,    -1,   137,   204,   231,   148,   231,   148,   231,   205,
      -1,   231,    -1,   232,   148,   231,    -1,    50,   203,   236,
      31,    -1,    50,   204,   234,   205,   203,   236,    31,    -1,
      38,    33,    82,    -1,    32,   203,   236,    31,    -1,   236,
     237,    -1,   237,    -1,    82,    33,   231,   203,    -1,    47,
     203,   239,    31,    -1,   239,   240,    -1,   240,    -1,    82,
     204,   276,   205,    33,   231,   203,    -1,   241,   148,   242,
      -1,   242,    -1,    57,    -1,    45,    -1,    83,    -1,   358,
      -1,   359,    -1,    -1,    76,   203,   244,   249,    31,    -1,
      -1,    76,   204,   346,   205,   203,   245,   249,    31,    -1,
      -1,    76,   204,   129,   205,   203,   246,   249,    31,    -1,
      -1,    76,   204,   119,   148,   241,   205,   247,   203,   249,
      31,    -1,    -1,    76,   204,   119,   205,   248,   203,   249,
      31,    -1,   249,   250,    -1,   249,   252,    -1,   250,    -1,
     252,    -1,   251,    33,   251,   203,    -1,   251,   203,    -1,
     204,   251,   205,    -1,   253,    -1,    42,    -1,    51,    -1,
     251,   150,   251,    -1,   251,   149,   251,    -1,   251,   151,
     251,    -1,   251,   152,   251,    -1,   251,   147,   251,    -1,
     251,   146,   251,    -1,   251,   145,   251,    -1,   251,   144,
     251,    -1,   251,   143,   251,    -1,   251,   141,   251,    -1,
     251,   154,   251,    -1,   149,   251,    -1,   150,   251,    -1,
     155,   204,   251,   205,    -1,   156,   204,   251,   205,    -1,
     157,   204,   251,   205,    -1,   158,   204,   251,   205,    -1,
     159,   204,   251,   205,    -1,   160,   204,   251,   205,    -1,
     161,   204,   251,   205,    -1,   162,   204,   251,   205,    -1,
     163,   204,   251,   205,    -1,   170,   204,   251,   205,    -1,
      64,   204,   251,   148,   251,   205,    -1,    72,   204,   251,
     148,   251,   205,    -1,   137,   204,   251,   148,   251,   148,
     251,   205,    -1,   206,    82,    33,   251,   203,    -1,    82,
      -1,    82,   204,   276,   205,    -1,   112,   203,   256,    31,
      -1,    78,   203,   256,    31,    -1,   256,   257,    -1,   257,
      -1,   131,    82,   203,   100,   258,   203,   130,   260,   203,
      -1,   131,    82,   203,   120,   231,   203,    -1,   131,    82,
      33,   231,   203,    -1,   131,    82,   148,    82,    33,   231,
     203,    -1,    22,    82,   148,    82,    33,   231,   203,    -1,
     258,   148,    51,    -1,   258,   148,    51,   207,    51,    -1,
      51,   207,    51,    -1,    51,    -1,   113,    33,   208,   261,
     209,   203,    -1,   260,   148,   231,    -1,   231,    -1,   261,
     203,   262,    -1,   262,    -1,   262,   148,   204,   231,   205,
      -1,   262,   148,    42,    -1,   262,   148,    51,    -1,   262,
     204,   231,   205,    -1,   262,    42,    -1,   262,    51,    -1,
     204,   231,   205,    -1,    42,    -1,    51,    -1,   121,   203,
      -1,   121,   204,   264,   205,   203,    -1,   264,   148,   265,
      -1,   265,    -1,   344,    -1,   440,    -1,   441,    -1,    19,
     203,    -1,    19,   204,   267,   205,   203,    -1,   267,   148,
     268,    -1,   268,    -1,   344,    -1,   114,   203,    -1,   114,
     204,   270,   205,   203,    -1,   270,   148,   271,    -1,   271,
      -1,   357,    -1,   363,    -1,   122,   203,    -1,   122,   204,
     273,   205,   203,    -1,   122,   275,   203,    -1,   122,   204,
     273,   205,   275,   203,    -1,   273,   148,   274,    -1,   274,
      -1,   343,    -1,   344,    -1,   345,    -1,   346,    -1,   347,
      -1,   348,    -1,   349,    -1,   350,    -1,   351,    -1,   352,
      -1,   353,    -1,   370,    -1,   354,    -1,   392,    -1,   355,
      -1,   356,    -1,   357,    -1,   358,    -1,   360,    -1,   361,
      -1,   362,    -1,   396,    -1,   397,    -1,   275,    82,    -1,
     275,    82,    33,    82,    -1,   275,   148,    82,    -1,   275,
     148,    82,    33,    82,    -1,    82,    -1,    82,    33,    82,
      -1,   150,    51,    -1,   149,    51,    -1,    51,    -1,   150,
      42,    -1,   149,    42,    -1,    42,    -1,    35,   203,   279,
      31,    -1,   279,   280,    -1,   280,    -1,   281,   148,   282,
     203,    -1,   120,    82,    -1,    82,    -1,    22,    82,   148,
      82,    -1,   290,   148,   283,    -1,   291,   148,   290,   148,
     283,    -1,   291,   148,   291,   148,   291,   148,   290,   148,
     283,    -1,   291,    -1,   291,   148,   291,   148,   291,    -1,
     291,   148,   291,    -1,   291,   148,   291,   148,   291,    -1,
     291,   148,   291,   148,   291,   148,   291,    -1,   291,   148,
     291,   148,   291,   148,   291,   148,   291,    -1,    37,   203,
     285,    31,    -1,   285,   286,    -1,   286,    -1,   120,    82,
     148,   291,   203,    -1,    22,    82,   148,    82,   148,   291,
     203,    -1,    82,   148,   291,   203,    -1,    36,   203,   288,
      31,    -1,   288,   289,    -1,   289,    -1,   120,    82,   148,
     291,   148,   291,   203,    -1,    22,    82,   148,    82,   148,
     291,   148,   291,   203,    -1,    82,   148,   291,   148,   291,
     203,    -1,     6,    -1,    44,    -1,    92,    -1,    52,    -1,
     127,    -1,    -1,    51,    -1,    42,    -1,    82,    -1,   149,
      51,    -1,   149,    42,    -1,    34,   203,    -1,    34,   204,
     293,   205,   203,    -1,    34,   275,   203,    -1,    34,   204,
     293,   205,   275,   203,    -1,   293,   148,   294,    -1,   294,
      -1,   363,    -1,   364,    -1,   365,    -1,   366,    -1,   367,
      -1,   368,    -1,   369,    -1,   370,    -1,   371,    -1,   372,
      -1,   373,    -1,   374,    -1,   375,    -1,   376,    -1,   377,
      -1,   378,    -1,   379,    -1,   380,    -1,   381,    -1,   382,
      -1,   383,    -1,   384,    -1,   385,    -1,   386,    -1,   354,
      -1,   387,    -1,   388,    -1,   389,    -1,   390,    -1,   391,
      -1,   393,    -1,   394,    -1,   398,    -1,   399,    -1,   400,
      -1,   344,    -1,   401,    -1,   402,    -1,   403,    -1,   106,
     204,   296,   205,   203,    -1,   106,   204,   296,   205,   275,
     203,    -1,   296,   148,   297,    -1,   297,    -1,   370,    -1,
     371,    -1,   380,    -1,   386,    -1,   354,    -1,   387,    -1,
     388,    -1,   389,    -1,   390,    -1,   391,    -1,   398,    -1,
     399,    -1,   400,    -1,   107,   204,   296,   205,   203,    -1,
     107,   204,   296,   205,   275,   203,    -1,   210,    82,   210,
     148,   210,    82,   210,    -1,   210,    82,   210,   148,   291,
      -1,   299,    -1,   300,   148,   299,    -1,   134,   275,   203,
      -1,    93,   203,   303,    31,    -1,   303,   304,    -1,   304,
      -1,    82,   204,   231,   205,   203,    -1,   128,   275,   203,
      -1,    95,   203,   307,    31,    -1,   307,    82,   231,   203,
      -1,   307,    82,   148,    82,   231,   203,    -1,    82,   231,
     203,    -1,    82,   148,    82,   231,   203,    -1,    98,   275,
     203,    -1,    97,   203,    -1,    97,   204,   274,   205,   203,
      -1,    97,   275,   203,    -1,    97,   204,   274,   205,   275,
     203,    -1,    18,   203,   311,    31,    -1,   311,   312,    -1,
     312,    -1,    82,   313,    33,   231,   203,    -1,    82,   148,
      82,   313,    33,   231,   203,    -1,     4,    82,   204,    51,
     205,   313,    33,   231,   203,    -1,    -1,   204,    51,   205,
      -1,   204,    42,   205,    -1,    17,   203,    -1,    17,   204,
      23,   205,   203,    -1,    30,   204,    82,   205,   203,    -1,
      30,   204,    82,   205,   275,   203,    -1,    30,    82,   203,
      -1,    30,   204,    82,   211,    82,   205,   203,    -1,    30,
     204,    82,   211,    82,   205,   275,   203,    -1,    30,    82,
     211,    82,   203,    -1,    29,   204,    82,   205,   203,    -1,
      29,   204,    82,   205,   275,   203,    -1,    29,    82,   203,
      -1,    29,   204,    82,   211,    82,   205,   203,    -1,    29,
     204,    82,   211,    82,   205,   275,   203,    -1,    29,    82,
     211,    82,   203,    -1,    77,   204,   318,   205,   320,   203,
      -1,   318,   148,   319,    -1,   319,    -1,   395,    -1,   396,
      -1,   397,    -1,   321,    -1,   320,   148,   321,    -1,   321,
     204,   291,   205,    -1,   320,   148,   321,   204,   291,   205,
      -1,   322,    -1,   321,   322,    -1,    82,    -1,   212,    -1,
     151,    -1,   207,    -1,   211,    -1,    -1,    -1,   101,   324,
     251,   325,   203,    -1,   124,   203,    -1,   124,   204,   327,
     205,   203,    -1,   124,   275,   203,    -1,   124,   204,   327,
     205,   275,   203,    -1,   327,   148,   328,    -1,   328,    -1,
     274,    -1,   404,    -1,   405,    -1,   406,    -1,   407,    -1,
     408,    -1,   409,    -1,   410,    -1,   411,    -1,   329,    -1,
     363,    -1,   398,    -1,   399,    -1,   365,    -1,   367,    -1,
     364,    -1,   366,    -1,   401,    -1,   402,    -1,   330,   148,
     331,    -1,   330,    -1,     7,    51,   203,    -1,     7,   204,
     331,   205,    51,   203,    -1,   330,    -1,   388,    -1,   371,
      -1,   412,    -1,   333,   148,   334,    -1,   333,    -1,     8,
      51,   203,    -1,     8,   204,   334,   205,    51,   203,    -1,
     171,   203,    -1,   171,   204,   337,   205,   203,    -1,   338,
     148,   337,    -1,   338,    -1,   413,    -1,   414,    -1,   415,
      -1,   416,    -1,   417,    -1,   418,    -1,   419,    -1,   420,
      -1,   421,    -1,   422,    -1,   423,    -1,   424,    -1,   425,
      -1,   426,    -1,   427,    -1,   428,    -1,   429,    -1,   430,
      -1,   432,    -1,   433,    -1,   434,    -1,   435,    -1,   436,
      -1,   437,    -1,   438,    -1,   439,    -1,   431,    -1,   138,
     203,   340,    31,    -1,   341,    -1,   340,   341,    -1,    82,
     148,   231,   148,   231,   203,    -1,    51,    -1,    42,    -1,
      26,    33,    51,    -1,   118,    33,    51,    -1,   115,    33,
      51,    -1,    60,    -1,    96,    33,    51,    -1,   110,    33,
      51,    -1,    27,    33,    51,    -1,     3,    33,    51,    -1,
      86,    -1,    88,    -1,    90,    -1,    53,    33,    51,    -1,
      48,    33,    51,    -1,    49,    33,    51,    -1,   100,    33,
      51,    -1,    24,    33,   342,    -1,    63,    33,   342,    -1,
     114,    -1,   116,    33,    51,    -1,   108,    33,   342,    -1,
      25,    33,    82,    -1,    84,    33,   445,    -1,    84,    33,
      51,    -1,    41,    33,    51,    -1,   102,    33,    51,    -1,
     103,    33,    51,    -1,    58,    33,    51,    -1,    59,    33,
      51,    -1,    89,    -1,    46,    -1,    20,    33,   342,    -1,
      70,    33,    51,    -1,    65,    33,   342,    -1,    67,    33,
     342,    -1,    94,    33,   204,   300,   205,    -1,    66,    33,
     342,    -1,    75,    33,    82,    -1,    74,    33,    51,    -1,
      73,    -1,   105,    33,   342,    -1,    68,    33,    51,    -1,
      69,    33,    51,    -1,    61,    -1,    62,    -1,    87,    -1,
       5,    -1,   123,    -1,    43,    33,    51,    -1,   117,    -1,
      81,    -1,    40,    -1,   109,    -1,    54,    33,    51,    -1,
      55,    33,    51,    -1,    79,    33,    56,    -1,    79,    33,
      80,    -1,   104,    -1,    91,    -1,   135,    33,    82,    -1,
     136,    33,   442,    -1,    39,    33,   445,    -1,    21,    -1,
      85,    -1,    71,    -1,   125,    33,   342,    -1,    14,    33,
     277,    -1,     9,    33,   342,    -1,    11,    33,   277,    -1,
      12,    33,   342,    -1,    13,    33,    51,    -1,    10,    -1,
      15,    33,    51,    -1,    16,    33,    51,    -1,   172,    33,
      51,    -1,   173,    33,    51,    -1,   174,    33,    51,    -1,
     175,    33,    51,    -1,   176,    33,    51,    -1,   177,    33,
      51,    -1,   178,    33,    51,    -1,   179,    33,    51,    -1,
     180,    33,    51,    -1,   181,    33,    51,    -1,   182,    33,
      51,    -1,   183,    33,    51,    -1,   184,    33,    51,    -1,
     185,    33,    51,    -1,   186,    33,    51,    -1,   187,    33,
     342,    -1,   188,    33,   342,    -1,   189,    33,    51,    -1,
     190,    33,   445,    -1,   191,    33,   342,    -1,   192,    33,
     342,    -1,   196,    33,    51,    -1,   197,    33,    51,    -1,
     199,    33,   342,    -1,   200,    33,    51,    -1,   201,    33,
     342,    -1,   202,    33,   342,    -1,   139,    33,    51,    -1,
     140,    33,    51,    -1,    82,   207,    82,    -1,    51,    -1,
      51,   207,    51,    -1,   208,   443,    -1,   444,   443,    -1,
     444,   209,    -1
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
     333,   338,   343,   348,   353,   360,   367,   372,   381,   383,
     387,   392,   400,   404,   409,   412,   414,   419,   424,   427,
     429,   437,   441,   443,   445,   447,   449,   451,   453,   454,
     460,   461,   470,   471,   480,   481,   492,   493,   502,   505,
     508,   510,   512,   517,   520,   524,   526,   528,   530,   534,
     538,   542,   546,   550,   554,   558,   562,   566,   570,   574,
     577,   580,   585,   590,   595,   600,   605,   610,   615,   620,
     625,   630,   637,   644,   653,   659,   661,   666,   671,   676,
     679,   681,   691,   698,   704,   712,   720,   724,   730,   734,
     736,   743,   747,   749,   753,   755,   761,   765,   769,   774,
     777,   780,   784,   786,   788,   791,   797,   801,   803,   805,
     807,   809,   812,   818,   822,   824,   826,   829,   835,   839,
     841,   843,   845,   848,   854,   858,   865,   869,   871,   873,
     875,   877,   879,   881,   883,   885,   887,   889,   891,   893,
     895,   897,   899,   901,   903,   905,   907,   909,   911,   913,
     915,   917,   920,   925,   929,   935,   937,   941,   944,   947,
     949,   952,   955,   957,   962,   965,   967,   972,   975,   977,
     982,   986,   992,  1002,  1004,  1010,  1014,  1020,  1028,  1038,
    1043,  1046,  1048,  1054,  1062,  1067,  1072,  1075,  1077,  1085,
    1095,  1102,  1104,  1106,  1108,  1110,  1112,  1113,  1115,  1117,
    1119,  1122,  1125,  1128,  1134,  1138,  1145,  1149,  1151,  1153,
    1155,  1157,  1159,  1161,  1163,  1165,  1167,  1169,  1171,  1173,
    1175,  1177,  1179,  1181,  1183,  1185,  1187,  1189,  1191,  1193,
    1195,  1197,  1199,  1201,  1203,  1205,  1207,  1209,  1211,  1213,
    1215,  1217,  1219,  1221,  1223,  1225,  1227,  1229,  1235,  1242,
    1246,  1248,  1250,  1252,  1254,  1256,  1258,  1260,  1262,  1264,
    1266,  1268,  1270,  1272,  1274,  1280,  1287,  1295,  1301,  1303,
    1307,  1311,  1316,  1319,  1321,  1327,  1331,  1336,  1341,  1348,
    1352,  1358,  1362,  1365,  1371,  1375,  1382,  1387,  1390,  1392,
    1398,  1406,  1416,  1417,  1421,  1425,  1428,  1434,  1440,  1447,
    1451,  1459,  1468,  1474,  1480,  1487,  1491,  1499,  1508,  1514,
    1521,  1525,  1527,  1529,  1531,  1533,  1535,  1539,  1544,  1551,
    1553,  1556,  1558,  1560,  1562,  1564,  1566,  1567,  1568,  1574,
    1577,  1583,  1587,  1594,  1598,  1600,  1602,  1604,  1606,  1608,
    1610,  1612,  1614,  1616,  1618,  1620,  1622,  1624,  1626,  1628,
    1630,  1632,  1634,  1636,  1638,  1642,  1644,  1648,  1655,  1657,
    1659,  1661,  1663,  1667,  1669,  1673,  1680,  1683,  1689,  1693,
    1695,  1697,  1699,  1701,  1703,  1705,  1707,  1709,  1711,  1713,
    1715,  1717,  1719,  1721,  1723,  1725,  1727,  1729,  1731,  1733,
    1735,  1737,  1739,  1741,  1743,  1745,  1747,  1749,  1754,  1756,
    1759,  1766,  1768,  1770,  1774,  1778,  1782,  1784,  1788,  1792,
    1796,  1800,  1802,  1804,  1806,  1810,  1814,  1818,  1822,  1826,
    1830,  1832,  1836,  1840,  1844,  1848,  1852,  1856,  1860,  1864,
    1868,  1872,  1874,  1876,  1880,  1884,  1888,  1892,  1898,  1902,
    1906,  1910,  1912,  1916,  1920,  1924,  1926,  1928,  1930,  1932,
    1934,  1938,  1940,  1942,  1944,  1946,  1950,  1954,  1958,  1962,
    1964,  1966,  1970,  1974,  1978,  1980,  1982,  1984,  1988,  1992,
    1996,  2000,  2004,  2008,  2010,  2014,  2018,  2022,  2026,  2030,
    2034,  2038,  2042,  2046,  2050,  2054,  2058,  2062,  2066,  2070,
    2074,  2078,  2082,  2086,  2090,  2094,  2098,  2102,  2106,  2110,
    2114,  2118,  2122,  2126,  2130,  2134,  2138,  2140,  2144,  2147,
    2150
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    99,    99,   100,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
     119,   120,   121,   122,   123,   124,   125,   126,   127,   128,
     129,   130,   131,   132,   133,   134,   135,   136,   137,   138,
     139,   140,   141,   142,   143,   146,   147,   148,   149,   153,
     155,   159,   161,   163,   165,   167,   169,   171,   173,   175,
     177,   179,   183,   185,   187,   189,   191,   193,   197,   199,
     201,   203,   205,   207,   211,   213,   215,   217,   219,   221,
     225,   227,   231,   233,   237,   239,   244,   246,   248,   250,
     252,   254,   256,   258,   260,   262,   264,   266,   268,   270,
     272,   274,   276,   278,   280,   282,   284,   286,   288,   290,
     292,   294,   296,   298,   300,   302,   304,   306,   310,   312,
     316,   318,   322,   324,   326,   327,   330,   332,   334,   335,
     338,   340,   341,   344,   346,   348,   350,   351,   354,   354,
     356,   356,   358,   358,   361,   360,   363,   363,   367,   368,
     369,   370,   373,   375,   379,   381,   382,   384,   386,   388,
     390,   392,   394,   396,   398,   400,   402,   404,   406,   408,
     410,   412,   414,   416,   418,   420,   422,   424,   426,   428,
     430,   432,   434,   436,   440,   443,   445,   449,   451,   453,
     454,   457,   459,   461,   463,   465,   469,   471,   473,   475,
     480,   483,   485,   489,   491,   495,   497,   499,   501,   503,
     505,   507,   509,   511,   515,   517,   521,   522,   525,   526,
     527,   530,   532,   536,   537,   540,   542,   544,   548,   549,
     552,   553,   556,   558,   560,   562,   566,   567,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
     582,   583,   584,   585,   586,   587,   588,   589,   590,   591,
     592,   595,   597,   599,   601,   603,   605,   609,   611,   613,
     617,   619,   621,   625,   627,   629,   633,   635,   641,   647,
     657,   662,   669,   680,   685,   696,   703,   712,   723,   738,
     741,   743,   747,   755,   765,   775,   778,   780,   784,   794,
     806,   818,   820,   822,   824,   826,   830,   831,   832,   833,
     834,   836,   840,   842,   844,   846,   850,   851,   854,   855,
     856,   857,   858,   859,   860,   861,   862,   863,   864,   865,
     866,   867,   868,   869,   870,   871,   872,   873,   874,   875,
     876,   877,   878,   879,   880,   881,   882,   883,   884,   885,
     886,   887,   888,   889,   890,   891,   892,   895,   897,   901,
     902,   905,   906,   907,   908,   909,   910,   911,   912,   913,
     914,   915,   916,   917,   920,   922,   926,   928,   932,   933,
     936,   938,   940,   941,   944,   946,   948,   950,   952,   954,
     956,   960,   962,   964,   966,   968,   972,   974,   975,   978,
     980,   982,   986,   987,   989,   993,   995,   999,  1001,  1003,
    1005,  1007,  1009,  1013,  1015,  1017,  1019,  1021,  1023,  1027,
    1030,  1031,  1034,  1035,  1036,  1039,  1041,  1043,  1045,  1049,
    1051,  1055,  1056,  1058,  1060,  1062,  1066,  1067,  1066,  1069,
    1071,  1073,  1075,  1079,  1080,  1083,  1084,  1087,  1088,  1089,
    1090,  1091,  1092,  1093,  1096,  1097,  1098,  1099,  1100,  1101,
    1102,  1103,  1104,  1105,  1108,  1109,  1112,  1114,  1118,  1119,
    1120,  1121,  1124,  1125,  1128,  1130,  1134,  1136,  1140,  1141,
    1144,  1145,  1146,  1147,  1148,  1149,  1150,  1151,  1152,  1153,
    1154,  1155,  1156,  1157,  1158,  1159,  1160,  1161,  1162,  1163,
    1164,  1165,  1166,  1167,  1168,  1169,  1170,  1173,  1176,  1177,
    1180,  1183,  1184,  1187,  1188,  1189,  1190,  1191,  1192,  1193,
    1194,  1195,  1196,  1197,  1198,  1199,  1200,  1201,  1203,  1204,
    1205,  1206,  1207,  1208,  1209,  1211,  1214,  1215,  1216,  1217,
    1218,  1219,  1221,  1224,  1225,  1226,  1227,  1228,  1229,  1230,
    1231,  1232,  1233,  1234,  1235,  1236,  1237,  1238,  1239,  1240,
    1241,  1242,  1243,  1244,  1245,  1246,  1247,  1248,  1250,  1253,
    1254,  1255,  1256,  1257,  1258,  1259,  1260,  1261,  1263,  1264,
    1265,  1266,  1267,  1268,  1269,  1270,  1272,  1273,  1274,  1275,
    1276,  1277,  1278,  1279,  1280,  1281,  1282,  1283,  1284,  1285,
    1286,  1287,  1288,  1289,  1290,  1292,  1293,  1299,  1300,  1304,
    1305,  1306,  1307,  1309,  1310,  1312,  1320,  1321,  1330,  1332,
    1341
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
       2,     2,     2,     2,     2,   206,     2,     2,     2,   210,
     204,   205,     2,     2,     2,     2,   211,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   207,   203,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   208,   212,   209,     2,     2,     2,     2,     2,     2,
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
     195,   196,   197,   198,   199,   200,   201,   202
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 2315;
  const int parser::yynnts_ = 233;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 163;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 213;

  const unsigned int parser::yyuser_token_number_max_ = 457;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1343 "DynareBison.yy"


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

