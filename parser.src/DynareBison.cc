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

  case 137:
#line 353 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 138:
#line 355 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 139:
#line 357 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 142:
#line 364 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 143:
#line 365 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 144:
#line 366 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 145:
#line 367 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 146:
#line 368 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 147:
#line 369 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 148:
#line 371 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 149:
#line 372 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 150:
#line 373 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 151:
#line 374 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 152:
#line 376 "DynareBison.yy"
    { driver.begin_model(); driver.sparse(); ;}
    break;

  case 153:
#line 377 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 154:
#line 378 "DynareBison.yy"
    { driver.begin_model(); driver.sparse(); ;}
    break;

  case 155:
#line 379 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 160:
#line 389 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 161:
#line 391 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val)); ;}
    break;

  case 162:
#line 395 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 164:
#line 398 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 165:
#line 400 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 166:
#line 402 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 167:
#line 404 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 168:
#line 406 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 169:
#line 408 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 170:
#line 410 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 171:
#line 412 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 172:
#line 414 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 173:
#line 416 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 174:
#line 418 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 175:
#line 420 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 176:
#line 422 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 177:
#line 424 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 178:
#line 426 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 179:
#line 428 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 180:
#line 430 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 181:
#line 432 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 182:
#line 434 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 183:
#line 436 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 184:
#line 438 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 185:
#line 440 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 186:
#line 442 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 187:
#line 444 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 188:
#line 446 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 189:
#line 448 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 190:
#line 450 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 191:
#line 452 "DynareBison.yy"
    { (yyval.node_val) = driver.add_normcdf((yysemantic_stack_[(8) - (3)].node_val),(yysemantic_stack_[(8) - (5)].node_val),(yysemantic_stack_[(8) - (7)].node_val));;}
    break;

  case 192:
#line 456 "DynareBison.yy"
    { driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 193:
#line 459 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 194:
#line 461 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 195:
#line 464 "DynareBison.yy"
    { driver.end_shocks(); ;}
    break;

  case 196:
#line 466 "DynareBison.yy"
    { driver.end_mshocks(); ;}
    break;

  case 199:
#line 473 "DynareBison.yy"
    { driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val)); ;}
    break;

  case 200:
#line 475 "DynareBison.yy"
    { driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 201:
#line 477 "DynareBison.yy"
    { driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 202:
#line 479 "DynareBison.yy"
    { driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 203:
#line 481 "DynareBison.yy"
    { driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 204:
#line 485 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 205:
#line 487 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 206:
#line 489 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 207:
#line 491 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 208:
#line 493 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 209:
#line 495 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 210:
#line 499 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 211:
#line 503 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 212:
#line 505 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 213:
#line 507 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].node_val));;}
    break;

  case 214:
#line 511 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 215:
#line 513 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 216:
#line 517 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 217:
#line 519 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 218:
#line 521 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 219:
#line 523 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 220:
#line 525 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 221:
#line 527 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 222:
#line 529 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 223:
#line 531 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 224:
#line 533 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 225:
#line 537 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 226:
#line 539 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 232:
#line 552 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 233:
#line 554 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 237:
#line 564 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 238:
#line 566 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 244:
#line 579 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 245:
#line 581 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 246:
#line 583 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 247:
#line 585 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 273:
#line 618 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 274:
#line 620 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 275:
#line 622 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 276:
#line 624 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 277:
#line 626 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 278:
#line 628 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 279:
#line 632 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 280:
#line 634 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 281:
#line 636 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 282:
#line 640 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 283:
#line 642 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 284:
#line 644 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 285:
#line 647 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 286:
#line 650 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 287:
#line 652 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 289:
#line 658 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 290:
#line 664 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 291:
#line 670 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 292:
#line 680 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 293:
#line 685 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 294:
#line 692 "DynareBison.yy"
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

  case 295:
#line 703 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 296:
#line 708 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 297:
#line 719 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 298:
#line 726 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 299:
#line 735 "DynareBison.yy"
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

  case 300:
#line 746 "DynareBison.yy"
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

  case 301:
#line 761 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 302:
#line 764 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 303:
#line 766 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 304:
#line 770 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 305:
#line 778 "DynareBison.yy"
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

  case 306:
#line 788 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 307:
#line 798 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 308:
#line 801 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 309:
#line 803 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 310:
#line 807 "DynareBison.yy"
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

  case 311:
#line 817 "DynareBison.yy"
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

  case 312:
#line 829 "DynareBison.yy"
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

  case 313:
#line 841 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 314:
#line 843 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 315:
#line 845 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 316:
#line 847 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 317:
#line 849 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 318:
#line 852 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 322:
#line 857 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 323:
#line 859 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 324:
#line 863 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 325:
#line 865 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 326:
#line 867 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 327:
#line 869 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 369:
#line 918 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 370:
#line 920 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 386:
#line 943 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 387:
#line 945 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 388:
#line 949 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 389:
#line 951 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 392:
#line 958 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 393:
#line 960 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 396:
#line 966 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 397:
#line 968 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 398:
#line 970 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 399:
#line 973 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 400:
#line 975 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 401:
#line 977 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 402:
#line 979 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 403:
#line 982 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 404:
#line 985 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 405:
#line 987 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 406:
#line 989 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 407:
#line 991 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 408:
#line 994 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 411:
#line 1001 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 412:
#line 1003 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 413:
#line 1005 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 414:
#line 1008 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 415:
#line 1010 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 416:
#line 1012 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 417:
#line 1016 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 418:
#line 1018 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 419:
#line 1022 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 420:
#line 1024 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 421:
#line 1026 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 422:
#line 1028 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 423:
#line 1030 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 424:
#line 1032 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 425:
#line 1036 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 426:
#line 1038 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 427:
#line 1040 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 428:
#line 1042 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 429:
#line 1044 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 430:
#line 1046 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 431:
#line 1050 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 437:
#line 1062 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 438:
#line 1064 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 439:
#line 1066 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 440:
#line 1068 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 441:
#line 1072 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 442:
#line 1074 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 444:
#line 1079 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 445:
#line 1081 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 446:
#line 1083 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 447:
#line 1085 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 448:
#line 1088 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 449:
#line 1089 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 451:
#line 1092 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 452:
#line 1094 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 453:
#line 1096 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 454:
#line 1098 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 478:
#line 1135 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 479:
#line 1137 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 486:
#line 1151 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 487:
#line 1153 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 488:
#line 1157 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 489:
#line 1159 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 519:
#line 1196 "DynareBison.yy"
    { driver.end_homotopy();;}
    break;

  case 522:
#line 1203 "DynareBison.yy"
    { driver.homotopy_val((yysemantic_stack_[(6) - (1)].string_val),(yysemantic_stack_[(6) - (3)].node_val),(yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 525:
#line 1209 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 526:
#line 1210 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 527:
#line 1211 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 528:
#line 1212 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 529:
#line 1213 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1214 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1215 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1216 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1217 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 534:
#line 1218 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 535:
#line 1219 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 536:
#line 1220 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 537:
#line 1221 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 538:
#line 1222 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1224 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 540:
#line 1225 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1226 "DynareBison.yy"
    { driver.option_num("simulation_method",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 542:
#line 1227 "DynareBison.yy"
    { driver.option_num("simulation_method", "0"); ;}
    break;

  case 543:
#line 1228 "DynareBison.yy"
    { driver.option_num("simulation_method", "1"); ;}
    break;

  case 544:
#line 1229 "DynareBison.yy"
    { driver.option_num("simulation_method", "2"); ;}
    break;

  case 545:
#line 1230 "DynareBison.yy"
    { driver.option_num("simulation_method", "3"); ;}
    break;

  case 546:
#line 1231 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 547:
#line 1232 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 548:
#line 1233 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 549:
#line 1234 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 550:
#line 1235 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 551:
#line 1237 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 552:
#line 1239 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 553:
#line 1241 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 554:
#line 1242 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 555:
#line 1243 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 556:
#line 1244 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 557:
#line 1245 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 558:
#line 1247 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 559:
#line 1249 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 560:
#line 1251 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 561:
#line 1252 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 562:
#line 1253 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 563:
#line 1254 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1256 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 566:
#line 1257 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1258 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 568:
#line 1259 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 569:
#line 1260 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 570:
#line 1261 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 571:
#line 1262 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1263 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 573:
#line 1264 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 574:
#line 1265 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 575:
#line 1266 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 576:
#line 1267 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 577:
#line 1268 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1269 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 579:
#line 1270 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 580:
#line 1271 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 581:
#line 1272 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 582:
#line 1273 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1274 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1276 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 585:
#line 1278 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 586:
#line 1280 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 587:
#line 1281 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 588:
#line 1282 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1283 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1284 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1285 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 592:
#line 1286 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 593:
#line 1287 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 594:
#line 1288 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 595:
#line 1290 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 596:
#line 1291 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 597:
#line 1292 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 598:
#line 1293 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 599:
#line 1294 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 600:
#line 1295 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 601:
#line 1296 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 602:
#line 1297 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 603:
#line 1299 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 604:
#line 1300 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 605:
#line 1301 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 606:
#line 1302 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 607:
#line 1303 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 608:
#line 1304 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 609:
#line 1305 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 610:
#line 1306 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 611:
#line 1307 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 612:
#line 1308 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 613:
#line 1309 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 614:
#line 1310 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 615:
#line 1311 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 616:
#line 1312 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 617:
#line 1313 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 618:
#line 1314 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 619:
#line 1315 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 620:
#line 1316 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 621:
#line 1317 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 622:
#line 1319 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 623:
#line 1320 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 624:
#line 1326 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 625:
#line 1327 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 626:
#line 1331 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 627:
#line 1332 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 628:
#line 1333 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 629:
#line 1334 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 630:
#line 1336 "DynareBison.yy"
    {driver.option_num("homotopy_mode",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 631:
#line 1337 "DynareBison.yy"
    {driver.option_num("homotopy_steps",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 632:
#line 1340 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 634:
#line 1349 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
               ;}
    break;

  case 635:
#line 1358 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 636:
#line 1360 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 637:
#line 1368 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2489 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1147;
  const short int
  parser::yypact_[] =
  {
      1129,    17,    19,   195,  -113,   420,   590,    66,   -11,    -2,
     -50,    13,   -36,   -22,   -17,   -15,   431,   592,   470,    -8,
      11,   305,   162,   169,    59,   329,   354,   323, -1147,   182,
     234,   329,   255,   432,   493,   495,    78,    80,   329,   389,
     406,   417,   329,   314,   498,   996, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147,   317,  1598,   324,  1562, -1147,   468,    99, -1147,
     437,   497,   355,    24,   -83,   486,   -51,   525,   547,   603,
   -1147,  1458,     0,    69,   391,   403,   564,   547,   614,   618,
     446, -1147,   282,    27,    51,  1261,   582,   584, -1147,  1705,
      28,   102,   552,   187,   632,   481,  1293,  1311,  1311,   209,
      51,   479, -1147,    79, -1147,   405, -1147,  1705,   256, -1147,
    1622,   274,   286,   567,   288,   578,   289,   579,   292,   293,
     626, -1147,  2198, -1147, -1147, -1147,   680, -1147,   682,   683,
     684,   685,   686, -1147,   690,   692,   693, -1147,   696,   698,
     699,   700, -1147,   581,   527, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147,   705,   706,   707, -1147,   588,   532, -1147, -1147,
   -1147,   534,   661,   -25,    84, -1147,   716,   -60, -1147, -1147,
     542, -1147,   543, -1147, -1147,   666,   283, -1147,   667,   414,
     721,    74, -1147,   670, -1147,   724, -1147, -1147,   732,   737,
     743,   746,   748, -1147, -1147,   749,   750,   751,   752,   753,
     754, -1147, -1147,   755,   756, -1147, -1147, -1147,   758,   761,
   -1147, -1147,   -47, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147,   762,   710, -1147,   711, -1147,   722,   530,
   -1147,   654,   723,   657,   726,   595, -1147,   729,   664,   733,
     596, -1147,   616,   363, -1147,   397,   788,   629,   620, -1147,
    1127, -1147,   -42,   -41,   631,   635,   793, -1147, -1147,   -33,
   -1147, -1147, -1147, -1147,   760,   764,   368, -1147, -1147, -1147,
     633,   639,   642,   644,  1261,  1261,   647,   649,   650,   651,
     653,   655,   658,   659,   660,   663,  1261,   612,   668,   445,
   -1147,   269,   559,   830,   833,   837,   840,   842,   843, -1147,
   -1147, -1147,   848,   849,   850, -1147,   851, -1147,   852,   853,
     677, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147,   759,   803, -1147,
     695, -1147, -1147, -1147,   697,   701,   702,   703,  1293,  1293,
     704,   708,   709,   712,   725,   727,   728,   734,   735,   772,
    1293,  2178, -1147,   -30, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,   -21, -1147,
     390,    55,   859,    -7, -1147, -1147, -1147, -1147,   872,   875,
     186, -1147, -1147, -1147, -1147,   216, -1147, -1147,   876, -1147,
     228, -1147, -1147, -1147, -1147, -1147,   784,   855, -1147, -1147,
     804,   856, -1147, -1147,   807,   871, -1147, -1147,   763,   572,
   -1147,   886,   925,   926,   928,   929,   930,   931,   932,   937,
     938,   949,   950,   952,   953,   955,   956,   957,   958,   959,
     960,   964,   965,   966,   967,   968,   969,   972,   796,   815,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147,   101,    41,   101,
     954,    41,   963,   922,   974,    18,   976,   981,   933,   936,
    1598,   982,   983,   101,   984,  1562,   985,   809,   802,   961,
     310,  1006, -1147, -1147,   987,   437,   810, -1147, -1147,   834,
      62,   970,   835,    76,   986,  1261, -1147, -1147, -1147,   828,
     991,   993,   997,   998,   999,   101,   101,   101,  1000,  1001,
    1017,  1018,   988,   864,   101,  1458,    89,   989,  1016,   927,
   -1147, -1147, -1147,   366,   943,    48,   951, -1147, -1147,   962,
      48,   971, -1147, -1147,   213, -1147, -1147, -1147,  1002,   894,
   -1147,  1019,    50, -1147,   126, -1147,    83, -1147,   597, -1147,
     895,   900,    72,    27,   137,   973,    46, -1147, -1147,  1261,
    1261,  1261,  1261,   920,   490,  1261,  1261,  1261,  1261,  1261,
    1261,  1261,  1261,  1261,  1261,   745,  1261,  1261,  1261,  1261,
    1261,  1261,  1261,  1261,  1261,  1261,  1261, -1147,  1261, -1147,
   -1147,  1023,  1238, -1147,  1210,  1057,   101,  1063,  1064,  1065,
    1067,  1075,  1077,   101,  1078,  1085,  1087,    90, -1147,   990,
   -1147,  1293,  1293,   213,  1293,   992,   508,  1293,  1293,  1293,
    1293,  1293,  1293,  1293,  1293,  1293,  1293,   774,  1293,  1293,
    1293,  1293,  1293,  1293,  1293,  1293,  1293,  1293,  1293,   934,
    1311,    92,    93, -1147, -1147, -1147,  1261,   423,    47,   513,
      79,   941,  1088,  1090,   405,   942,  1705,    96,   101,  1622,
     177, -1147,  1013, -1147,  1021, -1147,  1024,  1261, -1147, -1147,
    1101,  1103,  1107,  1109,  1114,  1115,  1117,  1118,  1120,  1121,
    1122,  1123,  1124,  1126,  1130,   101,   101,  1131,   828,   101,
     101,  1132,  1133,   101,  1134,   101,   101,   980,  2198, -1147,
   -1147, -1147, -1147,  1140,  1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147,  1137,    16, -1147, -1147, -1147, -1147,   979,
   -1147, -1147,   994, -1147, -1147, -1147, -1147,  1004, -1147,  1139,
     995,  1011,  1025,  1261, -1147, -1147, -1147, -1147, -1147,   294,
    1026, -1147, -1147,   297,  1027,  1336, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
     978, -1147, -1147, -1147,   298, -1147,  1110,  1111, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147,   345,  1030,  1045,  1046,
    1119,  1047,    48,  1143,  1033,    48, -1147,  1148,  1153,  1032,
   -1147,   547,  1174, -1147, -1147, -1147,  1293, -1147,  1175,   244,
   -1147, -1147, -1147,  1035, -1147, -1147, -1147,   247, -1147, -1147,
   -1147,  1040, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147,   467,   134, -1147,  1163,  1261,  1164,    14,  2196,
    2211,  2344,   248,  2260,   909,   935,  1068,  1170,  1479,  1508,
    1620,  1632,  1689,  1701, -1147,   517,   517,   517,   517,   517,
     517,   490,   490,   920,   920,  1092,  1713,  1261, -1147,  1167,
    1542, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147,   300, -1147,  2272,  2284,  1044,  2296,
    1725,  1737,  1758,  1770,  1782,  1794,  1806,  1827,  1839,  1851,
   -1147,   540,   540,   540,   540,   540,   540,   508,   508,   992,
     992,  1092, -1147, -1147, -1147,   302, -1147,   303,  1863,    55,
    1049, -1147, -1147,    56,  1261, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,   332, -1147,
   -1147, -1147,   398, -1147, -1147, -1147,  2308, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,  1048, -1147,
   -1147, -1147,  1173, -1147, -1147,  1051,  1231, -1147, -1147,  1989,
   -1147,   271, -1147,   278, -1147,  1184, -1147,   249, -1147, -1147,
   -1147, -1147, -1147, -1147,    48,   366,  1144,    48,  1145,  1146,
   -1147,  1066, -1147, -1147,  1240,   591,  1293,  2001,   101,    83,
   -1147,  1127,   597, -1147,  1127,  1127,  1127,   137, -1147,    48,
   -1147,  1242,  2013,  1243,  1218,  1261,  1261,  1261,  1261, -1147,
    1261, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147,  1071,  2034,  1261, -1147, -1147,  1293,  1293, -1147,  1293,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147,    47, -1147, -1147, -1147,  1261,  1875, -1147,
   -1147,  1261,  1227, -1147,   995,  1261, -1147, -1147,   401, -1147,
     404,  1069,   978, -1147, -1147,  1151,  1154,  1155,    48,  1093,
      48,    48, -1147,  1261, -1147,  2056, -1147, -1147, -1147,  1098,
      85, -1147,  1102,   306,   638,   669,   149,  1073,  1261, -1147,
    1261,  1097,    21,  2068,  1896,  1908,  2344,  2320, -1147, -1147,
    2080,  1920,  1932,  2332,  1944, -1147,  2095, -1147,  1252,  2111,
   -1147, -1147,  1158, -1147,    48,    48,    48,  1159, -1147,  1105,
    1125,  2123, -1147,  1127, -1147,  1127, -1147, -1147, -1147,    48,
   -1147,  2135,  2150,  1247,  1128,  1264,  1188, -1147, -1147, -1147,
    1261, -1147, -1147, -1147,  1293, -1147, -1147,  1261, -1147,    43,
    1177, -1147,  1179,    48, -1147, -1147, -1147,   785,   812,  1135,
   -1147, -1147, -1147,  1281,  1136,  1261,  1965,  1977,  2166,  1253,
   -1147,    48,    77,  1141, -1147, -1147, -1147, -1147,  1288,  2344,
      71, -1147, -1147, -1147,  1138,  1189,  1190, -1147, -1147,  1261,
   -1147, -1147, -1147,    48,    48,  2344,  1191, -1147,    48, -1147
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   448,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     2,     4,    29,    30,
      46,    47,    48,    45,     5,     6,     7,    12,     9,    10,
      11,     8,    13,    14,    15,    16,    17,    18,    19,    23,
      25,    24,    20,    21,    22,    26,    27,    28,    31,    32,
      33,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,     0,     0,     0,     0,   417,     0,     0,   232,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   277,
     324,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   142,     0,     0,     0,     0,     0,     0,   404,     0,
       0,     0,    76,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   237,     0,   225,     0,   244,     0,     0,   451,
       0,     0,     0,    58,     0,    64,     0,    70,     0,     0,
       0,   488,     0,     1,     3,   478,     0,   600,     0,     0,
       0,     0,     0,   591,     0,     0,     0,   592,     0,     0,
       0,     0,   466,   477,     0,   467,   472,   470,   473,   471,
     468,   469,   474,   475,   459,   460,   461,   462,   463,   464,
     465,   486,     0,     0,     0,   480,   485,     0,   482,   481,
     483,     0,     0,   414,     0,   410,     0,     0,   235,   236,
       0,    82,     0,    49,   427,     0,     0,   421,     0,     0,
       0,     0,   125,     0,   575,     0,   580,   559,     0,     0,
       0,     0,     0,   572,   573,     0,     0,     0,     0,     0,
       0,   593,   568,     0,     0,   579,   574,   558,     0,     0,
     578,   576,     0,   329,   365,   354,   330,   331,   332,   333,
     334,   335,   336,   337,   338,   339,   340,   341,   342,   343,
     344,   345,   346,   347,   348,   349,   350,   351,   352,   353,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     366,   367,   368,   273,     0,   326,     0,   290,     0,     0,
     287,     0,     0,     0,     0,     0,   309,     0,     0,     0,
       0,   303,     0,     0,   129,     0,     0,     0,     0,    84,
       0,   528,     0,     0,     0,     0,     0,   587,   586,     0,
     433,   434,   435,   436,     0,     0,     0,   198,    89,    90,
       0,     0,    88,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     395,     0,     0,     0,     0,     0,     0,     0,     0,   533,
     534,   535,     0,     0,     0,   581,     0,   547,     0,     0,
       0,   250,   251,   252,   253,   254,   255,   256,   257,   258,
     259,   260,   262,   264,   265,   266,   267,   268,   269,   270,
     261,   263,   271,   272,   406,   403,    79,    74,     0,    55,
       0,    80,   164,   165,     0,     0,   193,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   449,   163,     0,   372,   377,   373,   374,   375,   376,
     378,   379,   380,   381,   382,   383,   384,   385,     0,    51,
       0,     0,     0,     0,   240,   241,   243,   242,     0,     0,
       0,   228,   229,   230,   231,     0,   249,   246,     0,   457,
       0,   456,   458,   453,   397,    61,    56,     0,    52,    67,
      62,     0,    53,    73,    68,     0,    54,   392,     0,     0,
     520,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   491,
     492,   493,   494,   495,   496,   497,   498,   499,   500,   501,
     502,   503,   504,   505,   506,   507,   508,   509,   518,   510,
     511,   512,   513,   514,   515,   516,   517,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   408,   409,     0,     0,     0,    83,    50,     0,
       0,     0,     0,     0,     0,     0,   123,   124,   278,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   275,     0,
     289,   285,   286,   318,     0,   318,     0,   307,   308,     0,
     318,     0,   301,   302,     0,   127,   128,   120,     0,     0,
      85,     0,     0,   158,     0,   159,     0,   154,     0,   150,
       0,     0,     0,     0,     0,     0,     0,   196,   197,     0,
       0,     0,     0,   102,   103,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    86,     0,   393,
     394,     0,     0,   398,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    77,    75,
      81,     0,     0,     0,     0,   177,   178,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   195,   223,   224,     0,     0,   215,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    59,    57,    65,    63,    71,    69,     0,   519,   521,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   524,
     523,   596,   284,     0,     0,   597,   598,   599,   595,   601,
     550,   553,   552,     0,     0,   551,   554,   555,   588,     0,
     589,   476,     0,   602,   560,   577,   484,     0,   418,     0,
     414,     0,     0,     0,   526,   234,   233,   430,   425,     0,
       0,   424,   419,     0,     0,     0,   590,   536,   582,   583,
     556,   557,   562,   565,   563,   570,   571,   561,   567,   566,
       0,   569,   328,   325,     0,   274,     0,     0,   313,   320,
     314,   319,   316,   321,   315,   317,     0,     0,     0,   295,
       0,     0,   318,     0,     0,   318,   281,     0,     0,     0,
     122,     0,     0,   143,   156,   157,     0,   161,     0,     0,
     134,   140,   141,     0,   138,   137,   139,     0,   132,   135,
     136,     0,   146,   144,   584,   585,   432,   443,   445,   446,
     447,   444,     0,   437,   441,     0,     0,     0,     0,     0,
       0,   118,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    87,   101,   100,    99,    98,    97,
      96,    92,    91,    93,    94,    95,     0,     0,   401,     0,
       0,   532,   540,   525,   531,   537,   538,   529,   539,   549,
     530,   527,   548,   405,     0,    78,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     162,   175,   174,   173,   172,   171,   170,   167,   166,   168,
     169,   176,   450,   371,   369,     0,   386,     0,     0,     0,
       0,   220,   221,     0,     0,   545,   543,   544,   541,   542,
     239,   238,   630,   631,   227,   226,   248,   245,     0,   594,
     455,   452,     0,    60,    66,    72,     0,   603,   604,   605,
     606,   607,   608,   609,   610,   611,   612,   613,   614,   615,
     616,   617,   618,   619,   620,   621,   622,   623,   624,   625,
     626,   627,   628,   629,   489,   490,   283,   282,   633,   635,
     637,   636,     0,   479,   487,     0,     0,   416,   415,     0,
     426,     0,   420,     0,   126,     0,   390,     0,   327,   276,
     291,   323,   322,   288,   318,   318,     0,   318,     0,     0,
     306,     0,   280,   279,     0,     0,     0,     0,     0,     0,
     152,     0,     0,   148,     0,     0,     0,     0,   431,   318,
     442,     0,     0,     0,     0,     0,     0,     0,     0,   116,
       0,   104,   105,   106,   107,   108,   109,   110,   111,   112,
     113,     0,     0,     0,   399,   407,     0,     0,   194,     0,
     179,   180,   181,   182,   183,   184,   185,   186,   187,   188,
     370,   387,   222,   214,   210,   217,   218,     0,     0,   247,
     454,     0,     0,   632,   414,     0,   411,   428,     0,   422,
       0,     0,     0,   564,   292,     0,     0,     0,   318,     0,
     318,   318,   304,     0,   121,     0,   160,   546,   133,     0,
       0,   131,     0,     0,     0,     0,   438,     0,     0,   201,
       0,   209,     0,     0,     0,     0,   119,     0,   396,   402,
       0,     0,     0,     0,     0,   219,     0,   634,     0,     0,
     429,   423,     0,   391,   318,   318,   318,     0,   312,     0,
       0,     0,   192,     0,   155,     0,   151,   147,   145,   318,
     439,     0,     0,     0,   205,     0,     0,   200,   114,   115,
       0,   400,   189,   190,     0,   216,   522,     0,   412,   318,
     297,   293,   296,   318,   310,   305,   130,     0,     0,     0,
     203,   202,   208,     0,   204,     0,     0,     0,     0,     0,
     389,   318,     0,     0,   153,   149,   440,   207,     0,   213,
       0,   117,   191,   413,     0,   298,     0,   311,   206,     0,
     199,   212,   388,   318,   318,   211,   299,   294,   318,   300
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1147, -1147,  1312, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,  -352, -1147,
   -1147, -1147, -1147,  -112,  -230, -1147, -1147,  1036, -1147, -1147,
     246, -1147,  -613, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147,  -973,  -635,  -136,  -589, -1147, -1147, -1147,  1216,  -269,
   -1147, -1147, -1147, -1147,   353, -1147, -1147,   609, -1147, -1147,
     783, -1147, -1147,   619, -1147, -1147,  -120,   -13,   671,   821,
   -1147, -1147,  1074, -1147, -1147, -1146, -1147, -1147,  1072, -1147,
   -1147,  1084, -1029,  -582, -1147, -1147,   786, -1147,  1262,   662,
   -1147,   221, -1147, -1147, -1147, -1147,  1037, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147,  1193,  -766, -1147, -1147, -1147, -1147,
   -1147,   757, -1147,   291,  -858, -1147, -1147, -1147, -1147, -1147,
     652, -1147,   -53,   839, -1147, -1147,   838, -1147, -1147,   617,
   -1147, -1147, -1147,   913,  -551, -1147,   -96, -1147,  1292, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147,  -101, -1147, -1147,  -130,
    -567, -1147, -1147, -1147, -1147, -1147,  -105,   -90,   -85,   -81,
     -79, -1147, -1147,   -91,   -67, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147,   -76, -1147, -1147, -1147, -1147, -1147,   -72,
     -71,   -44,   -66,   -61,   -59, -1147, -1147, -1147, -1147,  -107,
    -106,   -80,   -77,   -52,   -49,   -43, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147, -1147,
   -1147, -1147, -1147, -1147, -1147, -1147,   604, -1147,  -530
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    45,    46,    47,    48,    49,    50,    51,    52,    53,
     154,   156,   158,   133,    54,    55,    56,    57,   367,   932,
      58,   327,    59,   231,   232,    60,   323,   324,   907,   899,
     908,   909,   910,    61,   330,  1116,  1115,  1202,   911,  1199,
     903,   642,   643,   644,   645,   442,    62,    63,   346,   347,
    1212,    64,  1300,   747,   748,    65,   470,   471,    66,   217,
     218,    67,   463,   464,    68,   475,   479,   112,   889,   805,
      69,   309,   310,   311,   877,  1184,    70,   320,   321,    71,
     315,   316,   878,  1185,    72,   262,   263,    73,   443,   444,
      74,  1086,  1087,    75,    76,   369,   370,    77,    78,   372,
      79,    80,    81,   214,   215,   581,    82,    83,    84,    85,
     339,   340,   922,   923,   924,    86,   136,   739,    87,   480,
     481,   182,   183,   184,    88,   206,   207,    89,    90,   528,
     529,    91,   499,   500,   801,   391,   392,   393,   394,   395,
     396,   397,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   466,   902,   407,   408,   409,   185,   186,   187,   188,
     189,   271,   272,   410,   447,   275,   276,   277,   278,   279,
     280,   281,   282,   448,   284,   285,   286,   287,   288,   449,
     450,   451,   452,   453,   454,   411,   295,   296,   341,   412,
     413,   190,   191,   457,   192,   193,   302,   482,   194,   195,
     196,   197,   198,   199,   200,   210,   530,   531,   532,   533,
     534,   535,   536,   537,   538,   539,   540,   541,   542,   543,
     544,   545,   546,   547,   548,   549,   550,   551,   552,   553,
     554,   555,   556,   473,   474,   820,  1069,   814,   815
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       441,   597,   663,   664,   219,   325,   266,   894,   806,   390,
     265,   130,   131,   465,   675,   264,   342,   343,   139,   692,
     273,   267,   824,   148,   151,   152,   268,   476,   208,   159,
     269,   297,   270,   900,   298,   283,   445,   445,   467,   289,
     290,   879,   205,   881,   274,   292,   446,   446,   884,   472,
     293,   209,   294,   895,   852,   853,   854,   455,   455,   299,
     456,   456,   300,   861,  1076,  1120,  1186,   291,   301,   846,
    1068,    92,   812,    94,   344,  1254,   104,   658,   222,   901,
     926,   901,   893,   868,   802,   106,   869,   303,   212,  1271,
    1011,   869,   306,   422,   585,   597,    98,   871,   744,  1165,
     109,  1012,   871,   212,   423,   174,   596,   615,   374,   745,
    1166,   336,   646,   648,   799,   303,   582,  1244,   424,  1124,
     103,   653,   870,   337,   740,   800,   224,   425,   422,   579,
     873,   914,   872,   740,   225,   873,   338,   426,  1200,   423,
    1125,  1203,  1204,  1205,   799,   962,   109,   750,   462,   109,
     898,   586,   969,   424,   304,   800,   307,   915,   227,   108,
     896,   230,   425,   109,   616,   109,   228,   109,  1317,   647,
     649,   213,   426,   113,   874,  1255,   109,   109,   654,   109,
     109,   741,   304,   109,   383,   580,   213,   114,   345,   303,
     742,   658,   115,   427,   116,   308,   803,   804,   876,   105,
     927,  1013,   123,   876,   751,   428,   429,  1029,   107,   305,
     875,   430,   431,   432,   433,   434,   435,   436,   437,   438,
     124,   917,   110,   111,   917,  1309,   439,    93,   427,    95,
    1256,  1070,   813,   223,  1052,  1053,   917,   414,  1056,  1057,
     428,   429,  1060,   845,  1062,  1063,   430,   431,   432,   433,
     434,   435,   436,   437,   438,   928,   304,  1014,  1055,  1289,
     440,   439,   641,  1306,   109,   746,  1167,   886,   128,   129,
    1277,   838,  1278,   728,   417,   729,   730,   731,   732,   733,
    1310,   734,   735,   736,   737,   842,   738,   146,   147,   149,
     150,   918,   715,   716,   918,   440,   303,   641,   863,   973,
    1098,  1004,  1006,  1101,   727,  1027,   918,   929,   930,   931,
     933,   415,   348,   934,   935,   936,   937,   938,   939,   940,
     941,   942,   943,   349,   945,   946,   947,   948,   949,   950,
     951,   952,   953,   954,   955,   897,   956,   350,  1246,   125,
     754,   418,   960,   303,  1119,   331,   351,   919,  1120,   422,
     919,   920,   921,   831,   920,   921,   352,   134,   109,  1249,
     423,   303,   919,   304,   832,   109,   920,   921,   887,   888,
     756,   126,   868,   303,   424,   486,   490,   135,   127,   494,
     303,   303,   759,   425,   303,   303,  1031,   303,  1091,   303,
     303,   344,   137,   426,  1008,   635,   419,   755,  1109,  1092,
     657,  1112,  1128,  1182,    96,    97,   332,   333,  1228,   869,
     304,   870,   353,   344,   312,  1036,   109,   334,   459,   303,
     871,   872,   743,   691,   354,   355,   317,   757,   304,   637,
     356,   357,   358,   359,   360,   361,   362,   363,   364,   760,
     304,   132,   487,   491,   138,   365,   495,   304,   304,   427,
     322,   304,   304,   873,   304,  1110,   304,   304,  1113,  1129,
    1183,   428,   429,   874,   140,   477,   141,   430,   431,   432,
     433,   434,   435,   436,   437,   438,   153,   689,   313,   366,
    1177,  1079,   439,   483,   230,   303,   304,  1179,   303,   219,
     318,   303,   211,   155,   590,   484,  1198,   488,   492,   875,
     591,   496,   497,  1080,   157,   345,  1082,  1088,   208,  1145,
     266,  1160,  1161,  1187,   265,  1189,   440,   314,   641,   264,
    1015,   876,   205,   160,   273,   267,   165,   345,   216,   319,
     268,   209,   368,   201,   269,   297,   270,  1207,   298,   283,
     220,  1169,   901,   289,   290,   901,   342,   343,   274,   292,
     468,   469,   304,   306,   293,   304,   294,  1197,   304,  1016,
     216,  1017,   621,   299,   221,   894,   300,  1018,   894,   894,
     894,   291,   301,   226,  1122,   976,   977,   839,   979,  1019,
     843,   980,   981,   982,   983,   984,   985,   986,   987,   988,
     989,   693,   991,   992,   993,   994,   995,   996,   997,   998,
     999,  1000,  1001,   864,   768,  1142,  1237,  1170,  1239,  1240,
    1230,   895,   229,  1231,   895,   895,   895,   307,   312,   317,
     465,  1117,   374,  1194,   101,   593,   119,   627,   632,    99,
     100,   594,  1009,   102,   230,   120,  1026,   233,  1010,   445,
     117,   118,   894,   894,   904,   467,   694,   684,   685,   446,
     686,   322,  1270,   326,  1272,   329,   308,   905,   472,   498,
     455,   328,  1168,   456,   898,   736,   737,  1279,   738,   368,
    1247,   371,   682,   683,   684,   685,  1118,   686,   230,   121,
     122,   422,   313,   318,   416,   906,   420,  1290,   895,   895,
     421,  1293,   423,   461,   974,   734,   735,   736,   737,   485,
     738,  1248,   142,   143,   144,   145,   424,   161,   162,  1305,
     489,   493,   422,   498,   557,   425,   558,   559,   560,   561,
     562,   314,   319,   423,   563,   426,   564,   565,  1005,  1007,
     566,  1316,   567,   568,   569,   570,  1319,   424,   571,   572,
     573,   574,   575,   576,  1028,   577,   425,  1032,   578,  1311,
     584,   587,   588,   589,   592,   595,   426,   598,   599,   676,
    1107,   677,   678,   679,   680,   681,   600,   682,   683,   684,
     685,   601,   686,  1213,  1214,  1215,  1216,   602,  1217,  1105,
     603,   427,   604,   605,   606,   607,   608,   609,   610,   611,
     612,  1220,   613,   428,   429,   614,   617,   618,   619,   430,
     431,   432,   433,   434,   435,   436,   437,   438,   623,   620,
     624,   625,   427,   626,   439,  1224,   629,  1294,   630,  1226,
     631,   687,   638,  1229,   428,   429,   634,   652,   422,   640,
     430,   431,   432,   433,   434,   435,   436,   437,   438,   423,
     639,  1241,   650,   659,  1295,   439,   651,   655,   440,   660,
     641,   656,   661,   424,   662,   422,  1251,   665,  1252,   666,
     667,   668,   425,   669,   695,   670,   423,   696,   671,   672,
     673,   697,   426,   674,   698,   597,   699,   700,   688,   440,
     424,   641,   701,   702,   703,   704,   705,   706,   707,   425,
     709,   708,   676,   749,   677,   678,   679,   680,   681,   426,
     682,   683,   684,   685,   710,   686,   752,   711,  1286,   753,
     758,   712,   713,   714,   717,  1288,   761,   767,   718,   719,
     770,   728,   720,   729,   730,   731,   732,   733,   427,   734,
     735,   736,   737,  1299,   738,   721,   763,   722,   723,   765,
     428,   429,   762,   764,   724,   725,   430,   431,   432,   433,
     434,   435,   436,   437,   438,   427,   944,  1315,   766,   771,
     772,   439,   773,   774,   775,   776,   777,   428,   429,   798,
    1195,   778,   779,   430,   431,   432,   433,   434,   435,   436,
     437,   438,   726,   780,   781,   990,   782,   783,   439,   784,
     785,   786,   787,   788,   789,   440,   163,   641,   790,   791,
     792,   793,   794,   795,     1,     2,   796,   797,   807,   810,
    1221,  1222,   829,  1223,     3,     4,     5,   809,   828,   836,
     818,     6,   440,   819,   641,     7,     8,     9,   811,    10,
     816,    11,    12,    13,    14,   817,   822,   823,   825,   827,
     833,   834,   813,   837,   841,   847,    15,   848,   830,    16,
     866,   849,   850,   851,   855,   856,   676,   840,   677,   678,
     679,   680,   681,    17,   682,   683,   684,   685,  1178,   686,
    1180,   857,   858,   844,   860,   859,   865,    18,    19,    20,
     686,   867,   676,    21,   677,   678,   679,   680,   681,   890,
     682,   683,   684,   685,    22,   686,    23,   880,    24,    25,
      26,    27,    28,   891,   912,   882,   892,    29,    30,   913,
     957,   961,    31,    32,    33,    34,   883,   963,   964,   965,
    1131,   966,   975,    35,    36,   885,    37,   925,  1287,   967,
      38,   968,   970,    39,    40,    41,    42,     1,     2,   971,
      43,   972,  1022,  1002,  1023,  1033,  1132,     3,     4,     5,
    1021,  1025,   738,  1034,     6,  1037,  1035,  1038,     7,     8,
       9,  1039,    10,  1040,    11,    12,    13,    14,  1041,  1042,
     422,  1043,  1044,    44,  1045,  1046,  1047,  1048,  1049,    15,
    1050,   423,    16,  1066,  1051,  1054,  1058,  1059,  1061,  1064,
    1067,  1068,  1072,  1075,  1085,   424,    17,  1089,  1090,  1094,
    1095,  1097,  1102,  1073,   425,   580,  1096,  1103,  1106,  1108,
      18,    19,    20,  1074,   426,   676,    21,   677,   678,   679,
     680,   681,  1077,   682,   683,   684,   685,    22,   686,    23,
    1099,    24,    25,    26,    27,    28,  1078,  1081,  1083,  1093,
      29,    30,  1100,  1104,  1111,    31,    32,    33,    34,  1114,
    1121,  1123,    -1,   348,  1143,  1148,    35,    36,  1164,    37,
    1173,  1172,  1174,    38,   349,  1175,    39,    40,    41,    42,
     427,  1181,  1211,    43,  1193,  1192,  1208,  1210,   350,  1133,
    1218,  1227,   428,   429,  1250,  1232,  1267,   351,   430,   431,
     432,   433,   434,   435,   436,   437,   438,   352,  1188,  1190,
    1191,  1282,  1238,   439,   348,  1234,    44,  1243,  1235,  1236,
    1253,  1245,  1269,  1273,  1274,   349,   234,   676,  1284,   677,
     678,   679,   680,   681,  1285,   682,   683,   684,   685,   350,
     686,  1291,   203,  1292,  1275,  1297,   422,   440,   351,   641,
    1304,  1283,  1308,  1313,  1314,  1318,  1296,   423,   352,  1298,
    1307,   235,   236,   353,  1312,   204,   460,   164,  1201,   636,
     237,   424,  1163,  1024,   959,   354,   355,   238,   835,  1020,
     425,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     426,  1134,   808,   622,   978,   676,   365,   677,   678,   679,
     680,   681,   633,   682,   683,   684,   685,   255,   686,   628,
     458,   862,  1003,  1233,   353,   257,   690,   583,  1206,   821,
     916,  1030,   769,   826,   335,  1065,   354,   355,  1071,     0,
     366,   259,   356,   357,   358,   359,   360,   361,   362,   363,
     364,     0,     0,   260,     0,     0,   427,   365,     0,     0,
     261,     0,     0,     0,     0,     0,     0,   958,   428,   429,
       0,     0,   180,   181,   430,   431,   432,   433,   434,   435,
     436,   437,   438,   234,     0,     0,     0,     0,     0,   439,
       0,   366,     0,     0,     0,     0,     0,     0,     0,   203,
     173,     0,     0,   676,   174,   677,   678,   679,   680,   681,
       0,   682,   683,   684,   685,     0,   686,     0,   235,   236,
     175,     0,   204,   440,     0,     0,     0,   237,     0,     0,
       0,     0,     0,     0,   238,   239,   240,     0,     0,   241,
     242,     0,   243,   244,     0,     0,     0,     0,   245,   246,
     247,   248,   249,   250,   251,     0,   252,   253,   254,     0,
       0,     0,     0,     0,   255,  1084,     0,   176,   177,     0,
     256,     0,   257,     0,     0,     0,     0,   258,     0,     0,
       0,     0,     0,     0,     0,   178,   179,     0,   259,     0,
       0,     0,   166,   167,   168,   169,   170,   171,   172,   202,
     260,   216,     0,   203,   173,     0,     0,   261,   174,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   180,
     181,     0,     0,     0,   175,     0,   204,     0,   166,   167,
     168,   169,   170,   171,   172,     0,     0,     0,     0,     0,
     173,     0,     0,     0,   174,   373,   676,     0,   677,   678,
     679,   680,   681,     0,   682,   683,   684,   685,     0,   686,
     175,     0,     0,     0,     0,     0,     0,   374,     0,   375,
     376,   176,   177,     0,     0,   676,     0,   677,   678,   679,
     680,   681,     0,   682,   683,   684,   685,     0,   686,   178,
     179,   237,     0,   377,   378,     0,     0,     0,   238,     0,
       0,     0,     0,     0,     0,   331,     0,   176,   177,   676,
    1135,   677,   678,   679,   680,   681,     0,   682,   683,   684,
     685,     0,   686,   180,   181,   178,   179,     0,   373,     0,
       0,     0,     0,   379,     0,   380,   257,   381,   337,  1136,
       0,     0,     0,   382,     0,     0,     0,   383,     0,     0,
     374,   338,   375,   376,     0,   384,   385,   386,     0,   180,
     181,   387,   388,   389,     0,   216,     0,     0,     0,     0,
       0,  1144,     0,   478,   237,     0,   377,   378,     0,     0,
       0,   238,     0,     0,     0,     0,     0,   676,   331,   677,
     678,   679,   680,   681,     0,   682,   683,   684,   685,   676,
     686,   677,   678,   679,   680,   681,     0,   682,   683,   684,
     685,     0,   686,     0,     0,     0,   379,     0,   380,   257,
     381,   337,     0,     0,     0,     0,   382,     0,     0,     0,
     383,     0,     0,     0,   338,     0,     0,     0,   384,   385,
     386,     0,     0,     0,   387,   388,   389,     0,   216,     0,
       0,  1137,     0,     0,     0,     0,   676,     0,   677,   678,
     679,   680,   681,  1138,   682,   683,   684,   685,   676,   686,
     677,   678,   679,   680,   681,     0,   682,   683,   684,   685,
     676,   686,   677,   678,   679,   680,   681,     0,   682,   683,
     684,   685,   728,   686,   729,   730,   731,   732,   733,     0,
     734,   735,   736,   737,   728,   738,   729,   730,   731,   732,
     733,     0,   734,   735,   736,   737,     0,   738,     0,     0,
    1139,     0,     0,     0,     0,   728,     0,   729,   730,   731,
     732,   733,  1140,   734,   735,   736,   737,   728,   738,   729,
     730,   731,   732,   733,  1141,   734,   735,   736,   737,   728,
     738,   729,   730,   731,   732,   733,  1150,   734,   735,   736,
     737,   728,   738,   729,   730,   731,   732,   733,  1151,   734,
     735,   736,   737,   728,   738,   729,   730,   731,   732,   733,
       0,   734,   735,   736,   737,     0,   738,     0,     0,  1152,
       0,     0,     0,     0,   728,     0,   729,   730,   731,   732,
     733,  1153,   734,   735,   736,   737,   728,   738,   729,   730,
     731,   732,   733,  1154,   734,   735,   736,   737,   728,   738,
     729,   730,   731,   732,   733,  1155,   734,   735,   736,   737,
     676,   738,   677,   678,   679,   680,   681,  1156,   682,   683,
     684,   685,   676,   686,   677,   678,   679,   680,   681,     0,
     682,   683,   684,   685,     0,   686,     0,     0,  1157,     0,
       0,     0,     0,   676,     0,   677,   678,   679,   680,   681,
    1158,   682,   683,   684,   685,   676,   686,   677,   678,   679,
     680,   681,  1159,   682,   683,   684,   685,   728,   686,   729,
     730,   731,   732,   733,  1162,   734,   735,   736,   737,   728,
     738,   729,   730,   731,   732,   733,  1225,   734,   735,   736,
     737,   676,   738,   677,   678,   679,   680,   681,     0,   682,
     683,   684,   685,     0,   686,     0,     0,  1258,     0,     0,
       0,     0,   676,     0,   677,   678,   679,   680,   681,  1259,
     682,   683,   684,   685,   728,   686,   729,   730,   731,   732,
     733,  1262,   734,   735,   736,   737,   676,   738,   677,   678,
     679,   680,   681,  1263,   682,   683,   684,   685,   728,   686,
     729,   730,   731,   732,   733,  1265,   734,   735,   736,   737,
     676,   738,   677,   678,   679,   680,   681,     0,   682,   683,
     684,   685,     0,   686,     0,     0,  1301,     0,     0,     0,
       0,   676,     0,   677,   678,   679,   680,   681,  1302,   682,
     683,   684,   685,     0,   686,     0,     0,     0,  1176,     0,
       0,     0,     0,   728,     0,   729,   730,   731,   732,   733,
    1196,   734,   735,   736,   737,   676,   738,   677,   678,   679,
     680,   681,  1209,   682,   683,   684,   685,   676,   686,   677,
     678,   679,   680,   681,     0,   682,   683,   684,   685,     0,
     686,     0,   676,  1219,   677,   678,   679,   680,   681,     0,
     682,   683,   684,   685,     0,   686,     0,     0,   676,     0,
     677,   678,   679,   680,   681,  1242,   682,   683,   684,   685,
     676,   686,   677,   678,   679,   680,   681,  1257,   682,   683,
     684,   685,   676,   686,   677,   678,   679,   680,   681,  1261,
     682,   683,   684,   685,     0,   686,     0,   676,     0,   677,
     678,   679,   680,   681,  1266,   682,   683,   684,   685,     0,
     686,     0,     0,   676,     0,   677,   678,   679,   680,   681,
    1268,   682,   683,   684,   685,   728,   686,   729,   730,   731,
     732,   733,  1276,   734,   735,   736,   737,     0,   738,     0,
       0,     0,     0,   676,  1280,   677,   678,   679,   680,   681,
    1126,   682,   683,   684,   685,     0,   686,     0,   676,  1281,
     677,   678,   679,   680,   681,  1127,   682,   683,   684,   685,
       0,   686,     0,     0,     0,  1303,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,     0,     0,     0,
     522,   523,     0,   524,   525,   526,   527,   676,     0,   677,
     678,   679,   680,   681,  1130,   682,   683,   684,   685,   728,
     686,   729,   730,   731,   732,   733,  1146,   734,   735,   736,
     737,   728,   738,   729,   730,   731,   732,   733,  1147,   734,
     735,   736,   737,   728,   738,   729,   730,   731,   732,   733,
    1149,   734,   735,   736,   737,   676,   738,   677,   678,   679,
     680,   681,  1171,   682,   683,   684,   685,   676,   686,   677,
     678,   679,   680,   681,  1260,   682,   683,   684,   685,   728,
     686,   729,   730,   731,   732,   733,  1264,   734,   735,   736,
     737,   676,   738,   677,   678,   679,   680,   681,     0,   682,
     683,   684,   685,     0,   686
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
       136,   231,   354,   355,   100,   117,   111,   642,   559,   129,
     111,    24,    25,   143,   366,   111,   123,   123,    31,   371,
     111,   111,   573,    36,    37,    38,   111,   147,    95,    42,
     111,   111,   111,   646,   111,   111,   137,   138,   143,   111,
     111,   623,    95,   625,   111,   111,   137,   138,   630,   145,
     111,    95,   111,   642,   605,   606,   607,   137,   138,   111,
     137,   138,   111,   614,   830,   923,  1095,   111,   111,   599,
      54,    54,    54,    54,    23,    54,    87,   346,    54,   646,
      34,   648,    32,     6,    43,    87,    43,    87,     4,  1235,
      43,    43,    23,    43,   154,   325,   209,    54,    43,    43,
      87,    54,    54,     4,    54,    26,    32,   154,    25,    54,
      54,    84,   154,   154,    43,    87,    32,    32,    68,   105,
      54,   154,    45,    96,   154,    54,   209,    77,    43,   154,
      87,    59,    55,   154,   217,    87,   109,    87,  1111,    54,
     126,  1114,  1115,  1116,    43,   696,    87,   154,    69,    87,
      67,   211,   703,    68,   154,    54,    87,    85,   209,   209,
      34,    87,    77,    87,   211,    87,   217,    87,  1314,   211,
     211,    87,    87,   209,    97,   154,    87,    87,   211,    87,
      87,   211,   154,    87,   105,   210,    87,   209,   137,    87,
     211,   460,   209,   143,   209,   126,   155,   156,   155,   210,
     154,   154,   210,   155,   211,   155,   156,   758,   210,   209,
     133,   161,   162,   163,   164,   165,   166,   167,   168,   169,
     209,    87,   209,   210,    87,   154,   176,   210,   143,   210,
     209,   215,   214,   209,   785,   786,    87,   209,   789,   790,
     155,   156,   793,   595,   795,   796,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   209,   154,   210,   788,   216,
     210,   176,   212,  1292,    87,   210,   210,    54,   209,   210,
    1243,   209,  1245,   147,    87,   149,   150,   151,   152,   153,
     209,   155,   156,   157,   158,   209,   160,   209,   210,   209,
     210,   157,   428,   429,   157,   210,    87,   212,   209,   209,
     882,   209,   209,   885,   440,   209,   157,   659,   660,   661,
     662,   209,    43,   665,   666,   667,   668,   669,   670,   671,
     672,   673,   674,    54,   676,   677,   678,   679,   680,   681,
     682,   683,   684,   685,   686,   209,   688,    68,    32,    34,
     154,   154,   694,    87,   210,    63,    77,   213,  1206,    43,
     213,   217,   218,    43,   217,   218,    87,    34,    87,   210,
      54,    87,   213,   154,    54,    87,   217,   218,   155,   156,
     154,   209,     6,    87,    68,    87,    87,    54,   209,    87,
      87,    87,   154,    77,    87,    87,   209,    87,    43,    87,
      87,    23,   210,    87,   746,    32,   209,   211,   154,    54,
      32,   154,   154,   154,   209,   210,   124,   125,  1174,    43,
     154,    45,   143,    23,    23,   767,    87,   135,   209,    87,
      54,    55,    32,   154,   155,   156,    23,   211,   154,    32,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   211,
     154,    87,   154,   154,   210,   176,   154,   154,   154,   143,
      87,   154,   154,    87,   154,   211,   154,   154,   211,   211,
     211,   155,   156,    97,   209,   209,    34,   161,   162,   163,
     164,   165,   166,   167,   168,   169,    87,    32,    87,   210,
     209,   833,   176,   209,    87,    87,   154,   209,    87,   585,
      87,    87,    24,    87,   211,   209,  1109,   209,   209,   133,
     217,   209,   209,   209,    87,   137,   209,   209,   575,   209,
     615,   209,   209,  1095,   615,  1097,   210,   126,   212,   615,
       7,   155,   575,   209,   615,   615,   209,   137,   123,   126,
     615,   575,    87,   209,   615,   615,   615,  1119,   615,   615,
      43,   209,  1109,   615,   615,  1112,   653,   653,   615,   615,
     145,   146,   154,    23,   615,   154,   615,  1108,   154,    46,
     123,    48,    32,   615,   209,  1200,   615,    54,  1203,  1204,
    1205,   615,   615,    87,   926,   711,   712,   590,   714,    66,
     593,   717,   718,   719,   720,   721,   722,   723,   724,   725,
     726,    32,   728,   729,   730,   731,   732,   733,   734,   735,
     736,   737,   738,   616,    32,   957,  1188,   209,  1190,  1191,
     209,  1200,    87,   209,  1203,  1204,  1205,    87,    23,    23,
     750,   154,    25,    32,    34,   211,    34,    32,    32,   209,
     210,   217,   209,    43,    87,    43,   756,    34,   215,   740,
     209,   210,  1277,  1278,    47,   750,    87,   157,   158,   740,
     160,    87,  1234,    39,  1236,   209,   126,    60,   754,    87,
     740,    43,  1014,   740,    67,   157,   158,  1249,   160,    87,
      32,    87,   155,   156,   157,   158,   209,   160,    87,   209,
     210,    43,    87,    87,   132,    88,    54,  1269,  1277,  1278,
     209,  1273,    54,   214,   707,   155,   156,   157,   158,   132,
     160,    32,   209,   210,   209,   210,    68,   209,   210,  1291,
     132,   132,    43,    87,    34,    77,    34,    34,    34,    34,
      34,   126,   126,    54,    34,    87,    34,    34,   741,   742,
      34,  1313,    34,    34,    34,   154,  1318,    68,   211,    34,
      34,    34,   154,   211,   757,   211,    77,   760,    87,  1300,
      34,   209,   209,    87,    87,    34,    87,    87,    34,   147,
     896,   149,   150,   151,   152,   153,    34,   155,   156,   157,
     158,    34,   160,  1125,  1126,  1127,  1128,    34,  1130,   891,
      34,   143,    34,    34,    34,    34,    34,    34,    34,    34,
      34,  1143,    34,   155,   156,    34,    34,    87,    87,   161,
     162,   163,   164,   165,   166,   167,   168,   169,   154,    87,
      87,   154,   143,    87,   176,  1167,    87,    32,   154,  1171,
      87,   209,    34,  1175,   155,   156,   210,    34,    43,   209,
     161,   162,   163,   164,   165,   166,   167,   168,   169,    54,
     211,  1193,   211,   210,    32,   176,   211,    87,   210,   210,
     212,    87,   210,    68,   210,    43,  1208,   210,  1210,   210,
     210,   210,    77,   210,    34,   210,    54,    34,   210,   210,
     210,    34,    87,   210,    34,  1105,    34,    34,   210,   210,
      68,   212,    34,    34,    34,    34,    34,    34,   211,    77,
      87,   132,   147,    34,   149,   150,   151,   152,   153,    87,
     155,   156,   157,   158,   209,   160,    34,   210,  1260,    34,
      34,   210,   210,   210,   210,  1267,   132,   154,   210,   210,
      34,   147,   210,   149,   150,   151,   152,   153,   143,   155,
     156,   157,   158,  1285,   160,   210,   132,   210,   210,   132,
     155,   156,    87,    87,   210,   210,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   143,   211,  1309,    87,    34,
      34,   176,    34,    34,    34,    34,    34,   155,   156,   154,
    1106,    34,    34,   161,   162,   163,   164,   165,   166,   167,
     168,   169,   210,    34,    34,   211,    34,    34,   176,    34,
      34,    34,    34,    34,    34,   210,     0,   212,    34,    34,
      34,    34,    34,    34,     8,     9,    34,   211,    54,    87,
    1146,  1147,   210,  1149,    18,    19,    20,    54,   209,   209,
      87,    25,   210,    87,   212,    29,    30,    31,    54,    33,
      54,    35,    36,    37,    38,    54,    54,    54,    54,    54,
      34,    54,   214,   209,   209,    54,    50,    54,    87,    53,
      34,    54,    54,    54,    54,    54,   147,    87,   149,   150,
     151,   152,   153,    67,   155,   156,   157,   158,  1081,   160,
    1083,    54,    54,    87,   210,    87,    87,    81,    82,    83,
     160,   154,   147,    87,   149,   150,   151,   152,   153,    87,
     155,   156,   157,   158,    98,   160,   100,   154,   102,   103,
     104,   105,   106,   209,   209,   154,    87,   111,   112,   209,
      87,    54,   116,   117,   118,   119,   154,    54,    54,    54,
     211,    54,   132,   127,   128,   154,   130,   154,  1264,    54,
     134,    54,    54,   137,   138,   139,   140,     8,     9,    54,
     144,    54,    54,   209,    54,   132,   211,    18,    19,    20,
     209,   209,   160,   132,    25,    54,   132,    54,    29,    30,
      31,    54,    33,    54,    35,    36,    37,    38,    54,    54,
      43,    54,    54,   177,    54,    54,    54,    54,    54,    50,
      54,    54,    53,    43,    54,    54,    54,    54,    54,   209,
      43,    54,   213,    54,   216,    68,    67,    87,    87,   154,
     154,   154,    54,   209,    77,   210,    87,    54,    34,    34,
      81,    82,    83,   209,    87,   147,    87,   149,   150,   151,
     152,   153,   211,   155,   156,   157,   158,    98,   160,   100,
      87,   102,   103,   104,   105,   106,   211,   211,   211,   209,
     111,   112,   209,   211,   209,   116,   117,   118,   119,   209,
      87,    87,   160,    43,    87,   211,   127,   128,   209,   130,
      87,   213,   211,   134,    54,    34,   137,   138,   139,   140,
     143,    87,    54,   144,    34,   209,    34,    34,    68,   211,
     209,    54,   155,   156,   211,   216,    34,    77,   161,   162,
     163,   164,   165,   166,   167,   168,   169,    87,   154,   154,
     154,    54,   209,   176,    43,   154,   177,   209,   154,   154,
     213,   209,   154,   154,   209,    54,     5,   147,    54,   149,
     150,   151,   152,   153,   136,   155,   156,   157,   158,    68,
     160,   154,    21,   154,   209,    54,    43,   210,    77,   212,
      87,   213,    54,   154,   154,   154,   211,    54,    87,   213,
     209,    40,    41,   143,   216,    44,   140,    45,  1112,   323,
      49,    68,  1009,   754,   154,   155,   156,    56,   585,   750,
      77,   161,   162,   163,   164,   165,   166,   167,   168,   169,
      87,   211,   561,   309,   713,   147,   176,   149,   150,   151,
     152,   153,   320,   155,   156,   157,   158,    86,   160,   315,
     138,   615,   740,  1182,   143,    94,   369,   214,  1117,   570,
     653,   759,   499,   575,   122,   798,   155,   156,   814,    -1,
     210,   110,   161,   162,   163,   164,   165,   166,   167,   168,
     169,    -1,    -1,   122,    -1,    -1,   143,   176,    -1,    -1,
     129,    -1,    -1,    -1,    -1,    -1,    -1,   209,   155,   156,
      -1,    -1,   141,   142,   161,   162,   163,   164,   165,   166,
     167,   168,   169,     5,    -1,    -1,    -1,    -1,    -1,   176,
      -1,   210,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    21,
      22,    -1,    -1,   147,    26,   149,   150,   151,   152,   153,
      -1,   155,   156,   157,   158,    -1,   160,    -1,    40,    41,
      42,    -1,    44,   210,    -1,    -1,    -1,    49,    -1,    -1,
      -1,    -1,    -1,    -1,    56,    57,    58,    -1,    -1,    61,
      62,    -1,    64,    65,    -1,    -1,    -1,    -1,    70,    71,
      72,    73,    74,    75,    76,    -1,    78,    79,    80,    -1,
      -1,    -1,    -1,    -1,    86,   209,    -1,    89,    90,    -1,
      92,    -1,    94,    -1,    -1,    -1,    -1,    99,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   107,   108,    -1,   110,    -1,
      -1,    -1,    10,    11,    12,    13,    14,    15,    16,    17,
     122,   123,    -1,    21,    22,    -1,    -1,   129,    26,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   141,
     142,    -1,    -1,    -1,    42,    -1,    44,    -1,    10,    11,
      12,    13,    14,    15,    16,    -1,    -1,    -1,    -1,    -1,
      22,    -1,    -1,    -1,    26,     3,   147,    -1,   149,   150,
     151,   152,   153,    -1,   155,   156,   157,   158,    -1,   160,
      42,    -1,    -1,    -1,    -1,    -1,    -1,    25,    -1,    27,
      28,    89,    90,    -1,    -1,   147,    -1,   149,   150,   151,
     152,   153,    -1,   155,   156,   157,   158,    -1,   160,   107,
     108,    49,    -1,    51,    52,    -1,    -1,    -1,    56,    -1,
      -1,    -1,    -1,    -1,    -1,    63,    -1,    89,    90,   147,
     211,   149,   150,   151,   152,   153,    -1,   155,   156,   157,
     158,    -1,   160,   141,   142,   107,   108,    -1,     3,    -1,
      -1,    -1,    -1,    91,    -1,    93,    94,    95,    96,   211,
      -1,    -1,    -1,   101,    -1,    -1,    -1,   105,    -1,    -1,
      25,   109,    27,    28,    -1,   113,   114,   115,    -1,   141,
     142,   119,   120,   121,    -1,   123,    -1,    -1,    -1,    -1,
      -1,   209,    -1,   131,    49,    -1,    51,    52,    -1,    -1,
      -1,    56,    -1,    -1,    -1,    -1,    -1,   147,    63,   149,
     150,   151,   152,   153,    -1,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,    -1,   155,   156,   157,
     158,    -1,   160,    -1,    -1,    -1,    91,    -1,    93,    94,
      95,    96,    -1,    -1,    -1,    -1,   101,    -1,    -1,    -1,
     105,    -1,    -1,    -1,   109,    -1,    -1,    -1,   113,   114,
     115,    -1,    -1,    -1,   119,   120,   121,    -1,   123,    -1,
      -1,   211,    -1,    -1,    -1,    -1,   147,    -1,   149,   150,
     151,   152,   153,   211,   155,   156,   157,   158,   147,   160,
     149,   150,   151,   152,   153,    -1,   155,   156,   157,   158,
     147,   160,   149,   150,   151,   152,   153,    -1,   155,   156,
     157,   158,   147,   160,   149,   150,   151,   152,   153,    -1,
     155,   156,   157,   158,   147,   160,   149,   150,   151,   152,
     153,    -1,   155,   156,   157,   158,    -1,   160,    -1,    -1,
     211,    -1,    -1,    -1,    -1,   147,    -1,   149,   150,   151,
     152,   153,   211,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,   211,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,   211,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,   211,   155,
     156,   157,   158,   147,   160,   149,   150,   151,   152,   153,
      -1,   155,   156,   157,   158,    -1,   160,    -1,    -1,   211,
      -1,    -1,    -1,    -1,   147,    -1,   149,   150,   151,   152,
     153,   211,   155,   156,   157,   158,   147,   160,   149,   150,
     151,   152,   153,   211,   155,   156,   157,   158,   147,   160,
     149,   150,   151,   152,   153,   211,   155,   156,   157,   158,
     147,   160,   149,   150,   151,   152,   153,   211,   155,   156,
     157,   158,   147,   160,   149,   150,   151,   152,   153,    -1,
     155,   156,   157,   158,    -1,   160,    -1,    -1,   211,    -1,
      -1,    -1,    -1,   147,    -1,   149,   150,   151,   152,   153,
     211,   155,   156,   157,   158,   147,   160,   149,   150,   151,
     152,   153,   211,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,   211,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,   211,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,    -1,   155,
     156,   157,   158,    -1,   160,    -1,    -1,   211,    -1,    -1,
      -1,    -1,   147,    -1,   149,   150,   151,   152,   153,   211,
     155,   156,   157,   158,   147,   160,   149,   150,   151,   152,
     153,   211,   155,   156,   157,   158,   147,   160,   149,   150,
     151,   152,   153,   211,   155,   156,   157,   158,   147,   160,
     149,   150,   151,   152,   153,   211,   155,   156,   157,   158,
     147,   160,   149,   150,   151,   152,   153,    -1,   155,   156,
     157,   158,    -1,   160,    -1,    -1,   211,    -1,    -1,    -1,
      -1,   147,    -1,   149,   150,   151,   152,   153,   211,   155,
     156,   157,   158,    -1,   160,    -1,    -1,    -1,   209,    -1,
      -1,    -1,    -1,   147,    -1,   149,   150,   151,   152,   153,
     209,   155,   156,   157,   158,   147,   160,   149,   150,   151,
     152,   153,   209,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,    -1,   155,   156,   157,   158,    -1,
     160,    -1,   147,   209,   149,   150,   151,   152,   153,    -1,
     155,   156,   157,   158,    -1,   160,    -1,    -1,   147,    -1,
     149,   150,   151,   152,   153,   209,   155,   156,   157,   158,
     147,   160,   149,   150,   151,   152,   153,   209,   155,   156,
     157,   158,   147,   160,   149,   150,   151,   152,   153,   209,
     155,   156,   157,   158,    -1,   160,    -1,   147,    -1,   149,
     150,   151,   152,   153,   209,   155,   156,   157,   158,    -1,
     160,    -1,    -1,   147,    -1,   149,   150,   151,   152,   153,
     209,   155,   156,   157,   158,   147,   160,   149,   150,   151,
     152,   153,   209,   155,   156,   157,   158,    -1,   160,    -1,
      -1,    -1,    -1,   147,   209,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,    -1,   160,    -1,   147,   209,
     149,   150,   151,   152,   153,   154,   155,   156,   157,   158,
      -1,   160,    -1,    -1,    -1,   209,   178,   179,   180,   181,
     182,   183,   184,   185,   186,   187,   188,   189,   190,   191,
     192,   193,   194,   195,   196,   197,   198,    -1,    -1,    -1,
     202,   203,    -1,   205,   206,   207,   208,   147,    -1,   149,
     150,   151,   152,   153,   154,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,   154,   155,
     156,   157,   158,   147,   160,   149,   150,   151,   152,   153,
     154,   155,   156,   157,   158,   147,   160,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,   147,   160,   149,
     150,   151,   152,   153,   154,   155,   156,   157,   158,   147,
     160,   149,   150,   151,   152,   153,   154,   155,   156,   157,
     158,   147,   160,   149,   150,   151,   152,   153,    -1,   155,
     156,   157,   158,    -1,   160
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
     138,   139,   140,   144,   177,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   233,   234,   235,   236,   239,   241,
     244,   252,   265,   266,   270,   274,   277,   280,   283,   289,
     295,   298,   303,   306,   309,   312,   313,   316,   317,   319,
     320,   321,   325,   326,   327,   328,   334,   337,   343,   346,
     347,   350,    54,   210,    54,   210,   209,   210,   209,   209,
     210,    34,    43,    54,    87,   210,    87,   210,   209,    87,
     209,   210,   286,   209,   209,   209,   209,   209,   210,    34,
      43,   209,   210,   210,   209,    34,   209,   209,   209,   210,
     286,   286,    87,   232,    34,    54,   335,   210,   210,   286,
     209,    34,   209,   210,   209,   210,   209,   210,   286,   209,
     210,   286,   286,    87,   229,    87,   230,    87,   231,   286,
     209,   209,   210,     0,   221,   209,    10,    11,    12,    13,
      14,    15,    16,    22,    26,    42,    89,    90,   107,   108,
     141,   142,   340,   341,   342,   375,   376,   377,   378,   379,
     410,   411,   413,   414,   417,   418,   419,   420,   421,   422,
     423,   209,    17,    21,    44,   341,   344,   345,   383,   400,
     424,    24,     4,    87,   322,   323,   123,   278,   279,   355,
      43,   209,    54,   209,   209,   217,    87,   209,   217,    87,
      87,   242,   243,    34,     5,    40,    41,    49,    56,    57,
      58,    61,    62,    64,    65,    70,    71,    72,    73,    74,
      75,    76,    78,    79,    80,    86,    92,    94,    99,   110,
     122,   129,   304,   305,   355,   365,   375,   376,   377,   378,
     379,   380,   381,   382,   383,   384,   385,   386,   387,   388,
     389,   390,   391,   392,   393,   394,   395,   396,   397,   398,
     399,   400,   401,   402,   403,   405,   406,   410,   411,   412,
     413,   414,   415,    87,   154,   209,    23,    87,   126,   290,
     291,   292,    23,    87,   126,   299,   300,    23,    87,   126,
     296,   297,    87,   245,   246,   242,    39,   240,    43,   209,
     253,    63,   124,   125,   135,   357,    84,    96,   109,   329,
     330,   407,   408,   409,    23,   137,   267,   268,    43,    54,
      68,    77,    87,   143,   155,   156,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   176,   210,   237,    87,   314,
     315,    87,   318,     3,    25,    27,    28,    51,    52,    91,
      93,    95,   101,   105,   113,   114,   115,   119,   120,   121,
     285,   354,   355,   356,   357,   358,   359,   360,   361,   362,
     363,   364,   365,   366,   367,   368,   369,   372,   373,   374,
     382,   404,   408,   409,   209,   209,   132,    87,   154,   209,
      54,   209,    43,    54,    68,    77,    87,   143,   155,   156,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   176,
     210,   262,   264,   307,   308,   365,   382,   383,   392,   398,
     399,   400,   401,   402,   403,   410,   411,   412,   307,   209,
     267,   214,    69,   281,   282,   368,   370,   375,   145,   146,
     275,   276,   355,   452,   453,   284,   285,   209,   131,   285,
     338,   339,   416,   209,   209,   132,    87,   154,   209,   132,
      87,   154,   209,   132,    87,   154,   209,   209,    87,   351,
     352,   178,   179,   180,   181,   182,   183,   184,   185,   186,
     187,   188,   189,   190,   191,   192,   193,   194,   195,   196,
     197,   198,   202,   203,   205,   206,   207,   208,   348,   349,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
     154,   211,    34,    34,    34,   154,   211,   211,    87,   154,
     210,   324,    32,   323,    34,   154,   211,   209,   209,    87,
     211,   217,    87,   211,   217,    34,    32,   243,    87,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,   154,   211,    34,    87,    87,
      87,    32,   291,   154,    87,   154,    87,    32,   300,    87,
     154,    87,    32,   297,   210,    32,   246,    32,    34,   211,
     209,   212,   260,   261,   262,   263,   154,   211,   154,   211,
     211,   211,    34,   154,   211,    87,    87,    32,   268,   210,
     210,   210,   210,   237,   237,   210,   210,   210,   210,   210,
     210,   210,   210,   210,   210,   237,   147,   149,   150,   151,
     152,   153,   155,   156,   157,   158,   160,   209,   210,    32,
     315,   154,   237,    32,    87,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,   211,   132,    87,
     209,   210,   210,   210,   210,   262,   262,   210,   210,   210,
     210,   210,   210,   210,   210,   210,   210,   262,   147,   149,
     150,   151,   152,   153,   155,   156,   157,   158,   160,   336,
     154,   211,   211,    32,    43,    54,   210,   272,   273,    34,
     154,   211,    34,    34,   154,   211,   154,   211,    34,   154,
     211,   132,    87,   132,    87,   132,    87,   154,    32,   352,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,   211,   154,    43,
      54,   353,    43,   155,   156,   288,   353,    54,   288,    54,
      87,    54,    54,   214,   456,   457,    54,    54,    87,    87,
     454,   342,    54,    54,   353,    54,   345,    54,   209,   210,
      87,    43,    54,    34,    54,   279,   209,   209,   209,   286,
      87,   209,   209,   286,    87,   237,   457,    54,    54,    54,
      54,    54,   353,   353,   353,    54,    54,    54,    54,    87,
     210,   353,   305,   209,   286,    87,    34,   154,     6,    43,
      45,    54,    55,    87,    97,   133,   155,   293,   301,   302,
     154,   302,   154,   154,   302,   154,    54,   155,   156,   287,
      87,   209,    87,    32,   261,   263,    34,   209,    67,   248,
     251,   369,   371,   259,    47,    60,    88,   247,   249,   250,
     251,   257,   209,   209,    59,    85,   330,    87,   157,   213,
     217,   218,   331,   332,   333,   154,    34,   154,   209,   237,
     237,   237,   238,   237,   237,   237,   237,   237,   237,   237,
     237,   237,   237,   237,   211,   237,   237,   237,   237,   237,
     237,   237,   237,   237,   237,   237,   237,    87,   209,   154,
     237,    54,   353,    54,    54,    54,    54,    54,    54,   353,
      54,    54,    54,   209,   286,   132,   262,   262,   287,   262,
     262,   262,   262,   262,   262,   262,   262,   262,   262,   262,
     211,   262,   262,   262,   262,   262,   262,   262,   262,   262,
     262,   262,   209,   308,   209,   286,   209,   286,   237,   209,
     215,    43,    54,   154,   210,     7,    46,    48,    54,    66,
     282,   209,    54,    54,   276,   209,   285,   209,   286,   353,
     339,   209,   286,   132,   132,   132,   237,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    54,
      54,    54,   353,   353,    54,   457,   353,   353,    54,    54,
     353,    54,   353,   353,   209,   348,    43,    43,    54,   455,
     215,   455,   213,   209,   209,    54,   324,   211,   211,   237,
     209,   211,   209,   211,   209,   216,   310,   311,   209,    87,
      87,    43,    54,   209,   154,   154,    87,   154,   302,    87,
     209,   302,    54,    54,   211,   242,    34,   262,    34,   154,
     211,   209,   154,   211,   209,   255,   254,   154,   209,   210,
     333,    87,   237,    87,   105,   126,   154,   154,   154,   211,
     154,   211,   211,   211,   211,   211,   211,   211,   211,   211,
     211,   211,   237,    87,   209,   209,   154,   154,   211,   154,
     211,   211,   211,   211,   211,   211,   211,   211,   211,   211,
     209,   209,   211,   273,   209,    43,    54,   210,   237,   209,
     209,   154,   213,    87,   211,    34,   209,   209,   286,   209,
     286,    87,   154,   211,   294,   302,   301,   302,   154,   302,
     154,   154,   209,    34,    32,   262,   209,   353,   251,   258,
     260,   249,   256,   260,   260,   260,   332,   302,    34,   209,
      34,    54,   269,   237,   237,   237,   237,   237,   209,   209,
     237,   262,   262,   262,   237,   211,   237,    54,   324,   237,
     209,   209,   216,   310,   154,   154,   154,   302,   209,   302,
     302,   237,   209,   209,    32,   209,    32,    32,    32,   210,
     211,   237,   237,   213,    54,   154,   209,   209,   211,   211,
     154,   209,   211,   211,   154,   211,   209,    34,   209,   154,
     302,   294,   302,   154,   209,   209,   209,   260,   260,   302,
     209,   209,    54,   213,    54,   136,   237,   262,   237,   216,
     302,   154,   154,   302,    32,    32,   211,    54,   213,   237,
     271,   211,   211,   209,    87,   302,   301,   209,    54,   154,
     209,   353,   216,   154,   154,   237,   302,   294,   154,   302
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
     455,   456,   457,   458,   459,   460,   461,   462,   463,    59,
      40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   219,   220,   220,   221,   221,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   221,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   221,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   221,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   222,   222,   222,   222,   223,
     223,   224,   225,   226,   227,   228,   229,   229,   229,   229,
     229,   229,   230,   230,   230,   230,   230,   230,   231,   231,
     231,   231,   231,   231,   232,   232,   232,   232,   232,   232,
     233,   233,   234,   234,   235,   235,   236,   237,   237,   237,
     237,   237,   237,   237,   237,   237,   237,   237,   237,   237,
     237,   237,   237,   237,   237,   237,   237,   237,   237,   237,
     237,   237,   237,   237,   237,   237,   237,   237,   238,   238,
     239,   239,   240,   241,   242,   242,   243,   244,   245,   245,
     246,   247,   247,   248,   248,   249,   249,   250,   250,   250,
     251,   251,   253,   252,   254,   252,   255,   252,   256,   252,
     257,   252,   258,   252,   259,   252,   260,   260,   260,   260,
     261,   261,   262,   262,   262,   262,   262,   262,   262,   262,
     262,   262,   262,   262,   262,   262,   262,   262,   262,   262,
     262,   262,   262,   262,   262,   262,   262,   262,   262,   262,
     262,   262,   263,   264,   264,   265,   266,   267,   267,   268,
     268,   268,   268,   268,   269,   269,   269,   269,   269,   269,
     270,   271,   271,   271,   272,   272,   273,   273,   273,   273,
     273,   273,   273,   273,   273,   274,   274,   275,   275,   276,
     276,   276,   277,   277,   278,   278,   279,   280,   280,   281,
     281,   282,   282,   282,   283,   283,   283,   283,   284,   284,
     285,   285,   285,   285,   285,   285,   285,   285,   285,   285,
     285,   285,   285,   285,   285,   285,   285,   285,   285,   285,
     285,   285,   285,   286,   286,   286,   286,   286,   286,   287,
     287,   287,   288,   288,   288,   289,   290,   290,   291,   292,
     292,   292,   293,   293,   293,   293,   293,   294,   294,   294,
     294,   295,   296,   296,   297,   297,   297,   298,   299,   299,
     300,   300,   300,   301,   301,   301,   301,   301,   302,   302,
     302,   302,   302,   302,   303,   303,   303,   303,   304,   304,
     305,   305,   305,   305,   305,   305,   305,   305,   305,   305,
     305,   305,   305,   305,   305,   305,   305,   305,   305,   305,
     305,   305,   305,   305,   305,   305,   305,   305,   305,   305,
     305,   305,   305,   305,   305,   305,   305,   305,   305,   306,
     306,   307,   307,   308,   308,   308,   308,   308,   308,   308,
     308,   308,   308,   308,   308,   308,   309,   309,   310,   310,
     311,   311,   312,   313,   314,   314,   315,   316,   317,   318,
     318,   318,   318,   319,   320,   320,   320,   320,   321,   322,
     322,   323,   323,   323,   324,   324,   324,   325,   325,   326,
     326,   326,   326,   326,   326,   327,   327,   327,   327,   327,
     327,   328,   329,   329,   330,   330,   330,   331,   331,   331,
     331,   332,   332,   333,   333,   333,   333,   333,   335,   336,
     334,   337,   337,   337,   337,   338,   338,   339,   339,   340,
     340,   340,   340,   340,   340,   340,   341,   341,   341,   341,
     341,   341,   341,   341,   341,   341,   342,   342,   343,   343,
     344,   344,   344,   344,   345,   345,   346,   346,   347,   347,
     348,   348,   349,   349,   349,   349,   349,   349,   349,   349,
     349,   349,   349,   349,   349,   349,   349,   349,   349,   349,
     349,   349,   349,   349,   349,   349,   349,   349,   349,   350,
     351,   351,   352,   353,   353,   354,   355,   356,   357,   358,
     359,   360,   361,   362,   363,   364,   365,   366,   367,   368,
     369,   370,   370,   370,   370,   370,   371,   372,   373,   374,
     375,   376,   376,   377,   378,   379,   380,   381,   382,   382,
     383,   384,   385,   386,   387,   388,   389,   390,   391,   392,
     393,   394,   395,   396,   397,   398,   399,   400,   401,   402,
     403,   404,   405,   406,   407,   407,   408,   409,   410,   411,
     412,   413,   414,   415,   416,   417,   418,   419,   420,   421,
     422,   423,   424,   425,   426,   427,   428,   429,   430,   431,
     432,   433,   434,   435,   436,   437,   438,   439,   440,   441,
     442,   443,   444,   445,   446,   447,   448,   449,   450,   451,
     452,   453,   454,   455,   455,   456,   456,   457
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
       7,     3,     1,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     0,     5,     0,     8,     0,     8,     0,    10,
       0,     8,     0,    10,     0,     8,     2,     2,     1,     1,
       4,     2,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     2,     2,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     6,
       6,     8,     5,     1,     4,     4,     4,     2,     1,     9,
       6,     5,     7,     7,     3,     2,     5,     4,     3,     1,
       6,     3,     2,     1,     3,     1,     5,     3,     3,     4,
       2,     2,     3,     1,     1,     2,     5,     3,     1,     1,
       1,     1,     2,     5,     3,     1,     1,     2,     5,     3,
       1,     1,     1,     1,     2,     5,     3,     6,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     2,     4,     3,     5,     1,     3,     2,
       2,     1,     2,     2,     1,     4,     2,     1,     4,     2,
       1,     4,     3,     5,     9,     1,     5,     3,     5,     7,
       9,     4,     2,     1,     5,     7,     4,     4,     2,     1,
       7,     9,     6,     1,     1,     1,     1,     1,     0,     1,
       1,     1,     2,     2,     2,     5,     3,     6,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     5,
       6,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     5,     6,     7,     5,
       1,     3,     3,     4,     2,     1,     5,     3,     4,     4,
       6,     3,     5,     3,     2,     5,     3,     6,     4,     2,
       1,     5,     7,     9,     0,     3,     3,     2,     5,     5,
       6,     3,     7,     8,     5,     5,     6,     3,     7,     8,
       5,     6,     3,     1,     1,     1,     1,     1,     3,     4,
       6,     1,     2,     1,     1,     1,     1,     1,     0,     0,
       5,     2,     5,     3,     6,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     1,     3,     6,
       1,     1,     1,     1,     3,     1,     3,     6,     2,     5,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     4,
       1,     2,     6,     1,     1,     3,     3,     3,     1,     3,
       3,     3,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     1,
       3,     3,     3,     3,     5,     3,     3,     3,     1,     3,
       3,     3,     1,     1,     1,     1,     1,     3,     1,     1,
       1,     1,     3,     3,     3,     3,     1,     1,     3,     3,
       3,     1,     1,     1,     3,     3,     3,     3,     3,     3,
       1,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     1,     3,     2,     2,     2
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
       220,     0,    -1,   221,    -1,   220,   221,    -1,   222,    -1,
     233,    -1,   234,    -1,   235,    -1,   252,    -1,   239,    -1,
     241,    -1,   244,    -1,   236,    -1,   265,    -1,   266,    -1,
     270,    -1,   274,    -1,   277,    -1,   280,    -1,   283,    -1,
     303,    -1,   306,    -1,   309,    -1,   289,    -1,   298,    -1,
     295,    -1,   312,    -1,   313,    -1,   316,    -1,   223,    -1,
     224,    -1,   317,    -1,   319,    -1,   320,    -1,   321,    -1,
     325,    -1,   326,    -1,   327,    -1,   328,    -1,   334,    -1,
     337,    -1,   343,    -1,   346,    -1,   347,    -1,   350,    -1,
     228,    -1,   225,    -1,   226,    -1,   227,    -1,    29,    54,
     209,    -1,    29,    54,    54,   209,    -1,   116,   286,   209,
      -1,   137,   229,   209,    -1,   138,   230,   209,    -1,   139,
     231,   209,    -1,   104,   232,   209,    -1,   229,    87,    -1,
     229,   154,    87,    -1,    87,    -1,   229,    87,   132,    -1,
     229,   154,    87,   132,    -1,    87,   132,    -1,   230,    87,
      -1,   230,   154,    87,    -1,    87,    -1,   230,    87,   132,
      -1,   230,   154,    87,   132,    -1,    87,   132,    -1,   231,
      87,    -1,   231,   154,    87,    -1,    87,    -1,   231,    87,
     132,    -1,   231,   154,    87,   132,    -1,    87,   132,    -1,
     232,    87,    -1,   232,   154,    87,    -1,    87,    -1,   232,
      87,   132,    -1,   232,   154,    87,   132,    -1,    87,   132,
      -1,   105,    54,   209,    -1,   105,    34,    54,   209,    -1,
      25,    43,   209,    -1,    25,    34,    43,   209,    -1,    67,
      43,   209,    -1,    67,    34,    43,   209,    -1,    87,    34,
     237,   209,    -1,   210,   237,   211,    -1,    87,    -1,    43,
      -1,    54,    -1,   237,   156,   237,    -1,   237,   155,   237,
      -1,   237,   157,   237,    -1,   237,   158,   237,    -1,   237,
     160,   237,    -1,   237,   153,   237,    -1,   237,   152,   237,
      -1,   237,   151,   237,    -1,   237,   150,   237,    -1,   237,
     149,   237,    -1,   237,   147,   237,    -1,   155,   237,    -1,
     156,   237,    -1,   161,   210,   237,   211,    -1,   162,   210,
     237,   211,    -1,   163,   210,   237,   211,    -1,   164,   210,
     237,   211,    -1,   165,   210,   237,   211,    -1,   166,   210,
     237,   211,    -1,   167,   210,   237,   211,    -1,   168,   210,
     237,   211,    -1,   169,   210,   237,   211,    -1,   176,   210,
     237,   211,    -1,    68,   210,   237,   154,   237,   211,    -1,
      77,   210,   237,   154,   237,   211,    -1,    87,   210,   238,
     211,    -1,   143,   210,   237,   154,   237,   154,   237,   211,
      -1,   237,    -1,   238,   154,   237,    -1,    53,   209,   242,
      32,    -1,    53,   210,   240,   211,   209,   242,    32,    -1,
      39,    34,    87,    -1,    33,   209,   242,    32,    -1,   242,
     243,    -1,   243,    -1,    87,    34,   237,   209,    -1,    50,
     209,   245,    32,    -1,   245,   246,    -1,   246,    -1,    87,
     210,   287,   211,    34,   237,   209,    -1,   247,   154,   249,
      -1,   249,    -1,   248,   154,   251,    -1,   251,    -1,   250,
      -1,   251,    -1,    60,    -1,    47,    -1,    88,    -1,   369,
      -1,   371,    -1,    -1,    81,   209,   253,   260,    32,    -1,
      -1,    81,   210,   357,   211,   209,   254,   260,    32,    -1,
      -1,    81,   210,   135,   211,   209,   255,   260,    32,    -1,
      -1,    81,   210,   125,   154,   247,   211,   256,   209,   260,
      32,    -1,    -1,    81,   210,   125,   211,   257,   209,   260,
      32,    -1,    -1,    81,   210,   124,   154,   248,   211,   258,
     209,   260,    32,    -1,    -1,    81,   210,   124,   211,   259,
     209,   260,    32,    -1,   260,   261,    -1,   260,   263,    -1,
     261,    -1,   263,    -1,   262,    34,   262,   209,    -1,   262,
     209,    -1,   210,   262,   211,    -1,   264,    -1,    43,    -1,
      54,    -1,   262,   156,   262,    -1,   262,   155,   262,    -1,
     262,   157,   262,    -1,   262,   158,   262,    -1,   262,   153,
     262,    -1,   262,   152,   262,    -1,   262,   151,   262,    -1,
     262,   150,   262,    -1,   262,   149,   262,    -1,   262,   147,
     262,    -1,   262,   160,   262,    -1,   155,   262,    -1,   156,
     262,    -1,   161,   210,   262,   211,    -1,   162,   210,   262,
     211,    -1,   163,   210,   262,   211,    -1,   164,   210,   262,
     211,    -1,   165,   210,   262,   211,    -1,   166,   210,   262,
     211,    -1,   167,   210,   262,   211,    -1,   168,   210,   262,
     211,    -1,   169,   210,   262,   211,    -1,   176,   210,   262,
     211,    -1,    68,   210,   262,   154,   262,   211,    -1,    77,
     210,   262,   154,   262,   211,    -1,   143,   210,   262,   154,
     262,   154,   262,   211,    -1,   212,    87,    34,   262,   209,
      -1,    87,    -1,    87,   210,   287,   211,    -1,   117,   209,
     267,    32,    -1,    83,   209,   267,    32,    -1,   267,   268,
      -1,   268,    -1,   137,    87,   209,   105,   269,   209,   136,
     271,   209,    -1,   137,    87,   209,   126,   237,   209,    -1,
     137,    87,    34,   237,   209,    -1,   137,    87,   154,    87,
      34,   237,   209,    -1,    23,    87,   154,    87,    34,   237,
     209,    -1,   269,   154,    54,    -1,   269,    54,    -1,   269,
     154,    54,   213,    54,    -1,   269,    54,   213,    54,    -1,
      54,   213,    54,    -1,    54,    -1,   118,    34,   214,   272,
     215,   209,    -1,   271,   154,   237,    -1,   271,   353,    -1,
     237,    -1,   272,   209,   273,    -1,   273,    -1,   273,   154,
     210,   237,   211,    -1,   273,   154,    43,    -1,   273,   154,
      54,    -1,   273,   210,   237,   211,    -1,   273,    43,    -1,
     273,    54,    -1,   210,   237,   211,    -1,    43,    -1,    54,
      -1,   127,   209,    -1,   127,   210,   275,   211,   209,    -1,
     275,   154,   276,    -1,   276,    -1,   355,    -1,   452,    -1,
     453,    -1,    20,   209,    -1,    20,   210,   278,   211,   209,
      -1,   278,   154,   279,    -1,   279,    -1,   355,    -1,   119,
     209,    -1,   119,   210,   281,   211,   209,    -1,   281,   154,
     282,    -1,   282,    -1,   368,    -1,   375,    -1,   370,    -1,
     128,   209,    -1,   128,   210,   284,   211,   209,    -1,   128,
     286,   209,    -1,   128,   210,   284,   211,   286,   209,    -1,
     284,   154,   285,    -1,   285,    -1,   354,    -1,   355,    -1,
     356,    -1,   357,    -1,   358,    -1,   359,    -1,   360,    -1,
     361,    -1,   362,    -1,   363,    -1,   364,    -1,   382,    -1,
     365,    -1,   404,    -1,   366,    -1,   367,    -1,   368,    -1,
     369,    -1,   372,    -1,   373,    -1,   374,    -1,   408,    -1,
     409,    -1,   286,    87,    -1,   286,    87,    34,    87,    -1,
     286,   154,    87,    -1,   286,   154,    87,    34,    87,    -1,
      87,    -1,    87,    34,    87,    -1,   156,    54,    -1,   155,
      54,    -1,    54,    -1,   156,    43,    -1,   155,    43,    -1,
      43,    -1,    36,   209,   290,    32,    -1,   290,   291,    -1,
     291,    -1,   292,   154,   293,   209,    -1,   126,    87,    -1,
      87,    -1,    23,    87,   154,    87,    -1,   301,   154,   294,
      -1,   302,   154,   301,   154,   294,    -1,   302,   154,   302,
     154,   302,   154,   301,   154,   294,    -1,   302,    -1,   302,
     154,   302,   154,   302,    -1,   302,   154,   302,    -1,   302,
     154,   302,   154,   302,    -1,   302,   154,   302,   154,   302,
     154,   302,    -1,   302,   154,   302,   154,   302,   154,   302,
     154,   302,    -1,    38,   209,   296,    32,    -1,   296,   297,
      -1,   297,    -1,   126,    87,   154,   302,   209,    -1,    23,
      87,   154,    87,   154,   302,   209,    -1,    87,   154,   302,
     209,    -1,    37,   209,   299,    32,    -1,   299,   300,    -1,
     300,    -1,   126,    87,   154,   302,   154,   302,   209,    -1,
      23,    87,   154,    87,   154,   302,   154,   302,   209,    -1,
      87,   154,   302,   154,   302,   209,    -1,     6,    -1,    45,
      -1,    97,    -1,    55,    -1,   133,    -1,    -1,    54,    -1,
      43,    -1,    87,    -1,   155,    54,    -1,   155,    43,    -1,
      35,   209,    -1,    35,   210,   304,   211,   209,    -1,    35,
     286,   209,    -1,    35,   210,   304,   211,   286,   209,    -1,
     304,   154,   305,    -1,   305,    -1,   375,    -1,   376,    -1,
     377,    -1,   378,    -1,   379,    -1,   380,    -1,   381,    -1,
     382,    -1,   383,    -1,   384,    -1,   385,    -1,   386,    -1,
     387,    -1,   388,    -1,   389,    -1,   390,    -1,   391,    -1,
     392,    -1,   393,    -1,   394,    -1,   395,    -1,   396,    -1,
     397,    -1,   398,    -1,   365,    -1,   399,    -1,   400,    -1,
     401,    -1,   402,    -1,   403,    -1,   405,    -1,   406,    -1,
     410,    -1,   411,    -1,   412,    -1,   355,    -1,   413,    -1,
     414,    -1,   415,    -1,   111,   210,   307,   211,   209,    -1,
     111,   210,   307,   211,   286,   209,    -1,   307,   154,   308,
      -1,   308,    -1,   382,    -1,   383,    -1,   392,    -1,   398,
      -1,   365,    -1,   399,    -1,   400,    -1,   401,    -1,   402,
      -1,   403,    -1,   410,    -1,   411,    -1,   412,    -1,   112,
     210,   307,   211,   209,    -1,   112,   210,   307,   211,   286,
     209,    -1,   216,    87,   216,   154,   216,    87,   216,    -1,
     216,    87,   216,   154,   302,    -1,   310,    -1,   311,   154,
     310,    -1,   140,   286,   209,    -1,    98,   209,   314,    32,
      -1,   314,   315,    -1,   315,    -1,    87,   210,   237,   211,
     209,    -1,   134,   286,   209,    -1,   100,   209,   318,    32,
      -1,   318,    87,   237,   209,    -1,   318,    87,   154,    87,
     237,   209,    -1,    87,   237,   209,    -1,    87,   154,    87,
     237,   209,    -1,   103,   286,   209,    -1,   102,   209,    -1,
     102,   210,   285,   211,   209,    -1,   102,   286,   209,    -1,
     102,   210,   285,   211,   286,   209,    -1,    19,   209,   322,
      32,    -1,   322,   323,    -1,   323,    -1,    87,   324,    34,
     237,   209,    -1,    87,   154,    87,   324,    34,   237,   209,
      -1,     4,    87,   210,    54,   211,   324,    34,   237,   209,
      -1,    -1,   210,    54,   211,    -1,   210,    43,   211,    -1,
      18,   209,    -1,    18,   210,    24,   211,   209,    -1,    31,
     210,    87,   211,   209,    -1,    31,   210,    87,   211,   286,
     209,    -1,    31,    87,   209,    -1,    31,   210,    87,   217,
      87,   211,   209,    -1,    31,   210,    87,   217,    87,   211,
     286,   209,    -1,    31,    87,   217,    87,   209,    -1,    30,
     210,    87,   211,   209,    -1,    30,   210,    87,   211,   286,
     209,    -1,    30,    87,   209,    -1,    30,   210,    87,   217,
      87,   211,   209,    -1,    30,   210,    87,   217,    87,   211,
     286,   209,    -1,    30,    87,   217,    87,   209,    -1,    82,
     210,   329,   211,   331,   209,    -1,   329,   154,   330,    -1,
     330,    -1,   407,    -1,   408,    -1,   409,    -1,   332,    -1,
     331,   154,   332,    -1,   332,   210,   302,   211,    -1,   331,
     154,   332,   210,   302,   211,    -1,   333,    -1,   332,   333,
      -1,    87,    -1,   218,    -1,   157,    -1,   213,    -1,   217,
      -1,    -1,    -1,   106,   335,   262,   336,   209,    -1,   130,
     209,    -1,   130,   210,   338,   211,   209,    -1,   130,   286,
     209,    -1,   130,   210,   338,   211,   286,   209,    -1,   338,
     154,   339,    -1,   339,    -1,   285,    -1,   416,    -1,   417,
      -1,   418,    -1,   419,    -1,   420,    -1,   421,    -1,   422,
      -1,   423,    -1,   340,    -1,   375,    -1,   410,    -1,   411,
      -1,   377,    -1,   379,    -1,   376,    -1,   378,    -1,   413,
      -1,   414,    -1,   341,   154,   342,    -1,   341,    -1,     8,
      54,   209,    -1,     8,   210,   342,   211,    54,   209,    -1,
     341,    -1,   400,    -1,   383,    -1,   424,    -1,   344,   154,
     345,    -1,   344,    -1,     9,    54,   209,    -1,     9,   210,
     345,   211,    54,   209,    -1,   177,   209,    -1,   177,   210,
     348,   211,   209,    -1,   349,   154,   348,    -1,   349,    -1,
     425,    -1,   426,    -1,   427,    -1,   428,    -1,   429,    -1,
     430,    -1,   431,    -1,   432,    -1,   433,    -1,   434,    -1,
     435,    -1,   436,    -1,   437,    -1,   438,    -1,   439,    -1,
     440,    -1,   441,    -1,   442,    -1,   444,    -1,   445,    -1,
     446,    -1,   447,    -1,   448,    -1,   449,    -1,   450,    -1,
     451,    -1,   443,    -1,   144,   209,   351,    32,    -1,   352,
      -1,   351,   352,    -1,    87,   154,   237,   154,   237,   209,
      -1,    54,    -1,    43,    -1,    27,    34,    54,    -1,   123,
      34,    54,    -1,   120,    34,    54,    -1,    63,    -1,   101,
      34,    54,    -1,   115,    34,    54,    -1,    28,    34,    54,
      -1,     3,    34,    54,    -1,    91,    -1,    93,    -1,    95,
      -1,    56,    34,    54,    -1,    51,    34,    54,    -1,    52,
      34,    54,    -1,   105,    34,    54,    -1,    25,    34,   353,
      -1,    69,    34,    54,    -1,    69,    34,    66,    -1,    69,
      34,    46,    -1,    69,    34,    48,    -1,    69,    34,     7,
      -1,    67,    34,   353,    -1,   119,    -1,   121,    34,    54,
      -1,   113,    34,   353,    -1,    26,    34,    87,    -1,    89,
      34,   457,    -1,    89,    34,    54,    -1,    42,    34,    54,
      -1,   107,    34,    54,    -1,   108,    34,    54,    -1,    61,
      34,    54,    -1,    62,    34,    54,    -1,    94,    -1,    49,
      -1,    21,    34,   353,    -1,    75,    34,    54,    -1,    70,
      34,   353,    -1,    72,    34,   353,    -1,    99,    34,   210,
     311,   211,    -1,    71,    34,   353,    -1,    80,    34,    87,
      -1,    79,    34,    54,    -1,    78,    -1,   110,    34,   353,
      -1,    73,    34,    54,    -1,    74,    34,    54,    -1,    64,
      -1,    65,    -1,    92,    -1,     5,    -1,   129,    -1,    44,
      34,    54,    -1,   122,    -1,    86,    -1,    41,    -1,   114,
      -1,    57,    34,    54,    -1,    58,    34,    54,    -1,    84,
      34,    59,    -1,    84,    34,    85,    -1,   109,    -1,    96,
      -1,   141,    34,    87,    -1,   142,    34,   454,    -1,    40,
      34,   457,    -1,    22,    -1,    90,    -1,    76,    -1,   131,
      34,   353,    -1,    15,    34,   288,    -1,    10,    34,   353,
      -1,    12,    34,   288,    -1,    13,    34,   353,    -1,    14,
      34,    54,    -1,    11,    -1,    16,    34,    54,    -1,    17,
      34,    54,    -1,   178,    34,    54,    -1,   179,    34,    54,
      -1,   180,    34,    54,    -1,   181,    34,    54,    -1,   182,
      34,    54,    -1,   183,    34,    54,    -1,   184,    34,    54,
      -1,   185,    34,    54,    -1,   186,    34,    54,    -1,   187,
      34,    54,    -1,   188,    34,    54,    -1,   189,    34,    54,
      -1,   190,    34,    54,    -1,   191,    34,    54,    -1,   192,
      34,    54,    -1,   193,    34,   353,    -1,   194,    34,   353,
      -1,   195,    34,    54,    -1,   196,    34,   457,    -1,   197,
      34,   353,    -1,   198,    34,   353,    -1,   202,    34,    54,
      -1,   203,    34,    54,    -1,   205,    34,   353,    -1,   206,
      34,    54,    -1,   207,    34,   353,    -1,   208,    34,   353,
      -1,   145,    34,    54,    -1,   146,    34,    54,    -1,    87,
     213,    87,    -1,    54,    -1,    54,   213,    54,    -1,   214,
     455,    -1,   456,   455,    -1,   456,   215,    -1
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
     429,   437,   441,   443,   447,   449,   451,   453,   455,   457,
     459,   461,   463,   464,   470,   471,   480,   481,   490,   491,
     502,   503,   512,   513,   524,   525,   534,   537,   540,   542,
     544,   549,   552,   556,   558,   560,   562,   566,   570,   574,
     578,   582,   586,   590,   594,   598,   602,   606,   609,   612,
     617,   622,   627,   632,   637,   642,   647,   652,   657,   662,
     669,   676,   685,   691,   693,   698,   703,   708,   711,   713,
     723,   730,   736,   744,   752,   756,   759,   765,   770,   774,
     776,   783,   787,   790,   792,   796,   798,   804,   808,   812,
     817,   820,   823,   827,   829,   831,   834,   840,   844,   846,
     848,   850,   852,   855,   861,   865,   867,   869,   872,   878,
     882,   884,   886,   888,   890,   893,   899,   903,   910,   914,
     916,   918,   920,   922,   924,   926,   928,   930,   932,   934,
     936,   938,   940,   942,   944,   946,   948,   950,   952,   954,
     956,   958,   960,   962,   965,   970,   974,   980,   982,   986,
     989,   992,   994,   997,  1000,  1002,  1007,  1010,  1012,  1017,
    1020,  1022,  1027,  1031,  1037,  1047,  1049,  1055,  1059,  1065,
    1073,  1083,  1088,  1091,  1093,  1099,  1107,  1112,  1117,  1120,
    1122,  1130,  1140,  1147,  1149,  1151,  1153,  1155,  1157,  1158,
    1160,  1162,  1164,  1167,  1170,  1173,  1179,  1183,  1190,  1194,
    1196,  1198,  1200,  1202,  1204,  1206,  1208,  1210,  1212,  1214,
    1216,  1218,  1220,  1222,  1224,  1226,  1228,  1230,  1232,  1234,
    1236,  1238,  1240,  1242,  1244,  1246,  1248,  1250,  1252,  1254,
    1256,  1258,  1260,  1262,  1264,  1266,  1268,  1270,  1272,  1274,
    1280,  1287,  1291,  1293,  1295,  1297,  1299,  1301,  1303,  1305,
    1307,  1309,  1311,  1313,  1315,  1317,  1319,  1325,  1332,  1340,
    1346,  1348,  1352,  1356,  1361,  1364,  1366,  1372,  1376,  1381,
    1386,  1393,  1397,  1403,  1407,  1410,  1416,  1420,  1427,  1432,
    1435,  1437,  1443,  1451,  1461,  1462,  1466,  1470,  1473,  1479,
    1485,  1492,  1496,  1504,  1513,  1519,  1525,  1532,  1536,  1544,
    1553,  1559,  1566,  1570,  1572,  1574,  1576,  1578,  1580,  1584,
    1589,  1596,  1598,  1601,  1603,  1605,  1607,  1609,  1611,  1612,
    1613,  1619,  1622,  1628,  1632,  1639,  1643,  1645,  1647,  1649,
    1651,  1653,  1655,  1657,  1659,  1661,  1663,  1665,  1667,  1669,
    1671,  1673,  1675,  1677,  1679,  1681,  1683,  1687,  1689,  1693,
    1700,  1702,  1704,  1706,  1708,  1712,  1714,  1718,  1725,  1728,
    1734,  1738,  1740,  1742,  1744,  1746,  1748,  1750,  1752,  1754,
    1756,  1758,  1760,  1762,  1764,  1766,  1768,  1770,  1772,  1774,
    1776,  1778,  1780,  1782,  1784,  1786,  1788,  1790,  1792,  1794,
    1799,  1801,  1804,  1811,  1813,  1815,  1819,  1823,  1827,  1829,
    1833,  1837,  1841,  1845,  1847,  1849,  1851,  1855,  1859,  1863,
    1867,  1871,  1875,  1879,  1883,  1887,  1891,  1895,  1897,  1901,
    1905,  1909,  1913,  1917,  1921,  1925,  1929,  1933,  1937,  1939,
    1941,  1945,  1949,  1953,  1957,  1963,  1967,  1971,  1975,  1977,
    1981,  1985,  1989,  1991,  1993,  1995,  1997,  1999,  2003,  2005,
    2007,  2009,  2011,  2015,  2019,  2023,  2027,  2029,  2031,  2035,
    2039,  2043,  2045,  2047,  2049,  2053,  2057,  2061,  2065,  2069,
    2073,  2075,  2079,  2083,  2087,  2091,  2095,  2099,  2103,  2107,
    2111,  2115,  2119,  2123,  2127,  2131,  2135,  2139,  2143,  2147,
    2151,  2155,  2159,  2163,  2167,  2171,  2175,  2179,  2183,  2187,
    2191,  2195,  2199,  2203,  2205,  2209,  2212,  2215
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
     338,   340,   341,   344,   345,   348,   349,   352,   354,   356,
     360,   361,   364,   364,   366,   366,   368,   368,   371,   370,
     373,   373,   376,   375,   378,   378,   382,   383,   384,   385,
     388,   390,   394,   396,   397,   399,   401,   403,   405,   407,
     409,   411,   413,   415,   417,   419,   421,   423,   425,   427,
     429,   431,   433,   435,   437,   439,   441,   443,   445,   447,
     449,   451,   455,   458,   460,   464,   466,   468,   469,   472,
     474,   476,   478,   480,   484,   486,   488,   490,   492,   494,
     499,   502,   504,   506,   510,   512,   516,   518,   520,   522,
     524,   526,   528,   530,   532,   536,   538,   542,   543,   546,
     547,   548,   551,   553,   557,   558,   561,   563,   565,   569,
     570,   573,   574,   575,   578,   580,   582,   584,   588,   589,
     592,   593,   594,   595,   596,   597,   598,   599,   600,   601,
     602,   603,   604,   605,   606,   607,   608,   609,   610,   611,
     612,   613,   614,   617,   619,   621,   623,   625,   627,   631,
     633,   635,   639,   641,   643,   647,   649,   651,   655,   657,
     663,   669,   679,   684,   691,   702,   707,   718,   725,   734,
     745,   760,   763,   765,   769,   777,   787,   797,   800,   802,
     806,   816,   828,   840,   842,   844,   846,   848,   852,   853,
     854,   855,   856,   858,   862,   864,   866,   868,   872,   873,
     876,   877,   878,   879,   880,   881,   882,   883,   884,   885,
     886,   887,   888,   889,   890,   891,   892,   893,   894,   895,
     896,   897,   898,   899,   900,   901,   902,   903,   904,   905,
     906,   907,   908,   909,   910,   911,   912,   913,   914,   917,
     919,   923,   924,   927,   928,   929,   930,   931,   932,   933,
     934,   935,   936,   937,   938,   939,   942,   944,   948,   950,
     954,   955,   958,   960,   962,   963,   966,   968,   970,   972,
     974,   976,   978,   982,   984,   986,   988,   990,   994,   996,
     997,  1000,  1002,  1004,  1008,  1009,  1011,  1015,  1017,  1021,
    1023,  1025,  1027,  1029,  1031,  1035,  1037,  1039,  1041,  1043,
    1045,  1049,  1052,  1053,  1056,  1057,  1058,  1061,  1063,  1065,
    1067,  1071,  1073,  1077,  1078,  1080,  1082,  1084,  1088,  1089,
    1088,  1091,  1093,  1095,  1097,  1101,  1102,  1105,  1106,  1109,
    1110,  1111,  1112,  1113,  1114,  1115,  1118,  1119,  1120,  1121,
    1122,  1123,  1124,  1125,  1126,  1127,  1130,  1131,  1134,  1136,
    1140,  1141,  1142,  1143,  1146,  1147,  1150,  1152,  1156,  1158,
    1162,  1163,  1166,  1167,  1168,  1169,  1170,  1171,  1172,  1173,
    1174,  1175,  1176,  1177,  1178,  1179,  1180,  1181,  1182,  1183,
    1184,  1185,  1186,  1187,  1188,  1189,  1190,  1191,  1192,  1195,
    1198,  1199,  1202,  1205,  1206,  1209,  1210,  1211,  1212,  1213,
    1214,  1215,  1216,  1217,  1218,  1219,  1220,  1221,  1222,  1223,
    1225,  1226,  1227,  1228,  1229,  1230,  1231,  1232,  1233,  1234,
    1235,  1236,  1238,  1241,  1242,  1243,  1244,  1245,  1246,  1248,
    1251,  1252,  1253,  1254,  1255,  1256,  1257,  1258,  1259,  1260,
    1261,  1262,  1263,  1264,  1265,  1266,  1267,  1268,  1269,  1270,
    1271,  1272,  1273,  1274,  1275,  1277,  1280,  1281,  1282,  1283,
    1284,  1285,  1286,  1287,  1288,  1290,  1291,  1292,  1293,  1294,
    1295,  1296,  1297,  1299,  1300,  1301,  1302,  1303,  1304,  1305,
    1306,  1307,  1308,  1309,  1310,  1311,  1312,  1313,  1314,  1315,
    1316,  1317,  1319,  1320,  1326,  1327,  1331,  1332,  1333,  1334,
    1336,  1337,  1339,  1347,  1348,  1357,  1359,  1368
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
       2,     2,     2,     2,     2,   212,     2,     2,     2,   216,
     210,   211,     2,     2,     2,     2,   217,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   213,   209,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   214,   218,   215,     2,     2,     2,     2,     2,     2,
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
     205,   206,   207,   208
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 2504;
  const int parser::yynnts_ = 239;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 163;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 219;

  const unsigned int parser::yyuser_token_number_max_ = 463;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1370 "DynareBison.yy"


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

