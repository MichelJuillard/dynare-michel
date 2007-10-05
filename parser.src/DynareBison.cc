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
#line 150 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 49:
#line 152 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 50:
#line 155 "DynareBison.yy"
    { driver.rplot(); ;}
    break;

  case 55:
#line 166 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 56:
#line 168 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 57:
#line 170 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 58:
#line 172 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 59:
#line 174 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 60:
#line 176 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 61:
#line 180 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 62:
#line 182 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 63:
#line 184 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 64:
#line 186 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 65:
#line 188 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 66:
#line 190 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 67:
#line 194 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 68:
#line 196 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 69:
#line 198 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 70:
#line 200 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 71:
#line 202 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 72:
#line 204 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 73:
#line 208 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 74:
#line 210 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 75:
#line 212 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 76:
#line 214 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 77:
#line 216 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 78:
#line 218 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 79:
#line 222 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 80:
#line 224 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 81:
#line 228 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 82:
#line 230 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 83:
#line 234 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 84:
#line 236 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 85:
#line 240 "DynareBison.yy"
    { driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 86:
#line 243 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 87:
#line 245 "DynareBison.yy"
    { (yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 88:
#line 247 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 89:
#line 249 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 90:
#line 251 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 91:
#line 253 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 92:
#line 255 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 93:
#line 257 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 94:
#line 259 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 95:
#line 261 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 96:
#line 263 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 97:
#line 265 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 98:
#line 267 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 99:
#line 269 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 100:
#line 271 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 101:
#line 273 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 102:
#line 275 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 103:
#line 277 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 104:
#line 279 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 105:
#line 281 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 106:
#line 283 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 107:
#line 285 "DynareBison.yy"
    { (yyval.node_val) = driver.add_dummy((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 108:
#line 287 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 109:
#line 289 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 110:
#line 291 "DynareBison.yy"
    { (yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val)); ;}
    break;

  case 111:
#line 295 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 112:
#line 297 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 113:
#line 301 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 114:
#line 303 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 115:
#line 306 "DynareBison.yy"
    { driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 116:
#line 308 "DynareBison.yy"
    { driver.end_endval(); ;}
    break;

  case 119:
#line 314 "DynareBison.yy"
    { driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 120:
#line 316 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 123:
#line 322 "DynareBison.yy"
    { driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 126:
#line 329 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 127:
#line 331 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 128:
#line 333 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 131:
#line 338 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 132:
#line 339 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 133:
#line 340 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 134:
#line 341 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 135:
#line 342 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 136:
#line 343 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 137:
#line 345 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 138:
#line 346 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 139:
#line 347 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 140:
#line 348 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 145:
#line 358 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 146:
#line 360 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val)); ;}
    break;

  case 147:
#line 364 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 149:
#line 367 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 150:
#line 369 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 151:
#line 371 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 152:
#line 373 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 153:
#line 375 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 154:
#line 377 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 155:
#line 379 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 156:
#line 381 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 157:
#line 383 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 158:
#line 385 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 159:
#line 387 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 160:
#line 389 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 161:
#line 391 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 162:
#line 393 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 163:
#line 395 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 164:
#line 397 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 165:
#line 399 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 166:
#line 401 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 167:
#line 403 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 168:
#line 405 "DynareBison.yy"
    { (yyval.node_val) = driver.add_dummy((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 169:
#line 407 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 170:
#line 409 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 171:
#line 413 "DynareBison.yy"
    { driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 172:
#line 416 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 173:
#line 418 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 174:
#line 421 "DynareBison.yy"
    { driver.end_shocks(); ;}
    break;

  case 175:
#line 423 "DynareBison.yy"
    { driver.end_mshocks(); ;}
    break;

  case 178:
#line 430 "DynareBison.yy"
    { driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val)); ;}
    break;

  case 179:
#line 432 "DynareBison.yy"
    { driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 180:
#line 434 "DynareBison.yy"
    { driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 181:
#line 436 "DynareBison.yy"
    { driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 182:
#line 438 "DynareBison.yy"
    { driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 183:
#line 442 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 184:
#line 444 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 185:
#line 446 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 186:
#line 448 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 187:
#line 450 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 188:
#line 452 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 189:
#line 456 "DynareBison.yy"
    { driver.add_value((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 190:
#line 458 "DynareBison.yy"
    { driver.add_value((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 191:
#line 461 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 192:
#line 464 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 193:
#line 466 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 194:
#line 470 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 195:
#line 472 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 196:
#line 474 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 197:
#line 476 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 198:
#line 478 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 199:
#line 480 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 200:
#line 482 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 201:
#line 484 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 202:
#line 486 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 203:
#line 490 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 204:
#line 492 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 208:
#line 502 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 209:
#line 504 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 213:
#line 514 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 214:
#line 516 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 219:
#line 528 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 220:
#line 530 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 221:
#line 532 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 222:
#line 534 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 248:
#line 567 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 249:
#line 569 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 250:
#line 571 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 251:
#line 573 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 252:
#line 575 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 253:
#line 577 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 254:
#line 581 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 255:
#line 583 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 256:
#line 585 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 257:
#line 589 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 258:
#line 591 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 259:
#line 593 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 260:
#line 596 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 261:
#line 599 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 262:
#line 601 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 264:
#line 607 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 265:
#line 613 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 266:
#line 619 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 267:
#line 629 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 268:
#line 634 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 269:
#line 641 "DynareBison.yy"
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

  case 270:
#line 652 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 271:
#line 657 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 272:
#line 668 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 273:
#line 675 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 274:
#line 684 "DynareBison.yy"
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

  case 275:
#line 695 "DynareBison.yy"
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

  case 276:
#line 710 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 277:
#line 713 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 278:
#line 715 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 279:
#line 719 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 280:
#line 727 "DynareBison.yy"
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

  case 281:
#line 737 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 282:
#line 747 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 283:
#line 750 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 284:
#line 752 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 285:
#line 756 "DynareBison.yy"
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

  case 286:
#line 766 "DynareBison.yy"
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

  case 287:
#line 778 "DynareBison.yy"
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

  case 288:
#line 790 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 289:
#line 792 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 290:
#line 794 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 291:
#line 796 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 292:
#line 798 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 293:
#line 801 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 297:
#line 806 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 298:
#line 808 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 299:
#line 812 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 300:
#line 814 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 301:
#line 816 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 302:
#line 818 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 344:
#line 867 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 345:
#line 869 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 361:
#line 892 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 362:
#line 894 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 363:
#line 898 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 364:
#line 900 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 367:
#line 907 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 368:
#line 909 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 371:
#line 915 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 372:
#line 917 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 373:
#line 919 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 374:
#line 922 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 375:
#line 924 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 376:
#line 926 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 377:
#line 928 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 378:
#line 931 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 379:
#line 934 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 380:
#line 936 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 381:
#line 938 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 382:
#line 940 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 383:
#line 943 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 386:
#line 950 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 387:
#line 952 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 388:
#line 954 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 389:
#line 957 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 390:
#line 959 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 391:
#line 961 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 392:
#line 965 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 393:
#line 967 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 394:
#line 971 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 395:
#line 973 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 396:
#line 975 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 397:
#line 977 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 398:
#line 979 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 399:
#line 981 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 400:
#line 985 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 401:
#line 987 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 402:
#line 989 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 403:
#line 991 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 404:
#line 993 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 405:
#line 995 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 406:
#line 999 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 412:
#line 1011 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 413:
#line 1013 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 414:
#line 1015 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 415:
#line 1017 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 416:
#line 1021 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 417:
#line 1023 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 419:
#line 1028 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 420:
#line 1030 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 421:
#line 1032 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 422:
#line 1034 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 423:
#line 1037 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 424:
#line 1038 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 426:
#line 1041 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 427:
#line 1043 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 428:
#line 1045 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 429:
#line 1047 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 453:
#line 1084 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 454:
#line 1086 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 461:
#line 1100 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 462:
#line 1102 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 463:
#line 1106 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 464:
#line 1108 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 496:
#line 1148 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 497:
#line 1149 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 498:
#line 1150 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 499:
#line 1151 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 500:
#line 1152 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 501:
#line 1153 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 502:
#line 1154 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 503:
#line 1155 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 504:
#line 1156 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 505:
#line 1157 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 506:
#line 1158 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 507:
#line 1159 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 508:
#line 1160 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 509:
#line 1161 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 510:
#line 1163 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 511:
#line 1164 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 512:
#line 1165 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 513:
#line 1166 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 514:
#line 1167 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 515:
#line 1168 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 516:
#line 1169 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 517:
#line 1171 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 518:
#line 1173 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 519:
#line 1175 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 520:
#line 1176 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 521:
#line 1177 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 522:
#line 1178 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 523:
#line 1179 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 524:
#line 1181 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 525:
#line 1183 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 526:
#line 1185 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 527:
#line 1186 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 528:
#line 1187 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1188 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1190 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1191 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1192 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 534:
#line 1193 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 535:
#line 1194 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 536:
#line 1195 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 537:
#line 1196 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 538:
#line 1197 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 539:
#line 1198 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 540:
#line 1199 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 541:
#line 1200 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 542:
#line 1201 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 543:
#line 1202 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 544:
#line 1203 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 545:
#line 1204 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 546:
#line 1205 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 547:
#line 1206 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 548:
#line 1207 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 549:
#line 1208 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 550:
#line 1210 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 551:
#line 1212 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 552:
#line 1214 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 553:
#line 1215 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 554:
#line 1216 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 555:
#line 1217 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 556:
#line 1218 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 557:
#line 1219 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 558:
#line 1220 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 559:
#line 1221 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 560:
#line 1222 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 561:
#line 1224 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 562:
#line 1225 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 563:
#line 1226 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 564:
#line 1227 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1228 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 566:
#line 1229 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 567:
#line 1230 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 568:
#line 1231 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 569:
#line 1233 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 570:
#line 1234 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 571:
#line 1235 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1236 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1237 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1238 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 575:
#line 1239 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 576:
#line 1240 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 577:
#line 1241 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1242 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 579:
#line 1243 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1244 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 581:
#line 1245 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 582:
#line 1246 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1247 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1248 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1249 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1250 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1251 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1253 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1254 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1260 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1261 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 592:
#line 1265 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 593:
#line 1266 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 594:
#line 1267 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 595:
#line 1268 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 596:
#line 1272 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 598:
#line 1281 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); 
               ;}
    break;

  case 599:
#line 1290 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 600:
#line 1292 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 601:
#line 1300 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2359 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1072;
  const short int
  parser::yypact_[] =
  {
      1308,    24,    38,   -68,   -78,   111,   106,    59,    62,    65,
     -64,    -4,   -57,   -52,   -31,    91,   226,   108,   323,    92,
     127,   410,   248,   269,    11,   368,   386,    68, -1072,   285,
     294,   368,   306,   469,   351,   357,    40,    45,   368,   438,
     441,   444,   368,   407,  1189, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
     343,  1544,   345,  1524, -1072,   523,    99, -1072,   439,   506,
     375,    30,   271,   490,   289,   495,   497,   570, -1072,  1424,
      22,   355,   356,   372,   524,   497,   569,   568,   419, -1072,
     381,    42,   137,   976,   533,   538, -1072,  1595,   160,   179,
     496,   206,   576,   436,  1035,   653,   653,   208,   137,   446,
   -1072,    82, -1072,   439, -1072,  1595,   209, -1072,  1043,   210,
     211,   510,   216,   518,   218,   525,   219,   220, -1072,  1565,
   -1072, -1072, -1072,   620, -1072,   621,   622,   625,   626,   627,
   -1072,   629,   631,   632, -1072,   633,   635,   636,   637, -1072,
     534,   483, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,   650,
     652,   658, -1072,   549,   493, -1072, -1072, -1072,   501,   616,
      -3,    74, -1072,   667,   -40, -1072, -1072,   512, -1072,   515,
   -1072, -1072,   619,   -91, -1072,   628,   239,   676,    87, -1072,
     645, -1072,   678, -1072, -1072,   681,   682,   683,   695,   699,
   -1072, -1072,   703,   704,   707,   716,   721,   723, -1072, -1072,
     724,   732, -1072, -1072, -1072,   735,   738, -1072, -1072,   -37,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
     745,   649, -1072,   657, -1072,   705,   376, -1072,   648,   709,
     659,   720,   383, -1072,   725,   672,   741,   401, -1072,   593,
     244, -1072,   400,   793,   630,   661, -1072,   791, -1072,   -30,
     646,   656,   794, -1072, -1072,    93, -1072, -1072, -1072, -1072,
     769,   770,   139, -1072,   663, -1072, -1072,   665,   669,   671,
     976,   976,   674,   675,   677,   688,   689,   690,   696,   697,
     700,   701,   976,   706,   712,   413, -1072,   815,   422,   827,
     832,   836,   841,   851,   858, -1072, -1072, -1072,   863,   866,
     874, -1072,   876, -1072,   877,   878,   729, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072,   792,   843, -1072,   736, -1072,   734, -1072,
   -1072,   751,   752,   753,  1035,  1035,   755,   756,   758,   762,
     765,   779,   780,   783,   785,   786,  1035,  1627, -1072,   258,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072,   276, -1072,   153,    41,   290, -1072,
   -1072, -1072,   292, -1072, -1072,   299, -1072, -1072,   901, -1072,
     305, -1072, -1072, -1072, -1072, -1072,   854,   903, -1072, -1072,
     873,   915, -1072, -1072,   875,   920, -1072, -1072,   972,   973,
     974,   980,   982,   983,   984,   990,   992,   993,   999,  1001,
    1002,  1003,  1005,  1006,  1008,  1009,  1013,  1014,  1026,  1027,
    1028,  1029,  1034,  1045,  1046,   893,   951, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072,    75,   321,    75,  1039,   321,  1042,
    1012,  1044,    21,  1049,  1050,  1020,  1022,  1544,  1054,  1055,
      75,  1061,  1524,  1065,   926,   937,  1037,   100,  1104, -1072,
   -1072,  1087,   439,   948, -1072, -1072,   949,    84,  1060,   952,
      97,  1063,   976, -1072, -1072, -1072,   953,  1095,  1097,  1098,
    1105,  1109,    75,    75,    75,  1111,  1112,  1114,  1115,  1072,
     962,    75,  1424,   103,  1085,  1137,  1038, -1072, -1072, -1072,
     307,  1040,   384,  1041, -1072, -1072,  1052,   384,  1053, -1072,
   -1072,   267, -1072, -1072, -1072,  1089,  1000, -1072,  1090,    48,
   -1072,    37, -1072,   428, -1072,  1011,  1023,    57,    42,    73,
    1056,    35, -1072, -1072,   976,   976,   976,   976,  1033,   478,
     976,   976,   976,   976,   976,   976,   976,   976,   976,   976,
     423,   976,   976,   976,   976,   976, -1072,   976, -1072, -1072,
    1119,  1235, -1072,   844,  1153,    75,  1163,  1169,  1171,  1176,
    1182,  1184,    75,  1186,  1187,  1194,   172, -1072,  1121, -1072,
    1035,  1035,  1035,   267,  1048,   499,  1035,  1035,  1035,  1035,
    1035,  1035,  1035,  1035,  1035,  1035,   430,  1035,  1035,  1035,
    1035,  1035,  1025,   653,   175,   180, -1072, -1072, -1072,   976,
     264,    31,    82,  1057,   439,  1058,  1595,   183,    75,  1043,
     196, -1072,  1122, -1072,  1125, -1072,  1128,  1204,  1209,  1211,
    1212,  1217,  1218,  1219,  1225,  1240,  1241,  1242,  1243,  1247,
    1248,  1257,    75,    75,  1265,   953,    75,    75,  1267,  1268,
      75,  1276,    75,    75,  1138,  1565, -1072, -1072, -1072, -1072,
    1198,  1287, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
    1281,    23, -1072, -1072, -1072, -1072,  1143, -1072, -1072,  1141,
   -1072, -1072, -1072, -1072,  1144, -1072,  1290,  1154,  1152,  1165,
     976, -1072, -1072, -1072, -1072, -1072,   223,  1183, -1072, -1072,
     224,  1185,  1253, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072,  1149, -1072, -1072,
   -1072,   228, -1072,  1275,  1283, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072,   404,  1188,  1244,  1245,  1301,  1250,   384,
    1306,  1197,   384, -1072,  1344,  1353,  1216, -1072,   497,  1378,
   -1072, -1072, -1072,  1035, -1072, -1072, -1072,  1379, -1072,   311,
   -1072, -1072, -1072,  1221, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072,   -29,     9, -1072,  1334,   976,  1336,
      16,   494,   623,   933,  1633,   312,   509,   535,   541,   612,
     634,   641,   654,   660,   691,   698, -1072,   478,   478,  1033,
    1033,  1282,   737,   976, -1072,  1342,  1259, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
     229, -1072,   763,  1213,  1620,  1232,   774,   781,   871,   881,
     890,   942,   970,  1059,  1070,  1091, -1072,   499,   499,  1048,
    1048,  1282, -1072, -1072, -1072,   231, -1072,   232,  1103,    41,
    1239, -1072, -1072,    43,   976, -1072, -1072, -1072, -1072, -1072,
   -1072,   233, -1072, -1072, -1072,   255, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
    1237, -1072, -1072, -1072,  1352, -1072, -1072,  1252,  1402, -1072,
   -1072,  1314, -1072,   198, -1072,   199, -1072,  1355, -1072,   320,
   -1072, -1072, -1072, -1072, -1072, -1072,   384,   307,  1310,   384,
    1312,  1313, -1072,  1246, -1072, -1072,  1423,   443,  1035,  1333,
      75,   428, -1072,   791,   791,   791,    73, -1072,   384, -1072,
    1425,  1381,  1426,  1409,   976, -1072,   976,   976,   976, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
    1269,  1431,   976, -1072, -1072, -1072,  1035,  1035, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072,    31, -1072, -1072, -1072,   976,  1118, -1072, -1072,
    1411, -1072,  1154,   976, -1072, -1072,   272, -1072,   280,  1270,
    1149, -1072, -1072,  1329,  1338,  1343,   384,  1289,   384,   384,
   -1072,   976, -1072,  1438, -1072, -1072, -1072,  1292,   101,   235,
     574,    71,  1293,   976, -1072,   976,  1300,    64,  1448,  1134,
    1140,  1633, -1072, -1072,  1460,  1166,  1223,  1229, -1072, -1072,
    1455,  1472, -1072, -1072,  1363, -1072,   384,   384,   384,  1364,
   -1072,  1311,  1315,  1495, -1072,   791, -1072, -1072, -1072,   384,
   -1072,  1527,  1533,  1451,  1316,  1453,  1380, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072,   976, -1072,    29,  1377, -1072,  1386,
     384, -1072, -1072, -1072,   666,  1321, -1072, -1072, -1072,  1465,
    1332,   976,  1577,  1435, -1072,   384,    79,  1339, -1072, -1072,
   -1072,  1479,  1633,   905, -1072,  1341,  1408,  1412, -1072, -1072,
   -1072,  1633, -1072,   384,   384,  1413, -1072,   384, -1072
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   423,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     4,    29,    30,    45,
      46,    47,    44,     5,     6,     7,    12,     9,    10,    11,
       8,    13,    14,    15,    16,    17,    18,    19,    23,    25,
      24,    20,    21,    22,    26,    27,    28,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,   392,     0,     0,   208,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   252,   299,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   131,
       0,     0,     0,     0,     0,     0,   379,     0,     0,     0,
      75,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     213,     0,   203,     0,   219,     0,     0,   426,     0,     0,
       0,    57,     0,    63,     0,    69,     0,     0,   463,     0,
       1,     3,   453,     0,   566,     0,     0,     0,     0,     0,
     557,     0,     0,     0,   558,     0,     0,     0,     0,   441,
     452,     0,   442,   447,   445,   448,   446,   443,   444,   449,
     450,   434,   435,   436,   437,   438,   439,   440,   461,     0,
       0,     0,   455,   460,     0,   457,   456,   458,     0,     0,
     389,     0,   385,     0,     0,   211,   212,     0,    81,     0,
      48,   402,     0,     0,   396,     0,     0,     0,     0,   118,
       0,   541,     0,   546,   525,     0,     0,     0,     0,     0,
     538,   539,     0,     0,     0,     0,     0,     0,   559,   534,
       0,     0,   545,   540,   524,     0,     0,   544,   542,     0,
     304,   340,   329,   305,   306,   307,   308,   309,   310,   311,
     312,   313,   314,   315,   316,   317,   318,   319,   320,   321,
     322,   323,   324,   325,   326,   327,   328,   330,   331,   332,
     333,   334,   335,   336,   337,   338,   339,   341,   342,   343,
     248,     0,   301,     0,   265,     0,     0,   262,     0,     0,
       0,     0,     0,   284,     0,     0,     0,     0,   278,     0,
       0,   122,     0,     0,     0,     0,    83,     0,   499,     0,
       0,     0,     0,   553,   552,     0,   408,   409,   410,   411,
       0,     0,     0,   177,     0,    88,    89,     0,     0,    87,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   370,     0,     0,     0,
       0,     0,     0,     0,     0,   504,   505,   506,     0,     0,
       0,   547,     0,   513,     0,     0,     0,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   237,   239,
     240,   241,   242,   243,   244,   245,   236,   238,   246,   247,
     381,   378,    78,    73,     0,    54,     0,    79,     0,   149,
     150,     0,     0,   172,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   424,   148,     0,
     347,   352,   348,   349,   350,   351,   353,   354,   355,   356,
     357,   358,   359,   360,     0,    50,     0,     0,     0,   216,
     217,   218,     0,   206,   207,     0,   224,   221,     0,   432,
       0,   431,   433,   428,   372,    60,    55,     0,    51,    66,
      61,     0,    52,    72,    67,     0,    53,   367,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   466,   467,   468,   469,
     470,   471,   472,   473,   474,   475,   476,   477,   478,   479,
     480,   481,   482,   483,   484,   493,   485,   486,   487,   488,
     489,   490,   491,   492,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   383,
     384,     0,     0,     0,    82,    49,     0,     0,     0,     0,
       0,     0,     0,   116,   117,   253,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   250,     0,   264,   260,   261,
     293,     0,   293,     0,   282,   283,     0,   293,     0,   276,
     277,     0,   120,   121,   113,     0,     0,    84,     0,     0,
     143,     0,   144,     0,   139,     0,     0,     0,     0,     0,
       0,     0,   175,   176,     0,     0,     0,     0,    95,    96,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    85,     0,   368,   369,
       0,     0,   373,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    76,    74,    80,
       0,     0,     0,     0,   156,   157,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   174,   201,   202,     0,
       0,   193,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    58,    56,    64,    62,    70,    68,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   495,   494,   562,   259,
       0,     0,   563,   564,   565,   561,   567,   516,   519,   518,
       0,     0,   517,   520,   521,   554,     0,   555,   451,     0,
     568,   526,   543,   459,     0,   393,     0,   389,     0,     0,
       0,   497,   210,   209,   405,   400,     0,     0,   399,   394,
       0,     0,     0,   556,   507,   548,   549,   522,   523,   528,
     531,   529,   536,   537,   527,   533,   532,     0,   535,   303,
     300,     0,   249,     0,     0,   288,   295,   289,   294,   291,
     296,   290,   292,     0,     0,     0,   270,     0,     0,   293,
       0,     0,   293,   256,     0,     0,     0,   115,     0,     0,
     132,   141,   142,     0,   146,   127,   126,     0,   128,     0,
     125,   129,   130,     0,   135,   133,   550,   551,   407,   418,
     420,   421,   422,   419,     0,   412,   416,     0,     0,     0,
       0,     0,     0,     0,   111,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    86,    91,    90,    92,
      93,    94,     0,     0,   376,     0,     0,   503,   511,   496,
     502,   508,   509,   500,   510,   515,   501,   498,   514,   380,
       0,    77,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   147,   152,   151,   153,
     154,   155,   425,   346,   344,     0,   361,     0,     0,     0,
       0,   198,   199,     0,     0,   215,   214,   205,   204,   223,
     220,     0,   560,   430,   427,     0,    59,    65,    71,   569,
     570,   571,   572,   573,   574,   575,   576,   577,   578,   579,
     580,   581,   582,   583,   584,   585,   586,   587,   588,   589,
     590,   591,   592,   593,   594,   595,   464,   465,   258,   257,
     597,   599,   601,   600,     0,   454,   462,     0,     0,   391,
     390,     0,   401,     0,   395,     0,   119,     0,   365,     0,
     302,   251,   266,   298,   297,   263,   293,   293,     0,   293,
       0,     0,   281,     0,   255,   254,     0,     0,     0,     0,
       0,     0,   137,     0,     0,     0,     0,   406,   293,   417,
       0,     0,     0,     0,     0,   107,     0,     0,     0,   110,
      97,    98,    99,   100,   101,   102,   103,   104,   105,   106,
       0,     0,     0,   374,   382,   168,     0,     0,   173,   158,
     159,   160,   161,   162,   163,   164,   165,   166,   167,   345,
     362,   200,   192,   191,   195,   196,     0,     0,   222,   429,
       0,   596,   389,     0,   386,   403,     0,   397,     0,     0,
       0,   530,   267,     0,     0,     0,   293,     0,   293,   293,
     279,     0,   114,     0,   145,   512,   124,     0,     0,     0,
       0,   413,     0,     0,   180,     0,   188,     0,     0,     0,
       0,   112,   371,   377,     0,     0,     0,     0,   197,   598,
       0,     0,   404,   398,     0,   366,   293,   293,   293,     0,
     287,     0,     0,     0,   171,     0,   140,   136,   134,   293,
     414,     0,     0,     0,   183,     0,     0,   179,   108,   109,
     375,   169,   170,   194,     0,   387,   293,   272,   268,   271,
     293,   285,   280,   123,     0,     0,   182,   181,   187,     0,
     185,     0,     0,     0,   364,   293,     0,     0,   138,   415,
     184,     0,   190,     0,   388,     0,   273,     0,   286,   186,
     178,   189,   363,   293,   293,   274,   269,   293,   275
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1072, -1072,  1503, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,  -322, -1072,
   -1072, -1072, -1072,  -111,  -225, -1072, -1072,  1256, -1072,   511,
   -1072, -1072, -1072, -1072, -1072, -1072,  -966,  -606,  -115,  -599,
   -1072, -1072, -1072,  1414,  -243, -1072, -1072, -1072, -1072,   604,
   -1072, -1072,   840, -1072, -1072,   995, -1072, -1072,   859, -1072,
   -1072,  -112,   -24,   891,  1047, -1072, -1072,  1277, -1072, -1072,
   -1071, -1072, -1072,  1274, -1072, -1072,  1273,  -971,  -569, -1072,
   -1072,   991, -1072,  1458,   883, -1072,   477, -1072, -1072, -1072,
   -1072,  1238, -1072, -1072, -1072, -1072, -1072, -1072, -1072,  1394,
    -733, -1072, -1072, -1072, -1072, -1072,   968, -1072,   551,  -815,
   -1072, -1072, -1072, -1072, -1072,   879, -1072,   -33,  1066, -1072,
   -1072,  1064, -1072, -1072,   850, -1072,  -536, -1072,   -92, -1072,
    1497, -1072, -1072, -1072, -1072, -1072, -1072, -1072,  -104, -1072,
   -1072,  -120,  -570, -1072, -1072, -1072, -1072,  -107,   -93,   -84,
     -82,   -72, -1072, -1072,  -100,   -71, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072,   -67, -1072, -1072, -1072, -1072, -1072,
     -65,   -60,   -54,   -59,   -55,   -50, -1072, -1072, -1072, -1072,
    -113,  -110,   -89,   -83,   -48,   -47,   -41, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072, -1072,
   -1072, -1072, -1072, -1072, -1072,   837, -1072,  -519
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    44,    45,    46,    47,    48,    49,    50,    51,    52,
     152,   154,   156,   131,    53,    54,    55,    56,   363,   895,
      57,   324,    58,   228,   229,    59,   320,   321,   869,   870,
      60,   327,  1055,  1054,  1137,   873,   629,   630,   631,   632,
     438,    61,    62,   342,   343,  1147,  1223,    63,   720,   721,
      64,   462,   463,    65,   214,   215,    66,   458,   459,    67,
     465,   469,   110,   856,   772,    68,   306,   307,   308,   844,
    1122,    69,   317,   318,    70,   312,   313,   845,  1123,    71,
     259,   260,    72,   439,   440,    73,  1028,  1029,    74,    75,
     365,   366,    76,    77,   368,    78,    79,    80,   211,   212,
     568,    81,    82,    83,    84,   335,   336,   884,   885,   886,
      85,   134,   712,    86,   470,   471,   179,   180,   181,    87,
     203,   204,    88,    89,   515,   516,   768,   387,   388,   389,
     390,   391,   392,   393,   394,   395,   396,   397,   398,   399,
     400,   401,   402,   872,   403,   404,   405,   182,   183,   184,
     185,   186,   268,   269,   406,   443,   272,   273,   274,   275,
     276,   277,   278,   279,   444,   281,   282,   283,   284,   285,
     445,   446,   447,   448,   449,   450,   407,   292,   293,   337,
     408,   409,   187,   188,   453,   189,   190,   299,   472,   191,
     192,   193,   194,   195,   196,   197,   207,   517,   518,   519,
     520,   521,   522,   523,   524,   525,   526,   527,   528,   529,
     530,   531,   532,   533,   534,   535,   536,   537,   538,   539,
     540,   541,   542,   543,   787,  1011,   781,   782
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       128,   129,   263,   584,   322,   262,   216,   137,   338,   270,
     773,   339,   146,   149,   150,   386,   264,   261,   157,   437,
     294,   460,   205,   861,   791,   265,   295,   266,   648,   649,
     862,   441,   441,   466,   461,   442,   442,   267,   271,   206,
     660,   846,   280,   848,   286,   671,   451,   451,   851,   287,
     289,   464,   452,   452,   290,   288,   819,   820,   821,   291,
     202,   296,   297,   871,  1018,   828,  1124,   813,   298,   888,
    1059,   863,   836,   779,   961,  1010,    90,   418,   209,   107,
     860,   838,   219,   962,   717,   835,  1104,  1138,  1139,  1140,
      92,   419,   879,   718,   107,  1105,  1198,   584,   572,   643,
     420,   602,   132,   209,   577,   300,   569,   171,   633,  1056,
     578,   101,   840,   421,   876,    96,  1184,  1063,   766,   583,
     133,   422,   332,   107,   837,    94,    95,   767,   107,   106,
     418,   423,   839,  1176,   333,   566,   111,  1064,   877,   918,
      99,   112,   117,   798,   419,   102,   925,   334,   104,   100,
     880,   118,   799,   420,   879,   573,   879,   210,   603,   340,
     301,   340,   113,  1236,  1057,   634,   421,   107,   843,   963,
     227,   642,   841,   889,   422,   340,   707,   708,   709,   710,
     107,   711,   210,   379,   423,   716,   107,   424,   425,   108,
     109,   567,   972,   426,   427,   428,   429,   430,   431,   432,
     433,   434,  1185,  1058,   126,   127,   881,   842,   435,  1204,
     882,   883,   880,   643,   880,   302,   994,   995,    91,   780,
     998,   999,  1012,   220,  1002,   964,  1004,  1005,   890,  1213,
     864,   638,    93,   144,   145,   719,   997,  1106,   147,   148,
     424,   425,   436,   300,   628,  1227,   426,   427,   428,   429,
     430,   431,   432,   433,   434,   107,   103,  1186,   107,   105,
     812,   435,   300,   107,   418,  1179,   107,  1177,   881,   341,
     881,   341,   882,   883,   882,   883,   622,   805,   419,   107,
    1040,   107,   107,  1043,   114,   341,   121,   420,   639,   413,
     809,   300,   300,   300,   300,   436,   830,   628,   301,   476,
     421,   480,   484,   300,    97,    98,   300,   300,   422,   694,
     695,   300,   300,   835,   300,   300,   300,   301,   423,   853,
     122,   706,   891,   892,   893,   894,  1059,   319,   896,   897,
     898,   899,   900,   901,   902,   903,   904,   905,   300,   907,
     908,   909,   910,   911,   414,   912,   301,   301,   301,   301,
     836,   916,   837,   410,   477,   300,   481,   485,   301,   838,
     839,   301,   301,   300,   769,   929,   301,   301,   954,   301,
     301,   301,   411,   956,   424,   425,   970,   303,   309,  1160,
     426,   427,   428,   429,   430,   431,   432,   433,   434,   974,
     840,  1115,  1117,   301,   314,   435,   713,   958,   303,   415,
     841,   455,   467,   473,   474,   309,   854,   855,   608,   478,
     301,   482,   486,   487,   713,   614,  1022,  1024,   301,   115,
     116,  1030,  1084,   314,  1099,  1100,  1108,   836,   722,   436,
     724,   628,   624,   619,   580,   842,   838,   726,   304,   310,
     581,   124,   328,   729,   123,   668,   843,  1033,  1109,  1051,
    1068,   107,   370,   714,   672,   315,  1034,   959,  1120,   304,
     770,   771,   125,   960,   221,  1162,   310,   840,  1125,   130,
    1127,   715,   222,  1163,   865,  1132,   305,   311,  1021,   135,
     216,   871,   224,   227,   315,   723,   866,   725,   136,  1142,
     225,   205,   867,   316,   727,   263,   364,   305,   262,   138,
     730,   329,   270,   139,   311,   673,  1052,  1069,   206,   264,
     261,   330,   868,   294,  1135,  1121,   119,   120,   265,   295,
     266,   151,   316,   843,   153,   338,   227,   155,   339,   202,
     267,   271,   861,   861,   861,   280,   162,   286,   198,   862,
     862,   862,   287,   289,   140,   141,   208,   290,   288,   217,
     142,   143,   291,   806,   296,   297,   810,  1169,   213,  1171,
    1172,   298,   661,   662,   663,   664,  1061,   665,   218,   707,
     708,   709,   710,   223,   711,   932,   933,   934,   226,   831,
     227,   936,   937,   938,   939,   940,   941,   942,   943,   944,
     945,  1081,   947,   948,   949,   950,   951,  1197,   861,  1199,
     158,   159,   460,   418,   230,   862,  1178,   319,   323,   441,
    1205,   325,   326,   442,   969,   461,   364,   419,   906,   663,
     664,   367,   665,   412,   451,   946,   420,  1214,   416,   417,
     452,  1217,   464,   661,   662,   663,   664,   475,   665,   421,
     709,   710,  1107,   711,   457,   479,  1226,   422,   661,   662,
     663,   664,   483,   665,   544,   545,   546,   423,   231,   547,
     548,   549,   930,   550,  1235,   551,   552,   553,  1238,   554,
     555,   556,   557,   200,   661,   662,   663,   664,   558,   665,
     661,   662,   663,   664,   559,   665,   560,   562,   563,  1065,
     955,   957,   561,   232,   233,   418,   564,   201,  1218,   565,
     234,   571,   576,   971,  1070,   574,   975,   235,   575,   419,
     582,   579,   586,   424,   425,   587,   588,   589,   420,   426,
     427,   428,   429,   430,   431,   432,   433,   434,   585,   590,
    1071,   421,   605,   591,   435,   252,  1072,   592,   593,   422,
     606,   594,  1148,   254,  1149,  1150,  1151,  1047,  1049,   423,
     595,   661,   662,   663,   664,   596,   665,   597,   598,   256,
    1154,  1066,   661,   662,   663,   664,   599,   665,   436,   600,
     628,   257,   601,   661,   662,   663,   664,   258,   665,   604,
     661,   662,   663,   664,  1157,   665,   610,   621,   607,   177,
     178,  1161,   611,   661,   662,   663,   664,   612,   665,   661,
     662,   663,   664,   613,   665,   424,   425,  1073,   616,  1173,
     617,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     418,  1181,   584,  1182,   618,   626,   435,   625,   637,  1074,
     661,   662,   663,   664,   419,   665,  1075,   661,   662,   663,
     664,   635,   665,   420,   344,   661,   662,   663,   664,  1076,
     665,   636,   640,   641,   627,  1077,   421,   644,   345,   645,
     436,   674,   628,   646,   422,   647,   675,   346,   650,   651,
     676,   652,  1212,   344,   423,   677,   661,   662,   663,   664,
     347,   665,   653,   654,   655,   678,  1078,   345,   348,  1222,
     656,   657,   679,  1079,   658,   659,   346,   680,   349,   666,
     681,  1231,   707,   708,   709,   710,   667,   711,   682,   347,
     683,   684,   685,   707,   708,   709,   710,   348,   711,   687,
     707,   708,   709,   710,   686,   711,   688,   349,   690,   689,
     424,   425,  1080,  1133,   344,   728,   426,   427,   428,   429,
     430,   431,   432,   433,   434,   691,   692,   693,   345,   696,
     697,   435,   698,   670,   350,   351,   699,   346,  1085,   700,
     352,   353,   354,   355,   356,   357,   358,   359,   360,  1089,
     347,  1155,  1156,   701,   702,   361,  1090,   703,   348,   704,
     705,   731,   915,   350,   351,   436,   732,   628,   349,   352,
     353,   354,   355,   356,   357,   358,   359,   360,   734,  1116,
     733,  1118,   735,   736,   361,   344,   737,   738,   739,   362,
     707,   708,   709,   710,   740,   711,   741,   742,   743,   345,
     707,   708,   709,   710,   744,   711,   745,   746,   346,   707,
     708,   709,   710,   747,   711,   748,   749,   750,   362,   751,
     752,   347,   753,   754,   350,   351,   369,   755,   756,   348,
     352,   353,   354,   355,   356,   357,   358,   359,   360,   349,
     757,   758,   759,   760,   418,   361,  1091,   370,   761,   371,
     372,  1067,   661,   662,   663,   664,  1092,   665,   419,   762,
     763,   707,   708,   709,   710,  1093,   711,   420,   764,   765,
     234,   774,   373,   374,   776,   777,   778,   235,  1230,   362,
     421,   783,   784,   785,   328,   786,   789,   790,   422,   707,
     708,   709,   710,   792,   711,   350,   351,   794,   423,   795,
     797,   352,   353,   354,   355,   356,   357,   358,   359,   360,
     375,   796,   376,   254,   377,   333,   361,  1094,   800,   801,
     378,   803,   804,   807,   379,   808,   811,   814,   334,   815,
     816,   780,   380,   381,   382,   826,   827,   817,   383,   384,
     385,   818,   213,   822,   823,  1095,   824,   825,   832,   468,
     362,   833,   857,   859,   424,   425,   834,   665,   847,   849,
     426,   427,   428,   429,   430,   431,   432,   433,   434,   160,
     850,   852,   711,   858,   887,   435,     1,     2,   707,   708,
     709,   710,   913,   711,   874,   917,     3,     4,     5,   707,
     708,   709,   710,     6,   711,   919,   875,     7,   952,     8,
       9,   920,    10,   921,    11,    12,    13,    14,   922,   436,
     707,   708,   709,   710,   923,   711,   924,    15,   926,   927,
      16,  1008,   661,   662,   663,   664,   928,   665,   931,   976,
     966,   968,   977,    17,  1096,   978,   979,   661,   662,   663,
     664,   980,   665,   981,   982,  1097,    18,    19,    20,   983,
     984,   985,    21,   661,   662,   663,   664,   986,   665,   661,
     662,   663,   664,    22,   665,    23,  1098,    24,    25,    26,
      27,    28,   987,   988,   989,   990,    29,    30,  1101,   991,
     992,    31,    32,    33,    34,   707,   708,   709,   710,   993,
     711,    35,    36,  1158,    37,     1,     2,   996,    38,  1000,
    1001,    39,    40,    41,    42,     3,     4,     5,  1003,  1188,
    1009,  1006,     6,  1010,  1015,  1189,     7,  1016,     8,     9,
    1014,    10,  1017,    11,    12,    13,    14,  1019,   567,  1027,
      43,  1086,   707,   708,   709,   710,    15,   711,  1031,    16,
    1020,  1191,   707,   708,   709,   710,  1032,   711,   661,   662,
     663,   664,    17,   665,   661,   662,   663,   664,  1023,   665,
    1025,  1035,  1036,  1037,  1038,    18,    19,    20,  1039,  1041,
    1042,    21,   661,   662,   663,   664,  1044,   665,   661,   662,
     663,   664,    22,   665,    23,  1045,    24,    25,    26,    27,
      28,  1046,  1048,  1050,  1053,    29,    30,  1060,  1192,  1062,
      31,    32,    33,    34,  1193,  1082,    -1,  1088,   914,   231,
      35,    36,  1103,    37,  1110,  1111,  1113,    38,  1119,  1130,
      39,    40,    41,    42,   200,   170,  1026,  1112,  1126,   171,
    1128,  1129,  1083,   661,   662,   663,   664,  1131,   665,  1143,
    1145,  1146,  1152,  1159,   232,   233,   172,  1166,   201,    43,
    1164,   234,   707,   708,   709,   710,  1167,   711,   235,   236,
     237,  1168,  1170,   238,   239,  1175,   240,   241,  1180,  1194,
     242,   243,   244,   245,   246,   247,   248,  1183,   249,   250,
     251,  1196,  1200,  1208,  1201,  1210,   252,  1114,  1202,   173,
     174,  1211,   253,  1209,   254,  1215,  1219,  1220,  1225,   255,
     661,   662,   663,   664,  1216,   665,  1134,   175,   176,  1221,
     256,  1229,  1228,   163,   164,   165,   166,   167,   168,   169,
     199,  1232,   257,   213,   200,   170,  1233,   161,   258,   171,
    1234,  1237,   456,   163,   164,   165,   166,   167,   168,   169,
     177,   178,  1136,  1102,   967,   170,   172,   802,   201,   171,
     661,   662,   663,   664,  1144,   665,   623,   707,   708,   709,
     710,   965,   711,   609,   935,   615,   172,   661,   662,   663,
     664,   620,   665,   829,   454,   775,   953,  1165,   369,   661,
     662,   663,   664,   669,   665,   570,   878,  1141,   973,   173,
     174,   661,   662,   663,   664,  1007,   665,   331,  1013,   370,
       0,   371,   372,   788,  1153,     0,   793,   175,   176,   173,
     174,  1174,     0,     0,   661,   662,   663,   664,     0,   665,
       0,  1187,   234,     0,   373,   374,     0,   175,   176,   235,
       0,     0,     0,  1190,     0,     0,   328,     0,     0,     0,
     177,   178,     0,     0,     0,  1195,   661,   662,   663,   664,
       0,   665,   661,   662,   663,   664,     0,   665,     0,     0,
     177,   178,   375,     0,   376,   254,   377,   333,  1203,     0,
       0,     0,   378,     0,     0,     0,   379,     0,     0,     0,
     334,     0,     0,     0,   380,   381,   382,     0,     0,     0,
     383,   384,   385,     0,   213,     0,   661,   662,   663,   664,
    1206,   665,     0,     0,     0,     0,  1207,   488,   489,   490,
     491,   492,   493,   494,   495,   496,   497,   498,   499,   500,
     501,   502,   503,   504,   505,   506,   507,   508,     0,     0,
       0,   509,   510,     0,   511,   512,   513,   514,  1087,   707,
     708,   709,   710,     0,   711,     0,   707,   708,   709,   710,
    1224,   711,   661,   662,   663,   664,     0,   665
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   109,   228,   115,   109,    98,    31,   121,   109,
     546,   121,    36,    37,    38,   127,   109,   109,    42,   134,
     109,   141,    93,   629,   560,   109,   109,   109,   350,   351,
     629,   135,   136,   145,   141,   135,   136,   109,   109,    93,
     362,   610,   109,   612,   109,   367,   135,   136,   617,   109,
     109,   143,   135,   136,   109,   109,   592,   593,   594,   109,
      93,   109,   109,   633,   797,   601,  1037,   586,   109,    34,
     885,    34,    43,    52,    43,    52,    52,    29,     4,    83,
      32,    52,    52,    52,    43,     6,    43,  1053,  1054,  1055,
      52,    43,    83,    52,    83,    52,  1167,   322,   138,   342,
      52,   138,    34,     4,   195,    83,    32,    25,   138,   138,
     201,    52,    83,    65,    57,   193,    52,   101,    43,    32,
      52,    73,    80,    83,    45,   193,   194,    52,    83,   193,
      29,    83,    53,    32,    92,   138,   193,   121,    81,   675,
      34,   193,    34,    43,    43,    83,   682,   105,    83,    43,
     141,    43,    52,    52,    83,   195,    83,    83,   195,    22,
     138,    22,   193,  1234,   193,   195,    65,    83,   139,   138,
      83,    32,    93,   138,    73,    22,   139,   140,   141,   142,
      83,   144,    83,   101,    83,    32,    83,   139,   140,   193,
     194,   194,   728,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   138,   194,   193,   194,   197,   128,   160,  1175,
     201,   202,   141,   456,   141,   193,   752,   753,   194,   198,
     756,   757,   199,   193,   760,   194,   762,   763,   193,   200,
     193,   138,   194,   193,   194,   194,   755,   194,   193,   194,
     139,   140,   194,    83,   196,  1216,   145,   146,   147,   148,
     149,   150,   151,   152,   153,    83,   194,   193,    83,   194,
     582,   160,    83,    83,    29,   194,    83,    32,   197,   132,
     197,   132,   201,   202,   201,   202,    32,   193,    43,    83,
     849,    83,    83,   852,   193,   132,   194,    52,   195,    83,
     193,    83,    83,    83,    83,   194,   193,   196,   138,    83,
      65,    83,    83,    83,   193,   194,    83,    83,    73,   424,
     425,    83,    83,     6,    83,    83,    83,   138,    83,    52,
     193,   436,   644,   645,   646,   647,  1141,    83,   650,   651,
     652,   653,   654,   655,   656,   657,   658,   659,    83,   661,
     662,   663,   664,   665,   138,   667,   138,   138,   138,   138,
      43,   673,    45,   193,   138,    83,   138,   138,   138,    52,
      53,   138,   138,    83,    43,   193,   138,   138,   193,   138,
     138,   138,   193,   193,   139,   140,   193,    22,    22,  1112,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   193,
      83,   193,   193,   138,    22,   160,   138,   719,    22,   193,
      93,   193,   193,   193,   193,    22,   139,   140,    32,   193,
     138,   193,   193,   193,   138,    32,   193,   193,   138,   193,
     194,   193,   193,    22,   193,   193,   193,    43,   138,   194,
     138,   196,    32,    32,   195,   128,    52,   138,    83,    83,
     201,   193,    61,   138,    34,    32,   139,    43,   193,   138,
     138,    83,    24,   195,    32,    83,    52,   193,   138,    83,
     139,   140,   193,   199,   193,   193,    83,    83,  1037,    83,
    1039,   195,   201,   193,    46,    32,   121,   121,   800,   194,
     572,  1051,   193,    83,    83,   195,    58,   195,   194,  1058,
     201,   562,    64,   121,   195,   602,    83,   121,   602,   193,
     195,   120,   602,    34,   121,    83,   195,   195,   562,   602,
     602,   130,    84,   602,  1050,   195,   193,   194,   602,   602,
     602,    83,   121,   139,    83,   638,    83,    83,   638,   562,
     602,   602,  1138,  1139,  1140,   602,   193,   602,   193,  1138,
    1139,  1140,   602,   602,   193,   194,    23,   602,   602,    43,
     193,   194,   602,   577,   602,   602,   580,  1126,   119,  1128,
    1129,   602,   139,   140,   141,   142,   888,   144,   193,   139,
     140,   141,   142,    83,   144,   690,   691,   692,    83,   603,
      83,   696,   697,   698,   699,   700,   701,   702,   703,   704,
     705,   913,   707,   708,   709,   710,   711,  1166,  1204,  1168,
     193,   194,   722,    29,    34,  1204,    32,    83,    39,   713,
    1179,    43,   193,   713,   726,   722,    83,    43,   195,   141,
     142,    83,   144,   127,   713,   195,    52,  1196,    52,   193,
     713,  1200,   724,   139,   140,   141,   142,   127,   144,    65,
     141,   142,   964,   144,   198,   127,  1215,    73,   139,   140,
     141,   142,   127,   144,    34,    34,    34,    83,     5,    34,
      34,    34,   686,    34,  1233,    34,    34,    34,  1237,    34,
      34,    34,   138,    20,   139,   140,   141,   142,   195,   144,
     139,   140,   141,   142,    34,   144,    34,   138,   195,   195,
     714,   715,    34,    40,    41,    29,   195,    44,    32,    83,
      47,    34,    83,   727,   195,   193,   730,    54,   193,    43,
      34,    83,    34,   139,   140,    34,    34,    34,    52,   145,
     146,   147,   148,   149,   150,   151,   152,   153,    83,    34,
     195,    65,    83,    34,   160,    82,   195,    34,    34,    73,
      83,    34,  1064,    90,  1066,  1067,  1068,   858,   863,    83,
      34,   139,   140,   141,   142,    34,   144,    34,    34,   106,
    1082,   138,   139,   140,   141,   142,    34,   144,   194,    34,
     196,   118,    34,   139,   140,   141,   142,   124,   144,    34,
     139,   140,   141,   142,  1106,   144,   138,   194,    83,   136,
     137,  1113,    83,   139,   140,   141,   142,   138,   144,   139,
     140,   141,   142,    83,   144,   139,   140,   195,    83,  1131,
     138,   145,   146,   147,   148,   149,   150,   151,   152,   153,
      29,  1143,  1047,  1145,    83,   195,   160,    34,    34,   195,
     139,   140,   141,   142,    43,   144,   195,   139,   140,   141,
     142,   195,   144,    52,    29,   139,   140,   141,   142,   195,
     144,   195,    83,    83,   193,   195,    65,   194,    43,   194,
     194,    34,   196,   194,    73,   194,    34,    52,   194,   194,
      34,   194,  1194,    29,    83,    34,   139,   140,   141,   142,
      65,   144,   194,   194,   194,    34,   195,    43,    73,  1211,
     194,   194,    34,   195,   194,   194,    52,    34,    83,   193,
      34,  1223,   139,   140,   141,   142,   194,   144,    34,    65,
      34,    34,    34,   139,   140,   141,   142,    73,   144,   127,
     139,   140,   141,   142,   195,   144,    83,    83,   194,   193,
     139,   140,   195,  1048,    29,    34,   145,   146,   147,   148,
     149,   150,   151,   152,   153,   194,   194,   194,    43,   194,
     194,   160,   194,   138,   139,   140,   194,    52,   195,   194,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   195,
      65,  1086,  1087,   194,   194,   160,   195,   194,    73,   194,
     194,   127,   138,   139,   140,   194,    83,   196,    83,   145,
     146,   147,   148,   149,   150,   151,   152,   153,    83,  1023,
     127,  1025,   127,    83,   160,    29,    34,    34,    34,   194,
     139,   140,   141,   142,    34,   144,    34,    34,    34,    43,
     139,   140,   141,   142,    34,   144,    34,    34,    52,   139,
     140,   141,   142,    34,   144,    34,    34,    34,   194,    34,
      34,    65,    34,    34,   139,   140,     3,    34,    34,    73,
     145,   146,   147,   148,   149,   150,   151,   152,   153,    83,
      34,    34,    34,    34,    29,   160,   195,    24,    34,    26,
      27,   138,   139,   140,   141,   142,   195,   144,    43,    34,
      34,   139,   140,   141,   142,   195,   144,    52,   195,   138,
      47,    52,    49,    50,    52,    83,    52,    54,   193,   194,
      65,    52,    52,    83,    61,    83,    52,    52,    73,   139,
     140,   141,   142,    52,   144,   139,   140,    52,    83,   193,
      83,   145,   146,   147,   148,   149,   150,   151,   152,   153,
      87,   194,    89,    90,    91,    92,   160,   195,    34,    52,
      97,   193,   193,    83,   101,   193,    83,    52,   105,    52,
      52,   198,   109,   110,   111,    83,   194,    52,   115,   116,
     117,    52,   119,    52,    52,   195,    52,    52,    83,   126,
     194,    34,    83,    83,   139,   140,   138,   144,   138,   138,
     145,   146,   147,   148,   149,   150,   151,   152,   153,     0,
     138,   138,   144,   193,   138,   160,     7,     8,   139,   140,
     141,   142,    83,   144,   193,    52,    17,    18,    19,   139,
     140,   141,   142,    24,   144,    52,   193,    28,   193,    30,
      31,    52,    33,    52,    35,    36,    37,    38,    52,   194,
     139,   140,   141,   142,    52,   144,    52,    48,    52,    52,
      51,    43,   139,   140,   141,   142,    52,   144,   127,   127,
     193,   193,   127,    64,   195,   127,    52,   139,   140,   141,
     142,    52,   144,    52,    52,   195,    77,    78,    79,    52,
      52,    52,    83,   139,   140,   141,   142,    52,   144,   139,
     140,   141,   142,    94,   144,    96,   195,    98,    99,   100,
     101,   102,    52,    52,    52,    52,   107,   108,   195,    52,
      52,   112,   113,   114,   115,   139,   140,   141,   142,    52,
     144,   122,   123,   195,   125,     7,     8,    52,   129,    52,
      52,   132,   133,   134,   135,    17,    18,    19,    52,   195,
      43,   193,    24,    52,   193,   195,    28,   193,    30,    31,
     197,    33,    52,    35,    36,    37,    38,   195,   194,   200,
     161,   138,   139,   140,   141,   142,    48,   144,    83,    51,
     195,   195,   139,   140,   141,   142,    83,   144,   139,   140,
     141,   142,    64,   144,   139,   140,   141,   142,   195,   144,
     195,   193,   138,   138,    83,    77,    78,    79,   138,    83,
     193,    83,   139,   140,   141,   142,    52,   144,   139,   140,
     141,   142,    94,   144,    96,    52,    98,    99,   100,   101,
     102,   195,    34,    34,   193,   107,   108,    83,   195,    83,
     112,   113,   114,   115,   195,    83,   144,   195,   193,     5,
     122,   123,   193,   125,   197,    83,    34,   129,    83,   193,
     132,   133,   134,   135,    20,    21,   193,   195,   138,    25,
     138,   138,   193,   139,   140,   141,   142,    34,   144,    34,
      34,    52,   193,    52,    40,    41,    42,   138,    44,   161,
     200,    47,   139,   140,   141,   142,   138,   144,    54,    55,
      56,   138,   193,    59,    60,   193,    62,    63,   195,    34,
      66,    67,    68,    69,    70,    71,    72,   197,    74,    75,
      76,   138,   138,    52,   193,    52,    82,   193,   193,    85,
      86,   131,    88,   197,    90,   138,   195,    52,    83,    95,
     139,   140,   141,   142,   138,   144,   193,   103,   104,   197,
     106,    52,   193,     9,    10,    11,    12,    13,    14,    15,
      16,   200,   118,   119,    20,    21,   138,    44,   124,    25,
     138,   138,   138,     9,    10,    11,    12,    13,    14,    15,
     136,   137,  1051,   959,   724,    21,    42,   572,    44,    25,
     139,   140,   141,   142,   193,   144,   320,   139,   140,   141,
     142,   722,   144,   306,   693,   312,    42,   139,   140,   141,
     142,   317,   144,   602,   136,   548,   713,  1120,     3,   139,
     140,   141,   142,   365,   144,   211,   638,  1056,   729,    85,
      86,   139,   140,   141,   142,   765,   144,   120,   781,    24,
      -1,    26,    27,   557,   193,    -1,   562,   103,   104,    85,
      86,   193,    -1,    -1,   139,   140,   141,   142,    -1,   144,
      -1,   193,    47,    -1,    49,    50,    -1,   103,   104,    54,
      -1,    -1,    -1,   193,    -1,    -1,    61,    -1,    -1,    -1,
     136,   137,    -1,    -1,    -1,   193,   139,   140,   141,   142,
      -1,   144,   139,   140,   141,   142,    -1,   144,    -1,    -1,
     136,   137,    87,    -1,    89,    90,    91,    92,   193,    -1,
      -1,    -1,    97,    -1,    -1,    -1,   101,    -1,    -1,    -1,
     105,    -1,    -1,    -1,   109,   110,   111,    -1,    -1,    -1,
     115,   116,   117,    -1,   119,    -1,   139,   140,   141,   142,
     193,   144,    -1,    -1,    -1,    -1,   193,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,    -1,    -1,
      -1,   186,   187,    -1,   189,   190,   191,   192,   138,   139,
     140,   141,   142,    -1,   144,    -1,   139,   140,   141,   142,
     193,   144,   139,   140,   141,   142,    -1,   144
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
     133,   134,   135,   161,   204,   205,   206,   207,   208,   209,
     210,   211,   212,   217,   218,   219,   220,   223,   225,   228,
     233,   244,   245,   250,   253,   256,   259,   262,   268,   274,
     277,   282,   285,   288,   291,   292,   295,   296,   298,   299,
     300,   304,   305,   306,   307,   313,   316,   322,   325,   326,
      52,   194,    52,   194,   193,   194,   193,   193,   194,    34,
      43,    52,    83,   194,    83,   194,   193,    83,   193,   194,
     265,   193,   193,   193,   193,   193,   194,    34,    43,   193,
     194,   194,   193,    34,   193,   193,   193,   194,   265,   265,
      83,   216,    34,    52,   314,   194,   194,   265,   193,    34,
     193,   194,   193,   194,   193,   194,   265,   193,   194,   265,
     265,    83,   213,    83,   214,    83,   215,   265,   193,   194,
       0,   205,   193,     9,    10,    11,    12,    13,    14,    15,
      21,    25,    42,    85,    86,   103,   104,   136,   137,   319,
     320,   321,   350,   351,   352,   353,   354,   385,   386,   388,
     389,   392,   393,   394,   395,   396,   397,   398,   193,    16,
      20,    44,   320,   323,   324,   358,   375,   399,    23,     4,
      83,   301,   302,   119,   257,   258,   331,    43,   193,    52,
     193,   193,   201,    83,   193,   201,    83,    83,   226,   227,
      34,     5,    40,    41,    47,    54,    55,    56,    59,    60,
      62,    63,    66,    67,    68,    69,    70,    71,    72,    74,
      75,    76,    82,    88,    90,    95,   106,   118,   124,   283,
     284,   331,   341,   350,   351,   352,   353,   354,   355,   356,
     357,   358,   359,   360,   361,   362,   363,   364,   365,   366,
     367,   368,   369,   370,   371,   372,   373,   374,   375,   376,
     377,   378,   380,   381,   385,   386,   387,   388,   389,   390,
      83,   138,   193,    22,    83,   121,   269,   270,   271,    22,
      83,   121,   278,   279,    22,    83,   121,   275,   276,    83,
     229,   230,   226,    39,   224,    43,   193,   234,    61,   120,
     130,   333,    80,    92,   105,   308,   309,   382,   383,   384,
      22,   132,   246,   247,    29,    43,    52,    65,    73,    83,
     139,   140,   145,   146,   147,   148,   149,   150,   151,   152,
     153,   160,   194,   221,    83,   293,   294,    83,   297,     3,
      24,    26,    27,    49,    50,    87,    89,    91,    97,   101,
     109,   110,   111,   115,   116,   117,   264,   330,   331,   332,
     333,   334,   335,   336,   337,   338,   339,   340,   341,   342,
     343,   344,   345,   347,   348,   349,   357,   379,   383,   384,
     193,   193,   127,    83,   138,   193,    52,   193,    29,    43,
      52,    65,    73,    83,   139,   140,   145,   146,   147,   148,
     149,   150,   151,   152,   153,   160,   194,   241,   243,   286,
     287,   341,   357,   358,   367,   373,   374,   375,   376,   377,
     378,   385,   386,   387,   286,   193,   246,   198,   260,   261,
     344,   350,   254,   255,   331,   263,   264,   193,   126,   264,
     317,   318,   391,   193,   193,   127,    83,   138,   193,   127,
      83,   138,   193,   127,    83,   138,   193,   193,   162,   163,
     164,   165,   166,   167,   168,   169,   170,   171,   172,   173,
     174,   175,   176,   177,   178,   179,   180,   181,   182,   186,
     187,   189,   190,   191,   192,   327,   328,   400,   401,   402,
     403,   404,   405,   406,   407,   408,   409,   410,   411,   412,
     413,   414,   415,   416,   417,   418,   419,   420,   421,   422,
     423,   424,   425,   426,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,   138,   195,    34,
      34,    34,   138,   195,   195,    83,   138,   194,   303,    32,
     302,    34,   138,   195,   193,   193,    83,   195,   201,    83,
     195,   201,    34,    32,   227,    83,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,   138,   195,    34,    83,    83,    83,    32,   270,
     138,    83,   138,    83,    32,   279,    83,   138,    83,    32,
     276,   194,    32,   230,    32,    34,   195,   193,   196,   239,
     240,   241,   242,   138,   195,   195,   195,    34,   138,   195,
      83,    83,    32,   247,   194,   194,   194,   194,   221,   221,
     194,   194,   194,   194,   194,   194,   194,   194,   194,   194,
     221,   139,   140,   141,   142,   144,   193,   194,    32,   294,
     138,   221,    32,    83,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,   195,   127,    83,   193,
     194,   194,   194,   194,   241,   241,   194,   194,   194,   194,
     194,   194,   194,   194,   194,   194,   241,   139,   140,   141,
     142,   144,   315,   138,   195,   195,    32,    43,    52,   194,
     251,   252,   138,   195,   138,   195,   138,   195,    34,   138,
     195,   127,    83,   127,    83,   127,    83,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,    34,    34,    34,    34,    34,    34,
      34,    34,    34,    34,   195,   138,    43,    52,   329,    43,
     139,   140,   267,   329,    52,   267,    52,    83,    52,    52,
     198,   429,   430,    52,    52,    83,    83,   427,   321,    52,
      52,   329,    52,   324,    52,   193,   194,    83,    43,    52,
      34,    52,   258,   193,   193,   193,   265,    83,   193,   193,
     265,    83,   221,   430,    52,    52,    52,    52,    52,   329,
     329,   329,    52,    52,    52,    52,    83,   194,   329,   284,
     193,   265,    83,    34,   138,     6,    43,    45,    52,    53,
      83,    93,   128,   139,   272,   280,   281,   138,   281,   138,
     138,   281,   138,    52,   139,   140,   266,    83,   193,    83,
      32,   240,   242,    34,   193,    46,    58,    64,    84,   231,
     232,   345,   346,   238,   193,   193,    57,    81,   309,    83,
     141,   197,   201,   202,   310,   311,   312,   138,    34,   138,
     193,   221,   221,   221,   221,   222,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   221,   195,   221,   221,   221,
     221,   221,   221,    83,   193,   138,   221,    52,   329,    52,
      52,    52,    52,    52,    52,   329,    52,    52,    52,   193,
     265,   127,   241,   241,   241,   266,   241,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   195,   241,   241,   241,
     241,   241,   193,   287,   193,   265,   193,   265,   221,   193,
     199,    43,    52,   138,   194,   261,   193,   255,   193,   264,
     193,   265,   329,   318,   193,   265,   127,   127,   127,    52,
      52,    52,    52,    52,    52,    52,    52,    52,    52,    52,
      52,    52,    52,    52,   329,   329,    52,   430,   329,   329,
      52,    52,   329,    52,   329,   329,   193,   327,    43,    43,
      52,   428,   199,   428,   197,   193,   193,    52,   303,   195,
     195,   221,   193,   195,   193,   195,   193,   200,   289,   290,
     193,    83,    83,    43,    52,   193,   138,   138,    83,   138,
     281,    83,   193,   281,    52,    52,   195,   226,    34,   241,
      34,   138,   195,   193,   236,   235,   138,   193,   194,   312,
      83,   221,    83,   101,   121,   195,   138,   138,   138,   195,
     195,   195,   195,   195,   195,   195,   195,   195,   195,   195,
     195,   221,    83,   193,   193,   195,   138,   138,   195,   195,
     195,   195,   195,   195,   195,   195,   195,   195,   195,   193,
     193,   195,   252,   193,    43,    52,   194,   221,   193,   193,
     197,    83,   195,    34,   193,   193,   265,   193,   265,    83,
     138,   195,   273,   281,   280,   281,   138,   281,   138,   138,
     193,    34,    32,   241,   193,   329,   232,   237,   239,   239,
     239,   311,   281,    34,   193,    34,    52,   248,   221,   221,
     221,   221,   193,   193,   221,   241,   241,   221,   195,    52,
     303,   221,   193,   193,   200,   289,   138,   138,   138,   281,
     193,   281,   281,   221,   193,   193,    32,    32,    32,   194,
     195,   221,   221,   197,    52,   138,   193,   193,   195,   195,
     193,   195,   195,   195,    34,   193,   138,   281,   273,   281,
     138,   193,   193,   193,   239,   281,   193,   193,    52,   197,
      52,   131,   221,   200,   281,   138,   138,   281,    32,   195,
      52,   197,   221,   249,   193,    83,   281,   280,   193,    52,
     193,   221,   200,   138,   138,   281,   273,   138,   281
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
     445,   446,   447,    59,    40,    41,    35,    58,    91,    93,
      39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   203,   204,   204,   205,   205,   205,   205,   205,   205,
     205,   205,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   205,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   205,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   205,   205,   205,   206,   206,   206,   206,   207,   207,
     208,   209,   210,   211,   212,   213,   213,   213,   213,   213,
     213,   214,   214,   214,   214,   214,   214,   215,   215,   215,
     215,   215,   215,   216,   216,   216,   216,   216,   216,   217,
     217,   218,   218,   219,   219,   220,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   221,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   221,   221,   221,   221,   221,
     221,   222,   222,   223,   223,   224,   225,   226,   226,   227,
     228,   229,   229,   230,   231,   231,   232,   232,   232,   232,
     232,   234,   233,   235,   233,   236,   233,   237,   233,   238,
     233,   239,   239,   239,   239,   240,   240,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   241,
     241,   242,   243,   243,   244,   245,   246,   246,   247,   247,
     247,   247,   247,   248,   248,   248,   248,   248,   248,   249,
     249,   250,   251,   251,   252,   252,   252,   252,   252,   252,
     252,   252,   252,   253,   253,   254,   254,   255,   256,   256,
     257,   257,   258,   259,   259,   260,   260,   261,   261,   262,
     262,   262,   262,   263,   263,   264,   264,   264,   264,   264,
     264,   264,   264,   264,   264,   264,   264,   264,   264,   264,
     264,   264,   264,   264,   264,   264,   264,   264,   265,   265,
     265,   265,   265,   265,   266,   266,   266,   267,   267,   267,
     268,   269,   269,   270,   271,   271,   271,   272,   272,   272,
     272,   272,   273,   273,   273,   273,   274,   275,   275,   276,
     276,   276,   277,   278,   278,   279,   279,   279,   280,   280,
     280,   280,   280,   281,   281,   281,   281,   281,   281,   282,
     282,   282,   282,   283,   283,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   285,   285,   286,   286,   287,   287,
     287,   287,   287,   287,   287,   287,   287,   287,   287,   287,
     287,   288,   288,   289,   289,   290,   290,   291,   292,   293,
     293,   294,   295,   296,   297,   297,   297,   297,   298,   299,
     299,   299,   299,   300,   301,   301,   302,   302,   302,   303,
     303,   303,   304,   304,   305,   305,   305,   305,   305,   305,
     306,   306,   306,   306,   306,   306,   307,   308,   308,   309,
     309,   309,   310,   310,   310,   310,   311,   311,   312,   312,
     312,   312,   312,   314,   315,   313,   316,   316,   316,   316,
     317,   317,   318,   318,   319,   319,   319,   319,   319,   319,
     319,   320,   320,   320,   320,   320,   320,   320,   320,   320,
     320,   321,   321,   322,   322,   323,   323,   323,   323,   324,
     324,   325,   325,   326,   326,   327,   327,   328,   328,   328,
     328,   328,   328,   328,   328,   328,   328,   328,   328,   328,
     328,   328,   328,   328,   328,   328,   328,   328,   328,   328,
     328,   328,   328,   328,   329,   329,   330,   331,   332,   333,
     334,   335,   336,   337,   338,   339,   340,   341,   342,   343,
     344,   345,   346,   347,   348,   349,   350,   351,   351,   352,
     353,   354,   355,   356,   357,   357,   358,   359,   360,   361,
     362,   363,   364,   365,   366,   367,   368,   369,   370,   371,
     372,   373,   374,   375,   376,   377,   378,   379,   380,   381,
     382,   382,   383,   384,   385,   386,   387,   388,   389,   390,
     391,   392,   393,   394,   395,   396,   397,   398,   399,   400,
     401,   402,   403,   404,   405,   406,   407,   408,   409,   410,
     411,   412,   413,   414,   415,   416,   417,   418,   419,   420,
     421,   422,   423,   424,   425,   426,   427,   428,   428,   429,
     429,   430
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
       3,     3,     3,     3,     3,     2,     2,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     6,     6,
       4,     1,     3,     4,     7,     3,     4,     2,     1,     4,
       4,     2,     1,     7,     3,     1,     1,     1,     1,     1,
       1,     0,     5,     0,     8,     0,     8,     0,    10,     0,
       8,     2,     2,     1,     1,     4,     2,     3,     1,     1,
       1,     3,     3,     3,     3,     3,     2,     2,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     6,
       6,     5,     1,     4,     4,     4,     2,     1,     9,     6,
       5,     7,     7,     2,     4,     3,     5,     3,     1,     2,
       1,     6,     3,     1,     5,     3,     3,     4,     2,     2,
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
  "VAROBS", "XLS_SHEET", "XLS_RANGE", "COMMA", "MINUS", "PLUS", "DIVIDE",
  "TIMES", "UMINUS", "POWER", "EXP", "LOG", "LOG10", "SIN", "COS", "TAN",
  "ASIN", "ACOS", "ATAN", "SINH", "COSH", "TANH", "ASINH", "ACOSH",
  "ATANH", "SQRT", "DYNARE_SENSITIVITY", "IDENTIFICATION", "MORRIS",
  "STAB", "REDFORM", "PPRIOR", "PRIOR_RANGE", "PPOST", "ILPTAU", "GLUE",
  "MORRIS_NLIV", "MORRIS_NTRA", "NSAM", "LOAD_REDFORM", "LOAD_RMSE",
  "LOAD_STAB", "ALPHA2_STAB", "KSSTAT", "LOGTRANS_REDFORM",
  "THRESHOLD_REDFORM", "KSSTAT_REDFORM", "ALPHA2_REDFORM", "NAMENDO",
  "NAMLAGENDO", "NAMEXO", "RMSE", "LIK_ONLY", "VAR_RMSE", "PFILT_RMSE",
  "ISTART_RMSE", "ALPHA_RMSE", "ALPHA2_RMSE", "';'", "'('", "')'", "'#'",
  "':'", "'['", "']'", "'''", "'.'", "'\\\\'", "$accept", "statement_list",
  "statement", "declaration", "dsample", "rplot", "var", "varexo",
  "varexo_det", "parameters", "var_list", "varexo_list", "varexo_det_list",
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
       204,     0,    -1,   205,    -1,   204,   205,    -1,   206,    -1,
     217,    -1,   218,    -1,   219,    -1,   233,    -1,   223,    -1,
     225,    -1,   228,    -1,   220,    -1,   244,    -1,   245,    -1,
     250,    -1,   253,    -1,   256,    -1,   259,    -1,   262,    -1,
     282,    -1,   285,    -1,   288,    -1,   268,    -1,   277,    -1,
     274,    -1,   291,    -1,   292,    -1,   295,    -1,   207,    -1,
     208,    -1,   296,    -1,   298,    -1,   299,    -1,   300,    -1,
     304,    -1,   305,    -1,   306,    -1,   307,    -1,   313,    -1,
     316,    -1,   322,    -1,   325,    -1,   326,    -1,   212,    -1,
     209,    -1,   210,    -1,   211,    -1,    28,    52,   193,    -1,
      28,    52,    52,   193,    -1,   112,   265,   193,    -1,   132,
     213,   193,    -1,   133,   214,   193,    -1,   134,   215,   193,
      -1,   100,   216,   193,    -1,   213,    83,    -1,   213,   138,
      83,    -1,    83,    -1,   213,    83,   127,    -1,   213,   138,
      83,   127,    -1,    83,   127,    -1,   214,    83,    -1,   214,
     138,    83,    -1,    83,    -1,   214,    83,   127,    -1,   214,
     138,    83,   127,    -1,    83,   127,    -1,   215,    83,    -1,
     215,   138,    83,    -1,    83,    -1,   215,    83,   127,    -1,
     215,   138,    83,   127,    -1,    83,   127,    -1,   216,    83,
      -1,   216,   138,    83,    -1,    83,    -1,   216,    83,   127,
      -1,   216,   138,    83,   127,    -1,    83,   127,    -1,   101,
      52,   193,    -1,   101,    34,    52,   193,    -1,    24,    43,
     193,    -1,    24,    34,    43,   193,    -1,    64,    43,   193,
      -1,    64,    34,    43,   193,    -1,    83,    34,   221,   193,
      -1,   194,   221,   195,    -1,    83,    -1,    43,    -1,    52,
      -1,   221,   140,   221,    -1,   221,   139,   221,    -1,   221,
     141,   221,    -1,   221,   142,   221,    -1,   221,   144,   221,
      -1,   139,   221,    -1,   140,   221,    -1,   145,   194,   221,
     195,    -1,   146,   194,   221,   195,    -1,   147,   194,   221,
     195,    -1,   148,   194,   221,   195,    -1,   149,   194,   221,
     195,    -1,   150,   194,   221,   195,    -1,   151,   194,   221,
     195,    -1,   152,   194,   221,   195,    -1,   153,   194,   221,
     195,    -1,   160,   194,   221,   195,    -1,    29,   194,   221,
     195,    -1,    65,   194,   221,   138,   221,   195,    -1,    73,
     194,   221,   138,   221,   195,    -1,    83,   194,   222,   195,
      -1,   221,    -1,   222,   138,   221,    -1,    51,   193,   226,
      32,    -1,    51,   194,   224,   195,   193,   226,    32,    -1,
      39,    34,    83,    -1,    33,   193,   226,    32,    -1,   226,
     227,    -1,   227,    -1,    83,    34,   221,   193,    -1,    48,
     193,   229,    32,    -1,   229,   230,    -1,   230,    -1,    83,
     194,   266,   195,    34,   221,   193,    -1,   231,   138,   232,
      -1,   232,    -1,    58,    -1,    46,    -1,    84,    -1,   345,
      -1,   346,    -1,    -1,    77,   193,   234,   239,    32,    -1,
      -1,    77,   194,   333,   195,   193,   235,   239,    32,    -1,
      -1,    77,   194,   130,   195,   193,   236,   239,    32,    -1,
      -1,    77,   194,   120,   138,   231,   195,   237,   193,   239,
      32,    -1,    -1,    77,   194,   120,   195,   238,   193,   239,
      32,    -1,   239,   240,    -1,   239,   242,    -1,   240,    -1,
     242,    -1,   241,    34,   241,   193,    -1,   241,   193,    -1,
     194,   241,   195,    -1,   243,    -1,    43,    -1,    52,    -1,
     241,   140,   241,    -1,   241,   139,   241,    -1,   241,   141,
     241,    -1,   241,   142,   241,    -1,   241,   144,   241,    -1,
     139,   241,    -1,   140,   241,    -1,   145,   194,   241,   195,
      -1,   146,   194,   241,   195,    -1,   147,   194,   241,   195,
      -1,   148,   194,   241,   195,    -1,   149,   194,   241,   195,
      -1,   150,   194,   241,   195,    -1,   151,   194,   241,   195,
      -1,   152,   194,   241,   195,    -1,   153,   194,   241,   195,
      -1,   160,   194,   241,   195,    -1,    29,   194,   241,   195,
      -1,    65,   194,   241,   138,   241,   195,    -1,    73,   194,
     241,   138,   241,   195,    -1,   196,    83,    34,   241,   193,
      -1,    83,    -1,    83,   194,   266,   195,    -1,   113,   193,
     246,    32,    -1,    79,   193,   246,    32,    -1,   246,   247,
      -1,   247,    -1,   132,    83,   193,   101,   248,   193,   131,
     249,   193,    -1,   132,    83,   193,   121,   221,   193,    -1,
     132,    83,    34,   221,   193,    -1,   132,    83,   138,    83,
      34,   221,   193,    -1,    22,    83,   138,    83,    34,   221,
     193,    -1,   248,    52,    -1,   248,    52,   197,    52,    -1,
     248,   138,    52,    -1,   248,   138,    52,   197,    52,    -1,
      52,   197,    52,    -1,    52,    -1,   249,   221,    -1,   221,
      -1,   114,    34,   198,   251,   199,   193,    -1,   251,   193,
     252,    -1,   252,    -1,   252,   138,   194,   221,   195,    -1,
     252,   138,    43,    -1,   252,   138,    52,    -1,   252,   194,
     221,   195,    -1,   252,    43,    -1,   252,    52,    -1,   194,
     221,   195,    -1,    43,    -1,    52,    -1,   122,   193,    -1,
     122,   194,   254,   195,   193,    -1,   254,   138,   255,    -1,
     255,    -1,   331,    -1,    19,   193,    -1,    19,   194,   257,
     195,   193,    -1,   257,   138,   258,    -1,   258,    -1,   331,
      -1,   115,   193,    -1,   115,   194,   260,   195,   193,    -1,
     260,   138,   261,    -1,   261,    -1,   344,    -1,   350,    -1,
     123,   193,    -1,   123,   194,   263,   195,   193,    -1,   123,
     265,   193,    -1,   123,   194,   263,   195,   265,   193,    -1,
     263,   138,   264,    -1,   264,    -1,   330,    -1,   331,    -1,
     332,    -1,   333,    -1,   334,    -1,   335,    -1,   336,    -1,
     337,    -1,   338,    -1,   339,    -1,   340,    -1,   357,    -1,
     341,    -1,   379,    -1,   342,    -1,   343,    -1,   344,    -1,
     345,    -1,   347,    -1,   348,    -1,   349,    -1,   383,    -1,
     384,    -1,   265,    83,    -1,   265,    83,    34,    83,    -1,
     265,   138,    83,    -1,   265,   138,    83,    34,    83,    -1,
      83,    -1,    83,    34,    83,    -1,   140,    52,    -1,   139,
      52,    -1,    52,    -1,   140,    43,    -1,   139,    43,    -1,
      43,    -1,    36,   193,   269,    32,    -1,   269,   270,    -1,
     270,    -1,   271,   138,   272,   193,    -1,   121,    83,    -1,
      83,    -1,    22,    83,   138,    83,    -1,   280,   138,   273,
      -1,   281,   138,   280,   138,   273,    -1,   281,   138,   281,
     138,   281,   138,   280,   138,   273,    -1,   281,    -1,   281,
     138,   281,   138,   281,    -1,   281,   138,   281,    -1,   281,
     138,   281,   138,   281,    -1,   281,   138,   281,   138,   281,
     138,   281,    -1,   281,   138,   281,   138,   281,   138,   281,
     138,   281,    -1,    38,   193,   275,    32,    -1,   275,   276,
      -1,   276,    -1,   121,    83,   138,   281,   193,    -1,    22,
      83,   138,    83,   138,   281,   193,    -1,    83,   138,   281,
     193,    -1,    37,   193,   278,    32,    -1,   278,   279,    -1,
     279,    -1,   121,    83,   138,   281,   138,   281,   193,    -1,
      22,    83,   138,    83,   138,   281,   138,   281,   193,    -1,
      83,   138,   281,   138,   281,   193,    -1,     6,    -1,    45,
      -1,    93,    -1,    53,    -1,   128,    -1,    -1,    52,    -1,
      43,    -1,    83,    -1,   139,    52,    -1,   139,    43,    -1,
      35,   193,    -1,    35,   194,   283,   195,   193,    -1,    35,
     265,   193,    -1,    35,   194,   283,   195,   265,   193,    -1,
     283,   138,   284,    -1,   284,    -1,   350,    -1,   351,    -1,
     352,    -1,   353,    -1,   354,    -1,   355,    -1,   356,    -1,
     357,    -1,   358,    -1,   359,    -1,   360,    -1,   361,    -1,
     362,    -1,   363,    -1,   364,    -1,   365,    -1,   366,    -1,
     367,    -1,   368,    -1,   369,    -1,   370,    -1,   371,    -1,
     372,    -1,   373,    -1,   341,    -1,   374,    -1,   375,    -1,
     376,    -1,   377,    -1,   378,    -1,   380,    -1,   381,    -1,
     385,    -1,   386,    -1,   387,    -1,   331,    -1,   388,    -1,
     389,    -1,   390,    -1,   107,   194,   286,   195,   193,    -1,
     107,   194,   286,   195,   265,   193,    -1,   286,   138,   287,
      -1,   287,    -1,   357,    -1,   358,    -1,   367,    -1,   373,
      -1,   341,    -1,   374,    -1,   375,    -1,   376,    -1,   377,
      -1,   378,    -1,   385,    -1,   386,    -1,   387,    -1,   108,
     194,   286,   195,   193,    -1,   108,   194,   286,   195,   265,
     193,    -1,   200,    83,   200,   138,   200,    83,   200,    -1,
     200,    83,   200,   138,   281,    -1,   289,    -1,   290,   138,
     289,    -1,   135,   265,   193,    -1,    94,   193,   293,    32,
      -1,   293,   294,    -1,   294,    -1,    83,   194,   221,   195,
     193,    -1,   129,   265,   193,    -1,    96,   193,   297,    32,
      -1,   297,    83,   221,   193,    -1,   297,    83,   138,    83,
     221,   193,    -1,    83,   221,   193,    -1,    83,   138,    83,
     221,   193,    -1,    99,   265,   193,    -1,    98,   193,    -1,
      98,   194,   264,   195,   193,    -1,    98,   265,   193,    -1,
      98,   194,   264,   195,   265,   193,    -1,    18,   193,   301,
      32,    -1,   301,   302,    -1,   302,    -1,    83,   303,    34,
     221,   193,    -1,    83,   138,    83,   303,    34,   221,   193,
      -1,     4,    83,   194,    52,   195,   303,    34,   221,   193,
      -1,    -1,   194,    52,   195,    -1,   194,    43,   195,    -1,
      17,   193,    -1,    17,   194,    23,   195,   193,    -1,    31,
     194,    83,   195,   193,    -1,    31,   194,    83,   195,   265,
     193,    -1,    31,    83,   193,    -1,    31,   194,    83,   201,
      83,   195,   193,    -1,    31,   194,    83,   201,    83,   195,
     265,   193,    -1,    31,    83,   201,    83,   193,    -1,    30,
     194,    83,   195,   193,    -1,    30,   194,    83,   195,   265,
     193,    -1,    30,    83,   193,    -1,    30,   194,    83,   201,
      83,   195,   193,    -1,    30,   194,    83,   201,    83,   195,
     265,   193,    -1,    30,    83,   201,    83,   193,    -1,    78,
     194,   308,   195,   310,   193,    -1,   308,   138,   309,    -1,
     309,    -1,   382,    -1,   383,    -1,   384,    -1,   311,    -1,
     310,   138,   311,    -1,   311,   194,   281,   195,    -1,   310,
     138,   311,   194,   281,   195,    -1,   312,    -1,   311,   312,
      -1,    83,    -1,   202,    -1,   141,    -1,   197,    -1,   201,
      -1,    -1,    -1,   102,   314,   241,   315,   193,    -1,   125,
     193,    -1,   125,   194,   317,   195,   193,    -1,   125,   265,
     193,    -1,   125,   194,   317,   195,   265,   193,    -1,   317,
     138,   318,    -1,   318,    -1,   264,    -1,   391,    -1,   392,
      -1,   393,    -1,   394,    -1,   395,    -1,   396,    -1,   397,
      -1,   398,    -1,   319,    -1,   350,    -1,   385,    -1,   386,
      -1,   352,    -1,   354,    -1,   351,    -1,   353,    -1,   388,
      -1,   389,    -1,   320,   138,   321,    -1,   320,    -1,     7,
      52,   193,    -1,     7,   194,   321,   195,    52,   193,    -1,
     320,    -1,   375,    -1,   358,    -1,   399,    -1,   323,   138,
     324,    -1,   323,    -1,     8,    52,   193,    -1,     8,   194,
     324,   195,    52,   193,    -1,   161,   193,    -1,   161,   194,
     327,   195,   193,    -1,   328,   138,   327,    -1,   328,    -1,
     400,    -1,   401,    -1,   402,    -1,   403,    -1,   404,    -1,
     405,    -1,   406,    -1,   407,    -1,   408,    -1,   409,    -1,
     410,    -1,   411,    -1,   412,    -1,   413,    -1,   414,    -1,
     415,    -1,   416,    -1,   417,    -1,   419,    -1,   420,    -1,
     421,    -1,   422,    -1,   423,    -1,   424,    -1,   425,    -1,
     426,    -1,   418,    -1,    52,    -1,    43,    -1,    26,    34,
      52,    -1,   119,    34,    52,    -1,   116,    34,    52,    -1,
      61,    -1,    97,    34,    52,    -1,   111,    34,    52,    -1,
      27,    34,    52,    -1,     3,    34,    52,    -1,    87,    -1,
      89,    -1,    91,    -1,    54,    34,    52,    -1,    49,    34,
      52,    -1,    50,    34,    52,    -1,   101,    34,    52,    -1,
      24,    34,   329,    -1,    64,    34,   329,    -1,   115,    -1,
     117,    34,    52,    -1,   109,    34,   329,    -1,    25,    34,
      83,    -1,    85,    34,   430,    -1,    85,    34,    52,    -1,
      42,    34,    52,    -1,   103,    34,    52,    -1,   104,    34,
      52,    -1,    59,    34,    52,    -1,    60,    34,    52,    -1,
      90,    -1,    47,    -1,    20,    34,   329,    -1,    71,    34,
      52,    -1,    66,    34,   329,    -1,    68,    34,   329,    -1,
      95,    34,   194,   290,   195,    -1,    67,    34,   329,    -1,
      76,    34,    83,    -1,    75,    34,    52,    -1,    74,    -1,
     106,    34,   329,    -1,    69,    34,    52,    -1,    70,    34,
      52,    -1,    62,    -1,    63,    -1,    88,    -1,     5,    -1,
     124,    -1,    44,    34,    52,    -1,   118,    -1,    82,    -1,
      41,    -1,   110,    -1,    55,    34,    52,    -1,    56,    34,
      52,    -1,    80,    34,    57,    -1,    80,    34,    81,    -1,
     105,    -1,    92,    -1,   136,    34,    83,    -1,   137,    34,
     427,    -1,    40,    34,   430,    -1,    21,    -1,    86,    -1,
      72,    -1,   126,    34,   329,    -1,    14,    34,   267,    -1,
       9,    34,   329,    -1,    11,    34,   267,    -1,    12,    34,
     329,    -1,    13,    34,    52,    -1,    10,    -1,    15,    34,
      52,    -1,    16,    34,    52,    -1,   162,    34,    52,    -1,
     163,    34,    52,    -1,   164,    34,    52,    -1,   165,    34,
      52,    -1,   166,    34,    52,    -1,   167,    34,    52,    -1,
     168,    34,    52,    -1,   169,    34,    52,    -1,   170,    34,
      52,    -1,   171,    34,    52,    -1,   172,    34,    52,    -1,
     173,    34,    52,    -1,   174,    34,    52,    -1,   175,    34,
      52,    -1,   176,    34,    52,    -1,   177,    34,   329,    -1,
     178,    34,   329,    -1,   179,    34,    52,    -1,   180,    34,
     430,    -1,   181,    34,   329,    -1,   182,    34,   329,    -1,
     186,    34,    52,    -1,   187,    34,    52,    -1,   189,    34,
     329,    -1,   190,    34,    52,    -1,   191,    34,   329,    -1,
     192,    34,   329,    -1,    83,   197,    83,    -1,    52,    -1,
      52,   197,    52,    -1,   198,   428,    -1,   429,   428,    -1,
     429,   199,    -1
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
     251,   255,   259,   263,   267,   271,   274,   277,   282,   287,
     292,   297,   302,   307,   312,   317,   322,   327,   332,   339,
     346,   351,   353,   357,   362,   370,   374,   379,   382,   384,
     389,   394,   397,   399,   407,   411,   413,   415,   417,   419,
     421,   423,   424,   430,   431,   440,   441,   450,   451,   462,
     463,   472,   475,   478,   480,   482,   487,   490,   494,   496,
     498,   500,   504,   508,   512,   516,   520,   523,   526,   531,
     536,   541,   546,   551,   556,   561,   566,   571,   576,   581,
     588,   595,   601,   603,   608,   613,   618,   621,   623,   633,
     640,   646,   654,   662,   665,   670,   674,   680,   684,   686,
     689,   691,   698,   702,   704,   710,   714,   718,   723,   726,
     729,   733,   735,   737,   740,   746,   750,   752,   754,   757,
     763,   767,   769,   771,   774,   780,   784,   786,   788,   790,
     793,   799,   803,   810,   814,   816,   818,   820,   822,   824,
     826,   828,   830,   832,   834,   836,   838,   840,   842,   844,
     846,   848,   850,   852,   854,   856,   858,   860,   862,   865,
     870,   874,   880,   882,   886,   889,   892,   894,   897,   900,
     902,   907,   910,   912,   917,   920,   922,   927,   931,   937,
     947,   949,   955,   959,   965,   973,   983,   988,   991,   993,
     999,  1007,  1012,  1017,  1020,  1022,  1030,  1040,  1047,  1049,
    1051,  1053,  1055,  1057,  1058,  1060,  1062,  1064,  1067,  1070,
    1073,  1079,  1083,  1090,  1094,  1096,  1098,  1100,  1102,  1104,
    1106,  1108,  1110,  1112,  1114,  1116,  1118,  1120,  1122,  1124,
    1126,  1128,  1130,  1132,  1134,  1136,  1138,  1140,  1142,  1144,
    1146,  1148,  1150,  1152,  1154,  1156,  1158,  1160,  1162,  1164,
    1166,  1168,  1170,  1172,  1174,  1180,  1187,  1191,  1193,  1195,
    1197,  1199,  1201,  1203,  1205,  1207,  1209,  1211,  1213,  1215,
    1217,  1219,  1225,  1232,  1240,  1246,  1248,  1252,  1256,  1261,
    1264,  1266,  1272,  1276,  1281,  1286,  1293,  1297,  1303,  1307,
    1310,  1316,  1320,  1327,  1332,  1335,  1337,  1343,  1351,  1361,
    1362,  1366,  1370,  1373,  1379,  1385,  1392,  1396,  1404,  1413,
    1419,  1425,  1432,  1436,  1444,  1453,  1459,  1466,  1470,  1472,
    1474,  1476,  1478,  1480,  1484,  1489,  1496,  1498,  1501,  1503,
    1505,  1507,  1509,  1511,  1512,  1513,  1519,  1522,  1528,  1532,
    1539,  1543,  1545,  1547,  1549,  1551,  1553,  1555,  1557,  1559,
    1561,  1563,  1565,  1567,  1569,  1571,  1573,  1575,  1577,  1579,
    1581,  1583,  1587,  1589,  1593,  1600,  1602,  1604,  1606,  1608,
    1612,  1614,  1618,  1625,  1628,  1634,  1638,  1640,  1642,  1644,
    1646,  1648,  1650,  1652,  1654,  1656,  1658,  1660,  1662,  1664,
    1666,  1668,  1670,  1672,  1674,  1676,  1678,  1680,  1682,  1684,
    1686,  1688,  1690,  1692,  1694,  1696,  1698,  1702,  1706,  1710,
    1712,  1716,  1720,  1724,  1728,  1730,  1732,  1734,  1738,  1742,
    1746,  1750,  1754,  1758,  1760,  1764,  1768,  1772,  1776,  1780,
    1784,  1788,  1792,  1796,  1800,  1802,  1804,  1808,  1812,  1816,
    1820,  1826,  1830,  1834,  1838,  1840,  1844,  1848,  1852,  1854,
    1856,  1858,  1860,  1862,  1866,  1868,  1870,  1872,  1874,  1878,
    1882,  1886,  1890,  1892,  1894,  1898,  1902,  1906,  1908,  1910,
    1912,  1916,  1920,  1924,  1928,  1932,  1936,  1938,  1942,  1946,
    1950,  1954,  1958,  1962,  1966,  1970,  1974,  1978,  1982,  1986,
    1990,  1994,  1998,  2002,  2006,  2010,  2014,  2018,  2022,  2026,
    2030,  2034,  2038,  2042,  2046,  2050,  2054,  2058,  2060,  2064,
    2067,  2070
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    96,    96,    97,   100,   101,   102,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   122,   123,   124,   125,
     126,   127,   128,   129,   130,   131,   132,   133,   134,   135,
     136,   137,   138,   139,   142,   143,   144,   145,   149,   151,
     155,   157,   159,   161,   163,   165,   167,   169,   171,   173,
     175,   179,   181,   183,   185,   187,   189,   193,   195,   197,
     199,   201,   203,   207,   209,   211,   213,   215,   217,   221,
     223,   227,   229,   233,   235,   240,   242,   244,   246,   248,
     250,   252,   254,   256,   258,   260,   262,   264,   266,   268,
     270,   272,   274,   276,   278,   280,   282,   284,   286,   288,
     290,   294,   296,   300,   302,   306,   308,   310,   311,   314,
     316,   318,   319,   322,   324,   325,   328,   330,   332,   334,
     335,   338,   338,   340,   340,   342,   342,   345,   344,   347,
     347,   351,   352,   353,   354,   357,   359,   363,   365,   366,
     368,   370,   372,   374,   376,   378,   380,   382,   384,   386,
     388,   390,   392,   394,   396,   398,   400,   402,   404,   406,
     408,   412,   415,   417,   421,   423,   425,   426,   429,   431,
     433,   435,   437,   441,   443,   445,   447,   449,   451,   455,
     457,   461,   463,   465,   469,   471,   473,   475,   477,   479,
     481,   483,   485,   489,   491,   495,   496,   499,   501,   503,
     507,   508,   511,   513,   515,   519,   520,   523,   524,   527,
     529,   531,   533,   537,   538,   541,   542,   543,   544,   545,
     546,   547,   548,   549,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   562,   563,   566,   568,
     570,   572,   574,   576,   580,   582,   584,   588,   590,   592,
     596,   598,   600,   604,   606,   612,   618,   628,   633,   640,
     651,   656,   667,   674,   683,   694,   709,   712,   714,   718,
     726,   736,   746,   749,   751,   755,   765,   777,   789,   791,
     793,   795,   797,   801,   802,   803,   804,   805,   807,   811,
     813,   815,   817,   821,   822,   825,   826,   827,   828,   829,
     830,   831,   832,   833,   834,   835,   836,   837,   838,   839,
     840,   841,   842,   843,   844,   845,   846,   847,   848,   849,
     850,   851,   852,   853,   854,   855,   856,   857,   858,   859,
     860,   861,   862,   863,   866,   868,   872,   873,   876,   877,
     878,   879,   880,   881,   882,   883,   884,   885,   886,   887,
     888,   891,   893,   897,   899,   903,   904,   907,   909,   911,
     912,   915,   917,   919,   921,   923,   925,   927,   931,   933,
     935,   937,   939,   943,   945,   946,   949,   951,   953,   957,
     958,   960,   964,   966,   970,   972,   974,   976,   978,   980,
     984,   986,   988,   990,   992,   994,   998,  1001,  1002,  1005,
    1006,  1007,  1010,  1012,  1014,  1016,  1020,  1022,  1026,  1027,
    1029,  1031,  1033,  1037,  1038,  1037,  1040,  1042,  1044,  1046,
    1050,  1051,  1054,  1055,  1058,  1059,  1060,  1061,  1062,  1063,
    1064,  1067,  1068,  1069,  1070,  1071,  1072,  1073,  1074,  1075,
    1076,  1079,  1080,  1083,  1085,  1089,  1090,  1091,  1092,  1095,
    1096,  1099,  1101,  1105,  1107,  1111,  1112,  1115,  1116,  1117,
    1118,  1119,  1120,  1121,  1122,  1123,  1124,  1125,  1126,  1127,
    1128,  1129,  1130,  1131,  1132,  1133,  1134,  1135,  1136,  1137,
    1138,  1139,  1140,  1141,  1144,  1145,  1148,  1149,  1150,  1151,
    1152,  1153,  1154,  1155,  1156,  1157,  1158,  1159,  1160,  1161,
    1162,  1164,  1165,  1166,  1167,  1168,  1169,  1170,  1172,  1175,
    1176,  1177,  1178,  1179,  1180,  1182,  1185,  1186,  1187,  1188,
    1189,  1190,  1191,  1192,  1193,  1194,  1195,  1196,  1197,  1198,
    1199,  1200,  1201,  1202,  1203,  1204,  1205,  1206,  1207,  1208,
    1209,  1211,  1214,  1215,  1216,  1217,  1218,  1219,  1220,  1221,
    1222,  1224,  1225,  1226,  1227,  1228,  1229,  1230,  1231,  1233,
    1234,  1235,  1236,  1237,  1238,  1239,  1240,  1241,  1242,  1243,
    1244,  1245,  1246,  1247,  1248,  1249,  1250,  1251,  1253,  1254,
    1260,  1261,  1265,  1266,  1267,  1268,  1271,  1279,  1280,  1289,
    1291,  1300
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
       2,     2,     2,     2,     2,   196,     2,     2,     2,   200,
     194,   195,     2,     2,     2,     2,   201,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   197,   193,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   198,   202,   199,     2,     2,     2,     2,     2,     2,
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
     185,   186,   187,   188,   189,   190,   191,   192
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1777;
  const int parser::yynnts_ = 228;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 160;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 203;

  const unsigned int parser::yyuser_token_number_max_ = 447;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1302 "DynareBison.yy"


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

