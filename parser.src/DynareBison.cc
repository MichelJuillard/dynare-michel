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
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 49:
#line 151 "DynareBison.yy"
    {driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 50:
#line 154 "DynareBison.yy"
    {driver.rplot();;}
    break;

  case 55:
#line 174 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 56:
#line 176 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 57:
#line 178 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 58:
#line 180 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 59:
#line 182 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 60:
#line 184 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 61:
#line 189 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 62:
#line 191 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 63:
#line 193 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 64:
#line 195 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 65:
#line 197 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 66:
#line 199 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 67:
#line 204 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 68:
#line 206 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 69:
#line 208 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 70:
#line 210 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 71:
#line 212 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 72:
#line 214 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 73:
#line 219 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 74:
#line 221 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 75:
#line 223 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 76:
#line 225 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 77:
#line 227 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 78:
#line 229 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 79:
#line 234 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 80:
#line 238 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 81:
#line 245 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 82:
#line 249 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 83:
#line 256 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 84:
#line 260 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 85:
#line 268 "DynareBison.yy"
    {driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 86:
#line 273 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 87:
#line 275 "DynareBison.yy"
    {(yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 88:
#line 277 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 89:
#line 279 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 90:
#line 281 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 91:
#line 283 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 92:
#line 285 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 93:
#line 287 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 94:
#line 289 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 95:
#line 291 "DynareBison.yy"
    {(yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 96:
#line 293 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 97:
#line 295 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 98:
#line 297 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 99:
#line 299 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 100:
#line 301 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 101:
#line 303 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 102:
#line 305 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 103:
#line 307 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 104:
#line 309 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 105:
#line 311 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 106:
#line 313 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 107:
#line 315 "DynareBison.yy"
    {(yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 108:
#line 320 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 109:
#line 322 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 110:
#line 326 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 111:
#line 328 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 112:
#line 332 "DynareBison.yy"
    {driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 113:
#line 337 "DynareBison.yy"
    {driver.end_endval();;}
    break;

  case 116:
#line 347 "DynareBison.yy"
    {driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 117:
#line 352 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 120:
#line 362 "DynareBison.yy"
    {driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 123:
#line 370 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 124:
#line 371 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 127:
#line 377 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 128:
#line 377 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 129:
#line 378 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 130:
#line 379 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 131:
#line 380 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 132:
#line 381 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 133:
#line 382 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 134:
#line 383 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 135:
#line 384 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 136:
#line 385 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 141:
#line 397 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 142:
#line 399 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val));;}
    break;

  case 143:
#line 403 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 145:
#line 406 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 146:
#line 408 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 147:
#line 410 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 148:
#line 412 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 149:
#line 414 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 150:
#line 416 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 151:
#line 418 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 152:
#line 420 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 153:
#line 422 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 154:
#line 424 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 155:
#line 426 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 156:
#line 428 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 157:
#line 430 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 158:
#line 432 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 159:
#line 434 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 160:
#line 436 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 161:
#line 438 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 162:
#line 440 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 163:
#line 442 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 164:
#line 446 "DynareBison.yy"
    {driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 165:
#line 450 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 166:
#line 452 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 167:
#line 456 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 168:
#line 460 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 171:
#line 470 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 172:
#line 472 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 173:
#line 474 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 174:
#line 476 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 175:
#line 478 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 176:
#line 483 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 177:
#line 485 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 178:
#line 487 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 179:
#line 489 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 180:
#line 491 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 181:
#line 493 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 182:
#line 499 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 183:
#line 501 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].node_val));;}
    break;

  case 184:
#line 506 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 185:
#line 511 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 186:
#line 513 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 187:
#line 518 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 188:
#line 520 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 189:
#line 522 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 190:
#line 524 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 191:
#line 526 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 192:
#line 528 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 193:
#line 530 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 194:
#line 532 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 195:
#line 534 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 196:
#line 539 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 197:
#line 543 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 201:
#line 555 "DynareBison.yy"
    {driver.check();;}
    break;

  case 202:
#line 557 "DynareBison.yy"
    {driver.check();;}
    break;

  case 206:
#line 569 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 207:
#line 571 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 212:
#line 584 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 213:
#line 586 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 214:
#line 588 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 215:
#line 590 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 241:
#line 624 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 242:
#line 626 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 243:
#line 628 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 244:
#line 630 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 245:
#line 632 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 246:
#line 634 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 247:
#line 639 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 248:
#line 641 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 249:
#line 643 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 250:
#line 648 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 251:
#line 650 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 252:
#line 652 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 253:
#line 657 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 254:
#line 662 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 255:
#line 664 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 257:
#line 673 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 258:
#line 678 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 259:
#line 683 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 260:
#line 693 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 261:
#line 698 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 262:
#line 704 "DynareBison.yy"
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

  case 263:
#line 714 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 264:
#line 719 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 265:
#line 730 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 266:
#line 736 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 267:
#line 744 "DynareBison.yy"
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

  case 268:
#line 754 "DynareBison.yy"
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

  case 269:
#line 768 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 270:
#line 772 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 271:
#line 774 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 272:
#line 778 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 273:
#line 785 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 274:
#line 794 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 275:
#line 803 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 276:
#line 807 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 277:
#line 809 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 278:
#line 813 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 279:
#line 822 "DynareBison.yy"
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

  case 280:
#line 833 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 281:
#line 845 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 282:
#line 847 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 283:
#line 849 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 284:
#line 851 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 285:
#line 853 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 286:
#line 857 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 290:
#line 862 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 291:
#line 864 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 292:
#line 871 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 293:
#line 873 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 294:
#line 875 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 295:
#line 877 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 337:
#line 928 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 338:
#line 930 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 354:
#line 956 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 355:
#line 958 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 356:
#line 962 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 357:
#line 963 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 360:
#line 973 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 361:
#line 978 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 364:
#line 987 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val));;}
    break;

  case 365:
#line 990 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 366:
#line 994 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 367:
#line 998 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 368:
#line 1000 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 369:
#line 1002 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 370:
#line 1004 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 371:
#line 1007 "DynareBison.yy"
    {driver.set_osr_params();;}
    break;

  case 372:
#line 1010 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 373:
#line 1011 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 374:
#line 1012 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 375:
#line 1013 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 376:
#line 1017 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 379:
#line 1024 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 380:
#line 1025 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 381:
#line 1026 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val));;}
    break;

  case 382:
#line 1029 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 383:
#line 1030 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 384:
#line 1031 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 385:
#line 1034 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 386:
#line 1035 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 387:
#line 1038 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 388:
#line 1039 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 389:
#line 1040 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 390:
#line 1041 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 391:
#line 1042 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 392:
#line 1043 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 393:
#line 1045 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 394:
#line 1046 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 395:
#line 1047 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 396:
#line 1048 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 397:
#line 1049 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 398:
#line 1050 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 399:
#line 1053 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 405:
#line 1065 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 406:
#line 1066 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 407:
#line 1067 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 408:
#line 1068 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 409:
#line 1071 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 410:
#line 1072 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 412:
#line 1076 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 413:
#line 1077 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 414:
#line 1078 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 415:
#line 1079 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 416:
#line 1082 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 417:
#line 1082 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 419:
#line 1086 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 420:
#line 1088 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 421:
#line 1090 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 422:
#line 1092 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 446:
#line 1130 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 447:
#line 1132 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 454:
#line 1146 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 455:
#line 1148 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 456:
#line 1152 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 457:
#line 1154 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 489:
#line 1192 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 490:
#line 1193 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 491:
#line 1194 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 492:
#line 1195 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 493:
#line 1196 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 494:
#line 1197 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1198 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 496:
#line 1199 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 497:
#line 1200 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 498:
#line 1201 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 499:
#line 1202 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 500:
#line 1203 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 501:
#line 1204 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 502:
#line 1205 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 503:
#line 1206 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 504:
#line 1207 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 505:
#line 1208 "DynareBison.yy"
    {driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 506:
#line 1209 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 507:
#line 1210 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 508:
#line 1211 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 509:
#line 1212 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 510:
#line 1213 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 511:
#line 1214 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 512:
#line 1216 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 513:
#line 1217 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 514:
#line 1218 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 515:
#line 1219 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 516:
#line 1220 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 517:
#line 1221 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 518:
#line 1222 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 519:
#line 1223 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 520:
#line 1224 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 521:
#line 1225 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 522:
#line 1226 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 524:
#line 1228 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 525:
#line 1229 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 526:
#line 1230 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 527:
#line 1231 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 528:
#line 1232 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 529:
#line 1233 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 530:
#line 1234 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 531:
#line 1235 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 532:
#line 1236 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 533:
#line 1237 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 534:
#line 1238 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 535:
#line 1239 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 536:
#line 1240 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 537:
#line 1241 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 538:
#line 1242 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 539:
#line 1243 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 540:
#line 1244 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 541:
#line 1245 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 542:
#line 1246 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 543:
#line 1249 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 544:
#line 1251 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 545:
#line 1253 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 546:
#line 1254 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 547:
#line 1255 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 548:
#line 1256 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 549:
#line 1257 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 550:
#line 1258 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 551:
#line 1259 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 552:
#line 1260 "DynareBison.yy"
    {driver.option_num("mh_recover", "1");;}
    break;

  case 553:
#line 1261 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 554:
#line 1263 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 555:
#line 1264 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 556:
#line 1265 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 557:
#line 1266 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 558:
#line 1267 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 559:
#line 1268 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 560:
#line 1269 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 561:
#line 1270 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 562:
#line 1272 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 563:
#line 1273 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 564:
#line 1274 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1275 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 566:
#line 1276 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1277 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 568:
#line 1278 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 569:
#line 1279 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 570:
#line 1280 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 571:
#line 1281 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1282 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1283 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1284 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 575:
#line 1285 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 576:
#line 1286 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 577:
#line 1287 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1288 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 579:
#line 1289 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1290 "DynareBison.yy"
    {driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 581:
#line 1292 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 582:
#line 1293 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1299 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1300 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1304 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1305 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1306 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1307 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1311 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 591:
#line 1320 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 592:
#line 1323 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 593:
#line 1325 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 594:
#line 1333 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2317 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1049;
  const short int
  parser::yypact_[] =
  {
       838,    21,    22,    88,   -78,   306,    74,    50,    -3,     2,
     -71,    49,   -53,   -46,   106,   116,   354,    97,   421,   -84,
     119,    90,   206,   231,    51,    58,    61,    81, -1049,   165,
     216,    58,   277,   425,   472,   499,    60,    64,    58,   385,
     389,   401,    58,   511,   722, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
     296,  1274,   301,  1239, -1049,   476,    94, -1049,   426,   450,
     336,    34,    82,   501,   120,   507,   509,   559, -1049,  1143,
      29,   322,   357,   377,   519,   509,   563,   562,   419, -1049,
      26,   529,   139,   982,   532,   536, -1049,  1366,   178,   181,
     506,   182,   583,   466,   998,   585,   585,   189,   139,   445,
   -1049,   131, -1049,   426, -1049,  1366,   214, -1049,  1298,   234,
     238,   540,   246,   551,   247,   553,   266,   267, -1049,  1426,
   -1049, -1049, -1049,   624, -1049,   651,   653,   660,   664,   669,
   -1049,   670,   671,   677, -1049,   679,   680,   681,   682, -1049,
     582,   543, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,   687,
     694,   698, -1049,   603,   547, -1049, -1049, -1049,   552,   662,
     -45,    79, -1049,   711,   -14, -1049, -1049,   556, -1049,   566,
   -1049, -1049,   668,   251, -1049,   683,   263,   716,    68, -1049,
     684, -1049,   720, -1049, -1049,   727,   728,   732,   733,   734,
   -1049, -1049,   737,   738,   740,   741,   743,   744, -1049, -1049,
     745,   747, -1049, -1049, -1049,   758,   761, -1049, -1049,   235,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
     762,   719, -1049,   721, -1049,   724,   399, -1049,   666,   726,
     688,   735,   418, -1049,   748,   689,   749,   560, -1049,   631,
      69, -1049,    71,   791,   636,   645, -1049,   500, -1049,   300,
     644,   647,   803, -1049, -1049,   333, -1049, -1049, -1049, -1049,
     759,   764,   410, -1049, -1049, -1049,   659,   982,   982,   686,
     690,   691,   692,   693,   696,   697,   705,   706,   708,   982,
     973,   709,   202, -1049,   842,   244,   809,   821,   836,   845,
     873,   875, -1049, -1049, -1049,   877,   878,   882, -1049,   884,
   -1049,   886,   887,   730, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
     800,   844, -1049,   742, -1049, -1049, -1049,   736,   998,   998,
     739,   750,   753,   760,   763,   765,   770,   771,   772,   774,
     998,   652, -1049,   348, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,   384, -1049,
     446,    27,   435, -1049, -1049, -1049,   459, -1049, -1049,   460,
   -1049, -1049,   894, -1049,   461, -1049, -1049, -1049, -1049, -1049,
     814,   859, -1049, -1049,   828,   864, -1049, -1049,   829,   885,
   -1049, -1049,   941,   946,   948,   963,   975,   978,   979,   980,
     992,   994,   995,   996,  1005,  1012,  1014,  1021,  1027,  1036,
    1038,  1040,  1041,  1042,  1043,  1044,  1050,  1052,  1053,   806,
     897, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,   271,   308,
     271,   820,   308,  1037,  1007,  1054,    14,  1056,  1063,  1009,
    1011,  1274,  1064,  1065,   271,  1068,  1239,  1069,   883,   902,
    1023,   420,  1089, -1049, -1049,  1081,   426,   947, -1049, -1049,
     961,    33,  1057,   962,    35,  1072,   982, -1049, -1049, -1049,
     959,  1105,  1106,  1107,  1108,  1109,   271,   271,   271,  1110,
    1114,  1115,  1116,  1090,   981,   271,  1143,    42,  1093,  1136,
    1046, -1049, -1049, -1049,   534,  1047,    70,  1051, -1049, -1049,
    1060,    70,  1066, -1049, -1049,    19, -1049, -1049, -1049,  1098,
    1010, -1049,  1123,    37, -1049,    31, -1049,   564, -1049,  1017,
    1028,    47,   529,    39,  1092,    41, -1049, -1049,   982,  1088,
     535,   982,   982,   982,   982,   982,   982,   982,   982,   982,
     982,   166,   982,   982,   982,   982,   982, -1049,   982, -1049,
   -1049,  1151,  1205, -1049,   858,  1181,   271,  1182,  1188,  1190,
    1193,  1196,  1210,   271,  1211,  1214,  1215,    57, -1049,  1133,
   -1049,    19,  1131,   558,   998,   998,   998,   998,   998,   998,
     998,   998,   998,   998,   530,   998,   998,   998,   998,   998,
    1084,   585,    66,    77, -1049, -1049, -1049,   982,   275,    24,
     131,  1101,   426,  1103,  1366,    85,   271,  1298,   110, -1049,
    1154, -1049,  1155, -1049,  1156,  1230,  1242,  1243,  1245,  1246,
    1247,  1249,  1251,  1256,  1262,  1265,  1266,  1267,  1268,  1279,
     271,   271,  1281,   959,   271,   271,  1282,  1283,   271,  1284,
     271,   271,  1147,  1426, -1049, -1049, -1049, -1049,  1145,  1306,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,  1286,    20,
   -1049, -1049, -1049, -1049,  1157, -1049, -1049,  1160, -1049, -1049,
   -1049, -1049,  1163, -1049,  1302,  1164,  1166,  1169,   982, -1049,
   -1049, -1049, -1049, -1049,   274,  1170, -1049, -1049,   294,  1171,
    1240, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049,  1167, -1049, -1049, -1049,   295,
   -1049,  1288,  1290, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049,   438,  1175,  1232,  1233,  1299,  1253,    70,  1308,  1200,
      70, -1049,  1345,  1346,  1207, -1049,   509,  1367, -1049, -1049,
   -1049,   998, -1049, -1049, -1049,  1368,   463, -1049, -1049, -1049,
    1213, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049,   100,    54, -1049,  1328,   982,  1336,   179,   723,   468,
     545,   571,   588,   646,   673,   754,   767,   835,   899,   906,
   -1049,   535,   535,  1088,  1088,  1277,   915,   982, -1049,  1338,
    1287, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049,   311, -1049,  1237,   921,   930,   944,
     958,   964,  1039,  1055,  1083,  1100,  1132, -1049,   558,   558,
    1131,  1131,  1277, -1049, -1049, -1049,   312, -1049,   313,  1168,
      27,  1257, -1049, -1049,    45,   982, -1049, -1049, -1049, -1049,
   -1049, -1049,   316, -1049, -1049, -1049,   318, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049,  1228, -1049, -1049, -1049,  1354, -1049, -1049,  1259,  1407,
   -1049, -1049,  1295, -1049,   130, -1049,   230, -1049,  1381, -1049,
     469, -1049, -1049, -1049, -1049, -1049, -1049,    70,   534,  1314,
      70,  1330,  1331, -1049,  1273, -1049, -1049,  1433,   276,   998,
    1301,   271,   564, -1049,   500,   500,   500,    39, -1049,    70,
   -1049,  1435,  1307,  1436,  1419,   982,   982, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,  1285,  1320,
     982, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049,    24, -1049, -1049,
   -1049,   982,  1174, -1049, -1049,  1424, -1049,  1164,   982, -1049,
   -1049,   342, -1049,   344,  1289,  1167, -1049, -1049,  1348,  1349,
    1387,    70,  1291,    70,    70, -1049,   982, -1049,  1351, -1049,
   -1049, -1049,  1313,    53,   217,   241,    67,  1324,   982, -1049,
     982,  1311,    30,  1357,   723, -1049, -1049,  1363,  1191, -1049,
   -1049,  1494,  1370, -1049, -1049,  1399, -1049,    70,    70,    70,
    1401, -1049,  1347,  1350,  1376, -1049,   500, -1049, -1049, -1049,
      70, -1049,  1382,  1388,  1486,  1352,  1487,  1414, -1049, -1049,
   -1049,   982, -1049,    25,  1408, -1049,  1409,    70, -1049, -1049,
   -1049,   269,  1353, -1049, -1049, -1049,  1496,  1355,   982,  1394,
    1469, -1049,    70,   118,  1361, -1049, -1049, -1049,  1500,   723,
     874, -1049,  1358,  1421,  1422, -1049, -1049, -1049,   723, -1049,
      70,    70,  1423, -1049,    70, -1049
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   416,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     4,    29,    30,    45,
      46,    47,    44,     5,     6,     7,    12,     9,    10,    11,
       8,    13,    14,    15,    16,    17,    18,    19,    23,    25,
      24,    20,    21,    22,    26,    27,    28,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,   385,     0,     0,   201,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   245,   292,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   127,
       0,     0,     0,     0,     0,     0,   372,     0,     0,     0,
      75,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     206,     0,   196,     0,   212,     0,     0,   419,     0,     0,
       0,    57,     0,    63,     0,    69,     0,     0,   456,     0,
       1,     3,   446,     0,   559,     0,     0,     0,     0,     0,
     550,     0,     0,     0,   551,     0,     0,     0,     0,   434,
     445,     0,   435,   440,   438,   441,   439,   436,   437,   442,
     443,   427,   428,   429,   430,   431,   432,   433,   454,     0,
       0,     0,   448,   453,     0,   450,   449,   451,     0,     0,
     382,     0,   378,     0,     0,   204,   205,     0,    81,     0,
      48,   395,     0,     0,   389,     0,     0,     0,     0,   115,
       0,   534,     0,   539,   518,     0,     0,     0,     0,     0,
     531,   532,     0,     0,     0,     0,     0,     0,   552,   527,
       0,     0,   538,   533,   517,     0,     0,   537,   535,     0,
     297,   333,   322,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   323,   324,   325,
     326,   327,   328,   329,   330,   331,   332,   334,   335,   336,
     241,     0,   294,     0,   258,     0,     0,   255,     0,     0,
       0,     0,     0,   277,     0,     0,     0,     0,   271,     0,
       0,   119,     0,     0,     0,     0,    83,     0,   492,     0,
       0,     0,     0,   546,   545,     0,   401,   402,   403,   404,
       0,     0,     0,   170,    88,    89,    87,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   363,     0,     0,     0,     0,     0,     0,
       0,     0,   497,   498,   499,     0,     0,     0,   540,     0,
     506,     0,     0,     0,   218,   219,   220,   221,   222,   223,
     224,   225,   226,   227,   228,   230,   232,   233,   234,   235,
     236,   237,   238,   229,   231,   239,   240,   374,   371,    78,
      73,     0,    54,     0,    79,   145,   146,   165,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   417,   144,     0,   340,   345,   341,   342,   343,   344,
     346,   347,   348,   349,   350,   351,   352,   353,     0,    50,
       0,     0,     0,   209,   210,   211,     0,   199,   200,     0,
     217,   214,     0,   425,     0,   424,   426,   421,   365,    60,
      55,     0,    51,    66,    61,     0,    52,    72,    67,     0,
      53,   360,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     459,   460,   461,   462,   463,   464,   465,   466,   467,   468,
     469,   470,   471,   472,   473,   474,   475,   476,   477,   486,
     478,   479,   480,   481,   482,   483,   484,   485,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   376,   377,     0,     0,     0,    82,    49,
       0,     0,     0,     0,     0,     0,     0,   113,   114,   246,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   243,
       0,   257,   253,   254,   286,     0,   286,     0,   275,   276,
       0,   286,     0,   269,   270,     0,   117,   118,   110,     0,
       0,    84,     0,     0,   139,     0,   140,     0,   135,     0,
       0,     0,     0,     0,     0,     0,   168,   169,     0,    95,
      96,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,     0,   361,
     362,     0,     0,   366,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    76,    74,
      80,     0,   152,   153,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   167,   194,   195,     0,     0,   186,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    58,
      56,    64,    62,    70,    68,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   488,   487,   555,   252,     0,     0,
     556,   557,   558,   554,   560,   509,   512,   511,     0,     0,
     510,   513,   514,   547,     0,   548,   444,     0,   561,   519,
     536,   452,     0,   386,     0,   382,     0,     0,     0,   490,
     203,   202,   398,   393,     0,     0,   392,   387,     0,     0,
       0,   549,   500,   541,   542,   515,   516,   521,   524,   522,
     529,   530,   520,   526,   525,     0,   528,   296,   293,     0,
     242,     0,     0,   281,   288,   282,   287,   284,   289,   283,
     285,     0,     0,     0,   263,     0,     0,   286,     0,     0,
     286,   249,     0,     0,     0,   112,     0,     0,   128,   137,
     138,     0,   142,   124,   123,     0,     0,   122,   125,   126,
       0,   131,   129,   543,   544,   400,   411,   413,   414,   415,
     412,     0,   405,   409,     0,     0,     0,     0,   108,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      86,    91,    90,    92,    93,    94,     0,     0,   369,     0,
       0,   496,   504,   489,   495,   501,   502,   493,   503,   508,
     494,   491,   507,   373,     0,    77,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   143,   148,   147,
     149,   150,   151,   418,   339,   337,     0,   354,     0,     0,
       0,     0,   191,   192,     0,     0,   208,   207,   198,   197,
     216,   213,     0,   553,   423,   420,     0,    59,    65,    71,
     562,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
     582,   583,   584,   585,   586,   587,   588,   457,   458,   251,
     250,   590,   592,   594,   593,     0,   447,   455,     0,     0,
     384,   383,     0,   394,     0,   388,     0,   116,     0,   358,
       0,   295,   244,   259,   291,   290,   256,   286,   286,     0,
     286,     0,     0,   274,     0,   248,   247,     0,     0,     0,
       0,     0,     0,   133,     0,     0,     0,     0,   399,   286,
     410,     0,     0,     0,     0,     0,     0,   107,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,     0,     0,
       0,   367,   375,   166,   154,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   338,   355,   193,   185,   184,   188,
     189,     0,     0,   215,   422,     0,   589,   382,     0,   379,
     396,     0,   390,     0,     0,     0,   523,   260,     0,     0,
       0,   286,     0,   286,   286,   272,     0,   111,     0,   141,
     505,   121,     0,     0,     0,     0,   406,     0,     0,   173,
       0,   181,     0,     0,   109,   364,   370,     0,     0,   190,
     591,     0,     0,   397,   391,     0,   359,   286,   286,   286,
       0,   280,     0,     0,     0,   164,     0,   136,   132,   130,
     286,   407,     0,     0,     0,   176,     0,     0,   172,   368,
     187,     0,   380,   286,   265,   261,   264,   286,   278,   273,
     120,     0,     0,   175,   174,   180,     0,   178,     0,     0,
       0,   357,   286,     0,     0,   134,   408,   177,     0,   183,
       0,   381,     0,   266,     0,   279,   179,   171,   182,   356,
     286,   286,   267,   262,   286,   268
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1049, -1049,  1509, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,  -313, -1049,
   -1049, -1049, -1049,  -111,  -225, -1049, -1049,  1238, -1049,   528,
   -1049, -1049, -1049, -1049, -1049, -1049,  -943,  -603,  -132,  -592,
   -1049, -1049, -1049,  1425,  -234, -1049, -1049, -1049, -1049,   621,
   -1049, -1049,   850, -1049, -1049,  1000, -1049, -1049,   854, -1049,
   -1049,  -116,   -24,   888,  1025, -1049, -1049,  1264, -1049, -1049,
   -1048, -1049, -1049,  1255, -1049, -1049,  1261,  -958,  -567, -1049,
   -1049,   972, -1049,  1438,   879, -1049,   480, -1049, -1049, -1049,
   -1049,  1216, -1049, -1049, -1049, -1049, -1049, -1049, -1049,  1365,
    -749, -1049, -1049, -1049, -1049, -1049,   949, -1049,   542,  -823,
   -1049, -1049, -1049, -1049, -1049,   865, -1049,   -70,  1059, -1049,
   -1049,  1049, -1049, -1049,   853, -1049,  -460, -1049,   -93, -1049,
    1495, -1049, -1049, -1049, -1049, -1049, -1049, -1049,   -88, -1049,
   -1049,  -133,  -594, -1049, -1049, -1049, -1049,  -103,   -99,   -92,
     -90,   -87, -1049, -1049,   -83,   -69, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049,   -82, -1049, -1049, -1049, -1049, -1049,
     -77,   -67,   -68,   -66,   -64,   -51, -1049, -1049, -1049, -1049,
    -112,  -106,   -81,   -79,   -50,   -48,   -47, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049, -1049,
   -1049, -1049, -1049, -1049, -1049,   847, -1049,  -517
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    44,    45,    46,    47,    48,    49,    50,    51,    52,
     152,   154,   156,   131,    53,    54,    55,    56,   360,   879,
      57,   324,    58,   228,   229,    59,   320,   321,   856,   857,
      60,   327,  1036,  1035,  1112,   860,   623,   624,   625,   626,
     432,    61,    62,   342,   343,  1122,  1190,    63,   708,   709,
      64,   456,   457,    65,   214,   215,    66,   452,   453,    67,
     459,   463,   110,   844,   760,    68,   306,   307,   308,   832,
    1097,    69,   317,   318,    70,   312,   313,   833,  1098,    71,
     259,   260,    72,   433,   434,    73,  1009,  1010,    74,    75,
     362,   363,    76,    77,   365,    78,    79,    80,   211,   212,
     562,    81,    82,    83,    84,   335,   336,   871,   872,   873,
      85,   134,   700,    86,   464,   465,   179,   180,   181,    87,
     203,   204,    88,    89,   509,   510,   756,   384,   385,   386,
     387,   388,   389,   390,   391,   392,   393,   394,   395,   396,
     397,   398,   399,   859,   400,   401,   402,   182,   183,   184,
     185,   186,   268,   269,   403,   437,   272,   273,   274,   275,
     276,   277,   278,   279,   438,   281,   282,   283,   284,   285,
     439,   440,   441,   442,   443,   444,   404,   292,   293,   337,
     405,   406,   187,   188,   447,   189,   190,   299,   466,   191,
     192,   193,   194,   195,   196,   197,   207,   511,   512,   513,
     514,   515,   516,   517,   518,   519,   520,   521,   522,   523,
     524,   525,   526,   527,   528,   529,   530,   531,   532,   533,
     534,   535,   536,   537,   775,   992,   769,   770
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       128,   129,   431,   578,   322,   216,   263,   137,   454,   338,
     264,   383,   146,   149,   150,   339,   261,   265,   157,   266,
     849,   262,   267,   202,   205,   206,   270,   280,   294,   460,
     295,   850,   286,   858,   639,   640,   999,   834,   455,   836,
     271,   288,   287,   289,   839,   290,   651,   435,   435,  1040,
     458,   662,   436,   436,   445,   445,   446,   446,   291,   296,
    1099,   297,   298,   801,   851,   767,   942,   824,   848,   705,
     841,   991,    90,    92,   875,   943,   826,   102,   706,   415,
     761,  1155,   104,   209,  1147,   219,   328,  1079,   416,   560,
    1165,  1113,  1114,  1115,   779,   415,  1080,   578,   209,   577,
     616,   101,   618,   863,   416,   828,   121,    99,   637,   300,
     563,    96,   824,   107,   132,   107,   100,   417,   106,   866,
     566,   826,   107,   123,   823,   864,   807,   808,   809,   107,
     117,   107,   133,   417,   866,   816,   111,   107,   107,   118,
     107,   130,   329,   112,   107,   561,   107,   866,   227,   319,
     828,   227,   330,  1203,   842,   843,   171,   107,   944,   210,
     831,   340,   825,   301,  1156,   107,   695,   696,   697,   698,
     827,   699,   418,   419,   210,   876,   867,   567,   420,   421,
     422,   423,   424,   425,   426,   427,   428,   103,   418,   419,
     107,   867,   105,   429,   420,   421,   422,   423,   424,   425,
     426,   427,   428,  1171,   867,   831,   902,   829,   768,   429,
     107,    91,    93,   909,   945,   993,   637,   707,   302,  1157,
     852,  1180,   793,   220,   797,  1194,   978,   430,   376,   622,
     877,   818,   868,   659,  1037,  1081,   869,   870,   108,   109,
     126,   127,   830,   430,  1039,   622,   913,   868,  1148,   144,
     145,   869,   870,   147,   148,   935,   953,  1150,   300,   415,
     868,   300,   410,   800,   869,   870,   937,   341,   416,   300,
    1021,   221,  1149,  1024,   951,   663,  1044,    94,    95,   222,
     975,   976,   361,   415,   979,   980,   682,   683,   983,  1038,
     985,   986,   416,  1040,   300,   113,  1045,   417,   694,   955,
    1185,   652,   653,   654,   655,   114,   656,  1107,   122,   224,
     107,   415,   301,   754,   300,   301,   411,   225,   300,  1090,
     416,   417,   755,   301,   664,   878,   470,   474,   880,   881,
     882,   883,   884,   885,   886,   887,   888,   889,  1131,   891,
     892,   893,   894,   895,   303,   896,   478,   300,   301,   417,
     757,   900,   418,   419,   300,   135,   227,   890,   420,   421,
     422,   423,   424,   425,   426,   427,   428,   407,   301,   596,
     408,   412,   301,   429,   300,   300,   418,   419,   449,   309,
     471,   475,   420,   421,   422,   423,   424,   425,   426,   427,
     428,   300,   300,   300,   939,   124,   300,   429,   300,   314,
     479,   301,   304,   461,   418,   419,   136,   430,   301,   622,
     420,   421,   422,   423,   424,   425,   426,   427,   428,  1092,
     125,   303,   300,   467,   300,   429,   597,   468,   301,   301,
     602,   430,   340,   622,   627,   472,   476,   310,   858,   305,
     309,   636,   571,   758,   759,   301,   301,   301,   572,   608,
     301,  1100,   301,  1102,   574,   480,   481,   315,   139,   430,
     575,   622,   786,  1003,   940,   151,   138,   632,   340,   153,
     941,   787,  1117,   216,   311,  1002,   301,   704,   301,   304,
    1014,   155,   701,  1005,  1011,   162,   202,   205,   206,  1015,
     198,   628,   217,   263,   316,    97,    98,   264,   310,   208,
    1062,  1074,  1075,   261,   265,  1083,   266,  1084,   262,   267,
     849,   849,   849,   270,   280,   294,   305,   295,   701,   286,
     338,   850,   850,   850,   633,   218,   339,   271,   288,   287,
     289,  1133,   290,  1134,  1140,   311,  1142,  1143,   341,   702,
     823,   213,   415,   115,   116,   291,   296,   794,   297,   298,
     798,   416,   917,   918,   919,   920,   921,   922,   923,   924,
     925,   926,  1042,   928,   929,   930,   931,   932,   849,   710,
    1164,  1110,  1166,   819,   341,   703,   824,   454,   825,   850,
     417,   223,   314,  1172,  1059,   826,   827,   226,   367,   227,
     231,   613,   230,   712,   714,   717,  1181,  1032,   950,   319,
    1184,   323,  1046,  1095,   325,   200,   332,   455,   326,   853,
     119,   120,   361,   435,   828,  1193,   364,   333,   436,   458,
     445,   854,   446,   829,   232,   233,   711,   855,   201,   409,
     334,   234,  1082,  1202,   413,   418,   419,  1205,   235,   451,
     315,   420,   421,   422,   423,   424,   425,   426,   427,   428,
     713,   715,   718,   914,  1033,   414,   429,   538,   830,  1047,
    1096,   140,   141,   469,   252,   695,   696,   697,   698,   831,
     699,   254,   654,   655,   473,   656,   477,   316,   936,   938,
     652,   653,   654,   655,   539,   656,   540,   256,   142,   143,
     430,   952,   622,   541,   956,   697,   698,   542,   699,   257,
     158,   159,   543,   544,   545,   258,   652,   653,   654,   655,
     546,   656,   547,   548,   549,   550,   551,   177,   178,  1030,
     553,   927,   160,   652,   653,   654,   655,   554,   656,     1,
       2,   555,  1123,  1124,   552,  1028,  1048,   556,   557,     3,
       4,     5,   559,   558,   565,   568,     6,  1127,   570,   576,
       7,     8,     9,   580,    10,   569,    11,    12,    13,    14,
     581,   582,  1049,   573,   579,   583,   584,   585,  1128,    15,
     586,   587,    16,   588,   589,  1132,   590,   591,   592,  1050,
     593,   652,   653,   654,   655,    17,   656,   695,   696,   697,
     698,   594,   699,  1144,   595,   598,    18,    19,    20,   599,
     604,   600,    21,   578,   601,  1152,   605,  1153,   652,   653,
     654,   655,    22,   656,    23,   607,    24,    25,    26,    27,
      28,   615,   606,   611,   619,    29,    30,   620,   610,   612,
      31,    32,    33,    34,   621,   629,   631,  1051,   630,   634,
      35,    36,   665,    37,   635,     1,     2,    38,  1179,   638,
      39,    40,    41,    42,   666,     3,     4,     5,   652,   653,
     654,   655,     6,   656,  1052,  1189,     7,     8,     9,   667,
      10,   762,    11,    12,    13,    14,   641,  1198,   668,    43,
     642,   643,   644,   645,   344,    15,   646,   647,    16,   652,
     653,   654,   655,   345,   656,   648,   649,  1108,   650,   658,
     344,    17,   652,   653,   654,   655,   669,   656,   670,   345,
     671,   672,    18,    19,    20,   673,   344,   674,    21,   675,
     676,   677,   346,   678,   679,   345,   681,   716,    22,   684,
      23,   680,    24,    25,    26,    27,    28,   719,   346,   720,
     685,    29,    30,   686,   722,  1053,    31,    32,    33,    34,
     687,   721,   723,   688,   346,   689,    35,    36,  1054,    37,
     690,   691,   692,    38,   693,   724,    39,    40,    41,    42,
     652,   653,   654,   655,   725,   656,   661,   347,   348,   726,
    1091,   727,  1093,   349,   350,   351,   352,   353,   354,   355,
     356,   357,   899,   347,   348,    43,   728,   752,   358,   349,
     350,   351,   352,   353,   354,   355,   356,   357,   729,   347,
     348,   730,   731,   732,   358,   349,   350,   351,   352,   353,
     354,   355,   356,   357,   344,   733,  1055,   734,   735,   736,
     358,   753,   359,   345,   652,   653,   654,   655,   737,   656,
     415,   652,   653,   654,   655,   738,   656,   739,   359,   416,
     652,   653,   654,   655,   740,   656,   695,   696,   697,   698,
     741,   699,   346,  1197,   359,   695,   696,   697,   698,   742,
     699,   743,   783,   744,   745,   746,   747,   748,   417,   695,
     696,   697,   698,   749,   699,   750,   751,   765,   764,   773,
    1056,   774,   784,   695,   696,   697,   698,  1057,   699,   695,
     696,   697,   698,   785,   699,   766,  1058,   771,   652,   653,
     654,   655,  1064,   656,   772,   777,   778,   347,   348,   780,
     782,  1065,   788,   349,   350,   351,   352,   353,   354,   355,
     356,   357,   789,   418,   419,  1066,   791,   795,   358,   420,
     421,   422,   423,   424,   425,   426,   427,   428,   231,  1067,
     792,   796,   799,   768,   429,  1068,   802,   803,   804,   805,
     806,   810,   657,   200,   170,   811,   812,   813,   171,   821,
     814,   815,   359,   820,   695,   696,   697,   698,   845,   699,
     822,   835,   232,   233,   172,   837,   201,   989,   430,   234,
     695,   696,   697,   698,   838,   699,   235,   236,   237,   846,
     840,   238,   239,   847,   240,   241,   861,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   862,   695,   696,
     697,   698,   252,   699,   173,   174,   874,   253,   656,   254,
    1069,   897,   901,   903,   255,   695,   696,   697,   698,   904,
     699,   905,   175,   176,   906,   256,  1070,   907,   163,   164,
     165,   166,   167,   168,   169,   199,   915,   257,   213,   200,
     170,   908,   910,   258,   171,   911,   912,   695,   696,   697,
     698,   699,   699,   933,  1071,   177,   178,   957,   958,   959,
     172,   960,   201,   163,   164,   165,   166,   167,   168,   169,
     947,  1072,   949,   961,   962,   170,   963,   964,   965,   171,
     966,   366,   967,   652,   653,   654,   655,   968,   656,   652,
     653,   654,   655,   969,   656,   172,   970,   971,   972,   973,
     173,   174,   367,  1073,   368,   369,   652,   653,   654,   655,
     974,   656,   977,   981,   982,   984,   987,   991,   175,   176,
     652,   653,   654,   655,   234,   656,   370,   371,   990,   996,
     995,   235,   997,   998,   561,   173,   174,  1000,   328,  1076,
    1001,  1004,  1006,  1008,  1016,  1129,  1017,  1018,  1012,   366,
    1013,   177,   178,   175,   176,   652,   653,   654,   655,  1019,
     656,   372,  1160,   373,   254,   374,   333,  1020,  1022,  1023,
     367,   375,   368,   369,   898,   376,  1025,  1026,  1027,   334,
    1029,  1031,  1034,   377,   378,   379,   177,   178,  1041,   380,
     381,   382,   234,   213,   370,   371,  1043,    -1,  1060,   235,
     462,  1085,   652,   653,   654,   655,   328,   656,  1063,  1007,
     652,   653,   654,   655,  1086,   656,   695,   696,   697,   698,
    1088,   699,   652,   653,   654,   655,  1078,   656,  1101,   372,
    1087,   373,   254,   374,   333,   652,   653,   654,   655,   375,
     656,  1094,  1105,   376,  1103,  1104,  1106,   334,  1118,  1120,
    1121,   377,   378,   379,  1125,  1130,  1061,   380,   381,   382,
    1141,   213,  1137,  1138,  1089,  1135,   695,   696,   697,   698,
    1109,   699,   652,   653,   654,   655,  1119,   656,   652,   653,
     654,   655,  1146,   656,  1154,   652,   653,   654,   655,  1126,
     656,   652,   653,   654,   655,  1151,   656,   652,   653,   654,
     655,  1139,   656,   652,   653,   654,   655,  1161,   656,   652,
     653,   654,   655,  1163,   656,  1167,  1168,  1175,  1177,  1169,
    1145,  1178,  1182,  1183,  1186,  1176,  1158,  1187,  1188,  1192,
    1195,  1196,  1159,   161,  1199,  1200,  1201,  1204,   617,  1162,
    1111,  1077,   948,   450,   946,  1170,   790,   763,   817,   916,
     603,  1173,   614,   609,   448,  1136,   564,  1174,   660,  1116,
     934,   865,   954,  1191,   482,   483,   484,   485,   486,   487,
     488,   489,   490,   491,   492,   493,   494,   495,   496,   497,
     498,   499,   500,   501,   502,   781,   988,     0,   503,   504,
     776,   505,   506,   507,   508,   331,   994
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   134,   228,   115,    98,   109,    31,   141,   121,
     109,   127,    36,    37,    38,   121,   109,   109,    42,   109,
     623,   109,   109,    93,    93,    93,   109,   109,   109,   145,
     109,   623,   109,   627,   347,   348,   785,   604,   141,   606,
     109,   109,   109,   109,   611,   109,   359,   135,   136,   872,
     143,   364,   135,   136,   135,   136,   135,   136,   109,   109,
    1018,   109,   109,   580,    33,    51,    42,    42,    31,    42,
      51,    51,    51,    51,    33,    51,    51,    80,    51,    42,
     540,    51,    80,     4,    31,    51,    60,    42,    51,   134,
    1138,  1034,  1035,  1036,   554,    42,    51,   322,     4,    31,
      31,    51,    31,    56,    51,    80,   190,    33,   342,    80,
      31,   189,    42,    80,    33,    80,    42,    80,   189,    80,
     134,    51,    80,    33,     6,    78,   586,   587,   588,    80,
      33,    80,    51,    80,    80,   595,   189,    80,    80,    42,
      80,    80,   116,   189,    80,   190,    80,    80,    80,    80,
      80,    80,   126,  1201,   135,   136,    25,    80,   134,    80,
     135,    22,    44,   134,   134,    80,   135,   136,   137,   138,
      52,   140,   135,   136,    80,   134,   137,   191,   141,   142,
     143,   144,   145,   146,   147,   148,   149,   190,   135,   136,
      80,   137,   190,   156,   141,   142,   143,   144,   145,   146,
     147,   148,   149,  1146,   137,   135,   666,    89,   194,   156,
      80,   190,   190,   673,   190,   195,   450,   190,   189,   189,
     189,   196,   189,   189,   189,  1183,   743,   190,    97,   192,
     189,   189,   193,    31,   134,   190,   197,   198,   189,   190,
     189,   190,   124,   190,   190,   192,   189,   193,    31,   189,
     190,   197,   198,   189,   190,   189,   716,   190,    80,    42,
     193,    80,    80,   576,   197,   198,   189,   128,    51,    80,
     837,   189,    31,   840,   189,    31,    97,   189,   190,   197,
     740,   741,    80,    42,   744,   745,   418,   419,   748,   189,
     750,   751,    51,  1116,    80,   189,   117,    80,   430,   189,
      31,   135,   136,   137,   138,   189,   140,    31,   189,   189,
      80,    42,   134,    42,    80,   134,   134,   197,    80,   189,
      51,    80,    51,   134,    80,   638,    80,    80,   641,   642,
     643,   644,   645,   646,   647,   648,   649,   650,  1087,   652,
     653,   654,   655,   656,    22,   658,    80,    80,   134,    80,
      42,   664,   135,   136,    80,   190,    80,   191,   141,   142,
     143,   144,   145,   146,   147,   148,   149,   189,   134,   134,
     189,   189,   134,   156,    80,    80,   135,   136,   189,    22,
     134,   134,   141,   142,   143,   144,   145,   146,   147,   148,
     149,    80,    80,    80,   707,   189,    80,   156,    80,    22,
     134,   134,    80,   189,   135,   136,   190,   190,   134,   192,
     141,   142,   143,   144,   145,   146,   147,   148,   149,   189,
     189,    22,    80,   189,    80,   156,   191,   189,   134,   134,
      31,   190,    22,   192,   134,   189,   189,    80,  1032,   117,
      22,    31,   191,   135,   136,   134,   134,   134,   197,    31,
     134,  1018,   134,  1020,   191,   189,   189,    80,    33,   190,
     197,   192,    42,   189,   189,    80,   189,   134,    22,    80,
     195,    51,  1039,   566,   117,   788,   134,    31,   134,    80,
      42,    80,   134,   189,   189,   189,   556,   556,   556,    51,
     189,   191,    42,   596,   117,   189,   190,   596,    80,    23,
     189,   189,   189,   596,   596,   189,   596,   189,   596,   596,
    1113,  1114,  1115,   596,   596,   596,   117,   596,   134,   596,
     632,  1113,  1114,  1115,   191,   189,   632,   596,   596,   596,
     596,   189,   596,   189,  1101,   117,  1103,  1104,   128,   191,
       6,   115,    42,   189,   190,   596,   596,   571,   596,   596,
     574,    51,   684,   685,   686,   687,   688,   689,   690,   691,
     692,   693,   875,   695,   696,   697,   698,   699,  1171,   134,
    1137,  1031,  1139,   597,   128,   191,    42,   710,    44,  1171,
      80,    80,    22,  1150,   897,    51,    52,    80,    24,    80,
       5,    31,    33,   134,   134,   134,  1163,   134,   714,    80,
    1167,    38,   134,   134,    42,    20,    77,   710,   189,    45,
     189,   190,    80,   701,    80,  1182,    80,    88,   701,   712,
     701,    57,   701,    89,    39,    40,   191,    63,    43,   123,
     101,    46,   945,  1200,    51,   135,   136,  1204,    53,   194,
      80,   141,   142,   143,   144,   145,   146,   147,   148,   149,
     191,   191,   191,   677,   191,   189,   156,    33,   124,   191,
     191,   189,   190,   123,    79,   135,   136,   137,   138,   135,
     140,    86,   137,   138,   123,   140,   123,   117,   702,   703,
     135,   136,   137,   138,    33,   140,    33,   102,   189,   190,
     190,   715,   192,    33,   718,   137,   138,    33,   140,   114,
     189,   190,    33,    33,    33,   120,   135,   136,   137,   138,
      33,   140,    33,    33,    33,    33,   134,   132,   133,   851,
      33,   191,     0,   135,   136,   137,   138,    33,   140,     7,
       8,    33,  1045,  1046,   191,   846,   191,   134,   191,    17,
      18,    19,    80,   191,    33,   189,    24,  1060,    80,    33,
      28,    29,    30,    33,    32,   189,    34,    35,    36,    37,
      33,    33,   191,    80,    80,    33,    33,    33,  1081,    47,
      33,    33,    50,    33,    33,  1088,    33,    33,    33,   191,
      33,   135,   136,   137,   138,    63,   140,   135,   136,   137,
     138,    33,   140,  1106,    33,    33,    74,    75,    76,    80,
     134,    80,    80,  1028,    80,  1118,    80,  1120,   135,   136,
     137,   138,    90,   140,    92,    80,    94,    95,    96,    97,
      98,   190,   134,   134,    33,   103,   104,   191,    80,    80,
     108,   109,   110,   111,   189,   191,    33,   191,   191,    80,
     118,   119,    33,   121,    80,     7,     8,   125,  1161,   190,
     128,   129,   130,   131,    33,    17,    18,    19,   135,   136,
     137,   138,    24,   140,   191,  1178,    28,    29,    30,    33,
      32,    51,    34,    35,    36,    37,   190,  1190,    33,   157,
     190,   190,   190,   190,    42,    47,   190,   190,    50,   135,
     136,   137,   138,    51,   140,   190,   190,  1029,   190,   190,
      42,    63,   135,   136,   137,   138,    33,   140,    33,    51,
      33,    33,    74,    75,    76,    33,    42,    33,    80,    33,
      33,   191,    80,   123,    80,    51,   190,    33,    90,   190,
      92,   189,    94,    95,    96,    97,    98,   123,    80,    80,
     190,   103,   104,   190,    80,   191,   108,   109,   110,   111,
     190,   123,   123,   190,    80,   190,   118,   119,   191,   121,
     190,   190,   190,   125,   190,    80,   128,   129,   130,   131,
     135,   136,   137,   138,    33,   140,   134,   135,   136,    33,
    1004,    33,  1006,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   134,   135,   136,   157,    33,   191,   156,   141,
     142,   143,   144,   145,   146,   147,   148,   149,    33,   135,
     136,    33,    33,    33,   156,   141,   142,   143,   144,   145,
     146,   147,   148,   149,    42,    33,   191,    33,    33,    33,
     156,   134,   190,    51,   135,   136,   137,   138,    33,   140,
      42,   135,   136,   137,   138,    33,   140,    33,   190,    51,
     135,   136,   137,   138,    33,   140,   135,   136,   137,   138,
      33,   140,    80,   189,   190,   135,   136,   137,   138,    33,
     140,    33,   189,    33,    33,    33,    33,    33,    80,   135,
     136,   137,   138,    33,   140,    33,    33,    80,    51,    80,
     191,    80,   190,   135,   136,   137,   138,   191,   140,   135,
     136,   137,   138,    80,   140,    51,   191,    51,   135,   136,
     137,   138,   191,   140,    51,    51,    51,   135,   136,    51,
      51,   191,    33,   141,   142,   143,   144,   145,   146,   147,
     148,   149,    51,   135,   136,   191,   189,    80,   156,   141,
     142,   143,   144,   145,   146,   147,   148,   149,     5,   191,
     189,   189,    80,   194,   156,   191,    51,    51,    51,    51,
      51,    51,   189,    20,    21,    51,    51,    51,    25,    33,
      80,   190,   190,    80,   135,   136,   137,   138,    80,   140,
     134,   134,    39,    40,    41,   134,    43,    42,   190,    46,
     135,   136,   137,   138,   134,   140,    53,    54,    55,   189,
     134,    58,    59,    80,    61,    62,   189,    64,    65,    66,
      67,    68,    69,    70,    71,    72,    73,   189,   135,   136,
     137,   138,    79,   140,    81,    82,   134,    84,   140,    86,
     191,    80,    51,    51,    91,   135,   136,   137,   138,    51,
     140,    51,    99,   100,    51,   102,   191,    51,     9,    10,
      11,    12,    13,    14,    15,    16,   123,   114,   115,    20,
      21,    51,    51,   120,    25,    51,    51,   135,   136,   137,
     138,   140,   140,   189,   191,   132,   133,   123,   123,   123,
      41,    51,    43,     9,    10,    11,    12,    13,    14,    15,
     189,   191,   189,    51,    51,    21,    51,    51,    51,    25,
      51,     3,    51,   135,   136,   137,   138,    51,   140,   135,
     136,   137,   138,    51,   140,    41,    51,    51,    51,    51,
      81,    82,    24,   191,    26,    27,   135,   136,   137,   138,
      51,   140,    51,    51,    51,    51,   189,    51,    99,   100,
     135,   136,   137,   138,    46,   140,    48,    49,    42,   189,
     193,    53,   189,    51,   190,    81,    82,   191,    60,   191,
     191,   191,   191,   196,   189,   191,   134,   134,    80,     3,
      80,   132,   133,    99,   100,   135,   136,   137,   138,    80,
     140,    83,   191,    85,    86,    87,    88,   134,    80,   189,
      24,    93,    26,    27,   189,    97,    51,    51,   191,   101,
      33,    33,   189,   105,   106,   107,   132,   133,    80,   111,
     112,   113,    46,   115,    48,    49,    80,   140,    80,    53,
     122,   193,   135,   136,   137,   138,    60,   140,   191,   189,
     135,   136,   137,   138,    80,   140,   135,   136,   137,   138,
      33,   140,   135,   136,   137,   138,   189,   140,   134,    83,
     191,    85,    86,    87,    88,   135,   136,   137,   138,    93,
     140,    80,   189,    97,   134,   134,    33,   101,    33,    33,
      51,   105,   106,   107,   189,    51,   189,   111,   112,   113,
     189,   115,   134,   134,   189,   196,   135,   136,   137,   138,
     189,   140,   135,   136,   137,   138,   189,   140,   135,   136,
     137,   138,   189,   140,   193,   135,   136,   137,   138,   189,
     140,   135,   136,   137,   138,   191,   140,   135,   136,   137,
     138,   134,   140,   135,   136,   137,   138,    33,   140,   135,
     136,   137,   138,   134,   140,   134,   189,    51,    51,   189,
     189,   127,   134,   134,   191,   193,   189,    51,   193,    80,
     189,    51,   189,    44,   196,   134,   134,   134,   320,   189,
    1032,   940,   712,   138,   710,   189,   566,   542,   596,   681,
     306,   189,   317,   312,   136,  1095,   211,   189,   362,  1037,
     701,   632,   717,   189,   158,   159,   160,   161,   162,   163,
     164,   165,   166,   167,   168,   169,   170,   171,   172,   173,
     174,   175,   176,   177,   178,   556,   753,    -1,   182,   183,
     551,   185,   186,   187,   188,   120,   769
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,    17,    18,    19,    24,    28,    29,    30,
      32,    34,    35,    36,    37,    47,    50,    63,    74,    75,
      76,    80,    90,    92,    94,    95,    96,    97,    98,   103,
     104,   108,   109,   110,   111,   118,   119,   121,   125,   128,
     129,   130,   131,   157,   200,   201,   202,   203,   204,   205,
     206,   207,   208,   213,   214,   215,   216,   219,   221,   224,
     229,   240,   241,   246,   249,   252,   255,   258,   264,   270,
     273,   278,   281,   284,   287,   288,   291,   292,   294,   295,
     296,   300,   301,   302,   303,   309,   312,   318,   321,   322,
      51,   190,    51,   190,   189,   190,   189,   189,   190,    33,
      42,    51,    80,   190,    80,   190,   189,    80,   189,   190,
     261,   189,   189,   189,   189,   189,   190,    33,    42,   189,
     190,   190,   189,    33,   189,   189,   189,   190,   261,   261,
      80,   212,    33,    51,   310,   190,   190,   261,   189,    33,
     189,   190,   189,   190,   189,   190,   261,   189,   190,   261,
     261,    80,   209,    80,   210,    80,   211,   261,   189,   190,
       0,   201,   189,     9,    10,    11,    12,    13,    14,    15,
      21,    25,    41,    81,    82,    99,   100,   132,   133,   315,
     316,   317,   346,   347,   348,   349,   350,   381,   382,   384,
     385,   388,   389,   390,   391,   392,   393,   394,   189,    16,
      20,    43,   316,   319,   320,   354,   371,   395,    23,     4,
      80,   297,   298,   115,   253,   254,   327,    42,   189,    51,
     189,   189,   197,    80,   189,   197,    80,    80,   222,   223,
      33,     5,    39,    40,    46,    53,    54,    55,    58,    59,
      61,    62,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    79,    84,    86,    91,   102,   114,   120,   279,
     280,   327,   337,   346,   347,   348,   349,   350,   351,   352,
     353,   354,   355,   356,   357,   358,   359,   360,   361,   362,
     363,   364,   365,   366,   367,   368,   369,   370,   371,   372,
     373,   374,   376,   377,   381,   382,   383,   384,   385,   386,
      80,   134,   189,    22,    80,   117,   265,   266,   267,    22,
      80,   117,   274,   275,    22,    80,   117,   271,   272,    80,
     225,   226,   222,    38,   220,    42,   189,   230,    60,   116,
     126,   329,    77,    88,   101,   304,   305,   378,   379,   380,
      22,   128,   242,   243,    42,    51,    80,   135,   136,   141,
     142,   143,   144,   145,   146,   147,   148,   149,   156,   190,
     217,    80,   289,   290,    80,   293,     3,    24,    26,    27,
      48,    49,    83,    85,    87,    93,    97,   105,   106,   107,
     111,   112,   113,   260,   326,   327,   328,   329,   330,   331,
     332,   333,   334,   335,   336,   337,   338,   339,   340,   341,
     343,   344,   345,   353,   375,   379,   380,   189,   189,   123,
      80,   134,   189,    51,   189,    42,    51,    80,   135,   136,
     141,   142,   143,   144,   145,   146,   147,   148,   149,   156,
     190,   237,   239,   282,   283,   337,   353,   354,   363,   369,
     370,   371,   372,   373,   374,   381,   382,   383,   282,   189,
     242,   194,   256,   257,   340,   346,   250,   251,   327,   259,
     260,   189,   122,   260,   313,   314,   387,   189,   189,   123,
      80,   134,   189,   123,    80,   134,   189,   123,    80,   134,
     189,   189,   158,   159,   160,   161,   162,   163,   164,   165,
     166,   167,   168,   169,   170,   171,   172,   173,   174,   175,
     176,   177,   178,   182,   183,   185,   186,   187,   188,   323,
     324,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,   134,   191,    33,    33,    33,   134,   191,   191,    80,
     134,   190,   299,    31,   298,    33,   134,   191,   189,   189,
      80,   191,   197,    80,   191,   197,    33,    31,   223,    80,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,   134,   191,    33,    80,
      80,    80,    31,   266,   134,    80,   134,    80,    31,   275,
      80,   134,    80,    31,   272,   190,    31,   226,    31,    33,
     191,   189,   192,   235,   236,   237,   238,   134,   191,   191,
     191,    33,   134,   191,    80,    80,    31,   243,   190,   217,
     217,   190,   190,   190,   190,   190,   190,   190,   190,   190,
     190,   217,   135,   136,   137,   138,   140,   189,   190,    31,
     290,   134,   217,    31,    80,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,   191,   123,    80,
     189,   190,   237,   237,   190,   190,   190,   190,   190,   190,
     190,   190,   190,   190,   237,   135,   136,   137,   138,   140,
     311,   134,   191,   191,    31,    42,    51,   190,   247,   248,
     134,   191,   134,   191,   134,   191,    33,   134,   191,   123,
      80,   123,    80,   123,    80,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   191,   134,    42,    51,   325,    42,   135,   136,
     263,   325,    51,   263,    51,    80,    51,    51,   194,   425,
     426,    51,    51,    80,    80,   423,   317,    51,    51,   325,
      51,   320,    51,   189,   190,    80,    42,    51,    33,    51,
     254,   189,   189,   189,   261,    80,   189,   189,   261,    80,
     217,   426,    51,    51,    51,    51,    51,   325,   325,   325,
      51,    51,    51,    51,    80,   190,   325,   280,   189,   261,
      80,    33,   134,     6,    42,    44,    51,    52,    80,    89,
     124,   135,   268,   276,   277,   134,   277,   134,   134,   277,
     134,    51,   135,   136,   262,    80,   189,    80,    31,   236,
     238,    33,   189,    45,    57,    63,   227,   228,   341,   342,
     234,   189,   189,    56,    78,   305,    80,   137,   193,   197,
     198,   306,   307,   308,   134,    33,   134,   189,   217,   218,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     191,   217,   217,   217,   217,   217,   217,    80,   189,   134,
     217,    51,   325,    51,    51,    51,    51,    51,    51,   325,
      51,    51,    51,   189,   261,   123,   262,   237,   237,   237,
     237,   237,   237,   237,   237,   237,   237,   191,   237,   237,
     237,   237,   237,   189,   283,   189,   261,   189,   261,   217,
     189,   195,    42,    51,   134,   190,   257,   189,   251,   189,
     260,   189,   261,   325,   314,   189,   261,   123,   123,   123,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,   325,   325,    51,   426,   325,
     325,    51,    51,   325,    51,   325,   325,   189,   323,    42,
      42,    51,   424,   195,   424,   193,   189,   189,    51,   299,
     191,   191,   217,   189,   191,   189,   191,   189,   196,   285,
     286,   189,    80,    80,    42,    51,   189,   134,   134,    80,
     134,   277,    80,   189,   277,    51,    51,   191,   222,    33,
     237,    33,   134,   191,   189,   232,   231,   134,   189,   190,
     308,    80,   217,    80,    97,   117,   134,   191,   191,   191,
     191,   191,   191,   191,   191,   191,   191,   191,   191,   217,
      80,   189,   189,   191,   191,   191,   191,   191,   191,   191,
     191,   191,   191,   191,   189,   189,   191,   248,   189,    42,
      51,   190,   217,   189,   189,   193,    80,   191,    33,   189,
     189,   261,   189,   261,    80,   134,   191,   269,   277,   276,
     277,   134,   277,   134,   134,   189,    33,    31,   237,   189,
     325,   228,   233,   235,   235,   235,   307,   277,    33,   189,
      33,    51,   244,   217,   217,   189,   189,   217,   217,   191,
      51,   299,   217,   189,   189,   196,   285,   134,   134,   134,
     277,   189,   277,   277,   217,   189,   189,    31,    31,    31,
     190,   191,   217,   217,   193,    51,   134,   189,   189,   189,
     191,    33,   189,   134,   277,   269,   277,   134,   189,   189,
     189,   235,   277,   189,   189,    51,   193,    51,   127,   217,
     196,   277,   134,   134,   277,    31,   191,    51,   193,   217,
     245,   189,    80,   277,   276,   189,    51,   189,   217,   196,
     134,   134,   277,   269,   134,   277
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
     435,   436,   437,   438,   439,   440,   441,   442,   443,    59,
      40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   199,   200,   200,   201,   201,   201,   201,   201,   201,
     201,   201,   201,   201,   201,   201,   201,   201,   201,   201,
     201,   201,   201,   201,   201,   201,   201,   201,   201,   201,
     201,   201,   201,   201,   201,   201,   201,   201,   201,   201,
     201,   201,   201,   201,   202,   202,   202,   202,   203,   203,
     204,   205,   206,   207,   208,   209,   209,   209,   209,   209,
     209,   210,   210,   210,   210,   210,   210,   211,   211,   211,
     211,   211,   211,   212,   212,   212,   212,   212,   212,   213,
     213,   214,   214,   215,   215,   216,   217,   217,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   218,   218,
     219,   219,   220,   221,   222,   222,   223,   224,   225,   225,
     226,   227,   227,   228,   228,   228,   228,   230,   229,   231,
     229,   232,   229,   233,   229,   234,   229,   235,   235,   235,
     235,   236,   236,   237,   237,   237,   237,   237,   237,   237,
     237,   237,   237,   237,   237,   237,   237,   237,   237,   237,
     237,   237,   237,   237,   238,   239,   239,   240,   241,   242,
     242,   243,   243,   243,   243,   243,   244,   244,   244,   244,
     244,   244,   245,   245,   246,   247,   247,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   249,   249,   250,   250,
     251,   252,   252,   253,   253,   254,   255,   255,   256,   256,
     257,   257,   258,   258,   258,   258,   259,   259,   260,   260,
     260,   260,   260,   260,   260,   260,   260,   260,   260,   260,
     260,   260,   260,   260,   260,   260,   260,   260,   260,   260,
     260,   261,   261,   261,   261,   261,   261,   262,   262,   262,
     263,   263,   263,   264,   265,   265,   266,   267,   267,   267,
     268,   268,   268,   268,   268,   269,   269,   269,   269,   270,
     271,   271,   272,   272,   272,   273,   274,   274,   275,   275,
     275,   276,   276,   276,   276,   276,   277,   277,   277,   277,
     277,   277,   278,   278,   278,   278,   279,   279,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   281,   281,   282,
     282,   283,   283,   283,   283,   283,   283,   283,   283,   283,
     283,   283,   283,   283,   284,   284,   285,   285,   286,   286,
     287,   288,   289,   289,   290,   291,   292,   293,   293,   293,
     293,   294,   295,   295,   295,   295,   296,   297,   297,   298,
     298,   298,   299,   299,   299,   300,   300,   301,   301,   301,
     301,   301,   301,   302,   302,   302,   302,   302,   302,   303,
     304,   304,   305,   305,   305,   306,   306,   306,   306,   307,
     307,   308,   308,   308,   308,   308,   310,   311,   309,   312,
     312,   312,   312,   313,   313,   314,   314,   315,   315,   315,
     315,   315,   315,   315,   316,   316,   316,   316,   316,   316,
     316,   316,   316,   316,   317,   317,   318,   318,   319,   319,
     319,   319,   320,   320,   321,   321,   322,   322,   323,   323,
     324,   324,   324,   324,   324,   324,   324,   324,   324,   324,
     324,   324,   324,   324,   324,   324,   324,   324,   324,   324,
     324,   324,   324,   324,   324,   324,   324,   325,   325,   326,
     327,   328,   329,   330,   331,   332,   333,   334,   335,   336,
     337,   338,   339,   340,   341,   342,   343,   344,   345,   346,
     347,   347,   348,   349,   350,   351,   352,   353,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   378,   379,   380,   381,   382,   383,
     384,   385,   386,   387,   388,   389,   390,   391,   392,   393,
     394,   395,   396,   397,   398,   399,   400,   401,   402,   403,
     404,   405,   406,   407,   408,   409,   410,   411,   412,   413,
     414,   415,   416,   417,   418,   419,   420,   421,   422,   423,
     424,   424,   425,   425,   426
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
       4,     4,     4,     4,     4,     4,     4,     4,     1,     3,
       4,     7,     3,     4,     2,     1,     4,     4,     2,     1,
       7,     3,     1,     1,     1,     1,     1,     0,     5,     0,
       8,     0,     8,     0,    10,     0,     8,     2,     2,     1,
       1,     4,     2,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     2,     2,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     5,     1,     4,     4,     4,     2,
       1,     9,     6,     5,     7,     7,     2,     4,     3,     5,
       3,     1,     2,     1,     6,     3,     1,     5,     3,     3,
       4,     2,     2,     3,     1,     1,     2,     5,     3,     1,
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
       1,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       3,     3,     1,     3,     3,     3,     3,     1,     1,     1,
       3,     3,     3,     3,     3,     3,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     3,
       3,     3,     3,     5,     3,     3,     3,     1,     3,     3,
       3,     1,     1,     1,     1,     1,     3,     1,     1,     1,
       1,     3,     3,     3,     3,     1,     1,     3,     3,     3,
       1,     1,     1,     3,     3,     3,     3,     3,     3,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       1,     3,     2,     2,     2
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
  "LIK_INIT", "LINEAR", "LOAD_MH_FILE", "LOGLINEAR", "MARKOWITZ",
  "MH_DROP", "MH_INIT_SCALE", "MH_JSCALE", "MH_MODE", "MH_NBLOCKS",
  "MH_REPLIC", "MH_RECOVER", "MODE_CHECK", "MODE_COMPUTE", "MODE_FILE",
  "MODEL", "MODEL_COMPARISON", "MSHOCKS", "MODEL_COMPARISON_APPROXIMATION",
  "MODIFIEDHARMONICMEAN", "MOMENTS_VARENDO", "NAME", "NOBS", "NOCONSTANT",
  "NOCORR", "NODIAGNOSTIC", "NOFUNCTIONS", "NOGRAPH", "NOMOMENTS",
  "NOPRINT", "NORMAL_PDF", "OBSERVATION_TRENDS", "OPTIM", "OPTIM_WEIGHTS",
  "ORDER", "OSR", "OSR_PARAMS", "PARAMETERS", "PERIODS",
  "PLANNER_OBJECTIVE", "PREFILTER", "PRESAMPLE", "PRINT", "PRIOR_TRUNC",
  "PRIOR_ANALYSIS", "POSTERIOR_ANALYSIS", "QZ_CRITERIUM", "RELATIVE_IRF",
  "REPLIC", "RPLOT", "SHOCKS", "SIGMA_E", "SIMUL", "SIMUL_ALGO",
  "SIMUL_SEED", "SMOOTHER", "SOLVE_ALGO", "SPARSE_DLL", "STDERR", "STEADY",
  "STOCH_SIMUL", "TEX", "RAMSEY_POLICY", "PLANNER_DISCOUNT", "TEX_NAME",
  "UNIFORM_PDF", "UNIT_ROOT_VARS", "USE_DLL", "VALUES", "VAR", "VAREXO",
  "VAREXO_DET", "VAROBS", "XLS_SHEET", "XLS_RANGE", "COMMA", "MINUS",
  "PLUS", "DIVIDE", "TIMES", "UMINUS", "POWER", "EXP", "LOG", "LOG10",
  "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN", "SINH", "COSH", "TANH",
  "ASINH", "ACOSH", "ATANH", "SQRT", "DYNARE_SENSITIVITY",
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
       200,     0,    -1,   201,    -1,   200,   201,    -1,   202,    -1,
     213,    -1,   214,    -1,   215,    -1,   229,    -1,   219,    -1,
     221,    -1,   224,    -1,   216,    -1,   240,    -1,   241,    -1,
     246,    -1,   249,    -1,   252,    -1,   255,    -1,   258,    -1,
     278,    -1,   281,    -1,   284,    -1,   264,    -1,   273,    -1,
     270,    -1,   287,    -1,   288,    -1,   291,    -1,   203,    -1,
     204,    -1,   292,    -1,   294,    -1,   295,    -1,   296,    -1,
     300,    -1,   301,    -1,   302,    -1,   303,    -1,   309,    -1,
     312,    -1,   318,    -1,   321,    -1,   322,    -1,   208,    -1,
     205,    -1,   206,    -1,   207,    -1,    28,    51,   189,    -1,
      28,    51,    51,   189,    -1,   108,   261,   189,    -1,   128,
     209,   189,    -1,   129,   210,   189,    -1,   130,   211,   189,
      -1,    96,   212,   189,    -1,   209,    80,    -1,   209,   134,
      80,    -1,    80,    -1,   209,    80,   123,    -1,   209,   134,
      80,   123,    -1,    80,   123,    -1,   210,    80,    -1,   210,
     134,    80,    -1,    80,    -1,   210,    80,   123,    -1,   210,
     134,    80,   123,    -1,    80,   123,    -1,   211,    80,    -1,
     211,   134,    80,    -1,    80,    -1,   211,    80,   123,    -1,
     211,   134,    80,   123,    -1,    80,   123,    -1,   212,    80,
      -1,   212,   134,    80,    -1,    80,    -1,   212,    80,   123,
      -1,   212,   134,    80,   123,    -1,    80,   123,    -1,    97,
      51,   189,    -1,    97,    33,    51,   189,    -1,    24,    42,
     189,    -1,    24,    33,    42,   189,    -1,    63,    42,   189,
      -1,    63,    33,    42,   189,    -1,    80,    33,   217,   189,
      -1,   190,   217,   191,    -1,    80,    -1,    42,    -1,    51,
      -1,   217,   136,   217,    -1,   217,   135,   217,    -1,   217,
     137,   217,    -1,   217,   138,   217,    -1,   217,   140,   217,
      -1,   135,   217,    -1,   136,   217,    -1,   141,   190,   217,
     191,    -1,   142,   190,   217,   191,    -1,   143,   190,   217,
     191,    -1,   144,   190,   217,   191,    -1,   145,   190,   217,
     191,    -1,   146,   190,   217,   191,    -1,   147,   190,   217,
     191,    -1,   148,   190,   217,   191,    -1,   149,   190,   217,
     191,    -1,   156,   190,   217,   191,    -1,    80,   190,   218,
     191,    -1,   217,    -1,   218,   134,   217,    -1,    50,   189,
     222,    31,    -1,    50,   190,   220,   191,   189,   222,    31,
      -1,    38,    33,    80,    -1,    32,   189,   222,    31,    -1,
     222,   223,    -1,   223,    -1,    80,    33,   217,   189,    -1,
      47,   189,   225,    31,    -1,   225,   226,    -1,   226,    -1,
      80,   190,   262,   191,    33,   217,   189,    -1,   227,   134,
     228,    -1,   228,    -1,    57,    -1,    45,    -1,   341,    -1,
     342,    -1,    -1,    74,   189,   230,   235,    31,    -1,    -1,
      74,   190,   329,   191,   189,   231,   235,    31,    -1,    -1,
      74,   190,   126,   191,   189,   232,   235,    31,    -1,    -1,
      74,   190,   116,   134,   227,   191,   233,   189,   235,    31,
      -1,    -1,    74,   190,   116,   191,   234,   189,   235,    31,
      -1,   235,   236,    -1,   235,   238,    -1,   236,    -1,   238,
      -1,   237,    33,   237,   189,    -1,   237,   189,    -1,   190,
     237,   191,    -1,   239,    -1,    42,    -1,    51,    -1,   237,
     136,   237,    -1,   237,   135,   237,    -1,   237,   137,   237,
      -1,   237,   138,   237,    -1,   237,   140,   237,    -1,   135,
     237,    -1,   136,   237,    -1,   141,   190,   237,   191,    -1,
     142,   190,   237,   191,    -1,   143,   190,   237,   191,    -1,
     144,   190,   237,   191,    -1,   145,   190,   237,   191,    -1,
     146,   190,   237,   191,    -1,   147,   190,   237,   191,    -1,
     148,   190,   237,   191,    -1,   149,   190,   237,   191,    -1,
     156,   190,   237,   191,    -1,   192,    80,    33,   237,   189,
      -1,    80,    -1,    80,   190,   262,   191,    -1,   109,   189,
     242,    31,    -1,    76,   189,   242,    31,    -1,   242,   243,
      -1,   243,    -1,   128,    80,   189,    97,   244,   189,   127,
     245,   189,    -1,   128,    80,   189,   117,   217,   189,    -1,
     128,    80,    33,   217,   189,    -1,   128,    80,   134,    80,
      33,   217,   189,    -1,    22,    80,   134,    80,    33,   217,
     189,    -1,   244,    51,    -1,   244,    51,   193,    51,    -1,
     244,   134,    51,    -1,   244,   134,    51,   193,    51,    -1,
      51,   193,    51,    -1,    51,    -1,   245,   217,    -1,   217,
      -1,   110,    33,   194,   247,   195,   189,    -1,   247,   189,
     248,    -1,   248,    -1,   248,   134,   190,   217,   191,    -1,
     248,   134,    42,    -1,   248,   134,    51,    -1,   248,   190,
     217,   191,    -1,   248,    42,    -1,   248,    51,    -1,   190,
     217,   191,    -1,    42,    -1,    51,    -1,   118,   189,    -1,
     118,   190,   250,   191,   189,    -1,   250,   134,   251,    -1,
     251,    -1,   327,    -1,    19,   189,    -1,    19,   190,   253,
     191,   189,    -1,   253,   134,   254,    -1,   254,    -1,   327,
      -1,   111,   189,    -1,   111,   190,   256,   191,   189,    -1,
     256,   134,   257,    -1,   257,    -1,   340,    -1,   346,    -1,
     119,   189,    -1,   119,   190,   259,   191,   189,    -1,   119,
     261,   189,    -1,   119,   190,   259,   191,   261,   189,    -1,
     259,   134,   260,    -1,   260,    -1,   326,    -1,   327,    -1,
     328,    -1,   329,    -1,   330,    -1,   331,    -1,   332,    -1,
     333,    -1,   334,    -1,   335,    -1,   336,    -1,   353,    -1,
     337,    -1,   375,    -1,   338,    -1,   339,    -1,   340,    -1,
     341,    -1,   343,    -1,   344,    -1,   345,    -1,   379,    -1,
     380,    -1,   261,    80,    -1,   261,    80,    33,    80,    -1,
     261,   134,    80,    -1,   261,   134,    80,    33,    80,    -1,
      80,    -1,    80,    33,    80,    -1,   136,    51,    -1,   135,
      51,    -1,    51,    -1,   136,    42,    -1,   135,    42,    -1,
      42,    -1,    35,   189,   265,    31,    -1,   265,   266,    -1,
     266,    -1,   267,   134,   268,   189,    -1,   117,    80,    -1,
      80,    -1,    22,    80,   134,    80,    -1,   276,   134,   269,
      -1,   277,   134,   276,   134,   269,    -1,   277,   134,   277,
     134,   277,   134,   276,   134,   269,    -1,   277,    -1,   277,
     134,   277,   134,   277,    -1,   277,   134,   277,    -1,   277,
     134,   277,   134,   277,    -1,   277,   134,   277,   134,   277,
     134,   277,    -1,   277,   134,   277,   134,   277,   134,   277,
     134,   277,    -1,    37,   189,   271,    31,    -1,   271,   272,
      -1,   272,    -1,   117,    80,   134,   277,   189,    -1,    22,
      80,   134,    80,   134,   277,   189,    -1,    80,   134,   277,
     189,    -1,    36,   189,   274,    31,    -1,   274,   275,    -1,
     275,    -1,   117,    80,   134,   277,   134,   277,   189,    -1,
      22,    80,   134,    80,   134,   277,   134,   277,   189,    -1,
      80,   134,   277,   134,   277,   189,    -1,     6,    -1,    44,
      -1,    89,    -1,    52,    -1,   124,    -1,    -1,    51,    -1,
      42,    -1,    80,    -1,   135,    51,    -1,   135,    42,    -1,
      34,   189,    -1,    34,   190,   279,   191,   189,    -1,    34,
     261,   189,    -1,    34,   190,   279,   191,   261,   189,    -1,
     279,   134,   280,    -1,   280,    -1,   346,    -1,   347,    -1,
     348,    -1,   349,    -1,   350,    -1,   351,    -1,   352,    -1,
     353,    -1,   354,    -1,   355,    -1,   356,    -1,   357,    -1,
     358,    -1,   359,    -1,   360,    -1,   361,    -1,   362,    -1,
     363,    -1,   364,    -1,   365,    -1,   366,    -1,   367,    -1,
     368,    -1,   369,    -1,   337,    -1,   370,    -1,   371,    -1,
     372,    -1,   373,    -1,   374,    -1,   376,    -1,   377,    -1,
     381,    -1,   382,    -1,   383,    -1,   327,    -1,   384,    -1,
     385,    -1,   386,    -1,   103,   190,   282,   191,   189,    -1,
     103,   190,   282,   191,   261,   189,    -1,   282,   134,   283,
      -1,   283,    -1,   353,    -1,   354,    -1,   363,    -1,   369,
      -1,   337,    -1,   370,    -1,   371,    -1,   372,    -1,   373,
      -1,   374,    -1,   381,    -1,   382,    -1,   383,    -1,   104,
     190,   282,   191,   189,    -1,   104,   190,   282,   191,   261,
     189,    -1,   196,    80,   196,   134,   196,    80,   196,    -1,
     196,    80,   196,   134,   277,    -1,   285,    -1,   286,   134,
     285,    -1,   131,   261,   189,    -1,    90,   189,   289,    31,
      -1,   289,   290,    -1,   290,    -1,    80,   190,   217,   191,
     189,    -1,   125,   261,   189,    -1,    92,   189,   293,    31,
      -1,   293,    80,   217,   189,    -1,   293,    80,   134,    80,
     217,   189,    -1,    80,   217,   189,    -1,    80,   134,    80,
     217,   189,    -1,    95,   261,   189,    -1,    94,   189,    -1,
      94,   190,   260,   191,   189,    -1,    94,   261,   189,    -1,
      94,   190,   260,   191,   261,   189,    -1,    18,   189,   297,
      31,    -1,   297,   298,    -1,   298,    -1,    80,   299,    33,
     217,   189,    -1,    80,   134,    80,   299,    33,   217,   189,
      -1,     4,    80,   190,    51,   191,   299,    33,   217,   189,
      -1,    -1,   190,    51,   191,    -1,   190,    42,   191,    -1,
      17,   189,    -1,    17,   190,    23,   191,   189,    -1,    30,
     190,    80,   191,   189,    -1,    30,   190,    80,   191,   261,
     189,    -1,    30,    80,   189,    -1,    30,   190,    80,   197,
      80,   191,   189,    -1,    30,   190,    80,   197,    80,   191,
     261,   189,    -1,    30,    80,   197,    80,   189,    -1,    29,
     190,    80,   191,   189,    -1,    29,   190,    80,   191,   261,
     189,    -1,    29,    80,   189,    -1,    29,   190,    80,   197,
      80,   191,   189,    -1,    29,   190,    80,   197,    80,   191,
     261,   189,    -1,    29,    80,   197,    80,   189,    -1,    75,
     190,   304,   191,   306,   189,    -1,   304,   134,   305,    -1,
     305,    -1,   378,    -1,   379,    -1,   380,    -1,   307,    -1,
     306,   134,   307,    -1,   307,   190,   277,   191,    -1,   306,
     134,   307,   190,   277,   191,    -1,   308,    -1,   307,   308,
      -1,    80,    -1,   198,    -1,   137,    -1,   193,    -1,   197,
      -1,    -1,    -1,    98,   310,   237,   311,   189,    -1,   121,
     189,    -1,   121,   190,   313,   191,   189,    -1,   121,   261,
     189,    -1,   121,   190,   313,   191,   261,   189,    -1,   313,
     134,   314,    -1,   314,    -1,   260,    -1,   387,    -1,   388,
      -1,   389,    -1,   390,    -1,   391,    -1,   392,    -1,   393,
      -1,   394,    -1,   315,    -1,   346,    -1,   381,    -1,   382,
      -1,   348,    -1,   350,    -1,   347,    -1,   349,    -1,   384,
      -1,   385,    -1,   316,   134,   317,    -1,   316,    -1,     7,
      51,   189,    -1,     7,   190,   317,   191,    51,   189,    -1,
     316,    -1,   371,    -1,   354,    -1,   395,    -1,   319,   134,
     320,    -1,   319,    -1,     8,    51,   189,    -1,     8,   190,
     320,   191,    51,   189,    -1,   157,   189,    -1,   157,   190,
     323,   191,   189,    -1,   324,   134,   323,    -1,   324,    -1,
     396,    -1,   397,    -1,   398,    -1,   399,    -1,   400,    -1,
     401,    -1,   402,    -1,   403,    -1,   404,    -1,   405,    -1,
     406,    -1,   407,    -1,   408,    -1,   409,    -1,   410,    -1,
     411,    -1,   412,    -1,   413,    -1,   415,    -1,   416,    -1,
     417,    -1,   418,    -1,   419,    -1,   420,    -1,   421,    -1,
     422,    -1,   414,    -1,    51,    -1,    42,    -1,    26,    33,
      51,    -1,   115,    33,    51,    -1,   112,    33,    51,    -1,
      60,    -1,    93,    33,    51,    -1,   107,    33,    51,    -1,
      27,    33,    51,    -1,     3,    33,    51,    -1,    83,    -1,
      85,    -1,    87,    -1,    53,    33,    51,    -1,    48,    33,
      51,    -1,    49,    33,    51,    -1,    97,    33,    51,    -1,
      24,    33,   325,    -1,    63,    33,   325,    -1,   111,    -1,
     113,    33,    51,    -1,   105,    33,   325,    -1,    25,    33,
      80,    -1,    81,    33,   426,    -1,    81,    33,    51,    -1,
      41,    33,    51,    -1,    99,    33,    51,    -1,   100,    33,
      51,    -1,    58,    33,    51,    -1,    59,    33,    51,    -1,
      86,    -1,    46,    -1,    20,    33,   325,    -1,    69,    33,
      51,    -1,    64,    33,   325,    -1,    66,    33,   325,    -1,
      91,    33,   190,   286,   191,    -1,    65,    33,   325,    -1,
      73,    33,    80,    -1,    72,    33,    51,    -1,    71,    -1,
     102,    33,   325,    -1,    67,    33,    51,    -1,    68,    33,
      51,    -1,    61,    -1,    62,    -1,    84,    -1,     5,    -1,
     120,    -1,    43,    33,    51,    -1,   114,    -1,    79,    -1,
      40,    -1,   106,    -1,    54,    33,    51,    -1,    55,    33,
      51,    -1,    77,    33,    56,    -1,    77,    33,    78,    -1,
     101,    -1,    88,    -1,   132,    33,    80,    -1,   133,    33,
     423,    -1,    39,    33,   426,    -1,    21,    -1,    82,    -1,
      70,    -1,   122,    33,   325,    -1,    14,    33,   263,    -1,
       9,    33,   325,    -1,    11,    33,   263,    -1,    12,    33,
     325,    -1,    13,    33,    51,    -1,    10,    -1,    15,    33,
      51,    -1,    16,    33,    51,    -1,   158,    33,    51,    -1,
     159,    33,    51,    -1,   160,    33,    51,    -1,   161,    33,
      51,    -1,   162,    33,    51,    -1,   163,    33,    51,    -1,
     164,    33,    51,    -1,   165,    33,    51,    -1,   166,    33,
      51,    -1,   167,    33,    51,    -1,   168,    33,    51,    -1,
     169,    33,    51,    -1,   170,    33,    51,    -1,   171,    33,
      51,    -1,   172,    33,    51,    -1,   173,    33,   325,    -1,
     174,    33,   325,    -1,   175,    33,    51,    -1,   176,    33,
     426,    -1,   177,    33,   325,    -1,   178,    33,   325,    -1,
     182,    33,    51,    -1,   183,    33,    51,    -1,   185,    33,
     325,    -1,   186,    33,    51,    -1,   187,    33,   325,    -1,
     188,    33,   325,    -1,    80,   193,    80,    -1,    51,    -1,
      51,   193,    51,    -1,   194,   424,    -1,   425,   424,    -1,
     425,   195,    -1
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
     292,   297,   302,   307,   312,   317,   322,   327,   332,   334,
     338,   343,   351,   355,   360,   363,   365,   370,   375,   378,
     380,   388,   392,   394,   396,   398,   400,   402,   403,   409,
     410,   419,   420,   429,   430,   441,   442,   451,   454,   457,
     459,   461,   466,   469,   473,   475,   477,   479,   483,   487,
     491,   495,   499,   502,   505,   510,   515,   520,   525,   530,
     535,   540,   545,   550,   555,   561,   563,   568,   573,   578,
     581,   583,   593,   600,   606,   614,   622,   625,   630,   634,
     640,   644,   646,   649,   651,   658,   662,   664,   670,   674,
     678,   683,   686,   689,   693,   695,   697,   700,   706,   710,
     712,   714,   717,   723,   727,   729,   731,   734,   740,   744,
     746,   748,   750,   753,   759,   763,   770,   774,   776,   778,
     780,   782,   784,   786,   788,   790,   792,   794,   796,   798,
     800,   802,   804,   806,   808,   810,   812,   814,   816,   818,
     820,   822,   825,   830,   834,   840,   842,   846,   849,   852,
     854,   857,   860,   862,   867,   870,   872,   877,   880,   882,
     887,   891,   897,   907,   909,   915,   919,   925,   933,   943,
     948,   951,   953,   959,   967,   972,   977,   980,   982,   990,
    1000,  1007,  1009,  1011,  1013,  1015,  1017,  1018,  1020,  1022,
    1024,  1027,  1030,  1033,  1039,  1043,  1050,  1054,  1056,  1058,
    1060,  1062,  1064,  1066,  1068,  1070,  1072,  1074,  1076,  1078,
    1080,  1082,  1084,  1086,  1088,  1090,  1092,  1094,  1096,  1098,
    1100,  1102,  1104,  1106,  1108,  1110,  1112,  1114,  1116,  1118,
    1120,  1122,  1124,  1126,  1128,  1130,  1132,  1134,  1140,  1147,
    1151,  1153,  1155,  1157,  1159,  1161,  1163,  1165,  1167,  1169,
    1171,  1173,  1175,  1177,  1179,  1185,  1192,  1200,  1206,  1208,
    1212,  1216,  1221,  1224,  1226,  1232,  1236,  1241,  1246,  1253,
    1257,  1263,  1267,  1270,  1276,  1280,  1287,  1292,  1295,  1297,
    1303,  1311,  1321,  1322,  1326,  1330,  1333,  1339,  1345,  1352,
    1356,  1364,  1373,  1379,  1385,  1392,  1396,  1404,  1413,  1419,
    1426,  1430,  1432,  1434,  1436,  1438,  1440,  1444,  1449,  1456,
    1458,  1461,  1463,  1465,  1467,  1469,  1471,  1472,  1473,  1479,
    1482,  1488,  1492,  1499,  1503,  1505,  1507,  1509,  1511,  1513,
    1515,  1517,  1519,  1521,  1523,  1525,  1527,  1529,  1531,  1533,
    1535,  1537,  1539,  1541,  1543,  1547,  1549,  1553,  1560,  1562,
    1564,  1566,  1568,  1572,  1574,  1578,  1585,  1588,  1594,  1598,
    1600,  1602,  1604,  1606,  1608,  1610,  1612,  1614,  1616,  1618,
    1620,  1622,  1624,  1626,  1628,  1630,  1632,  1634,  1636,  1638,
    1640,  1642,  1644,  1646,  1648,  1650,  1652,  1654,  1656,  1658,
    1662,  1666,  1670,  1672,  1676,  1680,  1684,  1688,  1690,  1692,
    1694,  1698,  1702,  1706,  1710,  1714,  1718,  1720,  1724,  1728,
    1732,  1736,  1740,  1744,  1748,  1752,  1756,  1760,  1762,  1764,
    1768,  1772,  1776,  1780,  1786,  1790,  1794,  1798,  1800,  1804,
    1808,  1812,  1814,  1816,  1818,  1820,  1822,  1826,  1828,  1830,
    1832,  1834,  1838,  1842,  1846,  1850,  1852,  1854,  1858,  1862,
    1866,  1868,  1870,  1872,  1876,  1880,  1884,  1888,  1892,  1896,
    1898,  1902,  1906,  1910,  1914,  1918,  1922,  1926,  1930,  1934,
    1938,  1942,  1946,  1950,  1954,  1958,  1962,  1966,  1970,  1974,
    1978,  1982,  1986,  1990,  1994,  1998,  2002,  2006,  2010,  2014,
    2018,  2020,  2024,  2027,  2030
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    94,    94,    95,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   143,   144,   145,   146,   150,   151,
     154,   157,   161,   165,   169,   173,   175,   177,   179,   181,
     183,   188,   190,   192,   194,   196,   198,   203,   205,   207,
     209,   211,   213,   218,   220,   222,   224,   226,   228,   233,
     237,   244,   248,   255,   259,   267,   272,   274,   276,   278,
     280,   282,   284,   286,   288,   290,   292,   294,   296,   298,
     300,   302,   304,   306,   308,   310,   312,   314,   319,   321,
     325,   327,   332,   336,   341,   342,   346,   351,   356,   357,
     361,   365,   366,   370,   371,   372,   373,   377,   377,   378,
     378,   380,   380,   382,   382,   384,   384,   389,   390,   391,
     392,   396,   398,   403,   404,   405,   407,   409,   411,   413,
     415,   417,   419,   421,   423,   425,   427,   429,   431,   433,
     435,   437,   439,   441,   445,   449,   451,   456,   460,   464,
     465,   469,   471,   473,   475,   477,   482,   484,   486,   488,
     490,   492,   498,   500,   505,   510,   512,   517,   519,   521,
     523,   525,   527,   529,   531,   533,   538,   542,   546,   547,
     550,   554,   556,   560,   561,   564,   568,   570,   574,   575,
     578,   579,   583,   585,   587,   589,   593,   594,   597,   598,
     599,   600,   601,   602,   603,   604,   605,   606,   607,   608,
     609,   610,   611,   612,   613,   614,   615,   616,   617,   618,
     619,   623,   625,   627,   629,   631,   633,   638,   640,   642,
     647,   649,   651,   656,   661,   663,   668,   672,   677,   682,
     692,   697,   703,   713,   718,   729,   735,   743,   753,   767,
     771,   773,   777,   784,   793,   802,   806,   808,   812,   821,
     832,   844,   846,   848,   850,   852,   857,   858,   859,   860,
     861,   863,   870,   872,   874,   876,   881,   882,   885,   886,
     887,   888,   889,   890,   891,   892,   893,   894,   895,   896,
     897,   898,   899,   900,   901,   902,   903,   904,   905,   906,
     907,   908,   909,   910,   911,   912,   913,   914,   915,   916,
     917,   918,   919,   920,   921,   922,   923,   927,   929,   934,
     935,   939,   940,   941,   942,   943,   944,   945,   946,   947,
     948,   949,   950,   951,   955,   957,   962,   963,   967,   968,
     972,   977,   982,   983,   986,   990,   993,   997,   999,  1001,
    1003,  1007,  1010,  1011,  1012,  1013,  1016,  1020,  1021,  1024,
    1025,  1026,  1029,  1030,  1031,  1034,  1035,  1038,  1039,  1040,
    1041,  1042,  1043,  1045,  1046,  1047,  1048,  1049,  1050,  1052,
    1056,  1057,  1060,  1061,  1062,  1065,  1066,  1067,  1068,  1071,
    1072,  1075,  1076,  1077,  1078,  1079,  1082,  1082,  1082,  1085,
    1087,  1089,  1091,  1096,  1097,  1100,  1101,  1104,  1105,  1106,
    1107,  1108,  1109,  1110,  1113,  1114,  1115,  1116,  1117,  1118,
    1119,  1120,  1121,  1122,  1125,  1126,  1129,  1131,  1135,  1136,
    1137,  1138,  1141,  1142,  1145,  1147,  1151,  1153,  1157,  1158,
    1161,  1162,  1163,  1164,  1165,  1166,  1167,  1168,  1169,  1170,
    1171,  1172,  1173,  1174,  1175,  1176,  1177,  1178,  1179,  1180,
    1181,  1182,  1183,  1184,  1185,  1186,  1187,  1190,  1190,  1192,
    1193,  1194,  1195,  1196,  1197,  1198,  1199,  1200,  1201,  1202,
    1203,  1204,  1205,  1206,  1207,  1208,  1209,  1210,  1211,  1212,
    1213,  1214,  1216,  1217,  1218,  1219,  1220,  1221,  1222,  1223,
    1224,  1225,  1226,  1227,  1228,  1229,  1230,  1231,  1232,  1233,
    1234,  1235,  1236,  1237,  1238,  1239,  1240,  1241,  1242,  1243,
    1244,  1245,  1246,  1248,  1250,  1253,  1254,  1255,  1256,  1257,
    1258,  1259,  1260,  1261,  1263,  1264,  1265,  1266,  1267,  1268,
    1269,  1270,  1272,  1273,  1274,  1275,  1276,  1277,  1278,  1279,
    1280,  1281,  1282,  1283,  1284,  1285,  1286,  1287,  1288,  1289,
    1290,  1292,  1293,  1299,  1300,  1304,  1305,  1306,  1307,  1310,
    1318,  1319,  1323,  1324,  1333
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
       2,     2,     2,     2,     2,   192,     2,     2,     2,   196,
     190,   191,     2,     2,     2,     2,   197,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   193,   189,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   194,   198,   195,     2,     2,     2,     2,     2,     2,
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
     185,   186,   187,   188
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1616;
  const int parser::yynnts_ = 228;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 160;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 199;

  const unsigned int parser::yyuser_token_number_max_ = 443;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1335 "DynareBison.yy"


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

