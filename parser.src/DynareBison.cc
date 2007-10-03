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

  case 125:
#line 372 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 128:
#line 378 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 129:
#line 378 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 130:
#line 379 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 131:
#line 380 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 132:
#line 381 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 133:
#line 382 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 134:
#line 383 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 135:
#line 384 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 136:
#line 385 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 137:
#line 386 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 142:
#line 398 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 143:
#line 400 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val));;}
    break;

  case 144:
#line 404 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 146:
#line 407 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 147:
#line 409 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 148:
#line 411 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 149:
#line 413 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 150:
#line 415 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 151:
#line 417 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 152:
#line 419 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 153:
#line 421 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 154:
#line 423 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 155:
#line 425 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 156:
#line 427 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 157:
#line 429 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 158:
#line 431 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 159:
#line 433 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 160:
#line 435 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 161:
#line 437 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 162:
#line 439 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 163:
#line 441 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 164:
#line 443 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 165:
#line 447 "DynareBison.yy"
    {driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 166:
#line 451 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 167:
#line 453 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 168:
#line 457 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 169:
#line 461 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 172:
#line 471 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 173:
#line 473 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 174:
#line 475 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 175:
#line 477 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 176:
#line 479 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 177:
#line 484 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 178:
#line 486 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 179:
#line 488 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 180:
#line 490 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 181:
#line 492 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 182:
#line 494 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 183:
#line 500 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 184:
#line 502 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].node_val));;}
    break;

  case 185:
#line 507 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 186:
#line 512 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 187:
#line 514 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 188:
#line 519 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 189:
#line 521 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 190:
#line 523 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 191:
#line 525 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 192:
#line 527 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 193:
#line 529 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 194:
#line 531 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 195:
#line 533 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 196:
#line 535 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 197:
#line 540 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 198:
#line 544 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 202:
#line 556 "DynareBison.yy"
    {driver.check();;}
    break;

  case 203:
#line 558 "DynareBison.yy"
    {driver.check();;}
    break;

  case 207:
#line 570 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 208:
#line 572 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 213:
#line 585 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 214:
#line 587 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 215:
#line 589 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 216:
#line 591 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 242:
#line 625 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 243:
#line 627 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 244:
#line 629 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 245:
#line 631 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 246:
#line 633 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 247:
#line 635 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 248:
#line 640 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 249:
#line 642 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 250:
#line 644 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 251:
#line 649 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 252:
#line 651 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 253:
#line 653 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 254:
#line 658 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 255:
#line 663 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 256:
#line 665 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 258:
#line 674 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 259:
#line 679 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 260:
#line 684 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 261:
#line 694 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 262:
#line 699 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 263:
#line 705 "DynareBison.yy"
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

  case 264:
#line 715 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 265:
#line 720 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 266:
#line 731 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 267:
#line 737 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 268:
#line 745 "DynareBison.yy"
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

  case 269:
#line 755 "DynareBison.yy"
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

  case 270:
#line 769 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 271:
#line 773 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 272:
#line 775 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 273:
#line 779 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 274:
#line 786 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 275:
#line 795 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 276:
#line 804 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 277:
#line 808 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 278:
#line 810 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 279:
#line 814 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 280:
#line 823 "DynareBison.yy"
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

  case 281:
#line 834 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 282:
#line 846 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 283:
#line 848 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 284:
#line 850 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 285:
#line 852 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 286:
#line 854 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 287:
#line 858 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 291:
#line 863 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 292:
#line 865 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 293:
#line 872 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 294:
#line 874 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 295:
#line 876 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 296:
#line 878 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 338:
#line 929 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 339:
#line 931 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 355:
#line 957 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 356:
#line 959 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 357:
#line 963 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 358:
#line 964 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 361:
#line 974 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 362:
#line 979 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 365:
#line 988 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val));;}
    break;

  case 366:
#line 991 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 367:
#line 995 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 368:
#line 999 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 369:
#line 1001 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 370:
#line 1003 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 371:
#line 1005 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 372:
#line 1008 "DynareBison.yy"
    {driver.set_osr_params();;}
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
#line 1014 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 377:
#line 1018 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 380:
#line 1025 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 381:
#line 1026 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 382:
#line 1027 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val));;}
    break;

  case 383:
#line 1030 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 384:
#line 1031 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 385:
#line 1032 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 386:
#line 1035 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 387:
#line 1036 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 388:
#line 1039 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 389:
#line 1040 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 390:
#line 1041 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 391:
#line 1042 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 392:
#line 1043 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 393:
#line 1044 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 394:
#line 1046 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 395:
#line 1047 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 396:
#line 1048 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 397:
#line 1049 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 398:
#line 1050 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 399:
#line 1051 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 400:
#line 1054 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 406:
#line 1066 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 407:
#line 1067 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 408:
#line 1068 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 409:
#line 1069 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 410:
#line 1072 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 411:
#line 1073 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 413:
#line 1077 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 414:
#line 1078 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 415:
#line 1079 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 416:
#line 1080 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 417:
#line 1083 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 418:
#line 1083 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 420:
#line 1087 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 421:
#line 1089 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 422:
#line 1091 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 423:
#line 1093 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 447:
#line 1131 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 448:
#line 1133 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 455:
#line 1147 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 456:
#line 1149 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 457:
#line 1153 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 458:
#line 1155 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 490:
#line 1193 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 491:
#line 1194 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 492:
#line 1195 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 493:
#line 1196 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 494:
#line 1197 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1198 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 496:
#line 1199 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 497:
#line 1200 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 498:
#line 1201 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 499:
#line 1202 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 500:
#line 1203 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 501:
#line 1204 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 502:
#line 1205 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 503:
#line 1206 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 504:
#line 1207 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 505:
#line 1208 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 506:
#line 1209 "DynareBison.yy"
    {driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 507:
#line 1210 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 508:
#line 1211 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 509:
#line 1212 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 510:
#line 1213 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 511:
#line 1214 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 512:
#line 1215 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 513:
#line 1217 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 514:
#line 1218 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 515:
#line 1219 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 516:
#line 1220 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 517:
#line 1221 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 518:
#line 1222 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 519:
#line 1223 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 520:
#line 1224 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 521:
#line 1225 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 522:
#line 1226 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 523:
#line 1227 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 525:
#line 1229 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 526:
#line 1230 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 527:
#line 1231 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 528:
#line 1232 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 529:
#line 1233 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 530:
#line 1234 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 531:
#line 1235 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 532:
#line 1236 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 533:
#line 1237 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 534:
#line 1238 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 535:
#line 1239 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 536:
#line 1240 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 537:
#line 1241 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 538:
#line 1242 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 539:
#line 1243 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 540:
#line 1244 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 541:
#line 1245 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 542:
#line 1246 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 543:
#line 1247 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 544:
#line 1250 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 545:
#line 1252 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 546:
#line 1254 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 547:
#line 1255 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 548:
#line 1256 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 549:
#line 1257 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 550:
#line 1258 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 551:
#line 1259 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 552:
#line 1260 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 553:
#line 1261 "DynareBison.yy"
    {driver.option_num("mh_recover", "1");;}
    break;

  case 554:
#line 1262 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 555:
#line 1264 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 556:
#line 1265 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 557:
#line 1266 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 558:
#line 1267 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 559:
#line 1268 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 560:
#line 1269 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 561:
#line 1270 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 562:
#line 1271 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 563:
#line 1273 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 564:
#line 1274 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1275 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 566:
#line 1276 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1277 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 568:
#line 1278 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 569:
#line 1279 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 570:
#line 1280 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 571:
#line 1281 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1282 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1283 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1284 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 575:
#line 1285 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 576:
#line 1286 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 577:
#line 1287 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1288 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 579:
#line 1289 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1290 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 581:
#line 1291 "DynareBison.yy"
    {driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 582:
#line 1293 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1294 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1300 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1301 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1305 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1306 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1307 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1308 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1312 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 592:
#line 1321 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 593:
#line 1324 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 594:
#line 1326 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 595:
#line 1334 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2322 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1051;
  const short int
  parser::yypact_[] =
  {
       889,    17,    26,   -74,   -72,   127,   261,    95,    -1,    20,
     -23,     5,   -21,   109,   120,   125,   226,   262,   246,   -38,
     214,   102,   219,   222,    48,   153,   353,    68, -1051,   256,
     264,   153,   259,   427,   253,   390,    54,    61,   153,   387,
     391,   415,   153,   398,   762, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
     320,  1079,   343,  1248, -1051,   474,    82, -1051,   397,   501,
     404,    40,   -76,   515,    98,   519,   524,   565, -1051,  1132,
      50,   344,   360,   410,   533,   524,   578,   576,   429, -1051,
      30,   375,   269,   965,   541,   543, -1051,  1352,    84,   178,
     500,   191,   575,   437,   997,  1203,  1203,   202,   269,   433,
   -1051,    62, -1051,   397, -1051,  1352,   210, -1051,  1308,   217,
     232,   506,   267,   514,   273,   516,   283,   299, -1051,  1453,
   -1051, -1051, -1051,   596, -1051,   602,   625,   627,   632,   634,
   -1051,   642,   644,   651, -1051,   653,   654,   655,   656, -1051,
     555,   459, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,   659,
     660,   665, -1051,   564,   508, -1051, -1051, -1051,   520,   631,
     -43,    74, -1051,   680,   -42, -1051, -1051,   528, -1051,   536,
   -1051, -1051,   635,   -60, -1051,   647,   -56,   681,    64, -1051,
     648, -1051,   696, -1051, -1051,   712,   713,   715,   716,   720,
   -1051, -1051,   725,   727,   728,   730,   731,   732, -1051, -1051,
     733,   734, -1051, -1051, -1051,   738,   739, -1051, -1051,   -37,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
     740,   694, -1051,   697, -1051,   698,   407, -1051,   649,   702,
     666,   703,   507, -1051,   708,   667,   709,   554, -1051,   604,
      79, -1051,   368,   767,   611,   616, -1051,   594, -1051,   288,
     619,   639,   775, -1051, -1051,   291, -1051, -1051, -1051, -1051,
     737,   743,   118, -1051, -1051, -1051,   643,   965,   965,   650,
     652,   661,   663,   671,   672,   673,   679,   685,   688,   965,
     544,   689,   399, -1051,   798,   423,   793,   800,   802,   806,
     811,   812, -1051, -1051, -1051,   813,   814,   817, -1051,   818,
   -1051,   823,   835,   693, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
     763,   809, -1051,   700, -1051, -1051, -1051,   695,   997,   997,
     704,   718,   719,   721,   723,   724,   736,   741,   746,   758,
     997,   184, -1051,   292, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,   293, -1051,
     254,    24,   296, -1051, -1051, -1051,   307, -1051, -1051,   323,
   -1051, -1051,   866, -1051,   347, -1051, -1051, -1051, -1051, -1051,
     780,   831, -1051, -1051,   792,   842, -1051, -1051,   804,   850,
   -1051, -1051,   898,   917,   918,   920,   921,   925,   926,   927,
     928,   929,   935,   941,   944,   945,   946,   957,   958,   959,
     963,   964,   969,   970,   971,   972,   973,   977,   979,   784,
     878, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,   258,   274,
     258,   966,   274,   974,   934,   976,    18,   987,   990,   948,
     960,  1079,   991,   992,   258,   995,  1248,   998,   860,   862,
     975,   263,  1021, -1051, -1051,  1005,   397,   876, -1051, -1051,
     877,    41,   994,   885,    59,   996,   965, -1051, -1051, -1051,
     888,  1034,  1035,  1036,  1044,  1045,   258,   258,   258,  1046,
    1047,  1048,  1052,  1043,   938,   258,  1132,    63,  1051,  1099,
    1001, -1051, -1051, -1051,   345,  1003,   540,  1013, -1051, -1051,
    1014,   540,  1015, -1051, -1051,   329, -1051, -1051, -1051,  1071,
     968, -1051,  1075,    29, -1051,    87, -1051,    70, -1051,   978,
     984,    51,   375,   108,  1024,    28, -1051, -1051,   965,  1019,
     464,   965,   965,   965,   965,   965,   965,   965,   965,   965,
     965,   505,   965,   965,   965,   965,   965, -1051,   965, -1051,
   -1051,  1083,   980, -1051,   887,  1125,   258,  1126,  1131,  1133,
    1141,  1144,  1155,   258,  1156,  1159,  1165,    65, -1051,  1059,
   -1051,   329,  1068,   470,   997,   997,   997,   997,   997,   997,
     997,   997,   997,   997,   518,   997,   997,   997,   997,   997,
     999,  1203,    78,    90, -1051, -1051, -1051,   965,    93,    21,
      62,  1031,   397,  1032,  1352,   103,   258,  1308,   129, -1051,
    1094, -1051,  1105, -1051,  1107,  1183,  1189,  1193,  1194,  1199,
    1200,  1201,  1204,  1216,  1219,  1220,  1221,  1227,  1229,  1242,
     258,   258,  1247,   888,   258,   258,  1249,  1257,   258,  1258,
     258,   258,  1064,  1453, -1051, -1051, -1051, -1051,  1245,  1259,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,  1261,    14,
   -1051, -1051, -1051, -1051,  1116, -1051, -1051,  1123, -1051, -1051,
   -1051, -1051,  1130, -1051,  1270,  1134,  1135,  1137,   965, -1051,
   -1051, -1051, -1051, -1051,   304,  1146, -1051, -1051,   306,  1148,
    1138, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051,  1129, -1051, -1051, -1051,   310,
   -1051,  1243,  1253, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051,   535,  1151,  1187,  1207,  1263,  1209,   540,  1265,  1157,
     540, -1051,  1295,  1300,  1160, -1051,   524,  1320, -1051, -1051,
   -1051,   997, -1051, -1051, -1051,  1326, -1051,   400, -1051, -1051,
   -1051,  1170, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051,   -29,    -7, -1051,  1286,   965,  1289,   183,   511,
     405,   525,   532,   584,   618,   677,   683,   691,   764,   834,
     932, -1051,   464,   464,  1019,  1019,  1230,   943,   965, -1051,
    1303,  1339, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051,   311, -1051,  1188,   989,  1028,
    1089,  1100,  1147,  1158,  1166,  1178,  1226,  1236, -1051,   470,
     470,  1068,  1068,  1230, -1051, -1051, -1051,   315, -1051,   321,
    1250,    24,  1195, -1051, -1051,    25,   965, -1051, -1051, -1051,
   -1051, -1051, -1051,   326, -1051, -1051, -1051,   333, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051,  1190, -1051, -1051, -1051,  1310, -1051, -1051,  1211,
    1360, -1051, -1051,  1345, -1051,   198, -1051,   204, -1051,  1319,
   -1051,   439, -1051, -1051, -1051, -1051, -1051, -1051,   540,   345,
    1269,   540,  1272,  1273, -1051,  1223, -1051, -1051,  1376,   428,
     997,  1353,   258,    70, -1051,   594,   594,   594,   108, -1051,
     540, -1051,  1378,  1359,  1384,  1368,   965,   965, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,  1233,
    1365,   965, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,    21, -1051,
   -1051, -1051,   965,  1296, -1051, -1051,  1374, -1051,  1134,   965,
   -1051, -1051,   334, -1051,   340,  1232,  1129, -1051, -1051,  1291,
    1292,  1309,   540,  1240,   540,   540, -1051,   965, -1051,  1371,
   -1051, -1051, -1051,  1255,    57,   212,   228,    43,  1251,   965,
   -1051,   965,  1254,    45,  1377,   511, -1051, -1051,  1383,  1333,
   -1051, -1051,  1414,  1395, -1051, -1051,  1314, -1051,   540,   540,
     540,  1316, -1051,  1262,  1266,  1401, -1051,   594, -1051, -1051,
   -1051,   540, -1051,  1409,  1415,  1402,  1267,  1404,  1329, -1051,
   -1051, -1051,   965, -1051,    32,  1327, -1051,  1328,   540, -1051,
   -1051, -1051,   559,  1275, -1051, -1051, -1051,  1422,  1285,   965,
    1421,  1405, -1051,   540,   221,  1297, -1051, -1051, -1051,  1442,
     511,   915, -1051,  1302,  1370,  1382, -1051, -1051, -1051,   511,
   -1051,   540,   540,  1388, -1051,   540, -1051
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   417,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     4,    29,    30,    45,
      46,    47,    44,     5,     6,     7,    12,     9,    10,    11,
       8,    13,    14,    15,    16,    17,    18,    19,    23,    25,
      24,    20,    21,    22,    26,    27,    28,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,   386,     0,     0,   202,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   246,   293,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   128,
       0,     0,     0,     0,     0,     0,   373,     0,     0,     0,
      75,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     207,     0,   197,     0,   213,     0,     0,   420,     0,     0,
       0,    57,     0,    63,     0,    69,     0,     0,   457,     0,
       1,     3,   447,     0,   560,     0,     0,     0,     0,     0,
     551,     0,     0,     0,   552,     0,     0,     0,     0,   435,
     446,     0,   436,   441,   439,   442,   440,   437,   438,   443,
     444,   428,   429,   430,   431,   432,   433,   434,   455,     0,
       0,     0,   449,   454,     0,   451,   450,   452,     0,     0,
     383,     0,   379,     0,     0,   205,   206,     0,    81,     0,
      48,   396,     0,     0,   390,     0,     0,     0,     0,   115,
       0,   535,     0,   540,   519,     0,     0,     0,     0,     0,
     532,   533,     0,     0,     0,     0,     0,     0,   553,   528,
       0,     0,   539,   534,   518,     0,     0,   538,   536,     0,
     298,   334,   323,   299,   300,   301,   302,   303,   304,   305,
     306,   307,   308,   309,   310,   311,   312,   313,   314,   315,
     316,   317,   318,   319,   320,   321,   322,   324,   325,   326,
     327,   328,   329,   330,   331,   332,   333,   335,   336,   337,
     242,     0,   295,     0,   259,     0,     0,   256,     0,     0,
       0,     0,     0,   278,     0,     0,     0,     0,   272,     0,
       0,   119,     0,     0,     0,     0,    83,     0,   493,     0,
       0,     0,     0,   547,   546,     0,   402,   403,   404,   405,
       0,     0,     0,   171,    88,    89,    87,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   364,     0,     0,     0,     0,     0,     0,
       0,     0,   498,   499,   500,     0,     0,     0,   541,     0,
     507,     0,     0,     0,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   231,   233,   234,   235,   236,
     237,   238,   239,   230,   232,   240,   241,   375,   372,    78,
      73,     0,    54,     0,    79,   146,   147,   166,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   418,   145,     0,   341,   346,   342,   343,   344,   345,
     347,   348,   349,   350,   351,   352,   353,   354,     0,    50,
       0,     0,     0,   210,   211,   212,     0,   200,   201,     0,
     218,   215,     0,   426,     0,   425,   427,   422,   366,    60,
      55,     0,    51,    66,    61,     0,    52,    72,    67,     0,
      53,   361,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     460,   461,   462,   463,   464,   465,   466,   467,   468,   469,
     470,   471,   472,   473,   474,   475,   476,   477,   478,   487,
     479,   480,   481,   482,   483,   484,   485,   486,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   377,   378,     0,     0,     0,    82,    49,
       0,     0,     0,     0,     0,     0,     0,   113,   114,   247,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   244,
       0,   258,   254,   255,   287,     0,   287,     0,   276,   277,
       0,   287,     0,   270,   271,     0,   117,   118,   110,     0,
       0,    84,     0,     0,   140,     0,   141,     0,   136,     0,
       0,     0,     0,     0,     0,     0,   169,   170,     0,    95,
      96,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,     0,   362,
     363,     0,     0,   367,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    76,    74,
      80,     0,   153,   154,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   168,   195,   196,     0,     0,   187,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    58,
      56,    64,    62,    70,    68,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   489,   488,   556,   253,     0,     0,
     557,   558,   559,   555,   561,   510,   513,   512,     0,     0,
     511,   514,   515,   548,     0,   549,   445,     0,   562,   520,
     537,   453,     0,   387,     0,   383,     0,     0,     0,   491,
     204,   203,   399,   394,     0,     0,   393,   388,     0,     0,
       0,   550,   501,   542,   543,   516,   517,   522,   525,   523,
     530,   531,   521,   527,   526,     0,   529,   297,   294,     0,
     243,     0,     0,   282,   289,   283,   288,   285,   290,   284,
     286,     0,     0,     0,   264,     0,     0,   287,     0,     0,
     287,   250,     0,     0,     0,   112,     0,     0,   129,   138,
     139,     0,   143,   124,   123,     0,   125,     0,   122,   126,
     127,     0,   132,   130,   544,   545,   401,   412,   414,   415,
     416,   413,     0,   406,   410,     0,     0,     0,     0,   108,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    86,    91,    90,    92,    93,    94,     0,     0,   370,
       0,     0,   497,   505,   490,   496,   502,   503,   494,   504,
     509,   495,   492,   508,   374,     0,    77,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   144,   149,
     148,   150,   151,   152,   419,   340,   338,     0,   355,     0,
       0,     0,     0,   192,   193,     0,     0,   209,   208,   199,
     198,   217,   214,     0,   554,   424,   421,     0,    59,    65,
      71,   563,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,   578,   579,   580,   581,
     582,   583,   584,   585,   586,   587,   588,   589,   458,   459,
     252,   251,   591,   593,   595,   594,     0,   448,   456,     0,
       0,   385,   384,     0,   395,     0,   389,     0,   116,     0,
     359,     0,   296,   245,   260,   292,   291,   257,   287,   287,
       0,   287,     0,     0,   275,     0,   249,   248,     0,     0,
       0,     0,     0,     0,   134,     0,     0,     0,     0,   400,
     287,   411,     0,     0,     0,     0,     0,     0,   107,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,     0,
       0,     0,   368,   376,   167,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   164,   339,   356,   194,   186,   185,
     189,   190,     0,     0,   216,   423,     0,   590,   383,     0,
     380,   397,     0,   391,     0,     0,     0,   524,   261,     0,
       0,     0,   287,     0,   287,   287,   273,     0,   111,     0,
     142,   506,   121,     0,     0,     0,     0,   407,     0,     0,
     174,     0,   182,     0,     0,   109,   365,   371,     0,     0,
     191,   592,     0,     0,   398,   392,     0,   360,   287,   287,
     287,     0,   281,     0,     0,     0,   165,     0,   137,   133,
     131,   287,   408,     0,     0,     0,   177,     0,     0,   173,
     369,   188,     0,   381,   287,   266,   262,   265,   287,   279,
     274,   120,     0,     0,   176,   175,   181,     0,   179,     0,
       0,     0,   358,   287,     0,     0,   135,   409,   178,     0,
     184,     0,   382,     0,   267,     0,   280,   180,   172,   183,
     357,   287,   287,   268,   263,   287,   269
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1051, -1051,  1467, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,  -314, -1051,
   -1051, -1051, -1051,  -111,  -225, -1051, -1051,  1206, -1051,   494,
   -1051, -1051, -1051, -1051, -1051, -1051,  -911,  -597,  -132,  -594,
   -1051, -1051, -1051,  1390,  -229, -1051, -1051, -1051, -1051,   589,
   -1051, -1051,   829, -1051, -1051,  1000, -1051, -1051,   853, -1051,
   -1051,  -121,   -24,   863,  1022, -1051, -1051,  1264, -1051, -1051,
   -1050, -1051, -1051,  1252, -1051, -1051,  1256,  -987,  -568, -1051,
   -1051,   981, -1051,  1429,   870, -1051,   476, -1051, -1051, -1051,
   -1051,  1212, -1051, -1051, -1051, -1051, -1051, -1051, -1051,  1364,
    -745, -1051, -1051, -1051, -1051, -1051,   947, -1051,   538,  -842,
   -1051, -1051, -1051, -1051, -1051,   861, -1051,   -63,  1029, -1051,
   -1051,  1025, -1051, -1051,   830, -1051,  -484, -1051,   -90, -1051,
    1462, -1051, -1051, -1051, -1051, -1051, -1051, -1051,   -94, -1051,
   -1051,  -131,  -618, -1051, -1051, -1051, -1051,  -104,   -98,   -92,
     -82,   -55, -1051, -1051,   -89,   -70, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051,   -54, -1051, -1051, -1051, -1051, -1051,
     -52,   -51,   -65,   -50,   -47,   -45, -1051, -1051, -1051, -1051,
    -105,  -100,   -87,   -84,   -28,   -27,   -25, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051, -1051,
   -1051, -1051, -1051, -1051, -1051,   815, -1051,  -545
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    44,    45,    46,    47,    48,    49,    50,    51,    52,
     152,   154,   156,   131,    53,    54,    55,    56,   360,   880,
      57,   324,    58,   228,   229,    59,   320,   321,   857,   858,
      60,   327,  1037,  1036,  1113,   861,   623,   624,   625,   626,
     432,    61,    62,   342,   343,  1123,  1191,    63,   708,   709,
      64,   456,   457,    65,   214,   215,    66,   452,   453,    67,
     459,   463,   110,   844,   760,    68,   306,   307,   308,   832,
    1098,    69,   317,   318,    70,   312,   313,   833,  1099,    71,
     259,   260,    72,   433,   434,    73,  1010,  1011,    74,    75,
     362,   363,    76,    77,   365,    78,    79,    80,   211,   212,
     562,    81,    82,    83,    84,   335,   336,   872,   873,   874,
      85,   134,   700,    86,   464,   465,   179,   180,   181,    87,
     203,   204,    88,    89,   509,   510,   756,   384,   385,   386,
     387,   388,   389,   390,   391,   392,   393,   394,   395,   396,
     397,   398,   399,   860,   400,   401,   402,   182,   183,   184,
     185,   186,   268,   269,   403,   437,   272,   273,   274,   275,
     276,   277,   278,   279,   438,   281,   282,   283,   284,   285,
     439,   440,   441,   442,   443,   444,   404,   292,   293,   337,
     405,   406,   187,   188,   447,   189,   190,   299,   466,   191,
     192,   193,   194,   195,   196,   197,   207,   511,   512,   513,
     514,   515,   516,   517,   518,   519,   520,   521,   522,   523,
     524,   525,   526,   527,   528,   529,   530,   531,   532,   533,
     534,   535,   536,   537,   775,   993,   769,   770
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       128,   129,   431,   578,   322,   263,   383,   137,   216,   859,
     454,   264,   146,   149,   150,   262,   338,   265,   157,   261,
     270,   339,   294,   205,   460,   295,   849,   266,   206,   850,
     202,  1041,  1100,   639,   640,   801,   834,   455,   836,   271,
    1000,   435,   435,   839,   288,   651,   436,   436,   445,   445,
     662,   446,   446,   458,   267,   280,   761,   286,   287,   289,
     848,   876,   290,   943,   291,   992,   705,  1080,    90,   767,
     779,   415,   944,   867,   824,   706,  1081,    92,   209,   102,
     416,   296,   297,   826,   298,   107,   209,   171,  1148,  1166,
     328,   219,   560,   566,   367,   577,  1156,   578,   596,   415,
     104,   132,   807,   808,   809,   563,  1038,   864,   416,   417,
     616,   816,   828,   637,   221,   853,    94,    95,    96,   133,
     851,   107,   222,   867,  1114,  1115,  1116,   854,   107,   865,
     300,   868,   571,   855,   107,   123,   574,   417,   572,   107,
     340,   107,   575,   107,   227,   107,   101,   329,   561,   636,
     567,   856,  1204,   121,   210,   597,   945,   330,   107,   319,
     376,  1039,   210,   877,   300,   418,   419,   106,   831,   111,
     107,   420,   421,   422,   423,   424,   425,   426,   427,   428,
    1157,   868,   903,   107,  1040,   301,   429,   869,   867,   910,
     103,   870,   871,   418,   419,   108,   109,  1195,   979,   420,
     421,   422,   423,   424,   425,   426,   427,   428,    91,   107,
     994,   105,   946,   768,   429,   707,  1082,    93,   878,   301,
     430,   637,   622,   695,   696,   697,   698,   823,   699,  1181,
     220,   793,   954,   107,  1151,  1158,  1172,   869,   126,   127,
     302,   870,   871,  1149,   144,   145,   868,   341,   430,   797,
     622,   147,   148,   818,   415,   914,   976,   977,   300,  1150,
     980,   981,   800,   416,   984,   825,   986,   987,   936,  1022,
     415,   410,  1025,   827,   407,  1041,   340,   852,   107,   416,
     938,  1045,   300,   941,   107,   704,   682,   683,   224,   942,
     300,   340,   417,   952,    99,   117,   225,   300,   694,   112,
     754,  1046,   869,   100,   118,   786,   870,   871,   417,   755,
     113,   829,   300,   301,   787,   114,   757,    97,    98,   956,
     695,   696,   697,   698,   879,   699,   411,   881,   882,   883,
     884,   885,   886,   887,   888,   889,   890,   301,   892,   893,
     894,   895,   896,  1132,   897,   301,   830,   470,   418,   419,
     901,   823,   301,   474,   420,   421,   422,   423,   424,   425,
     426,   427,   428,   478,   418,   419,   303,   301,   408,   429,
     420,   421,   422,   423,   424,   425,   426,   427,   428,   300,
     841,   412,   309,   341,   300,   429,   300,   824,  1091,   825,
     300,   300,   449,   940,  1093,   300,   826,   827,   341,   618,
     461,   300,   471,   430,   122,   622,   300,   467,   475,   124,
     758,   759,   125,   300,   300,   859,   115,   116,   479,   430,
     300,   622,   468,   627,   304,   828,   632,   701,   701,   303,
     659,   710,   314,   130,   301,   829,   119,   120,   602,   301,
     310,   301,   712,   140,   141,   301,   301,   135,   227,   138,
     301,  1101,   332,  1103,   663,   136,   301,   472,   714,  1108,
     139,   301,   305,   476,   333,   842,   843,   151,   301,   301,
     830,   153,  1118,   480,  1003,   301,   216,   334,   311,   361,
     628,   831,   717,   633,   702,   703,   205,   304,   711,   481,
     315,   206,   263,   202,  1004,   155,  1006,   208,   264,   713,
    1012,  1063,   262,   664,   265,  1075,   261,   270,   227,   294,
     162,  1076,   295,   213,   266,   715,  1084,   849,   849,   849,
     850,   850,   850,  1085,  1134,   305,   271,   338,   316,   309,
    1135,   288,   339,   198,  1141,  1033,  1143,  1144,   608,   718,
    1047,   267,   280,   217,   286,   287,   289,   794,  1111,   290,
     798,   291,   918,   919,   920,   921,   922,   923,   924,   925,
     926,   927,  1043,   929,   930,   931,   932,   933,   296,   297,
    1165,   298,  1167,   819,  1096,   849,   314,  1015,   850,   454,
     142,   143,   824,  1173,  1060,   613,  1016,   310,   158,   159,
    1186,   826,  1034,   951,   218,   223,  1182,  1048,   230,   226,
    1185,   415,   654,   655,   227,   656,   455,   435,   697,   698,
     416,   699,   436,   319,   445,  1194,   323,   446,   325,   326,
     828,   361,   458,   364,   409,   311,   413,   414,   451,   538,
     469,  1097,  1083,  1203,   315,   539,   415,  1206,   473,   417,
     477,   652,   653,   654,   655,   416,   656,   652,   653,   654,
     655,   552,   656,   915,   695,   696,   697,   698,   540,   699,
     541,   652,   653,   654,   655,   542,   656,   543,   652,   653,
     654,   655,   316,   656,   417,   544,   831,   545,   937,   939,
     652,   653,   654,   655,   546,   656,   547,   548,   549,   550,
     551,   953,   553,   554,   957,   418,   419,   891,   555,   556,
     557,   420,   421,   422,   423,   424,   425,   426,   427,   428,
     928,   559,   558,   565,   576,   570,   429,  1049,   568,  1031,
     652,   653,   654,   655,  1050,   656,   569,   573,   579,   580,
     418,   419,  1124,  1125,   657,  1029,   420,   421,   422,   423,
     424,   425,   426,   427,   428,   581,   582,  1128,   583,   584,
     430,   429,   622,   585,   652,   653,   654,   655,   586,   656,
     587,   588,   160,   589,   590,   591,   592,   593,  1129,     1,
       2,   594,   595,   598,   599,  1133,  1051,   600,   601,     3,
       4,     5,   605,   607,   604,   430,     6,   622,   610,   612,
       7,     8,     9,  1145,    10,   615,    11,    12,    13,    14,
     619,   606,   611,   620,   578,  1153,   621,  1154,   631,    15,
    1052,   629,    16,   652,   653,   654,   655,   634,   656,   652,
     653,   654,   655,   635,   656,    17,   665,   652,   653,   654,
     655,   630,   656,   666,   638,   667,    18,    19,    20,   668,
     344,   641,    21,   642,   669,   670,   671,   672,  1180,   345,
     673,   674,   643,    22,   644,    23,   675,    24,    25,    26,
      27,    28,   645,   646,   647,  1190,    29,    30,   676,  1053,
     648,    31,    32,    33,    34,  1054,   649,  1199,   346,   650,
     658,    35,    36,  1055,    37,   677,   681,   678,    38,   679,
     680,    39,    40,    41,    42,   684,     1,     2,  1109,   716,
     652,   653,   654,   655,   719,   656,     3,     4,     5,   685,
     686,   720,   687,     6,   688,   689,   721,     7,     8,     9,
      43,    10,   722,    11,    12,    13,    14,   690,   723,   344,
     724,   725,   691,   661,   347,   348,    15,   692,   345,    16,
     349,   350,   351,   352,   353,   354,   355,   356,   357,   693,
     726,   727,    17,   728,   729,   358,  1056,   344,   730,   731,
     732,   733,   734,    18,    19,    20,   345,   346,   735,    21,
     652,   653,   654,   655,   736,   656,   752,   737,   738,   739,
      22,  1092,    23,  1094,    24,    25,    26,    27,    28,   359,
     740,   741,   742,    29,    30,   346,   743,   744,    31,    32,
      33,    34,   745,   746,   747,   748,   749,   344,    35,    36,
     750,    37,   751,   753,   765,    38,   345,   762,    39,    40,
      41,    42,   900,   347,   348,   764,  1057,   766,   773,   349,
     350,   351,   352,   353,   354,   355,   356,   357,   771,   415,
     774,   772,   777,   778,   358,   346,   780,    43,   416,   782,
     783,   347,   348,   784,   788,   785,   789,   349,   350,   351,
     352,   353,   354,   355,   356,   357,   791,   792,   652,   653,
     654,   655,   358,   656,   795,   796,   799,   417,   359,   652,
     653,   654,   655,   768,   656,   802,   803,   804,   163,   164,
     165,   166,   167,   168,   169,   805,   806,   810,   811,   812,
     170,   347,   348,   813,   171,  1198,   359,   349,   350,   351,
     352,   353,   354,   355,   356,   357,   652,   653,   654,   655,
     172,   656,   358,   814,  1058,   695,   696,   697,   698,   815,
     699,   820,   821,   418,   419,  1059,   822,   231,   835,   420,
     421,   422,   423,   424,   425,   426,   427,   428,   837,   838,
     840,   845,   200,   170,   429,   847,   359,   171,   846,   875,
     656,   173,   174,   898,   695,   696,   697,   698,   862,   699,
     899,   232,   233,   172,   863,   201,   902,   904,   234,   175,
     176,  1065,   905,   916,   906,   235,   236,   237,   430,   934,
     238,   239,   907,   240,   241,   908,   242,   243,   244,   245,
     246,   247,   248,   249,   250,   251,   909,   911,   231,   699,
     912,   252,   177,   178,   173,   174,   913,   253,   958,   254,
    1066,   948,   950,   200,   255,   695,   696,   697,   698,   959,
     699,   960,   175,   176,   961,   256,   695,   696,   697,   698,
     962,   699,   232,   233,   963,   964,   201,   257,   213,   234,
     965,   966,   967,   258,   988,   968,   235,   163,   164,   165,
     166,   167,   168,   169,   199,   177,   178,   969,   200,   170,
     970,   971,   972,   171,   652,   653,   654,   655,   973,   656,
     974,  1067,   252,   695,   696,   697,   698,   990,   699,   172,
     254,   201,  1068,   975,   695,   696,   697,   698,   978,   699,
     982,   991,   695,   696,   697,   698,   256,   699,   983,   985,
     996,   366,   992,   997,   695,   696,   697,   698,   257,   699,
     998,   999,  1018,  1013,   258,   561,  1009,  1001,  1008,  1002,
     173,   174,   367,  1014,   368,   369,   177,   178,  1005,  1069,
    1007,  1017,  1019,  1020,  1021,  1023,  1026,  1024,   175,   176,
    1070,  1027,  1028,  1030,   234,   366,   370,   371,  1071,  1032,
    1035,   235,   695,   696,   697,   698,  1042,   699,   328,  1044,
    1072,    -1,   695,   696,   697,   698,   367,   699,   368,   369,
    1064,   177,   178,  1061,  1086,  1079,   652,   653,   654,   655,
    1087,   656,   372,  1089,   373,   254,   374,   333,   234,  1095,
     370,   371,   375,  1088,  1102,   235,   376,  1104,  1105,  1107,
     334,  1119,   328,  1106,   377,   378,   379,  1121,  1073,  1122,
     380,   381,   382,  1126,   213,  1131,  1138,  1139,  1074,  1136,
    1142,   462,   652,   653,   654,   655,   372,   656,   373,   254,
     374,   333,  1077,  1152,  1140,  1147,   375,  1162,  1155,  1164,
     376,  1168,  1169,  1176,   334,  1178,  1170,  1179,   377,   378,
     379,  1177,  1183,  1184,   380,   381,   382,  1187,   213,   652,
     653,   654,   655,  1188,   656,   652,   653,   654,   655,  1189,
     656,   652,   653,   654,   655,  1193,   656,  1196,  1130,   695,
     696,   697,   698,  1197,   699,   652,   653,   654,   655,  1200,
     656,   652,   653,   654,   655,  1201,   656,   695,   696,   697,
     698,   161,   699,   652,   653,   654,   655,  1202,   656,   652,
     653,   654,   655,  1205,   656,  1161,   617,  1112,   450,  1062,
    1078,   652,   653,   654,   655,  1090,   656,   652,   653,   654,
     655,   949,   656,  1110,   917,   652,   653,   654,   655,  1120,
     656,   652,   653,   654,   655,  1127,   656,   652,   653,   654,
     655,  1146,   656,   947,   763,   448,   790,  1159,   609,   614,
     603,   935,  1137,  1160,   660,   564,  1117,   817,   955,   866,
     776,   781,   331,   989,   995,  1163,     0,     0,     0,     0,
       0,  1171,     0,     0,     0,     0,     0,     0,     0,  1174,
       0,     0,     0,     0,     0,  1175,     0,     0,     0,     0,
       0,  1192,   482,   483,   484,   485,   486,   487,   488,   489,
     490,   491,   492,   493,   494,   495,   496,   497,   498,   499,
     500,   501,   502,     0,     0,     0,   503,   504,     0,   505,
     506,   507,   508
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   134,   228,   115,   109,   127,    31,    98,   627,
     141,   109,    36,    37,    38,   109,   121,   109,    42,   109,
     109,   121,   109,    93,   145,   109,   623,   109,    93,   623,
      93,   873,  1019,   347,   348,   580,   604,   141,   606,   109,
     785,   135,   136,   611,   109,   359,   135,   136,   135,   136,
     364,   135,   136,   143,   109,   109,   540,   109,   109,   109,
      31,    33,   109,    42,   109,    51,    42,    42,    51,    51,
     554,    42,    51,    80,    42,    51,    51,    51,     4,    80,
      51,   109,   109,    51,   109,    80,     4,    25,    31,  1139,
      60,    51,   135,   135,    24,    31,    51,   322,   135,    42,
      80,    33,   586,   587,   588,    31,   135,    56,    51,    80,
      31,   595,    80,   342,   190,    45,   190,   191,   190,    51,
      33,    80,   198,    80,  1035,  1036,  1037,    57,    80,    78,
      80,   138,   192,    63,    80,    33,   192,    80,   198,    80,
      22,    80,   198,    80,    80,    80,    51,   117,   191,    31,
     192,    81,  1202,   191,    80,   192,   135,   127,    80,    80,
      98,   190,    80,   135,    80,   136,   137,   190,   136,   190,
      80,   142,   143,   144,   145,   146,   147,   148,   149,   150,
     135,   138,   666,    80,   191,   135,   157,   194,    80,   673,
     191,   198,   199,   136,   137,   190,   191,  1184,   743,   142,
     143,   144,   145,   146,   147,   148,   149,   150,   191,    80,
     196,   191,   191,   195,   157,   191,   191,   191,   190,   135,
     191,   450,   193,   136,   137,   138,   139,     6,   141,   197,
     190,   190,   716,    80,   191,   190,  1147,   194,   190,   191,
     190,   198,   199,    31,   190,   191,   138,   129,   191,   190,
     193,   190,   191,   190,    42,   190,   740,   741,    80,    31,
     744,   745,   576,    51,   748,    44,   750,   751,   190,   837,
      42,    80,   840,    52,   190,  1117,    22,   190,    80,    51,
     190,    98,    80,   190,    80,    31,   418,   419,   190,   196,
      80,    22,    80,   190,    33,    33,   198,    80,   430,   190,
      42,   118,   194,    42,    42,    42,   198,   199,    80,    51,
     190,    90,    80,   135,    51,   190,    42,   190,   191,   190,
     136,   137,   138,   139,   638,   141,   135,   641,   642,   643,
     644,   645,   646,   647,   648,   649,   650,   135,   652,   653,
     654,   655,   656,  1088,   658,   135,   125,    80,   136,   137,
     664,     6,   135,    80,   142,   143,   144,   145,   146,   147,
     148,   149,   150,    80,   136,   137,    22,   135,   190,   157,
     142,   143,   144,   145,   146,   147,   148,   149,   150,    80,
      51,   190,    22,   129,    80,   157,    80,    42,   190,    44,
      80,    80,   190,   707,   190,    80,    51,    52,   129,    31,
     190,    80,   135,   191,   190,   193,    80,   190,   135,   190,
     136,   137,   190,    80,    80,  1033,   190,   191,   135,   191,
      80,   193,   190,   135,    80,    80,   135,   135,   135,    22,
      31,   135,    22,    80,   135,    90,   190,   191,    31,   135,
      80,   135,   135,   190,   191,   135,   135,   191,    80,   190,
     135,  1019,    77,  1021,    31,   191,   135,   190,   135,    31,
      33,   135,   118,   190,    89,   136,   137,    80,   135,   135,
     125,    80,  1040,   190,   788,   135,   566,   102,   118,    80,
     192,   136,   135,   192,   192,   192,   556,    80,   192,   190,
      80,   556,   596,   556,   190,    80,   190,    23,   596,   192,
     190,   190,   596,    80,   596,   190,   596,   596,    80,   596,
     190,   190,   596,   116,   596,   192,   190,  1114,  1115,  1116,
    1114,  1115,  1116,   190,   190,   118,   596,   632,   118,    22,
     190,   596,   632,   190,  1102,   135,  1104,  1105,    31,   192,
     135,   596,   596,    42,   596,   596,   596,   571,  1032,   596,
     574,   596,   684,   685,   686,   687,   688,   689,   690,   691,
     692,   693,   876,   695,   696,   697,   698,   699,   596,   596,
    1138,   596,  1140,   597,   135,  1172,    22,    42,  1172,   710,
     190,   191,    42,  1151,   898,    31,    51,    80,   190,   191,
      31,    51,   192,   714,   190,    80,  1164,   192,    33,    80,
    1168,    42,   138,   139,    80,   141,   710,   701,   138,   139,
      51,   141,   701,    80,   701,  1183,    38,   701,    42,   190,
      80,    80,   712,    80,   124,   118,    51,   190,   195,    33,
     124,   192,   946,  1201,    80,    33,    42,  1205,   124,    80,
     124,   136,   137,   138,   139,    51,   141,   136,   137,   138,
     139,   192,   141,   677,   136,   137,   138,   139,    33,   141,
      33,   136,   137,   138,   139,    33,   141,    33,   136,   137,
     138,   139,   118,   141,    80,    33,   136,    33,   702,   703,
     136,   137,   138,   139,    33,   141,    33,    33,    33,    33,
     135,   715,    33,    33,   718,   136,   137,   192,    33,   135,
     192,   142,   143,   144,   145,   146,   147,   148,   149,   150,
     192,    80,   192,    33,    33,    80,   157,   192,   190,   851,
     136,   137,   138,   139,   192,   141,   190,    80,    80,    33,
     136,   137,  1046,  1047,   190,   846,   142,   143,   144,   145,
     146,   147,   148,   149,   150,    33,    33,  1061,    33,    33,
     191,   157,   193,    33,   136,   137,   138,   139,    33,   141,
      33,    33,     0,    33,    33,    33,    33,    33,  1082,     7,
       8,    33,    33,    33,    80,  1089,   192,    80,    80,    17,
      18,    19,    80,    80,   135,   191,    24,   193,    80,    80,
      28,    29,    30,  1107,    32,   191,    34,    35,    36,    37,
      33,   135,   135,   192,  1029,  1119,   190,  1121,    33,    47,
     192,   192,    50,   136,   137,   138,   139,    80,   141,   136,
     137,   138,   139,    80,   141,    63,    33,   136,   137,   138,
     139,   192,   141,    33,   191,    33,    74,    75,    76,    33,
      42,   191,    80,   191,    33,    33,    33,    33,  1162,    51,
      33,    33,   191,    91,   191,    93,    33,    95,    96,    97,
      98,    99,   191,   191,   191,  1179,   104,   105,    33,   192,
     191,   109,   110,   111,   112,   192,   191,  1191,    80,   191,
     191,   119,   120,   192,   122,   192,   191,   124,   126,    80,
     190,   129,   130,   131,   132,   191,     7,     8,  1030,    33,
     136,   137,   138,   139,   124,   141,    17,    18,    19,   191,
     191,    80,   191,    24,   191,   191,   124,    28,    29,    30,
     158,    32,    80,    34,    35,    36,    37,   191,   124,    42,
      80,    33,   191,   135,   136,   137,    47,   191,    51,    50,
     142,   143,   144,   145,   146,   147,   148,   149,   150,   191,
      33,    33,    63,    33,    33,   157,   192,    42,    33,    33,
      33,    33,    33,    74,    75,    76,    51,    80,    33,    80,
     136,   137,   138,   139,    33,   141,   192,    33,    33,    33,
      91,  1005,    93,  1007,    95,    96,    97,    98,    99,   191,
      33,    33,    33,   104,   105,    80,    33,    33,   109,   110,
     111,   112,    33,    33,    33,    33,    33,    42,   119,   120,
      33,   122,    33,   135,    80,   126,    51,    51,   129,   130,
     131,   132,   135,   136,   137,    51,   192,    51,    80,   142,
     143,   144,   145,   146,   147,   148,   149,   150,    51,    42,
      80,    51,    51,    51,   157,    80,    51,   158,    51,    51,
     190,   136,   137,   191,    33,    80,    51,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   190,   190,   136,   137,
     138,   139,   157,   141,    80,   190,    80,    80,   191,   136,
     137,   138,   139,   195,   141,    51,    51,    51,     9,    10,
      11,    12,    13,    14,    15,    51,    51,    51,    51,    51,
      21,   136,   137,    51,    25,   190,   191,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   136,   137,   138,   139,
      41,   141,   157,    80,   192,   136,   137,   138,   139,   191,
     141,    80,    33,   136,   137,   192,   135,     5,   135,   142,
     143,   144,   145,   146,   147,   148,   149,   150,   135,   135,
     135,    80,    20,    21,   157,    80,   191,    25,   190,   135,
     141,    82,    83,    80,   136,   137,   138,   139,   190,   141,
     190,    39,    40,    41,   190,    43,    51,    51,    46,   100,
     101,   192,    51,   124,    51,    53,    54,    55,   191,   190,
      58,    59,    51,    61,    62,    51,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    51,    51,     5,   141,
      51,    79,   133,   134,    82,    83,    51,    85,   124,    87,
     192,   190,   190,    20,    92,   136,   137,   138,   139,   124,
     141,   124,   100,   101,    51,   103,   136,   137,   138,   139,
      51,   141,    39,    40,    51,    51,    43,   115,   116,    46,
      51,    51,    51,   121,   190,    51,    53,     9,    10,    11,
      12,    13,    14,    15,    16,   133,   134,    51,    20,    21,
      51,    51,    51,    25,   136,   137,   138,   139,    51,   141,
      51,   192,    79,   136,   137,   138,   139,    42,   141,    41,
      87,    43,   192,    51,   136,   137,   138,   139,    51,   141,
      51,    42,   136,   137,   138,   139,   103,   141,    51,    51,
     194,     3,    51,   190,   136,   137,   138,   139,   115,   141,
     190,    51,   135,    80,   121,   191,   197,   192,   190,   192,
      82,    83,    24,    80,    26,    27,   133,   134,   192,   192,
     192,   190,   135,    80,   135,    80,    51,   190,   100,   101,
     192,    51,   192,    33,    46,     3,    48,    49,   192,    33,
     190,    53,   136,   137,   138,   139,    80,   141,    60,    80,
     192,   141,   136,   137,   138,   139,    24,   141,    26,    27,
     192,   133,   134,    80,   194,   190,   136,   137,   138,   139,
      80,   141,    84,    33,    86,    87,    88,    89,    46,    80,
      48,    49,    94,   192,   135,    53,    98,   135,   135,    33,
     102,    33,    60,   190,   106,   107,   108,    33,   192,    51,
     112,   113,   114,   190,   116,    51,   135,   135,   192,   197,
     190,   123,   136,   137,   138,   139,    84,   141,    86,    87,
      88,    89,   192,   192,   135,   190,    94,    33,   194,   135,
      98,   135,   190,    51,   102,    51,   190,   128,   106,   107,
     108,   194,   135,   135,   112,   113,   114,   192,   116,   136,
     137,   138,   139,    51,   141,   136,   137,   138,   139,   194,
     141,   136,   137,   138,   139,    80,   141,   190,   192,   136,
     137,   138,   139,    51,   141,   136,   137,   138,   139,   197,
     141,   136,   137,   138,   139,   135,   141,   136,   137,   138,
     139,    44,   141,   136,   137,   138,   139,   135,   141,   136,
     137,   138,   139,   135,   141,   192,   320,  1033,   138,   190,
     941,   136,   137,   138,   139,   190,   141,   136,   137,   138,
     139,   712,   141,   190,   681,   136,   137,   138,   139,   190,
     141,   136,   137,   138,   139,   190,   141,   136,   137,   138,
     139,   190,   141,   710,   542,   136,   566,   190,   312,   317,
     306,   701,  1096,   190,   362,   211,  1038,   596,   717,   632,
     551,   556,   120,   753,   769,   190,    -1,    -1,    -1,    -1,
      -1,   190,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   190,
      -1,    -1,    -1,    -1,    -1,   190,    -1,    -1,    -1,    -1,
      -1,   190,   159,   160,   161,   162,   163,   164,   165,   166,
     167,   168,   169,   170,   171,   172,   173,   174,   175,   176,
     177,   178,   179,    -1,    -1,    -1,   183,   184,    -1,   186,
     187,   188,   189
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,    17,    18,    19,    24,    28,    29,    30,
      32,    34,    35,    36,    37,    47,    50,    63,    74,    75,
      76,    80,    91,    93,    95,    96,    97,    98,    99,   104,
     105,   109,   110,   111,   112,   119,   120,   122,   126,   129,
     130,   131,   132,   158,   201,   202,   203,   204,   205,   206,
     207,   208,   209,   214,   215,   216,   217,   220,   222,   225,
     230,   241,   242,   247,   250,   253,   256,   259,   265,   271,
     274,   279,   282,   285,   288,   289,   292,   293,   295,   296,
     297,   301,   302,   303,   304,   310,   313,   319,   322,   323,
      51,   191,    51,   191,   190,   191,   190,   190,   191,    33,
      42,    51,    80,   191,    80,   191,   190,    80,   190,   191,
     262,   190,   190,   190,   190,   190,   191,    33,    42,   190,
     191,   191,   190,    33,   190,   190,   190,   191,   262,   262,
      80,   213,    33,    51,   311,   191,   191,   262,   190,    33,
     190,   191,   190,   191,   190,   191,   262,   190,   191,   262,
     262,    80,   210,    80,   211,    80,   212,   262,   190,   191,
       0,   202,   190,     9,    10,    11,    12,    13,    14,    15,
      21,    25,    41,    82,    83,   100,   101,   133,   134,   316,
     317,   318,   347,   348,   349,   350,   351,   382,   383,   385,
     386,   389,   390,   391,   392,   393,   394,   395,   190,    16,
      20,    43,   317,   320,   321,   355,   372,   396,    23,     4,
      80,   298,   299,   116,   254,   255,   328,    42,   190,    51,
     190,   190,   198,    80,   190,   198,    80,    80,   223,   224,
      33,     5,    39,    40,    46,    53,    54,    55,    58,    59,
      61,    62,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    79,    85,    87,    92,   103,   115,   121,   280,
     281,   328,   338,   347,   348,   349,   350,   351,   352,   353,
     354,   355,   356,   357,   358,   359,   360,   361,   362,   363,
     364,   365,   366,   367,   368,   369,   370,   371,   372,   373,
     374,   375,   377,   378,   382,   383,   384,   385,   386,   387,
      80,   135,   190,    22,    80,   118,   266,   267,   268,    22,
      80,   118,   275,   276,    22,    80,   118,   272,   273,    80,
     226,   227,   223,    38,   221,    42,   190,   231,    60,   117,
     127,   330,    77,    89,   102,   305,   306,   379,   380,   381,
      22,   129,   243,   244,    42,    51,    80,   136,   137,   142,
     143,   144,   145,   146,   147,   148,   149,   150,   157,   191,
     218,    80,   290,   291,    80,   294,     3,    24,    26,    27,
      48,    49,    84,    86,    88,    94,    98,   106,   107,   108,
     112,   113,   114,   261,   327,   328,   329,   330,   331,   332,
     333,   334,   335,   336,   337,   338,   339,   340,   341,   342,
     344,   345,   346,   354,   376,   380,   381,   190,   190,   124,
      80,   135,   190,    51,   190,    42,    51,    80,   136,   137,
     142,   143,   144,   145,   146,   147,   148,   149,   150,   157,
     191,   238,   240,   283,   284,   338,   354,   355,   364,   370,
     371,   372,   373,   374,   375,   382,   383,   384,   283,   190,
     243,   195,   257,   258,   341,   347,   251,   252,   328,   260,
     261,   190,   123,   261,   314,   315,   388,   190,   190,   124,
      80,   135,   190,   124,    80,   135,   190,   124,    80,   135,
     190,   190,   159,   160,   161,   162,   163,   164,   165,   166,
     167,   168,   169,   170,   171,   172,   173,   174,   175,   176,
     177,   178,   179,   183,   184,   186,   187,   188,   189,   324,
     325,   397,   398,   399,   400,   401,   402,   403,   404,   405,
     406,   407,   408,   409,   410,   411,   412,   413,   414,   415,
     416,   417,   418,   419,   420,   421,   422,   423,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,   135,   192,    33,    33,    33,   135,   192,   192,    80,
     135,   191,   300,    31,   299,    33,   135,   192,   190,   190,
      80,   192,   198,    80,   192,   198,    33,    31,   224,    80,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,   135,   192,    33,    80,
      80,    80,    31,   267,   135,    80,   135,    80,    31,   276,
      80,   135,    80,    31,   273,   191,    31,   227,    31,    33,
     192,   190,   193,   236,   237,   238,   239,   135,   192,   192,
     192,    33,   135,   192,    80,    80,    31,   244,   191,   218,
     218,   191,   191,   191,   191,   191,   191,   191,   191,   191,
     191,   218,   136,   137,   138,   139,   141,   190,   191,    31,
     291,   135,   218,    31,    80,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,   192,   124,    80,
     190,   191,   238,   238,   191,   191,   191,   191,   191,   191,
     191,   191,   191,   191,   238,   136,   137,   138,   139,   141,
     312,   135,   192,   192,    31,    42,    51,   191,   248,   249,
     135,   192,   135,   192,   135,   192,    33,   135,   192,   124,
      80,   124,    80,   124,    80,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   192,   135,    42,    51,   326,    42,   136,   137,
     264,   326,    51,   264,    51,    80,    51,    51,   195,   426,
     427,    51,    51,    80,    80,   424,   318,    51,    51,   326,
      51,   321,    51,   190,   191,    80,    42,    51,    33,    51,
     255,   190,   190,   190,   262,    80,   190,   190,   262,    80,
     218,   427,    51,    51,    51,    51,    51,   326,   326,   326,
      51,    51,    51,    51,    80,   191,   326,   281,   190,   262,
      80,    33,   135,     6,    42,    44,    51,    52,    80,    90,
     125,   136,   269,   277,   278,   135,   278,   135,   135,   278,
     135,    51,   136,   137,   263,    80,   190,    80,    31,   237,
     239,    33,   190,    45,    57,    63,    81,   228,   229,   342,
     343,   235,   190,   190,    56,    78,   306,    80,   138,   194,
     198,   199,   307,   308,   309,   135,    33,   135,   190,   218,
     219,   218,   218,   218,   218,   218,   218,   218,   218,   218,
     218,   192,   218,   218,   218,   218,   218,   218,    80,   190,
     135,   218,    51,   326,    51,    51,    51,    51,    51,    51,
     326,    51,    51,    51,   190,   262,   124,   263,   238,   238,
     238,   238,   238,   238,   238,   238,   238,   238,   192,   238,
     238,   238,   238,   238,   190,   284,   190,   262,   190,   262,
     218,   190,   196,    42,    51,   135,   191,   258,   190,   252,
     190,   261,   190,   262,   326,   315,   190,   262,   124,   124,
     124,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    51,    51,    51,    51,   326,   326,    51,   427,
     326,   326,    51,    51,   326,    51,   326,   326,   190,   324,
      42,    42,    51,   425,   196,   425,   194,   190,   190,    51,
     300,   192,   192,   218,   190,   192,   190,   192,   190,   197,
     286,   287,   190,    80,    80,    42,    51,   190,   135,   135,
      80,   135,   278,    80,   190,   278,    51,    51,   192,   223,
      33,   238,    33,   135,   192,   190,   233,   232,   135,   190,
     191,   309,    80,   218,    80,    98,   118,   135,   192,   192,
     192,   192,   192,   192,   192,   192,   192,   192,   192,   192,
     218,    80,   190,   190,   192,   192,   192,   192,   192,   192,
     192,   192,   192,   192,   192,   190,   190,   192,   249,   190,
      42,    51,   191,   218,   190,   190,   194,    80,   192,    33,
     190,   190,   262,   190,   262,    80,   135,   192,   270,   278,
     277,   278,   135,   278,   135,   135,   190,    33,    31,   238,
     190,   326,   229,   234,   236,   236,   236,   308,   278,    33,
     190,    33,    51,   245,   218,   218,   190,   190,   218,   218,
     192,    51,   300,   218,   190,   190,   197,   286,   135,   135,
     135,   278,   190,   278,   278,   218,   190,   190,    31,    31,
      31,   191,   192,   218,   218,   194,    51,   135,   190,   190,
     190,   192,    33,   190,   135,   278,   270,   278,   135,   190,
     190,   190,   236,   278,   190,   190,    51,   194,    51,   128,
     218,   197,   278,   135,   135,   278,    31,   192,    51,   194,
     218,   246,   190,    80,   278,   277,   190,    51,   190,   218,
     197,   135,   135,   278,   270,   135,   278
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
      59,    40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   200,   201,   201,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   203,   203,   203,   203,   204,   204,
     205,   206,   207,   208,   209,   210,   210,   210,   210,   210,
     210,   211,   211,   211,   211,   211,   211,   212,   212,   212,
     212,   212,   212,   213,   213,   213,   213,   213,   213,   214,
     214,   215,   215,   216,   216,   217,   218,   218,   218,   218,
     218,   218,   218,   218,   218,   218,   218,   218,   218,   218,
     218,   218,   218,   218,   218,   218,   218,   218,   219,   219,
     220,   220,   221,   222,   223,   223,   224,   225,   226,   226,
     227,   228,   228,   229,   229,   229,   229,   229,   231,   230,
     232,   230,   233,   230,   234,   230,   235,   230,   236,   236,
     236,   236,   237,   237,   238,   238,   238,   238,   238,   238,
     238,   238,   238,   238,   238,   238,   238,   238,   238,   238,
     238,   238,   238,   238,   238,   239,   240,   240,   241,   242,
     243,   243,   244,   244,   244,   244,   244,   245,   245,   245,
     245,   245,   245,   246,   246,   247,   248,   248,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   250,   250,   251,
     251,   252,   253,   253,   254,   254,   255,   256,   256,   257,
     257,   258,   258,   259,   259,   259,   259,   260,   260,   261,
     261,   261,   261,   261,   261,   261,   261,   261,   261,   261,
     261,   261,   261,   261,   261,   261,   261,   261,   261,   261,
     261,   261,   262,   262,   262,   262,   262,   262,   263,   263,
     263,   264,   264,   264,   265,   266,   266,   267,   268,   268,
     268,   269,   269,   269,   269,   269,   270,   270,   270,   270,
     271,   272,   272,   273,   273,   273,   274,   275,   275,   276,
     276,   276,   277,   277,   277,   277,   277,   278,   278,   278,
     278,   278,   278,   279,   279,   279,   279,   280,   280,   281,
     281,   281,   281,   281,   281,   281,   281,   281,   281,   281,
     281,   281,   281,   281,   281,   281,   281,   281,   281,   281,
     281,   281,   281,   281,   281,   281,   281,   281,   281,   281,
     281,   281,   281,   281,   281,   281,   281,   281,   282,   282,
     283,   283,   284,   284,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   284,   285,   285,   286,   286,   287,
     287,   288,   289,   290,   290,   291,   292,   293,   294,   294,
     294,   294,   295,   296,   296,   296,   296,   297,   298,   298,
     299,   299,   299,   300,   300,   300,   301,   301,   302,   302,
     302,   302,   302,   302,   303,   303,   303,   303,   303,   303,
     304,   305,   305,   306,   306,   306,   307,   307,   307,   307,
     308,   308,   309,   309,   309,   309,   309,   311,   312,   310,
     313,   313,   313,   313,   314,   314,   315,   315,   316,   316,
     316,   316,   316,   316,   316,   317,   317,   317,   317,   317,
     317,   317,   317,   317,   317,   318,   318,   319,   319,   320,
     320,   320,   320,   321,   321,   322,   322,   323,   323,   324,
     324,   325,   325,   325,   325,   325,   325,   325,   325,   325,
     325,   325,   325,   325,   325,   325,   325,   325,   325,   325,
     325,   325,   325,   325,   325,   325,   325,   325,   326,   326,
     327,   328,   329,   330,   331,   332,   333,   334,   335,   336,
     337,   338,   339,   340,   341,   342,   343,   344,   345,   346,
     347,   348,   348,   349,   350,   351,   352,   353,   354,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   379,   380,   381,   382,   383,
     384,   385,   386,   387,   388,   389,   390,   391,   392,   393,
     394,   395,   396,   397,   398,   399,   400,   401,   402,   403,
     404,   405,   406,   407,   408,   409,   410,   411,   412,   413,
     414,   415,   416,   417,   418,   419,   420,   421,   422,   423,
     424,   425,   425,   426,   426,   427
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
       7,     3,     1,     1,     1,     1,     1,     1,     0,     5,
       0,     8,     0,     8,     0,    10,     0,     8,     2,     2,
       1,     1,     4,     2,     3,     1,     1,     1,     3,     3,
       3,     3,     3,     2,     2,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     5,     1,     4,     4,     4,
       2,     1,     9,     6,     5,     7,     7,     2,     4,     3,
       5,     3,     1,     2,     1,     6,     3,     1,     5,     3,
       3,     4,     2,     2,     3,     1,     1,     2,     5,     3,
       1,     1,     2,     5,     3,     1,     1,     2,     5,     3,
       1,     1,     1,     2,     5,     3,     6,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     2,     4,     3,     5,     1,     3,     2,     2,
       1,     2,     2,     1,     4,     2,     1,     4,     2,     1,
       4,     3,     5,     9,     1,     5,     3,     5,     7,     9,
       4,     2,     1,     5,     7,     4,     4,     2,     1,     7,
       9,     6,     1,     1,     1,     1,     1,     0,     1,     1,
       1,     2,     2,     2,     5,     3,     6,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     5,     6,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     5,     6,     7,     5,     1,
       3,     3,     4,     2,     1,     5,     3,     4,     4,     6,
       3,     5,     3,     2,     5,     3,     6,     4,     2,     1,
       5,     7,     9,     0,     3,     3,     2,     5,     5,     6,
       3,     7,     8,     5,     5,     6,     3,     7,     8,     5,
       6,     3,     1,     1,     1,     1,     1,     3,     4,     6,
       1,     2,     1,     1,     1,     1,     1,     0,     0,     5,
       2,     5,     3,     6,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     1,     3,     6,     1,
       1,     1,     1,     3,     1,     3,     6,     2,     5,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       3,     3,     3,     1,     3,     3,     3,     3,     1,     1,
       1,     3,     3,     3,     3,     3,     3,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     1,
       3,     3,     3,     3,     5,     3,     3,     3,     1,     3,
       3,     3,     1,     1,     1,     1,     1,     3,     1,     1,
       1,     1,     3,     3,     3,     3,     1,     1,     3,     3,
       3,     1,     1,     1,     3,     3,     3,     3,     3,     3,
       1,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     1,     3,     2,     2,     2
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
  "MODIFIEDHARMONICMEAN", "MOMENTS_VARENDO", "NAME", "NO_COMPILER", "NOBS",
  "NOCONSTANT", "NOCORR", "NODIAGNOSTIC", "NOFUNCTIONS", "NOGRAPH",
  "NOMOMENTS", "NOPRINT", "NORMAL_PDF", "OBSERVATION_TRENDS", "OPTIM",
  "OPTIM_WEIGHTS", "ORDER", "OSR", "OSR_PARAMS", "PARAMETERS", "PERIODS",
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
       201,     0,    -1,   202,    -1,   201,   202,    -1,   203,    -1,
     214,    -1,   215,    -1,   216,    -1,   230,    -1,   220,    -1,
     222,    -1,   225,    -1,   217,    -1,   241,    -1,   242,    -1,
     247,    -1,   250,    -1,   253,    -1,   256,    -1,   259,    -1,
     279,    -1,   282,    -1,   285,    -1,   265,    -1,   274,    -1,
     271,    -1,   288,    -1,   289,    -1,   292,    -1,   204,    -1,
     205,    -1,   293,    -1,   295,    -1,   296,    -1,   297,    -1,
     301,    -1,   302,    -1,   303,    -1,   304,    -1,   310,    -1,
     313,    -1,   319,    -1,   322,    -1,   323,    -1,   209,    -1,
     206,    -1,   207,    -1,   208,    -1,    28,    51,   190,    -1,
      28,    51,    51,   190,    -1,   109,   262,   190,    -1,   129,
     210,   190,    -1,   130,   211,   190,    -1,   131,   212,   190,
      -1,    97,   213,   190,    -1,   210,    80,    -1,   210,   135,
      80,    -1,    80,    -1,   210,    80,   124,    -1,   210,   135,
      80,   124,    -1,    80,   124,    -1,   211,    80,    -1,   211,
     135,    80,    -1,    80,    -1,   211,    80,   124,    -1,   211,
     135,    80,   124,    -1,    80,   124,    -1,   212,    80,    -1,
     212,   135,    80,    -1,    80,    -1,   212,    80,   124,    -1,
     212,   135,    80,   124,    -1,    80,   124,    -1,   213,    80,
      -1,   213,   135,    80,    -1,    80,    -1,   213,    80,   124,
      -1,   213,   135,    80,   124,    -1,    80,   124,    -1,    98,
      51,   190,    -1,    98,    33,    51,   190,    -1,    24,    42,
     190,    -1,    24,    33,    42,   190,    -1,    63,    42,   190,
      -1,    63,    33,    42,   190,    -1,    80,    33,   218,   190,
      -1,   191,   218,   192,    -1,    80,    -1,    42,    -1,    51,
      -1,   218,   137,   218,    -1,   218,   136,   218,    -1,   218,
     138,   218,    -1,   218,   139,   218,    -1,   218,   141,   218,
      -1,   136,   218,    -1,   137,   218,    -1,   142,   191,   218,
     192,    -1,   143,   191,   218,   192,    -1,   144,   191,   218,
     192,    -1,   145,   191,   218,   192,    -1,   146,   191,   218,
     192,    -1,   147,   191,   218,   192,    -1,   148,   191,   218,
     192,    -1,   149,   191,   218,   192,    -1,   150,   191,   218,
     192,    -1,   157,   191,   218,   192,    -1,    80,   191,   219,
     192,    -1,   218,    -1,   219,   135,   218,    -1,    50,   190,
     223,    31,    -1,    50,   191,   221,   192,   190,   223,    31,
      -1,    38,    33,    80,    -1,    32,   190,   223,    31,    -1,
     223,   224,    -1,   224,    -1,    80,    33,   218,   190,    -1,
      47,   190,   226,    31,    -1,   226,   227,    -1,   227,    -1,
      80,   191,   263,   192,    33,   218,   190,    -1,   228,   135,
     229,    -1,   229,    -1,    57,    -1,    45,    -1,    81,    -1,
     342,    -1,   343,    -1,    -1,    74,   190,   231,   236,    31,
      -1,    -1,    74,   191,   330,   192,   190,   232,   236,    31,
      -1,    -1,    74,   191,   127,   192,   190,   233,   236,    31,
      -1,    -1,    74,   191,   117,   135,   228,   192,   234,   190,
     236,    31,    -1,    -1,    74,   191,   117,   192,   235,   190,
     236,    31,    -1,   236,   237,    -1,   236,   239,    -1,   237,
      -1,   239,    -1,   238,    33,   238,   190,    -1,   238,   190,
      -1,   191,   238,   192,    -1,   240,    -1,    42,    -1,    51,
      -1,   238,   137,   238,    -1,   238,   136,   238,    -1,   238,
     138,   238,    -1,   238,   139,   238,    -1,   238,   141,   238,
      -1,   136,   238,    -1,   137,   238,    -1,   142,   191,   238,
     192,    -1,   143,   191,   238,   192,    -1,   144,   191,   238,
     192,    -1,   145,   191,   238,   192,    -1,   146,   191,   238,
     192,    -1,   147,   191,   238,   192,    -1,   148,   191,   238,
     192,    -1,   149,   191,   238,   192,    -1,   150,   191,   238,
     192,    -1,   157,   191,   238,   192,    -1,   193,    80,    33,
     238,   190,    -1,    80,    -1,    80,   191,   263,   192,    -1,
     110,   190,   243,    31,    -1,    76,   190,   243,    31,    -1,
     243,   244,    -1,   244,    -1,   129,    80,   190,    98,   245,
     190,   128,   246,   190,    -1,   129,    80,   190,   118,   218,
     190,    -1,   129,    80,    33,   218,   190,    -1,   129,    80,
     135,    80,    33,   218,   190,    -1,    22,    80,   135,    80,
      33,   218,   190,    -1,   245,    51,    -1,   245,    51,   194,
      51,    -1,   245,   135,    51,    -1,   245,   135,    51,   194,
      51,    -1,    51,   194,    51,    -1,    51,    -1,   246,   218,
      -1,   218,    -1,   111,    33,   195,   248,   196,   190,    -1,
     248,   190,   249,    -1,   249,    -1,   249,   135,   191,   218,
     192,    -1,   249,   135,    42,    -1,   249,   135,    51,    -1,
     249,   191,   218,   192,    -1,   249,    42,    -1,   249,    51,
      -1,   191,   218,   192,    -1,    42,    -1,    51,    -1,   119,
     190,    -1,   119,   191,   251,   192,   190,    -1,   251,   135,
     252,    -1,   252,    -1,   328,    -1,    19,   190,    -1,    19,
     191,   254,   192,   190,    -1,   254,   135,   255,    -1,   255,
      -1,   328,    -1,   112,   190,    -1,   112,   191,   257,   192,
     190,    -1,   257,   135,   258,    -1,   258,    -1,   341,    -1,
     347,    -1,   120,   190,    -1,   120,   191,   260,   192,   190,
      -1,   120,   262,   190,    -1,   120,   191,   260,   192,   262,
     190,    -1,   260,   135,   261,    -1,   261,    -1,   327,    -1,
     328,    -1,   329,    -1,   330,    -1,   331,    -1,   332,    -1,
     333,    -1,   334,    -1,   335,    -1,   336,    -1,   337,    -1,
     354,    -1,   338,    -1,   376,    -1,   339,    -1,   340,    -1,
     341,    -1,   342,    -1,   344,    -1,   345,    -1,   346,    -1,
     380,    -1,   381,    -1,   262,    80,    -1,   262,    80,    33,
      80,    -1,   262,   135,    80,    -1,   262,   135,    80,    33,
      80,    -1,    80,    -1,    80,    33,    80,    -1,   137,    51,
      -1,   136,    51,    -1,    51,    -1,   137,    42,    -1,   136,
      42,    -1,    42,    -1,    35,   190,   266,    31,    -1,   266,
     267,    -1,   267,    -1,   268,   135,   269,   190,    -1,   118,
      80,    -1,    80,    -1,    22,    80,   135,    80,    -1,   277,
     135,   270,    -1,   278,   135,   277,   135,   270,    -1,   278,
     135,   278,   135,   278,   135,   277,   135,   270,    -1,   278,
      -1,   278,   135,   278,   135,   278,    -1,   278,   135,   278,
      -1,   278,   135,   278,   135,   278,    -1,   278,   135,   278,
     135,   278,   135,   278,    -1,   278,   135,   278,   135,   278,
     135,   278,   135,   278,    -1,    37,   190,   272,    31,    -1,
     272,   273,    -1,   273,    -1,   118,    80,   135,   278,   190,
      -1,    22,    80,   135,    80,   135,   278,   190,    -1,    80,
     135,   278,   190,    -1,    36,   190,   275,    31,    -1,   275,
     276,    -1,   276,    -1,   118,    80,   135,   278,   135,   278,
     190,    -1,    22,    80,   135,    80,   135,   278,   135,   278,
     190,    -1,    80,   135,   278,   135,   278,   190,    -1,     6,
      -1,    44,    -1,    90,    -1,    52,    -1,   125,    -1,    -1,
      51,    -1,    42,    -1,    80,    -1,   136,    51,    -1,   136,
      42,    -1,    34,   190,    -1,    34,   191,   280,   192,   190,
      -1,    34,   262,   190,    -1,    34,   191,   280,   192,   262,
     190,    -1,   280,   135,   281,    -1,   281,    -1,   347,    -1,
     348,    -1,   349,    -1,   350,    -1,   351,    -1,   352,    -1,
     353,    -1,   354,    -1,   355,    -1,   356,    -1,   357,    -1,
     358,    -1,   359,    -1,   360,    -1,   361,    -1,   362,    -1,
     363,    -1,   364,    -1,   365,    -1,   366,    -1,   367,    -1,
     368,    -1,   369,    -1,   370,    -1,   338,    -1,   371,    -1,
     372,    -1,   373,    -1,   374,    -1,   375,    -1,   377,    -1,
     378,    -1,   382,    -1,   383,    -1,   384,    -1,   328,    -1,
     385,    -1,   386,    -1,   387,    -1,   104,   191,   283,   192,
     190,    -1,   104,   191,   283,   192,   262,   190,    -1,   283,
     135,   284,    -1,   284,    -1,   354,    -1,   355,    -1,   364,
      -1,   370,    -1,   338,    -1,   371,    -1,   372,    -1,   373,
      -1,   374,    -1,   375,    -1,   382,    -1,   383,    -1,   384,
      -1,   105,   191,   283,   192,   190,    -1,   105,   191,   283,
     192,   262,   190,    -1,   197,    80,   197,   135,   197,    80,
     197,    -1,   197,    80,   197,   135,   278,    -1,   286,    -1,
     287,   135,   286,    -1,   132,   262,   190,    -1,    91,   190,
     290,    31,    -1,   290,   291,    -1,   291,    -1,    80,   191,
     218,   192,   190,    -1,   126,   262,   190,    -1,    93,   190,
     294,    31,    -1,   294,    80,   218,   190,    -1,   294,    80,
     135,    80,   218,   190,    -1,    80,   218,   190,    -1,    80,
     135,    80,   218,   190,    -1,    96,   262,   190,    -1,    95,
     190,    -1,    95,   191,   261,   192,   190,    -1,    95,   262,
     190,    -1,    95,   191,   261,   192,   262,   190,    -1,    18,
     190,   298,    31,    -1,   298,   299,    -1,   299,    -1,    80,
     300,    33,   218,   190,    -1,    80,   135,    80,   300,    33,
     218,   190,    -1,     4,    80,   191,    51,   192,   300,    33,
     218,   190,    -1,    -1,   191,    51,   192,    -1,   191,    42,
     192,    -1,    17,   190,    -1,    17,   191,    23,   192,   190,
      -1,    30,   191,    80,   192,   190,    -1,    30,   191,    80,
     192,   262,   190,    -1,    30,    80,   190,    -1,    30,   191,
      80,   198,    80,   192,   190,    -1,    30,   191,    80,   198,
      80,   192,   262,   190,    -1,    30,    80,   198,    80,   190,
      -1,    29,   191,    80,   192,   190,    -1,    29,   191,    80,
     192,   262,   190,    -1,    29,    80,   190,    -1,    29,   191,
      80,   198,    80,   192,   190,    -1,    29,   191,    80,   198,
      80,   192,   262,   190,    -1,    29,    80,   198,    80,   190,
      -1,    75,   191,   305,   192,   307,   190,    -1,   305,   135,
     306,    -1,   306,    -1,   379,    -1,   380,    -1,   381,    -1,
     308,    -1,   307,   135,   308,    -1,   308,   191,   278,   192,
      -1,   307,   135,   308,   191,   278,   192,    -1,   309,    -1,
     308,   309,    -1,    80,    -1,   199,    -1,   138,    -1,   194,
      -1,   198,    -1,    -1,    -1,    99,   311,   238,   312,   190,
      -1,   122,   190,    -1,   122,   191,   314,   192,   190,    -1,
     122,   262,   190,    -1,   122,   191,   314,   192,   262,   190,
      -1,   314,   135,   315,    -1,   315,    -1,   261,    -1,   388,
      -1,   389,    -1,   390,    -1,   391,    -1,   392,    -1,   393,
      -1,   394,    -1,   395,    -1,   316,    -1,   347,    -1,   382,
      -1,   383,    -1,   349,    -1,   351,    -1,   348,    -1,   350,
      -1,   385,    -1,   386,    -1,   317,   135,   318,    -1,   317,
      -1,     7,    51,   190,    -1,     7,   191,   318,   192,    51,
     190,    -1,   317,    -1,   372,    -1,   355,    -1,   396,    -1,
     320,   135,   321,    -1,   320,    -1,     8,    51,   190,    -1,
       8,   191,   321,   192,    51,   190,    -1,   158,   190,    -1,
     158,   191,   324,   192,   190,    -1,   325,   135,   324,    -1,
     325,    -1,   397,    -1,   398,    -1,   399,    -1,   400,    -1,
     401,    -1,   402,    -1,   403,    -1,   404,    -1,   405,    -1,
     406,    -1,   407,    -1,   408,    -1,   409,    -1,   410,    -1,
     411,    -1,   412,    -1,   413,    -1,   414,    -1,   416,    -1,
     417,    -1,   418,    -1,   419,    -1,   420,    -1,   421,    -1,
     422,    -1,   423,    -1,   415,    -1,    51,    -1,    42,    -1,
      26,    33,    51,    -1,   116,    33,    51,    -1,   113,    33,
      51,    -1,    60,    -1,    94,    33,    51,    -1,   108,    33,
      51,    -1,    27,    33,    51,    -1,     3,    33,    51,    -1,
      84,    -1,    86,    -1,    88,    -1,    53,    33,    51,    -1,
      48,    33,    51,    -1,    49,    33,    51,    -1,    98,    33,
      51,    -1,    24,    33,   326,    -1,    63,    33,   326,    -1,
     112,    -1,   114,    33,    51,    -1,   106,    33,   326,    -1,
      25,    33,    80,    -1,    82,    33,   427,    -1,    82,    33,
      51,    -1,    41,    33,    51,    -1,   100,    33,    51,    -1,
     101,    33,    51,    -1,    58,    33,    51,    -1,    59,    33,
      51,    -1,    87,    -1,    46,    -1,    20,    33,   326,    -1,
      69,    33,    51,    -1,    64,    33,   326,    -1,    66,    33,
     326,    -1,    92,    33,   191,   287,   192,    -1,    65,    33,
     326,    -1,    73,    33,    80,    -1,    72,    33,    51,    -1,
      71,    -1,   103,    33,   326,    -1,    67,    33,    51,    -1,
      68,    33,    51,    -1,    61,    -1,    62,    -1,    85,    -1,
       5,    -1,   121,    -1,    43,    33,    51,    -1,   115,    -1,
      79,    -1,    40,    -1,   107,    -1,    54,    33,    51,    -1,
      55,    33,    51,    -1,    77,    33,    56,    -1,    77,    33,
      78,    -1,   102,    -1,    89,    -1,   133,    33,    80,    -1,
     134,    33,   424,    -1,    39,    33,   427,    -1,    21,    -1,
      83,    -1,    70,    -1,   123,    33,   326,    -1,    14,    33,
     264,    -1,     9,    33,   326,    -1,    11,    33,   264,    -1,
      12,    33,   326,    -1,    13,    33,    51,    -1,    10,    -1,
      15,    33,    51,    -1,    16,    33,    51,    -1,   159,    33,
      51,    -1,   160,    33,    51,    -1,   161,    33,    51,    -1,
     162,    33,    51,    -1,   163,    33,    51,    -1,   164,    33,
      51,    -1,   165,    33,    51,    -1,   166,    33,    51,    -1,
     167,    33,    51,    -1,   168,    33,    51,    -1,   169,    33,
      51,    -1,   170,    33,    51,    -1,   171,    33,    51,    -1,
     172,    33,    51,    -1,   173,    33,    51,    -1,   174,    33,
     326,    -1,   175,    33,   326,    -1,   176,    33,    51,    -1,
     177,    33,   427,    -1,   178,    33,   326,    -1,   179,    33,
     326,    -1,   183,    33,    51,    -1,   184,    33,    51,    -1,
     186,    33,   326,    -1,   187,    33,    51,    -1,   188,    33,
     326,    -1,   189,    33,   326,    -1,    80,   194,    80,    -1,
      51,    -1,    51,   194,    51,    -1,   195,   425,    -1,   426,
     425,    -1,   426,   196,    -1
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
     380,   388,   392,   394,   396,   398,   400,   402,   404,   405,
     411,   412,   421,   422,   431,   432,   443,   444,   453,   456,
     459,   461,   463,   468,   471,   475,   477,   479,   481,   485,
     489,   493,   497,   501,   504,   507,   512,   517,   522,   527,
     532,   537,   542,   547,   552,   557,   563,   565,   570,   575,
     580,   583,   585,   595,   602,   608,   616,   624,   627,   632,
     636,   642,   646,   648,   651,   653,   660,   664,   666,   672,
     676,   680,   685,   688,   691,   695,   697,   699,   702,   708,
     712,   714,   716,   719,   725,   729,   731,   733,   736,   742,
     746,   748,   750,   752,   755,   761,   765,   772,   776,   778,
     780,   782,   784,   786,   788,   790,   792,   794,   796,   798,
     800,   802,   804,   806,   808,   810,   812,   814,   816,   818,
     820,   822,   824,   827,   832,   836,   842,   844,   848,   851,
     854,   856,   859,   862,   864,   869,   872,   874,   879,   882,
     884,   889,   893,   899,   909,   911,   917,   921,   927,   935,
     945,   950,   953,   955,   961,   969,   974,   979,   982,   984,
     992,  1002,  1009,  1011,  1013,  1015,  1017,  1019,  1020,  1022,
    1024,  1026,  1029,  1032,  1035,  1041,  1045,  1052,  1056,  1058,
    1060,  1062,  1064,  1066,  1068,  1070,  1072,  1074,  1076,  1078,
    1080,  1082,  1084,  1086,  1088,  1090,  1092,  1094,  1096,  1098,
    1100,  1102,  1104,  1106,  1108,  1110,  1112,  1114,  1116,  1118,
    1120,  1122,  1124,  1126,  1128,  1130,  1132,  1134,  1136,  1142,
    1149,  1153,  1155,  1157,  1159,  1161,  1163,  1165,  1167,  1169,
    1171,  1173,  1175,  1177,  1179,  1181,  1187,  1194,  1202,  1208,
    1210,  1214,  1218,  1223,  1226,  1228,  1234,  1238,  1243,  1248,
    1255,  1259,  1265,  1269,  1272,  1278,  1282,  1289,  1294,  1297,
    1299,  1305,  1313,  1323,  1324,  1328,  1332,  1335,  1341,  1347,
    1354,  1358,  1366,  1375,  1381,  1387,  1394,  1398,  1406,  1415,
    1421,  1428,  1432,  1434,  1436,  1438,  1440,  1442,  1446,  1451,
    1458,  1460,  1463,  1465,  1467,  1469,  1471,  1473,  1474,  1475,
    1481,  1484,  1490,  1494,  1501,  1505,  1507,  1509,  1511,  1513,
    1515,  1517,  1519,  1521,  1523,  1525,  1527,  1529,  1531,  1533,
    1535,  1537,  1539,  1541,  1543,  1545,  1549,  1551,  1555,  1562,
    1564,  1566,  1568,  1570,  1574,  1576,  1580,  1587,  1590,  1596,
    1600,  1602,  1604,  1606,  1608,  1610,  1612,  1614,  1616,  1618,
    1620,  1622,  1624,  1626,  1628,  1630,  1632,  1634,  1636,  1638,
    1640,  1642,  1644,  1646,  1648,  1650,  1652,  1654,  1656,  1658,
    1660,  1664,  1668,  1672,  1674,  1678,  1682,  1686,  1690,  1692,
    1694,  1696,  1700,  1704,  1708,  1712,  1716,  1720,  1722,  1726,
    1730,  1734,  1738,  1742,  1746,  1750,  1754,  1758,  1762,  1764,
    1766,  1770,  1774,  1778,  1782,  1788,  1792,  1796,  1800,  1802,
    1806,  1810,  1814,  1816,  1818,  1820,  1822,  1824,  1828,  1830,
    1832,  1834,  1836,  1840,  1844,  1848,  1852,  1854,  1856,  1860,
    1864,  1868,  1870,  1872,  1874,  1878,  1882,  1886,  1890,  1894,
    1898,  1900,  1904,  1908,  1912,  1916,  1920,  1924,  1928,  1932,
    1936,  1940,  1944,  1948,  1952,  1956,  1960,  1964,  1968,  1972,
    1976,  1980,  1984,  1988,  1992,  1996,  2000,  2004,  2008,  2012,
    2016,  2020,  2022,  2026,  2029,  2032
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
     361,   365,   366,   370,   371,   372,   373,   374,   378,   378,
     379,   379,   381,   381,   383,   383,   385,   385,   390,   391,
     392,   393,   397,   399,   404,   405,   406,   408,   410,   412,
     414,   416,   418,   420,   422,   424,   426,   428,   430,   432,
     434,   436,   438,   440,   442,   446,   450,   452,   457,   461,
     465,   466,   470,   472,   474,   476,   478,   483,   485,   487,
     489,   491,   493,   499,   501,   506,   511,   513,   518,   520,
     522,   524,   526,   528,   530,   532,   534,   539,   543,   547,
     548,   551,   555,   557,   561,   562,   565,   569,   571,   575,
     576,   579,   580,   584,   586,   588,   590,   594,   595,   598,
     599,   600,   601,   602,   603,   604,   605,   606,   607,   608,
     609,   610,   611,   612,   613,   614,   615,   616,   617,   618,
     619,   620,   624,   626,   628,   630,   632,   634,   639,   641,
     643,   648,   650,   652,   657,   662,   664,   669,   673,   678,
     683,   693,   698,   704,   714,   719,   730,   736,   744,   754,
     768,   772,   774,   778,   785,   794,   803,   807,   809,   813,
     822,   833,   845,   847,   849,   851,   853,   858,   859,   860,
     861,   862,   864,   871,   873,   875,   877,   882,   883,   886,
     887,   888,   889,   890,   891,   892,   893,   894,   895,   896,
     897,   898,   899,   900,   901,   902,   903,   904,   905,   906,
     907,   908,   909,   910,   911,   912,   913,   914,   915,   916,
     917,   918,   919,   920,   921,   922,   923,   924,   928,   930,
     935,   936,   940,   941,   942,   943,   944,   945,   946,   947,
     948,   949,   950,   951,   952,   956,   958,   963,   964,   968,
     969,   973,   978,   983,   984,   987,   991,   994,   998,  1000,
    1002,  1004,  1008,  1011,  1012,  1013,  1014,  1017,  1021,  1022,
    1025,  1026,  1027,  1030,  1031,  1032,  1035,  1036,  1039,  1040,
    1041,  1042,  1043,  1044,  1046,  1047,  1048,  1049,  1050,  1051,
    1053,  1057,  1058,  1061,  1062,  1063,  1066,  1067,  1068,  1069,
    1072,  1073,  1076,  1077,  1078,  1079,  1080,  1083,  1083,  1083,
    1086,  1088,  1090,  1092,  1097,  1098,  1101,  1102,  1105,  1106,
    1107,  1108,  1109,  1110,  1111,  1114,  1115,  1116,  1117,  1118,
    1119,  1120,  1121,  1122,  1123,  1126,  1127,  1130,  1132,  1136,
    1137,  1138,  1139,  1142,  1143,  1146,  1148,  1152,  1154,  1158,
    1159,  1162,  1163,  1164,  1165,  1166,  1167,  1168,  1169,  1170,
    1171,  1172,  1173,  1174,  1175,  1176,  1177,  1178,  1179,  1180,
    1181,  1182,  1183,  1184,  1185,  1186,  1187,  1188,  1191,  1191,
    1193,  1194,  1195,  1196,  1197,  1198,  1199,  1200,  1201,  1202,
    1203,  1204,  1205,  1206,  1207,  1208,  1209,  1210,  1211,  1212,
    1213,  1214,  1215,  1217,  1218,  1219,  1220,  1221,  1222,  1223,
    1224,  1225,  1226,  1227,  1228,  1229,  1230,  1231,  1232,  1233,
    1234,  1235,  1236,  1237,  1238,  1239,  1240,  1241,  1242,  1243,
    1244,  1245,  1246,  1247,  1249,  1251,  1254,  1255,  1256,  1257,
    1258,  1259,  1260,  1261,  1262,  1264,  1265,  1266,  1267,  1268,
    1269,  1270,  1271,  1273,  1274,  1275,  1276,  1277,  1278,  1279,
    1280,  1281,  1282,  1283,  1284,  1285,  1286,  1287,  1288,  1289,
    1290,  1291,  1293,  1294,  1300,  1301,  1305,  1306,  1307,  1308,
    1311,  1319,  1320,  1324,  1325,  1334
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
       2,     2,     2,     2,     2,   193,     2,     2,     2,   197,
     191,   192,     2,     2,     2,     2,   198,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   194,   190,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   195,   199,   196,     2,     2,     2,     2,     2,     2,
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
     185,   186,   187,   188,   189
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1642;
  const int parser::yynnts_ = 228;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 160;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 200;

  const unsigned int parser::yyuser_token_number_max_ = 444;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1336 "DynareBison.yy"


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

