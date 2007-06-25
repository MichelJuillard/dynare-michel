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
#line 145 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 50:
#line 146 "DynareBison.yy"
    {driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 51:
#line 149 "DynareBison.yy"
    {driver.rplot();;}
    break;

  case 56:
#line 169 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 57:
#line 171 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 58:
#line 173 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 59:
#line 175 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 60:
#line 177 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 61:
#line 179 "DynareBison.yy"
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
#line 199 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 69:
#line 201 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 70:
#line 203 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 71:
#line 205 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 72:
#line 207 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 73:
#line 209 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 74:
#line 214 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 75:
#line 216 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 76:
#line 218 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 77:
#line 220 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 78:
#line 222 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 79:
#line 224 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 80:
#line 229 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 81:
#line 233 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 82:
#line 240 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 83:
#line 244 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 84:
#line 251 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 85:
#line 255 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 86:
#line 263 "DynareBison.yy"
    {driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 87:
#line 268 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 88:
#line 270 "DynareBison.yy"
    {(yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 89:
#line 272 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 90:
#line 274 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 91:
#line 276 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 92:
#line 278 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 93:
#line 280 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 94:
#line 282 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 95:
#line 284 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 96:
#line 286 "DynareBison.yy"
    {(yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 97:
#line 288 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 98:
#line 290 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 99:
#line 292 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 100:
#line 294 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 101:
#line 296 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 102:
#line 298 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 103:
#line 300 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 104:
#line 302 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 105:
#line 304 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 106:
#line 306 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 107:
#line 308 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 108:
#line 310 "DynareBison.yy"
    {(yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 109:
#line 315 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 110:
#line 317 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 111:
#line 321 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 112:
#line 323 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 113:
#line 327 "DynareBison.yy"
    {driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 114:
#line 332 "DynareBison.yy"
    {driver.end_endval();;}
    break;

  case 117:
#line 342 "DynareBison.yy"
    {driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 118:
#line 347 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 121:
#line 357 "DynareBison.yy"
    {driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 124:
#line 365 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 125:
#line 366 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 128:
#line 372 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 129:
#line 372 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 130:
#line 373 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 131:
#line 374 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 132:
#line 375 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 133:
#line 376 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 134:
#line 377 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 135:
#line 378 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 136:
#line 379 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 137:
#line 380 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 142:
#line 392 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 143:
#line 394 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val));;}
    break;

  case 144:
#line 398 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 146:
#line 401 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 147:
#line 403 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 148:
#line 405 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 149:
#line 407 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 150:
#line 409 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 151:
#line 411 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 152:
#line 413 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 153:
#line 415 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 154:
#line 417 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 155:
#line 419 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 156:
#line 421 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 157:
#line 423 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 158:
#line 425 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 159:
#line 427 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 160:
#line 429 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 161:
#line 431 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 162:
#line 433 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 163:
#line 435 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 164:
#line 437 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 165:
#line 441 "DynareBison.yy"
    {driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 166:
#line 445 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 167:
#line 447 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 168:
#line 451 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 169:
#line 455 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 172:
#line 465 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 173:
#line 467 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 174:
#line 469 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 175:
#line 471 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 176:
#line 473 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 177:
#line 478 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 178:
#line 480 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 179:
#line 482 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 180:
#line 484 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 181:
#line 486 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 182:
#line 488 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 183:
#line 494 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 184:
#line 496 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 185:
#line 498 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 186:
#line 500 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 187:
#line 502 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 188:
#line 504 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 189:
#line 506 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 190:
#line 508 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 191:
#line 513 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 192:
#line 518 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 193:
#line 520 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 194:
#line 525 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 195:
#line 527 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 196:
#line 529 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 197:
#line 531 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 198:
#line 533 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 199:
#line 535 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 200:
#line 537 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 201:
#line 539 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 202:
#line 541 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 203:
#line 546 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 204:
#line 550 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 208:
#line 562 "DynareBison.yy"
    {driver.check();;}
    break;

  case 209:
#line 564 "DynareBison.yy"
    {driver.check();;}
    break;

  case 213:
#line 576 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 214:
#line 578 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 219:
#line 591 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 220:
#line 593 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 221:
#line 595 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 222:
#line 597 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 246:
#line 629 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 247:
#line 631 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 248:
#line 633 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 249:
#line 635 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 250:
#line 637 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 251:
#line 639 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 252:
#line 644 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 253:
#line 646 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 254:
#line 648 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 255:
#line 653 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 256:
#line 655 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 257:
#line 657 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 258:
#line 662 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 259:
#line 667 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 260:
#line 669 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 262:
#line 678 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 263:
#line 683 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 264:
#line 688 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 265:
#line 698 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 266:
#line 703 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 267:
#line 709 "DynareBison.yy"
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

  case 268:
#line 719 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 269:
#line 724 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 270:
#line 735 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 271:
#line 741 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 272:
#line 749 "DynareBison.yy"
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

  case 273:
#line 759 "DynareBison.yy"
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

  case 274:
#line 773 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 275:
#line 777 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 276:
#line 779 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 277:
#line 783 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 278:
#line 790 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 279:
#line 799 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 280:
#line 808 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 281:
#line 812 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 282:
#line 814 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 283:
#line 818 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 284:
#line 827 "DynareBison.yy"
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

  case 285:
#line 838 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 286:
#line 850 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 287:
#line 852 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 288:
#line 854 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 289:
#line 856 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 290:
#line 858 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 291:
#line 862 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 295:
#line 867 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 296:
#line 869 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 297:
#line 876 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 298:
#line 878 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 299:
#line 880 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 300:
#line 882 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 342:
#line 933 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 343:
#line 935 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 359:
#line 961 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 360:
#line 963 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 361:
#line 967 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 362:
#line 968 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 365:
#line 978 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 366:
#line 983 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 369:
#line 992 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val));;}
    break;

  case 370:
#line 995 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 371:
#line 999 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 372:
#line 1003 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 373:
#line 1005 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 374:
#line 1007 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 375:
#line 1009 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 376:
#line 1012 "DynareBison.yy"
    {driver.set_osr_params();;}
    break;

  case 377:
#line 1015 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 378:
#line 1016 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 379:
#line 1017 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 380:
#line 1018 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 381:
#line 1021 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 382:
#line 1022 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 383:
#line 1023 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 384:
#line 1024 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 389:
#line 1035 "DynareBison.yy"
    {driver.set_olr_inst();;}
    break;

  case 390:
#line 1039 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 393:
#line 1046 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 394:
#line 1047 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 395:
#line 1048 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val));;}
    break;

  case 396:
#line 1051 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 397:
#line 1052 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 398:
#line 1053 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 399:
#line 1056 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 400:
#line 1057 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 401:
#line 1060 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 402:
#line 1061 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 403:
#line 1062 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 404:
#line 1063 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 405:
#line 1064 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 406:
#line 1065 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 407:
#line 1067 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 408:
#line 1068 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 409:
#line 1069 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 410:
#line 1070 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 411:
#line 1071 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 412:
#line 1072 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 413:
#line 1075 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 419:
#line 1087 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 420:
#line 1088 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 421:
#line 1089 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 422:
#line 1090 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 423:
#line 1093 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 424:
#line 1094 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 426:
#line 1098 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 427:
#line 1099 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 428:
#line 1100 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 429:
#line 1101 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 430:
#line 1104 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 431:
#line 1104 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 433:
#line 1108 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 434:
#line 1110 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 435:
#line 1112 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 436:
#line 1114 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 458:
#line 1150 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 459:
#line 1152 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 466:
#line 1166 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 467:
#line 1168 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 468:
#line 1171 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 469:
#line 1172 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 470:
#line 1173 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 471:
#line 1174 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 472:
#line 1175 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 473:
#line 1176 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 474:
#line 1177 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 475:
#line 1178 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 476:
#line 1179 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 477:
#line 1180 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 478:
#line 1181 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 479:
#line 1182 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 480:
#line 1183 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 481:
#line 1184 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 482:
#line 1185 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 483:
#line 1186 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 484:
#line 1187 "DynareBison.yy"
    {driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 485:
#line 1188 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 486:
#line 1189 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 487:
#line 1190 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 488:
#line 1191 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 489:
#line 1193 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 490:
#line 1194 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 491:
#line 1195 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 492:
#line 1197 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 493:
#line 1198 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 494:
#line 1199 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1200 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 496:
#line 1201 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 497:
#line 1202 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 498:
#line 1203 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 499:
#line 1204 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 500:
#line 1205 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 501:
#line 1206 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 502:
#line 1207 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 504:
#line 1209 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 505:
#line 1210 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 506:
#line 1211 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 507:
#line 1212 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 508:
#line 1213 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 509:
#line 1214 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 510:
#line 1215 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 511:
#line 1216 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 512:
#line 1217 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 513:
#line 1218 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 514:
#line 1219 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 515:
#line 1220 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 516:
#line 1221 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 517:
#line 1222 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 518:
#line 1223 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 519:
#line 1224 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 520:
#line 1225 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 521:
#line 1226 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 522:
#line 1227 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 523:
#line 1228 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 524:
#line 1229 "DynareBison.yy"
    {driver.option_num("olr_beta", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 525:
#line 1232 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 526:
#line 1234 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 527:
#line 1236 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 528:
#line 1237 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 529:
#line 1238 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 530:
#line 1239 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 531:
#line 1240 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 532:
#line 1241 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 533:
#line 1242 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 534:
#line 1243 "DynareBison.yy"
    {driver.option_num("load_mh_file", "-1");;}
    break;

  case 535:
#line 1244 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 536:
#line 1246 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 537:
#line 1247 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 538:
#line 1248 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1249 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 540:
#line 1250 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1251 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 542:
#line 1252 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 543:
#line 1253 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 544:
#line 1256 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 546:
#line 1265 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 547:
#line 1268 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 548:
#line 1270 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 549:
#line 1278 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2242 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -959;
  const short int
  parser::yypact_[] =
  {
       778,    34,    40,   -59,   -73,   292,   238,    39,   -12,     2,
     -46,   189,   -44,    -5,    55,    90,   465,   285,   535,   -19,
     122,   259,   186,   217,   268,   215,   220,   268,   305,    95,
    -959,   363,   366,   268,   258,   406,   582,   593,   237,   248,
     268,   401,   413,   482,   268,    76,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,   483,   762,   487,   820,  -959,   543,   275,  -959,
     546,   635,   521,    45,   422,   499,   461,   608,   623,   683,
    -959,   665,    -1,    52,   214,   245,   648,   623,   707,   719,
     604,  -959,   424,    50,    81,   946,   686,  -959,  1096,   242,
     267,   689,  -959,  1096,   289,   320,   654,   376,   730,   629,
     989,   556,   556,   392,    81,   634,  -959,    93,  -959,   546,
    -959,  1149,   394,  -959,   890,   409,   410,   678,   414,   679,
     440,   685,   455,   464,  -959,  -959,  -959,   776,  -959,   783,
     784,   789,   791,   793,   804,   805,   806,   811,   817,   822,
     823,  -959,   690,   695,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,   826,
     827,   829,  -959,   729,   705,  -959,  -959,  -959,   709,   792,
     124,   272,  -959,   841,   -22,  -959,  -959,   720,  -959,   722,
    -959,  -959,   807,   426,  -959,   808,   616,   861,    57,  -959,
     815,  -959,  -959,   863,  -959,  -959,   870,   871,   872,   874,
     875,  -959,  -959,   880,   885,   886,   887,   888,   896,  -959,
    -959,   898,   907,  -959,  -959,  -959,  -959,   908,   911,  -959,
    -959,   -18,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,   912,   866,  -959,   867,  -959,   868,    99,  -959,
     749,   881,   830,   891,   126,  -959,   892,   833,   894,   216,
    -959,   818,    61,  -959,   448,   919,   821,   824,  -959,   873,
    -959,   112,   839,   840,   949,  -959,  -959,   151,  -959,  -959,
    -959,  -959,   905,   909,   100,  -959,  -959,  -959,   834,   946,
     946,   842,   848,   849,   855,   869,   882,   899,   901,   902,
     906,   946,   794,   920,   492,  -959,   981,   994,   995,   996,
    1000,  1002,  -959,  -959,  -959,  1008,  1009,  1024,  1035,  -959,
    1037,  -959,  1047,  1049,  -959,  -959,   249,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,   900,   503,   278,  -959,  -959,  -959,   957,  1007,  -959,
     929,  -959,  -959,  -959,   939,   989,   989,   951,   953,   954,
     955,   956,   958,   960,   965,   968,   969,   989,   708,  -959,
     279,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,   300,  -959,   103,    35,   419,
    -959,  -959,  -959,   423,  -959,  -959,   431,  -959,  -959,  1085,
    -959,   469,  -959,  -959,  -959,  -959,  -959,   998,  1051,  -959,
    -959,  1017,  1066,  -959,  -959,  1021,  1071,  -959,  -959,  1111,
      75,  1112,  1104,    75,  1106,  1078,  1108,    16,  1113,  1119,
    1092,  1094,   762,  1126,  1127,  1138,  1134,   820,  1135,  1027,
    1030,  1114,   368,  1155,  -959,  -959,  1150,   546,  1033,  -959,
    -959,  1039,    84,  1120,  1043,   197,  1128,   946,  -959,  -959,
    -959,  1042,  1170,  1175,  1186,  1187,  1189,  1171,   572,  1177,
    1190,  1191,  1192,  1193,  1167,  1087,  1209,   665,   203,  1173,
    1217,  1117,  -959,  -959,  -959,   287,  1118,   236,  1123,  -959,
    -959,  1124,   236,  1125,  -959,  -959,    94,  -959,  -959,  -959,
    1176,  1133,  -959,  1198,    78,  -959,   202,  -959,   565,  -959,
    1159,  1185,   328,    50,   204,  1129,    32,  -959,  -959,   946,
    1141,   574,   946,   946,   946,   946,   946,   946,   946,   946,
     946,   946,   617,   946,   946,   946,   946,   946,  -959,   946,
    -959,  -959,  1221,  1248,  1265,  1291,  1317,  1320,   236,  1405,
    1406,   658,  1407,  1408,  1410,  1096,   253,  1382,   853,  -959,
     927,   257,  -959,  1337,  -959,    94,  1321,   650,   989,   989,
     989,   989,   989,   989,   989,   989,   989,   989,   680,   989,
     989,   989,   989,   989,  1305,   556,   286,   291,  -959,  -959,
    -959,   946,   632,    24,    93,  1307,   546,  1308,  1149,   332,
    1427,   890,   341,  -959,  1344,  -959,  1345,  -959,  1347,  -959,
    -959,  1432,  1433,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  1425,    12,  -959,  -959,  -959,  -959,  1313,  -959,  -959,
    1318,  -959,  -959,  -959,  -959,  1319,  -959,  1429,  1322,  1323,
    1324,   946,  -959,  -959,  -959,  -959,  -959,   470,  1325,  -959,
    -959,   479,  1326,  1234,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  1314,
    -959,  -959,  -959,   500,  -959,  1402,  1404,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,   659,  1329,  1353,  1354,  1412,
    1356,   236,  1414,  1335,   236,  -959,  1445,  1446,  1336,  -959,
     623,  1466,  -959,  -959,  -959,   989,  -959,  -959,  -959,  1467,
     474,  -959,  -959,  -959,  1341,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,   418,   190,  -959,  1422,   946,
    1423,   181,   915,   501,   787,   819,   825,   963,   970,  1022,
    1028,  1077,  1084,  1090,  -959,   574,   574,  1141,  1141,  1361,
    1130,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,   537,   946,  -959,
    1426,  1240,  -959,   538,  -959,  1343,  1136,  1142,  1148,  1156,
    1162,  1168,  1174,  1182,  1188,  1194,  -959,   650,   650,  1321,
    1321,  1361,  -959,  -959,  -959,   539,  -959,   547,  1200,    35,
    1348,  -959,  -959,    38,   946,  -959,  -959,  -959,  -959,  -959,
    -959,   552,  -959,  -959,  -959,   588,  -959,  -959,  -959,  -959,
    -959,  1346,  -959,  -959,  -959,  1431,  -959,  -959,  1350,  1474,
    -959,  -959,  1246,  -959,   380,  -959,   381,  -959,  1434,  -959,
     525,  -959,  -959,  -959,  -959,  -959,  -959,   236,   287,  1372,
     236,  1376,  1378,  -959,  1357,  -959,  -959,  1483,   507,   989,
    1252,  1476,   565,  -959,   873,   873,   873,   204,  -959,   236,
    -959,  1486,  1258,  1487,  1470,   946,   946,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  1362,  -959,
    1264,   946,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,    24,  -959,
    -959,  -959,   946,  1208,  -959,  -959,  1472,  -959,  1322,   946,
    -959,  -959,   602,  -959,   603,  1358,  1314,  -959,  -959,  1387,
    1389,  1390,   236,  1368,   236,   236,  -959,   946,  -959,  1270,
    -959,  -959,  -959,  1369,   160,   244,   284,   192,  1370,   946,
    -959,   946,  1366,    10,  1276,   915,  -959,  -959,  1282,  1214,
    -959,  -959,  1498,  1288,  -959,  -959,  1396,  -959,   236,   236,
     236,  1397,  -959,  1375,  1377,  1294,  -959,   873,  -959,  -959,
    -959,   236,  -959,  1300,  1306,  1485,  1374,  1488,  1411,  -959,
    -959,  -959,   946,  -959,    27,  1403,  -959,  1409,   236,  -959,
    -959,  -959,   506,  1380,  -959,  -959,  -959,  1492,  1381,    82,
    1312,  1464,  -959,   236,    91,  1388,  -959,  -959,  -959,  1496,
    -959,   671,   699,   946,   533,  -959,  -959,  -959,  1383,  1415,
    1416,  -959,  -959,  1220,  -959,  -959,   946,  -959,  -959,  -959,
     236,   236,  -959,  1226,  1417,  -959,  -959,   236,  -959
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     430,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     2,     4,    29,    30,
      46,    47,    48,    45,     5,     6,     7,    12,     9,    10,
      11,     8,    13,    14,    15,    16,    17,    18,    19,    23,
      25,    24,    20,    21,    22,    26,    27,    28,    31,    32,
      33,    38,    39,    34,    35,    36,    37,    40,    41,    42,
      43,    44,     0,     0,     0,     0,   399,     0,     0,   208,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   250,
     297,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   128,     0,     0,     0,     0,     0,   381,     0,     0,
       0,     0,   377,     0,     0,     0,    76,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   213,     0,   203,     0,
     219,     0,     0,   433,     0,     0,     0,    58,     0,    64,
       0,    70,     0,     0,     1,     3,   458,     0,   541,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   448,   457,     0,   449,   454,   452,   455,   453,   450,
     451,   441,   442,   443,   444,   445,   446,   447,   466,     0,
       0,     0,   460,   465,     0,   462,   461,   463,     0,     0,
     396,     0,   392,     0,     0,   211,   212,     0,    82,     0,
      49,   409,     0,     0,   403,     0,     0,     0,     0,   116,
       0,   515,   532,     0,   520,   498,     0,     0,     0,     0,
       0,   512,   513,     0,     0,     0,     0,     0,     0,   534,
     508,     0,     0,   519,   533,   514,   497,     0,     0,   518,
     516,     0,   302,   338,   327,   303,   304,   305,   306,   307,
     308,   309,   310,   311,   312,   313,   314,   315,   316,   317,
     318,   319,   320,   321,   322,   323,   324,   325,   326,   328,
     329,   330,   331,   332,   333,   334,   335,   336,   337,   339,
     340,   341,   246,     0,   299,     0,   263,     0,     0,   260,
       0,     0,     0,     0,     0,   282,     0,     0,     0,     0,
     276,     0,     0,   120,     0,     0,     0,     0,    84,     0,
     471,     0,     0,     0,     0,   528,   527,     0,   415,   416,
     417,   418,     0,     0,     0,   171,    89,    90,    88,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   368,     0,     0,     0,     0,
       0,     0,   476,   477,   478,     0,     0,     0,     0,   521,
       0,   485,     0,     0,   386,   387,     0,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   237,   239,
     240,   241,   242,   243,   244,   245,   236,   238,   385,   383,
     389,     0,     0,     0,   379,   376,    79,    74,     0,    55,
       0,    80,   146,   147,   166,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   431,   145,
       0,   345,   350,   346,   347,   348,   349,   351,   352,   353,
     354,   355,   356,   357,   358,     0,    51,     0,     0,     0,
     216,   217,   218,     0,   206,   207,     0,   224,   221,     0,
     439,     0,   438,   440,   435,   370,    61,    56,     0,    52,
      67,    62,     0,    53,    73,    68,     0,    54,   365,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   390,   391,     0,     0,     0,    83,
      50,     0,     0,     0,     0,     0,     0,     0,   114,   115,
     251,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     248,     0,   262,   258,   259,   291,     0,   291,     0,   280,
     281,     0,   291,     0,   274,   275,     0,   118,   119,   111,
       0,     0,    85,     0,     0,   140,     0,   141,     0,   136,
       0,     0,     0,     0,     0,     0,     0,   169,   170,     0,
      96,    97,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    86,     0,
     366,   367,     0,     0,     0,     0,     0,     0,   291,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   371,
       0,     0,    77,    75,    81,     0,   153,   154,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   168,   201,
     202,     0,     0,   193,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    59,    57,    65,    63,    71,    69,   537,
     257,     0,     0,   538,   539,   540,   536,   542,   489,   492,
     491,     0,     0,   490,   493,   494,   529,     0,   530,   456,
       0,   543,   499,   517,   464,     0,   400,     0,   396,     0,
       0,     0,   469,   210,   209,   412,   407,     0,     0,   406,
     401,     0,     0,     0,   531,   479,   522,   523,   495,   496,
     501,   504,   505,   502,   510,   511,   500,   507,   506,     0,
     509,   301,   298,     0,   247,     0,     0,   286,   293,   287,
     292,   289,   294,   288,   290,     0,     0,     0,   268,     0,
       0,   291,     0,     0,   291,   254,     0,     0,     0,   113,
       0,     0,   129,   138,   139,     0,   143,   125,   124,     0,
       0,   123,   126,   127,     0,   132,   130,   525,   526,   414,
     425,   427,   428,   429,   426,     0,   419,   423,     0,     0,
       0,     0,   109,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    87,    92,    91,    93,    94,    95,
       0,   475,   483,   468,   474,   480,   481,   524,   472,   482,
     488,   487,   473,   470,   486,   388,   382,     0,     0,   374,
       0,     0,   378,     0,    78,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   144,   149,   148,   150,
     151,   152,   432,   344,   342,     0,   359,     0,     0,     0,
       0,   198,   199,     0,     0,   215,   214,   205,   204,   223,
     220,     0,   535,   437,   434,     0,    60,    66,    72,   256,
     255,   545,   547,   549,   548,     0,   459,   467,     0,     0,
     398,   397,     0,   408,     0,   402,     0,   117,     0,   363,
       0,   300,   249,   264,   296,   295,   261,   291,   291,     0,
     291,     0,     0,   279,     0,   253,   252,     0,     0,     0,
       0,     0,     0,   134,     0,     0,     0,     0,   413,   291,
     424,     0,     0,     0,     0,     0,     0,   108,    98,    99,
     100,   101,   102,   103,   104,   105,   106,   107,     0,   384,
       0,     0,   372,   380,   167,   155,   156,   157,   158,   159,
     160,   161,   162,   163,   164,   343,   360,   200,   192,   191,
     195,   196,     0,     0,   222,   436,     0,   544,   396,     0,
     393,   410,     0,   404,     0,     0,     0,   503,   265,     0,
       0,     0,   291,     0,   291,   291,   277,     0,   112,     0,
     142,   484,   122,     0,     0,     0,     0,   420,     0,     0,
     174,     0,   182,     0,     0,   110,   369,   375,     0,     0,
     197,   546,     0,     0,   411,   405,     0,   364,   291,   291,
     291,     0,   285,     0,     0,     0,   165,     0,   137,   133,
     131,   291,   421,     0,     0,     0,   177,     0,     0,   173,
     373,   194,     0,   394,   291,   270,   266,   269,   291,   283,
     278,   121,     0,     0,   176,   175,   181,     0,   179,     0,
       0,     0,   362,   291,     0,     0,   135,   422,   178,     0,
     188,     0,     0,     0,     0,   187,   186,   395,     0,   271,
       0,   284,   180,     0,   185,   172,     0,   184,   183,   361,
     291,   291,   190,     0,   272,   267,   189,   291,   273
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -959,  -959,  1504,  -959,  -959,  -959,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -341,  -959,
    -959,  -959,  -959,  -111,  -226,  -959,  -959,  1229,  -959,   624,
    -959,  -959,  -959,  -959,  -959,  -959,  -863,  -550,  -135,  -546,
    -959,  -959,  -959,  1413,  -297,  -959,  -959,  -959,  -959,   691,
    -959,  -959,   889,  -959,  -959,  1041,  -959,  -959,   895,  -959,
    -959,  -129,   -23,  -606,  -482,  -959,  -959,  1253,  -959,  -959,
    -958,  -959,  -959,  1243,  -959,  -959,  1249,  -881,  -519,  -959,
    -959,  1018,  -959,  1424,   913,  -959,   568,  -959,  -959,  -959,
    -959,  1203,  -959,  -959,  -959,  -959,  -959,  -959,   944,  1437,
    -959,  -959,  -959,  1360,  -674,  -959,  -959,  -959,  -959,  -959,
     990,  -959,   637,  -744,  -959,  -959,  -959,  -959,  -959,   904,
    -959,   -63,  1070,  -959,  -959,  1069,  -959,  -959,   -93,  -959,
    1455,  -959,  -959,  -959,  -959,  -959,  -959,  -959,   -97,  -959,
    -959,  -134,  -545,  -959,  -959,  -959,  -959,   -99,   -80,   -76,
     -72,   -71,  -959,  -959,   -92,   -69,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,   -70,  -959,  -959,  -959,  -959,  -959,
     -60,   -56,   -65,   -52,   -51,   -49,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,   -88,   -84,   -47,  -959,  -959,  -959,  -959,
    -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,  -959,   893,
    -959,  1048
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    45,    46,    47,    48,    49,    50,    51,    52,    53,
     158,   160,   162,   137,    54,    55,    56,    57,   362,   803,
      58,   326,    59,   228,   229,    60,   322,   323,   780,   781,
      61,   329,   936,   935,  1013,   784,   574,   575,   576,   577,
     439,    62,    63,   344,   345,  1023,  1094,    64,   662,   663,
      65,   463,   464,    66,   214,   215,    67,   459,   460,    68,
     466,   384,   112,   768,   683,    69,   308,   309,   310,   756,
     998,    70,   319,   320,    71,   314,   315,   757,   999,    72,
     261,   262,    73,   440,   441,    74,   909,   910,    75,    76,
     364,   365,    77,    78,   412,    79,    80,    81,   385,   386,
      82,    83,   211,   212,   513,    84,    85,    86,    87,   337,
     338,   795,   796,   797,    88,   140,   654,    89,   471,   472,
     181,   182,   183,    90,   203,   204,    91,   387,   388,   389,
     390,   391,   392,   393,   394,   395,   396,   397,   398,   399,
     400,   401,   402,   783,   403,   404,   405,   184,   185,   186,
     187,   188,   270,   271,   406,   444,   274,   275,   276,   277,
     278,   279,   280,   281,   445,   283,   284,   285,   286,   287,
     446,   447,   448,   449,   450,   451,   407,   294,   295,   408,
     339,   340,   341,   189,   190,   454,   299,   300,   301,   473,
     191,   192,   193,   194,   195,   196,   197,   207,   698,   892,
     692,   693
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       129,   130,   529,   134,   135,   438,   324,   216,   590,   591,
     143,   686,   265,   461,   264,   152,   155,   156,   263,   272,
     602,   163,   467,   296,   773,   470,   205,   297,   774,   845,
     206,   266,   202,   782,   899,   267,   758,  1000,   760,   268,
     269,   282,   273,   763,   442,   442,   290,   588,   462,   443,
     443,   288,   940,   452,   452,   289,   465,   453,   453,   291,
     292,  1056,   293,   891,   298,   799,   871,   690,   104,   748,
     628,  1014,  1015,  1016,   305,   872,   164,   659,   750,   302,
     980,  1066,   106,     1,     2,    92,   660,    98,   528,   981,
     103,    94,   567,     3,     4,     5,   219,   747,   529,   827,
       6,    96,    97,   342,     7,     8,     9,   752,    10,   772,
      11,    12,    13,    14,   108,   517,   113,   680,   174,   547,
     422,   305,   342,    15,   680,   342,    16,   334,   138,   423,
     553,   587,   306,   765,   658,   749,   303,   227,   335,    17,
     518,   321,   123,   751,   548,   765,   139,  1057,   311,   105,
      18,    19,    20,  1115,   336,   114,    21,   559,   424,   304,
     588,   873,  1090,   107,   109,   755,    22,    23,    24,   800,
    1058,    25,   307,    26,    27,    28,    29,    30,   893,   306,
     753,   691,    31,    32,  1072,   874,   723,    33,    34,    35,
      36,  1048,   801,   377,  1081,    93,   661,    37,    38,   982,
      39,    95,   422,  1100,    40,   220,   312,    41,    42,    43,
      44,   423,   343,   681,   682,   115,   425,   426,   754,   307,
    1091,  1092,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   343,   766,   767,   343,   775,   311,   436,   316,   437,
     424,   573,   921,  1093,   716,   924,   313,   564,   802,   578,
     116,   804,   805,   806,   807,   808,   809,   810,   811,   812,
     813,   511,   815,   816,   817,   818,   819,   316,   820,   109,
     790,   101,   790,   940,   579,  1049,   209,   109,   748,   209,
     102,   944,   124,   109,   790,   512,   422,   750,   583,   841,
     636,   637,   125,   747,   312,   423,   317,   109,   425,   426,
     109,   945,   648,   514,   427,   428,   429,   430,   431,   432,
     433,   434,   435,   584,  1032,  1050,   752,   109,   119,   436,
     868,   437,   302,   573,   424,   317,   422,   120,   109,   748,
     791,   749,   791,   109,   313,   423,   318,   109,   750,   751,
     649,   650,   651,   652,   791,   653,   126,   302,   109,   110,
     111,   939,   210,  1051,   792,   210,   792,   720,   793,   794,
     793,   794,   776,   742,   424,   318,   109,   752,   792,   302,
     902,   109,   793,   794,   755,   131,   753,   127,   128,   303,
     132,   133,   425,   426,   787,   136,   625,   782,   427,   428,
     429,   430,   431,   432,   433,   434,   435,   150,   151,  1001,
     302,  1003,   409,   436,   303,   437,   788,   573,   153,   154,
     709,   626,   109,   836,   754,   625,   655,   842,   144,   710,
    1018,   109,   425,   426,   216,   755,   303,   410,   427,   428,
     429,   430,   431,   432,   433,   434,   435,   655,   205,   145,
     631,   656,   206,   436,   202,   437,   864,   573,   265,   414,
     264,   866,    99,   100,   263,   272,   417,   303,   942,   296,
     109,   109,   657,   297,   773,   773,   773,   266,   774,   774,
     774,   267,   302,  1095,   302,   268,   269,   282,   273,   569,
     415,   157,   290,  1041,   330,  1043,  1044,   288,  1107,   302,
     302,   289,   880,   159,   477,   291,   292,   960,   293,   717,
     298,   884,   721,   846,   847,   848,   849,   850,   851,   852,
     853,   854,   855,   418,   857,   858,   859,   860,   861,  1065,
     481,  1067,   773,   610,   141,   743,   774,   142,   227,   303,
     461,   303,  1073,   983,   629,   485,   419,  1086,  1008,   879,
     991,   993,   470,   331,   302,  1082,   303,   303,   422,  1085,
     302,   478,   456,   332,   468,   937,   664,   423,   442,   302,
     666,   231,   161,   443,  1099,   462,   208,   452,   668,   474,
     475,   453,   363,   465,   479,   680,   200,   482,   938,   223,
     302,   665,   221,   630,   765,   667,   424,   227,   522,   367,
     222,  1114,   486,   669,   523,   233,   234,  1096,  1118,   201,
     483,   303,   235,   837,  1024,  1025,   671,   303,   843,   236,
     777,   932,  1108,  1104,   731,   487,   303,   302,   302,   302,
    1028,   224,   778,   732,   488,   117,   118,   302,   779,   225,
     903,   672,   302,   865,   867,   253,   933,   303,   946,   905,
     930,  1029,   256,   166,   425,   426,   881,   198,  1033,   885,
     427,   428,   429,   430,   431,   432,   433,   434,   435,   928,
     911,   258,   996,   947,   213,   436,  1045,   437,   302,   573,
     231,  1091,  1092,   259,   303,   303,   303,   217,  1053,   260,
    1054,   218,   302,   302,   303,   200,   232,   997,   226,   303,
     174,   179,   180,  1105,  1106,   121,   122,   959,   963,   975,
     830,   914,   529,   227,   233,   234,   175,   976,   201,   831,
     915,   235,   984,   889,   605,   606,   230,   607,   236,   237,
     238,  1080,   925,   239,   240,   303,   241,   242,   321,   243,
     244,   245,   246,   247,   248,   249,   250,   251,   252,   303,
     303,   890,   146,   147,   253,   325,   176,   254,   985,   255,
     926,   256,  1103,   148,   149,   603,   604,   605,   606,   257,
     607,   327,  1034,  1035,   328,  1113,   363,   177,   178,   411,
     258,   167,   168,   169,   170,   171,   172,   173,   525,   814,
     416,   420,   259,   213,   526,     1,     2,   174,   260,   421,
     651,   652,   869,   653,  1009,     3,     4,     5,   870,   458,
     179,   180,     6,   175,   476,   480,     7,     8,     9,   489,
      10,   484,    11,    12,    13,    14,   490,   491,   649,   650,
     651,   652,   492,   653,   493,    15,   494,   502,    16,   167,
     168,   169,   170,   171,   172,   173,   199,   495,   496,   497,
     200,    17,   856,   176,   498,   174,   649,   650,   651,   652,
     499,   653,    18,    19,    20,   500,   501,   503,    21,   504,
     505,   175,   506,   201,   177,   178,   507,   508,    22,    23,
      24,   509,   510,    25,   516,    26,    27,    28,    29,    30,
     519,   992,   520,   994,    31,    32,   555,   521,   524,    33,
      34,    35,    36,   366,   527,   530,   531,   179,   180,    37,
      38,   176,    39,   532,   533,   534,    40,   535,   536,    41,
      42,    43,    44,   537,   367,   422,   368,   369,   538,   539,
     540,   541,   177,   178,   423,   603,   604,   605,   606,   542,
     607,   543,   603,   604,   605,   606,   235,   607,   370,   371,
     544,   545,   346,   236,   546,   549,   550,   551,   552,   948,
     330,   347,   570,   424,   608,   179,   180,   603,   604,   605,
     606,   556,   607,   603,   604,   605,   606,   557,   607,   346,
     562,   558,   561,   372,   563,   373,   256,   374,   347,   566,
     348,   949,   582,   571,   572,   585,   376,   950,   346,   586,
     377,   603,   604,   605,   606,   589,   607,   347,   378,   379,
     380,   580,   581,   592,   381,   382,   383,   348,   213,   593,
     594,   425,   426,   839,   612,   469,   595,   427,   428,   429,
     430,   431,   432,   433,   434,   435,   348,   613,   614,   615,
     596,   422,   436,   616,   437,   617,   573,   627,   349,   350,
     423,   618,   619,   597,   351,   352,   353,   354,   355,   356,
     357,   358,   359,   603,   604,   605,   606,   620,   607,   360,
     598,   361,   599,   600,   840,   349,   350,   601,   621,   424,
     622,   351,   352,   353,   354,   355,   356,   357,   358,   359,
     623,   609,   624,   632,   349,   350,   360,   633,   361,   634,
     351,   352,   353,   354,   355,   356,   357,   358,   359,   366,
     635,   603,   604,   605,   606,   360,   607,   361,   603,   604,
     605,   606,   638,   607,   639,   640,   641,   642,   670,   643,
     367,   644,   368,   369,   673,   951,   645,   425,   426,   646,
     647,   674,   952,   427,   428,   429,   430,   431,   432,   433,
     434,   435,   235,   675,   370,   371,   676,   677,   436,   236,
     437,   678,   366,   679,   684,   685,   330,   687,   688,   689,
     603,   604,   605,   606,   694,   607,   603,   604,   605,   606,
     695,   607,   696,   367,   697,   368,   369,   700,   701,   372,
     702,   373,   256,   374,   953,   703,   705,   706,   711,   375,
     954,   707,   376,   714,   708,   235,   377,   370,   371,   715,
     718,   712,   236,   719,   378,   379,   380,   691,   722,   330,
     381,   382,   383,   730,   213,   603,   604,   605,   606,   733,
     607,   725,   603,   604,   605,   606,   726,   607,   603,   604,
     605,   606,   372,   607,   373,   256,   374,   727,   728,   955,
     729,   734,   735,   736,   737,   376,   956,   738,   739,   377,
     745,   740,   957,   744,   746,   759,   769,   378,   379,   380,
     761,   762,   764,   381,   382,   383,   798,   213,   603,   604,
     605,   606,   821,   607,   649,   650,   651,   652,   771,   653,
     649,   650,   651,   652,   607,   653,   649,   650,   651,   652,
     822,   653,   958,   770,   649,   650,   651,   652,   965,   653,
     649,   650,   651,   652,   966,   653,   649,   650,   651,   652,
     967,   653,   649,   650,   651,   652,   823,   653,   968,   785,
     649,   650,   651,   652,   969,   653,   649,   650,   651,   652,
     970,   653,   649,   650,   651,   652,   971,   653,   603,   604,
     605,   606,   824,   607,   972,   786,   603,   604,   605,   606,
     973,   607,   603,   604,   605,   606,   974,   607,   603,   604,
     605,   606,   977,   607,   603,   604,   605,   606,   825,   607,
    1030,   826,   603,   604,   605,   606,  1061,   607,   603,   604,
     605,   606,  1112,   607,   603,   604,   605,   606,  1116,   607,
     649,   650,   651,   652,   907,   653,   603,   604,   605,   606,
     962,   607,   603,   604,   605,   606,   990,   607,   649,   650,
     651,   652,  1010,   653,   603,   604,   605,   606,  1020,   607,
     603,   604,   605,   606,  1027,   607,   603,   604,   605,   606,
    1046,   607,   603,   604,   605,   606,  1059,   607,   603,   604,
     605,   606,  1060,   607,   603,   604,   605,   606,  1063,   607,
     603,   604,   605,   606,  1071,   607,   828,   829,   832,   833,
    1074,   834,   838,   844,   653,   862,  1075,   876,   878,   882,
     886,   887,  1097,   888,   889,   890,   891,   895,   896,   897,
     898,   908,   912,   512,   913,   900,   901,   904,   906,   916,
     917,   918,   919,   920,   922,   923,   925,   926,   927,   929,
     931,   934,   941,   943,    -1,   964,   961,   989,   979,  1002,
     986,   987,   988,  1004,   995,  1005,  1007,  1006,  1011,  1019,
    1021,  1022,  1026,  1031,  1038,  1036,  1039,  1040,  1042,  1047,
    1055,  1062,  1052,  1064,  1068,  1069,  1076,  1070,  1077,  1078,
    1083,  1079,  1087,  1088,  1098,  1089,  1084,  1102,  1101,   165,
    1109,   568,  1110,  1111,  1117,   877,  1012,   457,   713,   875,
     978,   554,   565,   560,  1037,   741,   455,   611,   863,   835,
     413,   515,   699,   789,  1017,   883,   704,   333,     0,   724,
       0,     0,     0,     0,     0,   894
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        23,    24,   228,    26,    27,   140,   117,   100,   349,   350,
      33,   493,   111,   147,   111,    38,    39,    40,   111,   111,
     361,    44,   151,   111,   574,   154,    95,   111,   574,   635,
      95,   111,    95,   578,   708,   111,   555,   918,   557,   111,
     111,   111,   111,   562,   141,   142,   111,   344,   147,   141,
     142,   111,   796,   141,   142,   111,   149,   141,   142,   111,
     111,    51,   111,    51,   111,    33,    42,    51,    80,    42,
     411,   934,   935,   936,    22,    51,     0,    42,    51,    80,
      42,  1039,    80,     7,     8,    51,    51,   160,    31,    51,
      51,    51,    31,    17,    18,    19,    51,     6,   324,   618,
      24,   160,   161,    22,    28,    29,    30,    80,    32,    31,
      34,    35,    36,    37,   160,   137,   160,    42,    25,   137,
      42,    22,    22,    47,    42,    22,    50,    77,    33,    51,
      31,    31,    80,    51,    31,    44,   137,    80,    88,    63,
     162,    80,   161,    52,   162,    51,    51,   137,    22,   161,
      74,    75,    76,  1111,   104,   160,    80,    31,    80,   160,
     457,   137,    80,   161,    80,   138,    90,    91,    92,   137,
     160,    95,   120,    97,    98,    99,   100,   101,   166,    80,
      89,   165,   106,   107,  1047,   161,   527,   111,   112,   113,
     114,    31,   160,   100,   167,   161,   161,   121,   122,   161,
     124,   161,    42,  1084,   128,   160,    80,   131,   132,   133,
     134,    51,   131,   138,   139,   160,   138,   139,   127,   120,
     138,   139,   144,   145,   146,   147,   148,   149,   150,   151,
     152,   131,   138,   139,   131,    33,    22,   159,    22,   161,
      80,   163,   761,   161,   160,   764,   120,    31,   589,   137,
     160,   592,   593,   594,   595,   596,   597,   598,   599,   600,
     601,   137,   603,   604,   605,   606,   607,    22,   609,    80,
      80,    33,    80,  1017,   162,    31,     4,    80,    42,     4,
      42,   100,   160,    80,    80,   161,    42,    51,   137,   630,
     425,   426,    33,     6,    80,    51,    80,    80,   138,   139,
      80,   120,   437,    31,   144,   145,   146,   147,   148,   149,
     150,   151,   152,   162,   988,    31,    80,    80,    33,   159,
     661,   161,    80,   163,    80,    80,    42,    42,    80,    42,
     140,    44,   140,    80,   120,    51,   120,    80,    51,    52,
     138,   139,   140,   141,   140,   143,   160,    80,    80,   160,
     161,   161,    80,   161,   164,    80,   164,   160,   168,   169,
     168,   169,   160,   160,    80,   120,    80,    80,   164,    80,
     711,    80,   168,   169,   138,   160,    89,   160,   161,   137,
     160,   161,   138,   139,    56,    80,   137,   932,   144,   145,
     146,   147,   148,   149,   150,   151,   152,   160,   161,   918,
      80,   920,   160,   159,   137,   161,    78,   163,   160,   161,
      42,   162,    80,   160,   127,   137,   137,   160,   160,    51,
     939,    80,   138,   139,   517,   138,   137,   160,   144,   145,
     146,   147,   148,   149,   150,   151,   152,   137,   507,    33,
     162,   162,   507,   159,   507,   161,   160,   163,   547,   160,
     547,   160,   160,   161,   547,   547,    80,   137,   799,   547,
      80,    80,   162,   547,  1014,  1015,  1016,   547,  1014,  1015,
    1016,   547,    80,  1079,    80,   547,   547,   547,   547,    31,
     160,    80,   547,  1002,    60,  1004,  1005,   547,  1094,    80,
      80,   547,   160,    80,    80,   547,   547,   838,   547,   522,
     547,   160,   525,   638,   639,   640,   641,   642,   643,   644,
     645,   646,   647,   137,   649,   650,   651,   652,   653,  1038,
      80,  1040,  1072,    31,   161,   548,  1072,   161,    80,   137,
     664,   137,  1051,   874,    31,    80,   160,    31,    31,   668,
     160,   160,   671,   119,    80,  1064,   137,   137,    42,  1068,
      80,   137,   160,   129,   160,   137,   137,    51,   655,    80,
     137,     5,    80,   655,  1083,   664,    23,   655,   137,   160,
     160,   655,    80,   666,   160,    42,    20,   137,   160,    80,
      80,   162,   160,    80,    51,   162,    80,    80,   162,    24,
     168,  1110,   137,   162,   168,    39,    40,  1079,  1117,    43,
     160,   137,    46,   626,   945,   946,   137,   137,   631,    53,
      45,   137,  1094,    80,    42,   160,   137,    80,    80,    80,
     961,   160,    57,    51,   160,   160,   161,    80,    63,   168,
     160,   162,    80,   656,   657,    79,   162,   137,   137,   160,
     775,   982,    86,   160,   138,   139,   669,   160,   989,   672,
     144,   145,   146,   147,   148,   149,   150,   151,   152,   770,
     160,   105,   137,   162,   118,   159,  1007,   161,    80,   163,
       5,   138,   139,   117,   137,   137,   137,    42,  1019,   123,
    1021,   160,    80,    80,   137,    20,    21,   162,    80,   137,
      25,   135,   136,   160,   161,   160,   161,   160,   160,   160,
      42,    42,   928,    80,    39,    40,    41,   160,    43,    51,
      51,    46,   160,    42,   140,   141,    33,   143,    53,    54,
      55,  1062,    51,    58,    59,   137,    61,    62,    80,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,   137,
     137,    42,   160,   161,    79,    38,    81,    82,   160,    84,
      51,    86,  1093,   160,   161,   138,   139,   140,   141,    94,
     143,    42,   160,   160,   160,  1106,    80,   102,   103,    80,
     105,     9,    10,    11,    12,    13,    14,    15,   162,   162,
     126,    51,   117,   118,   168,     7,     8,    25,   123,   160,
     140,   141,   160,   143,   929,    17,    18,    19,   166,   165,
     135,   136,    24,    41,   126,   126,    28,    29,    30,    33,
      32,   126,    34,    35,    36,    37,    33,    33,   138,   139,
     140,   141,    33,   143,    33,    47,    33,   137,    50,     9,
      10,    11,    12,    13,    14,    15,    16,    33,    33,    33,
      20,    63,   162,    81,    33,    25,   138,   139,   140,   141,
      33,   143,    74,    75,    76,    33,    33,   162,    80,    33,
      33,    41,    33,    43,   102,   103,   137,   162,    90,    91,
      92,   162,    80,    95,    33,    97,    98,    99,   100,   101,
     160,   904,   160,   906,   106,   107,   137,    80,    80,   111,
     112,   113,   114,     3,    33,    80,    33,   135,   136,   121,
     122,    81,   124,    33,    33,    33,   128,    33,    33,   131,
     132,   133,   134,    33,    24,    42,    26,    27,    33,    33,
      33,    33,   102,   103,    51,   138,   139,   140,   141,    33,
     143,    33,   138,   139,   140,   141,    46,   143,    48,    49,
      33,    33,    42,    53,    33,    33,    80,    80,    80,   162,
      60,    51,    33,    80,   160,   135,   136,   138,   139,   140,
     141,    80,   143,   138,   139,   140,   141,   137,   143,    42,
     137,    80,    80,    83,    80,    85,    86,    87,    51,   161,
      80,   162,    33,   162,   160,    80,    96,   162,    42,    80,
     100,   138,   139,   140,   141,   161,   143,    51,   108,   109,
     110,   162,   162,   161,   114,   115,   116,    80,   118,   161,
     161,   138,   139,   160,    33,   125,   161,   144,   145,   146,
     147,   148,   149,   150,   151,   152,    80,    33,    33,    33,
     161,    42,   159,    33,   161,    33,   163,   137,   138,   139,
      51,    33,    33,   161,   144,   145,   146,   147,   148,   149,
     150,   151,   152,   138,   139,   140,   141,    33,   143,   159,
     161,   161,   161,   161,   137,   138,   139,   161,    33,    80,
      33,   144,   145,   146,   147,   148,   149,   150,   151,   152,
      33,   161,    33,   126,   138,   139,   159,    80,   161,   160,
     144,   145,   146,   147,   148,   149,   150,   151,   152,     3,
     161,   138,   139,   140,   141,   159,   143,   161,   138,   139,
     140,   141,   161,   143,   161,   161,   161,   161,    33,   161,
      24,   161,    26,    27,   126,   162,   161,   138,   139,   161,
     161,    80,   162,   144,   145,   146,   147,   148,   149,   150,
     151,   152,    46,   126,    48,    49,    80,   126,   159,    53,
     161,    80,     3,    42,    42,    51,    60,    51,    80,    51,
     138,   139,   140,   141,    51,   143,   138,   139,   140,   141,
      51,   143,    80,    24,    80,    26,    27,    51,    51,    83,
      42,    85,    86,    87,   162,    51,    51,   160,    33,    93,
     162,   161,    96,   160,    80,    46,   100,    48,    49,   160,
      80,    51,    53,   160,   108,   109,   110,   165,    80,    60,
     114,   115,   116,    42,   118,   138,   139,   140,   141,    42,
     143,    51,   138,   139,   140,   141,    51,   143,   138,   139,
     140,   141,    83,   143,    85,    86,    87,    51,    51,   162,
      51,    51,    51,    51,    51,    96,   162,    80,   161,   100,
      33,    42,   162,    80,   137,   137,    80,   108,   109,   110,
     137,   137,   137,   114,   115,   116,   137,   118,   138,   139,
     140,   141,    51,   143,   138,   139,   140,   141,    80,   143,
     138,   139,   140,   141,   143,   143,   138,   139,   140,   141,
      42,   143,   162,   160,   138,   139,   140,   141,   162,   143,
     138,   139,   140,   141,   162,   143,   138,   139,   140,   141,
     162,   143,   138,   139,   140,   141,    51,   143,   162,   160,
     138,   139,   140,   141,   162,   143,   138,   139,   140,   141,
     162,   143,   138,   139,   140,   141,   162,   143,   138,   139,
     140,   141,    51,   143,   162,   160,   138,   139,   140,   141,
     162,   143,   138,   139,   140,   141,   162,   143,   138,   139,
     140,   141,   162,   143,   138,   139,   140,   141,    51,   143,
     162,    51,   138,   139,   140,   141,   162,   143,   138,   139,
     140,   141,   162,   143,   138,   139,   140,   141,   162,   143,
     138,   139,   140,   141,   160,   143,   138,   139,   140,   141,
     160,   143,   138,   139,   140,   141,   160,   143,   138,   139,
     140,   141,   160,   143,   138,   139,   140,   141,   160,   143,
     138,   139,   140,   141,   160,   143,   138,   139,   140,   141,
     160,   143,   138,   139,   140,   141,   160,   143,   138,   139,
     140,   141,   160,   143,   138,   139,   140,   141,   160,   143,
     138,   139,   140,   141,   160,   143,    51,    51,    51,    51,
     160,    51,    80,   126,   143,   160,   160,   160,   160,    42,
     126,   126,   160,   126,    42,    42,    51,   164,   160,   160,
      51,   167,    80,   161,    80,   162,   162,   162,   162,   160,
     137,   137,    80,   137,    80,   160,    51,    51,   162,    33,
      33,   160,    80,    80,   143,   162,    80,    33,   160,   137,
     164,    80,   162,   137,    80,   137,    33,   160,    42,    33,
      33,    51,   160,    51,   137,   167,   137,   137,   160,   160,
     164,    33,   162,   137,   137,   160,    51,   160,   164,    51,
     137,   130,   162,    51,    80,   164,   137,    51,   160,    45,
     167,   322,   137,   137,   137,   666,   932,   144,   517,   664,
     869,   308,   319,   314,   996,   547,   142,   364,   655,   625,
     133,   211,   502,   583,   937,   671,   507,   122,    -1,   531,
      -1,    -1,    -1,    -1,    -1,   692
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,    17,    18,    19,    24,    28,    29,    30,
      32,    34,    35,    36,    37,    47,    50,    63,    74,    75,
      76,    80,    90,    91,    92,    95,    97,    98,    99,   100,
     101,   106,   107,   111,   112,   113,   114,   121,   122,   124,
     128,   131,   132,   133,   134,   171,   172,   173,   174,   175,
     176,   177,   178,   179,   184,   185,   186,   187,   190,   192,
     195,   200,   211,   212,   217,   220,   223,   226,   229,   235,
     241,   244,   249,   252,   255,   258,   259,   262,   263,   265,
     266,   267,   270,   271,   275,   276,   277,   278,   284,   287,
     293,   296,    51,   161,    51,   161,   160,   161,   160,   160,
     161,    33,    42,    51,    80,   161,    80,   161,   160,    80,
     160,   161,   232,   160,   160,   160,   160,   160,   161,    33,
      42,   160,   161,   161,   160,    33,   160,   160,   161,   232,
     232,   160,   160,   161,   232,   232,    80,   183,    33,    51,
     285,   161,   161,   232,   160,    33,   160,   161,   160,   161,
     160,   161,   232,   160,   161,   232,   232,    80,   180,    80,
     181,    80,   182,   232,     0,   172,   160,     9,    10,    11,
      12,    13,    14,    15,    25,    41,    81,   102,   103,   135,
     136,   290,   291,   292,   317,   318,   319,   320,   321,   353,
     354,   360,   361,   362,   363,   364,   365,   366,   160,    16,
      20,    43,   291,   294,   295,   325,   342,   367,    23,     4,
      80,   272,   273,   118,   224,   225,   298,    42,   160,    51,
     160,   160,   168,    80,   160,   168,    80,    80,   193,   194,
      33,     5,    21,    39,    40,    46,    53,    54,    55,    58,
      59,    61,    62,    64,    65,    66,    67,    68,    69,    70,
      71,    72,    73,    79,    82,    84,    86,    94,   105,   117,
     123,   250,   251,   298,   308,   317,   318,   319,   320,   321,
     322,   323,   324,   325,   326,   327,   328,   329,   330,   331,
     332,   333,   334,   335,   336,   337,   338,   339,   340,   341,
     342,   343,   344,   345,   347,   348,   353,   354,   355,   356,
     357,   358,    80,   137,   160,    22,    80,   120,   236,   237,
     238,    22,    80,   120,   245,   246,    22,    80,   120,   242,
     243,    80,   196,   197,   193,    38,   191,    42,   160,   201,
      60,   119,   129,   300,    77,    88,   104,   279,   280,   350,
     351,   352,    22,   131,   213,   214,    42,    51,    80,   138,
     139,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     159,   161,   188,    80,   260,   261,     3,    24,    26,    27,
      48,    49,    83,    85,    87,    93,    96,   100,   108,   109,
     110,   114,   115,   116,   231,   268,   269,   297,   298,   299,
     300,   301,   302,   303,   304,   305,   306,   307,   308,   309,
     310,   311,   312,   314,   315,   316,   324,   346,   349,   160,
     160,    80,   264,   269,   160,   160,   126,    80,   137,   160,
      51,   160,    42,    51,    80,   138,   139,   144,   145,   146,
     147,   148,   149,   150,   151,   152,   159,   161,   208,   210,
     253,   254,   308,   324,   325,   334,   340,   341,   342,   343,
     344,   345,   353,   354,   355,   253,   160,   213,   165,   227,
     228,   311,   317,   221,   222,   298,   230,   231,   160,   125,
     231,   288,   289,   359,   160,   160,   126,    80,   137,   160,
     126,    80,   137,   160,   126,    80,   137,   160,   160,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   137,   162,    33,    33,    33,   137,   162,   162,
      80,   137,   161,   274,    31,   273,    33,   137,   162,   160,
     160,    80,   162,   168,    80,   162,   168,    33,    31,   194,
      80,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,   137,   162,    33,
      80,    80,    80,    31,   237,   137,    80,   137,    80,    31,
     246,    80,   137,    80,    31,   243,   161,    31,   197,    31,
      33,   162,   160,   163,   206,   207,   208,   209,   137,   162,
     162,   162,    33,   137,   162,    80,    80,    31,   214,   161,
     188,   188,   161,   161,   161,   161,   161,   161,   161,   161,
     161,   161,   188,   138,   139,   140,   141,   143,   160,   161,
      31,   261,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,   137,   162,   137,   188,    31,
      80,   162,   126,    80,   160,   161,   208,   208,   161,   161,
     161,   161,   161,   161,   161,   161,   161,   161,   208,   138,
     139,   140,   141,   143,   286,   137,   162,   162,    31,    42,
      51,   161,   218,   219,   137,   162,   137,   162,   137,   162,
      33,   137,   162,   126,    80,   126,    80,   126,    80,    42,
      42,   138,   139,   234,    42,    51,   234,    51,    80,    51,
      51,   165,   370,   371,    51,    51,    80,    80,   368,   292,
      51,    51,    42,    51,   295,    51,   160,   161,    80,    42,
      51,    33,    51,   225,   160,   160,   160,   232,    80,   160,
     160,   232,    80,   188,   371,    51,    51,    51,    51,    51,
      42,    42,    51,    42,    51,    51,    51,    51,    80,   161,
      42,   251,   160,   232,    80,    33,   137,     6,    42,    44,
      51,    52,    80,    89,   127,   138,   239,   247,   248,   137,
     248,   137,   137,   248,   137,    51,   138,   139,   233,    80,
     160,    80,    31,   207,   209,    33,   160,    45,    57,    63,
     198,   199,   312,   313,   205,   160,   160,    56,    78,   280,
      80,   140,   164,   168,   169,   281,   282,   283,   137,    33,
     137,   160,   188,   189,   188,   188,   188,   188,   188,   188,
     188,   188,   188,   188,   162,   188,   188,   188,   188,   188,
     188,    51,    42,    51,    51,    51,    51,   248,    51,    51,
      42,    51,    51,    51,    51,   268,   160,   232,    80,   160,
     137,   188,   160,   232,   126,   233,   208,   208,   208,   208,
     208,   208,   208,   208,   208,   208,   162,   208,   208,   208,
     208,   208,   160,   254,   160,   232,   160,   232,   188,   160,
     166,    42,    51,   137,   161,   228,   160,   222,   160,   231,
     160,   232,    42,   289,   160,   232,   126,   126,   126,    42,
      42,    51,   369,   166,   369,   164,   160,   160,    51,   274,
     162,   162,   188,   160,   162,   160,   162,   160,   167,   256,
     257,   160,    80,    80,    42,    51,   160,   137,   137,    80,
     137,   248,    80,   160,   248,    51,    51,   162,   193,    33,
     208,    33,   137,   162,   160,   203,   202,   137,   160,   161,
     283,    80,   188,    80,   100,   120,   137,   162,   162,   162,
     162,   162,   162,   162,   162,   162,   162,   162,   162,   160,
     188,    80,   160,   160,   162,   162,   162,   162,   162,   162,
     162,   162,   162,   162,   162,   160,   160,   162,   219,   160,
      42,    51,   161,   188,   160,   160,   164,    80,   162,    33,
     160,   160,   232,   160,   232,    80,   137,   162,   240,   248,
     247,   248,   137,   248,   137,   137,   160,    33,    31,   208,
     160,    42,   199,   204,   206,   206,   206,   282,   248,    33,
     160,    33,    51,   215,   188,   188,   160,   160,   188,   188,
     162,    51,   274,   188,   160,   160,   167,   256,   137,   137,
     137,   248,   160,   248,   248,   188,   160,   160,    31,    31,
      31,   161,   162,   188,   188,   164,    51,   137,   160,   160,
     160,   162,    33,   160,   137,   248,   240,   248,   137,   160,
     160,   160,   206,   248,   160,   160,    51,   164,    51,   130,
     188,   167,   248,   137,   137,   248,    31,   162,    51,   164,
      80,   138,   139,   161,   216,   233,   234,   160,    80,   248,
     247,   160,    51,   188,    80,   160,   161,   233,   234,   167,
     137,   137,   162,   188,   248,   240,   162,   137,   248
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
      59,    40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   170,   171,   171,   172,   172,   172,   172,   172,   172,
     172,   172,   172,   172,   172,   172,   172,   172,   172,   172,
     172,   172,   172,   172,   172,   172,   172,   172,   172,   172,
     172,   172,   172,   172,   172,   172,   172,   172,   172,   172,
     172,   172,   172,   172,   172,   173,   173,   173,   173,   174,
     174,   175,   176,   177,   178,   179,   180,   180,   180,   180,
     180,   180,   181,   181,   181,   181,   181,   181,   182,   182,
     182,   182,   182,   182,   183,   183,   183,   183,   183,   183,
     184,   184,   185,   185,   186,   186,   187,   188,   188,   188,
     188,   188,   188,   188,   188,   188,   188,   188,   188,   188,
     188,   188,   188,   188,   188,   188,   188,   188,   188,   189,
     189,   190,   190,   191,   192,   193,   193,   194,   195,   196,
     196,   197,   198,   198,   199,   199,   199,   199,   201,   200,
     202,   200,   203,   200,   204,   200,   205,   200,   206,   206,
     206,   206,   207,   207,   208,   208,   208,   208,   208,   208,
     208,   208,   208,   208,   208,   208,   208,   208,   208,   208,
     208,   208,   208,   208,   208,   209,   210,   210,   211,   212,
     213,   213,   214,   214,   214,   214,   214,   215,   215,   215,
     215,   215,   215,   216,   216,   216,   216,   216,   216,   216,
     216,   217,   218,   218,   219,   219,   219,   219,   219,   219,
     219,   219,   219,   220,   220,   221,   221,   222,   223,   223,
     224,   224,   225,   226,   226,   227,   227,   228,   228,   229,
     229,   229,   229,   230,   230,   231,   231,   231,   231,   231,
     231,   231,   231,   231,   231,   231,   231,   231,   231,   231,
     231,   231,   231,   231,   231,   231,   232,   232,   232,   232,
     232,   232,   233,   233,   233,   234,   234,   234,   235,   236,
     236,   237,   238,   238,   238,   239,   239,   239,   239,   239,
     240,   240,   240,   240,   241,   242,   242,   243,   243,   243,
     244,   245,   245,   246,   246,   246,   247,   247,   247,   247,
     247,   248,   248,   248,   248,   248,   248,   249,   249,   249,
     249,   250,   250,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   252,   252,   253,   253,   254,   254,   254,   254,
     254,   254,   254,   254,   254,   254,   254,   254,   254,   255,
     255,   256,   256,   257,   257,   258,   259,   260,   260,   261,
     262,   263,   264,   264,   264,   264,   265,   266,   266,   266,
     266,   267,   267,   267,   267,   268,   268,   269,   269,   270,
     271,   272,   272,   273,   273,   273,   274,   274,   274,   275,
     275,   276,   276,   276,   276,   276,   276,   277,   277,   277,
     277,   277,   277,   278,   279,   279,   280,   280,   280,   281,
     281,   281,   281,   282,   282,   283,   283,   283,   283,   283,
     285,   286,   284,   287,   287,   287,   287,   288,   288,   289,
     289,   290,   290,   290,   290,   290,   290,   290,   291,   291,
     291,   291,   291,   291,   291,   291,   292,   292,   293,   293,
     294,   294,   294,   294,   295,   295,   296,   296,   297,   298,
     299,   300,   301,   302,   303,   304,   305,   306,   307,   308,
     309,   310,   311,   312,   313,   314,   315,   316,   316,   317,
     318,   318,   319,   320,   321,   322,   323,   324,   324,   325,
     326,   327,   328,   329,   330,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   350,   351,   352,   353,
     354,   355,   356,   357,   358,   359,   360,   361,   362,   363,
     364,   365,   366,   367,   368,   369,   369,   370,   370,   371
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
       1,     3,     3,     3,     3,     3,     2,     2,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     1,
       3,     4,     7,     3,     4,     2,     1,     4,     4,     2,
       1,     7,     3,     1,     1,     1,     1,     1,     0,     5,
       0,     8,     0,     8,     0,    10,     0,     8,     2,     2,
       1,     1,     4,     2,     3,     1,     1,     1,     3,     3,
       3,     3,     3,     2,     2,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     5,     1,     4,     4,     4,
       2,     1,     9,     6,     5,     7,     7,     2,     4,     3,
       5,     3,     1,     2,     2,     2,     1,     1,     1,     4,
       3,     6,     3,     1,     5,     3,     3,     4,     2,     2,
       3,     1,     1,     2,     5,     3,     1,     1,     2,     5,
       3,     1,     1,     2,     5,     3,     1,     1,     1,     2,
       5,     3,     6,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     2,     4,     3,     5,
       1,     3,     2,     2,     1,     2,     2,     1,     4,     2,
       1,     4,     2,     1,     4,     3,     5,     9,     1,     5,
       3,     5,     7,     9,     4,     2,     1,     5,     7,     4,
       4,     2,     1,     7,     9,     6,     1,     1,     1,     1,
       1,     0,     1,     1,     1,     2,     2,     2,     5,     3,
       6,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     5,     6,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     5,
       6,     7,     5,     1,     3,     3,     4,     2,     1,     5,
       3,     4,     4,     6,     3,     5,     3,     2,     5,     3,
       6,     2,     5,     3,     6,     1,     1,     1,     3,     3,
       4,     2,     1,     5,     7,     9,     0,     3,     3,     2,
       5,     5,     6,     3,     7,     8,     5,     5,     6,     3,
       7,     8,     5,     6,     3,     1,     1,     1,     1,     1,
       3,     4,     6,     1,     2,     1,     1,     1,     1,     1,
       0,     0,     5,     2,     5,     3,     6,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     1,     3,     6,
       1,     1,     1,     1,     3,     1,     3,     6,     3,     3,
       3,     1,     3,     3,     3,     3,     1,     1,     1,     3,
       3,     3,     3,     3,     3,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     3,
       3,     3,     3,     5,     3,     3,     3,     3,     1,     3,
       3,     3,     1,     1,     1,     1,     1,     3,     1,     1,
       1,     1,     3,     3,     3,     3,     3,     1,     1,     3,
       3,     3,     1,     1,     1,     3,     3,     3,     3,     3,
       3,     1,     3,     3,     3,     1,     3,     2,     2,     2
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
  "NOPRINT", "NORMAL_PDF", "OBSERVATION_TRENDS", "OLR", "OLR_INST",
  "OLR_BETA", "OPTIM", "OPTIM_WEIGHTS", "ORDER", "OSR", "OSR_PARAMS",
  "PARAMETERS", "PERIODS", "PLANNER_OBJECTIVE", "PREFILTER", "PRESAMPLE",
  "PRINT", "PRIOR_TRUNC", "PRIOR_ANALYSIS", "POSTERIOR_ANALYSIS",
  "QZ_CRITERIUM", "RELATIVE_IRF", "REPLIC", "RPLOT", "SHOCKS", "SIGMA_E",
  "SIMUL", "SIMUL_ALGO", "SIMUL_SEED", "SMOOTHER", "SOLVE_ALGO",
  "SPARSE_DLL", "STDERR", "STEADY", "STOCH_SIMUL", "TEX", "RAMSEY_POLICY",
  "PLANNER_DISCOUNT", "TEX_NAME", "UNIFORM_PDF", "UNIT_ROOT_VARS",
  "USE_DLL", "VALUES", "VAR", "VAREXO", "VAREXO_DET", "VAROBS",
  "XLS_SHEET", "XLS_RANGE", "COMMA", "MINUS", "PLUS", "DIVIDE", "TIMES",
  "UMINUS", "POWER", "EXP", "LOG", "LOG10", "SIN", "COS", "TAN", "ASIN",
  "ACOS", "ATAN", "SINH", "COSH", "TANH", "ASINH", "ACOSH", "ATANH",
  "SQRT", "';'", "'('", "')'", "'#'", "':'", "'['", "']'", "'''", "'.'",
  "'\\\\'", "$accept", "statement_list", "statement", "declaration",
  "dsample", "rplot", "var", "varexo", "varexo_det", "parameters",
  "var_list", "varexo_list", "varexo_det_list", "parameter_list",
  "periods", "cutoff", "markowitz", "init_param", "expression",
  "comma_expression", "initval", "initval_option", "endval",
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
  "osr_params", "osr", "olr", "olr_option", "olr_options", "olr_inst",
  "calib_var", "calib_var_list", "calib_arg1", "calib_arg2", "calib",
  "dynatype", "dynasave", "model_comparison", "model_comparison_options",
  "model_comparison_option", "filename_list", "filename", "filename_elem",
  "planner_objective", "@6", "@7", "ramsey_policy",
  "ramsey_policy_options_list", "ramsey_policy_options",
  "bvar_prior_option", "bvar_common_option", "bvar_density_options_list",
  "bvar_density", "bvar_forecast_option", "bvar_forecast_options_list",
  "bvar_forecast", "o_dr_algo", "o_solve_algo", "o_simul_algo", "o_linear",
  "o_order", "o_replic", "o_drop", "o_ar", "o_nocorr", "o_nofunctions",
  "o_nomoments", "o_irf", "o_hp_filter", "o_hp_ngrid", "o_periods",
  "o_cutoff", "o_markowitz", "o_simul", "o_simul_seed", "o_qz_criterium",
  "o_datafile", "o_nobs", "o_first_obs", "o_prefilter", "o_presample",
  "o_lik_algo", "o_lik_init", "o_nograph", "o_conf_sig", "o_mh_replic",
  "o_mh_drop", "o_mh_jscale", "o_optim", "o_mh_init_scale", "o_mode_file",
  "o_mode_compute", "o_mode_check", "o_prior_trunc", "o_mh_mode",
  "o_mh_nblcks", "o_load_mh_file", "o_loglinear", "o_nodiagnostic",
  "o_bayesian_irf", "o_tex", "o_forecast", "o_smoother",
  "o_moments_varendo", "o_filtered_vars", "o_relative_irf",
  "o_kalman_algo", "o_kalman_tol", "o_olr_beta",
  "o_model_comparison_approximation", "o_print", "o_noprint",
  "o_xls_sheet", "o_xls_range", "o_filter_step_ahead", "o_constant",
  "o_noconstant", "o_mh_recover", "o_planner_discount", "o_bvar_prior_tau",
  "o_bvar_prior_decay", "o_bvar_prior_lambda", "o_bvar_prior_mu",
  "o_bvar_prior_omega", "o_bvar_prior_flat", "o_bvar_prior_train",
  "o_bvar_replic", "range", "vec_int_elem", "vec_int_1", "vec_int", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const parser::rhs_number_type
  parser::yyrhs_[] =
  {
       171,     0,    -1,   172,    -1,   171,   172,    -1,   173,    -1,
     184,    -1,   185,    -1,   186,    -1,   200,    -1,   190,    -1,
     192,    -1,   195,    -1,   187,    -1,   211,    -1,   212,    -1,
     217,    -1,   220,    -1,   223,    -1,   226,    -1,   229,    -1,
     249,    -1,   252,    -1,   255,    -1,   235,    -1,   244,    -1,
     241,    -1,   258,    -1,   259,    -1,   262,    -1,   174,    -1,
     175,    -1,   263,    -1,   265,    -1,   266,    -1,   271,    -1,
     275,    -1,   276,    -1,   277,    -1,   267,    -1,   270,    -1,
     278,    -1,   284,    -1,   287,    -1,   293,    -1,   296,    -1,
     179,    -1,   176,    -1,   177,    -1,   178,    -1,    28,    51,
     160,    -1,    28,    51,    51,   160,    -1,   111,   232,   160,
      -1,   131,   180,   160,    -1,   132,   181,   160,    -1,   133,
     182,   160,    -1,    99,   183,   160,    -1,   180,    80,    -1,
     180,   137,    80,    -1,    80,    -1,   180,    80,   126,    -1,
     180,   137,    80,   126,    -1,    80,   126,    -1,   181,    80,
      -1,   181,   137,    80,    -1,    80,    -1,   181,    80,   126,
      -1,   181,   137,    80,   126,    -1,    80,   126,    -1,   182,
      80,    -1,   182,   137,    80,    -1,    80,    -1,   182,    80,
     126,    -1,   182,   137,    80,   126,    -1,    80,   126,    -1,
     183,    80,    -1,   183,   137,    80,    -1,    80,    -1,   183,
      80,   126,    -1,   183,   137,    80,   126,    -1,    80,   126,
      -1,   100,    51,   160,    -1,   100,    33,    51,   160,    -1,
      24,    42,   160,    -1,    24,    33,    42,   160,    -1,    63,
      42,   160,    -1,    63,    33,    42,   160,    -1,    80,    33,
     188,   160,    -1,   161,   188,   162,    -1,    80,    -1,    42,
      -1,    51,    -1,   188,   139,   188,    -1,   188,   138,   188,
      -1,   188,   140,   188,    -1,   188,   141,   188,    -1,   188,
     143,   188,    -1,   138,   188,    -1,   139,   188,    -1,   144,
     161,   188,   162,    -1,   145,   161,   188,   162,    -1,   146,
     161,   188,   162,    -1,   147,   161,   188,   162,    -1,   148,
     161,   188,   162,    -1,   149,   161,   188,   162,    -1,   150,
     161,   188,   162,    -1,   151,   161,   188,   162,    -1,   152,
     161,   188,   162,    -1,   159,   161,   188,   162,    -1,    80,
     161,   189,   162,    -1,   188,    -1,   189,   137,   188,    -1,
      50,   160,   193,    31,    -1,    50,   161,   191,   162,   160,
     193,    31,    -1,    38,    33,    80,    -1,    32,   160,   193,
      31,    -1,   193,   194,    -1,   194,    -1,    80,    33,   188,
     160,    -1,    47,   160,   196,    31,    -1,   196,   197,    -1,
     197,    -1,    80,   161,   233,   162,    33,   188,   160,    -1,
     198,   137,   199,    -1,   199,    -1,    57,    -1,    45,    -1,
     312,    -1,   313,    -1,    -1,    74,   160,   201,   206,    31,
      -1,    -1,    74,   161,   300,   162,   160,   202,   206,    31,
      -1,    -1,    74,   161,   129,   162,   160,   203,   206,    31,
      -1,    -1,    74,   161,   119,   137,   198,   162,   204,   160,
     206,    31,    -1,    -1,    74,   161,   119,   162,   205,   160,
     206,    31,    -1,   206,   207,    -1,   206,   209,    -1,   207,
      -1,   209,    -1,   208,    33,   208,   160,    -1,   208,   160,
      -1,   161,   208,   162,    -1,   210,    -1,    42,    -1,    51,
      -1,   208,   139,   208,    -1,   208,   138,   208,    -1,   208,
     140,   208,    -1,   208,   141,   208,    -1,   208,   143,   208,
      -1,   138,   208,    -1,   139,   208,    -1,   144,   161,   208,
     162,    -1,   145,   161,   208,   162,    -1,   146,   161,   208,
     162,    -1,   147,   161,   208,   162,    -1,   148,   161,   208,
     162,    -1,   149,   161,   208,   162,    -1,   150,   161,   208,
     162,    -1,   151,   161,   208,   162,    -1,   152,   161,   208,
     162,    -1,   159,   161,   208,   162,    -1,   163,    80,    33,
     208,   160,    -1,    80,    -1,    80,   161,   233,   162,    -1,
     112,   160,   213,    31,    -1,    76,   160,   213,    31,    -1,
     213,   214,    -1,   214,    -1,   131,    80,   160,   100,   215,
     160,   130,   216,   160,    -1,   131,    80,   160,   120,   188,
     160,    -1,   131,    80,    33,   188,   160,    -1,   131,    80,
     137,    80,    33,   188,   160,    -1,    22,    80,   137,    80,
      33,   188,   160,    -1,   215,    51,    -1,   215,    51,   164,
      51,    -1,   215,   137,    51,    -1,   215,   137,    51,   164,
      51,    -1,    51,   164,    51,    -1,    51,    -1,   216,   234,
      -1,   216,   233,    -1,   216,    80,    -1,   234,    -1,   233,
      -1,    80,    -1,   216,   161,   188,   162,    -1,   161,   188,
     162,    -1,   113,    33,   165,   218,   166,   160,    -1,   218,
     160,   219,    -1,   219,    -1,   219,   137,   161,   188,   162,
      -1,   219,   137,    42,    -1,   219,   137,    51,    -1,   219,
     161,   188,   162,    -1,   219,    42,    -1,   219,    51,    -1,
     161,   188,   162,    -1,    42,    -1,    51,    -1,   121,   160,
      -1,   121,   161,   221,   162,   160,    -1,   221,   137,   222,
      -1,   222,    -1,   298,    -1,    19,   160,    -1,    19,   161,
     224,   162,   160,    -1,   224,   137,   225,    -1,   225,    -1,
     298,    -1,   114,   160,    -1,   114,   161,   227,   162,   160,
      -1,   227,   137,   228,    -1,   228,    -1,   311,    -1,   317,
      -1,   122,   160,    -1,   122,   161,   230,   162,   160,    -1,
     122,   232,   160,    -1,   122,   161,   230,   162,   232,   160,
      -1,   230,   137,   231,    -1,   231,    -1,   297,    -1,   298,
      -1,   299,    -1,   300,    -1,   301,    -1,   302,    -1,   303,
      -1,   304,    -1,   305,    -1,   306,    -1,   307,    -1,   324,
      -1,   308,    -1,   346,    -1,   309,    -1,   310,    -1,   311,
      -1,   312,    -1,   314,    -1,   315,    -1,   316,    -1,   232,
      80,    -1,   232,    80,    33,    80,    -1,   232,   137,    80,
      -1,   232,   137,    80,    33,    80,    -1,    80,    -1,    80,
      33,    80,    -1,   139,    51,    -1,   138,    51,    -1,    51,
      -1,   139,    42,    -1,   138,    42,    -1,    42,    -1,    35,
     160,   236,    31,    -1,   236,   237,    -1,   237,    -1,   238,
     137,   239,   160,    -1,   120,    80,    -1,    80,    -1,    22,
      80,   137,    80,    -1,   247,   137,   240,    -1,   248,   137,
     247,   137,   240,    -1,   248,   137,   248,   137,   248,   137,
     247,   137,   240,    -1,   248,    -1,   248,   137,   248,   137,
     248,    -1,   248,   137,   248,    -1,   248,   137,   248,   137,
     248,    -1,   248,   137,   248,   137,   248,   137,   248,    -1,
     248,   137,   248,   137,   248,   137,   248,   137,   248,    -1,
      37,   160,   242,    31,    -1,   242,   243,    -1,   243,    -1,
     120,    80,   137,   248,   160,    -1,    22,    80,   137,    80,
     137,   248,   160,    -1,    80,   137,   248,   160,    -1,    36,
     160,   245,    31,    -1,   245,   246,    -1,   246,    -1,   120,
      80,   137,   248,   137,   248,   160,    -1,    22,    80,   137,
      80,   137,   248,   137,   248,   160,    -1,    80,   137,   248,
     137,   248,   160,    -1,     6,    -1,    44,    -1,    89,    -1,
      52,    -1,   127,    -1,    -1,    51,    -1,    42,    -1,    80,
      -1,   138,    51,    -1,   138,    42,    -1,    34,   160,    -1,
      34,   161,   250,   162,   160,    -1,    34,   232,   160,    -1,
      34,   161,   250,   162,   232,   160,    -1,   250,   137,   251,
      -1,   251,    -1,   317,    -1,   318,    -1,   319,    -1,   320,
      -1,   321,    -1,   322,    -1,   323,    -1,   324,    -1,   325,
      -1,   326,    -1,   327,    -1,   328,    -1,   329,    -1,   330,
      -1,   331,    -1,   332,    -1,   333,    -1,   334,    -1,   335,
      -1,   336,    -1,   337,    -1,   338,    -1,   339,    -1,   340,
      -1,   308,    -1,   341,    -1,   342,    -1,   343,    -1,   344,
      -1,   345,    -1,   347,    -1,   348,    -1,   353,    -1,   354,
      -1,   355,    -1,   298,    -1,   356,    -1,   357,    -1,   358,
      -1,   106,   161,   253,   162,   160,    -1,   106,   161,   253,
     162,   232,   160,    -1,   253,   137,   254,    -1,   254,    -1,
     324,    -1,   325,    -1,   334,    -1,   340,    -1,   308,    -1,
     341,    -1,   342,    -1,   343,    -1,   344,    -1,   345,    -1,
     353,    -1,   354,    -1,   355,    -1,   107,   161,   253,   162,
     160,    -1,   107,   161,   253,   162,   232,   160,    -1,   167,
      80,   167,   137,   167,    80,   167,    -1,   167,    80,   167,
     137,   248,    -1,   256,    -1,   257,   137,   256,    -1,   134,
     232,   160,    -1,    90,   160,   260,    31,    -1,   260,   261,
      -1,   261,    -1,    80,   161,   188,   162,   160,    -1,   128,
     232,   160,    -1,    95,   160,   264,    31,    -1,   264,    80,
     188,   160,    -1,   264,    80,   137,    80,   188,   160,    -1,
      80,   188,   160,    -1,    80,   137,    80,   188,   160,    -1,
      98,   232,   160,    -1,    97,   160,    -1,    97,   161,   269,
     162,   160,    -1,    97,   232,   160,    -1,    97,   161,   269,
     162,   232,   160,    -1,    91,   160,    -1,    91,   161,   269,
     162,   160,    -1,    91,   232,   160,    -1,    91,   161,   269,
     162,   232,   160,    -1,   349,    -1,   231,    -1,   268,    -1,
     269,   137,   268,    -1,    92,   232,   160,    -1,    18,   160,
     272,    31,    -1,   272,   273,    -1,   273,    -1,    80,   274,
      33,   188,   160,    -1,    80,   137,    80,   274,    33,   188,
     160,    -1,     4,    80,   161,    51,   162,   274,    33,   188,
     160,    -1,    -1,   161,    51,   162,    -1,   161,    42,   162,
      -1,    17,   160,    -1,    17,   161,    23,   162,   160,    -1,
      30,   161,    80,   162,   160,    -1,    30,   161,    80,   162,
     232,   160,    -1,    30,    80,   160,    -1,    30,   161,    80,
     168,    80,   162,   160,    -1,    30,   161,    80,   168,    80,
     162,   232,   160,    -1,    30,    80,   168,    80,   160,    -1,
      29,   161,    80,   162,   160,    -1,    29,   161,    80,   162,
     232,   160,    -1,    29,    80,   160,    -1,    29,   161,    80,
     168,    80,   162,   160,    -1,    29,   161,    80,   168,    80,
     162,   232,   160,    -1,    29,    80,   168,    80,   160,    -1,
      75,   161,   279,   162,   281,   160,    -1,   279,   137,   280,
      -1,   280,    -1,   350,    -1,   351,    -1,   352,    -1,   282,
      -1,   281,   137,   282,    -1,   282,   161,   248,   162,    -1,
     281,   137,   282,   161,   248,   162,    -1,   283,    -1,   282,
     283,    -1,    80,    -1,   169,    -1,   140,    -1,   164,    -1,
     168,    -1,    -1,    -1,   101,   285,   208,   286,   160,    -1,
     124,   160,    -1,   124,   161,   288,   162,   160,    -1,   124,
     232,   160,    -1,   124,   161,   288,   162,   232,   160,    -1,
     288,   137,   289,    -1,   289,    -1,   231,    -1,   359,    -1,
     360,    -1,   361,    -1,   362,    -1,   363,    -1,   364,    -1,
     365,    -1,   366,    -1,   290,    -1,   317,    -1,   353,    -1,
     354,    -1,   319,    -1,   321,    -1,   318,    -1,   320,    -1,
     291,   137,   292,    -1,   291,    -1,     7,    51,   160,    -1,
       7,   161,   292,   162,    51,   160,    -1,   291,    -1,   342,
      -1,   325,    -1,   367,    -1,   294,   137,   295,    -1,   294,
      -1,     8,    51,   160,    -1,     8,   161,   295,   162,    51,
     160,    -1,    26,    33,    51,    -1,   118,    33,    51,    -1,
     115,    33,    51,    -1,    60,    -1,    96,    33,    51,    -1,
     110,    33,    51,    -1,    27,    33,    51,    -1,     3,    33,
      51,    -1,    83,    -1,    85,    -1,    87,    -1,    53,    33,
      51,    -1,    48,    33,    51,    -1,    49,    33,    51,    -1,
     100,    33,    51,    -1,    24,    33,    42,    -1,    63,    33,
      42,    -1,   114,    -1,   116,    33,    51,    -1,   108,    33,
      51,    -1,   108,    33,    42,    -1,    25,    33,    80,    -1,
      81,    33,   371,    -1,    81,    33,    51,    -1,    41,    33,
      51,    -1,   102,    33,    51,    -1,   103,    33,    51,    -1,
      58,    33,    51,    -1,    59,    33,    51,    -1,    86,    -1,
      46,    -1,    20,    33,    42,    -1,    69,    33,    51,    -1,
      64,    33,    42,    -1,    66,    33,    42,    -1,    94,    33,
     161,   257,   162,    -1,    65,    33,    42,    -1,    65,    33,
      51,    -1,    73,    33,    80,    -1,    72,    33,    51,    -1,
      71,    -1,   105,    33,    42,    -1,    67,    33,    51,    -1,
      68,    33,    51,    -1,    61,    -1,    62,    -1,    84,    -1,
       5,    -1,   123,    -1,    43,    33,    51,    -1,   117,    -1,
      79,    -1,    40,    -1,   109,    -1,    54,    33,    51,    -1,
      55,    33,    51,    -1,    93,    33,   248,    -1,    77,    33,
      56,    -1,    77,    33,    78,    -1,   104,    -1,    88,    -1,
     135,    33,    80,    -1,   136,    33,   368,    -1,    39,    33,
     371,    -1,    21,    -1,    82,    -1,    70,    -1,   125,    33,
      42,    -1,    14,    33,   234,    -1,     9,    33,    42,    -1,
      11,    33,   234,    -1,    12,    33,    42,    -1,    13,    33,
      51,    -1,    10,    -1,    15,    33,    51,    -1,    16,    33,
      51,    -1,    80,   164,    80,    -1,    51,    -1,    51,   164,
      51,    -1,   165,   369,    -1,   370,   369,    -1,   370,   166,
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
      80,    82,    84,    86,    88,    90,    92,    94,    96,    98,
     102,   107,   111,   115,   119,   123,   127,   130,   134,   136,
     140,   145,   148,   151,   155,   157,   161,   166,   169,   172,
     176,   178,   182,   187,   190,   193,   197,   199,   203,   208,
     211,   215,   220,   224,   229,   233,   238,   243,   247,   249,
     251,   253,   257,   261,   265,   269,   273,   276,   279,   284,
     289,   294,   299,   304,   309,   314,   319,   324,   329,   334,
     336,   340,   345,   353,   357,   362,   365,   367,   372,   377,
     380,   382,   390,   394,   396,   398,   400,   402,   404,   405,
     411,   412,   421,   422,   431,   432,   443,   444,   453,   456,
     459,   461,   463,   468,   471,   475,   477,   479,   481,   485,
     489,   493,   497,   501,   504,   507,   512,   517,   522,   527,
     532,   537,   542,   547,   552,   557,   563,   565,   570,   575,
     580,   583,   585,   595,   602,   608,   616,   624,   627,   632,
     636,   642,   646,   648,   651,   654,   657,   659,   661,   663,
     668,   672,   679,   683,   685,   691,   695,   699,   704,   707,
     710,   714,   716,   718,   721,   727,   731,   733,   735,   738,
     744,   748,   750,   752,   755,   761,   765,   767,   769,   771,
     774,   780,   784,   791,   795,   797,   799,   801,   803,   805,
     807,   809,   811,   813,   815,   817,   819,   821,   823,   825,
     827,   829,   831,   833,   835,   837,   839,   842,   847,   851,
     857,   859,   863,   866,   869,   871,   874,   877,   879,   884,
     887,   889,   894,   897,   899,   904,   908,   914,   924,   926,
     932,   936,   942,   950,   960,   965,   968,   970,   976,   984,
     989,   994,   997,   999,  1007,  1017,  1024,  1026,  1028,  1030,
    1032,  1034,  1035,  1037,  1039,  1041,  1044,  1047,  1050,  1056,
    1060,  1067,  1071,  1073,  1075,  1077,  1079,  1081,  1083,  1085,
    1087,  1089,  1091,  1093,  1095,  1097,  1099,  1101,  1103,  1105,
    1107,  1109,  1111,  1113,  1115,  1117,  1119,  1121,  1123,  1125,
    1127,  1129,  1131,  1133,  1135,  1137,  1139,  1141,  1143,  1145,
    1147,  1149,  1151,  1157,  1164,  1168,  1170,  1172,  1174,  1176,
    1178,  1180,  1182,  1184,  1186,  1188,  1190,  1192,  1194,  1196,
    1202,  1209,  1217,  1223,  1225,  1229,  1233,  1238,  1241,  1243,
    1249,  1253,  1258,  1263,  1270,  1274,  1280,  1284,  1287,  1293,
    1297,  1304,  1307,  1313,  1317,  1324,  1326,  1328,  1330,  1334,
    1338,  1343,  1346,  1348,  1354,  1362,  1372,  1373,  1377,  1381,
    1384,  1390,  1396,  1403,  1407,  1415,  1424,  1430,  1436,  1443,
    1447,  1455,  1464,  1470,  1477,  1481,  1483,  1485,  1487,  1489,
    1491,  1495,  1500,  1507,  1509,  1512,  1514,  1516,  1518,  1520,
    1522,  1523,  1524,  1530,  1533,  1539,  1543,  1550,  1554,  1556,
    1558,  1560,  1562,  1564,  1566,  1568,  1570,  1572,  1574,  1576,
    1578,  1580,  1582,  1584,  1586,  1588,  1590,  1594,  1596,  1600,
    1607,  1609,  1611,  1613,  1615,  1619,  1621,  1625,  1632,  1636,
    1640,  1644,  1646,  1650,  1654,  1658,  1662,  1664,  1666,  1668,
    1672,  1676,  1680,  1684,  1688,  1692,  1694,  1698,  1702,  1706,
    1710,  1714,  1718,  1722,  1726,  1730,  1734,  1738,  1740,  1742,
    1746,  1750,  1754,  1758,  1764,  1768,  1772,  1776,  1780,  1782,
    1786,  1790,  1794,  1796,  1798,  1800,  1802,  1804,  1808,  1810,
    1812,  1814,  1816,  1820,  1824,  1828,  1832,  1836,  1838,  1840,
    1844,  1848,  1852,  1854,  1856,  1858,  1862,  1866,  1870,  1874,
    1878,  1882,  1884,  1888,  1892,  1896,  1898,  1902,  1905,  1908
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    88,    88,    89,    93,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
     119,   120,   121,   122,   123,   124,   125,   126,   127,   128,
     129,   130,   131,   132,   133,   138,   139,   140,   141,   145,
     146,   149,   152,   156,   160,   164,   168,   170,   172,   174,
     176,   178,   183,   185,   187,   189,   191,   193,   198,   200,
     202,   204,   206,   208,   213,   215,   217,   219,   221,   223,
     228,   232,   239,   243,   250,   254,   262,   267,   269,   271,
     273,   275,   277,   279,   281,   283,   285,   287,   289,   291,
     293,   295,   297,   299,   301,   303,   305,   307,   309,   314,
     316,   320,   322,   327,   331,   336,   337,   341,   346,   351,
     352,   356,   360,   361,   365,   366,   367,   368,   372,   372,
     373,   373,   375,   375,   377,   377,   379,   379,   384,   385,
     386,   387,   391,   393,   398,   399,   400,   402,   404,   406,
     408,   410,   412,   414,   416,   418,   420,   422,   424,   426,
     428,   430,   432,   434,   436,   440,   444,   446,   451,   455,
     459,   460,   464,   466,   468,   470,   472,   477,   479,   481,
     483,   485,   487,   493,   495,   497,   499,   501,   503,   505,
     507,   512,   517,   519,   524,   526,   528,   530,   532,   534,
     536,   538,   540,   545,   549,   553,   554,   557,   561,   563,
     567,   568,   571,   575,   577,   581,   582,   585,   586,   590,
     592,   594,   596,   600,   601,   604,   605,   606,   607,   608,
     609,   610,   611,   612,   613,   614,   615,   616,   617,   618,
     619,   620,   621,   622,   623,   624,   628,   630,   632,   634,
     636,   638,   643,   645,   647,   652,   654,   656,   661,   666,
     668,   673,   677,   682,   687,   697,   702,   708,   718,   723,
     734,   740,   748,   758,   772,   776,   778,   782,   789,   798,
     807,   811,   813,   817,   826,   837,   849,   851,   853,   855,
     857,   862,   863,   864,   865,   866,   868,   875,   877,   879,
     881,   886,   887,   890,   891,   892,   893,   894,   895,   896,
     897,   898,   899,   900,   901,   902,   903,   904,   905,   906,
     907,   908,   909,   910,   911,   912,   913,   914,   915,   916,
     917,   918,   919,   920,   921,   922,   923,   924,   925,   926,
     927,   928,   932,   934,   939,   940,   944,   945,   946,   947,
     948,   949,   950,   951,   952,   953,   954,   955,   956,   960,
     962,   967,   968,   972,   973,   977,   982,   987,   988,   991,
     995,   998,  1002,  1004,  1006,  1008,  1012,  1015,  1016,  1017,
    1018,  1021,  1022,  1023,  1024,  1027,  1028,  1031,  1032,  1035,
    1038,  1042,  1043,  1046,  1047,  1048,  1051,  1052,  1053,  1056,
    1057,  1060,  1061,  1062,  1063,  1064,  1065,  1067,  1068,  1069,
    1070,  1071,  1072,  1074,  1078,  1079,  1082,  1083,  1084,  1087,
    1088,  1089,  1090,  1093,  1094,  1097,  1098,  1099,  1100,  1101,
    1104,  1104,  1104,  1107,  1109,  1111,  1113,  1118,  1119,  1122,
    1123,  1126,  1127,  1128,  1129,  1130,  1131,  1132,  1135,  1136,
    1137,  1138,  1139,  1140,  1141,  1142,  1145,  1146,  1149,  1151,
    1155,  1156,  1157,  1158,  1161,  1162,  1165,  1167,  1171,  1172,
    1173,  1174,  1175,  1176,  1177,  1178,  1179,  1180,  1181,  1182,
    1183,  1184,  1185,  1186,  1187,  1188,  1189,  1190,  1191,  1193,
    1194,  1195,  1197,  1198,  1199,  1200,  1201,  1202,  1203,  1204,
    1205,  1206,  1207,  1208,  1209,  1210,  1211,  1212,  1213,  1214,
    1215,  1216,  1217,  1218,  1219,  1220,  1221,  1222,  1223,  1224,
    1225,  1226,  1227,  1228,  1229,  1231,  1233,  1236,  1237,  1238,
    1239,  1240,  1241,  1242,  1243,  1244,  1246,  1247,  1248,  1249,
    1250,  1251,  1252,  1253,  1255,  1263,  1264,  1268,  1269,  1278
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
       2,     2,     2,     2,     2,   163,     2,     2,     2,   167,
     161,   162,     2,     2,     2,     2,   168,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   164,   160,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   165,   169,   166,     2,     2,     2,     2,     2,     2,
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
     155,   156,   157,   158,   159
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1585;
  const int parser::yynnts_ = 202;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 164;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 170;

  const unsigned int parser::yyuser_token_number_max_ = 414;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1280 "DynareBison.yy"


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

