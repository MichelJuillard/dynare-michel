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

  case 248:
#line 631 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 249:
#line 633 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 250:
#line 635 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 251:
#line 637 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 252:
#line 639 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 253:
#line 641 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 254:
#line 646 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 255:
#line 648 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 256:
#line 650 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 257:
#line 655 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 258:
#line 657 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 259:
#line 659 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 260:
#line 664 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 261:
#line 669 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 262:
#line 671 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 264:
#line 680 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 265:
#line 685 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 266:
#line 690 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 267:
#line 700 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 268:
#line 705 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 269:
#line 711 "DynareBison.yy"
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

  case 270:
#line 721 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 271:
#line 726 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 272:
#line 737 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 273:
#line 743 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 274:
#line 751 "DynareBison.yy"
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

  case 275:
#line 761 "DynareBison.yy"
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

  case 276:
#line 775 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 277:
#line 779 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 278:
#line 781 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 279:
#line 785 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 280:
#line 792 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 281:
#line 801 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 282:
#line 810 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 283:
#line 814 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 284:
#line 816 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 285:
#line 820 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 286:
#line 829 "DynareBison.yy"
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

  case 287:
#line 840 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 288:
#line 852 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 289:
#line 854 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 290:
#line 856 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 291:
#line 858 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 292:
#line 860 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 293:
#line 864 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 297:
#line 869 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 298:
#line 871 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 299:
#line 878 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 300:
#line 880 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 301:
#line 882 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 302:
#line 884 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 344:
#line 935 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 345:
#line 937 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 361:
#line 963 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 362:
#line 965 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 363:
#line 969 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 364:
#line 970 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 367:
#line 980 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 368:
#line 985 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 371:
#line 994 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val));;}
    break;

  case 372:
#line 997 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 373:
#line 1001 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 374:
#line 1005 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 375:
#line 1007 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 376:
#line 1009 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 377:
#line 1011 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 378:
#line 1014 "DynareBison.yy"
    {driver.set_osr_params();;}
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
#line 1019 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 382:
#line 1020 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 383:
#line 1023 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 384:
#line 1024 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 385:
#line 1025 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 386:
#line 1026 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 391:
#line 1037 "DynareBison.yy"
    {driver.set_olr_inst();;}
    break;

  case 392:
#line 1041 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 395:
#line 1048 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 396:
#line 1049 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 397:
#line 1050 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val));;}
    break;

  case 398:
#line 1053 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 399:
#line 1054 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 400:
#line 1055 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 401:
#line 1058 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 402:
#line 1059 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 403:
#line 1062 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 404:
#line 1063 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 405:
#line 1064 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 406:
#line 1065 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 407:
#line 1066 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 408:
#line 1067 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 409:
#line 1069 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 410:
#line 1070 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 411:
#line 1071 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 412:
#line 1072 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 413:
#line 1073 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 414:
#line 1074 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 415:
#line 1077 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 421:
#line 1089 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 422:
#line 1090 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 423:
#line 1091 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 424:
#line 1092 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 425:
#line 1095 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 426:
#line 1096 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 428:
#line 1100 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 429:
#line 1101 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 430:
#line 1102 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 431:
#line 1103 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 432:
#line 1106 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 433:
#line 1106 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 435:
#line 1110 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 436:
#line 1112 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 437:
#line 1114 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 438:
#line 1116 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 460:
#line 1152 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 461:
#line 1154 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 468:
#line 1168 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 469:
#line 1170 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 470:
#line 1173 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 471:
#line 1174 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 472:
#line 1175 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 473:
#line 1176 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 474:
#line 1177 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 475:
#line 1178 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 476:
#line 1179 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 477:
#line 1180 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 478:
#line 1181 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 479:
#line 1182 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 480:
#line 1183 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 481:
#line 1184 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 482:
#line 1185 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 483:
#line 1186 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 484:
#line 1187 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 485:
#line 1188 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 486:
#line 1189 "DynareBison.yy"
    {driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 487:
#line 1190 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 488:
#line 1191 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 489:
#line 1192 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 490:
#line 1193 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 491:
#line 1195 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 492:
#line 1196 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 493:
#line 1197 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 494:
#line 1199 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1200 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 496:
#line 1201 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 497:
#line 1202 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 498:
#line 1203 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 499:
#line 1204 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 500:
#line 1205 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 501:
#line 1206 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 502:
#line 1207 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 503:
#line 1208 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 504:
#line 1209 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 506:
#line 1211 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 507:
#line 1212 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 508:
#line 1213 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 509:
#line 1214 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 510:
#line 1215 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 511:
#line 1216 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 512:
#line 1217 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 513:
#line 1218 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 514:
#line 1219 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 515:
#line 1220 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 516:
#line 1221 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 517:
#line 1222 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 518:
#line 1223 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 519:
#line 1224 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 520:
#line 1225 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 521:
#line 1226 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 522:
#line 1227 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 523:
#line 1228 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 524:
#line 1229 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 525:
#line 1230 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 526:
#line 1231 "DynareBison.yy"
    {driver.option_num("olr_beta", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 527:
#line 1234 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 528:
#line 1236 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 529:
#line 1238 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 530:
#line 1239 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 531:
#line 1240 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 532:
#line 1241 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 533:
#line 1242 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 534:
#line 1243 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 535:
#line 1244 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 536:
#line 1245 "DynareBison.yy"
    {driver.option_num("load_mh_file", "-1");;}
    break;

  case 537:
#line 1246 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 538:
#line 1248 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1249 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 540:
#line 1250 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1251 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 542:
#line 1252 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 543:
#line 1253 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 544:
#line 1254 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 545:
#line 1255 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 546:
#line 1258 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 548:
#line 1267 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 549:
#line 1270 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 550:
#line 1272 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 551:
#line 1280 "DynareBison.yy"
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
  const short int parser::yypact_ninf_ = -869;
  const short int
  parser::yypact_[] =
  {
       874,    36,    73,   115,   -77,   262,   202,    76,    -6,     4,
     -59,   118,   -12,    14,    77,   113,   586,   417,   589,   245,
     232,   501,   391,   121,   483,   406,   123,   483,   485,    49,
    -869,   416,   422,   483,   439,   594,   595,   617,   208,   249,
     483,   551,   566,   580,   483,    78,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,   521,    79,   577,   916,  -869,   643,   341,  -869,
     621,   653,   593,    48,   231,   668,   499,   695,   699,   747,
    -869,   761,   267,   210,   294,   310,   703,   699,   746,   743,
     627,  -869,   300,   359,   140,   979,   714,  -869,  1095,   275,
     281,   715,  -869,  1095,   318,   387,   671,   400,   748,   638,
    1014,   535,   535,   402,   140,   640,  -869,    50,  -869,   621,
    -869,  1166,   413,  -869,   271,   452,   453,   677,   463,   680,
     464,   691,   466,   474,  -869,  -869,  -869,   785,  -869,   788,
     791,   802,   811,   813,   819,   821,   829,   832,   834,   841,
     842,  -869,   739,   724,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,   844,
     854,   855,  -869,   752,   732,  -869,  -869,  -869,   733,   820,
     126,   548,  -869,   866,   270,  -869,  -869,   741,  -869,   745,
    -869,  -869,   838,   -25,  -869,   839,   360,   887,   331,  -869,
     843,  -869,  -869,   889,  -869,  -869,   900,   902,   906,   907,
     909,  -869,  -869,   910,   911,   913,   918,   919,   920,  -869,
    -869,   922,   925,  -869,  -869,  -869,  -869,   927,   928,  -869,
    -869,   283,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,   929,   876,  -869,   883,  -869,   888,   263,  -869,
     830,   890,   840,   898,   320,  -869,   899,   846,   904,   322,
    -869,   828,   333,  -869,   353,   957,   831,   862,  -869,   865,
    -869,   394,   837,   861,   959,  -869,  -869,   405,  -869,  -869,
    -869,  -869,   914,   921,   107,  -869,  -869,  -869,   864,   979,
     979,   868,   870,   871,   875,   877,   878,   892,   893,   897,
     912,   979,   730,   915,   357,  -869,   967,   994,  1004,  1016,
    1017,  1027,  -869,  -869,  -869,  1034,  1035,  1042,  1047,  -869,
    1048,  -869,  1049,  1062,  -869,  -869,   424,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,   896,   363,   457,  -869,  -869,  -869,   971,
    1020,  -869,   936,  -869,  -869,  -869,   945,  1014,  1014,   951,
     953,   954,   955,   972,   981,   984,   985,   986,   988,  1014,
     670,  -869,   467,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,   495,  -869,   108,
      26,   500,  -869,  -869,  -869,   532,  -869,  -869,   538,  -869,
    -869,  1087,  -869,   546,  -869,  -869,  -869,  -869,  -869,  1024,
    1071,  -869,  -869,  1028,  1077,  -869,  -869,  1041,  1088,  -869,
    -869,  1128,    28,  1129,  1121,    28,  1123,  1096,  1126,    18,
    1143,  1145,  1117,  1118,    79,  1149,  1150,  1137,  1151,   916,
    1155,  1056,  1057,  1127,   506,  1184,  -869,  -869,  1176,   621,
    1064,  -869,  -869,  1068,   -14,  1161,  1069,    -9,  1163,   979,
    -869,  -869,  -869,  1065,  1206,  1207,  1208,  1209,  1212,  1193,
     588,  1214,  1213,  1216,  1220,  1221,  1188,  1104,  1235,   761,
      -1,  1198,  1240,  1142,  -869,  -869,  -869,    67,  1146,   165,
    1152,  -869,  -869,  1158,   165,  1164,  -869,  -869,    89,  -869,
    -869,  -869,  1227,  1202,  -869,  1230,   162,  -869,    83,  -869,
     371,  -869,  1228,  1301,   221,   359,   317,  1196,    20,  -869,
    -869,   979,  1242,   444,   979,   979,   979,   979,   979,   979,
     979,   979,   979,   979,   630,   979,   979,   979,   979,   979,
    -869,   979,  -869,  -869,  1285,  1420,  1308,  1412,  1413,  1415,
     165,  1416,  1417,   682,  1418,  1419,  1421,  1095,    65,  1393,
     774,  -869,   940,   206,  -869,  1348,  -869,    89,  1332,   601,
    1014,  1014,  1014,  1014,  1014,  1014,  1014,  1014,  1014,  1014,
     650,  1014,  1014,  1014,  1014,  1014,  1316,   535,   216,   222,
    -869,  -869,  -869,   979,   355,    80,    50,  1318,   621,  1319,
    1166,   266,  1438,   271,   328,  -869,  1355,  -869,  1356,  -869,
    1357,  -869,  -869,  1442,  1443,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  1435,    29,  -869,  -869,  -869,  -869,  1323,
    -869,  -869,  1328,  -869,  -869,  -869,  -869,  1329,  -869,  1439,
    1330,  1331,  1333,   979,  -869,  -869,  -869,  -869,  -869,   478,
    1334,  -869,  -869,   496,  1335,   923,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  1325,  -869,  -869,  -869,   512,  -869,  1414,  1422,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,   685,  1338,  1362,
    1363,  1423,  1364,   165,  1424,  1345,   165,  -869,  1455,  1456,
    1346,  -869,   699,  1476,  -869,  -869,  -869,  1014,  -869,  -869,
    -869,  1477,   570,  -869,  -869,  -869,  1351,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,   -16,   209,  -869,
    1432,   979,  1433,   338,   964,   572,   698,   710,   718,   931,
     970,  1046,  1082,  1093,  1099,  1107,  -869,   444,   444,  1242,
    1242,  1371,  1147,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,   516,
     979,  -869,  1436,   996,  -869,   517,  -869,  1353,  1153,  1159,
    1165,  1173,  1179,  1185,  1191,  1199,  1205,  1211,  -869,   601,
     601,  1332,  1332,  1371,  -869,  -869,  -869,   536,  -869,   537,
    1217,    26,  1358,  -869,  -869,    75,   979,  -869,  -869,  -869,
    -869,  -869,  -869,   561,  -869,  -869,  -869,   565,  -869,  -869,
    -869,  -869,  -869,  1359,  -869,  -869,  -869,  1437,  -869,  -869,
    1360,  1486,  -869,  -869,  1251,  -869,   335,  -869,   465,  -869,
    1440,  -869,   573,  -869,  -869,  -869,  -869,  -869,  -869,   165,
      67,  1384,   165,  1387,  1388,  -869,  1366,  -869,  -869,  1494,
     386,  1014,  1257,  1487,   371,  -869,   865,   865,   865,   317,
    -869,   165,  -869,  1495,  1263,  1497,  1480,   979,   979,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    1372,  -869,  1269,   979,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
      80,  -869,  -869,  -869,   979,  1225,  -869,  -869,  1482,  -869,
    1330,   979,  -869,  -869,   583,  -869,   585,  1367,  1325,  -869,
    -869,  1398,  1399,  1400,   165,  1378,   165,   165,  -869,   979,
    -869,  1275,  -869,  -869,  -869,  1379,   189,   540,   567,   285,
    1380,   979,  -869,   979,  1376,    12,  1281,   964,  -869,  -869,
    1287,  1231,  -869,  -869,  1508,  1293,  -869,  -869,  1406,  -869,
     165,   165,   165,  1407,  -869,  1385,  1386,  1299,  -869,   865,
    -869,  -869,  -869,   165,  -869,  1305,  1311,  1496,  1389,  1498,
    1425,  -869,  -869,  -869,   979,  -869,    81,  1411,  -869,  1426,
     165,  -869,  -869,  -869,   613,  1390,  -869,  -869,  -869,  1499,
    1392,    25,  1317,  1471,  -869,   165,   274,  1394,  -869,  -869,
    -869,  1506,  -869,   687,   689,   979,    91,  -869,  -869,  -869,
    1391,  1427,  1428,  -869,  -869,  1237,  -869,  -869,   979,  -869,
    -869,  -869,   165,   165,  -869,  1243,  1429,  -869,  -869,   165,
    -869
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
     432,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     2,     4,    29,    30,
      46,    47,    48,    45,     5,     6,     7,    12,     9,    10,
      11,     8,    13,    14,    15,    16,    17,    18,    19,    23,
      25,    24,    20,    21,    22,    26,    27,    28,    31,    32,
      33,    38,    39,    34,    35,    36,    37,    40,    41,    42,
      43,    44,     0,     0,     0,     0,   401,     0,     0,   208,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   252,
     299,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   128,     0,     0,     0,     0,     0,   383,     0,     0,
       0,     0,   379,     0,     0,     0,    76,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   213,     0,   203,     0,
     219,     0,     0,   435,     0,     0,     0,    58,     0,    64,
       0,    70,     0,     0,     1,     3,   460,     0,   543,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   450,   459,     0,   451,   456,   454,   457,   455,   452,
     453,   443,   444,   445,   446,   447,   448,   449,   468,     0,
       0,     0,   462,   467,     0,   464,   463,   465,     0,     0,
     398,     0,   394,     0,     0,   211,   212,     0,    82,     0,
      49,   411,     0,     0,   405,     0,     0,     0,     0,   116,
       0,   517,   534,     0,   522,   500,     0,     0,     0,     0,
       0,   514,   515,     0,     0,     0,     0,     0,     0,   536,
     510,     0,     0,   521,   535,   516,   499,     0,     0,   520,
     518,     0,   304,   340,   329,   305,   306,   307,   308,   309,
     310,   311,   312,   313,   314,   315,   316,   317,   318,   319,
     320,   321,   322,   323,   324,   325,   326,   327,   328,   330,
     331,   332,   333,   334,   335,   336,   337,   338,   339,   341,
     342,   343,   248,     0,   301,     0,   265,     0,     0,   262,
       0,     0,     0,     0,     0,   284,     0,     0,     0,     0,
     278,     0,     0,   120,     0,     0,     0,     0,    84,     0,
     473,     0,     0,     0,     0,   530,   529,     0,   417,   418,
     419,   420,     0,     0,     0,   171,    89,    90,    88,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   370,     0,     0,     0,     0,
       0,     0,   478,   479,   480,     0,     0,     0,     0,   523,
       0,   487,     0,     0,   388,   389,     0,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   235,   237,   239,
     240,   241,   242,   243,   244,   245,   236,   238,   387,   246,
     247,   385,   391,     0,     0,     0,   381,   378,    79,    74,
       0,    55,     0,    80,   146,   147,   166,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     433,   145,     0,   347,   352,   348,   349,   350,   351,   353,
     354,   355,   356,   357,   358,   359,   360,     0,    51,     0,
       0,     0,   216,   217,   218,     0,   206,   207,     0,   224,
     221,     0,   441,     0,   440,   442,   437,   372,    61,    56,
       0,    52,    67,    62,     0,    53,    73,    68,     0,    54,
     367,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   392,   393,     0,     0,
       0,    83,    50,     0,     0,     0,     0,     0,     0,     0,
     114,   115,   253,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   250,     0,   264,   260,   261,   293,     0,   293,
       0,   282,   283,     0,   293,     0,   276,   277,     0,   118,
     119,   111,     0,     0,    85,     0,     0,   140,     0,   141,
       0,   136,     0,     0,     0,     0,     0,     0,     0,   169,
     170,     0,    96,    97,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      86,     0,   368,   369,     0,     0,     0,     0,     0,     0,
     293,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   373,     0,     0,    77,    75,    81,     0,   153,   154,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     168,   201,   202,     0,     0,   193,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    59,    57,    65,    63,    71,
      69,   539,   259,     0,     0,   540,   541,   542,   538,   544,
     491,   494,   493,     0,     0,   492,   495,   496,   531,     0,
     532,   458,     0,   545,   501,   519,   466,     0,   402,     0,
     398,     0,     0,     0,   471,   210,   209,   414,   409,     0,
       0,   408,   403,     0,     0,     0,   533,   481,   524,   525,
     497,   498,   503,   506,   507,   504,   512,   513,   502,   509,
     508,     0,   511,   303,   300,     0,   249,     0,     0,   288,
     295,   289,   294,   291,   296,   290,   292,     0,     0,     0,
     270,     0,     0,   293,     0,     0,   293,   256,     0,     0,
       0,   113,     0,     0,   129,   138,   139,     0,   143,   125,
     124,     0,     0,   123,   126,   127,     0,   132,   130,   527,
     528,   416,   427,   429,   430,   431,   428,     0,   421,   425,
       0,     0,     0,     0,   109,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    87,    92,    91,    93,
      94,    95,     0,   477,   485,   470,   476,   482,   483,   526,
     474,   484,   490,   489,   475,   472,   488,   390,   384,     0,
       0,   376,     0,     0,   380,     0,    78,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   144,   149,
     148,   150,   151,   152,   434,   346,   344,     0,   361,     0,
       0,     0,     0,   198,   199,     0,     0,   215,   214,   205,
     204,   223,   220,     0,   537,   439,   436,     0,    60,    66,
      72,   258,   257,   547,   549,   551,   550,     0,   461,   469,
       0,     0,   400,   399,     0,   410,     0,   404,     0,   117,
       0,   365,     0,   302,   251,   266,   298,   297,   263,   293,
     293,     0,   293,     0,     0,   281,     0,   255,   254,     0,
       0,     0,     0,     0,     0,   134,     0,     0,     0,     0,
     415,   293,   426,     0,     0,     0,     0,     0,     0,   108,
      98,    99,   100,   101,   102,   103,   104,   105,   106,   107,
       0,   386,     0,     0,   374,   382,   167,   155,   156,   157,
     158,   159,   160,   161,   162,   163,   164,   345,   362,   200,
     192,   191,   195,   196,     0,     0,   222,   438,     0,   546,
     398,     0,   395,   412,     0,   406,     0,     0,     0,   505,
     267,     0,     0,     0,   293,     0,   293,   293,   279,     0,
     112,     0,   142,   486,   122,     0,     0,     0,     0,   422,
       0,     0,   174,     0,   182,     0,     0,   110,   371,   377,
       0,     0,   197,   548,     0,     0,   413,   407,     0,   366,
     293,   293,   293,     0,   287,     0,     0,     0,   165,     0,
     137,   133,   131,   293,   423,     0,     0,     0,   177,     0,
       0,   173,   375,   194,     0,   396,   293,   272,   268,   271,
     293,   285,   280,   121,     0,     0,   176,   175,   181,     0,
     179,     0,     0,     0,   364,   293,     0,     0,   135,   424,
     178,     0,   188,     0,     0,     0,     0,   187,   186,   397,
       0,   273,     0,   286,   180,     0,   185,   172,     0,   184,
     183,   363,   293,   293,   190,     0,   274,   269,   189,   293,
     275
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -869,  -869,  1514,  -869,  -869,  -869,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -341,  -869,
    -869,  -869,  -869,  -111,  -226,  -869,  -869,  1238,  -869,   628,
    -869,  -869,  -869,  -869,  -869,  -869,  -802,  -548,  -135,  -545,
    -869,  -869,  -869,  1430,  -263,  -869,  -869,  -869,  -869,   690,
    -869,  -869,   901,  -869,  -869,  1051,  -869,  -869,   905,  -869,
    -869,  -132,   -23,  -604,  -488,  -869,  -869,  1259,  -869,  -869,
    -808,  -869,  -869,  1249,  -869,  -869,  1258,  -868,  -517,  -869,
    -869,  1026,  -869,  1431,   924,  -869,   578,  -869,  -869,  -869,
    -869,  1215,  -869,  -869,  -869,  -869,  -869,  -869,   950,  1445,
    -869,  -869,  -869,  1369,  -675,  -869,  -869,  -869,  -869,  -869,
     997,  -869,   644,  -747,  -869,  -869,  -869,  -869,  -869,   917,
    -869,   -61,  1080,  -869,  -869,  1076,  -869,  -869,   -88,  -869,
    1464,  -869,  -869,  -869,  -869,  -869,  -869,  -869,   -98,  -869,
    -869,  -136,  -541,  -869,  -869,  -869,  -869,   -97,   -84,   -75,
     -74,   -73,  -869,  -869,   -93,   -70,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,   -66,  -869,  -869,  -869,  -869,  -869,
     -53,   -52,   -65,   -51,   -49,   -47,  -869,  -869,  -869,  -869,
    -869,   -94,   -91,   -87,   -85,   -46,  -869,  -869,  -869,  -869,
    -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,  -869,   894,
    -869,  1054
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    45,    46,    47,    48,    49,    50,    51,    52,    53,
     158,   160,   162,   137,    54,    55,    56,    57,   362,   805,
      58,   326,    59,   228,   229,    60,   322,   323,   782,   783,
      61,   329,   938,   937,  1015,   786,   576,   577,   578,   579,
     441,    62,    63,   344,   345,  1025,  1096,    64,   664,   665,
      65,   465,   466,    66,   214,   215,    67,   461,   462,    68,
     468,   384,   112,   770,   685,    69,   308,   309,   310,   758,
    1000,    70,   319,   320,    71,   314,   315,   759,  1001,    72,
     261,   262,    73,   442,   443,    74,   911,   912,    75,    76,
     364,   365,    77,    78,   414,    79,    80,    81,   385,   386,
      82,    83,   211,   212,   515,    84,    85,    86,    87,   337,
     338,   797,   798,   799,    88,   140,   656,    89,   473,   474,
     181,   182,   183,    90,   203,   204,    91,   387,   388,   389,
     390,   391,   392,   393,   394,   395,   396,   397,   398,   399,
     400,   401,   402,   785,   403,   404,   405,   184,   185,   186,
     187,   188,   270,   271,   406,   446,   274,   275,   276,   277,
     278,   279,   280,   281,   447,   283,   284,   285,   286,   287,
     448,   449,   450,   451,   452,   453,   407,   294,   295,   408,
     339,   409,   410,   189,   190,   456,   299,   300,   301,   475,
     191,   192,   193,   194,   195,   196,   197,   207,   700,   894,
     694,   695
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       129,   130,   531,   134,   135,   440,   324,   688,   592,   593,
     143,   463,   216,   264,   265,   152,   155,   156,   272,   469,
     604,   163,   472,   263,   296,   205,   297,   266,   775,   340,
     206,   776,   341,   847,   202,   901,   267,   268,   269,   784,
     760,   273,   762,   444,   444,   282,   290,   765,   445,   445,
     464,   942,  1002,   801,   454,   454,   455,   455,   288,   289,
     291,   467,   292,  1058,   293,   298,   109,   682,   661,   692,
     682,   109,   630,   749,   104,   174,   767,   662,   164,   109,
     893,   590,   138,    98,   106,     1,     2,    92,   167,   168,
     169,   170,   171,   172,   173,     3,     4,     5,   531,   219,
     139,   108,     6,   829,   174,  1092,     7,     8,     9,   750,
      10,   751,    11,    12,    13,    14,   777,   982,   752,   753,
     175,   939,   873,   750,    94,    15,   983,   103,    16,   342,
     342,   874,   752,   682,  1016,  1017,  1018,   524,   589,   660,
     767,    17,   767,   525,   940,   109,   718,   754,   113,  1059,
     377,   722,    18,    19,    20,   105,   755,   802,    21,   744,
     176,   754,   342,  1093,  1094,   107,   683,   684,    22,    23,
      24,  1106,  1060,    25,   114,    26,    27,    28,    29,    30,
     803,   177,   178,   693,    31,    32,  1095,   663,   725,    33,
      34,    35,    36,   774,   756,   895,   590,    93,   109,    37,
      38,   109,    39,   109,   424,   757,    40,   750,   220,    41,
      42,    43,    44,   425,   179,   180,   752,   875,  1102,   757,
    1050,   651,   652,   653,   654,   838,   655,   768,   769,  1093,
    1094,   424,   305,  1068,    95,   101,   984,   115,   343,   343,
     425,   876,   426,   778,   102,   754,   923,  1074,  1083,   926,
     804,  1107,  1108,   806,   807,   808,   809,   810,   811,   812,
     813,   814,   815,   513,   817,   818,   819,   820,   821,   426,
     822,   343,   942,   116,   366,    96,    97,   789,   110,   111,
     749,   127,   128,   132,   133,   305,   109,   514,   109,   792,
     306,   843,   638,   639,   555,   367,   109,   368,   369,   790,
     427,   428,   109,   757,   650,  1117,   429,   430,   431,   432,
     433,   434,   435,   436,   437,  1034,   311,   235,   751,   370,
     371,   438,   870,   439,   236,   575,   753,   427,   428,   109,
     307,   330,   316,   429,   430,   431,   432,   433,   434,   435,
     436,   437,   311,   306,   316,   209,   109,   302,   438,   793,
     439,   561,   575,   566,   372,   302,   373,   256,   374,   335,
     330,   302,   530,   755,   569,   792,   844,   376,   150,   151,
     941,   377,   904,   794,   312,   336,   866,   795,   796,   378,
     379,   380,   868,   307,   571,   381,   382,   383,   612,   213,
     317,   221,   124,   784,   631,   367,   471,   792,   302,   222,
     312,   756,   317,  1003,   303,  1005,   123,   519,   109,   153,
     154,   227,   303,   321,   313,   109,   779,  1010,   303,   331,
     549,   210,    99,   100,  1020,   793,   882,   304,   780,   332,
     318,   216,   520,   227,   781,   411,   334,   363,   946,   205,
     313,   412,   318,   632,   206,   550,  1053,   335,   202,   794,
     119,   264,   265,   795,   796,   303,   272,   793,   947,   120,
     944,   263,   296,   336,   297,   266,   227,   302,   775,   775,
     775,   776,   776,   776,   267,   268,   269,  1097,   416,   273,
     419,   794,   302,   282,   290,   795,   796,  1043,   886,  1045,
    1046,   340,  1109,   302,   341,   993,   288,   289,   291,   962,
     292,   719,   293,   298,   723,   848,   849,   850,   851,   852,
     853,   854,   855,   856,   857,   871,   859,   860,   861,   862,
     863,   872,   527,  1067,   303,  1069,   775,   745,   528,   776,
     463,   580,   302,   302,   125,   985,  1075,   420,   881,   303,
     231,   472,   585,   479,   483,   109,   487,   417,   711,  1084,
     303,   126,   209,  1087,   302,   200,   581,   712,   302,   444,
     421,   627,   458,   109,   445,   136,   131,   586,  1101,   464,
     454,  1051,   455,   470,   233,   234,   302,   141,   201,   516,
     467,   235,   424,   142,   607,   608,   628,   609,   236,   303,
     303,   425,   302,  1098,   627,  1116,   302,   302,  1052,   144,
     480,   484,  1120,   488,   657,   839,  1026,  1027,  1110,   424,
     845,   303,   476,   477,   253,   303,   302,   302,   425,   633,
     426,   256,  1030,   481,   485,   995,   489,   145,   210,   658,
     733,   157,   657,   303,   490,   867,   869,   666,   905,   734,
     258,   302,   932,  1031,  1088,   302,   159,   426,   883,   303,
    1035,   887,   259,   303,   303,   424,   907,   659,   260,   224,
     161,   930,   667,   302,   425,   302,   208,   225,  1047,   668,
     179,   180,   913,   303,   303,   670,   961,   965,   427,   428,
    1055,   166,  1056,   673,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   426,   669,   217,   977,   978,   303,   438,
     671,   439,   303,   575,   531,   427,   428,   934,   674,   948,
     998,   429,   430,   431,   432,   433,   434,   435,   436,   437,
     303,   986,   303,  1082,   832,   987,   438,   916,   439,   891,
     575,   892,   935,   833,   949,   999,   917,   198,   927,   213,
     928,   653,   654,  1036,   655,  1037,   117,   118,   223,   121,
     122,   427,   428,   218,  1105,   146,   147,   429,   430,   431,
     432,   433,   434,   435,   436,   437,   231,  1115,   605,   606,
     607,   608,   438,   609,   439,   226,   575,   148,   149,   227,
     230,   200,   232,   321,   325,   327,   174,   328,   651,   652,
     653,   654,   816,   655,   363,   413,  1011,   418,   423,   422,
     233,   234,   175,   478,   201,   460,   482,   235,   651,   652,
     653,   654,   858,   655,   236,   237,   238,   486,   491,   239,
     240,   492,   241,   242,   493,   243,   244,   245,   246,   247,
     248,   249,   250,   251,   252,   494,   605,   606,   607,   608,
     253,   609,   176,   254,   495,   255,   496,   256,   605,   606,
     607,   608,   497,   609,   498,   257,   605,   606,   607,   608,
     950,   609,   499,   177,   178,   500,   258,   501,   605,   606,
     607,   608,   951,   609,   502,   503,   504,   506,   259,   213,
     952,     1,     2,   994,   260,   996,   505,   507,   508,   509,
     610,     3,     4,     5,   510,   511,   179,   180,     6,   518,
     512,   521,     7,     8,     9,   522,    10,   424,    11,    12,
      13,    14,   605,   606,   607,   608,   425,   609,   523,   526,
     529,    15,   533,   532,    16,   167,   168,   169,   170,   171,
     172,   173,   199,   534,   841,   535,   200,    17,   346,   536,
     537,   174,   538,   539,   540,   426,   541,   347,    18,    19,
      20,   542,   543,   544,    21,   545,   552,   175,   546,   201,
     547,   548,   551,   553,    22,    23,    24,   557,   554,    25,
     558,    26,    27,    28,    29,    30,   348,   559,   560,   563,
      31,    32,   346,   564,   565,    33,    34,    35,    36,   568,
     572,   347,   584,   573,   587,    37,    38,   176,    39,   582,
     614,   588,    40,   427,   428,    41,    42,    43,    44,   429,
     430,   431,   432,   433,   434,   435,   436,   437,   177,   178,
     348,   346,   574,   583,   438,   591,   439,   615,   575,   594,
     347,   595,   596,   629,   349,   350,   597,   616,   598,   599,
     351,   352,   353,   354,   355,   356,   357,   358,   359,   617,
     618,   179,   180,   600,   601,   360,   424,   361,   602,   348,
     619,   605,   606,   607,   608,   425,   609,   620,   621,   605,
     606,   607,   608,   603,   609,   622,   611,   842,   349,   350,
     623,   624,   625,   909,   351,   352,   353,   354,   355,   356,
     357,   358,   359,   953,   426,   626,   636,   634,   366,   360,
     635,   361,   605,   606,   607,   608,   637,   609,   605,   606,
     607,   608,   640,   609,   641,   642,   643,   349,   350,   367,
     672,   368,   369,   351,   352,   353,   354,   355,   356,   357,
     358,   359,   954,   644,   605,   606,   607,   608,   360,   609,
     361,   235,   645,   370,   371,   646,   647,   648,   236,   649,
     675,   676,   427,   428,   677,   330,   964,   678,   429,   430,
     431,   432,   433,   434,   435,   436,   437,   679,   680,   366,
     681,   686,   687,   438,   689,   439,   690,   691,   372,   704,
     373,   256,   374,   335,   605,   606,   607,   608,   375,   609,
     367,   376,   368,   369,   696,   377,   697,   698,   699,   336,
     702,   703,   705,   378,   379,   380,   707,   710,   955,   381,
     382,   383,   235,   213,   370,   371,   708,   713,   709,   236,
     605,   606,   607,   608,   716,   609,   330,   714,   717,   721,
     693,   605,   606,   607,   608,   732,   609,   605,   606,   607,
     608,   720,   609,   724,   956,   605,   606,   607,   608,   372,
     609,   373,   256,   374,   335,   957,   735,   727,   728,   729,
     730,   958,   376,   731,   736,   741,   377,   737,   740,   959,
     336,   738,   739,   747,   378,   379,   380,   742,   746,   748,
     381,   382,   383,   761,   213,   605,   606,   607,   608,   763,
     609,   651,   652,   653,   654,   764,   655,   651,   652,   653,
     654,   766,   655,   651,   652,   653,   654,   771,   655,   960,
     773,   651,   652,   653,   654,   967,   655,   651,   652,   653,
     654,   968,   655,   651,   652,   653,   654,   969,   655,   651,
     652,   653,   654,   800,   655,   970,   823,   651,   652,   653,
     654,   971,   655,   651,   652,   653,   654,   972,   655,   651,
     652,   653,   654,   973,   655,   605,   606,   607,   608,   825,
     609,   974,   772,   605,   606,   607,   608,   975,   609,   605,
     606,   607,   608,   976,   609,   605,   606,   607,   608,   979,
     609,   605,   606,   607,   608,   609,   609,  1032,   787,   605,
     606,   607,   608,  1063,   609,   651,   652,   653,   654,  1114,
     655,   605,   606,   607,   608,  1118,   609,   605,   606,   607,
     608,   992,   609,   651,   652,   653,   654,  1012,   655,   605,
     606,   607,   608,  1022,   609,   605,   606,   607,   608,  1029,
     609,   605,   606,   607,   608,  1048,   609,   605,   606,   607,
     608,  1061,   609,   605,   606,   607,   608,  1062,   609,   605,
     606,   607,   608,  1065,   609,   605,   606,   607,   608,  1073,
     609,   788,   824,   826,   827,  1076,   828,   830,   831,   834,
     835,  1077,   836,   840,   846,   655,   864,  1099,   878,   880,
     884,   888,   889,   890,   891,   892,   893,   897,   898,   899,
     900,   514,   910,   902,   914,   903,   906,   908,   918,   919,
     920,   922,   915,   921,   924,   925,   927,   928,   929,   931,
     933,   936,   943,   945,    -1,   966,   963,   989,   981,   991,
     997,  1004,   990,   988,  1006,  1007,  1008,  1009,  1021,  1013,
    1023,  1024,  1028,  1033,  1038,  1040,  1041,  1042,  1044,  1049,
    1057,  1064,  1054,  1066,  1070,  1071,  1072,  1078,  1085,  1080,
    1090,  1100,  1089,  1079,  1103,  1081,  1091,  1104,  1111,   165,
     570,   980,  1014,  1086,  1112,  1113,  1119,   556,   567,   879,
     715,   877,   562,   457,   459,   743,  1039,   837,   415,   613,
     517,   865,   791,  1019,   701,   706,   333,   726,   896,     0,
     885
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        23,    24,   228,    26,    27,   140,   117,   495,   349,   350,
      33,   147,   100,   111,   111,    38,    39,    40,   111,   151,
     361,    44,   154,   111,   111,    95,   111,   111,   576,   123,
      95,   576,   123,   637,    95,   710,   111,   111,   111,   580,
     557,   111,   559,   141,   142,   111,   111,   564,   141,   142,
     147,   798,   920,    33,   141,   142,   141,   142,   111,   111,
     111,   149,   111,    51,   111,   111,    80,    42,    42,    51,
      42,    80,   413,     6,    80,    25,    51,    51,     0,    80,
      51,   344,    33,   160,    80,     7,     8,    51,     9,    10,
      11,    12,    13,    14,    15,    17,    18,    19,   324,    51,
      51,   160,    24,   620,    25,    80,    28,    29,    30,    42,
      32,    44,    34,    35,    36,    37,    33,    42,    51,    52,
      41,   137,    42,    42,    51,    47,    51,    51,    50,    22,
      22,    51,    51,    42,   936,   937,   938,   162,    31,    31,
      51,    63,    51,   168,   160,    80,   160,    80,   160,   137,
     100,   160,    74,    75,    76,   161,    89,   137,    80,   160,
      81,    80,    22,   138,   139,   161,   138,   139,    90,    91,
      92,    80,   160,    95,   160,    97,    98,    99,   100,   101,
     160,   102,   103,   165,   106,   107,   161,   161,   529,   111,
     112,   113,   114,    31,   127,   166,   459,   161,    80,   121,
     122,    80,   124,    80,    42,   138,   128,    42,   160,   131,
     132,   133,   134,    51,   135,   136,    51,   137,  1086,   138,
      31,   138,   139,   140,   141,   160,   143,   138,   139,   138,
     139,    42,    22,  1041,   161,    33,   161,   160,   131,   131,
      51,   161,    80,   160,    42,    80,   763,  1049,   167,   766,
     591,   160,   161,   594,   595,   596,   597,   598,   599,   600,
     601,   602,   603,   137,   605,   606,   607,   608,   609,    80,
     611,   131,  1019,   160,     3,   160,   161,    56,   160,   161,
       6,   160,   161,   160,   161,    22,    80,   161,    80,    80,
      80,   632,   427,   428,    31,    24,    80,    26,    27,    78,
     138,   139,    80,   138,   439,  1113,   144,   145,   146,   147,
     148,   149,   150,   151,   152,   990,    22,    46,    44,    48,
      49,   159,   663,   161,    53,   163,    52,   138,   139,    80,
     120,    60,    22,   144,   145,   146,   147,   148,   149,   150,
     151,   152,    22,    80,    22,     4,    80,    80,   159,   140,
     161,    31,   163,    31,    83,    80,    85,    86,    87,    88,
      60,    80,    31,    89,    31,    80,   160,    96,   160,   161,
     161,   100,   713,   164,    80,   104,   160,   168,   169,   108,
     109,   110,   160,   120,    31,   114,   115,   116,    31,   118,
      80,   160,   160,   934,    31,    24,   125,    80,    80,   168,
      80,   127,    80,   920,   137,   922,   161,   137,    80,   160,
     161,    80,   137,    80,   120,    80,    45,    31,   137,   119,
     137,    80,   160,   161,   941,   140,   160,   160,    57,   129,
     120,   519,   162,    80,    63,   160,    77,    80,   100,   509,
     120,   160,   120,    80,   509,   162,   161,    88,   509,   164,
      33,   549,   549,   168,   169,   137,   549,   140,   120,    42,
     801,   549,   549,   104,   549,   549,    80,    80,  1016,  1017,
    1018,  1016,  1017,  1018,   549,   549,   549,  1081,   160,   549,
      80,   164,    80,   549,   549,   168,   169,  1004,   160,  1006,
    1007,   585,  1096,    80,   585,   160,   549,   549,   549,   840,
     549,   524,   549,   549,   527,   640,   641,   642,   643,   644,
     645,   646,   647,   648,   649,   160,   651,   652,   653,   654,
     655,   166,   162,  1040,   137,  1042,  1074,   550,   168,  1074,
     666,   137,    80,    80,    33,   876,  1053,   137,   670,   137,
       5,   673,   137,    80,    80,    80,    80,   160,    42,  1066,
     137,   160,     4,  1070,    80,    20,   162,    51,    80,   657,
     160,   137,   160,    80,   657,    80,   160,   162,  1085,   666,
     657,    31,   657,   160,    39,    40,    80,   161,    43,    31,
     668,    46,    42,   161,   140,   141,   162,   143,    53,   137,
     137,    51,    80,  1081,   137,  1112,    80,    80,    31,   160,
     137,   137,  1119,   137,   137,   628,   947,   948,  1096,    42,
     633,   137,   160,   160,    79,   137,    80,    80,    51,   162,
      80,    86,   963,   160,   160,   160,   160,    33,    80,   162,
      42,    80,   137,   137,   160,   658,   659,   137,   160,    51,
     105,    80,   777,   984,    31,    80,    80,    80,   671,   137,
     991,   674,   117,   137,   137,    42,   160,   162,   123,   160,
      80,   772,   162,    80,    51,    80,    23,   168,  1009,   137,
     135,   136,   160,   137,   137,   137,   160,   160,   138,   139,
    1021,   160,  1023,   137,   144,   145,   146,   147,   148,   149,
     150,   151,   152,    80,   162,    42,   160,   160,   137,   159,
     162,   161,   137,   163,   930,   138,   139,   137,   162,   137,
     137,   144,   145,   146,   147,   148,   149,   150,   151,   152,
     137,   160,   137,  1064,    42,   160,   159,    42,   161,    42,
     163,    42,   162,    51,   162,   162,    51,   160,    51,   118,
      51,   140,   141,   160,   143,   160,   160,   161,    80,   160,
     161,   138,   139,   160,  1095,   160,   161,   144,   145,   146,
     147,   148,   149,   150,   151,   152,     5,  1108,   138,   139,
     140,   141,   159,   143,   161,    80,   163,   160,   161,    80,
      33,    20,    21,    80,    38,    42,    25,   160,   138,   139,
     140,   141,   162,   143,    80,    80,   931,   126,   160,    51,
      39,    40,    41,   126,    43,   165,   126,    46,   138,   139,
     140,   141,   162,   143,    53,    54,    55,   126,    33,    58,
      59,    33,    61,    62,    33,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    33,   138,   139,   140,   141,
      79,   143,    81,    82,    33,    84,    33,    86,   138,   139,
     140,   141,    33,   143,    33,    94,   138,   139,   140,   141,
     162,   143,    33,   102,   103,    33,   105,    33,   138,   139,
     140,   141,   162,   143,    33,    33,   137,    33,   117,   118,
     162,     7,     8,   906,   123,   908,   162,    33,    33,   137,
     160,    17,    18,    19,   162,   162,   135,   136,    24,    33,
      80,   160,    28,    29,    30,   160,    32,    42,    34,    35,
      36,    37,   138,   139,   140,   141,    51,   143,    80,    80,
      33,    47,    33,    80,    50,     9,    10,    11,    12,    13,
      14,    15,    16,    33,   160,    33,    20,    63,    42,    33,
      33,    25,    33,    33,    33,    80,    33,    51,    74,    75,
      76,    33,    33,    33,    80,    33,    80,    41,    33,    43,
      33,    33,    33,    80,    90,    91,    92,   137,    80,    95,
      80,    97,    98,    99,   100,   101,    80,   137,    80,    80,
     106,   107,    42,   137,    80,   111,   112,   113,   114,   161,
      33,    51,    33,   162,    80,   121,   122,    81,   124,   162,
      33,    80,   128,   138,   139,   131,   132,   133,   134,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   102,   103,
      80,    42,   160,   162,   159,   161,   161,    33,   163,   161,
      51,   161,   161,   137,   138,   139,   161,    33,   161,   161,
     144,   145,   146,   147,   148,   149,   150,   151,   152,    33,
      33,   135,   136,   161,   161,   159,    42,   161,   161,    80,
      33,   138,   139,   140,   141,    51,   143,    33,    33,   138,
     139,   140,   141,   161,   143,    33,   161,   137,   138,   139,
      33,    33,    33,   160,   144,   145,   146,   147,   148,   149,
     150,   151,   152,   162,    80,    33,   160,   126,     3,   159,
      80,   161,   138,   139,   140,   141,   161,   143,   138,   139,
     140,   141,   161,   143,   161,   161,   161,   138,   139,    24,
      33,    26,    27,   144,   145,   146,   147,   148,   149,   150,
     151,   152,   162,   161,   138,   139,   140,   141,   159,   143,
     161,    46,   161,    48,    49,   161,   161,   161,    53,   161,
     126,    80,   138,   139,   126,    60,   160,    80,   144,   145,
     146,   147,   148,   149,   150,   151,   152,   126,    80,     3,
      42,    42,    51,   159,    51,   161,    80,    51,    83,    42,
      85,    86,    87,    88,   138,   139,   140,   141,    93,   143,
      24,    96,    26,    27,    51,   100,    51,    80,    80,   104,
      51,    51,    51,   108,   109,   110,    51,    80,   162,   114,
     115,   116,    46,   118,    48,    49,   160,    33,   161,    53,
     138,   139,   140,   141,   160,   143,    60,    51,   160,   160,
     165,   138,   139,   140,   141,    42,   143,   138,   139,   140,
     141,    80,   143,    80,   162,   138,   139,   140,   141,    83,
     143,    85,    86,    87,    88,   162,    42,    51,    51,    51,
      51,   162,    96,    51,    51,   161,   100,    51,    80,   162,
     104,    51,    51,    33,   108,   109,   110,    42,    80,   137,
     114,   115,   116,   137,   118,   138,   139,   140,   141,   137,
     143,   138,   139,   140,   141,   137,   143,   138,   139,   140,
     141,   137,   143,   138,   139,   140,   141,    80,   143,   162,
      80,   138,   139,   140,   141,   162,   143,   138,   139,   140,
     141,   162,   143,   138,   139,   140,   141,   162,   143,   138,
     139,   140,   141,   137,   143,   162,    51,   138,   139,   140,
     141,   162,   143,   138,   139,   140,   141,   162,   143,   138,
     139,   140,   141,   162,   143,   138,   139,   140,   141,    51,
     143,   162,   160,   138,   139,   140,   141,   162,   143,   138,
     139,   140,   141,   162,   143,   138,   139,   140,   141,   162,
     143,   138,   139,   140,   141,   143,   143,   162,   160,   138,
     139,   140,   141,   162,   143,   138,   139,   140,   141,   162,
     143,   138,   139,   140,   141,   162,   143,   138,   139,   140,
     141,   160,   143,   138,   139,   140,   141,   160,   143,   138,
     139,   140,   141,   160,   143,   138,   139,   140,   141,   160,
     143,   138,   139,   140,   141,   160,   143,   138,   139,   140,
     141,   160,   143,   138,   139,   140,   141,   160,   143,   138,
     139,   140,   141,   160,   143,   138,   139,   140,   141,   160,
     143,   160,    42,    51,    51,   160,    51,    51,    51,    51,
      51,   160,    51,    80,   126,   143,   160,   160,   160,   160,
      42,   126,   126,   126,    42,    42,    51,   164,   160,   160,
      51,   161,   167,   162,    80,   162,   162,   162,   160,   137,
     137,   137,    80,    80,    80,   160,    51,    51,   162,    33,
      33,   160,    80,    80,   143,   162,    80,    80,   160,    33,
      80,   137,   162,   164,   137,   137,   160,    33,    33,    42,
      33,    51,   160,    51,   167,   137,   137,   137,   160,   160,
     164,    33,   162,   137,   137,   160,   160,    51,   137,    51,
      51,    80,   162,   164,   160,   130,   164,    51,   167,    45,
     322,   871,   934,   137,   137,   137,   137,   308,   319,   668,
     519,   666,   314,   142,   144,   549,   998,   627,   133,   364,
     211,   657,   585,   939,   504,   509,   122,   533,   694,    -1,
     673
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
     310,   311,   312,   314,   315,   316,   324,   346,   349,   351,
     352,   160,   160,    80,   264,   269,   160,   160,   126,    80,
     137,   160,    51,   160,    42,    51,    80,   138,   139,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   159,   161,
     208,   210,   253,   254,   308,   324,   325,   334,   340,   341,
     342,   343,   344,   345,   353,   354,   355,   253,   160,   213,
     165,   227,   228,   311,   317,   221,   222,   298,   230,   231,
     160,   125,   231,   288,   289,   359,   160,   160,   126,    80,
     137,   160,   126,    80,   137,   160,   126,    80,   137,   160,
     160,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,   137,   162,    33,    33,    33,   137,
     162,   162,    80,   137,   161,   274,    31,   273,    33,   137,
     162,   160,   160,    80,   162,   168,    80,   162,   168,    33,
      31,   194,    80,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,   137,
     162,    33,    80,    80,    80,    31,   237,   137,    80,   137,
      80,    31,   246,    80,   137,    80,    31,   243,   161,    31,
     197,    31,    33,   162,   160,   163,   206,   207,   208,   209,
     137,   162,   162,   162,    33,   137,   162,    80,    80,    31,
     214,   161,   188,   188,   161,   161,   161,   161,   161,   161,
     161,   161,   161,   161,   188,   138,   139,   140,   141,   143,
     160,   161,    31,   261,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,   137,   162,   137,
     188,    31,    80,   162,   126,    80,   160,   161,   208,   208,
     161,   161,   161,   161,   161,   161,   161,   161,   161,   161,
     208,   138,   139,   140,   141,   143,   286,   137,   162,   162,
      31,    42,    51,   161,   218,   219,   137,   162,   137,   162,
     137,   162,    33,   137,   162,   126,    80,   126,    80,   126,
      80,    42,    42,   138,   139,   234,    42,    51,   234,    51,
      80,    51,    51,   165,   370,   371,    51,    51,    80,    80,
     368,   292,    51,    51,    42,    51,   295,    51,   160,   161,
      80,    42,    51,    33,    51,   225,   160,   160,   160,   232,
      80,   160,   160,   232,    80,   188,   371,    51,    51,    51,
      51,    51,    42,    42,    51,    42,    51,    51,    51,    51,
      80,   161,    42,   251,   160,   232,    80,    33,   137,     6,
      42,    44,    51,    52,    80,    89,   127,   138,   239,   247,
     248,   137,   248,   137,   137,   248,   137,    51,   138,   139,
     233,    80,   160,    80,    31,   207,   209,    33,   160,    45,
      57,    63,   198,   199,   312,   313,   205,   160,   160,    56,
      78,   280,    80,   140,   164,   168,   169,   281,   282,   283,
     137,    33,   137,   160,   188,   189,   188,   188,   188,   188,
     188,   188,   188,   188,   188,   188,   162,   188,   188,   188,
     188,   188,   188,    51,    42,    51,    51,    51,    51,   248,
      51,    51,    42,    51,    51,    51,    51,   268,   160,   232,
      80,   160,   137,   188,   160,   232,   126,   233,   208,   208,
     208,   208,   208,   208,   208,   208,   208,   208,   162,   208,
     208,   208,   208,   208,   160,   254,   160,   232,   160,   232,
     188,   160,   166,    42,    51,   137,   161,   228,   160,   222,
     160,   231,   160,   232,    42,   289,   160,   232,   126,   126,
     126,    42,    42,    51,   369,   166,   369,   164,   160,   160,
      51,   274,   162,   162,   188,   160,   162,   160,   162,   160,
     167,   256,   257,   160,    80,    80,    42,    51,   160,   137,
     137,    80,   137,   248,    80,   160,   248,    51,    51,   162,
     193,    33,   208,    33,   137,   162,   160,   203,   202,   137,
     160,   161,   283,    80,   188,    80,   100,   120,   137,   162,
     162,   162,   162,   162,   162,   162,   162,   162,   162,   162,
     162,   160,   188,    80,   160,   160,   162,   162,   162,   162,
     162,   162,   162,   162,   162,   162,   162,   160,   160,   162,
     219,   160,    42,    51,   161,   188,   160,   160,   164,    80,
     162,    33,   160,   160,   232,   160,   232,    80,   137,   162,
     240,   248,   247,   248,   137,   248,   137,   137,   160,    33,
      31,   208,   160,    42,   199,   204,   206,   206,   206,   282,
     248,    33,   160,    33,    51,   215,   188,   188,   160,   160,
     188,   188,   162,    51,   274,   188,   160,   160,   167,   256,
     137,   137,   137,   248,   160,   248,   248,   188,   160,   160,
      31,    31,    31,   161,   162,   188,   188,   164,    51,   137,
     160,   160,   160,   162,    33,   160,   137,   248,   240,   248,
     137,   160,   160,   160,   206,   248,   160,   160,    51,   164,
      51,   130,   188,   167,   248,   137,   137,   248,    31,   162,
      51,   164,    80,   138,   139,   161,   216,   233,   234,   160,
      80,   248,   247,   160,    51,   188,    80,   160,   161,   233,
     234,   167,   137,   137,   162,   188,   248,   240,   162,   137,
     248
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
     231,   231,   231,   231,   231,   231,   231,   231,   232,   232,
     232,   232,   232,   232,   233,   233,   233,   234,   234,   234,
     235,   236,   236,   237,   238,   238,   238,   239,   239,   239,
     239,   239,   240,   240,   240,   240,   241,   242,   242,   243,
     243,   243,   244,   245,   245,   246,   246,   246,   247,   247,
     247,   247,   247,   248,   248,   248,   248,   248,   248,   249,
     249,   249,   249,   250,   250,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   251,
     251,   251,   251,   251,   252,   252,   253,   253,   254,   254,
     254,   254,   254,   254,   254,   254,   254,   254,   254,   254,
     254,   255,   255,   256,   256,   257,   257,   258,   259,   260,
     260,   261,   262,   263,   264,   264,   264,   264,   265,   266,
     266,   266,   266,   267,   267,   267,   267,   268,   268,   269,
     269,   270,   271,   272,   272,   273,   273,   273,   274,   274,
     274,   275,   275,   276,   276,   276,   276,   276,   276,   277,
     277,   277,   277,   277,   277,   278,   279,   279,   280,   280,
     280,   281,   281,   281,   281,   282,   282,   283,   283,   283,
     283,   283,   285,   286,   284,   287,   287,   287,   287,   288,
     288,   289,   289,   290,   290,   290,   290,   290,   290,   290,
     291,   291,   291,   291,   291,   291,   291,   291,   292,   292,
     293,   293,   294,   294,   294,   294,   295,   295,   296,   296,
     297,   298,   299,   300,   301,   302,   303,   304,   305,   306,
     307,   308,   309,   310,   311,   312,   313,   314,   315,   316,
     316,   317,   318,   318,   319,   320,   321,   322,   323,   324,
     324,   325,   326,   327,   328,   329,   330,   330,   331,   332,
     333,   334,   335,   336,   337,   338,   339,   340,   341,   342,
     343,   344,   345,   346,   347,   348,   349,   350,   350,   351,
     352,   353,   354,   355,   356,   357,   358,   359,   360,   361,
     362,   363,   364,   365,   366,   367,   368,   369,   369,   370,
     370,   371
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
       5,     3,     6,     2,     5,     3,     6,     1,     1,     1,
       3,     3,     4,     2,     1,     5,     7,     9,     0,     3,
       3,     2,     5,     5,     6,     3,     7,     8,     5,     5,
       6,     3,     7,     8,     5,     6,     3,     1,     1,     1,
       1,     1,     3,     4,     6,     1,     2,     1,     1,     1,
       1,     1,     0,     0,     5,     2,     5,     3,     6,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     1,
       3,     6,     1,     1,     1,     1,     3,     1,     3,     6,
       3,     3,     3,     1,     3,     3,     3,     3,     1,     1,
       1,     3,     3,     3,     3,     3,     3,     1,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       1,     3,     3,     3,     3,     5,     3,     3,     3,     3,
       1,     3,     3,     3,     1,     1,     1,     1,     1,     3,
       1,     1,     1,     1,     3,     3,     3,     3,     3,     1,
       1,     3,     3,     3,     1,     1,     1,     3,     3,     3,
       3,     3,     3,     1,     3,     3,     3,     1,     3,     2,
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
      -1,   312,    -1,   314,    -1,   315,    -1,   316,    -1,   351,
      -1,   352,    -1,   232,    80,    -1,   232,    80,    33,    80,
      -1,   232,   137,    80,    -1,   232,   137,    80,    33,    80,
      -1,    80,    -1,    80,    33,    80,    -1,   139,    51,    -1,
     138,    51,    -1,    51,    -1,   139,    42,    -1,   138,    42,
      -1,    42,    -1,    35,   160,   236,    31,    -1,   236,   237,
      -1,   237,    -1,   238,   137,   239,   160,    -1,   120,    80,
      -1,    80,    -1,    22,    80,   137,    80,    -1,   247,   137,
     240,    -1,   248,   137,   247,   137,   240,    -1,   248,   137,
     248,   137,   248,   137,   247,   137,   240,    -1,   248,    -1,
     248,   137,   248,   137,   248,    -1,   248,   137,   248,    -1,
     248,   137,   248,   137,   248,    -1,   248,   137,   248,   137,
     248,   137,   248,    -1,   248,   137,   248,   137,   248,   137,
     248,   137,   248,    -1,    37,   160,   242,    31,    -1,   242,
     243,    -1,   243,    -1,   120,    80,   137,   248,   160,    -1,
      22,    80,   137,    80,   137,   248,   160,    -1,    80,   137,
     248,   160,    -1,    36,   160,   245,    31,    -1,   245,   246,
      -1,   246,    -1,   120,    80,   137,   248,   137,   248,   160,
      -1,    22,    80,   137,    80,   137,   248,   137,   248,   160,
      -1,    80,   137,   248,   137,   248,   160,    -1,     6,    -1,
      44,    -1,    89,    -1,    52,    -1,   127,    -1,    -1,    51,
      -1,    42,    -1,    80,    -1,   138,    51,    -1,   138,    42,
      -1,    34,   160,    -1,    34,   161,   250,   162,   160,    -1,
      34,   232,   160,    -1,    34,   161,   250,   162,   232,   160,
      -1,   250,   137,   251,    -1,   251,    -1,   317,    -1,   318,
      -1,   319,    -1,   320,    -1,   321,    -1,   322,    -1,   323,
      -1,   324,    -1,   325,    -1,   326,    -1,   327,    -1,   328,
      -1,   329,    -1,   330,    -1,   331,    -1,   332,    -1,   333,
      -1,   334,    -1,   335,    -1,   336,    -1,   337,    -1,   338,
      -1,   339,    -1,   340,    -1,   308,    -1,   341,    -1,   342,
      -1,   343,    -1,   344,    -1,   345,    -1,   347,    -1,   348,
      -1,   353,    -1,   354,    -1,   355,    -1,   298,    -1,   356,
      -1,   357,    -1,   358,    -1,   106,   161,   253,   162,   160,
      -1,   106,   161,   253,   162,   232,   160,    -1,   253,   137,
     254,    -1,   254,    -1,   324,    -1,   325,    -1,   334,    -1,
     340,    -1,   308,    -1,   341,    -1,   342,    -1,   343,    -1,
     344,    -1,   345,    -1,   353,    -1,   354,    -1,   355,    -1,
     107,   161,   253,   162,   160,    -1,   107,   161,   253,   162,
     232,   160,    -1,   167,    80,   167,   137,   167,    80,   167,
      -1,   167,    80,   167,   137,   248,    -1,   256,    -1,   257,
     137,   256,    -1,   134,   232,   160,    -1,    90,   160,   260,
      31,    -1,   260,   261,    -1,   261,    -1,    80,   161,   188,
     162,   160,    -1,   128,   232,   160,    -1,    95,   160,   264,
      31,    -1,   264,    80,   188,   160,    -1,   264,    80,   137,
      80,   188,   160,    -1,    80,   188,   160,    -1,    80,   137,
      80,   188,   160,    -1,    98,   232,   160,    -1,    97,   160,
      -1,    97,   161,   269,   162,   160,    -1,    97,   232,   160,
      -1,    97,   161,   269,   162,   232,   160,    -1,    91,   160,
      -1,    91,   161,   269,   162,   160,    -1,    91,   232,   160,
      -1,    91,   161,   269,   162,   232,   160,    -1,   349,    -1,
     231,    -1,   268,    -1,   269,   137,   268,    -1,    92,   232,
     160,    -1,    18,   160,   272,    31,    -1,   272,   273,    -1,
     273,    -1,    80,   274,    33,   188,   160,    -1,    80,   137,
      80,   274,    33,   188,   160,    -1,     4,    80,   161,    51,
     162,   274,    33,   188,   160,    -1,    -1,   161,    51,   162,
      -1,   161,    42,   162,    -1,    17,   160,    -1,    17,   161,
      23,   162,   160,    -1,    30,   161,    80,   162,   160,    -1,
      30,   161,    80,   162,   232,   160,    -1,    30,    80,   160,
      -1,    30,   161,    80,   168,    80,   162,   160,    -1,    30,
     161,    80,   168,    80,   162,   232,   160,    -1,    30,    80,
     168,    80,   160,    -1,    29,   161,    80,   162,   160,    -1,
      29,   161,    80,   162,   232,   160,    -1,    29,    80,   160,
      -1,    29,   161,    80,   168,    80,   162,   160,    -1,    29,
     161,    80,   168,    80,   162,   232,   160,    -1,    29,    80,
     168,    80,   160,    -1,    75,   161,   279,   162,   281,   160,
      -1,   279,   137,   280,    -1,   280,    -1,   350,    -1,   351,
      -1,   352,    -1,   282,    -1,   281,   137,   282,    -1,   282,
     161,   248,   162,    -1,   281,   137,   282,   161,   248,   162,
      -1,   283,    -1,   282,   283,    -1,    80,    -1,   169,    -1,
     140,    -1,   164,    -1,   168,    -1,    -1,    -1,   101,   285,
     208,   286,   160,    -1,   124,   160,    -1,   124,   161,   288,
     162,   160,    -1,   124,   232,   160,    -1,   124,   161,   288,
     162,   232,   160,    -1,   288,   137,   289,    -1,   289,    -1,
     231,    -1,   359,    -1,   360,    -1,   361,    -1,   362,    -1,
     363,    -1,   364,    -1,   365,    -1,   366,    -1,   290,    -1,
     317,    -1,   353,    -1,   354,    -1,   319,    -1,   321,    -1,
     318,    -1,   320,    -1,   291,   137,   292,    -1,   291,    -1,
       7,    51,   160,    -1,     7,   161,   292,   162,    51,   160,
      -1,   291,    -1,   342,    -1,   325,    -1,   367,    -1,   294,
     137,   295,    -1,   294,    -1,     8,    51,   160,    -1,     8,
     161,   295,   162,    51,   160,    -1,    26,    33,    51,    -1,
     118,    33,    51,    -1,   115,    33,    51,    -1,    60,    -1,
      96,    33,    51,    -1,   110,    33,    51,    -1,    27,    33,
      51,    -1,     3,    33,    51,    -1,    83,    -1,    85,    -1,
      87,    -1,    53,    33,    51,    -1,    48,    33,    51,    -1,
      49,    33,    51,    -1,   100,    33,    51,    -1,    24,    33,
      42,    -1,    63,    33,    42,    -1,   114,    -1,   116,    33,
      51,    -1,   108,    33,    51,    -1,   108,    33,    42,    -1,
      25,    33,    80,    -1,    81,    33,   371,    -1,    81,    33,
      51,    -1,    41,    33,    51,    -1,   102,    33,    51,    -1,
     103,    33,    51,    -1,    58,    33,    51,    -1,    59,    33,
      51,    -1,    86,    -1,    46,    -1,    20,    33,    42,    -1,
      69,    33,    51,    -1,    64,    33,    42,    -1,    66,    33,
      42,    -1,    94,    33,   161,   257,   162,    -1,    65,    33,
      42,    -1,    65,    33,    51,    -1,    73,    33,    80,    -1,
      72,    33,    51,    -1,    71,    -1,   105,    33,    42,    -1,
      67,    33,    51,    -1,    68,    33,    51,    -1,    61,    -1,
      62,    -1,    84,    -1,     5,    -1,   123,    -1,    43,    33,
      51,    -1,   117,    -1,    79,    -1,    40,    -1,   109,    -1,
      54,    33,    51,    -1,    55,    33,    51,    -1,    93,    33,
     248,    -1,    77,    33,    56,    -1,    77,    33,    78,    -1,
     104,    -1,    88,    -1,   135,    33,    80,    -1,   136,    33,
     368,    -1,    39,    33,   371,    -1,    21,    -1,    82,    -1,
      70,    -1,   125,    33,    42,    -1,    14,    33,   234,    -1,
       9,    33,    42,    -1,    11,    33,   234,    -1,    12,    33,
      42,    -1,    13,    33,    51,    -1,    10,    -1,    15,    33,
      51,    -1,    16,    33,    51,    -1,    80,   164,    80,    -1,
      51,    -1,    51,   164,    51,    -1,   165,   369,    -1,   370,
     369,    -1,   370,   166,    -1
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
     827,   829,   831,   833,   835,   837,   839,   841,   843,   846,
     851,   855,   861,   863,   867,   870,   873,   875,   878,   881,
     883,   888,   891,   893,   898,   901,   903,   908,   912,   918,
     928,   930,   936,   940,   946,   954,   964,   969,   972,   974,
     980,   988,   993,   998,  1001,  1003,  1011,  1021,  1028,  1030,
    1032,  1034,  1036,  1038,  1039,  1041,  1043,  1045,  1048,  1051,
    1054,  1060,  1064,  1071,  1075,  1077,  1079,  1081,  1083,  1085,
    1087,  1089,  1091,  1093,  1095,  1097,  1099,  1101,  1103,  1105,
    1107,  1109,  1111,  1113,  1115,  1117,  1119,  1121,  1123,  1125,
    1127,  1129,  1131,  1133,  1135,  1137,  1139,  1141,  1143,  1145,
    1147,  1149,  1151,  1153,  1155,  1161,  1168,  1172,  1174,  1176,
    1178,  1180,  1182,  1184,  1186,  1188,  1190,  1192,  1194,  1196,
    1198,  1200,  1206,  1213,  1221,  1227,  1229,  1233,  1237,  1242,
    1245,  1247,  1253,  1257,  1262,  1267,  1274,  1278,  1284,  1288,
    1291,  1297,  1301,  1308,  1311,  1317,  1321,  1328,  1330,  1332,
    1334,  1338,  1342,  1347,  1350,  1352,  1358,  1366,  1376,  1377,
    1381,  1385,  1388,  1394,  1400,  1407,  1411,  1419,  1428,  1434,
    1440,  1447,  1451,  1459,  1468,  1474,  1481,  1485,  1487,  1489,
    1491,  1493,  1495,  1499,  1504,  1511,  1513,  1516,  1518,  1520,
    1522,  1524,  1526,  1527,  1528,  1534,  1537,  1543,  1547,  1554,
    1558,  1560,  1562,  1564,  1566,  1568,  1570,  1572,  1574,  1576,
    1578,  1580,  1582,  1584,  1586,  1588,  1590,  1592,  1594,  1598,
    1600,  1604,  1611,  1613,  1615,  1617,  1619,  1623,  1625,  1629,
    1636,  1640,  1644,  1648,  1650,  1654,  1658,  1662,  1666,  1668,
    1670,  1672,  1676,  1680,  1684,  1688,  1692,  1696,  1698,  1702,
    1706,  1710,  1714,  1718,  1722,  1726,  1730,  1734,  1738,  1742,
    1744,  1746,  1750,  1754,  1758,  1762,  1768,  1772,  1776,  1780,
    1784,  1786,  1790,  1794,  1798,  1800,  1802,  1804,  1806,  1808,
    1812,  1814,  1816,  1818,  1820,  1824,  1828,  1832,  1836,  1840,
    1842,  1844,  1848,  1852,  1856,  1858,  1860,  1862,  1866,  1870,
    1874,  1878,  1882,  1886,  1888,  1892,  1896,  1900,  1902,  1906,
    1909,  1912
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
     619,   620,   621,   622,   623,   624,   625,   626,   630,   632,
     634,   636,   638,   640,   645,   647,   649,   654,   656,   658,
     663,   668,   670,   675,   679,   684,   689,   699,   704,   710,
     720,   725,   736,   742,   750,   760,   774,   778,   780,   784,
     791,   800,   809,   813,   815,   819,   828,   839,   851,   853,
     855,   857,   859,   864,   865,   866,   867,   868,   870,   877,
     879,   881,   883,   888,   889,   892,   893,   894,   895,   896,
     897,   898,   899,   900,   901,   902,   903,   904,   905,   906,
     907,   908,   909,   910,   911,   912,   913,   914,   915,   916,
     917,   918,   919,   920,   921,   922,   923,   924,   925,   926,
     927,   928,   929,   930,   934,   936,   941,   942,   946,   947,
     948,   949,   950,   951,   952,   953,   954,   955,   956,   957,
     958,   962,   964,   969,   970,   974,   975,   979,   984,   989,
     990,   993,   997,  1000,  1004,  1006,  1008,  1010,  1014,  1017,
    1018,  1019,  1020,  1023,  1024,  1025,  1026,  1029,  1030,  1033,
    1034,  1037,  1040,  1044,  1045,  1048,  1049,  1050,  1053,  1054,
    1055,  1058,  1059,  1062,  1063,  1064,  1065,  1066,  1067,  1069,
    1070,  1071,  1072,  1073,  1074,  1076,  1080,  1081,  1084,  1085,
    1086,  1089,  1090,  1091,  1092,  1095,  1096,  1099,  1100,  1101,
    1102,  1103,  1106,  1106,  1106,  1109,  1111,  1113,  1115,  1120,
    1121,  1124,  1125,  1128,  1129,  1130,  1131,  1132,  1133,  1134,
    1137,  1138,  1139,  1140,  1141,  1142,  1143,  1144,  1147,  1148,
    1151,  1153,  1157,  1158,  1159,  1160,  1163,  1164,  1167,  1169,
    1173,  1174,  1175,  1176,  1177,  1178,  1179,  1180,  1181,  1182,
    1183,  1184,  1185,  1186,  1187,  1188,  1189,  1190,  1191,  1192,
    1193,  1195,  1196,  1197,  1199,  1200,  1201,  1202,  1203,  1204,
    1205,  1206,  1207,  1208,  1209,  1210,  1211,  1212,  1213,  1214,
    1215,  1216,  1217,  1218,  1219,  1220,  1221,  1222,  1223,  1224,
    1225,  1226,  1227,  1228,  1229,  1230,  1231,  1233,  1235,  1238,
    1239,  1240,  1241,  1242,  1243,  1244,  1245,  1246,  1248,  1249,
    1250,  1251,  1252,  1253,  1254,  1255,  1257,  1265,  1266,  1270,
    1271,  1280
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
  const int parser::yylast_ = 1590;
  const int parser::yynnts_ = 202;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 164;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 170;

  const unsigned int parser::yyuser_token_number_max_ = 414;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1282 "DynareBison.yy"


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

