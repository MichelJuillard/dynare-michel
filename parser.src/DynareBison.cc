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
	  case 47:
#line 143 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 48:
#line 144 "DynareBison.yy"
    {driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 49:
#line 147 "DynareBison.yy"
    {driver.rplot();;}
    break;

  case 54:
#line 167 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 55:
#line 169 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 56:
#line 171 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 57:
#line 173 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 58:
#line 175 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 59:
#line 177 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 60:
#line 182 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 61:
#line 184 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 62:
#line 186 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 63:
#line 188 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 64:
#line 190 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 65:
#line 192 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 66:
#line 197 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 67:
#line 199 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 68:
#line 201 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 69:
#line 203 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 70:
#line 205 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 71:
#line 207 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 72:
#line 212 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 73:
#line 214 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 74:
#line 216 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 75:
#line 218 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 76:
#line 220 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 77:
#line 222 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 78:
#line 227 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 79:
#line 231 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 80:
#line 238 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 81:
#line 242 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 82:
#line 249 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 83:
#line 253 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 84:
#line 261 "DynareBison.yy"
    {driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 85:
#line 266 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 86:
#line 268 "DynareBison.yy"
    {(yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 87:
#line 270 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 88:
#line 272 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 89:
#line 274 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 90:
#line 276 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 91:
#line 278 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 92:
#line 280 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 93:
#line 282 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 94:
#line 284 "DynareBison.yy"
    {(yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 95:
#line 286 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 96:
#line 288 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 97:
#line 290 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 98:
#line 292 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 99:
#line 294 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 100:
#line 296 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 101:
#line 298 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 102:
#line 300 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 103:
#line 302 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 104:
#line 304 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 105:
#line 306 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 106:
#line 308 "DynareBison.yy"
    {(yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 107:
#line 313 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 108:
#line 315 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 109:
#line 319 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 110:
#line 321 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 111:
#line 325 "DynareBison.yy"
    {driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 112:
#line 330 "DynareBison.yy"
    {driver.end_endval();;}
    break;

  case 115:
#line 340 "DynareBison.yy"
    {driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 116:
#line 345 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 119:
#line 355 "DynareBison.yy"
    {driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 122:
#line 363 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 123:
#line 364 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 126:
#line 370 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 127:
#line 370 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 128:
#line 371 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 129:
#line 372 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 130:
#line 373 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 131:
#line 374 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 132:
#line 375 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
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

  case 140:
#line 390 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 141:
#line 392 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val));;}
    break;

  case 142:
#line 396 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 144:
#line 399 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 145:
#line 401 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 146:
#line 403 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 147:
#line 405 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 148:
#line 407 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 149:
#line 409 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 150:
#line 411 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 151:
#line 413 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 152:
#line 415 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 153:
#line 417 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 154:
#line 419 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 155:
#line 421 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 156:
#line 423 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 157:
#line 425 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 158:
#line 427 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 159:
#line 429 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 160:
#line 431 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 161:
#line 433 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 162:
#line 435 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 163:
#line 439 "DynareBison.yy"
    {driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 164:
#line 443 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 165:
#line 445 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 166:
#line 449 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 167:
#line 453 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 170:
#line 463 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 171:
#line 465 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 172:
#line 467 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 173:
#line 469 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 174:
#line 471 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 175:
#line 476 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 176:
#line 478 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 177:
#line 480 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 178:
#line 482 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 179:
#line 484 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 180:
#line 486 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 181:
#line 492 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 182:
#line 494 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 183:
#line 496 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 184:
#line 498 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 185:
#line 500 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 186:
#line 502 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 187:
#line 504 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 188:
#line 506 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 189:
#line 511 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 190:
#line 516 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 191:
#line 518 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 192:
#line 523 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 193:
#line 525 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 194:
#line 527 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 195:
#line 529 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 196:
#line 531 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 197:
#line 533 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 198:
#line 535 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 199:
#line 537 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 200:
#line 539 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 201:
#line 544 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 202:
#line 548 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 206:
#line 560 "DynareBison.yy"
    {driver.check();;}
    break;

  case 207:
#line 562 "DynareBison.yy"
    {driver.check();;}
    break;

  case 211:
#line 574 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 212:
#line 576 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 217:
#line 589 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 218:
#line 591 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 219:
#line 593 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 220:
#line 595 "DynareBison.yy"
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
#line 1022 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 384:
#line 1029 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 385:
#line 1030 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 386:
#line 1031 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val));;}
    break;

  case 387:
#line 1034 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 388:
#line 1035 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 389:
#line 1036 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 390:
#line 1039 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 391:
#line 1040 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 392:
#line 1043 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 393:
#line 1044 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 394:
#line 1045 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 395:
#line 1046 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 396:
#line 1047 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 397:
#line 1048 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 398:
#line 1050 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 399:
#line 1051 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 400:
#line 1052 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 401:
#line 1053 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 402:
#line 1054 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 403:
#line 1055 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 404:
#line 1058 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 410:
#line 1070 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 411:
#line 1071 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 412:
#line 1072 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 413:
#line 1073 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 414:
#line 1076 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 415:
#line 1077 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 417:
#line 1081 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 418:
#line 1082 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 419:
#line 1083 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 420:
#line 1084 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 421:
#line 1087 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 422:
#line 1087 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 424:
#line 1091 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 425:
#line 1093 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 426:
#line 1095 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 427:
#line 1097 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 449:
#line 1133 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 450:
#line 1135 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 457:
#line 1149 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 458:
#line 1151 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 459:
#line 1154 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 460:
#line 1155 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 461:
#line 1156 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 462:
#line 1157 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 463:
#line 1158 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 464:
#line 1159 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 465:
#line 1160 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 466:
#line 1161 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 467:
#line 1162 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 468:
#line 1163 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 469:
#line 1164 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 470:
#line 1165 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 471:
#line 1166 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 472:
#line 1167 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 473:
#line 1168 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 474:
#line 1169 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 475:
#line 1170 "DynareBison.yy"
    {driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 476:
#line 1171 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 477:
#line 1172 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 478:
#line 1173 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 479:
#line 1174 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 480:
#line 1176 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 481:
#line 1177 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 482:
#line 1178 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 483:
#line 1180 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 484:
#line 1181 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 485:
#line 1182 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 486:
#line 1183 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 487:
#line 1184 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 488:
#line 1185 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 489:
#line 1186 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 490:
#line 1187 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 491:
#line 1188 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 492:
#line 1189 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 493:
#line 1190 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1192 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 496:
#line 1193 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 497:
#line 1194 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 498:
#line 1195 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 499:
#line 1196 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 500:
#line 1197 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 501:
#line 1198 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 502:
#line 1199 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 503:
#line 1200 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 504:
#line 1201 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 505:
#line 1202 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 506:
#line 1203 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 507:
#line 1204 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 508:
#line 1205 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 509:
#line 1206 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 510:
#line 1207 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 511:
#line 1208 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 512:
#line 1209 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 513:
#line 1210 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 514:
#line 1211 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 515:
#line 1214 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 516:
#line 1216 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 517:
#line 1218 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 518:
#line 1219 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 519:
#line 1220 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 520:
#line 1221 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 521:
#line 1222 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 522:
#line 1223 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 523:
#line 1224 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 524:
#line 1225 "DynareBison.yy"
    {driver.option_num("mh_recover", "1");;}
    break;

  case 525:
#line 1226 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 526:
#line 1228 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 527:
#line 1229 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 528:
#line 1230 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1231 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1232 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1233 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 532:
#line 1234 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1235 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 534:
#line 1238 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 536:
#line 1247 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 537:
#line 1250 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 538:
#line 1252 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 539:
#line 1260 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2212 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -858;
  const short int
  parser::yypact_[] =
  {
       849,    35,    58,   -88,   -74,    77,   535,    22,    -6,    47,
     -44,     0,   -29,   -14,    -1,   101,   115,   585,   231,   129,
     162,   425,   321,   329,    53,    17,   342,   394,  -858,   355,
     369,    17,   380,   497,   421,   454,    62,   202,    17,   470,
     478,   500,    17,    71,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,   435,   265,
     451,   894,  -858,   611,    64,  -858,   508,   599,   487,    26,
     -26,   587,   283,   601,   603,   646,  -858,   731,   106,   190,
     259,   267,   604,   603,   647,   649,   549,  -858,   236,   302,
      69,   948,   630,   641,  -858,  1098,   292,   294,   595,   298,
     671,   567,   972,   768,   768,   389,    69,   577,  -858,   130,
    -858,   508,  -858,  1098,   391,  -858,   308,   399,   439,   618,
     442,   620,   449,   622,   464,   479,  -858,  -858,  -858,   716,
    -858,   717,   720,   721,   722,   729,   735,   746,   748,   750,
     758,   761,   772,  -858,   653,   650,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,   773,   783,   785,  -858,   685,   664,  -858,  -858,  -858,
     669,   740,    -4,    92,  -858,   799,   252,  -858,  -858,   681,
    -858,   687,  -858,  -858,   762,   218,  -858,   769,   367,   817,
      56,  -858,   775,  -858,  -858,   819,  -858,  -858,   820,   826,
     828,   832,   836,  -858,  -858,   838,   839,   841,   842,   843,
     847,  -858,  -858,   854,   862,  -858,  -858,  -858,  -858,   864,
     865,  -858,  -858,   253,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,   888,   813,  -858,   846,  -858,   848,
      98,  -858,   788,   850,   793,   858,   263,  -858,   868,   808,
     869,   322,  -858,   792,   289,  -858,   324,   918,   797,   804,
    -858,   860,  -858,   276,   803,   805,   932,  -858,  -858,   282,
    -858,  -858,  -858,  -858,   886,   889,    94,  -858,  -858,  -858,
     815,   948,   948,   818,   824,   825,   827,   829,   830,   831,
     833,   840,   853,   948,   536,   855,   344,  -858,   912,   351,
     938,   951,   953,   964,   967,   979,  -858,  -858,  -858,   982,
     984,   986,  -858,   988,  -858,   989,   991,   866,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,   906,   950,  -858,   874,  -858,  -858,
    -858,   875,   972,   972,   876,   877,   878,   879,   884,   887,
     891,   892,   893,   905,   972,   629,  -858,   381,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,   408,  -858,    95,    30,   411,  -858,  -858,  -858,
     413,  -858,  -858,   415,  -858,  -858,   999,  -858,   448,  -858,
    -858,  -858,  -858,  -858,   921,   987,  -858,  -858,   946,  1000,
    -858,  -858,   958,  1005,  -858,  -858,  1040,    93,  1045,  1047,
      93,  1048,  1020,  1051,    14,  1052,  1054,  1029,  1030,   265,
    1060,  1061,  1081,  1075,   894,  1076,   978,   971,  1063,   635,
    1108,  -858,  -858,  1094,   508,   992,  -858,  -858,   993,   158,
    1068,   995,   180,  1073,   948,  -858,  -858,  -858,   994,  1103,
    1106,  1109,  1116,  1122,  1120,   636,  1133,  1125,  1138,  1139,
    1141,  1114,  1001,  1146,   731,   206,  1118,  1164,  1062,  -858,
    -858,  -858,    60,  1066,    90,  1072,  -858,  -858,  1074,    90,
    1078,  -858,  -858,     2,  -858,  -858,  -858,  1127,  1067,  -858,
    1150,   156,  -858,   131,  -858,   530,  -858,  1082,  1134,   458,
     302,   184,  1084,    28,  -858,  -858,   948,  1096,   638,   948,
     948,   948,   948,   948,   948,   948,   948,   948,   948,   623,
     948,   948,   948,   948,   948,  -858,   948,  -858,  -858,  1182,
    1189,  -858,   930,  1151,  1223,  1263,  1266,  1289,  1292,  1293,
    1294,   661,  1296,  1313,  1369,   213,  -858,  1165,  -858,     2,
    1281,   703,   972,   972,   972,   972,   972,   972,   972,   972,
     972,   972,   689,   972,   972,   972,   972,   972,  1265,   768,
     261,   268,  -858,  -858,  -858,   948,   437,   175,   130,  1268,
     508,  1269,  1098,   286,  1381,   308,   287,  -858,  1304,  -858,
    1305,  -858,  1306,  -858,  -858,  1389,  1390,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  1382,    -3,  -858,  -858,  -858,
    -858,  1274,  -858,  -858,  1277,  -858,  -858,  -858,  -858,  1280,
    -858,  1387,  1282,  1283,  1284,   948,  -858,  -858,  -858,  -858,
    -858,   480,  1285,  -858,  -858,   485,  1286,  1195,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  1275,  -858,  -858,  -858,   486,  -858,  1361,
    1366,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,   665,
    1290,  1314,  1315,  1370,  1317,    90,  1372,  1297,    90,  -858,
    1402,  1404,  1298,  -858,   603,  1423,  -858,  -858,  -858,   972,
    -858,  -858,  -858,  1425,   463,  -858,  -858,  -858,  1302,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,    80,
     373,  -858,  1380,   948,  1383,   144,   780,   506,   699,   754,
     796,   903,   996,  1002,  1028,  1034,  1042,  1079,  -858,   638,
     638,  1096,  1096,  1321,  1085,   948,  -858,  1384,  1201,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,   491,  -858,  1303,  1091,  1097,  1105,  1111,
    1117,  1123,  1131,  1137,  1143,  1149,  -858,   703,   703,  1281,
    1281,  1321,  -858,  -858,  -858,   505,  -858,   548,  1157,    30,
    1308,  -858,  -858,    33,   948,  -858,  -858,  -858,  -858,  -858,
    -858,   558,  -858,  -858,  -858,   580,  -858,  -858,  -858,  -858,
    -858,  1307,  -858,  -858,  -858,  1386,  -858,  -858,  1310,  1434,
    -858,  -858,  1213,  -858,   319,  -858,   320,  -858,  1391,  -858,
     516,  -858,  -858,  -858,  -858,  -858,  -858,    90,    60,  1336,
      90,  1338,  1339,  -858,  1318,  -858,  -858,  1441,   354,   972,
    1219,  1435,   530,  -858,   860,   860,   860,   184,  -858,    90,
    -858,  1443,  1225,  1445,  1428,   948,   948,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  1323,  1231,
     948,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,   175,  -858,  -858,
    -858,   948,  1163,  -858,  -858,  1430,  -858,  1282,   948,  -858,
    -858,   589,  -858,   600,  1319,  1275,  -858,  -858,  1348,  1350,
    1351,    90,  1329,    90,    90,  -858,   948,  -858,  1237,  -858,
    -858,  -858,  1330,   182,   510,   553,   444,  1331,   948,  -858,
     948,  1327,    16,  1243,   780,  -858,  -858,  1249,  1169,  -858,
    -858,  1456,  1255,  -858,  -858,  1357,  -858,    90,    90,    90,
    1358,  -858,  1337,  1340,  1261,  -858,   860,  -858,  -858,  -858,
      90,  -858,  1267,  1273,  1442,  1334,  1447,  1373,  -858,  -858,
    -858,   948,  -858,    68,  1362,  -858,  1365,    90,  -858,  -858,
    -858,   584,  1342,  -858,  -858,  -858,  1451,  1343,   271,  1279,
    1426,  -858,    90,    70,  1346,  -858,  -858,  -858,  1454,  -858,
     666,   696,   948,    73,  -858,  -858,  -858,  1344,  1375,  1376,
    -858,  -858,  1175,  -858,  -858,   948,  -858,  -858,  -858,    90,
      90,  -858,  1183,  1377,  -858,  -858,    90,  -858
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   421,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     2,     4,    29,    30,    44,    45,
      46,    43,     5,     6,     7,    12,     9,    10,    11,     8,
      13,    14,    15,    16,    17,    18,    19,    23,    25,    24,
      20,    21,    22,    26,    27,    28,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,     0,     0,
       0,     0,   390,     0,     0,   206,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   250,   297,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   126,     0,     0,
       0,     0,     0,     0,   377,     0,     0,     0,    74,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   211,     0,
     201,     0,   217,     0,     0,   424,     0,     0,     0,    56,
       0,    62,     0,    68,     0,     0,     1,     3,   449,     0,
     531,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   439,   448,     0,   440,   445,   443,   446,
     444,   441,   442,   432,   433,   434,   435,   436,   437,   438,
     457,     0,     0,     0,   451,   456,     0,   453,   452,   454,
       0,     0,   387,     0,   383,     0,     0,   209,   210,     0,
      80,     0,    47,   400,     0,     0,   394,     0,     0,     0,
       0,   114,     0,   506,   522,     0,   511,   489,     0,     0,
       0,     0,     0,   503,   504,     0,     0,     0,     0,     0,
       0,   524,   499,     0,     0,   510,   523,   505,   488,     0,
       0,   509,   507,     0,   302,   338,   327,   303,   304,   305,
     306,   307,   308,   309,   310,   311,   312,   313,   314,   315,
     316,   317,   318,   319,   320,   321,   322,   323,   324,   325,
     326,   328,   329,   330,   331,   332,   333,   334,   335,   336,
     337,   339,   340,   341,   246,     0,   299,     0,   263,     0,
       0,   260,     0,     0,     0,     0,     0,   282,     0,     0,
       0,     0,   276,     0,     0,   118,     0,     0,     0,     0,
      82,     0,   462,     0,     0,     0,     0,   518,   517,     0,
     406,   407,   408,   409,     0,     0,     0,   169,    87,    88,
      86,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   368,     0,     0,
       0,     0,     0,     0,     0,     0,   467,   468,   469,     0,
       0,     0,   512,     0,   476,     0,     0,     0,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   235,
     237,   238,   239,   240,   241,   242,   243,   234,   236,   244,
     245,   379,   376,    77,    72,     0,    53,     0,    78,   144,
     145,   164,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   422,   143,     0,   345,   350,
     346,   347,   348,   349,   351,   352,   353,   354,   355,   356,
     357,   358,     0,    49,     0,     0,     0,   214,   215,   216,
       0,   204,   205,     0,   222,   219,     0,   430,     0,   429,
     431,   426,   370,    59,    54,     0,    50,    65,    60,     0,
      51,    71,    66,     0,    52,   365,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   381,   382,     0,     0,     0,    81,    48,     0,     0,
       0,     0,     0,     0,     0,   112,   113,   251,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   248,     0,   262,
     258,   259,   291,     0,   291,     0,   280,   281,     0,   291,
       0,   274,   275,     0,   116,   117,   109,     0,     0,    83,
       0,     0,   138,     0,   139,     0,   134,     0,     0,     0,
       0,     0,     0,     0,   167,   168,     0,    94,    95,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    84,     0,   366,   367,     0,
       0,   371,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    75,    73,    79,     0,
     151,   152,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   166,   199,   200,     0,     0,   191,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    57,    55,    63,
      61,    69,    67,   527,   257,     0,     0,   528,   529,   530,
     526,   532,   480,   483,   482,     0,     0,   481,   484,   485,
     519,     0,   520,   447,     0,   533,   490,   508,   455,     0,
     391,     0,   387,     0,     0,     0,   460,   208,   207,   403,
     398,     0,     0,   397,   392,     0,     0,     0,   521,   470,
     513,   514,   486,   487,   492,   495,   496,   493,   501,   502,
     491,   498,   497,     0,   500,   301,   298,     0,   247,     0,
       0,   286,   293,   287,   292,   289,   294,   288,   290,     0,
       0,     0,   268,     0,     0,   291,     0,     0,   291,   254,
       0,     0,     0,   111,     0,     0,   127,   136,   137,     0,
     141,   123,   122,     0,     0,   121,   124,   125,     0,   130,
     128,   515,   516,   405,   416,   418,   419,   420,   417,     0,
     410,   414,     0,     0,     0,     0,   107,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    85,    90,
      89,    91,    92,    93,     0,     0,   374,     0,     0,   466,
     474,   459,   465,   471,   472,   463,   473,   479,   478,   464,
     461,   477,   378,     0,    76,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   142,   147,   146,   148,
     149,   150,   423,   344,   342,     0,   359,     0,     0,     0,
       0,   196,   197,     0,     0,   213,   212,   203,   202,   221,
     218,     0,   525,   428,   425,     0,    58,    64,    70,   256,
     255,   535,   537,   539,   538,     0,   450,   458,     0,     0,
     389,   388,     0,   399,     0,   393,     0,   115,     0,   363,
       0,   300,   249,   264,   296,   295,   261,   291,   291,     0,
     291,     0,     0,   279,     0,   253,   252,     0,     0,     0,
       0,     0,     0,   132,     0,     0,     0,     0,   404,   291,
     415,     0,     0,     0,     0,     0,     0,   106,    96,    97,
      98,    99,   100,   101,   102,   103,   104,   105,     0,     0,
       0,   372,   380,   165,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,   343,   360,   198,   190,   189,   193,
     194,     0,     0,   220,   427,     0,   534,   387,     0,   384,
     401,     0,   395,     0,     0,     0,   494,   265,     0,     0,
       0,   291,     0,   291,   291,   277,     0,   110,     0,   140,
     475,   120,     0,     0,     0,     0,   411,     0,     0,   172,
       0,   180,     0,     0,   108,   369,   375,     0,     0,   195,
     536,     0,     0,   402,   396,     0,   364,   291,   291,   291,
       0,   285,     0,     0,     0,   163,     0,   135,   131,   129,
     291,   412,     0,     0,     0,   175,     0,     0,   171,   373,
     192,     0,   385,   291,   270,   266,   269,   291,   283,   278,
     119,     0,     0,   174,   173,   179,     0,   177,     0,     0,
       0,   362,   291,     0,     0,   133,   413,   176,     0,   186,
       0,     0,     0,     0,   185,   184,   386,     0,   271,     0,
     284,   178,     0,   183,   170,     0,   182,   181,   361,   291,
     291,   188,     0,   272,   267,   187,   291,   273
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -858,  -858,  1464,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -337,  -858,
    -858,  -858,  -858,  -107,  -218,  -858,  -858,  1198,  -858,   602,
    -858,  -858,  -858,  -858,  -858,  -858,  -822,  -533,  -129,  -530,
    -858,  -858,  -858,  1379,  -273,  -858,  -858,  -858,  -858,   667,
    -858,  -858,   863,  -858,  -858,  1013,  -858,  -858,   870,  -858,
    -858,   -99,   -24,  -592,  -472,  -858,  -858,  1220,  -858,  -858,
    -682,  -858,  -858,  1208,  -858,  -858,  1215,  -857,  -527,  -858,
    -858,   990,  -858,  1388,   890,  -858,   550,  -858,  -858,  -858,
    -858,  1167,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  1324,
    -657,  -858,  -858,  -858,  -858,  -858,   956,  -858,   613,  -744,
    -858,  -858,  -858,  -858,  -858,   873,  -858,   -58,  1043,  -858,
    -858,  1037,  -858,  -858,   -87,  -858,  1415,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,   -96,  -858,  -858,  -120,  -531,  -858,
    -858,  -858,  -858,   -97,   -67,   -64,   -62,   -60,  -858,  -858,
     -84,   -61,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
     -51,  -858,  -858,  -858,  -858,  -858,   -50,   -47,   -52,   -45,
     -43,   -25,  -858,  -858,  -858,  -858,   -95,   -90,   -82,   -75,
     -22,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,  -858,
    -858,  -858,  -858,  -858,   859,  -858,  1016
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    43,    44,    45,    46,    47,    48,    49,    50,    51,
     150,   152,   154,   129,    52,    53,    54,    55,   354,   787,
      56,   318,    57,   220,   221,    58,   314,   315,   764,   765,
      59,   321,   916,   915,   992,   768,   561,   562,   563,   564,
     426,    60,    61,   336,   337,  1002,  1073,    62,   646,   647,
      63,   450,   451,    64,   206,   207,    65,   446,   447,    66,
     453,   457,   108,   752,   667,    67,   300,   301,   302,   740,
     977,    68,   311,   312,    69,   306,   307,   741,   978,    70,
     253,   254,    71,   427,   428,    72,   889,   890,    73,    74,
     356,   357,    75,    76,   359,    77,    78,    79,   203,   204,
     500,    80,    81,    82,    83,   329,   330,   779,   780,   781,
      84,   132,   638,    85,   458,   459,   173,   174,   175,    86,
     195,   196,    87,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   767,
     394,   395,   396,   176,   177,   178,   179,   180,   262,   263,
     397,   431,   266,   267,   268,   269,   270,   271,   272,   273,
     432,   275,   276,   277,   278,   279,   433,   434,   435,   436,
     437,   438,   398,   286,   287,   331,   399,   400,   181,   182,
     441,   291,   292,   293,   460,   183,   184,   185,   186,   187,
     188,   189,   199,   682,   872,   676,   677
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       126,   127,   516,   425,   577,   578,   316,   135,   670,   208,
     257,   256,   144,   147,   148,   742,   589,   744,   155,   448,
     255,   600,   747,   264,   332,   288,   377,   825,   757,   333,
     197,   758,   289,   194,   766,   879,   920,   429,   429,   198,
     258,   979,   449,   259,   454,   260,   265,   261,   871,   430,
     430,   439,   439,   749,   452,   282,   274,   280,   440,   440,
     281,   783,   283,   575,   284,   674,   731,  1035,   201,    92,
      93,   156,   643,    99,   100,   959,   731,   211,     1,     2,
     105,   644,   285,    94,   960,   290,    88,   515,     3,     4,
       5,   334,   993,   994,   995,     6,   201,   105,   516,     7,
       8,     9,   732,    10,   733,    11,    12,    13,    14,    90,
     732,   734,   735,   104,   733,   664,   334,   334,    15,   734,
     297,    16,   735,   501,   749,   574,   642,   102,   109,   540,
     498,   213,   732,   105,    17,   664,   219,   750,   751,   214,
     736,   734,   105,   110,   202,    18,    19,    20,   736,   737,
    1036,    21,   101,  1083,   499,   166,   111,   106,   107,   737,
     873,    22,   784,    23,   759,    24,    25,    26,    27,    28,
     736,   575,   202,  1037,    29,    30,   675,   707,   298,    31,
      32,    33,    34,   212,   738,   785,   294,   756,   645,    35,
      36,   961,    37,    89,   738,   739,    38,   335,   409,    39,
      40,    41,    42,   739,  1051,   103,  1079,   410,  1070,  1071,
     124,   125,   297,  1027,   917,   299,    91,   851,   901,   142,
     143,   904,   335,   335,   409,   739,   852,   370,   665,   666,
    1084,  1085,  1060,   410,    95,    96,   411,   918,   105,   786,
     295,   924,   788,   789,   790,   791,   792,   793,   794,   795,
     796,   797,   920,   799,   800,   801,   802,   803,   112,   804,
     105,   925,   411,   296,   774,   808,   633,   634,   635,   636,
     298,   637,   113,   114,   159,   160,   161,   162,   163,   164,
     165,   303,   105,   620,   621,   303,   105,   119,   760,   308,
     166,   412,   413,   105,   546,   632,   322,   414,   415,   416,
     417,   418,   419,   420,   421,   422,   167,   299,   848,   853,
    1011,   360,   423,   664,   424,   700,   560,   412,   413,   120,
     554,   775,   749,   414,   415,   416,   417,   418,   419,   420,
     421,   422,   361,   854,   362,   363,  1045,   704,   423,   304,
     424,   105,   560,   304,   308,   776,   168,   309,   105,   777,
     778,  1069,   323,   551,   227,   556,   364,   365,   882,   145,
     146,   228,   324,   726,   169,   170,   105,   105,   322,   313,
     822,   980,   294,   982,   294,   597,   305,   509,   404,   326,
     305,   766,   601,   510,   310,   987,   504,   534,   117,   118,
     327,   366,   997,   367,   248,   368,   327,   171,   172,   105,
     105,   369,   309,   328,   219,   370,  1070,  1071,  1094,   328,
     565,   505,   535,   371,   372,   373,   570,   208,   844,   374,
     375,   376,   128,   205,   355,   846,   295,   130,   295,  1072,
     456,   602,   405,   197,   219,   566,   194,   257,   256,   310,
     216,   571,   198,   860,   864,   131,   922,   255,   217,   401,
     264,   402,   288,   774,  1020,   406,  1022,  1023,   121,   289,
     757,   757,   757,   758,   758,   758,  1074,   258,   939,   294,
     259,   294,   260,   265,   261,   332,   970,   972,   122,   294,
     333,  1086,   282,   274,   280,   701,   123,   281,   705,   283,
    1044,   284,  1046,   826,   827,   828,   829,   830,   831,   832,
     833,   834,   835,  1052,   837,   838,   839,   840,   841,   285,
     775,   727,   290,   133,   771,   639,  1061,   962,   757,   294,
    1064,   758,   464,   295,   774,   295,   512,   134,   448,   468,
     137,   919,   513,   295,   776,  1078,   772,   136,   777,   778,
     640,  1028,   639,   429,   472,   648,   443,   650,   455,   652,
     149,   449,   409,   859,   361,   430,   461,   439,   151,   294,
     294,   410,  1093,   452,   440,   294,   294,   641,    97,  1097,
     649,   294,   651,   295,   653,   761,   465,    98,   138,   139,
     153,   775,   655,   469,  1029,   294,  1075,   762,  1003,  1004,
     411,   823,   158,   763,   849,   409,   462,   912,   473,   466,
     850,  1087,  1030,  1007,   410,   776,   470,   656,   190,   777,
     778,   140,   141,   295,   295,  1065,   845,   847,   115,   295,
     295,   474,   913,   205,  1008,   295,   409,   116,   294,   861,
     910,  1012,   865,   411,   200,   410,   475,   883,   294,   295,
     926,   209,   885,   891,   210,   412,   413,   908,   942,  1024,
     975,   414,   415,   416,   417,   418,   419,   420,   421,   422,
     294,  1032,   954,  1033,   411,   927,   423,   215,   424,   294,
     560,   590,   591,   592,   593,   976,   594,   693,   715,   222,
     294,   218,   295,   219,   313,   317,   694,   716,   412,   413,
     516,   319,   295,   595,   414,   415,   416,   417,   418,   419,
     420,   421,   422,   817,  1059,   955,   320,   894,   869,   423,
     355,   424,   818,   560,   295,   963,   895,   905,   403,   412,
     413,   358,   407,   295,   408,   414,   415,   416,   417,   418,
     419,   420,   421,   422,   295,  1082,   223,   964,   870,   445,
     423,   463,   424,   467,   560,   471,  1013,   906,  1092,   476,
     477,   192,   224,   478,   479,   480,   166,  1014,   590,   591,
     592,   593,   481,   594,   633,   634,   635,   636,   482,   637,
     225,   226,   167,   223,   193,   592,   593,   227,   594,   483,
     988,   484,   798,   485,   228,   229,   230,   489,   192,   231,
     232,   486,   233,   234,   487,   235,   236,   237,   238,   239,
     240,   241,   242,   243,   244,   488,   491,   225,   226,   490,
     245,   193,   168,   246,   227,   247,   492,   248,   493,   494,
     497,   228,   249,   495,   633,   634,   635,   636,   496,   637,
     169,   170,   503,   250,   590,   591,   592,   593,   506,   594,
     635,   636,   508,   637,   507,   251,   205,   245,   836,   511,
     514,   252,   518,   519,   248,   517,     1,     2,   928,   520,
     971,   521,   973,   171,   172,   522,     3,     4,     5,   523,
     250,   524,   525,     6,   526,   527,   528,     7,     8,     9,
     529,    10,   251,    11,    12,    13,    14,   530,   252,   590,
     591,   592,   593,   537,   594,   531,    15,   532,   533,    16,
     171,   172,   409,   159,   160,   161,   162,   163,   164,   165,
     191,   410,    17,   929,   192,   590,   591,   592,   593,   166,
     594,   536,   542,    18,    19,    20,   538,   544,   539,    21,
     543,   590,   591,   592,   593,   167,   594,   193,   545,    22,
     411,    23,   549,    24,    25,    26,    27,    28,   548,   550,
     553,   557,    29,    30,   338,   930,   558,    31,    32,    33,
      34,   559,   567,   339,   568,   569,   572,    35,    36,   573,
      37,   603,   338,   576,    38,   168,   579,    39,    40,    41,
      42,   339,   580,   581,   604,   582,   605,   583,   584,   585,
     338,   586,   340,   169,   170,   412,   413,   606,   587,   339,
     607,   414,   415,   416,   417,   418,   419,   420,   421,   422,
     340,   588,   608,   596,   409,   609,   423,   610,   424,   611,
     560,   612,   613,   410,   614,   615,   171,   172,   340,   616,
     617,   618,   654,   619,   622,   623,   624,   625,   590,   591,
     592,   593,   626,   594,   657,   627,   599,   341,   342,   628,
     629,   630,   411,   343,   344,   345,   346,   347,   348,   349,
     350,   351,   931,   631,   807,   341,   342,   658,   352,   659,
     353,   343,   344,   345,   346,   347,   348,   349,   350,   351,
     660,   661,   663,   341,   342,   662,   352,   668,   353,   343,
     344,   345,   346,   347,   348,   349,   350,   351,   669,   671,
     672,   360,   673,   678,   352,   679,   353,   412,   413,   680,
     681,   684,   685,   414,   415,   416,   417,   418,   419,   420,
     421,   422,   361,   686,   362,   363,   687,   689,   423,   691,
     424,   590,   591,   592,   593,   690,   594,   590,   591,   592,
     593,   695,   594,   692,   227,   696,   364,   365,   702,   698,
     699,   228,   703,   706,   709,   932,   675,   710,   322,   723,
     711,   933,   714,   590,   591,   592,   593,   712,   594,   590,
     591,   592,   593,   713,   594,   717,   718,   590,   591,   592,
     593,   366,   594,   367,   248,   368,   327,   934,   724,   719,
     720,   369,   721,   935,   722,   370,   730,   729,   728,   328,
     743,   936,   809,   371,   372,   373,   745,   753,   746,   374,
     375,   376,   748,   205,   590,   591,   592,   593,   782,   594,
     590,   591,   592,   593,   754,   594,   633,   634,   635,   636,
     755,   637,   633,   634,   635,   636,   594,   637,   937,   769,
     633,   634,   635,   636,   938,   637,   633,   634,   635,   636,
     944,   637,   633,   634,   635,   636,   945,   637,   633,   634,
     635,   636,   805,   637,   946,   810,   633,   634,   635,   636,
     947,   637,   633,   634,   635,   636,   948,   637,   633,   634,
     635,   636,   949,   637,   633,   634,   635,   636,   824,   637,
     950,   770,   590,   591,   592,   593,   951,   594,   590,   591,
     592,   593,   952,   594,   590,   591,   592,   593,   953,   594,
     590,   591,   592,   593,   811,   594,   956,   812,   590,   591,
     592,   593,  1009,   594,   590,   591,   592,   593,  1040,   594,
     590,   591,   592,   593,  1091,   594,   590,   591,   592,   593,
     813,   594,  1095,   814,   815,   816,   806,   819,   590,   591,
     592,   593,   887,   594,   633,   634,   635,   636,   941,   637,
     590,   591,   592,   593,   820,   594,   590,   591,   592,   593,
     969,   594,   633,   634,   635,   636,   989,   637,   590,   591,
     592,   593,   999,   594,   590,   591,   592,   593,  1006,   594,
     590,   591,   592,   593,  1025,   594,   590,   591,   592,   593,
    1038,   594,   590,   591,   592,   593,  1039,   594,   590,   591,
     592,   593,  1042,   594,   590,   591,   592,   593,  1050,   594,
     821,   637,   842,   862,  1053,   856,   858,   866,   867,   868,
    1054,   869,   870,   871,   876,   875,  1076,   877,   878,   888,
     499,   892,   880,   881,   884,   886,   893,   896,   897,   898,
     899,   900,   902,   905,   903,   906,   909,   907,   911,   914,
     921,    -1,   943,   923,   940,   958,   966,   968,   965,   967,
     981,   974,   983,   984,   986,   985,   998,   990,  1000,  1001,
    1005,  1010,  1017,  1015,  1018,  1019,  1021,  1026,  1034,  1041,
    1031,  1043,  1047,  1055,  1048,  1056,  1062,  1049,  1057,  1063,
    1058,  1066,  1067,  1080,  1068,  1081,  1077,   157,  1088,  1089,
    1090,  1096,   555,   857,   991,   444,   957,   697,   855,   552,
     541,   547,   442,   598,   725,  1016,   773,   502,   863,   843,
     996,   688,   683,   325,   708,   874
  };

  /* YYCHECK.  */
  const unsigned short int
  parser::yycheck_[] =
  {
        24,    25,   220,   132,   341,   342,   113,    31,   480,    96,
     107,   107,    36,    37,    38,   542,   353,   544,    42,   139,
     107,   358,   549,   107,   119,   107,   125,   619,   561,   119,
      91,   561,   107,    91,   565,   692,   780,   133,   134,    91,
     107,   898,   139,   107,   143,   107,   107,   107,    51,   133,
     134,   133,   134,    51,   141,   107,   107,   107,   133,   134,
     107,    33,   107,   336,   107,    51,     6,    51,     4,   157,
     158,     0,    42,    51,    80,    42,     6,    51,     7,     8,
      80,    51,   107,   157,    51,   107,    51,    31,    17,    18,
      19,    22,   914,   915,   916,    24,     4,    80,   316,    28,
      29,    30,    42,    32,    44,    34,    35,    36,    37,    51,
      42,    51,    52,   157,    44,    42,    22,    22,    47,    51,
      22,    50,    52,    31,    51,    31,    31,    80,   157,    31,
     134,   157,    42,    80,    63,    42,    80,   135,   136,   165,
      80,    51,    80,   157,    80,    74,    75,    76,    80,    89,
     134,    80,   158,    80,   158,    25,   157,   157,   158,    89,
     163,    90,   134,    92,    33,    94,    95,    96,    97,    98,
      80,   444,    80,   157,   103,   104,   162,   514,    80,   108,
     109,   110,   111,   157,   124,   157,    80,    31,   158,   118,
     119,   158,   121,   158,   124,   135,   125,   128,    42,   128,
     129,   130,   131,   135,  1026,   158,  1063,    51,   135,   136,
     157,   158,    22,    31,   134,   117,   158,    42,   745,   157,
     158,   748,   128,   128,    42,   135,    51,    97,   135,   136,
     157,   158,   164,    51,   157,   158,    80,   157,    80,   576,
     134,    97,   579,   580,   581,   582,   583,   584,   585,   586,
     587,   588,   996,   590,   591,   592,   593,   594,   157,   596,
      80,   117,    80,   157,    80,   602,   135,   136,   137,   138,
      80,   140,   157,   158,     9,    10,    11,    12,    13,    14,
      15,    22,    80,   412,   413,    22,    80,   158,   157,    22,
      25,   135,   136,    80,    31,   424,    60,   141,   142,   143,
     144,   145,   146,   147,   148,   149,    41,   117,   645,   134,
     967,     3,   156,    42,   158,   157,   160,   135,   136,   157,
      31,   137,    51,   141,   142,   143,   144,   145,   146,   147,
     148,   149,    24,   158,    26,    27,  1018,   157,   156,    80,
     158,    80,   160,    80,    22,   161,    81,    80,    80,   165,
     166,    80,   116,    31,    46,    31,    48,    49,   695,   157,
     158,    53,   126,   157,    99,   100,    80,    80,    60,    80,
     157,   898,    80,   900,    80,    31,   117,   159,    80,    77,
     117,   912,    31,   165,   117,    31,   134,   134,   157,   158,
      88,    83,   919,    85,    86,    87,    88,   132,   133,    80,
      80,    93,    80,   101,    80,    97,   135,   136,  1090,   101,
     134,   159,   159,   105,   106,   107,   134,   504,   157,   111,
     112,   113,    80,   115,    80,   157,   134,    33,   134,   158,
     122,    80,   134,   494,    80,   159,   494,   534,   534,   117,
     157,   159,   494,   157,   157,    51,   783,   534,   165,   157,
     534,   157,   534,    80,   981,   157,   983,   984,    33,   534,
     993,   994,   995,   993,   994,   995,  1058,   534,   805,    80,
     534,    80,   534,   534,   534,   570,   157,   157,   157,    80,
     570,  1073,   534,   534,   534,   509,   157,   534,   512,   534,
    1017,   534,  1019,   622,   623,   624,   625,   626,   627,   628,
     629,   630,   631,  1030,   633,   634,   635,   636,   637,   534,
     137,   535,   534,   158,    56,   134,  1043,   854,  1051,    80,
    1047,  1051,    80,   134,    80,   134,   159,   158,   648,    80,
      33,   158,   165,   134,   161,  1062,    78,   157,   165,   166,
     159,    31,   134,   639,    80,   134,   157,   134,   157,   134,
      80,   648,    42,   652,    24,   639,   157,   639,    80,    80,
      80,    51,  1089,   650,   639,    80,    80,   159,    33,  1096,
     159,    80,   159,   134,   159,    45,   134,    42,   157,   158,
      80,   137,   134,   134,    31,    80,  1058,    57,   925,   926,
      80,   615,   157,    63,   157,    42,   157,   134,   134,   157,
     163,  1073,   158,   940,    51,   161,   157,   159,   157,   165,
     166,   157,   158,   134,   134,    31,   640,   641,    33,   134,
     134,   157,   159,   115,   961,   134,    42,    42,    80,   653,
     759,   968,   656,    80,    23,    51,   157,   157,    80,   134,
     134,    42,   157,   157,   157,   135,   136,   754,   157,   986,
     134,   141,   142,   143,   144,   145,   146,   147,   148,   149,
      80,   998,   157,  1000,    80,   159,   156,    80,   158,    80,
     160,   135,   136,   137,   138,   159,   140,    42,    42,    33,
      80,    80,   134,    80,    80,    38,    51,    51,   135,   136,
     908,    42,   134,   157,   141,   142,   143,   144,   145,   146,
     147,   148,   149,    42,  1041,   157,   157,    42,    42,   156,
      80,   158,    51,   160,   134,   157,    51,    51,   123,   135,
     136,    80,    51,   134,   157,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   134,  1072,     5,   157,    42,   162,
     156,   123,   158,   123,   160,   123,   157,    51,  1085,    33,
      33,    20,    21,    33,    33,    33,    25,   157,   135,   136,
     137,   138,    33,   140,   135,   136,   137,   138,    33,   140,
      39,    40,    41,     5,    43,   137,   138,    46,   140,    33,
     909,    33,   159,    33,    53,    54,    55,   134,    20,    58,
      59,    33,    61,    62,    33,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    33,    33,    39,    40,   159,
      79,    43,    81,    82,    46,    84,    33,    86,    33,   134,
      80,    53,    91,   159,   135,   136,   137,   138,   159,   140,
      99,   100,    33,   102,   135,   136,   137,   138,   157,   140,
     137,   138,    80,   140,   157,   114,   115,    79,   159,    80,
      33,   120,    33,    33,    86,    80,     7,     8,   159,    33,
     884,    33,   886,   132,   133,    33,    17,    18,    19,    33,
     102,    33,    33,    24,    33,    33,    33,    28,    29,    30,
      33,    32,   114,    34,    35,    36,    37,    33,   120,   135,
     136,   137,   138,    80,   140,    33,    47,    33,    33,    50,
     132,   133,    42,     9,    10,    11,    12,    13,    14,    15,
      16,    51,    63,   159,    20,   135,   136,   137,   138,    25,
     140,    33,   134,    74,    75,    76,    80,   134,    80,    80,
      80,   135,   136,   137,   138,    41,   140,    43,    80,    90,
      80,    92,   134,    94,    95,    96,    97,    98,    80,    80,
     158,    33,   103,   104,    42,   159,   159,   108,   109,   110,
     111,   157,   159,    51,   159,    33,    80,   118,   119,    80,
     121,    33,    42,   158,   125,    81,   158,   128,   129,   130,
     131,    51,   158,   158,    33,   158,    33,   158,   158,   158,
      42,   158,    80,    99,   100,   135,   136,    33,   158,    51,
      33,   141,   142,   143,   144,   145,   146,   147,   148,   149,
      80,   158,    33,   158,    42,    33,   156,    33,   158,    33,
     160,    33,    33,    51,    33,   159,   132,   133,    80,   123,
      80,   157,    33,   158,   158,   158,   158,   158,   135,   136,
     137,   138,   158,   140,   123,   158,   134,   135,   136,   158,
     158,   158,    80,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   159,   158,   134,   135,   136,    80,   156,   123,
     158,   141,   142,   143,   144,   145,   146,   147,   148,   149,
      80,   123,    42,   135,   136,    80,   156,    42,   158,   141,
     142,   143,   144,   145,   146,   147,   148,   149,    51,    51,
      80,     3,    51,    51,   156,    51,   158,   135,   136,    80,
      80,    51,    51,   141,   142,   143,   144,   145,   146,   147,
     148,   149,    24,    42,    26,    27,    51,    51,   156,   158,
     158,   135,   136,   137,   138,   157,   140,   135,   136,   137,
     138,    33,   140,    80,    46,    51,    48,    49,    80,   157,
     157,    53,   157,    80,    51,   159,   162,    51,    60,   158,
      51,   159,    42,   135,   136,   137,   138,    51,   140,   135,
     136,   137,   138,    51,   140,    42,    51,   135,   136,   137,
     138,    83,   140,    85,    86,    87,    88,   159,    42,    51,
      51,    93,    51,   159,    80,    97,   134,    33,    80,   101,
     134,   159,    51,   105,   106,   107,   134,    80,   134,   111,
     112,   113,   134,   115,   135,   136,   137,   138,   134,   140,
     135,   136,   137,   138,   157,   140,   135,   136,   137,   138,
      80,   140,   135,   136,   137,   138,   140,   140,   159,   157,
     135,   136,   137,   138,   159,   140,   135,   136,   137,   138,
     159,   140,   135,   136,   137,   138,   159,   140,   135,   136,
     137,   138,    80,   140,   159,    42,   135,   136,   137,   138,
     159,   140,   135,   136,   137,   138,   159,   140,   135,   136,
     137,   138,   159,   140,   135,   136,   137,   138,   123,   140,
     159,   157,   135,   136,   137,   138,   159,   140,   135,   136,
     137,   138,   159,   140,   135,   136,   137,   138,   159,   140,
     135,   136,   137,   138,    51,   140,   159,    51,   135,   136,
     137,   138,   159,   140,   135,   136,   137,   138,   159,   140,
     135,   136,   137,   138,   159,   140,   135,   136,   137,   138,
      51,   140,   159,    51,    51,    51,   157,    51,   135,   136,
     137,   138,   157,   140,   135,   136,   137,   138,   157,   140,
     135,   136,   137,   138,    51,   140,   135,   136,   137,   138,
     157,   140,   135,   136,   137,   138,   157,   140,   135,   136,
     137,   138,   157,   140,   135,   136,   137,   138,   157,   140,
     135,   136,   137,   138,   157,   140,   135,   136,   137,   138,
     157,   140,   135,   136,   137,   138,   157,   140,   135,   136,
     137,   138,   157,   140,   135,   136,   137,   138,   157,   140,
      51,   140,   157,    42,   157,   157,   157,   123,   123,   123,
     157,    42,    42,    51,   157,   161,   157,   157,    51,   164,
     158,    80,   159,   159,   159,   159,    80,   157,   134,   134,
      80,   134,    80,    51,   157,    51,    33,   159,    33,   157,
      80,   140,   159,    80,    80,   157,    80,    33,   161,   159,
     134,    80,   134,   134,    33,   157,    33,    42,    33,    51,
     157,    51,   134,   164,   134,   134,   157,   157,   161,    33,
     159,   134,   134,    51,   157,   161,   134,   157,    51,   134,
     127,   159,    51,   157,   161,    51,    80,    43,   164,   134,
     134,   134,   314,   650,   912,   136,   849,   504,   648,   311,
     300,   306,   134,   356,   534,   975,   570,   203,   655,   639,
     917,   494,   489,   118,   518,   676
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
     129,   130,   131,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   181,   182,   183,   184,   187,   189,   192,   197,
     208,   209,   214,   217,   220,   223,   226,   232,   238,   241,
     246,   249,   252,   255,   256,   259,   260,   262,   263,   264,
     268,   269,   270,   271,   277,   280,   286,   289,    51,   158,
      51,   158,   157,   158,   157,   157,   158,    33,    42,    51,
      80,   158,    80,   158,   157,    80,   157,   158,   229,   157,
     157,   157,   157,   157,   158,    33,    42,   157,   158,   158,
     157,    33,   157,   157,   157,   158,   229,   229,    80,   180,
      33,    51,   278,   158,   158,   229,   157,    33,   157,   158,
     157,   158,   157,   158,   229,   157,   158,   229,   229,    80,
     177,    80,   178,    80,   179,   229,     0,   169,   157,     9,
      10,    11,    12,    13,    14,    15,    25,    41,    81,    99,
     100,   132,   133,   283,   284,   285,   310,   311,   312,   313,
     314,   345,   346,   352,   353,   354,   355,   356,   357,   358,
     157,    16,    20,    43,   284,   287,   288,   318,   335,   359,
      23,     4,    80,   265,   266,   115,   221,   222,   291,    42,
     157,    51,   157,   157,   165,    80,   157,   165,    80,    80,
     190,   191,    33,     5,    21,    39,    40,    46,    53,    54,
      55,    58,    59,    61,    62,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    79,    82,    84,    86,    91,
     102,   114,   120,   247,   248,   291,   301,   310,   311,   312,
     313,   314,   315,   316,   317,   318,   319,   320,   321,   322,
     323,   324,   325,   326,   327,   328,   329,   330,   331,   332,
     333,   334,   335,   336,   337,   338,   340,   341,   345,   346,
     347,   348,   349,   350,    80,   134,   157,    22,    80,   117,
     233,   234,   235,    22,    80,   117,   242,   243,    22,    80,
     117,   239,   240,    80,   193,   194,   190,    38,   188,    42,
     157,   198,    60,   116,   126,   293,    77,    88,   101,   272,
     273,   342,   343,   344,    22,   128,   210,   211,    42,    51,
      80,   135,   136,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   156,   158,   185,    80,   257,   258,    80,   261,
       3,    24,    26,    27,    48,    49,    83,    85,    87,    93,
      97,   105,   106,   107,   111,   112,   113,   228,   290,   291,
     292,   293,   294,   295,   296,   297,   298,   299,   300,   301,
     302,   303,   304,   305,   307,   308,   309,   317,   339,   343,
     344,   157,   157,   123,    80,   134,   157,    51,   157,    42,
      51,    80,   135,   136,   141,   142,   143,   144,   145,   146,
     147,   148,   149,   156,   158,   205,   207,   250,   251,   301,
     317,   318,   327,   333,   334,   335,   336,   337,   338,   345,
     346,   347,   250,   157,   210,   162,   224,   225,   304,   310,
     218,   219,   291,   227,   228,   157,   122,   228,   281,   282,
     351,   157,   157,   123,    80,   134,   157,   123,    80,   134,
     157,   123,    80,   134,   157,   157,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,   134,
     159,    33,    33,    33,   134,   159,   159,    80,   134,   158,
     267,    31,   266,    33,   134,   159,   157,   157,    80,   159,
     165,    80,   159,   165,    33,    31,   191,    80,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,   134,   159,    33,    80,    80,    80,
      31,   234,   134,    80,   134,    80,    31,   243,    80,   134,
      80,    31,   240,   158,    31,   194,    31,    33,   159,   157,
     160,   203,   204,   205,   206,   134,   159,   159,   159,    33,
     134,   159,    80,    80,    31,   211,   158,   185,   185,   158,
     158,   158,   158,   158,   158,   158,   158,   158,   158,   185,
     135,   136,   137,   138,   140,   157,   158,    31,   258,   134,
     185,    31,    80,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,   159,   123,    80,   157,   158,
     205,   205,   158,   158,   158,   158,   158,   158,   158,   158,
     158,   158,   205,   135,   136,   137,   138,   140,   279,   134,
     159,   159,    31,    42,    51,   158,   215,   216,   134,   159,
     134,   159,   134,   159,    33,   134,   159,   123,    80,   123,
      80,   123,    80,    42,    42,   135,   136,   231,    42,    51,
     231,    51,    80,    51,    51,   162,   362,   363,    51,    51,
      80,    80,   360,   285,    51,    51,    42,    51,   288,    51,
     157,   158,    80,    42,    51,    33,    51,   222,   157,   157,
     157,   229,    80,   157,   157,   229,    80,   185,   363,    51,
      51,    51,    51,    51,    42,    42,    51,    42,    51,    51,
      51,    51,    80,   158,    42,   248,   157,   229,    80,    33,
     134,     6,    42,    44,    51,    52,    80,    89,   124,   135,
     236,   244,   245,   134,   245,   134,   134,   245,   134,    51,
     135,   136,   230,    80,   157,    80,    31,   204,   206,    33,
     157,    45,    57,    63,   195,   196,   305,   306,   202,   157,
     157,    56,    78,   273,    80,   137,   161,   165,   166,   274,
     275,   276,   134,    33,   134,   157,   185,   186,   185,   185,
     185,   185,   185,   185,   185,   185,   185,   185,   159,   185,
     185,   185,   185,   185,   185,    80,   157,   134,   185,    51,
      42,    51,    51,    51,    51,    51,    51,    42,    51,    51,
      51,    51,   157,   229,   123,   230,   205,   205,   205,   205,
     205,   205,   205,   205,   205,   205,   159,   205,   205,   205,
     205,   205,   157,   251,   157,   229,   157,   229,   185,   157,
     163,    42,    51,   134,   158,   225,   157,   219,   157,   228,
     157,   229,    42,   282,   157,   229,   123,   123,   123,    42,
      42,    51,   361,   163,   361,   161,   157,   157,    51,   267,
     159,   159,   185,   157,   159,   157,   159,   157,   164,   253,
     254,   157,    80,    80,    42,    51,   157,   134,   134,    80,
     134,   245,    80,   157,   245,    51,    51,   159,   190,    33,
     205,    33,   134,   159,   157,   200,   199,   134,   157,   158,
     276,    80,   185,    80,    97,   117,   134,   159,   159,   159,
     159,   159,   159,   159,   159,   159,   159,   159,   159,   185,
      80,   157,   157,   159,   159,   159,   159,   159,   159,   159,
     159,   159,   159,   159,   157,   157,   159,   216,   157,    42,
      51,   158,   185,   157,   157,   161,    80,   159,    33,   157,
     157,   229,   157,   229,    80,   134,   159,   237,   245,   244,
     245,   134,   245,   134,   134,   157,    33,    31,   205,   157,
      42,   196,   201,   203,   203,   203,   275,   245,    33,   157,
      33,    51,   212,   185,   185,   157,   157,   185,   185,   159,
      51,   267,   185,   157,   157,   164,   253,   134,   134,   134,
     245,   157,   245,   245,   185,   157,   157,    31,    31,    31,
     158,   159,   185,   185,   161,    51,   134,   157,   157,   157,
     159,    33,   157,   134,   245,   237,   245,   134,   157,   157,
     157,   203,   245,   157,   157,    51,   161,    51,   127,   185,
     164,   245,   134,   134,   245,    31,   159,    51,   161,    80,
     135,   136,   158,   213,   230,   231,   157,    80,   245,   244,
     157,    51,   185,    80,   157,   158,   230,   231,   164,   134,
     134,   159,   185,   245,   237,   159,   134,   245
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
     405,   406,   407,   408,   409,   410,   411,    59,    40,    41,
      35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   167,   168,   168,   169,   169,   169,   169,   169,   169,
     169,   169,   169,   169,   169,   169,   169,   169,   169,   169,
     169,   169,   169,   169,   169,   169,   169,   169,   169,   169,
     169,   169,   169,   169,   169,   169,   169,   169,   169,   169,
     169,   169,   169,   170,   170,   170,   170,   171,   171,   172,
     173,   174,   175,   176,   177,   177,   177,   177,   177,   177,
     178,   178,   178,   178,   178,   178,   179,   179,   179,   179,
     179,   179,   180,   180,   180,   180,   180,   180,   181,   181,
     182,   182,   183,   183,   184,   185,   185,   185,   185,   185,
     185,   185,   185,   185,   185,   185,   185,   185,   185,   185,
     185,   185,   185,   185,   185,   185,   185,   186,   186,   187,
     187,   188,   189,   190,   190,   191,   192,   193,   193,   194,
     195,   195,   196,   196,   196,   196,   198,   197,   199,   197,
     200,   197,   201,   197,   202,   197,   203,   203,   203,   203,
     204,   204,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   205,   205,   205,   205,   205,   205,   205,   205,   205,
     205,   205,   205,   206,   207,   207,   208,   209,   210,   210,
     211,   211,   211,   211,   211,   212,   212,   212,   212,   212,
     212,   213,   213,   213,   213,   213,   213,   213,   213,   214,
     215,   215,   216,   216,   216,   216,   216,   216,   216,   216,
     216,   217,   217,   218,   218,   219,   220,   220,   221,   221,
     222,   223,   223,   224,   224,   225,   225,   226,   226,   226,
     226,   227,   227,   228,   228,   228,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   228,   228,   228,   228,
     228,   228,   228,   228,   228,   228,   229,   229,   229,   229,
     229,   229,   230,   230,   230,   231,   231,   231,   232,   233,
     233,   234,   235,   235,   235,   236,   236,   236,   236,   236,
     237,   237,   237,   237,   238,   239,   239,   240,   240,   240,
     241,   242,   242,   243,   243,   243,   244,   244,   244,   244,
     244,   245,   245,   245,   245,   245,   245,   246,   246,   246,
     246,   247,   247,   248,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   248,   248,   248,   248,
     248,   248,   248,   248,   248,   248,   248,   248,   248,   248,
     248,   248,   249,   249,   250,   250,   251,   251,   251,   251,
     251,   251,   251,   251,   251,   251,   251,   251,   251,   252,
     252,   253,   253,   254,   254,   255,   256,   257,   257,   258,
     259,   260,   261,   261,   261,   261,   262,   263,   263,   263,
     263,   264,   265,   265,   266,   266,   266,   267,   267,   267,
     268,   268,   269,   269,   269,   269,   269,   269,   270,   270,
     270,   270,   270,   270,   271,   272,   272,   273,   273,   273,
     274,   274,   274,   274,   275,   275,   276,   276,   276,   276,
     276,   278,   279,   277,   280,   280,   280,   280,   281,   281,
     282,   282,   283,   283,   283,   283,   283,   283,   283,   284,
     284,   284,   284,   284,   284,   284,   284,   285,   285,   286,
     286,   287,   287,   287,   287,   288,   288,   289,   289,   290,
     291,   292,   293,   294,   295,   296,   297,   298,   299,   300,
     301,   302,   303,   304,   305,   306,   307,   308,   309,   309,
     310,   311,   311,   312,   313,   314,   315,   316,   317,   317,
     318,   319,   320,   321,   322,   323,   323,   324,   325,   326,
     327,   328,   329,   330,   331,   332,   333,   334,   335,   336,
     337,   338,   339,   340,   341,   342,   342,   343,   344,   345,
     346,   347,   348,   349,   350,   351,   352,   353,   354,   355,
     356,   357,   358,   359,   360,   361,   361,   362,   362,   363
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  parser::yyr2_[] =
  {
         0,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     3,     4,     3,
       3,     3,     3,     3,     2,     3,     1,     3,     4,     2,
       2,     3,     1,     3,     4,     2,     2,     3,     1,     3,
       4,     2,     2,     3,     1,     3,     4,     2,     3,     4,
       3,     4,     3,     4,     4,     3,     1,     1,     1,     3,
       3,     3,     3,     3,     2,     2,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     1,     3,     4,
       7,     3,     4,     2,     1,     4,     4,     2,     1,     7,
       3,     1,     1,     1,     1,     1,     0,     5,     0,     8,
       0,     8,     0,    10,     0,     8,     2,     2,     1,     1,
       4,     2,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     2,     2,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     5,     1,     4,     4,     4,     2,     1,
       9,     6,     5,     7,     7,     2,     4,     3,     5,     3,
       1,     2,     2,     2,     1,     1,     1,     4,     3,     6,
       3,     1,     5,     3,     3,     4,     2,     2,     3,     1,
       1,     2,     5,     3,     1,     1,     2,     5,     3,     1,
       1,     2,     5,     3,     1,     1,     1,     2,     5,     3,
       6,     3,     1,     1,     1,     1,     1,     1,     1,     1,
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
       6,     4,     2,     1,     5,     7,     9,     0,     3,     3,
       2,     5,     5,     6,     3,     7,     8,     5,     5,     6,
       3,     7,     8,     5,     6,     3,     1,     1,     1,     1,
       1,     3,     4,     6,     1,     2,     1,     1,     1,     1,
       1,     0,     0,     5,     2,     5,     3,     6,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     3,     1,     3,
       6,     1,     1,     1,     1,     3,     1,     3,     6,     3,
       3,     3,     1,     3,     3,     3,     3,     1,     1,     1,
       3,     3,     3,     3,     3,     3,     1,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     1,
       3,     3,     3,     3,     5,     3,     3,     3,     3,     1,
       3,     3,     3,     1,     1,     1,     1,     1,     3,     1,
       1,     1,     1,     3,     3,     3,     3,     1,     1,     3,
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
  "ASINH", "ACOSH", "ATANH", "SQRT", "';'", "'('", "')'", "'#'", "':'",
  "'['", "']'", "'''", "'.'", "'\\\\'", "$accept", "statement_list",
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
  "o_kalman_algo", "o_kalman_tol", "o_model_comparison_approximation",
  "o_print", "o_noprint", "o_xls_sheet", "o_xls_range",
  "o_filter_step_ahead", "o_constant", "o_noconstant", "o_mh_recover",
  "o_planner_discount", "o_bvar_prior_tau", "o_bvar_prior_decay",
  "o_bvar_prior_lambda", "o_bvar_prior_mu", "o_bvar_prior_omega",
  "o_bvar_prior_flat", "o_bvar_prior_train", "o_bvar_replic", "range",
  "vec_int_elem", "vec_int_1", "vec_int", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const parser::rhs_number_type
  parser::yyrhs_[] =
  {
       168,     0,    -1,   169,    -1,   168,   169,    -1,   170,    -1,
     181,    -1,   182,    -1,   183,    -1,   197,    -1,   187,    -1,
     189,    -1,   192,    -1,   184,    -1,   208,    -1,   209,    -1,
     214,    -1,   217,    -1,   220,    -1,   223,    -1,   226,    -1,
     246,    -1,   249,    -1,   252,    -1,   232,    -1,   241,    -1,
     238,    -1,   255,    -1,   256,    -1,   259,    -1,   171,    -1,
     172,    -1,   260,    -1,   262,    -1,   263,    -1,   264,    -1,
     268,    -1,   269,    -1,   270,    -1,   271,    -1,   277,    -1,
     280,    -1,   286,    -1,   289,    -1,   176,    -1,   173,    -1,
     174,    -1,   175,    -1,    28,    51,   157,    -1,    28,    51,
      51,   157,    -1,   108,   229,   157,    -1,   128,   177,   157,
      -1,   129,   178,   157,    -1,   130,   179,   157,    -1,    96,
     180,   157,    -1,   177,    80,    -1,   177,   134,    80,    -1,
      80,    -1,   177,    80,   123,    -1,   177,   134,    80,   123,
      -1,    80,   123,    -1,   178,    80,    -1,   178,   134,    80,
      -1,    80,    -1,   178,    80,   123,    -1,   178,   134,    80,
     123,    -1,    80,   123,    -1,   179,    80,    -1,   179,   134,
      80,    -1,    80,    -1,   179,    80,   123,    -1,   179,   134,
      80,   123,    -1,    80,   123,    -1,   180,    80,    -1,   180,
     134,    80,    -1,    80,    -1,   180,    80,   123,    -1,   180,
     134,    80,   123,    -1,    80,   123,    -1,    97,    51,   157,
      -1,    97,    33,    51,   157,    -1,    24,    42,   157,    -1,
      24,    33,    42,   157,    -1,    63,    42,   157,    -1,    63,
      33,    42,   157,    -1,    80,    33,   185,   157,    -1,   158,
     185,   159,    -1,    80,    -1,    42,    -1,    51,    -1,   185,
     136,   185,    -1,   185,   135,   185,    -1,   185,   137,   185,
      -1,   185,   138,   185,    -1,   185,   140,   185,    -1,   135,
     185,    -1,   136,   185,    -1,   141,   158,   185,   159,    -1,
     142,   158,   185,   159,    -1,   143,   158,   185,   159,    -1,
     144,   158,   185,   159,    -1,   145,   158,   185,   159,    -1,
     146,   158,   185,   159,    -1,   147,   158,   185,   159,    -1,
     148,   158,   185,   159,    -1,   149,   158,   185,   159,    -1,
     156,   158,   185,   159,    -1,    80,   158,   186,   159,    -1,
     185,    -1,   186,   134,   185,    -1,    50,   157,   190,    31,
      -1,    50,   158,   188,   159,   157,   190,    31,    -1,    38,
      33,    80,    -1,    32,   157,   190,    31,    -1,   190,   191,
      -1,   191,    -1,    80,    33,   185,   157,    -1,    47,   157,
     193,    31,    -1,   193,   194,    -1,   194,    -1,    80,   158,
     230,   159,    33,   185,   157,    -1,   195,   134,   196,    -1,
     196,    -1,    57,    -1,    45,    -1,   305,    -1,   306,    -1,
      -1,    74,   157,   198,   203,    31,    -1,    -1,    74,   158,
     293,   159,   157,   199,   203,    31,    -1,    -1,    74,   158,
     126,   159,   157,   200,   203,    31,    -1,    -1,    74,   158,
     116,   134,   195,   159,   201,   157,   203,    31,    -1,    -1,
      74,   158,   116,   159,   202,   157,   203,    31,    -1,   203,
     204,    -1,   203,   206,    -1,   204,    -1,   206,    -1,   205,
      33,   205,   157,    -1,   205,   157,    -1,   158,   205,   159,
      -1,   207,    -1,    42,    -1,    51,    -1,   205,   136,   205,
      -1,   205,   135,   205,    -1,   205,   137,   205,    -1,   205,
     138,   205,    -1,   205,   140,   205,    -1,   135,   205,    -1,
     136,   205,    -1,   141,   158,   205,   159,    -1,   142,   158,
     205,   159,    -1,   143,   158,   205,   159,    -1,   144,   158,
     205,   159,    -1,   145,   158,   205,   159,    -1,   146,   158,
     205,   159,    -1,   147,   158,   205,   159,    -1,   148,   158,
     205,   159,    -1,   149,   158,   205,   159,    -1,   156,   158,
     205,   159,    -1,   160,    80,    33,   205,   157,    -1,    80,
      -1,    80,   158,   230,   159,    -1,   109,   157,   210,    31,
      -1,    76,   157,   210,    31,    -1,   210,   211,    -1,   211,
      -1,   128,    80,   157,    97,   212,   157,   127,   213,   157,
      -1,   128,    80,   157,   117,   185,   157,    -1,   128,    80,
      33,   185,   157,    -1,   128,    80,   134,    80,    33,   185,
     157,    -1,    22,    80,   134,    80,    33,   185,   157,    -1,
     212,    51,    -1,   212,    51,   161,    51,    -1,   212,   134,
      51,    -1,   212,   134,    51,   161,    51,    -1,    51,   161,
      51,    -1,    51,    -1,   213,   231,    -1,   213,   230,    -1,
     213,    80,    -1,   231,    -1,   230,    -1,    80,    -1,   213,
     158,   185,   159,    -1,   158,   185,   159,    -1,   110,    33,
     162,   215,   163,   157,    -1,   215,   157,   216,    -1,   216,
      -1,   216,   134,   158,   185,   159,    -1,   216,   134,    42,
      -1,   216,   134,    51,    -1,   216,   158,   185,   159,    -1,
     216,    42,    -1,   216,    51,    -1,   158,   185,   159,    -1,
      42,    -1,    51,    -1,   118,   157,    -1,   118,   158,   218,
     159,   157,    -1,   218,   134,   219,    -1,   219,    -1,   291,
      -1,    19,   157,    -1,    19,   158,   221,   159,   157,    -1,
     221,   134,   222,    -1,   222,    -1,   291,    -1,   111,   157,
      -1,   111,   158,   224,   159,   157,    -1,   224,   134,   225,
      -1,   225,    -1,   304,    -1,   310,    -1,   119,   157,    -1,
     119,   158,   227,   159,   157,    -1,   119,   229,   157,    -1,
     119,   158,   227,   159,   229,   157,    -1,   227,   134,   228,
      -1,   228,    -1,   290,    -1,   291,    -1,   292,    -1,   293,
      -1,   294,    -1,   295,    -1,   296,    -1,   297,    -1,   298,
      -1,   299,    -1,   300,    -1,   317,    -1,   301,    -1,   339,
      -1,   302,    -1,   303,    -1,   304,    -1,   305,    -1,   307,
      -1,   308,    -1,   309,    -1,   343,    -1,   344,    -1,   229,
      80,    -1,   229,    80,    33,    80,    -1,   229,   134,    80,
      -1,   229,   134,    80,    33,    80,    -1,    80,    -1,    80,
      33,    80,    -1,   136,    51,    -1,   135,    51,    -1,    51,
      -1,   136,    42,    -1,   135,    42,    -1,    42,    -1,    35,
     157,   233,    31,    -1,   233,   234,    -1,   234,    -1,   235,
     134,   236,   157,    -1,   117,    80,    -1,    80,    -1,    22,
      80,   134,    80,    -1,   244,   134,   237,    -1,   245,   134,
     244,   134,   237,    -1,   245,   134,   245,   134,   245,   134,
     244,   134,   237,    -1,   245,    -1,   245,   134,   245,   134,
     245,    -1,   245,   134,   245,    -1,   245,   134,   245,   134,
     245,    -1,   245,   134,   245,   134,   245,   134,   245,    -1,
     245,   134,   245,   134,   245,   134,   245,   134,   245,    -1,
      37,   157,   239,    31,    -1,   239,   240,    -1,   240,    -1,
     117,    80,   134,   245,   157,    -1,    22,    80,   134,    80,
     134,   245,   157,    -1,    80,   134,   245,   157,    -1,    36,
     157,   242,    31,    -1,   242,   243,    -1,   243,    -1,   117,
      80,   134,   245,   134,   245,   157,    -1,    22,    80,   134,
      80,   134,   245,   134,   245,   157,    -1,    80,   134,   245,
     134,   245,   157,    -1,     6,    -1,    44,    -1,    89,    -1,
      52,    -1,   124,    -1,    -1,    51,    -1,    42,    -1,    80,
      -1,   135,    51,    -1,   135,    42,    -1,    34,   157,    -1,
      34,   158,   247,   159,   157,    -1,    34,   229,   157,    -1,
      34,   158,   247,   159,   229,   157,    -1,   247,   134,   248,
      -1,   248,    -1,   310,    -1,   311,    -1,   312,    -1,   313,
      -1,   314,    -1,   315,    -1,   316,    -1,   317,    -1,   318,
      -1,   319,    -1,   320,    -1,   321,    -1,   322,    -1,   323,
      -1,   324,    -1,   325,    -1,   326,    -1,   327,    -1,   328,
      -1,   329,    -1,   330,    -1,   331,    -1,   332,    -1,   333,
      -1,   301,    -1,   334,    -1,   335,    -1,   336,    -1,   337,
      -1,   338,    -1,   340,    -1,   341,    -1,   345,    -1,   346,
      -1,   347,    -1,   291,    -1,   348,    -1,   349,    -1,   350,
      -1,   103,   158,   250,   159,   157,    -1,   103,   158,   250,
     159,   229,   157,    -1,   250,   134,   251,    -1,   251,    -1,
     317,    -1,   318,    -1,   327,    -1,   333,    -1,   301,    -1,
     334,    -1,   335,    -1,   336,    -1,   337,    -1,   338,    -1,
     345,    -1,   346,    -1,   347,    -1,   104,   158,   250,   159,
     157,    -1,   104,   158,   250,   159,   229,   157,    -1,   164,
      80,   164,   134,   164,    80,   164,    -1,   164,    80,   164,
     134,   245,    -1,   253,    -1,   254,   134,   253,    -1,   131,
     229,   157,    -1,    90,   157,   257,    31,    -1,   257,   258,
      -1,   258,    -1,    80,   158,   185,   159,   157,    -1,   125,
     229,   157,    -1,    92,   157,   261,    31,    -1,   261,    80,
     185,   157,    -1,   261,    80,   134,    80,   185,   157,    -1,
      80,   185,   157,    -1,    80,   134,    80,   185,   157,    -1,
      95,   229,   157,    -1,    94,   157,    -1,    94,   158,   228,
     159,   157,    -1,    94,   229,   157,    -1,    94,   158,   228,
     159,   229,   157,    -1,    18,   157,   265,    31,    -1,   265,
     266,    -1,   266,    -1,    80,   267,    33,   185,   157,    -1,
      80,   134,    80,   267,    33,   185,   157,    -1,     4,    80,
     158,    51,   159,   267,    33,   185,   157,    -1,    -1,   158,
      51,   159,    -1,   158,    42,   159,    -1,    17,   157,    -1,
      17,   158,    23,   159,   157,    -1,    30,   158,    80,   159,
     157,    -1,    30,   158,    80,   159,   229,   157,    -1,    30,
      80,   157,    -1,    30,   158,    80,   165,    80,   159,   157,
      -1,    30,   158,    80,   165,    80,   159,   229,   157,    -1,
      30,    80,   165,    80,   157,    -1,    29,   158,    80,   159,
     157,    -1,    29,   158,    80,   159,   229,   157,    -1,    29,
      80,   157,    -1,    29,   158,    80,   165,    80,   159,   157,
      -1,    29,   158,    80,   165,    80,   159,   229,   157,    -1,
      29,    80,   165,    80,   157,    -1,    75,   158,   272,   159,
     274,   157,    -1,   272,   134,   273,    -1,   273,    -1,   342,
      -1,   343,    -1,   344,    -1,   275,    -1,   274,   134,   275,
      -1,   275,   158,   245,   159,    -1,   274,   134,   275,   158,
     245,   159,    -1,   276,    -1,   275,   276,    -1,    80,    -1,
     166,    -1,   137,    -1,   161,    -1,   165,    -1,    -1,    -1,
      98,   278,   205,   279,   157,    -1,   121,   157,    -1,   121,
     158,   281,   159,   157,    -1,   121,   229,   157,    -1,   121,
     158,   281,   159,   229,   157,    -1,   281,   134,   282,    -1,
     282,    -1,   228,    -1,   351,    -1,   352,    -1,   353,    -1,
     354,    -1,   355,    -1,   356,    -1,   357,    -1,   358,    -1,
     283,    -1,   310,    -1,   345,    -1,   346,    -1,   312,    -1,
     314,    -1,   311,    -1,   313,    -1,   284,   134,   285,    -1,
     284,    -1,     7,    51,   157,    -1,     7,   158,   285,   159,
      51,   157,    -1,   284,    -1,   335,    -1,   318,    -1,   359,
      -1,   287,   134,   288,    -1,   287,    -1,     8,    51,   157,
      -1,     8,   158,   288,   159,    51,   157,    -1,    26,    33,
      51,    -1,   115,    33,    51,    -1,   112,    33,    51,    -1,
      60,    -1,    93,    33,    51,    -1,   107,    33,    51,    -1,
      27,    33,    51,    -1,     3,    33,    51,    -1,    83,    -1,
      85,    -1,    87,    -1,    53,    33,    51,    -1,    48,    33,
      51,    -1,    49,    33,    51,    -1,    97,    33,    51,    -1,
      24,    33,    42,    -1,    63,    33,    42,    -1,   111,    -1,
     113,    33,    51,    -1,   105,    33,    51,    -1,   105,    33,
      42,    -1,    25,    33,    80,    -1,    81,    33,   363,    -1,
      81,    33,    51,    -1,    41,    33,    51,    -1,    99,    33,
      51,    -1,   100,    33,    51,    -1,    58,    33,    51,    -1,
      59,    33,    51,    -1,    86,    -1,    46,    -1,    20,    33,
      42,    -1,    69,    33,    51,    -1,    64,    33,    42,    -1,
      66,    33,    42,    -1,    91,    33,   158,   254,   159,    -1,
      65,    33,    42,    -1,    65,    33,    51,    -1,    73,    33,
      80,    -1,    72,    33,    51,    -1,    71,    -1,   102,    33,
      42,    -1,    67,    33,    51,    -1,    68,    33,    51,    -1,
      61,    -1,    62,    -1,    84,    -1,     5,    -1,   120,    -1,
      43,    33,    51,    -1,   114,    -1,    79,    -1,    40,    -1,
     106,    -1,    54,    33,    51,    -1,    55,    33,    51,    -1,
      77,    33,    56,    -1,    77,    33,    78,    -1,   101,    -1,
      88,    -1,   132,    33,    80,    -1,   133,    33,   360,    -1,
      39,    33,   363,    -1,    21,    -1,    82,    -1,    70,    -1,
     122,    33,    42,    -1,    14,    33,   231,    -1,     9,    33,
      42,    -1,    11,    33,   231,    -1,    12,    33,    42,    -1,
      13,    33,    51,    -1,    10,    -1,    15,    33,    51,    -1,
      16,    33,    51,    -1,    80,   161,    80,    -1,    51,    -1,
      51,   161,    51,    -1,   162,   361,    -1,   362,   361,    -1,
     362,   163,    -1
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
      80,    82,    84,    86,    88,    90,    92,    94,    98,   103,
     107,   111,   115,   119,   123,   126,   130,   132,   136,   141,
     144,   147,   151,   153,   157,   162,   165,   168,   172,   174,
     178,   183,   186,   189,   193,   195,   199,   204,   207,   211,
     216,   220,   225,   229,   234,   239,   243,   245,   247,   249,
     253,   257,   261,   265,   269,   272,   275,   280,   285,   290,
     295,   300,   305,   310,   315,   320,   325,   330,   332,   336,
     341,   349,   353,   358,   361,   363,   368,   373,   376,   378,
     386,   390,   392,   394,   396,   398,   400,   401,   407,   408,
     417,   418,   427,   428,   439,   440,   449,   452,   455,   457,
     459,   464,   467,   471,   473,   475,   477,   481,   485,   489,
     493,   497,   500,   503,   508,   513,   518,   523,   528,   533,
     538,   543,   548,   553,   559,   561,   566,   571,   576,   579,
     581,   591,   598,   604,   612,   620,   623,   628,   632,   638,
     642,   644,   647,   650,   653,   655,   657,   659,   664,   668,
     675,   679,   681,   687,   691,   695,   700,   703,   706,   710,
     712,   714,   717,   723,   727,   729,   731,   734,   740,   744,
     746,   748,   751,   757,   761,   763,   765,   767,   770,   776,
     780,   787,   791,   793,   795,   797,   799,   801,   803,   805,
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
    1297,  1304,  1309,  1312,  1314,  1320,  1328,  1338,  1339,  1343,
    1347,  1350,  1356,  1362,  1369,  1373,  1381,  1390,  1396,  1402,
    1409,  1413,  1421,  1430,  1436,  1443,  1447,  1449,  1451,  1453,
    1455,  1457,  1461,  1466,  1473,  1475,  1478,  1480,  1482,  1484,
    1486,  1488,  1489,  1490,  1496,  1499,  1505,  1509,  1516,  1520,
    1522,  1524,  1526,  1528,  1530,  1532,  1534,  1536,  1538,  1540,
    1542,  1544,  1546,  1548,  1550,  1552,  1554,  1556,  1560,  1562,
    1566,  1573,  1575,  1577,  1579,  1581,  1585,  1587,  1591,  1598,
    1602,  1606,  1610,  1612,  1616,  1620,  1624,  1628,  1630,  1632,
    1634,  1638,  1642,  1646,  1650,  1654,  1658,  1660,  1664,  1668,
    1672,  1676,  1680,  1684,  1688,  1692,  1696,  1700,  1704,  1706,
    1708,  1712,  1716,  1720,  1724,  1730,  1734,  1738,  1742,  1746,
    1748,  1752,  1756,  1760,  1762,  1764,  1766,  1768,  1770,  1774,
    1776,  1778,  1780,  1782,  1786,  1790,  1794,  1798,  1800,  1802,
    1806,  1810,  1814,  1816,  1818,  1820,  1824,  1828,  1832,  1836,
    1840,  1844,  1846,  1850,  1854,  1858,  1860,  1864,  1867,  1870
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    88,    88,    89,    93,    94,    95,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   118,
     119,   120,   121,   122,   123,   124,   125,   126,   127,   128,
     129,   130,   131,   136,   137,   138,   139,   143,   144,   147,
     150,   154,   158,   162,   166,   168,   170,   172,   174,   176,
     181,   183,   185,   187,   189,   191,   196,   198,   200,   202,
     204,   206,   211,   213,   215,   217,   219,   221,   226,   230,
     237,   241,   248,   252,   260,   265,   267,   269,   271,   273,
     275,   277,   279,   281,   283,   285,   287,   289,   291,   293,
     295,   297,   299,   301,   303,   305,   307,   312,   314,   318,
     320,   325,   329,   334,   335,   339,   344,   349,   350,   354,
     358,   359,   363,   364,   365,   366,   370,   370,   371,   371,
     373,   373,   375,   375,   377,   377,   382,   383,   384,   385,
     389,   391,   396,   397,   398,   400,   402,   404,   406,   408,
     410,   412,   414,   416,   418,   420,   422,   424,   426,   428,
     430,   432,   434,   438,   442,   444,   449,   453,   457,   458,
     462,   464,   466,   468,   470,   475,   477,   479,   481,   483,
     485,   491,   493,   495,   497,   499,   501,   503,   505,   510,
     515,   517,   522,   524,   526,   528,   530,   532,   534,   536,
     538,   543,   547,   551,   552,   555,   559,   561,   565,   566,
     569,   573,   575,   579,   580,   583,   584,   588,   590,   592,
     594,   598,   599,   602,   603,   604,   605,   606,   607,   608,
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
    1018,  1021,  1025,  1026,  1029,  1030,  1031,  1034,  1035,  1036,
    1039,  1040,  1043,  1044,  1045,  1046,  1047,  1048,  1050,  1051,
    1052,  1053,  1054,  1055,  1057,  1061,  1062,  1065,  1066,  1067,
    1070,  1071,  1072,  1073,  1076,  1077,  1080,  1081,  1082,  1083,
    1084,  1087,  1087,  1087,  1090,  1092,  1094,  1096,  1101,  1102,
    1105,  1106,  1109,  1110,  1111,  1112,  1113,  1114,  1115,  1118,
    1119,  1120,  1121,  1122,  1123,  1124,  1125,  1128,  1129,  1132,
    1134,  1138,  1139,  1140,  1141,  1144,  1145,  1148,  1150,  1154,
    1155,  1156,  1157,  1158,  1159,  1160,  1161,  1162,  1163,  1164,
    1165,  1166,  1167,  1168,  1169,  1170,  1171,  1172,  1173,  1174,
    1176,  1177,  1178,  1180,  1181,  1182,  1183,  1184,  1185,  1186,
    1187,  1188,  1189,  1190,  1191,  1192,  1193,  1194,  1195,  1196,
    1197,  1198,  1199,  1200,  1201,  1202,  1203,  1204,  1205,  1206,
    1207,  1208,  1209,  1210,  1211,  1213,  1215,  1218,  1219,  1220,
    1221,  1222,  1223,  1224,  1225,  1226,  1228,  1229,  1230,  1231,
    1232,  1233,  1234,  1235,  1237,  1245,  1246,  1250,  1251,  1260
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
       2,     2,     2,     2,     2,   160,     2,     2,     2,   164,
     158,   159,     2,     2,     2,     2,   165,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   161,   157,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   162,   166,   163,     2,     2,     2,     2,     2,     2,
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
     155,   156
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1535;
  const int parser::yynnts_ = 197;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 156;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 167;

  const unsigned int parser::yyuser_token_number_max_ = 411;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1262 "DynareBison.yy"


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

