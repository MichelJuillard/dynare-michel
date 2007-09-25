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

  case 451:
#line 1135 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 452:
#line 1137 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 459:
#line 1151 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 460:
#line 1153 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 461:
#line 1156 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 462:
#line 1157 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 463:
#line 1158 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 464:
#line 1159 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 465:
#line 1160 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 466:
#line 1161 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 467:
#line 1162 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 468:
#line 1163 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 469:
#line 1164 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 470:
#line 1165 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 471:
#line 1166 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 472:
#line 1167 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 473:
#line 1168 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 474:
#line 1169 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 475:
#line 1170 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 476:
#line 1171 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 477:
#line 1172 "DynareBison.yy"
    {driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 478:
#line 1173 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 479:
#line 1174 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 480:
#line 1175 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 481:
#line 1176 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 482:
#line 1178 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 483:
#line 1179 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 484:
#line 1180 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 485:
#line 1182 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 486:
#line 1183 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 487:
#line 1184 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 488:
#line 1185 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 489:
#line 1186 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 490:
#line 1187 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 491:
#line 1188 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 492:
#line 1189 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 493:
#line 1190 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 494:
#line 1191 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1192 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 497:
#line 1194 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 498:
#line 1195 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 499:
#line 1196 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 500:
#line 1197 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 501:
#line 1198 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 502:
#line 1199 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 503:
#line 1200 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 504:
#line 1201 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 505:
#line 1202 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 506:
#line 1203 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 507:
#line 1204 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 508:
#line 1205 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 509:
#line 1206 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 510:
#line 1207 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 511:
#line 1208 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 512:
#line 1209 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 513:
#line 1210 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 514:
#line 1211 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 515:
#line 1212 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 516:
#line 1213 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 517:
#line 1216 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 518:
#line 1218 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 519:
#line 1220 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 520:
#line 1221 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 521:
#line 1222 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 522:
#line 1223 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 523:
#line 1224 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 524:
#line 1225 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 525:
#line 1226 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 526:
#line 1227 "DynareBison.yy"
    {driver.option_num("mh_recover", "1");;}
    break;

  case 527:
#line 1228 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 528:
#line 1230 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1231 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1232 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1233 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1234 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1235 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 534:
#line 1236 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 535:
#line 1237 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 536:
#line 1240 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 538:
#line 1249 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 539:
#line 1252 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 540:
#line 1254 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 541:
#line 1262 "DynareBison.yy"
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
  const short int parser::yypact_ninf_ = -922;
  const short int
  parser::yypact_[] =
  {
       779,    10,    41,   174,   -73,   387,    93,    58,    -5,    49,
     -47,   167,    -9,    15,    32,    39,   488,   504,   521,    73,
      80,   104,   107,   135,   190,   219,   308,   438,  -922,   261,
     263,   219,   277,   403,   541,   560,   213,   218,   219,   377,
     388,   389,   219,    83,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,   316,   923,
     323,   902,  -922,   538,   120,  -922,   461,   543,   444,    70,
     383,   528,   485,   559,   563,   571,  -922,   661,   395,   203,
     286,   287,   573,   563,   632,   633,   523,  -922,   299,   352,
     105,   814,   593,   603,  -922,  1080,   440,   446,   569,   448,
     642,   552,  1002,   716,   716,   453,   105,   551,  -922,     9,
    -922,   461,  -922,  1080,   454,  -922,   964,   455,   456,   614,
     499,   625,   520,   630,   534,   537,  -922,  -922,  -922,   705,
    -922,   706,   725,   731,   733,   735,  -922,   737,   744,   746,
    -922,   750,   751,   752,   757,  -922,   657,   640,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,   772,   773,   784,  -922,   691,
     668,  -922,  -922,  -922,   669,   758,   -54,   166,  -922,   798,
     -20,  -922,  -922,   682,  -922,   683,  -922,  -922,   761,   -62,
    -922,   764,   -27,   812,    54,  -922,   766,  -922,   817,  -922,
    -922,   818,   819,   824,   825,   827,  -922,  -922,   828,   830,
     833,   834,   835,   837,  -922,  -922,   839,   845,  -922,  -922,
    -922,   846,   847,  -922,  -922,   162,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,   848,   767,  -922,   804,
    -922,   805,   246,  -922,   762,   806,   765,   811,   272,  -922,
     815,   769,   821,   280,  -922,   734,    56,  -922,    67,   872,
     760,   749,  -922,   480,  -922,   171,   770,   771,   887,  -922,
    -922,   267,  -922,  -922,  -922,  -922,   841,   844,    48,  -922,
    -922,  -922,   768,   814,   814,   781,   782,   783,   788,   789,
     793,   794,   795,   796,   807,   814,   891,   808,   211,  -922,
     851,   312,   892,   895,   909,   935,   936,   938,  -922,  -922,
    -922,   944,   946,   947,  -922,   948,  -922,   949,   956,   849,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,   880,   926,  -922,   854,
    -922,  -922,  -922,   856,  1002,  1002,   857,   860,   861,   862,
     863,   874,   875,   882,   884,   885,  1002,   697,  -922,   283,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,   284,  -922,   100,    27,   296,  -922,
    -922,  -922,   297,  -922,  -922,   305,  -922,  -922,   997,  -922,
     307,  -922,  -922,  -922,  -922,  -922,   922,   966,  -922,  -922,
     939,   978,  -922,  -922,   940,   979,  -922,  -922,  1022,    16,
    1024,  1016,    16,  1021,   988,  1023,     2,  1027,  1029,   993,
    1001,   923,  1033,  1034,  1049,  1042,   902,  1043,   941,   937,
    1019,   490,  1063,  -922,  -922,  1046,   461,   955,  -922,  -922,
     957,   -13,  1045,   970,    -7,  1051,   814,  -922,  -922,  -922,
     951,  1083,  1084,  1085,  1088,  1090,  1100,   508,  1114,  1108,
    1110,  1111,  1118,  1091,  1006,  1128,   661,    69,  1092,  1141,
    1041,  -922,  -922,  -922,   230,  1044,   237,  1048,  -922,  -922,
    1050,   237,  1054,  -922,  -922,    20,  -922,  -922,  -922,  1103,
    1032,  -922,  1120,    74,  -922,     5,  -922,   747,  -922,  1037,
    1055,   234,   352,   196,  1056,    33,  -922,  -922,   814,  1039,
     547,   814,   814,   814,   814,   814,   814,   814,   814,   814,
     814,   433,   814,   814,   814,   814,   814,  -922,   814,  -922,
    -922,  1126,  1203,  -922,   974,  1129,  1176,  1170,  1193,  1196,
    1219,  1222,  1245,   514,  1248,  1271,  1300,   159,  -922,  1202,
    -922,    20,  1208,   568,  1002,  1002,  1002,  1002,  1002,  1002,
    1002,  1002,  1002,  1002,   684,  1002,  1002,  1002,  1002,  1002,
    1201,   716,   220,   266,  -922,  -922,  -922,   814,   367,    26,
       9,  1204,   461,  1215,  1080,   270,  1317,   964,   275,  -922,
    1256,  -922,  1259,  -922,  1260,  -922,  -922,  1339,  1343,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  1341,     3,  -922,
    -922,  -922,  -922,  1242,  -922,  -922,  1249,  -922,  -922,  -922,
    -922,  1250,  -922,  1354,  1255,  1257,  1268,   814,  -922,  -922,
    -922,  -922,  -922,   555,  1270,  -922,  -922,   610,  1277,  1209,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  1266,  -922,  -922,  -922,   615,
    -922,  1358,  1360,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,   524,  1294,  1319,  1321,  1376,  1323,   237,  1378,  1302,
     237,  -922,  1409,  1410,  1303,  -922,   563,  1431,  -922,  -922,
    -922,  1002,  -922,  -922,  -922,  1432,   311,  -922,  -922,  -922,
    1309,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,   170,   441,  -922,  1387,   814,  1388,    11,   838,   313,
     901,   952,   965,  1017,  1061,  1067,  1073,  1079,  1087,  1093,
    -922,   547,   547,  1039,  1039,  1330,  1099,   814,  -922,  1391,
    1217,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,   616,  -922,  1313,  1105,  1113,
    1119,  1125,  1131,  1139,  1145,  1151,  1157,  1165,  -922,   568,
     568,  1208,  1208,  1330,  -922,  -922,  -922,   617,  -922,   623,
    1171,    27,  1316,  -922,  -922,    30,   814,  -922,  -922,  -922,
    -922,  -922,  -922,   631,  -922,  -922,  -922,   644,  -922,  -922,
    -922,  -922,  -922,  1314,  -922,  -922,  -922,  1394,  -922,  -922,
    1318,  1443,  -922,  -922,  1227,  -922,   276,  -922,   278,  -922,
    1398,  -922,   405,  -922,  -922,  -922,  -922,  -922,  -922,   237,
     230,  1345,   237,  1346,  1347,  -922,  1325,  -922,  -922,  1450,
     371,  1002,  1233,  1442,   747,  -922,   480,   480,   480,   196,
    -922,   237,  -922,  1452,  1240,  1453,  1436,   814,   814,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    1331,  1251,   814,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,    26,
    -922,  -922,  -922,   814,  1177,  -922,  -922,  1438,  -922,  1255,
     814,  -922,  -922,   655,  -922,   666,  1326,  1266,  -922,  -922,
    1357,  1359,  1361,   237,  1335,   237,   237,  -922,   814,  -922,
    1258,  -922,  -922,  -922,  1337,   193,   238,   264,   483,  1338,
     814,  -922,   814,  1340,   144,  1264,   838,  -922,  -922,  1274,
    1183,  -922,  -922,  1463,  1282,  -922,  -922,  1364,  -922,   237,
     237,   237,  1365,  -922,  1348,  1349,  1288,  -922,   480,  -922,
    -922,  -922,   237,  -922,  1297,  1306,  1449,  1342,  1451,  1377,
    -922,  -922,  -922,   814,  -922,   123,  1373,  -922,  1374,   237,
    -922,  -922,  -922,   516,  1350,  -922,  -922,  -922,  1459,  1351,
     187,  1312,  1433,  -922,   237,   117,  1362,  -922,  -922,  -922,
    1460,  -922,   530,   535,   814,   233,  -922,  -922,  -922,  1352,
    1380,  1381,  -922,  -922,  1191,  -922,  -922,   814,  -922,  -922,
    -922,   237,   237,  -922,  1197,  1383,  -922,  -922,   237,  -922
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
       0,    62,     0,    68,     0,     0,     1,     3,   451,     0,
     533,     0,     0,     0,     0,     0,   524,     0,     0,     0,
     525,     0,     0,     0,     0,   439,   450,     0,   440,   445,
     443,   446,   444,   441,   442,   447,   448,   432,   433,   434,
     435,   436,   437,   438,   459,     0,     0,     0,   453,   458,
       0,   455,   454,   456,     0,     0,   387,     0,   383,     0,
       0,   209,   210,     0,    80,     0,    47,   400,     0,     0,
     394,     0,     0,     0,     0,   114,     0,   508,     0,   513,
     491,     0,     0,     0,     0,     0,   505,   506,     0,     0,
       0,     0,     0,     0,   526,   501,     0,     0,   512,   507,
     490,     0,     0,   511,   509,     0,   302,   338,   327,   303,
     304,   305,   306,   307,   308,   309,   310,   311,   312,   313,
     314,   315,   316,   317,   318,   319,   320,   321,   322,   323,
     324,   325,   326,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   339,   340,   341,   246,     0,   299,     0,
     263,     0,     0,   260,     0,     0,     0,     0,     0,   282,
       0,     0,     0,     0,   276,     0,     0,   118,     0,     0,
       0,     0,    82,     0,   464,     0,     0,     0,     0,   520,
     519,     0,   406,   407,   408,   409,     0,     0,     0,   169,
      87,    88,    86,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   368,
       0,     0,     0,     0,     0,     0,     0,     0,   469,   470,
     471,     0,     0,     0,   514,     0,   478,     0,     0,     0,
     223,   224,   225,   226,   227,   228,   229,   230,   231,   232,
     233,   235,   237,   238,   239,   240,   241,   242,   243,   234,
     236,   244,   245,   379,   376,    77,    72,     0,    53,     0,
      78,   144,   145,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   422,   143,     0,
     345,   350,   346,   347,   348,   349,   351,   352,   353,   354,
     355,   356,   357,   358,     0,    49,     0,     0,     0,   214,
     215,   216,     0,   204,   205,     0,   222,   219,     0,   430,
       0,   429,   431,   426,   370,    59,    54,     0,    50,    65,
      60,     0,    51,    71,    66,     0,    52,   365,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   381,   382,     0,     0,     0,    81,    48,
       0,     0,     0,     0,     0,     0,     0,   112,   113,   251,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   248,
       0,   262,   258,   259,   291,     0,   291,     0,   280,   281,
       0,   291,     0,   274,   275,     0,   116,   117,   109,     0,
       0,    83,     0,     0,   138,     0,   139,     0,   134,     0,
       0,     0,     0,     0,     0,     0,   167,   168,     0,    94,
      95,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    84,     0,   366,
     367,     0,     0,   371,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    75,    73,
      79,     0,   151,   152,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   166,   199,   200,     0,     0,   191,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    57,
      55,    63,    61,    69,    67,   529,   257,     0,     0,   530,
     531,   532,   528,   534,   482,   485,   484,     0,     0,   483,
     486,   487,   521,     0,   522,   449,     0,   535,   492,   510,
     457,     0,   391,     0,   387,     0,     0,     0,   462,   208,
     207,   403,   398,     0,     0,   397,   392,     0,     0,     0,
     523,   472,   515,   516,   488,   489,   494,   497,   498,   495,
     503,   504,   493,   500,   499,     0,   502,   301,   298,     0,
     247,     0,     0,   286,   293,   287,   292,   289,   294,   288,
     290,     0,     0,     0,   268,     0,     0,   291,     0,     0,
     291,   254,     0,     0,     0,   111,     0,     0,   127,   136,
     137,     0,   141,   123,   122,     0,     0,   121,   124,   125,
       0,   130,   128,   517,   518,   405,   416,   418,   419,   420,
     417,     0,   410,   414,     0,     0,     0,     0,   107,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      85,    90,    89,    91,    92,    93,     0,     0,   374,     0,
       0,   468,   476,   461,   467,   473,   474,   465,   475,   481,
     480,   466,   463,   479,   378,     0,    76,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   142,   147,
     146,   148,   149,   150,   423,   344,   342,     0,   359,     0,
       0,     0,     0,   196,   197,     0,     0,   213,   212,   203,
     202,   221,   218,     0,   527,   428,   425,     0,    58,    64,
      70,   256,   255,   537,   539,   541,   540,     0,   452,   460,
       0,     0,   389,   388,     0,   399,     0,   393,     0,   115,
       0,   363,     0,   300,   249,   264,   296,   295,   261,   291,
     291,     0,   291,     0,     0,   279,     0,   253,   252,     0,
       0,     0,     0,     0,     0,   132,     0,     0,     0,     0,
     404,   291,   415,     0,     0,     0,     0,     0,     0,   106,
      96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
       0,     0,     0,   372,   380,   165,   153,   154,   155,   156,
     157,   158,   159,   160,   161,   162,   343,   360,   198,   190,
     189,   193,   194,     0,     0,   220,   427,     0,   536,   387,
       0,   384,   401,     0,   395,     0,     0,     0,   496,   265,
       0,     0,     0,   291,     0,   291,   291,   277,     0,   110,
       0,   140,   477,   120,     0,     0,     0,     0,   411,     0,
       0,   172,     0,   180,     0,     0,   108,   369,   375,     0,
       0,   195,   538,     0,     0,   402,   396,     0,   364,   291,
     291,   291,     0,   285,     0,     0,     0,   163,     0,   135,
     131,   129,   291,   412,     0,     0,     0,   175,     0,     0,
     171,   373,   192,     0,   385,   291,   270,   266,   269,   291,
     283,   278,   119,     0,     0,   174,   173,   179,     0,   177,
       0,     0,     0,   362,   291,     0,     0,   133,   413,   176,
       0,   186,     0,     0,     0,     0,   185,   184,   386,     0,
     271,     0,   284,   178,     0,   183,   170,     0,   182,   181,
     361,   291,   291,   188,     0,   272,   267,   187,   291,   273
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -922,  -922,  1475,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -333,  -922,
    -922,  -922,  -922,  -109,  -222,  -922,  -922,  1205,  -922,   606,
    -922,  -922,  -922,  -922,  -922,  -922,  -823,  -537,  -129,  -534,
    -922,  -922,  -922,  1386,  -264,  -922,  -922,  -922,  -922,   672,
    -922,  -922,   873,  -922,  -922,  1018,  -922,  -922,   876,  -922,
    -922,  -100,   -24,  -593,  -477,  -922,  -922,  1225,  -922,  -922,
    -921,  -922,  -922,  1216,  -922,  -922,  1220,  -867,  -507,  -922,
    -922,   994,  -922,  1397,   893,  -922,   556,  -922,  -922,  -922,
    -922,  1174,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  1328,
    -662,  -922,  -922,  -922,  -922,  -922,   967,  -922,   618,  -732,
    -922,  -922,  -922,  -922,  -922,   879,  -922,   -68,  1047,  -922,
    -922,  1052,  -922,  -922,   -90,  -922,  1422,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,   -98,  -922,  -922,  -123,  -536,  -922,
    -922,  -922,  -922,   -99,   -87,   -55,   -52,   -51,  -922,  -922,
     -92,   -42,  -922,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
     -50,  -922,  -922,  -922,  -922,  -922,   -48,   -45,   -31,   -44,
     -43,   -25,  -922,  -922,  -922,  -922,   -95,   -89,   -88,   -86,
     -21,   -19,   -18,  -922,  -922,  -922,  -922,  -922,  -922,  -922,
    -922,  -922,  -922,  -922,   864,  -922,  1025
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    43,    44,    45,    46,    47,    48,    49,    50,    51,
     150,   152,   154,   129,    52,    53,    54,    55,   356,   789,
      56,   320,    57,   224,   225,    58,   316,   317,   766,   767,
      59,   323,   918,   917,   994,   770,   563,   564,   565,   566,
     428,    60,    61,   338,   339,  1004,  1075,    62,   648,   649,
      63,   452,   453,    64,   210,   211,    65,   448,   449,    66,
     455,   459,   108,   754,   669,    67,   302,   303,   304,   742,
     979,    68,   313,   314,    69,   308,   309,   743,   980,    70,
     255,   256,    71,   429,   430,    72,   891,   892,    73,    74,
     358,   359,    75,    76,   361,    77,    78,    79,   207,   208,
     502,    80,    81,    82,    83,   331,   332,   781,   782,   783,
      84,   132,   640,    85,   460,   461,   175,   176,   177,    86,
     199,   200,    87,   380,   381,   382,   383,   384,   385,   386,
     387,   388,   389,   390,   391,   392,   393,   394,   395,   769,
     396,   397,   398,   178,   179,   180,   181,   182,   264,   265,
     399,   433,   268,   269,   270,   271,   272,   273,   274,   275,
     434,   277,   278,   279,   280,   281,   435,   436,   437,   438,
     439,   440,   400,   288,   289,   333,   401,   402,   183,   184,
     443,   185,   186,   295,   462,   187,   188,   189,   190,   191,
     192,   193,   203,   684,   874,   678,   679
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       126,   127,   518,   427,   318,   672,   212,   135,   259,   258,
     579,   580,   144,   147,   148,   266,   450,   257,   155,   290,
     260,   291,   591,   198,   334,   379,   759,   602,   827,   760,
     335,   768,   881,   981,   167,   431,   431,   744,   761,   746,
     451,   432,   432,   456,   749,   441,   441,   442,   442,   201,
     922,   454,   261,   676,   873,   262,   263,   276,   666,   282,
     202,    88,   283,   285,   286,   267,   785,   105,   853,   645,
     336,   751,   961,   105,   577,   100,   284,   854,   646,   576,
     500,   962,   287,   156,    94,   517,   292,   556,   293,   294,
       1,     2,    90,   995,   996,   997,   518,   511,   558,  1047,
       3,     4,     5,   512,   501,   758,   372,     6,   926,    99,
     104,     7,     8,     9,   506,    10,   411,    11,    12,    13,
      14,   215,   336,   733,   205,   412,    97,   336,   927,   102,
      15,   644,   514,    16,   223,    98,   315,   121,   515,   507,
     635,   636,   637,   638,   702,   639,    17,   223,   109,   105,
     706,   667,   668,   101,   413,   752,   753,    18,    19,    20,
     855,   735,   762,    21,   677,   734,   875,   786,    89,   737,
     205,  1096,   110,    22,   736,    23,   337,    24,    25,    26,
      27,    28,   577,   709,   856,   647,    29,    30,   963,   111,
     787,    31,    32,    33,    34,  1037,   112,   503,  1081,    91,
     206,    35,    36,   738,    37,  1053,   739,   103,    38,   414,
     415,    39,    40,    41,    42,   416,   417,   418,   419,   420,
     421,   422,   423,   424,  1029,   299,   728,   216,   337,   666,
     425,   119,   426,   337,   562,   411,   733,   120,   751,   105,
     903,   740,   599,   906,   412,   788,   206,   105,   790,   791,
     792,   793,   794,   795,   796,   797,   798,   799,   741,   801,
     802,   803,   804,   805,   122,   806,   922,  1071,   299,  1030,
     105,   810,   734,   413,   735,   666,   776,   542,  1038,   734,
     411,   736,   737,   300,   751,   622,   623,  1062,   736,   412,
     773,   357,   123,   105,   305,  1031,   536,   634,   105,   105,
     105,  1039,   310,   548,   919,   567,   411,  1013,   305,   310,
     738,   553,   774,  1085,   850,   412,   824,   738,   413,   739,
     301,   537,  1072,  1073,   106,   107,   300,   920,   414,   415,
     568,    92,    93,   777,   416,   417,   418,   419,   420,   421,
     422,   423,   424,   603,   413,  1074,   105,   124,   125,   425,
     105,   426,   306,   562,   740,   105,   105,   778,   105,   324,
     311,   779,   780,   301,   884,   741,   306,   311,  1072,  1073,
     142,   143,   741,   414,   415,   145,   146,   846,   768,   416,
     417,   418,   419,   420,   421,   422,   423,   424,   128,   307,
    1086,  1087,   604,   982,   425,   984,   426,   312,   562,   414,
     415,   572,   989,   307,   312,   416,   417,   418,   419,   420,
     421,   422,   423,   424,   999,   325,   212,   641,   641,   133,
     425,   134,   426,   848,   562,   326,   573,   862,   198,   328,
     650,   652,   866,   972,   136,   974,   137,   259,   258,   654,
     329,   657,   642,   643,   266,   914,   257,   928,   290,   260,
     291,   223,   924,   330,   201,   651,   653,   149,   759,   759,
     759,   760,   760,   760,   655,   202,   658,  1076,   151,   153,
     915,   130,   929,   158,   941,   296,  1022,   334,  1024,  1025,
     194,   261,  1088,   335,   262,   263,   276,   703,   282,   131,
     707,   283,   285,   286,   267,   828,   829,   830,   831,   832,
     833,   834,   835,   836,   837,   284,   839,   840,   841,   842,
     843,   287,  1046,   729,  1048,   292,   759,   293,   294,   760,
     296,   776,   411,   964,   851,  1054,   296,   450,   406,   297,
     852,   412,   695,   296,   296,   296,   296,   115,  1063,   977,
     217,   696,  1066,   431,    95,    96,   116,  1067,   218,   432,
     717,   451,   298,   441,   861,   442,   819,  1080,   411,   718,
     413,   204,   454,   776,   978,   820,   896,   412,   592,   593,
     594,   595,   871,   596,   297,   897,   209,   872,   777,   466,
     297,   907,   407,  1077,  1095,   213,   908,   297,   297,   297,
     297,  1099,   800,   825,  1005,  1006,   413,   403,  1089,   921,
     470,   214,   778,   404,   226,   408,   779,   780,   219,  1009,
     445,   457,   463,   464,   474,   414,   415,   296,   847,   849,
     777,   416,   417,   418,   419,   420,   421,   422,   423,   424,
    1010,   863,   912,   467,   867,   296,   425,  1014,   426,   222,
     562,  1032,   220,   223,   778,   113,   114,   910,   779,   780,
     221,   414,   415,   315,   471,  1026,   468,   416,   417,   418,
     419,   420,   421,   422,   423,   424,   227,  1034,   475,  1035,
     319,   297,   425,   357,   426,   321,   562,   472,   117,   118,
     322,   196,   166,   360,   594,   595,   167,   596,   518,   297,
     296,   476,   405,   409,   477,   296,   296,   296,   138,   139,
     228,   229,   168,   296,   197,   637,   638,   230,   639,   410,
    1061,   296,   885,   447,   231,   232,   233,   140,   141,   234,
     235,   227,   236,   237,   296,   238,   239,   240,   241,   242,
     243,   244,   245,   246,   247,   296,   196,   465,   478,   479,
     248,  1084,   169,   170,   297,   249,   296,   250,   469,   297,
     297,   297,   251,   473,  1094,   228,   229,   297,   480,   197,
     171,   172,   230,   252,   481,   297,   482,   887,   483,   231,
     484,   363,   893,   944,   956,   253,   209,   485,   297,   486,
     957,   254,   990,   487,   488,   489,     1,     2,   965,   297,
     490,   491,   763,   173,   174,   248,     3,     4,     5,   492,
     297,   966,   250,     6,   764,   493,   494,     7,     8,     9,
     765,    10,  1015,    11,    12,    13,    14,   495,   252,   635,
     636,   637,   638,  1016,   639,   496,    15,   497,   498,    16,
     253,   505,   635,   636,   637,   638,   254,   639,   499,   508,
     509,   510,    17,   838,   513,   516,   519,   539,   173,   174,
     520,   521,   522,    18,    19,    20,   340,   523,   524,    21,
     525,   526,   973,   527,   975,   341,   528,   529,   530,    22,
     531,    23,   532,    24,    25,    26,    27,    28,   533,   534,
     535,   538,    29,    30,   540,   541,   545,    31,    32,    33,
      34,   547,   555,   340,   342,   550,   544,    35,    36,   546,
      37,   552,   341,   551,    38,   559,   561,    39,    40,    41,
      42,   159,   160,   161,   162,   163,   164,   165,   195,   560,
     571,   574,   196,   166,   575,   605,   578,   167,   606,   569,
     570,   342,   159,   160,   161,   162,   163,   164,   165,   581,
     582,   583,   607,   168,   166,   197,   584,   585,   167,   343,
     344,   586,   587,   588,   589,   345,   346,   347,   348,   349,
     350,   351,   352,   353,   168,   590,   598,   362,   608,   609,
     354,   610,   355,   592,   593,   594,   595,   611,   596,   612,
     613,   614,   615,   169,   170,   601,   343,   344,   363,   616,
     364,   365,   345,   346,   347,   348,   349,   350,   351,   352,
     353,   171,   172,   618,   169,   170,   619,   354,   617,   355,
     230,   620,   366,   367,   621,   624,   340,   231,   625,   626,
     627,   628,   171,   172,   324,   341,   592,   593,   594,   595,
     656,   596,   629,   630,   173,   174,   592,   593,   594,   595,
     631,   596,   632,   633,   411,   659,   660,   368,   597,   369,
     250,   370,   329,   412,   342,   173,   174,   371,   662,   664,
     930,   372,   661,   663,   665,   330,   670,   671,   674,   373,
     374,   375,   673,   682,   675,   376,   377,   378,   680,   209,
     681,   683,   413,   362,   686,   687,   458,   592,   593,   594,
     595,   688,   596,   689,   691,   693,   697,   698,   692,   694,
     592,   593,   594,   595,   363,   596,   364,   365,   809,   343,
     344,   931,   700,   677,   701,   345,   346,   347,   348,   349,
     350,   351,   352,   353,   932,   704,   230,   705,   366,   367,
     354,   708,   355,   231,   711,   712,   713,   414,   415,   714,
     324,   715,   716,   416,   417,   418,   419,   420,   421,   422,
     423,   424,   592,   593,   594,   595,   719,   596,   425,   720,
     426,   721,   722,   368,   725,   369,   250,   370,   329,   723,
     726,   724,   730,   371,   731,   732,   933,   372,   745,   596,
     811,   330,   747,   755,   748,   373,   374,   375,   750,   756,
     784,   376,   377,   378,   771,   209,   592,   593,   594,   595,
     757,   596,   592,   593,   594,   595,   807,   596,   592,   593,
     594,   595,   772,   596,   592,   593,   594,   595,   812,   596,
     934,   813,   592,   593,   594,   595,   935,   596,   592,   593,
     594,   595,   936,   596,   592,   593,   594,   595,   937,   596,
     635,   636,   637,   638,   814,   639,   938,   815,   635,   636,
     637,   638,   939,   639,   635,   636,   637,   638,   940,   639,
     635,   636,   637,   638,   946,   639,   635,   636,   637,   638,
     816,   639,   947,   817,   635,   636,   637,   638,   948,   639,
     635,   636,   637,   638,   949,   639,   635,   636,   637,   638,
     950,   639,   635,   636,   637,   638,   818,   639,   951,   821,
     635,   636,   637,   638,   952,   639,   592,   593,   594,   595,
     953,   596,   592,   593,   594,   595,   954,   596,   592,   593,
     594,   595,   822,   596,   955,   826,   592,   593,   594,   595,
     958,   596,   592,   593,   594,   595,  1011,   596,   592,   593,
     594,   595,  1042,   596,   592,   593,   594,   595,   639,   596,
    1093,   823,   592,   593,   594,   595,  1097,   596,   844,   864,
     808,   858,   592,   593,   594,   595,   889,   596,   635,   636,
     637,   638,   860,   639,   943,   592,   593,   594,   595,   868,
     596,   871,   869,   870,   971,   872,   592,   593,   594,   595,
     991,   596,   873,   635,   636,   637,   638,  1001,   639,   592,
     593,   594,   595,   877,   596,   880,   878,   879,  1008,   592,
     593,   594,   595,   501,   596,  1027,   882,   592,   593,   594,
     595,  1040,   596,   592,   593,   594,   595,   883,   596,   886,
     890,  1041,   592,   593,   594,   595,   888,   596,   894,  1044,
     895,   592,   593,   594,   595,  1052,   596,   592,   593,   594,
     595,   898,   596,   899,  1055,   900,   901,   902,   904,   905,
     907,   908,   909,  1056,   911,   913,   916,   923,   925,  1078,
      -1,   942,   945,   960,   968,   967,   970,   969,   976,   983,
     985,   986,   987,   988,   992,  1000,  1002,  1003,  1007,  1012,
    1017,  1019,  1023,  1020,  1028,  1021,  1043,  1033,  1045,  1049,
    1057,  1036,  1059,  1058,  1060,  1050,  1051,  1064,  1065,  1068,
    1069,  1083,  1070,  1079,  1091,  1092,  1090,  1098,   157,  1082,
     993,   557,   446,   959,   699,   859,   857,   543,   549,   554,
     727,   444,   600,  1018,   845,   504,   865,   998,   685,   775,
     327,     0,   876,     0,     0,   710,     0,     0,   690
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   224,   132,   113,   482,    96,    31,   107,   107,
     343,   344,    36,    37,    38,   107,   139,   107,    42,   107,
     107,   107,   355,    91,   119,   125,   563,   360,   621,   563,
     119,   567,   694,   900,    25,   133,   134,   544,    33,   546,
     139,   133,   134,   143,   551,   133,   134,   133,   134,    91,
     782,   141,   107,    51,    51,   107,   107,   107,    42,   107,
      91,    51,   107,   107,   107,   107,    33,    80,    42,    42,
      22,    51,    42,    80,   338,    80,   107,    51,    51,    31,
     134,    51,   107,     0,   157,    31,   107,    31,   107,   107,
       7,     8,    51,   916,   917,   918,   318,   159,    31,  1020,
      17,    18,    19,   165,   158,    31,    97,    24,    97,    51,
     157,    28,    29,    30,   134,    32,    42,    34,    35,    36,
      37,    51,    22,     6,     4,    51,    33,    22,   117,    80,
      47,    31,   159,    50,    80,    42,    80,    33,   165,   159,
     135,   136,   137,   138,   157,   140,    63,    80,   157,    80,
     157,   135,   136,   158,    80,   135,   136,    74,    75,    76,
     134,    44,   157,    80,   162,    42,   163,   134,   158,    52,
       4,  1092,   157,    90,    51,    92,   128,    94,    95,    96,
      97,    98,   446,   516,   158,   158,   103,   104,   158,   157,
     157,   108,   109,   110,   111,    51,   157,    31,  1065,   158,
      80,   118,   119,    80,   121,  1028,    89,   158,   125,   135,
     136,   128,   129,   130,   131,   141,   142,   143,   144,   145,
     146,   147,   148,   149,    31,    22,   157,   157,   128,    42,
     156,   158,   158,   128,   160,    42,     6,   157,    51,    80,
     747,   124,    31,   750,    51,   578,    80,    80,   581,   582,
     583,   584,   585,   586,   587,   588,   589,   590,   135,   592,
     593,   594,   595,   596,   157,   598,   998,    80,    22,    31,
      80,   604,    42,    80,    44,    42,    80,    31,   134,    42,
      42,    51,    52,    80,    51,   414,   415,   164,    51,    51,
      56,    80,   157,    80,    22,    31,   134,   426,    80,    80,
      80,   157,    22,    31,   134,   134,    42,   969,    22,    22,
      80,    31,    78,    80,   647,    51,   157,    80,    80,    89,
     117,   159,   135,   136,   157,   158,    80,   157,   135,   136,
     159,   157,   158,   137,   141,   142,   143,   144,   145,   146,
     147,   148,   149,    31,    80,   158,    80,   157,   158,   156,
      80,   158,    80,   160,   124,    80,    80,   161,    80,    60,
      80,   165,   166,   117,   697,   135,    80,    80,   135,   136,
     157,   158,   135,   135,   136,   157,   158,   157,   914,   141,
     142,   143,   144,   145,   146,   147,   148,   149,    80,   117,
     157,   158,    80,   900,   156,   902,   158,   117,   160,   135,
     136,   134,    31,   117,   117,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   921,   116,   506,   134,   134,   158,
     156,   158,   158,   157,   160,   126,   159,   157,   496,    77,
     134,   134,   157,   157,   157,   157,    33,   536,   536,   134,
      88,   134,   159,   159,   536,   134,   536,   134,   536,   536,
     536,    80,   785,   101,   496,   159,   159,    80,   995,   996,
     997,   995,   996,   997,   159,   496,   159,  1060,    80,    80,
     159,    33,   159,   157,   807,    80,   983,   572,   985,   986,
     157,   536,  1075,   572,   536,   536,   536,   511,   536,    51,
     514,   536,   536,   536,   536,   624,   625,   626,   627,   628,
     629,   630,   631,   632,   633,   536,   635,   636,   637,   638,
     639,   536,  1019,   537,  1021,   536,  1053,   536,   536,  1053,
      80,    80,    42,   856,   157,  1032,    80,   650,    80,   134,
     163,    51,    42,    80,    80,    80,    80,    33,  1045,   134,
     157,    51,  1049,   641,   157,   158,    42,    31,   165,   641,
      42,   650,   157,   641,   654,   641,    42,  1064,    42,    51,
      80,    23,   652,    80,   159,    51,    42,    51,   135,   136,
     137,   138,    42,   140,   134,    51,   115,    42,   137,    80,
     134,    51,   134,  1060,  1091,    42,    51,   134,   134,   134,
     134,  1098,   159,   617,   927,   928,    80,   157,  1075,   158,
      80,   157,   161,   157,    33,   157,   165,   166,    80,   942,
     157,   157,   157,   157,    80,   135,   136,    80,   642,   643,
     137,   141,   142,   143,   144,   145,   146,   147,   148,   149,
     963,   655,   761,   134,   658,    80,   156,   970,   158,    80,
     160,   158,   157,    80,   161,   157,   158,   756,   165,   166,
     165,   135,   136,    80,   134,   988,   157,   141,   142,   143,
     144,   145,   146,   147,   148,   149,     5,  1000,   134,  1002,
      38,   134,   156,    80,   158,    42,   160,   157,   157,   158,
     157,    20,    21,    80,   137,   138,    25,   140,   910,   134,
      80,   157,   123,    51,   157,    80,    80,    80,   157,   158,
      39,    40,    41,    80,    43,   137,   138,    46,   140,   157,
    1043,    80,   157,   162,    53,    54,    55,   157,   158,    58,
      59,     5,    61,    62,    80,    64,    65,    66,    67,    68,
      69,    70,    71,    72,    73,    80,    20,   123,    33,    33,
      79,  1074,    81,    82,   134,    84,    80,    86,   123,   134,
     134,   134,    91,   123,  1087,    39,    40,   134,    33,    43,
      99,   100,    46,   102,    33,   134,    33,   157,    33,    53,
      33,    24,   157,   157,   157,   114,   115,    33,   134,    33,
     157,   120,   911,    33,    33,    33,     7,     8,   157,   134,
      33,   134,    45,   132,   133,    79,    17,    18,    19,   159,
     134,   157,    86,    24,    57,    33,    33,    28,    29,    30,
      63,    32,   157,    34,    35,    36,    37,    33,   102,   135,
     136,   137,   138,   157,   140,   134,    47,   159,   159,    50,
     114,    33,   135,   136,   137,   138,   120,   140,    80,   157,
     157,    80,    63,   159,    80,    33,    80,    80,   132,   133,
      33,    33,    33,    74,    75,    76,    42,    33,    33,    80,
      33,    33,   886,    33,   888,    51,    33,    33,    33,    90,
      33,    92,    33,    94,    95,    96,    97,    98,    33,    33,
      33,    33,   103,   104,    80,    80,    80,   108,   109,   110,
     111,    80,   158,    42,    80,    80,   134,   118,   119,   134,
     121,    80,    51,   134,   125,    33,   157,   128,   129,   130,
     131,     9,    10,    11,    12,    13,    14,    15,    16,   159,
      33,    80,    20,    21,    80,    33,   158,    25,    33,   159,
     159,    80,     9,    10,    11,    12,    13,    14,    15,   158,
     158,   158,    33,    41,    21,    43,   158,   158,    25,   135,
     136,   158,   158,   158,   158,   141,   142,   143,   144,   145,
     146,   147,   148,   149,    41,   158,   158,     3,    33,    33,
     156,    33,   158,   135,   136,   137,   138,    33,   140,    33,
      33,    33,    33,    81,    82,   134,   135,   136,    24,    33,
      26,    27,   141,   142,   143,   144,   145,   146,   147,   148,
     149,    99,   100,   123,    81,    82,    80,   156,   159,   158,
      46,   157,    48,    49,   158,   158,    42,    53,   158,   158,
     158,   158,    99,   100,    60,    51,   135,   136,   137,   138,
      33,   140,   158,   158,   132,   133,   135,   136,   137,   138,
     158,   140,   158,   158,    42,   123,    80,    83,   157,    85,
      86,    87,    88,    51,    80,   132,   133,    93,    80,    80,
     159,    97,   123,   123,    42,   101,    42,    51,    80,   105,
     106,   107,    51,    80,    51,   111,   112,   113,    51,   115,
      51,    80,    80,     3,    51,    51,   122,   135,   136,   137,
     138,    42,   140,    51,    51,   158,    33,    51,   157,    80,
     135,   136,   137,   138,    24,   140,    26,    27,   134,   135,
     136,   159,   157,   162,   157,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   159,    80,    46,   157,    48,    49,
     156,    80,   158,    53,    51,    51,    51,   135,   136,    51,
      60,    51,    42,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   135,   136,   137,   138,    42,   140,   156,    51,
     158,    51,    51,    83,   158,    85,    86,    87,    88,    51,
      42,    80,    80,    93,    33,   134,   159,    97,   134,   140,
      51,   101,   134,    80,   134,   105,   106,   107,   134,   157,
     134,   111,   112,   113,   157,   115,   135,   136,   137,   138,
      80,   140,   135,   136,   137,   138,    80,   140,   135,   136,
     137,   138,   157,   140,   135,   136,   137,   138,    42,   140,
     159,    51,   135,   136,   137,   138,   159,   140,   135,   136,
     137,   138,   159,   140,   135,   136,   137,   138,   159,   140,
     135,   136,   137,   138,    51,   140,   159,    51,   135,   136,
     137,   138,   159,   140,   135,   136,   137,   138,   159,   140,
     135,   136,   137,   138,   159,   140,   135,   136,   137,   138,
      51,   140,   159,    51,   135,   136,   137,   138,   159,   140,
     135,   136,   137,   138,   159,   140,   135,   136,   137,   138,
     159,   140,   135,   136,   137,   138,    51,   140,   159,    51,
     135,   136,   137,   138,   159,   140,   135,   136,   137,   138,
     159,   140,   135,   136,   137,   138,   159,   140,   135,   136,
     137,   138,    51,   140,   159,   123,   135,   136,   137,   138,
     159,   140,   135,   136,   137,   138,   159,   140,   135,   136,
     137,   138,   159,   140,   135,   136,   137,   138,   140,   140,
     159,    51,   135,   136,   137,   138,   159,   140,   157,    42,
     157,   157,   135,   136,   137,   138,   157,   140,   135,   136,
     137,   138,   157,   140,   157,   135,   136,   137,   138,   123,
     140,    42,   123,   123,   157,    42,   135,   136,   137,   138,
     157,   140,    51,   135,   136,   137,   138,   157,   140,   135,
     136,   137,   138,   161,   140,    51,   157,   157,   157,   135,
     136,   137,   138,   158,   140,   157,   159,   135,   136,   137,
     138,   157,   140,   135,   136,   137,   138,   159,   140,   159,
     164,   157,   135,   136,   137,   138,   159,   140,    80,   157,
      80,   135,   136,   137,   138,   157,   140,   135,   136,   137,
     138,   157,   140,   134,   157,   134,    80,   134,    80,   157,
      51,    51,   159,   157,    33,    33,   157,    80,    80,   157,
     140,    80,   159,   157,    80,   161,    33,   159,    80,   134,
     134,   134,   157,    33,    42,    33,    33,    51,   157,    51,
     164,   134,   157,   134,   157,   134,    33,   159,   134,   134,
      51,   161,    51,   161,   127,   157,   157,   134,   134,   159,
      51,    51,   161,    80,   134,   134,   164,   134,    43,   157,
     914,   316,   136,   851,   506,   652,   650,   302,   308,   313,
     536,   134,   358,   977,   641,   207,   657,   919,   491,   572,
     118,    -1,   678,    -1,    -1,   520,    -1,    -1,   496
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
      10,    11,    12,    13,    14,    15,    21,    25,    41,    81,
      82,    99,   100,   132,   133,   283,   284,   285,   310,   311,
     312,   313,   314,   345,   346,   348,   349,   352,   353,   354,
     355,   356,   357,   358,   157,    16,    20,    43,   284,   287,
     288,   318,   335,   359,    23,     4,    80,   265,   266,   115,
     221,   222,   291,    42,   157,    51,   157,   157,   165,    80,
     157,   165,    80,    80,   190,   191,    33,     5,    39,    40,
      46,    53,    54,    55,    58,    59,    61,    62,    64,    65,
      66,    67,    68,    69,    70,    71,    72,    73,    79,    84,
      86,    91,   102,   114,   120,   247,   248,   291,   301,   310,
     311,   312,   313,   314,   315,   316,   317,   318,   319,   320,
     321,   322,   323,   324,   325,   326,   327,   328,   329,   330,
     331,   332,   333,   334,   335,   336,   337,   338,   340,   341,
     345,   346,   347,   348,   349,   350,    80,   134,   157,    22,
      80,   117,   233,   234,   235,    22,    80,   117,   242,   243,
      22,    80,   117,   239,   240,    80,   193,   194,   190,    38,
     188,    42,   157,   198,    60,   116,   126,   293,    77,    88,
     101,   272,   273,   342,   343,   344,    22,   128,   210,   211,
      42,    51,    80,   135,   136,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   156,   158,   185,    80,   257,   258,
      80,   261,     3,    24,    26,    27,    48,    49,    83,    85,
      87,    93,    97,   105,   106,   107,   111,   112,   113,   228,
     290,   291,   292,   293,   294,   295,   296,   297,   298,   299,
     300,   301,   302,   303,   304,   305,   307,   308,   309,   317,
     339,   343,   344,   157,   157,   123,    80,   134,   157,    51,
     157,    42,    51,    80,   135,   136,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   156,   158,   205,   207,   250,
     251,   301,   317,   318,   327,   333,   334,   335,   336,   337,
     338,   345,   346,   347,   250,   157,   210,   162,   224,   225,
     304,   310,   218,   219,   291,   227,   228,   157,   122,   228,
     281,   282,   351,   157,   157,   123,    80,   134,   157,   123,
      80,   134,   157,   123,    80,   134,   157,   157,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,   134,   159,    33,    33,    33,   134,   159,   159,    80,
     134,   158,   267,    31,   266,    33,   134,   159,   157,   157,
      80,   159,   165,    80,   159,   165,    33,    31,   191,    80,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,   134,   159,    33,    80,
      80,    80,    31,   234,   134,    80,   134,    80,    31,   243,
      80,   134,    80,    31,   240,   158,    31,   194,    31,    33,
     159,   157,   160,   203,   204,   205,   206,   134,   159,   159,
     159,    33,   134,   159,    80,    80,    31,   211,   158,   185,
     185,   158,   158,   158,   158,   158,   158,   158,   158,   158,
     158,   185,   135,   136,   137,   138,   140,   157,   158,    31,
     258,   134,   185,    31,    80,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,   159,   123,    80,
     157,   158,   205,   205,   158,   158,   158,   158,   158,   158,
     158,   158,   158,   158,   205,   135,   136,   137,   138,   140,
     279,   134,   159,   159,    31,    42,    51,   158,   215,   216,
     134,   159,   134,   159,   134,   159,    33,   134,   159,   123,
      80,   123,    80,   123,    80,    42,    42,   135,   136,   231,
      42,    51,   231,    51,    80,    51,    51,   162,   362,   363,
      51,    51,    80,    80,   360,   285,    51,    51,    42,    51,
     288,    51,   157,   158,    80,    42,    51,    33,    51,   222,
     157,   157,   157,   229,    80,   157,   157,   229,    80,   185,
     363,    51,    51,    51,    51,    51,    42,    42,    51,    42,
      51,    51,    51,    51,    80,   158,    42,   248,   157,   229,
      80,    33,   134,     6,    42,    44,    51,    52,    80,    89,
     124,   135,   236,   244,   245,   134,   245,   134,   134,   245,
     134,    51,   135,   136,   230,    80,   157,    80,    31,   204,
     206,    33,   157,    45,    57,    63,   195,   196,   305,   306,
     202,   157,   157,    56,    78,   273,    80,   137,   161,   165,
     166,   274,   275,   276,   134,    33,   134,   157,   185,   186,
     185,   185,   185,   185,   185,   185,   185,   185,   185,   185,
     159,   185,   185,   185,   185,   185,   185,    80,   157,   134,
     185,    51,    42,    51,    51,    51,    51,    51,    51,    42,
      51,    51,    51,    51,   157,   229,   123,   230,   205,   205,
     205,   205,   205,   205,   205,   205,   205,   205,   159,   205,
     205,   205,   205,   205,   157,   251,   157,   229,   157,   229,
     185,   157,   163,    42,    51,   134,   158,   225,   157,   219,
     157,   228,   157,   229,    42,   282,   157,   229,   123,   123,
     123,    42,    42,    51,   361,   163,   361,   161,   157,   157,
      51,   267,   159,   159,   185,   157,   159,   157,   159,   157,
     164,   253,   254,   157,    80,    80,    42,    51,   157,   134,
     134,    80,   134,   245,    80,   157,   245,    51,    51,   159,
     190,    33,   205,    33,   134,   159,   157,   200,   199,   134,
     157,   158,   276,    80,   185,    80,    97,   117,   134,   159,
     159,   159,   159,   159,   159,   159,   159,   159,   159,   159,
     159,   185,    80,   157,   157,   159,   159,   159,   159,   159,
     159,   159,   159,   159,   159,   159,   157,   157,   159,   216,
     157,    42,    51,   158,   185,   157,   157,   161,    80,   159,
      33,   157,   157,   229,   157,   229,    80,   134,   159,   237,
     245,   244,   245,   134,   245,   134,   134,   157,    33,    31,
     205,   157,    42,   196,   201,   203,   203,   203,   275,   245,
      33,   157,    33,    51,   212,   185,   185,   157,   157,   185,
     185,   159,    51,   267,   185,   157,   157,   164,   253,   134,
     134,   134,   245,   157,   245,   245,   185,   157,   157,    31,
      31,    31,   158,   159,   185,   185,   161,    51,   134,   157,
     157,   157,   159,    33,   157,   134,   245,   237,   245,   134,
     157,   157,   157,   203,   245,   157,   157,    51,   161,    51,
     127,   185,   164,   245,   134,   134,   245,    31,   159,    51,
     161,    80,   135,   136,   158,   213,   230,   231,   157,    80,
     245,   244,   157,    51,   185,    80,   157,   158,   230,   231,
     164,   134,   134,   159,   185,   245,   237,   159,   134,   245
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
     284,   284,   284,   284,   284,   284,   284,   284,   284,   285,
     285,   286,   286,   287,   287,   287,   287,   288,   288,   289,
     289,   290,   291,   292,   293,   294,   295,   296,   297,   298,
     299,   300,   301,   302,   303,   304,   305,   306,   307,   308,
     309,   309,   310,   311,   311,   312,   313,   314,   315,   316,
     317,   317,   318,   319,   320,   321,   322,   323,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   342,   343,
     344,   345,   346,   347,   348,   349,   350,   351,   352,   353,
     354,   355,   356,   357,   358,   359,   360,   361,   361,   362,
     362,   363
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
       1,     1,     1,     1,     1,     1,     1,     1,     1,     3,
       1,     3,     6,     1,     1,     1,     1,     3,     1,     3,
       6,     3,     3,     3,     1,     3,     3,     3,     3,     1,
       1,     1,     3,     3,     3,     3,     3,     3,     1,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       1,     1,     3,     3,     3,     3,     5,     3,     3,     3,
       3,     1,     3,     3,     3,     1,     1,     1,     1,     1,
       3,     1,     1,     1,     1,     3,     3,     3,     3,     1,
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
     314,    -1,   311,    -1,   313,    -1,   348,    -1,   349,    -1,
     284,   134,   285,    -1,   284,    -1,     7,    51,   157,    -1,
       7,   158,   285,   159,    51,   157,    -1,   284,    -1,   335,
      -1,   318,    -1,   359,    -1,   287,   134,   288,    -1,   287,
      -1,     8,    51,   157,    -1,     8,   158,   288,   159,    51,
     157,    -1,    26,    33,    51,    -1,   115,    33,    51,    -1,
     112,    33,    51,    -1,    60,    -1,    93,    33,    51,    -1,
     107,    33,    51,    -1,    27,    33,    51,    -1,     3,    33,
      51,    -1,    83,    -1,    85,    -1,    87,    -1,    53,    33,
      51,    -1,    48,    33,    51,    -1,    49,    33,    51,    -1,
      97,    33,    51,    -1,    24,    33,    42,    -1,    63,    33,
      42,    -1,   111,    -1,   113,    33,    51,    -1,   105,    33,
      51,    -1,   105,    33,    42,    -1,    25,    33,    80,    -1,
      81,    33,   363,    -1,    81,    33,    51,    -1,    41,    33,
      51,    -1,    99,    33,    51,    -1,   100,    33,    51,    -1,
      58,    33,    51,    -1,    59,    33,    51,    -1,    86,    -1,
      46,    -1,    20,    33,    42,    -1,    69,    33,    51,    -1,
      64,    33,    42,    -1,    66,    33,    42,    -1,    91,    33,
     158,   254,   159,    -1,    65,    33,    42,    -1,    65,    33,
      51,    -1,    73,    33,    80,    -1,    72,    33,    51,    -1,
      71,    -1,   102,    33,    42,    -1,    67,    33,    51,    -1,
      68,    33,    51,    -1,    61,    -1,    62,    -1,    84,    -1,
       5,    -1,   120,    -1,    43,    33,    51,    -1,   114,    -1,
      79,    -1,    40,    -1,   106,    -1,    54,    33,    51,    -1,
      55,    33,    51,    -1,    77,    33,    56,    -1,    77,    33,
      78,    -1,   101,    -1,    88,    -1,   132,    33,    80,    -1,
     133,    33,   360,    -1,    39,    33,   363,    -1,    21,    -1,
      82,    -1,    70,    -1,   122,    33,    42,    -1,    14,    33,
     231,    -1,     9,    33,    42,    -1,    11,    33,   231,    -1,
      12,    33,    42,    -1,    13,    33,    51,    -1,    10,    -1,
      15,    33,    51,    -1,    16,    33,    51,    -1,    80,   161,
      80,    -1,    51,    -1,    51,   161,    51,    -1,   162,   361,
      -1,   362,   361,    -1,   362,   163,    -1
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
    1542,  1544,  1546,  1548,  1550,  1552,  1554,  1556,  1558,  1560,
    1564,  1566,  1570,  1577,  1579,  1581,  1583,  1585,  1589,  1591,
    1595,  1602,  1606,  1610,  1614,  1616,  1620,  1624,  1628,  1632,
    1634,  1636,  1638,  1642,  1646,  1650,  1654,  1658,  1662,  1664,
    1668,  1672,  1676,  1680,  1684,  1688,  1692,  1696,  1700,  1704,
    1708,  1710,  1712,  1716,  1720,  1724,  1728,  1734,  1738,  1742,
    1746,  1750,  1752,  1756,  1760,  1764,  1766,  1768,  1770,  1772,
    1774,  1778,  1780,  1782,  1784,  1786,  1790,  1794,  1798,  1802,
    1804,  1806,  1810,  1814,  1818,  1820,  1822,  1824,  1828,  1832,
    1836,  1840,  1844,  1848,  1850,  1854,  1858,  1862,  1864,  1868,
    1871,  1874
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
    1119,  1120,  1121,  1122,  1123,  1124,  1125,  1126,  1127,  1130,
    1131,  1134,  1136,  1140,  1141,  1142,  1143,  1146,  1147,  1150,
    1152,  1156,  1157,  1158,  1159,  1160,  1161,  1162,  1163,  1164,
    1165,  1166,  1167,  1168,  1169,  1170,  1171,  1172,  1173,  1174,
    1175,  1176,  1178,  1179,  1180,  1182,  1183,  1184,  1185,  1186,
    1187,  1188,  1189,  1190,  1191,  1192,  1193,  1194,  1195,  1196,
    1197,  1198,  1199,  1200,  1201,  1202,  1203,  1204,  1205,  1206,
    1207,  1208,  1209,  1210,  1211,  1212,  1213,  1215,  1217,  1220,
    1221,  1222,  1223,  1224,  1225,  1226,  1227,  1228,  1230,  1231,
    1232,  1233,  1234,  1235,  1236,  1237,  1239,  1247,  1248,  1252,
    1253,  1262
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
  const int parser::yylast_ = 1548;
  const int parser::yynnts_ = 197;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 156;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 167;

  const unsigned int parser::yyuser_token_number_max_ = 411;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1264 "DynareBison.yy"


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

