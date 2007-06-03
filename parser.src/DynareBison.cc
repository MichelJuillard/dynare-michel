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
#line 139 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 48:
#line 140 "DynareBison.yy"
    {driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 49:
#line 143 "DynareBison.yy"
    {driver.rplot();;}
    break;

  case 54:
#line 163 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 55:
#line 165 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 56:
#line 167 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 57:
#line 169 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 58:
#line 171 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 59:
#line 173 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 60:
#line 178 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 61:
#line 180 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 62:
#line 182 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 63:
#line 184 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 64:
#line 186 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 65:
#line 188 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 66:
#line 193 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 67:
#line 195 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 68:
#line 197 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 69:
#line 199 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 70:
#line 201 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 71:
#line 203 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 72:
#line 208 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 73:
#line 210 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 74:
#line 212 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 75:
#line 214 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 76:
#line 216 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 77:
#line 218 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 78:
#line 223 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 79:
#line 227 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 80:
#line 234 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 81:
#line 238 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 82:
#line 245 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 83:
#line 249 "DynareBison.yy"
    {
      driver.markowitz((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 84:
#line 257 "DynareBison.yy"
    {driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 85:
#line 262 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 86:
#line 264 "DynareBison.yy"
    {(yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 87:
#line 266 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 88:
#line 268 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 89:
#line 270 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 90:
#line 272 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 91:
#line 274 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 92:
#line 276 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 93:
#line 278 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 94:
#line 280 "DynareBison.yy"
    {(yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 95:
#line 282 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 96:
#line 284 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 97:
#line 286 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 98:
#line 288 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 99:
#line 290 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 100:
#line 292 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 101:
#line 294 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 102:
#line 296 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 103:
#line 298 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 104:
#line 300 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 105:
#line 302 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 106:
#line 304 "DynareBison.yy"
    {(yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 107:
#line 309 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 108:
#line 311 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 109:
#line 315 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 110:
#line 317 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 111:
#line 321 "DynareBison.yy"
    {driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 112:
#line 326 "DynareBison.yy"
    {driver.end_endval();;}
    break;

  case 115:
#line 336 "DynareBison.yy"
    {driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 116:
#line 341 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 119:
#line 351 "DynareBison.yy"
    {driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 122:
#line 359 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 123:
#line 360 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 126:
#line 366 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 127:
#line 366 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 128:
#line 367 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 129:
#line 368 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 130:
#line 369 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 131:
#line 370 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 132:
#line 371 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 133:
#line 372 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 134:
#line 373 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 135:
#line 374 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 140:
#line 386 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 141:
#line 388 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val));;}
    break;

  case 142:
#line 392 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 144:
#line 395 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 145:
#line 397 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 146:
#line 399 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 147:
#line 401 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 148:
#line 403 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 149:
#line 405 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 150:
#line 407 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 151:
#line 409 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 152:
#line 411 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 153:
#line 413 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 154:
#line 415 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 155:
#line 417 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 156:
#line 419 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 157:
#line 421 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 158:
#line 423 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 159:
#line 425 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 160:
#line 427 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 161:
#line 429 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 162:
#line 431 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 163:
#line 435 "DynareBison.yy"
    {driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 164:
#line 439 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 165:
#line 441 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 166:
#line 445 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 167:
#line 449 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 170:
#line 459 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 171:
#line 461 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 172:
#line 463 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 173:
#line 465 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 174:
#line 467 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 175:
#line 472 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 176:
#line 474 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 177:
#line 476 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 178:
#line 478 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 179:
#line 480 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 180:
#line 482 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 181:
#line 488 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 182:
#line 490 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 183:
#line 492 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 184:
#line 494 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 185:
#line 496 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 186:
#line 498 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 187:
#line 500 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 188:
#line 502 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 189:
#line 507 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 190:
#line 512 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 191:
#line 514 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 192:
#line 519 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 193:
#line 521 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 194:
#line 523 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 195:
#line 525 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 196:
#line 527 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 197:
#line 529 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 198:
#line 531 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 199:
#line 533 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 200:
#line 535 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 201:
#line 540 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 202:
#line 544 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 206:
#line 556 "DynareBison.yy"
    {driver.check();;}
    break;

  case 207:
#line 558 "DynareBison.yy"
    {driver.check();;}
    break;

  case 211:
#line 570 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 212:
#line 572 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 217:
#line 585 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 218:
#line 587 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 219:
#line 589 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 220:
#line 591 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 243:
#line 622 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 244:
#line 624 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 245:
#line 626 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 246:
#line 628 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 247:
#line 630 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 248:
#line 632 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 249:
#line 637 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 250:
#line 639 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 251:
#line 641 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 252:
#line 646 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 253:
#line 648 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 254:
#line 650 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 255:
#line 655 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 256:
#line 660 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 257:
#line 662 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 259:
#line 671 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 260:
#line 676 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 261:
#line 681 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 262:
#line 691 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 263:
#line 696 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 264:
#line 702 "DynareBison.yy"
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

  case 265:
#line 712 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 266:
#line 717 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 267:
#line 728 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 268:
#line 734 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 269:
#line 742 "DynareBison.yy"
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

  case 270:
#line 752 "DynareBison.yy"
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

  case 271:
#line 766 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 272:
#line 770 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 273:
#line 772 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 274:
#line 776 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 275:
#line 783 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 276:
#line 792 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 277:
#line 801 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 278:
#line 805 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 279:
#line 807 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 280:
#line 811 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 281:
#line 820 "DynareBison.yy"
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

  case 282:
#line 831 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 283:
#line 843 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 284:
#line 845 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 285:
#line 847 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 286:
#line 849 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 287:
#line 851 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 288:
#line 855 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 292:
#line 860 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 293:
#line 862 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 294:
#line 869 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 295:
#line 871 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 296:
#line 873 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 297:
#line 875 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 339:
#line 926 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 340:
#line 928 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 356:
#line 954 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 357:
#line 956 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 358:
#line 960 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 359:
#line 961 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 362:
#line 971 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 363:
#line 976 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 366:
#line 985 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val));;}
    break;

  case 367:
#line 988 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 368:
#line 992 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 369:
#line 996 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 370:
#line 998 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 371:
#line 1000 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 372:
#line 1002 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 373:
#line 1005 "DynareBison.yy"
    {driver.set_osr_params();;}
    break;

  case 374:
#line 1008 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 375:
#line 1009 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 376:
#line 1010 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 377:
#line 1011 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 378:
#line 1014 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 379:
#line 1015 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 380:
#line 1016 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 381:
#line 1017 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 386:
#line 1028 "DynareBison.yy"
    {driver.set_olr_inst();;}
    break;

  case 387:
#line 1032 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 390:
#line 1039 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 391:
#line 1040 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 392:
#line 1041 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val));;}
    break;

  case 393:
#line 1044 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 394:
#line 1045 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 395:
#line 1046 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 396:
#line 1049 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 397:
#line 1050 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 398:
#line 1053 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 399:
#line 1054 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 400:
#line 1055 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 401:
#line 1056 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 402:
#line 1057 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 403:
#line 1058 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 404:
#line 1060 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 405:
#line 1061 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 406:
#line 1062 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 407:
#line 1063 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 408:
#line 1064 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 409:
#line 1065 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 410:
#line 1068 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 416:
#line 1080 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 417:
#line 1081 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 418:
#line 1082 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 419:
#line 1083 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 420:
#line 1086 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 421:
#line 1087 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 423:
#line 1091 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 424:
#line 1092 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 425:
#line 1093 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 426:
#line 1094 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 427:
#line 1097 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 428:
#line 1097 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 430:
#line 1101 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 431:
#line 1103 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 432:
#line 1105 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 433:
#line 1107 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 438:
#line 1119 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 439:
#line 1120 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 440:
#line 1121 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 441:
#line 1122 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 442:
#line 1123 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 443:
#line 1124 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 444:
#line 1125 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 445:
#line 1126 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 446:
#line 1127 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 447:
#line 1128 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 448:
#line 1129 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 449:
#line 1130 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 450:
#line 1131 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 451:
#line 1132 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 452:
#line 1133 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 453:
#line 1134 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 454:
#line 1135 "DynareBison.yy"
    {driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 455:
#line 1136 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 456:
#line 1137 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 457:
#line 1138 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 458:
#line 1139 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 459:
#line 1141 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 460:
#line 1142 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 461:
#line 1143 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 462:
#line 1145 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 463:
#line 1146 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 464:
#line 1147 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 465:
#line 1148 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 466:
#line 1149 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 467:
#line 1150 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 468:
#line 1151 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 469:
#line 1152 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 470:
#line 1153 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 471:
#line 1154 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 472:
#line 1155 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 474:
#line 1157 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 475:
#line 1158 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 476:
#line 1159 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 477:
#line 1160 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 478:
#line 1161 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 479:
#line 1162 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 480:
#line 1163 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 481:
#line 1164 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 482:
#line 1165 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 483:
#line 1166 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 484:
#line 1167 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 485:
#line 1168 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 486:
#line 1169 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 487:
#line 1170 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 488:
#line 1171 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 489:
#line 1172 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 490:
#line 1173 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 491:
#line 1174 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 492:
#line 1175 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 493:
#line 1176 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 494:
#line 1177 "DynareBison.yy"
    {driver.option_num("olr_beta", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1180 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 496:
#line 1182 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 497:
#line 1184 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 498:
#line 1185 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 499:
#line 1186 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 500:
#line 1187 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 501:
#line 1188 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 502:
#line 1189 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 503:
#line 1190 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 504:
#line 1191 "DynareBison.yy"
    {driver.option_num("load_mh_file", "-1");;}
    break;

  case 505:
#line 1192 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 506:
#line 1195 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 508:
#line 1204 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 509:
#line 1207 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 510:
#line 1209 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 511:
#line 1217 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2182 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -914;
  const short int
  parser::yypact_[] =
  {
      1013,   394,   -64,   417,   515,    50,   -14,     4,   -55,    -4,
     -33,   -28,   -12,    80,   440,   517,   486,    84,   112,   263,
     117,   210,   235,   181,   234,   235,   237,    87,  -914,   184,
     194,   235,   217,   357,   541,   552,   239,   242,   235,   324,
     326,   328,   235,   642,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,   386,
      61,  -914,   301,   382,   290,    58,    70,   352,   457,   381,
     393,   480,  -914,   767,    -7,    38,   198,   200,   462,   393,
     535,   533,   429,  -914,    10,   339,    97,   853,   522,  -914,
    1139,     0,   302,   536,  -914,  1139,   307,   309,   487,   327,
     569,   467,   877,   106,   106,   338,    97,   465,  -914,    67,
    -914,   301,  -914,  1181,   348,  -914,  1068,   377,   385,   506,
     392,   531,   396,   547,   398,   403,  -914,  -914,   481,   575,
     -46,   212,  -914,   615,   -48,  -914,  -914,   508,  -914,   528,
    -914,  -914,   604,   428,  -914,   610,   507,   617,    43,  -914,
     614,  -914,   670,  -914,   677,   681,  -914,   682,   687,  -914,
     688,   694,   696,   711,   712,  -914,  -914,   717,   718,   719,
     732,   734,   739,  -914,  -914,   744,   745,  -914,   746,  -914,
    -914,  -914,   747,   750,   751,   752,  -914,  -914,   756,   757,
     -34,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,   758,   706,  -914,   713,  -914,   720,   248,  -914,   661,
     721,   665,   724,   272,  -914,   729,   668,   731,   277,  -914,
     653,    54,  -914,    55,   782,   654,   663,  -914,   761,  -914,
     -31,   662,   680,   794,  -914,  -914,   -29,  -914,  -914,  -914,
    -914,   764,   765,    40,  -914,  -914,  -914,   669,   853,   853,
     686,   689,   697,   699,   701,   710,   714,   715,   722,   727,
     853,  1031,   728,    62,  -914,   834,   840,   841,   846,   847,
     849,  -914,  -914,  -914,   864,   865,   868,   881,  -914,   883,
    -914,   884,   885,  -914,  -914,   151,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,   801,
      63,   221,  -914,  -914,  -914,   795,   843,  -914,   766,  -914,
    -914,  -914,   774,   877,   877,   775,   776,   780,   781,   783,
     793,   797,   798,   800,   802,   877,   753,  -914,   261,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,   264,  -914,    68,    27,   275,  -914,  -914,
    -914,   280,  -914,  -914,   304,  -914,  -914,   892,  -914,   389,
    -914,  -914,  -914,  -914,  -914,   817,   875,  -914,  -914,   830,
     888,  -914,  -914,   844,   891,  -914,  -914,   815,   816,   896,
     529,   894,  -914,  -914,   927,   301,   829,  -914,  -914,   833,
     -16,   910,   848,   183,   915,   853,  -914,  -914,  -914,   964,
     929,   842,   959,   960,   962,   966,   967,   968,   969,   991,
     543,   992,   984,   988,   989,   993,   975,    32,   897,   995,
    1006,  1017,   981,   982,   767,   189,   985,  1037,   935,  -914,
    -914,  -914,   258,   936,   231,   937,  -914,  -914,   946,   231,
     947,  -914,  -914,    12,  -914,  -914,  -914,  1005,   904,  -914,
    1011,    30,  -914,    23,  -914,   244,  -914,   926,   938,    44,
     339,   157,   963,    49,  -914,  -914,   853,   953,   595,   853,
     853,   853,   853,   853,   853,   853,   853,   853,   853,   524,
     853,   853,   853,   853,   853,  -914,   853,  -914,  -914,  1046,
    1060,  1048,  1055,  1056,  1058,   231,  1064,  1067,   546,  1072,
    1078,  1079,  1139,   196,  1042,  1077,  -914,   835,   213,  -914,
    1007,  -914,    12,   996,   628,   877,   877,   877,   877,   877,
     877,   877,   877,   877,   877,   542,   877,   877,   877,   877,
     877,   971,   106,   262,   350,  -914,  -914,  -914,   853,   527,
      53,    67,   972,   301,   976,  1181,   365,  1096,  1068,   368,
    -914,  1014,  -914,  1016,  -914,  1022,  -914,  1092,   999,   987,
    1000,   853,  -914,  -914,  -914,  -914,  -914,   406,  1019,  -914,
    -914,   427,  1023,  1214,  -914,  -914,  1099,     6,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  1008,  -914,  -914,
    -914,  -914,   990,  -914,  -914,  -914,   431,  -914,  1087,  1093,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,   564,  1029,
    1053,  1076,  1141,  1086,   231,  1147,  1080,   231,  -914,  1128,
    1177,  1069,  -914,   393,  1203,  -914,  -914,  -914,   877,  -914,
    -914,  -914,  1205,   399,  -914,  -914,  -914,  1082,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,   167,   211,
    -914,  1163,   853,  1164,   291,   791,   407,   557,   568,   585,
     619,   656,   716,   725,   826,   913,   928,  -914,   595,   595,
     953,   953,  1102,   939,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
     439,   853,  -914,  1166,  1222,  -914,   444,  -914,  1088,  1018,
    1057,  1063,  1071,  1120,  1131,  1145,  1162,  1170,  1176,  -914,
     628,   628,   996,   996,  1102,  -914,  -914,  -914,   447,  -914,
     458,  1182,    27,  1091,  -914,  -914,    37,   853,  -914,  -914,
    -914,  -914,  -914,  -914,   466,  -914,  -914,  -914,   473,  -914,
    -914,  -914,  1090,  1223,  -914,  -914,  1228,  -914,   371,  -914,
     400,  -914,  1098,  -914,  -914,  -914,  1185,  -914,   410,  1187,
    -914,  -914,  -914,  -914,  -914,  -914,   231,   258,  1136,   231,
    1138,  1142,  -914,  1116,  -914,  -914,  1245,   201,   877,  1237,
    1238,   244,  -914,   761,   761,   761,   157,  -914,   231,  -914,
    1254,  1246,  1259,  1243,   853,   853,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  1144,  -914,  1252,
     853,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,    53,  -914,  -914,
    -914,   853,  1188,  -914,  -914,   999,   853,  -914,  -914,   485,
    -914,   494,  1247,  1151,  1008,  -914,  -914,  -914,  1169,  1175,
    1193,   231,  1173,   231,   231,  -914,   853,  -914,  1260,  -914,
    -914,  -914,  1199,    56,   185,   215,   313,  1194,   853,  -914,
     853,  1209,    73,  1270,   791,  -914,  -914,  1276,  1196,  -914,
    1339,  1283,  -914,  -914,  -914,  1242,  -914,   231,   231,   231,
    1244,  -914,  1234,  1236,  1294,  -914,   761,  -914,  -914,  -914,
     231,  -914,  1301,  1307,  1332,  1240,  1351,  1275,  -914,  -914,
    -914,   853,  -914,    17,  1281,  -914,  1288,   231,  -914,  -914,
    -914,   490,  1265,  -914,  -914,  -914,  1356,  1264,    74,  1317,
    1349,  -914,   231,   223,  1271,  -914,  -914,  -914,  1387,  -914,
    -914,   570,   572,   853,   145,  -914,  -914,  -914,  1282,  1314,
    1315,  -914,  -914,  -914,  -914,  1202,  -914,  -914,   853,  -914,
    -914,  -914,   231,   231,  -914,  1208,  1316,  -914,  -914,   231,
    -914
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   427,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     2,     4,    29,    30,    44,    45,
      46,    43,     5,     6,     7,    12,     9,    10,    11,     8,
      13,    14,    15,    16,    17,    18,    19,    23,    25,    24,
      20,    21,    22,    26,    27,    28,    31,    32,    33,    38,
      39,    34,    35,    36,    37,    40,    41,    42,   396,     0,
       0,   206,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   247,   294,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   126,     0,     0,     0,     0,     0,   378,
       0,     0,     0,     0,   374,     0,     0,     0,    74,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   211,     0,
     201,     0,   217,     0,     0,   430,     0,     0,     0,    56,
       0,    62,     0,    68,     0,     0,     1,     3,     0,     0,
     393,     0,   389,     0,     0,   209,   210,     0,    80,     0,
      47,   406,     0,     0,   400,     0,     0,     0,     0,   114,
       0,   485,     0,   502,     0,     0,   490,     0,     0,   468,
       0,     0,     0,     0,     0,   482,   483,     0,     0,     0,
       0,     0,     0,   504,   478,     0,     0,   489,     0,   503,
     484,   467,     0,     0,     0,     0,   488,   486,     0,     0,
       0,   299,   335,   324,   300,   301,   302,   303,   304,   305,
     306,   307,   308,   309,   310,   311,   312,   313,   314,   315,
     316,   317,   318,   319,   320,   321,   322,   323,   325,   326,
     327,   328,   329,   330,   331,   332,   333,   334,   336,   337,
     338,   243,     0,   296,     0,   260,     0,     0,   257,     0,
       0,     0,     0,     0,   279,     0,     0,     0,     0,   273,
       0,     0,   118,     0,     0,     0,     0,    82,     0,   441,
       0,     0,     0,     0,   498,   497,     0,   412,   413,   414,
     415,     0,     0,     0,   169,    87,    88,    86,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   365,     0,     0,     0,     0,     0,
       0,   446,   447,   448,     0,     0,     0,     0,   491,     0,
     455,     0,     0,   383,   384,     0,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   236,   237,
     238,   239,   240,   241,   242,   235,   382,   380,   386,     0,
       0,     0,   376,   373,    77,    72,     0,    53,     0,    78,
     144,   145,   164,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   428,   143,     0,   342,
     347,   343,   344,   345,   346,   348,   349,   350,   351,   352,
     353,   354,   355,     0,    49,     0,     0,     0,   214,   215,
     216,     0,   204,   205,     0,   222,   219,     0,   436,     0,
     435,   437,   432,   367,    59,    54,     0,    50,    65,    60,
       0,    51,    71,    66,     0,    52,   362,     0,     0,     0,
       0,     0,   387,   388,     0,     0,     0,    81,    48,     0,
       0,     0,     0,     0,     0,     0,   112,   113,   248,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   245,     0,   259,
     255,   256,   288,     0,   288,     0,   277,   278,     0,   288,
       0,   271,   272,     0,   116,   117,   109,     0,     0,    83,
       0,     0,   138,     0,   139,     0,   134,     0,     0,     0,
       0,     0,     0,     0,   167,   168,     0,    94,    95,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    84,     0,   363,   364,     0,
       0,     0,     0,     0,     0,   288,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   368,     0,     0,    75,
      73,    79,     0,   151,   152,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   166,   199,   200,     0,     0,
     191,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      57,    55,    63,    61,    69,    67,   397,     0,   393,     0,
       0,     0,   439,   208,   207,   409,   404,     0,     0,   403,
     398,     0,     0,     0,   469,   459,     0,     0,   501,   462,
     487,   449,   492,   493,   465,   466,   471,   474,   475,   472,
     480,   481,   470,   477,   476,   461,   460,     0,   463,   464,
     479,   499,     0,   500,   298,   295,     0,   244,     0,     0,
     283,   290,   284,   289,   286,   291,   285,   287,     0,     0,
       0,   265,     0,     0,   288,     0,     0,   288,   251,     0,
       0,     0,   111,     0,     0,   127,   136,   137,     0,   141,
     123,   122,     0,     0,   121,   124,   125,     0,   130,   128,
     495,   496,   411,   422,   424,   425,   426,   423,     0,   416,
     420,     0,     0,     0,     0,   107,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    85,    90,    89,
      91,    92,    93,     0,   445,   453,   438,   444,   450,   451,
     494,   442,   452,   458,   457,   443,   440,   456,   385,   379,
       0,     0,   371,     0,     0,   375,     0,    76,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   142,
     147,   146,   148,   149,   150,   429,   341,   339,     0,   356,
       0,     0,     0,     0,   196,   197,     0,     0,   213,   212,
     203,   202,   221,   218,     0,   505,   434,   431,     0,    58,
      64,    70,     0,     0,   395,   394,     0,   405,     0,   399,
       0,   115,   507,   509,   511,   510,     0,   360,     0,     0,
     297,   246,   261,   293,   292,   258,   288,   288,     0,   288,
       0,     0,   276,     0,   250,   249,     0,     0,     0,     0,
       0,     0,   132,     0,     0,     0,     0,   410,   288,   421,
       0,     0,     0,     0,     0,     0,   106,    96,    97,    98,
      99,   100,   101,   102,   103,   104,   105,     0,   381,     0,
       0,   369,   377,   165,   153,   154,   155,   156,   157,   158,
     159,   160,   161,   162,   340,   357,   198,   190,   189,   193,
     194,     0,     0,   220,   433,   393,     0,   390,   407,     0,
     401,     0,     0,     0,     0,   473,   506,   262,     0,     0,
       0,   288,     0,   288,   288,   274,     0,   110,     0,   140,
     454,   120,     0,     0,     0,     0,   417,     0,     0,   172,
       0,   180,     0,     0,   108,   366,   372,     0,     0,   195,
       0,     0,   408,   402,   508,     0,   361,   288,   288,   288,
       0,   282,     0,     0,     0,   163,     0,   135,   131,   129,
     288,   418,     0,     0,     0,   175,     0,     0,   171,   370,
     192,     0,   391,   288,   267,   263,   266,   288,   280,   275,
     119,     0,     0,   174,   173,   179,     0,   177,     0,     0,
       0,   359,   288,     0,     0,   133,   419,   176,     0,   254,
     186,     0,     0,     0,     0,   185,   184,   392,     0,   268,
       0,   281,   178,   253,   252,     0,   183,   170,     0,   182,
     181,   358,   288,   288,   188,     0,   269,   264,   187,   288,
     270
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -914,  -914,  1406,  -914,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,  -914,  -302,  -914,
    -914,  -914,  -914,  -104,  -176,  -914,  -914,  1171,  -914,   592,
    -914,  -914,  -914,  -914,  -914,  -914,  -762,  -499,  -107,  -495,
    -914,  -914,  -914,  1318,  -255,  -914,  -914,  -914,  -914,   657,
    -914,  -914,   845,  -914,  -914,  1001,  -914,  -914,   850,  -914,
    -914,  -105,   -21,  -571,   436,  -914,  -914,  1195,  -914,  -914,
    -913,  -914,  -914,  1186,  -914,  -914,  1190,  -802,  -473,  -914,
    -914,   961,  -914,  1331,   866,  -914,   545,  -914,  -914,  -914,
    -914,  1143,  -914,  -914,  -914,  -914,  -914,  -914,   898,  1346,
    -914,  -914,  -914,  1311,  -585,  -914,  -914,  -914,  -914,  -914,
     943,  -914,   608,  -685,  -914,  -914,  -914,  -914,  -914,   857,
    -914,   -84,  -914,  1362,  -914,  -914,  -914,  -914,  -914,  -914,
    -914,   -94,  -914,  -914,  -116,  -483,  -914,  -914,  -914,  -914,
    -112,  -914,  -914,  -914,  -914,  -914,  -914,   -91,   -90,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,  -914,   -89,  -914,  -914,
    -914,  -914,  -914,   -83,   -79,   -75,   -73,   -71,   -70,  -914,
    -914,  -914,  -914,  -914,  -914,  -914,   -69,   -68,   -66,  -914,
    -914,  -914,  -914,  -914,   831,  -914,   994
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    43,    44,    45,    46,    47,    48,    49,    50,    51,
     150,   152,   154,   129,    52,    53,    54,    55,   321,   736,
      56,   285,    57,   178,   179,    58,   281,   282,   713,   714,
      59,   288,   865,   864,   942,   717,   521,   522,   523,   524,
     397,    60,    61,   303,   304,   952,  1024,    62,   609,   610,
      63,   421,   422,    64,   164,   165,    65,   417,   418,    66,
     424,   343,   104,   701,  1026,    67,   267,   268,   269,   689,
     927,    68,   278,   279,    69,   273,   274,   690,   928,    70,
     220,   221,    71,   398,   399,    72,   837,   838,    73,    74,
     323,   324,    75,    76,   370,    77,    78,    79,   344,   345,
      80,    81,   161,   162,   451,    82,    83,    84,    85,   296,
     297,   728,   729,   730,    86,   132,   601,    87,   429,   430,
     346,   347,   348,   349,   350,   351,   352,   353,   354,   355,
     356,   357,   358,   359,   360,   361,   716,   362,   363,   364,
     224,   225,   226,   227,   228,   229,   230,   401,   402,   233,
     234,   235,   236,   237,   238,   239,   240,   403,   242,   243,
     244,   245,   246,   404,   405,   406,   407,   408,   409,   365,
     253,   254,   366,   298,   299,   300,   410,   411,   412,   258,
     259,   260,   431,   673,   833,   647,   648
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       121,   122,   467,   126,   127,   283,   537,   538,   166,   223,
     135,   778,   231,   232,   241,   144,   147,   148,   549,   222,
     247,   155,   706,   419,   248,   396,   707,   420,   249,   691,
     250,   693,   251,   252,   255,   256,   696,   257,   425,   400,
     400,   428,   715,   823,   869,   929,   708,   832,   535,   681,
     264,   705,   301,   698,   101,   995,    96,   423,   683,   606,
     289,   534,   380,   261,   466,   159,   101,   575,   607,   909,
     261,   381,   732,   665,    98,   514,   516,   977,   910,   455,
     301,   449,   184,   557,   576,   804,    90,   685,   380,   605,
     720,    95,   760,   494,   805,   100,   525,   381,   530,   169,
     382,   943,   944,   945,   456,   450,  1019,   467,   265,   301,
     130,   181,   721,   177,   985,   698,   182,   105,   495,   290,
     262,   526,   106,   531,   280,   177,   382,   262,   131,   291,
    1047,   160,   322,   577,   636,   185,   186,    97,   107,   188,
     699,   700,   189,   263,  1020,   688,   102,   103,   266,   190,
     367,   596,   597,   598,   599,    99,   600,   336,   383,   384,
     535,   302,   834,   643,   385,   386,   387,   388,   389,   390,
     391,   392,   393,   709,  1010,   207,   733,  1019,   608,   394,
     806,   395,   211,   520,   383,   384,   698,   646,   911,   302,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   734,
     986,   215,  1021,  1022,   807,   394,   978,   395,   170,   520,
     270,  1030,   275,   216,  1001,  1036,   159,   380,   302,   217,
     171,   850,   937,   987,   853,  1023,   381,   723,   172,   680,
     108,   218,   219,   452,   735,   115,   979,   737,   738,   739,
     740,   741,   742,   743,   744,   745,   746,   380,   748,   749,
     750,   751,   752,   101,   753,   382,   381,   682,   326,   101,
     264,   869,   116,   681,   680,   684,   101,   118,   271,   500,
     276,   177,   683,  1021,  1022,   774,   583,   584,   572,   710,
     101,   723,   160,   101,   270,   382,   117,   724,   595,   275,
     681,   711,   682,   506,   866,  1037,  1038,   712,   511,   683,
     684,   685,   686,   573,   101,   101,   801,   128,   272,   101,
     277,   725,   101,   383,   384,   726,   727,   867,   265,   385,
     386,   387,   388,   389,   390,   391,   392,   393,   685,   826,
     960,   123,   101,   640,   394,   133,   395,   686,   520,   675,
     687,   724,   271,   383,   384,   134,   769,   276,   572,   385,
     386,   387,   388,   389,   390,   391,   392,   393,   266,   688,
     119,   120,   868,   775,   394,   725,   395,   136,   520,   726,
     727,   166,   261,   578,   930,   687,   932,   261,   715,   261,
     137,   873,   272,   723,   124,   125,   688,   277,   602,   142,
     143,   602,   145,   146,   149,   947,   151,   375,   153,   158,
     223,   874,   611,   231,   232,   241,   293,   613,   261,   163,
     222,   247,   797,   603,   167,   248,   604,   294,   261,   249,
     101,   250,   173,   251,   252,   255,   256,   612,   257,   262,
     871,   615,   614,   295,   262,   101,   262,  1025,   101,   637,
     168,   101,   641,   724,   706,   706,   706,   261,   707,   707,
     707,   176,   368,  1039,   376,   261,   616,   372,   970,   373,
     972,   973,   435,   177,   980,   262,   439,   725,   443,   889,
     101,   726,   727,   261,   676,   262,   261,   377,   779,   780,
     781,   782,   783,   784,   785,   786,   787,   788,   414,   790,
     791,   792,   793,   794,   994,   419,   996,   261,   426,   420,
     799,   261,   706,   180,   262,   912,   707,  1002,   400,   261,
     812,  1015,   262,   428,   261,   813,   618,   261,   817,   436,
    1011,   918,   380,   440,  1014,   444,   861,   432,   261,   423,
     262,   381,   280,   262,   875,   433,   261,   924,    93,  1029,
     111,   619,   437,   261,    88,    89,   441,    94,   445,   112,
     920,   862,   770,   446,   262,   261,   827,   776,   262,   876,
     382,   629,   925,   284,   261,   286,   262,    91,    92,  1046,
     630,   262,   953,   954,   262,   657,  1050,   829,   763,   287,
     460,   840,   798,   800,   658,   262,   461,   764,   957,   888,
     109,   110,   322,   262,   892,   814,   843,   904,   818,   857,
     262,   859,  1033,   374,  1034,   844,   369,   174,   905,   958,
     378,   854,   262,   855,   961,   175,   913,   379,   383,   384,
     416,   262,   434,   914,   385,   386,   387,   388,   389,   390,
     391,   392,   393,   447,   974,   962,   113,   114,   454,   394,
     465,   395,   156,   520,   963,   448,   982,   438,   983,     1,
       2,     3,   550,   551,   552,   553,     4,   554,   457,   463,
       5,     6,     7,   442,     8,   464,     9,    10,    11,    12,
     596,   597,   598,   599,   459,   600,   747,   802,   458,    13,
     462,   467,    14,   803,   468,   550,   551,   552,   553,  1009,
     554,   138,   139,   469,   789,    15,   550,   551,   552,   553,
     470,   554,   140,   141,   471,   472,    16,    17,    18,   877,
     473,   474,    19,   550,   551,   552,   553,   475,   554,   476,
     878,  1035,    20,    21,    22,   552,   553,    23,   554,    24,
      25,    26,    27,    28,   477,   478,  1045,   879,    29,    30,
     479,   480,   481,    31,    32,    33,    34,   550,   551,   552,
     553,   938,   554,    35,    36,   482,    37,   483,   598,   599,
      38,   600,   484,    39,    40,    41,    42,   485,   486,   487,
     488,   880,   181,   489,   490,   491,   497,   182,   183,   492,
     493,   496,   184,   498,   550,   551,   552,   553,   502,   554,
     499,   503,   504,   380,   505,   509,   185,   186,   187,   508,
     188,   510,   381,   189,   513,   517,   518,   919,   881,   921,
     190,   191,   192,   519,   527,   193,   194,   529,   195,   196,
     536,   197,   198,   199,   200,   201,   202,   203,   204,   205,
     206,   382,   528,   305,   532,   533,   207,   539,   208,   209,
     540,   210,   306,   211,   550,   551,   552,   553,   541,   554,
     542,   212,   543,   550,   551,   552,   553,   559,   554,   213,
     214,   544,   215,   560,   561,   545,   546,   305,   882,   562,
     563,   307,   564,   547,   216,   163,   306,   883,   548,   556,
     217,   596,   597,   598,   599,   305,   600,   565,   566,   383,
     384,   567,   218,   219,   306,   385,   386,   387,   388,   389,
     390,   391,   392,   393,   568,   307,   569,   570,   571,   380,
     394,   579,   395,   580,   520,   617,   581,   631,   381,   550,
     551,   552,   553,   307,   554,   582,   585,   586,   574,   308,
     309,   587,   588,   620,   589,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   590,   621,   622,   382,   591,   592,
     319,   593,   320,   594,   550,   551,   552,   553,   623,   554,
     624,   625,   773,   308,   309,   626,   628,   627,   632,   310,
     311,   312,   313,   314,   315,   316,   317,   318,   884,   634,
     638,   308,   309,   635,   319,   642,   320,   310,   311,   312,
     313,   314,   315,   316,   317,   318,   644,   646,   639,   645,
     649,   650,   319,   651,   320,   383,   384,   652,   653,   654,
     655,   385,   386,   387,   388,   389,   390,   391,   392,   393,
       1,     2,     3,   656,   659,   660,   394,     4,   395,   661,
     662,     5,     6,     7,   663,     8,   668,     9,    10,    11,
      12,   550,   551,   552,   553,   664,   554,   669,   667,   670,
      13,   671,   672,    14,   703,   677,   550,   551,   552,   553,
     678,   554,   679,   692,   694,   885,    15,   550,   551,   552,
     553,   325,   554,   695,   697,   702,   718,    16,    17,    18,
     886,   704,   326,    19,   327,   328,   554,   754,   719,   756,
     731,   887,   755,    20,    21,    22,   757,   758,    23,   759,
      24,    25,    26,    27,    28,   761,   329,   330,   762,    29,
      30,   190,   771,   765,    31,    32,    33,    34,   289,   766,
     767,   795,   809,   777,    35,    36,   811,    37,   815,   600,
     819,    38,   820,   822,    39,    40,    41,    42,   821,   824,
     832,   331,   325,   332,   839,   333,   596,   597,   598,   599,
     450,   600,   825,   326,   335,   327,   328,   841,   336,   550,
     551,   552,   553,   842,   554,   836,   337,   338,   339,   854,
     894,   828,   340,   341,   342,   830,   163,   329,   330,   845,
     846,   555,   190,   427,   325,   596,   597,   598,   599,   289,
     600,   596,   597,   598,   599,   326,   600,   327,   328,   596,
     597,   598,   599,   847,   600,   550,   551,   552,   553,   895,
     554,   848,   331,   849,   332,   896,   333,   851,   855,   329,
     330,   856,   334,   897,   190,   335,   858,   772,   860,   336,
     852,   289,   863,   870,   872,    -1,   890,   337,   338,   339,
     893,   908,   915,   340,   341,   342,   916,   163,   596,   597,
     598,   599,   922,   600,   331,   923,   332,   926,   333,   596,
     597,   598,   599,   931,   600,   933,   935,   335,   936,   934,
     940,   336,   898,   596,   597,   598,   599,   948,   600,   337,
     338,   339,   950,   899,   951,   340,   341,   342,   964,   163,
     596,   597,   598,   599,   955,   600,   967,   900,   596,   597,
     598,   599,   968,   600,   596,   597,   598,   599,   965,   600,
     550,   551,   552,   553,   901,   554,   550,   551,   552,   553,
     969,   554,   902,   971,   550,   551,   552,   553,   903,   554,
     550,   551,   552,   553,   906,   554,   550,   551,   552,   553,
     959,   554,   550,   551,   552,   553,   981,   554,   990,   976,
     550,   551,   552,   553,  1044,   554,   550,   551,   552,   553,
    1048,   554,   991,   984,   831,   596,   597,   598,   599,   993,
     600,   997,   891,  1005,   550,   551,   552,   553,   917,   554,
     550,   551,   552,   553,   998,   554,   999,   939,   596,   597,
     598,   599,  1007,   600,  1006,  1008,   949,  1017,   550,   551,
     552,   553,   956,   554,   550,   551,   552,   553,  1012,   554,
     975,   550,   551,   552,   553,  1013,   554,  1016,  1018,  1028,
     988,  1031,   550,   551,   552,   553,   989,   554,  1032,   550,
     551,   552,   553,   992,   554,   550,   551,   552,   553,  1041,
     554,  1042,  1043,  1049,  1000,   550,   551,   552,   553,   157,
     554,  1003,   515,   941,   415,   674,   633,  1004,   810,   907,
    1040,   808,   501,   507,   512,   413,   558,  1027,   796,   966,
     768,   371,   453,   722,   946,   816,   292,     0,   835,     0,
       0,   666
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        21,    22,   178,    24,    25,   109,   308,   309,    92,   103,
      31,   582,   103,   103,   103,    36,    37,    38,   320,   103,
     103,    42,   521,   139,   103,   132,   521,   139,   103,   502,
     103,   504,   103,   103,   103,   103,   509,   103,   143,   133,
     134,   146,   525,   628,   729,   847,    23,    41,   303,    32,
      12,    21,    12,    41,    70,   968,    70,   141,    41,    32,
      50,    21,    32,    70,    21,     4,    70,   369,    41,    32,
      70,    41,    23,    41,    70,    21,    21,    21,    41,   127,
      12,   127,    15,    21,    21,    32,   150,    70,    32,    21,
      46,    41,   565,   127,    41,   150,   127,    41,   127,    41,
      70,   863,   864,   865,   152,   151,    32,   283,    70,    12,
      23,     5,    68,    70,    41,    41,    10,   150,   152,   109,
     127,   152,   150,   152,    70,    70,    70,   127,    41,   119,
    1043,    70,    70,    70,   150,    29,    30,   151,   150,    33,
     128,   129,    36,   150,    70,   128,   150,   151,   110,    43,
     150,   128,   129,   130,   131,   151,   133,    90,   128,   129,
     415,   121,   156,   465,   134,   135,   136,   137,   138,   139,
     140,   141,   142,   150,   157,    69,   127,    32,   151,   149,
     127,   151,    76,   153,   128,   129,    41,   155,   151,   121,
     134,   135,   136,   137,   138,   139,   140,   141,   142,   150,
     127,    95,   128,   129,   151,   149,    21,   151,   150,   153,
      12,  1013,    12,   107,   976,    70,     4,    32,   121,   113,
     150,   694,    21,   150,   697,   151,    41,    70,   158,     6,
     150,   125,   126,    21,   536,   151,    21,   539,   540,   541,
     542,   543,   544,   545,   546,   547,   548,    32,   550,   551,
     552,   553,   554,    70,   556,    70,    41,    34,    14,    70,
      12,   946,   150,    32,     6,    42,    70,   150,    70,    21,
      70,    70,    41,   128,   129,   577,   383,   384,   127,    35,
      70,    70,    70,    70,    12,    70,    23,   130,   395,    12,
      32,    47,    34,    21,   127,   150,   151,    53,    21,    41,
      42,    70,    79,   152,    70,    70,   608,    70,   110,    70,
     110,   154,    70,   128,   129,   158,   159,   150,    70,   134,
     135,   136,   137,   138,   139,   140,   141,   142,    70,   631,
     915,   150,    70,   150,   149,   151,   151,    79,   153,   150,
     117,   130,    70,   128,   129,   151,   150,    70,   127,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   110,   128,
     150,   151,   151,   150,   149,   154,   151,   150,   153,   158,
     159,   455,    70,   152,   847,   117,   849,    70,   861,    70,
      23,    90,   110,    70,   150,   151,   128,   110,   127,   150,
     151,   127,   150,   151,    70,   868,    70,    70,    70,    13,
     494,   110,   127,   494,   494,   494,    67,   127,    70,   108,
     494,   494,   150,   152,    32,   494,   152,    78,    70,   494,
      70,   494,    70,   494,   494,   494,   494,   152,   494,   127,
     732,   127,   152,    94,   127,    70,   127,  1008,    70,   460,
     150,    70,   463,   130,   943,   944,   945,    70,   943,   944,
     945,    70,   150,  1024,   127,    70,   152,   150,   931,   150,
     933,   934,    70,    70,   151,   127,    70,   154,    70,   771,
      70,   158,   159,    70,   495,   127,    70,   150,   585,   586,
     587,   588,   589,   590,   591,   592,   593,   594,   150,   596,
     597,   598,   599,   600,   967,   611,   969,    70,   150,   611,
     150,    70,  1001,    23,   127,   807,  1001,   980,   602,    70,
     615,    21,   127,   618,    70,   150,   127,    70,   150,   127,
     993,   150,    32,   127,   997,   127,   127,   150,    70,   613,
     127,    41,    70,   127,   127,   150,    70,   127,    23,  1012,
      23,   152,   150,    70,   150,   151,   150,    32,   150,    32,
     150,   152,   573,   150,   127,    70,   150,   578,   127,   152,
      70,    32,   152,    28,    70,    32,   127,   150,   151,  1042,
      41,   127,   874,   875,   127,    32,  1049,   150,    32,   150,
     152,   150,   603,   604,    41,   127,   158,    41,   890,   150,
     150,   151,    70,   127,   150,   616,    32,   150,   619,   703,
     127,   708,    32,   116,    32,    41,    70,   150,   150,   911,
      41,    41,   127,    41,   916,   158,   150,   150,   128,   129,
     155,   127,   116,   150,   134,   135,   136,   137,   138,   139,
     140,   141,   142,   152,   936,   150,   150,   151,    23,   149,
      23,   151,     0,   153,   150,    70,   948,   116,   950,     7,
       8,     9,   128,   129,   130,   131,    14,   133,   150,   152,
      18,    19,    20,   116,    22,   158,    24,    25,    26,    27,
     128,   129,   130,   131,    70,   133,   152,   150,   150,    37,
      70,   857,    40,   156,    70,   128,   129,   130,   131,   991,
     133,   150,   151,    23,   152,    53,   128,   129,   130,   131,
      23,   133,   150,   151,    23,    23,    64,    65,    66,   152,
      23,    23,    70,   128,   129,   130,   131,    23,   133,    23,
     152,  1023,    80,    81,    82,   130,   131,    85,   133,    87,
      88,    89,    90,    91,    23,    23,  1038,   152,    96,    97,
      23,    23,    23,   101,   102,   103,   104,   128,   129,   130,
     131,   858,   133,   111,   112,    23,   114,    23,   130,   131,
     118,   133,    23,   121,   122,   123,   124,    23,    23,    23,
      23,   152,     5,    23,    23,    23,    70,    10,    11,    23,
      23,    23,    15,    70,   128,   129,   130,   131,   127,   133,
      70,    70,   127,    32,    70,   127,    29,    30,    31,    70,
      33,    70,    41,    36,   151,    23,   152,   828,   152,   830,
      43,    44,    45,   150,   152,    48,    49,    23,    51,    52,
     151,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      63,    70,   152,    32,    70,    70,    69,   151,    71,    72,
     151,    74,    41,    76,   128,   129,   130,   131,   151,   133,
     151,    84,   151,   128,   129,   130,   131,    23,   133,    92,
      93,   151,    95,    23,    23,   151,   151,    32,   152,    23,
      23,    70,    23,   151,   107,   108,    41,   152,   151,   151,
     113,   128,   129,   130,   131,    32,   133,    23,    23,   128,
     129,    23,   125,   126,    41,   134,   135,   136,   137,   138,
     139,   140,   141,   142,    23,    70,    23,    23,    23,    32,
     149,   116,   151,    70,   153,    23,   150,    23,    41,   128,
     129,   130,   131,    70,   133,   151,   151,   151,   127,   128,
     129,   151,   151,   116,   151,   134,   135,   136,   137,   138,
     139,   140,   141,   142,   151,    70,   116,    70,   151,   151,
     149,   151,   151,   151,   128,   129,   130,   131,    70,   133,
     116,    70,   127,   128,   129,   150,    70,   151,    41,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   152,   150,
      70,   128,   129,   150,   149,    70,   151,   134,   135,   136,
     137,   138,   139,   140,   141,   142,    32,   155,   150,    70,
      41,    41,   149,    41,   151,   128,   129,    41,    41,    41,
      41,   134,   135,   136,   137,   138,   139,   140,   141,   142,
       7,     8,     9,    32,    32,    41,   149,    14,   151,    41,
      41,    18,    19,    20,    41,    22,    41,    24,    25,    26,
      27,   128,   129,   130,   131,    70,   133,    41,   151,    32,
      37,    70,    70,    40,   150,    70,   128,   129,   130,   131,
      23,   133,   127,   127,   127,   152,    53,   128,   129,   130,
     131,     3,   133,   127,   127,    70,   150,    64,    65,    66,
     152,    70,    14,    70,    16,    17,   133,    41,   150,    41,
     127,   152,    32,    80,    81,    82,    41,    41,    85,    41,
      87,    88,    89,    90,    91,    41,    38,    39,    41,    96,
      97,    43,    70,    41,   101,   102,   103,   104,    50,    41,
      41,   150,   150,   116,   111,   112,   150,   114,    32,   133,
     116,   118,   116,    41,   121,   122,   123,   124,   116,   152,
      41,    73,     3,    75,   154,    77,   128,   129,   130,   131,
     151,   133,   152,    14,    86,    16,    17,    70,    90,   128,
     129,   130,   131,    70,   133,   157,    98,    99,   100,    41,
     152,   152,   104,   105,   106,   152,   108,    38,    39,   150,
     127,   150,    43,   115,     3,   128,   129,   130,   131,    50,
     133,   128,   129,   130,   131,    14,   133,    16,    17,   128,
     129,   130,   131,   127,   133,   128,   129,   130,   131,   152,
     133,    70,    73,   127,    75,   152,    77,    70,    41,    38,
      39,   152,    83,   152,    43,    86,    23,   150,    23,    90,
     150,    50,   150,    70,    70,   133,    70,    98,    99,   100,
     152,   150,   152,   104,   105,   106,    23,   108,   128,   129,
     130,   131,   154,   133,    73,    70,    75,    70,    77,   128,
     129,   130,   131,   127,   133,   127,   150,    86,    23,   127,
      32,    90,   152,   128,   129,   130,   131,    23,   133,    98,
      99,   100,    23,   152,    41,   104,   105,   106,    41,   108,
     128,   129,   130,   131,   150,   133,   127,   152,   128,   129,
     130,   131,   127,   133,   128,   129,   130,   131,   157,   133,
     128,   129,   130,   131,   152,   133,   128,   129,   130,   131,
     127,   133,   152,   150,   128,   129,   130,   131,   152,   133,
     128,   129,   130,   131,   152,   133,   128,   129,   130,   131,
     152,   133,   128,   129,   130,   131,   152,   133,   152,   150,
     128,   129,   130,   131,   152,   133,   128,   129,   130,   131,
     152,   133,    23,   154,   150,   128,   129,   130,   131,   127,
     133,   127,   150,    41,   128,   129,   130,   131,   150,   133,
     128,   129,   130,   131,   150,   133,   150,   150,   128,   129,
     130,   131,    41,   133,   154,   120,   150,    41,   128,   129,
     130,   131,   150,   133,   128,   129,   130,   131,   127,   133,
     150,   128,   129,   130,   131,   127,   133,   152,   154,    70,
     150,   150,   128,   129,   130,   131,   150,   133,    41,   128,
     129,   130,   131,   150,   133,   128,   129,   130,   131,   157,
     133,   127,   127,   127,   150,   128,   129,   130,   131,    43,
     133,   150,   281,   861,   136,   494,   455,   150,   613,   802,
    1024,   611,   267,   273,   278,   134,   323,   150,   602,   924,
     572,   125,   161,   530,   866,   618,   114,    -1,   647,    -1,
      -1,   487
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,     9,    14,    18,    19,    20,    22,    24,
      25,    26,    27,    37,    40,    53,    64,    65,    66,    70,
      80,    81,    82,    85,    87,    88,    89,    90,    91,    96,
      97,   101,   102,   103,   104,   111,   112,   114,   118,   121,
     122,   123,   124,   161,   162,   163,   164,   165,   166,   167,
     168,   169,   174,   175,   176,   177,   180,   182,   185,   190,
     201,   202,   207,   210,   213,   216,   219,   225,   231,   234,
     239,   242,   245,   248,   249,   252,   253,   255,   256,   257,
     260,   261,   265,   266,   267,   268,   274,   277,   150,   151,
     150,   150,   151,    23,    32,    41,    70,   151,    70,   151,
     150,    70,   150,   151,   222,   150,   150,   150,   150,   150,
     151,    23,    32,   150,   151,   151,   150,    23,   150,   150,
     151,   222,   222,   150,   150,   151,   222,   222,    70,   173,
      23,    41,   275,   151,   151,   222,   150,    23,   150,   151,
     150,   151,   150,   151,   222,   150,   151,   222,   222,    70,
     170,    70,   171,    70,   172,   222,     0,   162,    13,     4,
      70,   262,   263,   108,   214,   215,   281,    32,   150,    41,
     150,   150,   158,    70,   150,   158,    70,    70,   183,   184,
      23,     5,    10,    11,    15,    29,    30,    31,    33,    36,
      43,    44,    45,    48,    49,    51,    52,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,    69,    71,    72,
      74,    76,    84,    92,    93,    95,   107,   113,   125,   126,
     240,   241,   281,   291,   300,   301,   302,   303,   304,   305,
     306,   307,   308,   309,   310,   311,   312,   313,   314,   315,
     316,   317,   318,   319,   320,   321,   322,   323,   324,   325,
     326,   327,   328,   330,   331,   336,   337,   338,   339,   340,
     341,    70,   127,   150,    12,    70,   110,   226,   227,   228,
      12,    70,   110,   235,   236,    12,    70,   110,   232,   233,
      70,   186,   187,   183,    28,   181,    32,   150,   191,    50,
     109,   119,   283,    67,    78,    94,   269,   270,   333,   334,
     335,    12,   121,   203,   204,    32,    41,    70,   128,   129,
     134,   135,   136,   137,   138,   139,   140,   141,   142,   149,
     151,   178,    70,   250,   251,     3,    14,    16,    17,    38,
      39,    73,    75,    77,    83,    86,    90,    98,    99,   100,
     104,   105,   106,   221,   258,   259,   280,   281,   282,   283,
     284,   285,   286,   287,   288,   289,   290,   291,   292,   293,
     294,   295,   297,   298,   299,   329,   332,   150,   150,    70,
     254,   259,   150,   150,   116,    70,   127,   150,    41,   150,
      32,    41,    70,   128,   129,   134,   135,   136,   137,   138,
     139,   140,   141,   142,   149,   151,   198,   200,   243,   244,
     291,   307,   308,   317,   323,   324,   325,   326,   327,   328,
     336,   337,   338,   243,   150,   203,   155,   217,   218,   294,
     300,   211,   212,   281,   220,   221,   150,   115,   221,   278,
     279,   342,   150,   150,   116,    70,   127,   150,   116,    70,
     127,   150,   116,    70,   127,   150,   150,   152,    70,   127,
     151,   264,    21,   263,    23,   127,   152,   150,   150,    70,
     152,   158,    70,   152,   158,    23,    21,   184,    70,    23,
      23,    23,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,   127,   152,    23,    70,    70,    70,
      21,   227,   127,    70,   127,    70,    21,   236,    70,   127,
      70,    21,   233,   151,    21,   187,    21,    23,   152,   150,
     153,   196,   197,   198,   199,   127,   152,   152,   152,    23,
     127,   152,    70,    70,    21,   204,   151,   178,   178,   151,
     151,   151,   151,   151,   151,   151,   151,   151,   151,   178,
     128,   129,   130,   131,   133,   150,   151,    21,   251,    23,
      23,    23,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,   127,   152,   127,   178,    21,    70,   152,   116,
      70,   150,   151,   198,   198,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   198,   128,   129,   130,   131,
     133,   276,   127,   152,   152,    21,    32,    41,   151,   208,
     209,   127,   152,   127,   152,   127,   152,    23,   127,   152,
     116,    70,   116,    70,   116,    70,   150,   151,    70,    32,
      41,    23,    41,   215,   150,   150,   150,   222,    70,   150,
     150,   222,    70,   178,    32,    70,   155,   345,   346,    41,
      41,    41,    41,    41,    41,    41,    32,    32,    41,    32,
      41,    41,    41,    41,    70,    41,   346,   151,    41,    41,
      32,    70,    70,   343,   241,   150,   222,    70,    23,   127,
       6,    32,    34,    41,    42,    70,    79,   117,   128,   229,
     237,   238,   127,   238,   127,   127,   238,   127,    41,   128,
     129,   223,    70,   150,    70,    21,   197,   199,    23,   150,
      35,    47,    53,   188,   189,   295,   296,   195,   150,   150,
      46,    68,   270,    70,   130,   154,   158,   159,   271,   272,
     273,   127,    23,   127,   150,   178,   179,   178,   178,   178,
     178,   178,   178,   178,   178,   178,   178,   152,   178,   178,
     178,   178,   178,   178,    41,    32,    41,    41,    41,    41,
     238,    41,    41,    32,    41,    41,    41,    41,   258,   150,
     222,    70,   150,   127,   178,   150,   222,   116,   223,   198,
     198,   198,   198,   198,   198,   198,   198,   198,   198,   152,
     198,   198,   198,   198,   198,   150,   244,   150,   222,   150,
     222,   178,   150,   156,    32,    41,   127,   151,   218,   150,
     212,   150,   221,   150,   222,    32,   279,   150,   222,   116,
     116,   116,    41,   264,   152,   152,   178,   150,   152,   150,
     152,   150,    41,   344,   156,   344,   157,   246,   247,   154,
     150,    70,    70,    32,    41,   150,   127,   127,    70,   127,
     238,    70,   150,   238,    41,    41,   152,   183,    23,   198,
      23,   127,   152,   150,   193,   192,   127,   150,   151,   273,
      70,   178,    70,    90,   110,   127,   152,   152,   152,   152,
     152,   152,   152,   152,   152,   152,   152,   152,   150,   178,
      70,   150,   150,   152,   152,   152,   152,   152,   152,   152,
     152,   152,   152,   152,   150,   150,   152,   209,   150,    32,
      41,   151,   178,   150,   150,   152,    23,   150,   150,   222,
     150,   222,   154,    70,   127,   152,    70,   230,   238,   237,
     238,   127,   238,   127,   127,   150,    23,    21,   198,   150,
      32,   189,   194,   196,   196,   196,   272,   238,    23,   150,
      23,    41,   205,   178,   178,   150,   150,   178,   178,   152,
     264,   178,   150,   150,    41,   157,   246,   127,   127,   127,
     238,   150,   238,   238,   178,   150,   150,    21,    21,    21,
     151,   152,   178,   178,   154,    41,   127,   150,   150,   150,
     152,    23,   150,   127,   238,   230,   238,   127,   150,   150,
     150,   196,   238,   150,   150,    41,   154,    41,   120,   178,
     157,   238,   127,   127,   238,    21,   152,    41,   154,    32,
      70,   128,   129,   151,   206,   223,   224,   150,    70,   238,
     237,   150,    41,    32,    32,   178,    70,   150,   151,   223,
     224,   157,   127,   127,   152,   178,   238,   230,   152,   127,
     238
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
      59,    40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   160,   161,   161,   162,   162,   162,   162,   162,   162,
     162,   162,   162,   162,   162,   162,   162,   162,   162,   162,
     162,   162,   162,   162,   162,   162,   162,   162,   162,   162,
     162,   162,   162,   162,   162,   162,   162,   162,   162,   162,
     162,   162,   162,   163,   163,   163,   163,   164,   164,   165,
     166,   167,   168,   169,   170,   170,   170,   170,   170,   170,
     171,   171,   171,   171,   171,   171,   172,   172,   172,   172,
     172,   172,   173,   173,   173,   173,   173,   173,   174,   174,
     175,   175,   176,   176,   177,   178,   178,   178,   178,   178,
     178,   178,   178,   178,   178,   178,   178,   178,   178,   178,
     178,   178,   178,   178,   178,   178,   178,   179,   179,   180,
     180,   181,   182,   183,   183,   184,   185,   186,   186,   187,
     188,   188,   189,   189,   189,   189,   191,   190,   192,   190,
     193,   190,   194,   190,   195,   190,   196,   196,   196,   196,
     197,   197,   198,   198,   198,   198,   198,   198,   198,   198,
     198,   198,   198,   198,   198,   198,   198,   198,   198,   198,
     198,   198,   198,   199,   200,   200,   201,   202,   203,   203,
     204,   204,   204,   204,   204,   205,   205,   205,   205,   205,
     205,   206,   206,   206,   206,   206,   206,   206,   206,   207,
     208,   208,   209,   209,   209,   209,   209,   209,   209,   209,
     209,   210,   210,   211,   211,   212,   213,   213,   214,   214,
     215,   216,   216,   217,   217,   218,   218,   219,   219,   219,
     219,   220,   220,   221,   221,   221,   221,   221,   221,   221,
     221,   221,   221,   221,   221,   221,   221,   221,   221,   221,
     221,   221,   221,   222,   222,   222,   222,   222,   222,   223,
     223,   223,   224,   224,   224,   225,   226,   226,   227,   228,
     228,   228,   229,   229,   229,   229,   229,   230,   230,   230,
     230,   231,   232,   232,   233,   233,   233,   234,   235,   235,
     236,   236,   236,   237,   237,   237,   237,   237,   238,   238,
     238,   238,   238,   238,   239,   239,   239,   239,   240,   240,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   242,
     242,   243,   243,   244,   244,   244,   244,   244,   244,   244,
     244,   244,   244,   244,   244,   244,   245,   245,   246,   246,
     247,   247,   248,   249,   250,   250,   251,   252,   253,   254,
     254,   254,   254,   255,   256,   256,   256,   256,   257,   257,
     257,   257,   258,   258,   259,   259,   260,   261,   262,   262,
     263,   263,   263,   264,   264,   264,   265,   265,   266,   266,
     266,   266,   266,   266,   267,   267,   267,   267,   267,   267,
     268,   269,   269,   270,   270,   270,   271,   271,   271,   271,
     272,   272,   273,   273,   273,   273,   273,   275,   276,   274,
     277,   277,   277,   277,   278,   278,   279,   279,   280,   281,
     282,   283,   284,   285,   286,   287,   288,   289,   290,   291,
     292,   293,   294,   295,   296,   297,   298,   299,   299,   300,
     301,   301,   302,   303,   304,   305,   306,   307,   307,   308,
     309,   310,   311,   312,   313,   313,   314,   315,   316,   317,
     318,   319,   320,   321,   322,   323,   324,   325,   326,   327,
     328,   329,   330,   331,   332,   333,   333,   334,   335,   336,
     337,   338,   339,   340,   341,   342,   343,   344,   344,   345,
     345,   346
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
       6,     3,     5,     3,     2,     5,     3,     6,     2,     5,
       3,     6,     1,     1,     1,     3,     3,     4,     2,     1,
       5,     7,     9,     0,     3,     3,     2,     5,     5,     6,
       3,     7,     8,     5,     5,     6,     3,     7,     8,     5,
       6,     3,     1,     1,     1,     1,     1,     3,     4,     6,
       1,     2,     1,     1,     1,     1,     1,     0,     0,     5,
       2,     5,     3,     6,     3,     1,     1,     1,     3,     3,
       3,     1,     3,     3,     3,     3,     1,     1,     1,     3,
       3,     3,     3,     3,     3,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     3,
       3,     3,     3,     5,     3,     3,     3,     3,     1,     3,
       3,     3,     1,     1,     1,     1,     1,     3,     1,     1,
       1,     1,     3,     3,     3,     3,     3,     1,     1,     3,
       3,     3,     1,     1,     1,     3,     3,     1,     3,     2,
       2,     2
  };

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
  /* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
     First, the terminals, then, starting at \a yyntokens_, nonterminals.  */
  const char*
  const parser::yytname_[] =
  {
    "$end", "error", "$undefined", "AR", "AUTOCORR", "BAYESIAN_IRF",
  "BETA_PDF", "CALIB", "CALIB_VAR", "CHECK", "CONF_SIG", "CONSTANT",
  "CORR", "COVAR", "CUTOFF", "DATAFILE", "DR_ALGO", "DROP", "DSAMPLE",
  "DYNASAVE", "DYNATYPE", "END", "ENDVAL", "EQUAL", "ESTIMATION",
  "ESTIMATED_PARAMS", "ESTIMATED_PARAMS_BOUNDS", "ESTIMATED_PARAMS_INIT",
  "FILENAME", "FILTER_STEP_AHEAD", "FILTERED_VARS", "FIRST_OBS",
  "FLOAT_NUMBER", "FORECAST", "GAMMA_PDF", "GCC_COMPILER", "GRAPH",
  "HISTVAL", "HP_FILTER", "HP_NGRID", "INITVAL", "INT_NUMBER",
  "INV_GAMMA_PDF", "IRF", "KALMAN_ALGO", "KALMAN_TOL", "LAPLACE",
  "LCC_COMPILER", "LIK_ALGO", "LIK_INIT", "LINEAR", "LOAD_MH_FILE",
  "LOGLINEAR", "MARKOWITZ", "MH_DROP", "MH_INIT_SCALE", "MH_JSCALE",
  "MH_MODE", "MH_NBLOCKS", "MH_REPLIC", "MH_RECOVER", "MODE_CHECK",
  "MODE_COMPUTE", "MODE_FILE", "MODEL", "MODEL_COMPARISON", "MSHOCKS",
  "MODEL_COMPARISON_APPROXIMATION", "MODIFIEDHARMONICMEAN",
  "MOMENTS_VARENDO", "NAME", "NOBS", "NOCONSTANT", "NOCORR",
  "NODIAGNOSTIC", "NOFUNCTIONS", "NOGRAPH", "NOMOMENTS", "NOPRINT",
  "NORMAL_PDF", "OBSERVATION_TRENDS", "OLR", "OLR_INST", "OLR_BETA",
  "OPTIM", "OPTIM_WEIGHTS", "ORDER", "OSR", "OSR_PARAMS", "PARAMETERS",
  "PERIODS", "PLANNER_OBJECTIVE", "PREFILTER", "PRESAMPLE", "PRINT",
  "PRIOR_TRUNC", "PRIOR_ANALYSIS", "POSTERIOR_ANALYSIS", "QZ_CRITERIUM",
  "RELATIVE_IRF", "REPLIC", "RPLOT", "SHOCKS", "SIGMA_E", "SIMUL",
  "SIMUL_ALGO", "SIMUL_SEED", "SMOOTHER", "SOLVE_ALGO", "SPARSE_DLL",
  "STDERR", "STEADY", "STOCH_SIMUL", "TEX", "RAMSEY_POLICY",
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
  "ramsey_policy_options_list", "ramsey_policy_options", "o_dr_algo",
  "o_solve_algo", "o_simul_algo", "o_linear", "o_order", "o_replic",
  "o_drop", "o_ar", "o_nocorr", "o_nofunctions", "o_nomoments", "o_irf",
  "o_hp_filter", "o_hp_ngrid", "o_periods", "o_cutoff", "o_markowitz",
  "o_simul", "o_simul_seed", "o_qz_criterium", "o_datafile", "o_nobs",
  "o_first_obs", "o_prefilter", "o_presample", "o_lik_algo", "o_lik_init",
  "o_nograph", "o_conf_sig", "o_mh_replic", "o_mh_drop", "o_mh_jscale",
  "o_optim", "o_mh_init_scale", "o_mode_file", "o_mode_compute",
  "o_mode_check", "o_prior_trunc", "o_mh_mode", "o_mh_nblcks",
  "o_load_mh_file", "o_loglinear", "o_nodiagnostic", "o_bayesian_irf",
  "o_tex", "o_forecast", "o_smoother", "o_moments_varendo",
  "o_filtered_vars", "o_relative_irf", "o_kalman_algo", "o_kalman_tol",
  "o_olr_beta", "o_model_comparison_approximation", "o_print", "o_noprint",
  "o_xls_sheet", "o_xls_range", "o_filter_step_ahead", "o_constant",
  "o_noconstant", "o_mh_recover", "o_planner_discount", "range",
  "vec_int_elem", "vec_int_1", "vec_int", 0
  };
#endif

#if YYDEBUG
  /* YYRHS -- A `-1'-separated list of the rules' RHS.  */
  const parser::rhs_number_type
  parser::yyrhs_[] =
  {
       161,     0,    -1,   162,    -1,   161,   162,    -1,   163,    -1,
     174,    -1,   175,    -1,   176,    -1,   190,    -1,   180,    -1,
     182,    -1,   185,    -1,   177,    -1,   201,    -1,   202,    -1,
     207,    -1,   210,    -1,   213,    -1,   216,    -1,   219,    -1,
     239,    -1,   242,    -1,   245,    -1,   225,    -1,   234,    -1,
     231,    -1,   248,    -1,   249,    -1,   252,    -1,   164,    -1,
     165,    -1,   253,    -1,   255,    -1,   256,    -1,   261,    -1,
     265,    -1,   266,    -1,   267,    -1,   257,    -1,   260,    -1,
     268,    -1,   274,    -1,   277,    -1,   169,    -1,   166,    -1,
     167,    -1,   168,    -1,    18,    41,   150,    -1,    18,    41,
      41,   150,    -1,   101,   222,   150,    -1,   121,   170,   150,
      -1,   122,   171,   150,    -1,   123,   172,   150,    -1,    89,
     173,   150,    -1,   170,    70,    -1,   170,   127,    70,    -1,
      70,    -1,   170,    70,   116,    -1,   170,   127,    70,   116,
      -1,    70,   116,    -1,   171,    70,    -1,   171,   127,    70,
      -1,    70,    -1,   171,    70,   116,    -1,   171,   127,    70,
     116,    -1,    70,   116,    -1,   172,    70,    -1,   172,   127,
      70,    -1,    70,    -1,   172,    70,   116,    -1,   172,   127,
      70,   116,    -1,    70,   116,    -1,   173,    70,    -1,   173,
     127,    70,    -1,    70,    -1,   173,    70,   116,    -1,   173,
     127,    70,   116,    -1,    70,   116,    -1,    90,    41,   150,
      -1,    90,    23,    41,   150,    -1,    14,    32,   150,    -1,
      14,    23,    32,   150,    -1,    53,    32,   150,    -1,    53,
      23,    32,   150,    -1,    70,    23,   178,   150,    -1,   151,
     178,   152,    -1,    70,    -1,    32,    -1,    41,    -1,   178,
     129,   178,    -1,   178,   128,   178,    -1,   178,   130,   178,
      -1,   178,   131,   178,    -1,   178,   133,   178,    -1,   128,
     178,    -1,   129,   178,    -1,   134,   151,   178,   152,    -1,
     135,   151,   178,   152,    -1,   136,   151,   178,   152,    -1,
     137,   151,   178,   152,    -1,   138,   151,   178,   152,    -1,
     139,   151,   178,   152,    -1,   140,   151,   178,   152,    -1,
     141,   151,   178,   152,    -1,   142,   151,   178,   152,    -1,
     149,   151,   178,   152,    -1,    70,   151,   179,   152,    -1,
     178,    -1,   179,   127,   178,    -1,    40,   150,   183,    21,
      -1,    40,   151,   181,   152,   150,   183,    21,    -1,    28,
      23,    70,    -1,    22,   150,   183,    21,    -1,   183,   184,
      -1,   184,    -1,    70,    23,   178,   150,    -1,    37,   150,
     186,    21,    -1,   186,   187,    -1,   187,    -1,    70,   151,
     223,   152,    23,   178,   150,    -1,   188,   127,   189,    -1,
     189,    -1,    47,    -1,    35,    -1,   295,    -1,   296,    -1,
      -1,    64,   150,   191,   196,    21,    -1,    -1,    64,   151,
     283,   152,   150,   192,   196,    21,    -1,    -1,    64,   151,
     119,   152,   150,   193,   196,    21,    -1,    -1,    64,   151,
     109,   127,   188,   152,   194,   150,   196,    21,    -1,    -1,
      64,   151,   109,   152,   195,   150,   196,    21,    -1,   196,
     197,    -1,   196,   199,    -1,   197,    -1,   199,    -1,   198,
      23,   198,   150,    -1,   198,   150,    -1,   151,   198,   152,
      -1,   200,    -1,    32,    -1,    41,    -1,   198,   129,   198,
      -1,   198,   128,   198,    -1,   198,   130,   198,    -1,   198,
     131,   198,    -1,   198,   133,   198,    -1,   128,   198,    -1,
     129,   198,    -1,   134,   151,   198,   152,    -1,   135,   151,
     198,   152,    -1,   136,   151,   198,   152,    -1,   137,   151,
     198,   152,    -1,   138,   151,   198,   152,    -1,   139,   151,
     198,   152,    -1,   140,   151,   198,   152,    -1,   141,   151,
     198,   152,    -1,   142,   151,   198,   152,    -1,   149,   151,
     198,   152,    -1,   153,    70,    23,   198,   150,    -1,    70,
      -1,    70,   151,   223,   152,    -1,   102,   150,   203,    21,
      -1,    66,   150,   203,    21,    -1,   203,   204,    -1,   204,
      -1,   121,    70,   150,    90,   205,   150,   120,   206,   150,
      -1,   121,    70,   150,   110,   178,   150,    -1,   121,    70,
      23,   178,   150,    -1,   121,    70,   127,    70,    23,   178,
     150,    -1,    12,    70,   127,    70,    23,   178,   150,    -1,
     205,    41,    -1,   205,    41,   154,    41,    -1,   205,   127,
      41,    -1,   205,   127,    41,   154,    41,    -1,    41,   154,
      41,    -1,    41,    -1,   206,   224,    -1,   206,   223,    -1,
     206,    70,    -1,   224,    -1,   223,    -1,    70,    -1,   206,
     151,   178,   152,    -1,   151,   178,   152,    -1,   103,    23,
     155,   208,   156,   150,    -1,   208,   150,   209,    -1,   209,
      -1,   209,   127,   151,   178,   152,    -1,   209,   127,    32,
      -1,   209,   127,    41,    -1,   209,   151,   178,   152,    -1,
     209,    32,    -1,   209,    41,    -1,   151,   178,   152,    -1,
      32,    -1,    41,    -1,   111,   150,    -1,   111,   151,   211,
     152,   150,    -1,   211,   127,   212,    -1,   212,    -1,   281,
      -1,     9,   150,    -1,     9,   151,   214,   152,   150,    -1,
     214,   127,   215,    -1,   215,    -1,   281,    -1,   104,   150,
      -1,   104,   151,   217,   152,   150,    -1,   217,   127,   218,
      -1,   218,    -1,   294,    -1,   300,    -1,   112,   150,    -1,
     112,   151,   220,   152,   150,    -1,   112,   222,   150,    -1,
     112,   151,   220,   152,   222,   150,    -1,   220,   127,   221,
      -1,   221,    -1,   280,    -1,   281,    -1,   282,    -1,   283,
      -1,   284,    -1,   285,    -1,   286,    -1,   287,    -1,   288,
      -1,   289,    -1,   290,    -1,   291,    -1,   329,    -1,   292,
      -1,   293,    -1,   294,    -1,   295,    -1,   297,    -1,   298,
      -1,   299,    -1,   222,    70,    -1,   222,    70,    23,    70,
      -1,   222,   127,    70,    -1,   222,   127,    70,    23,    70,
      -1,    70,    -1,    70,    23,    70,    -1,   129,    41,    -1,
     128,    41,    -1,    41,    -1,   129,    32,    -1,   128,    32,
      -1,    32,    -1,    25,   150,   226,    21,    -1,   226,   227,
      -1,   227,    -1,   228,   127,   229,   150,    -1,   110,    70,
      -1,    70,    -1,    12,    70,   127,    70,    -1,   237,   127,
     230,    -1,   238,   127,   237,   127,   230,    -1,   238,   127,
     238,   127,   238,   127,   237,   127,   230,    -1,   238,    -1,
     238,   127,   238,   127,   238,    -1,   238,   127,   238,    -1,
     238,   127,   238,   127,   238,    -1,   238,   127,   238,   127,
     238,   127,   238,    -1,   238,   127,   238,   127,   238,   127,
     238,   127,   238,    -1,    27,   150,   232,    21,    -1,   232,
     233,    -1,   233,    -1,   110,    70,   127,   238,   150,    -1,
      12,    70,   127,    70,   127,   238,   150,    -1,    70,   127,
     238,   150,    -1,    26,   150,   235,    21,    -1,   235,   236,
      -1,   236,    -1,   110,    70,   127,   238,   127,   238,   150,
      -1,    12,    70,   127,    70,   127,   238,   127,   238,   150,
      -1,    70,   127,   238,   127,   238,   150,    -1,     6,    -1,
      34,    -1,    79,    -1,    42,    -1,   117,    -1,    -1,    41,
      -1,    32,    -1,    70,    -1,   128,    41,    -1,   128,    32,
      -1,    24,   150,    -1,    24,   151,   240,   152,   150,    -1,
      24,   222,   150,    -1,    24,   151,   240,   152,   222,   150,
      -1,   240,   127,   241,    -1,   241,    -1,   300,    -1,   301,
      -1,   302,    -1,   303,    -1,   304,    -1,   305,    -1,   306,
      -1,   307,    -1,   308,    -1,   309,    -1,   310,    -1,   311,
      -1,   312,    -1,   313,    -1,   314,    -1,   315,    -1,   316,
      -1,   317,    -1,   318,    -1,   319,    -1,   320,    -1,   321,
      -1,   322,    -1,   323,    -1,   291,    -1,   324,    -1,   325,
      -1,   326,    -1,   327,    -1,   328,    -1,   330,    -1,   331,
      -1,   336,    -1,   337,    -1,   338,    -1,   281,    -1,   339,
      -1,   340,    -1,   341,    -1,    96,   151,   243,   152,   150,
      -1,    96,   151,   243,   152,   222,   150,    -1,   243,   127,
     244,    -1,   244,    -1,   307,    -1,   308,    -1,   317,    -1,
     323,    -1,   291,    -1,   324,    -1,   325,    -1,   326,    -1,
     327,    -1,   328,    -1,   336,    -1,   337,    -1,   338,    -1,
      97,   151,   243,   152,   150,    -1,    97,   151,   243,   152,
     222,   150,    -1,   157,    70,   157,   127,   157,    70,   157,
      -1,   157,    70,   157,   127,   238,    -1,   246,    -1,   247,
     127,   246,    -1,   124,   222,   150,    -1,    80,   150,   250,
      21,    -1,   250,   251,    -1,   251,    -1,    70,   151,   178,
     152,   150,    -1,   118,   222,   150,    -1,    85,   150,   254,
      21,    -1,   254,    70,   178,   150,    -1,   254,    70,   127,
      70,   178,   150,    -1,    70,   178,   150,    -1,    70,   127,
      70,   178,   150,    -1,    88,   222,   150,    -1,    87,   150,
      -1,    87,   151,   259,   152,   150,    -1,    87,   222,   150,
      -1,    87,   151,   259,   152,   222,   150,    -1,    81,   150,
      -1,    81,   151,   259,   152,   150,    -1,    81,   222,   150,
      -1,    81,   151,   259,   152,   222,   150,    -1,   332,    -1,
     221,    -1,   258,    -1,   259,   127,   258,    -1,    82,   222,
     150,    -1,     8,   150,   262,    21,    -1,   262,   263,    -1,
     263,    -1,    70,   264,    23,   178,   150,    -1,    70,   127,
      70,   264,    23,   178,   150,    -1,     4,    70,   151,    41,
     152,   264,    23,   178,   150,    -1,    -1,   151,    41,   152,
      -1,   151,    32,   152,    -1,     7,   150,    -1,     7,   151,
      13,   152,   150,    -1,    20,   151,    70,   152,   150,    -1,
      20,   151,    70,   152,   222,   150,    -1,    20,    70,   150,
      -1,    20,   151,    70,   158,    70,   152,   150,    -1,    20,
     151,    70,   158,    70,   152,   222,   150,    -1,    20,    70,
     158,    70,   150,    -1,    19,   151,    70,   152,   150,    -1,
      19,   151,    70,   152,   222,   150,    -1,    19,    70,   150,
      -1,    19,   151,    70,   158,    70,   152,   150,    -1,    19,
     151,    70,   158,    70,   152,   222,   150,    -1,    19,    70,
     158,    70,   150,    -1,    65,   151,   269,   152,   271,   150,
      -1,   269,   127,   270,    -1,   270,    -1,   333,    -1,   334,
      -1,   335,    -1,   272,    -1,   271,   127,   272,    -1,   272,
     151,   238,   152,    -1,   271,   127,   272,   151,   238,   152,
      -1,   273,    -1,   272,   273,    -1,    70,    -1,   159,    -1,
     130,    -1,   154,    -1,   158,    -1,    -1,    -1,    91,   275,
     198,   276,   150,    -1,   114,   150,    -1,   114,   151,   278,
     152,   150,    -1,   114,   222,   150,    -1,   114,   151,   278,
     152,   222,   150,    -1,   278,   127,   279,    -1,   279,    -1,
     221,    -1,   342,    -1,    16,    23,    41,    -1,   108,    23,
      41,    -1,   105,    23,    41,    -1,    50,    -1,    86,    23,
      41,    -1,   100,    23,    41,    -1,    17,    23,    41,    -1,
       3,    23,    41,    -1,    73,    -1,    75,    -1,    77,    -1,
      43,    23,    41,    -1,    38,    23,    41,    -1,    39,    23,
      41,    -1,    90,    23,    41,    -1,    14,    23,    32,    -1,
      53,    23,    32,    -1,   104,    -1,   106,    23,    41,    -1,
      98,    23,    41,    -1,    98,    23,    32,    -1,    15,    23,
      70,    -1,    71,    23,   346,    -1,    71,    23,    41,    -1,
      31,    23,    41,    -1,    92,    23,    41,    -1,    93,    23,
      41,    -1,    48,    23,    41,    -1,    49,    23,    41,    -1,
      76,    -1,    36,    -1,    10,    23,    32,    -1,    59,    23,
      41,    -1,    54,    23,    32,    -1,    56,    23,    32,    -1,
      84,    23,   151,   247,   152,    -1,    55,    23,    32,    -1,
      55,    23,    41,    -1,    63,    23,    70,    -1,    62,    23,
      41,    -1,    61,    -1,    95,    23,    32,    -1,    57,    23,
      41,    -1,    58,    23,    41,    -1,    51,    -1,    52,    -1,
      74,    -1,     5,    -1,   113,    -1,    33,    23,    41,    -1,
     107,    -1,    69,    -1,    30,    -1,    99,    -1,    44,    23,
      41,    -1,    45,    23,    41,    -1,    83,    23,   238,    -1,
      67,    23,    46,    -1,    67,    23,    68,    -1,    94,    -1,
      78,    -1,   125,    23,    70,    -1,   126,    23,   343,    -1,
      29,    23,   346,    -1,    11,    -1,    72,    -1,    60,    -1,
     115,    23,    32,    -1,    70,   154,    70,    -1,    41,    -1,
      41,   154,    41,    -1,   155,   344,    -1,   345,   344,    -1,
     345,   156,    -1
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
     827,   829,   831,   833,   836,   841,   845,   851,   853,   857,
     860,   863,   865,   868,   871,   873,   878,   881,   883,   888,
     891,   893,   898,   902,   908,   918,   920,   926,   930,   936,
     944,   954,   959,   962,   964,   970,   978,   983,   988,   991,
     993,  1001,  1011,  1018,  1020,  1022,  1024,  1026,  1028,  1029,
    1031,  1033,  1035,  1038,  1041,  1044,  1050,  1054,  1061,  1065,
    1067,  1069,  1071,  1073,  1075,  1077,  1079,  1081,  1083,  1085,
    1087,  1089,  1091,  1093,  1095,  1097,  1099,  1101,  1103,  1105,
    1107,  1109,  1111,  1113,  1115,  1117,  1119,  1121,  1123,  1125,
    1127,  1129,  1131,  1133,  1135,  1137,  1139,  1141,  1143,  1145,
    1151,  1158,  1162,  1164,  1166,  1168,  1170,  1172,  1174,  1176,
    1178,  1180,  1182,  1184,  1186,  1188,  1190,  1196,  1203,  1211,
    1217,  1219,  1223,  1227,  1232,  1235,  1237,  1243,  1247,  1252,
    1257,  1264,  1268,  1274,  1278,  1281,  1287,  1291,  1298,  1301,
    1307,  1311,  1318,  1320,  1322,  1324,  1328,  1332,  1337,  1340,
    1342,  1348,  1356,  1366,  1367,  1371,  1375,  1378,  1384,  1390,
    1397,  1401,  1409,  1418,  1424,  1430,  1437,  1441,  1449,  1458,
    1464,  1471,  1475,  1477,  1479,  1481,  1483,  1485,  1489,  1494,
    1501,  1503,  1506,  1508,  1510,  1512,  1514,  1516,  1517,  1518,
    1524,  1527,  1533,  1537,  1544,  1548,  1550,  1552,  1554,  1558,
    1562,  1566,  1568,  1572,  1576,  1580,  1584,  1586,  1588,  1590,
    1594,  1598,  1602,  1606,  1610,  1614,  1616,  1620,  1624,  1628,
    1632,  1636,  1640,  1644,  1648,  1652,  1656,  1660,  1662,  1664,
    1668,  1672,  1676,  1680,  1686,  1690,  1694,  1698,  1702,  1704,
    1708,  1712,  1716,  1718,  1720,  1722,  1724,  1726,  1730,  1732,
    1734,  1736,  1738,  1742,  1746,  1750,  1754,  1758,  1760,  1762,
    1766,  1770,  1774,  1776,  1778,  1780,  1784,  1788,  1790,  1794,
    1797,  1800
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    84,    84,    85,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   132,   133,   134,   135,   139,   140,   143,
     146,   150,   154,   158,   162,   164,   166,   168,   170,   172,
     177,   179,   181,   183,   185,   187,   192,   194,   196,   198,
     200,   202,   207,   209,   211,   213,   215,   217,   222,   226,
     233,   237,   244,   248,   256,   261,   263,   265,   267,   269,
     271,   273,   275,   277,   279,   281,   283,   285,   287,   289,
     291,   293,   295,   297,   299,   301,   303,   308,   310,   314,
     316,   321,   325,   330,   331,   335,   340,   345,   346,   350,
     354,   355,   359,   360,   361,   362,   366,   366,   367,   367,
     369,   369,   371,   371,   373,   373,   378,   379,   380,   381,
     385,   387,   392,   393,   394,   396,   398,   400,   402,   404,
     406,   408,   410,   412,   414,   416,   418,   420,   422,   424,
     426,   428,   430,   434,   438,   440,   445,   449,   453,   454,
     458,   460,   462,   464,   466,   471,   473,   475,   477,   479,
     481,   487,   489,   491,   493,   495,   497,   499,   501,   506,
     511,   513,   518,   520,   522,   524,   526,   528,   530,   532,
     534,   539,   543,   547,   548,   551,   555,   557,   561,   562,
     565,   569,   571,   575,   576,   579,   580,   584,   586,   588,
     590,   594,   595,   598,   599,   600,   601,   602,   603,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   613,   614,
     615,   616,   617,   621,   623,   625,   627,   629,   631,   636,
     638,   640,   645,   647,   649,   654,   659,   661,   666,   670,
     675,   680,   690,   695,   701,   711,   716,   727,   733,   741,
     751,   765,   769,   771,   775,   782,   791,   800,   804,   806,
     810,   819,   830,   842,   844,   846,   848,   850,   855,   856,
     857,   858,   859,   861,   868,   870,   872,   874,   879,   880,
     883,   884,   885,   886,   887,   888,   889,   890,   891,   892,
     893,   894,   895,   896,   897,   898,   899,   900,   901,   902,
     903,   904,   905,   906,   907,   908,   909,   910,   911,   912,
     913,   914,   915,   916,   917,   918,   919,   920,   921,   925,
     927,   932,   933,   937,   938,   939,   940,   941,   942,   943,
     944,   945,   946,   947,   948,   949,   953,   955,   960,   961,
     965,   966,   970,   975,   980,   981,   984,   988,   991,   995,
     997,   999,  1001,  1005,  1008,  1009,  1010,  1011,  1014,  1015,
    1016,  1017,  1020,  1021,  1024,  1025,  1028,  1031,  1035,  1036,
    1039,  1040,  1041,  1044,  1045,  1046,  1049,  1050,  1053,  1054,
    1055,  1056,  1057,  1058,  1060,  1061,  1062,  1063,  1064,  1065,
    1067,  1071,  1072,  1075,  1076,  1077,  1080,  1081,  1082,  1083,
    1086,  1087,  1090,  1091,  1092,  1093,  1094,  1097,  1097,  1097,
    1100,  1102,  1104,  1106,  1111,  1112,  1115,  1116,  1119,  1120,
    1121,  1122,  1123,  1124,  1125,  1126,  1127,  1128,  1129,  1130,
    1131,  1132,  1133,  1134,  1135,  1136,  1137,  1138,  1139,  1141,
    1142,  1143,  1145,  1146,  1147,  1148,  1149,  1150,  1151,  1152,
    1153,  1154,  1155,  1156,  1157,  1158,  1159,  1160,  1161,  1162,
    1163,  1164,  1165,  1166,  1167,  1168,  1169,  1170,  1171,  1172,
    1173,  1174,  1175,  1176,  1177,  1179,  1181,  1184,  1185,  1186,
    1187,  1188,  1189,  1190,  1191,  1192,  1194,  1202,  1203,  1207,
    1208,  1217
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
       2,     2,     2,     2,     2,   153,     2,     2,     2,   157,
     151,   152,     2,     2,     2,     2,   158,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   154,   150,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   155,   159,   156,     2,     2,     2,     2,     2,     2,
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
     145,   146,   147,   148,   149
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1481;
  const int parser::yynnts_ = 187;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 156;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 160;

  const unsigned int parser::yyuser_token_number_max_ = 404;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1219 "DynareBison.yy"


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

