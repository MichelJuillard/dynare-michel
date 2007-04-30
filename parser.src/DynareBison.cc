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
	  case 46:
#line 138 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 47:
#line 139 "DynareBison.yy"
    {driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 48:
#line 142 "DynareBison.yy"
    {driver.rplot();;}
    break;

  case 53:
#line 162 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 54:
#line 164 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 55:
#line 166 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 56:
#line 168 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 57:
#line 170 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 58:
#line 172 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 59:
#line 177 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 60:
#line 179 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 61:
#line 181 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 62:
#line 183 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 63:
#line 185 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 64:
#line 187 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 65:
#line 192 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 66:
#line 194 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 67:
#line 196 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 68:
#line 198 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 69:
#line 200 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 70:
#line 202 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 71:
#line 207 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 72:
#line 209 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 73:
#line 211 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 74:
#line 213 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 75:
#line 215 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 76:
#line 217 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 77:
#line 222 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 78:
#line 226 "DynareBison.yy"
    {
      driver.periods((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 79:
#line 233 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(3) - (2)].string_val));
		;}
    break;

  case 80:
#line 237 "DynareBison.yy"
    {
      driver.cutoff((yysemantic_stack_[(4) - (3)].string_val));
		;}
    break;

  case 81:
#line 244 "DynareBison.yy"
    {driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 82:
#line 249 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 83:
#line 251 "DynareBison.yy"
    {(yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 84:
#line 253 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 85:
#line 255 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 86:
#line 257 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 87:
#line 259 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 88:
#line 261 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 89:
#line 263 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 90:
#line 265 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 91:
#line 267 "DynareBison.yy"
    {(yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 92:
#line 269 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 93:
#line 271 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 94:
#line 273 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 95:
#line 275 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 96:
#line 277 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 97:
#line 279 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 98:
#line 281 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 99:
#line 283 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 100:
#line 285 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 101:
#line 287 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 102:
#line 289 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 103:
#line 291 "DynareBison.yy"
    {(yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val));;}
    break;

  case 104:
#line 296 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 105:
#line 298 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 106:
#line 302 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 107:
#line 304 "DynareBison.yy"
    {driver.end_initval();;}
    break;

  case 108:
#line 308 "DynareBison.yy"
    {driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 109:
#line 313 "DynareBison.yy"
    {driver.end_endval();;}
    break;

  case 112:
#line 323 "DynareBison.yy"
    {driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 113:
#line 328 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 116:
#line 338 "DynareBison.yy"
    {driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 119:
#line 346 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 120:
#line 347 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 122:
#line 352 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 123:
#line 352 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 124:
#line 353 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 125:
#line 354 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 126:
#line 355 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 127:
#line 356 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 128:
#line 357 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 129:
#line 358 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 130:
#line 359 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 131:
#line 360 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 136:
#line 372 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 137:
#line 374 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val));;}
    break;

  case 138:
#line 378 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 140:
#line 381 "DynareBison.yy"
    {(yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 141:
#line 383 "DynareBison.yy"
    {(yysemantic_stack_[(1) - (1)].string_val)->append(".0"); (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 142:
#line 385 "DynareBison.yy"
    {(yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 143:
#line 387 "DynareBison.yy"
    {(yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 144:
#line 389 "DynareBison.yy"
    {(yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 145:
#line 391 "DynareBison.yy"
    {(yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 146:
#line 393 "DynareBison.yy"
    {(yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 147:
#line 395 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val));;}
    break;

  case 148:
#line 397 "DynareBison.yy"
    {(yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val);;}
    break;

  case 149:
#line 399 "DynareBison.yy"
    {(yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 150:
#line 401 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 151:
#line 403 "DynareBison.yy"
    {(yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 152:
#line 405 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 153:
#line 407 "DynareBison.yy"
    {(yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 154:
#line 409 "DynareBison.yy"
    {(yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 155:
#line 411 "DynareBison.yy"
    {(yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 156:
#line 413 "DynareBison.yy"
    {(yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 157:
#line 415 "DynareBison.yy"
    {(yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 158:
#line 417 "DynareBison.yy"
    {(yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 159:
#line 421 "DynareBison.yy"
    {driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 160:
#line 425 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 161:
#line 427 "DynareBison.yy"
    {(yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 162:
#line 431 "DynareBison.yy"
    {driver.end_shocks();;}
    break;

  case 163:
#line 435 "DynareBison.yy"
    {driver.end_mshocks();;}
    break;

  case 166:
#line 445 "DynareBison.yy"
    {driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val));;}
    break;

  case 167:
#line 447 "DynareBison.yy"
    {driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 168:
#line 449 "DynareBison.yy"
    {driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 169:
#line 451 "DynareBison.yy"
    {driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 170:
#line 453 "DynareBison.yy"
    {driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 171:
#line 458 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 172:
#line 460 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(4) - (2)].string_val),(yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 173:
#line 462 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 174:
#line 464 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 175:
#line 466 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 176:
#line 468 "DynareBison.yy"
    {driver.add_period((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 177:
#line 474 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 178:
#line 476 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 179:
#line 478 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 180:
#line 480 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 181:
#line 482 "DynareBison.yy"
    {driver.add_value_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 182:
#line 484 "DynareBison.yy"
    {driver.add_value_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 183:
#line 486 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 184:
#line 488 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 185:
#line 493 "DynareBison.yy"
    {driver.do_sigma_e();;}
    break;

  case 186:
#line 498 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 187:
#line 500 "DynareBison.yy"
    {driver.end_of_row();;}
    break;

  case 188:
#line 505 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 189:
#line 507 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 190:
#line 509 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 191:
#line 511 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 192:
#line 513 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 193:
#line 515 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 194:
#line 517 "DynareBison.yy"
    {driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 195:
#line 519 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 196:
#line 521 "DynareBison.yy"
    {driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 197:
#line 526 "DynareBison.yy"
    {
 			driver.steady();
 		;}
    break;

  case 198:
#line 530 "DynareBison.yy"
    {driver.steady();;}
    break;

  case 202:
#line 542 "DynareBison.yy"
    {driver.check();;}
    break;

  case 203:
#line 544 "DynareBison.yy"
    {driver.check();;}
    break;

  case 207:
#line 556 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 208:
#line 558 "DynareBison.yy"
    {driver.simulate();;}
    break;

  case 212:
#line 570 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 213:
#line 572 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 214:
#line 574 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 215:
#line 576 "DynareBison.yy"
    {driver.stoch_simul();;}
    break;

  case 238:
#line 607 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val));;}
    break;

  case 239:
#line 609 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val));;}
    break;

  case 240:
#line 611 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 241:
#line 613 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 242:
#line 615 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 243:
#line 617 "DynareBison.yy"
    {driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 244:
#line 622 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 245:
#line 624 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 246:
#line 626 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 247:
#line 631 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 248:
#line 633 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 249:
#line 635 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 250:
#line 640 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 251:
#line 645 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 252:
#line 647 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 254:
#line 656 "DynareBison.yy"
    {driver.estim_params.type = 1;
		 driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
     delete (yysemantic_stack_[(2) - (2)].string_val);
		;}
    break;

  case 255:
#line 661 "DynareBison.yy"
    {driver.estim_params.type = 2;
		 driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
     delete (yysemantic_stack_[(1) - (1)].string_val);
		;}
    break;

  case 256:
#line 666 "DynareBison.yy"
    {driver.estim_params.type = 3;
		 driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
		 driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
     delete (yysemantic_stack_[(4) - (2)].string_val);
     delete (yysemantic_stack_[(4) - (4)].string_val);
		;}
    break;

  case 257:
#line 676 "DynareBison.yy"
    {
      driver.estim_params.prior=*(yysemantic_stack_[(3) - (1)].string_val);
      delete (yysemantic_stack_[(3) - (1)].string_val);
    ;}
    break;

  case 258:
#line 681 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.prior=*(yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
		;}
    break;

  case 259:
#line 687 "DynareBison.yy"
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

  case 260:
#line 697 "DynareBison.yy"
    {
      driver.estim_params.init_val=*(yysemantic_stack_[(1) - (1)].string_val);
      delete (yysemantic_stack_[(1) - (1)].string_val);
    ;}
    break;

  case 261:
#line 702 "DynareBison.yy"
    {driver.estim_params.init_val=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.low_bound=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.up_bound=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 262:
#line 713 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(3) - (1)].string_val);
 		 driver.estim_params.std=*(yysemantic_stack_[(3) - (3)].string_val);
     delete (yysemantic_stack_[(3) - (1)].string_val);
     delete (yysemantic_stack_[(3) - (3)].string_val);
 		;}
    break;

  case 263:
#line 719 "DynareBison.yy"
    {driver.estim_params.mean=*(yysemantic_stack_[(5) - (1)].string_val);
		 driver.estim_params.std=*(yysemantic_stack_[(5) - (3)].string_val);
		 driver.estim_params.p3=*(yysemantic_stack_[(5) - (5)].string_val);
     delete (yysemantic_stack_[(5) - (1)].string_val);
     delete (yysemantic_stack_[(5) - (3)].string_val);
     delete (yysemantic_stack_[(5) - (5)].string_val);
		;}
    break;

  case 264:
#line 727 "DynareBison.yy"
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

  case 265:
#line 737 "DynareBison.yy"
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

  case 266:
#line 751 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 267:
#line 755 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 268:
#line 757 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 269:
#line 761 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(5) - (4)].string_val);
         delete (yysemantic_stack_[(5) - (2)].string_val);
         delete (yysemantic_stack_[(5) - (4)].string_val);
				;}
    break;

  case 270:
#line 768 "DynareBison.yy"
    {driver.estim_params.type = 3;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.name2 = *(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 271:
#line 777 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
				 driver.estim_params.init_val=*(yysemantic_stack_[(4) - (3)].string_val);
         delete (yysemantic_stack_[(4) - (1)].string_val);
         delete (yysemantic_stack_[(4) - (3)].string_val);
				;}
    break;

  case 272:
#line 786 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 273:
#line 790 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 274:
#line 792 "DynareBison.yy"
    {driver.add_estimated_params_element();;}
    break;

  case 275:
#line 796 "DynareBison.yy"
    {driver.estim_params.type = 1;
				 driver.estim_params.name = *(yysemantic_stack_[(7) - (2)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(7) - (4)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(7) - (6)].string_val);
         delete (yysemantic_stack_[(7) - (2)].string_val);
         delete (yysemantic_stack_[(7) - (4)].string_val);
         delete (yysemantic_stack_[(7) - (6)].string_val);
				;}
    break;

  case 276:
#line 805 "DynareBison.yy"
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

  case 277:
#line 816 "DynareBison.yy"
    {driver.estim_params.type = 2;
				 driver.estim_params.name = *(yysemantic_stack_[(6) - (1)].string_val);
				 driver.estim_params.low_bound=*(yysemantic_stack_[(6) - (3)].string_val);
				 driver.estim_params.up_bound=*(yysemantic_stack_[(6) - (5)].string_val);
         delete (yysemantic_stack_[(6) - (1)].string_val);
         delete (yysemantic_stack_[(6) - (3)].string_val);
         delete (yysemantic_stack_[(6) - (5)].string_val);
				;}
    break;

  case 278:
#line 828 "DynareBison.yy"
    {(yyval.string_val) = new string("1");;}
    break;

  case 279:
#line 830 "DynareBison.yy"
    {(yyval.string_val) = new string("2");;}
    break;

  case 280:
#line 832 "DynareBison.yy"
    {(yyval.string_val) = new string("3");;}
    break;

  case 281:
#line 834 "DynareBison.yy"
    {(yyval.string_val) = new string("4");;}
    break;

  case 282:
#line 836 "DynareBison.yy"
    {(yyval.string_val) = new string("5");;}
    break;

  case 283:
#line 840 "DynareBison.yy"
    {(yyval.string_val) = new string("NaN");;}
    break;

  case 287:
#line 845 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 288:
#line 847 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 289:
#line 854 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 290:
#line 856 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 291:
#line 858 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 292:
#line 860 "DynareBison.yy"
    {driver.run_estimation();;}
    break;

  case 334:
#line 911 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 335:
#line 913 "DynareBison.yy"
    {driver.run_prior_analysis();;}
    break;

  case 351:
#line 939 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 352:
#line 941 "DynareBison.yy"
    {driver.run_posterior_analysis();;}
    break;

  case 353:
#line 945 "DynareBison.yy"
    {driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val));;}
    break;

  case 354:
#line 946 "DynareBison.yy"
    {driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val));;}
    break;

  case 357:
#line 956 "DynareBison.yy"
    {driver.set_varobs();;}
    break;

  case 358:
#line 961 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 361:
#line 970 "DynareBison.yy"
    {driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val));;}
    break;

  case 362:
#line 973 "DynareBison.yy"
    {driver.set_unit_root_vars();;}
    break;

  case 363:
#line 977 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 364:
#line 981 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val));;}
    break;

  case 365:
#line 983 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val));;}
    break;

  case 366:
#line 985 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val));;}
    break;

  case 367:
#line 987 "DynareBison.yy"
    {driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 368:
#line 990 "DynareBison.yy"
    {driver.set_osr_params();;}
    break;

  case 369:
#line 993 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 370:
#line 994 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 371:
#line 995 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 372:
#line 996 "DynareBison.yy"
    {driver.run_osr();;}
    break;

  case 373:
#line 999 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 374:
#line 1000 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 375:
#line 1001 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 376:
#line 1002 "DynareBison.yy"
    {driver.run_olr();;}
    break;

  case 381:
#line 1013 "DynareBison.yy"
    {driver.set_olr_inst();;}
    break;

  case 382:
#line 1017 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 385:
#line 1024 "DynareBison.yy"
    {driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val));;}
    break;

  case 386:
#line 1025 "DynareBison.yy"
    {driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val));;}
    break;

  case 387:
#line 1026 "DynareBison.yy"
    {driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val));;}
    break;

  case 388:
#line 1029 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 389:
#line 1030 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 390:
#line 1031 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val);;}
    break;

  case 391:
#line 1034 "DynareBison.yy"
    {driver.run_calib(0);;}
    break;

  case 392:
#line 1035 "DynareBison.yy"
    {driver.run_calib(1);;}
    break;

  case 393:
#line 1038 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 394:
#line 1039 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 395:
#line 1040 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 396:
#line 1041 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 397:
#line 1042 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 398:
#line 1043 "DynareBison.yy"
    {driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 399:
#line 1045 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val));;}
    break;

  case 400:
#line 1046 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val));;}
    break;

  case 401:
#line 1047 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val));;}
    break;

  case 402:
#line 1048 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val));;}
    break;

  case 403:
#line 1049 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val));;}
    break;

  case 404:
#line 1050 "DynareBison.yy"
    {driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val));;}
    break;

  case 405:
#line 1053 "DynareBison.yy"
    {driver.run_model_comparison();;}
    break;

  case 411:
#line 1065 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val));;}
    break;

  case 412:
#line 1066 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 413:
#line 1067 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val));;}
    break;

  case 414:
#line 1068 "DynareBison.yy"
    {driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val));;}
    break;

  case 415:
#line 1071 "DynareBison.yy"
    {(yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val);;}
    break;

  case 416:
#line 1072 "DynareBison.yy"
    {(yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);;}
    break;

  case 418:
#line 1076 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 419:
#line 1077 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 420:
#line 1078 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 421:
#line 1079 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 422:
#line 1082 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 423:
#line 1082 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 425:
#line 1086 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 426:
#line 1088 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 427:
#line 1090 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 428:
#line 1092 "DynareBison.yy"
    {driver.ramsey_policy();;}
    break;

  case 433:
#line 1104 "DynareBison.yy"
    {driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 434:
#line 1105 "DynareBison.yy"
    {driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 435:
#line 1106 "DynareBison.yy"
    {driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 436:
#line 1107 "DynareBison.yy"
    {driver.linear();;}
    break;

  case 437:
#line 1108 "DynareBison.yy"
    {driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 438:
#line 1109 "DynareBison.yy"
    {driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 439:
#line 1110 "DynareBison.yy"
    {driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 440:
#line 1111 "DynareBison.yy"
    {driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 441:
#line 1112 "DynareBison.yy"
    {driver.option_num("nocorr", "1");;}
    break;

  case 442:
#line 1113 "DynareBison.yy"
    {driver.option_num("nofunctions", "1");;}
    break;

  case 443:
#line 1114 "DynareBison.yy"
    {driver.option_num("nomoments", "1");;}
    break;

  case 444:
#line 1115 "DynareBison.yy"
    {driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 445:
#line 1116 "DynareBison.yy"
    {driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 446:
#line 1117 "DynareBison.yy"
    {driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 447:
#line 1118 "DynareBison.yy"
    {driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1");;}
    break;

  case 448:
#line 1119 "DynareBison.yy"
    {driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 449:
#line 1120 "DynareBison.yy"
    {driver.option_num("simul", "1");;}
    break;

  case 450:
#line 1121 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 451:
#line 1122 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 452:
#line 1123 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 453:
#line 1125 "DynareBison.yy"
    {driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 454:
#line 1126 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 455:
#line 1127 "DynareBison.yy"
    {driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 456:
#line 1129 "DynareBison.yy"
    {driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 457:
#line 1130 "DynareBison.yy"
    {driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 458:
#line 1131 "DynareBison.yy"
    {driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 459:
#line 1132 "DynareBison.yy"
    {driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 460:
#line 1133 "DynareBison.yy"
    {driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 461:
#line 1134 "DynareBison.yy"
    {driver.option_num("nograph","1");;}
    break;

  case 462:
#line 1135 "DynareBison.yy"
    {driver.option_num("nograph", "0");;}
    break;

  case 463:
#line 1136 "DynareBison.yy"
    {driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 464:
#line 1137 "DynareBison.yy"
    {driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 465:
#line 1138 "DynareBison.yy"
    {driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 466:
#line 1139 "DynareBison.yy"
    {driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 468:
#line 1141 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 469:
#line 1142 "DynareBison.yy"
    {driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 470:
#line 1143 "DynareBison.yy"
    {driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 471:
#line 1144 "DynareBison.yy"
    {driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 472:
#line 1145 "DynareBison.yy"
    {driver.option_num("mode_check", "1");;}
    break;

  case 473:
#line 1146 "DynareBison.yy"
    {driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 474:
#line 1147 "DynareBison.yy"
    {driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 475:
#line 1148 "DynareBison.yy"
    {driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 476:
#line 1149 "DynareBison.yy"
    {driver.option_num("load_mh_file", "1");;}
    break;

  case 477:
#line 1150 "DynareBison.yy"
    {driver.option_num("loglinear", "1");;}
    break;

  case 478:
#line 1151 "DynareBison.yy"
    {driver.option_num("nodiagnostic", "1");;}
    break;

  case 479:
#line 1152 "DynareBison.yy"
    {driver.option_num("bayesian_irf", "1");;}
    break;

  case 480:
#line 1153 "DynareBison.yy"
    {driver.option_num("TeX", "1");;}
    break;

  case 481:
#line 1154 "DynareBison.yy"
    {driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 482:
#line 1155 "DynareBison.yy"
    {driver.option_num("smoother", "1");;}
    break;

  case 483:
#line 1156 "DynareBison.yy"
    {driver.option_num("moments_varendo", "1");;}
    break;

  case 484:
#line 1157 "DynareBison.yy"
    {driver.option_num("filtered_vars", "1");;}
    break;

  case 485:
#line 1158 "DynareBison.yy"
    {driver.option_num("relative_irf", "1");;}
    break;

  case 486:
#line 1159 "DynareBison.yy"
    {driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 487:
#line 1160 "DynareBison.yy"
    {driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 488:
#line 1161 "DynareBison.yy"
    {driver.option_num("olr_beta", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 489:
#line 1164 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 490:
#line 1166 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 491:
#line 1168 "DynareBison.yy"
    {driver.option_num("noprint", "0");;}
    break;

  case 492:
#line 1169 "DynareBison.yy"
    {driver.option_num("noprint", "1");;}
    break;

  case 493:
#line 1170 "DynareBison.yy"
    {driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 494:
#line 1171 "DynareBison.yy"
    {driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 495:
#line 1172 "DynareBison.yy"
    {driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 496:
#line 1173 "DynareBison.yy"
    {driver.option_num("noconstant", "0");;}
    break;

  case 497:
#line 1174 "DynareBison.yy"
    {driver.option_num("noconstant", "1");;}
    break;

  case 498:
#line 1175 "DynareBison.yy"
    {driver.option_num("load_mh_file", "-1");;}
    break;

  case 499:
#line 1176 "DynareBison.yy"
    {driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val));;}
    break;

  case 500:
#line 1179 "DynareBison.yy"
    {
    (yysemantic_stack_[(3) - (1)].string_val)->append(":");
    (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
    delete (yysemantic_stack_[(3) - (3)].string_val);
    (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
  ;}
    break;

  case 502:
#line 1188 "DynareBison.yy"
    { (yysemantic_stack_[(3) - (1)].string_val)->append(":"); (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val)); delete (yysemantic_stack_[(3) - (3)].string_val); (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); ;}
    break;

  case 503:
#line 1191 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 504:
#line 1193 "DynareBison.yy"
    {
               (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
               (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
               delete (yysemantic_stack_[(2) - (2)].string_val);
               (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
             ;}
    break;

  case 505:
#line 1201 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2163 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -887;
  const short int
  parser::yypact_[] =
  {
       990,   436,   -96,   452,   365,    46,   -17,   -15,   -29,   164,
       5,    81,   107,   142,   458,   482,   206,   210,   352,   242,
     186,   327,   258,   188,   327,   343,   237,  -887,   265,   270,
     327,   274,   430,   496,   514,   193,   195,   327,   374,   392,
     424,   327,   647,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,   501,   262,  -887,
     413,   490,   375,    36,   126,   472,   140,   474,   494,   541,
    -887,   761,    -9,   144,   203,   216,   505,   494,   551,  -887,
     348,    19,    96,   558,   512,  -887,  1114,   162,   232,   518,
    -887,  1114,   236,   260,   479,   264,   555,   454,   855,   797,
     797,   298,    96,   456,  -887,   535,  -887,   413,  -887,  1159,
     308,  -887,  1044,   311,   313,   522,   318,   524,   337,   526,
     359,   369,  -887,  -887,   492,   583,   -28,    60,  -887,   635,
     278,  -887,  -887,   519,  -887,   531,  -887,  -887,   613,   -32,
    -887,   619,   136,   666,    45,  -887,   621,  -887,   678,  -887,
     686,   690,  -887,   691,   698,  -887,   700,   702,   706,   709,
     715,  -887,  -887,   716,   721,   722,   723,   732,   736,  -887,
    -887,   738,   739,  -887,   740,  -887,  -887,  -887,   750,   758,
     765,   770,  -887,  -887,   772,   773,   385,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,   776,   696,  -887,
     734,  -887,   755,    67,  -887,   699,   759,   712,   766,   255,
    -887,   774,   713,   779,   259,  -887,   692,    54,  -887,    55,
     818,   694,   742,  -887,   386,   703,   733,   824,  -887,  -887,
     401,  -887,  -887,  -887,  -887,   780,   781,    71,  -887,  -887,
    -887,   701,   558,   558,   710,   714,   743,   745,   747,   748,
     751,   752,   754,   767,   558,   521,   768,    57,  -887,   840,
     843,   865,   866,   877,   892,  -887,  -887,  -887,   896,   897,
     900,   913,  -887,   914,  -887,   924,   925,  -887,  -887,   405,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,   805,    61,   407,  -887,  -887,  -887,   819,
     881,  -887,   802,  -887,  -887,  -887,   804,   855,   855,   809,
     810,   811,   812,   826,   829,   831,   834,   835,   836,   855,
     729,  -887,   408,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,   409,  -887,    90,
      24,   414,  -887,  -887,   420,  -887,  -887,   422,  -887,  -887,
     929,  -887,   427,  -887,  -887,  -887,  -887,  -887,   872,   931,
    -887,  -887,   886,   933,  -887,  -887,   898,   937,  -887,  -887,
     858,   861,   953,    52,  1001,  -887,  -887,   984,   413,   879,
    -887,  -887,   880,   -22,   957,   882,    32,   967,   558,  -887,
    -887,  -887,  1012,   976,   894,  1005,  1009,  1010,  1011,  1016,
    1022,  1023,  1017,   360,  1033,  1025,  1026,  1027,  1031,  1004,
       7,   934,  1034,  1040,  1056,  1020,  1028,   761,    73,  1029,
    1072,   970,  -887,  -887,  -887,    63,   973,   233,   978,  -887,
    -887,   979,   233,   980,  -887,  -887,    20,  -887,  -887,  -887,
    1039,   960,  1045,    30,  -887,    23,  -887,   102,  -887,   966,
     983,   217,    19,   226,   993,    34,  -887,  -887,   558,   995,
     468,   558,   558,   558,   558,   558,   558,   558,   558,   558,
     558,   493,   558,   558,   558,   558,   558,  -887,   558,  -887,
    -887,  1061,  1093,  1098,  1103,  1105,  1113,   233,  1115,  1119,
     361,  1120,  1128,  1130,  1114,   190,  1065,   575,  -887,   830,
     208,  -887,  1035,  -887,    20,  1042,   530,   855,   855,   855,
     855,   855,   855,   855,   855,   855,   855,   506,   855,   855,
     855,   855,   855,  1006,   797,   221,   223,  -887,  -887,  -887,
     558,   341,    62,   535,  1014,   413,  1032,  1159,   225,  1151,
    1044,   227,  -887,  1057,  -887,  1069,  -887,  1070,  -887,  1146,
    1041,  1043,  1049,   558,  -887,  -887,  -887,  -887,  -887,   377,
    1053,  -887,  -887,   379,  1063,  1187,  -887,  -887,  1152,     4,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  1036,
    -887,  -887,  -887,  -887,  1062,  -887,  -887,  -887,   380,  -887,
    1126,  1157,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
     362,  1067,  1094,  1102,  1161,  1106,   233,  1165,  1091,   233,
    -887,  1201,  1202,  1096,  -887,   494,  1222,  -887,  -887,  -887,
     855,  -887,  -887,  -887,   429,  -887,  -887,  1104,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,   -36,   -19,
    -887,  1186,   558,  1190,   197,   845,   431,   549,   590,   624,
     650,   657,   778,   784,   798,   891,   905,  -887,   468,   468,
     995,   995,  1129,   911,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
     381,   558,  -887,  1196,  1193,  -887,   389,  -887,  1127,   994,
    1008,  1038,  1050,  1078,  1095,  1109,  1122,  1140,  1147,  -887,
     530,   530,  1042,  1042,  1129,  -887,  -887,  -887,   390,  -887,
     395,  1153,    24,  1135,  -887,  -887,    27,   558,  -887,  -887,
    -887,  -887,  -887,  -887,   423,  -887,  -887,  -887,   428,  -887,
    -887,  -887,  1142,  1263,  -887,  -887,  1200,  -887,   229,  -887,
     376,  -887,  1118,  -887,  -887,  -887,  1237,  -887,   444,  1244,
    -887,  -887,  -887,  -887,  -887,  -887,   233,    63,  1207,   233,
    1208,  1209,  -887,  1177,  -887,  -887,  1314,   199,   855,  1211,
     102,  -887,   742,   742,   742,   226,  -887,   233,  -887,  1321,
    1218,  1332,  1316,   558,   558,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  1210,  -887,  1224,   558,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,    62,  -887,  -887,  -887,
     558,  1160,  -887,  -887,  1041,   558,  -887,  -887,   435,  -887,
     440,  1317,  1212,  1036,  -887,  -887,  -887,  1239,  1253,  1255,
     233,  1233,   233,   233,  -887,   558,  -887,  1234,  -887,  -887,
    1241,    59,   182,   213,   150,  1252,   558,  -887,   558,  1235,
      80,  1242,   845,  -887,  -887,  1248,  1167,  -887,  1369,  1257,
    -887,  -887,  -887,  1279,  -887,   233,   233,   233,  1285,  -887,
    1264,  1265,  1266,  -887,   742,  -887,  -887,  -887,   233,  -887,
    1272,  1280,  1375,  1270,  1383,  1306,  -887,  -887,  -887,   558,
    -887,    17,  1300,  -887,  1301,   233,  -887,  -887,  -887,   478,
    1277,  -887,  -887,  -887,  1389,  1278,    74,  1290,  1363,  -887,
     233,   183,  1284,  -887,  -887,  -887,  1393,  -887,  -887,   367,
     510,   558,    77,  -887,  -887,  -887,  1281,  1309,  1310,  -887,
    -887,  -887,  -887,  1173,  -887,  -887,   558,  -887,  -887,  -887,
     233,   233,  -887,  1180,  1312,  -887,  -887,   233,  -887
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   422,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     2,     4,    28,    29,    43,    44,    45,
      42,     5,     6,    11,     8,     9,    10,     7,    12,    13,
      14,    15,    16,    17,    18,    22,    24,    23,    19,    20,
      21,    25,    26,    27,    30,    31,    32,    37,    38,    33,
      34,    35,    36,    39,    40,    41,   391,     0,     0,   202,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   242,
     289,     0,     0,     0,     0,     0,     0,     0,     0,   122,
       0,     0,     0,     0,     0,   373,     0,     0,     0,     0,
     369,     0,     0,     0,    73,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   207,     0,   197,     0,   212,     0,
       0,   425,     0,     0,     0,    55,     0,    61,     0,    67,
       0,     0,     1,     3,     0,     0,   388,     0,   384,     0,
       0,   205,   206,     0,    79,     0,    46,   401,     0,     0,
     395,     0,     0,     0,     0,   111,     0,   479,     0,   496,
       0,     0,   484,     0,     0,   462,     0,     0,     0,     0,
       0,   476,   477,     0,     0,     0,     0,     0,     0,   498,
     472,     0,     0,   483,     0,   497,   478,   461,     0,     0,
       0,     0,   482,   480,     0,     0,     0,   294,   330,   319,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   320,   321,   322,   323,   324,   325,
     326,   327,   328,   329,   331,   332,   333,   238,     0,   291,
       0,   255,     0,     0,   252,     0,     0,     0,     0,     0,
     274,     0,     0,     0,     0,   268,     0,     0,   115,     0,
       0,     0,     0,   436,     0,     0,     0,     0,   492,   491,
       0,   407,   408,   409,   410,     0,     0,     0,   165,    84,
      85,    83,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   360,     0,
       0,     0,     0,     0,     0,   441,   442,   443,     0,     0,
       0,     0,   485,     0,   449,     0,     0,   378,   379,     0,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   231,   232,   233,   234,   235,   236,   237,   230,
     377,   375,   381,     0,     0,     0,   371,   368,    76,    71,
       0,    52,     0,    77,   140,   141,   160,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     423,   139,     0,   337,   342,   338,   339,   340,   341,   343,
     344,   345,   346,   347,   348,   349,   350,     0,    48,     0,
       0,     0,   210,   211,     0,   200,   201,     0,   217,   214,
       0,   431,     0,   430,   432,   427,   362,    58,    53,     0,
      49,    64,    59,     0,    50,    70,    65,     0,    51,   357,
       0,     0,     0,     0,     0,   382,   383,     0,     0,     0,
      80,    47,     0,     0,     0,     0,     0,     0,     0,   109,
     110,   243,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     240,     0,   254,   250,   251,   283,     0,   283,     0,   272,
     273,     0,   283,     0,   266,   267,     0,   113,   114,   106,
       0,     0,     0,     0,   134,     0,   135,     0,   130,     0,
       0,     0,     0,     0,     0,     0,   163,   164,     0,    91,
      92,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    81,     0,   358,
     359,     0,     0,     0,     0,     0,     0,   283,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   363,     0,
       0,    74,    72,    78,     0,   147,   148,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   162,   195,   196,
       0,     0,   187,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    56,    54,    62,    60,    68,    66,   392,     0,
     388,     0,     0,     0,   434,   204,   203,   404,   399,     0,
       0,   398,   393,     0,     0,     0,   463,   453,     0,     0,
     495,   456,   481,   444,   486,   487,   459,   460,   465,   468,
     469,   466,   474,   475,   464,   471,   470,   455,   454,     0,
     457,   458,   473,   493,     0,   494,   293,   290,     0,   239,
       0,     0,   278,   285,   279,   284,   281,   286,   280,   282,
       0,     0,     0,   260,     0,     0,   283,     0,     0,   283,
     246,     0,     0,     0,   108,     0,     0,   123,   132,   133,
       0,   137,   120,   119,     0,   118,   121,     0,   126,   124,
     489,   490,   406,   417,   419,   420,   421,   418,     0,   411,
     415,     0,     0,     0,     0,   104,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    82,    87,    86,
      88,    89,    90,     0,   440,   448,   433,   439,   445,   446,
     488,   437,   447,   452,   451,   438,   435,   450,   380,   374,
       0,     0,   366,     0,     0,   370,     0,    75,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   138,
     143,   142,   144,   145,   146,   424,   336,   334,     0,   351,
       0,     0,     0,     0,   192,   193,     0,     0,   209,   208,
     199,   198,   216,   213,     0,   499,   429,   426,     0,    57,
      63,    69,     0,     0,   390,   389,     0,   400,     0,   394,
       0,   112,   501,   503,   505,   504,     0,   355,     0,     0,
     292,   241,   256,   288,   287,   253,   283,   283,     0,   283,
       0,     0,   271,     0,   245,   244,     0,     0,     0,     0,
       0,   128,     0,     0,     0,     0,   405,   283,   416,     0,
       0,     0,     0,     0,     0,   103,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,     0,   376,     0,     0,
     364,   372,   161,   149,   150,   151,   152,   153,   154,   155,
     156,   157,   158,   335,   352,   194,   186,   185,   189,   190,
       0,     0,   215,   428,   388,     0,   385,   402,     0,   396,
       0,     0,     0,     0,   467,   500,   257,     0,     0,     0,
     283,     0,   283,   283,   269,     0,   107,     0,   136,   117,
       0,     0,     0,     0,   412,     0,     0,   168,     0,   176,
       0,     0,   105,   361,   367,     0,     0,   191,     0,     0,
     403,   397,   502,     0,   356,   283,   283,   283,     0,   277,
       0,     0,     0,   159,     0,   131,   127,   125,   283,   413,
       0,     0,     0,   171,     0,     0,   167,   365,   188,     0,
     386,   283,   262,   258,   261,   283,   275,   270,   116,     0,
       0,   170,   169,   175,     0,   173,     0,     0,     0,   354,
     283,     0,     0,   129,   414,   172,     0,   249,   182,     0,
       0,     0,     0,   181,   180,   387,     0,   263,     0,   276,
     174,   248,   247,     0,   179,   166,     0,   178,   177,   353,
     283,   283,   184,     0,   264,   259,   183,   283,   265
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
      -887,  -887,  1398,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,  -887,  -887,  -296,  -887,  -887,
    -887,  -887,  -102,  -172,  -887,  -887,  1164,  -887,   592,  -887,
    -887,  -887,  -887,  -887,  -887,  -780,  -501,  -108,  -491,  -887,
    -887,  -887,  1311,  -234,  -887,  -887,  -887,  -887,   652,  -887,
    -887,   841,  -887,  -887,   997,  -887,  -887,   844,  -887,  -887,
    -116,   -20,  -561,   437,  -887,  -887,  1185,  -887,  -887,  -886,
    -887,  -887,  1176,  -887,  -887,  1182,  -793,  -468,  -887,  -887,
     965,  -887,  1323,   860,  -887,   542,  -887,  -887,  -887,  -887,
    1139,  -887,  -887,  -887,  -887,  -887,  -887,   893,  1337,  -887,
    -887,  -887,  1302,  -578,  -887,  -887,  -887,  -887,  -887,   938,
    -887,   606,  -676,  -887,  -887,  -887,  -887,  -887,   852,  -887,
     -82,  -887,  1353,  -887,  -887,  -887,  -887,  -887,  -887,  -887,
     -92,  -887,  -887,  -124,  -477,  -887,  -887,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,   -87,   -77,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,  -887,   -76,  -887,  -887,  -887,  -887,
    -887,   -73,   -71,   -70,   -69,   -68,   -66,  -887,  -887,  -887,
    -887,  -887,  -887,  -887,   -65,   -62,   -60,  -887,  -887,  -887,
    -887,  -887,   825,  -887,   985
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    42,    43,    44,    45,    46,    47,    48,    49,    50,
     146,   148,   150,   125,    51,    52,    53,   315,   726,    54,
     281,    55,   174,   175,    56,   277,   278,   704,   705,    57,
     282,   854,   853,   930,   707,   513,   514,   515,   516,   391,
      58,    59,   297,   298,   940,  1012,    60,   601,   602,    61,
     414,   415,    62,   160,   161,    63,   411,   412,    64,   417,
     337,   102,   693,  1014,    65,   263,   264,   265,   681,   916,
      66,   274,   275,    67,   269,   270,   682,   917,    68,   216,
     217,    69,   392,   393,    70,   827,   828,    71,    72,   317,
     318,    73,    74,   364,    75,    76,    77,   338,   339,    78,
      79,   157,   158,   444,    80,    81,    82,    83,   290,   291,
     718,   719,   720,    84,   128,   593,    85,   422,   423,   340,
     341,   342,   343,   344,   345,   346,   347,   348,   349,   350,
     351,   352,   353,   354,   355,   356,   357,   358,   220,   221,
     222,   223,   224,   225,   226,   395,   396,   229,   230,   231,
     232,   233,   234,   235,   236,   397,   238,   239,   240,   241,
     242,   398,   399,   400,   401,   402,   403,   359,   249,   250,
     360,   292,   293,   294,   404,   405,   406,   254,   255,   256,
     424,   665,   823,   639,   640
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       117,   118,   460,   122,   123,   279,   529,   530,   162,   219,
     131,   413,   698,   768,   227,   140,   143,   144,   541,   218,
     390,   151,   699,   418,   228,   237,   421,   683,   243,   685,
     244,   245,   246,   247,   688,   248,   251,   394,   394,   252,
     706,   253,   813,   858,   918,   822,   700,    99,   657,   673,
     713,   697,    94,    88,    96,   416,   598,   722,   675,   898,
     257,   690,   374,   527,   155,   599,   459,   567,   899,   672,
     983,   375,   931,   932,   933,   507,   509,   165,   549,   260,
     965,   445,   568,   295,   621,   287,   677,    93,   493,   750,
     855,   374,   526,   622,   794,   673,   288,   674,   442,   376,
     375,    99,   295,   795,   675,   676,  1007,   460,   295,  1007,
     714,   597,   289,   856,   173,   690,   320,   258,   690,   453,
      98,   973,   443,   276,   173,   454,   316,   628,   376,   156,
     569,   857,   677,    95,   715,    97,   261,   702,   716,   717,
     259,   678,    99,  1008,   680,  1035,  1024,   691,   692,   703,
     588,   589,   590,   591,   103,   592,   260,   377,   378,   824,
     723,   638,   635,   379,   380,   381,   382,   383,   384,   385,
     386,   387,   701,   998,   600,   527,   262,   900,   388,   679,
     389,   632,   512,   724,   989,   166,   377,   378,   796,   672,
     680,   296,   379,   380,   381,   382,   383,   384,   385,   386,
     387,  1009,  1010,   966,  1009,  1010,   974,   388,  1018,   389,
     296,   512,   797,   261,   374,   266,   296,   674,   840,   713,
     926,   843,   667,   375,  1011,   676,  1025,  1026,   271,   975,
     104,   257,   725,    99,   967,   727,   728,   729,   730,   731,
     732,   733,   734,   735,   736,   374,   738,   739,   740,   741,
     742,   376,   743,   262,   375,    99,   105,    99,   858,    99,
     126,   678,    99,   710,    99,   673,   155,   266,   173,   575,
     576,   271,   267,   764,   675,   167,   499,    99,   127,   714,
     504,   587,   376,   168,   711,   272,   862,   456,   258,   170,
      99,   106,    99,   457,    99,   713,    99,   171,    99,   679,
     968,   257,   677,   715,   791,   257,   863,   716,   717,   377,
     378,   361,   268,   100,   101,   379,   380,   381,   382,   383,
     384,   385,   386,   387,   267,   273,   948,   816,   272,   257,
     388,   156,   389,   369,   512,   115,   116,   120,   121,   759,
     377,   378,   138,   139,   141,   142,   379,   380,   381,   382,
     383,   384,   385,   386,   387,   714,   111,   765,   258,   112,
     680,   388,   258,   389,   268,   512,   162,   257,   273,   919,
     787,   921,   789,   706,   803,   113,   807,   257,   907,   715,
     257,   362,   257,   716,   717,   366,   258,   428,    91,   935,
     370,   114,   649,   753,   833,   219,    99,    92,   283,  1021,
     227,   650,   754,   834,   448,   218,   432,   119,   844,   367,
     228,   237,   124,   371,   243,   129,   244,   245,   246,   247,
     130,   248,   251,   132,   258,   252,   860,   253,   436,   449,
     698,   698,   698,   629,   258,  1013,   633,   258,   257,   258,
     699,   699,   699,   145,   429,    99,   257,   408,   257,   257,
     257,  1027,   958,   133,   960,   961,   284,   419,   257,   257,
     425,   147,   426,   433,   257,   878,   285,   430,   668,   769,
     770,   771,   772,   773,   774,   775,   776,   777,   778,   413,
     780,   781,   782,   783,   784,   437,   434,   982,   698,   984,
     792,   802,   257,   149,   421,   258,   793,   257,   699,  1003,
     990,   901,   394,   258,   257,   258,   258,   258,   438,   257,
     374,   487,   517,   999,   154,   258,   258,  1002,   439,   375,
     159,   258,   163,   416,   164,   909,   817,   522,   819,   830,
     877,   564,  1017,   564,   594,   594,   488,   518,   881,   893,
     603,   169,  1022,   172,   894,   760,   605,   376,   607,   258,
     766,   845,   523,   610,   258,   850,   565,   864,   570,   595,
     596,   258,  1034,   173,   176,   604,   258,   941,   942,  1038,
     913,   606,   902,   608,   276,   788,   790,   903,   611,   280,
     851,   316,   865,   945,   950,    86,    87,   363,   804,   951,
     299,   808,   849,   847,   368,   914,   372,   544,   545,   300,
     546,    89,    90,   373,   946,   377,   378,   107,   108,   949,
     410,   379,   380,   381,   382,   383,   384,   385,   386,   387,
     542,   543,   544,   545,   330,   546,   388,   301,   389,   962,
     512,   109,   110,   588,   589,   590,   591,   427,   592,   431,
     970,   435,   971,   440,   737,   134,   135,   152,   542,   543,
     544,   545,   441,   546,     1,     2,     3,   779,   447,   590,
     591,     4,   592,   136,   137,     5,     6,     7,   450,     8,
     547,     9,    10,    11,    12,   460,   542,   543,   544,   545,
     451,   546,   452,   997,    13,   302,   303,    14,   455,   458,
     461,   304,   305,   306,   307,   308,   309,   310,   311,   312,
     866,   462,   542,   543,   544,   545,   313,   546,   314,   463,
      15,    16,    17,   464,   465,  1023,    18,   542,   543,   544,
     545,   466,   546,   467,   762,   468,    19,    20,    21,   469,
    1033,    22,   470,    23,    24,    25,    26,    27,   471,   472,
     927,   867,    28,    29,   473,   474,   475,    30,    31,    32,
      33,   542,   543,   544,   545,   476,   546,    34,    35,   477,
      36,   478,   479,   480,    37,   490,   177,    38,    39,    40,
      41,   178,   179,   481,   374,   868,   180,   542,   543,   544,
     545,   482,   546,   375,   542,   543,   544,   545,   483,   546,
     181,   182,   183,   484,   184,   485,   486,   185,   908,   489,
     910,   869,   177,   491,   186,   187,   188,   178,   870,   189,
     190,   376,   191,   192,   193,   194,   195,   196,   197,   198,
     199,   200,   201,   202,   492,   495,   181,   182,   496,   203,
     184,   204,   205,   185,   206,   498,   207,   299,   497,   502,
     186,   510,   506,   501,   208,   511,   300,   521,   503,   524,
     525,   528,   209,   210,   519,   211,   588,   589,   590,   591,
     531,   592,   299,   551,   532,   203,   552,   212,   159,   377,
     378,   300,   207,   213,   301,   379,   380,   381,   382,   383,
     384,   385,   386,   387,   520,   214,   215,   374,   553,   554,
     388,   211,   389,   533,   512,   534,   375,   535,   536,   301,
     555,   537,   538,   212,   539,   542,   543,   544,   545,   213,
     546,   542,   543,   544,   545,   556,   546,   540,   548,   557,
     558,   214,   215,   559,   376,   542,   543,   544,   545,   871,
     546,   566,   302,   303,   571,   872,   560,   561,   304,   305,
     306,   307,   308,   309,   310,   311,   312,   562,   563,   873,
     572,   573,   609,   313,   574,   314,   763,   302,   303,   577,
     578,   579,   580,   304,   305,   306,   307,   308,   309,   310,
     311,   312,   542,   543,   544,   545,   581,   546,   313,   582,
     314,   583,   377,   378,   584,   585,   586,   612,   379,   380,
     381,   382,   383,   384,   385,   386,   387,     1,     2,     3,
     613,   614,   615,   388,     4,   389,   617,   618,     5,     6,
       7,   619,     8,   616,     9,    10,    11,    12,   542,   543,
     544,   545,   620,   546,   623,   624,   630,    13,   626,   627,
      14,   631,   542,   543,   544,   545,   634,   546,   542,   543,
     544,   545,   874,   546,   636,   637,   641,   319,   638,   648,
     642,   643,   644,    15,    16,    17,   875,   645,   320,    18,
     321,   322,   876,   646,   647,   651,   652,   653,   654,    19,
      20,    21,   655,   656,    22,   660,    23,    24,    25,    26,
      27,   661,   323,   324,   659,    28,    29,   186,   662,   663,
      30,    31,    32,    33,   283,   670,   671,   664,   669,   684,
      34,    35,   744,    36,   686,   687,   689,    37,   694,   695,
      38,    39,    40,    41,   696,   708,   325,   319,   326,   721,
     327,   588,   589,   590,   591,   745,   592,   546,   320,   329,
     321,   322,   709,   330,   761,   588,   589,   590,   591,   746,
     592,   331,   332,   333,   747,   883,   748,   334,   335,   336,
     767,   159,   323,   324,   749,   785,   751,   186,   420,   884,
     752,   755,   319,   799,   283,   588,   589,   590,   591,   756,
     592,   757,   809,   320,   592,   321,   322,   588,   589,   590,
     591,   801,   592,   805,   810,   811,   325,   812,   326,   885,
     327,   443,   826,   822,   814,   831,   328,   323,   324,   329,
     815,   886,   186,   330,   818,   588,   589,   590,   591,   283,
     592,   331,   332,   333,   820,   829,   835,   334,   335,   336,
     836,   159,   588,   589,   590,   591,   832,   592,   837,   887,
     838,   325,   839,   326,   841,   327,   588,   589,   590,   591,
     842,   592,   844,   845,   329,   848,   888,   846,   330,   588,
     589,   590,   591,   852,   592,   859,   331,   332,   333,   861,
     889,    -1,   334,   335,   336,   879,   159,   588,   589,   590,
     591,   911,   592,   890,   588,   589,   590,   591,   882,   592,
     542,   543,   544,   545,   897,   546,   905,   542,   543,   544,
     545,   891,   546,   904,   542,   543,   544,   545,   892,   546,
     542,   543,   544,   545,   895,   546,   912,   542,   543,   544,
     545,   947,   546,   915,   542,   543,   544,   545,   978,   546,
     542,   543,   544,   545,  1032,   546,   924,   542,   543,   544,
     545,  1036,   546,   920,   922,   923,   821,   925,   588,   589,
     590,   591,   880,   592,   936,   542,   543,   544,   545,   906,
     546,   542,   543,   544,   545,   938,   546,   939,   952,   943,
     928,   588,   589,   590,   591,   955,   592,   937,   953,   542,
     543,   544,   545,   944,   546,   542,   543,   544,   545,   956,
     546,   957,   959,   963,   542,   543,   544,   545,   972,   546,
     964,   976,   979,   542,   543,   544,   545,   977,   546,   542,
     543,   544,   545,   969,   546,   981,   980,   542,   543,   544,
     545,   985,   546,   986,   987,   988,   993,   542,   543,   544,
     545,   991,   546,   994,   995,   996,  1000,  1001,  1004,   992,
    1005,  1006,  1016,  1019,  1020,  1030,  1031,  1029,  1037,  1015,
     153,   508,   929,   409,   896,   625,   800,   798,   494,  1028,
     505,   500,   666,   407,   786,   954,   550,   758,   365,   446,
     712,   934,   806,   286,   825,   658
  };

  /* YYCHECK.  */
  const unsigned short int
  parser::yycheck_[] =
  {
        20,    21,   174,    23,    24,   107,   302,   303,    90,   101,
      30,   135,   513,   574,   101,    35,    36,    37,   314,   101,
     128,    41,   513,   139,   101,   101,   142,   495,   101,   497,
     101,   101,   101,   101,   502,   101,   101,   129,   130,   101,
     517,   101,   620,   719,   837,    41,    23,    69,    41,    32,
      69,    21,    69,   149,    69,   137,    32,    23,    41,    32,
      69,    41,    32,   297,     4,    41,    21,   363,    41,     6,
     956,    41,   852,   853,   854,    21,    21,    41,    21,    12,
      21,    21,    21,    12,    32,    66,    69,    41,    21,   557,
     126,    32,    21,    41,    32,    32,    77,    34,   126,    69,
      41,    69,    12,    41,    41,    42,    32,   279,    12,    32,
     129,    21,    93,   149,    69,    41,    14,   126,    41,   151,
     149,    41,   150,    69,    69,   157,    69,   149,    69,    69,
      69,   150,    69,   150,   153,   150,    69,    35,   157,   158,
     149,    78,    69,    69,   127,  1031,    69,   127,   128,    47,
     127,   128,   129,   130,   149,   132,    12,   127,   128,   155,
     126,   154,   458,   133,   134,   135,   136,   137,   138,   139,
     140,   141,   149,   156,   150,   409,   109,   150,   148,   116,
     150,   149,   152,   149,   964,   149,   127,   128,   126,     6,
     127,   120,   133,   134,   135,   136,   137,   138,   139,   140,
     141,   127,   128,    21,   127,   128,   126,   148,  1001,   150,
     120,   152,   150,    69,    32,    12,   120,    34,   686,    69,
      21,   689,   149,    41,   150,    42,   149,   150,    12,   149,
     149,    69,   528,    69,    21,   531,   532,   533,   534,   535,
     536,   537,   538,   539,   540,    32,   542,   543,   544,   545,
     546,    69,   548,   109,    41,    69,   149,    69,   934,    69,
      23,    78,    69,    46,    69,    32,     4,    12,    69,   377,
     378,    12,    69,   569,    41,   149,    21,    69,    41,   129,
      21,   389,    69,   157,    67,    69,    89,   151,   126,   149,
      69,   149,    69,   157,    69,    69,    69,   157,    69,   116,
     150,    69,    69,   153,   600,    69,   109,   157,   158,   127,
     128,   149,   109,   149,   150,   133,   134,   135,   136,   137,
     138,   139,   140,   141,    69,   109,   904,   623,    69,    69,
     148,    69,   150,    69,   152,   149,   150,   149,   150,   149,
     127,   128,   149,   150,   149,   150,   133,   134,   135,   136,
     137,   138,   139,   140,   141,   129,   150,   149,   126,   149,
     127,   148,   126,   150,   109,   152,   448,    69,   109,   837,
     149,   839,   149,   850,   149,    23,   149,    69,   149,   153,
      69,   149,    69,   157,   158,   149,   126,    69,    23,   857,
     126,   149,    32,    32,    32,   487,    69,    32,    50,    32,
     487,    41,    41,    41,   126,   487,    69,   149,    41,   149,
     487,   487,    69,   149,   487,   150,   487,   487,   487,   487,
     150,   487,   487,   149,   126,   487,   722,   487,    69,   151,
     931,   932,   933,   453,   126,   996,   456,   126,    69,   126,
     931,   932,   933,    69,   126,    69,    69,   149,    69,    69,
      69,  1012,   920,    23,   922,   923,   108,   149,    69,    69,
     149,    69,   149,   126,    69,   761,   118,   149,   488,   577,
     578,   579,   580,   581,   582,   583,   584,   585,   586,   603,
     588,   589,   590,   591,   592,   126,   149,   955,   989,   957,
     149,   607,    69,    69,   610,   126,   155,    69,   989,    21,
     968,   797,   594,   126,    69,   126,   126,   126,   149,    69,
      32,   126,   126,   981,    13,   126,   126,   985,   149,    41,
     107,   126,    32,   605,   149,   149,   149,   126,   149,   149,
     149,   126,  1000,   126,   126,   126,   151,   151,   149,   149,
     126,    69,    32,    69,   149,   565,   126,    69,   126,   126,
     570,    41,   151,   126,   126,   126,   151,   126,   151,   151,
     151,   126,  1030,    69,    23,   151,   126,   863,   864,  1037,
     126,   151,   149,   151,    69,   595,   596,   149,   151,    28,
     151,    69,   151,   879,   149,   149,   150,    69,   608,   149,
      32,   611,   700,   695,   115,   151,    41,   129,   130,    41,
     132,   149,   150,   149,   900,   127,   128,   149,   150,   905,
     154,   133,   134,   135,   136,   137,   138,   139,   140,   141,
     127,   128,   129,   130,    89,   132,   148,    69,   150,   925,
     152,   149,   150,   127,   128,   129,   130,   115,   132,   115,
     936,   115,   938,   151,   151,   149,   150,     0,   127,   128,
     129,   130,    69,   132,     7,     8,     9,   151,    23,   129,
     130,    14,   132,   149,   150,    18,    19,    20,   149,    22,
     149,    24,    25,    26,    27,   847,   127,   128,   129,   130,
     149,   132,    69,   979,    37,   127,   128,    40,    69,    23,
      69,   133,   134,   135,   136,   137,   138,   139,   140,   141,
     151,    23,   127,   128,   129,   130,   148,   132,   150,    23,
      63,    64,    65,    23,    23,  1011,    69,   127,   128,   129,
     130,    23,   132,    23,   149,    23,    79,    80,    81,    23,
    1026,    84,    23,    86,    87,    88,    89,    90,    23,    23,
     848,   151,    95,    96,    23,    23,    23,   100,   101,   102,
     103,   127,   128,   129,   130,    23,   132,   110,   111,    23,
     113,    23,    23,    23,   117,    69,     5,   120,   121,   122,
     123,    10,    11,    23,    32,   151,    15,   127,   128,   129,
     130,    23,   132,    41,   127,   128,   129,   130,    23,   132,
      29,    30,    31,    23,    33,    23,    23,    36,   818,    23,
     820,   151,     5,    69,    43,    44,    45,    10,   151,    48,
      49,    69,    51,    52,    53,    54,    55,    56,    57,    58,
      59,    60,    61,    62,    69,   126,    29,    30,    69,    68,
      33,    70,    71,    36,    73,    69,    75,    32,   126,   126,
      43,    23,   150,    69,    83,   151,    41,    23,    69,    69,
      69,   150,    91,    92,   151,    94,   127,   128,   129,   130,
     150,   132,    32,    23,   150,    68,    23,   106,   107,   127,
     128,    41,    75,   112,    69,   133,   134,   135,   136,   137,
     138,   139,   140,   141,   151,   124,   125,    32,    23,    23,
     148,    94,   150,   150,   152,   150,    41,   150,   150,    69,
      23,   150,   150,   106,   150,   127,   128,   129,   130,   112,
     132,   127,   128,   129,   130,    23,   132,   150,   150,    23,
      23,   124,   125,    23,    69,   127,   128,   129,   130,   151,
     132,   126,   127,   128,   115,   151,    23,    23,   133,   134,
     135,   136,   137,   138,   139,   140,   141,    23,    23,   151,
      69,   149,    23,   148,   150,   150,   126,   127,   128,   150,
     150,   150,   150,   133,   134,   135,   136,   137,   138,   139,
     140,   141,   127,   128,   129,   130,   150,   132,   148,   150,
     150,   150,   127,   128,   150,   150,   150,   115,   133,   134,
     135,   136,   137,   138,   139,   140,   141,     7,     8,     9,
      69,   115,    69,   148,    14,   150,    69,   149,    18,    19,
      20,   150,    22,   115,    24,    25,    26,    27,   127,   128,
     129,   130,    69,   132,    23,    41,    69,    37,   149,   149,
      40,   149,   127,   128,   129,   130,    69,   132,   127,   128,
     129,   130,   151,   132,    32,    69,    41,     3,   154,    32,
      41,    41,    41,    63,    64,    65,   151,    41,    14,    69,
      16,    17,   151,    41,    41,    32,    41,    41,    41,    79,
      80,    81,    41,    69,    84,    41,    86,    87,    88,    89,
      90,    41,    38,    39,   150,    95,    96,    43,    32,    69,
     100,   101,   102,   103,    50,    23,   126,    69,    69,   126,
     110,   111,    41,   113,   126,   126,   126,   117,    69,   149,
     120,   121,   122,   123,    69,   149,    72,     3,    74,   126,
      76,   127,   128,   129,   130,    32,   132,   132,    14,    85,
      16,    17,   149,    89,    69,   127,   128,   129,   130,    41,
     132,    97,    98,    99,    41,   151,    41,   103,   104,   105,
     115,   107,    38,    39,    41,   149,    41,    43,   114,   151,
      41,    41,     3,   149,    50,   127,   128,   129,   130,    41,
     132,    41,   115,    14,   132,    16,    17,   127,   128,   129,
     130,   149,   132,    32,   115,   115,    72,    41,    74,   151,
      76,   150,   156,    41,   151,    69,    82,    38,    39,    85,
     151,   151,    43,    89,   151,   127,   128,   129,   130,    50,
     132,    97,    98,    99,   151,   153,   149,   103,   104,   105,
     126,   107,   127,   128,   129,   130,    69,   132,   126,   151,
      69,    72,   126,    74,    69,    76,   127,   128,   129,   130,
     149,   132,    41,    41,    85,    23,   151,   151,    89,   127,
     128,   129,   130,   149,   132,    69,    97,    98,    99,    69,
     151,   132,   103,   104,   105,    69,   107,   127,   128,   129,
     130,   153,   132,   151,   127,   128,   129,   130,   151,   132,
     127,   128,   129,   130,   149,   132,    23,   127,   128,   129,
     130,   151,   132,   151,   127,   128,   129,   130,   151,   132,
     127,   128,   129,   130,   151,   132,    69,   127,   128,   129,
     130,   151,   132,    69,   127,   128,   129,   130,   151,   132,
     127,   128,   129,   130,   151,   132,   149,   127,   128,   129,
     130,   151,   132,   126,   126,   126,   149,    23,   127,   128,
     129,   130,   149,   132,    23,   127,   128,   129,   130,   149,
     132,   127,   128,   129,   130,    23,   132,    41,    41,   149,
     149,   127,   128,   129,   130,   126,   132,   149,   156,   127,
     128,   129,   130,   149,   132,   127,   128,   129,   130,   126,
     132,   126,   149,   149,   127,   128,   129,   130,   153,   132,
     149,   149,    23,   127,   128,   129,   130,   149,   132,   127,
     128,   129,   130,   151,   132,   126,   149,   127,   128,   129,
     130,   126,   132,   149,   149,   149,    41,   127,   128,   129,
     130,   149,   132,   153,    41,   119,   126,   126,   151,   149,
      41,   153,    69,   149,    41,   126,   126,   156,   126,   149,
      42,   277,   850,   132,   792,   448,   605,   603,   263,  1012,
     274,   269,   487,   130,   594,   913,   317,   564,   121,   157,
     522,   855,   610,   110,   639,   480
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,     9,    14,    18,    19,    20,    22,    24,
      25,    26,    27,    37,    40,    63,    64,    65,    69,    79,
      80,    81,    84,    86,    87,    88,    89,    90,    95,    96,
     100,   101,   102,   103,   110,   111,   113,   117,   120,   121,
     122,   123,   160,   161,   162,   163,   164,   165,   166,   167,
     168,   173,   174,   175,   178,   180,   183,   188,   199,   200,
     205,   208,   211,   214,   217,   223,   229,   232,   237,   240,
     243,   246,   247,   250,   251,   253,   254,   255,   258,   259,
     263,   264,   265,   266,   272,   275,   149,   150,   149,   149,
     150,    23,    32,    41,    69,   150,    69,   150,   149,    69,
     149,   150,   220,   149,   149,   149,   149,   149,   150,   149,
     150,   150,   149,    23,   149,   149,   150,   220,   220,   149,
     149,   150,   220,   220,    69,   172,    23,    41,   273,   150,
     150,   220,   149,    23,   149,   150,   149,   150,   149,   150,
     220,   149,   150,   220,   220,    69,   169,    69,   170,    69,
     171,   220,     0,   161,    13,     4,    69,   260,   261,   107,
     212,   213,   279,    32,   149,    41,   149,   149,   157,    69,
     149,   157,    69,    69,   181,   182,    23,     5,    10,    11,
      15,    29,    30,    31,    33,    36,    43,    44,    45,    48,
      49,    51,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    68,    70,    71,    73,    75,    83,    91,
      92,    94,   106,   112,   124,   125,   238,   239,   279,   289,
     297,   298,   299,   300,   301,   302,   303,   304,   305,   306,
     307,   308,   309,   310,   311,   312,   313,   314,   315,   316,
     317,   318,   319,   320,   321,   322,   323,   324,   325,   327,
     328,   333,   334,   335,   336,   337,   338,    69,   126,   149,
      12,    69,   109,   224,   225,   226,    12,    69,   109,   233,
     234,    12,    69,   109,   230,   231,    69,   184,   185,   181,
      28,   179,   189,    50,   108,   118,   281,    66,    77,    93,
     267,   268,   330,   331,   332,    12,   120,   201,   202,    32,
      41,    69,   127,   128,   133,   134,   135,   136,   137,   138,
     139,   140,   141,   148,   150,   176,    69,   248,   249,     3,
      14,    16,    17,    38,    39,    72,    74,    76,    82,    85,
      89,    97,    98,    99,   103,   104,   105,   219,   256,   257,
     278,   279,   280,   281,   282,   283,   284,   285,   286,   287,
     288,   289,   290,   291,   292,   293,   294,   295,   296,   326,
     329,   149,   149,    69,   252,   257,   149,   149,   115,    69,
     126,   149,    41,   149,    32,    41,    69,   127,   128,   133,
     134,   135,   136,   137,   138,   139,   140,   141,   148,   150,
     196,   198,   241,   242,   289,   304,   305,   314,   320,   321,
     322,   323,   324,   325,   333,   334,   335,   241,   149,   201,
     154,   215,   216,   292,   209,   210,   279,   218,   219,   149,
     114,   219,   276,   277,   339,   149,   149,   115,    69,   126,
     149,   115,    69,   126,   149,   115,    69,   126,   149,   149,
     151,    69,   126,   150,   262,    21,   261,    23,   126,   151,
     149,   149,    69,   151,   157,    69,   151,   157,    23,    21,
     182,    69,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,    23,    23,    23,   126,   151,    23,
      69,    69,    69,    21,   225,   126,    69,   126,    69,    21,
     234,    69,   126,    69,    21,   231,   150,    21,   185,    21,
      23,   151,   152,   194,   195,   196,   197,   126,   151,   151,
     151,    23,   126,   151,    69,    69,    21,   202,   150,   176,
     176,   150,   150,   150,   150,   150,   150,   150,   150,   150,
     150,   176,   127,   128,   129,   130,   132,   149,   150,    21,
     249,    23,    23,    23,    23,    23,    23,    23,    23,    23,
      23,    23,    23,    23,   126,   151,   126,   176,    21,    69,
     151,   115,    69,   149,   150,   196,   196,   150,   150,   150,
     150,   150,   150,   150,   150,   150,   150,   196,   127,   128,
     129,   130,   132,   274,   126,   151,   151,    21,    32,    41,
     150,   206,   207,   126,   151,   126,   151,   126,   151,    23,
     126,   151,   115,    69,   115,    69,   115,    69,   149,   150,
      69,    32,    41,    23,    41,   213,   149,   149,   149,   220,
      69,   149,   149,   220,    69,   176,    32,    69,   154,   342,
     343,    41,    41,    41,    41,    41,    41,    41,    32,    32,
      41,    32,    41,    41,    41,    41,    69,    41,   343,   150,
      41,    41,    32,    69,    69,   340,   239,   149,   220,    69,
      23,   126,     6,    32,    34,    41,    42,    69,    78,   116,
     127,   227,   235,   236,   126,   236,   126,   126,   236,   126,
      41,   127,   128,   221,    69,   149,    69,    21,   195,   197,
      23,   149,    35,    47,   186,   187,   293,   193,   149,   149,
      46,    67,   268,    69,   129,   153,   157,   158,   269,   270,
     271,   126,    23,   126,   149,   176,   177,   176,   176,   176,
     176,   176,   176,   176,   176,   176,   176,   151,   176,   176,
     176,   176,   176,   176,    41,    32,    41,    41,    41,    41,
     236,    41,    41,    32,    41,    41,    41,    41,   256,   149,
     220,    69,   149,   126,   176,   149,   220,   115,   221,   196,
     196,   196,   196,   196,   196,   196,   196,   196,   196,   151,
     196,   196,   196,   196,   196,   149,   242,   149,   220,   149,
     220,   176,   149,   155,    32,    41,   126,   150,   216,   149,
     210,   149,   219,   149,   220,    32,   277,   149,   220,   115,
     115,   115,    41,   262,   151,   151,   176,   149,   151,   149,
     151,   149,    41,   341,   155,   341,   156,   244,   245,   153,
     149,    69,    69,    32,    41,   149,   126,   126,    69,   126,
     236,    69,   149,   236,    41,    41,   151,   181,    23,   196,
     126,   151,   149,   191,   190,   126,   149,   150,   271,    69,
     176,    69,    89,   109,   126,   151,   151,   151,   151,   151,
     151,   151,   151,   151,   151,   151,   151,   149,   176,    69,
     149,   149,   151,   151,   151,   151,   151,   151,   151,   151,
     151,   151,   151,   149,   149,   151,   207,   149,    32,    41,
     150,   176,   149,   149,   151,    23,   149,   149,   220,   149,
     220,   153,    69,   126,   151,    69,   228,   236,   235,   236,
     126,   236,   126,   126,   149,    23,    21,   196,   149,   187,
     192,   194,   194,   194,   270,   236,    23,   149,    23,    41,
     203,   176,   176,   149,   149,   176,   176,   151,   262,   176,
     149,   149,    41,   156,   244,   126,   126,   126,   236,   149,
     236,   236,   176,   149,   149,    21,    21,    21,   150,   151,
     176,   176,   153,    41,   126,   149,   149,   149,   151,    23,
     149,   126,   236,   228,   236,   126,   149,   149,   149,   194,
     236,   149,   149,    41,   153,    41,   119,   176,   156,   236,
     126,   126,   236,    21,   151,    41,   153,    32,    69,   127,
     128,   150,   204,   221,   222,   149,    69,   236,   235,   149,
      41,    32,    32,   176,    69,   149,   150,   221,   222,   156,
     126,   126,   151,   176,   236,   228,   151,   126,   236
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
     395,   396,   397,   398,   399,   400,   401,   402,   403,    59,
      40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   159,   160,   160,   161,   161,   161,   161,   161,   161,
     161,   161,   161,   161,   161,   161,   161,   161,   161,   161,
     161,   161,   161,   161,   161,   161,   161,   161,   161,   161,
     161,   161,   161,   161,   161,   161,   161,   161,   161,   161,
     161,   161,   162,   162,   162,   162,   163,   163,   164,   165,
     166,   167,   168,   169,   169,   169,   169,   169,   169,   170,
     170,   170,   170,   170,   170,   171,   171,   171,   171,   171,
     171,   172,   172,   172,   172,   172,   172,   173,   173,   174,
     174,   175,   176,   176,   176,   176,   176,   176,   176,   176,
     176,   176,   176,   176,   176,   176,   176,   176,   176,   176,
     176,   176,   176,   176,   177,   177,   178,   178,   179,   180,
     181,   181,   182,   183,   184,   184,   185,   186,   186,   187,
     187,   187,   189,   188,   190,   188,   191,   188,   192,   188,
     193,   188,   194,   194,   194,   194,   195,   195,   196,   196,
     196,   196,   196,   196,   196,   196,   196,   196,   196,   196,
     196,   196,   196,   196,   196,   196,   196,   196,   196,   197,
     198,   198,   199,   200,   201,   201,   202,   202,   202,   202,
     202,   203,   203,   203,   203,   203,   203,   204,   204,   204,
     204,   204,   204,   204,   204,   205,   206,   206,   207,   207,
     207,   207,   207,   207,   207,   207,   207,   208,   208,   209,
     209,   210,   211,   211,   212,   212,   213,   214,   214,   215,
     215,   216,   217,   217,   217,   217,   218,   218,   219,   219,
     219,   219,   219,   219,   219,   219,   219,   219,   219,   219,
     219,   219,   219,   219,   219,   219,   219,   219,   220,   220,
     220,   220,   220,   220,   221,   221,   221,   222,   222,   222,
     223,   224,   224,   225,   226,   226,   226,   227,   227,   227,
     227,   227,   228,   228,   228,   228,   229,   230,   230,   231,
     231,   231,   232,   233,   233,   234,   234,   234,   235,   235,
     235,   235,   235,   236,   236,   236,   236,   236,   236,   237,
     237,   237,   237,   238,   238,   239,   239,   239,   239,   239,
     239,   239,   239,   239,   239,   239,   239,   239,   239,   239,
     239,   239,   239,   239,   239,   239,   239,   239,   239,   239,
     239,   239,   239,   239,   239,   239,   239,   239,   239,   239,
     239,   239,   239,   239,   240,   240,   241,   241,   242,   242,
     242,   242,   242,   242,   242,   242,   242,   242,   242,   242,
     242,   243,   243,   244,   244,   245,   245,   246,   247,   248,
     248,   249,   250,   251,   252,   252,   252,   252,   253,   254,
     254,   254,   254,   255,   255,   255,   255,   256,   256,   257,
     257,   258,   259,   260,   260,   261,   261,   261,   262,   262,
     262,   263,   263,   264,   264,   264,   264,   264,   264,   265,
     265,   265,   265,   265,   265,   266,   267,   267,   268,   268,
     268,   269,   269,   269,   269,   270,   270,   271,   271,   271,
     271,   271,   273,   274,   272,   275,   275,   275,   275,   276,
     276,   277,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   296,   297,   298,   298,   299,   300,   301,   302,
     303,   304,   304,   305,   306,   307,   308,   309,   310,   310,
     311,   312,   313,   314,   315,   316,   317,   318,   319,   320,
     321,   322,   323,   324,   325,   326,   327,   328,   329,   330,
     330,   331,   332,   333,   334,   335,   336,   337,   338,   339,
     340,   341,   341,   342,   342,   343
  };

  /* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
  const unsigned char
  parser::yyr2_[] =
  {
         0,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     4,     3,     3,
       3,     3,     3,     2,     3,     1,     3,     4,     2,     2,
       3,     1,     3,     4,     2,     2,     3,     1,     3,     4,
       2,     2,     3,     1,     3,     4,     2,     3,     4,     3,
       4,     4,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     2,     2,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     1,     3,     4,     7,     3,     4,
       2,     1,     4,     4,     2,     1,     7,     3,     1,     1,
       1,     1,     0,     5,     0,     8,     0,     8,     0,    10,
       0,     8,     2,     2,     1,     1,     4,     2,     3,     1,
       1,     1,     3,     3,     3,     3,     3,     2,     2,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     5,
       1,     4,     4,     4,     2,     1,     9,     6,     5,     7,
       7,     2,     4,     3,     5,     3,     1,     2,     2,     2,
       1,     1,     1,     4,     3,     6,     3,     1,     5,     3,
       3,     4,     2,     2,     3,     1,     1,     2,     5,     3,
       1,     1,     2,     5,     3,     1,     1,     2,     5,     3,
       1,     1,     2,     5,     3,     6,     3,     1,     1,     1,
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
       1,     1,     1,     3,     3,     3,     1,     3,     3,     3,
       3,     1,     1,     1,     3,     3,     3,     3,     3,     1,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     1,     1,     3,     3,     3,     3,     5,     3,     3,
       3,     3,     1,     3,     3,     3,     1,     1,     1,     1,
       1,     3,     1,     1,     1,     1,     3,     3,     3,     3,
       3,     1,     1,     3,     3,     3,     1,     1,     1,     3,
       3,     1,     3,     2,     2,     2
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
  "LOGLINEAR", "MH_DROP", "MH_INIT_SCALE", "MH_JSCALE", "MH_MODE",
  "MH_NBLOCKS", "MH_REPLIC", "MH_RECOVER", "MODE_CHECK", "MODE_COMPUTE",
  "MODE_FILE", "MODEL", "MODEL_COMPARISON", "MSHOCKS",
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
  "periods", "cutoff", "init_param", "expression", "comma_expression",
  "initval", "initval_option", "endval", "initval_list", "initval_elem",
  "histval", "histval_list", "histval_elem", "model_sparse_options_list",
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
  "osr_params", "osr", "olr", "olr_option", "olr_options", "olr_inst",
  "calib_var", "calib_var_list", "calib_arg1", "calib_arg2", "calib",
  "dynatype", "dynasave", "model_comparison", "model_comparison_options",
  "model_comparison_option", "filename_list", "filename", "filename_elem",
  "planner_objective", "@6", "@7", "ramsey_policy",
  "ramsey_policy_options_list", "ramsey_policy_options", "o_dr_algo",
  "o_solve_algo", "o_simul_algo", "o_linear", "o_order", "o_replic",
  "o_drop", "o_ar", "o_nocorr", "o_nofunctions", "o_nomoments", "o_irf",
  "o_hp_filter", "o_hp_ngrid", "o_periods", "o_cutoff", "o_simul",
  "o_simul_seed", "o_qz_criterium", "o_datafile", "o_nobs", "o_first_obs",
  "o_prefilter", "o_presample", "o_lik_algo", "o_lik_init", "o_nograph",
  "o_conf_sig", "o_mh_replic", "o_mh_drop", "o_mh_jscale", "o_optim",
  "o_mh_init_scale", "o_mode_file", "o_mode_compute", "o_mode_check",
  "o_prior_trunc", "o_mh_mode", "o_mh_nblcks", "o_load_mh_file",
  "o_loglinear", "o_nodiagnostic", "o_bayesian_irf", "o_tex", "o_forecast",
  "o_smoother", "o_moments_varendo", "o_filtered_vars", "o_relative_irf",
  "o_kalman_algo", "o_kalman_tol", "o_olr_beta",
  "o_model_comparison_approximation", "o_print", "o_noprint",
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
       160,     0,    -1,   161,    -1,   160,   161,    -1,   162,    -1,
     173,    -1,   174,    -1,   188,    -1,   178,    -1,   180,    -1,
     183,    -1,   175,    -1,   199,    -1,   200,    -1,   205,    -1,
     208,    -1,   211,    -1,   214,    -1,   217,    -1,   237,    -1,
     240,    -1,   243,    -1,   223,    -1,   232,    -1,   229,    -1,
     246,    -1,   247,    -1,   250,    -1,   163,    -1,   164,    -1,
     251,    -1,   253,    -1,   254,    -1,   259,    -1,   263,    -1,
     264,    -1,   265,    -1,   255,    -1,   258,    -1,   266,    -1,
     272,    -1,   275,    -1,   168,    -1,   165,    -1,   166,    -1,
     167,    -1,    18,    41,   149,    -1,    18,    41,    41,   149,
      -1,   100,   220,   149,    -1,   120,   169,   149,    -1,   121,
     170,   149,    -1,   122,   171,   149,    -1,    88,   172,   149,
      -1,   169,    69,    -1,   169,   126,    69,    -1,    69,    -1,
     169,    69,   115,    -1,   169,   126,    69,   115,    -1,    69,
     115,    -1,   170,    69,    -1,   170,   126,    69,    -1,    69,
      -1,   170,    69,   115,    -1,   170,   126,    69,   115,    -1,
      69,   115,    -1,   171,    69,    -1,   171,   126,    69,    -1,
      69,    -1,   171,    69,   115,    -1,   171,   126,    69,   115,
      -1,    69,   115,    -1,   172,    69,    -1,   172,   126,    69,
      -1,    69,    -1,   172,    69,   115,    -1,   172,   126,    69,
     115,    -1,    69,   115,    -1,    89,    41,   149,    -1,    89,
      23,    41,   149,    -1,    14,    32,   149,    -1,    14,    23,
      32,   149,    -1,    69,    23,   176,   149,    -1,   150,   176,
     151,    -1,    69,    -1,    32,    -1,    41,    -1,   176,   128,
     176,    -1,   176,   127,   176,    -1,   176,   129,   176,    -1,
     176,   130,   176,    -1,   176,   132,   176,    -1,   127,   176,
      -1,   128,   176,    -1,   133,   150,   176,   151,    -1,   134,
     150,   176,   151,    -1,   135,   150,   176,   151,    -1,   136,
     150,   176,   151,    -1,   137,   150,   176,   151,    -1,   138,
     150,   176,   151,    -1,   139,   150,   176,   151,    -1,   140,
     150,   176,   151,    -1,   141,   150,   176,   151,    -1,   148,
     150,   176,   151,    -1,    69,   150,   177,   151,    -1,   176,
      -1,   177,   126,   176,    -1,    40,   149,   181,    21,    -1,
      40,   150,   179,   151,   149,   181,    21,    -1,    28,    23,
      69,    -1,    22,   149,   181,    21,    -1,   181,   182,    -1,
     182,    -1,    69,    23,   176,   149,    -1,    37,   149,   184,
      21,    -1,   184,   185,    -1,   185,    -1,    69,   150,   221,
     151,    23,   176,   149,    -1,   186,   126,   187,    -1,   187,
      -1,    47,    -1,    35,    -1,   293,    -1,    -1,    63,   149,
     189,   194,    21,    -1,    -1,    63,   150,   281,   151,   149,
     190,   194,    21,    -1,    -1,    63,   150,   118,   151,   149,
     191,   194,    21,    -1,    -1,    63,   150,   108,   126,   186,
     151,   192,   149,   194,    21,    -1,    -1,    63,   150,   108,
     151,   193,   149,   194,    21,    -1,   194,   195,    -1,   194,
     197,    -1,   195,    -1,   197,    -1,   196,    23,   196,   149,
      -1,   196,   149,    -1,   150,   196,   151,    -1,   198,    -1,
      32,    -1,    41,    -1,   196,   128,   196,    -1,   196,   127,
     196,    -1,   196,   129,   196,    -1,   196,   130,   196,    -1,
     196,   132,   196,    -1,   127,   196,    -1,   128,   196,    -1,
     133,   150,   196,   151,    -1,   134,   150,   196,   151,    -1,
     135,   150,   196,   151,    -1,   136,   150,   196,   151,    -1,
     137,   150,   196,   151,    -1,   138,   150,   196,   151,    -1,
     139,   150,   196,   151,    -1,   140,   150,   196,   151,    -1,
     141,   150,   196,   151,    -1,   148,   150,   196,   151,    -1,
     152,    69,    23,   196,   149,    -1,    69,    -1,    69,   150,
     221,   151,    -1,   101,   149,   201,    21,    -1,    65,   149,
     201,    21,    -1,   201,   202,    -1,   202,    -1,   120,    69,
     149,    89,   203,   149,   119,   204,   149,    -1,   120,    69,
     149,   109,   176,   149,    -1,   120,    69,    23,   176,   149,
      -1,   120,    69,   126,    69,    23,   176,   149,    -1,    12,
      69,   126,    69,    23,   176,   149,    -1,   203,    41,    -1,
     203,    41,   153,    41,    -1,   203,   126,    41,    -1,   203,
     126,    41,   153,    41,    -1,    41,   153,    41,    -1,    41,
      -1,   204,   222,    -1,   204,   221,    -1,   204,    69,    -1,
     222,    -1,   221,    -1,    69,    -1,   204,   150,   176,   151,
      -1,   150,   176,   151,    -1,   102,    23,   154,   206,   155,
     149,    -1,   206,   149,   207,    -1,   207,    -1,   207,   126,
     150,   176,   151,    -1,   207,   126,    32,    -1,   207,   126,
      41,    -1,   207,   150,   176,   151,    -1,   207,    32,    -1,
     207,    41,    -1,   150,   176,   151,    -1,    32,    -1,    41,
      -1,   110,   149,    -1,   110,   150,   209,   151,   149,    -1,
     209,   126,   210,    -1,   210,    -1,   279,    -1,     9,   149,
      -1,     9,   150,   212,   151,   149,    -1,   212,   126,   213,
      -1,   213,    -1,   279,    -1,   103,   149,    -1,   103,   150,
     215,   151,   149,    -1,   215,   126,   216,    -1,   216,    -1,
     292,    -1,   111,   149,    -1,   111,   150,   218,   151,   149,
      -1,   111,   220,   149,    -1,   111,   150,   218,   151,   220,
     149,    -1,   218,   126,   219,    -1,   219,    -1,   278,    -1,
     279,    -1,   280,    -1,   281,    -1,   282,    -1,   283,    -1,
     284,    -1,   285,    -1,   286,    -1,   287,    -1,   288,    -1,
     289,    -1,   326,    -1,   290,    -1,   291,    -1,   292,    -1,
     293,    -1,   294,    -1,   295,    -1,   296,    -1,   220,    69,
      -1,   220,    69,    23,    69,    -1,   220,   126,    69,    -1,
     220,   126,    69,    23,    69,    -1,    69,    -1,    69,    23,
      69,    -1,   128,    41,    -1,   127,    41,    -1,    41,    -1,
     128,    32,    -1,   127,    32,    -1,    32,    -1,    25,   149,
     224,    21,    -1,   224,   225,    -1,   225,    -1,   226,   126,
     227,   149,    -1,   109,    69,    -1,    69,    -1,    12,    69,
     126,    69,    -1,   235,   126,   228,    -1,   236,   126,   235,
     126,   228,    -1,   236,   126,   236,   126,   236,   126,   235,
     126,   228,    -1,   236,    -1,   236,   126,   236,   126,   236,
      -1,   236,   126,   236,    -1,   236,   126,   236,   126,   236,
      -1,   236,   126,   236,   126,   236,   126,   236,    -1,   236,
     126,   236,   126,   236,   126,   236,   126,   236,    -1,    27,
     149,   230,    21,    -1,   230,   231,    -1,   231,    -1,   109,
      69,   126,   236,   149,    -1,    12,    69,   126,    69,   126,
     236,   149,    -1,    69,   126,   236,   149,    -1,    26,   149,
     233,    21,    -1,   233,   234,    -1,   234,    -1,   109,    69,
     126,   236,   126,   236,   149,    -1,    12,    69,   126,    69,
     126,   236,   126,   236,   149,    -1,    69,   126,   236,   126,
     236,   149,    -1,     6,    -1,    34,    -1,    78,    -1,    42,
      -1,   116,    -1,    -1,    41,    -1,    32,    -1,    69,    -1,
     127,    41,    -1,   127,    32,    -1,    24,   149,    -1,    24,
     150,   238,   151,   149,    -1,    24,   220,   149,    -1,    24,
     150,   238,   151,   220,   149,    -1,   238,   126,   239,    -1,
     239,    -1,   297,    -1,   298,    -1,   299,    -1,   300,    -1,
     301,    -1,   302,    -1,   303,    -1,   304,    -1,   305,    -1,
     306,    -1,   307,    -1,   308,    -1,   309,    -1,   310,    -1,
     311,    -1,   312,    -1,   313,    -1,   314,    -1,   315,    -1,
     316,    -1,   317,    -1,   318,    -1,   319,    -1,   320,    -1,
     289,    -1,   321,    -1,   322,    -1,   323,    -1,   324,    -1,
     325,    -1,   327,    -1,   328,    -1,   333,    -1,   334,    -1,
     335,    -1,   279,    -1,   336,    -1,   337,    -1,   338,    -1,
      95,   150,   241,   151,   149,    -1,    95,   150,   241,   151,
     220,   149,    -1,   241,   126,   242,    -1,   242,    -1,   304,
      -1,   305,    -1,   314,    -1,   320,    -1,   289,    -1,   321,
      -1,   322,    -1,   323,    -1,   324,    -1,   325,    -1,   333,
      -1,   334,    -1,   335,    -1,    96,   150,   241,   151,   149,
      -1,    96,   150,   241,   151,   220,   149,    -1,   156,    69,
     156,   126,   156,    69,   156,    -1,   156,    69,   156,   126,
     236,    -1,   244,    -1,   245,   126,   244,    -1,   123,   220,
     149,    -1,    79,   149,   248,    21,    -1,   248,   249,    -1,
     249,    -1,    69,   150,   176,   151,   149,    -1,   117,   220,
     149,    -1,    84,   149,   252,    21,    -1,   252,    69,   176,
     149,    -1,   252,    69,   126,    69,   176,   149,    -1,    69,
     176,   149,    -1,    69,   126,    69,   176,   149,    -1,    87,
     220,   149,    -1,    86,   149,    -1,    86,   150,   257,   151,
     149,    -1,    86,   220,   149,    -1,    86,   150,   257,   151,
     220,   149,    -1,    80,   149,    -1,    80,   150,   257,   151,
     149,    -1,    80,   220,   149,    -1,    80,   150,   257,   151,
     220,   149,    -1,   329,    -1,   219,    -1,   256,    -1,   257,
     126,   256,    -1,    81,   220,   149,    -1,     8,   149,   260,
      21,    -1,   260,   261,    -1,   261,    -1,    69,   262,    23,
     176,   149,    -1,    69,   126,    69,   262,    23,   176,   149,
      -1,     4,    69,   150,    41,   151,   262,    23,   176,   149,
      -1,    -1,   150,    41,   151,    -1,   150,    32,   151,    -1,
       7,   149,    -1,     7,   150,    13,   151,   149,    -1,    20,
     150,    69,   151,   149,    -1,    20,   150,    69,   151,   220,
     149,    -1,    20,    69,   149,    -1,    20,   150,    69,   157,
      69,   151,   149,    -1,    20,   150,    69,   157,    69,   151,
     220,   149,    -1,    20,    69,   157,    69,   149,    -1,    19,
     150,    69,   151,   149,    -1,    19,   150,    69,   151,   220,
     149,    -1,    19,    69,   149,    -1,    19,   150,    69,   157,
      69,   151,   149,    -1,    19,   150,    69,   157,    69,   151,
     220,   149,    -1,    19,    69,   157,    69,   149,    -1,    64,
     150,   267,   151,   269,   149,    -1,   267,   126,   268,    -1,
     268,    -1,   330,    -1,   331,    -1,   332,    -1,   270,    -1,
     269,   126,   270,    -1,   270,   150,   236,   151,    -1,   269,
     126,   270,   150,   236,   151,    -1,   271,    -1,   270,   271,
      -1,    69,    -1,   158,    -1,   129,    -1,   153,    -1,   157,
      -1,    -1,    -1,    90,   273,   196,   274,   149,    -1,   113,
     149,    -1,   113,   150,   276,   151,   149,    -1,   113,   220,
     149,    -1,   113,   150,   276,   151,   220,   149,    -1,   276,
     126,   277,    -1,   277,    -1,   219,    -1,   339,    -1,    16,
      23,    41,    -1,   107,    23,    41,    -1,   104,    23,    41,
      -1,    50,    -1,    85,    23,    41,    -1,    99,    23,    41,
      -1,    17,    23,    41,    -1,     3,    23,    41,    -1,    72,
      -1,    74,    -1,    76,    -1,    43,    23,    41,    -1,    38,
      23,    41,    -1,    39,    23,    41,    -1,    89,    23,    41,
      -1,    14,    23,    32,    -1,   103,    -1,   105,    23,    41,
      -1,    97,    23,    41,    -1,    97,    23,    32,    -1,    15,
      23,    69,    -1,    70,    23,   343,    -1,    70,    23,    41,
      -1,    31,    23,    41,    -1,    91,    23,    41,    -1,    92,
      23,    41,    -1,    48,    23,    41,    -1,    49,    23,    41,
      -1,    75,    -1,    36,    -1,    10,    23,    32,    -1,    58,
      23,    41,    -1,    53,    23,    32,    -1,    55,    23,    32,
      -1,    83,    23,   150,   245,   151,    -1,    54,    23,    32,
      -1,    54,    23,    41,    -1,    62,    23,    69,    -1,    61,
      23,    41,    -1,    60,    -1,    94,    23,    32,    -1,    56,
      23,    41,    -1,    57,    23,    41,    -1,    51,    -1,    52,
      -1,    73,    -1,     5,    -1,   112,    -1,    33,    23,    41,
      -1,   106,    -1,    68,    -1,    30,    -1,    98,    -1,    44,
      23,    41,    -1,    45,    23,    41,    -1,    82,    23,   236,
      -1,    66,    23,    46,    -1,    66,    23,    67,    -1,    93,
      -1,    77,    -1,   124,    23,    69,    -1,   125,    23,   340,
      -1,    29,    23,   343,    -1,    11,    -1,    71,    -1,    59,
      -1,   114,    23,    32,    -1,    69,   153,    69,    -1,    41,
      -1,    41,   153,    41,    -1,   154,   341,    -1,   342,   341,
      -1,   342,   155,    -1
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
      80,    82,    84,    86,    88,    90,    92,    96,   101,   105,
     109,   113,   117,   121,   124,   128,   130,   134,   139,   142,
     145,   149,   151,   155,   160,   163,   166,   170,   172,   176,
     181,   184,   187,   191,   193,   197,   202,   205,   209,   214,
     218,   223,   228,   232,   234,   236,   238,   242,   246,   250,
     254,   258,   261,   264,   269,   274,   279,   284,   289,   294,
     299,   304,   309,   314,   319,   321,   325,   330,   338,   342,
     347,   350,   352,   357,   362,   365,   367,   375,   379,   381,
     383,   385,   387,   388,   394,   395,   404,   405,   414,   415,
     426,   427,   436,   439,   442,   444,   446,   451,   454,   458,
     460,   462,   464,   468,   472,   476,   480,   484,   487,   490,
     495,   500,   505,   510,   515,   520,   525,   530,   535,   540,
     546,   548,   553,   558,   563,   566,   568,   578,   585,   591,
     599,   607,   610,   615,   619,   625,   629,   631,   634,   637,
     640,   642,   644,   646,   651,   655,   662,   666,   668,   674,
     678,   682,   687,   690,   693,   697,   699,   701,   704,   710,
     714,   716,   718,   721,   727,   731,   733,   735,   738,   744,
     748,   750,   752,   755,   761,   765,   772,   776,   778,   780,
     782,   784,   786,   788,   790,   792,   794,   796,   798,   800,
     802,   804,   806,   808,   810,   812,   814,   816,   818,   821,
     826,   830,   836,   838,   842,   845,   848,   850,   853,   856,
     858,   863,   866,   868,   873,   876,   878,   883,   887,   893,
     903,   905,   911,   915,   921,   929,   939,   944,   947,   949,
     955,   963,   968,   973,   976,   978,   986,   996,  1003,  1005,
    1007,  1009,  1011,  1013,  1014,  1016,  1018,  1020,  1023,  1026,
    1029,  1035,  1039,  1046,  1050,  1052,  1054,  1056,  1058,  1060,
    1062,  1064,  1066,  1068,  1070,  1072,  1074,  1076,  1078,  1080,
    1082,  1084,  1086,  1088,  1090,  1092,  1094,  1096,  1098,  1100,
    1102,  1104,  1106,  1108,  1110,  1112,  1114,  1116,  1118,  1120,
    1122,  1124,  1126,  1128,  1130,  1136,  1143,  1147,  1149,  1151,
    1153,  1155,  1157,  1159,  1161,  1163,  1165,  1167,  1169,  1171,
    1173,  1175,  1181,  1188,  1196,  1202,  1204,  1208,  1212,  1217,
    1220,  1222,  1228,  1232,  1237,  1242,  1249,  1253,  1259,  1263,
    1266,  1272,  1276,  1283,  1286,  1292,  1296,  1303,  1305,  1307,
    1309,  1313,  1317,  1322,  1325,  1327,  1333,  1341,  1351,  1352,
    1356,  1360,  1363,  1369,  1375,  1382,  1386,  1394,  1403,  1409,
    1415,  1422,  1426,  1434,  1443,  1449,  1456,  1460,  1462,  1464,
    1466,  1468,  1470,  1474,  1479,  1486,  1488,  1491,  1493,  1495,
    1497,  1499,  1501,  1502,  1503,  1509,  1512,  1518,  1522,  1529,
    1533,  1535,  1537,  1539,  1543,  1547,  1551,  1553,  1557,  1561,
    1565,  1569,  1571,  1573,  1575,  1579,  1583,  1587,  1591,  1595,
    1597,  1601,  1605,  1609,  1613,  1617,  1621,  1625,  1629,  1633,
    1637,  1641,  1643,  1645,  1649,  1653,  1657,  1661,  1667,  1671,
    1675,  1679,  1683,  1685,  1689,  1693,  1697,  1699,  1701,  1703,
    1705,  1707,  1711,  1713,  1715,  1717,  1719,  1723,  1727,  1731,
    1735,  1739,  1741,  1743,  1747,  1751,  1755,  1757,  1759,  1761,
    1765,  1769,  1771,  1775,  1778,  1781
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    84,    84,    85,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   131,   132,   133,   134,   138,   139,   142,   145,
     149,   153,   157,   161,   163,   165,   167,   169,   171,   176,
     178,   180,   182,   184,   186,   191,   193,   195,   197,   199,
     201,   206,   208,   210,   212,   214,   216,   221,   225,   232,
     236,   243,   248,   250,   252,   254,   256,   258,   260,   262,
     264,   266,   268,   270,   272,   274,   276,   278,   280,   282,
     284,   286,   288,   290,   295,   297,   301,   303,   308,   312,
     317,   318,   322,   327,   332,   333,   337,   341,   342,   346,
     347,   348,   352,   352,   353,   353,   355,   355,   357,   357,
     359,   359,   364,   365,   366,   367,   371,   373,   378,   379,
     380,   382,   384,   386,   388,   390,   392,   394,   396,   398,
     400,   402,   404,   406,   408,   410,   412,   414,   416,   420,
     424,   426,   431,   435,   439,   440,   444,   446,   448,   450,
     452,   457,   459,   461,   463,   465,   467,   473,   475,   477,
     479,   481,   483,   485,   487,   492,   497,   499,   504,   506,
     508,   510,   512,   514,   516,   518,   520,   525,   529,   533,
     534,   537,   541,   543,   547,   548,   551,   555,   557,   561,
     562,   565,   569,   571,   573,   575,   579,   580,   583,   584,
     585,   586,   587,   588,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   599,   600,   601,   602,   606,   608,
     610,   612,   614,   616,   621,   623,   625,   630,   632,   634,
     639,   644,   646,   651,   655,   660,   665,   675,   680,   686,
     696,   701,   712,   718,   726,   736,   750,   754,   756,   760,
     767,   776,   785,   789,   791,   795,   804,   815,   827,   829,
     831,   833,   835,   840,   841,   842,   843,   844,   846,   853,
     855,   857,   859,   864,   865,   868,   869,   870,   871,   872,
     873,   874,   875,   876,   877,   878,   879,   880,   881,   882,
     883,   884,   885,   886,   887,   888,   889,   890,   891,   892,
     893,   894,   895,   896,   897,   898,   899,   900,   901,   902,
     903,   904,   905,   906,   910,   912,   917,   918,   922,   923,
     924,   925,   926,   927,   928,   929,   930,   931,   932,   933,
     934,   938,   940,   945,   946,   950,   951,   955,   960,   965,
     966,   969,   973,   976,   980,   982,   984,   986,   990,   993,
     994,   995,   996,   999,  1000,  1001,  1002,  1005,  1006,  1009,
    1010,  1013,  1016,  1020,  1021,  1024,  1025,  1026,  1029,  1030,
    1031,  1034,  1035,  1038,  1039,  1040,  1041,  1042,  1043,  1045,
    1046,  1047,  1048,  1049,  1050,  1052,  1056,  1057,  1060,  1061,
    1062,  1065,  1066,  1067,  1068,  1071,  1072,  1075,  1076,  1077,
    1078,  1079,  1082,  1082,  1082,  1085,  1087,  1089,  1091,  1096,
    1097,  1100,  1101,  1104,  1105,  1106,  1107,  1108,  1109,  1110,
    1111,  1112,  1113,  1114,  1115,  1116,  1117,  1118,  1119,  1120,
    1121,  1122,  1123,  1125,  1126,  1127,  1129,  1130,  1131,  1132,
    1133,  1134,  1135,  1136,  1137,  1138,  1139,  1140,  1141,  1142,
    1143,  1144,  1145,  1146,  1147,  1148,  1149,  1150,  1151,  1152,
    1153,  1154,  1155,  1156,  1157,  1158,  1159,  1160,  1161,  1163,
    1165,  1168,  1169,  1170,  1171,  1172,  1173,  1174,  1175,  1176,
    1178,  1186,  1187,  1191,  1192,  1201
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
       2,     2,     2,     2,     2,   152,     2,     2,     2,   156,
     150,   151,     2,     2,     2,     2,   157,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   153,   149,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   154,   158,   155,     2,     2,     2,     2,     2,     2,
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
     145,   146,   147,   148
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 1465;
  const int parser::yynnts_ = 185;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 152;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 159;

  const unsigned int parser::yyuser_token_number_max_ = 403;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1203 "DynareBison.yy"


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

