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
#line 301 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 113:
#line 303 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 114:
#line 306 "DynareBison.yy"
    { driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 115:
#line 308 "DynareBison.yy"
    { driver.end_endval(); ;}
    break;

  case 118:
#line 314 "DynareBison.yy"
    { driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 119:
#line 316 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 122:
#line 322 "DynareBison.yy"
    { driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 125:
#line 329 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 126:
#line 331 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 127:
#line 333 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 130:
#line 338 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 131:
#line 339 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 132:
#line 340 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 133:
#line 341 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 134:
#line 342 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 135:
#line 343 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 136:
#line 345 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 137:
#line 346 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 138:
#line 347 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 139:
#line 348 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 144:
#line 358 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 145:
#line 360 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val)); ;}
    break;

  case 146:
#line 364 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 148:
#line 367 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 149:
#line 369 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 150:
#line 371 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 151:
#line 373 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 152:
#line 375 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 153:
#line 377 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 154:
#line 379 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 155:
#line 381 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 156:
#line 383 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 157:
#line 385 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 158:
#line 387 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 159:
#line 389 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 160:
#line 391 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 161:
#line 393 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 162:
#line 395 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 163:
#line 397 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 164:
#line 399 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 165:
#line 401 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 166:
#line 403 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 167:
#line 405 "DynareBison.yy"
    { (yyval.node_val) = driver.add_dummy((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 168:
#line 407 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 169:
#line 409 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 170:
#line 413 "DynareBison.yy"
    { driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 171:
#line 416 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 172:
#line 418 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 173:
#line 421 "DynareBison.yy"
    { driver.end_shocks(); ;}
    break;

  case 174:
#line 423 "DynareBison.yy"
    { driver.end_mshocks(); ;}
    break;

  case 177:
#line 430 "DynareBison.yy"
    { driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val)); ;}
    break;

  case 178:
#line 432 "DynareBison.yy"
    { driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 179:
#line 434 "DynareBison.yy"
    { driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 180:
#line 436 "DynareBison.yy"
    { driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 181:
#line 438 "DynareBison.yy"
    { driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 182:
#line 442 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 183:
#line 444 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 184:
#line 446 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 185:
#line 448 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 186:
#line 450 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 187:
#line 452 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 188:
#line 456 "DynareBison.yy"
    { driver.add_value((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 189:
#line 458 "DynareBison.yy"
    { driver.add_value((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 190:
#line 461 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 191:
#line 464 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 192:
#line 466 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 193:
#line 470 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 194:
#line 472 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 195:
#line 474 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 196:
#line 476 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 197:
#line 478 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 198:
#line 480 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 199:
#line 482 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 200:
#line 484 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 201:
#line 486 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 202:
#line 490 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 203:
#line 492 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 207:
#line 502 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 208:
#line 504 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 212:
#line 514 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 213:
#line 516 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 218:
#line 528 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 219:
#line 530 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 220:
#line 532 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 221:
#line 534 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 247:
#line 567 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 248:
#line 569 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 249:
#line 571 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 250:
#line 573 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 251:
#line 575 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 252:
#line 577 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 253:
#line 581 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 254:
#line 583 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 255:
#line 585 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 256:
#line 589 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 257:
#line 591 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 258:
#line 593 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 259:
#line 596 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 260:
#line 599 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 261:
#line 601 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 263:
#line 607 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 264:
#line 613 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 265:
#line 619 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 266:
#line 629 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 267:
#line 634 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 268:
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

  case 269:
#line 652 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 270:
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

  case 271:
#line 668 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 272:
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

  case 273:
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

  case 274:
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

  case 275:
#line 710 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 276:
#line 713 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 277:
#line 715 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 278:
#line 719 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 279:
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

  case 280:
#line 737 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 281:
#line 747 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 282:
#line 750 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 283:
#line 752 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 284:
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

  case 285:
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

  case 286:
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

  case 287:
#line 790 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 288:
#line 792 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 289:
#line 794 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 290:
#line 796 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 291:
#line 798 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 292:
#line 801 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 296:
#line 806 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 297:
#line 808 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 298:
#line 812 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 299:
#line 814 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 300:
#line 816 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 301:
#line 818 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 343:
#line 867 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 344:
#line 869 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 360:
#line 892 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 361:
#line 894 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 362:
#line 898 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 363:
#line 900 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 366:
#line 907 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 367:
#line 909 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 370:
#line 915 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 371:
#line 917 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 372:
#line 919 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 373:
#line 922 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 374:
#line 924 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 375:
#line 926 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 376:
#line 928 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 377:
#line 931 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 378:
#line 934 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 379:
#line 936 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 380:
#line 938 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 381:
#line 940 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 382:
#line 943 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 385:
#line 950 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 386:
#line 952 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 387:
#line 954 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 388:
#line 957 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 389:
#line 959 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 390:
#line 961 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 391:
#line 965 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 392:
#line 967 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 393:
#line 971 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 394:
#line 973 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 395:
#line 975 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 396:
#line 977 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 397:
#line 979 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 398:
#line 981 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 399:
#line 985 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 400:
#line 987 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 401:
#line 989 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 402:
#line 991 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 403:
#line 993 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 404:
#line 995 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 405:
#line 999 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 411:
#line 1011 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 412:
#line 1013 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 413:
#line 1015 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 414:
#line 1017 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 415:
#line 1021 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 416:
#line 1023 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 418:
#line 1028 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 419:
#line 1030 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 420:
#line 1032 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 421:
#line 1034 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 422:
#line 1037 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 423:
#line 1038 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 425:
#line 1041 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 426:
#line 1043 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 427:
#line 1045 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 428:
#line 1047 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 452:
#line 1084 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 453:
#line 1086 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 460:
#line 1100 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 461:
#line 1102 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 462:
#line 1106 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 463:
#line 1108 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 495:
#line 1148 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 496:
#line 1149 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 497:
#line 1150 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 498:
#line 1151 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 499:
#line 1152 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 500:
#line 1153 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 501:
#line 1154 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 502:
#line 1155 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 503:
#line 1156 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 504:
#line 1157 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 505:
#line 1158 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 506:
#line 1159 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 507:
#line 1160 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 508:
#line 1161 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 509:
#line 1163 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 510:
#line 1164 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 511:
#line 1165 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 512:
#line 1166 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 513:
#line 1167 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 514:
#line 1168 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 515:
#line 1169 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 516:
#line 1171 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 517:
#line 1173 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 518:
#line 1175 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 519:
#line 1176 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 520:
#line 1177 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 521:
#line 1178 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 522:
#line 1179 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 523:
#line 1181 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 524:
#line 1183 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 525:
#line 1185 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 526:
#line 1186 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 527:
#line 1187 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 528:
#line 1188 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1190 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1191 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1192 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 533:
#line 1193 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 534:
#line 1194 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 535:
#line 1195 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 536:
#line 1196 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 537:
#line 1197 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 538:
#line 1198 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 539:
#line 1199 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 540:
#line 1200 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 541:
#line 1201 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 542:
#line 1202 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 543:
#line 1203 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 544:
#line 1204 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 545:
#line 1205 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 546:
#line 1206 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 547:
#line 1207 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 548:
#line 1208 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 549:
#line 1210 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 550:
#line 1212 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 551:
#line 1214 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 552:
#line 1215 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 553:
#line 1216 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 554:
#line 1217 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 555:
#line 1218 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 556:
#line 1219 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 557:
#line 1220 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 558:
#line 1221 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 559:
#line 1222 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 560:
#line 1224 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 561:
#line 1225 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 562:
#line 1226 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 563:
#line 1227 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 564:
#line 1228 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1229 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 566:
#line 1230 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 567:
#line 1231 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 568:
#line 1233 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 569:
#line 1234 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 570:
#line 1235 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 571:
#line 1236 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1237 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1238 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1239 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 575:
#line 1240 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 576:
#line 1241 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 577:
#line 1242 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1243 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 579:
#line 1244 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1245 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 581:
#line 1246 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 582:
#line 1247 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1248 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1249 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1250 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1251 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1253 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1254 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1260 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1261 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1265 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 592:
#line 1266 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 593:
#line 1267 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 594:
#line 1268 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 595:
#line 1272 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 597:
#line 1281 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); 
               ;}
    break;

  case 598:
#line 1290 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 599:
#line 1292 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 600:
#line 1300 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2354 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1068;
  const short int
  parser::yypact_[] =
  {
      1236,    26,    45,   -42,   -77,   288,    98,    69,    28,    64,
     -19,    34,    81,    85,   103,   119,   335,   106,   370,   -56,
     221,   111,   242,   246,    36,   374,   375,    78, -1068,   267,
     274,   374,   301,   479,   420,   433,    75,    77,   374,   437,
     454,   458,   374,   447,  1117, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
     337,  1472,   351,  1452, -1068,   535,    72, -1068,   436,   517,
     378,    48,   -79,   482,   259,   494,   497,   570, -1068,  1352,
     180,    52,   372,   386,   537,   497,   583,   581,   443, -1068,
     416,   538,    49,   937,   561,   562, -1068,  1576,   182,   206,
     525,   210,   602,   462,   963,   306,   306,   217,    49,   459,
   -1068,    58, -1068,   436, -1068,  1576,   218, -1068,  1523,   220,
     225,   531,   266,   533,   269,   534,   271,   281, -1068,  1562,
   -1068, -1068, -1068,   629, -1068,   630,   632,   633,   634,   636,
   -1068,   641,   643,   648, -1068,   658,   660,   661,   667, -1068,
     550,   510, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,   674,
     675,   676, -1068,   576,   516, -1068, -1068, -1068,   520,   645,
     -37,    71, -1068,   696,   -39, -1068, -1068,   542, -1068,   543,
   -1068, -1068,   649,   -62, -1068,   656,   118,   707,    59, -1068,
     679, -1068,   709, -1068, -1068,   731,   732,   734,   736,   738,
   -1068, -1068,   739,   744,   746,   747,   748,   750, -1068, -1068,
     758,   763, -1068, -1068, -1068,   765,   766, -1068, -1068,   -34,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
     767,   702, -1068,   719, -1068,   723,   400, -1068,   669,   735,
     681,   740,   401, -1068,   742,   693,   743,   405, -1068,   638,
     204, -1068,   384,   810,   650,   653, -1068,   688, -1068,    -4,
     652,   655,   817, -1068, -1068,     5, -1068, -1068, -1068, -1068,
     785,   786,   105, -1068,   677, -1068, -1068,   678,   680,   683,
     937,   937,   684,   685,   686,   687,   691,   692,   694,   700,
     704,   706,   937,   545,   711,   404, -1068,   800,   415,   841,
     862,   874,   875,   881,   884, -1068, -1068, -1068,   885,   886,
     887, -1068,   888, -1068,   893,   901,   741, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068,   814,   854, -1068,   761, -1068,   749, -1068,
   -1068,   764,   774,   787,   963,   963,   788,   789,   793,   794,
     796,   797,   799,   801,   803,   804,   963,   558, -1068,    12,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068,   166, -1068,   140,    30,   249, -1068,
   -1068, -1068,   268, -1068, -1068,   296, -1068, -1068,   921, -1068,
     300, -1068, -1068, -1068, -1068, -1068,   832,   878, -1068, -1068,
     873,   920, -1068, -1068,   877,   922, -1068, -1068,   973,   974,
     975,   983,   985,   987,   992,  1000,  1007,  1009,  1010,  1011,
    1013,  1014,  1015,  1016,  1017,  1020,  1021,  1022,  1023,  1041,
    1045,  1047,  1057,  1058,  1060,   900,   942, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068,   101,   314,   101,  1044,   314,  1048,
    1018,  1052,    20,  1053,  1066,  1036,  1037,  1472,  1069,  1070,
     101,  1074,  1452,  1081,   905,   949,  1061,   233,  1065, -1068,
   -1068,  1094,   436,   956, -1068, -1068,   958,    80,  1073,   969,
      84,  1083,   937, -1068, -1068, -1068,   966,  1115,  1122,  1124,
    1132,  1137,   101,   101,   101,  1139,  1140,  1145,  1146,  1086,
    1005,   101,  1352,    86,  1138,  1171,  1082, -1068, -1068, -1068,
     199,  1084,   346,  1085, -1068, -1068,  1089,   346,  1090, -1068,
   -1068,   350, -1068, -1068, -1068,  1151,  1042, -1068,  1153,    37,
   -1068,    54, -1068,   575, -1068,  1054,  1055,    50,   538,    65,
    1099,    33, -1068, -1068,   937,   937,   937,   937,  1101,   465,
     937,   937,   937,   937,   937,   937,   937,   937,   937,   937,
     428,   937,   937,   937,   937,   937, -1068,   937, -1068, -1068,
    1155,   751, -1068,   824,  1204,   101,  1206,  1207,  1209,  1210,
    1213,  1216,   101,  1218,  1223,  1224,    87, -1068,  1159, -1068,
     937,   937,   937,   350,  1133,   493,   963,   963,   963,   963,
     963,   963,   963,   963,   963,   963,   434,   963,   963,   963,
     963,   963,  1095,   306,    95,   109, -1068, -1068, -1068,   937,
     247,    41,    58,  1100,   436,  1102,  1576,   121,   101,  1523,
     198, -1068,  1175, -1068,  1176, -1068,  1183,  1256,  1259,  1260,
    1264,  1265,  1266,  1272,  1281,  1287,  1288,  1289,  1290,  1294,
    1295,  1301,   101,   101,  1302,   966,   101,   101,  1303,  1304,
     101,  1310,   101,   101,  1170,  1562, -1068, -1068, -1068, -1068,
    1321,  1323, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
    1315,    18, -1068, -1068, -1068, -1068,  1177, -1068, -1068,  1182,
   -1068, -1068, -1068, -1068,  1185, -1068,  1327,  1186,  1192,  1194,
     937, -1068, -1068, -1068, -1068, -1068,   282,  1195, -1068, -1068,
     293,  1196,  1261, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068,  1198, -1068, -1068,
   -1068,   303, -1068,  1298,  1312, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068,   264,  1211,  1271,  1275,  1333,  1279,   346,
    1342,  1217,   346, -1068,  1384,  1389,  1248, -1068,   497,  1410,
   -1068, -1068, -1068,   963, -1068, -1068, -1068,  1411, -1068,   311,
   -1068, -1068, -1068,  1253, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068,   279,    97, -1068,  1369,   937,  1374,
       8,   509,   891,   930,   616,  1274,   532,   539,   610,   635,
     647,   654,   715,   722,   762,   772, -1068,   465,   465,  1101,
    1101,  1316,   784,   937, -1068,  1391,  1309, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
     307, -1068,   790,  1141,  1291,  1280,   872,   883,   898,   988,
     998,  1019,  1031,  1038,  1046,  1062, -1068,   493,   493,  1133,
    1133,  1316, -1068, -1068, -1068,   310, -1068,   312,  1068,    30,
    1285, -1068, -1068,    38,   937, -1068, -1068, -1068, -1068, -1068,
   -1068,   318, -1068, -1068, -1068,   332, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
    1262, -1068, -1068, -1068,  1396, -1068, -1068,  1296,  1446, -1068,
   -1068,  1359, -1068,   207, -1068,   212, -1068,  1407, -1068,   328,
   -1068, -1068, -1068, -1068, -1068, -1068,   346,   199,  1354,   346,
    1357,  1371, -1068,  1311, -1068, -1068,  1477,   419,   963,  1366,
     101,   575, -1068,   688,   688,   688,    65, -1068,   346, -1068,
    1478,  1376,  1479,  1467,   937, -1068,   937,   937, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,  1328,
    1388,   937, -1068, -1068, -1068,   937,   937, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068,    41, -1068, -1068, -1068,   937,  1150, -1068, -1068,  1470,
   -1068,  1186,   937, -1068, -1068,   338, -1068,   342,  1324,  1198,
   -1068, -1068,  1385,  1387,  1393,   346,  1340,   346,   346, -1068,
     937, -1068,  1400, -1068, -1068, -1068,  1341,    63,   232,   573,
     165,  1348,   937, -1068,   937,  1338,    61,  1423,  1157,  1165,
   -1068, -1068,  1455,  1181,  1187,  1244, -1068, -1068,  1502,  1511,
   -1068, -1068,  1408, -1068,   346,   346,   346,  1413, -1068,  1355,
    1360,  1517, -1068,   688, -1068, -1068, -1068,   346, -1068,  1530,
    1557,  1493,  1363,  1509,  1435, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068,   937, -1068,    25,  1416, -1068,  1430,   346, -1068,
   -1068, -1068,   664,  1379, -1068, -1068, -1068,  1519,  1381,   937,
    1567,  1497, -1068,   346,    73,  1390, -1068, -1068, -1068,  1533,
     616,   913, -1068,  1382,  1448,  1449, -1068, -1068, -1068,   616,
   -1068,   346,   346,  1453, -1068,   346, -1068
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   422,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     4,    29,    30,    45,
      46,    47,    44,     5,     6,     7,    12,     9,    10,    11,
       8,    13,    14,    15,    16,    17,    18,    19,    23,    25,
      24,    20,    21,    22,    26,    27,    28,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,   391,     0,     0,   207,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   251,   298,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   130,
       0,     0,     0,     0,     0,     0,   378,     0,     0,     0,
      75,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     212,     0,   202,     0,   218,     0,     0,   425,     0,     0,
       0,    57,     0,    63,     0,    69,     0,     0,   462,     0,
       1,     3,   452,     0,   565,     0,     0,     0,     0,     0,
     556,     0,     0,     0,   557,     0,     0,     0,     0,   440,
     451,     0,   441,   446,   444,   447,   445,   442,   443,   448,
     449,   433,   434,   435,   436,   437,   438,   439,   460,     0,
       0,     0,   454,   459,     0,   456,   455,   457,     0,     0,
     388,     0,   384,     0,     0,   210,   211,     0,    81,     0,
      48,   401,     0,     0,   395,     0,     0,     0,     0,   117,
       0,   540,     0,   545,   524,     0,     0,     0,     0,     0,
     537,   538,     0,     0,     0,     0,     0,     0,   558,   533,
       0,     0,   544,   539,   523,     0,     0,   543,   541,     0,
     303,   339,   328,   304,   305,   306,   307,   308,   309,   310,
     311,   312,   313,   314,   315,   316,   317,   318,   319,   320,
     321,   322,   323,   324,   325,   326,   327,   329,   330,   331,
     332,   333,   334,   335,   336,   337,   338,   340,   341,   342,
     247,     0,   300,     0,   264,     0,     0,   261,     0,     0,
       0,     0,     0,   283,     0,     0,     0,     0,   277,     0,
       0,   121,     0,     0,     0,     0,    83,     0,   498,     0,
       0,     0,     0,   552,   551,     0,   407,   408,   409,   410,
       0,     0,     0,   176,     0,    88,    89,     0,     0,    87,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   369,     0,     0,     0,
       0,     0,     0,     0,     0,   503,   504,   505,     0,     0,
       0,   546,     0,   512,     0,     0,     0,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,   236,   238,
     239,   240,   241,   242,   243,   244,   235,   237,   245,   246,
     380,   377,    78,    73,     0,    54,     0,    79,     0,   148,
     149,     0,     0,   171,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   423,   147,     0,
     346,   351,   347,   348,   349,   350,   352,   353,   354,   355,
     356,   357,   358,   359,     0,    50,     0,     0,     0,   215,
     216,   217,     0,   205,   206,     0,   223,   220,     0,   431,
       0,   430,   432,   427,   371,    60,    55,     0,    51,    66,
      61,     0,    52,    72,    67,     0,    53,   366,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   465,   466,   467,   468,
     469,   470,   471,   472,   473,   474,   475,   476,   477,   478,
     479,   480,   481,   482,   483,   492,   484,   485,   486,   487,
     488,   489,   490,   491,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   382,
     383,     0,     0,     0,    82,    49,     0,     0,     0,     0,
       0,     0,     0,   115,   116,   252,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   249,     0,   263,   259,   260,
     292,     0,   292,     0,   281,   282,     0,   292,     0,   275,
     276,     0,   119,   120,   112,     0,     0,    84,     0,     0,
     142,     0,   143,     0,   138,     0,     0,     0,     0,     0,
       0,     0,   174,   175,     0,     0,     0,     0,    95,    96,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    85,     0,   367,   368,
       0,     0,   372,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    76,    74,    80,
       0,     0,     0,     0,   155,   156,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   173,   200,   201,     0,
       0,   192,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    58,    56,    64,    62,    70,    68,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   494,   493,   561,   258,
       0,     0,   562,   563,   564,   560,   566,   515,   518,   517,
       0,     0,   516,   519,   520,   553,     0,   554,   450,     0,
     567,   525,   542,   458,     0,   392,     0,   388,     0,     0,
       0,   496,   209,   208,   404,   399,     0,     0,   398,   393,
       0,     0,     0,   555,   506,   547,   548,   521,   522,   527,
     530,   528,   535,   536,   526,   532,   531,     0,   534,   302,
     299,     0,   248,     0,     0,   287,   294,   288,   293,   290,
     295,   289,   291,     0,     0,     0,   269,     0,     0,   292,
       0,     0,   292,   255,     0,     0,     0,   114,     0,     0,
     131,   140,   141,     0,   145,   126,   125,     0,   127,     0,
     124,   128,   129,     0,   134,   132,   549,   550,   406,   417,
     419,   420,   421,   418,     0,   411,   415,     0,     0,     0,
       0,     0,     0,     0,   111,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    86,    91,    90,    92,
      93,    94,     0,     0,   375,     0,     0,   502,   510,   495,
     501,   507,   508,   499,   509,   514,   500,   497,   513,   379,
       0,    77,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   146,   151,   150,   152,
     153,   154,   424,   345,   343,     0,   360,     0,     0,     0,
       0,   197,   198,     0,     0,   214,   213,   204,   203,   222,
     219,     0,   559,   429,   426,     0,    59,    65,    71,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,   578,
     579,   580,   581,   582,   583,   584,   585,   586,   587,   588,
     589,   590,   591,   592,   593,   594,   463,   464,   257,   256,
     596,   598,   600,   599,     0,   453,   461,     0,     0,   390,
     389,     0,   400,     0,   394,     0,   118,     0,   364,     0,
     301,   250,   265,   297,   296,   262,   292,   292,     0,   292,
       0,     0,   280,     0,   254,   253,     0,     0,     0,     0,
       0,     0,   136,     0,     0,     0,     0,   405,   292,   416,
       0,     0,     0,     0,     0,   107,     0,     0,   110,    97,
      98,    99,   100,   101,   102,   103,   104,   105,   106,     0,
       0,     0,   373,   381,   167,     0,     0,   172,   157,   158,
     159,   160,   161,   162,   163,   164,   165,   166,   344,   361,
     199,   191,   190,   194,   195,     0,     0,   221,   428,     0,
     595,   388,     0,   385,   402,     0,   396,     0,     0,     0,
     529,   266,     0,     0,     0,   292,     0,   292,   292,   278,
       0,   113,     0,   144,   511,   123,     0,     0,     0,     0,
     412,     0,     0,   179,     0,   187,     0,     0,     0,     0,
     370,   376,     0,     0,     0,     0,   196,   597,     0,     0,
     403,   397,     0,   365,   292,   292,   292,     0,   286,     0,
       0,     0,   170,     0,   139,   135,   133,   292,   413,     0,
       0,     0,   182,     0,     0,   178,   108,   109,   374,   168,
     169,   193,     0,   386,   292,   271,   267,   270,   292,   284,
     279,   122,     0,     0,   181,   180,   186,     0,   184,     0,
       0,     0,   363,   292,     0,     0,   137,   414,   183,     0,
     189,     0,   387,     0,   272,     0,   285,   185,   177,   188,
     362,   292,   292,   273,   268,   292,   274
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1068, -1068,  1546, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,  -322, -1068,
   -1068, -1068, -1068,  -112,  -226, -1068, -1068,  1278, -1068,   541,
   -1068, -1068, -1068, -1068, -1068, -1068,  -930,  -605,  -115,  -599,
   -1068, -1068, -1068,  1463,  -237, -1068, -1068, -1068, -1068,   646,
   -1068, -1068,   880, -1068, -1068,  1034, -1068, -1068,   889, -1068,
   -1068,  -118,   -24,   914,  1071, -1068, -1068,  1325, -1068, -1068,
   -1067, -1068, -1068,  1300, -1068, -1068,  1306,  -974,  -566, -1068,
   -1068,  1025, -1068,  1485,   909, -1068,   522, -1068, -1068, -1068,
   -1068,  1270, -1068, -1068, -1068, -1068, -1068, -1068, -1068,  1418,
    -737, -1068, -1068, -1068, -1068, -1068,  1006, -1068,   580,  -823,
   -1068, -1068, -1068, -1068, -1068,   916, -1068,   -44,  1097, -1068,
   -1068,  1098, -1068, -1068,   882, -1068,  -507, -1068,   -93, -1068,
    1526, -1068, -1068, -1068, -1068, -1068, -1068, -1068,  -103, -1068,
   -1068,  -121,  -586, -1068, -1068, -1068, -1068,  -105,   -92,   -88,
     -87,   -86, -1068, -1068,  -101,   -50, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068,   -83, -1068, -1068, -1068, -1068, -1068,
     -78,   -57,   -45,   -55,   -54,   -53, -1068, -1068, -1068, -1068,
    -111,   -96,   -98,   -94,   -52,   -51,   -48, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068, -1068,
   -1068, -1068, -1068, -1068, -1068,   894, -1068,  -521
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    44,    45,    46,    47,    48,    49,    50,    51,    52,
     152,   154,   156,   131,    53,    54,    55,    56,   363,   895,
      57,   324,    58,   228,   229,    59,   320,   321,   869,   870,
      60,   327,  1055,  1054,  1136,   873,   629,   630,   631,   632,
     438,    61,    62,   342,   343,  1146,  1221,    63,   720,   721,
      64,   462,   463,    65,   214,   215,    66,   458,   459,    67,
     465,   469,   110,   856,   772,    68,   306,   307,   308,   844,
    1121,    69,   317,   318,    70,   312,   313,   845,  1122,    71,
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
       128,   129,   584,   322,   263,   216,   262,   137,   270,   386,
     338,   294,   146,   149,   150,   295,   261,   264,   157,   437,
     460,   265,   266,   267,   861,   339,   280,   466,   648,   649,
     862,   286,   441,   441,   442,   442,   461,   451,   451,   773,
     660,   452,   452,   205,   846,   671,   848,   871,   206,   202,
     464,   851,   287,   791,   289,   290,   291,   296,   297,   271,
    1018,   298,  1059,  1123,   288,   813,   418,   888,   836,   860,
    1010,   340,   779,   717,   303,   209,   209,   838,    90,   835,
     419,  1103,   718,   171,   961,   819,   820,   821,   863,   420,
    1104,   583,   418,   962,   828,  1174,   584,    92,  1196,   572,
     219,   566,   421,   569,   602,   643,   419,   876,   840,  1063,
     422,   102,   132,  1182,   221,   420,    96,   107,   837,   107,
     423,   101,   222,  1137,  1138,  1139,   839,   340,   421,  1064,
     133,   877,    99,   577,   633,   304,   422,   642,   121,   578,
     117,   100,   227,   638,   766,   123,   423,   104,   879,   118,
     713,    94,    95,   767,   210,   210,   573,   567,   107,   379,
     107,   603,   340,   107,   843,  1234,   841,   107,   918,   107,
     107,   889,   716,   305,   106,   925,   424,   425,   107,   963,
     879,   341,   426,   427,   428,   429,   430,   431,   432,   433,
     434,   634,   107,   707,   708,   709,   710,   435,   711,  1183,
     639,   842,   424,   425,   107,   835,   880,   714,   426,   427,
     428,   429,   430,   431,   432,   433,   434,  1012,   780,   643,
      91,   972,   103,   435,   719,  1211,   890,   108,   109,   126,
     127,   436,  1105,   628,   997,   964,   622,   341,   880,    93,
    1225,   220,   836,  1202,   837,   994,   995,   864,   879,   998,
     999,   838,   839,  1002,  1184,  1004,  1005,   436,   105,   628,
     812,   418,   881,   300,  1175,   300,   882,   883,   144,   145,
     147,   148,   341,   805,   111,   419,   798,   809,   112,   830,
     929,   107,   840,  1040,   420,   799,  1043,   319,   954,   300,
     107,  1058,   841,   413,   881,   107,   113,   421,   882,   883,
     300,   300,   956,   300,   713,   422,   880,  1033,   300,   694,
     695,   231,   114,   580,   970,   423,  1034,  1059,   301,   581,
     301,   706,   891,   892,   893,   894,   200,   842,   896,   897,
     898,   899,   900,   901,   902,   903,   904,   905,   843,   907,
     908,   909,   910,   911,   301,   912,   232,   233,   414,   476,
     201,   916,   480,   234,   484,   301,   301,   769,   301,  1177,
     235,   715,   881,   301,   300,   300,   882,   883,   932,   933,
     934,   424,   425,   302,  1158,   410,   300,   426,   427,   428,
     429,   430,   431,   432,   433,   434,   300,   722,   252,   836,
     300,   974,   435,   300,   309,   300,   254,   958,   838,   411,
    1114,   300,   853,   415,   477,  1116,   724,   481,   314,   485,
     455,   467,   256,   473,   122,   300,   624,  1056,   474,   301,
     301,   300,   303,   309,   257,   300,   436,   314,   628,   840,
     258,   301,   608,   614,   726,   124,   668,   619,   729,   125,
     959,   301,   177,   178,   723,   301,   960,   672,   301,  1051,
     301,  1131,   224,   770,   771,   310,   301,   107,   130,   478,
     225,   135,   482,   725,   486,   871,  1119,   227,   136,   315,
     301,  1124,  1057,  1126,   487,  1022,   301,   328,  1021,   216,
     301,    97,    98,   304,   310,   843,  1024,   364,   315,   854,
     855,   727,  1141,   311,   138,   730,  1030,   263,   673,   262,
    1083,   270,   227,  1098,   294,  1099,  1052,   316,   295,   261,
     264,  1107,   205,   139,   265,   266,   267,   206,   202,   280,
     151,   305,   311,  1120,   286,  1108,   316,   338,   115,   116,
     162,  1160,   861,   861,   861,  1161,   329,   153,   862,   862,
     862,   155,   339,  1134,   198,   287,   330,   289,   290,   291,
     296,   297,   271,   806,   298,   213,   810,   288,   208,  1167,
     217,  1169,  1170,   119,   120,   223,  1061,   661,   662,   663,
     664,   218,   665,   707,   708,   709,   710,   226,   711,   831,
     227,   936,   937,   938,   939,   940,   941,   942,   943,   944,
     945,  1080,   947,   948,   949,   950,   951,   861,  1195,   370,
    1197,   460,   418,   862,   230,  1176,   663,   664,   969,   665,
     441,  1203,   442,   140,   141,   451,   419,   461,   332,   452,
     319,   865,   323,   906,   325,   420,   142,   143,  1212,   946,
     333,   464,  1215,   866,   709,   710,   326,   711,   421,   867,
     158,   159,  1106,   334,   364,   367,   422,  1224,   661,   662,
     663,   664,   412,   665,   416,   417,   423,   457,   475,   868,
     479,   483,   930,   544,   545,  1233,   546,   547,   548,  1236,
     549,   661,   662,   663,   664,   550,   665,   551,   661,   662,
     663,   664,   552,   665,   661,   662,   663,   664,   557,   665,
     955,   957,   553,   418,   554,   555,  1216,   707,   708,   709,
     710,   556,   711,   971,  1065,   558,   975,   419,   559,   560,
     561,   563,   424,   425,   562,   564,   420,   418,   426,   427,
     428,   429,   430,   431,   432,   433,   434,  1069,   565,   421,
     571,   419,   576,   435,  1070,   574,   575,   422,   666,   579,
     420,   582,  1147,   586,  1148,  1149,  1047,   423,  1049,   661,
     662,   663,   664,   421,   665,   661,   662,   663,   664,  1152,
     665,   422,   585,  1153,  1154,   587,   588,   436,   589,   628,
     590,   423,   591,   592,   661,   662,   663,   664,   593,   665,
     594,   595,   596,  1155,   597,   605,   661,   662,   663,   664,
    1159,   665,   598,   661,   662,   663,   664,   599,   665,   600,
     601,   604,   606,   424,   425,  1071,   607,   610,  1171,   426,
     427,   428,   429,   430,   431,   432,   433,   434,   611,   612,
    1179,   584,  1180,   613,   435,   616,   618,   424,   425,   344,
    1072,   617,   621,   426,   427,   428,   429,   430,   431,   432,
     433,   434,  1073,   345,   625,   626,   627,   635,   435,  1074,
     636,   637,   346,   344,   661,   662,   663,   664,   436,   665,
     628,   661,   662,   663,   664,   347,   665,   345,   640,   641,
    1210,   644,   645,   348,   646,   674,   346,   647,   650,   651,
     652,   653,   436,   349,   628,   654,   655,  1220,   656,   347,
     661,   662,   663,   664,   657,   665,   675,   348,   658,  1229,
     659,   661,   662,   663,   664,   667,   665,   349,   676,   677,
    1075,   661,   662,   663,   664,   678,   665,  1076,   679,   680,
     681,   682,   683,   661,   662,   663,   664,   684,   665,   661,
     662,   663,   664,  1132,   665,   685,   686,   688,   670,   350,
     351,   687,   344,   690,   914,   352,   353,   354,   355,   356,
     357,   358,   359,   360,   689,   728,   345,  1077,   691,   731,
     361,   732,   915,   350,   351,   346,   344,  1078,   692,   352,
     353,   354,   355,   356,   357,   358,   359,   360,   347,  1079,
     345,   693,   696,   697,   361,  1084,   348,   698,   699,   346,
     700,   701,   418,   702,   362,   703,   349,   704,   705,  1115,
     733,  1117,   347,   734,   735,   736,   419,   737,   738,   739,
     348,   707,   708,   709,   710,   420,   711,   740,   362,   741,
     349,   742,   707,   708,   709,   710,   743,   711,   421,  1066,
     661,   662,   663,   664,   744,   665,   422,   707,   708,   709,
     710,   745,   711,   746,   747,   748,   423,   749,   750,   751,
     752,   753,   350,   351,   754,   755,   756,   757,   352,   353,
     354,   355,   356,   357,   358,   359,   360,  1088,  1067,   661,
     662,   663,   664,   361,   665,   758,   350,   351,  1089,   759,
     765,   760,   352,   353,   354,   355,   356,   357,   358,   359,
     360,   761,   762,  1090,   763,   764,   774,   361,   795,   800,
     776,   777,   424,   425,   778,   783,  1228,   362,   426,   427,
     428,   429,   430,   431,   432,   433,   434,   160,   784,   785,
     786,   789,   790,   435,     1,     2,   792,   707,   708,   709,
     710,   362,   711,   794,     3,     4,     5,   707,   708,   709,
     710,     6,   711,   796,   797,     7,   801,     8,     9,   803,
      10,   804,    11,    12,    13,    14,   807,   436,   707,   708,
     709,   710,   808,   711,   780,    15,   811,   814,    16,   826,
     707,   708,   709,   710,   815,   711,   816,   707,   708,   709,
     710,    17,   711,  1091,   817,   707,   708,   709,   710,   818,
     711,   822,   823,  1092,    18,    19,    20,   824,   825,   827,
      21,   707,   708,   709,   710,   833,   711,   661,   662,   663,
     664,    22,   665,    23,  1093,    24,    25,    26,    27,    28,
     834,   832,   847,   849,    29,    30,  1094,   850,   852,    31,
      32,    33,    34,  1095,   857,   858,   859,   887,   913,    35,
      36,  1096,    37,     1,     2,   665,    38,   874,   875,    39,
      40,    41,    42,     3,     4,     5,   917,  1097,   919,   920,
       6,   921,   922,  1100,     7,   923,     8,     9,   924,    10,
     926,    11,    12,    13,    14,   927,   928,   711,    43,  1085,
     661,   662,   663,   664,    15,   665,   931,    16,   952,   661,
     662,   663,   664,   966,   665,   968,   661,   662,   663,   664,
      17,   665,   976,   977,   661,   662,   663,   664,   979,   665,
     978,   980,   981,    18,    19,    20,   982,   983,   984,    21,
     661,   662,   663,   664,   985,   665,   661,   662,   663,   664,
      22,   665,    23,   986,    24,    25,    26,    27,    28,   987,
     988,   989,   990,    29,    30,  1156,   991,   992,    31,    32,
      33,    34,  1186,   993,   996,  1000,  1001,   231,    35,    36,
    1187,    37,  1003,  1006,  1008,    38,  1009,  1010,    39,    40,
      41,    42,   200,   170,  1014,  1015,  1189,   171,  1016,  1017,
     567,  1031,  1190,   661,   662,   663,   664,  1019,   665,  1020,
    1023,  1025,   232,   233,   172,  1032,   201,    43,  1027,   234,
     661,   662,   663,   664,  1035,   665,   235,   236,   237,  1036,
    1042,   238,   239,  1037,   240,   241,  1038,  1039,   242,   243,
     244,   245,   246,   247,   248,  1041,   249,   250,   251,  1086,
     661,   662,   663,   664,   252,   665,  1044,   173,   174,  1191,
     253,  1045,   254,  1046,  1048,  1050,  1053,   255,   661,   662,
     663,   664,  1060,   665,  1026,   175,   176,  1062,   256,  1109,
      -1,   163,   164,   165,   166,   167,   168,   169,   199,  1068,
     257,   213,   200,   170,  1081,  1087,   258,   171,  1102,  1110,
    1112,   163,   164,   165,   166,   167,   168,   169,   177,   178,
    1118,  1111,  1125,   170,   172,  1127,   201,   171,   661,   662,
     663,   664,  1082,   665,  1129,   707,   708,   709,   710,  1128,
     711,  1130,  1142,  1144,   172,   661,   662,   663,   664,  1145,
     665,  1150,  1157,  1164,  1162,  1165,   369,   661,   662,   663,
     664,  1166,   665,  1168,  1173,  1181,  1192,   173,   174,   707,
     708,   709,   710,  1178,   711,  1206,  1194,   370,  1199,   371,
     372,  1198,  1113,  1200,  1213,   175,   176,   173,   174,  1133,
    1207,  1208,   661,   662,   663,   664,  1209,   665,  1214,  1143,
     234,  1218,   373,   374,  1217,   175,   176,   235,  1219,   369,
    1223,  1151,  1230,  1226,   328,  1227,  1231,  1232,   177,   178,
     161,  1235,  1135,  1172,   661,   662,   663,   664,   623,   665,
     370,   456,   371,   372,   967,  1101,   802,   935,   177,   178,
     375,   965,   376,   254,   377,   333,  1185,   620,   615,   775,
     378,   454,   953,   234,   379,   373,   374,   829,   334,   570,
     235,   609,   380,   381,   382,   669,  1140,   328,   383,   384,
     385,  1163,   213,     0,   878,   973,   331,  1007,  1188,   468,
     661,   662,   663,   664,   788,   665,   661,   662,   663,   664,
     793,   665,     0,   375,     0,   376,   254,   377,   333,   661,
     662,   663,   664,   378,   665,  1013,     0,   379,     0,     0,
       0,   334,     0,     0,     0,   380,   381,   382,     0,     0,
       0,   383,   384,   385,     0,   213,   661,   662,   663,   664,
       0,   665,     0,     0,  1193,     0,   661,   662,   663,   664,
    1201,   665,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1204,   488,   489,   490,   491,   492,   493,
     494,   495,   496,   497,   498,   499,   500,   501,   502,   503,
     504,   505,   506,   507,   508,     0,     0,     0,   509,   510,
    1205,   511,   512,   513,   514,     0,     0,     0,     0,     0,
    1222
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   228,   115,   109,    98,   109,    31,   109,   127,
     121,   109,    36,    37,    38,   109,   109,   109,    42,   134,
     141,   109,   109,   109,   629,   121,   109,   145,   350,   351,
     629,   109,   135,   136,   135,   136,   141,   135,   136,   546,
     362,   135,   136,    93,   610,   367,   612,   633,    93,    93,
     143,   617,   109,   560,   109,   109,   109,   109,   109,   109,
     797,   109,   885,  1037,   109,   586,    29,    34,    43,    32,
      52,    22,    52,    43,    22,     4,     4,    52,    52,     6,
      43,    43,    52,    25,    43,   592,   593,   594,    34,    52,
      52,    32,    29,    52,   601,    32,   322,    52,  1165,   138,
      52,   138,    65,    32,   138,   342,    43,    57,    83,   101,
      73,    83,    34,    52,   193,    52,   193,    83,    45,    83,
      83,    52,   201,  1053,  1054,  1055,    53,    22,    65,   121,
      52,    81,    34,   195,   138,    83,    73,    32,   194,   201,
      34,    43,    83,   138,    43,    34,    83,    83,    83,    43,
     138,   193,   194,    52,    83,    83,   195,   194,    83,   101,
      83,   195,    22,    83,   139,  1232,    93,    83,   675,    83,
      83,   138,    32,   121,   193,   682,   139,   140,    83,   138,
      83,   132,   145,   146,   147,   148,   149,   150,   151,   152,
     153,   195,    83,   139,   140,   141,   142,   160,   144,   138,
     195,   128,   139,   140,    83,     6,   141,   195,   145,   146,
     147,   148,   149,   150,   151,   152,   153,   199,   198,   456,
     194,   728,   194,   160,   194,   200,   193,   193,   194,   193,
     194,   194,   194,   196,   755,   194,    32,   132,   141,   194,
    1214,   193,    43,  1173,    45,   752,   753,   193,    83,   756,
     757,    52,    53,   760,   193,   762,   763,   194,   194,   196,
     582,    29,   197,    83,    32,    83,   201,   202,   193,   194,
     193,   194,   132,   193,   193,    43,    43,   193,   193,   193,
     193,    83,    83,   849,    52,    52,   852,    83,   193,    83,
      83,   194,    93,    83,   197,    83,   193,    65,   201,   202,
      83,    83,   193,    83,   138,    73,   141,    43,    83,   424,
     425,     5,   193,   195,   193,    83,    52,  1140,   138,   201,
     138,   436,   644,   645,   646,   647,    20,   128,   650,   651,
     652,   653,   654,   655,   656,   657,   658,   659,   139,   661,
     662,   663,   664,   665,   138,   667,    40,    41,   138,    83,
      44,   673,    83,    47,    83,   138,   138,    43,   138,   194,
      54,   195,   197,   138,    83,    83,   201,   202,   690,   691,
     692,   139,   140,   193,  1111,   193,    83,   145,   146,   147,
     148,   149,   150,   151,   152,   153,    83,   138,    82,    43,
      83,   193,   160,    83,    22,    83,    90,   719,    52,   193,
     193,    83,    52,   193,   138,   193,   138,   138,    22,   138,
     193,   193,   106,   193,   193,    83,    32,   138,   193,   138,
     138,    83,    22,    22,   118,    83,   194,    22,   196,    83,
     124,   138,    32,    32,   138,   193,    32,    32,   138,   193,
     193,   138,   136,   137,   195,   138,   199,    32,   138,   138,
     138,    32,   193,   139,   140,    83,   138,    83,    83,   193,
     201,   194,   193,   195,   193,  1051,   138,    83,   194,    83,
     138,  1037,   193,  1039,   193,   193,   138,    61,   800,   572,
     138,   193,   194,    83,    83,   139,   193,    83,    83,   139,
     140,   195,  1058,   121,   193,   195,   193,   602,    83,   602,
     193,   602,    83,   193,   602,   193,   195,   121,   602,   602,
     602,   193,   562,    34,   602,   602,   602,   562,   562,   602,
      83,   121,   121,   195,   602,   193,   121,   638,   193,   194,
     193,   193,  1137,  1138,  1139,   193,   120,    83,  1137,  1138,
    1139,    83,   638,  1050,   193,   602,   130,   602,   602,   602,
     602,   602,   602,   577,   602,   119,   580,   602,    23,  1125,
      43,  1127,  1128,   193,   194,    83,   888,   139,   140,   141,
     142,   193,   144,   139,   140,   141,   142,    83,   144,   603,
      83,   696,   697,   698,   699,   700,   701,   702,   703,   704,
     705,   913,   707,   708,   709,   710,   711,  1202,  1164,    24,
    1166,   722,    29,  1202,    34,    32,   141,   142,   726,   144,
     713,  1177,   713,   193,   194,   713,    43,   722,    80,   713,
      83,    46,    39,   195,    43,    52,   193,   194,  1194,   195,
      92,   724,  1198,    58,   141,   142,   193,   144,    65,    64,
     193,   194,   964,   105,    83,    83,    73,  1213,   139,   140,
     141,   142,   127,   144,    52,   193,    83,   198,   127,    84,
     127,   127,   686,    34,    34,  1231,    34,    34,    34,  1235,
      34,   139,   140,   141,   142,    34,   144,    34,   139,   140,
     141,   142,    34,   144,   139,   140,   141,   142,   138,   144,
     714,   715,    34,    29,    34,    34,    32,   139,   140,   141,
     142,    34,   144,   727,   195,   195,   730,    43,    34,    34,
      34,   195,   139,   140,   138,   195,    52,    29,   145,   146,
     147,   148,   149,   150,   151,   152,   153,   195,    83,    65,
      34,    43,    83,   160,   195,   193,   193,    73,   193,    83,
      52,    34,  1064,    34,  1066,  1067,   858,    83,   863,   139,
     140,   141,   142,    65,   144,   139,   140,   141,   142,  1081,
     144,    73,    83,  1085,  1086,    34,    34,   194,    34,   196,
      34,    83,    34,    34,   139,   140,   141,   142,    34,   144,
      34,    34,    34,  1105,    34,    83,   139,   140,   141,   142,
    1112,   144,    34,   139,   140,   141,   142,    34,   144,    34,
      34,    34,    83,   139,   140,   195,    83,   138,  1130,   145,
     146,   147,   148,   149,   150,   151,   152,   153,    83,   138,
    1142,  1047,  1144,    83,   160,    83,    83,   139,   140,    29,
     195,   138,   194,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   195,    43,    34,   195,   193,   195,   160,   195,
     195,    34,    52,    29,   139,   140,   141,   142,   194,   144,
     196,   139,   140,   141,   142,    65,   144,    43,    83,    83,
    1192,   194,   194,    73,   194,    34,    52,   194,   194,   194,
     194,   194,   194,    83,   196,   194,   194,  1209,   194,    65,
     139,   140,   141,   142,   194,   144,    34,    73,   194,  1221,
     194,   139,   140,   141,   142,   194,   144,    83,    34,    34,
     195,   139,   140,   141,   142,    34,   144,   195,    34,    34,
      34,    34,    34,   139,   140,   141,   142,    34,   144,   139,
     140,   141,   142,  1048,   144,    34,   195,    83,   138,   139,
     140,   127,    29,   194,   193,   145,   146,   147,   148,   149,
     150,   151,   152,   153,   193,    34,    43,   195,   194,   127,
     160,    83,   138,   139,   140,    52,    29,   195,   194,   145,
     146,   147,   148,   149,   150,   151,   152,   153,    65,   195,
      43,   194,   194,   194,   160,   195,    73,   194,   194,    52,
     194,   194,    29,   194,   194,   194,    83,   194,   194,  1023,
     127,  1025,    65,    83,   127,    83,    43,    34,    34,    34,
      73,   139,   140,   141,   142,    52,   144,    34,   194,    34,
      83,    34,   139,   140,   141,   142,    34,   144,    65,   138,
     139,   140,   141,   142,    34,   144,    73,   139,   140,   141,
     142,    34,   144,    34,    34,    34,    83,    34,    34,    34,
      34,    34,   139,   140,    34,    34,    34,    34,   145,   146,
     147,   148,   149,   150,   151,   152,   153,   195,   138,   139,
     140,   141,   142,   160,   144,    34,   139,   140,   195,    34,
     138,    34,   145,   146,   147,   148,   149,   150,   151,   152,
     153,    34,    34,   195,    34,   195,    52,   160,   193,    34,
      52,    83,   139,   140,    52,    52,   193,   194,   145,   146,
     147,   148,   149,   150,   151,   152,   153,     0,    52,    83,
      83,    52,    52,   160,     7,     8,    52,   139,   140,   141,
     142,   194,   144,    52,    17,    18,    19,   139,   140,   141,
     142,    24,   144,   194,    83,    28,    52,    30,    31,   193,
      33,   193,    35,    36,    37,    38,    83,   194,   139,   140,
     141,   142,   193,   144,   198,    48,    83,    52,    51,    83,
     139,   140,   141,   142,    52,   144,    52,   139,   140,   141,
     142,    64,   144,   195,    52,   139,   140,   141,   142,    52,
     144,    52,    52,   195,    77,    78,    79,    52,    52,   194,
      83,   139,   140,   141,   142,    34,   144,   139,   140,   141,
     142,    94,   144,    96,   195,    98,    99,   100,   101,   102,
     138,    83,   138,   138,   107,   108,   195,   138,   138,   112,
     113,   114,   115,   195,    83,   193,    83,   138,    83,   122,
     123,   195,   125,     7,     8,   144,   129,   193,   193,   132,
     133,   134,   135,    17,    18,    19,    52,   195,    52,    52,
      24,    52,    52,   195,    28,    52,    30,    31,    52,    33,
      52,    35,    36,    37,    38,    52,    52,   144,   161,   138,
     139,   140,   141,   142,    48,   144,   127,    51,   193,   139,
     140,   141,   142,   193,   144,   193,   139,   140,   141,   142,
      64,   144,   127,   127,   139,   140,   141,   142,    52,   144,
     127,    52,    52,    77,    78,    79,    52,    52,    52,    83,
     139,   140,   141,   142,    52,   144,   139,   140,   141,   142,
      94,   144,    96,    52,    98,    99,   100,   101,   102,    52,
      52,    52,    52,   107,   108,   195,    52,    52,   112,   113,
     114,   115,   195,    52,    52,    52,    52,     5,   122,   123,
     195,   125,    52,   193,    43,   129,    43,    52,   132,   133,
     134,   135,    20,    21,   197,   193,   195,    25,   193,    52,
     194,    83,   195,   139,   140,   141,   142,   195,   144,   195,
     195,   195,    40,    41,    42,    83,    44,   161,   200,    47,
     139,   140,   141,   142,   193,   144,    54,    55,    56,   138,
     193,    59,    60,   138,    62,    63,    83,   138,    66,    67,
      68,    69,    70,    71,    72,    83,    74,    75,    76,   138,
     139,   140,   141,   142,    82,   144,    52,    85,    86,   195,
      88,    52,    90,   195,    34,    34,   193,    95,   139,   140,
     141,   142,    83,   144,   193,   103,   104,    83,   106,   197,
     144,     9,    10,    11,    12,    13,    14,    15,    16,   195,
     118,   119,    20,    21,    83,   195,   124,    25,   193,    83,
      34,     9,    10,    11,    12,    13,    14,    15,   136,   137,
      83,   195,   138,    21,    42,   138,    44,    25,   139,   140,
     141,   142,   193,   144,   193,   139,   140,   141,   142,   138,
     144,    34,    34,    34,    42,   139,   140,   141,   142,    52,
     144,   193,    52,   138,   200,   138,     3,   139,   140,   141,
     142,   138,   144,   193,   193,   197,    34,    85,    86,   139,
     140,   141,   142,   195,   144,    52,   138,    24,   193,    26,
      27,   138,   193,   193,   138,   103,   104,    85,    86,   193,
     197,    52,   139,   140,   141,   142,   131,   144,   138,   193,
      47,    52,    49,    50,   195,   103,   104,    54,   197,     3,
      83,   193,   200,   193,    61,    52,   138,   138,   136,   137,
      44,   138,  1051,   193,   139,   140,   141,   142,   320,   144,
      24,   138,    26,    27,   724,   959,   572,   693,   136,   137,
      87,   722,    89,    90,    91,    92,   193,   317,   312,   548,
      97,   136,   713,    47,   101,    49,    50,   602,   105,   211,
      54,   306,   109,   110,   111,   365,  1056,    61,   115,   116,
     117,  1119,   119,    -1,   638,   729,   120,   765,   193,   126,
     139,   140,   141,   142,   557,   144,   139,   140,   141,   142,
     562,   144,    -1,    87,    -1,    89,    90,    91,    92,   139,
     140,   141,   142,    97,   144,   781,    -1,   101,    -1,    -1,
      -1,   105,    -1,    -1,    -1,   109,   110,   111,    -1,    -1,
      -1,   115,   116,   117,    -1,   119,   139,   140,   141,   142,
      -1,   144,    -1,    -1,   193,    -1,   139,   140,   141,   142,
     193,   144,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   193,   162,   163,   164,   165,   166,   167,
     168,   169,   170,   171,   172,   173,   174,   175,   176,   177,
     178,   179,   180,   181,   182,    -1,    -1,    -1,   186,   187,
     193,   189,   190,   191,   192,    -1,    -1,    -1,    -1,    -1,
     193
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
     265,   127,   221,   221,   221,   266,   241,   241,   241,   241,
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
      83,   221,    83,   101,   121,   195,   138,   138,   195,   195,
     195,   195,   195,   195,   195,   195,   195,   195,   195,   195,
     221,    83,   193,   193,   195,   138,   138,   195,   195,   195,
     195,   195,   195,   195,   195,   195,   195,   195,   193,   193,
     195,   252,   193,    43,    52,   194,   221,   193,   193,   197,
      83,   195,    34,   193,   193,   265,   193,   265,    83,   138,
     195,   273,   281,   280,   281,   138,   281,   138,   138,   193,
      34,    32,   241,   193,   329,   232,   237,   239,   239,   239,
     311,   281,    34,   193,    34,    52,   248,   221,   221,   221,
     193,   193,   221,   221,   221,   221,   195,    52,   303,   221,
     193,   193,   200,   289,   138,   138,   138,   281,   193,   281,
     281,   221,   193,   193,    32,    32,    32,   194,   195,   221,
     221,   197,    52,   138,   193,   193,   195,   195,   193,   195,
     195,   195,    34,   193,   138,   281,   273,   281,   138,   193,
     193,   193,   239,   281,   193,   193,    52,   197,    52,   131,
     221,   200,   281,   138,   138,   281,    32,   195,    52,   197,
     221,   249,   193,    83,   281,   280,   193,    52,   193,   221,
     200,   138,   138,   281,   273,   138,   281
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
     221,   222,   223,   223,   224,   225,   226,   226,   227,   228,
     229,   229,   230,   231,   231,   232,   232,   232,   232,   232,
     234,   233,   235,   233,   236,   233,   237,   233,   238,   233,
     239,   239,   239,   239,   240,   240,   241,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   241,
     241,   241,   241,   241,   241,   241,   241,   241,   241,   241,
     242,   243,   243,   244,   245,   246,   246,   247,   247,   247,
     247,   247,   248,   248,   248,   248,   248,   248,   249,   249,
     250,   251,   251,   252,   252,   252,   252,   252,   252,   252,
     252,   252,   253,   253,   254,   254,   255,   256,   256,   257,
     257,   258,   259,   259,   260,   260,   261,   261,   262,   262,
     262,   262,   263,   263,   264,   264,   264,   264,   264,   264,
     264,   264,   264,   264,   264,   264,   264,   264,   264,   264,
     264,   264,   264,   264,   264,   264,   264,   265,   265,   265,
     265,   265,   265,   266,   266,   266,   267,   267,   267,   268,
     269,   269,   270,   271,   271,   271,   272,   272,   272,   272,
     272,   273,   273,   273,   273,   274,   275,   275,   276,   276,
     276,   277,   278,   278,   279,   279,   279,   280,   280,   280,
     280,   280,   281,   281,   281,   281,   281,   281,   282,   282,
     282,   282,   283,   283,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   284,   284,   284,   284,   284,   284,   284,
     284,   284,   284,   285,   285,   286,   286,   287,   287,   287,
     287,   287,   287,   287,   287,   287,   287,   287,   287,   287,
     288,   288,   289,   289,   290,   290,   291,   292,   293,   293,
     294,   295,   296,   297,   297,   297,   297,   298,   299,   299,
     299,   299,   300,   301,   301,   302,   302,   302,   303,   303,
     303,   304,   304,   305,   305,   305,   305,   305,   305,   306,
     306,   306,   306,   306,   306,   307,   308,   308,   309,   309,
     309,   310,   310,   310,   310,   311,   311,   312,   312,   312,
     312,   312,   314,   315,   313,   316,   316,   316,   316,   317,
     317,   318,   318,   319,   319,   319,   319,   319,   319,   319,
     320,   320,   320,   320,   320,   320,   320,   320,   320,   320,
     321,   321,   322,   322,   323,   323,   323,   323,   324,   324,
     325,   325,   326,   326,   327,   327,   328,   328,   328,   328,
     328,   328,   328,   328,   328,   328,   328,   328,   328,   328,
     328,   328,   328,   328,   328,   328,   328,   328,   328,   328,
     328,   328,   328,   329,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   351,   352,   353,
     354,   355,   356,   357,   357,   358,   359,   360,   361,   362,
     363,   364,   365,   366,   367,   368,   369,   370,   371,   372,
     373,   374,   375,   376,   377,   378,   379,   380,   381,   382,
     382,   383,   384,   385,   386,   387,   388,   389,   390,   391,
     392,   393,   394,   395,   396,   397,   398,   399,   400,   401,
     402,   403,   404,   405,   406,   407,   408,   409,   410,   411,
     412,   413,   414,   415,   416,   417,   418,   419,   420,   421,
     422,   423,   424,   425,   426,   427,   428,   428,   429,   429,
     430
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
       4,     1,     4,     7,     3,     4,     2,     1,     4,     4,
       2,     1,     7,     3,     1,     1,     1,     1,     1,     1,
       0,     5,     0,     8,     0,     8,     0,    10,     0,     8,
       2,     2,     1,     1,     4,     2,     3,     1,     1,     1,
       3,     3,     3,     3,     3,     2,     2,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     6,     6,
       5,     1,     4,     4,     4,     2,     1,     9,     6,     5,
       7,     7,     2,     4,     3,     5,     3,     1,     2,     1,
       6,     3,     1,     5,     3,     3,     4,     2,     2,     3,
       1,     1,     2,     5,     3,     1,     1,     2,     5,     3,
       1,     1,     2,     5,     3,     1,     1,     1,     2,     5,
       3,     6,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     2,     4,     3,
       5,     1,     3,     2,     2,     1,     2,     2,     1,     4,
       2,     1,     4,     2,     1,     4,     3,     5,     9,     1,
       5,     3,     5,     7,     9,     4,     2,     1,     5,     7,
       4,     4,     2,     1,     7,     9,     6,     1,     1,     1,
       1,     1,     0,     1,     1,     1,     2,     2,     2,     5,
       3,     6,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     5,     6,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       5,     6,     7,     5,     1,     3,     3,     4,     2,     1,
       5,     3,     4,     4,     6,     3,     5,     3,     2,     5,
       3,     6,     4,     2,     1,     5,     7,     9,     0,     3,
       3,     2,     5,     5,     6,     3,     7,     8,     5,     5,
       6,     3,     7,     8,     5,     6,     3,     1,     1,     1,
       1,     1,     3,     4,     6,     1,     2,     1,     1,     1,
       1,     1,     0,     0,     5,     2,     5,     3,     6,     3,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       3,     1,     3,     6,     1,     1,     1,     1,     3,     1,
       3,     6,     2,     5,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     3,     3,     1,     3,
       3,     3,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     3,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     1,     1,     3,     3,     3,     3,     5,
       3,     3,     3,     1,     3,     3,     3,     1,     1,     1,
       1,     1,     3,     1,     1,     1,     1,     3,     3,     3,
       3,     1,     1,     3,     3,     3,     1,     1,     1,     3,
       3,     3,     3,     3,     3,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     3,     2,     2,
       2
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
      -1,   221,    -1,    51,   193,   226,    32,    -1,    51,   194,
     224,   195,   193,   226,    32,    -1,    39,    34,    83,    -1,
      33,   193,   226,    32,    -1,   226,   227,    -1,   227,    -1,
      83,    34,   221,   193,    -1,    48,   193,   229,    32,    -1,
     229,   230,    -1,   230,    -1,    83,   194,   266,   195,    34,
     221,   193,    -1,   231,   138,   232,    -1,   232,    -1,    58,
      -1,    46,    -1,    84,    -1,   345,    -1,   346,    -1,    -1,
      77,   193,   234,   239,    32,    -1,    -1,    77,   194,   333,
     195,   193,   235,   239,    32,    -1,    -1,    77,   194,   130,
     195,   193,   236,   239,    32,    -1,    -1,    77,   194,   120,
     138,   231,   195,   237,   193,   239,    32,    -1,    -1,    77,
     194,   120,   195,   238,   193,   239,    32,    -1,   239,   240,
      -1,   239,   242,    -1,   240,    -1,   242,    -1,   241,    34,
     241,   193,    -1,   241,   193,    -1,   194,   241,   195,    -1,
     243,    -1,    43,    -1,    52,    -1,   241,   140,   241,    -1,
     241,   139,   241,    -1,   241,   141,   241,    -1,   241,   142,
     241,    -1,   241,   144,   241,    -1,   139,   241,    -1,   140,
     241,    -1,   145,   194,   241,   195,    -1,   146,   194,   241,
     195,    -1,   147,   194,   241,   195,    -1,   148,   194,   241,
     195,    -1,   149,   194,   241,   195,    -1,   150,   194,   241,
     195,    -1,   151,   194,   241,   195,    -1,   152,   194,   241,
     195,    -1,   153,   194,   241,   195,    -1,   160,   194,   241,
     195,    -1,    29,   194,   221,   195,    -1,    65,   194,   221,
     138,   221,   195,    -1,    73,   194,   221,   138,   221,   195,
      -1,   196,    83,    34,   241,   193,    -1,    83,    -1,    83,
     194,   266,   195,    -1,   113,   193,   246,    32,    -1,    79,
     193,   246,    32,    -1,   246,   247,    -1,   247,    -1,   132,
      83,   193,   101,   248,   193,   131,   249,   193,    -1,   132,
      83,   193,   121,   221,   193,    -1,   132,    83,    34,   221,
     193,    -1,   132,    83,   138,    83,    34,   221,   193,    -1,
      22,    83,   138,    83,    34,   221,   193,    -1,   248,    52,
      -1,   248,    52,   197,    52,    -1,   248,   138,    52,    -1,
     248,   138,    52,   197,    52,    -1,    52,   197,    52,    -1,
      52,    -1,   249,   221,    -1,   221,    -1,   114,    34,   198,
     251,   199,   193,    -1,   251,   193,   252,    -1,   252,    -1,
     252,   138,   194,   221,   195,    -1,   252,   138,    43,    -1,
     252,   138,    52,    -1,   252,   194,   221,   195,    -1,   252,
      43,    -1,   252,    52,    -1,   194,   221,   195,    -1,    43,
      -1,    52,    -1,   122,   193,    -1,   122,   194,   254,   195,
     193,    -1,   254,   138,   255,    -1,   255,    -1,   331,    -1,
      19,   193,    -1,    19,   194,   257,   195,   193,    -1,   257,
     138,   258,    -1,   258,    -1,   331,    -1,   115,   193,    -1,
     115,   194,   260,   195,   193,    -1,   260,   138,   261,    -1,
     261,    -1,   344,    -1,   350,    -1,   123,   193,    -1,   123,
     194,   263,   195,   193,    -1,   123,   265,   193,    -1,   123,
     194,   263,   195,   265,   193,    -1,   263,   138,   264,    -1,
     264,    -1,   330,    -1,   331,    -1,   332,    -1,   333,    -1,
     334,    -1,   335,    -1,   336,    -1,   337,    -1,   338,    -1,
     339,    -1,   340,    -1,   357,    -1,   341,    -1,   379,    -1,
     342,    -1,   343,    -1,   344,    -1,   345,    -1,   347,    -1,
     348,    -1,   349,    -1,   383,    -1,   384,    -1,   265,    83,
      -1,   265,    83,    34,    83,    -1,   265,   138,    83,    -1,
     265,   138,    83,    34,    83,    -1,    83,    -1,    83,    34,
      83,    -1,   140,    52,    -1,   139,    52,    -1,    52,    -1,
     140,    43,    -1,   139,    43,    -1,    43,    -1,    36,   193,
     269,    32,    -1,   269,   270,    -1,   270,    -1,   271,   138,
     272,   193,    -1,   121,    83,    -1,    83,    -1,    22,    83,
     138,    83,    -1,   280,   138,   273,    -1,   281,   138,   280,
     138,   273,    -1,   281,   138,   281,   138,   281,   138,   280,
     138,   273,    -1,   281,    -1,   281,   138,   281,   138,   281,
      -1,   281,   138,   281,    -1,   281,   138,   281,   138,   281,
      -1,   281,   138,   281,   138,   281,   138,   281,    -1,   281,
     138,   281,   138,   281,   138,   281,   138,   281,    -1,    38,
     193,   275,    32,    -1,   275,   276,    -1,   276,    -1,   121,
      83,   138,   281,   193,    -1,    22,    83,   138,    83,   138,
     281,   193,    -1,    83,   138,   281,   193,    -1,    37,   193,
     278,    32,    -1,   278,   279,    -1,   279,    -1,   121,    83,
     138,   281,   138,   281,   193,    -1,    22,    83,   138,    83,
     138,   281,   138,   281,   193,    -1,    83,   138,   281,   138,
     281,   193,    -1,     6,    -1,    45,    -1,    93,    -1,    53,
      -1,   128,    -1,    -1,    52,    -1,    43,    -1,    83,    -1,
     139,    52,    -1,   139,    43,    -1,    35,   193,    -1,    35,
     194,   283,   195,   193,    -1,    35,   265,   193,    -1,    35,
     194,   283,   195,   265,   193,    -1,   283,   138,   284,    -1,
     284,    -1,   350,    -1,   351,    -1,   352,    -1,   353,    -1,
     354,    -1,   355,    -1,   356,    -1,   357,    -1,   358,    -1,
     359,    -1,   360,    -1,   361,    -1,   362,    -1,   363,    -1,
     364,    -1,   365,    -1,   366,    -1,   367,    -1,   368,    -1,
     369,    -1,   370,    -1,   371,    -1,   372,    -1,   373,    -1,
     341,    -1,   374,    -1,   375,    -1,   376,    -1,   377,    -1,
     378,    -1,   380,    -1,   381,    -1,   385,    -1,   386,    -1,
     387,    -1,   331,    -1,   388,    -1,   389,    -1,   390,    -1,
     107,   194,   286,   195,   193,    -1,   107,   194,   286,   195,
     265,   193,    -1,   286,   138,   287,    -1,   287,    -1,   357,
      -1,   358,    -1,   367,    -1,   373,    -1,   341,    -1,   374,
      -1,   375,    -1,   376,    -1,   377,    -1,   378,    -1,   385,
      -1,   386,    -1,   387,    -1,   108,   194,   286,   195,   193,
      -1,   108,   194,   286,   195,   265,   193,    -1,   200,    83,
     200,   138,   200,    83,   200,    -1,   200,    83,   200,   138,
     281,    -1,   289,    -1,   290,   138,   289,    -1,   135,   265,
     193,    -1,    94,   193,   293,    32,    -1,   293,   294,    -1,
     294,    -1,    83,   194,   221,   195,   193,    -1,   129,   265,
     193,    -1,    96,   193,   297,    32,    -1,   297,    83,   221,
     193,    -1,   297,    83,   138,    83,   221,   193,    -1,    83,
     221,   193,    -1,    83,   138,    83,   221,   193,    -1,    99,
     265,   193,    -1,    98,   193,    -1,    98,   194,   264,   195,
     193,    -1,    98,   265,   193,    -1,    98,   194,   264,   195,
     265,   193,    -1,    18,   193,   301,    32,    -1,   301,   302,
      -1,   302,    -1,    83,   303,    34,   221,   193,    -1,    83,
     138,    83,   303,    34,   221,   193,    -1,     4,    83,   194,
      52,   195,   303,    34,   221,   193,    -1,    -1,   194,    52,
     195,    -1,   194,    43,   195,    -1,    17,   193,    -1,    17,
     194,    23,   195,   193,    -1,    31,   194,    83,   195,   193,
      -1,    31,   194,    83,   195,   265,   193,    -1,    31,    83,
     193,    -1,    31,   194,    83,   201,    83,   195,   193,    -1,
      31,   194,    83,   201,    83,   195,   265,   193,    -1,    31,
      83,   201,    83,   193,    -1,    30,   194,    83,   195,   193,
      -1,    30,   194,    83,   195,   265,   193,    -1,    30,    83,
     193,    -1,    30,   194,    83,   201,    83,   195,   193,    -1,
      30,   194,    83,   201,    83,   195,   265,   193,    -1,    30,
      83,   201,    83,   193,    -1,    78,   194,   308,   195,   310,
     193,    -1,   308,   138,   309,    -1,   309,    -1,   382,    -1,
     383,    -1,   384,    -1,   311,    -1,   310,   138,   311,    -1,
     311,   194,   281,   195,    -1,   310,   138,   311,   194,   281,
     195,    -1,   312,    -1,   311,   312,    -1,    83,    -1,   202,
      -1,   141,    -1,   197,    -1,   201,    -1,    -1,    -1,   102,
     314,   241,   315,   193,    -1,   125,   193,    -1,   125,   194,
     317,   195,   193,    -1,   125,   265,   193,    -1,   125,   194,
     317,   195,   265,   193,    -1,   317,   138,   318,    -1,   318,
      -1,   264,    -1,   391,    -1,   392,    -1,   393,    -1,   394,
      -1,   395,    -1,   396,    -1,   397,    -1,   398,    -1,   319,
      -1,   350,    -1,   385,    -1,   386,    -1,   352,    -1,   354,
      -1,   351,    -1,   353,    -1,   388,    -1,   389,    -1,   320,
     138,   321,    -1,   320,    -1,     7,    52,   193,    -1,     7,
     194,   321,   195,    52,   193,    -1,   320,    -1,   375,    -1,
     358,    -1,   399,    -1,   323,   138,   324,    -1,   323,    -1,
       8,    52,   193,    -1,     8,   194,   324,   195,    52,   193,
      -1,   161,   193,    -1,   161,   194,   327,   195,   193,    -1,
     328,   138,   327,    -1,   328,    -1,   400,    -1,   401,    -1,
     402,    -1,   403,    -1,   404,    -1,   405,    -1,   406,    -1,
     407,    -1,   408,    -1,   409,    -1,   410,    -1,   411,    -1,
     412,    -1,   413,    -1,   414,    -1,   415,    -1,   416,    -1,
     417,    -1,   419,    -1,   420,    -1,   421,    -1,   422,    -1,
     423,    -1,   424,    -1,   425,    -1,   426,    -1,   418,    -1,
      52,    -1,    43,    -1,    26,    34,    52,    -1,   119,    34,
      52,    -1,   116,    34,    52,    -1,    61,    -1,    97,    34,
      52,    -1,   111,    34,    52,    -1,    27,    34,    52,    -1,
       3,    34,    52,    -1,    87,    -1,    89,    -1,    91,    -1,
      54,    34,    52,    -1,    49,    34,    52,    -1,    50,    34,
      52,    -1,   101,    34,    52,    -1,    24,    34,   329,    -1,
      64,    34,   329,    -1,   115,    -1,   117,    34,    52,    -1,
     109,    34,   329,    -1,    25,    34,    83,    -1,    85,    34,
     430,    -1,    85,    34,    52,    -1,    42,    34,    52,    -1,
     103,    34,    52,    -1,   104,    34,    52,    -1,    59,    34,
      52,    -1,    60,    34,    52,    -1,    90,    -1,    47,    -1,
      20,    34,   329,    -1,    71,    34,    52,    -1,    66,    34,
     329,    -1,    68,    34,   329,    -1,    95,    34,   194,   290,
     195,    -1,    67,    34,   329,    -1,    76,    34,    83,    -1,
      75,    34,    52,    -1,    74,    -1,   106,    34,   329,    -1,
      69,    34,    52,    -1,    70,    34,    52,    -1,    62,    -1,
      63,    -1,    88,    -1,     5,    -1,   124,    -1,    44,    34,
      52,    -1,   118,    -1,    82,    -1,    41,    -1,   110,    -1,
      55,    34,    52,    -1,    56,    34,    52,    -1,    80,    34,
      57,    -1,    80,    34,    81,    -1,   105,    -1,    92,    -1,
     136,    34,    83,    -1,   137,    34,   427,    -1,    40,    34,
     430,    -1,    21,    -1,    86,    -1,    72,    -1,   126,    34,
     329,    -1,    14,    34,   267,    -1,     9,    34,   329,    -1,
      11,    34,   267,    -1,    12,    34,   329,    -1,    13,    34,
      52,    -1,    10,    -1,    15,    34,    52,    -1,    16,    34,
      52,    -1,   162,    34,    52,    -1,   163,    34,    52,    -1,
     164,    34,    52,    -1,   165,    34,    52,    -1,   166,    34,
      52,    -1,   167,    34,    52,    -1,   168,    34,    52,    -1,
     169,    34,    52,    -1,   170,    34,    52,    -1,   171,    34,
      52,    -1,   172,    34,    52,    -1,   173,    34,    52,    -1,
     174,    34,    52,    -1,   175,    34,    52,    -1,   176,    34,
      52,    -1,   177,    34,   329,    -1,   178,    34,   329,    -1,
     179,    34,    52,    -1,   180,    34,   430,    -1,   181,    34,
     329,    -1,   182,    34,   329,    -1,   186,    34,    52,    -1,
     187,    34,    52,    -1,   189,    34,   329,    -1,   190,    34,
      52,    -1,   191,    34,   329,    -1,   192,    34,   329,    -1,
      83,   197,    83,    -1,    52,    -1,    52,   197,    52,    -1,
     198,   428,    -1,   429,   428,    -1,   429,   199,    -1
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
     346,   351,   353,   358,   366,   370,   375,   378,   380,   385,
     390,   393,   395,   403,   407,   409,   411,   413,   415,   417,
     419,   420,   426,   427,   436,   437,   446,   447,   458,   459,
     468,   471,   474,   476,   478,   483,   486,   490,   492,   494,
     496,   500,   504,   508,   512,   516,   519,   522,   527,   532,
     537,   542,   547,   552,   557,   562,   567,   572,   577,   584,
     591,   597,   599,   604,   609,   614,   617,   619,   629,   636,
     642,   650,   658,   661,   666,   670,   676,   680,   682,   685,
     687,   694,   698,   700,   706,   710,   714,   719,   722,   725,
     729,   731,   733,   736,   742,   746,   748,   750,   753,   759,
     763,   765,   767,   770,   776,   780,   782,   784,   786,   789,
     795,   799,   806,   810,   812,   814,   816,   818,   820,   822,
     824,   826,   828,   830,   832,   834,   836,   838,   840,   842,
     844,   846,   848,   850,   852,   854,   856,   858,   861,   866,
     870,   876,   878,   882,   885,   888,   890,   893,   896,   898,
     903,   906,   908,   913,   916,   918,   923,   927,   933,   943,
     945,   951,   955,   961,   969,   979,   984,   987,   989,   995,
    1003,  1008,  1013,  1016,  1018,  1026,  1036,  1043,  1045,  1047,
    1049,  1051,  1053,  1054,  1056,  1058,  1060,  1063,  1066,  1069,
    1075,  1079,  1086,  1090,  1092,  1094,  1096,  1098,  1100,  1102,
    1104,  1106,  1108,  1110,  1112,  1114,  1116,  1118,  1120,  1122,
    1124,  1126,  1128,  1130,  1132,  1134,  1136,  1138,  1140,  1142,
    1144,  1146,  1148,  1150,  1152,  1154,  1156,  1158,  1160,  1162,
    1164,  1166,  1168,  1170,  1176,  1183,  1187,  1189,  1191,  1193,
    1195,  1197,  1199,  1201,  1203,  1205,  1207,  1209,  1211,  1213,
    1215,  1221,  1228,  1236,  1242,  1244,  1248,  1252,  1257,  1260,
    1262,  1268,  1272,  1277,  1282,  1289,  1293,  1299,  1303,  1306,
    1312,  1316,  1323,  1328,  1331,  1333,  1339,  1347,  1357,  1358,
    1362,  1366,  1369,  1375,  1381,  1388,  1392,  1400,  1409,  1415,
    1421,  1428,  1432,  1440,  1449,  1455,  1462,  1466,  1468,  1470,
    1472,  1474,  1476,  1480,  1485,  1492,  1494,  1497,  1499,  1501,
    1503,  1505,  1507,  1508,  1509,  1515,  1518,  1524,  1528,  1535,
    1539,  1541,  1543,  1545,  1547,  1549,  1551,  1553,  1555,  1557,
    1559,  1561,  1563,  1565,  1567,  1569,  1571,  1573,  1575,  1577,
    1579,  1583,  1585,  1589,  1596,  1598,  1600,  1602,  1604,  1608,
    1610,  1614,  1621,  1624,  1630,  1634,  1636,  1638,  1640,  1642,
    1644,  1646,  1648,  1650,  1652,  1654,  1656,  1658,  1660,  1662,
    1664,  1666,  1668,  1670,  1672,  1674,  1676,  1678,  1680,  1682,
    1684,  1686,  1688,  1690,  1692,  1694,  1698,  1702,  1706,  1708,
    1712,  1716,  1720,  1724,  1726,  1728,  1730,  1734,  1738,  1742,
    1746,  1750,  1754,  1756,  1760,  1764,  1768,  1772,  1776,  1780,
    1784,  1788,  1792,  1796,  1798,  1800,  1804,  1808,  1812,  1816,
    1822,  1826,  1830,  1834,  1836,  1840,  1844,  1848,  1850,  1852,
    1854,  1856,  1858,  1862,  1864,  1866,  1868,  1870,  1874,  1878,
    1882,  1886,  1888,  1890,  1894,  1898,  1902,  1904,  1906,  1908,
    1912,  1916,  1920,  1924,  1928,  1932,  1934,  1938,  1942,  1946,
    1950,  1954,  1958,  1962,  1966,  1970,  1974,  1978,  1982,  1986,
    1990,  1994,  1998,  2002,  2006,  2010,  2014,  2018,  2022,  2026,
    2030,  2034,  2038,  2042,  2046,  2050,  2054,  2056,  2060,  2063,
    2066
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
     290,   294,   300,   302,   306,   308,   310,   311,   314,   316,
     318,   319,   322,   324,   325,   328,   330,   332,   334,   335,
     338,   338,   340,   340,   342,   342,   345,   344,   347,   347,
     351,   352,   353,   354,   357,   359,   363,   365,   366,   368,
     370,   372,   374,   376,   378,   380,   382,   384,   386,   388,
     390,   392,   394,   396,   398,   400,   402,   404,   406,   408,
     412,   415,   417,   421,   423,   425,   426,   429,   431,   433,
     435,   437,   441,   443,   445,   447,   449,   451,   455,   457,
     461,   463,   465,   469,   471,   473,   475,   477,   479,   481,
     483,   485,   489,   491,   495,   496,   499,   501,   503,   507,
     508,   511,   513,   515,   519,   520,   523,   524,   527,   529,
     531,   533,   537,   538,   541,   542,   543,   544,   545,   546,
     547,   548,   549,   550,   551,   552,   553,   554,   555,   556,
     557,   558,   559,   560,   561,   562,   563,   566,   568,   570,
     572,   574,   576,   580,   582,   584,   588,   590,   592,   596,
     598,   600,   604,   606,   612,   618,   628,   633,   640,   651,
     656,   667,   674,   683,   694,   709,   712,   714,   718,   726,
     736,   746,   749,   751,   755,   765,   777,   789,   791,   793,
     795,   797,   801,   802,   803,   804,   805,   807,   811,   813,
     815,   817,   821,   822,   825,   826,   827,   828,   829,   830,
     831,   832,   833,   834,   835,   836,   837,   838,   839,   840,
     841,   842,   843,   844,   845,   846,   847,   848,   849,   850,
     851,   852,   853,   854,   855,   856,   857,   858,   859,   860,
     861,   862,   863,   866,   868,   872,   873,   876,   877,   878,
     879,   880,   881,   882,   883,   884,   885,   886,   887,   888,
     891,   893,   897,   899,   903,   904,   907,   909,   911,   912,
     915,   917,   919,   921,   923,   925,   927,   931,   933,   935,
     937,   939,   943,   945,   946,   949,   951,   953,   957,   958,
     960,   964,   966,   970,   972,   974,   976,   978,   980,   984,
     986,   988,   990,   992,   994,   998,  1001,  1002,  1005,  1006,
    1007,  1010,  1012,  1014,  1016,  1020,  1022,  1026,  1027,  1029,
    1031,  1033,  1037,  1038,  1037,  1040,  1042,  1044,  1046,  1050,
    1051,  1054,  1055,  1058,  1059,  1060,  1061,  1062,  1063,  1064,
    1067,  1068,  1069,  1070,  1071,  1072,  1073,  1074,  1075,  1076,
    1079,  1080,  1083,  1085,  1089,  1090,  1091,  1092,  1095,  1096,
    1099,  1101,  1105,  1107,  1111,  1112,  1115,  1116,  1117,  1118,
    1119,  1120,  1121,  1122,  1123,  1124,  1125,  1126,  1127,  1128,
    1129,  1130,  1131,  1132,  1133,  1134,  1135,  1136,  1137,  1138,
    1139,  1140,  1141,  1144,  1145,  1148,  1149,  1150,  1151,  1152,
    1153,  1154,  1155,  1156,  1157,  1158,  1159,  1160,  1161,  1162,
    1164,  1165,  1166,  1167,  1168,  1169,  1170,  1172,  1175,  1176,
    1177,  1178,  1179,  1180,  1182,  1185,  1186,  1187,  1188,  1189,
    1190,  1191,  1192,  1193,  1194,  1195,  1196,  1197,  1198,  1199,
    1200,  1201,  1202,  1203,  1204,  1205,  1206,  1207,  1208,  1209,
    1211,  1214,  1215,  1216,  1217,  1218,  1219,  1220,  1221,  1222,
    1224,  1225,  1226,  1227,  1228,  1229,  1230,  1231,  1233,  1234,
    1235,  1236,  1237,  1238,  1239,  1240,  1241,  1242,  1243,  1244,
    1245,  1246,  1247,  1248,  1249,  1250,  1251,  1253,  1254,  1260,
    1261,  1265,  1266,  1267,  1268,  1271,  1279,  1280,  1289,  1291,
    1300
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
  const int parser::yylast_ = 1760;
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

