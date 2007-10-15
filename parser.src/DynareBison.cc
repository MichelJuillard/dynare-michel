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
#line 151 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 49:
#line 153 "DynareBison.yy"
    { driver.dsample((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 50:
#line 156 "DynareBison.yy"
    { driver.rplot(); ;}
    break;

  case 55:
#line 167 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 56:
#line 169 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 57:
#line 171 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 58:
#line 173 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 59:
#line 175 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 60:
#line 177 "DynareBison.yy"
    { driver.declare_endogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 61:
#line 181 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 62:
#line 183 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 63:
#line 185 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 64:
#line 187 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 65:
#line 189 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 66:
#line 191 "DynareBison.yy"
    { driver.declare_exogenous((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 67:
#line 195 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 68:
#line 197 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 69:
#line 199 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 70:
#line 201 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 71:
#line 203 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 72:
#line 205 "DynareBison.yy"
    { driver.declare_exogenous_det((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 73:
#line 209 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 74:
#line 211 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 75:
#line 213 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 76:
#line 215 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(3) - (2)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 77:
#line 217 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(4) - (3)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 78:
#line 219 "DynareBison.yy"
    { driver.declare_parameter((yysemantic_stack_[(2) - (1)].string_val), (yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 79:
#line 223 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 80:
#line 225 "DynareBison.yy"
    { driver.periods((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 81:
#line 229 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 82:
#line 231 "DynareBison.yy"
    { driver.cutoff((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 83:
#line 235 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 84:
#line 237 "DynareBison.yy"
    { driver.markowitz((yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 85:
#line 241 "DynareBison.yy"
    { driver.init_param((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 86:
#line 244 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 87:
#line 246 "DynareBison.yy"
    { (yyval.node_val) = driver.add_expression_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 88:
#line 248 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 89:
#line 250 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 90:
#line 252 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 91:
#line 254 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 92:
#line 256 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 93:
#line 258 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 94:
#line 260 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 95:
#line 262 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 96:
#line 264 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 97:
#line 266 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 98:
#line 268 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 99:
#line 270 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 100:
#line 272 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 101:
#line 274 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 102:
#line 276 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 103:
#line 278 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 104:
#line 280 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 105:
#line 282 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 106:
#line 284 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 107:
#line 286 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 108:
#line 288 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 109:
#line 290 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 110:
#line 292 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 111:
#line 294 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 112:
#line 296 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 113:
#line 298 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 114:
#line 300 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 115:
#line 302 "DynareBison.yy"
    { (yyval.node_val) = driver.add_unknown_function((yysemantic_stack_[(4) - (1)].string_val)); ;}
    break;

  case 116:
#line 306 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(1) - (1)].node_val)); ;}
    break;

  case 117:
#line 308 "DynareBison.yy"
    { driver.add_unknown_function_arg((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 118:
#line 312 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 119:
#line 314 "DynareBison.yy"
    { driver.end_initval(); ;}
    break;

  case 120:
#line 317 "DynareBison.yy"
    { driver.init_val_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 121:
#line 319 "DynareBison.yy"
    { driver.end_endval(); ;}
    break;

  case 124:
#line 325 "DynareBison.yy"
    { driver.init_val((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 125:
#line 327 "DynareBison.yy"
    { driver.end_histval(); ;}
    break;

  case 128:
#line 333 "DynareBison.yy"
    { driver.hist_val((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 131:
#line 340 "DynareBison.yy"
    { driver.init_compiler(0); ;}
    break;

  case 132:
#line 342 "DynareBison.yy"
    { driver.init_compiler(1); ;}
    break;

  case 133:
#line 344 "DynareBison.yy"
    { driver.init_compiler(2); ;}
    break;

  case 136:
#line 349 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 137:
#line 350 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 138:
#line 351 "DynareBison.yy"
    { driver.begin_model(); ;}
    break;

  case 139:
#line 352 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 140:
#line 353 "DynareBison.yy"
    { driver.begin_model(); driver.use_dll(); ;}
    break;

  case 141:
#line 354 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 142:
#line 356 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 143:
#line 357 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 144:
#line 358 "DynareBison.yy"
    { driver.begin_model(); driver.sparse_dll(); ;}
    break;

  case 145:
#line 359 "DynareBison.yy"
    { driver.reset_data_tree(); ;}
    break;

  case 150:
#line 369 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal((yysemantic_stack_[(4) - (1)].node_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 151:
#line 371 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_equal_with_zero_rhs((yysemantic_stack_[(2) - (1)].node_val)); ;}
    break;

  case 152:
#line 375 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(3) - (2)].node_val);;}
    break;

  case 154:
#line 378 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 155:
#line 380 "DynareBison.yy"
    { (yyval.node_val) = driver.add_constant((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 156:
#line 382 "DynareBison.yy"
    { (yyval.node_val) = driver.add_plus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 157:
#line 384 "DynareBison.yy"
    { (yyval.node_val) = driver.add_minus((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 158:
#line 386 "DynareBison.yy"
    { (yyval.node_val) = driver.add_divide((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 159:
#line 388 "DynareBison.yy"
    { (yyval.node_val) = driver.add_times((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 160:
#line 390 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 161:
#line 392 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 162:
#line 394 "DynareBison.yy"
    { (yyval.node_val) = driver.add_less_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 163:
#line 396 "DynareBison.yy"
    { (yyval.node_val) = driver.add_greater_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 164:
#line 398 "DynareBison.yy"
    { (yyval.node_val) = driver.add_equal_equal((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 165:
#line 400 "DynareBison.yy"
    { (yyval.node_val) = driver.add_different((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 166:
#line 402 "DynareBison.yy"
    { (yyval.node_val) = driver.add_power((yysemantic_stack_[(3) - (1)].node_val), (yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 167:
#line 404 "DynareBison.yy"
    { (yyval.node_val) = driver.add_uminus((yysemantic_stack_[(2) - (2)].node_val)); ;}
    break;

  case 168:
#line 406 "DynareBison.yy"
    { (yyval.node_val) = (yysemantic_stack_[(2) - (2)].node_val); ;}
    break;

  case 169:
#line 408 "DynareBison.yy"
    { (yyval.node_val) = driver.add_exp((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 170:
#line 410 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 171:
#line 412 "DynareBison.yy"
    { (yyval.node_val) = driver.add_log10((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 172:
#line 414 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 173:
#line 416 "DynareBison.yy"
    { (yyval.node_val) = driver.add_cos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 174:
#line 418 "DynareBison.yy"
    { (yyval.node_val) = driver.add_tan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 175:
#line 420 "DynareBison.yy"
    { (yyval.node_val) = driver.add_asin((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 176:
#line 422 "DynareBison.yy"
    { (yyval.node_val) = driver.add_acos((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 177:
#line 424 "DynareBison.yy"
    { (yyval.node_val) = driver.add_atan((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 178:
#line 426 "DynareBison.yy"
    { (yyval.node_val) = driver.add_sqrt((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 179:
#line 428 "DynareBison.yy"
    { (yyval.node_val) = driver.add_max((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 180:
#line 430 "DynareBison.yy"
    { (yyval.node_val) = driver.add_min((yysemantic_stack_[(6) - (3)].node_val) , (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 181:
#line 434 "DynareBison.yy"
    { driver.declare_and_init_model_local_variable((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 182:
#line 437 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 183:
#line 439 "DynareBison.yy"
    { (yyval.node_val) = driver.add_model_variable((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 184:
#line 442 "DynareBison.yy"
    { driver.end_shocks(); ;}
    break;

  case 185:
#line 444 "DynareBison.yy"
    { driver.end_mshocks(); ;}
    break;

  case 188:
#line 451 "DynareBison.yy"
    { driver.add_det_shock((yysemantic_stack_[(9) - (2)].string_val)); ;}
    break;

  case 189:
#line 453 "DynareBison.yy"
    { driver.add_stderr_shock((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 190:
#line 455 "DynareBison.yy"
    { driver.add_var_shock((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 191:
#line 457 "DynareBison.yy"
    { driver.add_covar_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 192:
#line 459 "DynareBison.yy"
    { driver.add_correl_shock((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 193:
#line 463 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 194:
#line 465 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 195:
#line 467 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 196:
#line 469 "DynareBison.yy"
    { driver.add_period((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 197:
#line 473 "DynareBison.yy"
    { driver.do_sigma_e(); ;}
    break;

  case 198:
#line 477 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(3) - (3)].node_val));;}
    break;

  case 199:
#line 479 "DynareBison.yy"
    {driver.add_value((yysemantic_stack_[(1) - (1)].node_val));;}
    break;

  case 200:
#line 483 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 201:
#line 485 "DynareBison.yy"
    { driver.end_of_row(); ;}
    break;

  case 202:
#line 489 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 203:
#line 491 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 204:
#line 493 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 205:
#line 495 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 206:
#line 497 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 207:
#line 499 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 208:
#line 501 "DynareBison.yy"
    { driver.add_to_row((yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 209:
#line 503 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 210:
#line 505 "DynareBison.yy"
    { driver.add_to_row_const((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 211:
#line 509 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 212:
#line 511 "DynareBison.yy"
    { driver.steady(); ;}
    break;

  case 216:
#line 521 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 217:
#line 523 "DynareBison.yy"
    { driver.check(); ;}
    break;

  case 221:
#line 533 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 222:
#line 535 "DynareBison.yy"
    { driver.simulate(); ;}
    break;

  case 227:
#line 547 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 228:
#line 549 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 229:
#line 551 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 230:
#line 553 "DynareBison.yy"
    { driver.stoch_simul(); ;}
    break;

  case 256:
#line 586 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(2) - (2)].string_val)); ;}
    break;

  case 257:
#line 588 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (4)].string_val)); ;}
    break;

  case 258:
#line 590 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 259:
#line 592 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 260:
#line 594 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 261:
#line 596 "DynareBison.yy"
    { driver.add_tmp_var((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 262:
#line 600 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 263:
#line 602 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 264:
#line 604 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 265:
#line 608 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 266:
#line 610 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 267:
#line 612 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 268:
#line 615 "DynareBison.yy"
    { driver.estimated_params(); ;}
    break;

  case 269:
#line 618 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 270:
#line 620 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 272:
#line 626 "DynareBison.yy"
    {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *(yysemantic_stack_[(2) - (2)].string_val);
                    delete (yysemantic_stack_[(2) - (2)].string_val);
                  ;}
    break;

  case 273:
#line 632 "DynareBison.yy"
    {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 274:
#line 638 "DynareBison.yy"
    {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *(yysemantic_stack_[(4) - (2)].string_val);
                    driver.estim_params.name2 = *(yysemantic_stack_[(4) - (4)].string_val);
                    delete (yysemantic_stack_[(4) - (2)].string_val);
                    delete (yysemantic_stack_[(4) - (4)].string_val);
                  ;}
    break;

  case 275:
#line 648 "DynareBison.yy"
    {
                    driver.estim_params.prior = *(yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                  ;}
    break;

  case 276:
#line 653 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.prior = *(yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                  ;}
    break;

  case 277:
#line 660 "DynareBison.yy"
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

  case 278:
#line 671 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(1) - (1)].string_val);
                    delete (yysemantic_stack_[(1) - (1)].string_val);
                  ;}
    break;

  case 279:
#line 676 "DynareBison.yy"
    {
                    driver.estim_params.init_val = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.low_bound = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.up_bound = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 280:
#line 687 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(3) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(3) - (3)].string_val);
                    delete (yysemantic_stack_[(3) - (1)].string_val);
                    delete (yysemantic_stack_[(3) - (3)].string_val);
                  ;}
    break;

  case 281:
#line 694 "DynareBison.yy"
    {
                    driver.estim_params.mean = *(yysemantic_stack_[(5) - (1)].string_val);
                    driver.estim_params.std = *(yysemantic_stack_[(5) - (3)].string_val);
                    driver.estim_params.p3 = *(yysemantic_stack_[(5) - (5)].string_val);
                    delete (yysemantic_stack_[(5) - (1)].string_val);
                    delete (yysemantic_stack_[(5) - (3)].string_val);
                    delete (yysemantic_stack_[(5) - (5)].string_val);
                  ;}
    break;

  case 282:
#line 703 "DynareBison.yy"
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

  case 283:
#line 714 "DynareBison.yy"
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

  case 284:
#line 729 "DynareBison.yy"
    { driver.estimated_params_init(); ;}
    break;

  case 285:
#line 732 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 286:
#line 734 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 287:
#line 738 "DynareBison.yy"
    {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *(yysemantic_stack_[(5) - (2)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(5) - (4)].string_val);
                        delete (yysemantic_stack_[(5) - (2)].string_val);
                        delete (yysemantic_stack_[(5) - (4)].string_val);
                      ;}
    break;

  case 288:
#line 746 "DynareBison.yy"
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

  case 289:
#line 756 "DynareBison.yy"
    {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *(yysemantic_stack_[(4) - (1)].string_val);
                        driver.estim_params.init_val = *(yysemantic_stack_[(4) - (3)].string_val);
                        delete (yysemantic_stack_[(4) - (1)].string_val);
                        delete (yysemantic_stack_[(4) - (3)].string_val);
                      ;}
    break;

  case 290:
#line 766 "DynareBison.yy"
    { driver.estimated_params_bounds(); ;}
    break;

  case 291:
#line 769 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 292:
#line 771 "DynareBison.yy"
    { driver.add_estimated_params_element(); ;}
    break;

  case 293:
#line 775 "DynareBison.yy"
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

  case 294:
#line 785 "DynareBison.yy"
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

  case 295:
#line 797 "DynareBison.yy"
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

  case 296:
#line 809 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 297:
#line 811 "DynareBison.yy"
    { (yyval.string_val) = new string("2"); ;}
    break;

  case 298:
#line 813 "DynareBison.yy"
    { (yyval.string_val) = new string("3"); ;}
    break;

  case 299:
#line 815 "DynareBison.yy"
    { (yyval.string_val) = new string("4"); ;}
    break;

  case 300:
#line 817 "DynareBison.yy"
    { (yyval.string_val) = new string("5"); ;}
    break;

  case 301:
#line 820 "DynareBison.yy"
    { (yyval.string_val) = new string("NaN"); ;}
    break;

  case 305:
#line 825 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 306:
#line 827 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "-"); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val); ;}
    break;

  case 307:
#line 831 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 308:
#line 833 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 309:
#line 835 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 310:
#line 837 "DynareBison.yy"
    { driver.run_estimation(); ;}
    break;

  case 352:
#line 886 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 353:
#line 888 "DynareBison.yy"
    { driver.run_prior_analysis(); ;}
    break;

  case 369:
#line 911 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 370:
#line 913 "DynareBison.yy"
    { driver.run_posterior_analysis(); ;}
    break;

  case 371:
#line 917 "DynareBison.yy"
    { driver.optim_options_string((yysemantic_stack_[(7) - (2)].string_val), (yysemantic_stack_[(7) - (6)].string_val)); ;}
    break;

  case 372:
#line 919 "DynareBison.yy"
    { driver.optim_options_num((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (5)].string_val)); ;}
    break;

  case 375:
#line 926 "DynareBison.yy"
    { driver.set_varobs(); ;}
    break;

  case 376:
#line 928 "DynareBison.yy"
    { driver.set_trends(); ;}
    break;

  case 379:
#line 934 "DynareBison.yy"
    { driver.set_trend_element((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].node_val)); ;}
    break;

  case 380:
#line 936 "DynareBison.yy"
    { driver.set_unit_root_vars(); ;}
    break;

  case 381:
#line 938 "DynareBison.yy"
    { driver.optim_weights(); ;}
    break;

  case 382:
#line 941 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(4) - (2)].string_val), (yysemantic_stack_[(4) - (3)].node_val)); ;}
    break;

  case 383:
#line 943 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(6) - (2)].string_val), (yysemantic_stack_[(6) - (4)].string_val), (yysemantic_stack_[(6) - (5)].node_val)); ;}
    break;

  case 384:
#line 945 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(3) - (1)].string_val), (yysemantic_stack_[(3) - (2)].node_val)); ;}
    break;

  case 385:
#line 947 "DynareBison.yy"
    { driver.set_optim_weights((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (3)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 386:
#line 950 "DynareBison.yy"
    { driver.set_osr_params(); ;}
    break;

  case 387:
#line 953 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 388:
#line 955 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 389:
#line 957 "DynareBison.yy"
    { driver.run_osr(); ;}
    break;

  case 390:
#line 959 "DynareBison.yy"
    {driver.run_osr(); ;}
    break;

  case 391:
#line 962 "DynareBison.yy"
    { driver.run_calib_var(); ;}
    break;

  case 394:
#line 969 "DynareBison.yy"
    { driver.set_calib_var((yysemantic_stack_[(5) - (1)].string_val), (yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].node_val)); ;}
    break;

  case 395:
#line 971 "DynareBison.yy"
    { driver.set_calib_covar((yysemantic_stack_[(7) - (1)].string_val), (yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (4)].string_val), (yysemantic_stack_[(7) - (6)].node_val)); ;}
    break;

  case 396:
#line 973 "DynareBison.yy"
    { driver.set_calib_ac((yysemantic_stack_[(9) - (2)].string_val), (yysemantic_stack_[(9) - (4)].string_val), (yysemantic_stack_[(9) - (6)].string_val), (yysemantic_stack_[(9) - (8)].node_val)); ;}
    break;

  case 397:
#line 976 "DynareBison.yy"
    { (yyval.string_val) = new string("1"); ;}
    break;

  case 398:
#line 978 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 399:
#line 980 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(3) - (2)].string_val); ;}
    break;

  case 400:
#line 984 "DynareBison.yy"
    { driver.run_calib(0); ;}
    break;

  case 401:
#line 986 "DynareBison.yy"
    { driver.run_calib(1); ;}
    break;

  case 402:
#line 990 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 403:
#line 992 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 404:
#line 994 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 405:
#line 996 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 406:
#line 998 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 407:
#line 1000 "DynareBison.yy"
    { driver.run_dynatype((yysemantic_stack_[(5) - (2)].string_val),(yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 408:
#line 1004 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (3)].string_val)); ;}
    break;

  case 409:
#line 1006 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(6) - (3)].string_val)); ;}
    break;

  case 410:
#line 1008 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 411:
#line 1010 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(7) - (3)].string_val), (yysemantic_stack_[(7) - (5)].string_val)); ;}
    break;

  case 412:
#line 1012 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(8) - (3)].string_val), (yysemantic_stack_[(8) - (5)].string_val)); ;}
    break;

  case 413:
#line 1014 "DynareBison.yy"
    { driver.run_dynasave((yysemantic_stack_[(5) - (2)].string_val), (yysemantic_stack_[(5) - (4)].string_val)); ;}
    break;

  case 414:
#line 1018 "DynareBison.yy"
    { driver.run_model_comparison(); ;}
    break;

  case 420:
#line 1030 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(1) - (1)].string_val)); ;}
    break;

  case 421:
#line 1032 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 422:
#line 1034 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(4) - (1)].string_val), (yysemantic_stack_[(4) - (3)].string_val)); ;}
    break;

  case 423:
#line 1036 "DynareBison.yy"
    { driver.add_mc_filename((yysemantic_stack_[(6) - (3)].string_val), (yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 424:
#line 1040 "DynareBison.yy"
    { (yyval.string_val) = (yysemantic_stack_[(1) - (1)].string_val); ;}
    break;

  case 425:
#line 1042 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val)); delete (yysemantic_stack_[(2) - (2)].string_val); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;

  case 427:
#line 1047 "DynareBison.yy"
    { (yyval.string_val) = new string("\\"); ;}
    break;

  case 428:
#line 1049 "DynareBison.yy"
    { (yyval.string_val) = new string("/"); ;}
    break;

  case 429:
#line 1051 "DynareBison.yy"
    { (yyval.string_val) = new string(":"); ;}
    break;

  case 430:
#line 1053 "DynareBison.yy"
    { (yyval.string_val) = new string("."); ;}
    break;

  case 431:
#line 1056 "DynareBison.yy"
    { driver.begin_planner_objective(); ;}
    break;

  case 432:
#line 1057 "DynareBison.yy"
    { driver.end_planner_objective((yysemantic_stack_[(3) - (3)].node_val)); ;}
    break;

  case 434:
#line 1060 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 435:
#line 1062 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 436:
#line 1064 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 437:
#line 1066 "DynareBison.yy"
    { driver.ramsey_policy(); ;}
    break;

  case 461:
#line 1103 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 462:
#line 1105 "DynareBison.yy"
    { driver.bvar_density((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 469:
#line 1119 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(3) - (2)].string_val)); ;}
    break;

  case 470:
#line 1121 "DynareBison.yy"
    { driver.bvar_forecast((yysemantic_stack_[(6) - (5)].string_val)); ;}
    break;

  case 471:
#line 1125 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 472:
#line 1127 "DynareBison.yy"
    { driver.dynare_sensitivity(); ;}
    break;

  case 504:
#line 1167 "DynareBison.yy"
    { driver.option_num("dr_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 505:
#line 1168 "DynareBison.yy"
    { driver.option_num("solve_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 506:
#line 1169 "DynareBison.yy"
    { driver.option_num("simul_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 507:
#line 1170 "DynareBison.yy"
    { driver.linear(); ;}
    break;

  case 508:
#line 1171 "DynareBison.yy"
    { driver.option_num("order", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 509:
#line 1172 "DynareBison.yy"
    { driver.option_num("replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 510:
#line 1173 "DynareBison.yy"
    { driver.option_num("drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 511:
#line 1174 "DynareBison.yy"
    { driver.option_num("ar", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 512:
#line 1175 "DynareBison.yy"
    { driver.option_num("nocorr", "1"); ;}
    break;

  case 513:
#line 1176 "DynareBison.yy"
    { driver.option_num("nofunctions", "1"); ;}
    break;

  case 514:
#line 1177 "DynareBison.yy"
    { driver.option_num("nomoments", "1"); ;}
    break;

  case 515:
#line 1178 "DynareBison.yy"
    { driver.option_num("irf", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 516:
#line 1179 "DynareBison.yy"
    { driver.option_num("hp_filter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 517:
#line 1180 "DynareBison.yy"
    { driver.option_num("hp_ngrid", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 518:
#line 1182 "DynareBison.yy"
    { driver.option_num("periods", (yysemantic_stack_[(3) - (3)].string_val)); driver.option_num("simul", "1"); ;}
    break;

  case 519:
#line 1183 "DynareBison.yy"
    { driver.option_num("cutoff", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 520:
#line 1184 "DynareBison.yy"
    { driver.option_num("markowitz", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 521:
#line 1185 "DynareBison.yy"
    { driver.option_num("simul", "1"); ;}
    break;

  case 522:
#line 1186 "DynareBison.yy"
    { driver.option_num("simul_seed", (yysemantic_stack_[(3) - (3)].string_val));}
    break;

  case 523:
#line 1187 "DynareBison.yy"
    { driver.option_num("qz_criterium", (yysemantic_stack_[(3) - (3)].string_val)) ;}
    break;

  case 524:
#line 1188 "DynareBison.yy"
    { driver.option_str("datafile", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 525:
#line 1190 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 526:
#line 1192 "DynareBison.yy"
    { driver.option_num("nobs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 527:
#line 1194 "DynareBison.yy"
    { driver.option_num("first_obs", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 528:
#line 1195 "DynareBison.yy"
    { driver.option_num("prefilter", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 529:
#line 1196 "DynareBison.yy"
    { driver.option_num("presample", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 530:
#line 1197 "DynareBison.yy"
    { driver.option_num("lik_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 531:
#line 1198 "DynareBison.yy"
    { driver.option_num("lik_init", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 532:
#line 1200 "DynareBison.yy"
    { driver.option_num("nograph","1"); ;}
    break;

  case 533:
#line 1202 "DynareBison.yy"
    { driver.option_num("nograph", "0"); ;}
    break;

  case 534:
#line 1204 "DynareBison.yy"
    { driver.option_num("conf_sig", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 535:
#line 1205 "DynareBison.yy"
    { driver.option_num("mh_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 536:
#line 1206 "DynareBison.yy"
    { driver.option_num("mh_drop", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 537:
#line 1207 "DynareBison.yy"
    { driver.option_num("mh_jscale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 539:
#line 1209 "DynareBison.yy"
    { driver.option_num("mh_init_scale", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 540:
#line 1210 "DynareBison.yy"
    { driver.option_str("mode_file", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 541:
#line 1211 "DynareBison.yy"
    { driver.option_num("mode_compute", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 542:
#line 1212 "DynareBison.yy"
    { driver.option_num("mode_check", "1"); ;}
    break;

  case 543:
#line 1213 "DynareBison.yy"
    { driver.option_num("prior_trunc", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 544:
#line 1214 "DynareBison.yy"
    { driver.option_num("mh_mode", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 545:
#line 1215 "DynareBison.yy"
    { driver.option_num("mh_nblck", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 546:
#line 1216 "DynareBison.yy"
    { driver.option_num("load_mh_file", "1"); ;}
    break;

  case 547:
#line 1217 "DynareBison.yy"
    { driver.option_num("loglinear", "1"); ;}
    break;

  case 548:
#line 1218 "DynareBison.yy"
    { driver.option_num("nodiagnostic", "1"); ;}
    break;

  case 549:
#line 1219 "DynareBison.yy"
    { driver.option_num("bayesian_irf", "1"); ;}
    break;

  case 550:
#line 1220 "DynareBison.yy"
    { driver.option_num("TeX", "1"); ;}
    break;

  case 551:
#line 1221 "DynareBison.yy"
    { driver.option_num("forecast", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 552:
#line 1222 "DynareBison.yy"
    { driver.option_num("smoother", "1"); ;}
    break;

  case 553:
#line 1223 "DynareBison.yy"
    { driver.option_num("moments_varendo", "1"); ;}
    break;

  case 554:
#line 1224 "DynareBison.yy"
    { driver.option_num("filtered_vars", "1"); ;}
    break;

  case 555:
#line 1225 "DynareBison.yy"
    { driver.option_num("relative_irf", "1"); ;}
    break;

  case 556:
#line 1226 "DynareBison.yy"
    { driver.option_num("kalman_algo", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 557:
#line 1227 "DynareBison.yy"
    { driver.option_num("kalman_tol", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 558:
#line 1229 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "Laplace"); ;}
    break;

  case 559:
#line 1231 "DynareBison.yy"
    { driver.option_str("model_comparison_approximation", "MODIFIEDHARMONICMEAN"); ;}
    break;

  case 560:
#line 1233 "DynareBison.yy"
    { driver.option_num("noprint", "0"); ;}
    break;

  case 561:
#line 1234 "DynareBison.yy"
    { driver.option_num("noprint", "1"); ;}
    break;

  case 562:
#line 1235 "DynareBison.yy"
    { driver.option_str("xls_sheet", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 563:
#line 1236 "DynareBison.yy"
    { driver.option_str("xls_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 564:
#line 1237 "DynareBison.yy"
    { driver.option_num("filter_step_ahead", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 565:
#line 1238 "DynareBison.yy"
    { driver.option_num("noconstant", "0"); ;}
    break;

  case 566:
#line 1239 "DynareBison.yy"
    { driver.option_num("noconstant", "1"); ;}
    break;

  case 567:
#line 1240 "DynareBison.yy"
    { driver.option_num("mh_recover", "1"); ;}
    break;

  case 568:
#line 1241 "DynareBison.yy"
    { driver.option_num("planner_discount",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 569:
#line 1243 "DynareBison.yy"
    { driver.option_num("bvar_prior_tau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 570:
#line 1244 "DynareBison.yy"
    { driver.option_num("bvar_prior_decay", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 571:
#line 1245 "DynareBison.yy"
    { driver.option_num("bvar_prior_lambda", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 572:
#line 1246 "DynareBison.yy"
    { driver.option_num("bvar_prior_mu", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 573:
#line 1247 "DynareBison.yy"
    { driver.option_num("bvar_prior_omega", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 574:
#line 1248 "DynareBison.yy"
    { driver.option_num("bvar_prior_flat", "1"); ;}
    break;

  case 575:
#line 1249 "DynareBison.yy"
    { driver.option_num("bvar_prior_train", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 576:
#line 1250 "DynareBison.yy"
    { driver.option_num("bvar_replic", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 577:
#line 1252 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 578:
#line 1253 "DynareBison.yy"
    { driver.option_num("morris", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 579:
#line 1254 "DynareBison.yy"
    { driver.option_num("stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 580:
#line 1255 "DynareBison.yy"
    { driver.option_num("redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 581:
#line 1256 "DynareBison.yy"
    { driver.option_num("pprior", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 582:
#line 1257 "DynareBison.yy"
    { driver.option_num("prior_range", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 583:
#line 1258 "DynareBison.yy"
    { driver.option_num("ppost", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 584:
#line 1259 "DynareBison.yy"
    { driver.option_num("ilptau", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 585:
#line 1260 "DynareBison.yy"
    { driver.option_num("glue", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 586:
#line 1261 "DynareBison.yy"
    { driver.option_num("morris_nliv", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 587:
#line 1262 "DynareBison.yy"
    { driver.option_num("morris_ntra", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 588:
#line 1263 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 589:
#line 1264 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 590:
#line 1265 "DynareBison.yy"
    { driver.option_num("load_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 591:
#line 1266 "DynareBison.yy"
    { driver.option_num("load_stab", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 592:
#line 1267 "DynareBison.yy"
    { driver.option_num("identification", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 593:
#line 1268 "DynareBison.yy"
    { driver.option_num("ksstat", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 594:
#line 1269 "DynareBison.yy"
    { driver.option_num("logtrans_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 595:
#line 1270 "DynareBison.yy"
    { driver.option_num("threshold_redfor",(yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 596:
#line 1272 "DynareBison.yy"
    { driver.option_num("ksstat_redfrom", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 597:
#line 1273 "DynareBison.yy"
    { driver.option_num("alpha2_redform", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 598:
#line 1279 "DynareBison.yy"
    { driver.option_num("rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 599:
#line 1280 "DynareBison.yy"
    { driver.option_num("lik_only", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 600:
#line 1284 "DynareBison.yy"
    { driver.option_num("pfilt_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 601:
#line 1285 "DynareBison.yy"
    { driver.option_num("istart_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 602:
#line 1286 "DynareBison.yy"
    { driver.option_num("alpha_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 603:
#line 1287 "DynareBison.yy"
    { driver.option_num("alpha2_rmse", (yysemantic_stack_[(3) - (3)].string_val)); ;}
    break;

  case 604:
#line 1291 "DynareBison.yy"
    {
          (yysemantic_stack_[(3) - (1)].string_val)->append(":");
          (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
          delete (yysemantic_stack_[(3) - (3)].string_val);
          (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val);
        ;}
    break;

  case 606:
#line 1300 "DynareBison.yy"
    {
                 (yysemantic_stack_[(3) - (1)].string_val)->append(":");
                 (yysemantic_stack_[(3) - (1)].string_val)->append(*(yysemantic_stack_[(3) - (3)].string_val));
                 delete (yysemantic_stack_[(3) - (3)].string_val);
                 (yyval.string_val) = (yysemantic_stack_[(3) - (1)].string_val); 
               ;}
    break;

  case 607:
#line 1309 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (2)].string_val)->insert(0, "["); (yyval.string_val) = (yysemantic_stack_[(2) - (2)].string_val);;}
    break;

  case 608:
#line 1311 "DynareBison.yy"
    {
              (yysemantic_stack_[(2) - (1)].string_val)->append(" ");
              (yysemantic_stack_[(2) - (1)].string_val)->append(*(yysemantic_stack_[(2) - (2)].string_val));
              delete (yysemantic_stack_[(2) - (2)].string_val);
              (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val);
            ;}
    break;

  case 609:
#line 1319 "DynareBison.yy"
    { (yysemantic_stack_[(2) - (1)].string_val)->append("]"); (yyval.string_val) = (yysemantic_stack_[(2) - (1)].string_val); ;}
    break;


    /* Line 675 of lalr1.cc.  */
#line 2399 "DynareBison.cc"
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
  const short int parser::yypact_ninf_ = -1090;
  const short int
  parser::yypact_[] =
  {
       914,    18,    28,   107,  -107,   393,    92,    44,    29,    34,
     -66,     4,   -64,   -38,    -4,   105,   415,   263,   417,   -58,
     109,   122,   179,   210,    24,   338,   345,    96, -1090,   247,
     250,   338,   235,   418,   426,   431,    82,    88,   338,   391,
     398,   412,   338,   434,   759, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
     312,  1248,   332,   662, -1090,   475,   100, -1090,   424,   523,
     422,    27,   -75,   541,   103,   555,   558,   608, -1090,  1225,
      -1,   385,   440,   537,   563,   558,   617,   614,   461, -1090,
     334,    35,   155,  1055,   579,   580, -1090,  1363,     6,   201,
     540,   212,   612,   468,  1082,   780,   780,   217,   155,   464,
   -1090,    68, -1090,   424, -1090,  1363,   219, -1090,  1310,   229,
     233,   553,   257,   554,   261,   559,   262,   275, -1090,  1944,
   -1090, -1090, -1090,   636, -1090,   648,   651,   653,   659,   661,
   -1090,   664,   667,   668, -1090,   669,   671,   673,   674, -1090,
     551,   507, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,   676,
     677,   679, -1090,   569,   514, -1090, -1090, -1090,   517,   637,
     -55,   408, -1090,   687,   -36, -1090, -1090,   522, -1090,   532,
   -1090, -1090,   650,   -80, -1090,   652,   251,   700,    69, -1090,
     654, -1090,   702, -1090, -1090,   705,   706,   708,   709,   710,
   -1090, -1090,   711,   712,   715,   716,   717,   718, -1090, -1090,
     724,   725, -1090, -1090, -1090,   727,   728, -1090, -1090,   -27,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
     729,   681, -1090,   686, -1090,   688,   444, -1090,   628,   692,
     631,   697,   479, -1090,   698,   638,   699,   576, -1090,   584,
      78, -1090,    81,   753,   589,   602, -1090,   911, -1090,   -26,
     601,   603,   770, -1090, -1090,   -25, -1090, -1090, -1090, -1090,
     723,   726,   140, -1090, -1090, -1090,   607,   610,   611,  1055,
    1055,   613,   615,   616,   618,   621,   624,   625,   627,   630,
     632,  1055,  1020,   634,   118, -1090,   934,   234,   779,   781,
     795,   805,   814,   822, -1090, -1090, -1090,   829,   830,   831,
   -1090,   834, -1090,   835,   841,   675, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090,   749,   796, -1090,   678, -1090, -1090, -1090,   682,
     684,   689,  1082,  1082,   694,   696,   701,   707,   714,   719,
     720,   730,   734,   736,  1082,  1928, -1090,   -21, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090,    -7, -1090,   386,    32,   281, -1090, -1090, -1090,
     282, -1090, -1090,   285, -1090, -1090,   853, -1090,   286, -1090,
   -1090, -1090, -1090, -1090,   762,   816, -1090, -1090,   773,   818,
   -1090, -1090,   776,   823, -1090, -1090,   873,   875,   876,   877,
     878,   879,   884,   885,   890,   891,   892,   894,   895,   896,
     904,   906,   907,   908,   912,   919,   921,   922,   923,   924,
     925,   926,   927,   746,   819, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090,   404,   357,   404,   915,   357,   916,   883,   917,
      26,   918,   928,   888,   898,  1248,   930,   931,   404,   933,
     662,   935,   774,   778,   905,   406,   941, -1090, -1090,   937,
     424,   790, -1090, -1090,   798,    49,   913,   800,    54,   920,
    1055, -1090, -1090, -1090,   797,   943,   949,   952,   953,   954,
     404,   404,   404,   957,   959,   967,   971,   942,   832,   404,
    1225,    58,   947,   990,   886, -1090, -1090, -1090,   389,   887,
      71,   889, -1090, -1090,   893,    71,   897, -1090, -1090,   384,
   -1090, -1090, -1090,   958,   840, -1090,   961,    33, -1090,   243,
   -1090,   587, -1090,   845,   850,   365,    35,    90,   909,    37,
   -1090, -1090,  1055,  1055,  1055,   900,   456,  1055,  1055,  1055,
    1055,  1055,  1055,  1055,  1055,  1055,  1055,   703,  1055,  1055,
    1055,  1055,  1055,  1055,  1055,  1055,  1055,  1055,  1055, -1090,
    1055, -1090, -1090,   969,  1702, -1090,  1032,   983,   404,  1001,
    1003,  1004,  1007,  1008,  1009,   404,  1010,  1021,  1022,    62,
   -1090,   945, -1090,  1082,  1082,   384,   932,   462,  1082,  1082,
    1082,  1082,  1082,  1082,  1082,  1082,  1082,  1082,  1002,  1082,
    1082,  1082,  1082,  1082,  1082,  1082,  1082,  1082,  1082,  1082,
     899,   780,    64,    70, -1090, -1090, -1090,  1055,   395,    25,
      68,   902,   424,   903,  1363,    75,   404,  1310,    85, -1090,
     950, -1090,   968, -1090,   973,  1024,  1033,  1044,  1052,  1054,
    1056,  1057,  1058,  1059,  1061,  1064,  1065,  1066,  1067,  1069,
     404,   404,  1070,   797,   404,   404,  1071,  1072,   404,  1074,
     404,   404,   929,  1944, -1090, -1090, -1090, -1090,  1084,  1087,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,  1079,    15,
   -1090, -1090, -1090, -1090,   948, -1090, -1090,   936, -1090, -1090,
   -1090, -1090,   939, -1090,  1080,   940,   955,   970,  1055, -1090,
   -1090, -1090, -1090, -1090,   278,   972, -1090, -1090,   292,   974,
    1721, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090,   963, -1090, -1090, -1090,   293,
   -1090,  1050,  1073, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090,   493,   975,   992,  1014,  1090,  1035,    71,  1098,   982,
      71, -1090,  1102,  1131,   991, -1090,   558,  1160, -1090, -1090,
   -1090,  1082, -1090, -1090, -1090,  1161, -1090,   299, -1090, -1090,
   -1090,   996, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090,   260,    59, -1090,  1114,  1055,  1115,   189,  1946,
    2006,  2042,   344,  1232,  1297,  1343,  1360,  1372,  1384,  1400,
    1412,  1429,  1441, -1090,   501,   501,   501,   501,   501,   501,
     456,   456,   900,   900,  1049,  1457,  1055, -1090,  1120,  1735,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090,   294, -1090,  2018,  2030,  1015,  1469,  1481,
    1498,  1514,  1526,  1538,  1550,  1564,  1583,  1595, -1090,   543,
     543,   543,   543,   543,   543,   462,   462,   932,   932,  1049,
   -1090, -1090, -1090,   305, -1090,   310,  1607,    32,  1005, -1090,
   -1090,    40,  1055, -1090, -1090, -1090, -1090, -1090, -1090,   320,
   -1090, -1090, -1090,   321, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,  1012, -1090,
   -1090, -1090,  1123, -1090, -1090,  1016,  1185, -1090, -1090,  1747,
   -1090,    91, -1090,   101, -1090,  1138, -1090,   352, -1090, -1090,
   -1090, -1090, -1090, -1090,    71,   389,  1078,    71,  1081,  1085,
   -1090,  1025, -1090, -1090,  1190,   413,  1082,  1763,   404,   587,
   -1090,   911,   911,   911,    90, -1090,    71, -1090,  1193,  1782,
    1198,  1191,  1055,  1055,  1055,  1055, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090,  1045,  1802,  1055,
   -1090, -1090,  1082,  1082, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,    25, -1090,
   -1090, -1090,  1055,  1621, -1090, -1090,  1192, -1090,   940,  1055,
   -1090, -1090,   328, -1090,   341,  1041,   963, -1090, -1090,  1105,
    1107,  1108,    71,  1068,    71,    71, -1090,  1055, -1090,  1818,
   -1090, -1090, -1090,  1075,    56,   196,   213,    72,  1053,  1055,
   -1090,  1055,  1094,   270,  1830,  1633,  1652,  2042, -1090, -1090,
    1845,  1664,  1678,  1690, -1090, -1090,  1220,  1857, -1090, -1090,
    1112, -1090,    71,    71,    71,  1126, -1090,  1076,  1077,  1873,
   -1090,   911, -1090, -1090, -1090,    71, -1090,  1885,  1900,  1221,
    1226,  1151, -1090, -1090, -1090, -1090, -1090, -1090, -1090,  1055,
   -1090,    48,  1141, -1090,  1144,    71, -1090, -1090, -1090,   571,
    1100, -1090, -1090, -1090,  1099,  1055,  1912,  1222, -1090,    71,
     104,  1104, -1090, -1090,  1254,  2042,   280, -1090,  1101,  1164,
    1167, -1090, -1090,  1055, -1090, -1090,    71,    71,  2042,  1171,
   -1090,    71, -1090
  };

  /* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
     doesn't specify something else to do.  Zero means the default is an
     error.  */
  const unsigned short int
  parser::yydefact_[] =
  {
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   431,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     2,     4,    29,    30,    45,
      46,    47,    44,     5,     6,     7,    12,     9,    10,    11,
       8,    13,    14,    15,    16,    17,    18,    19,    23,    25,
      24,    20,    21,    22,    26,    27,    28,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
       0,     0,     0,     0,   400,     0,     0,   216,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   260,   307,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   136,
       0,     0,     0,     0,     0,     0,   387,     0,     0,     0,
      75,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     221,     0,   211,     0,   227,     0,     0,   434,     0,     0,
       0,    57,     0,    63,     0,    69,     0,     0,   471,     0,
       1,     3,   461,     0,   574,     0,     0,     0,     0,     0,
     565,     0,     0,     0,   566,     0,     0,     0,     0,   449,
     460,     0,   450,   455,   453,   456,   454,   451,   452,   457,
     458,   442,   443,   444,   445,   446,   447,   448,   469,     0,
       0,     0,   463,   468,     0,   465,   464,   466,     0,     0,
     397,     0,   393,     0,     0,   219,   220,     0,    81,     0,
      48,   410,     0,     0,   404,     0,     0,     0,     0,   123,
       0,   549,     0,   554,   533,     0,     0,     0,     0,     0,
     546,   547,     0,     0,     0,     0,     0,     0,   567,   542,
       0,     0,   553,   548,   532,     0,     0,   552,   550,     0,
     312,   348,   337,   313,   314,   315,   316,   317,   318,   319,
     320,   321,   322,   323,   324,   325,   326,   327,   328,   329,
     330,   331,   332,   333,   334,   335,   336,   338,   339,   340,
     341,   342,   343,   344,   345,   346,   347,   349,   350,   351,
     256,     0,   309,     0,   273,     0,     0,   270,     0,     0,
       0,     0,     0,   292,     0,     0,     0,     0,   286,     0,
       0,   127,     0,     0,     0,     0,    83,     0,   507,     0,
       0,     0,     0,   561,   560,     0,   416,   417,   418,   419,
       0,     0,     0,   187,    88,    89,     0,     0,    87,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   378,     0,     0,     0,     0,
       0,     0,     0,     0,   512,   513,   514,     0,     0,     0,
     555,     0,   521,     0,     0,     0,   233,   234,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   245,   247,   248,
     249,   250,   251,   252,   253,   244,   246,   254,   255,   389,
     386,    78,    73,     0,    54,     0,    79,   154,   155,     0,
       0,   182,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   432,   153,     0,   355,   360,
     356,   357,   358,   359,   361,   362,   363,   364,   365,   366,
     367,   368,     0,    50,     0,     0,     0,   224,   225,   226,
       0,   214,   215,     0,   232,   229,     0,   440,     0,   439,
     441,   436,   380,    60,    55,     0,    51,    66,    61,     0,
      52,    72,    67,     0,    53,   375,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   474,   475,   476,   477,   478,   479,
     480,   481,   482,   483,   484,   485,   486,   487,   488,   489,
     490,   491,   492,   501,   493,   494,   495,   496,   497,   498,
     499,   500,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   391,   392,     0,
       0,     0,    82,    49,     0,     0,     0,     0,     0,     0,
       0,   121,   122,   261,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   258,     0,   272,   268,   269,   301,     0,
     301,     0,   290,   291,     0,   301,     0,   284,   285,     0,
     125,   126,   118,     0,     0,    84,     0,     0,   148,     0,
     149,     0,   144,     0,     0,     0,     0,     0,     0,     0,
     185,   186,     0,     0,     0,   101,   102,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    85,
       0,   376,   377,     0,     0,   381,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      76,    74,    80,     0,     0,     0,   167,   168,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   184,   209,   210,     0,     0,   201,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    58,
      56,    64,    62,    70,    68,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   503,   502,   570,   267,     0,     0,
     571,   572,   573,   569,   575,   524,   527,   526,     0,     0,
     525,   528,   529,   562,     0,   563,   459,     0,   576,   534,
     551,   467,     0,   401,     0,   397,     0,     0,     0,   505,
     218,   217,   413,   408,     0,     0,   407,   402,     0,     0,
       0,   564,   515,   556,   557,   530,   531,   536,   539,   537,
     544,   545,   535,   541,   540,     0,   543,   311,   308,     0,
     257,     0,     0,   296,   303,   297,   302,   299,   304,   298,
     300,     0,     0,     0,   278,     0,     0,   301,     0,     0,
     301,   264,     0,     0,     0,   120,     0,     0,   137,   146,
     147,     0,   151,   132,   131,     0,   133,     0,   130,   134,
     135,     0,   140,   138,   558,   559,   415,   426,   428,   429,
     430,   427,     0,   420,   424,     0,     0,     0,     0,     0,
       0,   116,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    86,   100,    99,    98,    97,    96,    95,
      91,    90,    92,    93,    94,     0,     0,   384,     0,     0,
     511,   519,   504,   510,   516,   517,   508,   518,   523,   509,
     506,   522,   388,     0,    77,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   152,   165,
     164,   163,   162,   161,   160,   157,   156,   158,   159,   166,
     433,   354,   352,     0,   369,     0,     0,     0,     0,   206,
     207,     0,     0,   223,   222,   213,   212,   231,   228,     0,
     568,   438,   435,     0,    59,    65,    71,   577,   578,   579,
     580,   581,   582,   583,   584,   585,   586,   587,   588,   589,
     590,   591,   592,   593,   594,   595,   596,   597,   598,   599,
     600,   601,   602,   603,   472,   473,   266,   265,   605,   607,
     609,   608,     0,   462,   470,     0,     0,   399,   398,     0,
     409,     0,   403,     0,   124,     0,   373,     0,   310,   259,
     274,   306,   305,   271,   301,   301,     0,   301,     0,     0,
     289,     0,   263,   262,     0,     0,     0,     0,     0,     0,
     142,     0,     0,     0,     0,   414,   301,   425,     0,     0,
       0,     0,     0,     0,     0,     0,   115,   103,   104,   105,
     106,   107,   108,   109,   110,   111,   112,     0,     0,     0,
     382,   390,     0,     0,   183,   169,   170,   171,   172,   173,
     174,   175,   176,   177,   178,   353,   370,   208,   200,   197,
     203,   204,     0,     0,   230,   437,     0,   604,   397,     0,
     394,   411,     0,   405,     0,     0,     0,   538,   275,     0,
       0,     0,   301,     0,   301,   301,   287,     0,   119,     0,
     150,   520,   129,     0,     0,     0,     0,   421,     0,     0,
     190,     0,   196,     0,     0,     0,     0,   117,   379,   385,
       0,     0,     0,     0,   205,   606,     0,     0,   412,   406,
       0,   374,   301,   301,   301,     0,   295,     0,     0,     0,
     181,     0,   145,   141,   139,   301,   422,     0,     0,     0,
       0,     0,   189,   113,   114,   383,   179,   180,   202,     0,
     395,   301,   280,   276,   279,   301,   293,   288,   128,     0,
       0,   192,   191,   195,   193,     0,     0,     0,   372,   301,
       0,     0,   143,   423,     0,   199,     0,   396,     0,   281,
       0,   294,   194,     0,   188,   371,   301,   301,   198,   282,
     277,   301,   283
  };

  /* YYPGOTO[NTERM-NUM].  */
  const short int
  parser::yypgoto_[] =
  {
     -1090, -1090,  1272, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,  -330, -1090,
   -1090, -1090, -1090,  -110,  -226, -1090, -1090,   997, -1090,   249,
   -1090, -1090, -1090, -1090, -1090, -1090,  -970,  -618,  -131,  -565,
   -1090, -1090, -1090,  1182,  -257, -1090, -1090, -1090, -1090,   346,
   -1090, -1090,   590, -1090, -1090,   751, -1090, -1090,   594, -1090,
   -1090,  -116,   -24,   640,   783, -1090, -1090,  1019, -1090, -1090,
   -1089, -1090, -1090,  1023, -1090, -1090,  1026,  -987,  -587, -1090,
   -1090,   731, -1090,  1203,   605, -1090,   205, -1090, -1090, -1090,
   -1090,   980, -1090, -1090, -1090, -1090, -1090, -1090, -1090,  1134,
    -747, -1090, -1090, -1090, -1090, -1090,   713, -1090,   272,  -830,
   -1090, -1090, -1090, -1090, -1090,   620, -1090,   -46,   792, -1090,
   -1090,   793, -1090, -1090,   581, -1090,  -519, -1090,   -94, -1090,
    1235, -1090, -1090, -1090, -1090, -1090, -1090, -1090,  -101, -1090,
   -1090,  -125,  -585, -1090, -1090, -1090, -1090,  -103,   -87,   -79,
     -77,   -76, -1090, -1090,   -92,   -53, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090,   -72, -1090, -1090, -1090, -1090, -1090,
     -68,   -67,   -48,   -61,   -57,   -54, -1090, -1090, -1090, -1090,
    -111,   -95,   -85,   -82,   -52,   -50,   -44, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090, -1090,
   -1090, -1090, -1090, -1090, -1090,   573, -1090,  -524
  };

  /* YYDEFGOTO[NTERM-NUM].  */
  const short int
  parser::yydefgoto_[] =
  {
        -1,    44,    45,    46,    47,    48,    49,    50,    51,    52,
     152,   154,   156,   131,    53,    54,    55,    56,   362,   902,
      57,   324,    58,   228,   229,    59,   320,   321,   877,   878,
      60,   327,  1073,  1072,  1153,   881,   627,   628,   629,   630,
     436,    61,    62,   342,   343,  1163,    63,  1236,   728,   729,
      64,   460,   461,    65,   214,   215,    66,   456,   457,    67,
     463,   467,   110,   864,   780,    68,   306,   307,   308,   852,
    1138,    69,   317,   318,    70,   312,   313,   853,  1139,    71,
     259,   260,    72,   437,   438,    73,  1046,  1047,    74,    75,
     364,   365,    76,    77,   367,    78,    79,    80,   211,   212,
     566,    81,    82,    83,    84,   335,   336,   892,   893,   894,
      85,   134,   720,    86,   468,   469,   179,   180,   181,    87,
     203,   204,    88,    89,   513,   514,   776,   386,   387,   388,
     389,   390,   391,   392,   393,   394,   395,   396,   397,   398,
     399,   400,   401,   880,   402,   403,   404,   182,   183,   184,
     185,   186,   268,   269,   405,   441,   272,   273,   274,   275,
     276,   277,   278,   279,   442,   281,   282,   283,   284,   285,
     443,   444,   445,   446,   447,   448,   406,   292,   293,   337,
     407,   408,   187,   188,   451,   189,   190,   299,   470,   191,
     192,   193,   194,   195,   196,   197,   207,   515,   516,   517,
     518,   519,   520,   521,   522,   523,   524,   525,   526,   527,
     528,   529,   530,   531,   532,   533,   534,   535,   536,   537,
     538,   539,   540,   541,   795,  1029,   789,   790
  };

  /* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule which
     number is the opposite.  If zero, do what YYDEFACT says.  */
  const signed char parser::yytable_ninf_ = -1;
  const short int
  parser::yytable_[] =
  {
       128,   129,   582,   435,   216,   322,   263,   137,   262,   869,
     338,   385,   146,   149,   150,   261,   458,   270,   157,   645,
     646,   854,   264,   856,   294,   781,   339,   295,   859,   464,
     265,   657,   266,   267,   439,   439,   674,   280,   459,   799,
     205,   286,   287,   440,   440,   206,   879,   202,   289,   462,
     449,   449,   290,   450,   450,   291,   271,   296,  1036,   297,
     821,   288,   870,  1077,   868,   298,  1028,   979,  1140,    90,
     896,   827,   828,   829,   725,   417,   980,   787,   219,    92,
     836,   300,  1120,   726,   418,   641,   107,  1192,   300,   564,
     844,  1121,    96,   171,  1213,   101,   582,   419,   417,   846,
     581,  1154,  1155,  1156,   209,   420,   107,   418,   570,   620,
     843,   102,   622,   844,   332,   421,   104,   600,   631,   636,
     419,   575,   846,   721,   221,    99,   333,   576,   420,   132,
     848,   107,   222,   106,   100,   111,   107,   721,   421,   334,
     107,   887,   121,   301,   107,   565,   107,   133,   845,   671,
     301,   227,   107,   848,   887,   123,   847,   107,  1250,   931,
     319,   112,   340,   227,   107,   571,   938,   107,   378,   981,
     107,   640,   887,   107,   601,   632,   637,   340,   422,   423,
     722,   897,   210,   107,   424,   425,   426,   427,   428,   429,
     430,   431,   432,   851,   723,   113,   849,   641,   302,   433,
     363,   422,   423,   108,   109,   409,   888,   424,   425,   426,
     427,   428,   429,   430,   431,   432,   851,   990,    91,   888,
    1030,  1219,   433,   126,   127,   982,   220,  1193,    93,   103,
     788,   850,   727,   434,   105,   626,   898,   888,   417,  1015,
    1122,  1012,  1013,  1240,  1194,  1016,  1017,   418,   813,  1020,
     820,  1022,  1023,   817,  1227,   417,   434,   838,   626,  1076,
     419,   942,   889,   972,   418,   675,   890,   891,   420,   974,
    1058,   341,  1195,  1061,   988,   889,   871,   419,   421,   890,
     891,   144,   145,   300,   992,   420,   341,   147,   148,  1081,
    1131,   696,   697,   889,   412,   421,   117,   890,   891,   300,
    1133,   300,   224,   708,   114,   118,    94,    95,   122,  1082,
     225,   300,   899,   900,   901,   300,   676,   903,   904,   905,
     906,   907,   908,   909,   910,   911,   912,  1077,   914,   915,
     916,   917,   918,   919,   920,   921,   922,   923,   924,   474,
     925,   422,   423,   478,   482,   301,   929,   424,   425,   426,
     427,   428,   429,   430,   431,   432,   413,   300,   422,   423,
     300,   301,   433,   301,   424,   425,   426,   427,   428,   429,
     430,   431,   432,   301,   300,   300,   300,   301,   124,   433,
     709,  1176,   710,   711,   712,   713,   714,   300,   715,   716,
     717,   718,   300,   719,   328,   843,   434,   976,   626,   777,
     410,   475,   300,   300,  1074,   479,   483,   303,   340,   125,
     300,   414,   209,   434,  1200,   626,   453,   724,   465,   301,
     107,   884,   301,   300,  1243,   730,   732,   130,   471,   734,
     737,   844,   472,   845,   138,   861,   301,   301,   301,   567,
     846,   847,   872,  1069,  1148,   885,   774,   135,   806,   301,
     136,   139,   578,   329,   301,   775,   476,   807,   579,  1075,
     480,   484,   309,   330,   301,   301,   303,   304,  1141,  1201,
    1143,   848,   301,   151,   485,   606,   216,  1040,  1039,  1244,
     153,   849,   731,   733,   879,   301,   735,   738,  1085,  1158,
     210,  1042,  1048,  1101,   155,   227,  1136,   263,   208,   262,
    1070,   309,   778,   779,  1115,   305,   261,   205,   270,  1116,
     612,   162,   206,   264,   202,   294,   850,   341,   295,  1124,
    1125,   265,   310,   266,   267,   338,   304,  1178,   280,   862,
     863,   198,   286,   287,   851,  1051,   869,   869,   869,   289,
    1179,   339,   213,   290,  1052,  1086,   291,   271,   296,  1151,
     297,   814,   288,  1137,   818,  1185,   298,  1187,  1188,   314,
     311,   310,   945,   946,   305,   217,  1079,   948,   949,   950,
     951,   952,   953,   954,   955,   956,   957,   839,   959,   960,
     961,   962,   963,   964,   965,   966,   967,   968,   969,   870,
     870,   870,    97,    98,   977,  1212,  1098,  1214,   314,   311,
     978,   869,  1232,   666,   667,   458,   668,   617,  1220,   717,
     718,   369,   719,   417,   115,   116,   119,   120,   987,   315,
     439,   218,   418,   223,  1228,   140,   141,   459,  1231,   440,
     142,   143,   873,   158,   159,   419,   449,   226,   462,   450,
     227,   230,  1239,   420,   874,   319,   664,   665,   666,   667,
     875,   668,  1123,   421,   870,   323,   325,   316,   315,  1249,
     326,   363,   366,   415,  1252,   943,   411,   416,   455,   542,
     876,   163,   164,   165,   166,   167,   168,   169,   199,   473,
     477,   543,   200,   170,   544,   481,   545,   171,   715,   716,
     717,   718,   546,   719,   547,   555,   316,   548,   973,   975,
     549,   550,   551,   172,   552,   201,   553,   554,   556,   557,
     558,   989,   559,   560,   993,   561,   422,   423,   562,   563,
     569,   572,   424,   425,   426,   427,   428,   429,   430,   431,
     432,   573,   574,   580,   577,   584,   583,   433,   585,   586,
    1067,   587,   588,   589,   590,   591,   173,   174,   592,   593,
     594,   595,  1164,  1165,  1166,  1167,  1065,   596,   597,   160,
     598,   599,   602,   603,   175,   176,     1,     2,   604,  1170,
     605,   434,   608,   626,   609,   610,     3,     4,     5,   611,
     614,   616,   615,     6,   619,   231,   623,     7,     8,     9,
     624,    10,  1173,    11,    12,    13,    14,   177,   178,  1177,
     200,   625,   633,   635,   634,   638,    15,   642,   639,    16,
     643,   644,   677,   647,   678,   648,   649,  1189,   650,   232,
     233,   651,    17,   201,   652,   653,   234,   654,   679,  1197,
     655,  1198,   656,   235,   670,    18,    19,    20,   680,   582,
     658,    21,   659,   660,   661,   662,   663,   681,   664,   665,
     666,   667,    22,   668,    23,   682,    24,    25,    26,    27,
      28,   252,   683,   684,   685,    29,    30,   686,   687,   254,
      31,    32,    33,    34,   688,   690,   689,   692,   691,  1226,
      35,    36,   693,    37,   694,   256,   736,    38,   739,   695,
      39,    40,    41,    42,   698,  1235,   699,   257,   740,   741,
     742,   700,   743,   258,   913,   744,   745,   701,   746,   747,
     748,   749,   750,  1248,   702,   177,   178,   751,   752,   703,
     704,     1,     2,   753,   754,   755,    43,   756,   757,   758,
     705,     3,     4,     5,   706,  1149,   707,   759,     6,   760,
     761,   762,     7,     8,     9,   763,    10,   772,    11,    12,
      13,    14,   764,   417,   765,   766,   767,   768,   769,   770,
     771,    15,   418,   773,    16,   785,   782,   784,   786,   791,
     793,  1171,  1172,   803,   808,   419,   344,    17,   804,   792,
     794,   797,   798,   420,   800,   345,   802,   805,   809,   811,
      18,    19,    20,   421,   822,   815,    21,   812,   346,   816,
     823,   788,   819,   824,   825,   826,   347,    22,   830,    23,
     831,    24,    25,    26,    27,    28,   348,  1132,   832,  1134,
      29,    30,   833,   841,   834,    31,    32,    33,    34,   840,
     842,   855,   835,   857,   930,    35,    36,   858,    37,   866,
     865,   860,    38,   867,   882,    39,    40,    41,    42,   883,
     668,   926,   932,   895,   933,   934,   422,   423,   935,   936,
     937,   939,   424,   425,   426,   427,   428,   429,   430,   431,
     432,   944,   940,   941,   344,   997,   994,   433,   673,   349,
     350,    43,   719,   345,   998,   351,   352,   353,   354,   355,
     356,   357,   358,   359,   995,   999,   346,   344,   970,   996,
     360,   984,   986,  1000,   347,  1001,   345,  1002,  1003,  1004,
    1005,   434,  1006,   626,   348,  1007,  1008,  1009,  1010,   346,
    1011,  1014,  1018,  1019,   417,  1021,  1026,   347,  1024,  1027,
    1028,  1035,  1049,   418,   361,  1033,  1054,   348,  1034,   709,
     565,   710,   711,   712,   713,   714,   419,   715,   716,   717,
     718,  1032,   719,  1062,   420,  1050,  1037,   658,  1055,   659,
     660,   661,   662,   663,   421,   664,   665,   666,   667,  1045,
     668,  1038,  1056,  1041,  1053,  1043,   928,   349,   350,  1057,
    1059,  1060,  1063,   351,   352,   353,   354,   355,   356,   357,
     358,   359,  1064,  1066,  1068,  1071,  1078,  1080,   360,    -1,
     349,   350,  1099,   958,  1119,  1127,   351,   352,   353,   354,
     355,   356,   357,   358,   359,  1126,  1104,  1128,  1129,   669,
    1135,   360,  1142,  1147,  1146,  1144,  1159,   422,   423,  1145,
     231,  1161,   361,   424,   425,   426,   427,   428,   429,   430,
     431,   432,  1162,  1175,  1168,   200,   170,  1180,   433,  1182,
     171,  1183,  1184,  1209,  1196,   361,  1211,   163,   164,   165,
     166,   167,   168,   169,   232,   233,   172,  1186,   201,   170,
    1215,   234,  1223,   171,  1191,  1216,  1217,  1224,   235,   236,
     237,  1225,   434,   238,   239,  1229,   240,   241,  1230,   172,
     242,   243,   244,   245,   246,   247,   248,  1199,   249,   250,
     251,  1233,  1234,  1241,  1238,  1242,   252,  1245,  1246,   173,
     174,  1247,   253,   368,   254,  1251,   161,   621,  1152,   255,
     454,   810,   985,  1118,   983,   607,   971,   175,   176,   783,
     256,   837,   173,   174,   369,   947,   370,   371,   613,   452,
     618,  1181,   257,   213,   672,   568,  1157,   796,   258,   886,
     175,   176,     0,   801,  1025,   331,   234,   991,   372,   373,
     177,   178,  1031,   235,     0,     0,   368,     0,     0,   658,
     328,   659,   660,   661,   662,   663,     0,   664,   665,   666,
     667,     0,   668,   177,   178,     0,     0,   369,     0,   370,
     371,     0,     0,     0,     0,     0,   374,     0,   375,   254,
     376,   333,     0,     0,     0,     0,   377,     0,     0,   234,
     378,   372,   373,     0,   334,     0,   235,     0,   379,   380,
     381,     0,     0,   328,   382,   383,   384,     0,   213,     0,
       0,     0,     0,  1087,   658,   466,   659,   660,   661,   662,
     663,     0,   664,   665,   666,   667,     0,   668,     0,   374,
       0,   375,   254,   376,   333,     0,     0,     0,     0,   377,
       0,     0,     0,   378,     0,     0,     0,   334,     0,     0,
       0,   379,   380,   381,     0,     0,     0,   382,   383,   384,
     658,   213,   659,   660,   661,   662,   663,     0,   664,   665,
     666,   667,     0,   668,     0,     0,     0,   658,  1088,   659,
     660,   661,   662,   663,     0,   664,   665,   666,   667,   658,
     668,   659,   660,   661,   662,   663,     0,   664,   665,   666,
     667,   658,   668,   659,   660,   661,   662,   663,     0,   664,
     665,   666,   667,     0,   668,     0,     0,   658,     0,   659,
     660,   661,   662,   663,  1089,   664,   665,   666,   667,   658,
     668,   659,   660,   661,   662,   663,     0,   664,   665,   666,
     667,  1090,   668,     0,     0,     0,   658,     0,   659,   660,
     661,   662,   663,  1091,   664,   665,   666,   667,   658,   668,
     659,   660,   661,   662,   663,  1092,   664,   665,   666,   667,
       0,   668,     0,     0,   658,     0,   659,   660,   661,   662,
     663,  1093,   664,   665,   666,   667,   709,   668,   710,   711,
     712,   713,   714,  1094,   715,   716,   717,   718,   709,   719,
     710,   711,   712,   713,   714,     0,   715,   716,   717,   718,
    1095,   719,     0,     0,     0,   709,     0,   710,   711,   712,
     713,   714,  1096,   715,   716,   717,   718,     0,   719,     0,
       0,   709,     0,   710,   711,   712,   713,   714,  1097,   715,
     716,   717,   718,   709,   719,   710,   711,   712,   713,   714,
    1105,   715,   716,   717,   718,   709,   719,   710,   711,   712,
     713,   714,  1106,   715,   716,   717,   718,   709,   719,   710,
     711,   712,   713,   714,     0,   715,   716,   717,   718,  1107,
     719,   709,     0,   710,   711,   712,   713,   714,     0,   715,
     716,   717,   718,     0,   719,  1108,     0,     0,     0,     0,
     709,     0,   710,   711,   712,   713,   714,  1109,   715,   716,
     717,   718,   709,   719,   710,   711,   712,   713,   714,  1110,
     715,   716,   717,   718,   658,   719,   659,   660,   661,   662,
     663,  1111,   664,   665,   666,   667,     0,   668,   658,     0,
     659,   660,   661,   662,   663,  1112,   664,   665,   666,   667,
     658,   668,   659,   660,   661,   662,   663,     0,   664,   665,
     666,   667,     0,   668,  1113,     0,     0,     0,     0,   658,
       0,   659,   660,   661,   662,   663,  1114,   664,   665,   666,
     667,   709,   668,   710,   711,   712,   713,   714,  1117,   715,
     716,   717,   718,     0,   719,   709,     0,   710,   711,   712,
     713,   714,  1174,   715,   716,   717,   718,   658,   719,   659,
     660,   661,   662,   663,  1203,   664,   665,   666,   667,   658,
     668,   659,   660,   661,   662,   663,     0,   664,   665,   666,
     667,     0,   668,  1204,     0,     0,     0,     0,   658,     0,
     659,   660,   661,   662,   663,  1206,   664,   665,   666,   667,
       0,   668,   658,     0,   659,   660,   661,   662,   663,  1207,
     664,   665,   666,   667,   658,   668,   659,   660,   661,   662,
     663,  1208,   664,   665,   666,   667,     0,   668,     0,     0,
     709,   927,   710,   711,   712,   713,   714,     0,   715,   716,
     717,   718,     0,   719,     0,     0,     0,     0,     0,   658,
    1044,   659,   660,   661,   662,   663,     0,   664,   665,   666,
     667,     0,   668,     0,  1100,     0,     0,     0,     0,   658,
       0,   659,   660,   661,   662,   663,  1130,   664,   665,   666,
     667,     0,   668,     0,     0,   709,     0,   710,   711,   712,
     713,   714,  1150,   715,   716,   717,   718,   658,   719,   659,
     660,   661,   662,   663,     0,   664,   665,   666,   667,     0,
     668,  1160,   658,     0,   659,   660,   661,   662,   663,     0,
     664,   665,   666,   667,   658,   668,   659,   660,   661,   662,
     663,  1169,   664,   665,   666,   667,     0,   668,     0,     0,
     658,     0,   659,   660,   661,   662,   663,  1190,   664,   665,
     666,   667,   658,   668,   659,   660,   661,   662,   663,  1202,
     664,   665,   666,   667,     0,   668,     0,   658,     0,   659,
     660,   661,   662,   663,  1205,   664,   665,   666,   667,   658,
     668,   659,   660,   661,   662,   663,  1210,   664,   665,   666,
     667,     0,   668,     0,     0,   709,     0,   710,   711,   712,
     713,   714,  1218,   715,   716,   717,   718,     0,   719,     0,
       0,     0,     0,   658,  1221,   659,   660,   661,   662,   663,
    1083,   664,   665,   666,   667,     0,   668,     0,     0,  1222,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1237,   486,   487,   488,   489,   490,   491,   492,   493,
     494,   495,   496,   497,   498,   499,   500,   501,   502,   503,
     504,   505,   506,     0,     0,     0,   507,   508,     0,   509,
     510,   511,   512,   658,     0,   659,   660,   661,   662,   663,
    1084,   664,   665,   666,   667,   709,   668,   710,   711,   712,
     713,   714,  1102,   715,   716,   717,   718,   709,   719,   710,
     711,   712,   713,   714,  1103,   715,   716,   717,   718,   658,
     719,   659,   660,   661,   662,   663,     0,   664,   665,   666,
     667,     0,   668
  };

  /* YYCHECK.  */
  const short int
  parser::yycheck_[] =
  {
        24,    25,   228,   134,    98,   115,   109,    31,   109,   627,
     121,   127,    36,    37,    38,   109,   141,   109,    42,   349,
     350,   608,   109,   610,   109,   544,   121,   109,   615,   145,
     109,   361,   109,   109,   135,   136,   366,   109,   141,   558,
      93,   109,   109,   135,   136,    93,   631,    93,   109,   143,
     135,   136,   109,   135,   136,   109,   109,   109,   805,   109,
     584,   109,   627,   893,    31,   109,    51,    42,  1055,    51,
      33,   590,   591,   592,    42,    42,    51,    51,    51,    51,
     599,    82,    42,    51,    51,   342,    82,    31,    82,   144,
      42,    51,   199,    25,  1183,    51,   322,    64,    42,    51,
      31,  1071,  1072,  1073,     4,    72,    82,    51,   144,    31,
       6,    82,    31,    42,    79,    82,    82,   144,   144,   144,
      64,   201,    51,   144,   199,    33,    91,   207,    72,    33,
      82,    82,   207,   199,    42,   199,    82,   144,    82,   104,
      82,    82,   200,   144,    82,   200,    82,    51,    44,    31,
     144,    82,    82,    82,    82,    33,    52,    82,  1247,   678,
      82,   199,    22,    82,    82,   201,   685,    82,   100,   144,
      82,    31,    82,    82,   201,   201,   201,    22,   145,   146,
     201,   144,    82,    82,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   145,   201,   199,    92,   454,   199,   166,
      82,   145,   146,   199,   200,   199,   147,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   145,   736,   200,   147,
     205,  1191,   166,   199,   200,   200,   199,    31,   200,   200,
     204,   127,   200,   200,   200,   202,   199,   147,    42,   763,
     200,   760,   761,  1230,    31,   764,   765,    51,   199,   768,
     580,   770,   771,   199,   206,    42,   200,   199,   202,   200,
      64,   199,   203,   199,    51,    31,   207,   208,    72,   199,
     857,   131,   200,   860,   199,   203,    33,    64,    82,   207,
     208,   199,   200,    82,   199,    72,   131,   199,   200,   100,
     199,   422,   423,   203,    82,    82,    33,   207,   208,    82,
     199,    82,   199,   434,   199,    42,   199,   200,   199,   120,
     207,    82,   642,   643,   644,    82,    82,   647,   648,   649,
     650,   651,   652,   653,   654,   655,   656,  1157,   658,   659,
     660,   661,   662,   663,   664,   665,   666,   667,   668,    82,
     670,   145,   146,    82,    82,   144,   676,   151,   152,   153,
     154,   155,   156,   157,   158,   159,   144,    82,   145,   146,
      82,   144,   166,   144,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   144,    82,    82,    82,   144,   199,   166,
     137,  1128,   139,   140,   141,   142,   143,    82,   145,   146,
     147,   148,    82,   150,    60,     6,   200,   727,   202,    42,
     199,   144,    82,    82,   144,   144,   144,    22,    22,   199,
      82,   199,     4,   200,   144,   202,   199,    31,   199,   144,
      82,    56,   144,    82,   144,   144,   144,    82,   199,   144,
     144,    42,   199,    44,   199,    51,   144,   144,   144,    31,
      51,    52,   199,   144,    31,    80,    42,   200,    42,   144,
     200,    33,   201,   119,   144,    51,   199,    51,   207,   199,
     199,   199,    22,   129,   144,   144,    22,    82,  1055,   199,
    1057,    82,   144,    82,   199,    31,   570,   199,   808,   199,
      82,    92,   201,   201,  1069,   144,   201,   201,   144,  1076,
      82,   199,   199,   199,    82,    82,   144,   600,    23,   600,
     201,    22,   145,   146,   199,   120,   600,   560,   600,   199,
      31,   199,   560,   600,   560,   600,   127,   131,   600,   199,
     199,   600,    82,   600,   600,   636,    82,   199,   600,   145,
     146,   199,   600,   600,   145,    42,  1154,  1155,  1156,   600,
     199,   636,   118,   600,    51,   201,   600,   600,   600,  1068,
     600,   575,   600,   201,   578,  1142,   600,  1144,  1145,    22,
     120,    82,   693,   694,   120,    42,   896,   698,   699,   700,
     701,   702,   703,   704,   705,   706,   707,   601,   709,   710,
     711,   712,   713,   714,   715,   716,   717,   718,   719,  1154,
    1155,  1156,   199,   200,   199,  1182,   926,  1184,    22,   120,
     205,  1219,    31,   147,   148,   730,   150,    31,  1195,   147,
     148,    24,   150,    42,   199,   200,   199,   200,   734,    82,
     721,   199,    51,    82,  1211,   199,   200,   730,  1215,   721,
     199,   200,    45,   199,   200,    64,   721,    82,   732,   721,
      82,    33,  1229,    72,    57,    82,   145,   146,   147,   148,
      63,   150,   982,    82,  1219,    38,    42,   120,    82,  1246,
     199,    82,    82,    51,  1251,   689,   126,   199,   204,    33,
      83,     9,    10,    11,    12,    13,    14,    15,    16,   126,
     126,    33,    20,    21,    33,   126,    33,    25,   145,   146,
     147,   148,    33,   150,    33,   144,   120,    33,   722,   723,
      33,    33,    33,    41,    33,    43,    33,    33,   201,    33,
      33,   735,    33,   144,   738,   201,   145,   146,   201,    82,
      33,   199,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   199,    82,    33,    82,    33,    82,   166,    33,    33,
     871,    33,    33,    33,    33,    33,    84,    85,    33,    33,
      33,    33,  1082,  1083,  1084,  1085,   866,    33,    33,     0,
      33,    33,    33,    82,   102,   103,     7,     8,    82,  1099,
      82,   200,   144,   202,    82,   144,    17,    18,    19,    82,
      82,    82,   144,    24,   200,     5,    33,    28,    29,    30,
     201,    32,  1122,    34,    35,    36,    37,   135,   136,  1129,
      20,   199,   201,    33,   201,    82,    47,   200,    82,    50,
     200,   200,    33,   200,    33,   200,   200,  1147,   200,    39,
      40,   200,    63,    43,   200,   200,    46,   200,    33,  1159,
     200,  1161,   200,    53,   200,    76,    77,    78,    33,  1065,
     137,    82,   139,   140,   141,   142,   143,    33,   145,   146,
     147,   148,    93,   150,    95,    33,    97,    98,    99,   100,
     101,    81,    33,    33,    33,   106,   107,    33,    33,    89,
     111,   112,   113,   114,    33,   126,   201,   199,    82,  1209,
     121,   122,   200,   124,   200,   105,    33,   128,   126,   200,
     131,   132,   133,   134,   200,  1225,   200,   117,    82,   126,
      82,   200,   126,   123,   201,    82,    33,   200,    33,    33,
      33,    33,    33,  1243,   200,   135,   136,    33,    33,   200,
     200,     7,     8,    33,    33,    33,   167,    33,    33,    33,
     200,    17,    18,    19,   200,  1066,   200,    33,    24,    33,
      33,    33,    28,    29,    30,    33,    32,   201,    34,    35,
      36,    37,    33,    42,    33,    33,    33,    33,    33,    33,
      33,    47,    51,   144,    50,    82,    51,    51,    51,    51,
      82,  1102,  1103,   199,    33,    64,    42,    63,   200,    51,
      82,    51,    51,    72,    51,    51,    51,    82,    51,   199,
      76,    77,    78,    82,    51,    82,    82,   199,    64,   199,
      51,   204,    82,    51,    51,    51,    72,    93,    51,    95,
      51,    97,    98,    99,   100,   101,    82,  1041,    51,  1043,
     106,   107,    51,    33,    82,   111,   112,   113,   114,    82,
     144,   144,   200,   144,    51,   121,   122,   144,   124,   199,
      82,   144,   128,    82,   199,   131,   132,   133,   134,   199,
     150,    82,    51,   144,    51,    51,   145,   146,    51,    51,
      51,    51,   151,   152,   153,   154,   155,   156,   157,   158,
     159,   126,    51,    51,    42,    51,   126,   166,   144,   145,
     146,   167,   150,    51,    51,   151,   152,   153,   154,   155,
     156,   157,   158,   159,   126,    51,    64,    42,   199,   126,
     166,   199,   199,    51,    72,    51,    51,    51,    51,    51,
      51,   200,    51,   202,    82,    51,    51,    51,    51,    64,
      51,    51,    51,    51,    42,    51,    42,    72,   199,    42,
      51,    51,    82,    51,   200,   199,   144,    82,   199,   137,
     200,   139,   140,   141,   142,   143,    64,   145,   146,   147,
     148,   203,   150,    51,    72,    82,   201,   137,   144,   139,
     140,   141,   142,   143,    82,   145,   146,   147,   148,   206,
     150,   201,    82,   201,   199,   201,   144,   145,   146,   144,
      82,   199,    51,   151,   152,   153,   154,   155,   156,   157,
     158,   159,   201,    33,    33,   199,    82,    82,   166,   150,
     145,   146,    82,   201,   199,    82,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   203,   201,   201,    33,   199,
      82,   166,   144,    33,   199,   144,    33,   145,   146,   144,
       5,    33,   200,   151,   152,   153,   154,   155,   156,   157,
     158,   159,    51,    51,   199,    20,    21,   206,   166,   144,
      25,   144,   144,    33,   201,   200,   144,     9,    10,    11,
      12,    13,    14,    15,    39,    40,    41,   199,    43,    21,
     144,    46,    51,    25,   199,   199,   199,    51,    53,    54,
      55,   130,   200,    58,    59,   144,    61,    62,   144,    41,
      65,    66,    67,    68,    69,    70,    71,   203,    73,    74,
      75,   201,   203,   199,    82,    51,    81,   206,   144,    84,
      85,   144,    87,     3,    89,   144,    44,   320,  1069,    94,
     138,   570,   732,   977,   730,   306,   721,   102,   103,   546,
     105,   600,    84,    85,    24,   695,    26,    27,   312,   136,
     317,  1136,   117,   118,   364,   211,  1074,   555,   123,   636,
     102,   103,    -1,   560,   773,   120,    46,   737,    48,    49,
     135,   136,   789,    53,    -1,    -1,     3,    -1,    -1,   137,
      60,   139,   140,   141,   142,   143,    -1,   145,   146,   147,
     148,    -1,   150,   135,   136,    -1,    -1,    24,    -1,    26,
      27,    -1,    -1,    -1,    -1,    -1,    86,    -1,    88,    89,
      90,    91,    -1,    -1,    -1,    -1,    96,    -1,    -1,    46,
     100,    48,    49,    -1,   104,    -1,    53,    -1,   108,   109,
     110,    -1,    -1,    60,   114,   115,   116,    -1,   118,    -1,
      -1,    -1,    -1,   201,   137,   125,   139,   140,   141,   142,
     143,    -1,   145,   146,   147,   148,    -1,   150,    -1,    86,
      -1,    88,    89,    90,    91,    -1,    -1,    -1,    -1,    96,
      -1,    -1,    -1,   100,    -1,    -1,    -1,   104,    -1,    -1,
      -1,   108,   109,   110,    -1,    -1,    -1,   114,   115,   116,
     137,   118,   139,   140,   141,   142,   143,    -1,   145,   146,
     147,   148,    -1,   150,    -1,    -1,    -1,   137,   201,   139,
     140,   141,   142,   143,    -1,   145,   146,   147,   148,   137,
     150,   139,   140,   141,   142,   143,    -1,   145,   146,   147,
     148,   137,   150,   139,   140,   141,   142,   143,    -1,   145,
     146,   147,   148,    -1,   150,    -1,    -1,   137,    -1,   139,
     140,   141,   142,   143,   201,   145,   146,   147,   148,   137,
     150,   139,   140,   141,   142,   143,    -1,   145,   146,   147,
     148,   201,   150,    -1,    -1,    -1,   137,    -1,   139,   140,
     141,   142,   143,   201,   145,   146,   147,   148,   137,   150,
     139,   140,   141,   142,   143,   201,   145,   146,   147,   148,
      -1,   150,    -1,    -1,   137,    -1,   139,   140,   141,   142,
     143,   201,   145,   146,   147,   148,   137,   150,   139,   140,
     141,   142,   143,   201,   145,   146,   147,   148,   137,   150,
     139,   140,   141,   142,   143,    -1,   145,   146,   147,   148,
     201,   150,    -1,    -1,    -1,   137,    -1,   139,   140,   141,
     142,   143,   201,   145,   146,   147,   148,    -1,   150,    -1,
      -1,   137,    -1,   139,   140,   141,   142,   143,   201,   145,
     146,   147,   148,   137,   150,   139,   140,   141,   142,   143,
     201,   145,   146,   147,   148,   137,   150,   139,   140,   141,
     142,   143,   201,   145,   146,   147,   148,   137,   150,   139,
     140,   141,   142,   143,    -1,   145,   146,   147,   148,   201,
     150,   137,    -1,   139,   140,   141,   142,   143,    -1,   145,
     146,   147,   148,    -1,   150,   201,    -1,    -1,    -1,    -1,
     137,    -1,   139,   140,   141,   142,   143,   201,   145,   146,
     147,   148,   137,   150,   139,   140,   141,   142,   143,   201,
     145,   146,   147,   148,   137,   150,   139,   140,   141,   142,
     143,   201,   145,   146,   147,   148,    -1,   150,   137,    -1,
     139,   140,   141,   142,   143,   201,   145,   146,   147,   148,
     137,   150,   139,   140,   141,   142,   143,    -1,   145,   146,
     147,   148,    -1,   150,   201,    -1,    -1,    -1,    -1,   137,
      -1,   139,   140,   141,   142,   143,   201,   145,   146,   147,
     148,   137,   150,   139,   140,   141,   142,   143,   201,   145,
     146,   147,   148,    -1,   150,   137,    -1,   139,   140,   141,
     142,   143,   201,   145,   146,   147,   148,   137,   150,   139,
     140,   141,   142,   143,   201,   145,   146,   147,   148,   137,
     150,   139,   140,   141,   142,   143,    -1,   145,   146,   147,
     148,    -1,   150,   201,    -1,    -1,    -1,    -1,   137,    -1,
     139,   140,   141,   142,   143,   201,   145,   146,   147,   148,
      -1,   150,   137,    -1,   139,   140,   141,   142,   143,   201,
     145,   146,   147,   148,   137,   150,   139,   140,   141,   142,
     143,   201,   145,   146,   147,   148,    -1,   150,    -1,    -1,
     137,   199,   139,   140,   141,   142,   143,    -1,   145,   146,
     147,   148,    -1,   150,    -1,    -1,    -1,    -1,    -1,   137,
     199,   139,   140,   141,   142,   143,    -1,   145,   146,   147,
     148,    -1,   150,    -1,   199,    -1,    -1,    -1,    -1,   137,
      -1,   139,   140,   141,   142,   143,   199,   145,   146,   147,
     148,    -1,   150,    -1,    -1,   137,    -1,   139,   140,   141,
     142,   143,   199,   145,   146,   147,   148,   137,   150,   139,
     140,   141,   142,   143,    -1,   145,   146,   147,   148,    -1,
     150,   199,   137,    -1,   139,   140,   141,   142,   143,    -1,
     145,   146,   147,   148,   137,   150,   139,   140,   141,   142,
     143,   199,   145,   146,   147,   148,    -1,   150,    -1,    -1,
     137,    -1,   139,   140,   141,   142,   143,   199,   145,   146,
     147,   148,   137,   150,   139,   140,   141,   142,   143,   199,
     145,   146,   147,   148,    -1,   150,    -1,   137,    -1,   139,
     140,   141,   142,   143,   199,   145,   146,   147,   148,   137,
     150,   139,   140,   141,   142,   143,   199,   145,   146,   147,
     148,    -1,   150,    -1,    -1,   137,    -1,   139,   140,   141,
     142,   143,   199,   145,   146,   147,   148,    -1,   150,    -1,
      -1,    -1,    -1,   137,   199,   139,   140,   141,   142,   143,
     144,   145,   146,   147,   148,    -1,   150,    -1,    -1,   199,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   199,   168,   169,   170,   171,   172,   173,   174,   175,
     176,   177,   178,   179,   180,   181,   182,   183,   184,   185,
     186,   187,   188,    -1,    -1,    -1,   192,   193,    -1,   195,
     196,   197,   198,   137,    -1,   139,   140,   141,   142,   143,
     144,   145,   146,   147,   148,   137,   150,   139,   140,   141,
     142,   143,   144,   145,   146,   147,   148,   137,   150,   139,
     140,   141,   142,   143,   144,   145,   146,   147,   148,   137,
     150,   139,   140,   141,   142,   143,    -1,   145,   146,   147,
     148,    -1,   150
  };

  /* STOS_[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
  const unsigned short int
  parser::yystos_[] =
  {
         0,     7,     8,    17,    18,    19,    24,    28,    29,    30,
      32,    34,    35,    36,    37,    47,    50,    63,    76,    77,
      78,    82,    93,    95,    97,    98,    99,   100,   101,   106,
     107,   111,   112,   113,   114,   121,   122,   124,   128,   131,
     132,   133,   134,   167,   210,   211,   212,   213,   214,   215,
     216,   217,   218,   223,   224,   225,   226,   229,   231,   234,
     239,   250,   251,   255,   259,   262,   265,   268,   274,   280,
     283,   288,   291,   294,   297,   298,   301,   302,   304,   305,
     306,   310,   311,   312,   313,   319,   322,   328,   331,   332,
      51,   200,    51,   200,   199,   200,   199,   199,   200,    33,
      42,    51,    82,   200,    82,   200,   199,    82,   199,   200,
     271,   199,   199,   199,   199,   199,   200,    33,    42,   199,
     200,   200,   199,    33,   199,   199,   199,   200,   271,   271,
      82,   222,    33,    51,   320,   200,   200,   271,   199,    33,
     199,   200,   199,   200,   199,   200,   271,   199,   200,   271,
     271,    82,   219,    82,   220,    82,   221,   271,   199,   200,
       0,   211,   199,     9,    10,    11,    12,    13,    14,    15,
      21,    25,    41,    84,    85,   102,   103,   135,   136,   325,
     326,   327,   356,   357,   358,   359,   360,   391,   392,   394,
     395,   398,   399,   400,   401,   402,   403,   404,   199,    16,
      20,    43,   326,   329,   330,   364,   381,   405,    23,     4,
      82,   307,   308,   118,   263,   264,   337,    42,   199,    51,
     199,   199,   207,    82,   199,   207,    82,    82,   232,   233,
      33,     5,    39,    40,    46,    53,    54,    55,    58,    59,
      61,    62,    65,    66,    67,    68,    69,    70,    71,    73,
      74,    75,    81,    87,    89,    94,   105,   117,   123,   289,
     290,   337,   347,   356,   357,   358,   359,   360,   361,   362,
     363,   364,   365,   366,   367,   368,   369,   370,   371,   372,
     373,   374,   375,   376,   377,   378,   379,   380,   381,   382,
     383,   384,   386,   387,   391,   392,   393,   394,   395,   396,
      82,   144,   199,    22,    82,   120,   275,   276,   277,    22,
      82,   120,   284,   285,    22,    82,   120,   281,   282,    82,
     235,   236,   232,    38,   230,    42,   199,   240,    60,   119,
     129,   339,    79,    91,   104,   314,   315,   388,   389,   390,
      22,   131,   252,   253,    42,    51,    64,    72,    82,   145,
     146,   151,   152,   153,   154,   155,   156,   157,   158,   159,
     166,   200,   227,    82,   299,   300,    82,   303,     3,    24,
      26,    27,    48,    49,    86,    88,    90,    96,   100,   108,
     109,   110,   114,   115,   116,   270,   336,   337,   338,   339,
     340,   341,   342,   343,   344,   345,   346,   347,   348,   349,
     350,   351,   353,   354,   355,   363,   385,   389,   390,   199,
     199,   126,    82,   144,   199,    51,   199,    42,    51,    64,
      72,    82,   145,   146,   151,   152,   153,   154,   155,   156,
     157,   158,   159,   166,   200,   247,   249,   292,   293,   347,
     363,   364,   373,   379,   380,   381,   382,   383,   384,   391,
     392,   393,   292,   199,   252,   204,   266,   267,   350,   356,
     260,   261,   337,   269,   270,   199,   125,   270,   323,   324,
     397,   199,   199,   126,    82,   144,   199,   126,    82,   144,
     199,   126,    82,   144,   199,   199,   168,   169,   170,   171,
     172,   173,   174,   175,   176,   177,   178,   179,   180,   181,
     182,   183,   184,   185,   186,   187,   188,   192,   193,   195,
     196,   197,   198,   333,   334,   406,   407,   408,   409,   410,
     411,   412,   413,   414,   415,   416,   417,   418,   419,   420,
     421,   422,   423,   424,   425,   426,   427,   428,   429,   430,
     431,   432,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,   144,   201,    33,    33,    33,
     144,   201,   201,    82,   144,   200,   309,    31,   308,    33,
     144,   201,   199,   199,    82,   201,   207,    82,   201,   207,
      33,    31,   233,    82,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
     144,   201,    33,    82,    82,    82,    31,   276,   144,    82,
     144,    82,    31,   285,    82,   144,    82,    31,   282,   200,
      31,   236,    31,    33,   201,   199,   202,   245,   246,   247,
     248,   144,   201,   201,   201,    33,   144,   201,    82,    82,
      31,   253,   200,   200,   200,   227,   227,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   227,   137,   139,
     140,   141,   142,   143,   145,   146,   147,   148,   150,   199,
     200,    31,   300,   144,   227,    31,    82,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,   201,
     126,    82,   199,   200,   200,   200,   247,   247,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   247,   137,
     139,   140,   141,   142,   143,   145,   146,   147,   148,   150,
     321,   144,   201,   201,    31,    42,    51,   200,   257,   258,
     144,   201,   144,   201,   144,   201,    33,   144,   201,   126,
      82,   126,    82,   126,    82,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,   201,   144,    42,    51,   335,    42,   145,   146,
     273,   335,    51,   273,    51,    82,    51,    51,   204,   435,
     436,    51,    51,    82,    82,   433,   327,    51,    51,   335,
      51,   330,    51,   199,   200,    82,    42,    51,    33,    51,
     264,   199,   199,   199,   271,    82,   199,   199,   271,    82,
     227,   436,    51,    51,    51,    51,    51,   335,   335,   335,
      51,    51,    51,    51,    82,   200,   335,   290,   199,   271,
      82,    33,   144,     6,    42,    44,    51,    52,    82,    92,
     127,   145,   278,   286,   287,   144,   287,   144,   144,   287,
     144,    51,   145,   146,   272,    82,   199,    82,    31,   246,
     248,    33,   199,    45,    57,    63,    83,   237,   238,   351,
     352,   244,   199,   199,    56,    80,   315,    82,   147,   203,
     207,   208,   316,   317,   318,   144,    33,   144,   199,   227,
     227,   227,   228,   227,   227,   227,   227,   227,   227,   227,
     227,   227,   227,   201,   227,   227,   227,   227,   227,   227,
     227,   227,   227,   227,   227,   227,    82,   199,   144,   227,
      51,   335,    51,    51,    51,    51,    51,    51,   335,    51,
      51,    51,   199,   271,   126,   247,   247,   272,   247,   247,
     247,   247,   247,   247,   247,   247,   247,   247,   201,   247,
     247,   247,   247,   247,   247,   247,   247,   247,   247,   247,
     199,   293,   199,   271,   199,   271,   227,   199,   205,    42,
      51,   144,   200,   267,   199,   261,   199,   270,   199,   271,
     335,   324,   199,   271,   126,   126,   126,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,   335,   335,    51,   436,   335,   335,    51,    51,
     335,    51,   335,   335,   199,   333,    42,    42,    51,   434,
     205,   434,   203,   199,   199,    51,   309,   201,   201,   227,
     199,   201,   199,   201,   199,   206,   295,   296,   199,    82,
      82,    42,    51,   199,   144,   144,    82,   144,   287,    82,
     199,   287,    51,    51,   201,   232,    33,   247,    33,   144,
     201,   199,   242,   241,   144,   199,   200,   318,    82,   227,
      82,   100,   120,   144,   144,   144,   201,   201,   201,   201,
     201,   201,   201,   201,   201,   201,   201,   201,   227,    82,
     199,   199,   144,   144,   201,   201,   201,   201,   201,   201,
     201,   201,   201,   201,   201,   199,   199,   201,   258,   199,
      42,    51,   200,   227,   199,   199,   203,    82,   201,    33,
     199,   199,   271,   199,   271,    82,   144,   201,   279,   287,
     286,   287,   144,   287,   144,   144,   199,    33,    31,   247,
     199,   335,   238,   243,   245,   245,   245,   317,   287,    33,
     199,    33,    51,   254,   227,   227,   227,   227,   199,   199,
     227,   247,   247,   227,   201,    51,   309,   227,   199,   199,
     206,   295,   144,   144,   144,   287,   199,   287,   287,   227,
     199,   199,    31,    31,    31,   200,   201,   227,   227,   203,
     144,   199,   199,   201,   201,   199,   201,   201,   201,    33,
     199,   144,   287,   279,   287,   144,   199,   199,   199,   245,
     287,   199,   199,    51,    51,   130,   227,   206,   287,   144,
     144,   287,    31,   201,   203,   227,   256,   199,    82,   287,
     286,   199,    51,   144,   199,   206,   144,   144,   227,   287,
     279,   144,   287
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
     445,   446,   447,   448,   449,   450,   451,   452,   453,    59,
      40,    41,    35,    58,    91,    93,    39,    46,    92
  };
#endif

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
  const unsigned short int
  parser::yyr1_[] =
  {
         0,   209,   210,   210,   211,   211,   211,   211,   211,   211,
     211,   211,   211,   211,   211,   211,   211,   211,   211,   211,
     211,   211,   211,   211,   211,   211,   211,   211,   211,   211,
     211,   211,   211,   211,   211,   211,   211,   211,   211,   211,
     211,   211,   211,   211,   212,   212,   212,   212,   213,   213,
     214,   215,   216,   217,   218,   219,   219,   219,   219,   219,
     219,   220,   220,   220,   220,   220,   220,   221,   221,   221,
     221,   221,   221,   222,   222,   222,   222,   222,   222,   223,
     223,   224,   224,   225,   225,   226,   227,   227,   227,   227,
     227,   227,   227,   227,   227,   227,   227,   227,   227,   227,
     227,   227,   227,   227,   227,   227,   227,   227,   227,   227,
     227,   227,   227,   227,   227,   227,   228,   228,   229,   229,
     230,   231,   232,   232,   233,   234,   235,   235,   236,   237,
     237,   238,   238,   238,   238,   238,   240,   239,   241,   239,
     242,   239,   243,   239,   244,   239,   245,   245,   245,   245,
     246,   246,   247,   247,   247,   247,   247,   247,   247,   247,
     247,   247,   247,   247,   247,   247,   247,   247,   247,   247,
     247,   247,   247,   247,   247,   247,   247,   247,   247,   247,
     247,   248,   249,   249,   250,   251,   252,   252,   253,   253,
     253,   253,   253,   254,   254,   254,   254,   255,   256,   256,
     257,   257,   258,   258,   258,   258,   258,   258,   258,   258,
     258,   259,   259,   260,   260,   261,   262,   262,   263,   263,
     264,   265,   265,   266,   266,   267,   267,   268,   268,   268,
     268,   269,   269,   270,   270,   270,   270,   270,   270,   270,
     270,   270,   270,   270,   270,   270,   270,   270,   270,   270,
     270,   270,   270,   270,   270,   270,   271,   271,   271,   271,
     271,   271,   272,   272,   272,   273,   273,   273,   274,   275,
     275,   276,   277,   277,   277,   278,   278,   278,   278,   278,
     279,   279,   279,   279,   280,   281,   281,   282,   282,   282,
     283,   284,   284,   285,   285,   285,   286,   286,   286,   286,
     286,   287,   287,   287,   287,   287,   287,   288,   288,   288,
     288,   289,   289,   290,   290,   290,   290,   290,   290,   290,
     290,   290,   290,   290,   290,   290,   290,   290,   290,   290,
     290,   290,   290,   290,   290,   290,   290,   290,   290,   290,
     290,   290,   290,   290,   290,   290,   290,   290,   290,   290,
     290,   290,   291,   291,   292,   292,   293,   293,   293,   293,
     293,   293,   293,   293,   293,   293,   293,   293,   293,   294,
     294,   295,   295,   296,   296,   297,   298,   299,   299,   300,
     301,   302,   303,   303,   303,   303,   304,   305,   305,   305,
     305,   306,   307,   307,   308,   308,   308,   309,   309,   309,
     310,   310,   311,   311,   311,   311,   311,   311,   312,   312,
     312,   312,   312,   312,   313,   314,   314,   315,   315,   315,
     316,   316,   316,   316,   317,   317,   318,   318,   318,   318,
     318,   320,   321,   319,   322,   322,   322,   322,   323,   323,
     324,   324,   325,   325,   325,   325,   325,   325,   325,   326,
     326,   326,   326,   326,   326,   326,   326,   326,   326,   327,
     327,   328,   328,   329,   329,   329,   329,   330,   330,   331,
     331,   332,   332,   333,   333,   334,   334,   334,   334,   334,
     334,   334,   334,   334,   334,   334,   334,   334,   334,   334,
     334,   334,   334,   334,   334,   334,   334,   334,   334,   334,
     334,   334,   335,   335,   336,   337,   338,   339,   340,   341,
     342,   343,   344,   345,   346,   347,   348,   349,   350,   351,
     352,   353,   354,   355,   356,   357,   357,   358,   359,   360,
     361,   362,   363,   363,   364,   365,   366,   367,   368,   369,
     370,   371,   372,   373,   374,   375,   376,   377,   378,   379,
     380,   381,   382,   383,   384,   385,   386,   387,   388,   388,
     389,   390,   391,   392,   393,   394,   395,   396,   397,   398,
     399,   400,   401,   402,   403,   404,   405,   406,   407,   408,
     409,   410,   411,   412,   413,   414,   415,   416,   417,   418,
     419,   420,   421,   422,   423,   424,   425,   426,   427,   428,
     429,   430,   431,   432,   433,   434,   434,   435,   435,   436
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
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     2,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     6,     6,     4,     1,     3,     4,     7,
       3,     4,     2,     1,     4,     4,     2,     1,     7,     3,
       1,     1,     1,     1,     1,     1,     0,     5,     0,     8,
       0,     8,     0,    10,     0,     8,     2,     2,     1,     1,
       4,     2,     3,     1,     1,     1,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     2,     2,     4,
       4,     4,     4,     4,     4,     4,     4,     4,     4,     6,
       6,     5,     1,     4,     4,     4,     2,     1,     9,     6,
       5,     7,     7,     3,     5,     3,     1,     6,     3,     1,
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
       6,     2,     5,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     3,     3,     1,     3,     3,
       3,     3,     1,     1,     1,     3,     3,     3,     3,     3,
       3,     1,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     3,     3,     3,     3,     5,     3,
       3,     3,     1,     3,     3,     3,     1,     1,     1,     1,
       1,     3,     1,     1,     1,     1,     3,     3,     3,     3,
       1,     1,     3,     3,     3,     1,     1,     1,     3,     3,
       3,     3,     3,     3,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     1,     3,     2,     2,     2
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
  "VAROBS", "XLS_SHEET", "XLS_RANGE", "EXCLAMATION_EQUAL", "EXCLAMATION",
  "EQUAL_EQUAL", "GREATER_EQUAL", "LESS_EQUAL", "GREATER", "LESS", "COMMA",
  "MINUS", "PLUS", "DIVIDE", "TIMES", "UMINUS", "POWER", "EXP", "LOG",
  "LOG10", "SIN", "COS", "TAN", "ASIN", "ACOS", "ATAN", "SINH", "COSH",
  "TANH", "ASINH", "ACOSH", "ATANH", "SQRT", "DYNARE_SENSITIVITY",
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
  "period_list", "sigma_e", "value_list", "triangular_matrix",
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
       210,     0,    -1,   211,    -1,   210,   211,    -1,   212,    -1,
     223,    -1,   224,    -1,   225,    -1,   239,    -1,   229,    -1,
     231,    -1,   234,    -1,   226,    -1,   250,    -1,   251,    -1,
     255,    -1,   259,    -1,   262,    -1,   265,    -1,   268,    -1,
     288,    -1,   291,    -1,   294,    -1,   274,    -1,   283,    -1,
     280,    -1,   297,    -1,   298,    -1,   301,    -1,   213,    -1,
     214,    -1,   302,    -1,   304,    -1,   305,    -1,   306,    -1,
     310,    -1,   311,    -1,   312,    -1,   313,    -1,   319,    -1,
     322,    -1,   328,    -1,   331,    -1,   332,    -1,   218,    -1,
     215,    -1,   216,    -1,   217,    -1,    28,    51,   199,    -1,
      28,    51,    51,   199,    -1,   111,   271,   199,    -1,   131,
     219,   199,    -1,   132,   220,   199,    -1,   133,   221,   199,
      -1,    99,   222,   199,    -1,   219,    82,    -1,   219,   144,
      82,    -1,    82,    -1,   219,    82,   126,    -1,   219,   144,
      82,   126,    -1,    82,   126,    -1,   220,    82,    -1,   220,
     144,    82,    -1,    82,    -1,   220,    82,   126,    -1,   220,
     144,    82,   126,    -1,    82,   126,    -1,   221,    82,    -1,
     221,   144,    82,    -1,    82,    -1,   221,    82,   126,    -1,
     221,   144,    82,   126,    -1,    82,   126,    -1,   222,    82,
      -1,   222,   144,    82,    -1,    82,    -1,   222,    82,   126,
      -1,   222,   144,    82,   126,    -1,    82,   126,    -1,   100,
      51,   199,    -1,   100,    33,    51,   199,    -1,    24,    42,
     199,    -1,    24,    33,    42,   199,    -1,    63,    42,   199,
      -1,    63,    33,    42,   199,    -1,    82,    33,   227,   199,
      -1,   200,   227,   201,    -1,    82,    -1,    42,    -1,    51,
      -1,   227,   146,   227,    -1,   227,   145,   227,    -1,   227,
     147,   227,    -1,   227,   148,   227,    -1,   227,   150,   227,
      -1,   227,   143,   227,    -1,   227,   142,   227,    -1,   227,
     141,   227,    -1,   227,   140,   227,    -1,   227,   139,   227,
      -1,   227,   137,   227,    -1,   145,   227,    -1,   146,   227,
      -1,   151,   200,   227,   201,    -1,   152,   200,   227,   201,
      -1,   153,   200,   227,   201,    -1,   154,   200,   227,   201,
      -1,   155,   200,   227,   201,    -1,   156,   200,   227,   201,
      -1,   157,   200,   227,   201,    -1,   158,   200,   227,   201,
      -1,   159,   200,   227,   201,    -1,   166,   200,   227,   201,
      -1,    64,   200,   227,   144,   227,   201,    -1,    72,   200,
     227,   144,   227,   201,    -1,    82,   200,   228,   201,    -1,
     227,    -1,   228,   144,   227,    -1,    50,   199,   232,    31,
      -1,    50,   200,   230,   201,   199,   232,    31,    -1,    38,
      33,    82,    -1,    32,   199,   232,    31,    -1,   232,   233,
      -1,   233,    -1,    82,    33,   227,   199,    -1,    47,   199,
     235,    31,    -1,   235,   236,    -1,   236,    -1,    82,   200,
     272,   201,    33,   227,   199,    -1,   237,   144,   238,    -1,
     238,    -1,    57,    -1,    45,    -1,    83,    -1,   351,    -1,
     352,    -1,    -1,    76,   199,   240,   245,    31,    -1,    -1,
      76,   200,   339,   201,   199,   241,   245,    31,    -1,    -1,
      76,   200,   129,   201,   199,   242,   245,    31,    -1,    -1,
      76,   200,   119,   144,   237,   201,   243,   199,   245,    31,
      -1,    -1,    76,   200,   119,   201,   244,   199,   245,    31,
      -1,   245,   246,    -1,   245,   248,    -1,   246,    -1,   248,
      -1,   247,    33,   247,   199,    -1,   247,   199,    -1,   200,
     247,   201,    -1,   249,    -1,    42,    -1,    51,    -1,   247,
     146,   247,    -1,   247,   145,   247,    -1,   247,   147,   247,
      -1,   247,   148,   247,    -1,   247,   143,   247,    -1,   247,
     142,   247,    -1,   247,   141,   247,    -1,   247,   140,   247,
      -1,   247,   139,   247,    -1,   247,   137,   247,    -1,   247,
     150,   247,    -1,   145,   247,    -1,   146,   247,    -1,   151,
     200,   247,   201,    -1,   152,   200,   247,   201,    -1,   153,
     200,   247,   201,    -1,   154,   200,   247,   201,    -1,   155,
     200,   247,   201,    -1,   156,   200,   247,   201,    -1,   157,
     200,   247,   201,    -1,   158,   200,   247,   201,    -1,   159,
     200,   247,   201,    -1,   166,   200,   247,   201,    -1,    64,
     200,   247,   144,   247,   201,    -1,    72,   200,   247,   144,
     247,   201,    -1,   202,    82,    33,   247,   199,    -1,    82,
      -1,    82,   200,   272,   201,    -1,   112,   199,   252,    31,
      -1,    78,   199,   252,    31,    -1,   252,   253,    -1,   253,
      -1,   131,    82,   199,   100,   254,   199,   130,   256,   199,
      -1,   131,    82,   199,   120,   227,   199,    -1,   131,    82,
      33,   227,   199,    -1,   131,    82,   144,    82,    33,   227,
     199,    -1,    22,    82,   144,    82,    33,   227,   199,    -1,
     254,   144,    51,    -1,   254,   144,    51,   203,    51,    -1,
      51,   203,    51,    -1,    51,    -1,   113,    33,   204,   257,
     205,   199,    -1,   256,   144,   227,    -1,   227,    -1,   257,
     199,   258,    -1,   258,    -1,   258,   144,   200,   227,   201,
      -1,   258,   144,    42,    -1,   258,   144,    51,    -1,   258,
     200,   227,   201,    -1,   258,    42,    -1,   258,    51,    -1,
     200,   227,   201,    -1,    42,    -1,    51,    -1,   121,   199,
      -1,   121,   200,   260,   201,   199,    -1,   260,   144,   261,
      -1,   261,    -1,   337,    -1,    19,   199,    -1,    19,   200,
     263,   201,   199,    -1,   263,   144,   264,    -1,   264,    -1,
     337,    -1,   114,   199,    -1,   114,   200,   266,   201,   199,
      -1,   266,   144,   267,    -1,   267,    -1,   350,    -1,   356,
      -1,   122,   199,    -1,   122,   200,   269,   201,   199,    -1,
     122,   271,   199,    -1,   122,   200,   269,   201,   271,   199,
      -1,   269,   144,   270,    -1,   270,    -1,   336,    -1,   337,
      -1,   338,    -1,   339,    -1,   340,    -1,   341,    -1,   342,
      -1,   343,    -1,   344,    -1,   345,    -1,   346,    -1,   363,
      -1,   347,    -1,   385,    -1,   348,    -1,   349,    -1,   350,
      -1,   351,    -1,   353,    -1,   354,    -1,   355,    -1,   389,
      -1,   390,    -1,   271,    82,    -1,   271,    82,    33,    82,
      -1,   271,   144,    82,    -1,   271,   144,    82,    33,    82,
      -1,    82,    -1,    82,    33,    82,    -1,   146,    51,    -1,
     145,    51,    -1,    51,    -1,   146,    42,    -1,   145,    42,
      -1,    42,    -1,    35,   199,   275,    31,    -1,   275,   276,
      -1,   276,    -1,   277,   144,   278,   199,    -1,   120,    82,
      -1,    82,    -1,    22,    82,   144,    82,    -1,   286,   144,
     279,    -1,   287,   144,   286,   144,   279,    -1,   287,   144,
     287,   144,   287,   144,   286,   144,   279,    -1,   287,    -1,
     287,   144,   287,   144,   287,    -1,   287,   144,   287,    -1,
     287,   144,   287,   144,   287,    -1,   287,   144,   287,   144,
     287,   144,   287,    -1,   287,   144,   287,   144,   287,   144,
     287,   144,   287,    -1,    37,   199,   281,    31,    -1,   281,
     282,    -1,   282,    -1,   120,    82,   144,   287,   199,    -1,
      22,    82,   144,    82,   144,   287,   199,    -1,    82,   144,
     287,   199,    -1,    36,   199,   284,    31,    -1,   284,   285,
      -1,   285,    -1,   120,    82,   144,   287,   144,   287,   199,
      -1,    22,    82,   144,    82,   144,   287,   144,   287,   199,
      -1,    82,   144,   287,   144,   287,   199,    -1,     6,    -1,
      44,    -1,    92,    -1,    52,    -1,   127,    -1,    -1,    51,
      -1,    42,    -1,    82,    -1,   145,    51,    -1,   145,    42,
      -1,    34,   199,    -1,    34,   200,   289,   201,   199,    -1,
      34,   271,   199,    -1,    34,   200,   289,   201,   271,   199,
      -1,   289,   144,   290,    -1,   290,    -1,   356,    -1,   357,
      -1,   358,    -1,   359,    -1,   360,    -1,   361,    -1,   362,
      -1,   363,    -1,   364,    -1,   365,    -1,   366,    -1,   367,
      -1,   368,    -1,   369,    -1,   370,    -1,   371,    -1,   372,
      -1,   373,    -1,   374,    -1,   375,    -1,   376,    -1,   377,
      -1,   378,    -1,   379,    -1,   347,    -1,   380,    -1,   381,
      -1,   382,    -1,   383,    -1,   384,    -1,   386,    -1,   387,
      -1,   391,    -1,   392,    -1,   393,    -1,   337,    -1,   394,
      -1,   395,    -1,   396,    -1,   106,   200,   292,   201,   199,
      -1,   106,   200,   292,   201,   271,   199,    -1,   292,   144,
     293,    -1,   293,    -1,   363,    -1,   364,    -1,   373,    -1,
     379,    -1,   347,    -1,   380,    -1,   381,    -1,   382,    -1,
     383,    -1,   384,    -1,   391,    -1,   392,    -1,   393,    -1,
     107,   200,   292,   201,   199,    -1,   107,   200,   292,   201,
     271,   199,    -1,   206,    82,   206,   144,   206,    82,   206,
      -1,   206,    82,   206,   144,   287,    -1,   295,    -1,   296,
     144,   295,    -1,   134,   271,   199,    -1,    93,   199,   299,
      31,    -1,   299,   300,    -1,   300,    -1,    82,   200,   227,
     201,   199,    -1,   128,   271,   199,    -1,    95,   199,   303,
      31,    -1,   303,    82,   227,   199,    -1,   303,    82,   144,
      82,   227,   199,    -1,    82,   227,   199,    -1,    82,   144,
      82,   227,   199,    -1,    98,   271,   199,    -1,    97,   199,
      -1,    97,   200,   270,   201,   199,    -1,    97,   271,   199,
      -1,    97,   200,   270,   201,   271,   199,    -1,    18,   199,
     307,    31,    -1,   307,   308,    -1,   308,    -1,    82,   309,
      33,   227,   199,    -1,    82,   144,    82,   309,    33,   227,
     199,    -1,     4,    82,   200,    51,   201,   309,    33,   227,
     199,    -1,    -1,   200,    51,   201,    -1,   200,    42,   201,
      -1,    17,   199,    -1,    17,   200,    23,   201,   199,    -1,
      30,   200,    82,   201,   199,    -1,    30,   200,    82,   201,
     271,   199,    -1,    30,    82,   199,    -1,    30,   200,    82,
     207,    82,   201,   199,    -1,    30,   200,    82,   207,    82,
     201,   271,   199,    -1,    30,    82,   207,    82,   199,    -1,
      29,   200,    82,   201,   199,    -1,    29,   200,    82,   201,
     271,   199,    -1,    29,    82,   199,    -1,    29,   200,    82,
     207,    82,   201,   199,    -1,    29,   200,    82,   207,    82,
     201,   271,   199,    -1,    29,    82,   207,    82,   199,    -1,
      77,   200,   314,   201,   316,   199,    -1,   314,   144,   315,
      -1,   315,    -1,   388,    -1,   389,    -1,   390,    -1,   317,
      -1,   316,   144,   317,    -1,   317,   200,   287,   201,    -1,
     316,   144,   317,   200,   287,   201,    -1,   318,    -1,   317,
     318,    -1,    82,    -1,   208,    -1,   147,    -1,   203,    -1,
     207,    -1,    -1,    -1,   101,   320,   247,   321,   199,    -1,
     124,   199,    -1,   124,   200,   323,   201,   199,    -1,   124,
     271,   199,    -1,   124,   200,   323,   201,   271,   199,    -1,
     323,   144,   324,    -1,   324,    -1,   270,    -1,   397,    -1,
     398,    -1,   399,    -1,   400,    -1,   401,    -1,   402,    -1,
     403,    -1,   404,    -1,   325,    -1,   356,    -1,   391,    -1,
     392,    -1,   358,    -1,   360,    -1,   357,    -1,   359,    -1,
     394,    -1,   395,    -1,   326,   144,   327,    -1,   326,    -1,
       7,    51,   199,    -1,     7,   200,   327,   201,    51,   199,
      -1,   326,    -1,   381,    -1,   364,    -1,   405,    -1,   329,
     144,   330,    -1,   329,    -1,     8,    51,   199,    -1,     8,
     200,   330,   201,    51,   199,    -1,   167,   199,    -1,   167,
     200,   333,   201,   199,    -1,   334,   144,   333,    -1,   334,
      -1,   406,    -1,   407,    -1,   408,    -1,   409,    -1,   410,
      -1,   411,    -1,   412,    -1,   413,    -1,   414,    -1,   415,
      -1,   416,    -1,   417,    -1,   418,    -1,   419,    -1,   420,
      -1,   421,    -1,   422,    -1,   423,    -1,   425,    -1,   426,
      -1,   427,    -1,   428,    -1,   429,    -1,   430,    -1,   431,
      -1,   432,    -1,   424,    -1,    51,    -1,    42,    -1,    26,
      33,    51,    -1,   118,    33,    51,    -1,   115,    33,    51,
      -1,    60,    -1,    96,    33,    51,    -1,   110,    33,    51,
      -1,    27,    33,    51,    -1,     3,    33,    51,    -1,    86,
      -1,    88,    -1,    90,    -1,    53,    33,    51,    -1,    48,
      33,    51,    -1,    49,    33,    51,    -1,   100,    33,    51,
      -1,    24,    33,   335,    -1,    63,    33,   335,    -1,   114,
      -1,   116,    33,    51,    -1,   108,    33,   335,    -1,    25,
      33,    82,    -1,    84,    33,   436,    -1,    84,    33,    51,
      -1,    41,    33,    51,    -1,   102,    33,    51,    -1,   103,
      33,    51,    -1,    58,    33,    51,    -1,    59,    33,    51,
      -1,    89,    -1,    46,    -1,    20,    33,   335,    -1,    70,
      33,    51,    -1,    65,    33,   335,    -1,    67,    33,   335,
      -1,    94,    33,   200,   296,   201,    -1,    66,    33,   335,
      -1,    75,    33,    82,    -1,    74,    33,    51,    -1,    73,
      -1,   105,    33,   335,    -1,    68,    33,    51,    -1,    69,
      33,    51,    -1,    61,    -1,    62,    -1,    87,    -1,     5,
      -1,   123,    -1,    43,    33,    51,    -1,   117,    -1,    81,
      -1,    40,    -1,   109,    -1,    54,    33,    51,    -1,    55,
      33,    51,    -1,    79,    33,    56,    -1,    79,    33,    80,
      -1,   104,    -1,    91,    -1,   135,    33,    82,    -1,   136,
      33,   433,    -1,    39,    33,   436,    -1,    21,    -1,    85,
      -1,    71,    -1,   125,    33,   335,    -1,    14,    33,   273,
      -1,     9,    33,   335,    -1,    11,    33,   273,    -1,    12,
      33,   335,    -1,    13,    33,    51,    -1,    10,    -1,    15,
      33,    51,    -1,    16,    33,    51,    -1,   168,    33,    51,
      -1,   169,    33,    51,    -1,   170,    33,    51,    -1,   171,
      33,    51,    -1,   172,    33,    51,    -1,   173,    33,    51,
      -1,   174,    33,    51,    -1,   175,    33,    51,    -1,   176,
      33,    51,    -1,   177,    33,    51,    -1,   178,    33,    51,
      -1,   179,    33,    51,    -1,   180,    33,    51,    -1,   181,
      33,    51,    -1,   182,    33,    51,    -1,   183,    33,   335,
      -1,   184,    33,   335,    -1,   185,    33,    51,    -1,   186,
      33,   436,    -1,   187,    33,   335,    -1,   188,    33,   335,
      -1,   192,    33,    51,    -1,   193,    33,    51,    -1,   195,
      33,   335,    -1,   196,    33,    51,    -1,   197,    33,   335,
      -1,   198,    33,   335,    -1,    82,   203,    82,    -1,    51,
      -1,    51,   203,    51,    -1,   204,   434,    -1,   435,   434,
      -1,   435,   205,    -1
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
     251,   255,   259,   263,   267,   271,   275,   279,   283,   287,
     291,   295,   298,   301,   306,   311,   316,   321,   326,   331,
     336,   341,   346,   351,   358,   365,   370,   372,   376,   381,
     389,   393,   398,   401,   403,   408,   413,   416,   418,   426,
     430,   432,   434,   436,   438,   440,   442,   443,   449,   450,
     459,   460,   469,   470,   481,   482,   491,   494,   497,   499,
     501,   506,   509,   513,   515,   517,   519,   523,   527,   531,
     535,   539,   543,   547,   551,   555,   559,   563,   566,   569,
     574,   579,   584,   589,   594,   599,   604,   609,   614,   619,
     626,   633,   639,   641,   646,   651,   656,   659,   661,   671,
     678,   684,   692,   700,   704,   710,   714,   716,   723,   727,
     729,   733,   735,   741,   745,   749,   754,   757,   760,   764,
     766,   768,   771,   777,   781,   783,   785,   788,   794,   798,
     800,   802,   805,   811,   815,   817,   819,   821,   824,   830,
     834,   841,   845,   847,   849,   851,   853,   855,   857,   859,
     861,   863,   865,   867,   869,   871,   873,   875,   877,   879,
     881,   883,   885,   887,   889,   891,   893,   896,   901,   905,
     911,   913,   917,   920,   923,   925,   928,   931,   933,   938,
     941,   943,   948,   951,   953,   958,   962,   968,   978,   980,
     986,   990,   996,  1004,  1014,  1019,  1022,  1024,  1030,  1038,
    1043,  1048,  1051,  1053,  1061,  1071,  1078,  1080,  1082,  1084,
    1086,  1088,  1089,  1091,  1093,  1095,  1098,  1101,  1104,  1110,
    1114,  1121,  1125,  1127,  1129,  1131,  1133,  1135,  1137,  1139,
    1141,  1143,  1145,  1147,  1149,  1151,  1153,  1155,  1157,  1159,
    1161,  1163,  1165,  1167,  1169,  1171,  1173,  1175,  1177,  1179,
    1181,  1183,  1185,  1187,  1189,  1191,  1193,  1195,  1197,  1199,
    1201,  1203,  1205,  1211,  1218,  1222,  1224,  1226,  1228,  1230,
    1232,  1234,  1236,  1238,  1240,  1242,  1244,  1246,  1248,  1250,
    1256,  1263,  1271,  1277,  1279,  1283,  1287,  1292,  1295,  1297,
    1303,  1307,  1312,  1317,  1324,  1328,  1334,  1338,  1341,  1347,
    1351,  1358,  1363,  1366,  1368,  1374,  1382,  1392,  1393,  1397,
    1401,  1404,  1410,  1416,  1423,  1427,  1435,  1444,  1450,  1456,
    1463,  1467,  1475,  1484,  1490,  1497,  1501,  1503,  1505,  1507,
    1509,  1511,  1515,  1520,  1527,  1529,  1532,  1534,  1536,  1538,
    1540,  1542,  1543,  1544,  1550,  1553,  1559,  1563,  1570,  1574,
    1576,  1578,  1580,  1582,  1584,  1586,  1588,  1590,  1592,  1594,
    1596,  1598,  1600,  1602,  1604,  1606,  1608,  1610,  1612,  1614,
    1618,  1620,  1624,  1631,  1633,  1635,  1637,  1639,  1643,  1645,
    1649,  1656,  1659,  1665,  1669,  1671,  1673,  1675,  1677,  1679,
    1681,  1683,  1685,  1687,  1689,  1691,  1693,  1695,  1697,  1699,
    1701,  1703,  1705,  1707,  1709,  1711,  1713,  1715,  1717,  1719,
    1721,  1723,  1725,  1727,  1729,  1733,  1737,  1741,  1743,  1747,
    1751,  1755,  1759,  1761,  1763,  1765,  1769,  1773,  1777,  1781,
    1785,  1789,  1791,  1795,  1799,  1803,  1807,  1811,  1815,  1819,
    1823,  1827,  1831,  1833,  1835,  1839,  1843,  1847,  1851,  1857,
    1861,  1865,  1869,  1871,  1875,  1879,  1883,  1885,  1887,  1889,
    1891,  1893,  1897,  1899,  1901,  1903,  1905,  1909,  1913,  1917,
    1921,  1923,  1925,  1929,  1933,  1937,  1939,  1941,  1943,  1947,
    1951,  1955,  1959,  1963,  1967,  1969,  1973,  1977,  1981,  1985,
    1989,  1993,  1997,  2001,  2005,  2009,  2013,  2017,  2021,  2025,
    2029,  2033,  2037,  2041,  2045,  2049,  2053,  2057,  2061,  2065,
    2069,  2073,  2077,  2081,  2085,  2089,  2091,  2095,  2098,  2101
  };

  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
  const unsigned short int
  parser::yyrline_[] =
  {
         0,    97,    97,    98,   101,   102,   103,   104,   105,   106,
     107,   108,   109,   110,   111,   112,   113,   114,   115,   116,
     117,   118,   119,   120,   121,   122,   123,   124,   125,   126,
     127,   128,   129,   130,   131,   132,   133,   134,   135,   136,
     137,   138,   139,   140,   143,   144,   145,   146,   150,   152,
     156,   158,   160,   162,   164,   166,   168,   170,   172,   174,
     176,   180,   182,   184,   186,   188,   190,   194,   196,   198,
     200,   202,   204,   208,   210,   212,   214,   216,   218,   222,
     224,   228,   230,   234,   236,   241,   243,   245,   247,   249,
     251,   253,   255,   257,   259,   261,   263,   265,   267,   269,
     271,   273,   275,   277,   279,   281,   283,   285,   287,   289,
     291,   293,   295,   297,   299,   301,   305,   307,   311,   313,
     317,   319,   321,   322,   325,   327,   329,   330,   333,   335,
     336,   339,   341,   343,   345,   346,   349,   349,   351,   351,
     353,   353,   356,   355,   358,   358,   362,   363,   364,   365,
     368,   370,   374,   376,   377,   379,   381,   383,   385,   387,
     389,   391,   393,   395,   397,   399,   401,   403,   405,   407,
     409,   411,   413,   415,   417,   419,   421,   423,   425,   427,
     429,   433,   436,   438,   442,   444,   446,   447,   450,   452,
     454,   456,   458,   462,   464,   466,   468,   473,   476,   478,
     482,   484,   488,   490,   492,   494,   496,   498,   500,   502,
     504,   508,   510,   514,   515,   518,   520,   522,   526,   527,
     530,   532,   534,   538,   539,   542,   543,   546,   548,   550,
     552,   556,   557,   560,   561,   562,   563,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,   578,   579,   580,   581,   582,   585,   587,   589,   591,
     593,   595,   599,   601,   603,   607,   609,   611,   615,   617,
     619,   623,   625,   631,   637,   647,   652,   659,   670,   675,
     686,   693,   702,   713,   728,   731,   733,   737,   745,   755,
     765,   768,   770,   774,   784,   796,   808,   810,   812,   814,
     816,   820,   821,   822,   823,   824,   826,   830,   832,   834,
     836,   840,   841,   844,   845,   846,   847,   848,   849,   850,
     851,   852,   853,   854,   855,   856,   857,   858,   859,   860,
     861,   862,   863,   864,   865,   866,   867,   868,   869,   870,
     871,   872,   873,   874,   875,   876,   877,   878,   879,   880,
     881,   882,   885,   887,   891,   892,   895,   896,   897,   898,
     899,   900,   901,   902,   903,   904,   905,   906,   907,   910,
     912,   916,   918,   922,   923,   926,   928,   930,   931,   934,
     936,   938,   940,   942,   944,   946,   950,   952,   954,   956,
     958,   962,   964,   965,   968,   970,   972,   976,   977,   979,
     983,   985,   989,   991,   993,   995,   997,   999,  1003,  1005,
    1007,  1009,  1011,  1013,  1017,  1020,  1021,  1024,  1025,  1026,
    1029,  1031,  1033,  1035,  1039,  1041,  1045,  1046,  1048,  1050,
    1052,  1056,  1057,  1056,  1059,  1061,  1063,  1065,  1069,  1070,
    1073,  1074,  1077,  1078,  1079,  1080,  1081,  1082,  1083,  1086,
    1087,  1088,  1089,  1090,  1091,  1092,  1093,  1094,  1095,  1098,
    1099,  1102,  1104,  1108,  1109,  1110,  1111,  1114,  1115,  1118,
    1120,  1124,  1126,  1130,  1131,  1134,  1135,  1136,  1137,  1138,
    1139,  1140,  1141,  1142,  1143,  1144,  1145,  1146,  1147,  1148,
    1149,  1150,  1151,  1152,  1153,  1154,  1155,  1156,  1157,  1158,
    1159,  1160,  1163,  1164,  1167,  1168,  1169,  1170,  1171,  1172,
    1173,  1174,  1175,  1176,  1177,  1178,  1179,  1180,  1181,  1183,
    1184,  1185,  1186,  1187,  1188,  1189,  1191,  1194,  1195,  1196,
    1197,  1198,  1199,  1201,  1204,  1205,  1206,  1207,  1208,  1209,
    1210,  1211,  1212,  1213,  1214,  1215,  1216,  1217,  1218,  1219,
    1220,  1221,  1222,  1223,  1224,  1225,  1226,  1227,  1228,  1230,
    1233,  1234,  1235,  1236,  1237,  1238,  1239,  1240,  1241,  1243,
    1244,  1245,  1246,  1247,  1248,  1249,  1250,  1252,  1253,  1254,
    1255,  1256,  1257,  1258,  1259,  1260,  1261,  1262,  1263,  1264,
    1265,  1266,  1267,  1268,  1269,  1270,  1272,  1273,  1279,  1280,
    1284,  1285,  1286,  1287,  1290,  1298,  1299,  1308,  1310,  1319
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
       2,     2,     2,     2,     2,   202,     2,     2,     2,   206,
     200,   201,     2,     2,     2,     2,   207,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   203,   199,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   204,   208,   205,     2,     2,     2,     2,     2,     2,
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
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198
    };
    if ((unsigned int) t <= yyuser_token_number_max_)
      return translate_table[t];
    else
      return yyundef_token_;
  }

  const int parser::yyeof_ = 0;
  const int parser::yylast_ = 2192;
  const int parser::yynnts_ = 228;
  const int parser::yyempty_ = -2;
  const int parser::yyfinal_ = 160;
  const int parser::yyterror_ = 1;
  const int parser::yyerrcode_ = 256;
  const int parser::yyntokens_ = 209;

  const unsigned int parser::yyuser_token_number_max_ = 453;
  const parser::token_number_type parser::yyundef_token_ = 2;

} // namespace yy

#line 1321 "DynareBison.yy"


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

