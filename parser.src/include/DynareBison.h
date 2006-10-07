/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

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

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     AR = 258,
     AUTOCORR = 259,
     BAYESIAN_IRF = 260,
     BETA_PDF = 261,
     CALIB = 262,
     CALIB_VAR = 263,
     CHECK = 264,
     CONF_SIG = 265,
     CORR = 266,
     COVAR = 267,
     DATAFILE = 268,
     DIAGNOSTIC = 269,
     DIFFUSE_D = 270,
     DOLLAR = 271,
     DR_ALGO = 272,
     DROP = 273,
     DSAMPLE = 274,
     DYN2VEC = 275,
     DYNASAVE = 276,
     DYNATYPE = 277,
     END = 278,
     ENDVAL = 279,
     EQUAL = 280,
     ESTIMATION = 281,
     ESTIMATED_PARAMS = 282,
     ESTIMATED_PARAMS_BOUNDS = 283,
     ESTIMATED_PARAMS_INIT = 284,
     FILTERED_VARS = 285,
     FIRST_OBS = 286,
     FLOAT_NUMBER = 287,
     FORECAST = 288,
     FUNCTIONS = 289,
     GAMMA_PDF = 290,
     GRAPH = 291,
     HISTVAL = 292,
     HP_FILTER = 293,
     HP_NGRID = 294,
     INITVAL = 295,
     INITVALF = 296,
     INT_NUMBER = 297,
     INV_GAMMA_PDF = 298,
     INV_GAMMA1_PDF = 299,
     INV_GAMMA2_PDF = 300,
     IRF = 301,
     KALMAN_ALGO = 302,
     KALMAN_TOL = 303,
     CONSTANT = 304,
     NOCONSTANT = 305,
     LAPLACE = 306,
     LIK_ALGO = 307,
     LIK_INIT = 308,
     LINEAR = 309,
     LOAD_MH_FILE = 310,
     LOGLINEAR = 311,
     MH_DROP = 312,
     MH_INIT_SCALE = 313,
     MH_JSCALE = 314,
     MH_MODE = 315,
     MH_NBLOCKS = 316,
     MH_REPLIC = 317,
     MODE_CHECK = 318,
     MODE_COMPUTE = 319,
     MODE_FILE = 320,
     MODEL = 321,
     MODEL_COMPARISON = 322,
     MODEL_COMPARISON_APPROXIMATION = 323,
     MODIFIEDHARMONICMEAN = 324,
     MOMENTS = 325,
     MOMENTS_VARENDO = 326,
     MSHOCKS = 327,
     NAME = 328,
     NOBS = 329,
     NOCORR = 330,
     NODIAGNOSTIC = 331,
     NOFUNCTIONS = 332,
     NOGRAPH = 333,
     XLS_SHEET = 334,
     XLS_RANGE = 335,
     NOMOMENTS = 336,
     NOPRINT = 337,
     NORMAL_PDF = 338,
     OBSERVATION_TRENDS = 339,
     OLR = 340,
     OLR_INST = 341,
     OLR_BETA = 342,
     OPTIM = 343,
     OPTIM_WEIGHTS = 344,
     ORDER = 345,
     OSR = 346,
     OSR_PARAMS = 347,
     PARAMETERS = 348,
     PERIODS = 349,
     PREFILTER = 350,
     PRESAMPLE = 351,
     PRINT = 352,
     PRIOR_TRUNC = 353,
     FILTER_STEP_AHEAD = 354,
     QZ_CRITERIUM = 355,
     RELATIVE_IRF = 356,
     REPLIC = 357,
     RESOL = 358,
     RPLOT = 359,
     SHOCKS = 360,
     SIGMA_E = 361,
     SIMUL = 362,
     SIMUL_ALGO = 363,
     SIMUL_SEED = 364,
     SMOOTHER = 365,
     SOLVE_ALGO = 366,
     STDERR = 367,
     STEADY = 368,
     STOCH_SIMUL = 369,
     TEX = 370,
     TEX_NAME = 371,
     UNIFORM_PDF = 372,
     UNIT_ROOT_VARS = 373,
     USE_DLL = 374,
     VALUES = 375,
     VAR = 376,
     VAREXO = 377,
     VAREXO_DET = 378,
     VAROBS = 379,
     XTICK = 380,
     XTICKLABEL = 381,
     COMMA = 382,
     MINUS = 383,
     PLUS = 384,
     DIVIDE = 385,
     TIMES = 386,
     UMINUS = 387,
     POWER = 388,
     FACTORIAL = 389,
     EXP = 390,
     LOG = 391,
     LOG10 = 392,
     LN = 393,
     SIN = 394,
     COS = 395,
     TAN = 396,
     ASIN = 397,
     ACOS = 398,
     ATAN = 399,
     SINH = 400,
     COSH = 401,
     TANH = 402,
     ASINH = 403,
     ACOSH = 404,
     ATANH = 405,
     SQRT = 406,
     ASSIGN = 407
   };
#endif
/* Tokens.  */
#define AR 258
#define AUTOCORR 259
#define BAYESIAN_IRF 260
#define BETA_PDF 261
#define CALIB 262
#define CALIB_VAR 263
#define CHECK 264
#define CONF_SIG 265
#define CORR 266
#define COVAR 267
#define DATAFILE 268
#define DIAGNOSTIC 269
#define DIFFUSE_D 270
#define DOLLAR 271
#define DR_ALGO 272
#define DROP 273
#define DSAMPLE 274
#define DYN2VEC 275
#define DYNASAVE 276
#define DYNATYPE 277
#define END 278
#define ENDVAL 279
#define EQUAL 280
#define ESTIMATION 281
#define ESTIMATED_PARAMS 282
#define ESTIMATED_PARAMS_BOUNDS 283
#define ESTIMATED_PARAMS_INIT 284
#define FILTERED_VARS 285
#define FIRST_OBS 286
#define FLOAT_NUMBER 287
#define FORECAST 288
#define FUNCTIONS 289
#define GAMMA_PDF 290
#define GRAPH 291
#define HISTVAL 292
#define HP_FILTER 293
#define HP_NGRID 294
#define INITVAL 295
#define INITVALF 296
#define INT_NUMBER 297
#define INV_GAMMA_PDF 298
#define INV_GAMMA1_PDF 299
#define INV_GAMMA2_PDF 300
#define IRF 301
#define KALMAN_ALGO 302
#define KALMAN_TOL 303
#define CONSTANT 304
#define NOCONSTANT 305
#define LAPLACE 306
#define LIK_ALGO 307
#define LIK_INIT 308
#define LINEAR 309
#define LOAD_MH_FILE 310
#define LOGLINEAR 311
#define MH_DROP 312
#define MH_INIT_SCALE 313
#define MH_JSCALE 314
#define MH_MODE 315
#define MH_NBLOCKS 316
#define MH_REPLIC 317
#define MODE_CHECK 318
#define MODE_COMPUTE 319
#define MODE_FILE 320
#define MODEL 321
#define MODEL_COMPARISON 322
#define MODEL_COMPARISON_APPROXIMATION 323
#define MODIFIEDHARMONICMEAN 324
#define MOMENTS 325
#define MOMENTS_VARENDO 326
#define MSHOCKS 327
#define NAME 328
#define NOBS 329
#define NOCORR 330
#define NODIAGNOSTIC 331
#define NOFUNCTIONS 332
#define NOGRAPH 333
#define XLS_SHEET 334
#define XLS_RANGE 335
#define NOMOMENTS 336
#define NOPRINT 337
#define NORMAL_PDF 338
#define OBSERVATION_TRENDS 339
#define OLR 340
#define OLR_INST 341
#define OLR_BETA 342
#define OPTIM 343
#define OPTIM_WEIGHTS 344
#define ORDER 345
#define OSR 346
#define OSR_PARAMS 347
#define PARAMETERS 348
#define PERIODS 349
#define PREFILTER 350
#define PRESAMPLE 351
#define PRINT 352
#define PRIOR_TRUNC 353
#define FILTER_STEP_AHEAD 354
#define QZ_CRITERIUM 355
#define RELATIVE_IRF 356
#define REPLIC 357
#define RESOL 358
#define RPLOT 359
#define SHOCKS 360
#define SIGMA_E 361
#define SIMUL 362
#define SIMUL_ALGO 363
#define SIMUL_SEED 364
#define SMOOTHER 365
#define SOLVE_ALGO 366
#define STDERR 367
#define STEADY 368
#define STOCH_SIMUL 369
#define TEX 370
#define TEX_NAME 371
#define UNIFORM_PDF 372
#define UNIT_ROOT_VARS 373
#define USE_DLL 374
#define VALUES 375
#define VAR 376
#define VAREXO 377
#define VAREXO_DET 378
#define VAROBS 379
#define XTICK 380
#define XTICKLABEL 381
#define COMMA 382
#define MINUS 383
#define PLUS 384
#define DIVIDE 385
#define TIMES 386
#define UMINUS 387
#define POWER 388
#define FACTORIAL 389
#define EXP 390
#define LOG 391
#define LOG10 392
#define LN 393
#define SIN 394
#define COS 395
#define TAN 396
#define ASIN 397
#define ACOS 398
#define ATAN 399
#define SINH 400
#define COSH 401
#define TANH 402
#define ASINH 403
#define ACOSH 404
#define ATANH 405
#define SQRT 406
#define ASSIGN 407




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

