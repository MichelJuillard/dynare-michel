/* A Bison parser, made by GNU Bison 2.1.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

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
     DLL = 271,
     DOLLAR = 272,
     DR_ALGO = 273,
     DROP = 274,
     DSAMPLE = 275,
     DYN2VEC = 276,
     DYNASAVE = 277,
     DYNATYPE = 278,
     END = 279,
     ENDVAL = 280,
     EQUAL = 281,
     ESTIMATION = 282,
     ESTIMATED_PARAMS = 283,
     ESTIMATED_PARAMS_BOUNDS = 284,
     ESTIMATED_PARAMS_INIT = 285,
     FILTERED_VARS = 286,
     FIRST_OBS = 287,
     FLOAT_NUMBER = 288,
     FORECAST = 289,
     FUNCTIONS = 290,
     GAMMA_PDF = 291,
     GRAPH = 292,
     HISTVAL = 293,
     HP_FILTER = 294,
     HP_NGRID = 295,
     INITVAL = 296,
     INITVALF = 297,
     INT_NUMBER = 298,
     INV_GAMMA_PDF = 299,
     INV_GAMMA1_PDF = 300,
     INV_GAMMA2_PDF = 301,
     IRF = 302,
     KALMAN_ALGO = 303,
     KALMAN_TOL = 304,
     CONSTANT = 305,
     NOCONSTANT = 306,
     LAPLACE = 307,
     LIK_ALGO = 308,
     LIK_INIT = 309,
     LINEAR = 310,
     LOAD_MH_FILE = 311,
     LOGLINEAR = 312,
     MH_DROP = 313,
     MH_INIT_SCALE = 314,
     MH_JSCALE = 315,
     MH_MODE = 316,
     MH_NBLOCKS = 317,
     MH_REPLIC = 318,
     MODE_CHECK = 319,
     MODE_COMPUTE = 320,
     MODE_FILE = 321,
     MODEL = 322,
     MODEL_COMPARISON = 323,
     MODEL_COMPARISON_APPROXIMATION = 324,
     MODIFIEDHARMONICMEAN = 325,
     MOMENTS = 326,
     MOMENTS_VARENDO = 327,
     MSHOCKS = 328,
     NAME = 329,
     NOBS = 330,
     NOCORR = 331,
     NODIAGNOSTIC = 332,
     NOFUNCTIONS = 333,
     NOGRAPH = 334,
     XLS_SHEET = 335,
     XLS_RANGE = 336,
     NOMOMENTS = 337,
     NOPRINT = 338,
     NORMAL_PDF = 339,
     OBSERVATION_TRENDS = 340,
     OLR = 341,
     OLR_INST = 342,
     OLR_BETA = 343,
     OPTIM = 344,
     OPTIM_WEIGHTS = 345,
     ORDER = 346,
     OSR = 347,
     OSR_PARAMS = 348,
     PARAMETERS = 349,
     PERIODS = 350,
     PREFILTER = 351,
     PRESAMPLE = 352,
     PRINT = 353,
     PRIOR_TRUNC = 354,
     FILTER_STEP_AHEAD = 355,
     QZ_CRITERIUM = 356,
     RELATIVE_IRF = 357,
     REPLIC = 358,
     RESOL = 359,
     RPLOT = 360,
     SHOCKS = 361,
     SIGMA_E = 362,
     SIMUL = 363,
     SIMUL_ALGO = 364,
     SIMUL_SEED = 365,
     SMOOTHER = 366,
     SOLVE_ALGO = 367,
     STDERR = 368,
     STEADY = 369,
     STOCH_SIMUL = 370,
     TEX = 371,
     TEX_NAME = 372,
     UNIFORM_PDF = 373,
     UNIT_ROOT_VARS = 374,
     USE_DLL = 375,
     VALUES = 376,
     VAR = 377,
     VAREXO = 378,
     VAREXO_DET = 379,
     VAROBS = 380,
     XTICK = 381,
     XTICKLABEL = 382,
     COMMA = 383,
     MINUS = 384,
     PLUS = 385,
     DIVIDE = 386,
     TIMES = 387,
     UMINUS = 388,
     POWER = 389,
     FACTORIAL = 390,
     EXP = 391,
     LOG = 392,
     LOG10 = 393,
     LN = 394,
     SIN = 395,
     COS = 396,
     TAN = 397,
     ASIN = 398,
     ACOS = 399,
     ATAN = 400,
     SINH = 401,
     COSH = 402,
     TANH = 403,
     ASINH = 404,
     ACOSH = 405,
     ATANH = 406,
     SQRT = 407,
     ASSIGN = 408
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
#define DLL 271
#define DOLLAR 272
#define DR_ALGO 273
#define DROP 274
#define DSAMPLE 275
#define DYN2VEC 276
#define DYNASAVE 277
#define DYNATYPE 278
#define END 279
#define ENDVAL 280
#define EQUAL 281
#define ESTIMATION 282
#define ESTIMATED_PARAMS 283
#define ESTIMATED_PARAMS_BOUNDS 284
#define ESTIMATED_PARAMS_INIT 285
#define FILTERED_VARS 286
#define FIRST_OBS 287
#define FLOAT_NUMBER 288
#define FORECAST 289
#define FUNCTIONS 290
#define GAMMA_PDF 291
#define GRAPH 292
#define HISTVAL 293
#define HP_FILTER 294
#define HP_NGRID 295
#define INITVAL 296
#define INITVALF 297
#define INT_NUMBER 298
#define INV_GAMMA_PDF 299
#define INV_GAMMA1_PDF 300
#define INV_GAMMA2_PDF 301
#define IRF 302
#define KALMAN_ALGO 303
#define KALMAN_TOL 304
#define CONSTANT 305
#define NOCONSTANT 306
#define LAPLACE 307
#define LIK_ALGO 308
#define LIK_INIT 309
#define LINEAR 310
#define LOAD_MH_FILE 311
#define LOGLINEAR 312
#define MH_DROP 313
#define MH_INIT_SCALE 314
#define MH_JSCALE 315
#define MH_MODE 316
#define MH_NBLOCKS 317
#define MH_REPLIC 318
#define MODE_CHECK 319
#define MODE_COMPUTE 320
#define MODE_FILE 321
#define MODEL 322
#define MODEL_COMPARISON 323
#define MODEL_COMPARISON_APPROXIMATION 324
#define MODIFIEDHARMONICMEAN 325
#define MOMENTS 326
#define MOMENTS_VARENDO 327
#define MSHOCKS 328
#define NAME 329
#define NOBS 330
#define NOCORR 331
#define NODIAGNOSTIC 332
#define NOFUNCTIONS 333
#define NOGRAPH 334
#define XLS_SHEET 335
#define XLS_RANGE 336
#define NOMOMENTS 337
#define NOPRINT 338
#define NORMAL_PDF 339
#define OBSERVATION_TRENDS 340
#define OLR 341
#define OLR_INST 342
#define OLR_BETA 343
#define OPTIM 344
#define OPTIM_WEIGHTS 345
#define ORDER 346
#define OSR 347
#define OSR_PARAMS 348
#define PARAMETERS 349
#define PERIODS 350
#define PREFILTER 351
#define PRESAMPLE 352
#define PRINT 353
#define PRIOR_TRUNC 354
#define FILTER_STEP_AHEAD 355
#define QZ_CRITERIUM 356
#define RELATIVE_IRF 357
#define REPLIC 358
#define RESOL 359
#define RPLOT 360
#define SHOCKS 361
#define SIGMA_E 362
#define SIMUL 363
#define SIMUL_ALGO 364
#define SIMUL_SEED 365
#define SMOOTHER 366
#define SOLVE_ALGO 367
#define STDERR 368
#define STEADY 369
#define STOCH_SIMUL 370
#define TEX 371
#define TEX_NAME 372
#define UNIFORM_PDF 373
#define UNIT_ROOT_VARS 374
#define USE_DLL 375
#define VALUES 376
#define VAR 377
#define VAREXO 378
#define VAREXO_DET 379
#define VAROBS 380
#define XTICK 381
#define XTICKLABEL 382
#define COMMA 383
#define MINUS 384
#define PLUS 385
#define DIVIDE 386
#define TIMES 387
#define UMINUS 388
#define POWER 389
#define FACTORIAL 390
#define EXP 391
#define LOG 392
#define LOG10 393
#define LN 394
#define SIN 395
#define COS 396
#define TAN 397
#define ASIN 398
#define ACOS 399
#define ATAN 400
#define SINH 401
#define COSH 402
#define TANH 403
#define ASINH 404
#define ACOSH 405
#define ATANH 406
#define SQRT 407
#define ASSIGN 408




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



