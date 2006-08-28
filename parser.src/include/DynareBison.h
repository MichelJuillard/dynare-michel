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
     LAPLACE = 304,
     LIK_ALGO = 305,
     LIK_INIT = 306,
     LINEAR = 307,
     LOAD_MH_FILE = 308,
     LOGLINEAR = 309,
     MH_DROP = 310,
     MH_INIT_SCALE = 311,
     MH_JSCALE = 312,
     MH_MODE = 313,
     MH_NBLOCKS = 314,
     MH_REPLIC = 315,
     MODE_CHECK = 316,
     MODE_COMPUTE = 317,
     MODE_FILE = 318,
     MODEL = 319,
     MODEL_COMPARISON = 320,
     MODEL_COMPARISON_APPROXIMATION = 321,
     MODIFIEDHARMONICMEAN = 322,
     MOMENTS = 323,
     MOMENTS_VARENDO = 324,
     MSHOCKS = 325,
     NAME = 326,
     NOBS = 327,
     NOCORR = 328,
     NODIAGNOSTIC = 329,
     NOFUNCTIONS = 330,
     NOGRAPH = 331,
     XLS_SHEET = 332,
     XLS_RANGE = 333,
     NOMOMENTS = 334,
     NOPRINT = 335,
     NORMAL_PDF = 336,
     OBSERVATION_TRENDS = 337,
     OLR = 338,
     OLR_INST = 339,
     OLR_BETA = 340,
     OPTIM = 341,
     OPTIM_WEIGHTS = 342,
     ORDER = 343,
     OSR = 344,
     OSR_PARAMS = 345,
     PARAMETERS = 346,
     PERIODS = 347,
     PREFILTER = 348,
     PRESAMPLE = 349,
     PRINT = 350,
     PRIOR_TRUNC = 351,
     FILTER_STEP_AHEAD = 352,
     QZ_CRITERIUM = 353,
     RELATIVE_IRF = 354,
     REPLIC = 355,
     RESOL = 356,
     RPLOT = 357,
     SHOCKS = 358,
     SIGMA_E = 359,
     SIMUL = 360,
     SIMUL_ALGO = 361,
     SIMUL_SEED = 362,
     SMOOTHER = 363,
     SOLVE_ALGO = 364,
     STDERR = 365,
     STEADY = 366,
     STOCH_SIMUL = 367,
     TEX = 368,
     TEX_NAME = 369,
     UNIFORM_PDF = 370,
     UNIT_ROOT_VARS = 371,
     USE_DLL = 372,
     VALUES = 373,
     VAR = 374,
     VAREXO = 375,
     VAREXO_DET = 376,
     VAROBS = 377,
     XTICK = 378,
     XTICKLABEL = 379,
     COMMA = 380,
     MINUS = 381,
     PLUS = 382,
     DIVIDE = 383,
     TIMES = 384,
     UMINUS = 385,
     POWER = 386,
     FACTORIAL = 387,
     EXP = 388,
     LOG = 389,
     LOG10 = 390,
     LN = 391,
     SIN = 392,
     COS = 393,
     TAN = 394,
     ASIN = 395,
     ACOS = 396,
     ATAN = 397,
     SINH = 398,
     COSH = 399,
     TANH = 400,
     ASINH = 401,
     ACOSH = 402,
     ATANH = 403,
     SQRT = 404,
     ASSIGN = 405
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
#define LAPLACE 304
#define LIK_ALGO 305
#define LIK_INIT 306
#define LINEAR 307
#define LOAD_MH_FILE 308
#define LOGLINEAR 309
#define MH_DROP 310
#define MH_INIT_SCALE 311
#define MH_JSCALE 312
#define MH_MODE 313
#define MH_NBLOCKS 314
#define MH_REPLIC 315
#define MODE_CHECK 316
#define MODE_COMPUTE 317
#define MODE_FILE 318
#define MODEL 319
#define MODEL_COMPARISON 320
#define MODEL_COMPARISON_APPROXIMATION 321
#define MODIFIEDHARMONICMEAN 322
#define MOMENTS 323
#define MOMENTS_VARENDO 324
#define MSHOCKS 325
#define NAME 326
#define NOBS 327
#define NOCORR 328
#define NODIAGNOSTIC 329
#define NOFUNCTIONS 330
#define NOGRAPH 331
#define XLS_SHEET 332
#define XLS_RANGE 333
#define NOMOMENTS 334
#define NOPRINT 335
#define NORMAL_PDF 336
#define OBSERVATION_TRENDS 337
#define OLR 338
#define OLR_INST 339
#define OLR_BETA 340
#define OPTIM 341
#define OPTIM_WEIGHTS 342
#define ORDER 343
#define OSR 344
#define OSR_PARAMS 345
#define PARAMETERS 346
#define PERIODS 347
#define PREFILTER 348
#define PRESAMPLE 349
#define PRINT 350
#define PRIOR_TRUNC 351
#define FILTER_STEP_AHEAD 352
#define QZ_CRITERIUM 353
#define RELATIVE_IRF 354
#define REPLIC 355
#define RESOL 356
#define RPLOT 357
#define SHOCKS 358
#define SIGMA_E 359
#define SIMUL 360
#define SIMUL_ALGO 361
#define SIMUL_SEED 362
#define SMOOTHER 363
#define SOLVE_ALGO 364
#define STDERR 365
#define STEADY 366
#define STOCH_SIMUL 367
#define TEX 368
#define TEX_NAME 369
#define UNIFORM_PDF 370
#define UNIT_ROOT_VARS 371
#define USE_DLL 372
#define VALUES 373
#define VAR 374
#define VAREXO 375
#define VAREXO_DET 376
#define VAROBS 377
#define XTICK 378
#define XTICKLABEL 379
#define COMMA 380
#define MINUS 381
#define PLUS 382
#define DIVIDE 383
#define TIMES 384
#define UMINUS 385
#define POWER 386
#define FACTORIAL 387
#define EXP 388
#define LOG 389
#define LOG10 390
#define LN 391
#define SIN 392
#define COS 393
#define TAN 394
#define ASIN 395
#define ACOS 396
#define ATAN 397
#define SINH 398
#define COSH 399
#define TANH 400
#define ASINH 401
#define ACOSH 402
#define ATANH 403
#define SQRT 404
#define ASSIGN 405




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



