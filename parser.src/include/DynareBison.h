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
     SHOCK_SIZE = 359,
     SIGMA_E = 360,
     SIMUL = 361,
     SIMUL_ALGO = 362,
     SIMUL_SEED = 363,
     SMOOTHER = 364,
     SOLVE_ALGO = 365,
     STDERR = 366,
     STEADY = 367,
     STOCH_SIMUL = 368,
     TEX = 369,
     TEX_NAME = 370,
     UNIFORM_PDF = 371,
     UNIT_ROOT_VARS = 372,
     USE_DLL = 373,
     VALUES = 374,
     VAR = 375,
     VAREXO = 376,
     VAREXO_DET = 377,
     VAROBS = 378,
     XTICK = 379,
     XTICKLABEL = 380,
     COMMA = 381,
     MINUS = 382,
     PLUS = 383,
     DIVIDE = 384,
     TIMES = 385,
     UMINUS = 386,
     POWER = 387,
     FACTORIAL = 388,
     EXP = 389,
     LOG = 390,
     LOG10 = 391,
     LN = 392,
     SIN = 393,
     COS = 394,
     TAN = 395,
     ASIN = 396,
     ACOS = 397,
     ATAN = 398,
     SINH = 399,
     COSH = 400,
     TANH = 401,
     ASINH = 402,
     ACOSH = 403,
     ATANH = 404,
     SQRT = 405
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
#define SHOCK_SIZE 359
#define SIGMA_E 360
#define SIMUL 361
#define SIMUL_ALGO 362
#define SIMUL_SEED 363
#define SMOOTHER 364
#define SOLVE_ALGO 365
#define STDERR 366
#define STEADY 367
#define STOCH_SIMUL 368
#define TEX 369
#define TEX_NAME 370
#define UNIFORM_PDF 371
#define UNIT_ROOT_VARS 372
#define USE_DLL 373
#define VALUES 374
#define VAR 375
#define VAREXO 376
#define VAREXO_DET 377
#define VAROBS 378
#define XTICK 379
#define XTICKLABEL 380
#define COMMA 381
#define MINUS 382
#define PLUS 383
#define DIVIDE 384
#define TIMES 385
#define UMINUS 386
#define POWER 387
#define FACTORIAL 388
#define EXP 389
#define LOG 390
#define LOG10 391
#define LN 392
#define SIN 393
#define COS 394
#define TAN 395
#define ASIN 396
#define ACOS 397
#define ATAN 398
#define SINH 399
#define COSH 400
#define TANH 401
#define ASINH 402
#define ACOSH 403
#define ATANH 404
#define SQRT 405




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



