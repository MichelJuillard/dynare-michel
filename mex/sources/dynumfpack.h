/*
 * Defines some prototypes for UMFPACK functions
 */

/*
 * Copyright (C) 2013 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _DYNUMFPACK_H
#define _DYNUMFPACK_H

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------------- */
/* size of Info and Control arrays */
/* -------------------------------------------------------------------------- */

/* These might be larger in future versions, since there are only 3 unused
 * entries in Info, and no unused entries in Control. */

#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20
/* used in all UMFPACK_report_* routines: */
#define UMFPACK_PRL 0			/* print level */
/* returned by all routines that use Info: */
#define UMFPACK_OK (0)
#define UMFPACK_STATUS 0	/* UMFPACK_OK, or other result */

#ifdef _WIN64
typedef long long int SuiteSparse_long;
#else
typedef long SuiteSparse_long;
#endif

void umfpack_dl_defaults(double Control[UMFPACK_CONTROL]);

SuiteSparse_long umfpack_dl_symbolic(SuiteSparse_long n_row, SuiteSparse_long n_col,
                                     const SuiteSparse_long Ap [ ], const SuiteSparse_long Ai [ ],
                                     const double Ax [ ], void **Symbolic,
                                     const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);

SuiteSparse_long umfpack_dl_numeric(const SuiteSparse_long Ap [ ], const SuiteSparse_long Ai [ ],
                                    const double Ax [ ], void *Symbolic, void **Numeric,
                                    const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);

SuiteSparse_long umfpack_dl_solve(SuiteSparse_long sys, const SuiteSparse_long Ap [ ],
                                  const SuiteSparse_long Ai [ ], const double Ax [ ],
                                  double X [ ], const double B [ ], void *Numeric,
                                  const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);

void umfpack_dl_report_info(const double Control [UMFPACK_CONTROL],
                            const double Info [UMFPACK_INFO]);

void umfpack_dl_report_status(const double Control [UMFPACK_CONTROL],
                              SuiteSparse_long status);

void umfpack_dl_free_symbolic(void **Symbolic);

void umfpack_dl_free_numeric(void **Numeric);

SuiteSparse_long umfpack_dl_load_symbolic (void **Symbolic, char *filename) ;

SuiteSparse_long umfpack_dl_load_numeric (void **Numeric, char *filename) ;

SuiteSparse_long umfpack_dl_save_symbolic (void *Symbolic, char *filename) ;

SuiteSparse_long umfpack_dl_save_numeric (void *Numeric, char *filename) ;

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* DYNUMFPACK */
