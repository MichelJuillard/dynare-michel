
#ifndef __DW_RANDOM__
#define __DW_RANDOM__

#ifdef __cplusplus
extern "C"
{
#endif

#include "prcsn.h"
#include <stdio.h>

void dw_initialize_generator(int init);

void* dw_get_generator_state(void);
int dw_get_generator_state_size(void);
void dw_set_generator_state(void *state);
void dw_print_generator_state(FILE *f);
void dw_read_generator_state(FILE *f);

PRECISION dw_uniform_rnd(void);
PRECISION dw_gaussian_rnd(void);
PRECISION dw_lognormal_rnd(PRECISION mean, PRECISION standard_deviation);
PRECISION dw_gamma_rnd(PRECISION a);

PRECISION dw_normal_cdf(PRECISION x);
PRECISION dw_chi_square_cdf(PRECISION x, int df);
PRECISION dw_chi_square_invcdf(PRECISION p, int df);
PRECISION dw_log_gamma(PRECISION x);

#ifdef __cplusplus
}
#endif

#endif
