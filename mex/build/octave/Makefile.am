ACLOCAL_AMFLAGS = -I ../../../m4

# libdynare++ must come before gensylv, k_order_perturbation, dynare_simul_
if DO_SOMETHING
SUBDIRS = mjdgges kronecker bytecode libdynare++ gensylv k_order_perturbation dynare_simul_ qzcomplex ordschur block_kalman_filter sobol local_state_space_iterations linsolve

if HAVE_GSL
if HAVE_MATIO
SUBDIRS += ms_sbvar
endif
endif

# libslicot must come before kalman_steady_state
if HAVE_FORT
SUBDIRS += libslicot kalman_steady_state
endif

if HAVE_MATIO
SUBDIRS += estimation
endif
endif
