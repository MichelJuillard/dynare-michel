ACLOCAL_AMFLAGS = -I ../../../m4

# libdynare++ must come before gensylv, k_order_perturbation, dynare_simul_
if DO_SOMETHING
SUBDIRS = mjdgges kronecker bytecode libdynare++ gensylv k_order_perturbation dynare_simul_ qzcomplex ordschur kalman_steady_state
if HAVE_FORT
SUBDIRS += libslicot
endif
if HAVE_GSL
SUBDIRS += swz
endif
if HAVE_MATIO
SUBDIRS += logposterior logMHMCMCposterior
endif
endif
