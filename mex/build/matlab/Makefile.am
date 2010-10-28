ACLOCAL_AMFLAGS = -I ../../../m4

# libdynare++ must come before gensylv, k_order_perturbation, dynare_simul_
if DO_SOMETHING
SUBDIRS = mjdgges kronecker bytecode libdynare++ gensylv k_order_perturbation dynare_simul_ logposterior logMHMCMCposterior
if HAVE_GSL
SUBDIRS += swz
endif
endif

EXTRA_DIST = mex.def mexFunction-MacOSX.map
