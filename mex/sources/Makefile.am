SUBDIRS = estimation

EXTRA_DIST = \
	dynblas.h \
	dynlapack.h \
	dynmex.h \
	mjdgges \
	kronecker \
	bytecode \
	k_order_perturbation

clean-local:
	rm -rf `find mex/sources -name *.o`

dist-hook:
	rm -rf `find $(distdir)/mex/sources -name *.o -or -name *~`
