check_PROGRAMS = tests

tests_SOURCES = MMMatrix.cpp MMMatrix.h tests.cpp
tests_LDADD = ../cc/libsylv.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
tests_CPPFLAGS = -I../cc

EXTRA_DIST = tdata.tgz

check-local: tdata.tgz
	tar xzf tdata.tgz
	./tests

CLEANFILES = *.mm
