check_PROGRAMS = tests

tests_SOURCES = tests.cpp
tests_CPPFLAGS = -I../cc -I../../tl/cc -I../../sylv/cc -I$(top_srcdir)/mex/sources
tests_CXXFLAGS = $(PTHREAD_CFLAGS)
tests_LDADD = ../../tl/cc/libtl.a ../../sylv/cc/libsylv.a ../cc/libinteg.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(PTHREAD_LIBS) $(LIBADD_MATIO)

check-local:
	./tests
