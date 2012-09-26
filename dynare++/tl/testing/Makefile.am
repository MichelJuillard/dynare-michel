check_PROGRAMS = tests

tests_SOURCES = factory.cpp factory.h monoms.cpp monoms.h tests.cpp
tests_CPPFLAGS = -I../cc -I../../sylv/cc
tests_CXXFLAGS = $(PTHREAD_CFLAGS)
tests_LDFLAGS = $(LDFLAGS_MATIO)
tests_LDADD = ../cc/libtl.a ../../sylv/cc/libsylv.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(PTHREAD_LIBS) $(LIBADD_MATIO)

check-local:
	./tests
