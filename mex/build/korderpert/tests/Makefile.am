vpath %.cpp $(top_srcdir)/mex/sources/korderpert/tests

check_PROGRAMS = k_order_test_main

CPPFLAGS += -I$(top_srcdir)/dynare++/src -I$(top_srcdir)/dynare++/kord -I$(top_srcdir)/dynare++/tl/cc -I$(top_srcdir)/dynare++/utils/cc -I$(top_srcdir)/dynare++/sylv/cc -I$(MATLAB)/extern/include -I$(top_srcdir)/mex/sources/korderpert -I$(top_srcdir)/mex/sources

CPPFLAGS += $(MATLAB_CPPFLAGS)

DEFS += -DPOSIX_THREADS -DMATLAB_MEX_FILE $(MATLAB_DEFS)

CXXFLAGS += $(PTHREAD_CFLAGS)

LDFLAGS += $(subst -shared,,$(MATLAB_LDFLAGS))

LDFLAGS += -Wl,-rpath,$(MATLAB)/bin/$(MATLAB_ARCH)

nodist_k_order_test_main_SOURCES = k_order_test_main.cpp

k_order_test_main_LDADD = ../matlab/libkorderpert.a $(MATLAB_LIBS) $(PTHREAD_LIBS) $(LIBADD_DLOPEN)
