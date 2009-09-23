VPATH = $(top_srcdir)/mex/sources/korderpert:$(top_srcdir)/mex/sources/korderpert/tests

check_PROGRAMS = k_order_test_main

CPPFLAGS += -DPOSIX_THREADS -I$(top_srcdir)/dynare++/src -I$(top_srcdir)/dynare++/kord -I$(top_srcdir)/dynare++/tl/cc -I$(top_srcdir)/dynare++/utils/cc -I$(top_srcdir)/dynare++/sylv/cc -I$(MATLAB)/extern/include -I$(top_srcdir)/mex/sources/korderpert

CPPFLAGS += -DMEXEXT=\"$(MEXEXT)\"

CXXFLAGS += $(PTHREAD_CFLAGS)

LIBS += $(LIBADD_DLOPEN) $(top_srcdir)/dynare++/kord/libkord.a $(top_srcdir)/dynare++/tl/cc/libtl.a $(top_srcdir)/dynare++/sylv/cc/libsylv.a $(top_srcdir)/dynare++/src/libnsolve.a $(MATLAB_LIBS)

LDFLAGS += $(subst -shared,,$(MATLAB_LDFLAGS))

LDFLAGS += -Wl,-rpath,$(MATLAB)/bin/$(MATLAB_ARCH)

nodist_k_order_test_main_SOURCES = \
	k_order_test_main.cpp \
	k_order_perturbation.cpp \
	k_ord_dynare.cpp \
	k_ord_dynare.h \
	dynamic_dll.cpp \
	dynamic_dll.h

