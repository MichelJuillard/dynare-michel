CWEBSRC = \
	faa_di_bruno.cweb \
	korder_stoch.cweb \
	journal.cweb \
	decision_rule.cweb \
	dynamic_model.cweb \
	random.cweb \
	first_order.cweb \
	normal_conjugate.cweb \
	approximation.cweb \
	global_check.cweb \
	korder.cweb \
	kord_exception.hweb \
	random.hweb \
	journal.hweb \
	approximation.hweb \
	korder_stoch.hweb \
	dynamic_model.hweb \
	decision_rule.hweb \
	korder.hweb \
	normal_conjugate.hweb \
	first_order.hweb \
	mersenne_twister.hweb \
	global_check.hweb \
	faa_di_bruno.hweb

GENERATED_FILES = \
	faa_di_bruno.cpp \
	korder_stoch.cpp \
	journal.cpp \
	decision_rule.cpp \
	dynamic_model.cpp \
	random.cpp \
	first_order.cpp \
	normal_conjugate.cpp \
	approximation.cpp \
	global_check.cpp \
	korder.cpp \
	kord_exception.h \
	random.h \
	journal.h \
	approximation.h \
	korder_stoch.h \
	dynamic_model.h \
	decision_rule.h \
	korder.h \
	normal_conjugate.h \
	first_order.h \
	mersenne_twister.h \
	global_check.h \
	faa_di_bruno.h

noinst_LIBRARIES = libkord.a

libkord_a_SOURCES = $(CWEBSRC) $(GENERATED_FILES)
libkord_a_CPPFLAGS = -I../sylv/cc -I../tl/cc -I../integ/cc -I$(top_srcdir)/mex/sources $(CPPFLAGS_MATIO)
libkord_a_CXXFLAGS = $(PTHREAD_CFLAGS)

BUILT_SOURCES = $(GENERATED_FILES)

EXTRA_DIST = main.web dummy.ch

check_PROGRAMS = tests

tests_SOURCES = tests.cpp
tests_CPPFLAGS = -I../sylv/cc -I../tl/cc -I../integ/cc -I$(top_srcdir)/mex/sources
tests_CXXFLAGS = $(PTHREAD_CFLAGS)
tests_LDADD = libkord.a ../tl/cc/libtl.a ../sylv/cc/libsylv.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(PTHREAD_LIBS) $(LIBADD_MATIO)

check-local:
	./tests

%.cpp: %.cweb dummy.ch
	$(CTANGLE) -bhp $< dummy.ch $@

%.h: %.hweb dummy.ch
	$(CTANGLE) -bhp $< dummy.ch $@

if HAVE_CWEAVE
if HAVE_PDFTEX
if HAVE_EPLAIN
pdf-local: kord.pdf

kord.pdf: main.web $(CWEBSRC)
	$(CWEAVE) -bhp main.web
	$(PDFTEX) main
	mv main.pdf kord.pdf
endif
endif
endif

CLEANFILES = kord.pdf main.idx main.log main.scn main.tex main.toc out.txt
