CWEBSRC = \
	product.cweb \
	quadrature.cweb \
	quasi_mcarlo.cweb \
	smolyak.cweb \
	vector_function.cweb \
	product.hweb \
	quadrature.hweb \
	quasi_mcarlo.hweb \
	smolyak.hweb \
	vector_function.hweb

GENERATED_FILES = \
	product.cpp \
	quadrature.cpp \
	quasi_mcarlo.cpp \
	smolyak.cpp \
	vector_function.cpp \
	product.h \
	quadrature.h \
	quasi_mcarlo.h \
	smolyak.h \
	vector_function.h

noinst_LIBRARIES = libinteg.a

libinteg_a_SOURCES = $(CWEBSRC) $(GENERATED_FILES) precalc_quadrature.dat
libinteg_a_CPPFLAGS = -I../../sylv/cc -I../../tl/cc -DPOSIX_THREADS
libinteg_a_CXXFLAGS = $(PTHREAD_CFLAGS)

BUILT_SOURCES = $(GENERATED_FILES)

EXTRA_DIST = main.web

dummy.ch:
	touch dummy.ch

# Don't depend on dummy.ch in the following rules, otherwise make will want to rebuild from CWEB even in a tarball, since dummy.ch has a recent timestamp

%.cpp: %.cweb
	$(CTANGLE) -bhp $< dummy.ch $@

%.h: %.hweb
	$(CTANGLE) -bhp $< dummy.ch $@

if HAVE_CWEAVE
if HAVE_PDFTEX
if HAVE_EPLAIN
pdf: main.web $(CWEBSRC)
	$(CWEAVE) -bhp main.web
	$(PDFTEX) main
	mv main.pdf integ.pdf
endif
endif
endif

CLEANFILES = dummy.ch integ.pdf main.idx main.log main.scn main.tex main.toc
