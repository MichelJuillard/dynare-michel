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

EXTRA_DIST = main.web dummy.ch

%.cpp: %.cweb dummy.ch
	$(CTANGLE) -bhp $< dummy.ch $@

%.h: %.hweb dummy.ch
	$(CTANGLE) -bhp $< dummy.ch $@

if HAVE_CWEAVE
if HAVE_PDFTEX
if HAVE_EPLAIN
pdf-local: integ.pdf

integ.pdf: main.web $(CWEBSRC)
	$(CWEAVE) -bhp main.web
	$(PDFTEX) main
	mv main.pdf integ.pdf
endif
endif
endif

CLEANFILES = integ.pdf main.idx main.log main.scn main.tex main.toc
