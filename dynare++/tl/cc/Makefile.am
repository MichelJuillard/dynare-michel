CWEBSRC = \
	normal_moments.cweb \
	int_sequence.cweb \
	tensor.cweb \
	ps_tensor.cweb \
	pyramid_prod2.cweb \
	equivalence.cweb \
	fine_container.cweb \
	kron_prod.cweb \
	ps_tensor.hweb \
	t_polynomial.cweb \
	symmetry.cweb \
	stack_container.cweb \
	sthread.hweb \
	twod_matrix.hweb \
	twod_matrix.cweb \
	symmetry.hweb \
	sparse_tensor.cweb \
	fine_container.hweb \
	sthread.cweb \
	int_sequence.hweb \
	tl_exception.hweb \
	pyramid_prod2.hweb \
	t_container.hweb \
	permutation.hweb \
	tensor.hweb \
	gs_tensor.cweb \
	rfs_tensor.hweb \
	pyramid_prod.hweb \
	t_polynomial.hweb \
	pyramid_prod.cweb \
	fs_tensor.cweb \
	sparse_tensor.hweb \
	permutation.cweb \
	equivalence.hweb \
	gs_tensor.hweb \
	normal_moments.hweb \
	tl_static.hweb \
	kron_prod.hweb \
	fs_tensor.hweb \
	stack_container.hweb \
	rfs_tensor.cweb \
	t_container.cweb \
	tl_static.cweb

GENERATED_FILES = \
	normal_moments.cpp \
	int_sequence.cpp \
	tensor.cpp \
	ps_tensor.cpp \
	pyramid_prod2.cpp \
	equivalence.cpp \
	fine_container.cpp \
	kron_prod.cpp \
	ps_tensor.h \
	t_polynomial.cpp \
	symmetry.cpp \
	stack_container.cpp \
	sthread.h \
	twod_matrix.h \
	twod_matrix.cpp \
	symmetry.h \
	sparse_tensor.cpp \
	fine_container.h \
	sthread.cpp \
	int_sequence.h \
	tl_exception.h \
	pyramid_prod2.h \
	t_container.h \
	permutation.h \
	tensor.h \
	gs_tensor.cpp \
	rfs_tensor.h \
	pyramid_prod.h \
	t_polynomial.h \
	pyramid_prod.cpp \
	fs_tensor.cpp \
	sparse_tensor.h \
	permutation.cpp \
	equivalence.h \
	gs_tensor.h \
	normal_moments.h \
	tl_static.h \
	kron_prod.h \
	fs_tensor.h \
	stack_container.h \
	rfs_tensor.cpp \
	t_container.cpp \
	tl_static.cpp

noinst_LIBRARIES = libtl.a

libtl_a_SOURCES = $(CWEBSRC) $(GENERATED_FILES)
libtl_a_CPPFLAGS = -I../../sylv/cc $(CPPFLAGS_MATIO)
libtl_a_CXXFLAGS = $(PTHREAD_CFLAGS)

BUILT_SOURCES = $(GENERATED_FILES)

EXTRA_DIST = main.web dummy.ch

%.cpp: %.cweb dummy.ch
	$(CTANGLE) -bhp $< dummy.ch $@

%.h: %.hweb dummy.ch
	$(CTANGLE) -bhp $< dummy.ch $@

if HAVE_CWEAVE
if HAVE_PDFTEX
if HAVE_EPLAIN
pdf-local: tl.pdf

tl.pdf: main.web $(CWEBSRC)
	$(CWEAVE) -bhp main.web
	$(PDFTEX) main
	mv main.pdf tl.pdf
endif
endif
endif

CLEANFILES = tl.pdf main.idx main.log main.scn main.tex main.toc
