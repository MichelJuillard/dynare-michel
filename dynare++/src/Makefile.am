bin_PROGRAMS = dynare++

GENERATED_FILES = dynglob_ll.cc dynglob_tab.cc dynglob_tab.hh

dynare___SOURCES = \
	main.cpp \
	dynare3.cpp \
	dynare_atoms.h \
	dynare_model.h \
	forw_subst_builder.h \
	planner_builder.cpp \
	dynare3.h \
	dynare_exception.h \
	dynare_params.cpp \
	nlsolve.cpp \
	planner_builder.h \
	dynare_atoms.cpp \
	dynare_model.cpp \
	dynare_params.h \
	forw_subst_builder.cpp \
	nlsolve.h \
	$(GENERATED_FILES)

dynare___CPPFLAGS = -I../sylv/cc -I../tl/cc -I../kord -I../integ/cc -I.. -DDYNVERSION=\"$(PACKAGE_VERSION)\" -DPOSIX_THREADS
dynare___LDADD = ../kord/libkord.a ../integ/cc/libinteg.a ../tl/cc/libtl.a ../parser/cc/libparser.a ../utils/cc/libutils.a ../sylv/cc/libsylv.a $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(PTHREAD_LIBS)
dynare___CXXFLAGS = $(PTHREAD_CFLAGS)

BUILT_SOURCES = $(GENERATED_FILES)
EXTRA_DIST = dynglob.lex dynglob.y

dynglob_tab.cc dynglob_tab.hh: dynglob.y
	$(YACC) -d -odynglob_tab.cc dynglob.y

dynglob_ll.cc: dynglob.lex
	$(LEX) -i -odynglob_ll.cc dynglob.lex
