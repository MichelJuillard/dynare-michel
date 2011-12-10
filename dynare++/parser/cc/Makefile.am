noinst_LIBRARIES = libparser.a

GENERATED_FILES = assign_tab.cc csv_tab.cc formula_tab.cc matrix_tab.cc namelist_tab.cc assign_tab.hh csv_tab.hh formula_tab.hh matrix_tab.hh namelist_tab.hh assign_ll.cc csv_ll.cc formula_ll.cc matrix_ll.cc namelist_ll.cc

libparser_a_SOURCES = \
	location.h \
	namelist.h \
	atom_assignings.cpp \
	atom_assignings.h \
	atom_substitutions.cpp \
	atom_substitutions.h \
	csv_parser.cpp \
	csv_parser.h \
	dynamic_atoms.cpp \
	dynamic_atoms.h \
	fine_atoms.cpp \
	fine_atoms.h \
	formula_parser.cpp \
	formula_parser.h \
	matrix_parser.cpp \
	matrix_parser.h \
	parser_exception.cpp \
	parser_exception.h \
	static_atoms.cpp \
	static_atoms.h \
	static_fine_atoms.cpp \
	static_fine_atoms.h \
	tree.cpp \
	tree.h \
	$(GENERATED_FILES)

libparser_a_CPPFLAGS = -I../.. $(BOOST_CPPFLAGS)

BUILT_SOURCES = $(GENERATED_FILES)

EXTRA_DIST = assign.y csv.y formula.y matrix.y namelist.y assign.lex csv.lex formula.lex matrix.lex namelist.lex

%_tab.cc %_tab.hh: %.y
	$(YACC) -d -o$*_tab.cc $<

%_ll.cc: %.lex
	$(LEX) -i -o$@ $<
