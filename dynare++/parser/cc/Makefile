# Copyright (C) 2005, Ondra Kamenik

# $Id: Makefile 1211 2007-03-19 21:43:42Z kamenik $

CC_FLAGS := -Wall -I../..

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g
else
	CC_FLAGS := $(CC_FLAGS) -O3
endif

objects := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
headers := $(wildcard *.h)
source  := $(wildcard *.cpp)
objects := $(objects) formula_tab.o formula_ll.o matrix_tab.o matrix_ll.o \
           assign_tab.o assign_ll.o namelist_tab.o namelist_ll.o \
           csv_tab.o csv_ll.o
headers := $(headers) formula_tab.hh matrix_tab.hh assign_tab.hh namelist_tab.hh csv_tab.hh
source  := $(source) formula_tab.cc formula_ll.cc formula.y formula.lex \
           matrix_tab.cc matrix_ll.cc matrix.y matrix.lex \
           assign_tab.cc assign_ll.cc assign.y assign.lex \
           namelist_tab.cc namelist_ll.cc namelist.y namelist.lex \
           csv_tab.cc csv_ll.cc csv.y csv.lex

formula_tab.cc: formula.y tree.h formula_parser.h
	bison  -d -t --verbose -oformula_tab.cc formula.y

formula_tab.hh: formula.y tree.h formula_parser.h
	bison  -d -t --verbose -oformula_tab.cc formula.y

formula_ll.cc: formula.lex formula_tab.hh
	flex -i -oformula_ll.cc formula.lex

formula_ll.o: formula_ll.cc $(headers)
	$(CC) -O3 -c formula_ll.cc

formula_tab.o: formula_tab.cc $(headers)
	$(CC) -O3 -c formula_tab.cc

matrix_tab.cc: matrix.y tree.h matrix_parser.h
	bison  -d -t --verbose -omatrix_tab.cc matrix.y

matrix_tab.hh: matrix.y tree.h matrix_parser.h
	bison  -d -t --verbose -omatrix_tab.cc matrix.y

matrix_ll.cc: matrix.lex matrix_tab.hh
	flex -i -omatrix_ll.cc matrix.lex

matrix_ll.o: matrix_ll.cc $(headers)
	$(CC) -O3 -c matrix_ll.cc

matrix_tab.o: matrix_tab.cc $(headers)
	$(CC) -O3 -c matrix_tab.cc

assign_tab.cc: assign.y static_atoms.h atom_assignings.h
	bison -d -t --verbose -oassign_tab.cc assign.y

assign_tab.hh: assign.y static_atoms.h atom_assignings.h
	bison -d -t --verbose -oassign_tab.cc assign.y

assign_ll.cc: assign.lex assign_tab.hh
	flex -i -oassign_ll.cc assign.lex

assign_ll.o: assign_ll.cc $(headers)
	$(CC) -O3 -c assign_ll.cc

assign_tab.o: assign_tab.cc $(hsource)
	$(CC) -O3 -c assign_tab.cc

namelist_tab.cc: namelist.y namelist.h
	bison -d -t --verbose -onamelist_tab.cc namelist.y

namelist_tab.hh: namelist.y namelist.h
	bison -d -t --verbose -onamelist_tab.cc namelist.y

namelist_ll.cc: namelist.lex namelist_tab.hh
	flex -i -onamelist_ll.cc namelist.lex

namelist_ll.o: namelist_ll.cc $(headers)
	$(CC) -O3 -c namelist_ll.cc

namelist_tab.o: namelist_tab.cc $(hsource)
	$(CC) -O3 -c namelist_tab.cc

csv_tab.cc: csv.y csv_parser.h
	bison -d -t --verbose -ocsv_tab.cc csv.y

csv_tab.hh: csv.y csv_parser.h
	bison -d -t --verbose -ocsv_tab.cc csv.y

csv_ll.cc: csv.lex csv_tab.hh
	flex -i -ocsv_ll.cc csv.lex

csv_ll.o: csv_ll.cc $(headers)
	$(CC) -O3 -c csv_ll.cc

csv_tab.o: csv_tab.cc $(hsource)
	$(CC) -O3 -c csv_tab.cc

%.o: %.cpp $(headers)
	$(CC) $(CC_FLAGS) -c $*.cpp

parser.a: $(objects) $(source) $(headers)
	ar cr parser.a $(objects) 

clear:
	rm -f *~
	rm -f *.o
	rm -f *.cc *.hh
	rm -f parser.a
	rm -rf html/
	rm -rf latex/
