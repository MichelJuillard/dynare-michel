# $Id: Makefile 843 2006-07-28 08:54:19Z tamas $
# Copyright 2004, Ondra Kamenik


LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c -lpthread
CC_FLAGS := -Wall -I../cc -I../../sylv/cc \
        $(CC_INCLUDE_PATH)
ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g -DTL_DEBUG=2
else
	CC_FLAGS := $(CC_FLAGS) -O2 -DPOSIX_THREADS
endif

matrix_interface := GeneralMatrix Vector SylvException 
matobjs := $(patsubst %, ../../sylv/cc/%.o, $(matrix_interface))
cwebsource := $(wildcard ../cc/*.cweb)
cppsource := $(patsubst %.cweb,%.cpp,$(cwebsource)) 
objects := $(patsubst %.cweb,%.o,$(cwebsource))
hwebsource := $(wildcard ../cc/*.hweb)
hsource := $(patsubst %.hweb,%.h,$(hwebsource))

../cc/dummy.ch:
	make -C ../cc dummy.ch

../cc/%.cpp: ../cc/%.cweb ../cc/dummy.ch
	make -C ../cc $*.cpp

../cc/%.h: ../cc/%.hweb ../cc/dummy.ch
	make -C ../cc $*.h

../cc/%.o: ../cc/%.cpp $(hsource)
	make -C ../cc $*.o

%.o: %.cpp factory.h monoms.h $(hwebsource) $(hsource)
	$(CC) $(CC_FLAGS) -c $*.cpp

tests: $(hwebsource) $(cwebsoure) $(hsource) $(cppsource) \
       tests.o factory.o monoms.o $(objects) 
	$(CC) $(CC_FLAGS) $(objects) $(matobjs) tests.o factory.o monoms.o -o tests $(LD_LIBS) 

clear:
	rm -f *.o
	rm -f tests
	make -C ../cc clear
