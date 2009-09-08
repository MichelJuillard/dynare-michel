# $Id: Makefile 843 2006-07-28 08:54:19Z tamas $
# Copyright 2005, Ondra Kamenik


LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c -lpthread
CC_FLAGS := -Wall -I../cc -I ../../tl/cc -I../../sylv/cc
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
tlcwebsource := $(wildcard ../../tl/cc/*.cweb)
tlcppsource := $(patsubst %.cweb,%.cpp,$(tlcwebsource)) 
tlobjects := $(patsubst %.cweb,%.o,$(tlcwebsource))
tlhwebsource := $(wildcard ../../tl/cc/*.hweb)
tlhsource := $(patsubst %.hweb,%.h,$(tlhwebsource))

../cc/dummy.ch:
	make -C ../cc dummy.ch

../cc/%.cpp: ../cc/%.cweb ../cc/dummy.ch
	make -C ../cc $*.cpp

../cc/%.h: ../cc/%.hweb ../cc/dummy.ch
	make -C ../cc $*.h

../cc/%.o: ../cc/%.cpp $(hsource)
	make -C ../cc $*.o

../../tl/cc/dummy.ch:
	make -C ../../tl/cc dummy.ch

../../tl/cc/%.cpp: ../../tl/cc/%.cweb ../../tl/cc/dummy.ch
	make -C ../../tl/cc $*.cpp

../../tl/cc/%.h: ../../tl/cc/%.hweb ../../tl/cc/dummy.ch
	make -C ../../tl/cc $*.h

../../tl/cc/%.o: ../../tl/cc/%.cpp $(tlhsource)
	make -C ../../tl/cc $*.o

%.o: %.cpp $(hwebsource) $(hsource) $(tlhwebsource) $(tlhsource) 
	$(CC) $(CC_FLAGS) -c $*.cpp

tests: $(hwebsource) $(cwebsoure) $(hsource) $(cppsource) \
       $(tlhwebsource) $(tlcwebsoure) $(tlhsource) $(tlcppsource) \
       tests.o $(objects) $(tlobjects) 
	$(CC) $(CC_FLAGS) $(objects) $(tlobjects) $(matobjs) tests.o -o tests $(LD_LIBS) 

clear:
	rm -f *.o
	rm -f tests
	make -C ../cc clear
	make -C ../../tl/cc clear
