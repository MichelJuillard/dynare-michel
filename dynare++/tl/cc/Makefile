# $Id: Makefile 2344 2009-02-09 20:36:08Z michel $
# Copyright 2004, Ondra Kamenik

include ../../Makefile.include

#LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c

CC_FLAGS := -I../../sylv/cc $(CC_FLAGS)

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g -DTL_DEBUG=2
else
	CC_FLAGS := $(CC_FLAGS) -O2 -DPOSIX_THREADS
endif

ifeq ($(OS),Windows_NT)
	CC_FLAGS := -mno-cygwin -mthreads $(CC_FLAGS)
	LD_LIBS := -mno-cygwin -mthreads $(LD_LIBS)  -lpthreadGC1
else
	LD_LIBS := $(LD_LIBS) -lpthread
endif


matrix_interface := GeneralMatrix Vector SylvException 
matobjs := $(patsubst %, ../../sylv/cc/%.o, $(matrix_interface))
cwebsource := $(wildcard *.cweb)
cppsource := $(patsubst %.cweb,%.cpp,$(cwebsource)) 
objects := $(patsubst %.cweb,%.o,$(cwebsource))
hwebsource := $(wildcard *.hweb)
hsource := $(patsubst %.hweb,%.h,$(hwebsource))


dummy.ch:
	touch dummy.ch

%.cpp: %.cweb dummy.ch
	ctangle -bhp $*.cweb dummy.ch $*.cpp

%.h: %.hweb dummy.ch
	ctangle -bhp $*.hweb dummy.ch $*.h

%.o : %.cpp $(hsource)
	$(CC) $(CC_FLAGS) $(EXTERN_DEFS) -c $*.cpp

all: $(objects) $(cppsource) $(hsource)

tl.pdf: doc

doc: main.web $(hwebsource) $(cwebsource)
	cweave -bhp main.web
	pdftex main
	mv main.pdf tl.pdf

clear:
	rm -f $(cppsource)
	rm -f $(hsource)
	rm -f *.o
	rm -f main.{dvi,idx,log,pdf,scn,tex,toc}
	rm -f dummy.ch
	rm -f *~
