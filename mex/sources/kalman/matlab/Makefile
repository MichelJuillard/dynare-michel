# $Id: Makefile 532 2005-11-30 13:51:33Z kamenik $

# Copyright 2005, Ondra Kamenik

CC = gcc
CC_FLAGS = -Wall -I../../sylv/cc/
LDFLAGS = -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g
else
	CC_FLAGS := $(CC_FLAGS) -O3
endif

ifeq ($(OS),Windows_NT)
	CC_FLAGS := -mno-cygwin $(CC_FLAGS)
	LDFLAGS := -mno-cygwin $(LDFLAGS)
	ARCH := w32
	MEX_SUFFIX = dll
else
	LDFLAGS := $(LDFLAGS)
	ARCH := linux
	MEX_SUFFIX = mexglx
endif

matrix_interface := GeneralMatrix Vector SylvException 
matobjs := $(patsubst %, ../../sylv/cc/%.o, $(matrix_interface))
mathsource := $(patsubst %, ../../sylv/cc/%.h, $(matrix_interface))
matcppsource := $(patsubst %, ../../sylv/cc/%.cpp, $(matrix_interface))

tscwebsource := $(wildcard ../cc/*.cweb)
tscppsource := $(patsubst %.cweb,%.cpp,$(tscwebsource))
tshwebsource := $(wildcard ../cc/*.hweb)
tshsource := $(patsubst %.hweb,%.h,$(tshwebsource))
tsobjects := $(patsubst %.cweb,%.o,$(tscwebsource))

cppsource := $(wildcard *.cpp)
mexobjects := $(patsubst %.cpp,%_.$(MEX_SUFFIX),$(cppsource))

all: $(mexobjects)

../cc/dummy.ch:
	make -C ../cc dummy.ch

../cc/%.cpp: ../cc/%.cweb ../cc/dummy.ch
	make -C ../cc $*.cpp

../cc/%.h: ../cc/%.hweb ../cc/dummy.ch
	make -C ../cc $*.h

../cc/%.o: ../cc/%.cpp $(tshsource)
	make -C ../cc $*.o

tslib.a: $(tshwebsource) $(tscwebsoure) $(tshsource) $(tscppsource) \
         $(mathsource) $(matcppsource) \
         $(tsobjects) $(matobjs)
	ar cr tslib.a $(tsobjects) $(matobjs)
	ranlib tslib.a

# to compile mex objects for Windows do:
# 1. install gnumex
# 2. create mexopts.bat via gnumex in this directory, specify MinGW compilation, and dll output
# 3. change -Ic:/MATLAB7/extern/include according your Matlab setup
# 4. pray it works
%_.$(MEX_SUFFIX): %.cpp $(tshwebsource) $(tscwebsoure) $(tshsource) $(tscppsource) \
         $(mathsource) $(matcppsource) \
         tslib.a
ifeq ($(OS),Windows_NT)
	mex.bat -I../../sylv/cc/ -I../cc COMPFLAGS\#"-c -DMATLAB_MEX_FILE" OPTIMFLAGS\#"-O3 -fexceptions -Ic:/MATLAB7/extern/include" GM_ADD_LIBS\#"$(LDFLAGS)" $*.cpp tslib.a
else
	mex -cxx -I../../sylv/cc/ -I../cc $*.cpp CFLAGS="$(CC_FLAGS) -fexceptions" tslib.a $(LDFLAGS)
endif
	mv $*.$(MEX_SUFFIX) $*_.$(MEX_SUFFIX)

clear:
	rm -f tslib.a
	rm -f *.mexglx
	rm -f *.dll
	make -C ../testing clear
	make -C ../cc clear
