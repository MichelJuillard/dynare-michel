# $Id: Makefile 534 2005-11-30 13:58:11Z kamenik $
# Copyright 2005, Ondra Kamenik

#LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c
#CC_FLAGS := -Wall -I../cc -I../../sylv/cc

CC_FLAGS := -Wall -I../cc -I../sylv/cc


DEBUG = yes
MATLAB = 1

# Added by GP
#    LDFLAGS := -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++ -lmingw32
        LDFLAGS := -Wl,--library-path $(LD_LIBRARY_PATH) \
             -Wl,-L'f:/MinGW/lib' \
	-Wl,-L"c:/Program Files"/MATLAB_SV71/extern/lib/win32/microsoft/ \
	-Wl,-llibmex -Wl,-llibmx  -Wl,-llibmwlapack -Wl,-llibdflapack \
	-lg2c -lmingw32  -lstdc++ $(LDFLAGS)

#                  -Wl,-L'f:/CygWin/usr/local/atlas/lib' \
#                  -Wl,-L'f:/CygWin/lib' \
#                 $(LDFLAGS)

    LD_LIBS :=$(LDFLAGS)
# end add

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g
else
	CC_FLAGS := $(CC_FLAGS) -O2
endif

matrix_interface := GeneralMatrix Vector SylvException 
matobjs := $(patsubst %, ../sylv/cc/%.o, $(matrix_interface))
mathsource := $(patsubst %, ../sylv/cc/%.h, $(matrix_interface))
#cwebsource := $(wildcard ../cc/*.cweb)
cppsource := $(wildcard ../cc/*.cpp)
#cppsource := $(patsubst %.cweb,%.cpp,$(cwebsource)) 
mainobjects := $(patsubst %.cpp,%.o,$(cppsource))
#objects := $(patsubst %.cpp,%.o, *.cpp)
objects := $(wildcard *.o)
#hwebsource := $(wildcard ../cc/*.hweb)
hsource := $(wildcard ../cc/*.h)
#hsource := $(patsubst %.hweb,%.h,$(hwebsource))

../cc/dummy.ch:
	make -C ../cc dummy.ch

../cc/%.o: ../cc/%.cpp $(hsource) $(mathsource)
	make -C ../cc $*.o

%.o: %.cpp $(hsource) ascii_matrix.h
	c++ $(CC_FLAGS) -c $*.cpp

#tests: $(hmainsource) $(cppmainsource) $(hsource) $(cppsource) objects) \

tests: $(hsource) $(cppsource) $(mainobjects) \
       tests.o ascii_matrix.o
	c++ $(CC_FLAGS) $(mainobjects) $(matobjs) ascii_matrix.o tests.o -o tests $(LD_LIBS)

clear:
	make -C ../cc clear
	rm -f *.o
	rm -f tests
