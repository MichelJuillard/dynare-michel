# $Id: Makefile 531 2005-11-30 13:49:48Z kamenik $
# Copyright 2005, Ondra Kamenik

DEBUG = yes

LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c
CC_FLAGS :=  -DMATLAB -DWINDOWS  -DNO_BLAS_H -DNO_LAPACK_H \
	-Wall -I../../gensylv/cc -I../cc \
	-Ic:/"Program Files"/MATLAB_SV71/extern/include 

ifeq ($(DEBUG),yes)
	CC_FLAGS := -DDEBUG $(CC_FLAGS) -g
else
	CC_FLAGS := $(CC_FLAGS) -O2
endif

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

matrix_interface := GeneralMatrix Vector SylvException 
matobjs := $(patsubst %, ../../gensylv/cc/%.o, $(matrix_interface))
mathsource := $(patsubst %,  ../../gensylv/cc/%.h, $(matrix_interface))
matcppsource := $(patsubst %,  ../../gensylv/cc/%.cpp, $(matrix_interface))
# cppsource := $(patsubst %.cweb,%.cpp,$(cwebsource)) 
kalmancppsource := $(wildcard ../cc/*.cpp)
kalmanhsource := $(wildcard ../cc/*.h)
kalmanobjects := $(patsubst %.cpp,%.o,$(kalmancppsource))
cppsource := $(wildcard *.cpp)
hsource := $(wildcard *.h)
objects := $(patsubst %.cpp,%.o,$(cppsource))
hwebsource := $(wildcard *.hweb)
cwebsource := $(wildcard *.cweb)

dummy.ch:
	touch dummy.ch

# %.cpp: %.cweb dummy.ch
# 	ctangle -bhp $*.cweb dummy.ch $*.cpp

# %.h: %.hweb dummy.ch
# 	ctangle -bhp $*.hweb dummy.ch $*.h

#%.o: %.cpp $(hsource) $(cppsource) #$(kalmanhsource) $(mathsource)

%.o: %.cpp $(hsource) $(cppsource)
	c++ $(CC_FLAGS) -c $*.cpp


kalmanlib.a: $(objects) #$(kalmanobjects) $(matobjs) 
	ar cr kalmanlib.a $(objects) $(kalmanobjects) $(matobjs) 
	ranlib kalmanlib.a

kalman_smoother_dll.dll: kalman_smoother.o kalmanlib.a #$(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o kalman_smoother_dll.dll kalman_smoother.o \
	kalmanlib.a $(LD_LIBS) 

kalman_filter_dll.dll: kalman_filters.o kalmanlib.a # $(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o kalman_filter_dll.dll kalman_filters.o \
	kalmanlib.a $(LD_LIBS) 

all: $(objects) kalmanlib.a kalman_smoother_dll.dll kalman_filter_dll.dll # $(cppsource) $(hsource) $(kalmanhsource) $(kalmancppsource)

doc: main.web $(hwebsource) $(cwebsource)
	cweave -bhp main.web
	pdftex main
	mv main.pdf ts.pdf 

clear:
	rm -f *.o
	rm -f *.{pdf,dvi,log,scn,idx,toc}
#	rm -f *.cpp
#	rm -f *.h
