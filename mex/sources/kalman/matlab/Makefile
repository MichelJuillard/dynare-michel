# $Id: Makefile 531 2005-11-30 13:49:48Z kamenik $
# Copyright 2005, Ondra Kamenik


#DEBUG = yes

#LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c

#CC_FLAGS :=  -DWINDOWS -mno-cygwin 

CC_FLAGS :=  -DMATLAB -DNO_BLAS_H -DNO_LAPACK_H \
	-Wall -fpic -I../qt/cc -I../sylv/cc -I../cc \
	-I$(MATLAB_PATH)/extern/include  #-pg

ifeq ($(DEBUG),yes)
	CC_FLAGS := -DDEBUG $(CC_FLAGS) -g
#	CC_FLAGS := -DTIMING_LOOP -DDEBUG $(CC_FLAGS) -g #-pg #-Wl,-pg
	KALMANLIB := kalmanlib_dbg.a
else
#	CC_FLAGS := $(CC_FLAGS) -O3
	CC_FLAGS := -DTIMING_LOOP $(CC_FLAGS) -O3
	KALMANLIB := kalmanlib.a
endif

# Added by GP
#    LDFLAGS := -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++ -lmingw32
        #LDFLAG := -Wl,--library-path  $(LD_LIBRARY_PATH) 
        #-Wl,-L'C:/MinGW/lib/gcc-lib/i686-pc-mingw32/4.0.4' 
        #-Wl,-llibmex -Wl,-llibmx -Wl,-llibmwlapack -Wl,-llibdflapack -lf95 
	#-lg2c -lmingw32  kalmanlib.def -Wl,-lmwm_ir 
                  #-Wl,-L'f:/CygWin/lib' 
        #LD_LIBS := -Wl,--library-path
        #          -Wl,-L'/usr/lib' 

        LD_LIBS := -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 \
	-Wl,-L$(MATLAB_PATH)/bin/glnxa64 \
        -Wl,-lmex -lmx -lmwlapack -lmwblas -lmat -lm \
        -Wl,-lstdc++ $(LDFLAGS) 

	#-Wl,-L'/usr/lib'  

#                  -Wl,-L'f:/CygWin/usr/local/atlas/lib' 
#                  -Wl,-L'f:/CygWin/lib' 
#              -Wl,-L'f:/MinGW/lib' 
#                 $(LDFLAGS)
#    LD_LIBS :=$(LDFLAGS)
# end add

matrix_interface := GeneralMatrix Vector SylvException 
matobjs := $(patsubst %, ../sylv/cc/%.o, $(matrix_interface))
mathsource := $(patsubst %,  ../sylv/cc/%.h, $(matrix_interface))
matcppsource := $(patsubst %,  ../sylv/cc/%.cpp, $(matrix_interface))
qtf90source := $(wildcard ../qt/f90/*.f90)
qtobjs := $(patsubst %.f90,%.o,$(qtf90source))
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


$(KALMANLIB): $(objects) # $(matobjs) $(qtobjs) #$(kalmanobjects)
	ar cr $(KALMANLIB) $(kalmanobjects) $(matobjs) $(qtobjs)  
	ranlib $(KALMANLIB)

kalman_smoother_dll.dll: kalman_smoother.o $(KALMANLIB) #$(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o kalman_smoother_dll.dll kalman_smoother.o \
	kalmanlib.a $(LD_LIBS) 

minv.dll: minv.o $(KALMANLIB) # $(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o minv.dll minv.o \
	kalmanlib.a $(LD_LIBS) 

gmvm.dll: gmvm.o $(KALMANLIB) # $(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o gmvm.dll gmvm.o \
	kalmanlib.a $(LD_LIBS) 

qtamvm.dll: qtamvm.o $(KALMANLIB) # $(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o qtamvm.dll qtamvm.o \
	kalmanlib.a $(LD_LIBS) 

qtmvm.dll: qtmvm.o $(KALMANLIB) # $(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o qtmvm.dll qtmvm.o \
	kalmanlib.a $(LD_LIBS) 

disclyap_fast_dll.dll: disclyap_fast_dll.o $(KALMANLIB) # $(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o disclyap_fast_dll.dll disclyap_fast_dll.o \
	$(KALMANLIB) $(LD_LIBS) kalmanlib.def

kalman_filter_dll.dll: kalman_filters.o $(KALMANLIB) # $(hsource) $(cppsource) 
	gcc -shared $(CC_FLAGS) -o kalman_filter_dll.dll kalman_filters.o \
	$(KALMANLIB) $(LD_LIBS)

kalman_filters_testx.exe: kalman_filters_testx.o $(KALMANLIB) # $(hsource) $(cppsource) 
	gcc $(CC_FLAGS) -pg -o kalman_filters_testx.exe kalman_filters_testx.o \
	$(KALMANLIB) $(LD_LIBS) 

all: $(objects) $(KALMANLIB) kalman_smoother_dll.dll kalman_filter_dll.dll # $(cppsource) $(hsource) $(kalmanhsource) $(kalmancppsource)

#kalman_filter_loop.o: kalman_filters.cpp
#    c++ -DTIMING_LOOP $(CC_FLAGS) -o kalman_filter_loop.o kalman_filters.cpp

#kalman_filter_loop.dll: kalman_filter_loop.o kalmanlib.a # $(hsource) $(cppsource) 
#    gcc -shared -DTIMING_LOOP $(CC_FLAGS) -o kalman_filter_loop.dll kalman_filter_loop.o \
#    kalmanlib.a $(LD_LIBS) 

doc: main.web $(hwebsource) $(cwebsource)
	cweave -bhp main.web
	pdftex main
	mv main.pdf ts.pdf 

clear:
	rm -f *.o
	rm -f *.{pdf,dvi,log,scn,idx,toc}
#	rm -f *.cpp
#	rm -f *.h
