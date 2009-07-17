# $Id: Makefile 531 2005-11-30 13:49:48Z kamenik $
# Copyright 2005, Ondra Kamenik

#DEBUG = yes

#LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c
#-mno-cygwin -DWINDOWS  
CC_FLAGS :=  -DMATLAB -DNO_BLAS_H -DNO_LAPACK_H \
	-Wall -I../cc -I../../sylv/cc -I../cc \
	-I$(MATLAB_PATH)/extern/include  #-pg

ifeq ($(DEBUG),yes)
	CC_FLAGS := -DDEBUG $(CC_FLAGS) -g -pg
#	CC_FLAGS := -DTIMING_LOOP -DDEBUG $(CC_FLAGS) -g  
	KALMANLIB := kalmanlib_dbg.a
else
#	CC_FLAGS := $(CC_FLAGS) -O3
	CC_FLAGS := -DTIMING_LOOP $(CC_FLAGS) -O3
	KALMANLIB := kalmanlib.a
endif

# Added by GP
#    LDFLAGS := -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++ -lmingw32
        #LDFLAGS := -Wl,--library-path  $(LD_LIBRARY_PATH) 
        #      -Wl,-L'f:/MinGW/lib' 
        #-Wl,-L'C:/MinGW/lib/gcc-lib/i686-pc-mingw32/4.0.4' -Wl,-L'C:/MinGW/lib'

        LD_LIBS := -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 \
        -Wl,-L$(MATLAB_PATH)/bin/glnxa64 \
        -Wl,-lmex -lmx -lmwlapack -lmwblas -lmat -lm \
        -Wl,-lstdc++ -lgfortran $(LDFLAGS)

#                  -Wl,-L'f:/CygWin/usr/local/atlas/lib' 
#                  -Wl,-L'f:/CygWin/lib' 
#                 $(LDFLAGS)
#    LD_LIBS :=$(LDFLAGS)
# end add

#matrix_interface := GeneralMatrix Vector SylvException 
#matobjs := $(patsubst %, ../sylv/cc/%.o, $(matrix_interface))
#mathsource := $(patsubst %,  ../sylv/cc/%.h, $(matrix_interface))
#matcppsource := $(patsubst %,  ../sylv/cc/%.cpp, $(matrix_interface))
qtf90source := $(wildcard ../f90/*.f90)
qtobjs := $(patsubst %.f90,%.o,$(qtf90source))
cppsource := $(wildcard *.cpp)
hsource := $(wildcard *.h)
objects := $(patsubst %.cpp,%.o,$(cppsource))

dummy.ch:
	touch dummy.ch


%.o: %.cpp $(hsource) $(cppsource)
	c++ $(CC_FLAGS) -c $*.cpp

dgemmmtm_exe.exe: dgemmmtm_exe.o $(objects) $(qtobjs) $(hsource) $(cppsource) 
	gcc $(CC_FLAGS) -o dgemmmtm_exe.exe dgemmmtm_exe.o ascii_array.o \
	$(qtobjs) $(LD_LIBS) 

qtmmmtm_exe.exe: qtmmmtm_exe.o $(objects) $(qtobjs) $(hsource) $(cppsource) 
	gcc $(CC_FLAGS) -o qtmmmtm_exe.exe qtmmmtm_exe.o ascii_array.o \
	$(qtobjs) $(LD_LIBS) 

dgemvm_exe.exe: dgemvm_exe.o $(qtobjs) $(hsource) $(cppsource) 
	gcc $(CC_FLAGS) -o dgemvm_exe.exe dgemvm_exe.o ascii_array.o \
	$(qtobjs) $(LD_LIBS) 

qtamvm_exe.exe: qtamvm_exe.o $(qtobjs) $(hsource) $(cppsource) 
	gcc $(CC_FLAGS) -o qtamvm_exe.exe qtamvm_exe.o ascii_array.o \
	$(qtobjs) $(LD_LIBS) 

qtvmvm_exe.exe: qtvmvm_exe.o $(qtobjs) $(hsource) $(cppsource) 
	gcc $(CC_FLAGS) -o qtvmvm_exe.exe qtvmvm_exe.o ascii_array.o \
	$(qtobjs) $(LD_LIBS) 

qtv1mvm_exe.exe: qtv1mvm_exe.o $(qtobjs) $(hsource) $(cppsource) 
	gcc $(CC_FLAGS) -o qtv1mvm_exe.exe qtv1mvm_exe.o ascii_array.o \
	$(qtobjs) $(LD_LIBS) 

all: $(objects) dgemmmtm_exe.exe qtmmmtm_exe.exe qtv1mvm_exe.exe dgemvm_exe.exe qtvmvm_exe.exe qtamvm_exe.exe # $(cppsource) $(hsource) $(kalmanhsource) $(kalmancppsource)

clear:
	rm -f *.o
	rm -f *.a
	rm -f *.{exe,dll}
