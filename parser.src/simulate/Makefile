CPP = g++
CPPFLAGS = -Wall
MATLABPATH = C:/MATLAB7

ifeq ($(shell uname -o), Cygwin)
	# Detection of uninitialized variables is buggy in Cygwin and generates spurious warnings
	CPPFLAGS += -Wno-uninitialized
	CPPFLAGS += -mno-cygwin
endif

ifeq ($(CROSS_WIN32), yes)
	CPP = i586-mingw32msvc-g++
	# Detection of uninitialized variables is buggy in MinGW and generates spurious warnings
	CPPFLAGS += -Wno-uninitialized
endif

ifeq ($(DEBUG),yes)
	CPPFLAGS += -ggdb
else
	CPPFLAGS += -O3 -DMATLAB_MEX_FILE -malign-double -fno-exceptions -mtune=pentium
endif

ifeq ($(VALGRIND), yes)
	CPPFLAGS = -Wall -O -g -fno-inline
endif

################################################################################
### Build ######################################################################
################################################################################

all: simulate.o linbcg.o sparsematrix.o mem_mngr.o interpreter.o
	dlltool --def c:/matlab7/extern/include/libmex.def --output-lib temp1.a
	dlltool --def c:/matlab7/extern/include/libmx.def --output-lib temp2.a
	$(CPP) -mno-cygwin -shared simulate.o linbcg.o sparsematrix.o mem_mngr.o interpreter.o -o ../../matlab/simulate.dll temp1.a temp2.a -lstdc++
	rm -f temp1.a temp2.a simulate.d simulate.P

################################################################################
### Compile ####################################################################
################################################################################
mem_mngr.o : mem_mngr.cc
	$(CPP) $(CPPFLAGS) -MD -I include -I $(MATLABPATH)/extern/include -c mem_mngr.cc
linbcg.o : ../linbcg.cc
	$(CPP) $(CPPFLAGS) -MD -I ../include -I $(MATLABPATH)/extern/include -c ../linbcg.cc
sparsematrix.o : sparsematrix.cc
	$(CPP) $(CPPFLAGS) -MD -I include -I $(MATLABPATH)/extern/include -c sparsematrix.cc
interpreter.o : interpreter.cc
	$(CPP) $(CPPFLAGS) -MD -I include -I $(MATLABPATH)/extern/include -c interpreter.cc
simulate.o : ../simulate.cc
	$(CPP) $(CPPFLAGS) -MD -I ../include -I $(MATLABPATH)/extern/include -c ../simulate.cc
	@cp $*.d $*.P;

