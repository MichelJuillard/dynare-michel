CPP = g++

CPPFLAGS = -Wall

ifeq ($(shell uname -o), Cygwin)
	# Detection of uninitialized variables is buggy in Cygwin and generates spurious warnings
	CPPFLAGS += -Wno-uninitialized
	CPPFLAGS += -mno-cygwin
	DYNARE_M = dynare_m.exe
	DYNARE_S = dynare_s.exe
else 
	DYNARE_M = dynare_m
	DYNARE_S = dynare_s
endif

ifeq ($(CROSS_WIN32), yes)
	CPP = i586-mingw32msvc-g++
	# Detection of uninitialized variables is buggy in MinGW and generates spurious warnings
	CPPFLAGS += -Wno-uninitialized
	DYNARE_M = dynare_m.exe
	DYNARE_S = dynare_s.exe
endif

ifeq ($(DEBUG),yes)
	CPPFLAGS += -ggdb
else
	CPPFLAGS += -O3
endif

ifeq ($(VALGRIND), yes)
	CPPFLAGS = -Wall -O -g -fno-inline
endif


COMMON_OBJ=\
	DynareFlex.o\
	DynareBison.o\
	ComputingTasks.o\
	DynareMain.o\
	ModelTree.o\
	NumericalConstants.o\
	NumericalInitialization.o\
	Shocks.o\
	SigmaeInitialization.o\
	SymbolTable.o\
	TmpSymbolTable.o\
	VariableTable.o\
	ParsingDriver.o\
	DataTree.o \
	ModFile.o \
	Statement.o \
	ExprNode.o \
	ModelNormalization.o \
	ModelBlocks.o \
	BlockTriangular.o \
	Model_Graph.o \
	SymbolGaussElim.o

MATLAB_OBJ = InterfaceMatlab.o

SCILAB_OBJ = InterfaceScilab.o


################################################################################
### Build ######################################################################
################################################################################

all: $(DYNARE_M) $(DYNARE_S)

$(DYNARE_M): $(COMMON_OBJ) $(MATLAB_OBJ)
	$(CPP) $(CPPFLAGS) -o $(DYNARE_M) $(COMMON_OBJ) $(MATLAB_OBJ)
	cp $(DYNARE_M) ../matlab/

$(DYNARE_S): $(COMMON_OBJ) $(SCILAB_OBJ)
	$(CPP) $(CPPFLAGS) -o $(DYNARE_S) $(COMMON_OBJ) $(SCILAB_OBJ)
	cp $(DYNARE_S) ../scilab/


################################################################################
### Compile ####################################################################
################################################################################

%.o : %.cc
	$(CPP) $(CPPFLAGS) -MD -I include -c $<
	@cp $*.d $*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	  rm -f $*.d

-include $(COMMON_OBJ:.o=.P) $(MATLAB_OBJ:.o=.P) $(SCILAB_OBJ:.o=.P)

DynareFlex.cc: DynareFlex.ll include/DynareBison.hh include/ParsingDriver.hh
	flex -oDynareFlex.cc DynareFlex.ll

DynareBison.cc include/DynareBison.hh: DynareBison.yy include/ParsingDriver.hh
	bison --verbose -o DynareBison.cc DynareBison.yy
	mv DynareBison.hh location.hh stack.hh position.hh include/


################################################################################
### Clean ######################################################################
################################################################################

clean:
	rm -f *.o *.P \
		*~ include/*~ \
		DynareBison.output \
		DynareFlex.cc \
		$(DYNARE_M) \
		$(DYNARE_S)

distclean: clean
	rm -f DynareBison.cc \
		include/position.hh \
		include/stack.hh \
		include/location.hh \
		include/DynareBison.hh

