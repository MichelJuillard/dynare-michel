CPP = c++
ifeq ($(DEBUG),yes)
	CPPFLAGS = -mno-cygwin -ggdb -pg -Wall
	FLEXFLAGS = -i -d
else
ifeq ($(MINGW),yes)
	CPPFLAGS = -O2
	FLEXFLAGS = -i
else
	CPPFLAGS = -mno-cygwin -O2
	FLEXFLAGS = -i
endif
endif
ifeq ($(PROFILE),yes)
	CPPFLAGS = -mno-cygwin -O2 -pg
endif

COMMON_OBJ=\
	DynareFlex.o\
	DynareBison.o\
	ComputingTasks.o\
	DynareMain.o\
	Expression.o\
	ModelParameters.o\
	ModelTree.o\
	NumericalConstants.o\
	NumericalInitialization.o\
	OperatorTable.o\
	OutputFile.o\
	Shocks.o\
	SigmaeInitialization.o\
	SymbolTable.o\
	TmpSymbolTable.o\
	VariableTable.o\
	DynareParser.o\
	DataTree.o

MATLAB_OBJ = InterfaceMatlab.o

SCILAB_OBJ = InterfaceScilab.o


################################################################################
### Build ######################################################################
################################################################################

all: dynare.exe dynare_s.exe

dynare.exe: $(COMMON_OBJ) $(MATLAB_OBJ)
	$(CPP) $(CPPFLAGS) -o  dynare.exe $(COMMON_OBJ) $(MATLAB_OBJ);\
	cp dynare.exe ../matlab/

dynare_s.exe: $(COMMON_OBJ) $(SCILAB_OBJ)
	$(CPP) $(CPPFLAGS) -o  dynare_s.exe $(COMMON_OBJ) $(SCILAB_OBJ);\
	cp dynare_s.exe ../scilab/



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

DynareFlex.cc: DynareFlex.ll include/DynareScanner.h
	flex $(FLEXFLAGS) -oDynareFlex.cc DynareFlex.ll

DynareBison.cc include/DynareBison.h: DynareBison.yy include/DynareParser.h
	 (bison  -v -d -b --verbose -o DynareBison.cc DynareBison.yy; mv DynareBison.hh include/DynareBison.h)


################################################################################
### Clean ######################################################################
################################################################################

clean:
	rm -f *.o *.P \
		DynareBison.cc \
		include/DynareBison.h \
		DynareBison.output \
		DynareFlex.cc \
		dynare.exe \
		dynare_s.exe
