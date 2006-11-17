CPP = c++
ifeq ($(DEBUG),yes)
	CPPFLAGS = -ggdb -pg -Wall
else
ifeq ($(MINGW),yes)
	CPPFLAGS = -O2
else
	CPPFLAGS = -mno-cygwin -O2
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
	ParsingDriver.o\
	DataTree.o

MATLAB_OBJ = InterfaceMatlab.o

SCILAB_OBJ = InterfaceScilab.o


################################################################################
### Build ######################################################################
################################################################################

all: dynare_m.exe dynare_s.exe

dynare_m.exe: $(COMMON_OBJ) $(MATLAB_OBJ)
	$(CPP) $(CPPFLAGS) -o  dynare_m.exe $(COMMON_OBJ) $(MATLAB_OBJ)
	cp dynare_m.exe ../matlab/

dynare_s.exe: $(COMMON_OBJ) $(SCILAB_OBJ)
	$(CPP) $(CPPFLAGS) -o  dynare_s.exe $(COMMON_OBJ) $(SCILAB_OBJ)
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
		DynareBison.cc \
		include/DynareBison.hh \
		DynareBison.output \
		DynareFlex.cc \
		dynare_m.exe \
		dynare_s.exe
