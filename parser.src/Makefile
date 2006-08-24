CPP = c++
ifeq ($(DEBUG),yes)
	CPPFLAGS = -mno-cygwin -ggdb -pg -Wall
	FLEXFLAGS = -i
else
ifeq ($(MINGW),yes)
	CPPFLAGS = -O2
	FLEXFLAGS = -i -d
else
	CPPFLAGS = -mno-cygwin -O2
	FLEXFLAGS = -i
endif
endif
ifeq ($(PROFILE),yes)
	CPPFLAGS = -mno-cygwin -O2 -pg
endif
OBJ=\
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
	DataTree.o\
	Interface.o

################################################################################
### Build ######################################################################
################################################################################

all: dynare.exe

dynare.exe: $(OBJ) 
	$(CPP) $(CPPFLAGS) -o  dynare.exe $(OBJ);\
	cp dynare.exe ../matlab/dynare_m.exe;\

dynare_s.exe: $(OBJ) 
	$(CPP) $(CPPFLAGS) -Iinclude -DSCILAB -o DynareMain.o -c DynareMain.cc
	$(CPP) $(CPPFLAGS) -o  dynare_s.exe $(OBJ);\
	cp dynare_s.exe ../scilab/dynare_s.exe;\



################################################################################
### Compile ####################################################################
################################################################################

%.o : %.cc
	$(CPP) $(CPPFLAGS) -MD -I include -c $<
	@cp $*.d $*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	  rm -f $*.d

-include $(OBJ:.o=.P)

DynareFlex.cc: DynareFlex.ll include/DynareScanner.h
	flex $(FLEXFLAGS) -oDynareFlex.cc DynareFlex.ll

DynareBison.cc include/DynareBison.h: DynareBison.yy include/DynareParser.h
	 (bison  -v -d -b --verbose -o DynareBison.cc DynareBison.yy; mv DynareBison.hh include/DynareBison.h)


################################################################################
### Clean ######################################################################
################################################################################

clean:
	rm "ComputingTasks.o" \
	"DynareParser.o" \
	"DynareBison.o" \
	"DynareFlex.o" \
	"DynareMain.o" \
	"Expression.o" \
	"ModelParameters.o" \
	"ModelTree.o" \
	"NumericalConstants.o" \
	"NumericalInitialization.o" \
	"OperatorTable.o" \
	"OutputFile.o" \
	"Shocks.o" \
	"SigmaeInitialization.o" \
	"SymbolTable.o" \
	"TmpSymbolTable.o" \
	"VariableTable.o" \
	"DataTree.o" \
	"dynare.exe" \
	"DynareBison.cc" \
	"include/DynareBison.h" \
	"DynareFlex.cc"

