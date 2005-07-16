CPP = c++
ifeq ($(DEBUG),yes)
	CPPFLAGS = -mno-cygwin -ggdb -Wall
	CPPFLAGS = -ggdb -Wall
else
	CPPFLAGS = -mno-cygwin -O1
endif
OBJ=\
	ComputingTasks.o\
	DynareBison.o\
	DynareFlex.o\
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


SRC=../DynareFlex.cc ../DynareBison.cc ../include/DynareBison.h


################################################################################
### Build ######################################################################
################################################################################

all: dynare.exe

dynare.exe: $(SRC) $(OBJ) 
	pwd;\
	$(CPP) $(CPPFLAGS) -o  dynare.exe $(OBJ);\
	cp dynare.exe ../../matlab/dynare_m.exe;\



################################################################################
### Compile ####################################################################
################################################################################

ComputingTasks.o: ../ComputingTasks.cc ModelTree.o ../include/ComputingTasks.h
	$(CPP) $(CPPFLAGS) -c   -I../include  -o ComputingTasks.o ../ComputingTasks.cc


DynareFlex.o: ../DynareFlex.cc 
	$(CPP) $(CPPFLAGS) -c   -I../include -o DynareFlex.o ../DynareFlex.cc


DynareMain.o: ../DynareMain.cc ../include/DynareParser.h ../include/DynareScanner.h ../include/OutputFile.h
	$(CPP) $(CPPFLAGS) -c   -I../include -o DynareMain.o ../DynareMain.cc


Expression.o: ../Expression.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o Expression.o ../Expression.cc


ModelParameters.o: ../ModelParameters.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o ModelParameters.o ../ModelParameters.cc


ModelTree.o: ../ModelTree.cc ../include/ModelTree.h ../include/ModelTypes.h ../include/DataTree.h
	$(CPP) $(CPPFLAGS) -c   -I../include -o ModelTree.o ../ModelTree.cc

DataTree.o: ../DataTree.cc ../include/DataTree.h ../include/ModelTypes.h
	$(CPP) $(CPPFLAGS) -c   -I../include -o DataTree.o ../DataTree.cc

NumericalConstants.o: ../NumericalConstants.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o NumericalConstants.o ../NumericalConstants.cc


NumericalInitialization.o: ../NumericalInitialization.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o NumericalInitialization.o ../NumericalInitialization.cc


OperatorTable.o: ../OperatorTable.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o OperatorTable.o ../OperatorTable.cc


OutputFile.o: ../OutputFile.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o OutputFile.o ../OutputFile.cc


Shocks.o: ../Shocks.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o Shocks.o ../Shocks.cc


SigmaeInitialization.o: ../SigmaeInitialization.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o SigmaeInitialization.o ../SigmaeInitialization.cc


SymbolTable.o: ../SymbolTable.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o SymbolTable.o ../SymbolTable.cc


TmpSymbolTable.o: ../TmpSymbolTable.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o TmpSymbolTable.o ../TmpSymbolTable.cc


VariableTable.o: ../VariableTable.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o VariableTable.o ../VariableTable.cc

DynareParser.o: ../include/DynareParser.h ../DynareParser.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o DynareParser.o ../DynareParser.cc

DynareBison.o: ../DynareBison.cc
	$(CPP) $(CPPFLAGS) -c   -I../include -o DynareBison.o ../DynareBison.cc

../DynareFlex.cc: ../DynareFlex.ll ../include/DynareScanner.h
	flex -i -o../DynareFlex.cc ../DynareFlex.ll

../DynareBison.cc ../include/DynareBison.h: ../DynareBison.yy ../include/DynareParser.h
	 (bison  -v -d -b --verbose -o../DynareBison.cc ../DynareBison.yy; mv ../DynareBison.hh ../include/DynareBison.h)


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
	"../DynareBison.cc" \
	"../include/DynareBison.h" \
	"../DynareFlex.cc"

