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
	flex $(FLEXFLAGS) -o../DynareFlex.cc ../DynareFlex.ll

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

# DO NOT DELETE

ComputingTasks.o: ../include/ComputingTasks.h ../include/TmpSymbolTable.h
ComputingTasks.o: ../include/SymbolTable.h ../include/ModelParameters.h
ComputingTasks.o: ../include/SymbolTableTypes.h
DataTree.o: /usr/include/time.h /usr/include/_ansi.h /usr/include/newlib.h
DataTree.o: /usr/include/sys/config.h /usr/include/machine/ieeefp.h
DataTree.o: /usr/include/sys/reent.h /usr/include/sys/_types.h
DataTree.o: /usr/include/sys/lock.h /usr/include/machine/time.h
DataTree.o: /usr/include/sys/types.h /usr/include/machine/_types.h
DataTree.o: /usr/include/limits.h /usr/include/features.h
DataTree.o: /usr/include/sys/cdefs.h /usr/include/sys/features.h
DataTree.o: /usr/include/machine/types.h ../include/DynareBison.h
DataTree.o: ../include/VariableTable.h ../include/SymbolTable.h
DataTree.o: ../include/ModelParameters.h ../include/SymbolTableTypes.h
DataTree.o: ../include/NumericalConstants.h ../include/DataTree.h
DataTree.o: /usr/include/stdio.h /usr/include/sys/stdio.h
DataTree.o: ../include/OperatorTable.h ../include/ModelTypes.h
DynareBison.o: ../include/DynareParser.h ../include/ModelParameters.h
DynareBison.o: ../include/SymbolTable.h ../include/SymbolTableTypes.h
DynareBison.o: ../include/Expression.h ../include/OperatorTable.h
DynareBison.o: ../include/DynareBison.h ../include/NumericalConstants.h
DynareBison.o: ../include/NumericalInitialization.h ../include/ModelTree.h
DynareBison.o: ../include/ModelTypes.h ../include/DataTree.h
DynareBison.o: /usr/include/stdio.h /usr/include/_ansi.h
DynareBison.o: /usr/include/newlib.h /usr/include/sys/config.h
DynareBison.o: /usr/include/machine/ieeefp.h /usr/include/sys/reent.h
DynareBison.o: /usr/include/sys/_types.h /usr/include/sys/lock.h
DynareBison.o: /usr/include/sys/types.h /usr/include/machine/_types.h
DynareBison.o: /usr/include/limits.h /usr/include/features.h
DynareBison.o: /usr/include/sys/cdefs.h /usr/include/sys/features.h
DynareBison.o: /usr/include/machine/types.h /usr/include/sys/stdio.h
DynareBison.o: ../include/VariableTable.h ../include/Shocks.h
DynareBison.o: ../include/SigmaeInitialization.h
DynareBison.o: ../include/ComputingTasks.h ../include/TmpSymbolTable.h
DynareBison.o: ../include/Objects.h ../include/ylmm/basic_parser.hh
DynareBison.o: ../include/ylmm/basic_messenger.hh
DynareBison.o: ../include/ylmm/basic_location.hh
DynareBison.o: ../include/ylmm/basic_scanner.hh
DynareBison.o: ../include/ylmm/basic_buffer.hh ../include/ylmm/yaccmm.hh
DynareFlex.o: /usr/include/stdio.h /usr/include/_ansi.h
DynareFlex.o: /usr/include/newlib.h /usr/include/sys/config.h
DynareFlex.o: /usr/include/machine/ieeefp.h /usr/include/sys/reent.h
DynareFlex.o: /usr/include/sys/_types.h /usr/include/sys/lock.h
DynareFlex.o: /usr/include/sys/types.h /usr/include/machine/_types.h
DynareFlex.o: /usr/include/limits.h /usr/include/features.h
DynareFlex.o: /usr/include/sys/cdefs.h /usr/include/sys/features.h
DynareFlex.o: /usr/include/machine/types.h /usr/include/sys/stdio.h
DynareFlex.o: /usr/include/errno.h /usr/include/sys/errno.h
DynareFlex.o: /usr/include/unistd.h /usr/include/sys/unistd.h
DynareFlex.o: /usr/include/string.h /usr/include/sys/string.h
DynareFlex.o: ../include/DynareScanner.h ../include/DynareBison.h
DynareFlex.o: ../include/ylmm/basic_scanner.hh
DynareFlex.o: ../include/ylmm/basic_buffer.hh ../include/Objects.h
DynareFlex.o: ../include/SymbolTable.h ../include/ModelParameters.h
DynareFlex.o: ../include/SymbolTableTypes.h ../include/ModelTypes.h
DynareFlex.o: ../include/ylmm/lexmm.hh /usr/include/stdlib.h
DynareFlex.o: /usr/include/machine/stdlib.h /usr/include/alloca.h
DynareMain.o: /usr/include/stdio.h /usr/include/_ansi.h
DynareMain.o: /usr/include/newlib.h /usr/include/sys/config.h
DynareMain.o: /usr/include/machine/ieeefp.h /usr/include/sys/reent.h
DynareMain.o: /usr/include/sys/_types.h /usr/include/sys/lock.h
DynareMain.o: /usr/include/sys/types.h /usr/include/machine/_types.h
DynareMain.o: /usr/include/limits.h /usr/include/features.h
DynareMain.o: /usr/include/sys/cdefs.h /usr/include/sys/features.h
DynareMain.o: /usr/include/machine/types.h /usr/include/sys/stdio.h
DynareMain.o: ../include/DynareParser.h ../include/ModelParameters.h
DynareMain.o: ../include/SymbolTable.h ../include/SymbolTableTypes.h
DynareMain.o: ../include/Expression.h ../include/OperatorTable.h
DynareMain.o: ../include/DynareBison.h ../include/NumericalConstants.h
DynareMain.o: ../include/NumericalInitialization.h ../include/ModelTree.h
DynareMain.o: ../include/ModelTypes.h ../include/DataTree.h
DynareMain.o: ../include/VariableTable.h ../include/Shocks.h
DynareMain.o: ../include/SigmaeInitialization.h
DynareMain.o: ../include/ComputingTasks.h ../include/TmpSymbolTable.h
DynareMain.o: ../include/Objects.h ../include/ylmm/basic_parser.hh
DynareMain.o: ../include/ylmm/basic_messenger.hh
DynareMain.o: ../include/ylmm/basic_location.hh
DynareMain.o: ../include/ylmm/basic_scanner.hh
DynareMain.o: ../include/ylmm/basic_buffer.hh ../include/DynareScanner.h
DynareMain.o: ../include/OutputFile.h
DynareParser.o: ../include/ModelParameters.h ../include/SymbolTable.h
DynareParser.o: ../include/SymbolTableTypes.h ../include/Expression.h
DynareParser.o: ../include/OperatorTable.h ../include/DynareBison.h
DynareParser.o: ../include/NumericalConstants.h
DynareParser.o: ../include/NumericalInitialization.h
DynareParser.o: ../include/ModelTree.h ../include/ModelTypes.h
DynareParser.o: ../include/DataTree.h /usr/include/stdio.h
DynareParser.o: /usr/include/_ansi.h /usr/include/newlib.h
DynareParser.o: /usr/include/sys/config.h /usr/include/machine/ieeefp.h
DynareParser.o: /usr/include/sys/reent.h /usr/include/sys/_types.h
DynareParser.o: /usr/include/sys/lock.h /usr/include/sys/types.h
DynareParser.o: /usr/include/machine/_types.h /usr/include/limits.h
DynareParser.o: /usr/include/features.h /usr/include/sys/cdefs.h
DynareParser.o: /usr/include/sys/features.h /usr/include/machine/types.h
DynareParser.o: /usr/include/sys/stdio.h ../include/VariableTable.h
DynareParser.o: ../include/Shocks.h ../include/SigmaeInitialization.h
DynareParser.o: ../include/ComputingTasks.h ../include/TmpSymbolTable.h
DynareParser.o: ../include/DynareParser.h ../include/Objects.h
DynareParser.o: ../include/ylmm/basic_parser.hh
DynareParser.o: ../include/ylmm/basic_messenger.hh
DynareParser.o: ../include/ylmm/basic_location.hh
DynareParser.o: ../include/ylmm/basic_scanner.hh
DynareParser.o: ../include/ylmm/basic_buffer.hh
Expression.o: ../include/Expression.h ../include/SymbolTable.h
Expression.o: ../include/ModelParameters.h ../include/SymbolTableTypes.h
Expression.o: ../include/OperatorTable.h ../include/DynareBison.h
Expression.o: ../include/NumericalConstants.h
ModelParameters.o: ../include/ModelParameters.h /usr/include/limits.h
ModelParameters.o: /usr/include/features.h /usr/include/sys/cdefs.h
ModelParameters.o: /usr/include/sys/features.h
ModelTree.o: /usr/include/stdio.h /usr/include/_ansi.h
ModelTree.o: /usr/include/newlib.h /usr/include/sys/config.h
ModelTree.o: /usr/include/machine/ieeefp.h /usr/include/sys/reent.h
ModelTree.o: /usr/include/sys/_types.h /usr/include/sys/lock.h
ModelTree.o: /usr/include/sys/types.h /usr/include/machine/_types.h
ModelTree.o: /usr/include/limits.h /usr/include/features.h
ModelTree.o: /usr/include/sys/cdefs.h /usr/include/sys/features.h
ModelTree.o: /usr/include/machine/types.h /usr/include/sys/stdio.h
ModelTree.o: /usr/include/time.h /usr/include/machine/time.h
ModelTree.o: ../include/VariableTable.h ../include/SymbolTable.h
ModelTree.o: ../include/ModelParameters.h ../include/SymbolTableTypes.h
ModelTree.o: ../include/d_tab.h ../include/DynareBison.h
ModelTree.o: ../include/NumericalConstants.h ../include/ModelTree.h
ModelTree.o: ../include/OperatorTable.h ../include/ModelTypes.h
ModelTree.o: ../include/DataTree.h
NumericalConstants.o: ../include/NumericalConstants.h
NumericalInitialization.o: ../include/NumericalInitialization.h
NumericalInitialization.o: ../include/SymbolTable.h
NumericalInitialization.o: ../include/ModelParameters.h
NumericalInitialization.o: ../include/SymbolTableTypes.h
OperatorTable.o: ../include/OperatorTable.h ../include/DynareBison.h
OutputFile.o: ../include/OutputFile.h ../include/SymbolTable.h
OutputFile.o: ../include/ModelParameters.h ../include/SymbolTableTypes.h
OutputFile.o: ../include/ModelTree.h ../include/OperatorTable.h
OutputFile.o: ../include/DynareBison.h ../include/NumericalConstants.h
OutputFile.o: ../include/ModelTypes.h ../include/DataTree.h
OutputFile.o: /usr/include/stdio.h /usr/include/_ansi.h
OutputFile.o: /usr/include/newlib.h /usr/include/sys/config.h
OutputFile.o: /usr/include/machine/ieeefp.h /usr/include/sys/reent.h
OutputFile.o: /usr/include/sys/_types.h /usr/include/sys/lock.h
OutputFile.o: /usr/include/sys/types.h /usr/include/machine/_types.h
OutputFile.o: /usr/include/limits.h /usr/include/features.h
OutputFile.o: /usr/include/sys/cdefs.h /usr/include/sys/features.h
OutputFile.o: /usr/include/machine/types.h /usr/include/sys/stdio.h
OutputFile.o: ../include/VariableTable.h
Shocks.o: ../include/Shocks.h ../include/ModelParameters.h
SigmaeInitialization.o: ../include/SigmaeInitialization.h
SymbolTable.o: ../include/SymbolTable.h ../include/ModelParameters.h
SymbolTable.o: ../include/SymbolTableTypes.h
TmpSymbolTable.o: ../include/SymbolTable.h ../include/ModelParameters.h
TmpSymbolTable.o: ../include/SymbolTableTypes.h
TmpSymbolTable.o: ../include/TmpSymbolTable.h
VariableTable.o: ../include/VariableTable.h ../include/SymbolTable.h
VariableTable.o: ../include/ModelParameters.h
VariableTable.o: ../include/SymbolTableTypes.h
