DYNARE=../../../matlab/dynare_m
am_test_dr_OBJECTS = test_dr-Matrix.$(OBJEXT) test_dr-Vector.$(OBJEXT) \
	test_dr-QRDecomposition.$(OBJEXT) \
	test_dr-GeneralizedSchurDecomposition.$(OBJEXT) \
	test_dr-LUSolver.$(OBJEXT) test_dr-DecisionRules.$(OBJEXT) \
	test_dr-test-dr.$(OBJEXT)

all: test1
Matrix.o : ../../../mex/sources/estimation/libmat/Matrix.cc
	gcc -g -c ../../../mex/sources/estimation/libmat/Matrix.cc -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
Vector.o : ../../../mex/sources/estimation/libmat/Vector.cc
	gcc -g -c ../../../mex/sources/estimation/libmat/Vector.cc -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
QRDecomposition.o : ../../../mex/sources/estimation/libmat/QRDecomposition.cc
	gcc -g -c ../../../mex/sources/estimation/libmat/QRDecomposition.cc -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
GeneralizedSchurDecomposition.o : ../../../mex/sources/estimation/libmat/GeneralizedSchurDecomposition.cc
	gcc -g -c ../../../mex/sources/estimation/libmat/GeneralizedSchurDecomposition.cc -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
LUSolver.o : ../../../mex/sources/estimation/libmat/LUSolver.cc
	gcc -g -c ../../../mex/sources/estimation/libmat/LUSolver.cc -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
DecisionRules.o : ../../../mex/sources/estimation/DecisionRules.cc
	gcc -g -c ../../../mex/sources/estimation/DecisionRules.cc -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat

test1.o : test1.cc ../dynare_cpp_driver.hh ../dynare_cpp_driver.cc
	gcc -g -c test1.cc -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
dynare_cpp_driver.o: ../dynare_cpp_driver.cc ../dynare_cpp_driver.hh
	gcc -g -c ../dynare_cpp_driver.cc -I..
example1.cc example1_steadystate.cc: example1.mod
	$(DYNARE) example1.mod output=first
example1.o: example1.cc
	gcc -g -c example1.cc -I..
example1_steadystate.o: example1_steadystate.cc
	gcc -g -c example1_steadystate.cc
example1_first_derivatives.o: example1_first_derivatives.cc
	gcc -g -c example1_first_derivatives.cc

test1 : test1.o example1.o example1_steadystate.o example1_first_derivatives.o dynare_cpp_driver.o Matrix.o Vector.o QRDecomposition.o GeneralizedSchurDecomposition.o LUSolver.o DecisionRules.o
	gcc -g -o test1 test1.o example1.o example1_steadystate.o example1_first_derivatives.o dynare_cpp_driver.o Matrix.o Vector.o QRDecomposition.o GeneralizedSchurDecomposition.o LUSolver.o DecisionRules.o -llapack -lblas -lm -lstdc++

.cc.o:
	gcc -g -c -o $@ $< -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
