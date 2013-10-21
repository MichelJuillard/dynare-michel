DYNARE=../../../matlab/dynare_m
all: test1
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
test1 : test1.o example1.o example1_steadystate.o dynare_cpp_driver.o
	gcc -g -o test1 test1.o example1.o example1_steadystate.o dynare_cpp_driver.o -lm -lstdc++
