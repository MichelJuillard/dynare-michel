DYNARE=../../../matlab/dynare_m
all: test1
test1.o : test1.cc ../dynare_cpp_driver.hh ../dynare_cpp_driver.cc
	gcc -g -c test1.cc -I.. -I../../../mex/sources -I../../../mex/sources/estimation -I../../../mex/sources/estimation/libmat
dynare_cpp_driver.o: ../dynare_cpp_driver.cc ../dynare_cpp_driver.hh
	gcc -g -c ../dynare_cpp_driver.cc -I..
example1.o: example1.mod
	$(DYNARE) example1.mod output=first
	gcc -g -c example1.cc -I..
test1 : test1.o example1.o dynare_cpp_driver.o
	gcc -g -o test1 test1.o example1.o dynare_cpp_driver.o -lstdc++
