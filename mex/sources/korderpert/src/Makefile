include ../Makefile.include
all: k_orddbgtest
#all: k_order_perturbation.mexa64

k_order_test_main.o: k_order_test_main.cpp
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I$(MATLAB_PATH)/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -c k_order_test_main.cpp 

k_ord_dynare.o: k_ord_dynare.cpp k_ord_dynare.h
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I$(MATLAB_PATH)/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -fPIC -c k_ord_dynare.cpp

k_order_perturbation.o: k_order_perturbation.cpp k_order_perturbation.h	
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I$(MATLAB_PATH)/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -fPIC -c k_order_perturbation.cpp

nlsolve.o: $(DYNAREPPDIR)/src/nlsolve.cpp
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I$(MATLAB_PATH)/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -fPIC -c k_order_test_main.cpp -c k_ord_dynare.cpp -c k_order_perturbation.cpp -c $(DYNAREPPDIR)/src/nlsolve.cpp


fs2000a_dynamic.o: fs2000a_dynamic.c
	gcc -g -c  -I$(MATLAB_PATH)/extern/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "fs2000a_dynamic.c"

mexversion.o: $(MATLAB_PATH)/extern/src/mexversion.c
	gcc -c  -I$(MATLAB_PATH)/extern/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "$(MATLAB_PATH)/extern/src/mexversion.c"

fs2000a_dynamic.so: fs2000a_dynamic.o mexversion.o dynamic.map
#	gcc -O -pthread -shared -Wl,--version-script,dynamic.map -Wl,--no-undefined -o  "fs2000a_dynamic.so"  fs2000a_dynamic.o mexversion.o  -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 -L$(MATLAB_PATH)/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
	gcc -O -pthread -shared -Wl,--no-undefined -o  "fs2000a_dynamic.so"  fs2000a_dynamic.o mexversion.o  -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 -L$(MATLAB_PATH)/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++

fs2000k_dynamic.o: fs2000k_dynamic.c
	gcc -g -c  -I$(MATLAB_PATH)/extern/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "fs2000k_dynamic.c"

fs2000k_dynamic.so: fs2000k_dynamic.o mexversion.o dynamic.map
#	gcc -O -pthread -shared -Wl,--version-script,dynamic.map -Wl,--no-undefined -o  "fs2000k_dynamic.so"  fs2000k_dynamic.o mexversion.o  -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 -L$(MATLAB_PATH)/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
	gcc -O -pthread -shared -Wl,--no-undefined -o  "fs2000k_dynamic.so"  fs2000k_dynamic.o mexversion.o  -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 -L$(MATLAB_PATH)/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++

ramst1_dynamic.o: ramst1_dynamic.c
	gcc -g -c  -I$(MATLAB_PATH)/extern/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "ramst1_dynamic.c"

ramst1_dynamic.mexa64: ramst1_dynamic.o mexversion.o dynamic.map
	gcc -O -pthread -shared -Wl,--no-undefined -o  "ramst1_dynamic.mexa64"  ramst1_dynamic.o mexversion.o  -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 -L$(MATLAB_PATH)/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++


k_orddbgtest: k_order_test_main.o k_ord_dynare.o k_order_perturbation.o nlsolve.o  dynarelib.a
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS  -o k_orddbgtest k_order_test_main.o k_ord_dynare.o k_order_perturbation.o nlsolve.o dynarelib.a -Wl,-rpath,$(MATLAB_PATH)/bin/glnxa64 -Wl,-rpath,$(MATLAB_PATH)/sys/os/glnxa64 -Wl,-rpath,. -L$(MATLAB_PATH)/bin/glnxa64 -lmex -lmx -lmwlapack -lmwblas -lmat -lm -lstdc++

k_orddbgtest_ramst1: k_order_test_main_ramst1.o k_ord_dynare.o k_order_perturbation.o nlsolve.o  dynarelib.a
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS  -o k_orddbgtest_ramst1 k_order_test_main_ramst1.o k_ord_dynare.o k_order_perturbation.o nlsolve.o dynarelib.a -Wl,-rpath,$(MATLAB_PATH)/bin/glnxa64 -Wl,-rpath,$(MATLAB_PATH)/sys/os/glnxa64 -Wl,-rpath,. -L$(MATLAB_PATH)/bin/glnxa64 -lmex -lmx -lmwlapack -lmwblas -lmat -lm -lstdc++

k_order_perturbation.mexa64: k_order_perturbation.o k_ord_dynare.o nlsolve.o dynarelib.a
	$(MEX) -v -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I$(MATLAB_PATH)/extern/include -I$(DYNAREPPDIR)/src -I$(DYNAREPPDIR)/kord -I$(DYNAREPPDIR)/tl/cc -I$(DYNAREPPDIR)/utils/cc/ -I$(DYNAREPPDIR)/sylv/cc/ k_order_perturbation.cpp k_ord_dynare.o nlsolve.o dynarelib.a -lmwlapack -lmwblas -ldl -lstdc++ LDFLAGS='-pthread -shared -Wl,--version-script,$(MATLAB_PATH)/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -Wl,-rpath,.'