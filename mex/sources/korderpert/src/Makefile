include ../Makefile.include
all: k_orddbgtest
#all: k_order_perturbation.mexa64

k_order_test_main.o: k_order_test_main.cpp
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I/home/michel/usr/bin/matlab/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -c k_order_test_main.cpp 

k_ord_dynare.o: k_ord_dynare.cpp k_ord_dynare.h
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I/home/michel/usr/bin/matlab/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -fPIC -c k_ord_dynare.cpp

k_order_perturbation.o: k_order_perturbation.cpp k_order_perturbation.h	
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I/home/michel/usr/bin/matlab/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -fPIC -c k_order_perturbation.cpp

nlsolve.o: $(DYNAREPPDIR)/src/nlsolve.cpp
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I/home/michel/usr/bin/matlab/extern/include -I $(DYNAREPPDIR)/src -I $(DYNAREPPDIR)/kord -I $(DYNAREPPDIR)/tl/cc -I $(DYNAREPPDIR)/utils/cc/ -I $(DYNAREPPDIR)/sylv/cc/ -fPIC -c k_order_test_main.cpp -c k_ord_dynare.cpp -c k_order_perturbation.cpp -c $(DYNAREPPDIR)/src/nlsolve.cpp


fs2000a_dynamic.o: fs2000a_dynamic.c
	gcc -g -c  -I/home/michel/usr/bin/matlab/extern/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "fs2000a_dynamic.c"

mexversion.o: /home/michel/usr/bin/matlab/extern/src/mexversion.c
	gcc -c  -I/home/michel/usr/bin/matlab/extern/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "/home/michel/usr/bin/matlab/extern/src/mexversion.c"

fs2000a_dynamic.so: fs2000a_dynamic.o mexversion.o dynamic.map
#	gcc -O -pthread -shared -Wl,--version-script,dynamic.map -Wl,--no-undefined -o  "fs2000a_dynamic.so"  fs2000a_dynamic.o mexversion.o  -Wl,-rpath-link,/home/michel/usr/bin/matlab/bin/glnxa64 -L/home/michel/usr/bin/matlab/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
	gcc -O -pthread -shared -Wl,--no-undefined -o  "fs2000a_dynamic.so"  fs2000a_dynamic.o mexversion.o  -Wl,-rpath-link,/home/michel/usr/bin/matlab/bin/glnxa64 -L/home/michel/usr/bin/matlab/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++

fs2000k_dynamic.o: fs2000k_dynamic.c
	gcc -g -c  -I/home/michel/usr/bin/matlab/extern/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread  -DMX_COMPAT_32 -O -DNDEBUG  "fs2000k_dynamic.c"

fs2000k_dynamic.so: fs2000k_dynamic.o mexversion.o dynamic.map
#	gcc -O -pthread -shared -Wl,--version-script,dynamic.map -Wl,--no-undefined -o  "fs2000k_dynamic.so"  fs2000k_dynamic.o mexversion.o  -Wl,-rpath-link,/home/michel/usr/bin/matlab/bin/glnxa64 -L/home/michel/usr/bin/matlab/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++
	gcc -O -pthread -shared -Wl,--no-undefined -o  "fs2000k_dynamic.so"  fs2000k_dynamic.o mexversion.o  -Wl,-rpath-link,/home/michel/usr/bin/matlab/bin/glnxa64 -L/home/michel/usr/bin/matlab/bin/glnxa64 -lmx -lmex -lmat -lm -lstdc++

k_orddbgtest: k_order_test_main.o k_ord_dynare.o k_order_perturbation.o nlsolve.o  dynarelib.a
	gcc -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS  -o k_orddbgtest k_order_test_main.o k_ord_dynare.o k_order_perturbation.o nlsolve.o dynarelib.a -Wl,-rpath,/home/michel/usr/bin/matlab/bin/glnxa64 -Wl,-rpath,/home/michel/usr/bin/matlab/sys/os/glnxa64 -Wl,-rpath,. -L/home/michel/usr/bin/matlab/bin/glnxa64 -lmex -lmx -lmwlapack -lmwblas -lmat -lm -lstdc++

k_order_perturbation.mexa64: k_order_perturbation.o k_ord_dynare.o nlsolve.o dynarelib.a
	mex -v -DMATLAB_MEX_FILE -g  -DDEBUG -DPOSIX_THREADS -I/home/michel/usr/bin/matlab/extern/include -I$(DYNAREPPDIR)/src -I$(DYNAREPPDIR)/kord -I$(DYNAREPPDIR)/tl/cc -I$(DYNAREPPDIR)/utils/cc/ -I$(DYNAREPPDIR)/sylv/cc/ k_order_perturbation.cpp k_ord_dynare.o nlsolve.o dynarelib.a -lmwlapack -lmwblas -ldl -lstdc++ 