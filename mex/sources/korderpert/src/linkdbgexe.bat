gcc -DMATLAB_MEX_FILE -mthreads -g -DWINDOWS -DPOSIX_THREADS -I"c:/Program Files"/MATLAB_SV71/extern/include -I Dynare_pp/src -I Dynare_pp/kord -I Dynare_pp/tl/cc -I Dynare_pp/utils/cc/ -I Dynare_pp/sylv/cc/ -I"f:/Pthreads/Pre-built.2/include" -I"f:/mingw/include" -o k_orddbgtest.exe k_order_test_main.o k_ord_dynare.o k_order_perturbation.o nlsolve.o "dynare_pp/extern/matlab"/dynarelib.a -Wl,-L"c:/Program Files"/MATLAB_SV71/extern/lib/win32/microsoft/ -Wl,-llibmex -Wl,-llibmx -Wl,-llibmwlapack -Wl,-llibdflapack -lg2c -lmingw32  -lstdc++ -L"f:/Pthreads/Pre-built.2/lib" -lpthreadGC2
