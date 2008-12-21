k_order_perturbation project status

Last upate: 21/12/08

Makefiles are not complete and in final verison yet 
but they and C++ files are ready for creation of a library 
that can be linked with the set of k_order_perturbation 
Dynare extensions:

k_order_perturbation.cpp (and .h)
k_ord_dynare.cpp (and .h)

NOTE: at the moment the comp.bat compiles also:
k_order_test_main.cpp - exe test driver main function
nlsolv.cpp - i.e. not covered by the lib make files.

to create debug, POSIX THREAD and Matlab Lapack verion of k_ord_dynarelib.a
k_ord_dynarelibML_PTRD_DB.a, 

1) Compile and create dynarelib.a with MingW and MATLAB LAPACK:
go to Dynare_pp/extern/matlab and run make

NOTE: At the momment the library still uses integ files - 
those will be removed once makefiles are comopleted

2) move/copy dynarelib.a to root src and rename it to 
k_ord_dynarelibML_PTRD_DB.a

Then:
3) Compile extensions with MingW using the temporary comp.bat:
go to the root src directory and run comp.bat 


4) to link statically linked test exe driver k_orddbgtest.exe run 
linkdbgexe.bat
 
5) to link Matlab callable DLL 
linkDLL.bat


NOTE: The dll is called from Matlab Dynare dr1, after set_state_space as:

[ysteady, gx, gu]=k_ord_dynare_perturbation(dr,task,M_,options_, oo_ )

