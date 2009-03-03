k_order_perturbation project status

Last upate: 3/3/09

Makefiles are not complete and in final verison yet 
but they and C++ files are ready for creation of a library 
that can be linked with the set of k_order_perturbation 
Dynare extensions:

k_order_perturbation.cpp (and .h)
k_ord_dynare.cpp (and .h)

NOTE: at the moment the comp.bat compiles also:
k_order_test_main.cpp - exe test driver main function
nlsolv.cpp - i.e. not covered by the lib make files.

to create a debug, POSIX THREAD and Matlab Lapack verion of k_ord_dynarelib.a
k_ord_dynarelibML_PTRD_DB.a, 

1) use cygwin shell to compile and create k_ord_lib.a under Windows 
with MingW and MATLAB LAPACK:

go to Dynare_pp/extern/matlab assuming MATLAB is already defined with directory
and in your cygwin shell set 

>WINDOWS=1
>DEBUG=1
>export WINDOWS DEBUG

and run 

> make dynarelib.a

NOTE: At the momment the library still uses integ files too - 
those may be removed at a later stage


Then:
2) Compile extensions with MingW using the temporary comp.bat:
go to the root src directory and run comp.bat 


3) to link statically linked test exe driver k_orddbgtest.exe run 
linkdbgexe.bat
 
4) to link Matlab callable DLL 
linkDLL.bat


NOTE: The dll  (mexw32 or so) is called from new Matlab Dynare function dr1_k_order
derived from dr1, after set_state_space as:

[ysteady, ghx_u]=k_ord_dynare_perturbation(dr,task,M_,options_, oo_, ['.' mexext])

where last term is optional but it will default to .dll on windows and .so on linux.

dr1_k_order is called by amended resol.m:

elseif(options_.use_k_order==1)&& (check_flag == 0)
    [dr,info,M_,options_,oo_] = dr1_k_order(dr,check_flag,M_,options_,oo_);
else

and requirese options to be set

options_.use_k_order=1;

==================
Tests:

first_order.m is matlab emulation of Dynare++ c++ first_order.cpp for testing pruposes


==================
ToDO:
==================
1) amend <model>.m to use Dynamic_mexopts.bat

   mex -f Dynamic_mexopts.bat -O fs2000k_no_both_UR_dynamic.c

or amend preprocessor to make mex to export Dynamic() function as well as mexFunction() as the Dynamic_mexopts.bat does, e.g.:

set LINKFLAGS=/dll /export:Dynamic /export:%ENTRYPOINT% /MAP ....

2) make k_order_perturbation handle models  which have the "both" variables (i.e. variables that appear both as lag and as lead)