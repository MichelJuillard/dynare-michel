k_order_perturbation project status

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
