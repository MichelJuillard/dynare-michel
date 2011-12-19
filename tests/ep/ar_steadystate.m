function [ys, info] = ar_steadystate(ys, exogenous)
% Steady state routine for ar.mod (First order autoregressive process)
    
global M_
    
info = 0;

ys(1)=M_.params(2);
ys(2)=0;