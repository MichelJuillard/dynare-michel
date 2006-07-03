function []=steady_()
// Copyright (C) 2001 Michel Juillard
// 
 
global ys_
 
x = ys_;
xlen = xkmin_+xkmax_+1;
nn = size(iy_,2);
it_ = ykmin_+1;
temp = ex_;
ex_ = ones(xlen,1)*exe_';
 
 
[ys_,%check] = solve(fname_+'_fff',x);
 
if %check==0 then
  dyn_disp(' ');
  dyn_disp('STEADY: Stationary state reached.');
else
  dyn_disp('STEADY: convergence problems');
end
 
ex_ = temp;
 
// 06/24/01 MJ: steady_ no results printer; steady with printed results
// 09/26/01 MJ: translated to Scilab, changed function to <logname_>_fff 
 
 
 
 
 
 

