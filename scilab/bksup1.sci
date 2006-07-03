function [d]=bksup1(ny,jcf)
d=[];
// Copyright (C) 2001 Michel Juillard
// 

global c

ir = (iter_-2)*ny+1:ny+(iter_-2)*ny;
%irf = iyf+(iter_-1)*ny;
icf = 1:size(iyf,2);
 
for i = 2:iter_
  c(ir,jcf) = c(ir,jcf)-c(ir,icf)*c(%irf,jcf);
  ir = ir-ny;
  %irf = %irf-ny;
end
 
d = c(:,jcf);
 
