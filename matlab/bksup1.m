% Copyright (C) 2001 Michel Juillard
%
function d = bksup1(ny,jcf)

global options_ iyf c 

ir = [(options_.periods-2)*ny+1:ny+(options_.periods-2)*ny] ;
irf = iyf+(options_.periods-1)*ny ;
icf = [1:size(iyf,2)] ;

for i = 2:options_.periods
	c(ir,jcf) = c(ir,jcf)-c(ir,icf)*c(irf,jcf) ;
	ir = ir-ny ;
	irf = irf-ny ;
end

d = c(:,jcf) ;

return ;
