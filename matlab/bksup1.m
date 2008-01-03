function d = bksup1(ny,jcf)

% function d = bksup1(ny,jcf)
% Solves deterministic models recursively by backsubstitution for one lead/lag
%
% INPUTS
%    ny:             number of endogenous variables
%    jcf:            variables index forward
%    
% OUTPUTS
%    d:              vector of backsubstitution results
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2007)
% Gnu Public License.


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

