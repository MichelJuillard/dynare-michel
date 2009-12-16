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

% Copyright (C) 2003-2007 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

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

