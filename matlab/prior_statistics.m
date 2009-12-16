function pstat = prior_statistics(info,M_,bayestopt_,options_)
% This function computes various statistics associated with the prior distribution.
%
% INPUTS 
%   info        [integer]    1*p vector.  
%   M_          [structure]  Model description.
%   bayestopt_  [structure]  Prior distribution description.  
%   options_    [structure]  Global options of Dynare.
%    
% OUTPUTS:
%   pstat       [structure]  Various statistics 
%
% SPECIAL REQUIREMENTS
%   none 

% Copyright (C) 2009 Dynare Team
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

if ~info(0)
    results = prior_sampler(1,M_,bayestopt_,options_);
    results.prior.mass
end