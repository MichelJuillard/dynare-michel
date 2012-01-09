function [hasLicense] = user_has_matlab_license(toolbox)
%[hasLicense] = user_has_matlab_license(toolbox)
% checks for license using the appropriate function call
%
% INPUTS
%   toolbox: string for toolbox name
%
% OUTPUTS
%   hasLicense: bool indicating whether or not the user has the license
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2012 Dynare Team
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

if matlab_ver_less_than('7.12')
    hasLicense = license('test', toolbox);
else
    [hasLicense ~] = license('checkout',toolbox);
end
end
