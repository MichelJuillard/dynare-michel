function [m,s,p1,p2] = uniform_specification(m,s,p3,p4)

% function [m,s,p1,p2] = uniform_specification(m,s,p3,p4)
% Specification of the uniform density function parameters
%
% INPUTS
%    m:      mean
%    s:      standard deviation 
%    p3:     lower bound 
%    p4:     upper bound 

% OUTPUTS
%    m:      mean
%    s:      standard deviation 
%    p1:     lower bound 
%    p2:     upper bound 
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2004-2008 Dynare Team
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

    if ~(isnan(p3) | isnan(p4))
      p1 = p3;
      p2 = p4;
      m = (p3+p4)/2;
      s = (p4-p3)/(sqrt(12));
    else
      p1 = m-s*sqrt(3);
      p2 = m+s*sqrt(3);
    end
