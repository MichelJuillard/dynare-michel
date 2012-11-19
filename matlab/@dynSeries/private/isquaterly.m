function b = isquaterly(date)
    
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

[year, remain] = strtok(date,'Q');
if ~isint(str2num(year))
    b = 0;
    return 
end
[quarter, remain] = strtok(remain,'Q');
if ~isempty(remain)
    b = 0;
    return
end
quarter = str2num(quarter);
if ~isint(quarter) || quarter<1 || quarter>4
    b = 0;
    return
end
b = 1;