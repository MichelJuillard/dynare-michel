function b = isyearly(date)

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

year = str2num(date);
if ~isempty(year)
    if isint(year)
        b = 1;
    else
        b = 0;
    end
    return
else
    [year, remain] = strtok(date,'Y');
    if ~isequal(remain,'Y')
        b = 0;
        return
    end
    if ~isint(str2num(year))
        b = 0;
        return
    end
    b = 2;
end