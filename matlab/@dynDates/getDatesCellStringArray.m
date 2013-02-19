function m = getDatesCellStringArray(dd)
%function m = getDatesCellStringArray(dd)
% Returns a cell array of strings containing the dates
%
% INPUTS
%   dd - dynDates object
%
% OUTPUTS
%   m  - cell array of strings
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013 Dynare Team
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

dateSep = '';
switch dd.freq
    case 1
    case 4
        dateSep = 'q';
    case 12
        dateSep = 'm';
    case 52
        dateSep = 'w';
    otherwise
        error('Unknown frequency %d', dd.freq);
end

m = cell(0);
for i = 1:dd.ndat
    if isempty(dateSep)
        newdate = num2str(dd.time(i,1));
    else
        newdate = [num2str(dd.time(i,1)) dateSep num2str(dd.time(i,2))];
    end
    m = { m{:} newdate };
end
end