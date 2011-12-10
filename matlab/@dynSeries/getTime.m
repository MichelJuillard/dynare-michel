function time = getTime(ts)

% Copyright (C) 2011 Dynare Team
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

if ts.freq==1
    time = ts.Time(:,1);
    return
end

switch ts.freq
  case 4
    time = [num2str(ts.Time(1,1)) 'Q' num2str(ts.Time(1,2))];
    for i=2:ts.nobs
        time = char(time,[num2str(ts.Time(i,1)) 'Q' num2str(ts.Time(i,2))]);
    end
  case 12
    time = [num2str(ts.Time(1,1)) 'M' num2str(ts.Time(1,2))];
    for i=2:ts.nobs
        time = char(time,[num2str(ts.Time(i,1)) 'M' num2str(ts.Time(i,2))]);
    end
  case 52
    time = [num2str(ts.Time(1,1)) 'W' num2str(ts.Time(1,2))];
    for i=2:ts.nobs
        time = char(time,[num2str(ts.Time(i,1)) 'W' num2str(ts.Time(i,2))]);
    end
  otherwise
    error('dynSeries::getTime: Unknown type of frequency!')
end