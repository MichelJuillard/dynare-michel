function display(d)

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

    fprintf('%s = <dynDate: %d', inputname(1), d.time(1));
    switch d.freq
        case 1
        case 4
            fprintf('Q')
        case 12
            fprintf('M')
        case 52
            fprintf('W')
        otherwise
            error('Unknown frequency %d', d.freq)
    end
    if d.freq ~= 1
        fprintf('%d', d.time(2))
    end
    fprintf('>\n')
end
