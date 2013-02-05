function display(dd)

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

    max_displayed = 5;

    fprintf('%s = <dynDates: ', inputname(1));
    
    for i = 1:min(max_displayed, dd.ndat)

        fprintf('%d', dd.time(i,1))
        switch dd.freq
          case 1
          case 4
            fprintf('Q')
          case 12
            fprintf('M')
          case 52
            fprintf('W')
          otherwise
            error('Unknown frequency %d', dd.freq)
        end
        if dd.freq ~= 1
            fprintf('%d', dd.time(i,2))
        end
        if i ~= min(dd.ndat, max_displayed)
            fprintf(', ')
        end
    end
    if dd.ndat > max_displayed
        fprintf(', ...')
    end
    fprintf('>\n')
end
