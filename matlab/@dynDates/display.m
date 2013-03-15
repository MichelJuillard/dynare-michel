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
first_displayed = 2;
    
fprintf('%s = <dynDates: ', inputname(1));
    
if dd.ndat<=max_displayed
    for i=1:ndat
        fprintf(format(dynDate(dd.time(i,:),dd.freq)))
        if i<ndat
            fprintf(', ')
        else
            fprintf('>\n')
        end
    end
else
    for i=1:first_displayed
        fprintf(format(dynDate(dd.time(i,:),dd.freq)))
        fprintf(', ')
    end
    fprintf(' ..., ')
    fprintf(format(dynDate(dd.time(dd.ndat-1,:),dd.freq)))
    fprintf(', ')
    fprintf(format(dynDate(dd.time(dd.ndat,:),dd.freq)))
    fprintf('>\n')
end