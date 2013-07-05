function info = is_unitary_test_available(fun)

%@info:
%! @deftypefn {Function File} {@var{info} =} is_unitary_test_available (@var{fun})
%! @anchor{is_unitary_test_available}
%! @sp 1
%! Tests if matlab/octave routine @var{fun} has unitary tests.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item fun
%! string, name of the matlab/octave routine to be tested.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item info
%! Integer scalar equal to one if unitary tests are available, zero otherwise.
%! @end table
%! @sp 2
%! @end deftypefn
%@eod:

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
    
info = 0;

fid = fopen(fun,'r');
first_line = fgetl(fid);

if strcmp(first_line,'% --*-- Unitary tests --*--')
    info = 1;
end

fclose(fid);