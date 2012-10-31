function a = isfile(b)

%@info:
%! @deftypefn {Function File} {@var{a} =} isfile (@var{b})
%! @anchor{isfile}
%! @sp 1
%! Test if @var{b} is a file.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @var
%! @item b
%! A matlab/octave string.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item a
%! Integer scalar, equal to 1 if @var{b} is a file, zero otherwise.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

[base,ext] = strtok(b,'.');

if isempty(ext)
    % File has no extension.
    [status, c] = fileattrib(b);
    if status
        a = ~c.directory;
    else
        a = 0;
    end
else
    a = isequal(exist(b,'file'),2);
end