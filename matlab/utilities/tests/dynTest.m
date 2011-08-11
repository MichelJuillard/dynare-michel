function dynTest(fun)

%@info:
%! @deftypefn {Function File} dynTest (@var{fun})
%! @anchor{dynTest}
%! Tests matlab/octave routine @var{fun.m}.
%!
%! @strong{Inputs}
%! @table @ @var
%! @item fun
%! string, name of the matlab/octave routine to be tested.
%! @end table
%!
%! @strong{Outputs}
%! None
%!
%! @strong{This function is called by:}
%! @ref{dynare}, @ref{mroutines}
%!
%! @strong{This function calls:}
%! @ref{mtest}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011 Dynare Team
% stephane DOT adjemian AT univ DASH lemans DOT fr
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

[pathstr, name, ext] = fileparts(which(fun));
if ~( isempty(pathstr) || isempty(name) || isempty(ext) ) && strcmp(ext(2:end),'m')
    check = mtest(name,pathstr);
    if check
        disp(['Succesfull test(s) for ' fun ' routine!'])
    end
else
    disp([fun  'is not a known matlab/octave routine!'])
end