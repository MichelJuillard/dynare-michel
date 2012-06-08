function dynTest(fun,dynare_path)

%@info:
%! @deftypefn {Function File} dynTest (@var{fun})
%! @anchor{dynTest}
%! @sp 1
%! Tests matlab/octave routine @var{fun}.m.
%! @sp 2
%!
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item fun
%! string, name of the matlab/octave routine to be tested.
%! @end table
%! @sp 2
%!
%! @strong{Outputs}
%! @sp 1
%! None
%! @sp 2
%!
%! @strong{This function is called by:}
%! @sp 1
%! @ref{internals}, @ref{mroutines}
%! @sp 2
%!
%! @strong{This function calls:}
%! @sp 1
%! @ref{mtest}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2012 Dynare Team
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

original_directory = pwd();

[pathstr1, name1, ext1] = fileparts(fun);

pathstr1 = [original_directory filesep pathstr1];

cd([dynare_path filesep '..' filesep 'tests']);

mex_flag = 0;
if exist(name1)==3
    mex_flag = 1;
end

class_flag = 0;
if ~isempty(strfind(fun,'@')) || ~isempty(strfind(which(name1),'@'))
    class_flag = 1;
end

check = mtest(name1,pathstr1);

if check
    if mex_flag
        disp(['Succesfull test(s) for ' name1 ' mex file!'])
    elseif class_flag
        disp(['Succesfull test(s) for ' name1 ' method!'])
    else
        disp(['Succesfull test(s) for ' name1 ' routine!'])
    end
end

cd(original_directory);