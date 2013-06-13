function b = isempty(a)

%@info:
%! @deftypefn {Function File} {@var{b} =} isempty (@var{a})
%! @anchor{@dynDate/isempty}
%! @sp 1
%! Overloads the isempty function for the @ref{dynDate} class.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Dynare date object instantiated by @ref{dynDate}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item b
%! Integer scalar (equal to zero if @var{a} is not empty).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2012-2013 Dynare Team
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

if ~isa(a,'dynDate')
    error(['dynDate::isempty: Input argument ' inputname(1) ' have to be a dynDate object!'])
end

b = all(isnan(a.time)) && isnan(a.freq);

%@test:1
%$ % Instantiate an empty dynDate object
%$ d = dynDate();
%$ % Test if this object is empty
%$ t(1) = isempty(d);
%$ T = all(t);
%@eof:1
