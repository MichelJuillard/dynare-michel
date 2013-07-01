function a = horzcat(varargin)

%@info:
%! @deftypefn {Function file} {@var{a} =} horzcat (@var{b},@var{c}, ...)
%! @anchor{horzcat}
%! @sp 1
%! Method of the dynDate class.
%! @sp 1
%! Concatenate dynDate objects to form a dynDates object. This method overloads the horizontal concatenation operator, so that
%! two (or more) dynDate objects an be concatenated :
%!
%!     a = [b, c, d];
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item b
%! dynDate object, instantiated by @ref{dynDate}.
%! @item c
%! dynDate object, instantiated by @ref{dynDate}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item a
%! dynDates object.
%! @end table
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

a = dynDates(varargin{:});