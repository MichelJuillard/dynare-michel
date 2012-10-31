function B = subsref(A,S)

%@info:
%! @deftypefn {Function File} {@var{us} =} subsref (@var{ts},S)
%! @anchor{@dynDate/subsref}
%! @sp 1
%! Overloads the subsref method for the Dynare dates class (@ref{dynDate}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare date object instantiated by @ref{dynDate}.
%! @item S
%! Matlab's structure array S with two fields, type and subs. The type field is string containing '()', '@{@}', or '.', where '()' specifies
%! integer subscripts, '@{@}' specifies cell array subscripts, and '.' specifies subscripted structure fields. The subs field is a cell array
%! or a string containing the actual subscripts (see matlab's documentation).
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item B
%! A matlab object (public member of the @ref{dynDate} object).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp2
%!
%! @end deftypefn
%@eod:

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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

B = builtin('subsref', A, S);