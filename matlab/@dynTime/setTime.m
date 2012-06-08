function sp = setTime(sp,i,date)
% Update time member of a dynTime object.

%@info:
%! @deftypefn {Function File} {@var{sp} =} setTime (@var{sp}, @var{i}, @var{freq})
%! @anchor{@dynTime/setTime}
%! @sp 1
%! Update time member of a dynTime object.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item sp
%! dynTime object instantiated by @ref{dynTime}
%! @item i
%! (if nargin==3) scalar integer, targeting a row of the matrix @var{sp}.time (see @ref{dynTime} for a description of this matrix).
%! @sp 1
%! (if nargin==2) nobs*2 matrix of integers (see @ref{dynTime} for a description of this array).
%! @item date
%! row vector of two integers (see @ref{dynTime} for a description of the rows of @var{sp}).
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item sp
%! Updated @ref{dynTime} object.
%! @end table
%! @sp 2
%! @strong{Example}
%! @sp 1
%! Let @var{sp} be an object instantiated by @ref{dynTime}, both following syntaxes can be used to update the time member:
%! @sp 1
%! @example
%! sp = setTime(sp,2,[1950,4]);
%! @end example
%! or
%! @example
%! sp = sp.setFreq(2,[1950,4]);
%! @end example
%! @sp 1
%! Note that the second syntax is probably slower than the first one, and should not be used in a large loop.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%!
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

if nargin==3
    sp.time(i,:) = date;
elseif nargin==2
    if isa(i,'dynTime')
        sp.time=i.time;
    else
        sp.time=i;
    end
end