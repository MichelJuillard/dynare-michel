function sp = setFreq(sp,freq)
% Set frequency of a dynTime object.

%@info:
%! @deftypefn {Function File} {@var{sp} =} setFreq (@var{sp}, @var{freq})
%! @anchor{@dynTime/setFreq}
%! @sp 1
%! Set frequency of a dynTime object.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item sp
%! dynTime object instantiated by @ref{dynTime}
%! @item freq
%! scalar integer equal to 1 (yearly data), 4 (quaterly data), 12 (monthly data) or 52 (weekly data).
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
%! Let @var{sp} be an object instantiated by @ref{dynTime}, both following syntaxes can be used to define the frequency:
%! @sp 1
%! @example
%! sp = setFreq(sp,4);
%! @end example
%! or
%! @example
%! sp = sp.setFreq(4);
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

sp.freq = freq;