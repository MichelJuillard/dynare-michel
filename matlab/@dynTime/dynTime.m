function sp = dynTime(a)
%@info:
%! @deftypefn {Function File} {@var{sp} =} dynTime (@var{a})
%! @anchor{dynTime}
%! @sp 1
%! Constructor for the Dynare dynTime class.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! dynTime object instantiated by @ref{dynTime}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item sp
%! dynTime object.
%! @end table
%! @sp 2
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item time
%! Array of integers (nobs*2). The first column defines the years associated to each observation. The second column,
%! depending on the frequency, indicates the week, month or quarter numbers. For yearly data or unspecified frequency
%! the second column is filled by ones.
%! @end table
%! @sp 1
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
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

sp = struct;

sp.freq = [];
sp.time = [];

sp = class(sp,'dynTime');

switch nargin
  case 0
    return
  case 1
    if isa(a,'dynTime')
        sp = a;
    else
        error(['dynTime::dynTime: Input argument ' inputname(1) ' must be a dynTime object!'])
    end
  otherwise
    error(['dynTime::dynTime: Too many input arguments!'])
end