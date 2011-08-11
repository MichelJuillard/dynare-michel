function ts = set_time(ts)

%@info:
%! @deftypefn {Function File} {@var{ts} =} set_time(@var{ts})
%! @anchor{dynSeries}
%! @sp 1
%! Fill the time field of a Dynare time series object instantiated by @ref{dynSeries}.
%! @sp 2
%! @strong{Inputs}
%! @sp 2
%! @table @ @var
%! @item ts
%! Dynare time series object, instantiated by @ref{dynSeries}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item t
%! Dynare time series object.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @table @ @ref
%! @item dynSeries
%! @item horzcat
%! @end table
%! @sp 2
%! @strong{This function calls:}
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
ts.time = zeros(ts.nobs,2);
ts.time(1,:) = ts.init;
for obs=2:ts.nobs
    if ts.time(obs-1,2)<ts.freq
        ts.time(obs,1) = ts.time(obs-1,1);
        ts.time(obs,2) = ts.time(obs-1,2)+1;
    else
        ts.time(obs,1) = ts.time(obs-1,1)+1;
        ts.time(obs,2) = 1;
    end
end