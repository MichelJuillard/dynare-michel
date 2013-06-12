function a = horzcat2(b,c)

%@info:
%! @deftypefn {Function file} {@var{a} =} horzcat2 (@var{b},@var{c}, ...)
%! @anchor{private/horzcat2}
%! @sp 1
%! Private method of the dynSeries class.
%! @sp 1
%! Merge two Dynare time series objects. This method overloads the horizontal concatenation operator, so that
%! two time series objects can be merged using the following syntax
%!
%!     a = [b, c];
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item b
%! Dynare time series object, instantiated by @ref{dynSeries}.
%! @item c
%! Dynare time series object, instantiated by @ref{dynSeries}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item a
%! Dynare time series object.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @ref{descriptive_statistics}
%!
%! @strong{This function calls:}
%! @ref{dynSeries}
%!
%! @strong{Remark 1.} It is assumed that the two time series objects have the same frequencies. The two time series objects can cover
%! different time ranges.
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2013 Dynare Team
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

[n,message] = common_strings_in_cell_arrays(b.name,c.name);
if n
    error(['dynSeries::horzcat: I cannot concatenate dynSeries objects with common variable names (' message ')!'])
end
    
if b.freq ~= c.freq
    error('dynSeries::horzcat: All time series objects must have common frequency!')
else
    a = dynSeries();
    a.freq = b.freq;
end

d_nobs_flag = 0;
if b.nobs ~= c.nobs
    d_nobs_flag = 1;
else
    a.nobs = b.nobs;
end

d_init_flag = 0;
if ~isequal(b.init,c.init)
    d_init_flag = 1;
end

a.vobs = b.vobs+c.vobs;
a.name = vertcat(b.name,c.name);
a.tex  = vertcat(b.tex,c.tex);

if ~( d_nobs_flag(1) || d_init_flag(1) )
    a.init = b.init;
    a.data = [b.data,c.data];
    a.time = b.time;
else
    if b.init<=c.init
        a.init = b.init;
        if b.init<c.init
            c.data = [NaN(c.init-b.init,c.vobs); c.data];
        end
    else
        a.init = c.init;
        b_first_lines = b.init-c.init;
        b.data = [NaN(b.init-c.init,b.vobs); b.data];
    end
    b_last_date = b.init+b.nobs;
    c_last_date = c.init+c.nobs;
    if b_last_date<c_last_date
        b.data = [b.data; NaN(c_last_date-b_last_date,b.vobs)];
    elseif b_last_date>c_last_date
        c.data = [c.data; NaN(b_last_date-c_last_date,c.vobs)];
    end
    a.data = [b.data, c.data];
    a.time = unique(b.time.append(c.time));
end
a.nobs = size(a.data,1);
