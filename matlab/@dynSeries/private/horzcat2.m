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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

if ~(isa(b,'dynSeries') && isa(c,'dynSeries'))
    error('dynSeries::horzcat: All input arguments have to be Dynare time series objects!')
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
if ~isequal(b.Time(1),c.Time(1))
    d_init_flag = 1;
end

d_last_flag = 0;
if ~isequal(b.Time(end),c.Time(end))
    d_last_flag = 1;
end

a.vobs = b.vobs+c.vobs;
a.name = char(b.name,c.name);
a.tex  = char(b.tex,c.tex);

if ~( d_nobs_flag(1) || d_init_flag(1) || d_last_flag(1) )
    a.Time = b.Time;
    a.data = [b.data,c.data];
else
    [JUNK,IB] = setdiff(b.Time(:),c.Time(:),'rows');
    [JUNK,IC] = setdiff(c.Time(:),b.Time(:),'rows');
    [JUNK,JB,JC] = intersect(b.Time(:),c.Time(:),'rows');
    a.Time = a.Time.setTime(sortrows([b.Time(IB); b.Time(JB); c.Time(IC)],[1 2]));
    a.nobs = rows(a.Time(:));
    a.data = NaN(a.nobs,a.vobs);
    [junk,ia,ib] = intersect(a.Time(:),b.Time(:),'rows');
    a.data(ia,1:b.vobs) = b.data;
    [junk,ia,ic] = intersect(a.Time(:),c.Time(:),'rows');
    a.data(ia,b.vobs+(1:c.vobs)) = c.data;
end