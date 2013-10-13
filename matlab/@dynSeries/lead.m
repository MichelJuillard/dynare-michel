function us = lead(ts,p) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{us} =} lead (@var{ts})
%! @anchor{lag}
%! @sp 1
%! Computes leaded time series.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @var
%! @item ts
%! Dynare time series object, instantiated by @ref{dynSeries}
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item us
%! Dynare time series object with transformed data field.
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

% Set default number of leads
if nargin<2
    p = 1;
end

if p<=0
    error('dynSeries::lead: Second input argument must be strictly positive! Use lag method instead.')
end

% Copy of ts dynSeries object
us = ts;

% Update data member
us.data = [  ts.data(p+1:end,:); NaN(p,ts.vobs);];

for i=1:ts.vobs
    us.name(i) = {[ 'lead(' ts.name{i} ',' int2str(p) ')']};
    us.tex(i) = {[ ts.tex{i} '_{+' int2str(p) '}']};
end

%@test:1
%$ t = zeros(4,1);
%$
%$ try
%$     data = transpose(0:1:50);
%$     ts = dynSeries(data,'1950Q1');
%$     a = ts.lead;
%$     b = ts.lead.lead;
%$     c = lead(ts,2);
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     DATA = [transpose(1:50); NaN(1,ts.vobs)];
%$     t(2) = dyn_assert(a.data,DATA,1e-15);
%$     DATA = [transpose(2:50); NaN(2,ts.vobs)];
%$     t(3) = dyn_assert(b.data,DATA,1e-15);
%$     t(4) = dyn_assert(b.data,c.data,1e-15);
%$ end
%$
%$ T = all(t);
%@eof:1