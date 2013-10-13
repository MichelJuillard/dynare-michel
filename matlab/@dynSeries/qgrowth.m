function us = qgrowth(ts) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{us} =} qgrowth (@var{ts})
%! @anchor{qgrowth}
%! @sp 1
%! Computes quaterly growth rates.
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

us = ts;

switch ts.freq
  case 1
    error('dynSeries::qgrowth: I cannot compute quaterly growth rates from yearly data!')
  case 4
    us.data(2:end,:) = ts.data(2:end,:)./ts.data(1:end-1,:) - 1;
    us.data(1,:) = NaN;
    for i = 1:ts.vobs
        us.name(i) = {['QGROWTH_' us.name{i}]};
        us.tex(i) = {['\delta ' us.tex{i}]};
    end
  case 12
    us.data(4:end,:) = ts.data(4:end,:)./ts.data(1:end-3,:) - 1;
    us.data(1:3,:) = NaN;
    for i = 1:ts.vobs
        us.name(i) = {['QGROWTH_' us.name{i}]};
        us.tex(i) = {['\delta_3 ' us.tex{i}]};
    end
  case 52
    error('dynSeries::qgrowth: I do not know yet how to compute quaterly growth rates from weekly data!')
  otherwise
    error(['dynSeries::ygrowth: object ' inputname(1) ' has unknown frequency']);
end

%@test:1
%$ t = zeros(2,1);
%$
%$ try
%$     data = (1+.01).^transpose(0:1:50);
%$     ts = dynSeries(data,'1950Q1');
%$     ts = ts.qgrowth;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     DATA = NaN(1,ts.vobs);
%$     DATA = [DATA; .01*ones(ts.nobs-1,ts.vobs)];
%$     t(2) = dyn_assert(ts.data,DATA,1e-15);
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(2,1);
%$
%$ try
%$     data = (1+.01).^transpose(0:1:80);
%$     ts = dynSeries(data,'1950M1');
%$     ts = ts.qgrowth;
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     DATA = NaN(3,ts.vobs);
%$     DATA = [DATA; (1.01^3-1)*ones(ts.nobs-3,ts.vobs)];
%$     t(2) = dyn_assert(ts.data,DATA,1e-15);
%$ end
%$
%$ T = all(t);
%@eof:2
