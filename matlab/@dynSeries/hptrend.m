function ts = hptrend(ts, lambda) % --*-- Unitary tests --*--

% ts = hptrend(ts, lambda)
%
% Extracts the trend component form a dynSeries object using Hodrick Prescott filter.
%
% INPUTS 
%  o ts         dynSeries object.
%  o lambda     positive scalar, trend smoothness parameter.
%
% OUTPUTS 
%  o ts         dynSeries object, with time series replaced by the trend component of the original time series.

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

if nargin>1 
    if lambda<=0
        error(['dynSeries::hptrend: Lambda must be a positive integer!'])
    end
else
    lambda = [];
end

for i=1:ts.vobs
    ts.name(i) = {['hptrend(' ts.name{i} ')']};
    ts.tex(i) = {['\text{hptrend}(' ts.tex{i} ')']};
end

ts.data = sample_hp_filter(ts.data,lambda);

%@test:1
%$ plot_flag = 0; 
%$
%$ % Create a dataset.
%$ e = .2*randn(200,1);
%$ u = randn(200,1);
%$ stochastic_trend = cumsum(e); 
%$ deterministic_trend = .1*transpose(1:200);
%$ x = zeros(200,1);
%$ for i=2:200
%$    x(i) = .75*x(i-1) + e(i);
%$ end
%$ y = x + stochastic_trend + deterministic_trend;
%$
%$ % Test the routine.
%$ try
%$     ts = dynSeries(y,'1950Q1');
%$     ts = ts.hptrend(1600);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1950, 1]);
%$     t(5) = dyn_assert(ts.vobs,1);
%$     t(6) = dyn_assert(ts.nobs,200);
%$ end
%$
%$ % Show results
%$ if plot_flag
%$     plot(stochastic_trend+deterministic_trend,'-k');
%$     hold on
%$     plot(ts.data,'--r');
%$     hold off
%$     axis tight
%$ end
%$ T = all(t);
%@eof:1