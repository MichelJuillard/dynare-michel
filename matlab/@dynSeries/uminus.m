function A = uminus(B)
%@info:
%! @deftypefn {Function File} {@var{A} =} plus (@var{B},@var{C})
%! @anchor{@dynSeries/uminus}
%! @sp 1
%! Overloads the uminus method for the Dynare time series class (@ref{dynSeries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item B
%! Dynare time series object instantiated by @ref{dynSeries}.
%! @item C
%! Dynare time series object instantiated by @ref{dynSeries}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare time series object.
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

A = dynSeries();

A.freq = B.freq;
A.nobs = B.nobs;
A.vobs = B.vobs;
A.init = B.init;
A.time = B.time;
A.name = cell(A.vobs,1);
A.tex = cell(A.vobs,1);
for i = 1:A.vobs
    A.name(i) = {[ '-' B.name{i}]};
    A.tex(i) = {[ '-' B.tex{i}]};
end
A.data = -(B.data);

%@test:1
%$ % Define a datasets.
%$ A = rand(10,2);
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ A_tex = {'A_1';'A_2'};
%$ t = zeros(6,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,A_tex);
%$    ts2 = -ts1;
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dyn_assert(ts2.vobs,2);
%$    t(3) = dyn_assert(ts2.nobs,10);
%$    t(4) = dyn_assert(ts2.data,-A,1e-15);
%$    t(5) = dyn_assert(ts2.name,{'-A1';'-A2'});
%$    t(6) = dyn_assert(ts2.tex,{'-A_1';'-A_2'});
%$ end
%$ T = all(t);
%@eof:1
