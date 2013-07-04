% --*-- Unitary tests --*--
function C = eq(A,B)

%@info:
%! @deftypefn {Function File} {@var{C} =} eq (@var{A},@var{B})
%! @anchor{@dynSeries/eq}
%! @sp 1
%! Overloads the eq (equal) operator for the @ref{dynSeries} class.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! @ref{dynSeries} object.
%! @item B
%! @ref{dynSeries} object.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item C
%! scalar integer equal to one if a==b, 0 otherwise.
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

if nargin~=2
    error('dynDates::eq: I need exactly two input arguments!')
end

C = isequal(A,B);


%@test:1
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Define names
%$ A_name = {'A1';'A2';'A3'};
%$
%$ t = zeros(2,1);
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    ts2 = ts1;
%$    a = isequal(ts1,ts2);
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dyn_assert(a,1);
%$ end
%$ T = all(t);
%@eof:1