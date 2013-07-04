% --*-- Unitary tests --*--
function a = vertcat(varargin)

%@info:
%! @deftypefn {Function file} {@var{a} =} vertcat (@var{b},@var{c}, ...)
%! @anchor{horzcat}
%! @sp 1
%! Method of the dynSeries class.
%! @sp 1
%! This method overloads the vertical concatenation operator, so that
%! two (or more) time series objects containing the same variables 
%! can be merged using the following syntax:
%!
%!     a = [b; c; d];
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

    if nargin==0 
        a = DynSeries();
    elseif nargin == 1
        a = varargin{1};
    elseif nargin>1
        a = varargin{1};
        for i=2:nargin
            a = vertcat_(a,varargin{i});
        end
    end

end

function d = vertcat_(b, c)
    d = NaN;
    if ~isequal(b.freq, c.freq)
        error('dynSeries::vertcat: Frequencies must be common!')
    end
    if ~isequal(b.vobs, c.vobs)
        error('dynSeries::vertcat: Number of variables must be common!')
    end
    if ~isequal(b.name, c.name)
        error('dynSeries::vertcat: Variables must be common!')
    end
    d = b;
    d.data = [b.data; c.data];
    d.nobs = b.nobs+c.nobs;
end

%@test:1
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'A1';'A2'};
%$
%$ % Define expected results.
%$ e.time = dynDate(1);
%$ e.freq = 1;
%$ e.name = {'A1';'A2'};
%$ e.data = [A;B];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(B,[],B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1;ts2];
%$
%$ % Check the results.
%$
%$ t(1) = dyn_assert(ts3.init,e.time);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ t(5) = dyn_assert(ts3.nobs,20);
%$ T = all(t);
%@eof:1


%@test:2
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),2*transpose(1:10)];
%$ C = [transpose(1:10),3*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'A1';'A2'};
%$ C_name = {'A1';'A2'};
%$
%$ % Define expected results.
%$ e.time = dynDate(1);
%$ e.freq = 1;
%$ e.name = {'A1';'A2'};
%$ e.data = [A;B;C];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(B,[],B_name,[]);
%$ ts3 = dynSeries(C,[],C_name,[]);
%$
%$ % Call the tested method.
%$ ts4 = [ts1; ts2; ts3];
%$
%$ % Check the results.
%$
%$ t(1) = dyn_assert(ts4.init,e.time);
%$ t(2) = dyn_assert(ts4.freq,e.freq);
%$ t(3) = dyn_assert(ts4.data,e.data);
%$ t(4) = dyn_assert(ts4.name,e.name);
%$ t(5) = dyn_assert(ts4.nobs,30);
%$ T = all(t);
%@eof:2
