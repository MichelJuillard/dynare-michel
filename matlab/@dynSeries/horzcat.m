function a = horzcat(varargin)

%@info:
%! @deftypefn {Function file} {@var{a} =} horzcat (@var{b},@var{c}, ...)
%! @anchor{horzcat}
%! @sp 1
%! Method of the dynSeries class.
%! @sp 1
%! Merge Dynare time series objects. This method overloads the horizontal concatenation operator, so that
%! two (or more) time series objects can be merged using the following syntax:
%!
%!     a = [b, c, d];
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
%! @ref{dynSeries}, @ref{private/horzcat2}
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

if nargin==0 
    a = DynSeries();
elseif nargin == 1
    a = varargin{1};
elseif nargin==2
    a = horzcat2(varargin{1},varargin{2});
else
    a = horzcat2(varargin{1},varargin{2});
    for i=3:nargin
        a = horzcat2(a,varargin{i});
    end
end

%@test:1
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define expected results.
%$ e.time = dynDate(1);
%$ e.freq = 1;
%$ e.name = {'A1';'A2';'B1';'B2'};
%$ e.data = [A,B];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(B,[],B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$
%$ t(1) = dyn_assert(ts3.init,e.time);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(5:12),2*transpose(5:12)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial date
%$ A_init = 2001;
%$ B_init = 2005;
%$
%$ % Define expected results.
%$ e.time = [transpose(2000+(1:12)), ones(12,1)];
%$ e.init = dynDate(2001);
%$ e.freq = 1;
%$ e.name = {'A1';'A2';'B1';'B2'};
%$ e.data = [ [A; NaN(2,2)], [NaN(4,2); B]];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,A_init,A_name,[]);
%$ ts2 = dynSeries(B,B_init,B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts3.init,e.init);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a data set.
%$ A = [transpose(1:7),2*transpose(1:7)];
%$ B = [transpose(5:11),2*transpose(5:11)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial date
%$ A_init = '1950Q1';
%$ B_init = '1950Q3';
%$
%$ % Define expected results.
%$ e.freq = 4;
%$ e.init = dynDate('1950Q1');
%$ e.name = {'A1';'A2';'B1';'B2'};
%$ e.data = [ [A; NaN(2,2)], [NaN(2,2); B]];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,A_init,A_name,[]);
%$ ts2 = dynSeries(B,B_init,B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts3.init,e.init);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a data set.
%$ A = [transpose(1:7),2*transpose(1:7)];
%$ B = [transpose(5:9),2*transpose(5:9)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'B1';'B2'};
%$
%$ % Define initial date
%$ A_init = '1950Q1';
%$ B_init = '1950Q3';
%$
%$ % Define expected results.
%$ e.init = dynDate(A_init);
%$ e.freq = 4;
%$ e.name = {'A1';'A2';'B1';'B2'};
%$ e.data = [ A, [NaN(2,2); B]];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,A_init,A_name,[]);
%$ ts2 = dynSeries(B,B_init,B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1,ts2];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts3.init,e.init);
%$ t(2) = dyn_assert(ts3.freq,e.freq);
%$ t(3) = dyn_assert(ts3.data,e.data);
%$ t(4) = dyn_assert(ts3.name,e.name);
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),3*transpose(1:10)];
%$ C = [transpose(1:10),4*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'B1';'B2'};
%$ C_name = {'C1';'C2'};
%$
%$ % Define expected results.
%$ e.init = dynDate(1);
%$ e.freq = 1;
%$ e.name = {'A1';'A2';'B1';'B2';'C1';'C2'};
%$ e.data = [A,B,C];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(B,[],B_name,[]);
%$ ts3 = dynSeries(C,[],C_name,[]);
%$
%$ % Call the tested method.
%$ ts4 = [ts1,ts2,ts3];
%$
%$ % Check the results.
%$ t(1) = dyn_assert(ts4.init,e.init);
%$ t(2) = dyn_assert(ts4.freq,e.freq);
%$ t(3) = dyn_assert(ts4.data,e.data);
%$ t(4) = dyn_assert(ts4.name,e.name);
%$ T = all(t);
%@eof:5

%@test:6
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'B1';'A2'};
%$
%$ % Instantiate two time series objects.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(B,[],B_name,[]);
%$
%$ % Call the tested method.
%$ try
%$   ts3 = [ts1,ts2];
%$   t = 0;
%$ catch
%$   t = 1;
%$ end
%$
%$ T = t;
%@eof:6
