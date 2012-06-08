function B = subsref(A, S)
%@info:
%! @deftypefn {Function File} {@var{B} =} subsref (@var{A},@var{S})
%! @anchor{@dynTime/subsref}
%! @sp 1
%! Overloads the subsref method for the Dynare time series class (@ref{dynTime}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare time series object instantiated by @ref{dynSeries}.
%! @item S
%! Matlab's structure array S with two fields, type and subs. The type field is string containing '()', '@{@}', or '.', where '()' specifies
%! integer subscripts, '@{@}' specifies cell array subscripts, and '.' specifies subscripted structure fields. The subs field is a cell array
%! or a string containing the actual subscripts (see matlab's documentation).
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item B
%! Dynare time series object. Depending on the calling sequence @var{us} is a transformation of @var{ts} obtained by applying a public method on @var{ts},
%! or a dynSeries object built by extracting a variable from @var{ts}, or a dynSeries object containing a subsample of the all the variable in @var{ts}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{dynTime}, @ref{@@dynTime/setFreq}, @ref{@@dynTime/setSize}, @ref{@@dynTime/setTime}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2012 Dynare Team
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

if isequal(S(1).type,'.')
    switch S(1).subs
      case {'time','freq','init','last'}                                   % Public members.
        B = builtin('subsref', A, S(1));
      case {'setFreq','setSize','setTime'}                                 % Give "dot access" to public methods.
        if length(S)==1
            B = feval(S(1).subs,A);
        else
            if isequal(S(2).type,'()')
                B = feval(S(1).subs,A,S(2).subs{:});
            else
                error('dynTime::subsref: Something is wrong in your syntax!')
            end
        end
      otherwise
        error('dynTime::subsref: Unknown public method or member!')
    end
elseif isequal(S.type,'()')                                                    % Extract a sub-sample.
    if length(S.subs)==1
        S.subs = [S.subs, ':'];
    end
    B = builtin('subsref', A.time, S);
else
    error('dynTime::subsref: Something is wrong in your syntax!')
end

%@test:1
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1(2:9);
%$
%$ % Expected results.
%$ e.data = [transpose(2:9),2*transpose(2:9)];
%$ e.nobs = 8;g
%$ e.vobs = 2;
%$ e.name = char('A1','A2');
%$ e.freq = 1;
%$ tmp = ts1.time; e.time = tmp(2:9,:);
%$ e.init = e.time(1,:);
%$ e.last = e.time(end,:);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.time,e.time);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ t(6) = dyn_assert(a.init,e.init);
%$ t(7) = dyn_assert(a.last,e.last);
%$ T = all(t);
%@eof:1

%@test:2
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1.A1;
%$
%$ % Expected results.
%$ e.data = transpose(1:10);
%$ e.nobs = 10;
%$ e.vobs = 1;
%$ e.name = char('A1');
%$ e.freq = 1;
%$ e.time = [transpose(1:10),ones(10,1)];
%$ e.init = e.time(1,:);
%$ e.last = e.time(end,:);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.time,e.time);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ t(6) = dyn_assert(a.init,e.init);
%$ t(7) = dyn_assert(a.last,e.last);
%$ T = all(t);
%@eof:2

%@test:3
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1.log;
%$
%$ % Expected results.
%$ e.data = log(A);
%$ e.nobs = 10;
%$ e.vobs = 2;
%$ e.name = char('A1','A2');
%$ e.freq = 1;
%$ tmp = ts1.time; e.time = tmp(1:10,:);
%$ e.init = e.time(1,:);
%$ e.last = e.time(end,:);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.time,e.time);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ t(6) = dyn_assert(a.init,e.init);
%$ t(7) = dyn_assert(a.last,e.last);
%$ T = all(t);
%@eof:3
