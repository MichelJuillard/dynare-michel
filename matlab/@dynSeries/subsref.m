function us = subsref(ts, S)
%@info:
%! @deftypefn {Function File} {@var{us} =} subsref (@var{ts},S)
%! @anchor{@dynSeries/subsref}
%! @sp 1
%! Overloads the subsref method for the Dynare time series class (@ref{dynSeries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ts
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
%! @item us
%! Dynare time series object. Depending on the calling sequence @var{us} is a transformation of @var{ts} obtained by applying a public method on @var{ts},
%! or a dynSeries object built by extracting a variable from @var{ts}, or a dynSeries object containing a subsample of the all the variable in @var{ts}.
%! @end table
%! @sp 2
%! @strong{Example 1.} Let @var{ts} be a dynSeries object containing three variables named 'A1', 'A2' and 'A3'. Then the following syntax:
%! @example
%!   us = ts.A1;
%! @end example
%!will create a new dynSeries object @var{us} containing the variable 'A1'.
%! @sp 1
%! @strong{Example 2.} Let @var{ts} be a dynSeries object. Then the following syntax:
%! @example
%!   us = ts.log;
%! @end example
%!will create a new dynSeries object @var{us} containing all the variables of @var{ts} transformed by the neperian logarithm.
%! @sp 1
%! @strong{Example 3.} Let @var{ts} be a dynSeries object. The following syntax:
%! @example
%!   us = ts(3:50);
%! @end example
%!will create a new dynSeries object @var{us} by selecting a subsample out of @var{ts}.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{dynSeries}, @ref{log}, @ref{exp}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011, 2012 Dynare Team
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

if length(S)==1 && isequal(S.type,'.')
    switch S.subs
      case {'data','nobs','vobs','name','tex','freq','time','init'}        % Public members.
        us = builtin('subsref', ts, S);
      case {'log','exp'}                                                   % Give "dot access" to public methods.
        us = feval(S.subs,ts);
      otherwise                                                            % Extract a sub-object by selecting one variable.
        ndx = strmatch(S.subs,ts.name);
        if ~isempty(ndx)
            us = dynSeries();
            us.data = ts.data(:,ndx);
            us.name = deblank(ts.name(ndx,:));
            us.tex  = deblank(ts.tex(ndx,:));
            us.nobs = ts.nobs;
            us.vobs = 1;
            us.freq = ts.freq;
            us.init = ts.init;
            return
        else
            error('dynSeries::subsref: Unknown public method, public member or variable!')
        end
    end
    return
end

if length(S)==1 && isequal(S.type,'()')                                                    % Extract a sub-object by selecting a sub-sample.
    us = dynSeries();
    if size(ts.data,2)>1
        S.subs = [S.subs, ':'];
    end
    us.data = builtin('subsref', ts.data, S);
    us.nobs = size(us.data,1);
    us.vobs = ts.vobs;
    us.freq = ts.freq;
    us.time = builtin('subsref', ts.time, S);
    us.init = ts.init+S.subs{1}(1);
    us.name = ts.name;
    us.tex  = ts.tex;
    return
end

if (length(S)==2) && (isequal(S(1).subs,'init'))
    if isequal(S(2).type,'.') && ( isequal(S(2).subs,'freq') || isequal(S(2).subs,'time') )
        us = builtin('subsref', ts.init, S(2));
    else
        error('dynSeries:subsref:: I don''t understand what you are trying to do!')
    end
end

%@test:1
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1(2:9);
%$
%$ % Expected results.
%$ e.data = [transpose(2:9),2*transpose(2:9)];
%$ e.nobs = 8;
%$ e.vobs = 2;
%$ e.name = {'A1';'A2'};
%$ e.freq = 1;
%$ e.init = dynDate(2);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.nobs,e.nobs);
%$ t(3) = dyn_assert(a.vobs,e.vobs);
%$ t(4) = dyn_assert(a.freq,e.freq);
%$ t(5) = dyn_assert(a.init,e.init);
%$ T = all(t);
%@eof:1

%@test:2
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
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
%$ e.name = {'A1'};
%$ e.freq = 1;
%$ e.init = dynDate(1);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.init,e.init);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ T = all(t);
%@eof:2

%@test:3
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
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
%$ e.name = {'A1';'A2'};
%$ e.freq = 1;
%$ e.init = dynDate(1);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.nobs,e.nobs);
%$ t(3) = dyn_assert(a.vobs,e.vobs);
%$ t(4) = dyn_assert(a.freq,e.freq);
%$ t(5) = dyn_assert(a.init,e.init);
%$ T = all(t);
%@eof:3
