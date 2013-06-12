function B = subsref(A,S)

%@info:
%! @deftypefn {Function File} {@var{us} =} subsref (@var{ts},S)
%! @anchor{@dynDate/subsref}
%! @sp 1
%! Overloads the subsref method for the Dynare dates class (@ref{dynDate}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare date object instantiated by @ref{dynDate}.
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
%! A matlab object (public member of the @ref{dynDate} object).
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp2
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

switch S(1).type
  case '.'
    switch S(1).subs
      case 'format'
        B = format(A);
      case {'time', 'freq'}
        B = builtin('subsref', A, S(1));
      otherwise
        error('dynDate::subsref: Unknown public member of method!')
    end
  case '()'
    switch length(S(1).subs)
      case 1
        if ischar(S(1).subs{1})
            if numel(S(1).subs{1})==1 && isempty(strmatch(S(1).subs{1},{'W','M','Q','Y'},'exact'))
                error(['dynDate::subsref: To set the frequency, the input argument of dynDate object ''' inputname(1) ''' should be ''W'', ''M'', ''Q'' or ''Y''.'])
            end
            % Set the frequency (if numel==1) of an empty dynDate object or set the date (if numel>1).
            B = dynDate(S(1).subs{1});
        elseif isnumeric(S(1).subs{1}) && isscalar(S(1).subs{1}) && isint(S(1).subs{1})
            if (~isnan(A.freq) && A.freq==1) || isnan(A.freq)
                B = dynDate(S(1).subs{1});
            else
                error(['dynDate::subsref: dynDate object ''' inputname(1) ''' was not instantiated for years.'])
            end
        else
            error('dynDate::subsref: Something is wrong in your syntax!')
        end
      case 2% Populate an empty dynDate object
        if isnan(A.freq)
            error(['dynDate::subsref: I cannot interpret the two inputs of dynDate object ''' inputname(1) ''' because frequency is not set.'])
        else
            tmp = [];
            switch A.freq
              case 4
                % Quaterly data
                if isint(S(1).subs{2}) && isint(S(1).subs{1}) && S(1).subs{2}<5 && S(1).subs{2}>0
                    tmp = [int2str(S(1).subs{1}), 'Q' int2str(S(1).subs{2})];
                else
                    if ~isint(S(1).subs{2}) || ~(S(1).subs{2}<5 && S(1).subs{2}>0)
                        error(['dynDate::subsref: The second input argument of dynDate object ''' inputname(1) ''' (' num2str(S(1).subs{2})  ') should be a positive integer less than or equal to 4.'])
                    end
                    if ~isint(S(1).subs{2})
                        error(['dynDate::subsref: The first input argument of dynDate object ''' inputname(1) ''' (' num2str(S(1).subs{1})  ') should be an integer.'])
                    end
                end
              case 12
                % Monthly data
                if isint(S(1).subs{2}) && isint(S(1).subs{1}) && S(1).subs{2}<13 && S(1).subs{2}>0
                    tmp = [num2str(S(1).subs{1}), 'M' num2str(S(1).subs{2})];
                else
                    if ~isint(S(1).subs{2}) || ~(S(1).subs{2}<13 && S(1).subs{2}>0)
                        error(['dynDate::subsref: The second input argument of dynDate object ''' inputname(1) ''' (' num2str(S(1).subs{2})  ') should be a positive integer less than or equal to 12.'])
                    end
                    if ~isint(S(1).subs{2})
                        error(['dynDate::subsref: The first input argument of dynDate object ''' inputname(1) ''' (' num2str(S(1).subs{1})  ') should be an integer.'])
                    end
                end
              case 52
                % Weekly data
                if isint(S(1).subs{2}) && isint(S(1).subs{1}) && S(1).subs{2}<53 && S(1).subs{2}>0
                    tmp = [num2str(S(1).subs{1}), 'W' num2str(S(1).subs{2})];
                else
                    if ~isint(S(1).subs{2}) || ~(S(1).subs{2}<53 && S(1).subs{2}>0)
                        error(['dynDate::subsref: The second input argument of dynDate object ''' inputname(1) ''' (' num2str(S(1).subs{2})  ') should be a positive integer less than or equal to 52.'])
                    end
                    if ~isint(S(1).subs{2})
                        error(['dynDate::subsref: The first input argument of dynDate object ''' inputname(1) ''' (' num2str(S(1).subs{1})  ') should be an integer.'])
                    end                    
                end
              case 1
                % Yearly data
                error('dynDate::subsref: Frequency is set for years. You should not provide more than one integer input argument (to set the year)!')
              otherwise
                error('dynDate::subsref: Unknown frequency!')
            end
            if ~isempty(tmp)
                B = dynDate(tmp);
            end
        end
      otherwise
        error(['dynDate::subsref: dynDate object ''' inputname(1) ''' cannot have more than two inputs.'])
    end
  otherwise
    error('dynDate::subsref: Something is wrong in your syntax!')
end

S = shiftS(S);
if ~isempty(S)
    B = subsref(B, S);
end

%@test:1
%$ t = zeros(3,1);
%$
%$ % Instantiate an empty dynDate object
%$ a = dynDate();
%$ if all(isnan(a.time)) && isnan(a.freq)
%$     t(1) = 1;
%$ end
%$
%$ % Populate the empty dynDate object
%$ try
%$     a = a('1950Q1');
%$     if isequal(a.time,[1950 1]) && isequal(a.freq,4)
%$         t(2) = 1;
%$     end
%$ catch
%$     % Nothing to do here...
%$ end
%$
%$ % "Overwrite" a dynDate object
%$ try
%$     a = a('1945Q3');
%$     if isequal(a.time,[1945 3]) && isequal(a.freq,4)
%$         t(3) = 1;
%$     end
%$ catch
%$     % Nothing to do here...
%$ end
%$
%$ % Check the results.
%$ T = all(t);
%@eof:1

%@test:2
%$ % Instantiate a dynDate object
%$ a = dynDate('1938Q4');
%$
%$ % Try to access a non existent (or forbidden) property
%$ try
%$     a.Time;
%$     t = 0;
%$ catch
%$     t = 1;
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ % Try more complex call to overloaded subsref
%$ t = zeros(3,1);
%$ try
%$     a = dynDate('M');
%$     time = a('1973M1').time;
%$     t(1) = 1;
%$     t(2) = dyn_assert(time,[1973,1]);
%$     t(3) = dyn_assert(a.freq,12);
%$ catch
%$     % Nothing to do here.
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ t = NaN(3,1);
%$ % Instantiate an empty object for quaterly date
%$ qq = dynDate('Q');
%$ % Populate this object
%$ a = qq(1938,4);
%$ try
%$    a = dynDate(1938,11);
%$    t(3) = 0;
%$ catch
%$    t(3) = 1;
%$ end
%$
%$ % Check the results
%$ t(1) = dyn_assert(a.freq,4);
%$ t(2) = dyn_assert(a.time,[1938,4]);
%$ T = all(t);
%@eof:4

%@test:5
%$ t = NaN(2,1);
%$ % Instantiate an empty object for quaterly date
%$ qq = dynDate('Q');
%$ % Populate this object and get the time member
%$ time = qq(1938,4).time;
%$ 
%$ % Check the results
%$ t(1) = dyn_assert(qq.freq,4);
%$ t(2) = dyn_assert(time,[1938,4]);
%$ T = all(t);
%@eof:5
