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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

% Allow to populate an empty dynDate object or update a dynDate object
if isequal(length(S),1) && isequal(S.type,'()')
    if isequal(length(S.subs),1) && ischar(S.subs{1})
        B = dynDate(S.subs{1});
        return
    elseif isequal(length(S.subs),1) && isnumeric(S.subs{1})
        % Yearly data are assumed.
        if isequal(A.freq,1)
            B = dynDate(S.subs{1});
            return
        end
    elseif isequal(length(S.subs),2) && isequal(length(S.subs{1}),1) && isequal(length(S.subs{2}),1)
        tmp = [];
        switch A.freq
          case 4
            % Quaterly data
            if S.subs{2}<5 && S.subs{2}>0
                tmp = [num2str(S.subs{1}), 'Q' num2str(S.subs{2})];
            end
          case 12
            % Monthly data
            if S.subs{2}<13 && S.subs{2}>0
                tmp = [num2str(S.subs{1}), 'M' num2str(S.subs{2})];
            end
          case 52
            % Weekly data
            if S.subs{2}<53 && S.subs{2}>0
                tmp = [num2str(S.subs{1}), 'W' num2str(S.subs{2})];
            end
          otherwise
            %
        end
        if ~isempty(tmp)
            B = dynDate(tmp);
            return
        end
    end
end

% Give access to dynDate methods (format).
if isequal(length(S),1) && isequal(S.type,'.') && ( strcmp(S.subs,'format') )
    B = format(A);
    return
end


% Give access to dynDate properties (time and freq).
if isequal(length(S),1) && isequal(S.type,'.') && ( strcmp(S.subs,'time') || strcmp(S.subs,'freq') )
    B = builtin('subsref', A, S);
    return
end

% Allow more complex call to subsref such that:
%
% a = dynDate();
% a('2009M4').time
%
% should return a row vector [2009 4]. Note that the object name should not match any function name
% declared in the matlab's path.
if length(S)>1 && isequal(S(1).type,'()') && isequal(S(2).type,'.')
    tmp = dynDate(S(1).subs{1});
    B = builtin('subsref', tmp, S(2));
    return
end

error('dynDate::subsref: You''re trying to do something wrong!')

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
%$ t = zeros(1,1);
%$ try
%$     a = dynDate();
%$     time = a('1973M1').time;
%$     t(1) = 1;
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
%@eof:3