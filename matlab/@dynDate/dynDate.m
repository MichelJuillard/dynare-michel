function date = dynDate(a,b) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{date} =} dynDate (@var{a})
%! @anchor{dynDate}
%! @sp 1
%! Constructor for the Dynare dates class.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Date. For Quaterly, Monthly or Weekly data, a must be a string. For yearly data or if the frequence is not
%! defined  must be an integer. If @var{a} is a dynDate object, then date will be a copy of this object. If
%! the constructor is called without input argument, it will return an empty dynDate object.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item date
%! Dynare date object.
%! @end table
%! @sp 1
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item time
%! Row vector of integers (1*2) indicating the year and the week, month or quarter of the first observation.
%! @end table
%! @sp 1
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @ref{set_time}
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

date = struct;

date.freq = NaN;
date.time = NaN(1,2);

date = class(date,'dynDate');

switch nargin
  case 0
    % Return an empty dynDate object (without specified frequency).
    return
  case 1
    if ischar(a)% Weekly, Monthly or Quaterly data.
        a = upper(a);
        if length(a)>1
            yearly = findstr('Y',a);
            if isempty(yearly)
                yearly = findstr('A',a);
            end
            quaterly = findstr('Q',a);
            monthly  = findstr('M',a);
            weekly   = findstr('W',a);
            if ~isempty(yearly)
                date.freq = 1;
                date.time(1) = str2num(a(1:yearly-1));
                date.time(2) = 1;
            end
            if ~isempty(quaterly)
                date.freq = 4;
                date.time(1) = str2num(a(1:quaterly-1));
                date.time(2) = str2num(a(quaterly+1:end));
            end
            if ~isempty(monthly)
                date.freq = 12;
                date.time(1) = str2num(a(1:monthly-1));
                date.time(2) = str2num(a(monthly+1:end));
            end
            if ~isempty(weekly)
                date.freq = 52;
                date.time(1) = str2num(a(1:weekly-1));
                date.time(2) = str2num(a(weekly+1:end));
            end
            if isempty(yearly) && isempty(quaterly) && isempty(monthly) && isempty(weekly)
                if isaletter(a)
                    error('dynDate:: Using a string as an input argument, I can only handle weekly (W), monthly (M), quaterly (Q) or yearly (Y, A), data!');
                else
                    % Yearly data declared with a string
                    date.freq = 1;
                    date.time(1) = str2num(a);
                    date.time(2) = 1;
                end
            end
        else
            % Return an empty dynDate object (with a specified frequency).
            switch a
              case {'Y','A'}
                date.freq = 1;
              case 'Q'
                date.freq = 4;
              case 'M'
                date.freq = 12;
              case 'W'
                date.freq = 52;
              otherwise
                error(['dynDate:: With one string argument of length one, you must provide one of weekly (''W''), monthly (''M''), quaterly (''Q'') or yearly (''Y'').']);
            end
        end
    elseif isa(a,'dynDate') % If input argument is a dynDate object then do a copy.
        date = a;
    else
        if isequal(length(a),1) && isnumeric(a) && isint(a)
            % If b is not a string then yearly data are assumed.
            date.freq = 1;
            date.time(1) = a;
            date.time(2) = 1;
        else
            error('dynDate:: Can''t instantiate the class, wrong calling sequence!')
        end            
    end
  case 2 % provide time and freq to instantiate a dynDate object
    date = dynDate();
    if isnumeric(b) && isscalar(b) && ismember(b,[1,4,12,52])
        date.freq = b;
    elseif ischar(b) && isequal(length(b),1) && ismember(upper(b),{'Y','A','W','M','Q'})
        b = upper(b);
        switch b
          case {'Y','A'}
            date.freq = 1;
          case 'W'
            date.freq = 52;
          case 'M'
            date.freq = 12;
          case 'Q'
            date.freq = 4;
          otherwise
            error(['dynDate:: Can''t instantiate the class! The second argument ' inputname(2) ' must be an integer equal to 1, 4, 12 or 52, or a character equal to ''Y'', ''A'', ''W'', ''M'' or ''Q''!'])
        end
    else
        error(['dynDate:: Can''t instantiate the class! The second argument ' inputname(2) ' must be an integer equal to 1, 4, 12 or 52, or a character equal to ''Y'', ''A'', ''W'', ''M'' or ''Q''!'])
    end
    if ~isnumeric(a)
        error(['dynDate:: Can''t instantiate the class! The first argument ' inputname(1) ' must be numeric!'])
    end
    if ~all(isint(a))
        error(['dynDate:: Can''t instantiate the class! The first argument ' inputname(1) ' must be a scalar or a 1*2 vector of integers!'])
    end
    if ~isequal(size(a),[1,2])
        if date.freq>1
            error(['dynDate:: Can''t instantiate the class! The first argument ' inputname(1) ' must be a 1*2 vector of integers.'])
        end
    else
        if isequal(date.freq,1) && ~isequal(a(2),1)
            error(['dynDate:: Can''t instantiate the class! The second element of the first input argument ' inputname(1) ' must be equal to one (because freq==1)!'])
        end
    end
    if a(2)<=0 || a(2)>date.freq
        error(['dynDate:: Can''t instantiate the class! The second element of the first argument ' inputname(1) ' must be a positive integer be <=' int2str(2) '!' ])
    end
    if length(a)==1
        date.time = [a, 1];
    else
        date.time = a;
    end
  otherwise
    error('dynDate:: Can''t instantiate the class, wrong calling sequence!')
end

%@test:1
%$ % Define some dates
%$ date_1 = 1950;
%$ date_2 = '1950Q2';
%$ date_3 = '1950m10';
%$ date_4 = '1950w50';
%$ date_5 = '1950';
%$ date_6 = '1967y';
%$ date_7 = '2009A';
%$
%$ % Define expected results.
%$ e_date_1 = [1950 1];
%$ e_freq_1 = 1;
%$ e_date_2 = [1950 2];
%$ e_freq_2 = 4;
%$ e_date_3 = [1950 10];
%$ e_freq_3 = 12;
%$ e_date_4 = [1950 50];
%$ e_freq_4 = 52;
%$ e_date_5 = [1950 1];
%$ e_freq_5 = 1;
%$ e_date_6 = [1967 1];
%$ e_freq_6 = 1;
%$ e_date_7 = [2009 1];
%$ e_freq_7 = 1;
%$
%$ % Call the tested routine.
%$ d1 = dynDate(date_1);
%$ d2 = dynDate(date_2);
%$ d3 = dynDate(date_3);
%$ d4 = dynDate(date_4);
%$ d5 = dynDate(date_5);
%$ d6 = dynDate(date_6);
%$ d7 = dynDate(date_7);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(d1.time,e_date_1);
%$ t(2) = dyn_assert(d2.time,e_date_2);
%$ t(3) = dyn_assert(d3.time,e_date_3);
%$ t(4) = dyn_assert(d4.time,e_date_4);
%$ t(5) = dyn_assert(d5.time,e_date_5);
%$ t(6) = dyn_assert(d1.freq,e_freq_1);
%$ t(7) = dyn_assert(d2.freq,e_freq_2);
%$ t(8) = dyn_assert(d3.freq,e_freq_3);
%$ t(9) = dyn_assert(d4.freq,e_freq_4);
%$ t(10)= dyn_assert(d5.freq,e_freq_5);
%$ t(11)= dyn_assert(d6.freq,e_freq_6);
%$ t(12)= dyn_assert(d7.freq,e_freq_7);
%$ t(13)= dyn_assert(d6.time,e_date_6);
%$ t(14)= dyn_assert(d7.time,e_date_7);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Instatiate an empty objects for quaterly, monthly and weekly dates.
%$ qq = dynDate('Q');
%$ mm = dynDate('M');
%$ ww = dynDate('W');
%$ yy = dynDate('Y');
%$
%$ % Check the results.
%$ t(1) = dyn_assert(qq.freq,4);
%$ t(2) = dyn_assert(mm.freq,12);
%$ t(3) = dyn_assert(ww.freq,52);
%$ t(4) = dyn_assert(yy.freq,1);
%$ t(5) = dyn_assert(all(isnan(qq.time)),1);
%$ t(6) = dyn_assert(all(isnan(mm.time)),1);
%$ t(7) = dyn_assert(all(isnan(ww.time)),1);
%$ t(8) = dyn_assert(all(isnan(yy.time)),1);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Try to instatiate dynDate objects.
%$ try
%$    a = dynDate([1950 1],4);
%$    t(1) = 1;
%$ catch
%$    t(1) = 0;
%$ end
%$ try
%$    a = dynDate([1950 5],4);
%$    t(2) = 0;
%$ catch
%$    t(2) = 1;
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ % Instatiate an empty objects for quaterly, monthly and weekly dates.
%$ qq = dynDate('q');
%$ mm = dynDate('m');
%$ ww = dynDate('w');
%$
%$ % Check the results.
%$ t(1) = dyn_assert(qq.freq,4);
%$ t(2) = dyn_assert(mm.freq,12);
%$ t(3) = dyn_assert(ww.freq,52);
%$ t(4) = dyn_assert(all(isnan(qq.time)),1);
%$ t(5) = dyn_assert(all(isnan(mm.time)),1);
%$ t(6) = dyn_assert(all(isnan(ww.time)),1);
%$ T = all(t);
%@eof:4

%@test:5
%$ % Try to instatiate dynDate objects.
%$ try
%$    a = dynDate([1950 1],'Q');
%$    t(1) = 1;
%$ catch
%$    t(1) = 0;
%$ end
%$ try
%$    a = dynDate([1950 5],'Q');
%$    t(2) = 0;
%$ catch
%$    t(2) = 1;
%$ end
%$ T = all(t);
%@eof:5
