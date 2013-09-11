function ts = dynSeries(varargin) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{ts} =} dynSeries (@var{a},@var{b},@var{c},@var{d})
%! @anchor{dynSeries}
%! @sp 1
%! Constructor for the Dynare time series class.
%! @sp 2
%! @strong{Inputs}
%! @sp 2
%! If @code{nargin==0} then an empty dynSeries object is created. The object can be populated with data subsequently using the overloaded subsref method.
%! @sp 2
%! If @code{nargin==1} and if the input argument is a @ref{dynDate} object, then a dynSeries object without data is created. This object can be populated with the overload subsref method.
%! @sp 2
%! If @code{nargin==1} and if the input argument is a string for the name of a csv, m or mat file containing data, then a dynSeries object is created from these data.
%! @sp 2
%! If @code{nargin>1}:
%! @sp 1
%! @table @ @var
%! @item a
%! T*1 vector or T*N matrix of data.
%! @item b
%! Initial date. For Quaterly, Monthly or Weekly data, b must be a string. For yearly data or if the frequence is not defined b must be an integer.
%! @item c
%! N*1 cell array of strings. Names of the N time series.
%! @item d
%! N*p array of characters. TeX names of the N time series.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Dynare time series object.
%! @end table
%! @sp 2
%! @strong{Properties}
%! @sp 1
%! The constructor defines the following properties:
%! @sp 1
%! @table @ @var
%! @item data
%! Array of doubles (nobs*vobs).
%! @item nobs
%! Scalar integer, the number of observations.
%! @item vobs
%! Scalar integer, the number of variables.
%! @item name
%! Cell array of strings, names of the variables.
%! @item tex
%! Cell array of strings, tex names of the variables.
%! @item freq
%! Scalar integer, the frequency of the time series. @var{freq} is equal to 1 if data are on a yearly basis or if
%! frequency is unspecified. @var{freq} is equal to 4 if data are on a quaterly basis. @var{freq} is equal to
%! 12 if data are on a monthly basis. @var{freq} is equal to 52 if data are on a weekly basis.
%! @item init
%! @ref{dynDate} object, initial date of the dataset.
%! @end table
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

ts = struct;

ts.data = [];
ts.nobs = 0;
ts.vobs = 0;
ts.name = {};
ts.tex  = {};
ts.freq = [];
ts.init = dynDate();
ts.time = dynDates();

ts = class(ts,'dynSeries');

switch nargin
  case 0
    %  Create an empty dynSeries object.
    return
  case 1
    if isa(varargin{1},'dynDate')
        if isempty(varargin{1})
            error(['dynSeries:: ' inputname(1) ' (identified as a dynDate object) must be non empty!'])
        else
            % Create an empty dynSeries object with an initial date.
            ts.init = varargin{1};
            ts.freq = varargin{1}.freq;
        end
        return
    elseif ischar(varargin{1})
        % Create a dynSeries object loading data in a file (*.csv, *.m, *.mat).
        if isempty(varargin{1})
            error('dynSeries:: Wrong calling sequence! Input argument cannot be an empty string.')
        elseif check_file_extension(varargin{1},'m')
            [freq,init,data,varlist,tex] = load_m_file_data(varargin{1});
        elseif check_file_extension(varargin{1},'mat')
            [freq,init,data,varlist,tex] = load_mat_file_data(varargin{1});
        elseif check_file_extension(varargin{1},'csv')
            [freq,init,data,varlist] = load_csv_file_data(varargin{1});
            tex = [];
        elseif check_file_extension(varargin{1},'xls') || check_file_extension(varargin{1},'xlsx')
            if ~isempty(who('global','options_'));
                % Check that the object is instantiated within a dynare session so that options_ global structure exists.
                % Should provide latter a mechanism to pass range and sheet to dynSeries constructor...
                range = evalin('base','options_.xls_range');
                sheet = evalin('base','options_.xls_sheet');
            else
                % By default only the (whole) first sheet is loaded.
                range = [];
                sheet = [];
            end
            [freq,init,data,varlist] = load_xls_file_data(varargin{1}, sheet, range);
            tex = [];
        else
            error(['dynSeries:: I''m not able to load data from ' inputname(1) '!'])
        end
        ts.init = init;
        ts.freq = freq;
        ts.data = data;
        ts.name = varlist;
        ts.vobs = length(varlist);
        ts.nobs = size(data,1);
        if isempty(tex)
            ts.tex = name2tex(varlist);
        else
            ts.tex = tex;
        end
    elseif isnumeric(varargin{1}) && isequal(ndims(varargin{1}),2)
        ts.data = varargin{1};
        [ts.nobs, ts.vobs] = size(ts.data);
        ts.freq = 1;
        ts.init = dynDate(1);
        ts.time = ts.init:ts.init+ts.nobs;
        ts.name = default_name(ts.vobs);
        ts.tex = name2tex(ts.name);
    end
  case {2,3,4}
    a = varargin{1};
    b = varargin{2};
    if nargin<4
        d = {};
    else
        d = varargin{4};
    end
    if nargin<3
        c = {};
    else
        c = varargin{3};
    end
    % Get data, number of observations and number of variables.
    ts.data = a;
    ts.nobs = size(a,1);
    ts.vobs = size(a,2);
    % Get the first date and set the frequency.
    if ~isempty(b)
        if ischar(b)% Weekly, Monthly or Quaterly data.
            ts.init = dynDate(b);
            ts.freq = ts.init.freq;
        elseif isa(b,'dynDate') && ~isempty(b)
            ts.freq = b.freq;
            ts.init = b;
        elseif isnumeric(b) && isreal(b) && isint(b)
            ts.freq = 1;
            ts.init = dynDate(b);
        else
            error('dynSeries::dynSeries: Wrong calling sequence!');
        end
    else% If b is empty.
        ts.freq = 1;
        ts.init = dynDate(1);
    end
    % Get the names of the variables.
    if ~isempty(c)
        if ts.vobs==length(c)
            for i=1:ts.vobs
                ts.name = vertcat(ts.name, c(i) );
            end
        else
            error('dynSeries::dynSeries: The number of declared names does not match the number of variables!')
        end
    else
        ts.name = default_name(ts.vobs);
    end
    if ~isempty(d)
        if ts.vobs==length(d)
            for i=1:ts.vobs
                ts.tex = vertcat(ts.tex, d(i));
            end
        else
            error('dynSeries::dynSeries: The number of declared tex names does not match the number of variables!')
        end
    else
        ts.tex = name2tex(ts.name);
    end
  otherwise
    error('dynSeries::dynSeries: Can''t instantiate the class, wrong calling sequence!')
end

ts.time = ts.init:(ts.init+ts.nobs);

%@test:1
%$ % Test if we can instantiate an empty dynSeries object.
%$ try
%$     ts = dynSeries();
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(4,1);
%$
%$ try
%$     aa = dynDate('1938M11');
%$     ts = dynSeries(aa);
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,12);
%$     t(3) = dyn_assert(ts.init.freq,12);
%$     t(4) = dyn_assert(ts.init.time,[1938, 11]);
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ t = zeros(6,1);
%$
%$ try
%$     ts = dynSeries('dynseries_test_data.m');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1994, 3]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,100);
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ t = zeros(6,1);
%$
%$ try
%$     ts = dynSeries('dynseries_test_data.mat');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1994, 3]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,100);
%$ end
%$
%$ T = all(t);
%@eof:4

%@test:5
%$ t = zeros(8,1);
%$
%$ try
%$     ts = dynSeries('dynseries_test_data.csv');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1990, 1]);
%$     t(5) = dyn_assert(ts.vobs,4);
%$     t(6) = dyn_assert(ts.nobs,4);
%$     t(7) = dyn_assert(ts.name,{'azert';'yuiop';'qsdfg';'jklm'});
%$     t(8) = dyn_assert(ts.tex,{'azert';'yuiop';'qsdfg';'jklm'});
%$ end
%$
%$ T = all(t);
%@eof:5

%@test:6
%$ t = zeros(8,1);
%$
%$ %try
%$     ts = dynSeries(transpose(1:5),[]);
%$     t(1) = 1;
%$ %catch
%$ %    t = 0;
%$ %end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,1);
%$     t(3) = dyn_assert(ts.init.freq,1);
%$     t(4) = dyn_assert(ts.init.time,[1, 1]);
%$     t(5) = dyn_assert(ts.vobs,1);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Variable_1'});
%$     t(8) = dyn_assert(ts.tex,{'Variable\\_1'});
%$ end
%$
%$ T = all(t);
%@eof:6

%@test:7
%$ t = zeros(8,1);
%$
%$ try
%$     ts = dynSeries(transpose(1:5),'1950Q1');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1950, 1]);
%$     t(5) = dyn_assert(ts.vobs,1);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Variable_1'});
%$     t(8) = dyn_assert(ts.tex,{'Variable\\_1'});
%$ end
%$
%$ T = all(t);
%@eof:7

%@test:8
%$ t = zeros(8,1);
%$
%$ try
%$     ts = dynSeries([transpose(1:5), transpose(6:10)],'1950q1',{'Output'; 'Consumption'}, {'Y_t'; 'C_t'});
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1950, 1]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Output'; 'Consumption'});
%$     t(8) = dyn_assert(ts.tex,{'Y_t'; 'C_t'});
%$ end
%$
%$ T = all(t);
%@eof:8

%@test:9
%$ try
%$     ts = dynSeries('dynseries_test_data-1.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1990, 1]);
%$     t(5) = dyn_assert(ts.vobs,3);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'GDP';'Consumption';'CPI'});
%$     t(8) = dyn_assert(ts.tex,{'GDP';'Consumption';'CPI'});
%$ end
%$
%$ T = all(t);
%@eof:9

%@test:10
%$ try
%$     ts = dynSeries('dynseries_test_data-2.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1990, 1]);
%$     t(5) = dyn_assert(ts.vobs,3);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Variable_1';'Variable_2';'Variable_3'});
%$     t(8) = dyn_assert(ts.tex,{'Variable\\_1';'Variable\\_2';'Variable\\_3'});
%$ end
%$
%$ T = all(t);
%@eof:10

%@test:11
%$ try
%$     ts = dynSeries('dynseries_test_data-3.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(ts.freq,1);
%$     t(3) = dyn_assert(ts.init.freq,1);
%$     t(4) = dyn_assert(ts.init.time,[1, 1]);
%$     t(5) = dyn_assert(ts.vobs,3);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Variable_1';'Variable_2';'Variable_3'});
%$     t(8) = dyn_assert(ts.tex,{'Variable\\_1';'Variable\\_2';'Variable\\_3'});
%$ end
%$
%$ T = all(t);
%@eof:11

%@test:12
%$ try
%$     ts = dynSeries('dynseries_test_data-4.xls');
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ if t(1)
%$     t(2) = dyn_assert(ts.freq,1);
%$     t(3) = dyn_assert(ts.init.freq,1);
%$     t(4) = dyn_assert(ts.init.time,[1, 1]);
%$     t(5) = dyn_assert(ts.vobs,3);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'GDP';'Consumption';'CPI'});
%$     t(8) = dyn_assert(ts.tex,{'GDP';'Consumption';'CPI'});
%$ end
%$
%$ T = all(t);
%@eof:12