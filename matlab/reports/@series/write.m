function o = write(o, fid, dates, precision)
%function o = write(o, fid, dates, precision)
% Write Table Row
%
% INPUTS
%   o       [series]    series object
%   xrange  [dynDates]  range of x values for line
%
% OUTPUTS
%   o       [series]    series object
%
% SPECIAL REQUIREMENTS
%   none

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

%% Validate options passed to function
assert(fid ~= -1);
assert(isa(dates, 'dynDates'));
assert(isint(precision));

%% Validate options provided by user
assert(~isempty(o.data) && isa(o.data, 'dynSeries'), ...
       '@series.write: must provide data as a dynSeries');

assert(ischar(o.color), '@series.write: color must be a string');
assert(ischar(o.table_neg_color), '@series.write: table_neg_color must be a string');
assert(ischar(o.table_pos_color), '@series.write: table_pos_color must be a string');
assert(islogical(o.table_markers), '@series.write: table_markers must be a string');

%% Write Output
dataString = ['%.' num2str(precision) 'f'];
precision  = 10^precision;

fprintf(fid, '%% Table Row (series)\n');
fprintf(fid, '%s', o.data.name{:});
data = o.data(dates);
data = data.data;
for i=1:size(data,1)
    thisCellData = round(data(i)*precision)/precision;

    fprintf(fid, ' &');
    if o.table_markers
        if thisCellData < 0
            fprintf(fid, '\\color{%s}', o.table_neg_color);
        elseif thisCellData > 0
            fprintf(fid, '\\color{%s}', o.table_pos_color);
        end
        fprintf(fid, '[');
    end

    fprintf(fid, dataString, thisCellData);

    if o.table_markers
        fprintf(fid, ']');
    end
end
fprintf(fid, ' \\\\\n\n');
end
