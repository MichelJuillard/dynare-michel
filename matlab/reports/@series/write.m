function o = write(o, fid, dates, precision, yrsForAvgs)
%function o = write(o, fid, dates, precision, yrsForAvgs)
% Write Table Row
%
% INPUTS
%   o            [series]    series object
%   fid          [int]       file id
%   dates        [dynDates]  dates for series slice
%   precision    [float]     precision with which to print the data
%   yrsForAvgs   [bool]      the years for which to compute averages
%
%
% OUTPUTS
%   o            [series]    series object
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
assert(ischar(o.tableSubSectionHeader), '@series.write: tableSubSectionHeader must be a string');
if isempty(o.tableSubSectionHeader)
    assert(~isempty(o.data) && isa(o.data, 'dynSeries'), ...
           '@series.write: must provide data as a dynSeries');
end

assert(ischar(o.tableNegColor), '@series.write: tableNegColor must be a string');
assert(ischar(o.tablePosColor), '@series.write: tablePosColor must be a string');
assert(ischar(o.tableRowColor), '@series.write: tableRowColor must be a string');
assert(islogical(o.tableShowMarkers), '@series.write: tableShowMarkers must be true or false');
assert(islogical(o.tableAlignRight), '@series.write: tableAlignRight must be true or false');
assert(isfloat(o.tableMarkerLimit), '@series,write: tableMarkerLimit must be a float');

%% Write Output
dataString = ['%.' num2str(precision) 'f'];
precision  = 10^precision;

fprintf(fid, '%% Table Row (series)\n');
if ~isempty(o.tableRowColor)
    fprintf(fid, '\\rowcolor{%s}', o.tableRowColor);
end
if ~isempty(o.tableSubSectionHeader)
    fprintf(fid, '%s', o.tableSubSectionHeader);
    for i=1:size(dates)+length(yrsForAvgs)
        fprintf(fid, '&');
    end
    fprintf(fid, '\\\\%%\n');
    return;
end
if o.tableAlignRight
    fprintf(fid, '\\multicolumn{1}{r}{');
end
fprintf(fid, '%s', o.data.tex{:});
if o.tableAlignRight
    fprintf(fid, '}');
end
data = o.data(dates);
data = data.data;
for i=1:size(data,1)
    fprintf(fid, '&');
    if o.tableShowMarkers
        if data(i) < -o.tableMarkerLimit
            fprintf(fid, '\\color{%s}', o.tableNegColor);
        elseif data(i) > o.tableMarkerLimit
            fprintf(fid, '\\color{%s}', o.tablePosColor);
        end
        fprintf(fid, '[');
    end

    fprintf(fid, dataString, round(data(i)*precision)/precision);

    if o.tableShowMarkers
        fprintf(fid, ']');
    end
end

% Calculate and display annual averages
for i=1:length(yrsForAvgs)
    slice = o.data(dynDate([num2str(yrsForAvgs(i)) 'q1']):dynDate([num2str(yrsForAvgs(i)) 'q4']));
    avg = sum(slice.data)/size(slice);
    fprintf(fid, '&');
    if o.tableShowMarkers
        if avg < -o.tableMarkerLimit
            fprintf(fid, '\\color{%s}', o.tableNegColor);
        elseif avg > o.tableMarkerLimit
            fprintf(fid, '\\color{%s}', o.tablePosColor);
        end
        fprintf(fid, '[');
    end

    fprintf(fid, dataString, round(avg*precision)/precision);

    if o.tableShowMarkers
        fprintf(fid, ']');
    end
end
fprintf(fid, '\\\\%%\n');
end
