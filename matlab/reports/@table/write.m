function o = write(o, fid)
%function o = write(o, fid)
% Write a Table object
%
% INPUTS
%   o   - Table Object
%   fid - int, file id
%
% OUTPUTS
%   o   - Table Object
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

assert(fid ~= -1);
if isempty(o.data)
    return
end

if isempty(o.seriestouse)
    ds = o.data;
else
    ds = o.data{o.seriestouse{:}};
end

if ~isempty(o.range)
    ds  = ds(o.range);
end

%number of left-hand columns, 1 until we allow the user to group data,
% e.g.: GDP Europe
%         GDP France
%         GDP Germany
% this example would be two lh columns, with GDP Europe spanning both
nlhc = 1;

disp('creating table.........');
fprintf(fid, '%% Table Object\n');
fprintf(fid, '\\begin{tabular}{l');

dates = ds.time;
ndates = dates.ndat;

for i=1:ndates
    if o.vlines
        fprintf(fid, '|');
    end
    fprintf(fid, 'r');
end
fprintf(fid, '}%%\n');
if ~isempty(o.title)
    fprintf(fid, '\\multicolumn{%d}{c}{%s} \\\\\n', ndates+nlhc, o.title);
end
fprintf(fid, '\\toprule%%\n');

% Column Headers
datedata = dates.time;
years = unique(datedata(:, 1));
thdr = num2cell(years, size(years, 1));
lind = nlhc;
switch dates.freq
    case 1
        for i=1:size(thdr, 1)
            fprintf(fid, ' & %d', thdr{i, 1});
        end
    case 4
        thdr{1, 2} = datedata(:, 2)';
        if size(thdr, 1) > 1
            for i=2:size(thdr, 1)
                split = find(thdr{i-1, 2} == 4, 1, 'first');
                if isempty(split)
                    error('@table.write: Shouldn''t arrive here');
                else
                    thdr{i, 2} = thdr{i-1, 2}(split+1:end);
                    thdr{i-1, 2} = thdr{i-1, 2}(1:split);
                end
            end
        end

        for i=1:size(thdr, 1)
            fprintf(fid, ' & \\multicolumn{%d}{c}{%d}', size(thdr{i,2}, 2), thdr{i,1});
        end
        fprintf(fid, '\\\\[-10pt]%%\n');
        for i=1:size(thdr, 1)
            fprintf(fid, ' & \\multicolumn{%d}{c}{\\hrulefill}', size(thdr{i,2}, 2));
        end
        fprintf(fid, '\\\\%%\n');
        for i=1:size(thdr, 1)
            quarters = thdr{i, 2};
            for j=1:size(quarters, 2)
                fprintf(fid, ' & \\multicolumn{1}{c}{Q%d}', quarters(j));
            end
        end
    case 12
        error('@table.write: weekly dates not yet implemented');
    otherwise
        error('@table.write: invalid dynSeries frequency');
end
fprintf(fid, '\\\\[-10pt]%%\n');
for i=1:ndates
    fprintf(fid, ' & \\hrulefill');
end
fprintf(fid, '\\\\%%\n');
fprintf(fid, '%%\n');

% Table Data
vars = ds.name;
nvars = size(vars);
data = ds.data;
assert(isint(o.precision));
precision = 10^o.precision;
dataString = [' & %.' num2str(o.precision) 'f'];
for i=1:nvars
    fprintf(fid, '%% Table Row %d\n', i);
    fprintf(fid, '%s', vars{i});
    for j=1:ndates
        fprintf(fid, dataString, round(data(j,i)*precision)/precision);
    end
    fprintf(fid, ' \\\\\n\n');
end

fprintf(fid, '\\bottomrule%%\n');
fprintf(fid, '\\end{tabular}%%\n\n');
fprintf(fid, '\\medskip\\par\n\n');
fprintf(fid, '%% End Table Object\n');
end
