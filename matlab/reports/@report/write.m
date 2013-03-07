function o = write(o)
%function o = write(o)
% Write Report object
%
% INPUTS
%   o   - Report Object
%
% OUTPUTS
%   o   - Report Object
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

[fid, msg] = fopen(o.filename, 'w');
if fid == -1
    error(msg);
end

fprintf(fid, '%% Report Object\n');
fprintf(fid, '\\documentclass[11pt]{article}\n');

fprintf(fid, '\\usepackage[%spaper,margin=%s', o.paper, o.margin);
if strcmpi(o.orientation, 'landscape')
    fprintf(fid, ',landscape');
end
fprintf(fid, ']{geometry}\n');
fprintf(fid, '\\usepackage{graphicx, pdflscape, pgf, pgfplots}\n');
fprintf(fid, ['\\makeatletter\n' ...
              '\\def\\blfootnote{\\gdef\\@thefnmark{}\\@footnotetext}\n' ...
              '\\makeatother\n']);
if o.showdate
    fprintf(fid, '\\usepackage{fancyhdr, datetime}\n');
    fprintf(fid, '\\newdateformat{reportdate}{\\THEDAY\\ \\shortmonthname\\ \\THEYEAR}\n');
    fprintf(fid, '\\pagestyle{fancy}\n');
    fprintf(fid, '\\renewcommand{\\headrulewidth}{0pt}\n');
    fprintf(fid, '\\renewcommand{\\footrulewidth}{0.5pt}\n');
    fprintf(fid, '\\rfoot{\\scriptsize\\reportdate\\today\\ -- \\currenttime}\n');
end
fprintf(fid, '\\begin{document}\n');

o.pages.write(fid);

fprintf(fid, '\\end{document}\n');
fprintf(fid, '%% End Report Object\n');
status = fclose(fid);
if status == -1
    error('Error closing %s\n', o.filename);
end
end