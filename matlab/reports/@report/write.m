function write(o)
%function write(o)
% Write Report object
%
% INPUTS
%   none
%
% OUTPUTS
%   none
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

fprintf(fid, '\\usepackage[%spaper,margin=2.5cm', o.paper);
if strcmpi(o.orientation, 'landscape')
    fprintf(fid, ',landscape');
end
fprintf(fid, ']{geometry}\n');
fprintf(fid, '\\usepackage{graphicx}\n');
fprintf(fid, '\\usepackage{pdflscape}\n')
fprintf(fid, '\\begin{document}\n');

o.pages.write(fid, addIndentation(''));

fprintf(fid, '\\end{document}\n');
fprintf(fid, '%% End Report Object\n');
status = fclose(fid);
if status == -1
    error('Error closing %s\n', o.filename);
end
end