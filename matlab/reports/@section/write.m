function o = write(o, fid)
%function o = write(o, fid)
% Write Section object
%
% INPUTS
%   fid - int, file id
%
% OUTPUTS
%   o   - this
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

fprintf(fid, '%% Section Object\n');
fprintf(fid, '\\maxsizebox{\\textwidth}{!}{%%\n');
fprintf(fid, '\\begin{tabular}[t]{');
for i=1:o.cols
    fprintf(fid, 'c');
end
fprintf(fid, '}\n');
ne = numElements(o);
for i=1:ne
    o.elements(i).write(fid);
    if rem(i, o.cols)
        fprintf(fid, ' & ');
    else
        fprintf(fid, '\\\\\n');
    end
end
fprintf(fid, '\\end{tabular}}\n');
fprintf(fid, '%% End Section Object\n\n');
end