function write(o, fid, indent)
%function write(o, fid, indent)
% Write Section object
%
% INPUTS
%   fid - int, file id
%   indent - char, number of spaces to indent tex code
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

assert(fid ~= -1);

fprintf(fid, '\n%s%% Section Object\n', indent);
fprintf(fid, '%s\\noindent\\begin{minipage}[%s]{0.32\\hsize}\n', indent, o.align);

ne = numElements(o);
for i=1:ne
    o.elements(i).write(fid, addIndentation(indent));
end

fprintf(fid, '%s\\end{minipage}\n', indent);
fprintf(fid, '%s%% End Section Object\n\n', indent);
end