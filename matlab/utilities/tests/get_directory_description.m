function flist = get_directory_description(basedir)

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

dd = dir(basedir);
flist = {};
file = 1;

for f=1:length(dd)
    if ~(isequal(dd(f).name,'.') || isequal(dd(f).name,'..'))
        if dd(f).isdir
            flist(file) = { get_directory_description([ basedir filesep dd(f).name]) };
        else
            flist(file) = { [basedir filesep dd(f).name] };
        end
        file = file + 1; 
    end
end