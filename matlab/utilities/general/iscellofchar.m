function a = iscellofchar(b)

% Copyright (C) 2012-2013 Dynare Team
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

a = 1;

if ndims(b)>2
    error(['iscellofchar: Input argument ''' inputname(1) ''' has to be a two dimensional cell array!'])
end

[n,m] = size(b);
p = numel(b);
q = 0;

for i=1:m
    for j=1:n
        if ischar(b{j,i})
            q = q + 1;
        end
    end
end

if ~isequal(q,p)
    a = 0;
end