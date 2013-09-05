function B = subsref(A,S) % --*-- Unitary tests --*--

% B = subsref(A,S)
%
% Overloads the subsref method for dynTimeIndex class. This method only allows to get
% the value of the field `index`.
    
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

if length(S)>1
    error('dynTimeIndex::subsref: Something is wrong in your syntax!')
end

if isequal(S.type,'.')
    if isequal(S.subs,'index')
        B = builtin('subsref', A, S(1));
    else
        error(['dynTimeIndex::subsref: ' S.subs  ' is not a known member!'])
    end
else
    error('dynTimeIndex::subsref: Something is wrong in your syntax!')
end

%@test:1
%$ % Instantiate a dynTimeIndex object
%$ u = dynTimeIndex();
%$ try
%$    v = u.index;
%$    t(1) = 1;
%$ catch
%$    t(1) = 0;
%$ end
%$
%$ if t(1)
%$    t(2) = isequal(v,int8(0));
%$ end
%$
%$ T = all(t);
%@eof:1