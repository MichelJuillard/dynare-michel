function t = dyn_assert(A,B,tol)
% This function tests the equality of two objects.

% Copyright (C) 2011-2012 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

if ( (nargin<3) || isempty(tol) )
    use_isequal_matlab_builtin = 1;
else
    use_isequal_matlab_builtin = 0;
end

[nA,mA] = size(A);
[nB,mB] = size(B);

if nA-nB
    error('assert:: Dimensions of compared objects A and B don''t match!')
end

if mA-mB
    error('assert:: Dimensions of compared objects A and B don''t match!')
end

if isstruct(B) && ~isstruct(A)
    error('assert:: Compared objects are not of the same type!')
end

if iscell(B) && ~iscell(A)
    error('assert:: Compared objects are not of the same type!')
end

if use_isequal_matlab_builtin
    t = isequal(A,B);
    if ~t
        t = isequalwithequalnans(A,B);
    end
else
    t = 1;
    if ~(isstruct(B) || iscell(B))
        if max(abs(A(:)-B(:)))>tol
            t = 0;
        end
    else
        % not yet implemented
        t = NaN;
    end
end