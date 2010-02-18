function C = bsxfun(fun,A,B)
% (Imperfect) Clone of matlab's bsxfun built-in function.
    
% Copyright (C) 2010 Dynare Team
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

dA = size(A);
dB = size(B);

if length(dB)~=length(dA)
    % Note that this function crashes if A and B are k-dimensional with k>2 and if the size of the last dimension is one.
    % It is a bug, but we don't need this feature, so I do not investigate further... 
    error(['A is a ' int2str(length(dA)) '-dimensional array whereas B is a ' int2str(length(dB)) '-dimensional array!'])
end

if all(dA==dB)
    C = fun(A,B);
else
    tB = dB<dA;
    tA = dA<dB;
    if any(tB)
        iB = find(tB);
        if all(dB(iB)==1)
            B = repmat(B,tB.*dA+~tB);
        else
            error('Non-singleton dimensions of the two input arrays must match each other.')
        end
    end
    if any(tA)
        iA = find(tA);
        if all(dA(iA)==1)
            A = repmat(A,tA.*dB+~tA);
        else
            error('Non-singleton dimensions of the two input arrays must match each other.')
        end
    end
    C = fun(A,B);
end