function D = A_times_B_kronecker_C(A,B,C)
% Computes A * kron(B,C). 
%
% INPUTS
%   A  [double] mA*nA matrix.
%   B  [double] mB*nB matrix.
%   C  [double] mC*nC matrix.
%  
% OUTPUTS
%   D  [double] mA*(nC*nB) or mA*(nB*nB) matrix.
%  
% ALGORITHM
%   none.    
%
% SPECIAL REQUIREMENTS
%   none.
%  
%  
% part of DYNARE, copyright Dynare Team (1996-2008)
% Gnu Public License.
    
% Chek number of inputs and outputs.
if nargin>3 | nargin<2
    error('Two or Three input arguments required!')
end
if nargout>1
    error('Too many output arguments!')
end
% Get & check dimensions. Initialization of the output matrix.
[mA,nA] = size(A);
[mB,nB] = size(B);
if nargin == 3
    [mC,nC] = size(C);
    if mB*mc ~= nA
        error('Input dimension error!')
    end
    D = zeros(mA,nB*nC);
    loop = (mB*nB*mC*nC > 1e7);
else
    if mB*mB ~= nA
        error('Input dimension error!')
    end
    D = D = zeros(mA,nB*nB);
    loop = (mB*nB*mB*nB > 1e7);
end
% Computational part.
if loop
    if nargin == 3
        k1 = 1; 
        for i1=1:nB
            for i2=1:nC
                D(:,k1) = A * kron(B(:,i1),C(:,i2));
                k1 = k1 + 1; 
            end
        end
    else
        k1 = 1;
        for i1=1:nB
            for i2=1:nB
                D(:,k1) = A * kron(B(:,i1),B(:,i2));
                k1 = k1 + 1; 
            end
        end
    end
else
    if nargin == 3
        D = A * kron(B,C);
    else
        D = A * kron(B,B);
    end
end