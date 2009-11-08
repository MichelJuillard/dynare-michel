function T = reduced_rank_cholesky(X)
% Computes the cholesky decomposition of a symetric semidefinite matrix or of a definite positive matrix.
%
% INPUTS:
%    X   [double]   n*n matrix to be factorized.    
%
% OUTPUTS
%    T   [double]   q*n matrix such that T'*T = X, where q is the number of positive eigenvalues in X.
%
% NOTES:
%    If X is not positive definite, then X has to be a symetric semidefinite matrix.
%    The matrix T is upper triangular iff X is positive definite.     

% Copyright (C) 2009 Dynare Team
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
    
    
    [T,X_is_not_positive_definite] = chol(X);
    
    if X_is_not_positive_definite
        n = length(X);
        [U,D] = eig(.5*(X+X'));
        [tmp,max_elements_indices] = max(abs(U),[],1);
        negloc = (U(max_elements_indices+(0:n:(n-1)*n))<0);
        U(:,negloc) = -U(:,negloc);
        D = diag(D);
        tol = eps(max(D)) * length(D)*100;
        t = (abs(D) > tol);
        D = D(t);
        if ~(sum(D<0))
            T = diag(sqrt(D)) * U(:,t)';
        else
            disp('reduced_rank_cholesky:: Input matrix is not semidefinite positive!')
            T = NaN;
        end
    end