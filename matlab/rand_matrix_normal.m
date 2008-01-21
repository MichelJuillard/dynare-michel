function B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)

% function B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)
% Pseudo random matrices drawn from a matrix-normal distribution
% B ~ MN_n*p(M, Omega, Sigma) 
% Equivalent to vec(B) ~ N(vec(Mu), kron(Omega, Sigma))
%
% INPUTS
%    n:                 row
%    p:                 column
%    M:                 (n*p) matrix, mean
%    Omega_lower_chol:  (p*p), lower Cholesky decomposition of Omega,
%                       (Omega_lower_chol = chol(Omega, 'lower'))
%    Sigma_lower_chol:  (n*n), lower Cholesky decomposition of Sigma,
%                       (Sigma_lower_chol = chol(Sigma, 'lower'))
%    
% OUTPUTS
%    B:                 (n*p) matrix drawn from a Matrix-normal distribution
%        
% SPECIAL REQUIREMENTS
%    Same notations than: http://en.wikipedia.org/wiki/Matrix_normal_distribution
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

    B1 = randn(n * p, 1);
    B2 = kron(Omega_lower_chol, Sigma_lower_chol) * B1;
    B3 = reshape(B2, n, p);
    B = B3 + M;
