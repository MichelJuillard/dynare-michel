function B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)
% rand_matrix_normal  Pseudo random matrices drawn from a
%                     matrix-normal distribution
%
% B = rand_matrix_normal(n, p, M, Omega_lower_chol, Sigma_lower_chol)
%
% Returns an n-by-p matrix drawn from a Matrix-normal distribution
%
% B ~ MN_n*p(M, Omega, Sigma)   
%
% Equivalent to vec(B) ~ N(vec(Mu), kron(Omega, Sigma))
%
% Same notations than: http://en.wikipedia.org/wiki/Matrix_normal_distribution
%
% M is the mean, n-by-p matrix
% Omega_lower_chol is p-by-p, lower Cholesky decomposition of Omega
% (Omega_lower_chol = chol(Omega, 'lower'))
% Sigma_lower_chol is n-by-n, lower Cholesky decomposition of Sigma
% (Sigma_lower_chol = chol(Sigma, 'lower'))
    
    B1 = randn(n * p, 1);
    B2 = kron(Omega_lower_chol, Sigma_lower_chol) * B1;
    B3 = reshape(B2, n, p);
    B = B3 + M;
