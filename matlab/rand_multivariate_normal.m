function draw = rand_multivariate_normal(Mean,Sigma_upper_chol,n)
% Pseudo random draws from a multivariate normal distribution,
% \mathcal N_n(Mean,Sigma), with expectation Mean and variance Sigma.
%
% INPUTS 
%
%    Mean               [double]    1*n vector, expectation of the multivariate random variable.
%    Sigma_upper_chol   [double]    n*n matrix, upper triangular Cholesky decomposition of Sigma (the covariance matrix).
%    n                  [integer]   dimension.
%    
% OUTPUTS 
%    draw               [double]    1*n vector drawn from a multivariate normal distribution with expectation Mean and
%                                   covariance Sigma 
%        
% SPECIAL REQUIREMENTS
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.
    draw = Mean + randn(1,n) * Sigma_upper_chol;