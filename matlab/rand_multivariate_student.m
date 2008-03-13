function draw = rand_multivariate_student(Mean,Sigma_upper_chol,df)
% Pseudo random draws from a multivariate student distribution,
% with expectation Mean, variance Sigma*df/(df-2) and degrees of freedom df>0.
%
% INPUTS 
%
%    Mean               [double]    1*n vector, expectation of the multivariate random variable.
%    Sigma_upper_chol   [double]    n*n matrix, upper triangular Cholesky decomposition of Sigma 
%                                   (the covariance matrix up to a factor df/(df-2)).
%    df                 [integer]   degrees of freedom.
%    
% OUTPUTS 
%    draw               [double]    1*n vector drawn from a multivariate normal distribution with expectation Mean and
%                                   covariance Sigma.
%        
% REMARK This is certainly not the most efficient way...
%
% NOTE See Zellner (appendix B.2, 1971) for a definition.     
%    
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.
    n = length(Mean);
    draw = Mean + randn(1,n) * Sigma_upper_chol * sqrt(df/sum(randn(df,1).^2));