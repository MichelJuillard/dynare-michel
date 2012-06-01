function m = mean_preserving_spread(autoregressive_parameter,sigma)
% Computes the mean preserving spread for first order autoregressive process.
%
% The mean preserving spread m is a constant such that the mean of the process
%
%   X_t = X^{\star} * e^{x_t - m}
%   x_t = \rho x_{t-1} + \varepsilon_t
%   \varepsilon_t \sim N(0,\sigma^2)
%
% is X^{\star}. This constant is such that the unconditional expectation of X_t is equal
% to the deterministic steady state of X_t
%
% AUTHOR(S) 
%  stephane DOT adjemian AT univ DASH lemans DOT fr
%  frederic DOT karame AT univ DASH evry DOT fr

m = sigma/(1-autoregressive_parameter*autoregressive_parameter);