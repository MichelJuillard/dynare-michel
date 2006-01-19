function [H,prob,d] = smirnov(x1 , x2 , alpha )
% Smirnov test for 2 distributions
%   [H,prob,d] = smirnov(x1 , x2 , alpha )
%
% Copyright (C) 2005 Marco Ratto
%



if nargin<3
    alpha  =  0.05;
end


% empirical cdfs.

bin    =  [-inf ; sort([x1;x2]) ; inf];

ncount1  =  histc (x1 , bin);
ncount2  =  histc (x2 , bin);

cum1  =  cumsum(ncount1)./sum(ncount1);
cum2  =  cumsum(ncount2)./sum(ncount2);

cum1  =  cum1(1:end-1);
cum2  =  cum2(1:end-1);

n1     =  length(x1);
n2     =  length(x2);
n      =  n1 * n2 /(n1 + n2);

% Compute the d(n1,n2) statistic.

d  =  max(abs(cum1 - cum2));

%
% Compute P-value check H0 hypothesis
%

lam =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * d , 0);

j       =  [1:101]';
prob  =  2 * sum((-1).^(j-1).*exp(-2*lam*lam*j.^2));

if prob < 0 , prob = 0; end
if prob > 1 , prob = 1; end


H  =  (alpha >= prob);
