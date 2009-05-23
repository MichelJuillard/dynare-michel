% [y,epsilon,alpha,eta] = gendata(T, ssf, a0)
%
% generates random data of the length T for the given state space form
% and initial state

% $Id: gendata.m 532 2005-11-30 13:51:33Z kamenik $
% Copyright 2005, Ondra Kamenik

function [y,epsilon,alpha,eta] = gendata(T, ssf, a0)

  m = size(ssf.T, 1);
  p = size(ssf.Z, 1);
  r = size(ssf.R, 2);
  
  cholH = chol(ssf.H);
  cholQ = chol(ssf.Q);

  epsilon = cholH*randn(p,T);
  eta = cholQ*randn(r, T);

  y = zeros(p, T);
  alpha = zeros(m,T);
  alpha(:,1) = a0;
  
  for t = 1:T
    y(:,t) = ssf.Z*alpha(:,t) + epsilon(:,t);
    if t ~= T
      alpha(:,t+1) = ssf.T*alpha(:,t) + ssf.R*eta(:,t);
    end
  end
