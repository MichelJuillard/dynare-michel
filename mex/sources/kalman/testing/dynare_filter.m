%
% loglik = dynare_filter(Z,H,T,R,Q,Y,Pstar,Pinf)
% 
% This is just an interface to DiffuseLikelihoodH1 of Dynare. It
% takes state space in the form
%       y_t         = Z*alpha_t + epsilon_t
%       alpha_{t+1} = T*alpha_t + R*eta_t
% where epsilon covariance is H, eta covariance is Q
% 
% It returns log likelihood.
%
% Copyright 2005, Ondra Kamenik

% $Id: dynare_filter.m 534 2005-11-30 13:58:11Z kamenik $

function lik = dynare_filter(Z,H,T,R,Q,Y,Pstar,Pinf)
  global bayestopt_ options_
  
  pp = size(Z,1);
  mm = size(T,1);
  rr = size(R,2);
  dT = [zeros(pp,pp) Z*T; zeros(mm,pp) T];
  dR = [eye(pp) Z*R; zeros(mm,pp) R];
  dQ = [zeros(pp,pp) zeros(pp,rr); zeros(rr,pp) Q];
  dPinf = [zeros(pp,pp) zeros(pp,mm); zeros(mm,pp) Pinf];
  dPstar = [Z*Pstar*Z' Z*Pstar; Pstar*Z' Pstar];
  
  bayestopt_.mf = [1:pp];
  options_.kalman_tol = 1e-10;
  
  lik = DiffuseLikelihoodH1(dT,dR,dQ,H,dPinf,dPstar,Y,zeros(pp,size(Y,2)),1);
  
  
