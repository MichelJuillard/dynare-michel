%
% [alpha,epsilon,eta] = dynare_smoother(Z,H,T,R,Q,Y,Pstar,Pinf)
% 
% This is just an interface to DiffuseKalmanSmootherH1 of Dynare. It
% takes state space in the form
%       y_t         = Z*alpha_t + epsilon_t
%       alpha_{t+1} = T*alpha_t + R*eta_t
% where epsilon covariance is H, eta covariance is Q
% 
% It returns smoothed alpha, epsilon and eta.
%
% Copyright 2005, Ondra Kamenik

% $Id: dynare_smoother.m 534 2005-11-30 13:58:11Z kamenik $

function [alpha,epsilon,eta] = dynare_smoother(Z,H,T,R,Q,Y,Pstar,Pinf)
  global options_
  
  pp = size(Z,1);
  mm = size(T,1);
  rr = size(R,2);
  dT = [zeros(pp,pp) Z*T; zeros(mm,pp) T];
  dR = [eye(pp) Z*R; zeros(mm,pp) R];
  dQ = [zeros(pp,pp) zeros(pp,rr); zeros(rr,pp) Q];
  dPinf = [Z*Pinf*Z' Z*Pinf; Pinf*Z' Pinf];
  dPstar = [Z*Pstar*Z' Z*Pstar; Pstar*Z' Pstar];
  
  mf = [1:pp];
  options_.kalman_tol = 1e-10;
  options_.nk = 0;
% if you want DiffuseKalmanSmootherH3, uncomment the following and set
% diffuse_d (possibly empty [])
% options_.diffuse_d = [7];
  
  [alpha,epsilon,eta] = DiffuseKalmanSmootherH1(dT,dR,dQ,H,dPinf,dPstar,Y,zeros(pp,size(Y,2)),pp,mm+pp,size(Y,2),mf);
  alpha = alpha(pp+1:end,:);
  eta = eta(pp+1:end,:);
