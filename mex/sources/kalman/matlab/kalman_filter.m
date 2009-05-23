%
% SYNOPSIS
%
% [loglik,per,d] = kalman_filter(Z,H,T,R,Q,Y,a,P)
% [loglik,per,d] = kalman_filter(Z,H,T,R,Q,Y,a,P,flag)
% [loglik,per,d] = kalman_filter(Z,H,T,R,Q,Y,a,Pstar,Pinf)
% [loglik,per,d] = kalman_filter(Z,H,T,R,Q,Y,a,Pstar,Pinf,flag)
%
% SEMANTICS
%
% The first two commands run a Kalman filter for non-diffuse
% initial conditions, the other two for diffuse initial conditions.
%
% Input:
%        Z,H,T,R,Q   gives a state space form
%        Y           observed data (columns correspond to periods)
%        a           mean of initial state
%        P           covariance of initial non-diffuse state
%        Pstar       finite part of covariance of initial diffuse state
%        Pinf        infinite part of covariance of initial diffuse state
%        flag        string starting with 'u', or 'U' runs a univariate
%                    form of the filter; if omitted, a multivariate version
%                    is run by default
% 
% Output:
%        loglik      data log likelihood
%        per         number of succesfully filtered periods; if no error
%                    then per equals to the number of columns of Y
%        d           number of initial periods for which the state is
%                    still diffuse (d is always 0 for non-diffuse case)
%
% Copyright 2005, Ondra Kamenik
% 

function [loglik, per, d] = kalman_filter(varargin)

  [loglik, per, d] = kalman_filter_(varargin{:});
