% function [LIK per d lik] = kalman_filters(T,R,Q,H,Y,start,Z,a,P,[Pinf |u/U flag])
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    R                      [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                      [double]    rr*rr covariance matrix of the structural innovations.
%    H                      [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    Z                      [double]    pp*mm matrix mapping state to pp observations 
%    a                      [vector]    mm vector of mean initial state, usually of 0s
%    P                      [double]    mm*mm variance-covariance matrix with stationary variables
%    Pinf   [optional]      [double]    mm*mm variance-covariance matrix with stationary variables
%    [u/U]flag   [optional] [char]      u/U univariate falg
%
% SYNOPSIS
%
% [LIK,per,d,lik] = kalman_filters(T,R,Q,H,Y,start,a, Z,P)
% [LIK,per,d,lik] = kalman_filters(T,R,Q,H,Y,start,a, Z,P,flag)
% [LIK,per,d,lik] = kalman_filters(T,R,Q,H,Y,start,a, Z,Pstar,Pinf)
% [LIK,per,d,lik] = kalman_filters(T,R,Q,H,Y,start,a, Z, Pstar, Pinf, flag)
%
% SEMANTICS
%
% The first two commands run a Kalman filter for non-diffuse initial conditions, 
% univariate or multivariate, the other two for diffuse initial conditions.
%
% 
% Output:
%        LIK      data log likelihood
%        per         number of succesfully filtered periods; if no error
%                    then per equals to the number of columns of Y
%        d           number of initial periods for which the state is
%                    still diffuse (d is always 0 for non-diffuse case)
%
% Copyright 2005, Ondra Kamenik
% 

%function [LIK per d lik] = kalman_filters(varargin)
function [LIK per d lik] = kalman_filters(T,R,Q,H,Y,start,Z,a,P,varargin)
if isempty(H)
    H=zeros(size(Y,1), size(Y,1))
elseif H==0
    H=zeros(size(Y,1), size(Y,1))
end
if isempty(a)
    a=zeros(size(T,1),1)
elseif a==0
    a=zeros(size(T,1),1)
end
if size(Z,1)== 1 && size(Z,2)==size(Y,1) && size(Y,1)> 1  
    ZM=zeros(size(Y,1), size(T,2))
    for i = 1:size(Y,1)
    ZM(i,Z(i))=1
    end
else
    ZM=Z;
end
% run basic multivariate kalman filter
[LIK per d lik] = kalman_filter_dll(T,R,Q,H,Y,start,ZM,a, P); 

