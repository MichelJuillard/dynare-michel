function [YtY,XtY,YtX,XtX,Y,X] = var_sample_moments(FirstObservation,LastObservation,qlag,var_trend_order)
% Computes the sample moments of a VAR model.  
%
% The VAR(p) model is defined by:
%
%   y_t = \sum_{k=1}^p y_{t-k} A_k + z_t C + e_t  for t = 1,...,T  
%
% where y_t is a 1*m vector of observed endogenous variables, p is the
% number of lags, A_k is an m*m real matrix, z_t is a 1*q vector of
% exogenous (deterministic) variables, C is a q*m real matrix and
% e_t is a vector of exogenous stochastic shocks. T is the number
% of observations. The deterministic exogenous variables are assumed to 
% be a polynomial trend of order q = "var_trend_ordre".  
%
% We define: 
%
%  <>  Y = (y_1',y_2',...,y_T')' a T*m matrix,
%
%  <>  x_t = (y_{t-1},y_{t-2},...,y_{t-p},z_t) a 1*(mp+q) row vector, 
%
%  <>  X = (x_1',x_2',...,x_T')' a T*(mp+q) matrix, 
%
%  <>  E = (e_1',e_2',...,e_T')' a T*m matrix and
%
%  <>  A = (A_1',A_2',...,A_p',C')' an (mp+q)*m matrix of coefficients.   
%
% So that we can equivalently write the VAR(p) model using the following
% matrix representation:
%
%   Y = X * A +E
%
%
% INPUTS 
%   o FirstObservation    [integer] First observation.
%   o LastObservation     [integer] Last observation.
%   o qlag                [integer] Number of lags in the VAR model.
%   o var_trend_order     [integer] Order of the polynomial exogenous trend: 
%                                       = -1 no constant and no linear trend,
%                                       =  0 constant and no linear trend,
%                                       =  1 constant and linear trend.
%
% OUTPUTS 
%   o YtY                 [double]  Y'*Y an m*m matrix.
%   o XtY                 [double]  X'*Y an (mp+q)*m matrix. 
%   o YtX                 [double]  Y'*X an m*(mp+q) matrix.
%   o XtX                 [double]  X'*X an (mp+q)*(mp+q) matrix.
%   o Y                   [double]  Y a T*m matrix.
%   o X                   [double]  X a T*(mp+q) matrix.
%
% ALGORITHM
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% part of DYNARE, copyright Dynare Team (2007)
% Gnu Public License.
global options_

X = [];
Y = [];
YtY = [];
YtX = [];
XtY = [];
XtX = [];

if exist(options_.datafile)
  eval(options_.datafile);
else
  eval(['load ' options_.datafile]);
end

data = [ ];
for i=1:size(options_.varobs,1)% m is equal to options_.varobs
  data = [data eval(deblank(options_.varobs(i,:)))];
end

if qlag > FirstObservation
  disp('VarSampleMoments :: not enough data to initialize! Try to increase FirstObservation.')
  return
end

NumberOfObservations = LastObservation-FirstObservation+1;% This is T.
NumberOfVariables = size(options_.varobs,1);% This is m.
if var_trend_order == -1% No constant no linear trend case.
    X = zeros(NumberOfObservations,NumberOfVariables*qlag);
elseif var_trend_order == 0% Constant and no linear trend case.
    X = zeros(NumberOfObservations,NumberOfVariables*qlag+1);
    indx = NumberOfVariables*qlag+1;
elseif var_trend_order == 1;% Constant and linear trend case.
    X = zeros(NumberOfObservations,NumberOfVariables*qlag+2);
    indx = NumberOfVariables*qlag+1:NumberOfVariables*qlag+2;
else
    disp('var_sample_moments :: trend must be equal to -1,0 or 1!')
    return
end

% I build matrices Y and X  
Y = data(FirstObservation:LastObservation,:);
for t=1:NumberOfObservations
  line = t + FirstObservation-1;
  for lag = 1:qlag
      X(t,(lag-1)*NumberOfVariables+1:lag*NumberOfVariables) = data(line-lag,:);
  end
  if var_trend_order == 0
      X(t,indx) = 1;
  elseif var_trend_order == 1
      X(t,indx) = [ 1 , t ];
  end
end

YtY = Y'*Y;
YtX = Y'*X;
XtY = X'*Y;
XtX = X'*X;