function [fit,SE,par,V,COMP,o1,y0]= dlr0(y,z,IRW,nvr,alpha,P0,x0,smooth,ALG,p)
% Version for use with SDP01 for fast single parameter estimation

% DLR  Dynamic Linear Regression analysis
%
% [fit,fitse,par,parse,comp,e,y0]=dlr(y,z,TVP,nvr,alpha,P0,x0,sm,ALG)
%                                     1 2  3   4    5   6  7  8   9
%
% y: Time series (*)
% z: Regressors (*)
% TVP: Model type for each TVP (0-RW/AR(1), 1-IRW/SRW) (0)
% nvr: NVR hyper-parameters (0)
% alpha: alpha hyper-parameters for SRW model (1)
% P0: Initial P matrix (1e5)
% x0: Initial state vector (0)
% sm: Smoothing on (1-default) or off (0-saves memory)
% ALG: Smoothing algorithm: P (0) or Q (1-default)
%
% fit: Model fit
% fitse: Standard error of the fit
% par: Parameter estimates
% parse: Standard errors of parameters
% comp: Linear components
% e: Normalised innovations; use e=e(~isnan(e)) to remove NaNs
% y0: Interpolated data
%
% Example: dlr(y, [ones(size(u)) u], [0 1], 0.001, [0.95 1])
%   regression type model y = c1 + c2(t)*u, with an AR(1) model for
%   c1 (alpha=0.95) and an IRW model for c2 (NVR=0.001 in both cases)
%
% See also DLROPT, FCAST, STAND, DHR, DHROPT, SDP

% Copyright (c) 2004 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Revision history :
%   04/03/2003, JT, cosmetic
%   22/08/2002, JT, interpolated data returned
%   18/01/1999, JT, revised output arguments
%   02/11/1998, JT, fit output argument
%   23/09/1998, JT, separate y and regressors input arguments, Q/P algorithm
%   22/09/1998, JT, SRW option, initial P0 not overwritten
%   01/10/1997, DP, univariate toolbox

%      Auxiliary functions:
%           dlr_opt:
%           dlrfun:

y=y(:);  % ensure y is a column, 28/6/99, JT

Z=[y, z];  %23/9/98, JT
[n, m]= size(Z); 

% Checking inputs

if nargin<10, p=[]; end

if nargin<9, ALG= 1; end  %23/9/98, JT
if isempty(ALG), ALG= 1; end


if nargin<8, smooth= 1; end
if isempty(smooth), smooth=1; end

if nargin<7, x0= []; end

if nargin<6, P0= []; end
if ~isempty(P0), lP0= length(P0); 
  P0= P0(:)'; P0= [P0 ones(1, m-1-lP0)*P0(lP0)]; 
end

if nargin<5, alpha= 1; end  %JT, 22/9/98, SRW option
if isempty(alpha), alpha= 1; end

if any(size(P0))==1
P0= diag(P0); 
end


% Checking missing values
I= zeros(n, 1); maxp= m-1; 
Miss= I; 

% SS F,  Q and P0 matrices
if IRW
   F=[1 1 ; 0 1]; Q=[ 0 0 ; 0 nvr];
else
   F=1; Q=nvr;
end


[row, col]= size(F); 

if isempty(P0)  %JT, 22/09/98, don't overwrite input argument
    P0= eye(row)*1e5;
end

% Filtering and Smoothing
[Xkk, Pkk, er, ykk1, Xkk1, yp, ep, ykk, XkN, ys, PkN, ers]=ksirw2c(Z, IRW, Miss, I, F, P0, Q, 0, x0, [], [], ALG); 
Xkk= XkN; Pkk= PkN; er= ers;

% Building outputs (linear components,  parameters,  ...)
i= cumsum([1 IRW(1:maxp-1)+1]); 

% TVP parameters
par= Xkk(i, :)'; 

% Linear components
COMP= par.*Z(:, 2:m); 

% Normalised innovations and SE of predictions
o1= [ones(row,1)*nan; (Z(row+1:n, 1)-ykk1(row+1:n))./sqrt(er(row+1:n))];   % Normalised innovations
o11= o1(~isnan(o1));
sigma= o11'*o11/(length(o11));                     % Innovations variance
SE= sqrt(sigma.*er);                               % SE of innovations

% Standard error of parameters
Pkk= Pkk'; 
W= Pkk*sigma;  
V= zeros(n, m-1); [rP, cP]= size(Pkk);
for j= 1:length(i)
  V(:, j)= sqrt(W(i(j):cP:rP, i(j))); 
end

% Fit, 2/10/98, JT
if size(COMP, 2)>1
   fit=sum(COMP')';
else
   fit=COMP;
end

% interpolated data, 22/08/02, JT
ii=find(isnan(y));
y0=y;
y0(ii)=fit(ii);

% end of m-file
