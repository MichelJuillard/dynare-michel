function [out,er,o1]= dlrfun(b0, X, nvr, RWtype, Miss, F, P0, steps, ALG, p, stop_info)
%DLRFUN   Likelihood function for DLR hyper-parameter estimation.
%         out= dlrfun(b0, X, nvr, RWtype, Miss, F, P0, steps, ALG, p, stop_info)
%
%      out: Value of the likelihood function.
%
%      b0: Alphas and/or NVR (transformed) parameters for evaluation of the likelihood.
%      X: Time series (output,   inputs).
%      nvr: Constrained estimation (number of inputs vector). Options:
%              -2: Free estimation of nvr parameters (default).
%              -1: Constrained estimation (all the nvrs at locations
%                  where nvr is -1 will be the same).
%              any value > 0: Nvr constrained to such value (not estimated).
%      RWtype: RW type of TVP parameter (0-RW,  1-IRW).
%      Miss: Location of missing values
%      F: SS transition matrix. When SRW models required, the diagonal of F may contain:
%              -2: Free estimation of alpha parameters.
%              -1: Constrained estimation (all the alphas at locations
%                  where F(i, i) is -1 will be the same).
%              any value > 0: Alpha constrained to such value (default: 1).
%      P0: Initial P matrix.
%      steps: Number of steps ahead forecasts for objective (0: ML)
%      ALG: Optimisation algorithm: 0=fmins, 1=fminu, 2=leastsq (not ML) (0).
%      p: extra parameter for DAR modelling
%      stop_info: graphics handle for optimisation window

% Copyright (c) 2004 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

[n, cX]= size(X); 

NVR= 10.^(b0);

% Constructing nvr vector
x= NVR; 
i1= min(find(nvr== -1)); i2= sum(nvr(1:i1)== -2); 
iii= find(nvr<0); jjj= find(nvr>-1); 
if ~isempty(i1), rest= i2+1; else, rest= 0; end
j= 1; 
for i= 1:cX-1
    if nvr(i)== -2
        if j== rest, j= j+1; end
        nvr(i)= x(j); j= j+1; 
    end
    if nvr(i)== -1, nvr(i)= x(rest); end
end

if ~isempty(iii), x(iii)= 10.^(nvr(iii)); end
if ~isempty(jjj), x(jjj)= nvr(jjj); end
x= x(:); 

if RWtype; 
    Q=[0 0;0 NVR]; 
else
    Q=NVR;
end

[Xkk, Pkk, er, ykk1, Xkk1, yp, ep, ykk]= ...
    ksirw2c(X, RWtype, Miss, zeros(size(Miss)), F, P0, Q, 0, [], steps); 
ind= ((length(F)+2):length(er))';
nd=length(ind);
out= (X(ind, 1)-ykk1(ind))./sqrt(er(ind)); 
out= sum(log(er(ind)))+nd*log((out'*out)/nd);

if nargout>1, o1= X(ind, 1)-ykk1(ind); er=er(ind); end

%JT, 9/4/99
if stop_info
    if steps==0  % ML Log-Likelihood
        set(stop_info, 'String', num2str(-0.5*out-size(X, 1)/2*(log(2*pi)+1)))
    elseif ALG==2
        set(stop_info, 'String', num2str(out'*out))  % display scalar output
    else
        set(stop_info, 'String', num2str(out))    
    end  
    drawnow
end

% end of m-file