function rootsinv = fn_varoots(Bhat,nvar,lags)
%
%    Using eigenvalues to find the inverse of all roots associated with the VAR proceess:
%          y_t' = C + y_{t-1}'*B_1 + ... + Y_{t-p}'*B_p + u_t'.
%    where columns correspond to equations.  See also Judge (1), pp.753-755 where rows correspond to equations.
% Bhat:  ncoef-by-nvar where ncoef=nvar*lags+nexo and nvar is the number of endogenous variables.
%    Columns corresponds to equations with
%    ncoef=[nvar for 1st lag, ..., nvar for last lag, other exogenous terms, const term]
%                       ..., nvar coef in the last lag, and nexo coefficients.
%    Note that entries in the rows of Bhat that > nvar*lags are irrelevant.
% nvar: number of endogenous variables.
% lags: number of lags.
%-------
% rootsinv:  a vector of nvar*lags inverse roots.  When > 1, explosive.  When all < 1, stationary.
%
% Tao Zha, September 2000


if size(Bhat,1)<nvar*lags
   disp(' ')
   warning('Make sure that Bhat has at least nvar*lags rows')
   return
end

%--------- Strack the VAR(p) to the VAR(1) with z_t = Az_{t-1}.
%
A1 = diag(ones(nvar*(lags-1),1));
A2 = [A1 zeros(nvar*(lags-1),nvar)];
A = [Bhat(1:nvar*lags,:)'; A2];
rootsinv=eig(A);
