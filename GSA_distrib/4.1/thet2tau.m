function tau = thet2tau(params, indx, indexo, flagmoments,mf,nlags,useautocorr)
global M_ oo_ options_

if nargin==1,
    indx = [1:M_.param_nbr];
    indexo = [];
end

if nargin<4,
  flagmoments=0;
end
if nargin<7 | isempty(useautocorr),
  useautocorr=0;
end

M_.params(indx) = params(length(indexo)+1:end);
if ~isempty(indexo)
  M_.Sigma_e(indexo,indexo) = diag(params(1:length(indexo)).^2);
end
% [A(oo_.dr.order_var,oo_.dr.order_var),B(oo_.dr.order_var,:)]=dynare_resolve;
[A,B]=dynare_resolve;
if flagmoments==0,
tau = [A(:); vech(B*M_.Sigma_e*B')];
else
GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.qz_criterium,options_.lyapunov_complex_threshold);
k = find(abs(GAM) < 1e-12);
GAM(k) = 0;
if useautocorr,
  sy = sqrt(diag(GAM));
  sy = sy*sy';
  sy0 = sy-diag(diag(sy))+eye(length(sy));
  dum = GAM./sy0;
  tau = vech(dum(mf,mf));
else
  tau = vech(GAM(mf,mf));
end
for ii = 1:nlags
  dum = A^(ii)*GAM;
  if useautocorr,
    dum = dum./sy;
  end
  tau = [tau;vec(dum(mf,mf))];
end
end