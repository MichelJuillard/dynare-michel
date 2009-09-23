function [JJ, H, A, B, GAM] = getJJ(M_,oo_,options_,kronflag,indx,indexo,mf,nlags,useautocorr)

if nargin<5 | isempty(indx), indx = [1:M_.param_nbr];, end,
if nargin<6 | isempty(indexo), indexo = [];, end,
if nargin<8 | isempty(nlags), nlags=3; end,
if nargin<9 | isempty(useautocorr), useautocorr=0; end,

  if useautocorr,
    warning('off','MATLAB:divideByZero')
  end
if kronflag == -1,
  fun = 'thet2tau';
  params0 = M_.params;
  JJ = fdjac(fun,[sqrt(diag(M_.Sigma_e(indexo,indexo))); M_.params(indx)],indx,indexo,1,mf,nlags,useautocorr);
  M_.params = params0;
  assignin('base','M_', M_);
  assignin('base','oo_', oo_);
else
  [H, A, B, dA, dOm, info] = getH(M_,oo_,kronflag,indx,indexo);
  if info(1) > 0
    JJ = [];
    GAM = [];
    return
  end
  m = length(A);

  GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.qz_criterium,options_.lyapunov_complex_threshold,1);
  k = find(abs(GAM) < 1e-12);
  GAM(k) = 0;
  if useautocorr,
  sdy = sqrt(diag(GAM));
  sy = sdy*sdy';
  end
  
%   BB = dOm*0;
%   for j=1:length(indx),
%     BB(:,:,j)= dA(:,:,j)*GAM*A'+A*GAM*dA(:,:,j)'+dOm(:,:,j);
%   end
%   XX =  lyapunov_symm_mr(A,BB,options_.qz_criterium,options_.lyapunov_complex_threshold,0);
  for j=1:length(indexo),
    dum =  lyapunov_symm(A,dOm(:,:,j),options_.qz_criterium,options_.lyapunov_complex_threshold,2);
%     dum =  XX(:,:,j);
    k = find(abs(dum) < 1e-12);
    dum(k) = 0;
    if useautocorr
      dsy = 1/2./sdy.*diag(dum);
      dsy = dsy*sdy'+sdy*dsy';
      dum1=dum;
      dum1 = (dum1.*sy-dsy.*GAM)./(sy.*sy);
      dum1 = dum1-diag(diag(dum1))+diag(diag(dum));
      dumm = vech(dum1(mf,mf));
    else
      dumm = vech(dum(mf,mf));
    end
    for i=1:nlags,
      dum1 = A^i*dum;
      if useautocorr
        dum1 = (dum1.*sy-dsy.*(A^i*GAM))./(sy.*sy);
      end
      dumm = [dumm; vec(dum1(mf,mf))];
    end
    JJ(:,j) = dumm;
  end
  nexo = length(indexo);
  for j=1:length(indx),
    dum =  lyapunov_symm(A,dA(:,:,j+nexo)*GAM*A'+A*GAM*dA(:,:,j+nexo)'+dOm(:,:,j+nexo),options_.qz_criterium,options_.lyapunov_complex_threshold,2);
%     dum =  XX(:,:,j);
    k = find(abs(dum) < 1e-12);
    dum(k) = 0;
    if useautocorr
      dsy = 1/2./sdy.*diag(dum);
      dsy = dsy*sdy'+sdy*dsy';
      dum1=dum;
      dum1 = (dum1.*sy-dsy.*GAM)./(sy.*sy);
      dum1 = dum1-diag(diag(dum1))+diag(diag(dum));
      dumm = vech(dum1(mf,mf));
    else
      dumm = vech(dum(mf,mf));
    end
    for i=1:nlags,
      dum1 = A^i*dum;
      for ii=1:i,
        dum1 = dum1 + A^(ii-1)*dA(:,:,j+nexo)*A^(i-ii)*GAM;
      end
      if useautocorr
        dum1 = (dum1.*sy-dsy.*(A^i*GAM))./(sy.*sy);
      end
      dumm = [dumm; vec(dum1(mf,mf))];
    end
    JJ(:,j+nexo) = dumm;
  end
  
  if nargout >4,
    [GAM,stationary_vars] = th_autocovariances(oo_.dr,oo_.dr.order_var(mf),M_,options_);
  end
end

  if useautocorr,
    warning('on','MATLAB:divideByZero')
  end