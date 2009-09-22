function JJ = getJJ(M_,oo_,options_,kronflag,indx,mf,nlags)

if nargin<5 | isempty(indx), indx = [1:M_.param_nbr];, end,
if nargin<7 | isempty(nlags), nlags=3; end,

if kronflag == -1,
  fun = 'thet2tau';
  params0 = M_.params;
  JJ = fdjac(fun,M_.params(indx),indx,1,mf,nlags);
  assignin('base','M_', M_);
  assignin('base','oo_', oo_);
else
  [H, A, B, dA, dOm] = getH(M_,oo_,kronflag,indx);
  m = length(A);
  Dm = duplication(m);
  dA = reshape(H(1:m^2,:),[m m length(indx)]);
  dOm = Dm*H(m^2+1:end,:);
  dOm = reshape(dOm,[m m length(indx)]);

  GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.qz_criterium,options_.lyapunov_complex_threshold);
  
  for j=1:length(indx),
    dum =  lyapunov_symm(A,dA(:,:,j)*GAM*A'+A*GAM*dA(:,:,j)'+dOm(:,:,j),options_.qz_criterium,options_.lyapunov_complex_threshold);
    dumm = vech(dum(mf,mf));
    for i=1:nlags,
      dum1 = A^i*dum;
      for ii=1:i,
        dum1 = dum1 + A^(ii-1)*dA(:,:,j)*A^(i-ii)*GAM;
      end
      dumm = [dumm; vec(dum1(mf,mf))];
    end
    JJ(:,j) = dumm;
  end
end