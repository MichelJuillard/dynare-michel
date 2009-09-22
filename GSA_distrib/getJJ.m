function JJ = getJJ(M_,oo_,options_,kronflag,indx,mf)

if kronflag == -1,
  fun = 'thet2tau';
  params0 = M_.params;
  JJ = fdjac(fun,M_.params(indx),indx,1,mf);
  assignin('base','M_', M_);
  assignin('base','oo_', oo_);
else
  [H, A, B] = getH(M_,oo_,kronflag,indx);
  m = length(A);
  Dm = duplication(m);
  dA = reshape(H(1:m^2,:),[m m length(indx)]);
  dOm = Dm*H(m^2+1:end,:);
  dOm = reshape(dOm,[m m length(indx)]);

  GAM =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.qz_criterium,options_.lyapunov_complex_threshold);

  for j=1:length(indx),
    dum =  lyapunov_symm(A,(squeeze(dA(:,:,j))*GAM*A'+A*GAM*squeeze(dA(:,:,j))'+squeeze(dOm(:,:,j))),options_.qz_criterium,options_.lyapunov_complex_threshold);
    JJ(:,j) = vech(dum(mf,mf));
  end
end