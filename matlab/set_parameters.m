function set_parameters(xparam1)
  global estim_params_ M_
  
  nvx = estim_params_.nvx;
  ncx = estim_params_.ncx;
  np = estim_params_.np;
  Sigma_e = M_.Sigma_e;
  offset = 0;
  if nvx
    offset = offset + nvx;
    var_exo = estim_params_.var_exo;
    for i=1:nvx
      k = var_exo(i,1);
      Sigma_e(k,k) = xparam1(i)^2;
    end
  end
  
  if ncx
    offset = offset + estim_params_.nvn;
    corrx = estim_params_.corrx;
    for i=1:ncx
      k1 = corrx(i,1);
      k2 = corrx(i,2);
      Sigma_e(k1,k2) = xparam1(i+offset)*sqrt(Sigma_e(k1,k1)*Sigma_e(k2,k2));
      Sigma_e(k2,k1) = Sigma_e(k1,k2);
    end
  end
  
  if np
    offset = offset+estim_params_.ncx+estim_params_.ncn;
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
  end
  
  M_.Sigma_e = Sigma_e;