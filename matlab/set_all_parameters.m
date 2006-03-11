function set_all_parameters(xparam1)
  global estim_params_ M_
  
  nvx = estim_params_.nvx;
  ncx = estim_params_.ncx;
  nvn = estim_params_.nvn;
  ncn = estim_params_.ncn;
  np = estim_params_.np;
  Sigma_e = M_.Sigma_e;
  H = M_.H;

  if nvx
    var_exo = estim_params_.var_exo;
    for i=1:nvx
      k =var_exo(i,1);
      Sigma_e(k,k) = xparam1(i)^2;
    end
  end
  
  if ncx
    offset = nvx+nvn;
    corrx = estim_params_.corrx;
    for i=1:ncx
      k1 = corrx(i,1);
      k2 = corrx(i,2);
      Sigma_e(k1,k2) = xparam1(i+offset)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
      Sigma_e(k2,k1) = Sigma_e_(k1,k2);
    end
  end
  
  if nvn
    offset = nvx;
    var_endo = estim_params_.var_endo;
    for i=1:nvn
      k = var_endo(i,1);
      H(k,k) = xparam1(i+offset)^2;
    end
  end
  
  if ncx
    offset = nvx+nvn+ncx;
    corrn = estim_params_.corrn;
    for i=1:ncx
      k1 = corr(i,1);
      k2 = corr(i,2);
      H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
      H(k2,k1) = H(k1,k2);
    end
  end
  
  if np
    offset = nvx+ncx+nvn+ncn;
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
  end
  
  if nvx
    M_.Sigma_e = Sigma_e;
  end

  M_.H = H;
