function y_=irf(dr, e1, long_, drop_, replic, iorder)
  global M_ oo_ options_

  old_iter = options_.periods;
  options_.periods = long_;
  
  temps = repmat(dr.ys,1,M_.maximum_lag);

  y_	= 0;
  if iorder == 1
    options_.periods = long_;
    y1_ = repmat(dr.ys,1,options_.periods);
    ex2_ = zeros(options_.periods,M_.exo_nbr);
    ex2_(1,:) = e1';
    y2_ = simult_(repmat(dr.ys,1,M_.maximum_lag),dr,ex2_,iorder);
    y_ = y2_(:,M_.maximum_lag:end)-y1_;% <-- y2_(:,M_.maximum_lag+1:end)-y1_
  
  else
    % eliminate shocks with 0 variance
    i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 ));
    nxs = length(i_exo_var);
    ex1_ = zeros(long_+drop_+M_.maximum_lag,M_.exo_nbr);
    ex2_ = ex1_;
    chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));

    for j = 1: replic
      randn('seed',j);
      ex1_(:,i_exo_var) = randn(long_+drop_+M_.maximum_lag,nxs)*chol_S;
      ex2_ = ex1_;
      ex2_(drop_+1,:) = ex2_(drop_+1,:)+e1';   
      y1_ = simult_(repmat(dr.ys,1,M_.maximum_lag),dr,ex1_,iorder);
      y2_ = simult_(repmat(dr.ys,1,M_.maximum_lag),dr,ex2_,iorder);
      y_ = y_+(y2_(:,M_.maximum_lag+drop_+1:end)-y1_(:,M_.maximum_lag+drop_+1:end));
    end
    y_=y_/replic;
  end
  options_.periods = old_iter;
% 01/18/02 MJ corrected for many lags
% 03/11/22 MJ input is now entire shock vector e1 (for orthogonalized IRFs)


