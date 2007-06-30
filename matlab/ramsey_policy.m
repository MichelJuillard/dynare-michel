function info = ramsey_policy(var_list)
  global options_ oo_ M_
  
  oldoptions = options_;
  options_.ramsey_policy = 1;
  options_.order = 1;
  info = stoch_simul(var_list);
  
  return
  if info
    return
  end
  
  options_.ramsey_policy = 2;
  info = stoch_simul(var_list);
  dr = oo_.dr;
  
  orig_model = M_.orig_model;
  endo_nbr1 = orig_model.endo_nbr;
  exo_nbr1 = M_.exo_nbr;
  
  npred = dr.npred;
  nspred = dr.nspred;
  inv_order_var1 = dr.inv_order_var(1:endo_nbr1);
  nstatic = dr.nstatic;
  
  ghx = dr.ghx;
  ghx1 = ghx(inv_order_var1,:);
  [ghx2,ghu2] = transition_matrix(dr);
  ghu = dr.ghu;
  ghu1 = ghu(inv_order_var1,:);
  ghxx = dr.ghxx;
  ghxx1 = ghxx(inv_order_var1,:);
  ghxx2 = [ghxx(nstatic+(1:npred),:); zeros(nspred-npred,size(ghxx,2))];
  ghxu = dr.ghxu;
  ghxu1 = ghxu(inv_order_var1,:);
  ghxu2 = [ghxu(nstatic+(1:npred),:); zeros(nspred-npred,size(ghxu,2))];
  ghuu = dr.ghuu;
  ghuu1 = ghuu(inv_order_var1,:);
  ghuu2 = [ghuu(nstatic+(1:npred),:); zeros(nspred-npred,size(ghuu,2))];
  ghs2 = dr.ghs2;  
  ghs21 = ghs2(inv_order_var1,:);
  ghs22 = [ghs2(nstatic+(1:npred),:); zeros(nspred-npred,size(ghs2,2))];
  
  fname = [M_.fname '_objective_static'];
  [ubar,uj,uh] = feval(fname,dr.ys(1:endo_nbr1),zeros(exo_nbr1,1));
  
  bet = options_.planner_discount;
  wbar = ubar/(1-bet);
  % wx = ux*gx + bet*wx*hx
  % wx*(I-bet*hx) = ux*gx
  wx = uj*ghx1/(eye(nspred)-bet*ghx2);
  % wu = ux*gu + bet*wx*hu
  wu = uj*ghu1 + bet*wx*ghu2;
  % wxx = uxx*kron(gx,gx)+ux*ghxx + bet(wxx*kron(hx,hx) +
  %       wx*hxx)
  wxx = (uh*kron(ghx1,ghx1) + uj*ghxx1 + bet*wx*ghxx2)/(eye(nspred^2)-bet*kron(ghx2,ghx2));
  wxu = uh*kron(ghx1,ghu1) + uj*ghxu1 + bet*(wx*ghxu2+wxx*kron(ghx2,ghu2));
  wuu = uh*kron(ghu1,ghu1) + uj*ghuu1 + bet*(wx*ghuu2+wxx*kron(ghu2,ghu2));
  % ws2 = ux*ghs2 + bet*(ws2+wx*hs2+wuu*Sigma_u)
  ws2 = (uj*ghs21 + bet*(wx*ghs22+wuu*M_.Sigma_e(:)))/(1-bet);
  oo_.welfare.wbar = wbar;
  oo_.welfare.wx = wx;
  oo_.welfare.wu = wu;
  oo_.welfare.wxx = wxx;
  oo_.welfare.wxu = wxu;
  oo_.welfare.wuu = wuu;
  oo_.welfare.ws2 = ws2;

  disp(' ')
  disp(' ')
  disp(['Welfare at the deterministic steady state is ' num2str(wbar+0.5* ...
						  ws2)]);
  
  options_ = oldoptions;