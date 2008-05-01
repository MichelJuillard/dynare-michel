function yf=forcst2a(y0,dr,e)
  global M_ options_
  
  Sigma_e_ = M_.Sigma_e;
  endo_nbr = M_.endo_nbr;
  exo_nbr = M_.exo_nbr;
  ykmin_ = M_.maximum_endo_lag;
  
  horizon = size(e,1);
  options_ = set_default_option(options_,'simul_seed',0);
  order = options_.order;

  k1 = [ykmin_:-1:1];
  k2 = dr.kstate(find(dr.kstate(:,2) <= ykmin_+1),[1 2]);
  k2 = k2(:,1)+(ykmin_+1-k2(:,2))*endo_nbr;

  yf = zeros(horizon+ykmin_,endo_nbr);
  yf(1:ykmin_,:) = y0';
  
  j = ykmin_*endo_nbr;
  for i=ykmin_+(1:horizon)
    tempx = yf(k1,:)';
    yf(i,:) = tempx(k2)'*dr.ghx';
    k1 = k1+1;
  end
  
   yf(:,dr.order_var) = yf;
   
