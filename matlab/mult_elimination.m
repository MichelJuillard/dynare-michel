function dr=mult_elimination(void)
% function mult_elimination()
% replaces Lagrange multipliers in Ramsey policy by lagged value of state
% and shock variables
% INPUT
%   none  
% OUTPUT
%   dr: a structure with the new decision rule
%
  
  global M_ options_ oo_

  dr = oo_.dr;
  
  nstatic = dr.nstatic;
  npred = dr.npred;
  order_var = dr.order_var;
  nstates = M_.endo_names(order_var(nstatic+(1:npred)),:);
  
  il = strmatch('mult_',nstates);
  nil = setdiff(1:dr.npred,il);
  m_nbr = length(il);
  nm_nbr = length(nil);
  
  AA1 = dr.ghx(:,nil);
  AA2 = dr.ghx(:,il);
  A1 = dr.ghx(nstatic+(1:npred),nil);
  A2 = dr.ghx(nstatic+(1:npred),il);
  B = dr.ghu(nstatic+(1:npred),:);
  A11 = A1(nil,:);
  A21 = A1(il,:);
  A12 = A2(nil,:);
  A22 = A2(il,:)

  [Q1,R1,E1] = qr(A2);
  n1 = sum(abs(diag(R1)) > 1e-8);
  
  Q1_12 = Q1(1:nm_nbr,n1+1:end);
  Q1_22 = Q1(nm_nbr+1:end,n1+1:end);
  [Q2,R2,E2] = qr(Q1_22');
  n2 = sum(abs(diag(R2)) > 1e-8);
  
  R2_1 = inv(R2(1:n2,1:n2));

  M1(order_var,:) = AA1 - AA2*E2*[R2_1*Q2(:,1:n2)'*Q1_12'; zeros(m_nbr-n2,m_nbr)];
  M2(order_var,:) = AA2*E2*[R2_1*Q2(:,1:n2)'*[Q1_12' Q1_22']*A1; zeros(m_nbr-n2,nil)];
  M3(order_var,:) = dr.ghu;
  M4(order_var,:) = AA2*E2*[R2_1*Q2(:,1:n2)'*[Q1_12' Q1_22']*B; zeros(m_nbr-n2,size(B,2))];

  endo_nbr = M_.orig_model.endo_nbr;
  exo_nbr = M_.exo_nbr;

  lead_lag_incidence = M_.lead_lag_incidence(:,1:endo_nbr+exo_nbr);
  lead_lag_incidence1 = lead_lag_incidence(1,:) > 0;
  maximum_lag = M_.maximum_lag;
  for i=1:maximum_lag-1
    lead_lag_incidence1 = [lead_lag_incidence1; lead_lag_incidence(i,:)| ...
		    lead_lag_incidence(i+1,:)];
  end
  lead_lag_incidence1 = [lead_lag_incidence1; ...
		    lead_lag_incidence(M_.maximum_lag,:) > 0];
  k = find(lead_lag_incidence1');
  lead_lag_incidence1 = zeros(size(lead_lag_incidence1'));
  lead_lag_incidence1(k) = 1:length(k);
  lead_lag_incidence1 = lead_lag_incidence1';
  
  kstate = zeros(0,2);
  for i=maximum_lag:-1:1
    k = find(lead_lag_incidence(i,:));
    kstate = [kstate; [k repmat(i+1,length(k),1)]];
  end
  
  dr.M1 = M1;
  dr.M2 = M2;
  dr.M3 = M3;
  dr.M4 = M4;