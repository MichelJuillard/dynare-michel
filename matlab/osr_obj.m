% the beginning and the end of this function may be adapted by the userx
function [z,vx]=osr_obj(x,params,weights);
  global M_ oo_ optimal_Q_ it_
  
  % set parameters of the policiy rule
  np = size(params,1);
  for i=1:np
    assignin('base',deblank(params(i,:)),x(i))
  end
  
  % don't change below until the part where the loss function is computed
  it_ = M_.maximum_lag+1;
  oo_.dr = resol(oo_.steady_state,1,0,1);
  nstatic = oo_.dr.nstatic;
  npred = oo_.dr.npred;
  ghx = oo_.dr.ghx;
  ghu = oo_.dr.ghu;
  order = oo_.dr.order_var;
  k=[nstatic+1:nstatic+npred]';
  vx1 = ghu(k,:)*M_.Sigma_e*ghu(k,:)';
  
  % compute variance of predetermined variables
  vx = (eye(npred*npred)-kron(ghx(k,:),oo_.dr.ghx(k,:)))\vx1(:);
  vx=reshape(vx,npred,npred);

  % compute variance of all variables
  if M_.endo_nbr > npred
    qx = eye(npred);
    qu = zeros(npred,M_.exo_nbr);
    if nstatic > 0
      qx = [ghx(1:nstatic,:);qx];
      qu = [ghu(1:nstatic,:);qu];
    end
    if M_.endo_nbr > nstatic+npred
      qx = [qx;ghx(nstatic+npred+1:end,:)];
      qu = [qu;ghu(nstatic+npred+1:end,:)];
    end
    vx = qx*vx*qx'+qu*M_.Sigma_e*qu';
  end
  % end of the non touch region of the program
  
  % computes the loss function
  weights = weights(order,order);
  z = weights(:)'*vx(:);
  














