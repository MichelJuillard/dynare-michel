% homotopy1 implements a homotopy method with a fixed number of intervals

function homotopy1(params,exo,exodet, step_nbr)
  global M_ oo_ options_
  
  options_.jacobian_flag = 1;

  np = size(params,1);
  ip = zeros(np,1);
  vp = zeros(np,step_nbr+1);

  nx = size(exo,1);
  ix = zeros(nx,1);
  vx = zeros(nx,step_nbr+1);
 
  nxd = size(exodet,1);
  ixd = zeros(nxd,1);
  vxd = zeros(nxd,step_nbr+1);
  
  for i = 1:np
    temp1 = strmatch(params{i,1},M_.param_names,'exact');
    if isempty(temp1)
      error(['HOMOTOPY: unknown parameter name: ' params{i,1}])
    end
    ip(i) = temp1;
    vp(i,:) = params{i,2}:(params{i,3}-params{i,2})/step_nbr:params{i,3};
  end
  
  for i = 1:nx
    temp1 = strmatch(exo{i,1},M_.exo_names,'exact');
    if isempty(temp1)
      error(['HOMOTOPY: unknown exogenous variable name: ' exo{i,1}])
    end
    ix(i) = temp1;
    vp(i,:) = exo{i,2}:(exo{i,3}-exo{i,2})/step_nbr:exo{i,3};
  end
  
  for i = 1:nxd
    temp1 = strmatch(exodet{i,1},M_.exodet_names,'exact');
    if isempty(temp1)
      error(['HOMOTOPY: unknown deterministic exogenous name: ' exodet{i,1}])
    end
    ixd(i) = temp1;
    vxd(i,:) = exodet{i,2}:(exodet{i,3}-exodet{i,2})/step_nbr:exodet{i,3};
  end
  
  for i=1:step_nbr+1
    M_.params(ip) = vp(:,i);
    for j=1:np
      assignin('base',params{j,1},vp(j,i));
    end
    
    oo_.exo_steady_state(ix) = vx(:,i);
    oo_.exo_det_steady_state(ixd) = vxd(:,i);

    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
					    oo_.steady_state,...
					    options_.jacobian_flag, ...	    
					    [oo_.exo_steady_state; ...
		    oo_.exo_det_steady_state]);
  
    if check
      error('HOMOTOPY didn''t succeed')
    end
    steady_;
  end