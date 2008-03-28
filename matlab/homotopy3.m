% homotopy3 implements a homotopy method that reduces the step as much as necessary

function homotopy3(params,exo,exodet, step_nbr)
  global M_ oo_ options_
  
  options_.jacobian_flag = 1;
  np = length(param_names);
  ip = zeros(np,1);
  oldvalues = zeros(np,1);
  iplus = [];
  iminus = [];
  for i = 1:np
    temp1 = strmatch(param_names{i},M_.param_names,'exact');
    if isempty(temp1)
      error(['HOMOTOPY: unknown parameter name: ' param_names{i}])
    end
    ip(i) = temp1;
    oldvalues(i) = param_values{i}(1);
    targetvalues(i) = param_values{i}(2);
    if targetvalues(i) > oldvalues(i)
      iplus = [iplus i];
    else
      iminus = [iminus i];
    end
  end
  
  iter = 1
  maxiter = 500;
  values = oldvalues;
  inc = (targetvalues-oldvalues)/2;
  k = [];
  old_ss = oo_.steady_state;
  while iter < maxiter
    for j=1:np
      M_.params(ip(j)) = values(j,i);
      assignin('base',param_names{1},values(j,1));
    end
    
    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
					    oo_.steady_state,...
					    options_.jacobian_flag, ...	    
					    [oo_.exo_steady_state; ...
		    oo_.exo_det_steady_state]);
  
    if check
      inc = inc/2;
      oo_.steady_state = old_ss;
    else
      if length(k) == np
	return
      end
      oldvalues = values;
      inc = 2*inc;
    end
    values = oldvalues + inc
    k = find(values(iplus) > targetvalues(iplus));
    values(k) = targetvalues(k);
    k = find(values(iminus) < targetvalues(iminus));
    values(k) = targetvalues(k);
    values
    if max(abs(inc)) < 1e-8
        error('HOMOTOPY didn''t succeed')
    end
  end
  error('HOMOTOPY didn''t succeed')