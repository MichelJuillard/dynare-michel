function homotopy3(values, step_nbr)
% function homotopy3(values, step_nbr)
%
% Implements homotopy (mode 3) for steady-state computation.
% Tries first the most extreme values. If it fails to compute the steady
% state, the interval between initial and desired values is divided by two
% for each parameter. Every time that it is impossible to find a steady
% state, the previous interval is divided by two. When one succeed to find
% a steady state, the previous interval is multiplied by two.
%
% INPUTS
%    values:        a matrix with 4 columns, representing the content of
%                   homotopy_setup block, with one variable per line.
%                   Column 1 is variable type (1 for exogenous, 2 for
%                   exogenous deterministic, 4 for parameters)
%                   Column 2 is symbol integer identifier.
%                   Column 3 is initial value, and column 4 is final value.
%                   Column 3 can contain NaNs, in which case previous
%                   initialization of variable will be used as initial value.
%    step_nbr:      maximum number of steps to try before aborting
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.
  global M_ oo_ options_
  
  tol = 1e-8;
  
  nv = size(values,1);

  ip = find(values(:,1) == 4); % Parameters
  ix = find(values(:,1) == 1); % Exogenous
  ixd = find(values(:,1) == 2); % Exogenous deterministic

  if length([ip, ix, ixd]) ~= nv
    error('HOMOTOPY: incorrect variable types specified')
  end

  % Construct vector of starting values, using previously initialized values
  % when initial value has not been given in homotopy_setup block
  oldvalues = values(:,3);
  ipn = find(values(:,1) == 4 & isnan(oldvalues));
  oldvalues(ipn) = M_.params(values(ipn, 2));
  ixn = find(values(:,1) == 1 & isnan(oldvalues));
  oldvalues(ixn) = oo_.exo_steady_state(values(ixn, 2));
  ixdn = find(values(:,1) == 2 & isnan(oldvalues));
  oldvalues(ixdn) = oo_.exo_det_steady_state(values(ixdn, 2));

  targetvalues = values(:,4);

  if min(abs(targetvalues - oldvalues)) < tol
    error('HOMOTOPY: distance between initial and final values should be at least %e for all variables', tol)
  end
  iplus = find(targetvalues > oldvalues);
  iminus = find(targetvalues < oldvalues);
  
  curvalues = oldvalues;
  inc = (targetvalues-oldvalues)/2;
  kplus = [];
  kminus = [];
  old_ss = oo_.steady_state;

  iter = 1;
  while iter < step_nbr
    M_.params(values(ip,2)) = curvalues(ip);
    oo_.exo_steady_state(values(ix,2)) = curvalues(ix);
    oo_.exo_det_steady_state(values(ixd,2)) = curvalues(ixd);
    
    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
                                            oo_.steady_state,...
                                            options_.jacobian_flag, ...	    
                                            [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state], M_.params);
  
    if check
      inc = inc/2;
      oo_.steady_state = old_ss;
    else
      if length([kplus; kminus]) == nv
        return
      end
      oldvalues = curvalues;
      inc = 2*inc;
    end
    curvalues = oldvalues + inc;
    kplus = find(curvalues(iplus) >= targetvalues(iplus));
    curvalues(kplus) = targetvalues(kplus);
    kminus = find(curvalues(iminus) <= targetvalues(iminus));
    curvalues(kminus) = targetvalues(kminus);

    if max(abs(inc)) < tol
        error('HOMOTOPY didn''t succeed')
    end
    
    iter = iter + 1;
  end
  error('HOMOTOPY didn''t succeed')
