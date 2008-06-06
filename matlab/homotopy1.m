function homotopy1(values, step_nbr)
% function homotopy1(values, step_nbr)
%
% Implements homotopy (mode 1) for steady-state computation.
% The multi-dimensional vector going from the set of initial values
% to the set of final values is divided in as many sub-vectors as
% there are steps, and the problem is solved as many times.
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
%    step_nbr:      number of steps for homotopy
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

  nv = size(values, 1);
  
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
  
  if any(oldvalues == values(:,4))
    error('HOMOTOPY: initial and final values should be different')
  end
  
  points = zeros(nv, step_nbr+1);
  for i = 1:nv
    points(i,:) = oldvalues(i):(values(i,4)-oldvalues(i))/step_nbr:values(i,4);
  end
  
  for i=1:step_nbr+1
    M_.params(values(ip,2)) = points(ip,i);
    oo_.exo_steady_state(values(ix,2)) = points(ix,i);
    oo_.exo_det_steady_state(values(ixd,2)) = points(ixd,i);

    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
                                            oo_.steady_state,...
                                            options_.jacobian_flag, ...	    
                                            [oo_.exo_steady_state; ...
                        oo_.exo_det_steady_state], M_.params);
  
    if check
      error('HOMOTOPY didn''t succeed')
    end
  end
