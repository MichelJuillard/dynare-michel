function homotopy2(values, step_nbr)
% function homotopy2(values, step_nbr)
%
% Implements homotopy (mode 2) for steady-state computation.
% Only one parameter/exogenous is changed at a time.
% Computation jumps to next variable only when current variable has been
% brought to its final value.
% Variables are processed in the order in which they appear in "values".
% The problem is solved var_nbr*step_nbr times.
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
  
  oldvalues = values(:,3);
  
  % Initialize all variables with initial value, or the reverse...
  for i = 1:nv
    switch values(i,1)
     case 1
      if isnan(oldvalues(i))
        oldvalues(i) = oo_.exo_steady_state(values(i,2));
      else
        oo_.exo_steady_state(values(i,2)) = oldvalues(i);
      end
     case 2
      if isnan(oldvalues(i))
        oldvalues(i) = oo_.exo_det_steady_state(values(i,2));
      else
        oo_.exo_det_steady_state(values(i,2)) = oldvalues(i);
      end
     case 4
      if isnan(oldvalues(i))
        oldvalues(i) = M_.params(values(i,2));
      else
        M_.params(values(i,2)) = oldvalues(i);
      end
     otherwise
      error('HOMOTOPY: incorrect variable types specified')
    end
  end

  if any(oldvalues == values(:,4))
    error('HOMOTOPY: initial and final values should be different')
  end
  
  % Actually do the homotopy
  for i = 1:nv
    for v = oldvalues(i):(values(i,4)-oldvalues(i))/step_nbr:values(i,4)
      switch values(i,1)
       case 1
        oo_.exo_steady_state(values(i,2)) = v;
       case 2
        oo_.exo_det_steady_state(values(i,2)) = v;
       case 4
        M_.params(values(i,2)) = v;
      end
      
      [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
                                              oo_.steady_state,...
                                              options_.jacobian_flag, ...	    
                                              [oo_.exo_steady_state; ...
                          oo_.exo_det_steady_state], M_.params);
      
      if check
        error('HOMOTOPY didn''t succeed')
      end
    end
  end
  