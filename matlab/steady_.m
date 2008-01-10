
function steady_()

% function steady_()
% computes the steady state 
%  
% INPUTS
%   none
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none
%  
% part of DYNARE, copyright Dynare Team (2001-2007)
% Gnu Public License.

  global M_ oo_ it_ options_
  
  
  if options_.steadystate_flag
    [oo_.steady_state,check] = feval([M_.fname '_steadystate'],...
				     oo_.steady_state,...
				     [oo_.exo_steady_state; ...
		                      oo_.exo_det_steady_state]);
    % Check if the steady state obtained from the _steadystate file is a 
    % steady state.
    check1 = 0;
    if isfield(options_,'unit_root_vars')
      if isempty(options_.unit_root_vars)
	check1 = max(abs(feval([M_.fname '_static'],...
			       oo_.steady_state,...
			       [oo_.exo_steady_state; ...
		    oo_.exo_det_steady_state]))) > options_.dynatol ;
	if check1
	  error(['The seadystate values returned by ' M_.fname ...
		'_steadystate.m don''t solve the static model!' ])
	end
      end
    end
    if ~isempty(options_.steadystate_partial)
      ssvar = options_.steadystate_partial.ssvar;
      nov   = length(ssvar);
      indv  = zeros(nov,1);
      for i = 1:nov
	indv(i) = strmatch(ssvar(i),M_.endo_names,'exact');
      end
      [oo_.steady_state,check] = dynare_solve('restricted_steadystate',...
					      oo_.steady_state(indv),...
					      options_.jacobian_flag, ...	    
					      [oo_.exo_steady_state;oo_.exo_det_steady_state],indv);
    end
  else
    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
				     oo_.steady_state,...
				     options_.jacobian_flag, ...	    
			             [oo_.exo_steady_state; ...
		                      oo_.exo_det_steady_state]);
  end

  if check ~= 0
    error('STEADY: convergence problems')
  end