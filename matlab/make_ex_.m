
function make_ex_
% function make_ex_
% forms oo_.exo_simul and oo_.exo_det_simul
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
% SPECIAL REQUIREMENTS
%   none
%  
%  
% part of DYNARE, copyright Dynare Team (1996-2007)
% Gnu Public License.
  
  global M_ options_ oo_ ex0_ ex_det0_
  
  options_ = set_default_option(options_,'periods',0);
  
  if isempty(oo_.exo_steady_state)
    oo_.exo_steady_state = zeros(M_.exo_nbr,1);
  end
  if M_.exo_det_nbr > 1 & isempty(oo_.exo_det_steady_state)
    oo_.exo_det_steady_state = zeros(M_.exo_det_nbr,1);
  end
  if isempty(oo_.exo_simul)
    if isempty(ex0_)
      oo_.exo_simul = [ones(M_.maximum_lag+options_.periods+M_.maximum_lead,1)*oo_.exo_steady_state'];
    else
      oo_.exo_simul = [ones(M_.maximum_lag,1)*ex0_';ones(options_.periods+M_.maximum_lead,1)*oo_.exo_steady_state'];
    end
  elseif size(oo_.exo_simul,1) < M_.maximum_lag+M_.maximum_lead+options_.periods
    oo_.exo_simul = [oo_.exo_simul; ones(M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.exo_simul,1),1)*oo_.exo_steady_state'];
  end

  if M_.exo_det_nbr > 0
    if isempty(oo_.exo_det_simul)
      if isempty(ex_det0_)
	oo_.exo_det_simul = [ones(M_.maximum_lag+options_.periods+M_.maximum_lead,1)*oo_.exo_det_steady_state'];
      else
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*ex_det0_';ones(options_.periods+M_.maximum_lead,1)*oo_.exo_det_steady_state'];
      end
    elseif size(oo_.exo_det_simul,1) < M_.maximum_lag+M_.maximum_lead+options_.periods
      oo_.exo_det_simul = [oo_.exo_det_simul; ones(M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.exo_det_simul,1),1)*oo_.exo_det_steady_state'];
    end
  end
    
	     