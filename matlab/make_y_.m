% Copyright (C) 2001 Michel Juillard
%
function make_y_
  global M_ options_ oo_ ys0_ 
  
  options_ = set_default_option(options_,'periods',0);
  
  if isempty(oo_.steady_state)
    oo_.steady_state = ones(M_.endo_nbr,1);
  end
  
  
  if isempty(oo_.y_simul)
    if isempty(ys0_)
      oo_.y_simul = [oo_.steady_state*ones(1,M_.maximum_lag+options_.periods+M_.maximum_lead)];
    else
      oo_.y_simul = [ys0_*ones(1,M_.maximum_lag);oo_.steady_state*ones(1,options_.periods+M_.maximum_lead)];
    end
  elseif size(oo_.y_simul,1) < length(oo_.steady_state)
    k = size(oo_.y_simul,1)+1:length(oo_.steady_state)
    if isempty(ys0_)
      oo_.y_simul = [oo_.y_simul; oo_.steady_state(k)*ones(1,M_.maximum_lag+size(oo_.y_simul,1)+M_.maximum_lead)];
    else
      oo_.y_simul = [oo_.y_simul; [ys0_(k)*ones(1,M_.maximum_lag); oo_.steady_state(k)*ones(1,size(oo_.y_simul,2)-M_.maximum_lag+ ...
						     M_.maximum_lead)]];
    end
  elseif size(oo_.y_simul,2) < M_.maximum_lag+M_.maximum_lead+options_.periods
    if isempty(ys0_)
      oo_.y_simul = [oo_.y_simul oo_.steady_state*ones(1,M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.y_simul,2),1)];
    else
      oo_.y_simul = [ys0_*ones(1,M_.maximum_lag) oo_.y_simul  oo_.steady_state*ones(1,options_.periods+M_.maximum_lead-size(oo_.y_simul, ...
						  2))];
    end
  end
    
	     