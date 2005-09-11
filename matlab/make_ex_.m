% Copyright (C) 2001 Michel Juillard
%
function make_ex_
  global M_ options_ oo_ ex0_ ex_det0_
  
  options_ = set_default_option(options_,'periods',0);
  
  if isempty(oo_.exo_steady_state)
    oo_.exo_steady_state = zeros(M_.exo_nbr,1);
  end
  if M_.exo_det_nbr > 1 & isempty(oo_.exo_det_steadystate)
    oo_.exo_det_steadystate = zeros(M_.exo_det_nbr,1);
  end
  if isempty(oo_.exo_simul)
    if isempty(ex0_)
      oo_.exo_simul = [ones(M_.maximum_lag+options_.periods+M_.maximum_lead,1)*oo_.exo_steady_state'];
    else
      oo_.exo_simul = [ones(M_.maximum_lag,1)*ex0_';ones(options_.periods+M_.maximum_lead,1)*oo_.exo_steady_state'];
    end
  elseif size(oo_.exo_simul,2) < length(oo_.exo_steady_state)
    k = size(oo_.exo_simul,2)+1:length(oo_.exo_steady_state);
    if isempty(ex0_)
      oo_.exo_simul = [oo_.exo_simul ones(M_.maximum_lag+size(oo_.exo_simul,1)+M_.maximum_lead,1)*oo_.exo_steady_state(k)'];
    else
      oo_.exo_simul = [oo_.exo_simul [ones(M_.maximum_lag,1)*ex0_(k)'; ones(size(oo_.exo_simul,1)-M_.maximum_lag+M_.maximum_lead, ...
						1)*oo_.exo_steady_state(k)']];
    end
  elseif size(oo_.exo_simul,1) < M_.maximum_lag+M_.maximum_lead+options_.periods
    if isempty(ex0_)
      oo_.exo_simul = [oo_.exo_simul; ones(M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.exo_simul,1),1)*oo_.exo_steady_state'];
    else
      oo_.exo_simul = [ones(M_.maximum_lag,1)*ex0_'; oo_.exo_simul; ones(options_.periods+M_.maximum_lead-size(oo_.exo_simul, ...
						  1),1)*oo_.exo_steady_state'];
    end
  end
  if M_.exo_det_nbr > 0
    if isempty(oo_.exo_det_simul)
      if isempty(ex_det0_)
	oo_.exo_det_simul = [ones(ykmin_+options_.periods+ykmax_,1)*M_.exo_det_steadystate'];
      else
	oo_.exo_det_simul = [ones(ykmin_,1)*ex_det0_';ones(options_.periods+ykmax_,1)*M_.exo_det_steadystate'];
      end
    elseif size(oo_.exo_det_simul,2) < length(M_.exo_det_steadystate)
      k = size(oo_.exo_det_simul,2)+1:length(M_.exo_det_steadystate);
      if isempty(ex_det0_)
	oo_.exo_det_simul = [oo_.exo_det_simul ones(ykmin_+size(oo_.exo_det_simul,1)+ykmax_,1)*M_.exo_det_steadystate(k)'];
      else
	oo_.exo_det_simul = [oo_.exo_det_simul [ones(ykmin_,1)*ex_det0_(k)'; ones(size(oo_.exo_det_simul,1)-ykmin_+ykmax_, ...
						  1)*M_.exo_det_steadystate(k)']];
      end
    elseif size(oo_.exo_det_simul,1) < ykmin_+ykmax_+options_.periods
      if isempty(ex_det0_)
	oo_.exo_det_simul = [oo_.exo_det_simul; ones(ykmin_+options_.periods+ykmax_-size(oo_.exo_det_simul,1),1)*M_.exo_det_steadystate'];
      else
	oo_.exo_det_simul = [ones(ykmin_,1)*ex_det0_'; oo_.exo_det_simul; ones(options_.periods+ykmax_-size(oo_.exo_det_simul, ...
						  1),1)*M_.exo_det_steadystate'];
      end
    end
  end
    
	     