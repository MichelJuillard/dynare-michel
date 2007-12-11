% Copyright (C) 2005 Michel Juillard
%
function forecast(var_list)
  global options_ dr_ oo_ M_
  
  old_options = options_;
  options_ = set_default_option(options_,'periods',40);
  if options_.periods == 0
    options_.periods = 40;
  end
  options_ = set_default_option(options_,'conf_sig',0.9);
  
  if size(oo_.endo_simul,2) < M_.maximum_lag
    y0 = repmat(oo_.steady_state,M_.maximum_lag);
  else
    y0 = oo_.endo_simul(:,1:M_.maximum_lag);
  end
  
  if M_.exo_det_nbr == 0
    [yf,int_width] = forcst(oo_.dr,y0,options_.periods,var_list);
  else
    exo_det_length = size(oo_.exo_det_simul,1);
    if options_.periods > exo_det_length
      ex = zeros(options_.periods,M_.exo_nbr);
      oo_.exo_det_simul = [ oo_.exo_det_simul;...
		    repmat(oo_.exo_det_steady_state',...
			   options_.periods- ... 
			   exo_det_length,1)];
			   %ex_det_length,1),1)];
    elseif options_.periods < exo_det_length 
      ex = zeros(exo_det_length,M_.exo_nbr); 
    end
    [yf,int_width] = simultxdet(y0,dr_,ex,oo_.exo_det_simul,...
				options_.order,var_list);
  end
  
  for i=1:M_.endo_nbr
    eval(['oo_.forecast.Mean.' M_.endo_names(i,:) '= yf(' int2str(i) ',M_.maximum_lag+1:end)'';']);
    eval(['oo_.forecast.HPDinf.' M_.endo_names(i,:) '= yf(' int2str(i) ',M_.maximum_lag+1:end)''-' ...
		    ' int_width(:,' int2str(i) ');']);
    eval(['oo_.forecast.HPDsup.' M_.endo_names(i,:) '= yf(' int2str(i) ',M_.maximum_lag+1:end)''+' ...
		    ' int_width(:,' int2str(i) ');']);
  end

  for i=1:M_.exo_det_nbr
    eval(['oo_.forecast.Exogenous.' lgx_det_(i,:) '= M_.exo_det_simul(:,' int2str(i) ');']);
  end
  
  options_ = old_options;