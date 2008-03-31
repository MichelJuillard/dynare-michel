function steady()

% function steady()
% computes and prints the steady state calculations
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


  global M_ oo_ options_ ys0_ 

  options_ = set_default_option(options_,'jacobian_flag',1);
  options_ = set_default_option(options_,'steadystate_flag',0);
  if exist([M_.fname '_steadystate.m'])
    options_.steadystate_flag = 1;
  end 

  if options_.steadystate_flag && options_.homotopy_mode
    error('STEADY: Can''t use homotopy when providing a steady state external file');
  end
  
  switch options_.homotopy_mode
    case 1
     homotopy1(options_.homotopy_values, options_.homotopy_steps);
    case 2
     homotopy2(options_.homotopy_values, options_.homotopy_steps);
    case 3
     homotopy3(options_.homotopy_values, options_.homotopy_steps);
  end
  
  steady_;
  
  disp(' ')
  disp('STEADY-STATE RESULTS:')
  disp(' ')
  endo_names = M_.endo_names;
  steady_state = oo_.steady_state;
  for i=1:size(oo_.steady_state,1)
    disp(sprintf('%s \t\t %g',endo_names(i,:),steady_state(i)));
  end
  
  %%% Unless I'm wrong, this is (should be?) done in make_y_.m
% $$$   if isempty(ys0_)
% $$$     oo_.endo_simul(:,1:M_.maximum_lag) = oo_.steady_state * ones(1, M_.maximum_lag);
% $$$   else
% $$$     options_ =set_default_option(options_,'periods',1);
% $$$     oo_.endo_simul(:,M_.maximum_lag+1:M_.maximum_lag+options_.periods+ ...
% $$$                    M_.maximum_lead) = oo_.steady_state * ones(1,options_.periods+M_.maximum_lead);
% $$$   end