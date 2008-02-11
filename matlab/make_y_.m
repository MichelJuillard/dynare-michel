
function make_y_
% function make_y_
% forms oo_.endo_simul as guess values for deterministic simulations
%  
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
  global M_ options_ oo_ ys0_ 
  
  options_ = set_default_option(options_,'periods',0);
  
  if isempty(oo_.steady_state)
    oo_.steady_state = zeros(M_.endo_nbr,1);
  end
  
  
  if isempty(oo_.endo_simul)
    if isempty(ys0_)
      oo_.endo_simul = [oo_.steady_state*ones(1,M_.maximum_lag+options_.periods+M_.maximum_lead)];
    else
      oo_.endo_simul = [ys0_*ones(1,M_.maximum_lag) oo_.steady_state*ones(1,options_.periods+M_.maximum_lead)];
    end
  elseif size(oo_.endo_simul,2) < M_.maximum_lag+M_.maximum_lead+options_.periods
      if size(oo_.exo_simul)
          oo_.endo_simul = [oo_.endo_simul ...
                        oo_.steady_state*ones(1,M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.endo_simul,2),1)];
      else
          %% A linear approximation is used to initiate the solution.
          oldopt = options_;
          options_.order = 1;
          dr = oo_.dr;
          dr.ys = oo_.steady_state;
          [dr,info]=dr1(dr,0);
          exogenous_variables = zeros(M_.maximum_lag+options_.periods+M_.maximum_lead-size(oo_.endo_simul,2)+1,0);
          y0 = oo_.endo_simul(:,1:M_.maximum_lag);
          oo_.endo_simul=simult_(y0,dr,exogenous_variables,1);
          options_ = oldopt;
      end
  end