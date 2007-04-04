% Copyright (C) 2001 Michel Juillard
%
function [dr,info]=resol(ys,check_flag)
global M_ options_ oo_ bayestopt_
% info: same as dr1 
global it_
% plus: 
% 11 .... same as dr1 for dr_algo = 2
% 20: can't find steady state info(2) contains sum of sqare residuals

 
%unfinished
jacobian_flag = 0; 

options_ = set_default_option(options_,'jacobian_flag',1);
info = 0;

it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end

% check if ys is steady state
tempex = oo_.exo_simul;
oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
if M_.exo_det_nbr > 0 
  tempexdet = oo_.exo_det_simul;
  oo_.exo_det_simul = repmat(oo_.exo_det_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
end
dr.ys = ys;
check1 = 0;
% testing for steadystate file
fh = str2func([M_.fname '_static']);
if options_.steadystate_flag
  [dr.ys,check1] = feval([M_.fname '_steadystate'],dr.ys,...
			 [oo_.exo_steady_state; oo_.exo_det_steady_state]);
else
  % testing if ys isn't a steady state or if we aren't computing Ramsey policy
  if max(abs(feval(fh,dr.ys,[oo_.exo_steady_state; oo_.exo_det_steady_state]))) ...
	> options_.dynatol & options_.ramsey_policy == 0
    if options_.linear == 0
      % nonlinear models
      [dr.ys,check1] = dynare_solve(fh,dr.ys,options_.jacobian_flag,...
				    [oo_.exo_steady_state; oo_.exo_det_steady_state]);
    else
      % linear models
      dr.ys = dr.ys-jacob\fvec;
    end
  end
end
% testing for problem
if check1
  info(1) = 20;
  resid = feval(fh,ys,oo_.exo_steady_state);
  info(2) = resid'*resid; % penalty...
  return
end

dr.fbias = zeros(M_.endo_nbr,1);
[dr,info] = dr1(dr,check_flag);

if info(1)
  return
end

if options_.dr_algo == 1 & options_.order > 1
  dr.ys = dynare_solve('dr2',ys,0,dr);
  dr.fbias = 2*feval([M_.fname '_static'],dr.ys,oo_.exo_steady_state);
  [dr, info1] = dr1(dr,check_flag);
  if info1(1)
    info(1) = info(1)+10;
    return
  end
end
if M_.exo_det_nbr > 0
  oo_.exo_det_simul = tempexdet;
end
oo_.exo_simul = tempex;
tempex = [];

% 01/01/2003 MJ added dr_algo == 1
% 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
%               in dr.ghs2 
% 05/26/2003 MJ added temporary values for oo_.exo_simul
