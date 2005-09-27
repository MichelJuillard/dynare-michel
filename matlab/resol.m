% Copyright (C) 2001 Michel Juillard
%
function [dr,info]=resol(ys,check_flag)
global M_ options_ oo_
% info: same as dr1 
global it_
% plus: 
% 11 .... same as dr1 for dr_algo = 2
% 20: can't find steady state info(2) contains sum of sqare residuals

  

options_ = set_default_option(options_,'olr',0);
info = 0;

it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end

if M_.maximum_exo_lag > 0 || M_.maximum_exo_lead > 0
  error (['RESOL: stochastic exogenous variables must appear only at the' ...
	  ' current period. Use additional endogenous variables']) ;
end


% check if ys is steady state
tempex = oo_.exo_simul;
oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
if M_.exo_det_nbr > 0 
  tempexdet = oo_.exo_det_simul;
  oo_.exo_det_simul = ones(M_.maximum_lag+1,1)*oo_.exo_steady_statedet_';
end
dr.ys = ys;
fh = str2func([M_.fname '_static']);
if options_.linear == 0
  if max(abs(feval(fh,dr.ys,oo_.exo_steady_state))) > options_.dynatol & options_.olr == 0
    if exist([M_.fname '_steadystate'])
      [dr.ys,check1] = feval([M_.fname '_steadystate'],dr.ys,oo_.exo_steady_state);
    else
      [dr.ys,check1] = dynare_solve(fh,dr.ys,oo_.exo_steady_state);
    end
    if check1
      info(1) = 20;
      resid = feval(fh,ys,oo_.exo_steady_state);
      info(2) = resid'*resid; % penalty...
      return
    end
  end
else
  [fvec,jacob] = feval(fh,dr.ys,oo_.exo_steady_state);
  if max(abs(fvec)) > options_.dynatol & options_.olr == 0
    dr.ys = dr.ys-jacob\fvec;
  end
end

dr.fbias = zeros(M_.endo_nbr,1);
[dr,info] = dr1(dr,check_flag);

if info(1)
  return
end

if options_.dr_algo == 1 & options_.order > 1
  dr.ys = dynare_solve('dr2',ys,dr);
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