% Copyright (C) 2001 Michel Juillard
%
function dr=resol(ys,algo,linear,iorder)

global M_ options_ oo_
global  it_  means_ stderrs_


if linear == 1
  iorder =1;
end

xlen = M_.maximum_lead + M_.maximum_lag + 1;
klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv = M_.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = M_.maximum_lag + 1 ;

if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end

if ~ M_.lead_lag_incidence(M_.maximum_lag+1,:) > 0
  error ('RESOL: Error in model specification: some variables don"t appear as current') ;
end

%if xlen > 1
%  error (['RESOL: stochastic exogenous variables must appear only at the' ...
%	  ' current period. Use additional endogenous variables']) ;
%end


% check if ys is steady state
tempex = oo_.exo_simul;
oo_.exo_simul = repmat(oo_.exo_steady_state',options_.periods + M_.maximum_lag + M_.maximum_lead,1);
fh = str2func([M_.fname '_static']);
if max(abs(feval(fh,ys,oo_.exo_simul))) > options_.dynatol
  if exist([M_.fname '_steadystate'])
    [dr.ys,check] = feval([M_.fname '_steadystate'],ys,oo_.exo_simul);
  else
    [dr.ys, check] = dynare_solve([M_.fname '_static'],ys,oo_.exo_simul);
  end
  if check
    error('RESOL: convergence problem in DYNARE_SOLVE')
  end
else
  dr.ys = ys;
end

dr.fbias = zeros(M_.endo_nbr,1);
dr = dr1(iorder,dr,0);

if algo == 1 & iorder > 1
  dr.ys = dynare_solve('dr2',ys,dr);
  dr.fbias = 2*feval([M_.fname '_static'],dr.ys,oo_.exo_simul);
  dr = dr1(iorder,dr,0);
end
oo_.exo_simul = tempex;
tempex = [];

% 01/01/2003 MJ added dr_algo == 1
% 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
%               in dr.ghs2 
% 05/26/2003 MJ added temporary values for oo_.exo_simul



