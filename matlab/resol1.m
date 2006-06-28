% Copyright (C) 2001 Michel Juillard
%
function dr=resol1(ys,algo,linear,iorder)

global M_  options_ oo_ bayestopt_
global  it_  means_ stderrs_ dr1_test_

dr1_test_ = zeros(2,1);

if linear == 1
  iorder =1;
end

xlen = M_.maximum_lead + M_.maximum_lag + 1;
klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv  = M_.lead_lag_incidence';
iyv  = iyv(:);
iyr0 = find(iyv) ;
it_  = M_.maximum_lag + 1 ;

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
oo_.exo_simul = repmat(transpose(oo_.exo_steady_state), M_.maximum_lag + ...
		       M_.maximum_lead+1,1);
fh = str2func([M_.fname '_static']);
if max(abs(feval(fh,ys,oo_.exo_simul))) > options_.dynatol
  % dirty trick to call either user function or dynare_solve
  [dr.ys, cheik] = feval(bayestopt_.static_solve,[M_.fname '_static'],ys,oo_.exo_simul);
  if cheik
    dr1_test_(1) = 1; % dynare_solve did not converge to the steady state.  
    resid = feval([M_.fname '_static'],dr.ys,oo_.exo_simul);
    dr1_test_(2) = resid'*resid; % penalty...
    disp('dynare_solve is unable to find the steady state.')
    return
  end
else
  dr.ys = ys;
end

dr.fbias = zeros(M_.endo_nbr,1);
dr = dr11(iorder,dr,0);

if algo == 1 & iorder > 1
  dr.ys = dynare_solve('dr2',ys,dr);
  dr.fbias = 2*feval([M_.fname '_static'],dr.ys,oo_.exo_simul);
  dr = dr11(iorder,dr,0);
end
oo_.exo_simul = tempex;
tempex = [];

% 01/01/2003 MJ added dr_algo == 1
% 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
%               in dr.ghs2 
% 05/26/2003 MJ added temporary values for oo_.exo_simul
% 01/22/2005 SA check --> cheik (to avoid confusion with the function check.m)
