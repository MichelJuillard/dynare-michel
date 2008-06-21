function [A,B,ys,info] = dynare_resolve(iv,ic,aux)

% function [A,B,ys,info] = dynare_resolve(iv,ic,aux)
% Computes the linear approximation and the matrices A and B of the
% transition equation
%
% INPUTS
%    iv:             selected variables (observed and state variables)
%    ic:             state variables position in the transition matrix columns
%    aux:            indices for auxiliary equations
%
% OUTPUTS
%    A:              matrix of predetermined variables effects in linear solution (ghx)
%    B:              matrix of shocks effects in linear solution (ghu)
%    ys:             steady state of original endogenous variables
%    info=1:         the model doesn't determine the current variables '...' uniquely
%    info=2:         MJDGGES returns the following error code'
%    info=3:         Blanchard Kahn conditions are not satisfied: no stable '...' equilibrium
%    info=4:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy
%    info=5:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy due to rank failure
%    info=11:        same as dr1 for dr_algo = 2
%    info=20:        can't find steady state info(2) contains sum of sqare residuals
%    info=30:        variance can't be computed
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2007)
% Gnu Public License.


  global oo_ M_
  
  [oo_.dr,info] = resol(oo_.steady_state,0);

  if info(1) > 0
    A = [];
    B = [];
    ys = [];
    return
  end
  
  if nargin == 0
    endo_nbr = M_.endo_nbr;
    nstatic = oo_.dr.nstatic;
    npred = oo_.dr.npred;
    iv = (1:endo_nbr)';
    ic = [ nstatic+(1:npred) endo_nbr+(1:size(oo_.dr.ghx,2)-npred) ]';
    aux = oo_.dr.transition_auxiliary_variables;
    k = find(aux(:,2) > npred);
    aux(:,2) = aux(:,2) + nstatic;
    aux(k,2) = aux(k,2) + oo_.dr.nfwrd;
  end
  
  [A,B] = kalman_transition_matrix(oo_.dr,iv,ic,aux,M_.exo_nbr);
  ys = oo_.dr.ys;