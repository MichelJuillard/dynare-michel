function [A,B] = ghx2transition(mm,iv,ic,aux)
% [A,B] = ghx2transition(mm,iv,ic,aux)
%
% Adapted by M. Ratto from kalman_transition_matrix.m 
% (kalman_transition_matrix.m is part of DYNARE, copyright M. Juillard)
%
% Part of the Sensitivity Analysis Toolbox for DYNARE
%
% Written by Marco Ratto, 2006
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it 
%
% Disclaimer: This software is not subject to copyright protection and is in the public domain. 
% It is an experimental system. The Joint Research Centre of European Commission 
% assumes no responsibility whatsoever for its use by other parties
% and makes no guarantees, expressed or implied, about its quality, reliability, or any other
% characteristic. We would appreciate acknowledgement if the software is used.
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.
%

global oo_ M_

  [nr1, nc1] = size(mm);
  ghx = mm(:, [1:(nc1-M_.exo_nbr)]);
  ghu = mm(:, [(nc1-M_.exo_nbr+1):end] );
  if nargin == 1
    oo_.dr.ghx = ghx;
    oo_.dr.ghu = ghu;
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
  n_iv = length(iv);
  n_ir1 = size(aux,1);
  nr = n_iv + n_ir1;
  
  A = zeros(nr,nr);
  B = zeros(nr,M_.exo_nbr);
  
  i_n_iv = 1:n_iv;
  A(i_n_iv,ic) = ghx(iv,:);
  if n_ir1 > 0
    A(n_iv+1:end,:) = sparse(aux(:,1),aux(:,2),ones(n_ir1,1),n_ir1,nr);
  end
  
  B(i_n_iv,:) = ghu(iv,:);
