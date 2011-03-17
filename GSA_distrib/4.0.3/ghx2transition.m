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

global M_

  [nr1, nc1] = size(mm);
  ghx = mm(:, [1:(nc1-M_.exo_nbr)]);
  ghu = mm(:, [(nc1-M_.exo_nbr+1):end] );
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
