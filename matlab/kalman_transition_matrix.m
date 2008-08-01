function [A,B] = kalman_transition_matrix(dr,iv,ic,aux,exo_nbr)
% Builds the transition equation of the state space representation out of ghx and ghu for Kalman filter
% 
% INPUTS
%    dr:      structure of decisions rules for stochastic simulations
%    iv:      selected variables
%    ic:      state variables position in the transition matrix columns
%    aux:     indices for auxiliary equations
%    exo_nbr: number of exogenous variables
%    
% OUTPUTS
%    A:       matrix of predetermined variables effects in linear solution (ghx)
%    B:       matrix of shocks effects in linear solution (ghu)
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2008 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
  
  n_iv = length(iv);
  n_ir1 = size(aux,1);
  nr = n_iv + n_ir1;
  
  A = zeros(nr,nr);
  B = zeros(nr,exo_nbr);
  
  i_n_iv = 1:n_iv;
  A(i_n_iv,ic) = dr.ghx(iv,:);
  if n_ir1 > 0
    A(n_iv+1:end,:) = sparse(aux(:,1),aux(:,2),ones(n_ir1,1),n_ir1,nr);
  end
  
  B(i_n_iv,:) = dr.ghu(iv,:);


% $$$ function [A,B] = kalman_transition_matrix(dr)
% $$$   global M_
% $$$   nx = size(dr.ghx,2)+dr.nfwrd+dr.nstatic;
% $$$   kstate = dr.kstate;
% $$$   ikx = [dr.nstatic+1:dr.nstatic+dr.npred];
% $$$   
% $$$   A = zeros(nx,nx);
% $$$   k0 = kstate(find(kstate(:,2) <= M_.maximum_lag+1),:);
% $$$   i0 = find(k0(:,2) == M_.maximum_lag+1);
% $$$   n0 = size(dr.ghx,1);
% $$$   A(1:n0,dr.nstatic+1:dr.nstatic+dr.npred) = dr.ghx(:,1:dr.npred);
% $$$   A(1:n0,dr.nstatic+dr.npred+dr.nfwrd+1:end) = dr.ghx(:,dr.npred+1:end);
% $$$   B = zeros(nx,M_.exo_nbr);
% $$$   B(1:n0,:) = dr.ghu;
% $$$   offset_col = dr.nstatic;
% $$$   for i=M_.maximum_lag:-1:2
% $$$     i1 = find(k0(:,2) == i);
% $$$     n1 = size(i1,1);
% $$$     j = zeros(n1,1);
% $$$     for j1 = 1:n1
% $$$       j(j1) = find(k0(i0,1)==k0(i1(j1),1));
% $$$     end
% $$$     if i == M_.maximum_lag-1
% $$$       offset_col = dr.nstatic+dr.nfwrd;
% $$$     end
% $$$     A(n0+i1-dr.npred,offset_col+i0(j))=eye(n1);
% $$$     i0 = i1;
% $$$   end
