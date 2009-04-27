function [resids, rJ,mult] = dyn_ramsey_static_(x,M_,options_,oo_,it_)

% function [resids, rJ,mult] = dyn_ramsey_static_(x)
% Computes the static first order conditions for optimal policy
%
% INPUTS
%    x:         vector of endogenous variables
%
% OUTPUTS
%    resids:    residuals of non linear equations
%    rJ:        Jacobian
%    mult:      Lagrangian multipliers
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2007 Dynare Team
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

  
  % recovering usefull fields
  endo_nbr = M_.endo_nbr;
  exo_nbr = M_.exo_nbr;
  fname = M_.fname;
  % inst_nbr = M_.inst_nbr;
  % i_endo_no_inst = M_.endogenous_variables_without_instruments;
  max_lead = M_.maximum_lead;
  max_lag = M_.maximum_lag;
  beta =  options_.planner_discount;
  
  % indices of all endogenous variables
  i_endo = [1:endo_nbr]';
  % indices of endogenous variable except instruments
  % i_inst = M_.instruments;
  % lead_lag incidence matrix for endogenous variables
  i_lag = M_.lead_lag_incidence;
  
  % value and Jacobian of objective function
  ex = zeros(1,M_.exo_nbr);
  [U,Uy,Uyy] = feval([fname '_objective_static'],x(i_endo),ex, M_.params);
  Uy = Uy';
  Uyy = reshape(Uyy,endo_nbr,endo_nbr);
  
  % value and Jacobian of dynamic function
  y = repmat(x(i_endo),1,max_lag+max_lead+1);
  k = find(i_lag');
  it_ = 1;
%  [f,fJ,fH] = feval([fname '_dynamic'],y(k),ex);
  [f,fJ] = feval([fname '_dynamic'],y(k),[oo_.exo_simul oo_.exo_det_simul], M_.params, it_);
  % indices of Lagrange multipliers
  inst_nbr = endo_nbr - size(f,1);
  i_mult = [endo_nbr+1:2*endo_nbr-inst_nbr]';
  
  % derivatives of Lagrangian with respect to endogenous variables
%  res1 = Uy;
  A = zeros(endo_nbr,endo_nbr-inst_nbr);
  for i=1:max_lag+max_lead+1
    % select variables present in the model at a given lag
    [junk,k1,k2] = find(i_lag(i,:));
%    res1(k1) = res1(k1) + beta^(max_lag-i+1)*fJ(:,k2)'*x(i_mult); 
    A(k1,:) = A(k1,:) + beta^(max_lag-i+1)*fJ(:,k2)';
  end
  
%  i_inst = var_index(options_.olr_inst);
%  k = setdiff(1:size(A,1),i_inst);
%  mult = -A(k,:)\Uy(k);
  mult = -A\Uy;
%  resids = [f; Uy(i_inst)+A(i_inst,:)*mult];
  resids1 = Uy+A*mult;
%  resids = [f; sqrt(resids1'*resids1/endo_nbr)]; 
  [q,r,e] = qr([A Uy]');
  resids = [f; r(end,(endo_nbr-inst_nbr+1:end))'];
  rJ = [];
  return;
  
  % Jacobian of first order conditions
  n = nnz(i_lag)+exo_nbr;
  iH = reshape(1:n^2,n,n);
  rJ = zeros(2*endo_nbr-inst_nbr,2*endo_nbr-inst_nbr);
  
  rJ(i_endo,i_endo) = Uyy;
  for i=1:max_lag+max_lead+1
    % select variables present in the model at a given lag
    [junk,k1,k2] = find(i_lag(i,:));
    k3 = length(k2);
    rJ(k1,k1) = rJ(k1,k1) + beta^(max_lag-i+1)*reshape(fH(:,iH(k2,k2))'*x(i_mult),k3,k3); 
    rJ(k1,i_mult) = rJ(k1,i_mult) + beta^(max_lag-1+1)*fJ(:,k2)';
    rJ(i_mult,k1) = rJ(i_mult,k1) + fJ(:,k2);
  end
  
%  rJ = 1e-3*rJ;
%  rJ(209,210) = rJ(209,210)+1-1e-3;


  