function y = irf(dr, e1, long, drop, replic, iorder)

% function y = irf(dr, e1, long, drop, replic, iorder)
% Computes impulse response functions
% 
% INPUTS
%    dr:     structure of decisions rules for stochastic simulations
%    e1:     exogenous variables value in time 1 after one shock
%    long
%    drop:   truncation (in order 2)
%    replic: number of replications (in order 2)
%    iorder: first or second order approximation
%
% OUTPUTS
%    y:      impulse response matrix
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  global M_ oo_ options_


  temps = repmat(dr.ys,1,M_.maximum_lag);
  y	= 0;
  
  if iorder == 1
    y1 = repmat(dr.ys,1,long);
    ex2 = zeros(long,M_.exo_nbr);
    ex2(1,:) = e1';
    y2 = simult_(temps,dr,ex2,iorder);
    y = y2(:,M_.maximum_lag+1:end)-y1;
  else
    % eliminate shocks with 0 variance
    i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 ));
    nxs = length(i_exo_var);
    ex1 = zeros(long+drop,M_.exo_nbr);
    ex2 = ex1;
    chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));
    for j = 1: replic
      randn('seed',j);
      ex1(:,i_exo_var) = randn(long+drop,nxs)*chol_S;
      ex2 = ex1;
      ex2(drop+1,:) = ex2(drop+1,:)+e1';   
      y1 = simult_(temps,dr,ex1,iorder);
      y2 = simult_(temps,dr,ex2,iorder);
      y = y+(y2(:,M_.maximum_lag+drop+1:end)-y1(:,M_.maximum_lag+drop+1:end));
    end
    y=y/replic;
  end
