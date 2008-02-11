% function y_=simult_(y0,dr,ex_,iorder)
% Monte Carlo simulations
% Draws random shocks
%
% INPUTS
%    y0:       starting values
%    dr:       structure of decisions rules for stochastic simulations
%    ex_:      matrix of shocks
%    iorder=0: first-order approximation
%    iorder=1: second-order approximation
%
% OUTPUTS
%    y_:       stochastic simulations results
%
% SPECIAL REQUIREMENTS
%    none
%  
%  
% part of DYNARE, copyright Dynare Team (2001-2007)
% Gnu Public License.

function y_=simult_(y0,dr,ex_,iorder)
global M_ options_ it_
  iter = size(ex_,1);
  if ~isempty(dr.ghu)
      nx = size(dr.ghu,2);
  end
  y_ = zeros(size(y0,1),iter+M_.maximum_lag);
  
  y_(:,1:M_.maximum_lag) = y0;
  k1 = [M_.maximum_lag:-1:1];
  k2 = dr.kstate(find(dr.kstate(:,2) <= M_.maximum_lag+1),[1 2]);
  k2 = k2(:,1)+(M_.maximum_lag+1-k2(:,2))*M_.endo_nbr;
  k3 = M_.lead_lag_incidence(1:M_.maximum_lag,:)';
  k3 = find(k3(:));
  k4 = dr.kstate(find(dr.kstate(:,2) < M_.maximum_lag+1),[1 2]);
  k4 = k4(:,1)+(M_.maximum_lag+1-k4(:,2))*M_.endo_nbr;
  
  options_ = set_default_option(options_,'simul_algo',0);
  if options_.simul_algo == 1
    o1 = dr.nstatic+1;
    o2 = dr.nstatic+dr.npred;
    o3 = o2-dr.nboth+1;
    [junk, k5] = sort(dr.order_var(o1:o2));
    [junk, k6] = sort(dr.order_var(o3:end));
  end

  if iorder == 1    
      if ~isempty(dr.ghu)
          for i = M_.maximum_lag+1: iter+M_.maximum_lag
              tempx1 = y_(dr.order_var,k1);
              tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,M_.maximum_lag);
              tempx = tempx2(k2);
              if options_.simul_algo == 0
                  y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx+dr.ghu* ...
                      ex_(i-M_.maximum_lag,:)';
              elseif options_.simul_algo == 1
                  it_ = i;
                  m = dr.ys(dr.order_var);
                  [y_(:,i), check] = dynare_solve('ff_simul1',y_(:,i-1),tempx1(k3), ...
                                                  m(o3:end),tempx(k4),o1,o2,o3,k6);
              end
              k1 = k1+1;
          end
      else
          for i = M_.maximum_lag+1: iter+M_.maximum_lag
              tempx1 = y_(dr.order_var,k1);
              tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,M_.maximum_lag);
              tempx = tempx2(k2);
              if options_.simul_algo == 0
                  y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx;
              elseif options_.simul_algo == 1
                  it_ = i;
                  m = dr.ys(dr.order_var);
                  [y_(:,i), check] = dynare_solve('ff_simul1',y_(:,i-1),tempx1(k3), ...
                                                  m(o3:end),tempx(k4),o1,o2,o3,k6);
              end
              k1 = k1+1;
          end
      end
  elseif iorder == 2
    for i = M_.maximum_lag+1: iter+M_.maximum_lag
      tempx1 = y_(dr.order_var,k1);
      tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,M_.maximum_lag);
      tempx = tempx2(k2);
      tempu = ex_(i-M_.maximum_lag,:)';
      tempuu = kron(tempu,tempu);
      if options_.simul_algo == 0
	tempxx = kron(tempx,tempx);
	tempxu = kron(tempx,tempu);
	y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+ ...
	    dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu*tempxu;
      elseif options_.simul_algo == 1
	it_ = i;
	m = dr.ys(dr.order_var)+dr.ghs2/2;
	tempx1 = y_(:,k1);
	[y_(:,i), check] = dynare_solve('ff_simul2',y_(:,i-1),tempx1(k3), ...
					m(o3:end),tempx(k4),o1,o2,o3,k6);
      end
      k1 = k1+1;
    end
  end

% MJ 08/30/02 corrected bug at order 2