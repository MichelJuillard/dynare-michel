function [y_,int_width]=simultxdet(y0,ex,ex_det, iorder,var_list,M_,oo_,options_)
%function [y_,int_width]=simultxdet(y0,ex,ex_det, iorder,var_list,M_,oo_,options_)

% Copyright (C) 2008 Dynare Team
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

%  global endo_nbr ykmin_ xkmin_ it_ options_ iy_ M_ exe_det_ Sigma_e_ lgy_

  dr = oo_.dr;
  ykmin = M_.maximum_lag;
  endo_nbr = M_.endo_nbr;
  exo_det_steady_state = oo_.exo_det_steady_state;
  nstatic = dr.nstatic;
  npred =dr.npred;
  nc = size(dr.ghx,2);
  iter = size(ex,1);
  nx = size(dr.ghu,2);
  y_ = zeros(size(y0,1),iter+ykmin);
  y_(:,1:ykmin) = y0;
  k1 = [ykmin:-1:1];
  k2 = dr.kstate(find(dr.kstate(:,2) <= ykmin+1),[1 2]);
  k2 = k2(:,1)+(ykmin+1-k2(:,2))*endo_nbr;
  k3 = M_.lead_lag_incidence(1:ykmin,:)';
  k3 = find(k3(:));
  k4 = dr.kstate(find(dr.kstate(:,2) < ykmin+1),[1 2]);
  k4 = k4(:,1)+(ykmin+1-k4(:,2))*endo_nbr;
  
  if options_.simul_algo == 1
    o1 = dr.nstatic+1;
    o2 = dr.nstatic+dr.npred;
    o3 = o2-dr.nboth+1;
    [junk, k5] = sort(dr.order_var(o1:o2));
    [junk, k6] = sort(dr.order_var(o3:end));
  end
  
  nvar = size(var_list,1);
  if nvar == 0
    nvar = endo_nbr;
    ivar = [1:nvar];
  else
    ivar=zeros(nvar,1);
    for i=1:nvar
      i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
      if isempty(i_tmp)
	disp(var_list(i,:));
	error (['One of the variable specified does not exist']) ;
      else
	ivar(i) = i_tmp;
      end
    end
  end

  if iorder == 1
    for i = ykmin+1: iter+ykmin
      tempx1 = y_(dr.order_var,k1);
      tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin);
      tempx = tempx2(k2);
      if options_.simul_algo == 0
	y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx+dr.ghu* ...
	    ex(i,:)';
	for j=1:min(ykmin+M_.exo_det_length+1-i,M_.exo_det_length)
	  y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*(ex_det(i+j-1,:)'-exo_det_steady_state');
	end
      elseif options_.simul_algo == 1
	it_ = i;
	m = dr.ys(dr.order_var);
	[y_(:,i), check] = dynare_solve('ff_simul1',y_(:,i-1),tempx1(k3), ...
					m(o3:end),tempx(k4),o1,o2,o3,k6);
      end
	
      k1 = k1+1;
    end
  elseif iorder == 2
    for i = ykmin+1: iter+ykmin
      tempx1 = y_(dr.order_var,k1);
      tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin);
      tempx = tempx2(k2);
      tempu = ex(i-ykmin,:)';
      tempuu = kron(tempu,tempu);
      if options_.simul_algo == 0
	tempxx = kron(tempx,tempx);
	tempxu = kron(tempx,tempu);
	y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+ ...
	    dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu* ...
	    tempxu;
	for j=1:min(ykmin+M_.exo_det_length+1-i,M_.exo_det_length)
	  tempud = ex_det(i+j-1,:)'-exo_det_steady_state;
	  tempudud = kron(tempud,tempud);
	  tempxud = kron(tempx,tempud);
	  tempuud = kron(tempu,tempud);
	  y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*tempud + ...
	      dr.ghxud{j}*tempxud + dr.ghuud{j}*tempuud + ...
	      0.5*dr.ghudud{j,j}*tempudud;
	  for k=1:j-1
	    tempudk = ex_det(i+k-1,:)'-exo_det_steady_state;
	    tempududk = kron(tempudk,tempud);
	    y_(dr.order_var,i) = y_(dr.order_var,i) + ...
		dr.ghudud{k,j}*tempududk;
	  end
	end
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

  [A,B] = kalman_transition_matrix(dr,nstatic+(1:npred),1:nc,dr.transition_auxiliary_variables,M_.exo_nbr);

  inv_order_var = dr.inv_order_var;
  ghx1 = dr.ghx(inv_order_var(ivar),:);
  ghu1 = dr.ghu(inv_order_var(ivar),:);

  sigma_u = B*M_.Sigma_e*B';
  sigma_u1 = ghu1*M_.Sigma_e*ghu1';
  sigma_y = 0;
  
  for i=1:iter
      sigma_y1 = ghx1*sigma_y*ghx1'+sigma_u1;
      var_yf(i,:) = diag(sigma_y1)';
      if i == iter
          break
      end
      sigma_u = A*sigma_u*A';
      sigma_y = sigma_y+sigma_u;
  end

  fact = norminv((1-options_.conf_sig)/2,0,1);

  int_width = zeros(iter,endo_nbr);
  for i=1:endo_nbr
    int_width(:,i) = fact*sqrt(var_yf(:,i));
  end
  
  for i=1:nvar
    my_subplot(i,nvar,2,3,'Forecasts');
    
    plot([-ykmin+1:0],y0(ivar(i),1:ykmin),'b-',...
	 [1:iter],y_(ivar(i),ykmin+1:end),'g-',...
	 [1:iter],y_(ivar(i),ykmin+1:end)'+int_width(:,ivar(i)),'g:',...
	 [1:iter],y_(ivar(i),ykmin+1:end)'-int_width(:,ivar(i)),'g:',...
	 [1 iter],repmat(dr.ys(ivar(i)),1,2),'r-');
    title(M_.endo_names(ivar(i),:));
  end

