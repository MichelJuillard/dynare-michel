function [dr,info,M_,options_,oo_] = dr1_sparse(dr,task,M_,options_,oo_)
% Computes the reduced form solution of a rational expectation model (first or second order
% approximation of the stochastic model around the deterministic steady state). 
%
% INPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   task       [integer]          if task = 0 then dr1 computes decision rules.
%                                 if task = 1 then dr1 computes eigenvalues.
%   M_         [matlab structure] Definition of the model.           
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results 
%    
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   info       [integer]          info=1: the model doesn't define current variables uniquely
%                                 info=2: problem in mjdgges.dll info(2) contains error code. 
%                                 info=3: BK order condition not satisfied info(2) contains "distance"
%                                         absence of stable trajectory.
%                                 info=4: BK order condition not satisfied info(2) contains "distance"
%                                         indeterminacy.
%                                 info=5: BK rank condition not satisfied.
%   M_         [matlab structure]            
%   options_   [matlab structure]
%   oo_        [matlab structure]
%  
% ALGORITHM
%   ...
%    
% SPECIAL REQUIREMENTS
%   none.
%  

% Copyright (C) 1996-2008 Dynare Team
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

    info = 0;
  
    options_ = set_default_option(options_,'loglinear',0);
    options_ = set_default_option(options_,'noprint',0);
    options_ = set_default_option(options_,'olr',0);
    options_ = set_default_option(options_,'olr_beta',1);
    options_ = set_default_option(options_,'qz_criterium',1.000001);
    
    xlen = M_.maximum_endo_lead + M_.maximum_endo_lag + 1;
    klen = M_.maximum_endo_lag + M_.maximum_endo_lead + 1;
    iyv = M_.lead_lag_incidence';
    iyv = iyv(:);
    iyr0 = find(iyv) ;
    it_ = M_.maximum_lag + 1 ;
    
    if M_.exo_nbr == 0
        oo_.exo_steady_state = [] ;
    end
    
    % expanding system for Optimal Linear Regulator
    if options_.ramsey_policy
        if isfield(M_,'orig_model')
            orig_model = M_.orig_model;
            M_.endo_nbr = orig_model.endo_nbr;
            M_.endo_names = orig_model.endo_names;
            M_.lead_lag_incidence = orig_model.lead_lag_incidence;
            M_.maximum_lead = orig_model.maximum_lead;
            M_.maximum_endo_lead = orig_model.maximum_endo_lead;
            M_.maximum_lag = orig_model.maximum_lag;
            M_.maximum_endo_lag = orig_model.maximum_endo_lag;
        end
        old_solve_algo = options_.solve_algo;
        %  options_.solve_algo = 1;
        oo_.steady_state = dynare_solve('ramsey_static',oo_.steady_state,0,M_,options_,oo_,it_);
        options_.solve_algo = old_solve_algo;
        [junk,junk,multbar] = ramsey_static(oo_.steady_state,M_,options_,oo_,it_);
        [jacobia_,M_] = ramsey_dynamic(oo_.steady_state,multbar,M_,options_,oo_,it_);
        klen = M_.maximum_lag + M_.maximum_lead + 1;
        dr.ys = [oo_.steady_state;zeros(M_.exo_nbr,1);multbar];
% $$$         if options_.ramsey_policy == 2
% $$$             mask = M_.orig_model.lead_lag_incidence ~= 0;
% $$$             incidence_submatrix = M_.lead_lag_incidence(M_.orig_model.maximum_lead+(1:size(mask,1)),1:M_.orig_model.endo_nbr); 
% $$$             k = nonzeros((incidence_submatrix.*mask)');
% $$$             nl = nnz(M_.lead_lag_incidence);
% $$$             k = [k; nl+(1:M_.exo_nbr)'];
% $$$             kk = reshape(1:(nl+M_.exo_nbr)^2,nl+M_.exo_nbr,nl+M_.exo_nbr);
% $$$             kk2 = kk(k,k);
% $$$             
% $$$             k1 = find(M_.orig_model.lead_lag_incidence');
% $$$             y = repmat(oo_.dr.ys(1:M_.orig_model.endo_nbr),1,M_.orig_model.maximum_lag+M_.orig_model.maximum_lead+1);
% $$$             [f,fJ,fh] = feval([M_.fname '_dynamic'],y(k1),zeros(1,M_.exo_nbr), M_.params, it_);
% $$$             
% $$$             % looking for dynamic variables that are both in the original model
% $$$             % and in the optimal policy model
% $$$             k1 = k1+nnz(M_.lead_lag_incidence(1:M_.orig_model.maximum_lead,1:M_.orig_model.endo_nbr));
% $$$             hessian = sparse([],[],[],size(jacobia_,1),(nl+M_.exo_nbr)^2,nnz(fh));
% $$$             hessian(M_.orig_model.endo_nbr+(1:size(fh,1)),kk2) = fh;
% $$$             options_.order = 2;
% $$$         elseif options_.ramsey_policy == 3
% $$$             maxlag1 = M_.orig_model.maximum_lag;
% $$$             maxlead1 = M_.orig_model.maximum_lead;
% $$$             endo_nbr1 = M_.orig_model.endo_nbr;
% $$$             lead_lag_incidence1 = M_.orig_model.lead_lag_incidence;
% $$$             y = repmat(oo_.dr.ys(1:M_.orig_model.endo_nbr),1,M_.orig_model.maximum_lag+M_.orig_model.maximum_lead+1);
% $$$             k1 = find(M_.orig_model.lead_lag_incidence');
% $$$             [f,fj,fh] = feval([M_.fname '_dynamic'],y(k1),zeros(1,M_.exo_nbr), M_.params, it_);
% $$$             nrj = size(fj,1); 
% $$$             
% $$$             iy = M_.lead_lag_incidence;
% $$$             kstate = oo_.dr.kstate;
% $$$             inv_order_var = oo_.dr.inv_order_var;
% $$$             offset = 0;
% $$$             i3 = zeros(0,1);
% $$$             i4 = find(kstate(:,2) <= M_.maximum_lag+1);
% $$$             kstate1 = kstate(i4,:);
% $$$             kk2 = zeros(0,1);
% $$$             % lagged variables
% $$$             for i=2:M_.maximum_lag + 1
% $$$                 i1 = find(kstate(:,2) == i);
% $$$                 k1 = kstate(i1,:);
% $$$                 i2 = find(oo_.dr.order_var(k1(:,1)) <= M_.orig_model.endo_nbr);
% $$$                 i3 = [i3; i2+offset]; 
% $$$                 offset = offset + size(k1,1);
% $$$                 i4 = find(kstate1(:,2) == i);
% $$$                 kk2 = [kk2; i4];
% $$$             end
% $$$             i2 = find(oo_.dr.order_var(k1(:,1)) > M_.orig_model.endo_nbr);
% $$$             j2 = k1(i2,1);
% $$$             nj2 = length(j2);
% $$$             k2 = offset+(1:nj2)';
% $$$             offset = offset + length(i2);
% $$$             i3 = [i3; ...
% $$$                   find(M_.orig_model.lead_lag_incidence(M_.orig_model.maximum_lag+1:end,:)')+offset];
% $$$             i3 = [i3; (1:M_.exo_nbr)'+length(i3)];
% $$$             ni3 = length(i3);
% $$$             nrfj = size(fj,1);
% $$$             jacobia_ = zeros(nrfj+length(j2),ni3);
% $$$             jacobia_(1:nrfj,i3) = fj;
% $$$             jacobia_(nrfj+(1:nj2),1:size(oo_.dr.ghx,2)) = oo_.dr.ghx(j2,:);
% $$$             jacobia_(nrfj+(1:nj2),k2) = eye(nj2);
% $$$             kk1 = reshape(1:ni3^2,ni3,ni3);
% $$$             hessian =  zeros(nrfj+length(j2),ni3^2);
% $$$             hessian(1:nrfj,kk1(i3,i3)) = fh;
% $$$             
% $$$             k = find(any(M_.lead_lag_incidence(1:M_.maximum_lag, ...
% $$$                                                M_.orig_model.endo_nbr+1:end)));
% $$$             if maxlead1 > maxlag1
% $$$                 M_.lead_lag_incidence = [ [zeros(maxlead1-maxlag1,endo_nbr1); ...
% $$$                                     lead_lag_incidence1] ...
% $$$                                     [M_.lead_lag_incidence(M_.maximum_lag+(1:maxlead1), ...
% $$$                                                            k); zeros(maxlead1,length(k))]];
% $$$             elseif maxlag1 > maxlead1
% $$$                 M_.lead_lag_incidence = [ [lead_lag_incidence1; zeros(maxlag1- ...
% $$$                                                                   maxlead1,endo_nbr1);] ...
% $$$                                     [M_.lead_lag_incidence(M_.maximum_lag+(1:maxlead1), ...
% $$$                                                            k); zeros(maxlead1,length(k))]];
% $$$             else % maxlag1 == maxlead1
% $$$                 M_.lead_lag_incidence = [ lead_lag_incidence1
% $$$                                     [M_.lead_lag_incidence(M_.maximum_lag+(1:maxlead1), ...
% $$$                                                            k); zeros(maxlead1,length(k))]];
% $$$             end
% $$$             M_.maximum_lag = max(maxlead1,maxlag1);
% $$$             M_.maximum_endo_lag = M_.maximum_lag;
% $$$             M_.maximum_lead = M_.maximum_lag;
% $$$             M_.maximum_endo_lead = M_.maximum_lag;
% $$$             
% $$$             M_.endo_names = strvcat(M_.orig_model.endo_names, M_.endo_names(endo_nbr1+k,:));
% $$$             M_.endo_nbr = endo_nbr1+length(k);  
% $$$         end
    else
        klen = M_.maximum_lag + M_.maximum_lead + 1;
        iyv = M_.lead_lag_incidence';
        iyv = iyv(:);
        iyr0 = find(iyv) ;
        it_ = M_.maximum_lag + 1 ;
        
        if M_.exo_nbr == 0
            oo_.exo_steady_state = [] ;
        end
        
        it_ = M_.maximum_lag + 1;
        z = repmat(dr.ys,1,klen);
        z = z(iyr0) ;
        if options_.model_mode==0
          if options_.order == 1
              [junk,jacobia_] = feval([M_.fname '_dynamic'],z,[oo_.exo_simul ...
                                  oo_.exo_det_simul], M_.params, it_);
              hessian = 0;
            elseif options_.order == 2
              [junk,jacobia_,hessian] = feval([M_.fname '_dynamic'],z,...
                                              [oo_.exo_simul ...
                                  oo_.exo_det_simul], M_.params, it_);
          end
          dr=set_state_space(dr,M_);
          if options_.debug
            save([M_.fname '_debug.mat'],'jacobia_')
          end
          [dr,info,M_,options_,oo_] = dr11_sparse(dr,task,M_,options_,oo_, jacobia_, hessian);
          dr.nyf = nnz(dr.kstate(:,2)>M_.maximum_lag+1);
        elseif options_.model_mode==1
            if options_.order == 1
                [junk,jacobia_] = feval([M_.fname '_dynamic'],ones(M_.maximum_lag+M_.maximum_lead+1,1)*dr.ys',[oo_.exo_simul ...
                    oo_.exo_det_simul], it_);
                 dr.eigval = [];
                 dr.nyf = 0;
                 dr.rank = 0;
                 for i=1:length(M_.block_structure.block)
                     %disp(['block ' num2str(i)]);
                     M_.block_structure.block(i).dr.Null=0;
                     M_.block_structure.block(i).dr=set_state_space(M_.block_structure.block(i).dr,M_.block_structure.block(i));
                     jcb_=jacobia_(M_.block_structure.block(i).equation,repmat(M_.block_structure.block(i).variable,1,M_.block_structure.block(i).maximum_endo_lag+M_.block_structure.block(i).maximum_endo_lead+1)+kron([M_.maximum_endo_lag-M_.block_structure.block(i).maximum_endo_lag:M_.maximum_endo_lag+M_.block_structure.block(i).maximum_endo_lead],M_.endo_nbr*ones(1,M_.block_structure.block(i).endo_nbr)));
                     jcb_=jcb_(:,find(any(jcb_,1)));
                     hss_=0; %hessian(M_.block_structure.block(i).equation,M_.block_structure.block(i).variable);
                     dra = M_.block_structure.block(i).dr;
                     M_.block_structure.block(i).exo_nbr=M_.exo_nbr;
                     [dra ,info,M_.block_structure.block(i),options_,oo_] = dr11_sparse(dra ,task,M_.block_structure.block(i),options_,oo_, jcb_, hss_);
                     M_.block_structure.block(i).dr = dra;
                     dr.eigval = [dr.eigval; dra.eigval];
                     dr.nyf = dr.nyf + nnz(dra.kstate(:,2)>M_.block_structure.block(i).maximum_endo_lag+1);
                     dr.rank = dr.rank + dra.rank;
                 end;
            end
        end
    end
    
