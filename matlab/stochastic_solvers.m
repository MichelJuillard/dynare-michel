function [dr,info,M_,options_,oo_] = stochastic_solvers(dr,task,M_,options_,oo_)
% function [dr,info,M_,options_,oo_] = stochastic_solvers(dr,task,M_,options_,oo_)
% computes the reduced form solution of a rational expectation model (first or second order
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
%                                 info=6: The jacobian matrix evaluated at the steady state is complex.        
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

% Copyright (C) 1996-2009 Dynare Team
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

if (options_.aim_solver == 1) && (options_.order > 1)
        error('Option "aim_solver" is incompatible with order >= 2')
end

if options_.k_order_solver;
    if options_.risky_steadystate
        [dr,info] = dyn_risky_steadystate_solver(oo_.steady_state,M_,dr, ...
                                             options_,oo_);
    else
        dr = set_state_space(dr,M_);
        [dr,info] = k_order_pert(dr,M_,options_,oo_);
    end
    return;
end

if options_.ramsey_policy
    % expanding system for Optimal Linear Regulator
    [jacobia_,dr,info,M_,oo_] = dyn_ramsey_linearized_foc(dr,M_,options_,oo_);
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
    if options_.order == 1
        [junk,jacobia_] = feval([M_.fname '_dynamic'],z,[oo_.exo_simul ...
                            oo_.exo_det_simul], M_.params, dr.ys, it_);
    elseif options_.order == 2
        [junk,jacobia_,hessian1] = feval([M_.fname '_dynamic'],z,...
                                         [oo_.exo_simul ...
                            oo_.exo_det_simul], M_.params, dr.ys, it_);
        if options_.use_dll
            % In USE_DLL mode, the hessian is in the 3-column sparse representation
            hessian1 = sparse(hessian1(:,1), hessian1(:,2), hessian1(:,3), ...
                              size(jacobia_, 1), size(jacobia_, 2)*size(jacobia_, 2));
        end
    end
end

if options_.debug
    save([M_.fname '_debug.mat'],'jacobia_')
end

if ~isreal(jacobia_)
    if max(max(abs(imag(jacobia_)))) < 1e-15
        jacobia_ = real(jacobia_);
    else
        info(1) = 6;
        info(2) = sum(sum(imag(jacobia_).^2));
        return
    end
end

kstate = dr.kstate;
kad = dr.kad;
kae = dr.kae;
nstatic = dr.nstatic;
nfwrd = dr.nfwrd;
npred = dr.npred;
nboth = dr.nboth;
order_var = dr.order_var;
nd = size(kstate,1);
nz = nnz(M_.lead_lag_incidence);

sdyn = M_.endo_nbr - nstatic;

[junk,cols_b,cols_j] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+1, ...
                                                  order_var));
b = zeros(M_.endo_nbr,M_.endo_nbr);
b(:,cols_b) = jacobia_(:,cols_j);

if M_.maximum_endo_lead == 0
    % backward models: simplified code exist only at order == 1
    % If required, use AIM solver if not check only
    if options_.order == 1
        [k1,junk,k2] = find(kstate(:,4));
        temp = -b\jacobia_(:,[k2 nz+1:end]);
        dr.ghx = temp(:,1:npred);
        if M_.exo_nbr
            dr.ghu = temp(:,npred+1:end);
        end
        dr.eigval = eig(transition_matrix(dr));
        dr.rank = 0;
        if any(abs(dr.eigval) > options_.qz_criterium)
            temp = sort(abs(dr.eigval));
            nba = nnz(abs(dr.eigval) > options_.qz_criterium);
            temp = temp(nd-nba+1:nd)-1-options_.qz_criterium;
            info(1) = 3;
            info(2) = temp'*temp;
        end
    else
        error(['2nd and 3rd order approximation not implemented for purely ' ...
               'backward models'])
    end
elseif M_.maximum_endo_lag == 0
    % purely forward model
    dr.ghx = [];
    dr.ghu = -b\jacobia_(:,nz+1:end);
elseif options_.risky_steadystate
    [dr,info] = dyn_risky_steadystate_solver(oo_.steady_state,M_,dr, ...
                                             options_,oo_);
else
    % If required, use AIM solver if not check only
    if (options_.aim_solver == 1) && (task == 0)
        [dr,info] = AIM_first_order_solver(jacobia_,M_,dr,sdim);

    else  % use original Dynare solver
        [dr,info] = dyn_first_order_solver(jacobia_,b,M_,dr,options_,task);
        if info
            return;
        end
    end

    if options_.loglinear == 1
        k = find(dr.kstate(:,2) <= M_.maximum_endo_lag+1);
        klag = dr.kstate(k,[1 2]);
        k1 = dr.order_var;
        
        dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
                 repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
        dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
    end

    %exogenous deterministic variables
    if M_.exo_det_nbr > 0
        f1 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+2:end,order_var))));
        f0 = sparse(jacobia_(:,nonzeros(M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var))));
        fudet = sparse(jacobia_(:,nz+M_.exo_nbr+1:end));
        M1 = inv(f0+[zeros(M_.endo_nbr,nstatic) f1*gx zeros(M_.endo_nbr,nyf-nboth)]);
        M2 = M1*f1;
        dr.ghud = cell(M_.exo_det_length,1);
        dr.ghud{1} = -M1*fudet;
        for i = 2:M_.exo_det_length
            dr.ghud{i} = -M2*dr.ghud{i-1}(end-nyf+1:end,:);
        end
    end

    if options_.order > 1
        % Second order
        dr = dyn_second_order_solver(jacobia_,hessian1,dr,M_,...
                                     options_.threads.kronecker.A_times_B_kronecker_C,...
                                     options_.threads.kronecker.sparse_hessian_times_B_kronecker_C);
    end
end
oo.dr = dr;