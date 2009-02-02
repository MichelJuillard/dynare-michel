function [dr,info,M_,options_,oo_] = dr1(dr,task,M_,options_,oo_)
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

    end
    
    if options_.debug
        save([M_.fname '_debug.mat'],'jacobia_')
    end
    
    dr=set_state_space(dr,M_);
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
    
    k0 = M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var);
    k1 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_endo_lag+1),:);
   
        
        
    if options_.order == 1
        M_.var_order_endo_names=M_.endo_names(dr.order_var,:);
%            z = repmat(dr.ys,1,klen);
%            z = z(iyr0) ;
%            oo_.dyn_ys=z;  % extended ys 
        try
            [ysteady, gx, gu]=k_order_perturbation(dr,task,M_,options_, oo_ );
            load(M_.fname);
            ghxu = eval([M_.fname '_g_1']);
            sss= size(ghxu,2);
            dr.ghx= ghxu(:,1:sss-M_.exo_nbr); 
            dr.ghu= ghxu(:,sss-M_.exo_nbr+1:end); 
            dr.ys=eval([M_.fname '_ss']);
        catch
            disp('*************************************************************************************');
%           disp('Problem with using k_order perturbation solver - Using Dynare solver instead');
%            warning('Problem with using k_order perturbation solver - Using Dynare solver instead');
            error('Problem with using k_order perturbation solver ');
            disp('*****************************************************************************');
            options_.use_k_order=0; % and then try mjdgges instead
            info(1) = 4;
            info(2) = 1000;
            return
        end

    elseif options_.order > 1
        error(' can not use order > 1 with K-Order yet!')
        % or ???
        disp('********************************************************************');
        disp(' can not use order > 1 with K-Order yet - Using Dynare solver instead');
        disp('********************************************************************');
        options_.use_k_order= 0; % and then try mjdgges instead
        info(1) = 4;
        info(2) = 1000;
        return
    end

    
    if M_.maximum_endo_lead == 0;  % backward models
        % If required, try Gary Anderson and G Moore AIM solver if not
        % check only and if 1st order (added by GP July'08)

        dr.eigval = eig(transition_matrix(dr));
        dr.rank = 0;
        if any(abs(dr.eigval) > options_.qz_criterium)
            temp = sort(abs(dr.eigval));
            nba = nnz(abs(dr.eigval) > options_.qz_criterium);
            temp = temp(nd-nba+1:nd)-1-options_.qz_criterium;
            info(1) = 3;
            info(2) = temp'*temp;
        end
        return;
    end
    
    %forward--looking models
        [A,B] =transition_matrix(dr);
        dr.eigval = eig(A);
%            if any(abs(dr.eigval) > options_.qz_criterium)
%                temp = sort(abs(dr.eigval));
%                nba = nnz(abs(dr.eigval) > options_.qz_criterium);
%                temp = temp(nd-nba+1:nd)-1-options_.qz_criterium;
%                info(1) = 3;
%                info(2) = temp'*temp;
%                return
%            end
        sdim = sum( abs(dr.eigval) < options_.qz_criterium );
        nba = nd-sdim;

        nyf = sum(kstate(:,2) > M_.maximum_endo_lag+1);
        if nba ~= nyf
            temp = sort(abs(dr.eigval));
            if nba > nyf
                temp = temp(nd-nba+1:nd-nyf)-1-options_.qz_criterium;
                info(1) = 3;
            elseif nba < nyf;
                temp = temp(nd-nyf+1:nd-nba)-1-options_.qz_criterium;
                info(1) = 4;
            end
            info(2) = temp'*temp;
            return
        end

    
    if options_.loglinear == 1
        k = find(dr.kstate(:,2) <= M_.maximum_endo_lag+1);
        klag = dr.kstate(k,[1 2]);
        k1 = dr.order_var;
        
        dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
                 repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
        dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
    end
    
    dr.ghx = real(dr.ghx);
    dr.ghu = real(dr.ghu);
        
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    if options_.order == 1
        return
    end
    
    % Second order
    %tempex = oo_.exo_simul ;
        [junk,jacobia_,hessian] = feval([M_.fname '_dynamic'],z,...
                [oo_.exo_simul ...
                oo_.exo_det_simul], M_.params, it_);
    
    %hessian = real(hessext('ff1_',[z; oo_.exo_steady_state]))' ;
    kk = flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_endo_lag+1:end,order_var)),1));
    if M_.maximum_endo_lag > 0
        kk = [cumsum(M_.lead_lag_incidence(1:M_.maximum_endo_lag,order_var),1); kk];
    end
    kk = kk';
    kk = find(kk(:));
    nk = size(kk,1) + M_.exo_nbr + M_.exo_det_nbr;
    k1 = M_.lead_lag_incidence(:,order_var);
    k1 = k1';
    k1 = k1(:);
    k1 = k1(kk);
    k2 = find(k1);
    kk1(k1(k2)) = k2;
    kk1 = [kk1 length(k1)+1:length(k1)+M_.exo_nbr+M_.exo_det_nbr];
    kk = reshape([1:nk^2],nk,nk);
    kk1 = kk(kk1,kk1);
    %[junk,junk,hessian] = feval([M_.fname '_dynamic'],z, oo_.exo_steady_state);
    hessian(:,kk1(:)) = hessian;
    
    %oo_.exo_simul = tempex ;
    %clear tempex
    
    n1 = 0;
    n2 = np;
    zx = zeros(np,np);
    zu=zeros(np,M_.exo_nbr);
    for i=2:M_.maximum_endo_lag+1
        k1 = sum(kstate(:,2) == i);
        zx(n1+1:n1+k1,n2-k1+1:n2)=eye(k1);
        n1 = n1+k1;
        n2 = n2-k1;
    end
    kk = flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_endo_lag+1:end,order_var)),1));
    k0 = [1:M_.endo_nbr];
    gx1 = dr.ghx;
    hu = dr.ghu(nstatic+[1:npred],:);
    zx = [zx; gx1];
    zu = [zu; dr.ghu];
    for i=1:M_.maximum_endo_lead
        k1 = find(kk(i+1,k0) > 0);
        zu = [zu; gx1(k1,1:npred)*hu];
        gx1 = gx1(k1,:)*hx;
        zx = [zx; gx1];
        kk = kk(:,k0);
        k0 = k1;
    end
    zx=[zx; zeros(M_.exo_nbr,np);zeros(M_.exo_det_nbr,np)];
    zu=[zu; eye(M_.exo_nbr);zeros(M_.exo_det_nbr,M_.exo_nbr)];
    [nrzx,nczx] = size(zx);
    
    rhs = -sparse_hessian_times_B_kronecker_C(hessian,zx);
    
    %lhs
    n = M_.endo_nbr+sum(kstate(:,2) > M_.maximum_endo_lag+1 & kstate(:,2) < M_.maximum_endo_lag+M_.maximum_endo_lead+1);
    A = zeros(n,n);
    B = zeros(n,n);
    A(1:M_.endo_nbr,1:M_.endo_nbr) = jacobia_(:,M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var));
    % variables with the highest lead
    k1 = find(kstate(:,2) == M_.maximum_endo_lag+M_.maximum_endo_lead+1);
    if M_.maximum_endo_lead > 1
        k2 = find(kstate(:,2) == M_.maximum_endo_lag+M_.maximum_endo_lead);
        [junk,junk,k3] = intersect(kstate(k1,1),kstate(k2,1));
    else
        k2 = [1:M_.endo_nbr];
        k3 = kstate(k1,1);
    end
    % Jacobian with respect to the variables with the highest lead
    B(1:M_.endo_nbr,end-length(k2)+k3) = jacobia_(:,kstate(k1,3)+M_.endo_nbr);
    offset = M_.endo_nbr;
    k0 = [1:M_.endo_nbr];
    gx1 = dr.ghx;
    for i=1:M_.maximum_endo_lead-1
        k1 = find(kstate(:,2) == M_.maximum_endo_lag+i+1);
        [k2,junk,k3] = find(kstate(k1,3));
        A(1:M_.endo_nbr,offset+k2) = jacobia_(:,k3+M_.endo_nbr);
        n1 = length(k1);
        A(offset+[1:n1],nstatic+[1:npred]) = -gx1(kstate(k1,1),1:npred);
        gx1 = gx1*hx;
        A(offset+[1:n1],offset+[1:n1]) = eye(n1);
        n0 = length(k0);
        E = eye(n0);
        if i == 1
            [junk,junk,k4]=intersect(kstate(k1,1),[1:M_.endo_nbr]);
        else
            [junk,junk,k4]=intersect(kstate(k1,1),kstate(k0,1));
        end
        i1 = offset-n0+n1;
        B(offset+[1:n1],offset-n0+[1:n0]) = -E(k4,:);
        k0 = k1;
        offset = offset + n1;
    end
    [junk,k1,k2] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+M_.maximum_endo_lead+1,order_var));
    A(1:M_.endo_nbr,nstatic+1:nstatic+npred)=...
        A(1:M_.endo_nbr,nstatic+[1:npred])+jacobia_(:,k2)*gx1(k1,1:npred);
    C = hx;
    D = [rhs; zeros(n-M_.endo_nbr,size(rhs,2))];
    
    
    dr.ghxx = gensylv(2,A,B,C,D);
    
    %ghxu
    %rhs
    hu = dr.ghu(nstatic+1:nstatic+npred,:);
    %kk = reshape([1:np*np],np,np);
    %kk = kk(1:npred,1:npred);
    %rhs = -hessian*kron(zx,zu)-f1*dr.ghxx(end-nyf+1:end,kk(:))*kron(hx(1:npred,:),hu(1:npred,:));
    
    rhs = sparse_hessian_times_B_kronecker_C(hessian,zx,zu);
    
    nyf1 = sum(kstate(:,2) == M_.maximum_endo_lag+2);
    hu1 = [hu;zeros(np-npred,M_.exo_nbr)];
    %B1 = [B(1:M_.endo_nbr,:);zeros(size(A,1)-M_.endo_nbr,size(B,2))];
    [nrhx,nchx] = size(hx);
    [nrhu1,nchu1] = size(hu1);
    
    B1 = B*A_times_B_kronecker_C(dr.ghxx,hx,hu1);
    rhs = -[rhs; zeros(n-M_.endo_nbr,size(rhs,2))]-B1;
    
    
    %lhs
    dr.ghxu = A\rhs;
    
    %ghuu
    %rhs
    kk = reshape([1:np*np],np,np);
    kk = kk(1:npred,1:npred);
    
    rhs = sparse_hessian_times_B_kronecker_C(hessian,zu);
    
    
    B1 = A_times_B_kronecker_C(B*dr.ghxx,hu1);
    rhs = -[rhs; zeros(n-M_.endo_nbr,size(rhs,2))]-B1;
    
    %lhs
    dr.ghuu = A\rhs;
    
    dr.ghxx = dr.ghxx(1:M_.endo_nbr,:);
    dr.ghxu = dr.ghxu(1:M_.endo_nbr,:);
    dr.ghuu = dr.ghuu(1:M_.endo_nbr,:);
    
    
    % dr.ghs2
    % derivatives of F with respect to forward variables
    % reordering predetermined variables in diminishing lag order
    O1 = zeros(M_.endo_nbr,nstatic);
    O2 = zeros(M_.endo_nbr,M_.endo_nbr-nstatic-npred);
    LHS = jacobia_(:,M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var));
    RHS = zeros(M_.endo_nbr,M_.exo_nbr^2);
    kk = find(kstate(:,2) == M_.maximum_endo_lag+2);
    gu = dr.ghu; 
    guu = dr.ghuu; 
    Gu = [dr.ghu(nstatic+[1:npred],:); zeros(np-npred,M_.exo_nbr)];
    Guu = [dr.ghuu(nstatic+[1:npred],:); zeros(np-npred,M_.exo_nbr*M_.exo_nbr)];
    E = eye(M_.endo_nbr);
    M_.lead_lag_incidenceordered = flipud(cumsum(flipud(M_.lead_lag_incidence(M_.maximum_endo_lag+1:end,order_var)),1));
    if M_.maximum_endo_lag > 0
        M_.lead_lag_incidenceordered = [cumsum(M_.lead_lag_incidence(1:M_.maximum_endo_lag,order_var),1); M_.lead_lag_incidenceordered];
    end
    M_.lead_lag_incidenceordered = M_.lead_lag_incidenceordered';
    M_.lead_lag_incidenceordered = M_.lead_lag_incidenceordered(:);
    k1 = find(M_.lead_lag_incidenceordered);
    M_.lead_lag_incidenceordered(k1) = [1:length(k1)]';
    M_.lead_lag_incidenceordered =reshape(M_.lead_lag_incidenceordered,M_.endo_nbr,M_.maximum_endo_lag+M_.maximum_endo_lead+1)';
    kh = reshape([1:nk^2],nk,nk);
    kp = sum(kstate(:,2) <= M_.maximum_endo_lag+1);
    E1 = [eye(npred); zeros(kp-npred,npred)];
    H = E1;
    hxx = dr.ghxx(nstatic+[1:npred],:);
    for i=1:M_.maximum_endo_lead
        for j=i:M_.maximum_endo_lead
            [junk,k2a,k2] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+j+1,order_var));
            [junk,k3a,k3] = ...
                find(M_.lead_lag_incidenceordered(M_.maximum_endo_lag+j+1,:));
            nk3a = length(k3a);
            B1 = sparse_hessian_times_B_kronecker_C(hessian(:,kh(k3,k3)),gu(k3a,:));
            RHS = RHS + jacobia_(:,k2)*guu(k2a,:)+B1;
        end
        % LHS
        [junk,k2a,k2] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+i+1,order_var));
        LHS = LHS + jacobia_(:,k2)*(E(k2a,:)+[O1(k2a,:) dr.ghx(k2a,:)*H O2(k2a,:)]);
        
        if i == M_.maximum_endo_lead 
            break
        end
        
        kk = find(kstate(:,2) == M_.maximum_endo_lag+i+1);
        gu = dr.ghx*Gu;
        [nrGu,ncGu] = size(Gu);
        G1 = A_times_B_kronecker_C(dr.ghxx,Gu);
        G2 = A_times_B_kronecker_C(hxx,Gu);
        guu = dr.ghx*Guu+G1;
        Gu = hx*Gu;
        Guu = hx*Guu;
        Guu(end-npred+1:end,:) = Guu(end-npred+1:end,:) + G2;
        H = E1 + hx*H;
    end
    RHS = RHS*M_.Sigma_e(:);
    dr.fuu = RHS;
    %RHS = -RHS-dr.fbias;
    RHS = -RHS;
    dr.ghs2 = LHS\RHS;
    
    % deterministic exogenous variables
    if M_.exo_det_nbr > 0
        hud = dr.ghud{1}(nstatic+1:nstatic+npred,:);
        zud=[zeros(np,M_.exo_det_nbr);dr.ghud{1};gx(:,1:npred)*hud;zeros(M_.exo_nbr,M_.exo_det_nbr);eye(M_.exo_det_nbr)];
        R1 = hessian*kron(zx,zud);
        dr.ghxud = cell(M_.exo_det_length,1);
        kf = [M_.endo_nbr-nyf+1:M_.endo_nbr];
        kp = nstatic+[1:npred];
        dr.ghxud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{1}(kp,:)));
        Eud = eye(M_.exo_det_nbr);
        for i = 2:M_.exo_det_length
            hudi = dr.ghud{i}(kp,:);
            zudi=[zeros(np,M_.exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
            R2 = hessian*kron(zx,zudi);
            dr.ghxud{i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(hx,Eud)+dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{i}(kp,:)))-M1*R2;
        end
        R1 = hessian*kron(zu,zud);
        dr.ghudud = cell(M_.exo_det_length,1);
        kf = [M_.endo_nbr-nyf+1:M_.endo_nbr];
        
        dr.ghuud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghu(kp,:),dr.ghud{1}(kp,:)));
        Eud = eye(M_.exo_det_nbr);
        for i = 2:M_.exo_det_length
            hudi = dr.ghud{i}(kp,:);
            zudi=[zeros(np,M_.exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
            R2 = hessian*kron(zu,zudi);
            dr.ghuud{i} = -M2*dr.ghxud{i-1}(kf,:)*kron(hu,Eud)-M1*R2;
        end
        R1 = hessian*kron(zud,zud);
        dr.ghudud = cell(M_.exo_det_length,M_.exo_det_length);
        dr.ghudud{1,1} = -M1*R1-M2*dr.ghxx(kf,:)*kron(hud,hud);
        for i = 2:M_.exo_det_length
            hudi = dr.ghud{i}(nstatic+1:nstatic+npred,:);
            zudi=[zeros(np,M_.exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi+dr.ghud{i-1}(kf,:);zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
            R2 = hessian*kron(zudi,zudi);
            dr.ghudud{i,i} = -M2*(dr.ghudud{i-1,i-1}(kf,:)+...
                                  2*dr.ghxud{i-1}(kf,:)*kron(hudi,Eud) ...
                                  +dr.ghxx(kf,:)*kron(hudi,hudi))-M1*R2;
            R2 = hessian*kron(zud,zudi);
            dr.ghudud{1,i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(hud,Eud)+...
                                  dr.ghxx(kf,:)*kron(hud,hudi))...
                -M1*R2;
            for j=2:i-1
                hudj = dr.ghud{j}(kp,:);
                zudj=[zeros(np,M_.exo_det_nbr);dr.ghud{j};gx(:,1:npred)*hudj;zeros(M_.exo_nbr+M_.exo_det_nbr,M_.exo_det_nbr)];
                R2 = hessian*kron(zudj,zudi);
                dr.ghudud{j,i} = -M2*(dr.ghudud{j-1,i-1}(kf,:)+dr.ghxud{j-1}(kf,:)* ...
                                      kron(hudi,Eud)+dr.ghxud{i-1}(kf,:)* ...
                                      kron(hudj,Eud)+dr.ghxx(kf,:)*kron(hudj,hudi))-M1*R2;
            end
            
        end
    end