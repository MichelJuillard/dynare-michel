function [dr,info] = dyn_risky_steadystate_solver(ys0,M, ...
                                                  dr,options,oo)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info}] =} dyn_risky_steadystate_solver (@var{ys0},@var{M},@var{dr},@var{options},@var{oo})
%! @anchor{dyn_risky_steadystate_solver}
%! @sp 1
%! Computes the second order risky steady state and first and second order reduced form of the DSGE model.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ys0
%! Vector containing a guess value for the risky steady state
%! @item M
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item options
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @item oo
%! Matlab's structure gathering the results (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @end table
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2012 Dynare Team
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
    lead_lag_incidence = M.lead_lag_incidence;
    order_var = dr.order_var;
    exo_nbr = M.exo_nbr;
    
    M.var_order_endo_names = M.endo_names(dr.order_var,:);
 
    [junk,dr.i_fwrd_g,i_fwrd_f] = find(lead_lag_incidence(3,order_var));
    dr.i_fwrd_f = i_fwrd_f;
    nd = nnz(lead_lag_incidence) + M.exo_nbr;
    dr.nd = nd;
    kk = reshape(1:nd^2,nd,nd);
    kkk = reshape(1:nd^3,nd^2,nd);
    dr.i_fwrd2_f = kk(i_fwrd_f,i_fwrd_f);
    dr.i_fwrd2a_f = kk(i_fwrd_f,:);
    dr.i_fwrd3_f = kkk(dr.i_fwrd2_f,:);
    dr.i_uu = kk(end-exo_nbr+1:end,end-exo_nbr+1:end);
    if options.k_order_solver
        func = @risky_residuals_k_order;
    else
        func = @risky_residuals;
    end
    
    if isfield(options,'portfolio') && options.portfolio == 1
        eq_tags = M.equations_tags;
        n_tags = size(eq_tags,1);
        l_var = zeros(n_tags,1);
        for i=1:n_tags
            l_var(i) = find(strncmp(eq_tags(i,3),M.endo_names, ...
                                    length(cell2mat(eq_tags(i,3)))));
        end
        dr.ys = ys0;
        x0 = ys0(l_var);
        %        dr = first_step_ds(x0,M,dr,options,oo);
        n = size(ys0);
        %x0 = ys0;
        [x, info] = solve1(@risky_residuals_ds,x0,1:n_tags,1:n_tags,0,1,M,dr,options,oo);
        %[x, info] = solve1(@risky_residuals,x0,1:n,1:n,0,1,M,dr,options,oo);
        %        ys0(l_var) = x;
        ys0(l_var) = x;
        dr.ys = ys0;
        oo.dr = dr;
        oo.steady_state = ys0;
        disp_steady_state(M,oo);
    end
        
    [ys, info] = csolve(func,ys0,[],1e-10,100,M,dr,options,oo);
    
    if options.k_order_solver
        [resid,dr] = risky_residuals_k_order(ys,M,dr,options,oo);
    else
        [resid,dr] = risky_residuals(ys,M,dr,options,oo);
    end
    
    dr.ys = ys;
    for i=1:M.endo_nbr
        disp(sprintf('%16s %12.6f %12.6f',M.endo_names(i,:),ys0(i), ys(i)))
    end
    
    dr.ghs2 = zeros(size(dr.ghs2));

    k_var = setdiff(1:M.endo_nbr,l_var);
    dr.ghx(k_var,:) = dr.ghx;
    dr.ghu(k_var,:) = dr.ghu;
    dr.ghxx(k_var,:) = dr.ghxx;
    dr.ghxu(k_var,:) = dr.ghxu;
    dr.ghuu(k_var,:) = dr.ghuu;
    dr.ghs2(k_var,:) = dr.ghs2;
end

function [resid,dr] = risky_residuals(ys,M,dr,options,oo)
    persistent old_ys old_resid
    
    lead_lag_incidence = M.lead_lag_incidence;
    iyv = lead_lag_incidence';
    iyv = iyv(:);
    iyr0 = find(iyv) ;
    
    if M.exo_nbr == 0
        oo.exo_steady_state = [] ;
    end
    
    z = repmat(ys,1,3);
    z = z(iyr0) ;
    [resid1,d1,d2] = feval([M.fname '_dynamic'],z,...
                                     [oo.exo_simul ...
                        oo.exo_det_simul], M.params, dr.ys, 2);
    if ~isreal(d1) || ~isreal(d2)
        pause
    end
    
    if options.use_dll
        % In USE_DLL mode, the hessian is in the 3-column sparse representation
        d2 = sparse(d2(:,1), d2(:,2), d2(:,3), ...
                         size(d1, 1), size(d1, 2)*size(d1, 2));
    end

    if isfield(options,'portfolio') && options.portfolio == 1
        eq_tags = M.equations_tags;
        n_tags = size(eq_tags,1);
        portfolios_eq = cell2mat(eq_tags(strcmp(eq_tags(:,2), ...
                                                'portfolio'),1));
        eq = setdiff(1:M.endo_nbr,portfolios_eq);
        l_var = zeros(n_tags,1);
        for i=1:n_tags
            l_var(i) = find(strncmp(eq_tags(i,3),M.endo_names, ...
                                    length(cell2mat(eq_tags(i,3)))));
        end
        k_var = setdiff(1:M.endo_nbr,l_var);
        lli1 = lead_lag_incidence(:,k_var);
        lead_incidence = lli1(3,:)';
        k = find(lli1');
        lli2 = lli1';
        lli2(k) = 1:nnz(lli1);
        lead_lag_incidence = lli2';
        x = ys(l_var);
        dr = first_step_ds(x,M,dr,options,oo);

        
        M.lead_lag_incidence = lead_lag_incidence;
        lli1a = [nonzeros(lli1'); size(d1,2)+(-M.exo_nbr+1:0)'];
        d1a = d1(eq,lli1a);
        ih = 1:size(d2,2);
        ih = reshape(ih,size(d1,2),size(d1,2));
        ih1 = ih(lli1a,lli1a);
        d2a = d2(eq,ih1);
        
        M.endo_nbr = M.endo_nbr-n_tags;
        dr = set_state_space(dr,M,options);
    
        [junk,dr.i_fwrd_g] = find(lead_lag_incidence(3,dr.order_var));
        i_fwrd_f = nonzeros(lead_incidence(dr.order_var));
        i_fwrd2_f = ih(i_fwrd_f,i_fwrd_f);
        dr.i_fwrd_f = i_fwrd_f;
        dr.i_fwrd2_f = i_fwrd2_f;
        nd = nnz(lead_lag_incidence) + M.exo_nbr;
        dr.nd = nd;
        kk = reshape(1:nd^2,nd,nd);
        kkk = reshape(1:nd^3,nd^2,nd);
        dr.i_fwrd2a_f = kk(i_fwrd_f,:);
        %        dr.i_fwrd3_f = kkk(i_fwrd2_f,:);
        dr.i_uu = kk(end-M.exo_nbr+1:end,end-M.exo_nbr+1:end);
    else
        d1a = d1;
        d2a = d2;
    end
    
% $$$     [junk,cols_b,cols_j] = find(lead_lag_incidence(2,dr.order_var));
% $$$     b = zeros(M.endo_nbr,M.endo_nbr);
% $$$     b(:,cols_b) = d1a(:,cols_j);
% $$$ 
% $$$     [dr,info] = dyn_first_order_solver(d1a,b,M,dr,options,0);
% $$$     if info
% $$$         [m1,m2]=max(abs(ys-old_ys));
% $$$         disp([m1 m2])
% $$$         %        print_info(info,options.noprint);
% $$$         resid = old_resid+info(2)/40;
% $$$         return
% $$$     end
% $$$     
% $$$     dr = dyn_second_order_solver(d1a,d2a,dr,M);
    
    gu1 = dr.ghu(dr.i_fwrd_g,:);

    resid = resid1+0.5*(d1(:,dr.i_fwrd_f)*dr.ghuu(dr.i_fwrd_g,:)+ ...
                        d2(:,dr.i_fwrd2_f)*kron(gu1,gu1))*vec(M.Sigma_e);
    disp(d1(:,dr.i_fwrd_f)*dr.ghuu(dr.i_fwrd_g,:)*vec(M.Sigma_e));
    old_ys = ys;
    disp(max(abs(resid)))
    old_resid = resid;
end

function [resid,dr] = risky_residuals_ds(x,M,dr,options,oo)
    persistent old_ys old_resid old_resid1 old_d1 old_d2
    
    dr = first_step_ds(x,M,dr,options,oo);

    lead_lag_incidence = M.lead_lag_incidence;
    iyv = lead_lag_incidence';
    iyv = iyv(:);
    iyr0 = find(iyv) ;
    
    if M.exo_nbr == 0
        oo.exo_steady_state = [] ;
    end
    
    eq_tags = M.equations_tags;
    n_tags = size(eq_tags,1);
    portfolios_eq = cell2mat(eq_tags(strcmp(eq_tags(:,2), ...
                                            'portfolio'),1));
    eq = setdiff(1:M.endo_nbr,portfolios_eq);
    l_var = zeros(n_tags,1);
    for i=1:n_tags
        l_var(i) = find(strncmp(eq_tags(i,3),M.endo_names, ...
                                length(cell2mat(eq_tags(i,3)))));
    end
    k_var = setdiff(1:M.endo_nbr,l_var);
    lli1 = lead_lag_incidence(:,k_var);
    k = find(lli1');
    lli2 = lli1';
    lli2(k) = 1:nnz(lli1);
    lead_lag_incidence = lli2';

    ys = dr.ys;
    ys(l_var) = x;
    
    z = repmat(ys,1,3);
    z = z(iyr0) ;
    [resid1,d1,d2] = feval([M.fname '_dynamic'],z,...
                                     [oo.exo_simul ...
                        oo.exo_det_simul], M.params, dr.ys, 2);
% $$$     if isempty(old_resid)
% $$$         old_resid1 = resid1;
% $$$         old_d1 = d1;
% $$$         old_d2 = d2;
% $$$         old_ys = ys;
% $$$     else
% $$$         if ~isequal(resid1,old_resid)
% $$$             disp('ys')
% $$$             disp((ys-old_ys)');
% $$$             disp('resids1')
% $$$             disp((resid1-old_resid1)')
% $$$             old_resid1 = resid1;
% $$$             pause
% $$$         end
% $$$         if ~isequal(d1,old_d1)
% $$$             disp('d1')
% $$$             disp(d1-old_d1);
% $$$             old_d1 = d1;
% $$$             pause
% $$$         end
% $$$         if ~isequal(d2,old_d2)
% $$$             disp('d2')
% $$$             disp(d2-old_d2);
% $$$             old_d2 = d2;
% $$$             pause
% $$$         end
% $$$     end
    if ~isreal(d1) || ~isreal(d2)
        pause
    end
    
    if options.use_dll
        % In USE_DLL mode, the hessian is in the 3-column sparse representation
        d2 = sparse(d2(:,1), d2(:,2), d2(:,3), ...
                         size(d1, 1), size(d1, 2)*size(d1, 2));
    end

% $$$     if isfield(options,'portfolio') && options.portfolio == 1
% $$$         lli1a = [nonzeros(lli1'); size(d1,2)+(-M.exo_nbr+1:0)'];
% $$$         d1a = d1(eq,lli1a);
% $$$         ih = 1:size(d2,2);
% $$$         ih = reshape(ih,size(d1,2),size(d1,2));
% $$$         ih1 = ih(lli1a,lli1a);
% $$$         d2a = d2(eq,ih1);
% $$$         
% $$$         M.endo_nbr = M.endo_nbr-n_tags;
% $$$         dr = set_state_space(dr,M);
% $$$     
% $$$         dr.i_fwrd_g = find(lead_lag_incidence(3,dr.order_var)');
% $$$     else
% $$$         d1a = d1;
% $$$         d2a = d2;
% $$$     end
% $$$     
% $$$     [junk,cols_b,cols_j] = find(lead_lag_incidence(2,dr.order_var));
% $$$     b = zeros(M.endo_nbr,M.endo_nbr);
% $$$     b(:,cols_b) = d1a(:,cols_j);
% $$$ 
% $$$     [dr,info] = dyn_first_order_solver(d1a,b,M,dr,options,0);
% $$$     if info
% $$$         [m1,m2]=max(abs(ys-old_ys));
% $$$         disp([m1 m2])
% $$$         %        print_info(info,options.noprint);
% $$$         resid = old_resid+info(2)/40;
% $$$         return
% $$$     end
% $$$     
% $$$     dr = dyn_second_order_solver(d1a,d2a,dr,M);
    
    gu1 = dr.ghu(dr.i_fwrd_g,:);

    %    resid = resid1+0.5*(d1(:,dr.i_fwrd_f)*dr.ghuu(dr.i_fwrd_g,:)+ ...
    %                    d2(:,dr.i_fwrd2_f)*kron(gu1,gu1))*vec(M.Sigma_e);
    resid = resid1+0.5*(d2(:,dr.i_fwrd2_f)*kron(gu1,gu1))*vec(M.Sigma_e);

% $$$     if isempty(old_resid)
% $$$         old_resid = resid;
% $$$     else
% $$$         disp('resid')
% $$$         dr = (resid-old_resid)';
% $$$         %        disp(dr)
% $$$         %        disp(dr(portfolios_eq))
% $$$         old_resid = resid;
% $$$     end
    resid = resid(portfolios_eq)
end

function [dr] = first_step_ds(x,M,dr,options,oo)
    
    lead_lag_incidence = M.lead_lag_incidence;
    iyv = lead_lag_incidence';
    iyv = iyv(:);
    iyr0 = find(iyv) ;
    
    if M.exo_nbr == 0
        oo.exo_steady_state = [] ;
    end
    
    eq_tags = M.equations_tags;
    n_tags = size(eq_tags,1);
    portfolios_eq = cell2mat(eq_tags(strcmp(eq_tags(:,2), ...
                                            'portfolio'),1));
    eq = setdiff(1:M.endo_nbr,portfolios_eq);
    l_var = zeros(n_tags,1);
    for i=1:n_tags
        l_var(i) = find(strncmp(eq_tags(i,3),M.endo_names, ...
                                length(cell2mat(eq_tags(i,3)))));
    end
    k_var = setdiff(1:M.endo_nbr,l_var);
    lli1 = lead_lag_incidence(:,k_var);
    k = find(lli1');
    lli2 = lli1';
    lli2(k) = 1:nnz(lli1);
    lead_lag_incidence = lli2';
    M.lead_lag_incidence = lead_lag_incidence;

    ys = dr.ys;
    ys(l_var) = x;
    
    z = repmat(ys,1,3);
    z = z(iyr0) ;
    [resid1,d1,d2] = feval([M.fname '_dynamic'],z,...
                                     [oo.exo_simul ...
                        oo.exo_det_simul], M.params, dr.ys, 2);
    if ~isreal(d1) || ~isreal(d2)
        pause
    end
    
    if options.use_dll
        % In USE_DLL mode, the hessian is in the 3-column sparse representation
        d2 = sparse(d2(:,1), d2(:,2), d2(:,3), ...
                         size(d1, 1), size(d1, 2)*size(d1, 2));
    end

    if isfield(options,'portfolio') && options.portfolio == 1
        lli1a = [nonzeros(lli1'); size(d1,2)+(-M.exo_nbr+1:0)'];
        d1a = d1(eq,lli1a);
        ih = 1:size(d2,2);
        ih = reshape(ih,size(d1,2),size(d1,2));
        ih1 = ih(lli1a,lli1a);
        d2a = d2(eq,ih1);
        
        M.endo_nbr = M.endo_nbr-n_tags;
        dr = set_state_space(dr,M,options);
    
        dr.i_fwrd_g = find(lead_lag_incidence(3,dr.order_var)');
    else
        d1a = d1;
        d2a = d2;
    end
    
    [junk,cols_b,cols_j] = find(lead_lag_incidence(2,dr.order_var));
    b = zeros(M.endo_nbr,M.endo_nbr);
    b(:,cols_b) = d1a(:,cols_j);

    [dr,info] = dyn_first_order_solver(d1a,M,dr,options,0);
    if info
        [m1,m2]=max(abs(ys-old_ys));
        disp([m1 m2])
        %        print_info(info,options.noprint);
        resid = old_resid+info(2)/40;
        return
    end
    
    dr = dyn_second_order_solver(d1a,d2a,dr,M,...
                                 options.threads.kronecker.A_times_B_kronecker_C,...
                                 options.threads.kronecker.sparse_hessian_times_B_kronecker_C);
end

function [resid,dr] = risky_residuals_k_order(ys,M,dr,options,oo)
    
    lead_lag_incidence = M.lead_lag_incidence;
    npred = dr.npred;
    exo_nbr = M.exo_nbr;
    vSigma_e = vec(M.Sigma_e);
    
    iyv = lead_lag_incidence';
    iyv = iyv(:);
    iyr0 = find(iyv) ;
    
    if M.exo_nbr == 0
        oo.exo_steady_state = [] ;
    end
    
    z = repmat(ys,1,3);
    z = z(iyr0) ;
    [resid1,d1,d2,d3] = feval([M.fname '_dynamic'],z,...
                                     [oo.exo_simul ...
                        oo.exo_det_simul], M.params, dr.ys, 2);

% $$$     hessian = sparse(d2(:,1), d2(:,2), d2(:,3), ...
% $$$                      size(d1, 1), size(d1, 2)*size(d1, 2));
% $$$     fy3 = sparse(d2(:,1), d2(:,2), d2(:,3), ...
% $$$                      size(d1, 1), size(d1, 2)^3);

    options.order = 3;
    
    nu2 = exo_nbr*(exo_nbr+1)/2;
% $$$     d1_0 = d1;
% $$$     gu1 = dr.ghu(dr.i_fwrd_g,:);
% $$$     guu = dr.ghuu;
% $$$     for i=1:2
% $$$         d1 = d1_0 + 0.5*(hessian(:,dr.i_fwrd2a_f)*kron(eye(dr.nd),guu(dr.i_fwrd_g,:)*vSigma_e)+ ...
% $$$                        fy3(:,dr.i_fwrd3_f)*kron(eye(dr.nd),kron(gu1,gu1)*vSigma_e));
% $$$     [junk,cols_b,cols_j] = find(lead_lag_incidence(2,dr.order_var));
% $$$     b = zeros(M.endo_nbr,M.endo_nbr);
% $$$     b(:,cols_b) = d1(:,cols_j);

% $$$     [dr,info] = dyn_first_order_solver(d1,b,M,dr,options,0);
        [err,g_0, g_1, g_2, g_3] = k_order_perturbation(dr,M,options);
        mexErrCheck('k_order_perturbation', err);
        gu1 = g_1(dr.i_fwrd_g,end-exo_nbr+1:end);
        guu = unfold(g_2(:,end-nu2+1:end),exo_nbr);
        d1old = d1;
        %        disp(max(max(abs(d1-d1old))));
        %    end
    
    [junk,cols_b,cols_j] = find(lead_lag_incidence(2,dr.order_var));
    
    resid = resid1+0.5*(d1(:,dr.i_fwrd_f)*guu(dr.i_fwrd_g,:)+hessian(:,dr.i_fwrd2_f)*kron(gu1,gu1))*vec(M.Sigma_e);

    if nargout > 1
        [dr,info] = k_order_pert(dr,M,options,oo);
    end
end

function y=unfold(x,n)
    y = zeros(size(x,1),n*n);
    k = 1;
    for i=1:n
        for j=i:n
            y(:,(i-1)*n+j) = x(:,k);
            if i ~= j
                y(:,(j-1)*n+i) = x(:,k);
            end
        end
    end
end
