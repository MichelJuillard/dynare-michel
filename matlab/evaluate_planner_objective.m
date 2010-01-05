function oo1 = evaluate_planner_objective(dr,M,oo,options)
    
    oo1 = oo;
    endo_nbr = M.endo_nbr;
    exo_nbr = M.exo_nbr;
    nstatic = dr.nstatic;
    npred = dr.npred;
    lead_lag_incidence = M.lead_lag_incidence;
    if (size(lead_lag_incidence,1) == 5) ...
            && (nnz(lead_lag_incidence(1,:)) > 0 ...
                || nnz(lead_lag_incidence(5,:)) > 0)
        error(['procedure not yet implemented for leads and lags on more ' ...
               'than one period'])
    end
    if options.ramsey_policy
        i_org = (1:M.orig_model.endo_nbr)';
    else
        i_org = (1:M.endo_nbr)';
    end
    ipred = find(lead_lag_incidence(M.maximum_lag,:))';
    order_var = dr.order_var;
    ys = oo.dr.ys;
    [U,Uy,Uyy] = feval([M.fname '_objective_static'],ys,zeros(1,exo_nbr), ...
                       M.params);
    %K = reshape(1:endo_nbr^2,endo_nbr,endo_nbr);
    %K = K(order_var,order_var);
    %Uyy = Uyy(:,K(:));
    beta = options.planner_discount;
    gy(dr.order_var,:) = dr.ghx;
    gu(dr.order_var,:) = dr.ghu;
    Gy = gy(ipred,:);
    Gu = gu(ipred,:);
    gy = gy(i_org,:);
    gu = gu(i_org,:);
    
    Wbar = U/(1-beta);
    Wy = Uy*gy/(eye(npred)-beta*Gy);
    Wu = Uy*gu + beta*Wy*Gu;
    Wyy = Uyy*kron(gy,gy)/(eye(npred^2)-beta*kron(Gy,Gy));
    Wyu = Uyy*kron(gy,gu)+beta*Wyy*kron(Gy,Gu);
    Wuu = Uyy*kron(gu,gu)+beta*Wyy*kron(Gu,Gu);
    Wss = beta*Wuu*vec(M.Sigma_e)/(1-beta);
    
    if options.ramsey_policy
        yhat = [oo.endo_simul; ...
                zeros(M.endo_nbr-size(oo.endo_simul,1),1)];
    else
        yhat = oo.endo_simul;
    end
    yhat = yhat(dr.order_var(nstatic+(1:npred)),1)-dr.ys(dr.order_var(nstatic+(1:npred)));
    u = oo.exo_simul(1,:)';

    planner_objective_value = Wbar+Wy*yhat+Wu*u+Wyu*kron(yhat,u) ...
        + 0.5*(Wyy*kron(yhat,yhat) + Wuu*kron(u,u)+Wss);
    oo1.planner_objective_value = planner_objective_value;
    if ~options.noprint
        disp(' ')
        if options.ramsey_policy
            disp(['Value of planner objective function under Ramsey policy: ' ...
                  num2str(planner_objective_value)])
        else
            disp(['Value of planner objective function: ' ...
                  num2str(planner_objective_value)])
            disp(' ')
        end
    end