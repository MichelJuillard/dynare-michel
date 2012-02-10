function [flag,endo_simul,err] = solve_stochastic_perfect_foresight_model(endo_simul,exo_simul,pfm,nnodes)

    flag = 0;
    err = 0;
    stop = 0;

    number_of_shocks = size(exo_simul,2);

    [nodes,weights] = gauss_hermite_weights_and_nodes(nnodes);

    if number_of_shocks>1
        for i=1:number_of_shocks
            rr(i) = {nodes};
            ww(i) = {weights};
        end
        nodes = cartesian_product_of_sets(rr{:});
        weights = prod(cartesian_product_of_sets(ww{:}),2);
    end

    innovations = zeros(pfm.periods+2,number_of_shocks);

    model_dynamic = pfm.dynamic_model;

    dimension = (2+pfm.periods)*pfm.ny; % First n are given, dimension-n is the number of unknowns.

    Y = repmat(endo_simul(:),dimension/pfm.ny,1);

    if pfm.verbose
        disp ([' -----------------------------------------------------']);
        disp (['MODEL SIMULATION :']);
        fprintf('\n');
    end

    z = Y(find(pfm.lead_lag_incidence'));
    [d1,jacobian] = model_dynamic(z,exo_simul,pfm.params,pfm.steady_state,2);

    A = sparse([],[],[],dimension,dimension,dimension/pfm.ny*nnz(jacobian));
    res = zeros(dimension,1);

    h1 = clock;
    for iter = 1:pfm.maxit_
        h2 = clock;
        i_rows = 1:pfm.ny;
        i_cols = find(pfm.lead_lag_incidence');
        i_cols_A = i_cols;
        for it = 2:(pfm.periods+1)
            i_cols
            if it == 2
                y = Y(i_cols);
                expectations = zeros(pfm.nyf,1);
                tt = pfm.ny+pfm.ny;
                for n=1:nnodes
                    expectations = expectations+weights(n)*Y(tt+(n-1)*pfm.ny+pfm.iyf);
                end
                y(it*pfm.ny+pfm.iyf) = expectations;
                [d1,jacobian] = model_dynamic(y,exo_simul,pfm.params,pfm.steady_state,it);
                A(i_rows,pfm.i_cols_A1) = jacobian(:,pfm.i_cols_1);
            elseif it == pfm.periods+1
                A(i_rows,i_cols_A(pfm.i_cols_T)) = jacobian(:,pfm.i_cols_T);
            else
                for n=1:nnodes
                    innovations(3,:) = nodes(n,:);
                    [d1,jacobian] = model_dynamic(Y(i_cols),innovations,pfm.params,pfm.steady_state,it);
                    A(i_rows,i_cols_A) = jacobian(:,pfm.i_cols_j);
                end
                return
            end
            res(i_rows) = d1;
            i_rows = i_rows + pfm.ny;
            i_cols = i_cols + pfm.ny;
            if it > 2
                i_cols_A = i_cols_A + pfm.ny;
            end
        end
        err = max(abs(res));
        if err < pfm.tolerance
            stop = 1 ;
            if pfm.verbose
                fprintf('\n') ;
                disp([' Total time of simulation        :' num2str(etime(clock,h1))]) ;
                fprintf('\n') ;
                disp([' Convergency obtained.']) ;
                fprintf('\n') ;
            end
            flag = 0;% Convergency obtained.
            endo_simul = reshape(Y,pfm.ny,pfm.periods+2);
            break
        end
        dy = -A\res;
        Y(pfm.i_upd) =   Y(pfm.i_upd) + dy;
    end

    if ~stop
        if pfm.verbose
            fprintf('\n') ;
            disp(['     Total time of simulation        :' num2str(etime(clock,h1))]) ;
            fprintf('\n') ;
            disp(['WARNING : maximum number of iterations is reached (modify options_.maxit_).']) ;
            fprintf('\n') ;
        end
        flag = 1;% more iterations are needed.
        endo_simul = 1;
    end
    if pfm.verbose
        disp (['-----------------------------------------------------']) ;
    end